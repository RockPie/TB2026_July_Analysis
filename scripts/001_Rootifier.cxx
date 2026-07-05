#include <spdlog/spdlog.h>

#include <TFile.h>
#include <TNamed.h>

#include "cxxopts.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

constexpr std::size_t kRowSize = 32;
constexpr std::size_t kWordSize = 4;
constexpr std::size_t kFrameWordOffset = 8;
constexpr std::size_t kFrameWordCount = 6;
constexpr std::size_t kTimestampHexWidth = 12;

constexpr std::array<std::uint8_t, kWordSize> kZeroWord{0x00, 0x00, 0x00, 0x00};
constexpr std::array<std::uint8_t, kWordSize> kFrameIdleWord{0xcc, 0xcc, 0xcc, 0xac};
constexpr std::array<std::uint8_t, kWordSize> kFmcIdleWord{0x36, 0x36, 0x36, 0x36};
constexpr std::array<std::uint8_t, kWordSize> kFmcInjectionWord{0x2d, 0x2d, 0x2d, 0x2d};
constexpr std::array<std::uint8_t, kWordSize> kFmcL1aWord{0x4b, 0x4b, 0x4b, 0x4b};

enum class LineType {
    heartbeat,
    trigger,
    data_idle,
    data_daq,
    dummy,
    other,
    partial_tail,
};

enum class DisplayMode {
    raw,
    non_dummy_non_heartbeat,
    data_idle,
    data_daq,
    heartbeat,
    trigger,
    dummy,
    other,
};

struct LineStats {
    std::size_t heartbeat = 0;
    std::size_t trigger = 0;
    std::size_t data_idle = 0;
    std::size_t data_daq = 0;
    std::size_t dummy = 0;
    std::size_t other = 0;
    std::size_t partial_tail = 0;
};

struct ParseSummary {
    std::string input_path;
    std::uintmax_t file_size = 0;
    std::size_t total_lines = 0;
    std::size_t total_bytes_parsed = 0;
    std::size_t partial_tail_bytes = 0;
    std::size_t printed_lines = 0;
    double elapsed_sec = 0.0;
    LineStats stats;
};

using RawLine = std::array<std::uint8_t, kRowSize>;
using Word = std::array<std::uint8_t, kWordSize>;
using LineRange = std::pair<std::size_t, std::size_t>;

struct PrintableLine {
    std::size_t line_index = 0;
    RawLine row{};
    LineType type = LineType::other;
    std::optional<std::uint64_t> last_trigger_timestamp;
    std::optional<RawLine> hb_tail;
};

bool data_lanes_all_match(const RawLine& row, const Word& expected);

struct DisplayOptions {
    std::vector<LineRange> selected_ranges;
    DisplayMode mode = DisplayMode::non_dummy_non_heartbeat;
    std::size_t limit = 40;
    bool use_selected_ranges = false;
};

std::string_view line_type_name(LineType type)
{
    switch (type) {
    case LineType::heartbeat:
        return "heartbeat line";
    case LineType::trigger:
        return "phy trigger line";
    case LineType::data_idle:
        return "data idle line";
    case LineType::data_daq:
        return "data daq line";
    case LineType::dummy:
        return "dummy zero line";
    case LineType::other:
        return "other line";
    case LineType::partial_tail:
        return "partial tail";
    }
    return "unknown";
}

std::string_view line_type_code(LineType type)
{
    switch (type) {
    case LineType::heartbeat:
        return "[HB]";
    case LineType::trigger:
        return "[PT]";
    case LineType::data_idle:
        return "[AC]";
    case LineType::data_daq:
        return "[DQ]";
    case LineType::dummy:
        return "[00]";
    case LineType::other:
        return "[XX]";
    case LineType::partial_tail:
        return "[tail]";
    }
    return "[??]";
}

std::string color_line_type_code(LineType type)
{
    constexpr std::string_view reset = "\033[0m";
    std::string_view color;
    switch (type) {
    case LineType::heartbeat:
        color = "\033[1;36m";
        break;
    case LineType::trigger:
        color = "\033[1;33m";
        break;
    case LineType::data_idle:
        color = "\033[90m";
        break;
    case LineType::data_daq:
        color = "\033[1;32m";
        break;
    case LineType::dummy:
        color = "\033[1;34m";
        break;
    case LineType::other:
        color = "\033[1;31m";
        break;
    case LineType::partial_tail:
        color = "\033[1;35m";
        break;
    }
    return std::string(color) + std::string(line_type_code(type)) + std::string(reset);
}

std::string color_section(std::string_view text, std::string_view color)
{
    constexpr std::string_view reset = "\033[0m";
    return std::string(color) + std::string(text) + std::string(reset);
}

std::string_view gbt_color(std::uint8_t gbt_id)
{
    constexpr std::array<std::string_view, 6> colors{
        "\033[1;36m",
        "\033[1;32m",
        "\033[1;33m",
        "\033[1;35m",
        "\033[1;34m",
        "\033[1;31m",
    };
    return colors[gbt_id % colors.size()];
}

std::string_view fcmd_color(std::string_view fcmd)
{
    if (fcmd == "idle") {
        return "\033[90m";
    }
    if (fcmd == "L1A") {
        return "\033[1;33m";
    }
    if (fcmd == "inj") {
        return "\033[1;35m";
    }
    if (fcmd == "zero") {
        return "\033[1;34m";
    }
    return "\033[1;31m";
}

std::string_view word_color(std::size_t word_index)
{
    constexpr std::array<std::string_view, 4> colors{
        "\033[1;36m",
        "\033[1;32m",
        "\033[1;33m",
        "\033[1;35m",
    };
    return colors[word_index % colors.size()];
}

std::string_view version_color(std::uint8_t version)
{
    constexpr std::array<std::string_view, 6> colors{
        "\033[1;36m",
        "\033[1;32m",
        "\033[1;33m",
        "\033[1;35m",
        "\033[1;34m",
        "\033[1;31m",
    };
    return colors[version % colors.size()];
}

void increment_stats(LineStats& stats, LineType type)
{
    switch (type) {
    case LineType::heartbeat:
        ++stats.heartbeat;
        break;
    case LineType::trigger:
        ++stats.trigger;
        break;
    case LineType::data_idle:
        ++stats.data_idle;
        break;
    case LineType::data_daq:
        ++stats.data_daq;
        break;
    case LineType::dummy:
        ++stats.dummy;
        break;
    case LineType::other:
        ++stats.other;
        break;
    case LineType::partial_tail:
        ++stats.partial_tail;
        break;
    }
}

bool is_zero_row(const RawLine& row)
{
    return std::all_of(row.begin(), row.end(), [](std::uint8_t value) { return value == 0x00; });
}

LineType classify_line(const RawLine& row)
{
    if (is_zero_row(row)) {
        return LineType::dummy;
    }
    if (row[0] == 0x07) {
        return LineType::heartbeat;
    }
    if (row[0] == 0xbb) {
        return LineType::trigger;
    }
    if (row[0] == 0xac) {
        return data_lanes_all_match(row, kFrameIdleWord) ? LineType::data_idle : LineType::data_daq;
    }
    return LineType::other;
}

std::uint64_t read_u48_le(const RawLine& row, std::size_t offset)
{
    std::uint64_t value = 0;
    for (std::size_t index = 0; index < 6; ++index) {
        value |= static_cast<std::uint64_t>(row[offset + index]) << (8 * index);
    }
    return value;
}

std::uint16_t read_u16_le(const std::array<std::uint8_t, 64>& bytes, std::size_t offset)
{
    return static_cast<std::uint16_t>(bytes[offset]) |
        static_cast<std::uint16_t>(bytes[offset + 1]) << 8;
}

std::uint32_t read_u32_le(const std::array<std::uint8_t, 64>& bytes, std::size_t offset)
{
    return static_cast<std::uint32_t>(bytes[offset]) |
        (static_cast<std::uint32_t>(bytes[offset + 1]) << 8) |
        (static_cast<std::uint32_t>(bytes[offset + 2]) << 16) |
        (static_cast<std::uint32_t>(bytes[offset + 3]) << 24);
}

std::uint32_t extract_trigger_timestamp(const RawLine& row)
{
    return static_cast<std::uint32_t>(row[4]) |
        (static_cast<std::uint32_t>(row[5]) << 8) |
        (static_cast<std::uint32_t>(row[12]) << 16) |
        (static_cast<std::uint32_t>(row[13]) << 24);
}

std::uint64_t extract_stream_timestamp(const RawLine& row)
{
    return read_u48_le(row, 2);
}

Word frame_word(const RawLine& row, std::size_t word_index)
{
    Word word{};
    const std::size_t offset = kFrameWordOffset + word_index * kWordSize;
    std::copy_n(row.begin() + static_cast<std::ptrdiff_t>(offset), kWordSize, word.begin());
    return word;
}

std::uint32_t reordered_word_value(const Word& word)
{
    return static_cast<std::uint32_t>(word[3]) |
        (static_cast<std::uint32_t>(word[2]) << 8) |
        (static_cast<std::uint32_t>(word[1]) << 16) |
        (static_cast<std::uint32_t>(word[0]) << 24);
}

bool data_lanes_all_match(const RawLine& row, const Word& expected)
{
    for (std::size_t word_index = 0; word_index < 4; ++word_index) {
        if (frame_word(row, word_index) != expected) {
            return false;
        }
    }
    return true;
}

bool data_lanes_have_activity(const RawLine& row)
{
    for (std::size_t word_index = 0; word_index < 4; ++word_index) {
        const auto word = frame_word(row, word_index);
        if (word != kZeroWord && word != kFrameIdleWord) {
            return true;
        }
    }
    return false;
}

std::string format_count(std::size_t value)
{
    std::string digits = std::to_string(value);
    for (int insert_position = static_cast<int>(digits.size()) - 3; insert_position > 0; insert_position -= 3) {
        digits.insert(static_cast<std::size_t>(insert_position), ",");
    }
    return digits;
}

std::string format_byte_size(std::uintmax_t bytes)
{
    constexpr std::array<std::string_view, 6> units{"B", "KB", "MB", "GB", "TB", "PB"};
    double value = static_cast<double>(bytes);
    std::size_t unit_index = 0;
    while (value >= 1000.0 && unit_index + 1 < units.size()) {
        value /= 1000.0;
        ++unit_index;
    }

    std::ostringstream out;
    if (unit_index == 0) {
        out << format_count(static_cast<std::size_t>(bytes)) << ' ' << units[unit_index];
    } else {
        out << std::fixed << std::setprecision(3) << value << ' ' << units[unit_index];
    }
    return out.str();
}

std::string format_byte_rate(std::uintmax_t bytes, double elapsed_sec)
{
    if (elapsed_sec <= 0.0) {
        return "0 B/s";
    }
    return format_byte_size(static_cast<std::uintmax_t>(static_cast<double>(bytes) / elapsed_sec)) + "/s";
}

std::string format_seconds(double seconds)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << seconds << " s";
    return out.str();
}

std::string file_name_only(const std::string& path)
{
    const auto separator = path.find_last_of("/\\");
    if (separator == std::string::npos) {
        return path;
    }
    return path.substr(separator + 1);
}

void print_box(std::string_view title, const std::vector<std::pair<std::string, std::string>>& rows)
{
    std::size_t label_width = title.size();
    std::size_t value_width = 0;
    for (const auto& [label, value] : rows) {
        label_width = std::max(label_width, label.size());
        value_width = std::max(value_width, value.size());
    }

    const std::size_t inner_width = label_width + value_width + 5;
    spdlog::info("+{}+", std::string(inner_width + 2, '-'));
    spdlog::info("| {:<{}} |", title, inner_width);
    spdlog::info("+{}+", std::string(inner_width + 2, '-'));
    for (const auto& [label, value] : rows) {
        spdlog::info("|   {:<{}} : {:>{}} |", label, label_width, value, value_width);
    }
    spdlog::info("+{}+", std::string(inner_width + 2, '-'));
}

std::string format_hex(std::uint64_t value, std::size_t width)
{
    std::ostringstream out;
    out << "0x" << std::uppercase << std::hex << std::setw(static_cast<int>(width)) << std::setfill('0') << value;
    return out.str();
}

std::string format_compact_hex(std::uint64_t value)
{
    std::ostringstream out;
    out << "0x" << std::uppercase << std::hex << value;
    return out.str();
}

std::string format_delta(std::optional<std::int64_t> value)
{
    if (!value.has_value()) {
        return "     -";
    }
    const auto delta = *value;
    const auto magnitude = delta < 0
        ? static_cast<std::uint64_t>(-(delta + 1)) + 1
        : static_cast<std::uint64_t>(delta);
    std::ostringstream out;
    if (delta < 0) {
        out << '-';
    }
    out << "0x" << std::uppercase << std::hex << std::setw(4) << std::setfill('0') << magnitude;
    return out.str();
}

std::string format_bytes(const std::uint8_t* bytes, std::size_t size)
{
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (std::size_t index = 0; index < size; ++index) {
        if (index != 0) {
            out << ' ';
        }
        out << std::setw(2) << static_cast<int>(bytes[index]);
    }
    return out.str();
}

std::string format_data_bytes(const RawLine& row)
{
    constexpr std::string_view reset = "\033[0m";
    constexpr std::string_view white = "\033[1;37m";
    constexpr std::string_view gray = "\033[90m";

    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (std::size_t index = 0; index < row.size(); ++index) {
        if (index != 0) {
            out << ' ';
        }
        if (index % kWordSize == 0) {
            if (index != 0) {
                out << reset;
            }
            out << (((index / kWordSize) % 2 == 0) ? white : gray);
        }
        out << std::setw(2) << static_cast<int>(row[index]);
    }
    out << reset;
    return out.str();
}

std::string format_word_bytes(const Word& word)
{
    return format_bytes(word.data(), word.size());
}

std::string describe_fmc_word(const Word& fmc_word)
{
    if (fmc_word == kFmcInjectionWord) {
        return "inj";
    }
    if (fmc_word == kFmcL1aWord) {
        return "L1A";
    }
    if (fmc_word == kFmcIdleWord) {
        return "idle";
    }
    if (fmc_word == kZeroWord) {
        return "zero";
    }
    return "other";
}

std::string describe_data_lanes(const RawLine& row)
{
    if (data_lanes_all_match(row, kZeroWord)) {
        return "all0";
    }
    if (data_lanes_all_match(row, kFrameIdleWord)) {
        return "allAC";
    }
    if (!data_lanes_have_activity(row)) {
        return "idle/mixed";
    }
    return "active";
}

std::string describe_frame_words(const RawLine& row)
{
    std::ostringstream out;
    for (std::size_t word_index = 0; word_index < kFrameWordCount; ++word_index) {
        if (word_index != 0) {
            out << " | ";
        }
        const auto word = frame_word(row, word_index);
        out << "w" << word_index << '=' << format_word_bytes(word)
            << "/" << format_hex(reordered_word_value(word), 8);
    }
    return out.str();
}

std::string format_byte_hex(std::uint8_t value)
{
    return format_compact_hex(value);
}

std::string describe_heartbeat(const RawLine& row, const std::optional<RawLine>& hb_tail)
{
    if (!hb_tail.has_value()) {
        return "heartbeat [HB] line; missing continuation row";
    }

    std::array<std::uint8_t, 64> bytes{};
    std::copy(row.begin(), row.end(), bytes.begin());
    std::copy(hb_tail->begin(), hb_tail->end(), bytes.begin() + static_cast<std::ptrdiff_t>(row.size()));

    std::ostringstream out;
    out << color_section("HB", "\033[1;36m")
        << ' ' << color_section("v" + format_byte_hex(bytes[0]), version_color(bytes[0]))
        << ' ' << color_section("hdrSize", "\033[1;36m") << ' ' << format_byte_hex(bytes[1])
        << ' ' << color_section("feeId", "\033[1;32m") << ' ' << format_compact_hex(read_u16_le(bytes, 2))
        << ' ' << color_section("priority", "\033[90m") << ' ' << format_byte_hex(bytes[4])
        << ' ' << color_section("sourceID", "\033[1;33m") << ' ' << format_byte_hex(bytes[5])
        << ' ' << color_section("offset", "\033[1;35m") << ' ' << format_compact_hex(read_u16_le(bytes, 8))
        << ' ' << color_section("memSize", "\033[1;35m") << ' ' << format_compact_hex(read_u16_le(bytes, 10))
        << ' ' << color_section("linkID", "\033[1;34m") << ' ' << format_byte_hex(bytes[12])
        << ' ' << color_section("pktCnt", "\033[1;34m") << ' ' << format_byte_hex(bytes[13])
        << ' ' << color_section("cru/ep", "\033[1;34m") << ' ' << format_compact_hex(read_u16_le(bytes, 14)) << '\n'
        << color_section("BC", "\033[1;31m") << ' ' << format_compact_hex(read_u16_le(bytes, 16))
        << ' ' << color_section("orbit", "\033[1;31m") << ' ' << format_compact_hex(read_u32_le(bytes, 20))
        << ' ' << color_section("dataFmt", "\033[90m") << ' ' << format_byte_hex(bytes[24])
        << ' ' << color_section("triggerType", "\033[1;33m") << ' ' << format_compact_hex(read_u32_le(bytes, 32))
        << ' ' << color_section("pageCnt", "\033[1;33m") << ' ' << format_compact_hex(read_u16_le(bytes, 36))
        << ' ' << color_section("stop", "\033[1;33m") << ' ' << format_byte_hex(bytes[38])
        << ' ' << color_section("detectorField", "\033[1;32m") << ' ' << format_compact_hex(read_u32_le(bytes, 48))
        << ' ' << color_section("detectorPAR", "\033[1;32m") << ' ' << format_compact_hex(read_u16_le(bytes, 52));
    return out.str();
}

std::string describe_line(const RawLine& row, LineType type, std::optional<std::uint64_t> last_trigger_timestamp, const std::optional<RawLine>& hb_tail)
{
    std::ostringstream out;
    switch (type) {
    case LineType::trigger: {
        const auto timestamp = extract_trigger_timestamp(row);
        out << "ts=" << format_compact_hex(timestamp);
        break;
    }
    case LineType::data_idle:
    case LineType::data_daq: {
        const auto fmc = frame_word(row, 4);
        const auto fcmd = describe_fmc_word(fmc);
        out << color_section("gbt", gbt_color(row[1])) << ' ' << color_section(std::to_string(row[1]), gbt_color(row[1]))
            << ' ' << color_section("FCMD", fcmd_color(fcmd)) << ' ' << color_section(fcmd, fcmd_color(fcmd))
            << ' ' << color_section("cnt", "\033[1;31m") << ' ' << color_section(format_compact_hex(reordered_word_value(frame_word(row, 5))), "\033[1;31m");
        for (std::size_t word_index = 0; word_index < 4; ++word_index) {
            out << ' ' << color_section("w" + std::to_string(word_index), word_color(word_index))
                << ' ' << color_section(format_word_bytes(frame_word(row, word_index)), word_color(word_index));
        }
        break;
    }
    case LineType::heartbeat:
        out << describe_heartbeat(row, hb_tail);
        break;
    case LineType::dummy:
        break;
    case LineType::other:
        out << "other [XX] line first_byte=" << format_compact_hex(row[0]);
        break;
    case LineType::partial_tail:
        break;
    }
    return out.str();
}

std::string trim(std::string_view text)
{
    std::size_t first = 0;
    while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first]))) {
        ++first;
    }
    std::size_t last = text.size();
    while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1]))) {
        --last;
    }
    return std::string(text.substr(first, last - first));
}

std::size_t parse_positive_index(const std::string& text)
{
    if (text.empty()) {
        throw std::invalid_argument("empty line index");
    }
    std::size_t consumed = 0;
    const auto value = std::stoull(text, &consumed, 0);
    if (consumed != text.size() || value == 0) {
        throw std::invalid_argument("line indices are 1-based positive integers: " + text);
    }
    return static_cast<std::size_t>(value);
}

std::uintmax_t parse_size_limit(const std::string& text)
{
    const auto value_text = trim(text);
    if (value_text.empty()) {
        throw std::invalid_argument("empty size limit");
    }

    std::size_t consumed = 0;
    const double value = std::stod(value_text, &consumed);
    if (value <= 0.0) {
        throw std::invalid_argument("size limit must be positive: " + text);
    }

    const auto unit = trim(std::string_view(value_text).substr(consumed));
    std::uintmax_t multiplier = 1;
    if (unit.empty() || unit == "B" || unit == "b") {
        multiplier = 1;
    } else if (unit == "K" || unit == "KB" || unit == "Kb" || unit == "k" || unit == "kb") {
        multiplier = 1000;
    } else if (unit == "M" || unit == "MB" || unit == "Mb" || unit == "m" || unit == "mb") {
        multiplier = 1000 * 1000;
    } else if (unit == "G" || unit == "GB" || unit == "Gb" || unit == "g" || unit == "gb") {
        multiplier = 1000 * 1000 * 1000;
    } else if (unit == "T" || unit == "TB" || unit == "Tb" || unit == "t" || unit == "tb") {
        multiplier = 1000ULL * 1000ULL * 1000ULL * 1000ULL;
    } else {
        throw std::invalid_argument("unsupported size unit: " + unit);
    }

    return static_cast<std::uintmax_t>(value * static_cast<double>(multiplier));
}

bool looks_like_line_range_selector(const std::string& text)
{
    return text.find(',') != std::string::npos || text.find('-') != std::string::npos;
}

std::vector<LineRange> parse_line_ranges(const std::string& text)
{
    std::vector<LineRange> ranges;
    std::size_t start = 0;
    while (start < text.size()) {
        const std::size_t comma = text.find(',', start);
        const auto token = trim(std::string_view(text).substr(start, comma == std::string::npos ? std::string::npos : comma - start));
        if (!token.empty()) {
            const auto dash = token.find('-');
            if (dash == std::string::npos) {
                const auto line_index = parse_positive_index(token);
                ranges.emplace_back(line_index, line_index);
            } else {
                const auto first = parse_positive_index(trim(std::string_view(token).substr(0, dash)));
                const auto last = parse_positive_index(trim(std::string_view(token).substr(dash + 1)));
                if (last < first) {
                    throw std::invalid_argument("line range ends before it starts: " + token);
                }
                ranges.emplace_back(first, last);
            }
        }
        if (comma == std::string::npos) {
            break;
        }
        start = comma + 1;
    }
    return ranges;
}

DisplayOptions parse_display_selector(const std::string& text, DisplayMode mode)
{
    DisplayOptions options;
    options.mode = mode;
    if (looks_like_line_range_selector(text)) {
        options.selected_ranges = parse_line_ranges(text);
        options.use_selected_ranges = true;
        return options;
    }
    options.limit = parse_positive_index(trim(text));
    return options;
}

bool line_matches_display_mode(const RawLine& row, LineType type, DisplayMode mode)
{
    switch (mode) {
    case DisplayMode::raw:
        return true;
    case DisplayMode::non_dummy_non_heartbeat:
        return type != LineType::dummy && type != LineType::heartbeat;
    case DisplayMode::data_idle:
        return type == LineType::data_idle;
    case DisplayMode::data_daq:
        return type == LineType::data_daq;
    case DisplayMode::heartbeat:
        return type == LineType::heartbeat;
    case DisplayMode::trigger:
        return type == LineType::trigger;
    case DisplayMode::dummy:
        return type == LineType::dummy;
    case DisplayMode::other:
        return type == LineType::other;
    }
    return false;
}

bool should_print_line(std::size_t line_index, const RawLine& row, LineType type, const DisplayOptions& display_options, std::size_t matched_print_lines)
{
    if (display_options.use_selected_ranges) {
        const bool selected = std::any_of(display_options.selected_ranges.begin(), display_options.selected_ranges.end(), [line_index](const LineRange& range) {
            return range.first <= line_index && line_index <= range.second;
        });
        return selected && line_matches_display_mode(row, type, display_options.mode);
    }
    return matched_print_lines < display_options.limit && line_matches_display_mode(row, type, display_options.mode);
}

void print_lines(const std::vector<PrintableLine>& lines)
{
    std::size_t index_width = 1;
    for (const auto& line : lines) {
        index_width = std::max(index_width, std::to_string(line.line_index).size());
    }

    for (const auto& line : lines) {
        const auto detail = describe_line(line.row, line.type, line.last_trigger_timestamp, line.hb_tail);
        if (detail.empty()) {
            spdlog::info("L{:>{}} {} {}", line.line_index, index_width, color_line_type_code(line.type), format_bytes(line.row.data(), line.row.size()));
            continue;
        }
        if (line.type == LineType::data_idle || line.type == LineType::data_daq) {
            spdlog::info("L{:>{}} {} {}", line.line_index, index_width, color_line_type_code(line.type), format_data_bytes(line.row));
            spdlog::info("{}{}", std::string(index_width + 7, ' '), detail);
            continue;
        }
        if (line.type == LineType::heartbeat && line.hb_tail.has_value()) {
            spdlog::info("L{:>{}} {} {}", line.line_index, index_width, color_line_type_code(line.type), format_bytes(line.row.data(), line.row.size()));
            spdlog::info("L{:>{}} {} {}", line.line_index + 1, index_width, color_line_type_code(line.type), format_bytes(line.hb_tail->data(), line.hb_tail->size()));
            const auto detail_indent = std::string(index_width + 7, ' ');
            const auto line_break = detail.find('\n');
            if (line_break == std::string::npos) {
                spdlog::info("{}{}", detail_indent, detail);
            } else {
                spdlog::info("{}{}", detail_indent, detail.substr(0, line_break));
                spdlog::info("{}{}", detail_indent, detail.substr(line_break + 1));
            }
            continue;
        }
        spdlog::info("L{:>{}} {} {}    || {}", line.line_index, index_width, color_line_type_code(line.type), format_bytes(line.row.data(), line.row.size()), detail);
    }
}

ParseSummary scan_raw_file(const std::string& input_path, const DisplayOptions& display_options, std::optional<std::size_t> max_rows, std::optional<std::uintmax_t> max_bytes)
{
    const auto start_time = std::chrono::steady_clock::now();

    std::ifstream input(input_path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("failed to open input file: " + input_path);
    }

    ParseSummary summary;
    summary.input_path = input_path;

    std::optional<std::uint64_t> last_trigger_timestamp;
    RawLine row{};
    std::size_t matched_print_lines = 0;
    bool incomplete_read = false;
    std::vector<PrintableLine> printable_lines;
    while ((!max_rows.has_value() || summary.total_lines < *max_rows) && (!max_bytes.has_value() || summary.total_bytes_parsed + row.size() <= *max_bytes)) {
        if (!input.read(reinterpret_cast<char*>(row.data()), static_cast<std::streamsize>(row.size()))) {
            incomplete_read = true;
            break;
        }
        ++summary.total_lines;
        summary.total_bytes_parsed += row.size();
        const auto line_index = summary.total_lines;
        const auto type = classify_line(row);
        increment_stats(summary.stats, type);

        std::optional<RawLine> hb_tail;
        if (type == LineType::heartbeat && (!max_rows.has_value() || summary.total_lines < *max_rows) && (!max_bytes.has_value() || summary.total_bytes_parsed + row.size() <= *max_bytes)) {
            RawLine tail{};
            if (input.read(reinterpret_cast<char*>(tail.data()), static_cast<std::streamsize>(tail.size()))) {
                ++summary.total_lines;
                summary.total_bytes_parsed += tail.size();
                increment_stats(summary.stats, LineType::heartbeat);
                hb_tail = tail;
            } else {
                incomplete_read = true;
            }
        }

        if (should_print_line(line_index, row, type, display_options, matched_print_lines)) {
            printable_lines.push_back(PrintableLine{line_index, row, type, last_trigger_timestamp, hb_tail});
            ++summary.printed_lines;
            if (!display_options.use_selected_ranges) {
                ++matched_print_lines;
            }
        }

        if (type == LineType::trigger) {
            last_trigger_timestamp = extract_trigger_timestamp(row);
        }
    }

    const bool stopped_by_size = max_bytes.has_value() && summary.total_bytes_parsed < *max_bytes && summary.total_bytes_parsed + row.size() > *max_bytes;
    if (stopped_by_size) {
        const auto remaining_bytes = static_cast<std::size_t>(*max_bytes - summary.total_bytes_parsed);
        std::vector<char> tail(remaining_bytes);
        if (remaining_bytes != 0) {
            input.read(tail.data(), static_cast<std::streamsize>(tail.size()));
            summary.partial_tail_bytes = static_cast<std::size_t>(input.gcount());
        }
    } else if (incomplete_read) {
        summary.partial_tail_bytes = static_cast<std::size_t>(input.gcount());
    }
    if (input.bad()) {
        throw std::runtime_error("failed while reading input file: " + input_path);
    }
    if (summary.partial_tail_bytes != 0) {
        increment_stats(summary.stats, LineType::partial_tail);
    }

    std::ifstream size_input(input_path, std::ios::binary | std::ios::ate);
    if (size_input) {
        summary.file_size = static_cast<std::uintmax_t>(size_input.tellg());
    }

    const auto end_time = std::chrono::steady_clock::now();
    summary.elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();
    print_lines(printable_lines);
    return summary;
}

void print_summary(const ParseSummary& summary)
{
    const auto parsed_bytes = static_cast<std::uintmax_t>(summary.total_bytes_parsed + summary.partial_tail_bytes);

    print_box("File summary", {
        {"File", file_name_only(summary.input_path)},
        {"Size", format_byte_size(parsed_bytes)},
        {"Parsed lines", format_count(summary.total_lines) + " x " + std::to_string(kRowSize) + " bytes"},
        {"Speed", format_byte_rate(summary.total_bytes_parsed, summary.elapsed_sec) + " in " + format_seconds(summary.elapsed_sec)},
        {"Printed lines", format_count(summary.printed_lines)},
    });

    const auto line_stat_value = [&summary](std::size_t count) {
        const double percentage = summary.total_lines > 0
            ? 100.0 * static_cast<double>(count) / static_cast<double>(summary.total_lines)
            : 0.0;
        std::ostringstream out;
        out << std::setw(12) << format_count(count) << " (" << std::setw(6) << std::fixed << std::setprecision(2) << percentage << "%)";
        return out.str();
    };

    std::vector<std::pair<std::string, std::string>> stat_rows{
        {std::string(line_type_name(LineType::heartbeat)), line_stat_value(summary.stats.heartbeat)},
        {std::string(line_type_name(LineType::trigger)), line_stat_value(summary.stats.trigger)},
        {std::string(line_type_name(LineType::data_idle)), line_stat_value(summary.stats.data_idle)},
        {std::string(line_type_name(LineType::data_daq)), line_stat_value(summary.stats.data_daq)},
        {std::string(line_type_name(LineType::dummy)), line_stat_value(summary.stats.dummy)},
        {std::string(line_type_name(LineType::other)), line_stat_value(summary.stats.other)},
    };
    if (summary.partial_tail_bytes != 0) {
        stat_rows.emplace_back("partial tail", format_byte_size(summary.partial_tail_bytes));
    }
    print_box("Line type statistics", stat_rows);
}

void write_root_summary(const std::string& output_path, const ParseSummary& summary)
{
    TFile output(output_path.c_str(), "RECREATE");
    if (output.IsZombie()) {
        throw std::runtime_error("failed to create ROOT file: " + output_path);
    }

    TNamed("Rootifier_version", "001 raw line browser").Write();
    TNamed("Rootifier_input_file", summary.input_path.c_str()).Write();
    TNamed("Rootifier_total_lines", std::to_string(summary.total_lines).c_str()).Write();
    TNamed("Rootifier_HB_heartbeat_lines", std::to_string(summary.stats.heartbeat).c_str()).Write();
    TNamed("Rootifier_PT_phy_trigger_lines", std::to_string(summary.stats.trigger).c_str()).Write();
    TNamed("Rootifier_AC_data_idle_lines", std::to_string(summary.stats.data_idle).c_str()).Write();
    TNamed("Rootifier_DQ_data_daq_lines", std::to_string(summary.stats.data_daq).c_str()).Write();
    TNamed("Rootifier_00_dummy_zero_lines", std::to_string(summary.stats.dummy).c_str()).Write();
    TNamed("Rootifier_XX_other_lines", std::to_string(summary.stats.other).c_str()).Write();
    output.Write();
    output.Close();
}

} // namespace

int main(int argc, char** argv)
{
    std::string script_full_name = argv[0];
    auto script_name = script_full_name.substr(script_full_name.find_last_of("/\\") + 1);
    spdlog::info("Running script: {}", script_name);
    spdlog::info("----------------------------------------");

    cxxopts::Options options(script_name, "Browse raw 32-byte CH2G lines and print selected line indices");
    options.add_options()
        ("i,input", "Input raw data file", cxxopts::value<std::string>())
        ("o,output", "Optional ROOT summary output file", cxxopts::value<std::string>())
        ("line-ac", "Print data idle [AC] lines by count or 1-based indices/ranges, e.g. 20 or 1,2,40-48", cxxopts::value<std::string>())
        ("line-dq", "Print data daq [DQ] lines by count or 1-based indices/ranges, e.g. 20 or 1,2,40-48", cxxopts::value<std::string>())
        ("line-hb", "Print heartbeat [HB] lines by count or 1-based indices/ranges", cxxopts::value<std::string>())
        ("line-bb", "Print phy trigger [PT] lines by count or 1-based indices/ranges", cxxopts::value<std::string>())
        ("line-00", "Print dummy zero [00] lines by count or 1-based indices/ranges", cxxopts::value<std::string>())
        ("line-xx", "Print other [XX] lines by count or 1-based indices/ranges", cxxopts::value<std::string>())
        ("line", "Print non-[00] and non-[HB] lines by count or 1-based indices/ranges", cxxopts::value<std::string>())
        ("line-raw", "Print raw lines by count or 1-based indices/ranges", cxxopts::value<std::string>())
        ("max-rows", "Stop after this many full 32-byte rows", cxxopts::value<std::size_t>())
        ("size", "Stop after parsing this much input, e.g. 20KB, 20MB, 20Mb", cxxopts::value<std::string>())
        ("v,verbose", "Verbose mode kept for compatibility with 001", cxxopts::value<int>()->default_value("2"))
        ("h,help", "Print help");

    cxxopts::ParseResult parsed_opts;
    try {
        parsed_opts = options.parse(argc, argv);
    } catch (const std::exception& error) {
        spdlog::error("Argument error: {}", error.what());
        std::cout << options.help() << std::endl;
        return 2;
    }

    if (parsed_opts.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (!parsed_opts.count("input")) {
        spdlog::error("Input file is required. Use --input <path>.");
        std::cout << options.help() << std::endl;
        return 2;
    }

    try {
        const auto input_path = parsed_opts["input"].as<std::string>();
        DisplayOptions display_options;

        std::size_t line_mode_count = 0;
        const auto set_display_mode = [&](std::string_view option_name, DisplayMode mode) {
            if (!parsed_opts.count(std::string(option_name))) {
                return;
            }
            display_options = parse_display_selector(parsed_opts[std::string(option_name)].as<std::string>(), mode);
            ++line_mode_count;
        };
        set_display_mode("line-ac", DisplayMode::data_idle);
        set_display_mode("line-dq", DisplayMode::data_daq);
        set_display_mode("line-hb", DisplayMode::heartbeat);
        set_display_mode("line-bb", DisplayMode::trigger);
        set_display_mode("line-00", DisplayMode::dummy);
        set_display_mode("line-xx", DisplayMode::other);
        set_display_mode("line", DisplayMode::non_dummy_non_heartbeat);
        set_display_mode("line-raw", DisplayMode::raw);
        if (line_mode_count > 1) {
            throw std::invalid_argument("Use only one of --line, --line-ac, --line-dq, --line-hb, --line-bb, --line-00, --line-xx, or --line-raw at a time");
        }
        const std::optional<std::size_t> max_rows = parsed_opts.count("max-rows")
            ? std::optional<std::size_t>(parsed_opts["max-rows"].as<std::size_t>())
            : std::nullopt;
        const std::optional<std::uintmax_t> max_bytes = parsed_opts.count("size")
            ? std::optional<std::uintmax_t>(parse_size_limit(parsed_opts["size"].as<std::string>()))
            : std::nullopt;

        const auto summary = scan_raw_file(input_path, display_options, max_rows, max_bytes);
        print_summary(summary);

        if (parsed_opts.count("output")) {
            const auto output_path = parsed_opts["output"].as<std::string>();
            write_root_summary(output_path, summary);
            spdlog::info("Wrote ROOT summary to {}", output_path);
        }
    } catch (const std::exception& error) {
        spdlog::error("{}", error.what());
        return 1;
    }

    return 0;
}