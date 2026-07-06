#include <spdlog/spdlog.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TTree.h>

#include <binparse/parser.hpp>
#include <binparse/sample.hpp>

#include "cxxopts.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstring>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <ctime>
#include <vector>

namespace {

constexpr std::size_t kRowSize = 32;
constexpr std::size_t kWordSize = 4;
constexpr std::size_t kFrameWordOffset = 8;
constexpr std::size_t kFrameWordCount = 6;
constexpr std::size_t kTimestampHexWidth = 12;
constexpr std::size_t kInputBufferSize = 16 * 1024 * 1024;

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

struct ReconEvent {
    std::size_t trigger_line = 0;
    LineStats counts;
    std::vector<bp::DataLine> dq_lines;
    std::vector<std::size_t> dq_line_indices;
    std::vector<bp::Sample> samples;
};

struct ReconSummary {
    std::string input_path;
    std::filesystem::path output_dir;
    std::filesystem::path pdf_path;
    std::filesystem::path sample_root_path;
    std::size_t total_lines = 0;
    std::size_t total_bytes_parsed = 0;
    std::size_t partial_tail_bytes = 0;
    double elapsed_sec = 0.0;
    LineStats stats;
    std::uintmax_t cached_event_memory_bytes = 0;
    std::optional<std::uintmax_t> memory_limit_bytes;
    std::vector<ReconEvent> events;
};

struct ReconStatistics {
    std::string input_path;
    std::filesystem::path output_dir;
    std::filesystem::path pdf_path;
    std::filesystem::path sample_root_path;
    std::size_t total_lines = 0;
    std::size_t total_bytes_parsed = 0;
    std::size_t partial_tail_bytes = 0;
    double elapsed_sec = 0.0;
    LineStats stats;
    std::size_t pt_events = 0;
    std::size_t total_dq_lines = 0;
    std::size_t total_samples = 0;
    std::size_t total_sample_crc_ok_elinks = 0;
    std::size_t total_samples_all_elinks_crc_ok = 0;
    std::array<std::vector<std::size_t>, 7> event_counts_by_type;
    std::array<std::vector<std::size_t>, 256> dq_cluster_lengths_by_gbt;
    std::array<std::vector<std::size_t>, 256> sample_cluster_lengths_by_gbt_all;
    std::array<std::vector<std::size_t>, 256> sample_cluster_lengths_by_gbt_crc;
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
    for (const auto value : row) {
        if (value != 0x00) {
            return false;
        }
    }
    return true;
}

bool frame_word_matches(const RawLine& row, std::size_t word_index, const Word& expected)
{
    const auto offset = kFrameWordOffset + word_index * kWordSize;
    for (std::size_t byte_index = 0; byte_index < kWordSize; ++byte_index) {
        if (row[offset + byte_index] != expected[byte_index]) {
            return false;
        }
    }
    return true;
}

std::uint32_t reordered_word_value_from_row(const RawLine& row, std::size_t word_index)
{
    const auto offset = kFrameWordOffset + word_index * kWordSize;
    return static_cast<std::uint32_t>(row[offset + 3]) |
        (static_cast<std::uint32_t>(row[offset + 2]) << 8) |
        (static_cast<std::uint32_t>(row[offset + 1]) << 16) |
        (static_cast<std::uint32_t>(row[offset]) << 24);
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

bp::DataLine parse_data_line_struct(const RawLine& row)
{
    const auto stream_timestamp = extract_stream_timestamp(row);
    bp::DataLine line{};
    line.header_vldb_id = row[1];
    line.bx_cnt = static_cast<std::uint16_t>(stream_timestamp & 0xFFFFu);
    line.ob_cnt = static_cast<std::uint32_t>(stream_timestamp >> 16);
    line.data_word0 = reordered_word_value_from_row(row, 0);
    line.data_word1 = reordered_word_value_from_row(row, 1);
    line.data_word2 = reordered_word_value_from_row(row, 2);
    line.data_word3 = reordered_word_value_from_row(row, 3);
    line.data_word4 = reordered_word_value_from_row(row, 4);
    line.data_word5 = reordered_word_value_from_row(row, 5);
    return line;
}

bool data_lanes_all_match(const RawLine& row, const Word& expected)
{
    for (std::size_t word_index = 0; word_index < 4; ++word_index) {
        if (!frame_word_matches(row, word_index, expected)) {
            return false;
        }
    }
    return true;
}

bool data_lanes_have_activity(const RawLine& row)
{
    for (std::size_t word_index = 0; word_index < 4; ++word_index) {
        if (!frame_word_matches(row, word_index, kZeroWord) && !frame_word_matches(row, word_index, kFrameIdleWord)) {
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

std::string format_count_with_percentage(std::size_t count, std::size_t total)
{
    const double percentage = total == 0 ? 0.0 : 100.0 * static_cast<double>(count) / static_cast<double>(total);
    std::ostringstream out;
    out << format_count(count) << " (" << std::fixed << std::setprecision(2) << percentage << "%)";
    return out.str();
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

std::string read_failure_message(const std::string& input_path, std::size_t total_lines, std::uintmax_t total_bytes_parsed, std::streamsize last_gcount)
{
    std::ostringstream message;
    message << "failed while reading input file: " << input_path
        << " after " << format_count(total_lines) << " full rows ("
        << format_byte_size(total_bytes_parsed) << ", offset " << format_count(static_cast<std::size_t>(total_bytes_parsed)) << " bytes)"
        << ", last gcount " << last_gcount;
    if (errno != 0) {
        message << ", errno " << errno << " (" << std::strerror(errno) << ")";
    }
    return message.str();
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

std::string join_arguments(int argc, char** argv)
{
    std::ostringstream out;
    bool first = true;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--input" || argument == "-i") {
            ++index;
            continue;
        }
        if (argument.rfind("--input=", 0) == 0 || argument.rfind("-i=", 0) == 0) {
            continue;
        }
        if (!first) {
            out << ' ';
        }
        out << argument;
        first = false;
    }
    return out.str();
}

std::string current_time_minute()
{
    const auto now = std::time(nullptr);
    std::tm local_time{};
    localtime_r(&now, &local_time);
    std::ostringstream out;
    out << std::put_time(&local_time, "%Y-%m-%d %H:%M");
    return out.str();
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
        if (index == 0) {
            out << white;
        } else if (index == 2) {
            out << reset << gray;
        } else if (index == kFrameWordOffset) {
            out << reset << white;
        } else if (index > kFrameWordOffset && (index - kFrameWordOffset) % kWordSize == 0) {
            if (index != 0) {
                out << reset;
            }
            out << ((((index - kFrameWordOffset) / kWordSize) % 2 == 0) ? white : gray);
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

    std::vector<char> input_buffer(kInputBufferSize);
    std::ifstream input;
    input.rdbuf()->pubsetbuf(input_buffer.data(), static_cast<std::streamsize>(input_buffer.size()));
    input.open(input_path, std::ios::binary);
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
        throw std::runtime_error(read_failure_message(input_path, summary.total_lines, summary.total_bytes_parsed, input.gcount()));
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

std::filesystem::path recon_output_dir(const std::string& input_path)
{
    return std::filesystem::path("dump") / "001" / file_name_only(input_path);
}

void fill_event_count(ReconEvent& event, LineType type)
{
    increment_stats(event.counts, type);
}

std::size_t total_recon_dq_lines(const std::vector<ReconEvent>& events)
{
    std::size_t total = 0;
    for (const auto& event : events) {
        total += event.dq_lines.size();
    }
    return total;
}

std::size_t total_recon_samples(const std::vector<ReconEvent>& events)
{
    std::size_t total = 0;
    for (const auto& event : events) {
        total += event.samples.size();
    }
    return total;
}

std::size_t total_sample_crc_ok_elinks(const std::vector<ReconEvent>& events)
{
    std::size_t total = 0;
    for (const auto& event : events) {
        for (const auto& sample : event.samples) {
            for (const auto& elink : sample.elinks) {
                if (elink.crc_ok) {
                    ++total;
                }
            }
        }
    }
    return total;
}

bool sample_all_elinks_crc_ok(const bp::Sample& sample)
{
    return std::all_of(sample.elinks.begin(), sample.elinks.end(), [](const auto& elink) { return elink.crc_ok; });
}

std::size_t total_samples_all_elinks_crc_ok(const std::vector<ReconEvent>& events)
{
    std::size_t total = 0;
    for (const auto& event : events) {
        for (const auto& sample : event.samples) {
            if (sample_all_elinks_crc_ok(sample)) {
                ++total;
            }
        }
    }
    return total;
}

std::size_t total_cluster_count(const std::array<std::vector<std::size_t>, 256>& lengths_by_gbt)
{
    std::size_t total = 0;
    for (const auto& lengths : lengths_by_gbt) {
        total += lengths.size();
    }
    return total;
}

std::optional<std::size_t> peak_cluster_length(const std::array<std::vector<std::size_t>, 256>& lengths_by_gbt)
{
    std::map<std::size_t, std::size_t> count_frequencies;
    for (const auto& lengths : lengths_by_gbt) {
        for (const auto length : lengths) {
            ++count_frequencies[length];
        }
    }
    if (count_frequencies.empty()) {
        return std::nullopt;
    }
    const auto peak_iter = std::max_element(count_frequencies.begin(), count_frequencies.end(), [](const auto& left, const auto& right) {
        return left.second < right.second;
    });
    return peak_iter->first;
}

std::array<bool, 256> active_gbts_from_cluster_lengths(const std::array<std::vector<std::size_t>, 256>& lengths_by_gbt)
{
    std::array<bool, 256> active_gbts{};
    for (std::size_t gbt_index = 0; gbt_index < lengths_by_gbt.size(); ++gbt_index) {
        active_gbts[gbt_index] = !lengths_by_gbt[gbt_index].empty();
    }
    return active_gbts;
}

bool has_active_gbt(const std::array<bool, 256>& active_gbts)
{
    return std::any_of(active_gbts.begin(), active_gbts.end(), [](bool active) { return active; });
}

std::uint32_t compact_word_value(const std::array<std::uint8_t, kWordSize>& word)
{
    return static_cast<std::uint32_t>(word[0]) |
        (static_cast<std::uint32_t>(word[1]) << 8) |
        (static_cast<std::uint32_t>(word[2]) << 16) |
        (static_cast<std::uint32_t>(word[3]) << 24);
}

bool data_line_has_l1a_fcmd(const bp::DataLine& line)
{
    return line.data_word4 == compact_word_value(kFmcL1aWord);
}

void collect_samples_from_cluster(ReconEvent& event, const std::vector<bp::DataLine>& cluster_lines)
{
    for (std::size_t first_line = 0; first_line + bp::kSampleLegalLineCount <= cluster_lines.size();) {
        if (auto sample = bp::make_sample(cluster_lines, first_line)) {
            event.samples.push_back(*sample);
            first_line += bp::kSampleLegalLineCount;
            continue;
        }
        ++first_line;
    }
}

void build_event_samples(ReconEvent& event)
{
    std::array<std::vector<bp::DataLine>, 256> lines_by_gbt;
    for (const auto& line : event.dq_lines) {
        lines_by_gbt[line.header_vldb_id].push_back(line);
    }

    for (auto& lines : lines_by_gbt) {
        if (lines.empty()) {
            continue;
        }
        std::sort(lines.begin(), lines.end(), [](const auto& left, const auto& right) {
            return bp::data_line_timestamp(left) < bp::data_line_timestamp(right);
        });

        std::vector<bp::DataLine> cluster_lines;
        cluster_lines.push_back(lines.front());
        for (std::size_t index = 1; index < lines.size(); ++index) {
            if (bp::data_line_timestamp(lines[index]) == bp::data_line_timestamp(lines[index - 1]) + 1) {
                cluster_lines.push_back(lines[index]);
                continue;
            }
            collect_samples_from_cluster(event, cluster_lines);
            cluster_lines.clear();
            cluster_lines.push_back(lines[index]);
        }
        collect_samples_from_cluster(event, cluster_lines);
    }
}

std::uintmax_t estimated_recon_event_memory_bytes(const ReconEvent& event)
{
    return sizeof(ReconEvent) +
        event.dq_lines.capacity() * sizeof(bp::DataLine) +
        event.dq_line_indices.capacity() * sizeof(std::size_t) +
        event.samples.capacity() * sizeof(bp::Sample);
}

void check_recon_memory_budget(const ReconSummary& summary, const ReconEvent& pending_event)
{
    if (!summary.memory_limit_bytes.has_value()) {
        return;
    }
    const auto estimated_bytes = summary.cached_event_memory_bytes + estimated_recon_event_memory_bytes(pending_event);
    if (estimated_bytes <= *summary.memory_limit_bytes) {
        return;
    }
    std::ostringstream message;
    message << "001 --recon stopped before exceeding the in-memory event cache budget: estimated "
            << format_byte_size(estimated_bytes) << " > limit " << format_byte_size(*summary.memory_limit_bytes)
            << ". Reduce input with --size, process a smaller raw file, or raise/disable --memory-limit.";
    throw std::runtime_error(message.str());
}

void push_recon_event(ReconSummary& summary, ReconEvent event)
{
    build_event_samples(event);
    check_recon_memory_budget(summary, event);
    summary.cached_event_memory_bytes += estimated_recon_event_memory_bytes(event);
    summary.events.push_back(std::move(event));
}

ReconSummary reconstruct_events(const std::string& input_path, std::optional<std::size_t> max_rows, std::optional<std::uintmax_t> max_bytes,
    std::optional<std::uintmax_t> memory_limit_bytes)
{
    const auto start_time = std::chrono::steady_clock::now();

    std::vector<char> input_buffer(kInputBufferSize);
    std::ifstream input;
    input.rdbuf()->pubsetbuf(input_buffer.data(), static_cast<std::streamsize>(input_buffer.size()));
    input.open(input_path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("failed to open input file: " + input_path);
    }

    ReconSummary summary;
    summary.input_path = input_path;
    summary.output_dir = recon_output_dir(input_path);
    summary.pdf_path = summary.output_dir / "pt_line_type_distributions.pdf";
    summary.sample_root_path = summary.output_dir / "sample_events.root";
    summary.memory_limit_bytes = memory_limit_bytes;

    std::optional<ReconEvent> current_event;
    RawLine row{};
    bool incomplete_read = false;
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
        if (type == LineType::other) {
            continue;
        }

        if (type == LineType::trigger) {
            if (current_event.has_value()) {
                push_recon_event(summary, std::move(*current_event));
            }
            current_event = ReconEvent{};
            current_event->trigger_line = line_index;
            fill_event_count(*current_event, type);
            continue;
        }

        if (type == LineType::heartbeat && (!max_rows.has_value() || summary.total_lines < *max_rows) && (!max_bytes.has_value() || summary.total_bytes_parsed + row.size() <= *max_bytes)) {
            RawLine tail{};
            if (input.read(reinterpret_cast<char*>(tail.data()), static_cast<std::streamsize>(tail.size()))) {
                ++summary.total_lines;
                summary.total_bytes_parsed += tail.size();
                increment_stats(summary.stats, LineType::heartbeat);
                if (current_event.has_value()) {
                    fill_event_count(*current_event, LineType::heartbeat);
                    fill_event_count(*current_event, LineType::heartbeat);
                }
                continue;
            }
            incomplete_read = true;
        }

        if (current_event.has_value()) {
            fill_event_count(*current_event, type);
            if (type == LineType::data_daq) {
                current_event->dq_lines.push_back(parse_data_line_struct(row));
                current_event->dq_line_indices.push_back(line_index);
                check_recon_memory_budget(summary, *current_event);
            }
        }
    }

    if (current_event.has_value()) {
        push_recon_event(summary, std::move(*current_event));
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
        throw std::runtime_error(read_failure_message(input_path, summary.total_lines, summary.total_bytes_parsed, input.gcount()));
    }
    if (summary.partial_tail_bytes != 0) {
        increment_stats(summary.stats, LineType::partial_tail);
    }

    const auto end_time = std::chrono::steady_clock::now();
    summary.elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();

    return summary;
}

std::vector<std::size_t> event_counts_for_type(const std::vector<ReconEvent>& events, LineType type)
{
    std::vector<std::size_t> values;
    values.reserve(events.size());
    for (const auto& event : events) {
        switch (type) {
        case LineType::heartbeat:
            values.push_back(event.counts.heartbeat);
            break;
        case LineType::trigger:
            values.push_back(event.counts.trigger);
            break;
        case LineType::data_idle:
            values.push_back(event.counts.data_idle);
            break;
        case LineType::data_daq:
            values.push_back(event.counts.data_daq);
            break;
        case LineType::dummy:
            values.push_back(event.counts.dummy);
            break;
        case LineType::other:
            values.push_back(event.counts.other);
            break;
        case LineType::partial_tail:
            values.push_back(event.counts.partial_tail);
            break;
        }
    }
    return values;
}

std::vector<std::size_t> dq_cluster_lengths(const std::vector<ReconEvent>& events)
{
    std::vector<std::size_t> lengths;
    for (const auto& event : events) {
        std::array<std::vector<std::uint64_t>, 256> timestamps_by_gbt;
        for (const auto& line : event.dq_lines) {
            timestamps_by_gbt[line.header_vldb_id].push_back(bp::data_line_timestamp(line));
        }

        for (auto& timestamps : timestamps_by_gbt) {
            if (timestamps.empty()) {
                continue;
            }
            std::sort(timestamps.begin(), timestamps.end());
            std::size_t cluster_length = 1;
            for (std::size_t index = 1; index < timestamps.size(); ++index) {
                if (timestamps[index] == timestamps[index - 1] + 1) {
                    ++cluster_length;
                    continue;
                }
                lengths.push_back(cluster_length);
                cluster_length = 1;
            }
            lengths.push_back(cluster_length);
        }
    }
    return lengths;
}

std::array<std::vector<std::size_t>, 256> dq_cluster_lengths_by_gbt(const std::vector<ReconEvent>& events)
{
    std::array<std::vector<std::size_t>, 256> lengths_by_gbt;
    for (const auto& event : events) {
        std::array<std::vector<std::uint64_t>, 256> timestamps_by_gbt;
        for (const auto& line : event.dq_lines) {
            timestamps_by_gbt[line.header_vldb_id].push_back(bp::data_line_timestamp(line));
        }

        for (std::size_t gbt_index = 0; gbt_index < timestamps_by_gbt.size(); ++gbt_index) {
            auto& timestamps = timestamps_by_gbt[gbt_index];
            if (timestamps.empty()) {
                continue;
            }
            std::sort(timestamps.begin(), timestamps.end());
            std::size_t cluster_length = 1;
            for (std::size_t index = 1; index < timestamps.size(); ++index) {
                if (timestamps[index] == timestamps[index - 1] + 1) {
                    ++cluster_length;
                    continue;
                }
                lengths_by_gbt[gbt_index].push_back(cluster_length);
                cluster_length = 1;
            }
            lengths_by_gbt[gbt_index].push_back(cluster_length);
        }
    }
    return lengths_by_gbt;
}

std::array<std::vector<std::size_t>, 256> sample_cluster_lengths_by_gbt(const std::vector<ReconEvent>& events, bool crc_only)
{
    std::array<std::vector<std::size_t>, 256> lengths_by_gbt;
    for (const auto& event : events) {
        std::array<std::vector<std::uint64_t>, 256> first_timestamps_by_gbt;
        for (const auto& sample : event.samples) {
            if (crc_only && !sample_all_elinks_crc_ok(sample)) {
                continue;
            }
            first_timestamps_by_gbt[sample.gbt_index].push_back(sample.first_timestamp);
        }

        for (std::size_t gbt_index = 0; gbt_index < first_timestamps_by_gbt.size(); ++gbt_index) {
            auto& timestamps = first_timestamps_by_gbt[gbt_index];
            if (timestamps.empty()) {
                continue;
            }
            std::sort(timestamps.begin(), timestamps.end());
            std::size_t cluster_length = 1;
            for (std::size_t index = 1; index < timestamps.size(); ++index) {
                if (timestamps[index] == timestamps[index - 1] + bp::kSampleLegalLineCount + 1) {
                    ++cluster_length;
                    continue;
                }
                lengths_by_gbt[gbt_index].push_back(cluster_length);
                cluster_length = 1;
            }
            lengths_by_gbt[gbt_index].push_back(cluster_length);
        }
    }
    return lengths_by_gbt;
}

std::array<std::vector<std::size_t>, 256> sample_cluster_lengths_by_gbt(const ReconEvent& event, bool crc_only)
{
    std::array<std::vector<std::size_t>, 256> lengths_by_gbt;
    std::array<std::vector<std::uint64_t>, 256> first_timestamps_by_gbt;
    for (const auto& sample : event.samples) {
        if (crc_only && !sample_all_elinks_crc_ok(sample)) {
            continue;
        }
        first_timestamps_by_gbt[sample.gbt_index].push_back(sample.first_timestamp);
    }

    for (std::size_t gbt_index = 0; gbt_index < first_timestamps_by_gbt.size(); ++gbt_index) {
        auto& timestamps = first_timestamps_by_gbt[gbt_index];
        if (timestamps.empty()) {
            continue;
        }
        std::sort(timestamps.begin(), timestamps.end());
        std::size_t cluster_length = 1;
        for (std::size_t index = 1; index < timestamps.size(); ++index) {
            if (timestamps[index] == timestamps[index - 1] + bp::kSampleLegalLineCount + 1) {
                ++cluster_length;
                continue;
            }
            lengths_by_gbt[gbt_index].push_back(cluster_length);
            cluster_length = 1;
        }
        lengths_by_gbt[gbt_index].push_back(cluster_length);
    }
    return lengths_by_gbt;
}

std::array<std::vector<std::size_t>, 256> dq_cluster_lengths_by_gbt(const ReconEvent& event)
{
    std::array<std::vector<std::size_t>, 256> lengths_by_gbt;
    std::array<std::vector<std::uint64_t>, 256> timestamps_by_gbt;
    for (const auto& line : event.dq_lines) {
        timestamps_by_gbt[line.header_vldb_id].push_back(bp::data_line_timestamp(line));
    }

    for (std::size_t gbt_index = 0; gbt_index < timestamps_by_gbt.size(); ++gbt_index) {
        auto& timestamps = timestamps_by_gbt[gbt_index];
        if (timestamps.empty()) {
            continue;
        }
        std::sort(timestamps.begin(), timestamps.end());
        std::size_t cluster_length = 1;
        for (std::size_t index = 1; index < timestamps.size(); ++index) {
            if (timestamps[index] == timestamps[index - 1] + 1) {
                ++cluster_length;
                continue;
            }
            lengths_by_gbt[gbt_index].push_back(cluster_length);
            cluster_length = 1;
        }
        lengths_by_gbt[gbt_index].push_back(cluster_length);
    }
    return lengths_by_gbt;
}

std::size_t line_count_for_type(const LineStats& counts, LineType type)
{
    switch (type) {
    case LineType::heartbeat:
        return counts.heartbeat;
    case LineType::trigger:
        return counts.trigger;
    case LineType::data_idle:
        return counts.data_idle;
    case LineType::data_daq:
        return counts.data_daq;
    case LineType::dummy:
        return counts.dummy;
    case LineType::other:
        return counts.other;
    case LineType::partial_tail:
        return counts.partial_tail;
    }
    return 0;
}

void append_lengths_by_gbt(std::array<std::vector<std::size_t>, 256>& destination,
    const std::array<std::vector<std::size_t>, 256>& source)
{
    for (std::size_t index = 0; index < destination.size(); ++index) {
        destination[index].insert(destination[index].end(), source[index].begin(), source[index].end());
    }
}

void update_recon_statistics_from_event(ReconStatistics& statistics, ReconEvent& event)
{
    build_event_samples(event);
    ++statistics.pt_events;
    statistics.total_dq_lines += event.dq_lines.size();
    statistics.total_samples += event.samples.size();
    for (const auto& sample : event.samples) {
        for (const auto& elink : sample.elinks) {
            if (elink.crc_ok) {
                ++statistics.total_sample_crc_ok_elinks;
            }
        }
        if (sample_all_elinks_crc_ok(sample)) {
            ++statistics.total_samples_all_elinks_crc_ok;
        }
    }

    const std::array<LineType, 7> types{
        LineType::heartbeat,
        LineType::trigger,
        LineType::data_idle,
        LineType::data_daq,
        LineType::dummy,
        LineType::other,
        LineType::partial_tail,
    };
    for (const auto type : types) {
        statistics.event_counts_by_type[static_cast<std::size_t>(type)].push_back(line_count_for_type(event.counts, type));
    }

    append_lengths_by_gbt(statistics.dq_cluster_lengths_by_gbt, dq_cluster_lengths_by_gbt(event));
    append_lengths_by_gbt(statistics.sample_cluster_lengths_by_gbt_all, sample_cluster_lengths_by_gbt(event, false));
    append_lengths_by_gbt(statistics.sample_cluster_lengths_by_gbt_crc, sample_cluster_lengths_by_gbt(event, true));
}

ReconStatistics scan_recon_statistics(const std::string& input_path, std::optional<std::size_t> max_rows, std::optional<std::uintmax_t> max_bytes)
{
    const auto start_time = std::chrono::steady_clock::now();

    std::vector<char> input_buffer(kInputBufferSize);
    std::ifstream input;
    input.rdbuf()->pubsetbuf(input_buffer.data(), static_cast<std::streamsize>(input_buffer.size()));
    input.open(input_path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("failed to open input file: " + input_path);
    }

    ReconStatistics statistics;
    statistics.input_path = input_path;
    statistics.output_dir = recon_output_dir(input_path);
    statistics.pdf_path = statistics.output_dir / "pt_line_type_distributions.pdf";
    statistics.sample_root_path = statistics.output_dir / "sample_events.root";

    std::optional<ReconEvent> current_event;
    RawLine row{};
    bool incomplete_read = false;
    while ((!max_rows.has_value() || statistics.total_lines < *max_rows) && (!max_bytes.has_value() || statistics.total_bytes_parsed + row.size() <= *max_bytes)) {
        if (!input.read(reinterpret_cast<char*>(row.data()), static_cast<std::streamsize>(row.size()))) {
            incomplete_read = true;
            break;
        }

        ++statistics.total_lines;
        statistics.total_bytes_parsed += row.size();
        const auto line_index = statistics.total_lines;
        const auto type = classify_line(row);
        increment_stats(statistics.stats, type);
        if (type == LineType::other) {
            continue;
        }

        if (type == LineType::trigger) {
            if (current_event.has_value()) {
                update_recon_statistics_from_event(statistics, *current_event);
            }
            current_event = ReconEvent{};
            current_event->trigger_line = line_index;
            fill_event_count(*current_event, type);
            continue;
        }

        if (type == LineType::heartbeat && (!max_rows.has_value() || statistics.total_lines < *max_rows) && (!max_bytes.has_value() || statistics.total_bytes_parsed + row.size() <= *max_bytes)) {
            RawLine tail{};
            if (input.read(reinterpret_cast<char*>(tail.data()), static_cast<std::streamsize>(tail.size()))) {
                ++statistics.total_lines;
                statistics.total_bytes_parsed += tail.size();
                increment_stats(statistics.stats, LineType::heartbeat);
                if (current_event.has_value()) {
                    fill_event_count(*current_event, LineType::heartbeat);
                    fill_event_count(*current_event, LineType::heartbeat);
                }
                continue;
            }
            incomplete_read = true;
        }

        if (current_event.has_value()) {
            fill_event_count(*current_event, type);
            if (type == LineType::data_daq) {
                current_event->dq_lines.push_back(parse_data_line_struct(row));
                current_event->dq_line_indices.push_back(line_index);
            }
        }
    }

    if (current_event.has_value()) {
        update_recon_statistics_from_event(statistics, *current_event);
    }

    const bool stopped_by_size = max_bytes.has_value() && statistics.total_bytes_parsed < *max_bytes && statistics.total_bytes_parsed + row.size() > *max_bytes;
    if (stopped_by_size) {
        const auto remaining_bytes = static_cast<std::size_t>(*max_bytes - statistics.total_bytes_parsed);
        std::vector<char> tail(remaining_bytes);
        if (remaining_bytes != 0) {
            input.read(tail.data(), static_cast<std::streamsize>(tail.size()));
            statistics.partial_tail_bytes = static_cast<std::size_t>(input.gcount());
        }
    } else if (incomplete_read) {
        statistics.partial_tail_bytes = static_cast<std::size_t>(input.gcount());
    }
    if (input.bad()) {
        throw std::runtime_error(read_failure_message(input_path, statistics.total_lines, statistics.total_bytes_parsed, input.gcount()));
    }
    if (statistics.partial_tail_bytes != 0) {
        increment_stats(statistics.stats, LineType::partial_tail);
    }

    const auto end_time = std::chrono::steady_clock::now();
    statistics.elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();
    return statistics;
}

std::array<std::vector<std::size_t>, 256> sample_indices_by_gbt(const ReconEvent& event, bool crc_only)
{
    std::array<std::vector<std::size_t>, 256> indices_by_gbt;
    for (std::size_t sample_index = 0; sample_index < event.samples.size(); ++sample_index) {
        const auto& sample = event.samples[sample_index];
        if (crc_only && !sample_all_elinks_crc_ok(sample)) {
            continue;
        }
        indices_by_gbt[sample.gbt_index].push_back(sample_index);
    }
    for (auto& indices : indices_by_gbt) {
        std::sort(indices.begin(), indices.end(), [&](std::size_t left, std::size_t right) {
            return event.samples[left].first_timestamp < event.samples[right].first_timestamp;
        });
    }
    return indices_by_gbt;
}

std::array<std::vector<std::vector<std::size_t>>, 256> sample_clusters_by_gbt(const ReconEvent& event, bool crc_only)
{
    std::array<std::vector<std::vector<std::size_t>>, 256> clusters_by_gbt;
    const auto indices_by_gbt = sample_indices_by_gbt(event, crc_only);
    for (std::size_t gbt_index = 0; gbt_index < indices_by_gbt.size(); ++gbt_index) {
        const auto& indices = indices_by_gbt[gbt_index];
        if (indices.empty()) {
            continue;
        }
        std::vector<std::size_t> cluster{indices.front()};
        for (std::size_t index = 1; index < indices.size(); ++index) {
            const auto previous_timestamp = event.samples[indices[index - 1]].first_timestamp;
            const auto timestamp = event.samples[indices[index]].first_timestamp;
            if (timestamp == previous_timestamp + bp::kSampleLegalLineCount + 1) {
                cluster.push_back(indices[index]);
                continue;
            }
            clusters_by_gbt[gbt_index].push_back(cluster);
            cluster.clear();
            cluster.push_back(indices[index]);
        }
        clusters_by_gbt[gbt_index].push_back(cluster);
    }
    return clusters_by_gbt;
}

std::size_t count_events_all_required_gbts_have_sample_cluster_length(const std::vector<ReconEvent>& events, bool crc_only,
    const std::array<bool, 256>& required_gbts, std::size_t target_length)
{
    std::size_t total = 0;
    for (const auto& event : events) {
        const auto lengths_by_gbt = sample_cluster_lengths_by_gbt(event, crc_only);
        bool all_required_gbts_match = true;
        for (std::size_t gbt_index = 0; gbt_index < required_gbts.size(); ++gbt_index) {
            if (!required_gbts[gbt_index]) {
                continue;
            }
            const auto& lengths = lengths_by_gbt[gbt_index];
            if (std::find(lengths.begin(), lengths.end(), target_length) == lengths.end()) {
                all_required_gbts_match = false;
                break;
            }
        }
        if (all_required_gbts_match) {
            ++total;
        }
    }
    return total;
}

std::string format_all_lpgbt_sample_peak_events(const std::vector<ReconEvent>& events, bool crc_only,
    const std::array<std::vector<std::size_t>, 256>& sample_cluster_lengths, std::size_t raw_pt_line_count)
{
    const auto required_gbts = active_gbts_from_cluster_lengths(sample_cluster_lengths_by_gbt(events, false));
    if (!has_active_gbt(required_gbts)) {
        return format_count_with_percentage(0, raw_pt_line_count) + " @ L=none";
    }
    const auto peak_length = peak_cluster_length(sample_cluster_lengths);
    if (!peak_length.has_value()) {
        return format_count_with_percentage(0, raw_pt_line_count) + " @ L=none";
    }
    const auto event_count = count_events_all_required_gbts_have_sample_cluster_length(events, crc_only, required_gbts, *peak_length);
    return format_count_with_percentage(event_count, raw_pt_line_count) + " @ L=" + format_count(*peak_length);
}

struct SampleEventTreeBranches {
    std::uint32_t event_index{};
    std::uint64_t trigger_line{};
    std::uint32_t peak_cluster_length{};
    bool crc_required{};
    std::uint32_t n_lpgbt{};
    std::vector<std::uint8_t> lpgbt_id;
    std::vector<std::uint32_t> n_samples;
    std::vector<std::uint32_t> payload_offset;
    std::vector<std::uint32_t> payload_count;
    std::vector<std::uint32_t> sample_offset;
    std::vector<std::uint64_t> sample_timestamp;
    std::vector<std::uint32_t> first_l1a_fcmd_line_delta;
    std::vector<std::uint32_t> l1a_fcmd_line_count;
    std::vector<std::uint8_t> tc;
    std::vector<std::uint8_t> tp;
    std::vector<std::uint16_t> val0;
    std::vector<std::uint16_t> val1;
    std::vector<std::uint16_t> val2;
    std::vector<std::uint32_t> header_word;
    std::vector<std::uint32_t> crc_word;
    std::vector<std::uint32_t> computed_crc;
    std::vector<std::uint8_t> crc_ok;

    void clear_vectors()
    {
        lpgbt_id.clear();
        n_samples.clear();
        payload_offset.clear();
        payload_count.clear();
        sample_offset.clear();
        sample_timestamp.clear();
        first_l1a_fcmd_line_delta.clear();
        l1a_fcmd_line_count.clear();
        tc.clear();
        tp.clear();
        val0.clear();
        val1.clear();
        val2.clear();
        header_word.clear();
        crc_word.clear();
        computed_crc.clear();
        crc_ok.clear();
    }
};

void create_sample_event_tree_branches(TTree& tree, SampleEventTreeBranches& branches)
{
    tree.Branch("event_index", &branches.event_index);
    tree.Branch("trigger_line", &branches.trigger_line);
    tree.Branch("peak_cluster_length", &branches.peak_cluster_length);
    tree.Branch("crc_required", &branches.crc_required);
    tree.Branch("n_lpgbt", &branches.n_lpgbt);
    tree.Branch("lpgbt_id", &branches.lpgbt_id);
    tree.Branch("n_samples", &branches.n_samples);
    tree.Branch("payload_offset", &branches.payload_offset);
    tree.Branch("payload_count", &branches.payload_count);
    tree.Branch("sample_offset", &branches.sample_offset);
    tree.Branch("sample_timestamp", &branches.sample_timestamp);
    tree.Branch("first_l1a_fcmd_line_delta", &branches.first_l1a_fcmd_line_delta);
    tree.Branch("l1a_fcmd_line_count", &branches.l1a_fcmd_line_count);
    tree.Branch("tc", &branches.tc);
    tree.Branch("tp", &branches.tp);
    tree.Branch("val0", &branches.val0);
    tree.Branch("val1", &branches.val1);
    tree.Branch("val2", &branches.val2);
    tree.Branch("header_word", &branches.header_word);
    tree.Branch("crc_word", &branches.crc_word);
    tree.Branch("computed_crc", &branches.computed_crc);
    tree.Branch("crc_ok", &branches.crc_ok);
}

std::pair<std::uint32_t, std::uint32_t> l1a_fcmd_stats_for_gbt(const ReconEvent& event, std::uint8_t gbt_index)
{
    std::uint32_t count = 0;
    std::uint32_t first_delta = std::numeric_limits<std::uint32_t>::max();
    for (std::size_t index = 0; index < event.dq_lines.size() && index < event.dq_line_indices.size(); ++index) {
        const auto& line = event.dq_lines[index];
        if (line.header_vldb_id != gbt_index || !data_line_has_l1a_fcmd(line)) {
            continue;
        }
        ++count;
        if (first_delta == std::numeric_limits<std::uint32_t>::max()) {
            first_delta = static_cast<std::uint32_t>(event.dq_line_indices[index] - event.trigger_line);
        }
    }
    return {first_delta, count};
}

void append_sample_to_tree_vectors(SampleEventTreeBranches& branches, const bp::Sample& sample)
{
    branches.sample_timestamp.push_back(sample.first_timestamp);
    for (const auto& elink : sample.elinks) {
        branches.header_word.push_back(elink.header.has_value() ? elink.header->raw_word : 0u);
        branches.crc_word.push_back(elink.crc_word);
        branches.computed_crc.push_back(elink.computed_crc);
        branches.crc_ok.push_back(static_cast<std::uint8_t>(elink.crc_ok ? 1 : 0));
    }
    for (std::size_t row_index = 0; row_index < bp::kSamplePayloadLineCount; ++row_index) {
        for (std::size_t elink_index = 0; elink_index < bp::kSampleElinkCount; ++elink_index) {
            const auto& elink = sample.elinks[elink_index];
            branches.tc.push_back(elink.tc[row_index]);
            branches.tp.push_back(elink.tp[row_index]);
            branches.val0.push_back(elink.val0[row_index]);
            branches.val1.push_back(elink.val1[row_index]);
            branches.val2.push_back(elink.val2[row_index]);
        }
    }
}

bool fill_sample_event_tree_entry(SampleEventTreeBranches& branches, const ReconEvent& event, std::size_t event_index,
    bool crc_only, const std::array<bool, 256>& required_gbts, std::size_t peak_length)
{
    const auto clusters_by_gbt = sample_clusters_by_gbt(event, crc_only);
    branches.clear_vectors();
    branches.event_index = static_cast<std::uint32_t>(event_index);
    branches.trigger_line = event.trigger_line;
    branches.peak_cluster_length = static_cast<std::uint32_t>(peak_length);
    branches.crc_required = crc_only;

    for (std::size_t gbt_index = 0; gbt_index < required_gbts.size(); ++gbt_index) {
        if (!required_gbts[gbt_index]) {
            continue;
        }
        const auto& clusters = clusters_by_gbt[gbt_index];
        const auto cluster_iter = std::find_if(clusters.begin(), clusters.end(), [&](const auto& cluster) {
            return cluster.size() == peak_length;
        });
        if (cluster_iter == clusters.end()) {
            return false;
        }

        branches.lpgbt_id.push_back(static_cast<std::uint8_t>(gbt_index));
        branches.n_samples.push_back(static_cast<std::uint32_t>(cluster_iter->size()));
        branches.sample_offset.push_back(static_cast<std::uint32_t>(branches.sample_timestamp.size()));
        branches.payload_offset.push_back(static_cast<std::uint32_t>(branches.tc.size()));
        branches.payload_count.push_back(static_cast<std::uint32_t>(cluster_iter->size() * bp::kSamplePayloadLineCount * bp::kSampleElinkCount));
        const auto [first_delta, fcmd_count] = l1a_fcmd_stats_for_gbt(event, static_cast<std::uint8_t>(gbt_index));
        branches.first_l1a_fcmd_line_delta.push_back(first_delta);
        branches.l1a_fcmd_line_count.push_back(fcmd_count);
        for (const auto sample_index : *cluster_iter) {
            append_sample_to_tree_vectors(branches, event.samples[sample_index]);
        }
    }
    branches.n_lpgbt = static_cast<std::uint32_t>(branches.lpgbt_id.size());
    return branches.n_lpgbt != 0;
}

std::size_t write_sample_events_root(const ReconSummary& summary, bool crc_only_sample_clusters)
{
    const auto sample_cluster_lengths = sample_cluster_lengths_by_gbt(summary.events, crc_only_sample_clusters);
    const auto peak_length = peak_cluster_length(sample_cluster_lengths);
    const auto required_gbts = active_gbts_from_cluster_lengths(sample_cluster_lengths_by_gbt(summary.events, false));
    TFile output(summary.sample_root_path.string().c_str(), "RECREATE");
    if (output.IsZombie()) {
        throw std::runtime_error("failed to create sample ROOT file: " + summary.sample_root_path.string());
    }
    TNamed("Rootifier_input_file", summary.input_path.c_str()).Write();
    TNamed("Rootifier_tree_layout", "payload index order: lpGBT-major, sample-major, row-major, elink-major").Write();

    TTree tree("sample_events", "Peak sample cluster events");
    SampleEventTreeBranches branches;
    create_sample_event_tree_branches(tree, branches);

    if (!peak_length.has_value() || !has_active_gbt(required_gbts)) {
        TNamed("Rootifier_sample_events_note", "No peak sample clusters available").Write();
        tree.Write();
        output.Write();
        output.Close();
        return 0;
    }

    std::size_t entries = 0;
    for (std::size_t event_index = 0; event_index < summary.events.size(); ++event_index) {
        if (!fill_sample_event_tree_entry(branches, summary.events[event_index], event_index, crc_only_sample_clusters, required_gbts, *peak_length)) {
            continue;
        }
        tree.Fill();
        ++entries;
    }
    tree.Write();
    output.Write();
    output.Close();
    return entries;
}

std::size_t write_sample_events_root(const ReconStatistics& statistics, bool crc_only_sample_clusters,
    std::optional<std::size_t> max_rows, std::optional<std::uintmax_t> max_bytes)
{
    const auto& sample_cluster_lengths = crc_only_sample_clusters
        ? statistics.sample_cluster_lengths_by_gbt_crc
        : statistics.sample_cluster_lengths_by_gbt_all;
    const auto peak_length = peak_cluster_length(sample_cluster_lengths);
    const auto required_gbts = active_gbts_from_cluster_lengths(statistics.sample_cluster_lengths_by_gbt_all);

    TFile output(statistics.sample_root_path.string().c_str(), "RECREATE");
    if (output.IsZombie()) {
        throw std::runtime_error("failed to create sample ROOT file: " + statistics.sample_root_path.string());
    }
    TNamed("Rootifier_input_file", statistics.input_path.c_str()).Write();
    TNamed("Rootifier_tree_layout", "payload index order: lpGBT-major, sample-major, row-major, elink-major").Write();

    TTree tree("sample_events", "Peak sample cluster events");
    SampleEventTreeBranches branches;
    create_sample_event_tree_branches(tree, branches);

    if (!peak_length.has_value() || !has_active_gbt(required_gbts)) {
        TNamed("Rootifier_sample_events_note", "No peak sample clusters available").Write();
        tree.Write();
        output.Write();
        output.Close();
        return 0;
    }

    std::vector<char> input_buffer(kInputBufferSize);
    std::ifstream input;
    input.rdbuf()->pubsetbuf(input_buffer.data(), static_cast<std::streamsize>(input_buffer.size()));
    input.open(statistics.input_path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("failed to open input file: " + statistics.input_path);
    }

    std::optional<ReconEvent> current_event;
    RawLine row{};
    std::size_t total_lines = 0;
    std::size_t total_bytes_parsed = 0;
    std::size_t event_index = 0;
    std::size_t entries = 0;

    const auto flush_current_event = [&]() {
        if (!current_event.has_value()) {
            return;
        }
        build_event_samples(*current_event);
        if (fill_sample_event_tree_entry(branches, *current_event, event_index, crc_only_sample_clusters, required_gbts, *peak_length)) {
            tree.Fill();
            ++entries;
        }
        ++event_index;
        current_event.reset();
    };

    while ((!max_rows.has_value() || total_lines < *max_rows) && (!max_bytes.has_value() || total_bytes_parsed + row.size() <= *max_bytes)) {
        if (!input.read(reinterpret_cast<char*>(row.data()), static_cast<std::streamsize>(row.size()))) {
            break;
        }

        ++total_lines;
        total_bytes_parsed += row.size();
        const auto line_index = total_lines;
        const auto type = classify_line(row);
        if (type == LineType::other) {
            continue;
        }

        if (type == LineType::trigger) {
            flush_current_event();
            current_event = ReconEvent{};
            current_event->trigger_line = line_index;
            fill_event_count(*current_event, type);
            continue;
        }

        if (type == LineType::heartbeat && (!max_rows.has_value() || total_lines < *max_rows) && (!max_bytes.has_value() || total_bytes_parsed + row.size() <= *max_bytes)) {
            RawLine tail{};
            if (input.read(reinterpret_cast<char*>(tail.data()), static_cast<std::streamsize>(tail.size()))) {
                ++total_lines;
                total_bytes_parsed += tail.size();
                if (current_event.has_value()) {
                    fill_event_count(*current_event, LineType::heartbeat);
                    fill_event_count(*current_event, LineType::heartbeat);
                }
                continue;
            }
        }

        if (current_event.has_value()) {
            fill_event_count(*current_event, type);
            if (type == LineType::data_daq) {
                current_event->dq_lines.push_back(parse_data_line_struct(row));
                current_event->dq_line_indices.push_back(line_index);
            }
        }
    }
    flush_current_event();

    if (input.bad()) {
        throw std::runtime_error(read_failure_message(statistics.input_path, total_lines, total_bytes_parsed, input.gcount()));
    }

    tree.Write();
    output.Write();
    output.Close();
    return entries;
}

void draw_recon_hist_page(TCanvas& canvas, const std::string& pdf_path, const std::vector<std::size_t>& values,
    const std::string& plot_title, const std::string& x_title, const std::string& y_title,
    const std::string& input_path, std::size_t pt_events, const std::string& run_arguments, const std::string& run_time, int color)
{
    if (values.empty()) {
        return;
    }

    const auto [min_value_iter, max_value_iter] = std::minmax_element(values.begin(), values.end());
    const auto min_value = *min_value_iter;
    const auto max_value = *max_value_iter;
    const auto total_count = std::accumulate(values.begin(), values.end(), std::size_t{0});
    const double mean_count = static_cast<double>(total_count) / static_cast<double>(values.size());
    double x_min = 0.0;
    double x_max = 2.0;
    if (min_value == max_value) {
        x_min = max_value == 0 ? 0.0 : static_cast<double>(max_value) - 1.0;
        x_max = static_cast<double>(max_value) + 2.0;
    } else {
        const auto span = static_cast<double>(max_value - min_value);
        const auto padding = std::max(1.0, span * 0.1);
        x_min = static_cast<double>(min_value) - padding;
        x_max = static_cast<double>(max_value) + padding;
    }
    const int bins = std::max(2, static_cast<int>(std::ceil(x_max - x_min)));
    const auto title = plot_title + ";" + x_title + ";" + y_title;
    TH1D hist(("h_" + plot_title).c_str(), title.c_str(), bins, x_min, x_max);
    for (const auto value : values) {
        hist.Fill(static_cast<double>(value));
    }
    hist.SetMinimum(0.5);
    hist.SetLineWidth(2);
    hist.SetFillColorAlpha(color, 0.35);
    hist.Draw("hist");

    std::map<std::size_t, std::size_t> count_frequencies;
    for (const auto value : values) {
        ++count_frequencies[value];
    }
    const auto peak_iter = std::max_element(count_frequencies.begin(), count_frequencies.end(), [](const auto& left, const auto& right) {
        return left.second < right.second;
    });
    const auto peak_bin_value = peak_iter->first;
    const auto peak_bin_count = peak_iter->second;
    const double peak_bin_percentage = 100.0 * static_cast<double>(peak_bin_count) / static_cast<double>(values.size());

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextAlign(13);
    latex.SetTextSize(0.042);
    latex.DrawLatex(0.115, 0.88, "#bf{FoCal TB2026 July}");
    latex.SetTextSize(0.032);
    latex.DrawLatex(0.115, 0.835, plot_title.c_str());
    latex.DrawLatex(0.115, 0.795, file_name_only(input_path).c_str());
    latex.DrawLatex(0.115, 0.755, run_arguments.c_str());
    latex.DrawLatex(0.115, 0.715, run_time.c_str());

    std::ostringstream stat_text;
    stat_text << "mean count " << std::fixed << std::setprecision(2) << mean_count;
    latex.SetTextAlign(33);
    latex.DrawLatex(0.88, 0.88, stat_text.str().c_str());
    latex.DrawLatex(0.88, 0.835, ("total PT " + format_count(pt_events)).c_str());
    std::ostringstream peak_text;
    peak_text << "peak bin " << peak_bin_value;
    latex.DrawLatex(0.88, 0.79, peak_text.str().c_str());
    std::ostringstream peak_fraction_text;
    peak_fraction_text << "peak fraction " << std::fixed << std::setprecision(2) << peak_bin_percentage << "%";
    latex.DrawLatex(0.88, 0.745, peak_fraction_text.str().c_str());
    canvas.Print(pdf_path.c_str());
}

void draw_cluster_hist_page(TCanvas& canvas, const std::string& pdf_path,
    const std::array<std::vector<std::size_t>, 256>& lengths_by_gbt,
    const std::string& plot_title, const std::string& x_title, const std::string& y_title,
    const std::string& hist_prefix,
    const std::string& input_path, std::size_t pt_events, const std::string& run_arguments, const std::string& run_time)
{
    std::vector<std::size_t> active_gbts;
    std::size_t min_value = std::numeric_limits<std::size_t>::max();
    std::size_t max_value = 0;
    std::size_t total_clusters = 0;
    std::size_t total_lengths = 0;
    std::map<std::size_t, std::size_t> count_frequencies;
    for (std::size_t gbt_index = 0; gbt_index < lengths_by_gbt.size(); ++gbt_index) {
        const auto& values = lengths_by_gbt[gbt_index];
        if (values.empty()) {
            continue;
        }
        active_gbts.push_back(gbt_index);
        for (const auto value : values) {
            min_value = std::min(min_value, value);
            max_value = std::max(max_value, value);
            total_lengths += value;
            ++total_clusters;
            ++count_frequencies[value];
        }
    }
    if (active_gbts.empty()) {
        return;
    }
    canvas.SetLogy(1);

    double x_min = 0.0;
    double x_max = 2.0;
    if (min_value == max_value) {
        x_min = max_value == 0 ? 0.0 : static_cast<double>(max_value) - 1.0;
        x_max = static_cast<double>(max_value) + 2.0;
    } else {
        const auto span = static_cast<double>(max_value - min_value);
        const auto padding = std::max(1.0, span * 0.1);
        x_min = static_cast<double>(min_value) - padding;
        x_max = static_cast<double>(max_value) + padding;
    }
    const int bins = std::max(2, static_cast<int>(std::ceil(x_max - x_min)));
    const std::array<int, 10> colors = {kAzure + 1, kOrange + 7, kGreen + 2, kViolet + 1, kRed + 1, kTeal + 3, kPink + 1, kBlue + 2, kSpring + 5, kMagenta + 1};

    std::vector<std::unique_ptr<TH1D>> histograms;
    double max_bin_content = 0.0;
    for (const auto gbt_index : active_gbts) {
        const auto title = plot_title + ";" + x_title + ";" + y_title;
        auto hist = std::make_unique<TH1D>((hist_prefix + "_gbt_" + std::to_string(gbt_index)).c_str(), title.c_str(), bins, x_min, x_max);
        for (const auto value : lengths_by_gbt[gbt_index]) {
            hist->Fill(static_cast<double>(value));
        }
        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
            max_bin_content = std::max(max_bin_content, hist->GetBinContent(bin));
        }
        const auto color = colors[histograms.size() % colors.size()];
        hist->SetLineColor(color);
        hist->SetLineWidth(2);
        hist->SetFillStyle(0);
        histograms.push_back(std::move(hist));
    }

    TLegend legend(0.74, 0.73, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.028);
    legend.SetTextAlign(32);
    legend.SetMargin(0.18);

    for (std::size_t index = 0; index < histograms.size(); ++index) {
        auto& hist = histograms[index];
        hist->SetMinimum(0.5);
        hist->SetMaximum(std::max(1.0, max_bin_content * 1.25));
        const auto draw_option = index == 0 ? "hist" : "hist same";
        hist->Draw(draw_option);
        legend.AddEntry(hist.get(), ("lpGBT " + std::to_string(active_gbts[index])).c_str(), "l");
    }
    legend.Draw();

    const double mean_count = static_cast<double>(total_lengths) / static_cast<double>(total_clusters);
    const auto peak_iter = std::max_element(count_frequencies.begin(), count_frequencies.end(), [](const auto& left, const auto& right) {
        return left.second < right.second;
    });
    const auto peak_bin_value = peak_iter->first;
    const auto peak_bin_count = peak_iter->second;
    const double peak_bin_percentage = 100.0 * static_cast<double>(peak_bin_count) / static_cast<double>(total_clusters);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextAlign(13);
    latex.SetTextSize(0.042);
    latex.DrawLatex(0.115, 0.88, "#bf{FoCal TB2026 July}");
    latex.SetTextSize(0.032);
    latex.DrawLatex(0.115, 0.835, plot_title.c_str());
    latex.DrawLatex(0.115, 0.795, file_name_only(input_path).c_str());
    latex.DrawLatex(0.115, 0.755, run_arguments.c_str());
    latex.DrawLatex(0.115, 0.715, run_time.c_str());

    std::ostringstream stat_text;
    stat_text << "mean count " << std::fixed << std::setprecision(2) << mean_count;
    latex.SetTextAlign(33);
    latex.DrawLatex(0.88, 0.69, stat_text.str().c_str());
    latex.DrawLatex(0.88, 0.645, ("total PT " + format_count(pt_events)).c_str());
    latex.DrawLatex(0.88, 0.60, ("active lpGBT " + format_count(active_gbts.size())).c_str());
    std::ostringstream peak_text;
    peak_text << "peak bin " << peak_bin_value;
    latex.DrawLatex(0.88, 0.555, peak_text.str().c_str());
    std::ostringstream peak_fraction_text;
    peak_fraction_text << "peak fraction " << std::fixed << std::setprecision(2) << peak_bin_percentage << "%";
    latex.DrawLatex(0.88, 0.51, peak_fraction_text.str().c_str());
    canvas.Print(pdf_path.c_str());
}

void write_recon_pdf(const ReconSummary& summary, const std::string& run_arguments, const std::string& run_time, bool crc_only_sample_clusters)
{
    if (summary.events.empty()) {
        throw std::runtime_error("reconstruction found no PT events; no PDF was written");
    }

    std::filesystem::create_directories(summary.output_dir);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas canvas("recon_canvas", "PT event line count distributions", 1100, 800);
    canvas.SetLogy();
    const auto pdf_path = summary.pdf_path.string();
    canvas.Print((pdf_path + "[").c_str());

    const std::array<LineType, 6> types{
        LineType::heartbeat,
        LineType::trigger,
        LineType::data_idle,
        LineType::data_daq,
        LineType::dummy,
        LineType::other,
    };

    for (const auto type : types) {
        const auto values = event_counts_for_type(summary.events, type);
        const auto plot_title = std::string(line_type_name(type)) + " per PT event";
        draw_recon_hist_page(canvas, pdf_path, values, plot_title, "line count per physical trigger", "physical triggers", summary.input_path, summary.events.size(), run_arguments, run_time, kAzure + static_cast<int>(type));
    }

    draw_cluster_hist_page(canvas, pdf_path, dq_cluster_lengths_by_gbt(summary.events),
        "DQ timestamp cluster length", "consecutive DQ lines per cluster", "DQ line cluster count", "h_dq_cluster",
        summary.input_path, summary.events.size(), run_arguments, run_time);
    draw_cluster_hist_page(canvas, pdf_path, sample_cluster_lengths_by_gbt(summary.events, crc_only_sample_clusters),
        crc_only_sample_clusters ? "Sample timestamp cluster length (CRC OK)" : "Sample timestamp cluster length", "consecutive 40-line samples per cluster", "sample cluster count", "h_sample_cluster",
        summary.input_path, summary.events.size(), run_arguments, run_time);
    canvas.Print((pdf_path + "]").c_str());
}

void write_recon_pdf(const ReconStatistics& statistics, const std::string& run_arguments, const std::string& run_time, bool crc_only_sample_clusters)
{
    if (statistics.pt_events == 0) {
        throw std::runtime_error("reconstruction found no PT events; no PDF was written");
    }

    std::filesystem::create_directories(statistics.output_dir);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas canvas("recon_canvas", "PT event line count distributions", 1100, 800);
    canvas.SetLogy();
    const auto pdf_path = statistics.pdf_path.string();
    canvas.Print((pdf_path + "[").c_str());

    const std::array<LineType, 6> types{
        LineType::heartbeat,
        LineType::trigger,
        LineType::data_idle,
        LineType::data_daq,
        LineType::dummy,
        LineType::other,
    };

    for (const auto type : types) {
        const auto& values = statistics.event_counts_by_type[static_cast<std::size_t>(type)];
        const auto plot_title = std::string(line_type_name(type)) + " per PT event";
        draw_recon_hist_page(canvas, pdf_path, values, plot_title, "line count per physical trigger", "physical triggers", statistics.input_path, statistics.pt_events, run_arguments, run_time, kAzure + static_cast<int>(type));
    }

    draw_cluster_hist_page(canvas, pdf_path, statistics.dq_cluster_lengths_by_gbt,
        "DQ timestamp cluster length", "consecutive DQ lines per cluster", "DQ line cluster count", "h_dq_cluster",
        statistics.input_path, statistics.pt_events, run_arguments, run_time);
    draw_cluster_hist_page(canvas, pdf_path, crc_only_sample_clusters ? statistics.sample_cluster_lengths_by_gbt_crc : statistics.sample_cluster_lengths_by_gbt_all,
        crc_only_sample_clusters ? "Sample timestamp cluster length (CRC OK)" : "Sample timestamp cluster length", "consecutive 40-line samples per cluster", "sample cluster count", "h_sample_cluster",
        statistics.input_path, statistics.pt_events, run_arguments, run_time);
    canvas.Print((pdf_path + "]").c_str());
}

void print_recon_summary(const ReconSummary& summary, bool crc_only_sample_clusters, std::size_t sample_root_entries)
{
    const auto total_samples = total_recon_samples(summary.events);
    const auto sample_cluster_lengths = sample_cluster_lengths_by_gbt(summary.events, crc_only_sample_clusters);
    print_box("Reconstruction summary", {
        {"File", file_name_only(summary.input_path)},
        {"PDF", summary.pdf_path.filename().string()},
        {"Sample ROOT", summary.sample_root_path.filename().string()},
        {"Sample ROOT entries", format_count(sample_root_entries)},
        {"Parsed lines", format_count(summary.total_lines) + " x " + std::to_string(kRowSize) + " bytes"},
        {"PT events", format_count(summary.events.size())},
        {"DQ DataLine", format_count(total_recon_dq_lines(summary.events))},
        {"DQ clusters", format_count(dq_cluster_lengths(summary.events).size())},
        {"Samples", format_count(total_samples)},
        {"Estimated event cache", format_byte_size(summary.cached_event_memory_bytes)},
        {"Event cache limit", summary.memory_limit_bytes.has_value() ? format_byte_size(*summary.memory_limit_bytes) : std::string("disabled")},
        {crc_only_sample_clusters ? "Sample clusters (CRC OK)" : "Sample clusters", format_count(total_cluster_count(sample_cluster_lengths))},
        {crc_only_sample_clusters ? "All lpGBT sample peak events (CRC OK)" : "All lpGBT sample peak events", format_all_lpgbt_sample_peak_events(summary.events, crc_only_sample_clusters, sample_cluster_lengths, summary.stats.trigger)},
        {"Sample CRC OK elinks", format_count(total_sample_crc_ok_elinks(summary.events))},
        {"40-line all-elink CRC OK", format_count_with_percentage(total_samples_all_elinks_crc_ok(summary.events), total_samples)},
    });
}

void print_recon_summary(const ReconStatistics& statistics, bool crc_only_sample_clusters, std::size_t sample_root_entries)
{
    const auto& sample_cluster_lengths = crc_only_sample_clusters
        ? statistics.sample_cluster_lengths_by_gbt_crc
        : statistics.sample_cluster_lengths_by_gbt_all;
    print_box("Reconstruction summary", {
        {"File", file_name_only(statistics.input_path)},
        {"PDF", statistics.pdf_path.filename().string()},
        {"Sample ROOT", statistics.sample_root_path.filename().string()},
        {"Sample ROOT entries", format_count(sample_root_entries)},
        {"Parsed lines", format_count(statistics.total_lines) + " x " + std::to_string(kRowSize) + " bytes"},
        {"PT events", format_count(statistics.pt_events)},
        {"DQ DataLine", format_count(statistics.total_dq_lines)},
        {"DQ clusters", format_count(total_cluster_count(statistics.dq_cluster_lengths_by_gbt))},
        {"Samples", format_count(statistics.total_samples)},
        {crc_only_sample_clusters ? "Sample clusters (CRC OK)" : "Sample clusters", format_count(total_cluster_count(sample_cluster_lengths))},
        {crc_only_sample_clusters ? "All lpGBT sample peak events (CRC OK)" : "All lpGBT sample peak events", format_count_with_percentage(sample_root_entries, statistics.stats.trigger)},
        {"Sample CRC OK elinks", format_count(statistics.total_sample_crc_ok_elinks)},
        {"40-line all-elink CRC OK", format_count_with_percentage(statistics.total_samples_all_elinks_crc_ok, statistics.total_samples)},
    });
}

ParseSummary parse_summary_from_recon(const ReconSummary& recon_summary)
{
    ParseSummary summary;
    summary.input_path = recon_summary.input_path;
    summary.total_lines = recon_summary.total_lines;
    summary.total_bytes_parsed = recon_summary.total_bytes_parsed;
    summary.partial_tail_bytes = recon_summary.partial_tail_bytes;
    summary.printed_lines = 0;
    summary.elapsed_sec = recon_summary.elapsed_sec;
    summary.stats = recon_summary.stats;
    return summary;
}

ParseSummary parse_summary_from_recon(const ReconStatistics& recon_statistics)
{
    ParseSummary summary;
    summary.input_path = recon_statistics.input_path;
    summary.total_lines = recon_statistics.total_lines;
    summary.total_bytes_parsed = recon_statistics.total_bytes_parsed;
    summary.partial_tail_bytes = recon_statistics.partial_tail_bytes;
    summary.printed_lines = 0;
    summary.elapsed_sec = recon_statistics.elapsed_sec;
    summary.stats = recon_statistics.stats;
    return summary;
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
        ("recon", "Run PT-event reconstruction and write ROOT histogram PDF under dump/001/<input filename>")
        ("crc", "For --recon sample clusters, only use samples where all four elinks pass CRC")
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

        if (parsed_opts.count("recon")) {
            const bool crc_only_sample_clusters = parsed_opts.count("crc") != 0;
            const auto recon_statistics = scan_recon_statistics(input_path, max_rows, max_bytes);
            write_recon_pdf(recon_statistics, join_arguments(argc, argv), current_time_minute(), crc_only_sample_clusters);
            const auto sample_root_entries = write_sample_events_root(recon_statistics, crc_only_sample_clusters, max_rows, max_bytes);
            print_summary(parse_summary_from_recon(recon_statistics));
            print_recon_summary(recon_statistics, crc_only_sample_clusters, sample_root_entries);
            spdlog::info("Wrote reconstruction PDF to {}", recon_statistics.pdf_path.string());
            spdlog::info("Wrote sample event ROOT file to {}", recon_statistics.sample_root_path.string());
            return 0;
        }

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