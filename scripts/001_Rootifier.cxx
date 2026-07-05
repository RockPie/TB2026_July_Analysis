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
    stream,
    dummy,
    other,
    partial_tail,
};

struct LineStats {
    std::size_t heartbeat = 0;
    std::size_t trigger = 0;
    std::size_t stream = 0;
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

std::string_view line_type_name(LineType type)
{
    switch (type) {
    case LineType::heartbeat:
        return "heartbeat";
    case LineType::trigger:
        return "trigger";
    case LineType::stream:
        return "stream";
    case LineType::dummy:
        return "dummy";
    case LineType::other:
        return "other";
    case LineType::partial_tail:
        return "partial_tail";
    }
    return "unknown";
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
    case LineType::stream:
        ++stats.stream;
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
        return LineType::stream;
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

std::string format_hex(std::uint64_t value, std::size_t width)
{
    std::ostringstream out;
    out << "0x" << std::uppercase << std::hex << std::setw(static_cast<int>(width)) << std::setfill('0') << value;
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

std::string describe_line(const RawLine& row, LineType type, std::optional<std::uint64_t> last_trigger_timestamp)
{
    std::ostringstream out;
    switch (type) {
    case LineType::trigger: {
        const auto timestamp = extract_trigger_timestamp(row);
        out << "ts=" << format_hex(timestamp, 8);
        break;
    }
    case LineType::stream: {
        const auto stream_timestamp = extract_stream_timestamp(row);
        std::optional<std::int64_t> trigger_delta;
        if (last_trigger_timestamp.has_value()) {
            trigger_delta = static_cast<std::int64_t>(stream_timestamp) - static_cast<std::int64_t>(*last_trigger_timestamp);
        }
        const auto fmc = frame_word(row, 4);
        out << "lpGBT=" << std::setw(3) << static_cast<int>(row[1])
            << " ts=" << format_hex(stream_timestamp, kTimestampHexWidth)
            << " dBB=" << format_delta(trigger_delta)
            << " lanes=" << describe_data_lanes(row)
            << " FMC=" << describe_fmc_word(fmc)
            << " words=[" << describe_frame_words(row) << ']';
        break;
    }
    case LineType::heartbeat:
        out << "heartbeat continuation follows in raw stream if present";
        break;
    case LineType::dummy:
        out << "all-zero firmware dummy line";
        break;
    case LineType::other:
        out << "first_byte=" << format_hex(row[0], 2);
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

bool should_print_line(std::size_t line_index, const std::vector<LineRange>& selected_ranges, std::size_t preview_limit)
{
    if (selected_ranges.empty()) {
        return line_index <= preview_limit;
    }
    return std::any_of(selected_ranges.begin(), selected_ranges.end(), [line_index](const LineRange& range) {
        return range.first <= line_index && line_index <= range.second;
    });
}

void print_line(std::size_t line_index, const RawLine& row, LineType type, std::optional<std::uint64_t> last_trigger_timestamp)
{
    const auto detail = describe_line(row, type, last_trigger_timestamp);
    if (detail.empty()) {
        spdlog::info("L{:>10}  {:<9}  {}", line_index, line_type_name(type), format_bytes(row.data(), row.size()));
        return;
    }
    spdlog::info("L{:>10}  {:<9}  {}    || {}", line_index, line_type_name(type), format_bytes(row.data(), row.size()), detail);
}

ParseSummary scan_raw_file(const std::string& input_path, const std::vector<LineRange>& selected_ranges, std::size_t preview_limit, std::optional<std::size_t> max_rows)
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
    while (!max_rows.has_value() || summary.total_lines < *max_rows) {
        if (!input.read(reinterpret_cast<char*>(row.data()), static_cast<std::streamsize>(row.size()))) {
            break;
        }
        ++summary.total_lines;
        summary.total_bytes_parsed += row.size();
        const auto type = classify_line(row);
        increment_stats(summary.stats, type);

        if (should_print_line(summary.total_lines, selected_ranges, preview_limit)) {
            print_line(summary.total_lines, row, type, last_trigger_timestamp);
            ++summary.printed_lines;
        }

        if (type == LineType::trigger) {
            last_trigger_timestamp = extract_trigger_timestamp(row);
        }
    }

    if (input.bad()) {
        throw std::runtime_error("failed while reading input file: " + input_path);
    }

    if (!max_rows.has_value() || summary.total_lines < *max_rows) {
        summary.partial_tail_bytes = static_cast<std::size_t>(input.gcount());
        if (summary.partial_tail_bytes != 0) {
            increment_stats(summary.stats, LineType::partial_tail);
        }
    }

    std::ifstream size_input(input_path, std::ios::binary | std::ios::ate);
    if (size_input) {
        summary.file_size = static_cast<std::uintmax_t>(size_input.tellg());
    }

    const auto end_time = std::chrono::steady_clock::now();
    summary.elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();
    return summary;
}

void print_summary(const ParseSummary& summary)
{
    const double mib_per_sec = summary.elapsed_sec > 0.0
        ? static_cast<double>(summary.total_bytes_parsed) / summary.elapsed_sec / (1024.0 * 1024.0)
        : 0.0;

    spdlog::info("File: {}", summary.input_path);
    spdlog::info("Size: {} bytes", format_count(static_cast<std::size_t>(summary.file_size)));
    spdlog::info("Parsed lines: {} x {} bytes", format_count(summary.total_lines), kRowSize);
    spdlog::info("Speed: {:.3f} s, {:.1f} MiB/s", summary.elapsed_sec, mib_per_sec);
    spdlog::info("Printed lines: {}", format_count(summary.printed_lines));
    spdlog::info("");
    spdlog::info("Line type statistics:");
    spdlog::info("  heartbeat lines : {}", format_count(summary.stats.heartbeat));
    spdlog::info("  physical trigger BB lines : {}", format_count(summary.stats.trigger));
    spdlog::info("  data stream AC lines : {}", format_count(summary.stats.stream));
    spdlog::info("  firmware dummy 00 lines : {}", format_count(summary.stats.dummy));
    spdlog::info("  other lines : {}", format_count(summary.stats.other));
    if (summary.partial_tail_bytes != 0) {
        spdlog::info("  partial tail bytes : {}", format_count(summary.partial_tail_bytes));
    }
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
    TNamed("Rootifier_stream_lines", std::to_string(summary.stats.stream).c_str()).Write();
    TNamed("Rootifier_trigger_lines", std::to_string(summary.stats.trigger).c_str()).Write();
    TNamed("Rootifier_heartbeat_lines", std::to_string(summary.stats.heartbeat).c_str()).Write();
    TNamed("Rootifier_dummy_lines", std::to_string(summary.stats.dummy).c_str()).Write();
    TNamed("Rootifier_other_lines", std::to_string(summary.stats.other).c_str()).Write();
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
        ("l,lines", "1-based line indices/ranges to print, e.g. 1,2,40-48", cxxopts::value<std::string>()->default_value(""))
        ("p,preview", "Number of first lines to print when --lines is not set", cxxopts::value<std::size_t>()->default_value("40"))
        ("max-rows", "Stop after this many full 32-byte rows", cxxopts::value<std::size_t>())
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
        const auto line_ranges = parse_line_ranges(parsed_opts["lines"].as<std::string>());
        const auto preview_limit = parsed_opts["preview"].as<std::size_t>();
        const std::optional<std::size_t> max_rows = parsed_opts.count("max-rows")
            ? std::optional<std::size_t>(parsed_opts["max-rows"].as<std::size_t>())
            : std::nullopt;

        const auto summary = scan_raw_file(input_path, line_ranges, preview_limit, max_rows);
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