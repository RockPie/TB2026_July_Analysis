#pragma once

#include <binparse/parser.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace bp {

constexpr std::size_t kSampleLegalLineCount = 40;
constexpr std::size_t kSamplePayloadLineCount = 38;
constexpr std::size_t kSampleElinkCount = 4;
constexpr std::uint32_t kCrc32Mpeg2Polynomial = 0x04C11DB7;

struct SampleHeader {
	std::uint32_t raw_word{};
	std::uint8_t header{};
	std::uint16_t bx_counter{};
	std::uint8_t event_counter{};
	std::uint8_t orbit_counter{};
	std::array<std::uint8_t, 3> hamming_flags{};
	std::uint8_t trailer{};
	std::optional<bool> first_event;
};

struct SampleElink {
	std::optional<SampleHeader> header;
	std::uint32_t crc_word{};
	std::uint32_t computed_crc{};
	bool crc_ok = false;
	std::array<std::uint8_t, kSamplePayloadLineCount> tc{};
	std::array<std::uint8_t, kSamplePayloadLineCount> tp{};
	std::array<std::uint16_t, kSamplePayloadLineCount> val0{};
	std::array<std::uint16_t, kSamplePayloadLineCount> val1{};
	std::array<std::uint16_t, kSamplePayloadLineCount> val2{};
};

struct Sample {
	std::uint8_t gbt_index{};
	std::uint64_t first_timestamp{};
	std::array<SampleElink, kSampleElinkCount> elinks{};
};

inline std::uint64_t data_line_timestamp(const DataLine& line)
{
	return (static_cast<std::uint64_t>(line.ob_cnt) << 16) | line.bx_cnt;
}

inline std::uint32_t sample_word(std::uint32_t displayed_word)
{
	return ((displayed_word & 0x000000FFu) << 24) |
		((displayed_word & 0x0000FF00u) << 8) |
		((displayed_word & 0x00FF0000u) >> 8) |
		((displayed_word & 0xFF000000u) >> 24);
}

inline bool is_sample_header_word(std::uint32_t raw_word)
{
	return ((raw_word >> 28) & 0xFu) == 0xFu && ((raw_word & 0xFu) == 0x2u || (raw_word & 0xFu) == 0x5u);
}

inline SampleHeader parse_sample_header(std::uint32_t raw_word)
{
	const auto trailer = static_cast<std::uint8_t>(raw_word & 0xFu);
	SampleHeader header{};
	header.raw_word = raw_word;
	header.header = static_cast<std::uint8_t>((raw_word >> 28) & 0xFu);
	header.bx_counter = static_cast<std::uint16_t>((raw_word >> 16) & 0x0FFFu);
	header.event_counter = static_cast<std::uint8_t>((raw_word >> 10) & 0x3Fu);
	header.orbit_counter = static_cast<std::uint8_t>((raw_word >> 7) & 0x7u);
	header.hamming_flags = {
		static_cast<std::uint8_t>((raw_word >> 6) & 0x1u),
		static_cast<std::uint8_t>((raw_word >> 5) & 0x1u),
		static_cast<std::uint8_t>((raw_word >> 4) & 0x1u),
	};
	header.trailer = trailer;
	if (trailer == 0x2u) {
		header.first_event = true;
	} else if (trailer == 0x5u) {
		header.first_event = false;
	}
	return header;
}

inline std::uint32_t crc32_mpeg2_word_update(std::uint32_t crc_value, std::uint32_t word)
{
	for (int byte_index = 3; byte_index >= 0; --byte_index) {
		const auto byte_value = static_cast<std::uint8_t>((word >> (8 * byte_index)) & 0xFFu);
		crc_value ^= static_cast<std::uint32_t>(byte_value) << 24;
		for (int bit = 0; bit < 8; ++bit) {
			crc_value = (crc_value & 0x80000000u) != 0
				? ((crc_value << 1) ^ kCrc32Mpeg2Polynomial) & 0xFFFFFFFFu
				: (crc_value << 1) & 0xFFFFFFFFu;
		}
	}
	return crc_value;
}

inline std::uint32_t crc32_mpeg2_words(const std::array<std::uint32_t, kSampleLegalLineCount>& words)
{
	std::uint32_t crc_value = 0;
	for (std::size_t index = 0; index + 1 < words.size(); ++index) {
		crc_value = crc32_mpeg2_word_update(crc_value, words[index]);
	}
	return crc_value;
}

inline std::uint32_t data_word_at(const DataLine& line, std::size_t elink_index)
{
	switch (elink_index) {
	case 0:
		return line.data_word0;
	case 1:
		return line.data_word1;
	case 2:
		return line.data_word2;
	case 3:
		return line.data_word3;
	default:
		return 0;
	}
}

inline void fill_sample_word(SampleElink& elink, std::size_t sample_index, std::uint32_t word)
{
	elink.tc[sample_index] = static_cast<std::uint8_t>((word >> 31) & 0x1u);
	elink.tp[sample_index] = static_cast<std::uint8_t>((word >> 30) & 0x1u);
	elink.val0[sample_index] = static_cast<std::uint16_t>((word >> 20) & 0x03FFu);
	elink.val1[sample_index] = static_cast<std::uint16_t>((word >> 10) & 0x03FFu);
	elink.val2[sample_index] = static_cast<std::uint16_t>(word & 0x03FFu);
}

inline std::optional<Sample> make_sample(const std::vector<DataLine>& lines, std::size_t first_line)
{
	if (first_line + kSampleLegalLineCount > lines.size()) {
		return std::nullopt;
	}

	Sample sample{};
	sample.gbt_index = lines[first_line].header_vldb_id;
	sample.first_timestamp = data_line_timestamp(lines[first_line]);

	for (std::size_t line_offset = 0; line_offset < kSampleLegalLineCount; ++line_offset) {
		const auto& line = lines[first_line + line_offset];
		if (line.header_vldb_id != sample.gbt_index || data_line_timestamp(line) != sample.first_timestamp + line_offset) {
			return std::nullopt;
		}
	}

	for (std::size_t elink_index = 0; elink_index < kSampleElinkCount; ++elink_index) {
		std::array<std::uint32_t, kSampleLegalLineCount> words{};
		for (std::size_t line_offset = 0; line_offset < kSampleLegalLineCount; ++line_offset) {
			words[line_offset] = sample_word(data_word_at(lines[first_line + line_offset], elink_index));
		}
		if (!is_sample_header_word(words.front())) {
			return std::nullopt;
		}
		sample.elinks[elink_index].header = parse_sample_header(words.front());
		sample.elinks[elink_index].crc_word = words.back();
		sample.elinks[elink_index].computed_crc = crc32_mpeg2_words(words);
		sample.elinks[elink_index].crc_ok = sample.elinks[elink_index].computed_crc == sample.elinks[elink_index].crc_word;
	}

	for (std::size_t payload_index = 0; payload_index < kSamplePayloadLineCount; ++payload_index) {
		const auto& line = lines[first_line + payload_index + 1];
		for (std::size_t elink_index = 0; elink_index < kSampleElinkCount; ++elink_index) {
			fill_sample_word(sample.elinks[elink_index], payload_index, sample_word(data_word_at(line, elink_index)));
		}
	}

	return sample;
}

} // namespace bp
