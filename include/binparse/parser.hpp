#pragma once

#include <binparse/bytecursor.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

namespace bp {

struct Packet {};
struct Heartbeat {};
struct RDH_L0 {};
struct RDH_L1 {};

struct DataLine {
	uint8_t header_vldb_id{};
	uint16_t bx_cnt{};
	uint32_t ob_cnt{};
	uint32_t data_word0{};
	uint32_t data_word1{};
	uint32_t data_word2{};
	uint32_t data_word3{};
	uint32_t data_word4{};
	uint32_t data_word5{};
};

struct TrgLine {
	uint16_t bx_cnt{};
	uint32_t ob_cnt{};
};

class StreamParser {
public:
	template <typename OnPacket, typename OnHeartbeat, typename OnSync, typename OnRdhL0,
			  typename OnRdhL1, typename OnDataLine, typename OnTrgLine>
	StreamParser(OnPacket&& on_packet,
				 OnHeartbeat&& on_heartbeat,
				 OnSync&& on_sync,
				 OnRdhL0&& on_rdh_l0,
				 OnRdhL1&& on_rdh_l1,
				 OnDataLine&& on_data_line,
				 OnTrgLine&& on_trg_line)
		: on_packet_(std::forward<OnPacket>(on_packet)),
		  on_heartbeat_(std::forward<OnHeartbeat>(on_heartbeat)),
		  on_sync_(std::forward<OnSync>(on_sync)),
		  on_rdh_l0_(std::forward<OnRdhL0>(on_rdh_l0)),
		  on_rdh_l1_(std::forward<OnRdhL1>(on_rdh_l1)),
		  on_data_line_(std::forward<OnDataLine>(on_data_line)),
		  on_trg_line_(std::forward<OnTrgLine>(on_trg_line))
	{
	}

	void feed(const std::vector<std::byte>& bytes)
	{
		feed(ByteSpan(bytes.data(), bytes.size()));
	}

	void feed(ByteSpan bytes)
	{
		for (std::size_t offset = 0; offset + ByteCursor::kLineSize <= bytes.size(); offset += ByteCursor::kLineSize) {
			parse_line(bytes.subspan(offset, ByteCursor::kLineSize));
		}
	}

private:
	static bool looks_like_data_payload(uint32_t word)
	{
		const uint32_t high_nibble = word >> 28;
		const uint32_t low_nibble = word & 0xFu;
		return word == 0xACCCCCCCu || word == 0x4B00004Bu ||
			(high_nibble == 0xFu && (low_nibble == 0x5u || low_nibble == 0x2u));
	}

	static uint16_t bx_from_words(const std::array<uint32_t, 8>& words)
	{
		return static_cast<uint16_t>(words[6] & 0x0FFFu);
	}

	static uint32_t orbit_from_words(const std::array<uint32_t, 8>& words)
	{
		return (words[6] >> 12) | ((words[7] & 0xFFFFFu) << 20);
	}

	static uint8_t vldb_from_words(const std::array<uint32_t, 8>& words)
	{
		return static_cast<uint8_t>((words[7] >> 24) & 0xFFu);
	}

	void parse_line(ByteSpan raw)
	{
		const auto words = ByteCursor::read_words_le(raw);

		if (words[0] == 0xBEEF4007u) {
			on_rdh_l0_(RDH_L0{}, raw);
			return;
		}

		if ((words[0] & 0xFFFFu) == 0x4883u || words[0] == 0x00004003u) {
			on_rdh_l1_(RDH_L1{}, raw);
			return;
		}

		if (words[0] == 0xACCCCCCCu && words[1] == 0xACCCCCCCu && words[2] == 0xACCCCCCCu) {
			on_heartbeat_(Heartbeat{});
			return;
		}

		if (looks_like_data_payload(words[0]) || looks_like_data_payload(words[1]) ||
			looks_like_data_payload(words[2]) || looks_like_data_payload(words[3]) ||
			looks_like_data_payload(words[4])) {
			DataLine line{};
			line.data_word0 = words[0];
			line.data_word1 = words[1];
			line.data_word2 = words[2];
			line.data_word3 = words[3];
			line.data_word4 = words[4];
			line.data_word5 = words[5];
			line.bx_cnt = bx_from_words(words);
			line.ob_cnt = orbit_from_words(words);
			line.header_vldb_id = vldb_from_words(words);
			on_data_line_(line, raw);
			return;
		}

		TrgLine line{};
		line.bx_cnt = bx_from_words(words);
		line.ob_cnt = orbit_from_words(words);
		on_trg_line_(line, raw);
	}

	std::function<void(const Packet&)> on_packet_;
	std::function<void(const Heartbeat&)> on_heartbeat_;
	std::function<void(ByteSpan)> on_sync_;
	std::function<void(const RDH_L0&, ByteSpan)> on_rdh_l0_;
	std::function<void(const RDH_L1&, ByteSpan)> on_rdh_l1_;
	std::function<void(const DataLine&, ByteSpan)> on_data_line_;
	std::function<void(const TrgLine&, ByteSpan)> on_trg_line_;
};

} // namespace bp