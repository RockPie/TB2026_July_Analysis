#pragma once

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

#include <binparse/bytecursor.hpp>

namespace bp {

struct TailOptions {
	int poll_ms = 1;
	std::size_t read_chunk = 1u << 20;
	int inactivity_timeout_ms = 10000;
};

template <typename OnChunk>
void tail_growing_file(const std::string& path, const TailOptions& opts, OnChunk&& on_chunk)
{
	(void)opts;
	std::ifstream input(path, std::ios::binary);
	if (!input) {
		return;
	}

	std::vector<std::byte> buffer(opts.read_chunk);
	while (input) {
		input.read(reinterpret_cast<char*>(buffer.data()), static_cast<std::streamsize>(buffer.size()));
		const auto count = input.gcount();
		if (count > 0) {
			on_chunk(ByteSpan(buffer.data(), static_cast<std::size_t>(count)));
		}
	}
}

} // namespace bp