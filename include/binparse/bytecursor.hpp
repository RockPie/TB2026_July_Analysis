#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

namespace bp {

class ByteSpan {
public:
	ByteSpan() = default;
	ByteSpan(const std::byte* data, std::size_t size) : data_(data), size_(size) {}

	const std::byte* begin() const { return data_; }
	const std::byte* end() const { return data_ + size_; }
	const std::byte* data() const { return data_; }
	std::size_t size() const { return size_; }
	bool empty() const { return size_ == 0; }
	const std::byte& operator[](std::size_t index) const { return data_[index]; }

	ByteSpan first(std::size_t count) const { return ByteSpan(data_, count); }
	ByteSpan subspan(std::size_t offset) const { return ByteSpan(data_ + offset, size_ - offset); }
	ByteSpan subspan(std::size_t offset, std::size_t count) const { return ByteSpan(data_ + offset, count); }

private:
	const std::byte* data_ = nullptr;
	std::size_t size_ = 0;
};

struct ByteCursor {
	static constexpr std::size_t kLineSize = 32;

	static uint32_t read_u32_le(ByteSpan bytes, std::size_t offset)
	{
		return static_cast<uint32_t>(bytes[offset]) |
			(static_cast<uint32_t>(bytes[offset + 1]) << 8) |
			(static_cast<uint32_t>(bytes[offset + 2]) << 16) |
			(static_cast<uint32_t>(bytes[offset + 3]) << 24);
	}

	static std::array<uint32_t, 8> read_words_le(ByteSpan line)
	{
		std::array<uint32_t, 8> words{};
		for (std::size_t i = 0; i < words.size(); ++i) {
			words[i] = read_u32_le(line, i * sizeof(uint32_t));
		}
		return words;
	}
};

} // namespace bp