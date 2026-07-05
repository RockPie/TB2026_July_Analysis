#!/usr/bin/env python3
from __future__ import annotations

import argparse
import signal
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence


COMMON_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = COMMON_DIR / "temp"
ROW_SIZE = 32
WORD_SIZE = 4
FRAME_WORD_OFFSET = 8
FRAME_WORD_COUNT = 6
SAMPLE_TIMESTAMP_WINDOW = 40
SAMPLE_WORD_COUNT = 40
FRAME_EVENT_PREVIEW_COUNT = 64
FRAME_SAMPLE_PREVIEW_INTERVALS = 2
FRAME_SAMPLE_PREVIEW_LINKS = {0, 1}
FRAME_FULL_PREVIEW_LINKS = {0, 1}
DATA_FRAME_WORD_COUNT = 40
DATA_FRAME_PAYLOAD_WORD_COUNT = 39
CRC32_POLYNOMIAL = 0x04C11DB7
ZERO_ROW = b"\x00" * ROW_SIZE
ZERO_WORD = b"\x00" * WORD_SIZE
FRAME_IDLE_WORD = b"\xcc\xcc\xcc\xac"
FMC_IDLE_WORD = b"\x36\x36\x36\x36"
FMC_INJECTION_WORD = b"\x2d\x2d\x2d\x2d"
FMC_L1A_WORD = b"\x4b\x4b\x4b\x4b"

TYPE_HEARTBEAT = "heartbeat"
TYPE_TRIGGER = "trigger"
TYPE_STREAM = "stream"
TYPE_DUMMY = "dummy"
TYPE_OTHER = "other"
TYPE_PARTIAL = "partial_tail"
TYPE_IDLE_AC_SUMMARY = "idle_ac_summary"
ANSI_GRAY = "\033[90m"
ANSI_BLUE = "\033[1;34m"
ANSI_INDEX_COLORS = ("\033[1;36m", "\033[1;33m", "\033[1;35m", "\033[1;32m", "\033[1;31m", "\033[1;37m")
ANSI_GROUP_COLORS = ("\033[38;5;208m", "\033[38;5;45m", "\033[38;5;220m", "\033[38;5;141m", "\033[38;5;82m", "\033[38;5;203m")
ANSI_RESET = "\033[0m"
RICH_INDEX_COLORS = ("cyan", "yellow", "magenta", "green", "red", "white")
RICH_GROUP_COLORS = ("orange1", "deep_sky_blue1", "gold1", "medium_purple1", "spring_green2", "indian_red1")
TIMESTAMP_HEX_WIDTH = 12


@dataclass
class StreamPreview:
	line_no: int
	row_type: str
	lpgbt_index: int | None
	timestamp: int | None
	trigger_delta: int | None
	hex_text: str
	parsed_text: str = ""


@dataclass(frozen=True)
class FmcInjectionSummary:
	interval_no: int
	lpgbt_index: int
	injection_no: int
	idle_before_l1a: int
	l1a_count: int


@dataclass
class FrameWordPreview:
	interval_no: int
	line_no: int
	link_id: int
	byte_offset: int
	row_type: str
	lpgbt_index: int | None
	timestamp: int | None
	bb_delta: int | None
	word: bytes
	reordered_word: bytes


@dataclass
class OpenFrameSample:
	interval_no: int
	link_id: int
	header_line_no: int
	header_timestamp: int
	header_word: bytes
	words: list[bytes] = field(default_factory=list)
	last_line_no: int = 0
	last_timestamp: int = 0


@dataclass
class CompleteFrameSample:
	interval_no: int
	link_id: int
	header_line_no: int
	header_timestamp: int
	header_word: bytes
	last_line_no: int
	last_timestamp: int
	words: list[bytes]


@dataclass
class DataFrameHeader:
	raw_word: int
	header: int
	bx_counter: int
	event_counter: int
	orbit_counter: int
	hamming_flags: tuple[int, int, int]
	trailer: int
	first_event: bool | None


@dataclass
class DataFrameChannel:
	name: str
	word_index: int
	raw_word: int
	tc: int
	tp: int
	adc: int
	tot: int
	toa: int


@dataclass
class DataFrameParseResult:
	ok: bool
	errors: list[str]
	words: tuple[int, ...]
	header: DataFrameHeader | None = None
	channels: list[DataFrameChannel] = field(default_factory=list)
	crc_word: int | None = None
	computed_crc: int | None = None
	crc_ok: bool = False


@dataclass
class FrameIntervalSummary:
	interval_no: int
	sample_counts: tuple[int, ...]


@dataclass
class FrameEventPreview:
	interval_no: int
	event_index: int
	link_words: tuple[tuple[bytes, ...], ...]


@dataclass
class FrameCrcErrorSample:
	interval_no: int
	event_index: int
	link_id: int
	parsed: DataFrameParseResult
	previews: list[FrameWordPreview]


@dataclass
class FrameParseResult:
	path: Path
	file_size: int
	elapsed_sec: float
	total_lines: int = 0
	total_bytes_parsed: int = 0
	partial_tail_bytes: int = 0
	trigger_lines: int = 0
	trigger_intervals: int = 0
	frame_words: int = 0
	idle_words: int = 0
	zero_words: int = 0
	heartbeat_payload_words: int = 0
	pre_trigger_words: int = 0
	trailing_open_words: int = 0
	header_words: int = 0
	complete_samples: int = 0
	incomplete_samples: int = 0
	expired_samples: int = 0
	preview: list[FrameWordPreview] = field(default_factory=list)
	interval_summaries: list[FrameIntervalSummary] = field(default_factory=list)
	event_previews: list[FrameEventPreview] = field(default_factory=list)
	crc_error_samples: list[FrameCrcErrorSample] = field(default_factory=list)
	sample_preview: list[CompleteFrameSample] = field(default_factory=list)

	@property
	def mib_per_sec(self) -> float:
		if self.elapsed_sec <= 0:
			return 0.0
		return self.total_bytes_parsed / self.elapsed_sec / (1024 * 1024)


@dataclass
class ParseResult:
	path: Path
	file_size: int
	elapsed_sec: float
	total_lines: int = 0
	total_bytes_parsed: int = 0
	partial_tail_bytes: int = 0
	counts: Counter[str] = field(default_factory=Counter)
	heartbeat_frames: int = 0
	heartbeat_unpaired_lines: int = 0
	lpgbt_counts: Counter[int] = field(default_factory=Counter)
	lpgbt_all_zero_counts: Counter[int] = field(default_factory=Counter)
	lpgbt_all_idle_counts: Counter[int] = field(default_factory=Counter)
	trigger_interval_count: int = 0
	trigger_interval_lpgbt_totals: Counter[int] = field(default_factory=Counter)
	trigger_interval_dummy_total: int = 0
	other_first_byte_counts: Counter[int] = field(default_factory=Counter)
	bb_first_l1a_lpgbt_line_deltas: dict[int, Counter[int]] = field(default_factory=dict)
	bb_first_l1a_lpgbt_missing_intervals: Counter[int] = field(default_factory=Counter)
	bb_first_l1a_lpgbt_pair_deltas: Counter[tuple[int | None, int | None]] = field(default_factory=Counter)
	bb_l1a_lpgbt_line_counts: dict[int, Counter[int]] = field(default_factory=dict)
	fmc_summaries: list[FmcInjectionSummary] = field(default_factory=list)
	fmc_pending_injections: int = 0
	stream_preview: list[StreamPreview] = field(default_factory=list)

	@property
	def mib_per_sec(self) -> float:
		if self.elapsed_sec <= 0:
			return 0.0
		return self.total_bytes_parsed / self.elapsed_sec / (1024 * 1024)


def latest_nonempty_ch2g(data_dir: Path = DEFAULT_DATA_DIR) -> Path:
	files = [path for path in data_dir.glob("data*.ch2g") if path.is_file() and path.stat().st_size > 0]
	if not files:
		raise FileNotFoundError(f"No non-empty data*.ch2g files found in {data_dir}")
	return max(files, key=lambda path: path.stat().st_mtime)


def classify_row(row: bytes) -> str:
	if row == ZERO_ROW:
		return TYPE_DUMMY
	first = row[0]
	if first == 0x07:
		return TYPE_HEARTBEAT
	if first == 0xBB:
		return TYPE_TRIGGER
	if first == 0xAC:
		return TYPE_STREAM
	return TYPE_OTHER


def extract_stream_timestamp(row: bytes) -> int:
	return int.from_bytes(row[2:8], "little")


def extract_trigger_timestamp(row: bytes) -> int:
	return row[4] | (row[5] << 8) | (row[12] << 16) | (row[13] << 24)


def reorder_frame_word(word: bytes) -> bytes:
	return word[::-1]


def is_sample_header_word(reordered_word: bytes) -> bool:
	if reordered_word == reorder_frame_word(FRAME_IDLE_WORD):
		return False
	return (reordered_word[0] & 0xF0) == 0xF0 and (reordered_word[-1] & 0x0F) in {0x2, 0x5}


def crc32_mpeg2(data: bytes) -> int:
	crc_value = 0
	for byte_value in data:
		crc_value ^= byte_value << 24
		for _ in range(8):
			if crc_value & 0x80000000:
				crc_value = ((crc_value << 1) ^ CRC32_POLYNOMIAL) & 0xFFFFFFFF
			else:
				crc_value = (crc_value << 1) & 0xFFFFFFFF
	return crc_value


def data_frame_channel_name(word_index: int) -> str:
	if word_index == 1:
		return "CM"
	if 2 <= word_index <= 19:
		return f"ch{word_index - 2}"
	if word_index == 20:
		return "Calib"
	if 21 <= word_index <= 38:
		return f"ch{word_index - 3}"
	return f"word{word_index}"


def normalize_data_frame_words(words: Sequence[int | bytes]) -> tuple[list[int], list[str]]:
	normalized_words: list[int] = []
	errors: list[str] = []
	for index, word in enumerate(words):
		if isinstance(word, int):
			if word < 0 or word > 0xFFFFFFFF:
				errors.append(f"word {index} is outside 32-bit range: {word}")
				continue
			normalized_words.append(word)
		elif isinstance(word, bytes):
			if len(word) != WORD_SIZE:
				errors.append(f"word {index} has {len(word)} bytes, expected {WORD_SIZE}")
				continue
			normalized_words.append(int.from_bytes(word, "big"))
		else:
			errors.append(f"word {index} has unsupported type {type(word).__name__}")
	return normalized_words, errors


def parse_data_frame_header(raw_word: int) -> DataFrameHeader:
	trailer = raw_word & 0xF
	first_event = True if trailer == 0x2 else False if trailer == 0x5 else None
	return DataFrameHeader(
		raw_word=raw_word,
		header=(raw_word >> 28) & 0xF,
		bx_counter=(raw_word >> 16) & 0xFFF,
		event_counter=(raw_word >> 10) & 0x3F,
		orbit_counter=(raw_word >> 7) & 0x7,
		hamming_flags=((raw_word >> 6) & 0x1, (raw_word >> 5) & 0x1, (raw_word >> 4) & 0x1),
		trailer=trailer,
		first_event=first_event,
	)


def parse_data_frame_channel(word_index: int, raw_word: int) -> DataFrameChannel:
	return DataFrameChannel(
		name=data_frame_channel_name(word_index),
		word_index=word_index,
		raw_word=raw_word,
		tc=(raw_word >> 31) & 0x1,
		tp=(raw_word >> 30) & 0x1,
		adc=(raw_word >> 20) & 0x3FF,
		tot=(raw_word >> 10) & 0x3FF,
		toa=raw_word & 0x3FF,
	)


def parse_h2gcroc_data_frame(words: Sequence[int | bytes]) -> DataFrameParseResult:
	normalized_words, errors = normalize_data_frame_words(words)
	if len(normalized_words) != DATA_FRAME_WORD_COUNT:
		errors.append(f"expected {DATA_FRAME_WORD_COUNT} data words, got {len(normalized_words)}")

	result = DataFrameParseResult(ok=False, errors=errors, words=tuple(normalized_words))
	if len(normalized_words) < DATA_FRAME_WORD_COUNT:
		return result

	frame_words = normalized_words[:DATA_FRAME_WORD_COUNT]
	result.words = tuple(frame_words)
	result.header = parse_data_frame_header(frame_words[0])
	if result.header.header != 0xF:
		result.errors.append(f"invalid DAQ header nibble 0x{result.header.header:X}, expected 0xF")
	if result.header.first_event is None:
		result.errors.append(f"unknown trigger flag 0x{result.header.trailer:X}, expected 0x2 or 0x5")

	result.channels = [parse_data_frame_channel(word_index, frame_words[word_index]) for word_index in range(1, DATA_FRAME_PAYLOAD_WORD_COUNT)]
	result.crc_word = frame_words[-1]
	payload = b"".join(word.to_bytes(WORD_SIZE, "big") for word in frame_words[:DATA_FRAME_PAYLOAD_WORD_COUNT])
	result.computed_crc = crc32_mpeg2(payload)
	result.crc_ok = result.computed_crc == result.crc_word
	if not result.crc_ok:
		result.errors.append(f"CRC mismatch: expected 0x{result.crc_word:08X}, computed 0x{result.computed_crc:08X}")
	result.ok = not result.errors
	return result


def format_hex_value(value: int | None) -> str:
	if value is None:
		return "-"
	width = max(1, ((value.bit_length() + 7) // 8) * 2)
	return f"0x{value:0{width}X}"


def format_delta_value(value: int | None) -> str:
	if value is None:
		return "     -"
	return f"0x{value:04X}"


def format_timestamp(value: int | None) -> str:
	if value is None:
		return "-"
	return f"0x{value:0{TIMESTAMP_HEX_WIDTH}X}"


def format_group_summary(count: int) -> str:
	full_groups, remainder = divmod(count, SAMPLE_WORD_COUNT)
	if remainder:
		return f"{full_groups}*{SAMPLE_WORD_COUNT}+{remainder}"
	return f"{full_groups}*{SAMPLE_WORD_COUNT}"


def parse_ch2g(path: Path, preview_limit: int = 40, chunk_rows: int = 262144, max_rows: int | None = None, max_bb_intervals: int | None = None, preview_idle: bool = False) -> ParseResult:
	start_time = time.perf_counter()
	file_size = path.stat().st_size
	result = ParseResult(path=path, file_size=file_size, elapsed_sec=0.0)
	chunk_size = ROW_SIZE * chunk_rows
	heartbeat_continuation_pending = False
	seen_trigger = False
	last_trigger_timestamp: int | None = None
	current_trigger_interval_lpgbt_counts: Counter[int] = Counter()
	current_trigger_interval_dummy_count = 0
	current_fmc_states: dict[int, dict[str, int | bool]] = {}
	current_interval_no = 0
	current_first_l1a_lpgbt_line_deltas: dict[int, int] = {}
	current_l1a_lpgbt_line_counts: Counter[int] = Counter()
	current_word_counts_by_lpgbt: dict[int, list[int]] = {}
	line_no = 0
	stop_parse = False
	idle_ac_preview_counts: Counter[str] = Counter()

	def preview_has_space() -> bool:
		return len(result.stream_preview) < preview_limit

	def flush_idle_ac_preview() -> None:
		nonlocal idle_ac_preview_counts
		total_count = sum(idle_ac_preview_counts.values())
		if total_count > 1 and preview_has_space():
			line_word = "line" if total_count == 1 else "lines"
			details = ", ".join(
				f"{name}={format_count(count)}"
				for name, count in (("all0", idle_ac_preview_counts["all0"]), ("allAC", idle_ac_preview_counts["allAC"]), ("mixed", idle_ac_preview_counts["mixed"]))
				if count
			)
			result.stream_preview.append(StreamPreview(line_no, TYPE_IDLE_AC_SUMMARY, None, None, None, "", f"{format_count(total_count)} idle AC {line_word} ({details})"))
		idle_ac_preview_counts = Counter()

	def reset_fmc_states() -> None:
		nonlocal current_fmc_states
		current_fmc_states = {}

	def get_fmc_state(lpgbt_index: int) -> dict[str, int | bool]:
		return current_fmc_states.setdefault(lpgbt_index, {"pending": False, "idle": 0, "l1a": 0, "injection_no": 0})

	def flush_fmc_state(lpgbt_index: int, state: dict[str, int | bool]) -> None:
		if not state.get("pending"):
			return
		result.fmc_summaries.append(FmcInjectionSummary(
			interval_no=current_interval_no,
			lpgbt_index=lpgbt_index,
			injection_no=int(state["injection_no"]),
			idle_before_l1a=int(state["idle"]),
			l1a_count=int(state["l1a"]),
		))
		state["pending"] = False
		state["idle"] = 0
		state["l1a"] = 0

	def flush_fmc_states() -> None:
		for lpgbt_index, state in sorted(current_fmc_states.items()):
			flush_fmc_state(lpgbt_index, state)

	def update_fmc_state(lpgbt_index: int, row: bytes) -> None:
		fmc_word_start = FRAME_WORD_OFFSET + 4 * WORD_SIZE
		fmc_word = row[fmc_word_start:fmc_word_start + WORD_SIZE]
		state = get_fmc_state(lpgbt_index)
		if fmc_word == FMC_INJECTION_WORD:
			flush_fmc_state(lpgbt_index, state)
			state["pending"] = True
			state["idle"] = 0
			state["l1a"] = 0
			state["injection_no"] = int(state["injection_no"]) + 1
			result.fmc_pending_injections += 1
		elif fmc_word == FMC_IDLE_WORD and state.get("pending") and int(state["l1a"]) == 0:
			state["idle"] = int(state["idle"]) + 1
		elif fmc_word == FMC_L1A_WORD and state.get("pending"):
			state["l1a"] = int(state["l1a"]) + 1

	def has_fmc_l1a_word(row: bytes) -> bool:
		fmc_word_start = FRAME_WORD_OFFSET + 4 * WORD_SIZE
		return row[fmc_word_start:fmc_word_start + WORD_SIZE] == FMC_L1A_WORD

	def finish_trigger_interval() -> None:
		for lpgbt_index in sorted(current_trigger_interval_lpgbt_counts):
			line_delta = current_first_l1a_lpgbt_line_deltas.get(lpgbt_index)
			if line_delta is None:
				result.bb_first_l1a_lpgbt_missing_intervals[lpgbt_index] += 1
			else:
				result.bb_first_l1a_lpgbt_line_deltas.setdefault(lpgbt_index, Counter())[line_delta] += 1
			result.bb_l1a_lpgbt_line_counts.setdefault(lpgbt_index, Counter())[current_l1a_lpgbt_line_counts[lpgbt_index]] += 1
		if 0 in current_trigger_interval_lpgbt_counts or 1 in current_trigger_interval_lpgbt_counts:
			result.bb_first_l1a_lpgbt_pair_deltas[(
				current_first_l1a_lpgbt_line_deltas.get(0),
				current_first_l1a_lpgbt_line_deltas.get(1),
			)] += 1

	def is_idle_stream_row(row: bytes) -> bool:
		if row[0] != 0xAC:
			return False
		fmc_word_start = FRAME_WORD_OFFSET + 4 * WORD_SIZE
		fmc_word = row[fmc_word_start:fmc_word_start + WORD_SIZE]
		if fmc_word == FMC_IDLE_WORD:
			state = get_fmc_state(row[1])
			if state.get("pending") and int(state["l1a"]) == 0:
				return False
		return fmc_word in {FMC_IDLE_WORD, ZERO_WORD} and not data_lanes_have_activity(row)

	def data_lanes_have_activity(row: bytes) -> bool:
		return any(
			row[FRAME_WORD_OFFSET + frame_word_index * WORD_SIZE:FRAME_WORD_OFFSET + (frame_word_index + 1) * WORD_SIZE] not in {ZERO_WORD, FRAME_IDLE_WORD}
			for frame_word_index in range(4)
		)

	def data_lanes_all_zero(row: bytes) -> bool:
		return all(
			row[FRAME_WORD_OFFSET + frame_word_index * WORD_SIZE:FRAME_WORD_OFFSET + (frame_word_index + 1) * WORD_SIZE] == ZERO_WORD
			for frame_word_index in range(4)
		)

	def data_lanes_all_idle(row: bytes) -> bool:
		return all(
			row[FRAME_WORD_OFFSET + frame_word_index * WORD_SIZE:FRAME_WORD_OFFSET + (frame_word_index + 1) * WORD_SIZE] == FRAME_IDLE_WORD
			for frame_word_index in range(4)
		)

	def idle_ac_preview_category(row: bytes) -> str:
		if data_lanes_all_zero(row):
			return "all0"
		if data_lanes_all_idle(row):
			return "allAC"
		return "mixed"

	def describe_stream_row(row: bytes) -> str:
		parts: list[str] = []
		lpgbt_index = row[1]
		current_word_counts = current_word_counts_by_lpgbt.setdefault(lpgbt_index, [0] * 4)
		all_zero = data_lanes_all_zero(row)
		fmc_word_start = FRAME_WORD_OFFSET + 4 * WORD_SIZE
		fmc_word = row[fmc_word_start:fmc_word_start + WORD_SIZE]
		if not data_lanes_have_activity(row):
			if fmc_word == FMC_INJECTION_WORD:
				return "FMC=inj"
			if fmc_word == FMC_L1A_WORD:
				return "FMC=L1A"
			if fmc_word == FMC_IDLE_WORD:
				state = get_fmc_state(lpgbt_index)
				if state.get("pending") and int(state["l1a"]) == 0:
					return "FMC=idle"
				return ""
			if fmc_word in {FMC_IDLE_WORD, ZERO_WORD}:
				return ""
			if fmc_word != FMC_IDLE_WORD:
				return "FMC=other"
			return ""
		for frame_word_index in range(4):
			word_start = FRAME_WORD_OFFSET + frame_word_index * WORD_SIZE
			word = row[word_start:word_start + WORD_SIZE]
			if word == ZERO_WORD and all_zero:
				parts.append("zero")
				continue
			if word == FRAME_IDLE_WORD:
				parts.append("idle")
				continue
			word_index = current_word_counts[frame_word_index] % DATA_FRAME_WORD_COUNT
			current_word_counts[frame_word_index] += 1
			reordered_word = reorder_frame_word(word)
			raw_word = int.from_bytes(reordered_word, "big")
			if word_index == 0:
				parts.append("header")
			elif word_index == DATA_FRAME_WORD_COUNT - 1:
				parts.append("crc")
			else:
				channel = parse_data_frame_channel(word_index, raw_word)
				parts.append(f"{channel.name:<5s} ADC={channel.adc:4d} ToT={channel.tot:4d} ToA={channel.toa:4d}")
		if fmc_word == FMC_INJECTION_WORD:
			parts.append("FMC=inj")
		elif fmc_word == FMC_L1A_WORD:
			parts.append("FMC=L1A")
		return " | ".join(parts)

	with path.open("rb", buffering=1024 * 1024) as data_file:
		while not stop_parse:
			chunk = data_file.read(chunk_size)
			if not chunk:
				break

			full_rows = len(chunk) // ROW_SIZE
			full_bytes = full_rows * ROW_SIZE
			result.total_bytes_parsed += full_bytes

			for offset in range(0, full_bytes, ROW_SIZE):
				if max_rows is not None and line_no >= max_rows:
					stop_parse = True
					break
				line_no += 1
				row = chunk[offset:offset + ROW_SIZE]

				if heartbeat_continuation_pending:
					result.counts[TYPE_HEARTBEAT] += 1
					result.heartbeat_frames += 1
					heartbeat_continuation_pending = False
					continue

				if row[0] == 0x07:
					result.counts[TYPE_HEARTBEAT] += 1
					heartbeat_continuation_pending = True
					continue

				row_type = classify_row(row)
				result.counts[row_type] += 1

				if row_type == TYPE_TRIGGER:
					flush_fmc_states()
					flush_idle_ac_preview()
					trigger_timestamp = extract_trigger_timestamp(row)
					last_trigger_timestamp = trigger_timestamp
					current_word_counts_by_lpgbt = {}
					if seen_trigger:
						finish_trigger_interval()
						result.trigger_interval_count += 1
						result.trigger_interval_lpgbt_totals.update(current_trigger_interval_lpgbt_counts)
						result.trigger_interval_dummy_total += current_trigger_interval_dummy_count
						current_trigger_interval_lpgbt_counts.clear()
						current_trigger_interval_dummy_count = 0
						if max_bb_intervals is not None and result.trigger_interval_count >= max_bb_intervals:
							stop_parse = True
					else:
						seen_trigger = True
					current_interval_no += 1
					current_first_l1a_lpgbt_line_deltas = {}
					current_l1a_lpgbt_line_counts = Counter()
					reset_fmc_states()
					if preview_has_space():
						result.stream_preview.append(StreamPreview(line_no, row_type, None, trigger_timestamp, None, row.hex(" ")))
				elif row_type == TYPE_STREAM:
					lpgbt_index = row[1]
					stream_timestamp = extract_stream_timestamp(row)
					trigger_delta = None if last_trigger_timestamp is None else stream_timestamp - last_trigger_timestamp
					parsed_text = describe_stream_row(row) if seen_trigger else ""
					result.lpgbt_counts[lpgbt_index] += 1
					if data_lanes_all_zero(row):
						result.lpgbt_all_zero_counts[lpgbt_index] += 1
					elif data_lanes_all_idle(row):
						result.lpgbt_all_idle_counts[lpgbt_index] += 1
					if seen_trigger:
						current_trigger_interval_lpgbt_counts[lpgbt_index] += 1
						update_fmc_state(lpgbt_index, row)
						if has_fmc_l1a_word(row):
							current_l1a_lpgbt_line_counts[lpgbt_index] += 1
							if lpgbt_index not in current_first_l1a_lpgbt_line_deltas:
								current_first_l1a_lpgbt_line_deltas[lpgbt_index] = current_trigger_interval_lpgbt_counts[lpgbt_index]
					if not preview_idle and is_idle_stream_row(row):
						idle_ac_preview_counts[idle_ac_preview_category(row)] += 1
					elif preview_has_space():
						flush_idle_ac_preview()
						result.stream_preview.append(StreamPreview(line_no, row_type, lpgbt_index, stream_timestamp, trigger_delta, row.hex(" "), parsed_text))
				elif row_type == TYPE_DUMMY:
					if seen_trigger:
						current_trigger_interval_dummy_count += 1
				elif row_type == TYPE_OTHER:
					result.other_first_byte_counts[row[0]] += 1
				if stop_parse:
					break

			tail_bytes = len(chunk) - full_bytes
			if tail_bytes:
				result.partial_tail_bytes += tail_bytes
				result.counts[TYPE_PARTIAL] += 1

	if heartbeat_continuation_pending:
		result.heartbeat_unpaired_lines += 1
	flush_fmc_states()
	flush_idle_ac_preview()

	result.total_lines = line_no
	result.elapsed_sec = time.perf_counter() - start_time
	return result


def parse_ch2g_frames(path: Path, preview_limit: int = 40, chunk_rows: int = 262144, max_rows: int | None = None, max_bb_intervals: int | None = None, crc_error_word_limit: int | None = None) -> FrameParseResult:
	start_time = time.perf_counter()
	file_size = path.stat().st_size
	result = FrameParseResult(path=path, file_size=file_size, elapsed_sec=0.0)
	chunk_size = ROW_SIZE * chunk_rows
	line_no = 0
	seen_trigger = False
	interval_no = 0
	current_trigger_timestamp: int | None = None
	current_frame_words = 0
	current_idle_words = 0
	current_zero_words = 0
	current_heartbeat_payload_words = 0
	current_word_counts = [0] * FRAME_WORD_COUNT
	current_preview: list[FrameWordPreview] = []
	current_event_words = [[[] for _ in range(FRAME_WORD_COUNT)] for _ in range(FRAME_EVENT_PREVIEW_COUNT)]
	current_crc_sample_counts: Counter[tuple[int, int]] = Counter()
	current_crc_words: dict[tuple[int, int], list[bytes]] = {}
	current_crc_previews: dict[tuple[int, int], list[FrameWordPreview]] = {}
	crc_error_preview_words = 0
	open_samples: dict[int, OpenFrameSample] = {}
	stop_parse = False

	def data_lanes_all_zero(row: bytes) -> bool:
		return all(
			row[FRAME_WORD_OFFSET + frame_word_index * WORD_SIZE:FRAME_WORD_OFFSET + (frame_word_index + 1) * WORD_SIZE] == ZERO_WORD
			for frame_word_index in range(4)
		)

	def flush_event_previews() -> None:
		if interval_no > FRAME_SAMPLE_PREVIEW_INTERVALS:
			return
		for event_index, link_words in enumerate(current_event_words):
			if not any(len(words) >= SAMPLE_WORD_COUNT for words in link_words):
				continue
			result.event_previews.append(FrameEventPreview(
				interval_no=interval_no,
				event_index=event_index,
				link_words=tuple(tuple(words[:SAMPLE_WORD_COUNT]) for words in link_words),
			))

	def append_crc_error_sample(link_id: int, event_index: int, parsed: DataFrameParseResult, previews: list[FrameWordPreview]) -> None:
		nonlocal crc_error_preview_words
		if parsed.crc_ok:
			return
		if crc_error_word_limit is not None and crc_error_preview_words + len(previews) > crc_error_word_limit:
			return
		result.crc_error_samples.append(FrameCrcErrorSample(
			interval_no=interval_no,
			event_index=event_index,
			link_id=link_id,
			parsed=parsed,
			previews=list(previews),
		))
		crc_error_preview_words += len(previews)

	with path.open("rb", buffering=1024 * 1024) as data_file:
		while not stop_parse:
			chunk = data_file.read(chunk_size)
			if not chunk:
				break

			full_rows = len(chunk) // ROW_SIZE
			full_bytes = full_rows * ROW_SIZE
			result.total_bytes_parsed += full_bytes

			for offset in range(0, full_bytes, ROW_SIZE):
				if max_rows is not None and line_no >= max_rows:
					stop_parse = True
					break
				line_no += 1
				row = chunk[offset:offset + ROW_SIZE]
				row_type = classify_row(row)
				all_data_lanes_zero = data_lanes_all_zero(row)

				if row_type == TYPE_TRIGGER:
					current_trigger_timestamp = extract_trigger_timestamp(row)
					result.trigger_lines += 1
					if seen_trigger:
						result.trigger_intervals += 1
						if len(result.interval_summaries) < preview_limit:
							result.interval_summaries.append(
								FrameIntervalSummary(interval_no, tuple(count // SAMPLE_WORD_COUNT for count in current_word_counts))
							)
						flush_event_previews()
						result.incomplete_samples += len(open_samples)
						result.frame_words += current_frame_words
						result.idle_words += current_idle_words
						result.zero_words += current_zero_words
						result.heartbeat_payload_words += current_heartbeat_payload_words
						result.preview.extend(current_preview)
						current_frame_words = 0
						current_idle_words = 0
						current_zero_words = 0
						current_heartbeat_payload_words = 0
						current_word_counts = [0] * FRAME_WORD_COUNT
						current_preview = []
						current_event_words = [[[] for _ in range(FRAME_WORD_COUNT)] for _ in range(FRAME_EVENT_PREVIEW_COUNT)]
						current_crc_sample_counts = Counter()
						current_crc_words = {}
						current_crc_previews = {}
						open_samples = {}
						if max_bb_intervals is not None and result.trigger_intervals >= max_bb_intervals:
							stop_parse = True
					else:
						seen_trigger = True
					interval_no += 1
					if stop_parse:
						break
					continue

				for frame_word_index in range(FRAME_WORD_COUNT):
					word_start = FRAME_WORD_OFFSET + frame_word_index * WORD_SIZE
					word = row[word_start:word_start + WORD_SIZE]
					if not seen_trigger:
						if word != ZERO_WORD and word != FRAME_IDLE_WORD:
							result.pre_trigger_words += 1
						continue
					if row_type == TYPE_HEARTBEAT:
						if word != ZERO_WORD and word != FRAME_IDLE_WORD:
							current_heartbeat_payload_words += 1
						continue
					if word == ZERO_WORD and (frame_word_index >= 4 or all_data_lanes_zero):
						current_zero_words += 1
						continue
					if word == FRAME_IDLE_WORD:
						current_idle_words += 1
						continue
					reordered_word = reorder_frame_word(word)
					timestamp = extract_stream_timestamp(row) if row_type == TYPE_STREAM else None
					bb_delta = None if timestamp is None or current_trigger_timestamp is None else timestamp - current_trigger_timestamp
					current_frame_words += 1
					word_count = current_word_counts[frame_word_index]
					current_word_counts[frame_word_index] += 1
					event_index = word_count // SAMPLE_WORD_COUNT
					word_index = word_count % SAMPLE_WORD_COUNT
					word_preview = FrameWordPreview(
						interval_no,
						line_no,
						frame_word_index,
						word_start,
						row_type,
						row[1] if row_type == TYPE_STREAM else None,
						timestamp,
						bb_delta,
						word,
						reordered_word,
					)
					if len(result.event_previews) < FRAME_EVENT_PREVIEW_COUNT:
						if event_index < FRAME_EVENT_PREVIEW_COUNT:
							current_event_words[event_index][frame_word_index].append(reordered_word)
					if frame_word_index == 0 and len(result.preview) + len(current_preview) < preview_limit:
						current_preview.append(word_preview)
					if row_type == TYPE_STREAM and frame_word_index < 4:
						lpgbt_index = row[1]
						crc_key = (lpgbt_index, frame_word_index)
						if crc_key not in current_crc_words:
							if not is_sample_header_word(reordered_word):
								continue
							current_crc_words[crc_key] = []
							current_crc_previews[crc_key] = []
						current_crc_words[crc_key].append(reordered_word)
						current_crc_previews[crc_key].append(word_preview)
						if len(current_crc_words[crc_key]) == DATA_FRAME_WORD_COUNT:
							append_crc_error_sample(
								frame_word_index,
								current_crc_sample_counts[crc_key],
								parse_h2gcroc_data_frame(current_crc_words[crc_key]),
								current_crc_previews[crc_key],
							)
							current_crc_sample_counts[crc_key] += 1
							current_crc_words.pop(crc_key, None)
							current_crc_previews.pop(crc_key, None)

					if timestamp is None:
						continue

					open_sample = open_samples.get(frame_word_index)
					if open_sample is not None and timestamp - open_sample.header_timestamp >= SAMPLE_TIMESTAMP_WINDOW:
						result.expired_samples += 1
						open_samples.pop(frame_word_index)
						open_sample = None

					if open_sample is None and is_sample_header_word(reordered_word):
						result.header_words += 1
						open_samples[frame_word_index] = OpenFrameSample(
							interval_no=interval_no,
							link_id=frame_word_index,
							header_line_no=line_no,
							header_timestamp=timestamp,
							header_word=reordered_word,
							last_line_no=line_no,
							last_timestamp=timestamp,
						)
						continue

					if open_sample is None:
						continue
					open_sample.words.append(reordered_word)
					open_sample.last_line_no = line_no
					open_sample.last_timestamp = timestamp
					if len(open_sample.words) == SAMPLE_WORD_COUNT:
						result.complete_samples += 1
						if open_sample.interval_no <= FRAME_SAMPLE_PREVIEW_INTERVALS and open_sample.link_id in FRAME_SAMPLE_PREVIEW_LINKS:
							result.sample_preview.append(CompleteFrameSample(
								interval_no=open_sample.interval_no,
								link_id=open_sample.link_id,
								header_line_no=open_sample.header_line_no,
								header_timestamp=open_sample.header_timestamp,
								header_word=open_sample.header_word,
								last_line_no=open_sample.last_line_no,
								last_timestamp=open_sample.last_timestamp,
								words=list(open_sample.words),
							))
						open_samples.pop(frame_word_index)

			tail_bytes = len(chunk) - full_bytes
			if tail_bytes:
				result.partial_tail_bytes += tail_bytes

	result.total_lines = line_no
	result.trailing_open_words = current_frame_words
	result.incomplete_samples += len(open_samples)
	result.elapsed_sec = time.perf_counter() - start_time
	return result


def format_count(value: int) -> str:
	return f"{value:,}"


def format_word(word: bytes) -> str:
	return word.hex(" ")


def format_word_value(word: bytes) -> str:
	return f"0x{int.from_bytes(word, 'little'):08X}"


def format_reordered_word_value(word: bytes) -> str:
	return f"0x{int.from_bytes(word, 'big'):08X}"


def format_data_frame_parse(parsed: DataFrameParseResult, indent: str = "") -> list[str]:
	lines: list[str] = []
	status = "OK" if parsed.ok else "FAIL"
	lines.append(f"{indent}parse {status}")
	if parsed.errors:
		for error in parsed.errors:
			lines.append(f"{indent}  error: {error}")
	if parsed.header is not None:
		first_event_text = "yes" if parsed.header.first_event is True else "no" if parsed.header.first_event is False else "unknown"
		lines.append(
			f"{indent}  header  Bx {parsed.header.bx_counter}  Event {parsed.header.event_counter}  "
			f"Orbit {parsed.header.orbit_counter}  "
			f"H {parsed.header.hamming_flags[0]},{parsed.header.hamming_flags[1]},{parsed.header.hamming_flags[2]}  "
			f"Tr 0x{parsed.header.trailer:X}  first_event {first_event_text}"
		)
	if parsed.crc_word is not None and parsed.computed_crc is not None:
		crc_status = "OK" if parsed.crc_ok else "FAIL"
		lines.append(f"{indent}  crc {crc_status}  word 0x{parsed.crc_word:08X}  computed 0x{parsed.computed_crc:08X}")
	if parsed.channels:
		lines.append(f"{indent}  channels: name Tc Tp ADC ToT ToA")
		for channel in parsed.channels:
			lines.append(
				f"{indent}    {channel.name:5s}  {channel.tc}  {channel.tp}  "
				f"{channel.adc:4d} {channel.tot:4d} {channel.toa:4d}"
			)
	return lines


def format_compact_frame_sample(interval_no: int, sample_index: int, link_id: int, words: Sequence[bytes], crc_errors_only: bool = False) -> list[str]:
	parsed = parse_h2gcroc_data_frame(words)
	if crc_errors_only and parsed.crc_ok:
		return []
	status = "OK" if parsed.ok else "FAIL"
	if parsed.header is None:
		return [
			f"    sample {sample_index + 1}: parse {status}",
			"      no decoded header",
		]
	first_event_text = "yes" if parsed.header.first_event is True else "no" if parsed.header.first_event is False else "unknown"
	lines = [
		f"    sample {sample_index + 1}: "
		f"parse {status}  Bx {parsed.header.bx_counter}  Event {parsed.header.event_counter}  "
		f"Orbit {parsed.header.orbit_counter}  first_event {first_event_text}",
	]
	if crc_errors_only and parsed.crc_word is not None and parsed.computed_crc is not None:
		crc_status = "OK" if parsed.crc_ok else "FAIL"
		lines.append(f"      crc {crc_status}  word 0x{parsed.crc_word:08X}  computed 0x{parsed.computed_crc:08X}")
	if parsed.errors:
		displayed_errors = [error for error in parsed.errors if not error.startswith("CRC mismatch")]
		lines.append("      errors: " + "; ".join(displayed_errors[:2]))
	channel_text = "  ".join(
		f"{channel.name}:{channel.adc}"
		for channel in parsed.channels
		if channel.name.startswith("ch")
	)
	lines.append(f"      ADC {channel_text}")
	return lines


def describe_crc_sample_word(word_index: int, raw_word: int, parsed: DataFrameParseResult) -> str:
	if word_index == 0 and parsed.header is not None:
		first_event_text = "yes" if parsed.header.first_event is True else "no" if parsed.header.first_event is False else "unknown"
		return (
			f"header Bx={parsed.header.bx_counter} Event={parsed.header.event_counter} "
			f"Orbit={parsed.header.orbit_counter} Tr=0x{parsed.header.trailer:X} first_event={first_event_text}"
		)
	if word_index == DATA_FRAME_WORD_COUNT - 1 and parsed.crc_word is not None and parsed.computed_crc is not None:
		crc_status = "OK" if parsed.crc_ok else "FAIL"
		return f"crc {crc_status} word=0x{parsed.crc_word:08X} computed=0x{parsed.computed_crc:08X}"
	channel = parse_data_frame_channel(word_index, raw_word)
	return f"{channel.name:<5s} ADC={channel.adc:4d} ToT={channel.tot:4d} ToA={channel.toa:4d}"


def format_crc_error_sample(sample: FrameCrcErrorSample) -> list[str]:
	crc_text = "crc FAIL"
	if sample.parsed.crc_word is not None and sample.parsed.computed_crc is not None:
		crc_text = f"crc FAIL word 0x{sample.parsed.crc_word:08X} computed 0x{sample.parsed.computed_crc:08X}"
	lines = [f"  BB#{sample.interval_no:6d}  link {sample.link_id}  sample {sample.event_index + 1}: {crc_text}"]
	for word_index, preview in enumerate(sample.previews):
		raw_word = int.from_bytes(preview.reordered_word, "big")
		lpgbt_text = "lpGBT  -" if preview.lpgbt_index is None else f"lpGBT{preview.lpgbt_index:3d}"
		line = (
			f"    {lpgbt_text}  link {preview.link_id}  w{word_index:02d}  "
			f"dBB {format_delta_value(preview.bb_delta)}  {format_word(preview.word)}"
		)
		line += f"    || {describe_crc_sample_word(word_index, raw_word, sample.parsed)}"
		lines.append(line)
	return lines


def compact_preview_by_interval_and_link(event_previews: Sequence[FrameEventPreview]) -> dict[int, dict[int, list[tuple[int, tuple[bytes, ...]]]]]:
	preview: dict[int, dict[int, list[tuple[int, tuple[bytes, ...]]]]] = {}
	for event in event_previews:
		for link_id, words in enumerate(event.link_words):
			if link_id not in FRAME_FULL_PREVIEW_LINKS or len(words) < SAMPLE_WORD_COUNT:
				continue
			preview.setdefault(event.interval_no, {}).setdefault(link_id, []).append((event.event_index, words))
	return preview


def format_average_count(value: float) -> str:
	if value.is_integer():
		return format_count(int(value))
	return f"{value:.3f}"


def color_byte(byte_text: str, color_mode: str, ansi_color: str, rich_color: str) -> str:
	if color_mode == "ansi":
		return f"{ansi_color}{byte_text}{ANSI_RESET}"
	if color_mode == "rich":
		return f"[{rich_color}]{byte_text}[/]"
	return byte_text


def color_stream_hex(hex_text: str, color_mode: str) -> str:
	if color_mode == "plain":
		return hex_text

	bytes_text = hex_text.split()
	colored = []
	for position, byte_text in enumerate(bytes_text):
		if position == 0:
			colored.append(color_byte(byte_text, color_mode, ANSI_GRAY, "grey50"))
		elif position == 1:
			index_value = int(byte_text, 16)
			palette_index = index_value % len(ANSI_INDEX_COLORS)
			colored.append(color_byte(byte_text, color_mode, ANSI_INDEX_COLORS[palette_index], RICH_INDEX_COLORS[palette_index]))
		elif position < 8:
			colored.append(color_byte(byte_text, color_mode, ANSI_BLUE, "blue"))
		else:
			palette_index = ((position - 8) // 4) % len(ANSI_GROUP_COLORS)
			colored.append(color_byte(byte_text, color_mode, ANSI_GROUP_COLORS[palette_index], RICH_GROUP_COLORS[palette_index]))
	return " ".join(colored)


def color_trigger_hex(hex_text: str, color_mode: str) -> str:
	if color_mode == "plain":
		return hex_text

	bytes_text = hex_text.split()
	colored = []
	for position, byte_text in enumerate(bytes_text):
		if position < 4:
			colored.append(color_byte(byte_text, color_mode, ANSI_GRAY, "grey50"))
		else:
			colored.append(byte_text)
	return " ".join(colored)


def color_preview_hex(row_type: str, hex_text: str, color_mode: str) -> str:
	if row_type == TYPE_TRIGGER:
		return color_trigger_hex(hex_text, color_mode)
	return color_stream_hex(hex_text, color_mode)


def result_lines(result: ParseResult, color_mode: str = "plain") -> list[str]:
	all_zero_ac_lines = sum(result.lpgbt_all_zero_counts.values())
	all_idle_ac_lines = sum(result.lpgbt_all_idle_counts.values())
	other_ac_lines = result.counts[TYPE_STREAM] - all_zero_ac_lines - all_idle_ac_lines
	lines = [
		f"File: {result.path}",
		f"Size: {format_count(result.file_size)} bytes",
		f"Parsed lines: {format_count(result.total_lines)} x {ROW_SIZE} bytes",
		f"Speed: {result.elapsed_sec:.3f} s, {result.mib_per_sec:.1f} MiB/s",
		"",
		"Line type statistics:",
		f"  heartbeat lines : {format_count(result.counts[TYPE_HEARTBEAT])} ({format_count(result.heartbeat_frames)} frames, 2 lines/frame)",
		f"  heartbeat unpaired lines : {format_count(result.heartbeat_unpaired_lines)}",
		f"  physical trigger BB lines : {format_count(result.counts[TYPE_TRIGGER])}",
		f"  data stream AC lines : {format_count(result.counts[TYPE_STREAM])} (all0={format_count(all_zero_ac_lines)}, allAC={format_count(all_idle_ac_lines)}, other={format_count(other_ac_lines)})",
		f"  firmware dummy 00 lines : {format_count(result.counts[TYPE_DUMMY])}",
		f"  other lines : {format_count(result.counts[TYPE_OTHER])}",
	]

	if result.partial_tail_bytes:
		lines.append(f"  partial tail bytes : {format_count(result.partial_tail_bytes)}")

	lines.append("")
	lines.append("Data stream lines by lpGBT index:")
	if result.lpgbt_counts:
		for index, count in sorted(result.lpgbt_counts.items()):
			all_zero = result.lpgbt_all_zero_counts[index]
			all_idle = result.lpgbt_all_idle_counts[index]
			other = count - all_zero - all_idle
			lines.append(f"  lpGBT {index:3d} : {format_count(count)} (all0={format_count(all_zero)}, allAC={format_count(all_idle)}, other={format_count(other)})")
	else:
		lines.append("  none")

	lines.append("")
	lines.append("Average AC lines between two BB trigger lines:")
	if result.trigger_interval_count and result.lpgbt_counts:
		lines.append(f"  BB intervals : {format_count(result.trigger_interval_count)}")
		average_dummy = result.trigger_interval_dummy_total / result.trigger_interval_count
		for index in sorted(result.lpgbt_counts):
			average = result.trigger_interval_lpgbt_totals[index] / result.trigger_interval_count
			lines.append(f"  lpGBT {index:3d} : {format_average_count(average)} (dummy lines: {format_average_count(average_dummy)})")
	else:
		lines.append("  none")

	lines.append("")
	lines.append("First FMC 0x4b4b4b4b AC-line distance from previous BB line by lpGBT:")
	lpgbt_indices = sorted(set(result.bb_first_l1a_lpgbt_line_deltas) | set(result.bb_first_l1a_lpgbt_missing_intervals))
	if lpgbt_indices:
		for lpgbt_index in lpgbt_indices:
			delta_counts = result.bb_first_l1a_lpgbt_line_deltas.get(lpgbt_index, Counter())
			measured_intervals = sum(delta_counts.values())
			missing_intervals = result.bb_first_l1a_lpgbt_missing_intervals[lpgbt_index]
			lines.append(f"  lpGBT {lpgbt_index:3d}:")
			lines.append(f"    measured BB intervals : {format_count(measured_intervals)}")
			lines.append(f"    missing intervals : {format_count(missing_intervals)}")
			for line_delta, count in sorted(delta_counts.items()):
				lines.append(f"    {format_count(line_delta)} AC lines : {format_count(count)}")
	else:
		lines.append("  none")

	lines.append("")
	lines.append("Paired first FMC 0x4b4b4b4b AC-line distances per BB interval:")
	if result.bb_first_l1a_lpgbt_pair_deltas:
		for (lpgbt0_delta, lpgbt1_delta), interval_count in sorted(
			result.bb_first_l1a_lpgbt_pair_deltas.items(),
			key=lambda item: ((-1 if item[0][0] is None else item[0][0]), (-1 if item[0][1] is None else item[0][1])),
		):
			lpgbt0_text = "missing" if lpgbt0_delta is None else f"{format_count(lpgbt0_delta)} AC"
			lpgbt1_text = "missing" if lpgbt1_delta is None else f"{format_count(lpgbt1_delta)} AC"
			lines.append(f"  lpGBT0 {lpgbt0_text} + lpGBT1 {lpgbt1_text} : {format_count(interval_count)}")
	else:
		lines.append("  none")

	lines.append("")
	lines.append("FMC 0x4b4b4b4b AC-line counts per BB interval by lpGBT:")
	count_lpgbt_indices = sorted(result.bb_l1a_lpgbt_line_counts)
	if count_lpgbt_indices:
		for lpgbt_index in count_lpgbt_indices:
			lines.append(f"  lpGBT {lpgbt_index:3d}:")
			for l1a_count, interval_count in sorted(result.bb_l1a_lpgbt_line_counts[lpgbt_index].items()):
				lines.append(f"    {format_count(l1a_count)} L1A AC lines : {format_count(interval_count)}")
	else:
		lines.append("  none")

	if result.other_first_byte_counts:
		lines.append("")
		lines.append("Other first-byte statistics:")
		for first_byte, count in result.other_first_byte_counts.most_common():
			lines.append(f"  0x{first_byte:02X} : {format_count(count)}")

	lines.append("")
	lines.append("FMC injection/L1A timing by BB interval and lpGBT:")
	if result.fmc_summaries:
		for summary in result.fmc_summaries:
			lines.append(
				f"  BB#{summary.interval_no:6d}  lpGBT {summary.lpgbt_index:3d}  "
				f"inj#{summary.injection_no}: idle36={format_count(summary.idle_before_l1a)}  "
				f"L1A4b={format_count(summary.l1a_count)}"
			)
	else:
		lines.append("  none")

	lines.append("")
	lines.append(f"First {len(result.stream_preview)} AC/BB lines:")
	if result.stream_preview:
		for item in result.stream_preview:
			if item.row_type == TYPE_IDLE_AC_SUMMARY:
				lines.append(f"  {item.parsed_text}")
				continue
			hex_text = color_preview_hex(item.row_type, item.hex_text, color_mode)
			if item.row_type == TYPE_TRIGGER:
				lines.append(f"  trigger    {hex_text}")
			else:
				line = f"  lpGBT{item.lpgbt_index:3d}  dBB {format_delta_value(item.trigger_delta)}  {hex_text}"
				if item.parsed_text:
					line += f"    || {item.parsed_text}"
				lines.append(line)
	else:
		lines.append("  none")

	return lines


def frame_result_lines(result: FrameParseResult, crc_errors_only: bool = False) -> list[str]:
	lines = [
		f"File: {result.path}",
		f"Size: {format_count(result.file_size)} bytes",
		f"Parsed lines: {format_count(result.total_lines)} x {ROW_SIZE} bytes",
		f"Speed: {result.elapsed_sec:.3f} s, {result.mib_per_sec:.1f} MiB/s",
		"",
		"Frame mode statistics:",
		f"  BB trigger lines : {format_count(result.trigger_lines)}",
		f"  BB intervals : {format_count(result.trigger_intervals)}",
		f"  non-idle frame words : {format_count(result.frame_words)}",
		f"  idle cc cc cc ac words : {format_count(result.idle_words)}",
		f"  zero 00 00 00 00 words : {format_count(result.zero_words)}",
		f"  skipped heartbeat payload words : {format_count(result.heartbeat_payload_words)}",
		f"  non-idle words before first BB : {format_count(result.pre_trigger_words)}",
		f"  unclosed non-idle words after last BB : {format_count(result.trailing_open_words)}",
		f"  sample header words : {format_count(result.header_words)}",
		f"  complete samples : {format_count(result.complete_samples)}",
		f"  incomplete samples : {format_count(result.incomplete_samples)}",
		f"  expired samples : {format_count(result.expired_samples)}",
	]

	if result.partial_tail_bytes:
		lines.append(f"  partial tail bytes : {format_count(result.partial_tail_bytes)}")

	if crc_errors_only:
		lines.append("")
		lines.append(f"CRC-error samples, {DATA_FRAME_WORD_COUNT} word lines per sample:")
		if result.crc_error_samples:
			for sample in result.crc_error_samples:
				lines.extend(format_crc_error_sample(sample))
		else:
			lines.append("  none")
		return lines

	lines.append("")
	lines.append(f"First {len(result.interval_summaries)} BB interval sample counts:")
	if result.interval_summaries:
		for summary in result.interval_summaries:
			total_samples = sum(summary.sample_counts)
			word_summary = "  ".join(
				f"w{word_index} {format_count(sample_count)}"
				for word_index, sample_count in enumerate(summary.sample_counts)
			)
			lines.append(f"  BB#{summary.interval_no:6d}  {word_summary}  total {format_count(total_samples)}")
	else:
		lines.append("  none")

	lines.append("")
	lines.append(f"Compact samples from first {FRAME_SAMPLE_PREVIEW_INTERVALS} BB intervals for links 0 and 1:")
	compact_preview = compact_preview_by_interval_and_link(result.event_previews)
	if compact_preview:
		displayed_sample = False
		for interval_no, links in sorted(compact_preview.items()):
			interval_lines: list[str] = []
			for link_id in sorted(FRAME_FULL_PREVIEW_LINKS):
				samples = sorted(links.get(link_id, []))
				if not samples:
					if not crc_errors_only:
						interval_lines.append(f"    link {link_id}: no complete preview samples")
					continue
				sample_lines: list[str] = []
				displayed_samples = 0
				for sample_index, words in samples:
					formatted_sample = format_compact_frame_sample(interval_no, sample_index, link_id, words, crc_errors_only=crc_errors_only)
					if formatted_sample:
						displayed_samples += 1
						sample_lines.extend(formatted_sample)
				if sample_lines:
					displayed_sample = True
					interval_lines.append(f"    link {link_id}: {format_count(displayed_samples)} samples")
					interval_lines.extend(sample_lines)
				elif not crc_errors_only:
					interval_lines.append(f"    link {link_id}: {len(samples)} samples")
			if interval_lines:
				lines.append(f"  BB#{interval_no:6d}")
				lines.extend(interval_lines)
		if crc_errors_only and not displayed_sample:
			lines.append("  none")
	else:
		lines.append("  none")

	return lines


def print_result(result: ParseResult) -> None:
	print("\n".join(result_lines(result, color_mode="ansi")))


def print_frame_result(result: FrameParseResult, crc_errors_only: bool = False) -> None:
	print("\n".join(frame_result_lines(result, crc_errors_only=crc_errors_only)))


def run_textual(path: Path, preview_limit: int, chunk_rows: int) -> None:
	try:
		from textual.app import App, ComposeResult
		from textual.containers import Horizontal, Vertical
		from textual.widgets import Button, Footer, Header, Input, RichLog, Static
	except ModuleNotFoundError:
		print("Textual is not installed in this Python environment; using CLI output instead.\n")
		print_result(parse_ch2g(path, preview_limit=preview_limit, chunk_rows=chunk_rows))
		return

	class RawBrowserApp(App):
		TITLE = "CH2G Raw Browser"
		SUB_TITLE = "32-byte line parser"
		BINDINGS = [("q", "quit", "Quit"), ("r", "parse_file", "Parse")]

		CSS = """
		Screen { layout: vertical; }
		#top_bar { height: 3; padding: 0 1; }
		#path_input { width: 1fr; }
		#parse_button { width: 14; }
		#status { height: 1; padding: 0 1; color: $accent; }
		#output { height: 1fr; border: solid $accent 50%; }
		"""

		def compose(self) -> ComposeResult:
			yield Header()
			with Vertical():
				with Horizontal(id="top_bar"):
					yield Input(value=str(path), id="path_input")
					yield Button("Parse", id="parse_button", variant="primary")
				yield Static("Ready", id="status")
				yield RichLog(id="output", wrap=False, highlight=False)
			yield Footer()

		def on_mount(self) -> None:
			self.action_parse_file()

		def on_button_pressed(self, event: Button.Pressed) -> None:
			if event.button.id == "parse_button":
				self.action_parse_file()

		def action_parse_file(self) -> None:
			input_path = Path(self.query_one("#path_input", Input).value).expanduser()
			output = self.query_one("#output", RichLog)
			status = self.query_one("#status", Static)
			output.clear()
			status.update("Parsing...")
			try:
				parsed = parse_ch2g(input_path, preview_limit=preview_limit, chunk_rows=chunk_rows)
			except Exception as exc:
				status.update("Parse failed")
				output.write(f"ERROR: {exc}")
				return
			status.update(f"Parsed {format_count(parsed.total_lines)} lines in {parsed.elapsed_sec:.3f} s")
			for line in result_lines(parsed, color_mode="rich"):
				output.write(line)

	RawBrowserApp().run()


def parse_and_print(path: Path, mode: str, preview_limit: int, chunk_rows: int, max_rows: int | None = None, max_bb_intervals: int | None = None, preview_idle: bool = False, crc_errors_only: bool = False) -> None:
	if mode == "frames":
		crc_error_word_limit = preview_limit if crc_errors_only else None
		print_frame_result(parse_ch2g_frames(path, preview_limit=preview_limit, chunk_rows=chunk_rows, max_rows=max_rows, max_bb_intervals=max_bb_intervals, crc_error_word_limit=crc_error_word_limit), crc_errors_only=crc_errors_only)
	else:
		print_result(parse_ch2g(path, preview_limit=preview_limit, chunk_rows=chunk_rows, max_rows=max_rows, max_bb_intervals=max_bb_intervals, preview_idle=preview_idle))


def build_arg_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(description="Fast 32-byte line browser/statistics for dataXXX.ch2g files.")
	parser.add_argument("file", nargs="?", type=Path, help="dataXXX.ch2g file, default: latest non-empty COMMON/temp/data*.ch2g")
	mode = parser.add_mutually_exclusive_group()
	mode.add_argument("--words", action="store_const", const="words", dest="mode", help="word/line parser mode (default)")
	mode.add_argument("--frames", action="store_const", const="frames", dest="mode", help="frame parser mode: print sample counts and the first two full events")
	parser.set_defaults(mode="words")
	parser.add_argument("--preview", type=int, default=40, help="number of data stream lines to print")
	parser.add_argument("--chunk-rows", type=int, default=262144, help="rows read per chunk; larger is faster but uses more memory")
	parser.add_argument("--max-rows", type=int, help="stop after parsing this many 32-byte rows")
	parser.add_argument("--max-bb-intervals", type=int, help="stop after this many complete BB-to-BB intervals")
	parser.add_argument("--show-idle-ac", action="store_true", help="include AC rows whose frame-word payload is only idle/zero words in --words preview")
	parser.add_argument("--crc", action="store_true", help="use frame parser mode and show only compact samples with CRC mismatch")
	parser.add_argument("--ui", action="store_true", help="start Textual UI for --words mode; falls back to CLI if Textual is unavailable")
	return parser


def main() -> None:
	if hasattr(signal, "SIGPIPE"):
		signal.signal(signal.SIGPIPE, signal.SIG_DFL)

	args = build_arg_parser().parse_args()
	path = args.file.expanduser() if args.file else latest_nonempty_ch2g()
	mode = "frames" if args.crc else args.mode

	if args.ui and mode == "words":
		run_textual(path, args.preview, args.chunk_rows)
	else:
		parse_and_print(path, mode, args.preview, args.chunk_rows, max_rows=args.max_rows, max_bb_intervals=args.max_bb_intervals, preview_idle=args.show_idle_ac, crc_errors_only=args.crc)


if __name__ == "__main__":
	main()
