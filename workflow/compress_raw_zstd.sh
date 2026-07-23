#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="data"
OUT_DIR="data/zstd"
LEVEL="3"
THREADS="0"
FORCE=0
INPUT_FILE=""

usage() {
	cat <<'USAGE'
Usage: workflow/compress_raw_zstd.sh [options]

Compress all data/*.raw files to data/zstd/*.raw.zst and verify byte-for-byte
consistency with SHA-256 of the original stream and decompressed stream.

Options:
	-i, --input FILE     Compress only this raw file
  -d, --data-dir DIR    Input data directory (default: data)
  -o, --out-dir DIR     Output zstd directory (default: data/zstd)
  -l, --level N         zstd compression level (default: 3)
  -T, --threads N       zstd thread count, 0 means all cores (default: 0)
  -f, --force           Recompress even if a verified .zst already exists
  -h, --help            Show this help

Outputs per raw file:
  OUT_DIR/<name>.raw.zst
  OUT_DIR/<name>.raw.zst.sha256

The .sha256 file records both original and decompressed SHA-256 plus sizes.
USAGE
}

while (($#)); do
	case "$1" in
		-i|--input)
			INPUT_FILE="${2:?missing value for $1}"
			shift 2
			;;
		-d|--data-dir)
			DATA_DIR="${2:?missing value for $1}"
			shift 2
			;;
		-o|--out-dir)
			OUT_DIR="${2:?missing value for $1}"
			shift 2
			;;
		-l|--level)
			LEVEL="${2:?missing value for $1}"
			shift 2
			;;
		-T|--threads)
			THREADS="${2:?missing value for $1}"
			shift 2
			;;
		-f|--force)
			FORCE=1
			shift
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			echo "Unknown option: $1" >&2
			usage >&2
			exit 2
			;;
	esac
done

command -v zstd >/dev/null || { echo "missing required command: zstd" >&2; exit 127; }
command -v sha256sum >/dev/null || { echo "missing required command: sha256sum" >&2; exit 127; }
command -v stat >/dev/null || { echo "missing required command: stat" >&2; exit 127; }

mkdir -p "$OUT_DIR"

shopt -s nullglob
if [[ -n "$INPUT_FILE" ]]; then
	if [[ ! -f "$INPUT_FILE" ]]; then
		echo "Input raw file does not exist: $INPUT_FILE" >&2
		exit 1
	fi
	raw_files=("$INPUT_FILE")
else
	raw_files=("$DATA_DIR"/*.raw)
fi
if ((${#raw_files[@]} == 0)); then
	echo "No raw files found under $DATA_DIR" >&2
	exit 1
fi

total=${#raw_files[@]}
index=0
for raw in "${raw_files[@]}"; do
	base=$(basename "$raw")
	case "$base" in
		._*)
			continue
			;;
	esac
	((++index))
	zst="$OUT_DIR/$base.zst"
	meta="$zst.sha256"
	tmp="$zst.tmp.$$"

	echo "[$index/$total] $raw -> $zst"

	original_size=$(stat -c '%s' "$raw")
	if [[ -s "$zst" && -s "$meta" && "$FORCE" -eq 0 ]]; then
		stored_original=$(awk -F'  ' '$2 == "original" {print $1}' "$meta" || true)
		stored_decompressed=$(awk -F'  ' '$2 == "decompressed" {print $1}' "$meta" || true)
		stored_size=$(awk -F'  ' '$2 == "original_size" {print $1}' "$meta" || true)
		if [[ -n "$stored_original" && "$stored_original" == "$stored_decompressed" && "$stored_size" == "$original_size" ]]; then
			echo "  verified output already exists; skipping"
			continue
		fi
		echo "  existing metadata is missing or stale; recompressing"
	fi

	rm -f "$tmp"
	zstd -T"$THREADS" -"$LEVEL" -o "$tmp" -- "$raw"
	zstd -t -- "$tmp"

	original_hash=$(sha256sum "$raw" | awk '{print $1}')
	decompressed_hash=$(zstd -dc -- "$tmp" | sha256sum | awk '{print $1}')
	if [[ "$original_hash" != "$decompressed_hash" ]]; then
		rm -f "$tmp"
		echo "  ERROR: SHA-256 mismatch for $raw" >&2
		echo "  original:     $original_hash" >&2
		echo "  decompressed: $decompressed_hash" >&2
		exit 1
	fi

	mv -f "$tmp" "$zst"
	compressed_size=$(stat -c '%s' "$zst")
	{
		printf '%s  original\n' "$original_hash"
		printf '%s  decompressed\n' "$decompressed_hash"
		printf '%s  original_size\n' "$original_size"
		printf '%s  compressed_size\n' "$compressed_size"
	} > "$meta"

	ratio=$(awk -v compressed="$compressed_size" -v original="$original_size" 'BEGIN { if (original == 0) print "0.00"; else printf "%.2f", 100.0 * compressed / original }')
	echo "  ok: $compressed_size / $original_size bytes (${ratio}%)"
done

echo "Done. Compressed files are under $OUT_DIR"
