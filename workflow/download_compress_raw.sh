#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="sjia@lxplus.cern.ch"
REMOTE_DIR="/eos/experiment/alice/focal/TB_2026_Wk27_H4/data/focalh-nas-1/Data"
DATA_DIR="data"
OUT_DIR="data/zstd"
LEVEL="3"
THREADS="0"
IDENTITY_FILE=""
JOBS="1"
CONTROL_PATH=""

usage() {
	cat <<'USAGE'
Usage: workflow/download_compress_raw.sh [options]

Discover all remote .raw files and process them in bounded parallel batches: download,
compress with workflow/compress_raw_zstd.sh, verify SHA-256, then delete the
local raw file. Existing verified .raw.zst files are skipped.

Options:
  -H, --host HOST       SSH host (default: sjia@lxplus.cern.ch)
	-i, --identity FILE   OpenSSH private key used by ssh and scp
	-S, --control-path P  Reuse an authenticated SSH control socket
  -r, --remote-dir DIR  Remote directory containing .raw files
  -d, --data-dir DIR    Temporary raw download directory (default: data)
  -o, --out-dir DIR     Compressed output directory (default: data/zstd)
	-j, --jobs N          Files to download/compress in parallel (default: 1)
  -l, --level N         zstd compression level (default: 3)
  -T, --threads N       zstd thread count, 0 means all cores (default: 0)
  -h, --help            Show this help

Interrupted downloads remain as DATA_DIR/<name>.part and resume on the next
run. Complete local .raw files are reused. Raw files are deleted only
after the compressed stream and its checksum metadata have been verified.
USAGE
}

while (($#)); do
	case "$1" in
		-H|--host)
			REMOTE_HOST="${2:?missing value for $1}"
			shift 2
			;;
		-i|--identity)
			IDENTITY_FILE="${2:?missing value for $1}"
			shift 2
			;;
		-S|--control-path)
			CONTROL_PATH="${2:?missing value for $1}"
			shift 2
			;;
		-r|--remote-dir)
			REMOTE_DIR="${2:?missing value for $1}"
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
		-j|--jobs)
			JOBS="${2:?missing value for $1}"
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

if [[ ! "$JOBS" =~ ^[1-9][0-9]*$ ]]; then
	echo "jobs must be a positive integer: $JOBS" >&2
	exit 2
fi

ssh_options=()
if [[ -n "$IDENTITY_FILE" ]]; then
	if [[ ! -f "$IDENTITY_FILE" ]]; then
		echo "SSH private key does not exist: $IDENTITY_FILE" >&2
		exit 1
	fi
	key_permissions=$(stat -c '%a' "$IDENTITY_FILE")
	if [[ "$key_permissions" != "400" && "$key_permissions" != "600" ]]; then
		echo "SSH private key must have permissions 400 or 600: $IDENTITY_FILE ($key_permissions)" >&2
		exit 1
	fi
	ssh_options+=(-i "$IDENTITY_FILE" -o IdentitiesOnly=yes)
fi

for command_name in ssh rsync zstd sha256sum stat awk sort; do
	command -v "$command_name" >/dev/null || {
		echo "missing required command: $command_name" >&2
		exit 127
	}
done

owns_control_socket=0
control_dir=""
if [[ -n "$CONTROL_PATH" ]]; then
	ssh_options+=(-o BatchMode=yes -o ControlMaster=no -o "ControlPath=$CONTROL_PATH")
else
	control_dir=$(mktemp -d)
	CONTROL_PATH="$control_dir/ssh-%C"
	ssh_options+=(-o GSSAPIAuthentication=yes -o ControlMaster=auto -o ControlPersist=yes -o "ControlPath=$CONTROL_PATH")
	owns_control_socket=1
fi
printf -v rsync_rsh '%q ' ssh "${ssh_options[@]}"

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
COMPRESS_SCRIPT="$SCRIPT_DIR/compress_raw_zstd.sh"
if [[ ! -x "$COMPRESS_SCRIPT" ]]; then
	echo "compression script is not executable: $COMPRESS_SCRIPT" >&2
	exit 1
fi

mkdir -p "$DATA_DIR" "$OUT_DIR"
manifest=$(mktemp)
cleanup() {
	if ((owns_control_socket)); then
		ssh "${ssh_options[@]}" -O exit "$REMOTE_HOST" >/dev/null 2>&1 || true
	fi
	rm -f "$manifest"
	if [[ -n "$control_dir" ]]; then
		rmdir "$control_dir" >/dev/null 2>&1 || true
	fi
}
trap cleanup EXIT
printf -v remote_dir_quoted '%q' "$REMOTE_DIR"
echo "Discovering .raw files on $REMOTE_HOST:$REMOTE_DIR"
ssh "${ssh_options[@]}" "$REMOTE_HOST" "find $remote_dir_quoted -maxdepth 1 -type f -name '*.raw' -printf '%f\\n' | sort" > "$manifest"

mapfile -t raw_names < "$manifest"
if ((${#raw_names[@]} == 0)); then
	echo "No remote .raw files found"
	exit 0
fi

metadata_value() {
	local metadata=$1
	local field=$2
	awk -F'  ' -v field="$field" '$2 == field {print $1}' "$metadata"
}

verify_compressed_output() {
	local compressed=$1
	local metadata="$compressed.sha256"
	[[ -s "$compressed" && -s "$metadata" ]] || return 1

	local original_hash decompressed_hash stored_compressed_size actual_compressed_size actual_decompressed_hash
	original_hash=$(metadata_value "$metadata" original)
	decompressed_hash=$(metadata_value "$metadata" decompressed)
	stored_compressed_size=$(metadata_value "$metadata" compressed_size)
	[[ -n "$original_hash" && "$original_hash" == "$decompressed_hash" ]] || return 1

	actual_compressed_size=$(stat -c '%s' "$compressed")
	[[ "$stored_compressed_size" == "$actual_compressed_size" ]] || return 1
	zstd -t -- "$compressed" >/dev/null
	actual_decompressed_hash=$(zstd -dc -- "$compressed" | sha256sum | awk '{print $1}')
	[[ "$actual_decompressed_hash" == "$decompressed_hash" ]]
}

total=${#raw_names[@]}
process_file() {
	local index=$1
	local base=$2
	if [[ -z "$base" || "$base" == */* || "$base" != *.raw || "$base" == ._* ]]; then
		echo "Ignoring unexpected remote filename: $base" >&2
		return 0
	fi

	local raw="$DATA_DIR/$base"
	local partial="$raw.part"
	local compressed="$OUT_DIR/$base.zst"
	local metadata="$compressed.sha256"
	local stored_hash local_hash
	echo "[$((index + 1))/$total] $base"

	if verify_compressed_output "$compressed"; then
		if [[ -f "$raw" ]]; then
			stored_hash=$(metadata_value "$metadata" original)
			local_hash=$(sha256sum "$raw" | awk '{print $1}')
			if [[ "$local_hash" == "$stored_hash" ]]; then
				rm -f -- "$raw" "$partial"
				echo "  verified compressed output; removed matching local raw"
				return 0
			fi
			echo "  local raw differs from existing output; recompressing local raw"
		else
			rm -f -- "$partial"
			echo "  verified compressed output already exists; skipping download"
			return 0
		fi
	fi

	if [[ ! -f "$raw" ]]; then
		if [[ -f "$partial" ]]; then
			echo "  resuming incomplete download: $partial ($(stat -c '%s' "$partial") bytes)"
		fi
		echo "  downloading $REMOTE_HOST:$REMOTE_DIR/$base"
		local attempt download_ok=0
		for attempt in 1 2 3; do
			if rsync --partial --append-verify --compress --protect-args -e "$rsync_rsh" -- "$REMOTE_HOST:$REMOTE_DIR/$base" "$partial"; then
				download_ok=1
				break
			fi
			echo "  WARNING: download attempt $attempt/3 failed for $base" >&2
		done
		if ((download_ok == 0)); then
			echo "  ERROR: download failed after 3 attempts for $base; keeping $partial" >&2
			return 1
		fi
		mv -f -- "$partial" "$raw"
	else
		echo "  using existing local raw: $raw"
	fi

	"$COMPRESS_SCRIPT" --input "$raw" --out-dir "$OUT_DIR" --level "$LEVEL" --threads "$THREADS" --force
	if ! verify_compressed_output "$compressed"; then
		echo "  ERROR: compressed output verification failed; keeping $raw" >&2
		exit 1
	fi
	stored_hash=$(metadata_value "$metadata" original)
	local_hash=$(sha256sum "$raw" | awk '{print $1}')
	if [[ "$local_hash" != "$stored_hash" ]]; then
		echo "  ERROR: raw checksum does not match metadata; keeping $raw" >&2
		exit 1
	fi

	rm -f -- "$raw"
	echo "  verified and removed local raw"
}

pids=()
wait_for_workers() {
	local failed=0
	local pid
	for pid in "${pids[@]}"; do
		if ! wait "$pid"; then
			failed=1
		fi
	done
	pids=()
	return "$failed"
}

for index in "${!raw_names[@]}"; do
	process_file "$index" "${raw_names[index]}" &
	pids+=("$!")
	if ((${#pids[@]} >= JOBS)); then
		if ! wait_for_workers; then
			echo "ERROR: one or more workers failed; unverified raw files were kept" >&2
			exit 1
		fi
	fi
done
if ! wait_for_workers; then
	echo "ERROR: one or more workers failed; unverified raw files were kept" >&2
	exit 1
fi

echo "Done. Processed or verified $total remote file(s) with up to $JOBS worker(s)."
echo "Compressed files are under $OUT_DIR"
