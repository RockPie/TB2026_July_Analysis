# Snakemake workflow

Run from the repository root:

```sh
snakemake -s workflow/Snakefile --cores 1 --rerun-incomplete
```

The workflow discovers all real `data/*.raw` files and ignores macOS `._*.raw` resource files.

Rootify every nonempty `data/zstd/*.raw.zst` file that has a checksum sidecar:

```sh
snakemake -s workflow/Snakefile --cores 6 --rerun-incomplete rootify_all_verified
```

Before each rootifier job, `verify_zstd_raw` checks zstd integrity, the stored
compressed size, nonempty original size, and SHA-256 of the decompressed stream. Failed files do not
reach `001_Rootifier`; verification logs are written under
`dump/logs/verify_zstd/`. Start a new Snakemake invocation to discover sidecars
created after the workflow began.

Outputs per raw file:

- `dump/001/<raw filename>/pt_line_type_distributions.pdf`
- `dump/001/<raw filename>/sample_events.root`
- `dump/002/<raw filename>/adc_val0_accumulated_waveform.pdf`
- `dump/002/<raw filename>/event_adc_sum_histograms.root`
- `dump/003/adc_sum_comparison.pdf`
- logs under `dump/logs/001/`, `dump/logs/002/`, and `dump/logs/003/`

Optional config values:

```sh
snakemake -s workflow/Snakefile --cores 1 --rerun-incomplete --config size_limit=1GB crc=true
```

By default, `size_limit` is empty, so the full raw file is processed, and `crc=true` adds `--crc` to `001_Rootifier`.

`--rerun-incomplete` is recommended because long ROOT jobs may leave Snakemake incomplete metadata if a run is interrupted.

003 can also compare selected runs from a JSON config, including parameter-trend pages:

```sh
build/bin/003_ADC_Sum_Comparison --config config/temperature_adc_sum_comparison_2026_07_05.json
```

This writes `dump/003/temperature_adc_sum_comparison.pdf`.

The workflow also has a config-driven 003 rule that reads the same JSON and tracks each listed raw file through the downstream 002 histogram output:

```sh
snakemake -s workflow/Snakefile --cores 1 --rerun-incomplete dump/003/temperature_adc_sum_comparison.pdf
```

Compress raw data files with zstd and verify byte-for-byte decompression:

```sh
workflow/compress_raw_zstd.sh
```

By default this compresses `data/*.raw` to `data/zstd/*.raw.zst` with `zstd -T0 -3`, skips already verified outputs, and writes a `.sha256` sidecar for each compressed file. Use `--force` to recompress or `--level N` to choose another zstd level.

Download all CERN EOS raw files one at a time, compress and verify each file,
then remove its local raw copy before downloading the next file:

```sh
workflow/download_compress_raw.sh
```

The default source is
`sjia@lxplus.cern.ch:/eos/experiment/alice/focal/TB_2026_Wk27_H4/data/focalh-nas-1/Data`.
Existing verified files under `data/zstd` are skipped, complete `data/*.raw`
files are reused, and interrupted `.part` downloads are resumed on the next
run. Override connection or output settings with `--host`, `--remote-dir`,
`--data-dir`, and `--out-dir`; run with `--help` for all options.

For unattended lxplus authentication, obtain a CERN Kerberos ticket first:

```sh
kinit sjia@CERN.CH
klist
workflow/download_compress_raw.sh
```

`kinit` asks for the CERN password once and stores a temporary ticket, not the
password. If Kerberos is unavailable, run the downloader directly and enter
the lxplus password once when prompted. The downloader keeps that first SSH
connection open and reuses it for every resumable `rsync`, then closes it on exit. Do not
store the account password in a text file or pass it through `sshpass`.

To run six downloads/compressions in parallel under `nohup`, first create and
authenticate a persistent SSH connection in the foreground, then launch the
downloader using its control socket:

```sh
mkdir -p dump/logs
CONTROL_DIR=$(mktemp -d /tmp/lxplus-download.XXXXXX)
CONTROL_SOCKET="$CONTROL_DIR/control"
printf '%s\n' "$CONTROL_SOCKET" > dump/logs/download_compress_raw.control
ssh -M -S "$CONTROL_SOCKET" -o ControlPersist=yes -fN sjia@lxplus.cern.ch
nohup workflow/download_compress_raw.sh --control-path "$CONTROL_SOCKET" --jobs 6 --threads 1 > dump/logs/download_compress_raw.log 2>&1 &
echo $! > dump/logs/download_compress_raw.pid
```

The `ssh` command prompts for the lxplus password once before returning. Check
progress with `tail -f dump/logs/download_compress_raw.log`. The downloader
discovers the complete remote file list and skips locally verified outputs.
After it finishes, close the persistent connection and remove its directory:

```sh
CONTROL_SOCKET=$(cat dump/logs/download_compress_raw.control)
ssh -S "$CONTROL_SOCKET" -O exit sjia@lxplus.cern.ch
rm -rf -- "$(dirname -- "$CONTROL_SOCKET")"
```
