# Snakemake workflow

Run from the repository root:

```sh
snakemake -s workflow/Snakefile --cores 1 --rerun-incomplete
```

The workflow discovers all real `data/*.raw` files and ignores macOS `._*.raw` resource files.

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
