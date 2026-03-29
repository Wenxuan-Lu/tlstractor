# Convert Local-Ancestry Tracts Format from GDS to Text

Converts ancestry-specific dosage and hapcount matrices stored in a
tracts GDS file into tab-delimited text files (plain or gzipped).

## Usage

``` r
convert_gds_to_txt(
  gds_path,
  output_prefix = NULL,
  output_format = "txt.gz",
  chunk_size = 1024L
)
```

## Arguments

- gds_path:

  Character; input tracts GDS path.

- output_prefix:

  Character or `NULL`; output file prefix. If `NULL` (default), uses the
  input filename prefix in the same directory as `gds_path`.

- output_format:

  Character; output format, one of `"txt"` or `"txt.gz"` (default).

- chunk_size:

  Integer; number of variants processed per chunk. Default is 1024;
  higher values can improve speed but increase memory usage.

## Value

Invisibly returns `NULL`. For each ancestry `k`, writes two files:
`<output_prefix>.anc{k}.dosage.txt(.gz)` and
`<output_prefix>.anc{k}.hapcount.txt(.gz)`.

## Details

Output tables contain columns `CHROM`, `POS`, `ID`, `REF`, `ALT`,
followed by one column per sample from `sample.id` in the input GDS.
