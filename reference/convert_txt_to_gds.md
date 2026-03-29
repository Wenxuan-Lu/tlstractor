# Convert Local-Ancestry Tracts Format from Text to GDS

Converts ancestry-specific dosage/hapcount text files into a tracts GDS
file.

## Usage

``` r
convert_txt_to_gds(
  input_prefix,
  num_ancs,
  output_prefix = NULL,
  chunk_size = 1024L
)
```

## Arguments

- input_prefix:

  Character; input file prefix used to locate ancestry files. Expected
  files per ancestry `k`: `<input_prefix>.anc{k}.dosage.txt(.gz)` and
  `<input_prefix>.anc{k}.hapcount.txt(.gz)`.

- num_ancs:

  Integer; number of ancestral populations (must be \> 1).

- output_prefix:

  Character or `NULL`; output GDS prefix. If `NULL` (default), uses
  `input_prefix`.

- chunk_size:

  Integer; number of variants processed per chunk. Default is 1024;
  higher values can improve speed but increase memory usage.

## Value

Invisibly returns `NULL`. Writes `<output_prefix>.gds`.

## Details

The generated GDS contains `sample.id`, variant metadata
(`snp.chromosome`, `snp.position`, `snp.id`, `snp.ref`, `snp.alt`), and
ancestry-specific nodes `dosage/anc0..ancK` and `hapcount/anc0..ancK`.
