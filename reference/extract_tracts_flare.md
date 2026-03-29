# Extract Local-Ancestry Tracts from a FLARE VCF file

Parses a FLARE VCF and writes ancestry-specific outputs in one or more
formats: `gds`, `txt`, `txt.gz`, `vcf`, `vcf.gz`.

## Usage

``` r
extract_tracts_flare(
  vcf_path,
  num_ancs,
  output_dir = NULL,
  output_formats = "gds",
  chunk_size = 1024L
)
```

## Arguments

- vcf_path:

  Character; input FLARE VCF path (`.vcf` or `.vcf.gz`).

- num_ancs:

  Integer; number of ancestral populations (must be \> 1).

- output_dir:

  Character or `NULL`; output directory (must exist). If `NULL`
  (default), uses the directory of `vcf_path`.

- output_formats:

  Character vector; output format specification (for example `"gds"` or
  `c("gds", "txt.gz")`). Default is `"gds"`. Supported values: `gds`,
  `txt`, `txt.gz`, `vcf`, `vcf.gz`. If both compressed and uncompressed
  versions are requested for the same type, compressed output is
  written: `txt.gz` over `txt`, and `vcf.gz` over `vcf`.

- chunk_size:

  Integer; number of variants processed per chunk. Default is 1024;
  higher values can improve speed but increase memory usage.

## Value

Invisibly returns `NULL`. Output files are written with filename prefix
derived from `vcf_path` in `output_dir`.

## Details

Output files (prefix = basename of `vcf_path` without `.vcf` or
`.vcf.gz`):

- `gds`: one file `<prefix>.gds` containing `sample.id`, variant
  metadata (`snp.chromosome`, `snp.position`, `snp.id`, `snp.ref`,
  `snp.alt`), and ancestry-specific nodes `dosage/anc0..ancK` and
  `hapcount/anc0..ancK`.

- `txt` / `txt.gz`: for each ancestry `k`, two files
  `<prefix>.anc{k}.dosage.txt(.gz)` and
  `<prefix>.anc{k}.hapcount.txt(.gz)`. Columns are
  `CHROM POS ID REF ALT` followed by one column per sample.

- `vcf` / `vcf.gz`: for each ancestry `k`, one file
  `<prefix>.anc{k}.vcf(.gz)` with `FORMAT=GT`; haplotypes not assigned
  to ancestry `k` are written as `.`.

Assumptions for the input FLARE VCF:

- Biallelic variants.

- No missing genotype/local-ancestry fields.

- No duplicate variants.

- Autosomal chromosomes only (1-22, 01-22, `chr1`-`chr22`,
  `chr01`-`chr22`; case-insensitive).

- Alleles are uppercase `A/C/G/T`.

- Sample FORMAT is exactly `GT:AN1:AN2`.

- `GT` is phased (`0|0`, `0|1`, `1|0`, `1|1`) with single-digit allele
  codes.
