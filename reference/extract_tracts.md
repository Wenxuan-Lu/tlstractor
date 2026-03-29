# Extract Local-Ancestry Tracts from Phased VCF and MSP (RFMix2/Gnomix)

Parses a phased VCF together with an MSP local-ancestry file and writes
ancestry-specific outputs in one or more formats: `gds`, `txt`,
`txt.gz`, `vcf`, `vcf.gz`.

## Usage

``` r
extract_tracts(
  vcf_path,
  msp_path,
  num_ancs,
  output_dir = NULL,
  output_formats = "gds",
  chunk_size = 1024L
)
```

## Arguments

- vcf_path:

  Character; input phased VCF path (`.vcf` or `.vcf.gz`).

- msp_path:

  Character; input MSP path.

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

Assumptions for inputs:

- Biallelic variants.

- No missing genotype/local-ancestry fields.

- No duplicate variants.

- Autosomal chromosomes only (1-22, 01-22, `chr1`-`chr22`,
  `chr01`-`chr22`; case-insensitive).

- Alleles are uppercase `A/C/G/T`.

- `GT` is the first sample subfield (`GT:...`).

- `GT` is phased (`0|0`, `0|1`, `1|0`, `1|1`) with single-digit allele
  codes.

- One chromosome per VCF and per MSP file.

- Sample order is consistent between VCF and MSP files.
