# Munge GWAS Summary Statistics to Align with GDS File

Aligns GWAS summary statistics with a GDS file using chromosome-position
or variant ID matching. Performs quality control filtering, allele
flipping for discordant variants, and outputs filtered results.

## Usage

``` r
munge_sumstats(
  gds_path,
  sumstats_path,
  match_by = "CHR-POS",
  chr_col = "CHR",
  pos_col = "POS",
  id_col = "ID",
  ref_col = "REF",
  alt_col = "ALT",
  beta_col = "BETA",
  se_col = "SE",
  af_col = "AF",
  remove_ambiguous = TRUE,
  output_path = NULL
)
```

## Arguments

- gds_path:

  Character; path to the input GDS file.

- sumstats_path:

  Character; path to the GWAS summary statistics file.

- match_by:

  Character; `"CHR-POS"` (default) or `"ID"`. `"CHR-POS"`: match by
  chromosome and position (autosomes 1-22 only). `"ID"`: match by
  variant ID.

- chr_col:

  Character; chromosome column name in the GWAS summary statistics file
  (default `"CHR"`).

- pos_col:

  Character; position column name in the GWAS summary statistics file
  (default `"POS"`).

- id_col:

  Character; variant ID column name in the GWAS summary statistics file
  (default `"ID"`).

- ref_col:

  Character; reference allele column name in the GWAS summary statistics
  file (default `"REF"`).

- alt_col:

  Character; alternative allele column name in the GWAS summary
  statistics file (default `"ALT"`).

- beta_col:

  Character; effect size column name in the GWAS summary statistics file
  (default `"BETA"`).

- se_col:

  Character; standard error column name in the GWAS summary statistics
  file (default `"SE"`).

- af_col:

  Character; allele frequency column name in the GWAS summary statistics
  file (default `"AF"`). Optional; omitted from output if absent in
  input GWAS summary statistics.

- remove_ambiguous:

  Logical; if `TRUE` (default), exclude ambiguous SNPs.

- output_path:

  Character or `NULL`; output file path. If `NULL` (default),
  auto-generated as `<sumstats_prefix>_munged.txt` in the same directory
  as `sumstats_path`, where `<sumstats_prefix>` is the filename prefix
  of `sumstats_path` without extension.

## Value

Invisibly returns `NULL`. Filtered results written to `output_path` as
tab-delimited text with columns: `CHR`, `POS`, `ID`, `REF`, `ALT`,
`BETA`, `SE`, `AF` (optional), `GDS_ID`. `GDS_ID` is the 1-based variant
index in the input GDS file to match variants in the summary statistics
with variants in the GDS file.

## Details

**Required columns:** The input summary statistics must contain columns
for `REF`, `ALT`, `BETA`, and `SE`, along with either (`CHR` and `POS`)
when `match_by = "CHR-POS"` or `ID` when `match_by = "ID"`. The `AF`
column is optional. Additional columns are allowed but ignored. When
`match_by = "CHR-POS"`, accepted chromosome formats are autosomes only:
`1-22`, `01-22`, `chr1-chr22`, or `chr01-chr22` (case-insensitive). The
`BETA` column should represent the effect size estimate for the
alternative allele (log odds ratio for logistic regression or linear
coefficient for linear regression), and the `SE` column should contain
the corresponding standard error.

**Quality control filters for summary statistics (applied in order):**

1.  REF/ALT: single A/C/G/T nucleotides, biallelic

2.  Ambiguous SNPs: optionally remove A/T, T/A, C/G, G/C

3.  Autosomes: restrict to chromosomes 1-22 (when
    `match_by = "CHR-POS"`)

4.  POS: `POS > 0` (when `match_by = "CHR-POS"`)

5.  ID: non-missing and non-empty (when `match_by = "ID"`)

6.  Duplicates: retain variants with unique `(CHR, POS)` pairs (when
    `match_by = "CHR-POS"`) or unique `ID` (when `match_by = "ID"`)

7.  Effect/SE: `BETA != NA`, `SE > 0`

8.  AF range: `0 < AF < 1` (if present)
