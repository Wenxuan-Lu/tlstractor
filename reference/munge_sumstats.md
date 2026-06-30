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
  n_col = "N",
  af_col = "AF",
  n_case_col = "N_case",
  n_control_col = "N_control",
  af_case_col = "AF_case",
  af_control_col = "AF_control",
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

- n_col:

  Character; total sample size column name in the GWAS summary
  statistics file (default `"N"`). Optional; omitted from output if
  absent in input GWAS summary statistics.

- af_col:

  Character; allele frequency column name in the GWAS summary statistics
  file (default `"AF"`). Optional; omitted from output if absent in
  input GWAS summary statistics.

- n_case_col:

  Character; case sample size column name in the GWAS summary statistics
  file (default `"N_case"`). Optional; omitted from output if absent in
  input GWAS summary statistics.

- n_control_col:

  Character; control sample size column name in the GWAS summary
  statistics file (default `"N_control"`). Optional; omitted from output
  if absent in input GWAS summary statistics.

- af_case_col:

  Character; case allele frequency column name in the GWAS summary
  statistics file (default `"AF_case"`). Optional; omitted from output
  if absent in input GWAS summary statistics.

- af_control_col:

  Character; control allele frequency column name in the GWAS summary
  statistics file (default `"AF_control"`). Optional; omitted from
  output if absent in input GWAS summary statistics.

- remove_ambiguous:

  Logical; if `TRUE` (default), exclude ambiguous SNPs.

- output_path:

  Character or `NULL`; output file path. If `NULL` (default),
  auto-generated as `<sumstats_prefix>_munged.txt` in the same directory
  as `sumstats_path`, where `<sumstats_prefix>` is the filename prefix
  of `sumstats_path` without extension.

## Value

Invisibly returns `NULL`. Writes filtered results to `output_path` as
tab-delimited text. Output always includes `CHR`, `POS`, `ID`, `REF`,
`ALT`, `BETA`, `SE`, and `GDS_ID`. Optional input columns `N`, `AF`,
`N_case`, `N_control`, `AF_case`, and `AF_control` are preserved if
present. `GDS_ID` is the 1-based variant index in the input GDS file to
match variants in the summary statistics with variants in the GDS file.

## Details

**Required columns:** The input summary statistics must contain columns
for `REF`, `ALT`, `BETA`, and `SE`, plus either (`CHR` and `POS`) when
`match_by = "CHR-POS"` or `ID` when `match_by = "ID"`. Optional columns
are `N`, `AF`, `N_case`, `N_control`, `AF_case`, and `AF_control`.
Additional columns are allowed but ignored. When `match_by = "CHR-POS"`,
accepted chromosome formats are autosomes only: `1-22`, `01-22`,
`chr1-chr22`, or `chr01-chr22` (case-insensitive). The `BETA` column
should represent the effect size estimate for the alternative allele
(log odds ratio for logistic regression or linear coefficient for linear
regression), and the `SE` column should contain the corresponding
standard error.

**Quality control filters for summary statistics (applied in order):**

1.  REF/ALT: single A/C/G/T nucleotides, biallelic

2.  Ambiguous SNPs: optionally remove A/T, T/A, C/G, and G/C

3.  Autosomes: restrict to chromosomes 1-22 when `match_by = "CHR-POS"`

4.  POS: require `POS > 0` when `match_by = "CHR-POS"`

5.  ID: require non-missing, non-empty IDs when `match_by = "ID"`

6.  Duplicates: keep variants with unique `(CHR, POS)` or unique `ID`
    and remove all duplicates

7.  Effect/SE: require `BETA` not missing and `SE > 0`

8.  Sample sizes and AF ranges: require `N`, `N_case`, and `N_control`
    to be positive when present, and require `AF`, `AF_case`, and
    `AF_control` to lie in `(0, 1)` when present
