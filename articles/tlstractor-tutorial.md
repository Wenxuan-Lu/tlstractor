# TLS-Tractor Tutorial

``` r
library(tlstractor)
```

## Overview

This tutorial goes through the end-to-end pipeline to:

1.  extract tracts and ancestry dosages
2.  munge external GWAS summary statistics
3.  run TLS-Tractor

Below, we walk through a step-by-step example for a continuous
phenotype, then provide a runnable short script for a binary phenotype.

## Simulated data in `inst/extdata/`

The examples in this tutorial use simulated files bundled with the
package under `inst/extdata/`.

The simulated data include genotype, covariates, and phenotype for
unrelated two-way admixed individuals. Individuals are split into
non-overlapping internal (`n = 100`) and external (`n = 8000`) cohorts.
In the external cohort, standard GWAS models
(`y ~ total genotype dosage + covariates`) are used to generate external
GWAS summary statistics. The internal cohort provides individual-level
covariates, phenotype data, and local ancestry inference (LAI) results
in both RFMix2/Gnomix and FLARE formats.

``` r
# Get data directory path
data_dir <- system.file("extdata", package = "tlstractor")

# To see all available simulated files, use:
list.files(data_dir, full.names = FALSE)
#>  [1] "simulate_binary.covar.txt"                 
#>  [2] "simulate_binary.external_gwas_sumstats.txt"
#>  [3] "simulate_binary.flare.vcf"                 
#>  [4] "simulate_binary.internal_data.txt"         
#>  [5] "simulate_binary.meta.txt"                  
#>  [6] "simulate_binary.msp.tsv"                   
#>  [7] "simulate_binary.phased.vcf"                
#>  [8] "simulate_binary.pheno.txt"                 
#>  [9] "simulate_linear.covar.txt"                 
#> [10] "simulate_linear.external_gwas_sumstats.txt"
#> [11] "simulate_linear.flare.vcf"                 
#> [12] "simulate_linear.internal_data.txt"         
#> [13] "simulate_linear.meta.txt"                  
#> [14] "simulate_linear.msp.tsv"                   
#> [15] "simulate_linear.phased.vcf"                
#> [16] "simulate_linear.pheno.txt"
```

### Linear case inputs

The linear example (`simulate_linear.*`) includes:

1.  Phased genotype inputs and local ancestry inference outputs for the
    RFMix2/Gnomix style pipeline for the internal cohort:
    - `simulate_linear.phased.vcf`
    - `simulate_linear.msp.tsv`
    - Phased VCF stores genotypes in the `GT` field; MSP stores local
      ancestry labels for each haplotype window.
2.  Local ancestry inference output from FLARE for the internal cohort:
    - `simulate_linear.flare.vcf`
    - VCF with `GT:AN1:AN2` format where `GT` denotes phased genotypes
      and `AN1` and `AN2` encode local ancestry labels for each
      haplotype.
3.  Individual-level phenotype and covariates for the internal cohort:
    - `simulate_linear.pheno.txt`
    - `simulate_linear.covar.txt`
4.  External GWAS summary statistics:
    - `simulate_linear.external_gwas_sumstats.txt`
    - Columns include `CHR POS ID REF ALT BETA SE AF N`.

### Binary case inputs

The binary simulation has the same file layout with prefix
`simulate_binary`:

1.  `simulate_binary.phased.vcf`
2.  `simulate_binary.msp.tsv`
3.  `simulate_binary.flare.vcf`
4.  `simulate_binary.pheno.txt`
5.  `simulate_binary.covar.txt`
6.  `simulate_binary.external_gwas_sumstats.txt`

## Step 1: Extract tracts and convert formats

This step extracts local ancestry tracts from LAI results and converts
them into formats that are more flexible for downstream analysis. For
downstream steps in TLS-Tractor, **GDS is the only required format**.

First, create a directory to store tract outputs.

``` r
output_dir_tracts <- file.path(tempdir(), 
                               "tlstractor_tutorial", 
                               "linear", "tracts")
dir.create(output_dir_tracts, recursive = TRUE, showWarnings = FALSE)
```

Second, use **one** extractor based on your LAI format. If you have
RFMix2/Gnomix style phased VCF + MSP, use
[`extract_tracts()`](https://wenxuan-lu.github.io/tlstractor/reference/extract_tracts.md).
If you have FLARE VCF, use
[`extract_tracts_flare()`](https://wenxuan-lu.github.io/tlstractor/reference/extract_tracts_flare.md).

### Generate GDS file (required for downstream steps)

#### Option A: RFMix2/Gnomix style (phased VCF + MSP)

Use
**[`extract_tracts()`](https://wenxuan-lu.github.io/tlstractor/reference/extract_tracts.md)**
when your local ancestry inference output is a phased VCF plus an MSP
file.

- `vcf_path`: phased VCF input with phased genotypes (`GT` field
  required to be the first FORMAT field)
- `msp_path`: MSP input with local ancestry labels
- `num_ancs`: number of ancestral populations in the dataset
- `output_dir`: directory for output files (`NULL` uses the input VCF
  directory)
- `output_formats = "gds"`: write GDS output required for downstream
  TLS-Tractor steps
- `chunk_size`: number of variants processed per chunk (default: 1024L).
  Larger values can improve speed but increase memory usage.

The GDS output path is `<prefix>.gds` within `output_dir`, where
`prefix` is the basename of the input VCF file `vcf_path` without `.vcf`
or `.vcf.gz`.

``` r
extract_tracts(
  vcf_path = file.path(data_dir, "simulate_linear.phased.vcf"),
  msp_path = file.path(data_dir, "simulate_linear.msp.tsv"),
  num_ancs = 2L,
  output_dir = output_dir_tracts,
  output_formats = "gds"
)

# After running, the GDS file will be at:
gds_path <- file.path(output_dir_tracts, "simulate_linear.phased.gds")
```

#### Option B: FLARE VCF

Use
**[`extract_tracts_flare()`](https://wenxuan-lu.github.io/tlstractor/reference/extract_tracts_flare.md)**
when your local ancestry inference output is a FLARE VCF.

- `vcf_path`: FLARE VCF input with `GT:AN1:AN2` format
- `num_ancs`: number of ancestral populations in the dataset
- `output_dir`: directory for output files (`NULL` uses the input VCF
  directory)
- `output_formats = "gds"`: write GDS output required for downstream
  TLS-Tractor steps
- `chunk_size`: number of variants processed per chunk (default: 1024L).
  Larger values can improve speed but increase memory usage.

The GDS output path is `<prefix>.gds` within `output_dir`, where
`prefix` is the basename of the input VCF file `vcf_path` without `.vcf`
or `.vcf.gz`.

``` r
extract_tracts_flare(
  vcf_path = file.path(data_dir, "simulate_linear.flare.vcf"),
  num_ancs = 2L,
  output_dir = output_dir_tracts,
  output_formats = "gds"
)

# After running, the GDS file will be at:
gds_path <- file.path(output_dir_tracts, "simulate_linear.flare.gds")
```

#### Understanding the GDS output

After extraction, the GDS file stores local ancestry and risk allele
information in a hierarchical structure optimized for analysis. You can
think of it as a compact container with two types of content:

1.  sample IDs and SNP metadata that describe samples and variants
2.  ancestry-specific arrays for risk allele dosages and haplotype
    counts

In this structure,

- `dosage` contains the number of copies of the risk allele from each
  ancestry
- `hapcount` contains the number of haplotypes from each ancestry

``` r
# Open the GDS file in read-only mode
gds <- gdsfmt::openfn.gds(gds_path, readonly = TRUE)

# Display the GDS structure
print(gds)
#> File: /tmp/RtmphmsPgt/tlstractor_tutorial/linear/tracts/simulate_linear.flare.gds (3.5K)
#> +    [  ]
#> |--+ sample.id   { Str8 100 ZIP_ra(35.0%), 287B }
#> |--+ snp.chromosome   { Str8 5 ZIP_ra(104.0%), 33B }
#> |--+ snp.position   { Int32 5 ZIP_ra(170.0%), 41B }
#> |--+ snp.id   { Str8 5 ZIP_ra(72.0%), 43B }
#> |--+ snp.ref   { Str8 5 ZIP_ra(220.0%), 29B }
#> |--+ snp.alt   { Str8 5 ZIP_ra(220.0%), 29B }
#> |--+ dosage   [  ]
#> |  |--+ anc0   { Bit2 100x5 LZ4_ra(119.2%), 156B }
#> |  \--+ anc1   { Bit2 100x5 LZ4_ra(119.2%), 156B }
#> \--+ hapcount   [  ]
#>    |--+ anc0   { Bit2 100x5 LZ4_ra(119.2%), 156B }
#>    \--+ anc1   { Bit2 100x5 LZ4_ra(119.2%), 156B }

# Key variables in the GDS:
# - sample.id: individual identifiers
# - snp.chromosome / snp.position / snp.id / snp.ref / snp.alt: SNP metadata
# - dosage/anc0, dosage/anc1, ...: ancestry-specific risk-allele dosage
# - hapcount/anc0, hapcount/anc1, ...: ancestry-specific haplotype count

# List top-level nodes
gdsfmt::ls.gdsn(gds)
#> [1] "sample.id"      "snp.chromosome" "snp.position"   "snp.id"        
#> [5] "snp.ref"        "snp.alt"        "dosage"         "hapcount"

# Close the GDS connection
gdsfmt::closefn.gds(gds)
```

The chunk below demonstrates basic usage of GDS files: reading sample
IDs and SNP metadata, then extracting ancestry-specific dosage and
hapcount for selected SNPs or individuals.

``` r
# Open the GDS file in read-only mode
gds <- gdsfmt::openfn.gds(gds_path, readonly = TRUE)

# Read sample IDs
sample_id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id"))

# Read SNP metadata
snp_meta <- data.frame(
  snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id")),
  snp.chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.chromosome")),
  snp.position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.position")),
  snp.ref = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.ref")),
  snp.alt = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.alt")),
  stringsAsFactors = FALSE
)

# Access ancestry-specific dosage/hapcount nodes
dosage_anc0_node <- gdsfmt::index.gdsn(gds, "dosage/anc0")
hapcount_anc0_node <- gdsfmt::index.gdsn(gds, "hapcount/anc0")

# Check dimensions (should be n_samples x n_snps)
dosage_anc0_dim <- gdsfmt::objdesp.gdsn(dosage_anc0_node)$dim
hapcount_anc0_dim <- gdsfmt::objdesp.gdsn(hapcount_anc0_node)$dim

# Set safe example indices
n_sample <- length(sample_id)
n_snp <- nrow(snp_meta)
snp_idx <- min(2L, n_snp)
sample_idx <- min(5L, n_sample)

# 1) All samples for one SNP in ancestry 0
dosage_anc0_snp <- gdsfmt::read.gdsn(dosage_anc0_node, 
                                     start = c(1L, snp_idx), 
                                     count = c(n_sample, 1L)
                                     )
hapcount_anc0_snp <- gdsfmt::read.gdsn(hapcount_anc0_node, 
                                       start = c(1L, snp_idx), 
                                       count = c(n_sample, 1L)
                                       )

# 2) All SNPs for one individual in ancestry 0
dosage_anc0_sample <- gdsfmt::read.gdsn(dosage_anc0_node,
                                        start = c(sample_idx, 1L), 
                                        count = c(1L, n_snp)
                                        )
hapcount_anc0_sample <- gdsfmt::read.gdsn(hapcount_anc0_node, 
                                          start = c(sample_idx, 1L), 
                                          count = c(1L, n_snp)
                                          )

# 3) One sample x one SNP in ancestry 0
dosage_anc0_cell <- gdsfmt::read.gdsn(dosage_anc0_node, 
                                      start = c(sample_idx, snp_idx),
                                      count = c(1L, 1L)
                                      )
hapcount_anc0_cell <- gdsfmt::read.gdsn(hapcount_anc0_node, 
                                        start = c(sample_idx, snp_idx), 
                                        count = c(1L, 1L)
                                        )

# Access another ancestry by changing the node path, e.g. "dosage/anc1"
dosage_anc1_node <- gdsfmt::index.gdsn(gds, "dosage/anc1")

# Close the GDS connection
gdsfmt::closefn.gds(gds)
```

**Benefits of GDS format:**

- **Memory efficient**: Uses a compressed binary representation that is
  typically much smaller than `.vcf.gz` or `.txt.gz`, while enabling
  on-demand access without loading the full dataset into memory
- **Fast access**: Supports efficient sequential and indexed random
  access across both variants and samples
- **Portable**: Stores genotypes, local ancestry, and meta data in a
  single, self-contained file
- **Scalable**: Designed for large-scale genomic datasets, enabling
  efficient analysis of genome-wide data in large cohorts

### Additional output format options

The
[`extract_tracts()`](https://wenxuan-lu.github.io/tlstractor/reference/extract_tracts.md)
and
[`extract_tracts_flare()`](https://wenxuan-lu.github.io/tlstractor/reference/extract_tracts_flare.md)
functions support multiple output formats via the `output_formats`
argument. Though **GDS is the only format required for downstream
TLS-Tractor steps**, you can optionally export to text or VCF formats
for inspection, sharing, or use in external tools.

- **`"gds"`** (default): Genomic Data Structure format—binary,
  compressed, optimized for fast access and scalability (required for
  downstream TLS-Tractor steps).
- **`"txt"`** and **`"txt.gz"`**: Plain-text and gzip-compressed text
  format. For each ancestry, two files are created:
  `<prefix>.anc{k}.dosage.txt(.gz)` and
  `<prefix>.anc{k}.hapcount.txt(.gz)`, where `<prefix>` is the basename
  of the input VCF file `vcf_path` without `.vcf` or `.vcf.gz`. The
  former contains the number of copies of the risk allele from each
  ancestry, and the latter contains the number of haplotypes from each
  ancestry. Output files have columns `CHROM POS ID REF ALT` followed by
  one column per sample.
- **`"vcf"`** and **`"vcf.gz"`**: VCF format. For each ancestry, one
  file is created: `<prefix>.anc{k}.vcf(.gz)` with `FORMAT=GT`, where
  `<prefix>` is the basename of the input VCF file `vcf_path` without
  `.vcf` or `.vcf.gz`. Haplotypes not assigned to that ancestry are
  written as `.`.

You can request multiple formats in a single run by passing a vector
(e.g., `output_formats = c("gds", "txt.gz", "vcf.gz")`). Note that if
you request both compressed and uncompressed versions of the same format
(e.g., both `"txt"` and `"txt.gz"`), only the compressed version will be
written to disk to save space.

Example: extract to both GDS and compressed text in one call using the
FLARE VCF input.

``` r
extract_tracts_flare(
  vcf_path = file.path(data_dir, "simulate_linear.flare.vcf"),
  num_ancs = 2L,
  output_dir = output_dir_tracts,
  output_formats = c("gds", "txt.gz")  # Multiple formats in one call
)

# Outputs will be in output_dir_tracts/ with the filenames:
# - simulate_linear.flare.gds (GDS format)
# - simulate_linear.flare.anc0.dosage.txt.gz (dosage for ancestry 0)
# - simulate_linear.flare.anc0.hapcount.txt.gz (hapcount for ancestry 0)
# - simulate_linear.flare.anc1.dosage.txt.gz (dosage for ancestry 1)
# - simulate_linear.flare.anc1.hapcount.txt.gz (hapcount for ancestry 1)
```

### Format conversion utilities

After extracting tracts, you can convert between GDS and text formats
using the following functions:

- **[`convert_txt_to_gds()`](https://wenxuan-lu.github.io/tlstractor/reference/convert_txt_to_gds.md)**:
  Converts text format (from a previous extraction or external source)
  to GDS format. Useful for reprocessing data or constructing a GDS from
  manually formatted tract files.
  - `input_prefix`: prefix of input text files, without
    `.anc{k}.dosage.txt(.gz)` or `.anc{k}.hapcount.txt(.gz)`
  - `num_ancs`: number of ancestral populations in the dataset
  - `output_prefix`: prefix for the output `.gds` file, without
    extension (`NULL` uses `input_prefix`)
  - `chunk_size`: number of variants processed per chunk (default:
    1024L). Larger values can improve speed but increase memory usage.
  - The GDS output will be `<output_prefix>.gds`.
- **[`convert_gds_to_txt()`](https://wenxuan-lu.github.io/tlstractor/reference/convert_gds_to_txt.md)**:
  Exports GDS data to plain-text format for inspection, sharing, or use
  in external tools.
  - `gds_path`: input GDS file path
  - `output_prefix`: prefix for output text files (`NULL` uses the input
    GDS filename prefix in the same directory)
  - `output_format`: output text format, either `"txt"` or `"txt.gz"`
    (default: `"txt.gz"`)
  - `chunk_size`: number of variants processed per chunk (default:
    1024L). Larger values can improve speed but increase memory usage.
  - The output will include two files per ancestry:
    `<output_prefix>.anc{k}.dosage.txt(.gz)` and
    `<output_prefix>.anc{k}.hapcount.txt(.gz)`.

``` r
# TXT (or txt.gz) -> GDS
convert_txt_to_gds(
  input_prefix = file.path(output_dir_tracts, "simulate_linear.flare"),
  num_ancs = 2L,
  output_prefix = file.path(output_dir_tracts, "simulate_linear_from_txt")
)
# Output: output_dir_tracts/simulate_linear_from_txt.gds

# GDS -> TXT (or txt.gz)
convert_gds_to_txt(
  gds_path = file.path(output_dir_tracts, "simulate_linear.flare.gds"),
  output_prefix = file.path(output_dir_tracts, "simulate_linear_from_gds"),
  output_format = "txt.gz"
)
# Outputs will be in output_dir_tracts/ with the filenames:
# - simulate_linear_from_gds.anc0.dosage.txt.gz
# - simulate_linear_from_gds.anc0.hapcount.txt.gz
# - simulate_linear_from_gds.anc1.dosage.txt.gz
# - simulate_linear_from_gds.anc1.hapcount.txt.gz
```

### Read function documentation

Use `?` to inspect documentation quickly.

``` r
?tlstractor::extract_tracts
?tlstractor::extract_tracts_flare
?tlstractor::convert_txt_to_gds
?tlstractor::convert_gds_to_txt
```

## Step 2: Munge external GWAS summary statistics

After creating the GDS file from individual-level data in the
main/internal study, use
[`munge_sumstats()`](https://wenxuan-lu.github.io/tlstractor/reference/munge_sumstats.md)
to harmonize external GWAS summary statistics with the GDS file. The
function performs quality control (QC), aligns variants between the
summary statistics and the GDS file using either `CHR-POS` or `ID`, and
outputs munged summary statistics restricted to SNPs present in both
datasets.

Expected input format for GWAS summary statistics:

- The input file should be tab-delimited.
- Required columns for `match_by = "CHR-POS"`: chromosome, position,
  reference allele, alternative allele, effect size, and standard error
  (default column names: `CHR`, `POS`, `REF`, `ALT`, `BETA`, `SE`).
- Required columns for `match_by = "ID"`: variant ID, reference allele,
  alternative allele, effect size, and standard error (default column
  names: `ID`, `REF`, `ALT`, `BETA`, `SE`).
- Optional column: allele frequency (default `AF`), if available.
- Additional columns are allowed but will be ignored in the munging
  process.
- For `match_by = "CHR-POS"`, accepted chromosome formats are autosomes
  only: `1-22`, `01-22`, `chr1-chr22`, or `chr01-chr22`
  (case-insensitive).
- The `BETA` column should represent the effect size estimate for the
  alternative allele (log odds ratio for logistic regression or linear
  coefficient for linear regression), and the `SE` column should contain
  the corresponding standard error.

If the column names in your summary statistics file differ from the
defaults, you can specify the column names using the function arguments:
`chr_col`, `pos_col`, `id_col`, `ref_col`, `alt_col`, `beta_col`,
`se_col`, and `af_col`.

QC performed by
[`munge_sumstats()`](https://wenxuan-lu.github.io/tlstractor/reference/munge_sumstats.md)
on GWAS summary statistics (applied in order):

1.  REF/ALT checks: keep biallelic SNPs with single-nucleotide `A/C/G/T`
    alleles and remove the rest
2.  optional removal of ambiguous SNPs (`A/T`, `T/A`, `C/G`, `G/C`)
3.  autosome restriction to chromosomes `1-22` when
    `match_by = "CHR-POS"`
4.  `POS > 0` when `match_by = "CHR-POS"`
5.  non-missing and non-empty `ID` when `match_by = "ID"`
6.  duplicate removal by key: unique `(CHR, POS)` when
    `match_by = "CHR-POS"`, or unique `ID` when `match_by = "ID"`
7.  effect/SE checks: `BETA != NA` and `SE > 0`
8.  optional `AF` range check (`0 < AF < 1`) when `AF` is provided
9.  allele flipping when needed to align with GDS reference allele
    definition

Output format for munged summary statistics:

- The output file is tab-delimited.
- Output columns are standardized as `CHR`, `POS`, `ID`, `REF`, `ALT`,
  `BETA`, `SE`, `AF` (optional), and `GDS_ID`.
- The output is aligned to GDS metadata (chromosome, position, variant
  ID, and alleles) so variant definitions are consistent with the GDS
  file.
- `GDS_ID` is the 1-based variant index in the GDS file to match
  variants in the summary statistics with variants in the GDS file. It
  is used by
  [`tlstractor()`](https://wenxuan-lu.github.io/tlstractor/reference/tlstractor.md)
  in Step 3.

In the example call below:

- `gds_path`: path to the GDS file generated in Step 1
- `sumstats_path`: path to the external GWAS summary statistics file
- `match_by`: method to align variants between summary statistics and
  GDS file, either `"CHR-POS"` (default) or `"ID"`
- `remove_ambiguous`: whether to remove ambiguous SNPs (default: `TRUE`)
- `output_path`: output path for munged summary statistics (`NULL` uses
  `<sumstats_prefix>_munged.txt` in the same directory as
  `sumstats_path`, where `<sumstats_prefix>` is the input summary
  statistics filename without extension)

``` r
# Create output directory for munged summary statistics
output_dir_munge <- file.path(tempdir(), 
                              "tlstractor_tutorial",
                              "linear", "munged_sumstats")
dir.create(output_dir_munge, recursive = TRUE, showWarnings = FALSE)

# Define output path for munged summary statistics
munged_sumstats_path <- file.path(output_dir_munge,
                                  "external_gwas_sumstats_munged.txt")

# Munge summary statistics
munge_sumstats(
  gds_path = gds_path,
  sumstats_path = file.path(data_dir,
                            "simulate_linear.external_gwas_sumstats.txt"),
  match_by = "CHR-POS",
  remove_ambiguous = TRUE,
  output_path = munged_sumstats_path
)
#> 0 rows excluded: invalid alleles (not biallelic SNPs).
#> 0 rows excluded: ambiguous SNPs (A/T, T/A, C/G, G/C).
#> 0 rows excluded: non-autosomal or invalid chromosomes.
#> 0 rows excluded: invalid POS values.
#> 0 rows excluded: duplicate (CHR, POS) pairs.
#> 0 rows excluded: invalid BETA or SE values.
#> 0 rows excluded: invalid AF values.
#> Variants retained after QC and matching: 5
#> Munged summary statistics written to: /tmp/RtmphmsPgt/tlstractor_tutorial/linear/munged_sumstats/external_gwas_sumstats_munged.txt
```

Preview the munged summary statistics.

``` r
munged_sumstats_preview <- utils::read.delim(
  munged_sumstats_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
head(munged_sumstats_preview)
#>    CHR  POS        ID REF ALT         BETA         SE        AF GDS_ID
#> 1 chr1 1000 rs0000001   A   G  0.065739730 0.01338672 0.2983750      1
#> 2 chr1 1010 rs0000002   A   G  0.051836263 0.01287591 0.3506250      2
#> 3 chr1 1020 rs0000003   A   G -0.004204495 0.01285285 0.3522500      3
#> 4 chr1 1030 rs0000004   A   G  0.038744761 0.01246723 0.3762500      4
#> 5 chr1 1040 rs0000005   A   G  0.059339484 0.01267735 0.3783125      5
```

Use `?` to inspect full documentation.

``` r
?tlstractor::munge_sumstats
```

## Step 3: Run TLS-Tractor

With the individual-level GDS file and phenotype/covariates from the
internal study, and the munged external GWAS summary statistics, we are
ready to run TLS-Tractor. TLS-Tractor estimates ancestry-specific
genetic effects with optional local-ancestry adjustment. The key
advantage is that although external summary statistics come from
standard GWAS models, TLS-Tractor strategically leverages them via
transfer learning to obtain unbiased and efficient estimates for
ancestry-specific parameters.

Model assumptions for TLS-Tractor: - Unrelated individuals in the
internal cohort - No sample overlap between internal and external
cohorts - **Transportability assumption**: comparable admixture profiles
across cohorts (e.g., similar ancestral populations and proportions).
This assumption is testable in practice, and a concrete diagnostic
workflow is provided in the [Assumption checking](#assumption-checking)
subsection below. Looking ahead, to relax this assumption, the next
version of TLS-Tractor is planned to support settings where internal and
external cohorts share the same ancestral populations but have different
ancestry proportions.

### Input parameters

- **`gds_path`**: Path to the GDS file generated in Step 1, containing
  local ancestry and risk allele dosage information for the internal
  cohort.
- **`sumstats_path`**: Path to the munged external GWAS summary
  statistics file generated in Step 2, which must be harmonized and
  aligned with the GDS file.
- **`method`**: Choose `"linear"` for continuous phenotypes or
  `"logistic"` for binary phenotypes.
- **`cond_local`**: Set to `TRUE` to condition on local-ancestry terms,
  which is recommended for most admixed population studies. This keeps
  interpretation of ancestry-specific SNP effects consistent with
  interpretation of the marginal effect in single-ancestry GWAS. Set to
  `FALSE` to estimate ancestry-specific SNP effects without
  local-ancestry adjustment, which can increase power but may be more
  susceptible to confounding by local ancestry.
- **`pheno_path`, `pheno_id_col`, `pheno_col`**: Specify the phenotype
  file, along with the column containing sample IDs and the column
  containing the phenotype values. Sample IDs must match between the GDS
  file and phenotype file, and only the intersecting samples will be
  used.
- **`covar_path`, `covar_id_col`, `covar_cols`**: Optionally provide a
  covariate file, the column containing sample IDs, and a vector of
  columns to include as covariates (e.g., PCs, age, sex). If covariates
  are supplied, sample IDs must match between the GDS file and covariate
  files, and only the intersecting samples will be used. If no
  covariates are needed, set `covar_path = NULL`.
- **`output_prefix`**: Output path prefix for result files. TLS-Tractor
  writes the main result to `<output_prefix>.txt.gz` and may write
  `<output_prefix>.excluded_samples.txt` when samples are dropped from
  the GDS file during intersection/filtering.
- **`scratch_dir`**: Optional directory for temporary per-task files
  before merge. If `NULL` (default), a run-specific temporary directory
  is created under `dirname(output_prefix)` with name
  `<basename(output_prefix)>_<pid>_<YYYYmmdd_HHMMSS>_tmp`. Temporary
  files are deleted upon successful completion. The directory is removed
  only if it is created by the function.
- **`snp_start`**: 1-based starting SNP index in the GDS file (default
  `1L`). Useful when analyzing a subset of SNPs.
- **`snp_count`**: Number of SNPs to process from `snp_start`. If `NULL`
  (default), analysis continues to the end of the GDS file.
- **`n_cores`**: Number of CPU cores to use for parallelization (default
  `1L`). Increase to speed up runtime but require more memory.
- **`chunk_size`**: Number of variants processed per chunk per core
  (default `1024L`). Larger values can improve speed but increase memory
  usage.
- **`local_ancestry_mac_threshold`**: Ancestry-specific minor-allele
  count threshold (default `20L`). SNPs with any ancestry-specific MAC
  below this threshold are skipped and returned with `NA` for
  inferential columns.
- **`use_fast_version`**: Defaults to `TRUE` (recommended, especially
  for large datasets). In fast mode, TLS-Tractor assumes that the
  estimated non-genetic covariate effects from the null model
  (`phenotype ~ covariates`) are close to those from the full standard
  GWAS model (`phenotype ~ SNP dosage + covariates`) for any single SNP,
  so that estimated covariate effects from the null model can be reused
  across variants to reduce computation. The fast mode is substantially
  faster with little loss of accuracy. Set to `FALSE` to force full
  per-variant fitting when you want maximum robustness. If no covariates
  are provided, fast mode is automatically disabled.

Run TLS-Tractor using internal individual-level data and the munged
external summary statistics:

``` r
# Create output directory for TLS-Tractor results
output_dir_results <- file.path(tempdir(), "tlstractor_tutorial",
                                "linear", "results")
dir.create(output_dir_results, recursive = TRUE, showWarnings = FALSE)

# Run TLS-Tractor
tlstractor(
  gds_path = gds_path,
  sumstats_path = munged_sumstats_path,
  method = "linear",
  cond_local = TRUE,
  pheno_path = file.path(data_dir, "simulate_linear.pheno.txt"),
  pheno_id_col = "id",
  pheno_col = "y",
  covar_path = file.path(data_dir, "simulate_linear.covar.txt"),
  covar_id_col = "id",
  covar_cols = c("cov1", "cov2", "cov3"),
  output_prefix = file.path(output_dir_results, "tlstractor_linear")
)
#> Reading GDS metadata...
#> Total samples in GDS: 100
#> Total SNPs in GDS: 5
#> Loading phenotype data...
#> Loading covariate data...
#> Filtering samples...
#> Samples after filtering: 100
#> Preprocessing data...
#> Loading summary statistics...
#> Adjusting chunk_size from 1024 to 5 to ensure enough tasks for 1 cores.
#> Parallel setup: using 1 of 4 available cores. Planned 1 tasks to process SNPs [1, 5] (5 total). Each task handles up to 5 SNPs (last task may be smaller). Genotypes are read in chunks of 5 SNPs within each task.
#> Scratch directory: /tmp/RtmphmsPgt/tlstractor_tutorial/linear/results/tlstractor_linear_14286_20260412_180149_tmp
#> Output filepath: /tmp/RtmphmsPgt/tlstractor_tutorial/linear/results/tlstractor_linear.txt.gz
#> Initializing parallel cluster...
#> Running TLS-Tractor analysis in parallel...
#> Merging results...
#> Cleaning up temporary files...
#> TLS-Tractor analysis complete!
#> Results written to: /tmp/RtmphmsPgt/tlstractor_tutorial/linear/results/tlstractor_linear.txt.gz
```

### Output interpretation

TLS-Tractor outputs a summary statistics file `<output_prefix>.txt.gz`
with the following columns:

1.  **Variant metadata**:
    - `CHROM`: Chromosome
    - `POS`: Base-pair position
    - `ID`: Variant ID
    - `REF`: Reference allele
    - `ALT`: Alternate allele
    - `main_N`: Number of samples analyzed in the internal/main study
2.  **Frequency and ancestry summaries**:
    - `AF`: Overall allele frequency (from internal study)
    - `AF_anc*`: Ancestry-specific allele frequency for each ancestry
      (e.g., `AF_anc0`, `AF_anc1` for two-way admixture)
    - `LAprop_anc*`: Local ancestry proportion for each ancestry
3.  **Analysis metadata**:
    - `has_sumstats`: Logical indicator of whether external summary
      statistics are available for this variant. Variants without
      summary statistics are analyzed using internal data alone.
    - `fallback_used`: Logical indicator of whether QR decomposition
      fallback was used during estimation; when `TRUE`, results may have
      reduced numerical stability and should be interpreted with caution
4.  **Association results** (main results):
    - `joint_pval`: Joint p-value for testing all ancestry-specific SNP
      dosage effects
    - `beta_anc*`: Effect size (log odds for logistic, linear
      coefficient for linear) for each ancestry
    - `se_anc*`: Standard error of the effect size for each ancestry
    - `pval_anc*`: P-value for each ancestry-specific effect
5.  **Local ancestry effects** (when `cond_local = TRUE`):
    - `LAeff_anc*`: Effect size (log odds for logistic, linear
      coefficient for linear) of the local ancestry term for each
      ancestry
    - `LAse_anc*`: Standard error of the local ancestry effect for each
      ancestry
    - `LApval_anc*`: P-value for the local ancestry effect for each
      ancestry

``` r
# Path to the output file
tlstractor_output_path <- file.path(output_dir_results,
                                    "tlstractor_linear.txt.gz")

# Read the results
results <- utils::read.delim(
  tlstractor_output_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Preview the first few rows
head(results)
#>   CHROM  POS        ID REF ALT main_N    AF joint_pval has_sumstats   AF_anc0
#> 1  chr1 1000 rs0000001   A   G    100 0.300 0.06567632         TRUE 0.3076923
#> 2  chr1 1010 rs0000002   A   G    100 0.315 0.30876716         TRUE 0.2056075
#> 3  chr1 1020 rs0000003   A   G    100 0.370 0.97100185         TRUE 0.4019608
#> 4  chr1 1030 rs0000004   A   G    100 0.380 0.45949422         TRUE 0.5154639
#> 5  chr1 1040 rs0000005   A   G    100 0.380 0.01970145         TRUE 0.2747253
#>     AF_anc1 LAprop_anc0 LAprop_anc1   beta_anc0    beta_anc1    se_anc0
#> 1 0.2916667       0.520       0.480  0.15165327 -0.074057144 0.09800820
#> 2 0.4408602       0.535       0.465  0.11108844 -0.008669883 0.13797109
#> 3 0.3367347       0.510       0.490 -0.01886665  0.023525459 0.11604381
#> 4 0.2524272       0.485       0.515  0.06819441 -0.002498296 0.09589733
#> 5 0.4678899       0.455       0.545  0.08477358  0.040311881 0.13228691
#>      se_anc1 pval_anc0 pval_anc1  LAeff_anc0 LAse_anc0 LApval_anc0
#> 1 0.12161558 0.1217782 0.5425612 -0.07223941 0.1200606   0.5473798
#> 2 0.09677068 0.4207289 0.9286114 -0.08770070 0.1318639   0.5059961
#> 3 0.10977918 0.8708475 0.8303147 -0.04317887 0.1396011   0.7570919
#> 4 0.16867080 0.4770106 0.9881824 -0.05379375 0.1253552   0.6678283
#> 5 0.09816028 0.5216323 0.6813115 -0.03541323 0.1354502   0.7937468
#>   fallback_used
#> 1         FALSE
#> 2         FALSE
#> 3         FALSE
#> 4         FALSE
#> 5         FALSE

# Check output dimensions
dim(results)
#> [1]  5 23

# List all column names
colnames(results)
#>  [1] "CHROM"         "POS"           "ID"            "REF"          
#>  [5] "ALT"           "main_N"        "AF"            "joint_pval"   
#>  [9] "has_sumstats"  "AF_anc0"       "AF_anc1"       "LAprop_anc0"  
#> [13] "LAprop_anc1"   "beta_anc0"     "beta_anc1"     "se_anc0"      
#> [17] "se_anc1"       "pval_anc0"     "pval_anc1"     "LAeff_anc0"   
#> [21] "LAse_anc0"     "LApval_anc0"   "fallback_used"
```

### Assumption checking

A practical way to assess TLS-Tractor’s transportability assumption is
to test for effect heterogeneity between the internal and external
cohorts. Specifically, - from the internal cohort, obtain effect
estimates $\beta_{1}$ and standard errors $SE_{1}$ for each SNP from a
standard GWAS model (fitted using tools such as PLINK2) - from the
external GWAS summary statistics, obtain effect estimates $\beta_{2}$
and standard errors $SE_{2}$ for the same SNPs (these should already be
available from the munged summary statistics file)

For each SNP, test whether the effect sizes are consistent across
cohorts using:
$$z = \frac{\beta_{1} - \beta_{2}}{\sqrt{SE_{1}^{2} + SE_{2}^{2}}}$$

Under the null hypothesis $H_{0}:\beta_{1} = \beta_{2}$, the $z$
statistic approximately follows a standard normal distribution, assuming
the two estimates are independent. A two-sided p-value can be computed
accordingly. This is a Wald test for heterogeneity between two
independent estimates, and is equivalent to Cochran’s Q test in the
special case of two studies (i.e., $Q = z^{2}$).

After computing p-values genome-wide, a QQ plot can be used to assess
calibration of the heterogeneity test. If the transportability
assumption is satisfied, most SNPs should exhibit no evidence of
heterogeneity, and the observed p-values should follow the expected null
distribution (i.e., align with the diagonal). Systematic deviations may
indicate violations of the assumption.

``` r
# Example workflow:
# 1) Run standard GWAS in the internal cohort and obtain per-SNP effect size
#    (beta_1) and standard error (se_1).
# 2) Merge with effect size (beta_2) and standard error (se_2) from external
#    GWAS summary statistics by SNP ID (or CHR/POS/alleles).

assump_dt <- merge(
  internal_gwas[, c("ID", "beta_1", "se_1")],
  external_sumstats[, c("ID", "beta_2", "se_2")],
  by = "ID"
)

# Wald z-test for equality of two independent estimates
assump_dt$z <- (assump_dt$beta_1 - assump_dt$beta_2) /
  sqrt(assump_dt$se_1^2 + assump_dt$se_2^2)
assump_dt$pval <- 2 * stats::pnorm(-abs(assump_dt$z))

# QQ plot of p-values
p_sorted <- sort(assump_dt$pval)
n <- length(p_sorted)
exp_logp <- -log10(stats::ppoints(n))
obs_logp <- -log10(p_sorted)

# 95% CI under the null (order-statistic based)
idx <- seq_len(n)
ci_low <- stats::qbeta(0.025, idx, n - idx + 1)
ci_high <- stats::qbeta(0.975, idx, n - idx + 1)
ci_lower_logp <- -log10(pmax(ci_high, .Machine$double.xmin))
ci_upper_logp <- -log10(pmax(ci_low, .Machine$double.xmin))

plot(
  exp_logp,
  obs_logp,
  type = "n",
  xlab = expression(Expected~ -log[10](p)),
  ylab = expression(Observed~ -log[10](p)),
  main = "QQ Plot of Heterogeneity Test P-values"
)

polygon(
  c(exp_logp, rev(exp_logp)),
  c(ci_lower_logp, rev(ci_upper_logp)),
  col = grDevices::adjustcolor("grey70", alpha.f = 0.35),
  border = NA
)

points(exp_logp, obs_logp, pch = 20, cex = 0.6)
abline(0, 1, col = "red", lwd = 2)
```

### Read function documentation

Use `?` to inspect full documentation.

``` r
?tlstractor::tlstractor
```

## Full binary pipeline code

The code chunk below illustrates the full workflow for simulated binary
data, from tract extraction to the final TLS-Tractor output.

``` r
# Get path to tutorial data included in the package
data_dir <- system.file("extdata", package = "tlstractor")

# Store all temporary outputs under tempdir()/tlstractor_tutorial/binary/
work_dir <- file.path(tempdir(), "tlstractor_tutorial", "binary")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Extract tracts to GDS (choose one extractor; this example uses FLARE input)
extract_tracts_flare(
  vcf_path = file.path(data_dir, "simulate_binary.flare.vcf"),
  num_ancs = 2L,
  output_dir = work_dir,
  output_formats = "gds"
)
gds_path <- file.path(work_dir, "simulate_binary.flare.gds")

# If using RFMix2/Gnomix inputs instead, use:
# extract_tracts(
#   vcf_path = file.path(data_dir, "simulate_binary.phased.vcf"),
#   msp_path = file.path(data_dir, "simulate_binary.msp.tsv"),
#   num_ancs = 2L,
#   output_dir = work_dir,
#   output_formats = "gds"
# )
# gds_path <- file.path(work_dir, "simulate_binary.phased.gds")

# 2) Munge external GWAS summary statistics
munged_sumstats_path <- file.path(work_dir, "external_gwas_sumstats_munged.txt")
munge_sumstats(
  gds_path = gds_path,
  sumstats_path = file.path(data_dir,
                            "simulate_binary.external_gwas_sumstats.txt"),
  match_by = "CHR-POS",
  remove_ambiguous = TRUE,
  output_path = munged_sumstats_path
)

# 3) Run TLS-Tractor
tlstractor(
  gds_path = gds_path,
  sumstats_path = munged_sumstats_path,
  method = "logistic",
  cond_local = TRUE,
  pheno_path = file.path(data_dir, "simulate_binary.pheno.txt"),
  pheno_id_col = "id",
  pheno_col = "y",
  covar_path = file.path(data_dir, "simulate_binary.covar.txt"),
  covar_id_col = "id",
  covar_cols = c("cov1", "cov2", "cov3"),
  output_prefix = file.path(work_dir, "tlstractor_binary")
)
```

## Cleanup

Remove all temporary tutorial outputs created under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html).

``` r
unlink(file.path(tempdir(), "tlstractor_tutorial"),
       recursive = TRUE, force = TRUE)
```
