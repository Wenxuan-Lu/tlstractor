# tlstractor

The R package **tlstractor** implements **TLS-Tractor** (**T**ransfer **L**earning of **S**ummary Statistics to **Tractor**), a transfer learning framework to improve the power of local ancestry-aware genome-wide association studies (GWAS) in admixed populations. A common practical challenge in local ancestry-aware GWAS is limited sample size. Standard GWAS can scale by meta-analyzing summary statistics across cohorts, but local ancestry-aware analyses typically require individual-level data and cannot directly use existing summary statistics from other cohorts. TLS-Tractor addresses this gap by combining internal individual-level data with external GWAS summary statistics to estimate ancestry-specific genetic effects, with optional adjustment for local ancestry. Core computational steps are implemented in C++ to keep analyses efficient and scalable for large datasets. The method assumes unrelated individuals in the internal cohort, no sample overlap between internal and external cohorts, and comparable admixture profiles across cohorts.

Last updated: March 29, 2026

Current version: 0.0.0.9000 (in development)

## Contents

-   [Installation](#installation)
-   [Quick start](#quick-start)
-   [Installation troubleshooting](#installation-troubleshooting)
-   [Bug reports](#bug-reports)
-   [License](#license)

<!-- After installation notes,
- [Citation](#citation)
-->

## Installation

### Optional pre-setup (Conda-based R only)

If you're not using conda, skip this subsection and continue to the install options below.

Create and use a dedicated conda environment for tlstractor before installation:

``` bash
conda create -n r-tlstractor -c conda-forge r-base r-essentials -y
conda activate r-tlstractor
conda install -c conda-forge compilers make pkg-config zlib
```

### Install tlstractor

#### Option 1) Install the development version (GitHub)

Use this option for the latest features and fixes. Installation may take a few minutes.

``` r
install.packages("pak")
pak::pkg_install("Wenxuan-Lu/tlstractor")
```

#### Option 2) Install a specific release (source tarball)

Use this option for a fixed, reproducible version.

``` r
install.packages("pak")
pak::pkg_install("https://github.com/Wenxuan-Lu/tlstractor/releases/download/vX.Y.Z/tlstractor_X.Y.Z.tar.gz")
```

If installation fails or you want to build vignettes locally, see [Installation troubleshooting](#installation-troubleshooting).

[Back to Contents](#contents)

## Quick start

TLS-Tractor extends [Tractor](https://github.com/Atkinson-Lab/Tractor) by incorporating external GWAS summary statistics into local ancestry-aware association analysis. For background on GWAS in admixed populations and local ancestry inference, see the [Tractor tutorial](https://atkinson-lab.github.io/Tractor-tutorial/).

This package follows a 3-step workflow: extract local-ancestry tracts, harmonize external GWAS summary statistics, then combine the two data sources and run local ancestry-aware association analysis.

Our implementation of the first step (`extract_tracts()` / `extract_tracts_flare()`) is workflow-compatible with Tractor. It is designed as a drop-in replacement and implemented in C++ for substantially faster tract extraction. For workflows that need tract files in different formats, we provide functions `convert_gds_to_txt()` and `convert_txt_to_gds()` for quick format conversion.

We provide a minimal runnable script below for illustration. For a full pipeline walkthrough, see the [TLS-Tractor tutorial](https://wenxuan-lu.github.io/tlstractor/).

**Required inputs:**

Individual-level data from the main study / internal cohort:

-   Phased genotypes with local ancestry inference results (choose one):

    -   If using RFMix2/Gnomix for LAI: `cohort.vcf(.gz)` + `cohort.msp.tsv(.gz)`

        -   `cohort.vcf(.gz)`: phased genotypes (`GT`, e.g. `0|1`)

        -   `cohort.msp.tsv(.gz)`: local ancestry labels

    -   If using FLARE for LAI: `cohort.vcf(.gz)`

        -   FLARE output VCF containing both phased genotypes and local ancestry labels

-   Phenotype file: `pheno.txt` with columns `id`, `pheno`

-   Optional covariates: `covariates.txt` with columns `id`, `cov1`, `cov2`, `cov3`

Summary-level data from an external cohort:

-   External GWAS summary statistics (from a standard GWAS model): `external_sumstats.txt` with columns `CHR`, `POS`, `ID`, `REF`, `ALT`, `BETA`, `SE`, and an optional column `AF`.

Filenames and column names are illustrative; map them to your own data schema.

``` r
library(tlstractor)

# 1) Extract tracts (choose one option)

# Option A: phased VCF + MSP (RFMix2/Gnomix)
extract_tracts(
  vcf_path = "cohort.vcf.gz",
  msp_path = "cohort.msp.tsv",
  num_ancs = 2L,
  output_formats = "gds"
)

# Option B: FLARE VCF
# extract_tracts_flare(
#   vcf_path = "cohort.vcf.gz",
#   num_ancs = 2L,
#   output_formats = "gds"
# )

# 2) QC + align external summary statistics
munge_sumstats(
  gds_path = "cohort.gds",
  sumstats_path = "external_sumstats.txt",
  match_by = "CHR-POS",
  output_path = "external_sumstats_munged.txt"
)

# 3) Run TLS-Tractor
tlstractor(
  gds_path = "cohort.gds",
  sumstats_path = "external_sumstats_munged.txt",
  method = "linear", # or "logistic"
  cond_local = TRUE, # whether to condition on local ancestry terms
  pheno_path = "pheno.txt",
  pheno_id_col = "id",
  pheno_col = "pheno",
  covar_path = "covariates.txt",
  covar_id_col = "id",
  covar_cols = c("cov1", "cov2", "cov3"),
  output_prefix = "tlstractor",
  n_cores = 4L
)
```

**Main outputs:**

1.  **Tracts GDS** (`cohort.gds`) — HDF5-based GDS file containing:

    -   `sample.id` — sample identifiers
    -   `snp.chromosome`, `snp.position`, `snp.id`, `snp.ref`, `snp.alt` — variant metadata
    -   `dosage/anc0..ancK`, `hapcount/anc0..ancK` — ancestry-specific dosage (0, 1, 2) and local ancestry haplotype count (0, 1, 2) matrices

2.  **Munged summary statistics** (`external_sumstats_munged.txt`) — Tab-delimited file with columns:

    -   `CHR`, `POS`, `ID`, `REF`, `ALT` — variant information
    -   `BETA`, `SE` — effect size and standard error
    -   `AF` (optional) — allele frequency if available
    -   `GDS_ID` — internal index linking to GDS file variant

3.  **TLS-Tractor results** (`tlstractor.txt.gz`) — Gzip-compressed tab-delimited file with output columns organized as follows:

    **Variant metadata:**

    -   `CHROM`, `POS`, `ID`, `REF`, `ALT` — Variant information
    -   `main_N` — Number of samples analyzed from the internal/main study

    **Allele frequency and ancestry composition:**

    -   `AF` — Overall allele frequency (calculated from internal cohort)
    -   `AF_anc*` — Ancestry-specific allele frequency for each ancestry (e.g., `AF_anc0`, `AF_anc1` for two-way admixture)
    -   `LAprop_anc*` — Proportion of local ancestry for each ancestry (reflects the ancestry makeup at this locus)

    **Analysis metadata:**

    -   `has_sumstats` — Logical indicator (`TRUE`/`FALSE`) of whether external summary statistics are available for this variant. When `FALSE`, the variant was analyzed using internal data alone
    -   `fallback_used` — Logical indicator (`TRUE`/`FALSE`) of whether QR decomposition fallback was used during statistical estimation. When `TRUE`, results may have reduced numerical stability and should be interpreted cautiously

    **Association results (main output):**

    -   `joint_pval` — Joint p-value for testing all ancestry-specific SNP dosage effects
    -   `beta_anc*` — Effect size (linear coefficient for linear phenotype, log odds for logistic phenotype) for each ancestry
    -   `se_anc*` — Standard error of the ancestry-specific effect size
    -   `pval_anc*` — P-value for each ancestry-specific effect

    **Local ancestry effects (only present when `cond_local=TRUE`):**

    -   `LAeff_anc*` — Effect size of the local ancestry term for each ancestry
    -   `LAse_anc*` — Standard error of the local ancestry effect
    -   `LApval_anc*` — P-value for the local ancestry effect

**Next:**

-   Full pipeline walkthrough: [TLS-Tractor tutorial](https://wenxuan-lu.github.io/tlstractor/)

-   Tutorial on statistical phasing and local ancestry inference: [Tractor tutorial](https://atkinson-lab.github.io/Tractor-tutorial/)

-   Function help in R: `?extract_tracts`, `?extract_tracts_flare`, `?convert_gds_to_txt`, `?convert_txt_to_gds`, `?munge_sumstats`, `?tlstractor`

[Back to Contents](#contents)

## Installation troubleshooting

Read this section only if installation fails.

### 1) Check dependencies first

If installation fails, the most common issue is that `gdsfmt` may not install correctly during package installation.

Install it explicitly from Bioconductor:

```r
install.packages("BiocManager")
BiocManager::install("gdsfmt")
```

Then retry installing `tlstractor`.

### 2) Check system requirements

This package includes C++ code and links against zlib. Most users will not need to do anything.

#### 2.1) Quick check

Run `Sys.which(c("make", "g++", "pkg-config"))` in R. All entries should be non-empty.

To check for zlib, run `system("pkg-config --exists zlib") == 0` (skip this if `pkg-config` is not installed). This should return `TRUE` if zlib is installed.

#### 2.2) Install required tools

-   **Windows**: install Rtools (matching your R version).
-   **macOS**: run `xcode-select --install`, then `brew install zlib pkg-config`
-   **Debian/Ubuntu**: `sudo apt-get install -y build-essential zlib1g-dev pkg-config`
-   **Fedora/RHEL/CentOS**: `sudo dnf install -y gcc gcc-c++ make zlib-devel pkgconf-pkg-config`
-   **Alpine**: `apk add --no-cache build-base zlib-dev pkgconf`

#### 2.3) Retry installation

Retry with `pak::pkg_install("Wenxuan-Lu/tlstractor")` or install a release tarball with\
`pak::pkg_install("https://github.com/Wenxuan-Lu/tlstractor/releases/download/vX.Y.Z/tlstractor_X.Y.Z.tar.gz")`.

#### 2.4) Advanced: manual zlib configuration

Only needed if installation still fails after installing zlib.

Set zlib paths explicitly (replace the paths with your system’s locations):

`Sys.setenv(ZLIB_CPPFLAGS = "-I/path/include", ZLIB_LIBS = "-L/path/lib -lz")`

Then retry installation.

[Back to Contents](#contents)

<!-- ### 3) Optional: build vignettes

We recommend installing without building vignettes; instead, read tutorials on the package website.

To install vignettes locally, ensure Pandoc is installed by checking `rmarkdown::pandoc_available()`.
Also ensure `gdsfmt` package is installed; if not, run `BiocManager::install("gdsfmt")`.

Then, install with vignettes building:

-   `remotes::install_github("Wenxuan-Lu/tlstractor", build_vignettes = TRUE)`

-   or `remotes::install_url("https://github.com/Wenxuan-Lu/tlstractor/releases/download/vX.Y.Z/tlstractor_X.Y.Z.tar.gz", build_vignettes = TRUE)`

View vignettes with `browseVignettes("tlstractor")`. -->

<!--
## Citation

If you use this package, please cite our manuscript.

TODO: add usethis::use_citation() later to include inst/CITATION

[Back to Contents](#contents) -->

## Bug reports

If you encounter a bug or have a feature request, please open an issue at <https://github.com/Wenxuan-Lu/tlstractor/issues>.

For other inquiries or feedback, please contact Wenxuan Lu at [wlu15@jhu.edu](mailto:[wlu15@jhu.edu]).

[Back to Contents](#contents)

## License

This package is licensed under the MIT License.

[Back to Contents](#contents)
