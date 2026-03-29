# Run TLS-Tractor

## Usage

``` r
tlstractor(
  gds_path,
  sumstats_path,
  method,
  cond_local,
  pheno_path,
  pheno_id_col,
  pheno_col,
  covar_path = NULL,
  covar_id_col = NULL,
  covar_cols = NULL,
  output_prefix,
  scratch_dir = NULL,
  snp_start = 1L,
  snp_count = NULL,
  n_cores = 1L,
  chunk_size = 1024L,
  local_ancestry_mac_threshold = 20L,
  use_fast_version = TRUE
)
```

## Arguments

- gds_path:

  Character; path to the input GDS file.

- sumstats_path:

  Character; path to the munged summary statistics file.

- method:

  Character; one of `"linear"` or `"logistic"`.

- cond_local:

  Logical; whether to condition on local-ancestry dosage terms.

- pheno_path:

  Character; path to the phenotype file.

- pheno_id_col:

  Character; sample-ID column name in the phenotype file.

- pheno_col:

  Character; phenotype column name in the phenotype file.

- covar_path:

  Character or `NULL`; path to an optional covariate file. Default is
  `NULL`.

- covar_id_col:

  Character or `NULL`; sample-ID column name in the covariate file
  (required when `covar_path` is provided). Default is `NULL`.

- covar_cols:

  Character vector or `NULL`; covariate column names in the covariate
  file (required when `covar_path` is provided). Default is `NULL`.

- output_prefix:

  Character; output path prefix.

- scratch_dir:

  Character or `NULL`; directory used to temporarily store per-task
  result files before merge. If `NULL`, a run-specific directory is
  created under `dirname(output_prefix)` with name
  `<basename(output_prefix)>_<pid>_<YYYYmmdd_HHMMSS>_tmp`. Default is
  `NULL`.

- snp_start:

  Integer; 1-based starting SNP index in the GDS file. Default is `1L`.

- snp_count:

  Integer or `NULL`; number of SNPs to process. If `NULL`, processing
  starts at `snp_start` and continues to the end of the GDS file.
  Default is `NULL`.

- n_cores:

  Integer; number of CPU cores to use. Default is `1L`.

- chunk_size:

  Integer; number of variants processed per chunk within each CPU core.
  Default is 1024; larger values may improve speed but increase memory
  usage.

- local_ancestry_mac_threshold:

  Integer; ancestry-specific minor-allele count threshold. SNPs with any
  ancestry-specific MAC below this threshold are skipped. Skipped SNPs
  are returned in the output with `NA` values for inferential columns.
  Default is `20L`.

- use_fast_version:

  Logical; whether to use the fast TLS-Tractor mode. The fast mode
  assumes that the estimated non-genetic covariate effects from the null
  model (phenotype ~ covariates) are close to those from the full
  standard GWAS model (phenotype ~ SNP dosage + covariates) for any
  single SNP, so that estimated covariate effects from the null model
  can be reused across variants to reduce computation. The fast mode can
  provide a substantial speedup with minimal loss of accuracy. Default
  is `TRUE` (recommended for large datasets). Setting `FALSE` fources
  the full per-variant fitting path, which is generally more robust but
  substantially slower when covariates are included. If no covariates
  are provided, fast mode is automatically disabled.

## Value

Invisibly returns `NULL`. Writes gzipped GWAS results to
`<output_prefix>.txt.gz`. The output includes:

- Variant metadata: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `main_N` (sample
  size in the main/internal study)

- Frequency and ancestry summaries: `AF` (overall allele frequency),
  `AF_anc*` (ancestry-specific allele frequency), `LAprop_anc*`
  (local-ancestry-specific haplotype proportion)

- Analysis metadata: `has_sumstats` (indicates whether external summary
  statistics were available), `fallback_used` (indicates whether QR
  decomposition fallback was used; when TRUE, results may have reduced
  numerical stability and should be interpreted with caution)

- Association results: `joint_pval` (p-value for joint testing all
  ancestry-specific SNP dosage effects), `beta_anc*`, `se_anc*`,
  `pval_anc*` (effect size, standard error, and p-value for each
  ancestry-specific SNP dosage effect)

- When `cond_local = TRUE`: `LAeff_anc*`, `LAse_anc*`, `LApval_anc*`
  (effect size, standard error, and p-value for each local ancestry
  term) }

Performs local ancestry-aware GWAS by integrating individual-level data
with external GWAS summary statistics via transfer learning. If present,
output the file `<output_prefix>.excluded_samples.txt` containing sample
IDs present in the GDS file but excluded from analysis after sample
intersection/filtering.During execution, temporary per-task result files
are written to `scratch_dir/<run_tag>_task_<id>.txt.gz`, where `run_tag`
is a run-specific identifier with format
`tlstractor_<pid>_<YYYYmmdd_HHMMSS>`. These temporary files are removed
on successful cleanup. The directory `scratch_dir` is removed only if it
was created by the function.SNPs present in the GDS file but absent from
the summary statistics file are still analyzed using the Tractor model
with individual-level data only.
