# Simulate multi-SNP two-way admixed data with internal/external split.
# The output includes:
# - RFMix2-style phased VCF (GT only) + MSP for internal samples
# - FLARE-style VCF (GT:AN1:AN2) for internal samples
# - covariate file for internal samples
# - phenotype file for internal samples
# - external GWAS summary statistics

normalize_to_sd1 <- function(x) {
  x <- as.numeric(x)
  x <- x - mean(x)
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx == 0) {
    return(rep(0, length(x)))
  }
  x / sx
}

ensure_length_n <- function(x, n, name) {
  if (length(x) == n) {
    return(x)
  }
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  stop(sprintf("%s must have length 1 or n_snps (%d).", name, n), call. = FALSE)
}

as_logical_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

# Helper 1: simulate local ancestry and genotype for one SNP.
simulate_one_snp_genotype_local_ancestry <- function(
  n_individuals,
  prop_anc1,
  maf_anc0,
  maf_anc1
) {
  n_individuals <- as.integer(n_individuals)
  stopifnot(n_individuals > 1L)
  stopifnot(is.finite(prop_anc1), prop_anc1 > 0, prop_anc1 < 1)
  stopifnot(is.finite(maf_anc0), maf_anc0 > 0, maf_anc0 < 1)
  stopifnot(is.finite(maf_anc1), maf_anc1 > 0, maf_anc1 < 1)

  anc1_h1 <- stats::rbinom(n_individuals, size = 1L, prob = prop_anc1)
  anc1_h2 <- stats::rbinom(n_individuals, size = 1L, prob = prop_anc1)
  local_anc1_count <- anc1_h1 + anc1_h2

  p_h1 <- ifelse(anc1_h1 == 1L, maf_anc1, maf_anc0)
  p_h2 <- ifelse(anc1_h2 == 1L, maf_anc1, maf_anc0)

  alt_h1 <- stats::rbinom(n_individuals, size = 1L, prob = p_h1)
  alt_h2 <- stats::rbinom(n_individuals, size = 1L, prob = p_h2)

  dosage_anc1 <- alt_h1 * anc1_h1 + alt_h2 * anc1_h2
  dosage_anc0 <- alt_h1 * (1L - anc1_h1) + alt_h2 * (1L - anc1_h2)
  dosage_total <- alt_h1 + alt_h2

  list(
    anc1_h1 = anc1_h1,
    anc1_h2 = anc1_h2,
    local_anc1_count = local_anc1_count,
    alt_h1 = alt_h1,
    alt_h2 = alt_h2,
    dosage_anc0 = dosage_anc0,
    dosage_anc1 = dosage_anc1,
    dosage_total = dosage_total
  )
}

# Helper 2: generate covariates and phenotype from multi-SNP genotype/local ancestry.
simulate_covariates_and_phenotype <- function(
  genotype_list,
  n_covariates = 3L,
  beta_covar = c(0.9, -0.7, 0.5),
  beta_dosage_anc0 = 0.10,
  beta_dosage_anc1 = 0.15,
  beta_local_anc1 = 0.08,
  covariate_r2 = 0.35,
  genetic_r2 = 0.03,
  local_ancestry_r2 = 0.02,
  phenotype_type = c("linear", "binary")
) {
  phenotype_type <- match.arg(phenotype_type)
  n_covariates <- as.integer(n_covariates)
  n_snps <- length(genotype_list)
  stopifnot(n_snps >= 1L, n_covariates >= 1L)

  beta_dosage_anc0 <- ensure_length_n(beta_dosage_anc0, n_snps, "beta_dosage_anc0")
  beta_dosage_anc1 <- ensure_length_n(beta_dosage_anc1, n_snps, "beta_dosage_anc1")
  beta_local_anc1 <- ensure_length_n(beta_local_anc1, n_snps, "beta_local_anc1")

  if (length(beta_covar) != n_covariates) {
    stop("length(beta_covar) must equal n_covariates.", call. = FALSE)
  }

  noise_r2 <- 1 - (covariate_r2 + genetic_r2 + local_ancestry_r2)
  if (noise_r2 <= 0) {
    stop("covariate_r2 + genetic_r2 + local_ancestry_r2 must be < 1.", call. = FALSE)
  }

  n <- length(genotype_list[[1]]$dosage_total)
  covariates <- matrix(stats::rnorm(n * n_covariates), nrow = n, ncol = n_covariates)
  colnames(covariates) <- paste0("cov", seq_len(n_covariates))

  cov_raw <- as.numeric(covariates %*% beta_covar)
  gen_raw <- rep(0, n)
  la_raw <- rep(0, n)

  for (j in seq_len(n_snps)) {
    g <- genotype_list[[j]]
    gen_raw <- gen_raw + beta_dosage_anc0[j] * g$dosage_anc0 + beta_dosage_anc1[j] * g$dosage_anc1
    la_raw <- la_raw + beta_local_anc1[j] * g$local_anc1_count
  }

  cov_component <- sqrt(covariate_r2) * normalize_to_sd1(cov_raw)
  gen_component <- sqrt(genetic_r2) * normalize_to_sd1(gen_raw)
  la_component <- sqrt(local_ancestry_r2) * normalize_to_sd1(la_raw)
  noise_component <- stats::rnorm(n, mean = 0, sd = sqrt(noise_r2))

  latent <- cov_component + gen_component + la_component + noise_component

  if (phenotype_type == "linear") {
    y <- latent
  } else {
    pr <- stats::plogis(latent)
    y <- stats::rbinom(n, size = 1L, prob = pr)
  }

  list(
    covariates = covariates,
    phenotype = y,
    phenotype_type = phenotype_type
  )
}

# Helper 3: split into internal/external and compute external GWAS summary stats for one SNP.
make_internal_external_one_snp <- function(
  genotype,
  covariates,
  phenotype,
  sample_ids,
  internal_idx,
  chrom,
  pos,
  snp_id,
  ref,
  alt
) {
  n <- length(sample_ids)
  stopifnot(length(internal_idx) >= 2L, length(internal_idx) < n)
  external_idx <- setdiff(seq_len(n), internal_idx)

  ext_df <- data.frame(
    y = phenotype[external_idx],
    dosage = genotype$dosage_total[external_idx],
    stringsAsFactors = FALSE
  )

  cov_df_ext <- as.data.frame(covariates[external_idx, , drop = FALSE], stringsAsFactors = FALSE)
  names(cov_df_ext) <- colnames(covariates)
  ext_df <- cbind(ext_df, cov_df_ext)

  all_binary <- all(ext_df$y %in% c(0, 1))
  rhs <- paste(c("dosage", names(cov_df_ext)), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  if (all_binary) {
    fit <- stats::glm(form, data = ext_df, family = stats::binomial())
  } else {
    fit <- stats::lm(form, data = ext_df)
  }

  fit_sum <- summary(fit)$coefficients
  if (!"dosage" %in% rownames(fit_sum)) {
    stop("Failed to estimate dosage coefficient in external model.", call. = FALSE)
  }

  af_external <- mean(genotype$dosage_total[external_idx]) / 2

  gwas_summary <- data.frame(
    CHR = chrom,
    POS = as.integer(pos),
    ID = snp_id,
    REF = ref,
    ALT = alt,
    BETA = unname(fit_sum["dosage", "Estimate"]),
    SE = unname(fit_sum["dosage", "Std. Error"]),
    AF = af_external,
    N = length(external_idx),
    stringsAsFactors = FALSE
  )

  internal <- data.frame(
    id = sample_ids[internal_idx],
    CHR = chrom,
    POS = as.integer(pos),
    ID = snp_id,
    REF = ref,
    ALT = alt,
    anc1_h1 = genotype$anc1_h1[internal_idx],
    anc1_h2 = genotype$anc1_h2[internal_idx],
    local_anc1_count = genotype$local_anc1_count[internal_idx],
    alt_h1 = genotype$alt_h1[internal_idx],
    alt_h2 = genotype$alt_h2[internal_idx],
    dosage_anc0 = genotype$dosage_anc0[internal_idx],
    dosage_anc1 = genotype$dosage_anc1[internal_idx],
    dosage_total = genotype$dosage_total[internal_idx],
    stringsAsFactors = FALSE
  )

  list(internal = internal, gwas_summary = gwas_summary)
}

write_multi_snp_internal_outputs <- function(
  internal_by_snp,
  gwas_summary,
  pheno_df,
  covar_df,
  output_dir,
  prefix,
  write_internal_data = TRUE,
  pos_step = 10L
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  write_internal_data <- as_logical_flag(write_internal_data, "write_internal_data")

  internal_all <- do.call(rbind, internal_by_snp)
  internal_all <- internal_all[order(internal_all$POS, internal_all$id), , drop = FALSE]

  sample_ids <- unique(internal_all$id)
  pheno_df <- pheno_df[order(pheno_df$id), , drop = FALSE]
  covar_df <- covar_df[order(covar_df$id), , drop = FALSE]

  pheno_path <- file.path(output_dir, sprintf("%s.pheno.txt", prefix))
  covar_path <- file.path(output_dir, sprintf("%s.covar.txt", prefix))
  internal_path <- file.path(output_dir, sprintf("%s.internal_data.txt", prefix))
  gwas_path <- file.path(output_dir, sprintf("%s.external_gwas_sumstats.txt", prefix))
  flare_vcf_path <- file.path(output_dir, sprintf("%s.flare.vcf", prefix))
  phased_vcf_path <- file.path(output_dir, sprintf("%s.phased.vcf", prefix))
  msp_path <- file.path(output_dir, sprintf("%s.msp.tsv", prefix))

  utils::write.table(pheno_df, file = pheno_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(covar_df, file = covar_path, sep = "\t", quote = FALSE, row.names = FALSE)
  if (isTRUE(write_internal_data)) {
    utils::write.table(internal_all, file = internal_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  utils::write.table(gwas_summary, file = gwas_path, sep = "\t", quote = FALSE, row.names = FALSE)

  # FLARE VCF
  con_flare <- file(flare_vcf_path, open = "wt")
  on.exit(try(close(con_flare), silent = TRUE), add = TRUE)
  writeLines("##fileformat=VCFv4.2", con_flare)
  writeLines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", con_flare)
  writeLines("##FORMAT=<ID=AN1,Number=1,Type=Integer,Description=\"Ancestry of first haplotype\">", con_flare)
  writeLines("##FORMAT=<ID=AN2,Number=1,Type=Integer,Description=\"Ancestry of second haplotype\">", con_flare)
  header <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_ids)
  writeLines(paste(header, collapse = "\t"), con_flare)

  for (j in seq_len(nrow(gwas_summary))) {
    snp <- gwas_summary[j, , drop = FALSE]
    d <- internal_all[internal_all$ID == snp$ID, , drop = FALSE]
    d <- d[match(sample_ids, d$id), , drop = FALSE]
    flare_fields <- paste0(d$alt_h1, "|", d$alt_h2, ":", d$anc1_h1, ":", d$anc1_h2)
    line <- c(
      as.character(snp$CHR), as.character(snp$POS), as.character(snp$ID),
      as.character(snp$REF), as.character(snp$ALT), ".", ".", ".", "GT:AN1:AN2", flare_fields
    )
    writeLines(paste(line, collapse = "\t"), con_flare)
  }

  # RFMix2-style phased VCF
  con_phased <- file(phased_vcf_path, open = "wt")
  on.exit(try(close(con_phased), silent = TRUE), add = TRUE)
  writeLines("##fileformat=VCFv4.2", con_phased)
  writeLines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", con_phased)
  writeLines(paste(header, collapse = "\t"), con_phased)

  for (j in seq_len(nrow(gwas_summary))) {
    snp <- gwas_summary[j, , drop = FALSE]
    d <- internal_all[internal_all$ID == snp$ID, , drop = FALSE]
    d <- d[match(sample_ids, d$id), , drop = FALSE]
    gt_fields <- paste0(d$alt_h1, "|", d$alt_h2)
    line <- c(
      as.character(snp$CHR), as.character(snp$POS), as.character(snp$ID),
      as.character(snp$REF), as.character(snp$ALT), ".", ".", ".", "GT", gt_fields
    )
    writeLines(paste(line, collapse = "\t"), con_phased)
  }

  # MSP: one window per SNP.
  con_msp <- file(msp_path, open = "wt")
  on.exit(try(close(con_msp), silent = TRUE), add = TRUE)
  writeLines("# Simulated MSP file", con_msp)
  sample_hap_labels <- as.vector(rbind(paste0(sample_ids, ".0"), paste0(sample_ids, ".1")))
  msp_header <- c("#chm", "spos", "epos", "sgpos", "egpos", "n_snps", sample_hap_labels)
  writeLines(paste(msp_header, collapse = "\t"), con_msp)

  for (j in seq_len(nrow(gwas_summary))) {
    snp <- gwas_summary[j, , drop = FALSE]
    d <- internal_all[internal_all$ID == snp$ID, , drop = FALSE]
    d <- d[match(sample_ids, d$id), , drop = FALSE]
    msp_calls <- as.vector(rbind(d$anc1_h1, d$anc1_h2))
    pos <- as.integer(snp$POS)
    line <- c(
      as.character(snp$CHR),
      as.character(pos),
      as.character(pos + pos_step),
      "0",
      "0",
      "1",
      as.character(msp_calls)
    )
    writeLines(paste(line, collapse = "\t"), con_msp)
  }

  invisible(list(
    flare_vcf_path = flare_vcf_path,
    phased_vcf_path = phased_vcf_path,
    msp_path = msp_path,
    pheno_path = pheno_path,
    covar_path = covar_path,
    internal_data_path = if (isTRUE(write_internal_data)) internal_path else NA_character_,
    gwas_summary_path = gwas_path
  ))
}

# Final helper: run full multi-SNP pipeline and write outputs.
simulate_and_write_multi_snp <- function(
  n_individuals = 1000L,
  main_n = 500L,
  n_snps = 5L,
  prop_anc1 = c(0.20, 0.30, 0.35, 0.45, 0.55),
  maf_anc0 = c(0.08, 0.10, 0.12, 0.06, 0.14),
  maf_anc1 = c(0.18, 0.24, 0.30, 0.20, 0.28),
  n_covariates = 3L,
  beta_covar = c(0.9, -0.7, 0.5),
  beta_dosage_anc0 = c(0.10, 0.12, 0.08, 0.11, 0.09),
  beta_dosage_anc1 = c(0.15, 0.20, 0.13, 0.18, 0.14),
  beta_local_anc1 = c(0.05, 0.07, 0.06, 0.08, 0.05),
  covariate_r2 = 0.35,
  genetic_r2 = 0.03,
  local_ancestry_r2 = 0.02,
  phenotype_type = c("linear", "binary"),
  output_dir = "inst/extdata",
  prefix = "sim_multi_snp",
  write_meta = TRUE,
  write_internal_data = TRUE,
  chrom = "1",
  pos_start = 1000L,
  pos_step = 10L,
  ref = "A",
  alt = "G",
  seed = 1L
) {
  phenotype_type <- match.arg(phenotype_type)
  n_snps <- as.integer(n_snps)
  n_individuals <- as.integer(n_individuals)
  main_n <- as.integer(main_n)

  stopifnot(n_snps >= 1L, n_individuals > 2L, main_n >= 2L, main_n < n_individuals)
  write_meta <- as_logical_flag(write_meta, "write_meta")
  write_internal_data <- as_logical_flag(write_internal_data, "write_internal_data")

  set.seed(as.integer(seed))

  prop_anc1 <- ensure_length_n(prop_anc1, n_snps, "prop_anc1")
  maf_anc0 <- ensure_length_n(maf_anc0, n_snps, "maf_anc0")
  maf_anc1 <- ensure_length_n(maf_anc1, n_snps, "maf_anc1")

  sample_ids <- sprintf("id%05d", seq_len(n_individuals))
  positions <- pos_start + (seq_len(n_snps) - 1L) * as.integer(pos_step)
  snp_ids <- sprintf("rs%07d", seq_len(n_snps))

  genotype_list <- vector("list", n_snps)
  for (j in seq_len(n_snps)) {
    genotype_list[[j]] <- simulate_one_snp_genotype_local_ancestry(
      n_individuals = n_individuals,
      prop_anc1 = prop_anc1[j],
      maf_anc0 = maf_anc0[j],
      maf_anc1 = maf_anc1[j]
    )
  }

  pheno <- simulate_covariates_and_phenotype(
    genotype_list = genotype_list,
    n_covariates = n_covariates,
    beta_covar = beta_covar,
    beta_dosage_anc0 = beta_dosage_anc0,
    beta_dosage_anc1 = beta_dosage_anc1,
    beta_local_anc1 = beta_local_anc1,
    covariate_r2 = covariate_r2,
    genetic_r2 = genetic_r2,
    local_ancestry_r2 = local_ancestry_r2,
    phenotype_type = phenotype_type
  )

  internal_idx <- sort(sample.int(n_individuals, size = main_n, replace = FALSE))

  internal_by_snp <- vector("list", n_snps)
  gwas_rows <- vector("list", n_snps)
  for (j in seq_len(n_snps)) {
    res_j <- make_internal_external_one_snp(
      genotype = genotype_list[[j]],
      covariates = pheno$covariates,
      phenotype = pheno$phenotype,
      sample_ids = sample_ids,
      internal_idx = internal_idx,
      chrom = chrom,
      pos = positions[j],
      snp_id = snp_ids[j],
      ref = ref,
      alt = alt
    )
    internal_by_snp[[j]] <- res_j$internal
    gwas_rows[[j]] <- res_j$gwas_summary
  }

  gwas_summary <- do.call(rbind, gwas_rows)
  gwas_summary <- gwas_summary[order(gwas_summary$POS), , drop = FALSE]

  outputs <- write_multi_snp_internal_outputs(
    internal_by_snp = internal_by_snp,
    gwas_summary = gwas_summary,
    pheno_df = data.frame(id = sample_ids[internal_idx], y = pheno$phenotype[internal_idx], stringsAsFactors = FALSE),
    covar_df = data.frame(id = sample_ids[internal_idx], pheno$covariates[internal_idx, , drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE),
    output_dir = output_dir,
    prefix = prefix,
    write_internal_data = write_internal_data,
    pos_step = pos_step
  )

  meta_path <- NA_character_
  if (isTRUE(write_meta)) {
    meta_path <- file.path(output_dir, sprintf("%s.meta.txt", prefix))
    meta <- data.frame(
      key = c(
        "n_individuals", "main_n", "external_n", "n_snps", "phenotype_type",
        "covariate_r2", "genetic_r2", "local_ancestry_r2", "seed"
      ),
      value = c(
        n_individuals, main_n, n_individuals - main_n, n_snps, phenotype_type,
        covariate_r2, genetic_r2, local_ancestry_r2, seed
      ),
      stringsAsFactors = FALSE
    )
    utils::write.table(meta, file = meta_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  invisible(c(outputs, list(meta_path = meta_path)))
}

# ---------------------------
# Simulation parameters
# ---------------------------

n_individuals <- 8100L
main_n <- 100L
n_snps <- 5L

phenotype_type <- "linear"
prefix <- "simulate_linear"
seed <- 123L

simulate_and_write_multi_snp(
  n_individuals = n_individuals,
  main_n = main_n,
  n_snps = n_snps,
  prop_anc1 = c(0.50, 0.50, 0.50, 0.50, 0.50),
  maf_anc0 = c(0.30, 0.30, 0.40, 0.45, 0.30),
  maf_anc1 = c(0.30, 0.40, 0.30, 0.30, 0.45),
  n_covariates = 3L,
  beta_covar = c(0.9, -0.7, 0.5),
  beta_dosage_anc0 = c(0.15, 0.20, 0.13, 0.18, 0.14),
  beta_dosage_anc1 = c(0.10, 0.12, 0.08, 0.11, 0.09),
  beta_local_anc1 = c(0.05, 0.07, 0.06, 0.08, 0.05),
  covariate_r2 = 0.40,
  genetic_r2 = 0.005,
  local_ancestry_r2 = 0.003,
  phenotype_type = phenotype_type,
  output_dir = "inst/extdata",
  prefix = prefix,
  write_meta = TRUE,
  write_internal_data = FALSE,
  chrom = "chr1",
  pos_start = 1000L,
  pos_step = 10L,
  ref = "A",
  alt = "G",
  seed = seed
)

phenotype_type <- "binary"
prefix <- "simulate_binary"
seed <- 456L

simulate_and_write_multi_snp(
  n_individuals = n_individuals,
  main_n = main_n,
  n_snps = n_snps,
  prop_anc1 = c(0.50, 0.50, 0.50, 0.50, 0.50),
  maf_anc0 = c(0.30, 0.30, 0.40, 0.45, 0.30),
  maf_anc1 = c(0.30, 0.40, 0.30, 0.30, 0.45),
  n_covariates = 3L,
  beta_covar = c(0.9, -0.7, 0.5),
  beta_dosage_anc0 = c(0.15, 0.20, 0.13, 0.18, 0.14),
  beta_dosage_anc1 = c(0.10, 0.12, 0.08, 0.11, 0.09),
  beta_local_anc1 = c(0.05, 0.07, 0.06, 0.08, 0.05),
  covariate_r2 = 0.40,
  genetic_r2 = 0.005,
  local_ancestry_r2 = 0.003,
  phenotype_type = phenotype_type,
  output_dir = "inst/extdata",
  prefix = prefix,
  write_meta = TRUE,
  write_internal_data = FALSE,
  chrom = "chr1",
  pos_start = 1000L,
  pos_step = 10L,
  ref = "A",
  alt = "G",
  seed = seed
)