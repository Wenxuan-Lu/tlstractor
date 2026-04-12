#' Run TLS-Tractor
#'
#' Performs local ancestry-aware GWAS by integrating individual-level data with
#' external GWAS summary statistics via transfer learning.
#'
#' @param gds_path Character; path to the input GDS file.
#' @param sumstats_path Character; path to the munged summary statistics file.
#' @param method Character; one of `"linear"` or `"logistic"`.
#' @param cond_local Logical; whether to condition on local-ancestry dosage
#'   terms.
#' @param pheno_path Character; path to the phenotype file.
#' @param pheno_id_col Character; sample-ID column name in the phenotype file.
#' @param pheno_col Character; phenotype column name in the phenotype file.
#' @param covar_path Character or `NULL`; path to an optional covariate file.
#'   Default is `NULL`.
#' @param covar_id_col Character or `NULL`; sample-ID column name in the
#'   covariate file (required when `covar_path` is provided). Default is `NULL`.
#' @param covar_cols Character vector or `NULL`; covariate column names in the
#'   covariate file (required when `covar_path` is provided). Default is `NULL`.
#' @param output_prefix Character; output path prefix.
#' @param scratch_dir Character or `NULL`; directory used to temporarily store
#'   per-task result files before merge. If `NULL`, a run-specific directory is
#'   created under `dirname(output_prefix)` with name
#'   `<basename(output_prefix)>_<pid>_<YYYYmmdd_HHMMSS>_tmp`. Default is `NULL`.
#' @param snp_start Integer; 1-based starting SNP index in the GDS file.
#'   Default is `1L`.
#' @param snp_count Integer or `NULL`; number of SNPs to process. If `NULL`,
#'   processing starts at `snp_start` and continues to the end of the GDS file.
#'   Default is `NULL`.
#' @param n_cores Integer; number of CPU cores to use. Default is `1L`.
#' @param chunk_size Integer; number of variants processed per chunk within each
#'   CPU core. Default is 1024; larger values may improve speed but increase
#'   memory usage.
#' @param local_ancestry_mac_threshold Integer; ancestry-specific minor-allele
#'   count threshold. SNPs with any ancestry-specific MAC below this threshold
#'   are skipped. Skipped SNPs are returned in the output with `NA` values for
#'   inferential columns. Default is `20L`.
#' @param use_fast_version Logical; whether to use the fast TLS-Tractor mode.
#'   The fast mode assumes that the estimated non-genetic covariate effects from
#'   the null model (phenotype ~ covariates) are close to those from the full
#'   standard GWAS model (phenotype ~ SNP dosage + covariates) for any single
#'   SNP, so that estimated covariate effects from the null model can be reused
#'   across variants to reduce computation. The fast mode can provide a
#'   substantial speedup with minimal loss of accuracy. Default is `TRUE`
#'   (recommended for large datasets). Setting `FALSE` fources the full
#'   per-variant fitting path, which is generally more robust but substantially
#'   slower when covariates are included. If no covariates are provided, fast
#'   mode is automatically disabled.
#'
#' @return Invisibly returns `NULL`. Writes gzipped GWAS results to
#'   `<output_prefix>.txt.gz`. The output includes:
#'   \itemize{
#'     \item Variant metadata: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `main_N`
#'       (sample size in the main/internal study)
#'     \item Frequency and ancestry summaries: `AF` (overall allele frequency),
#'       `AF_anc*` (ancestry-specific allele frequency), `LAprop_anc*`
#'       (local-ancestry-specific haplotype proportion)
#'     \item Analysis metadata: `has_sumstats` (indicates whether external
#'       summary statistics were available), `fallback_used` (indicates whether
#'       QR decomposition fallback was used; when TRUE, results may have reduced
#'       numerical stability and should be interpreted with caution)
#'     \item Association results: `joint_pval` (p-value for joint testing all
#'       ancestry-specific SNP dosage effects), `beta_anc*`, `se_anc*`,
#'       `pval_anc*` (effect size, standard error, and p-value for each
#'       ancestry-specific SNP dosage effect)
#'     \item When `cond_local = TRUE`: `LAeff_anc*`, `LAse_anc*`, `LApval_anc*`
#'       (effect size, standard error, and p-value for each local ancestry term)
#'   }
#'
#' @details
#' If present, output the file `<output_prefix>.excluded_samples.txt` containing
#' sample IDs present in the GDS file but excluded from analysis after sample
#' intersection/filtering.
#'
#' During execution, temporary per-task result files are written to
#' `scratch_dir/<run_tag>_task_<id>.txt.gz`, where `run_tag` is a run-specific
#' identifier with format `tlstractor_<pid>_<YYYYmmdd_HHMMSS>`. These temporary
#' files are removed on successful cleanup. The directory `scratch_dir` is
#' removed only if it was created by the function.
#'
#'
#' SNPs present in the GDS file but absent from the summary statistics file are
#' still analyzed using the Tractor model with individual-level data only.
#' 
#' To perform Tractor anlysis on all SNPs using only the main study,
#' set `sumstats_path` to a file with one tab-delimited header line containing:
#' `CHR`, `POS`, `ID`, `REF`, `ALT`, `BETA`, `SE`, and `GDS_ID`.
#'
#' @export
tlstractor <- function(gds_path, sumstats_path, method, cond_local,
                       pheno_path, pheno_id_col, pheno_col,
                       covar_path=NULL, covar_id_col=NULL, covar_cols=NULL,
                       output_prefix, scratch_dir=NULL,
                       snp_start=1L, snp_count=NULL, n_cores=1L, chunk_size=1024L,
                       local_ancestry_mac_threshold=20L,
                       use_fast_version=TRUE) {
    # Validate parameters
    validate_param_type(gds_path, "character", "gds_path", check_length_one = TRUE)
    if (!file.exists(gds_path)) stop("Input GDS file not found.", call. = FALSE)

    validate_param_type(sumstats_path, "character", "sumstats_path", check_length_one = TRUE)
    if (!file.exists(sumstats_path)) stop("Summary statistics file not found.", call. = FALSE)
    validate_col_names(sumstats_path, "GDS_ID", c("CHR", "POS", "ID", "REF", "ALT", "BETA", "SE"), "summary statistics file")

    validate_param_type(method, "character", "method", check_length_one = TRUE)
    method <- tolower(method)
    if (!method %in% c("linear", "logistic")) {
      stop("Error: method must be 'linear' or 'logistic'.", call. = FALSE)
    }

    validate_param_type(cond_local, "logical", "cond_local", check_length_one = TRUE)

    validate_param_type(pheno_path, "character", "pheno_path", check_length_one = TRUE)
    validate_param_type(pheno_id_col, "character", "pheno_id_col", check_length_one = TRUE)
    validate_param_type(pheno_col, "character", "pheno_col", check_length_one = TRUE)
    if (!file.exists(pheno_path)) stop("Phenotype file not found.", call. = FALSE)
    validate_col_names(pheno_path, pheno_id_col, pheno_col, "phenotype file")
    
    has_covar <- !is.null(covar_path)
    if (has_covar) {
        validate_param_type(covar_path, "character", "covar_path", check_length_one = TRUE)
        validate_param_type(covar_id_col, "character", "covar_id_col", check_length_one = TRUE)
        validate_param_type(covar_cols, "character", "covar_cols", check_length_one = FALSE)
        if (!file.exists(covar_path)) stop("Covariate file not found.", call. = FALSE)
        validate_col_names(covar_path, covar_id_col, covar_cols, "covariate file")
    }

    # TODO: weights will be added in the next update
    weight_path <- NULL # path to an optional sample weights file; if NULL, no weights are used
    weight_id_col <- NULL # sample ID column name in the weights file; required if weight_path is not NULL
    weight_col <- NULL # weights column name in the weights file; required if weight_path is not NULL
    has_weight <- !is.null(weight_path)
    if (has_weight) {
        validate_param_type(weight_path, "character", "weight_path", check_length_one = TRUE)
        validate_param_type(weight_id_col, "character", "weight_id_col", check_length_one = TRUE)
        validate_param_type(weight_col, "character", "weight_col", check_length_one = TRUE)
        if (!file.exists(weight_path)) stop("Weights file not found.", call. = FALSE)
        validate_col_names(weight_path, weight_id_col, weight_col, "weight file")
    }

    validate_param_type(output_prefix, "character", "output_prefix", check_length_one = TRUE)
    output_prefix <- normalizePath(trimws(output_prefix), mustWork=FALSE)
    out_dir <- dirname(output_prefix)
    if (!dir.exists(out_dir)) stop(sprintf("Error: Output directory does not exist:\n  %s", out_dir), call. = FALSE)

    pkg <- "tlstractor"
    pid <- Sys.getpid()
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    run_tag <- sprintf("%s_%d_%s", pkg, pid, timestamp)
    if (!is.null(scratch_dir)) {
        validate_param_type(scratch_dir, "character", "scratch_dir", check_length_one = TRUE)
        scratch_dir <- normalizePath(trimws(scratch_dir), mustWork = FALSE)
    } else {
        scratch_dir <- file.path(out_dir, sprintf("%s_%d_%s_tmp", basename(output_prefix), pid, timestamp))
    }

    validate_param_type(snp_start, "integer", "snp_start", check_length_one = TRUE)
    snp_start <- as.integer(snp_start)
    if (snp_start <= 0L) stop("Error: snp_start must be a positive integer.", call. = FALSE)
    if (!is.null(snp_count)) {
        validate_param_type(snp_count, "integer", "snp_count", check_length_one = TRUE)
        snp_count <- as.integer(snp_count)
        if (snp_count <= 0L) stop("Error: snp_count must be a positive integer.", call. = FALSE)
    }
    validate_param_type(n_cores, "integer", "n_cores", check_length_one = TRUE)
    n_cores <- as.integer(n_cores)
    if (n_cores <= 0L) stop("Error: n_cores must be a positive integer.", call. = FALSE)
    validate_param_type(chunk_size, "integer", "chunk_size", check_length_one = TRUE)
    chunk_size <- as.integer(chunk_size)
    if (chunk_size <= 0L) stop("Error: chunk_size must be a positive integer.", call. = FALSE)

    validate_param_type(local_ancestry_mac_threshold, "integer", "local_ancestry_mac_threshold", check_length_one = TRUE)
    local_ancestry_mac_threshold <- as.integer(local_ancestry_mac_threshold)
    if (local_ancestry_mac_threshold < 0L) stop("Error: local_ancestry_mac_threshold must be non-negative.", call. = FALSE)

    validate_param_type(use_fast_version, "logical", "use_fast_version", check_length_one = TRUE)

    # Internal parameters
    use_offset <- TRUE
    refine_C <- FALSE
    sd_tol <- 1e-8
    use_qr_fallback <- TRUE
    refine_W <- FALSE

    local_ancestry_mac_soft_threshold <- 100L # Placeholder; overwritten later.
    validate_param_type(local_ancestry_mac_soft_threshold, "integer", "local_ancestry_mac_soft_threshold", check_length_one = TRUE)
    local_ancestry_mac_soft_threshold <- as.integer(local_ancestry_mac_soft_threshold)
    if (local_ancestry_mac_soft_threshold < 0L) stop("Error: local_ancestry_mac_soft_threshold must be non-negative.", call. = FALSE)

    control <- list(
        local_ancestry_mac_hard_threshold = local_ancestry_mac_threshold,
        local_ancestry_mac_soft_threshold = local_ancestry_mac_soft_threshold,
        use_offset = use_offset,
        refine_C = refine_C,
        use_qr_fallback = use_qr_fallback,
        refine_W = refine_W
    )
    
    # Load GDS metadata and determine SNP range
    message("Reading GDS metadata...")
    in_gds <- gdsfmt::openfn.gds(gds_path, readonly = TRUE)
    on.exit({
        try({
            if (inherits(in_gds, "gds.class")) {
                gdsfmt::closefn.gds(in_gds)
                in_gds <- NULL
            }
        }, silent = TRUE)
    }, add = TRUE)

    sample_ids_gds <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(in_gds, "sample.id"))
    nsamples_total <- length(sample_ids_gds)
    if (nsamples_total <= 1L) stop("Error: at least two samples need to be present.", call. = FALSE)
    message("Total samples in GDS: ", nsamples_total)

    chrom_node <- gdsfmt::index.gdsn(in_gds, "snp.chromosome")
    nsnp <- gdsfmt::objdesp.gdsn(chrom_node)$dim
    if (is.null(nsnp)) nsnp <- gdsfmt::objdesp.gdsn(chrom_node)$size
    nsnp <- as.integer(nsnp)
    if (nsnp <= 0L) stop("Error: no SNPs found in GDS file.", call. = FALSE)
    message("Total SNPs in GDS: ", nsnp)

    num_ancs <- length(gdsfmt::ls.gdsn(gdsfmt::index.gdsn(in_gds, "dosage")))
    if (is.null(num_ancs) || is.na(num_ancs) || num_ancs < 2L) {
        stop("Error: at least two ancestries are required.", call. = FALSE)
    }

    for (k in 0:(num_ancs - 1L)) {
        dos_node <- gdsfmt::index.gdsn(in_gds, paste0("dosage/anc", k))
        dos_dim <- gdsfmt::objdesp.gdsn(dos_node)$dim
        if (is.null(dos_dim) || length(dos_dim) != 2L || dos_dim[1] != nsamples_total || dos_dim[2] != nsnp) {
            stop(sprintf("Error: dosage/anc%d has incorrect dimensions. Expected [%d, %d], got [%s]",
                       k, nsamples_total, nsnp, paste(dos_dim, collapse = ", ")), call. = FALSE)
        }
        
        hap_node <- gdsfmt::index.gdsn(in_gds, paste0("hapcount/anc", k))
        hap_dim <- gdsfmt::objdesp.gdsn(hap_node)$dim
        if (is.null(hap_dim) || length(hap_dim) != 2L || hap_dim[1] != nsamples_total || hap_dim[2] != nsnp) {
            stop(sprintf("Error: hapcount/anc%d has incorrect dimensions. Expected [%d, %d], got [%s]",
                       k, nsamples_total, nsnp, paste(hap_dim, collapse = ", ")), call. = FALSE)
        }
    }

    if (snp_start > nsnp) stop("Error: snp_start exceeds number of SNPs in GDS file.", call. = FALSE)
    if (is.null(snp_count)) {
        snp_end <- nsnp
    } else {
        snp_end <- min(nsnp, snp_start + snp_count - 1L)
    }
    
    # Load phenotype, covariates, and weights
    message("Loading phenotype data...")
    pheno_dt <- data.table::fread(pheno_path, select=c(pheno_id_col, pheno_col),
                    colClasses=setNames(c("character", "numeric"), c(pheno_id_col, pheno_col)),
                    showProgress=FALSE)
    data.table::setnames(pheno_dt, c(pheno_id_col, pheno_col), c("id","pheno"))
    if (anyDuplicated(pheno_dt$id)) stop("Error: Duplicate sample IDs found in phenotype file", call. = FALSE)
    n_orig <- nrow(pheno_dt)
    pheno_dt <- pheno_dt[is.finite(pheno)]
    if (nrow(pheno_dt) == 0L) stop("Error: All phenotype values are NA or Inf.", call. = FALSE)
    n_removed <- n_orig - nrow(pheno_dt)
    if (n_removed > 0L) message(sprintf("Excluding %d individuals with NA or Inf phenotype values", n_removed))
    pheno_range <- range(pheno_dt$pheno)
    if (pheno_range[1] == pheno_range[2]) stop("Error: phenotype has zero variance.", call. = FALSE)

    if (has_covar) {
        message("Loading covariate data...")
        covar_dt <- data.table::fread(covar_path, select=c(covar_id_col, covar_cols),
                        colClasses=setNames(c("character", rep("numeric", length(covar_cols))), c(covar_id_col, covar_cols)),
                        showProgress=FALSE)
        data.table::setnames(covar_dt, covar_id_col, "id")
        if (anyDuplicated(covar_dt$id)) stop("Error: Duplicate sample IDs found in covariate file.", call. = FALSE)
        covar_mat <- as.matrix(covar_dt[, ..covar_cols])
        storage.mode(covar_mat) <- "double"
        ok <- rowSums(!is.finite(covar_mat)) == 0L
        covar_ids <- covar_dt$id[ok]
        covar_mat <- covar_mat[ok, , drop = FALSE]
        dimnames(covar_mat) <- NULL
        if (nrow(covar_mat) == 0L) stop("Error: All covariate rows contain NA or Inf values.", call. = FALSE)
        n_removed <- sum(!ok)
        if (n_removed > 0L) message(sprintf("Excluding %d individuals with NA or Inf covariate values", n_removed))
        covar_mins <- matrixStats::colMins(covar_mat)
        covar_maxs <- matrixStats::colMaxs(covar_mat)
        zero_var <- (covar_mins == covar_maxs)
        if (any(zero_var)) {
            message(sprintf("Removing %d covariate column(s) with zero variance: %s",
                          sum(zero_var), paste(covar_cols[zero_var], collapse = ", ")))
            covar_cols <- covar_cols[!zero_var]
            covar_mat <- covar_mat[, !zero_var, drop = FALSE]
            if (length(covar_cols) == 0L) {
                has_covar <- FALSE
                covar_mat <- NULL
                covar_ids <- NULL
                warning("All covariate columns have zero variance. Proceeding without covariates.", call. = FALSE)
            }
        }
        rm(covar_dt, ok, covar_mins, covar_maxs, zero_var)
    }

    if (has_weight) {
        message("Loading sample weights...")
        weight_dt <- data.table::fread(weight_path, select=c(weight_id_col, weight_col),
                        colClasses=setNames(c("character", "numeric"), c(weight_id_col, weight_col)),
                        showProgress=FALSE)
        data.table::setnames(weight_dt, c(weight_id_col, weight_col), c("id","weight"))
        if (anyDuplicated(weight_dt$id)) stop("Error: Duplicate sample IDs found in weights file.", call.=FALSE)
        if (any(!is.finite(weight_dt$weight)) || any(weight_dt$weight < 0)) {
            stop("Error: Weights contain invalid values. All weights must be non-negative finite numbers.", call.=FALSE)
        }
    }

    # Filter samples to those present in GDS
    message("Filtering samples...")
    keep <- data.table::`%chin%`(sample_ids_gds, pheno_dt$id)
    if (has_covar)  keep <- keep & data.table::`%chin%`(sample_ids_gds, covar_ids)
    if (has_weight) keep <- keep & data.table::`%chin%`(sample_ids_gds, weight_dt$id)
    ids_in_gds <- which(keep)
    common_ids <- sample_ids_gds[keep]
    nsamples <- length(common_ids)
    if (nsamples <= 1L) stop("Error: at least two samples need to be present after intersection.", call. = FALSE)
    message("Samples after filtering: ", nsamples)

    if (nsamples_total - nsamples > 0L) {
        excluded_ids <- sample_ids_gds[!keep]
        exclude_sample_ids_path <- paste0(output_prefix, ".excluded_samples.txt")
        writeLines(excluded_ids, exclude_sample_ids_path)
        message(sprintf("Excluded %d samples from the GDS file. Sample IDs written to: %s", length(excluded_ids), exclude_sample_ids_path))
        rm(excluded_ids)
    }

    idx <- match(common_ids, pheno_dt$id)
    pheno_vec <- pheno_dt$pheno[idx]
    pheno_vec <- as.double(pheno_vec)

    if (has_covar) {
        idx <- match(common_ids, covar_ids)
        covar_mat <- covar_mat[idx, , drop = FALSE]
        rm(covar_ids)
    } else {
        covar_mat <- NULL # has_covar == FALSE
    }

    if (has_weight) {
        idx <- match(common_ids, weight_dt$id)
        weights_vec <- weight_dt$weight[idx]
        weights_vec <- as.double(weights_vec)
        rm(weight_dt)
    } else {
        weights_vec <- NULL # has_weight == FALSE
    }

    rm(sample_ids_gds, pheno_dt, idx, common_ids, keep)

    # Sanity checks and preprocessing
    message("Preprocessing data...")
    if (has_covar) {
        covar_mins <- matrixStats::colMins(covar_mat)
        covar_maxs <- matrixStats::colMaxs(covar_mat)
        zero_var <- (covar_mins == covar_maxs)
        if (any(zero_var)) {
            message(sprintf("Removing %d covariate column(s) with zero variance after filtering: %s",
                          sum(zero_var), paste(covar_cols[zero_var], collapse = ", ")))
            covar_cols <- covar_cols[!zero_var]
            covar_mat <- covar_mat[, !zero_var, drop = FALSE]
            if (length(covar_cols) == 0L) {
                has_covar <- FALSE
                covar_mat <- NULL
                warning("All covariates have zero variance after filtering. Proceeding without covariates.", call. = FALSE)
            }
        }
        rm(covar_mins, covar_maxs, zero_var)
    }

    if (has_covar) {
        covar_mat <- scale(covar_mat, center = TRUE, scale = TRUE)
        covar_sds <- attr(covar_mat, "scaled:scale")
        small_sd_cols <- which(!is.finite(covar_sds) | covar_sds < sd_tol)
        if (length(small_sd_cols) > 0L) {
            warning(sprintf("Covariate column(s) with very small/invalid SD after filtering. Proceeding with standardization, but results may be unstable. Consider removing these covariates: %s", 
                            paste(small_sd_cols, collapse = ", ")), call. = FALSE)
        }
        if (!all(is.finite(covar_mat))) {
            stop("Error: Covariate matrix contains non-finite values after standardization.", call. = FALSE)
        }
        attr(covar_mat, "scaled:center") <- NULL
        attr(covar_mat, "scaled:scale") <- NULL
        rm(covar_sds, small_sd_cols)
    }

    covar_cols <- NULL
    if (method == "logistic") {
        if (has_covar) {
            covar_mat <- cbind(1.0, covar_mat)
        } else {
            covar_mat <- matrix(1.0, nrow = nsamples, ncol = 1L)
            has_covar <- TRUE
        }
    }

    min_samples_required <- 2L + num_ancs + if (has_covar) ncol(covar_mat) else 0L + if (cond_local) (num_ancs - 1L) else 0L
    if (nsamples < min_samples_required) {
        stop(sprintf("Error: Not enough samples after filtering to fit the model. Need at least %d samples, but got %d.", min_samples_required, nsamples), call. = FALSE)
    }
    if (has_weight) {
        positive_weight <- weights_vec > 0
        n_pos_weight <- sum(positive_weight)
        if (n_pos_weight < min_samples_required) {
            stop(sprintf("Error: Not enough samples with positive weight to fit the model. Need at least %d samples with positive weight, but got %d.", min_samples_required, n_pos_weight), call. = FALSE)
        }
    }

    if (!has_covar && use_fast_version) {
        use_fast_version <- FALSE
        message("No covariates provided. Disabling fast mode and falling back to default mode.")
    }
    if (use_fast_version) {
        local_ancestry_mac_soft_threshold <- 0L
        control$local_ancestry_mac_soft_threshold <- local_ancestry_mac_soft_threshold
    } else {
        local_ancestry_mac_soft_threshold <- nsamples * 2 + 1L
        control$local_ancestry_mac_soft_threshold <- local_ancestry_mac_soft_threshold
    }

    if (method == "linear") {
        pheno_range <- range(pheno_vec)
        if (pheno_range[1] == pheno_range[2]) stop("Error: phenotype has zero variance.", call. = FALSE)
        pheno_mean <- mean(pheno_vec)
        pheno_sd <- stats::sd(pheno_vec)
        pheno_vec <- (pheno_vec - pheno_mean) / pheno_sd
        if (pheno_sd < sd_tol) {
            warning("Phenotype has near-zero standard deviation after filtering. Proceeding with standardization, but results may be unstable.", call. = FALSE)
        }
        if (!all(is.finite(pheno_vec))) {
            stop("Error: Phenotype vector contains non-finite values after standardization.", call. = FALSE)
        }
        if (has_weight) {
            rm(positive_weight)
        }
        linear_precomp_error <- "Error: Linear precomputation (pheno ~ covariates) failed due to high covariate collinearity."
        if (has_covar) {
            precomp <- tlstractor_linear_precompute(pheno_vec, covar_mat)
            if (isTRUE(precomp$rank < 1L) || isTRUE(precomp$R1_not_ok)) {
                stop(linear_precomp_error, call. = FALSE)
            }
            validate_precomp_fields(precomp, linear_precomp_error)
        } else {
            precomp <- list(y = pheno_vec)
        }
    } else {
        pheno_unique <- unique(pheno_vec)
        if (length(pheno_unique) != 2L || !all(pheno_unique %in% c(0, 1))) {
            stop("Error: logistic method requires phenotype values 0/1, and both values must be present.", call. = FALSE)
        }
        if (has_weight) {
            pheno_0_weight <- any(pheno_vec == 0 & positive_weight)
            pheno_1_weight <- any(pheno_vec == 1 & positive_weight)
            if (!pheno_0_weight || !pheno_1_weight) {
                stop("Error: logistic method requires both 0 and 1 among samples with positive weight.", call. = FALSE)
            }
            rm(positive_weight, pheno_0_weight, pheno_1_weight)
        }
        logistic_precomp_error <- "Error: Logistic precomputation (pheno ~ covariates) failed due to near-separation, covariate collinearity, or non-convergence."
        precomp <- tlstractor_logistic_precompute(pheno_vec, covar_mat)
        if (!isTRUE(precomp$converged) || isTRUE(precomp$error) || isTRUE(precomp$L_AtWA_not_ok) || isTRUE(precomp$rank < 1L)) {
            stop(logistic_precomp_error, call. = FALSE)
        }
        validate_precomp_fields(precomp, logistic_precomp_error)
    }

    # For weighted analyses, a weighted version of precomputation is needed.

    if (has_covar) rm(covar_mat)
    if (has_weight) rm(weights_vec)
    rm(pheno_vec)
    
    # Load summary statistics
    message("Loading summary statistics...")
    sumstats <- data.table::fread(sumstats_path, select = c("GDS_ID", "BETA", "SE"),
                    colClasses = c(GDS_ID = "integer", BETA = "numeric", SE = "numeric"),
                    showProgress = FALSE)

    # Set number of cores
    avail_cores <- tryCatch(
        parallelly::availableCores(),
        error = function(e) NA_integer_
    )
    if (is.na(avail_cores) || is.null(avail_cores) || length(avail_cores) != 1L || avail_cores <= 0L) {
        avail_cores <- parallel::detectCores(logical = FALSE)
    }
    if (is.null(avail_cores) || is.na(avail_cores) || length(avail_cores) != 1L || avail_cores <= 0L) {
        stop("Error: failed to detect available CPU cores.", call. = FALSE)
    }
    max_cores <- max(1L, avail_cores - 1L)
    if (n_cores > avail_cores) stop(sprintf("Error: n_cores (%d) exceeds max usable cores (%d).", n_cores, avail_cores), call. = FALSE)
    n_cores <- min(n_cores, max_cores)

    # Determine task chunking
    nsnp_to_process <- snp_end - snp_start + 1L
    min_snp_per_chunk <- 4096L
    max_snp_per_chunk <- ceiling(nsnp_to_process / n_cores)
    max_tasks <- max(n_cores * 10L, 1000L)

    snp_per_task <- chunk_size
    if (snp_per_task > max_snp_per_chunk) {
        message(sprintf("Adjusting chunk_size from %d to %d to ensure enough tasks for %d cores.", chunk_size, max_snp_per_chunk, n_cores))
        chunk_size <- max_snp_per_chunk
        snp_per_task <- max_snp_per_chunk
        num_tasks <- ceiling(nsnp_to_process / snp_per_task)
    } else if (min_snp_per_chunk >= max_snp_per_chunk) {
        snp_per_task <- max_snp_per_chunk
        num_tasks <- ceiling(nsnp_to_process / snp_per_task)       
    } else {
        if (snp_per_task < min_snp_per_chunk) {
            snp_per_task <- min_snp_per_chunk
        }
        num_tasks <- ceiling(nsnp_to_process / snp_per_task)
        if (num_tasks > max_tasks) {
            snp_per_task <- ceiling(nsnp_to_process / max_tasks)
            num_tasks <- ceiling(nsnp_to_process / snp_per_task)
        }
    }

    # Create task list
    starts <- seq.int(from = snp_start, to = snp_end, by = snp_per_task)
    ends <- pmin.int(starts + snp_per_task - 1L, snp_end)
    sumstats_gds_ids <- sumstats$GDS_ID
    n_sumstats <- nrow(sumstats)
    sumstats_idx_start <- integer(num_tasks)
    sumstats_idx_end <- integer(num_tasks)
    curr_sumstats_idx <- 1L
    for (i in seq_len(num_tasks)) {
        curr_start <- starts[i]
        curr_end <- ends[i]
        while (curr_sumstats_idx <= n_sumstats && sumstats_gds_ids[curr_sumstats_idx] < curr_start) {
            curr_sumstats_idx <- curr_sumstats_idx + 1L
        }
        if (curr_sumstats_idx > n_sumstats) {
            sumstats_idx_start[i] <- -1L
            sumstats_idx_end[i] <- -1L
            next
        }
        sumstats_idx_start[i] <- curr_sumstats_idx
        while (curr_sumstats_idx <= n_sumstats && sumstats_gds_ids[curr_sumstats_idx] <= curr_end) {
            curr_sumstats_idx <- curr_sumstats_idx + 1L
        }
        sumstats_idx_end[i] <- curr_sumstats_idx - 1L
        if (sumstats_idx_end[i] < sumstats_idx_start[i]) {
            sumstats_idx_start[i] <- -1L
            sumstats_idx_end[i] <- -1L
        }
    }
    tasks <- vector("list", num_tasks)
    for (i in seq_len(num_tasks)) {
        if (sumstats_idx_start[i] != -1L) {
            idx <- sumstats_idx_start[i]:sumstats_idx_end[i]
            tasks[[i]] <- list(
                id = i,
                start = starts[i],
                end = ends[i],
                sumstats_gds_id = sumstats_gds_ids[idx],
                sumstats_beta = sumstats$BETA[idx],
                sumstats_se = sumstats$SE[idx]
            )
        } else {
        tasks[[i]] <- list(
            id = i,
            start = starts[i],
            end = ends[i],
            sumstats_gds_id = integer(0),
            sumstats_beta = numeric(0),
            sumstats_se = numeric(0)
        )
        }
    }
    idx <- NULL
    rm(sumstats, starts, ends, sumstats_gds_ids, sumstats_idx_start, sumstats_idx_end)
    
    # Prepare scratch directory and output path
    created_scratch_dir <- FALSE
    if (!dir.exists(scratch_dir)) {
        created_scratch_dir <- TRUE
        ok <- dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
        if (!ok && !dir.exists(scratch_dir)) stop(sprintf("Failed to create the scratch directory:\n  %s", scratch_dir), call. = FALSE)
    }
    output_path <- paste0(output_prefix, ".txt.gz")
    task_file_prefix <- file.path(scratch_dir, paste0(run_tag, "_task_"))
    task_file_suffix <- ".txt.gz"
    
    # Log parallel setup information
    message(sprintf(
        paste0(
            "Parallel setup: using %d of %d available cores. ",
            "Planned %d tasks to process SNPs [%d, %d] (%d total). ",
            "Each task handles up to %d SNPs (last task may be smaller). ",
            "Genotypes are read in chunks of %d SNPs within each task.\n",
            "Scratch directory: %s\n",
            "Output filepath: %s"
        ),
        n_cores, avail_cores,
        num_tasks, snp_start, snp_end, nsnp_to_process,
        snp_per_task, chunk_size,
        scratch_dir, output_path
    ))

    # Decide on genotype reading strategy
    sample_keep_frac <- nsamples / nsamples_total
    runs <- 1L + sum(diff(ids_in_gds) != 1L)
    holeyness <- runs / nsamples
    keep_frac_low <- 0.20
    keep_frac_high <- 0.60
    holeyness_mid_threshold <- 0.12

    if (sample_keep_frac <= keep_frac_low) {
        read_gds_all_then_subset <- FALSE
    } else if (sample_keep_frac >= keep_frac_high) {
        read_gds_all_then_subset <- TRUE
    } else {
        read_gds_all_then_subset <- (holeyness >= holeyness_mid_threshold)
    }

    # Initialize parallel cluster
    message("Initializing parallel cluster...")
    thread_env_names <- c(
        "OMP_NUM_THREADS",
        "OMP_THREAD_LIMIT",
        "OPENBLAS_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "BLIS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS"
    )
    thread_env_vals <- rep.int("1", length(thread_env_names))
    names(thread_env_vals) <- thread_env_names

    cl <- parallelly::makeClusterPSOCK(
        n_cores,
        outfile = "",
        rscript_envs = as.list(thread_env_vals)
    ) # outfile=NULL
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    
    # With load balancing, tasks are assigned dynamically, so worker-based RNG
    # streams may not give reproducible task-level results across runs.
    # If random seed is needed later, seed per task instead of per worker.
    # parallel::clusterSetRNGStream(cl, iseed = 1)

    # Initialize per-worker state once and persist it on each worker.

    close_worker_gds_handles <- function(cl, pkg) {
        try(
            parallel::clusterCall(
                cl,
                function(pkg) {
                    ns <- asNamespace(pkg)
                    pkg_env <- get(".tlstractor_env", envir = ns, inherits = FALSE)
                    # No thread/env restoration is needed here: workers are
                    # short-lived PSOCK processes and are terminated at
                    # cluster shutdown, so process-local settings do not leak.

                    st <- pkg_env$worker_state
                    if (is.list(st) && !is.null(st$g) && inherits(st$g, "gds.class")) {
                        gdsfmt::closefn.gds(st$g)
                        st$g <- NULL
                    }

                    pkg_env$worker_state <- NULL
                    pkg_env$run_task_fn <- NULL
                    NULL
                },
                pkg = pkg
            ),
            silent = TRUE
        )
        invisible(NULL)
    }
    on.exit(close_worker_gds_handles(cl, pkg), add = TRUE)

    parallel::clusterCall(
        cl,
        function(pkg, gds_path, ids_in_gds, precomp, control, num_ancs, has_covar, method, cond_local, chunk_size, task_file_prefix, task_file_suffix, read_gds_all_then_subset) {
            # Load package namespace on workers
            # loadNamespace() automatically loads all Imports dependencies from DESCRIPTION
            ns <- loadNamespace(pkg)
            pkg_env <- get(".tlstractor_env", envir = ns, inherits = FALSE)

            # Env thread limits are set at worker launch via rscript_envs.

            data.table::setDTthreads(1L)

            if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
                try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
                try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
            }

            # Access worker_init via pkg::: (internal function)
            worker_init_fn <- getFromNamespace("worker_init", pkg)

            # Save worker state in package-private worker storage so it persists
            # across all parLapplyLB task invocations handled by that worker.
            pkg_env$worker_state <- worker_init_fn(
                gds_path = gds_path,
                ids_in_gds = ids_in_gds,
                precomp = precomp,
                control = control,
                num_ancs = num_ancs,
                has_covar = has_covar,
                method = method,
                cond_local = cond_local,
                chunk_size = chunk_size,
                task_file_prefix = task_file_prefix,
                task_file_suffix = task_file_suffix,
                read_gds_all_then_subset = read_gds_all_then_subset
            )

            # Cache run_task function once per worker to avoid repeated lookups.
            pkg_env$run_task_fn <- get("run_task", envir = ns, inherits = FALSE)
            NULL
        },
        pkg = pkg,
        gds_path = gds_path,
        ids_in_gds = ids_in_gds,
        precomp = precomp,
        control = control,
        num_ancs = num_ancs,
        has_covar = has_covar,
        method = method,
        cond_local = cond_local,
        chunk_size = chunk_size,
        task_file_prefix = task_file_prefix,
        task_file_suffix = task_file_suffix,
        read_gds_all_then_subset = read_gds_all_then_subset
    )
    
    # Run tasks in parallel with dynamic scheduling
    message("Running TLS-Tractor analysis in parallel...")
    results <- parallel::parLapply(
        cl,
        X = tasks,
        fun = function(task, pkg) {
            tryCatch({
                ns <- asNamespace(pkg)
                pkg_env <- get(".tlstractor_env", envir = ns, inherits = FALSE)

                st <- pkg_env$worker_state
                if (!is.list(st)) {
                    stop("Worker state is missing; worker initialization may have failed")
                }

                run_task_fn <- pkg_env$run_task_fn
                if (!is.function(run_task_fn)) {
                    stop("run_task function cache is missing on worker; worker initialization may have failed")
                }

                run_task_fn(task, st)
                # Return NULL on success
                NULL
            }, error = function(e) {
                msg <- sprintf("Worker %d (task %d): %s", Sys.getpid(), task$id, conditionMessage(e))
                warning(msg, call. = FALSE)
                # Return error object on failure
                structure(list(error = TRUE, message = msg, task_id = task$id), class = "task-error")
            })
        },
        pkg = pkg
    )

    # Close worker GDS handles as soon as SNP tasks finish.
    close_worker_gds_handles(cl, pkg)

    # Check for task errors
    error_results <- Filter(function(r) inherits(r, "task-error"), results)
    if (length(error_results) > 0L) {
        error_ids <- sort(vapply(error_results, `[[`, integer(1), "task_id"))
        warning(
            sprintf("%d task(s) failed. Failed task ID(s): %s. See warnings() for details.",
                    length(error_ids), paste(error_ids, collapse = ", ")),
            call. = FALSE
        )
    }

    # Merge results
    message("Merging results...")
    merge_task_results_cpp(output_path, task_file_prefix, task_file_suffix, num_tasks)

    # Cleanup
    message("Cleaning up temporary files...")
    if (created_scratch_dir) {
        suppressWarnings({
            ok <- unlink(scratch_dir, recursive = TRUE, force = TRUE)
            if (ok != 0) message("Note: could not completely remove scratch directory: ", scratch_dir)
        })
    } else {
        pattern <- paste0("^", run_tag)
        files <- list.files(scratch_dir, pattern = pattern, full.names = TRUE, all.files = TRUE)
        if (length(files)) {
            suppressWarnings({
                ok <- unlink(files, recursive = TRUE, force = TRUE)
                if (any(ok != 0)) message("Note: could not remove some temporary files in: ", scratch_dir)
            })
        }
    }

    message("TLS-Tractor analysis complete!")
    message("Results written to: ", output_path)
}


#' @keywords internal
#' @noRd
validate_col_names <- function(file_path, id_col, value_cols, file_context) {
    if (!identical(id_col, make.names(id_col))) {
    stop(sprintf("Error: Invalid ID column name in %s: '%s'. Use a valid R column name (letters, numbers, dots, or underscores; no leading digit).", file_context, id_col), call. = FALSE)
    }

    for (col in value_cols) {
        if (!identical(col, make.names(col))) {
            stop(sprintf("Error: Invalid column name in %s: '%s'. Use a valid R column name (letters, numbers, dots, or underscores; no leading digit).", file_context, col), call. = FALSE)
        }
    }

    header <- tryCatch(
        data.table::fread(file_path, nrows = 0L, data.table = FALSE),
    error = function(e) {
        stop(sprintf("Error: Failed to read %s header: %s", file_context, conditionMessage(e)), call. = FALSE)
    }
    )
    header_names <- names(header)

    if (is.null(header_names) || length(header_names) == 0L) {
    stop(sprintf("%s header is empty or invalid.", file_context), call. = FALSE)
    }

    if (!id_col %in% header_names) {
        stop(sprintf("%s missing required ID column: '%s'", file_context, id_col), call. = FALSE)
    }

    missing_cols <- setdiff(value_cols, header_names)
    if (length(missing_cols) > 0L) {
        stop(sprintf("%s missing required column(s): %s", file_context, paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
}


#' @keywords internal
#' @noRd
validate_precomp_fields <- function(precomp_obj, error_message) {
    if (!is.list(precomp_obj) || length(precomp_obj) == 0L) {
        stop(error_message, call. = FALSE)
    }
    nm <- names(precomp_obj)
    if (is.null(nm) || anyNA(nm) || any(nm == "")) {
        stop(error_message, call. = FALSE)
    }
    for (name_i in nm) {
        value_i <- precomp_obj[[name_i]]
        if (is.null(value_i) || length(value_i) == 0L) {
            stop(error_message, call. = FALSE)
        }
        if (is.matrix(value_i) || is.array(value_i) || is.numeric(value_i) || is.integer(value_i)) {
            if (any(!is.finite(value_i))) {
                stop(error_message, call. = FALSE)
            }
        } else if (is.logical(value_i)) {
            if (anyNA(value_i)) {
                stop(error_message, call. = FALSE)
            }
        } else if (is.character(value_i)) {
            if (anyNA(value_i)) {
                stop(error_message, call. = FALSE)
            }
        }
    }
}


#' @keywords internal
#' @noRd
worker_init <- function(gds_path, ids_in_gds, precomp, control, num_ancs,
                        has_covar,
                        method, cond_local, chunk_size, task_file_prefix, task_file_suffix,
                        read_gds_all_then_subset) {
    g <- gdsfmt::openfn.gds(gds_path, readonly = TRUE)

    chrom_node <- gdsfmt::index.gdsn(g, "snp.chromosome")
    pos_node <- gdsfmt::index.gdsn(g, "snp.position")
    id_node <- gdsfmt::index.gdsn(g, "snp.id")
    ref_node <- gdsfmt::index.gdsn(g, "snp.ref")
    alt_node <- gdsfmt::index.gdsn(g, "snp.alt")
    dosage_nodes <- lapply(0:(num_ancs - 1L), function(k) gdsfmt::index.gdsn(g, paste0("dosage/anc", k)))
    hapcount_nodes <- lapply(0:(num_ancs - 1L), function(k) gdsfmt::index.gdsn(g, paste0("hapcount/anc", k)))

    list(
        g = g,
        chrom_node = chrom_node,
        pos_node = pos_node,
        id_node = id_node,
        ref_node = ref_node,
        alt_node = alt_node,
        dosage_nodes = dosage_nodes,
        hapcount_nodes = hapcount_nodes,
        ids_in_gds = ids_in_gds,
        precomp = precomp,
        control = control,
        num_ancs = num_ancs,
        has_covar = has_covar,
        method = method,
        cond_local = cond_local,
        chunk_size = chunk_size,
        task_file_prefix = task_file_prefix,
        task_file_suffix = task_file_suffix,
        read_gds_all_then_subset = read_gds_all_then_subset
    )
}


#' @keywords internal
#' @noRd
run_task <- function(task, worker_env) {
    g <- worker_env$g
    chrom_node <- worker_env$chrom_node
    pos_node <- worker_env$pos_node
    id_node <- worker_env$id_node
    ref_node <- worker_env$ref_node
    alt_node <- worker_env$alt_node
    dosage_nodes <- worker_env$dosage_nodes
    hapcount_nodes <- worker_env$hapcount_nodes
    ids_in_gds <- worker_env$ids_in_gds
    precomp <- worker_env$precomp
    control <- worker_env$control
    num_ancs <- worker_env$num_ancs
    has_covar <- worker_env$has_covar
    method <- worker_env$method
    cond_local <- worker_env$cond_local
    chunk_size <- worker_env$chunk_size
    read_gds_all_then_subset <- worker_env$read_gds_all_then_subset

    output_path <- paste0(worker_env$task_file_prefix, task$id, worker_env$task_file_suffix)

    wrote_header <- FALSE
    local_ancestry_mac_hard_threshold <- as.double(control$local_ancestry_mac_hard_threshold)
    local_ancestry_mac_soft_threshold <- as.double(control$local_ancestry_mac_soft_threshold)
    dos_list <- vector("list", num_ancs)
    hap_list <- vector("list", num_ancs)

    anc_names <- paste0("anc", 0:(num_ancs - 1L))
    anc_names_la <- paste0("anc", 0:(num_ancs - 2L))
    af_col_names <- paste0("AF_", anc_names)
    laprop_col_names <- paste0("LAprop_", anc_names)
    beta_col_names <- paste0("beta_", anc_names)
    se_col_names <- paste0("se_", anc_names)
    pval_col_names <- paste0("pval_", anc_names)
    laeff_col_names <- paste0("LAeff_", anc_names_la)
    lase_col_names <- paste0("LAse_", anc_names_la)
    lapval_col_names <- paste0("LApval_", anc_names_la)
    fallback_used_col_name <- "fallback_used"
    model_col_names <- c(beta_col_names, se_col_names, pval_col_names)
    if (cond_local) {
        model_col_names <- c(model_col_names, laeff_col_names, lase_col_names, lapval_col_names)
    }
    model_col_names <- c(model_col_names, fallback_used_col_name)

    n_sumstats_task <- length(task$sumstats_gds_id)
    sumstats_ptr_start <- 1L
    sumstats_ptr_end <- 0L

    start <- task$start
    idx_len <- chunk_size

    # Keep singleton chunks (idx_len == 1) as 2D matrices to avoid
    # "incorrect number of dimensions" when indexing with [, k].
    ensure_matrix_cols <- function(x, ncol_expected, label = "unnamed") {
        if (is.null(x)) {
            stop(sprintf(
                "run_task: NULL encountered for '%s' while shaping matrix columns (task=%d, chunk_start=%d, chunk_end=%d).",
                label, task$id, start, min(task$end, start + idx_len - 1L)
            ), call. = FALSE)
        }
        if (is.matrix(x)) return(x)
        if (is.null(dim(x))) {
            return(matrix(x, ncol = ncol_expected))
        }
        as.matrix(x)
    }

    ensure_matrix_shape <- function(x, nrow_expected, ncol_expected, label = "unnamed") {
        if (is.null(x)) {
            stop(sprintf(
                "run_task: NULL encountered for '%s' while shaping matrix (expected %d x %d; task=%d, chunk_start=%d, chunk_end=%d).",
                label, nrow_expected, ncol_expected, task$id, start, min(task$end, start + idx_len - 1L)
            ), call. = FALSE)
        }
        if (is.matrix(x)) return(x)
        if (is.null(dim(x))) {
            return(matrix(x, nrow = nrow_expected, ncol = ncol_expected))
        }
        mat <- as.matrix(x)
        if (ncol(mat) != ncol_expected || nrow(mat) != nrow_expected) {
            matrix(as.vector(mat), nrow = nrow_expected, ncol = ncol_expected)
        } else {
            mat
        }
    }

    while (start <= task$end) {
        end <- start + chunk_size - 1L
        if (end > task$end) {
            end <- task$end
            idx_len <- end - start + 1L
        }
        gds_ids <- start:end

        # Obtain sumstats for this chunk
        has_sumstats <- rep.int(FALSE, idx_len)
        sumstats_beta_chunk <- rep(NA_real_, idx_len)
        sumstats_se_chunk <- rep(NA_real_, idx_len)
        if (n_sumstats_task > 0L && sumstats_ptr_start <= n_sumstats_task) {
            while (sumstats_ptr_start <= n_sumstats_task && task$sumstats_gds_id[sumstats_ptr_start] < start) {
                sumstats_ptr_start <- sumstats_ptr_start + 1L
            }
            sumstats_ptr_end <- sumstats_ptr_start - 1L
            while ((sumstats_ptr_end + 1L) <= n_sumstats_task && task$sumstats_gds_id[sumstats_ptr_end + 1L] <= end) {
                sumstats_ptr_end <- sumstats_ptr_end + 1L
            }

            if (sumstats_ptr_end >= sumstats_ptr_start) {
                sum_idx <- sumstats_ptr_start:sumstats_ptr_end
                rel_local <- task$sumstats_gds_id[sum_idx] - start + 1L
                has_sumstats[rel_local] <- TRUE
                sumstats_beta_chunk[rel_local] <- task$sumstats_beta[sum_idx]
                sumstats_se_chunk[rel_local] <- task$sumstats_se[sum_idx]
                sumstats_ptr_start <- sumstats_ptr_end + 1L
            }
        }

        # Read variant information
        chroms <- gdsfmt::read.gdsn(chrom_node, start = start, count = idx_len)
        pos <- gdsfmt::read.gdsn(pos_node, start = start, count = idx_len)
        ids <- gdsfmt::read.gdsn(id_node, start = start, count = idx_len)
        refs <- gdsfmt::read.gdsn(ref_node, start = start, count = idx_len)
        alts <- gdsfmt::read.gdsn(alt_node, start = start, count = idx_len)

        # Read dosage and hapcount data
        if (read_gds_all_then_subset) {
            for (k in seq_len(num_ancs)) {
                    dmat <- gdsfmt::read.gdsn(dosage_nodes[[k]], start = c(1L, start), count = c(-1L, idx_len))
                    dmat <- ensure_matrix_cols(dmat, idx_len, sprintf("dosage read.gdsn anc%d", k - 1L))
                    dmat <- dmat[ids_in_gds, , drop = FALSE]
                    dos_list[[k]] <- dmat
                    hmat <- gdsfmt::read.gdsn(hapcount_nodes[[k]], start = c(1L, start), count = c(-1L, idx_len))
                    hmat <- ensure_matrix_cols(hmat, idx_len, sprintf("hapcount read.gdsn anc%d", k - 1L))
                    hmat <- hmat[ids_in_gds, , drop = FALSE]
                    hap_list[[k]] <- hmat
            }
        } else {
            for (k in seq_len(num_ancs)) {
                dmat <- gdsfmt::readex.gdsn(dosage_nodes[[k]], sel = list(ids_in_gds, gds_ids))
                dmat <- ensure_matrix_shape(dmat, length(ids_in_gds), idx_len, sprintf("dosage readex.gdsn anc%d", k - 1L))
                dos_list[[k]] <- dmat
                hmat <- gdsfmt::readex.gdsn(hapcount_nodes[[k]], sel = list(ids_in_gds, gds_ids))
                hmat <- ensure_matrix_shape(hmat, length(ids_in_gds), idx_len, sprintf("hapcount readex.gdsn anc%d", k - 1L))
                hap_list[[k]] <- hmat
            }
        }

        main_n <- nrow(dos_list[[1L]])
        main_hap_n <- 2 * main_n
        dosage_counts <- vapply(dos_list, matrixStats::colSums2,
                        FUN.VALUE = numeric(idx_len))
        hap_counts <- vapply(hap_list, matrixStats::colSums2,
                     FUN.VALUE = numeric(idx_len))
        dosage_counts <- ensure_matrix_shape(dosage_counts, idx_len, num_ancs)
        hap_counts <- ensure_matrix_shape(hap_counts, idx_len, num_ancs)
        total_dosage_counts <- rowSums(dosage_counts)

        chunk_dt <- data.table::data.table(
            CHROM = chroms,
            POS = pos,
            ID = ids,
            REF = refs,
            ALT = alts,
            main_N = main_n,
            AF = total_dosage_counts / main_hap_n,
            joint_pval = NA_real_,
            has_sumstats = has_sumstats
        )
        for (k in seq_len(num_ancs)) {
            hapk <- hap_counts[, k]
            afk <- dosage_counts[, k] / hapk
            afk[hapk <= 0] <- NA_real_
            data.table::set(chunk_dt, j = af_col_names[k], value = afk)
            data.table::set(chunk_dt, j = laprop_col_names[k], value = hapk / main_hap_n)
        }
        data.table::setcolorder(chunk_dt, c("CHROM", "POS", "ID", "REF", "ALT", "main_N", "AF", "joint_pval", "has_sumstats", af_col_names, laprop_col_names))
        model_init_values <- setNames(as.list(rep(NA_real_, length(model_col_names))), model_col_names)
        model_init_values[[fallback_used_col_name]] <- NA
        for (col_name in model_col_names) {
            data.table::set(chunk_dt, j = col_name, value = model_init_values[[col_name]])
        }

        min_dosage_count <- matrixStats::rowMins(dosage_counts)
        hard_pass <- min_dosage_count > local_ancestry_mac_hard_threshold
        soft_pass <- min_dosage_count > local_ancestry_mac_soft_threshold

        idx_softpass <- which(hard_pass & has_sumstats & soft_pass)
        idx_softfail <- which(hard_pass & has_sumstats & !soft_pass)
        idx_nosum <- which(hard_pass & !has_sumstats)

        if (length(idx_softpass) > 0L) {
            res_softpass <- if (method == "linear") {
                fill_chunk_with_sumstats_softpass_linear(
                    dos_list = dos_list,
                    hap_list = hap_list,
                    idx = idx_softpass,
                    sumstats_beta = sumstats_beta_chunk[idx_softpass],
                    sumstats_se = sumstats_se_chunk[idx_softpass],
                    precomp = precomp,
                    control = control,
                    has_covar = has_covar,
                    cond_local = cond_local,
                    num_ancs = num_ancs
                )
            } else {
                fill_chunk_with_sumstats_softpass_logistic(
                    dos_list = dos_list,
                    hap_list = hap_list,
                    idx = idx_softpass,
                    sumstats_beta = sumstats_beta_chunk[idx_softpass],
                    sumstats_se = sumstats_se_chunk[idx_softpass],
                    precomp = precomp,
                    control = control,
                    cond_local = cond_local,
                    num_ancs = num_ancs
                )
            }
            z_mat <- ensure_matrix_shape(res_softpass$z, length(idx_softpass), num_ancs, "res_softpass$z")
            res_softpass$beta <- ensure_matrix_shape(res_softpass$beta, length(idx_softpass), num_ancs, "res_softpass$beta")
            res_softpass$se <- ensure_matrix_shape(res_softpass$se, length(idx_softpass), num_ancs, "res_softpass$se")
            pval_mat <- 2.0 * stats::pnorm(abs(z_mat), lower.tail = FALSE)
            if (cond_local) {
                laz_mat <- ensure_matrix_shape(res_softpass$laz, length(idx_softpass), num_ancs - 1L, "res_softpass$laz")
                res_softpass$laeff <- ensure_matrix_shape(res_softpass$laeff, length(idx_softpass), num_ancs - 1L, "res_softpass$laeff")
                res_softpass$lase <- ensure_matrix_shape(res_softpass$lase, length(idx_softpass), num_ancs - 1L, "res_softpass$lase")
                lapval_mat <- 2.0 * stats::pnorm(abs(laz_mat), lower.tail = FALSE)
            }
            wald_vec <- res_softpass$wald
            valid_wald <- is.finite(wald_vec) & (wald_vec >= 0)
            joint_pval_vec <- stats::pchisq(wald_vec, df = num_ancs, lower.tail = FALSE)
            joint_pval_vec[!valid_wald] <- NA_real_

            data.table::set(chunk_dt, i = idx_softpass, j = "joint_pval", value = joint_pval_vec)
            data.table::set(chunk_dt, i = idx_softpass, j = fallback_used_col_name, value = as.logical(res_softpass$used_qr_fallback))
            for (k in seq_len(num_ancs)) {
                data.table::set(chunk_dt, i = idx_softpass, j = beta_col_names[k], value = res_softpass$beta[, k])
                data.table::set(chunk_dt, i = idx_softpass, j = se_col_names[k], value = res_softpass$se[, k])
                data.table::set(chunk_dt, i = idx_softpass, j = pval_col_names[k], value = pval_mat[, k])
            }
            if (cond_local) {
                for (k in seq_len(num_ancs - 1L)) {
                    data.table::set(chunk_dt, i = idx_softpass, j = laeff_col_names[k], value = res_softpass$laeff[, k])
                    data.table::set(chunk_dt, i = idx_softpass, j = lase_col_names[k], value = res_softpass$lase[, k])
                    data.table::set(chunk_dt, i = idx_softpass, j = lapval_col_names[k], value = lapval_mat[, k])
                }
            }
        }

        if (length(idx_softfail) > 0L) {
            res_softfail <- if (method == "linear") {
                fill_chunk_with_sumstats_softfail_linear(
                    dos_list = dos_list,
                    hap_list = hap_list,
                    idx = idx_softfail,
                    sumstats_beta = sumstats_beta_chunk[idx_softfail],
                    sumstats_se = sumstats_se_chunk[idx_softfail],
                    precomp = precomp,
                    control = control,
                    has_covar = has_covar,
                    cond_local = cond_local,
                    num_ancs = num_ancs
                )
            } else {
                fill_chunk_with_sumstats_softfail_logistic(
                    dos_list = dos_list,
                    hap_list = hap_list,
                    idx = idx_softfail,
                    sumstats_beta = sumstats_beta_chunk[idx_softfail],
                    sumstats_se = sumstats_se_chunk[idx_softfail],
                    precomp = precomp,
                    control = control,
                    cond_local = cond_local,
                    num_ancs = num_ancs
                )
            }
            z_mat <- ensure_matrix_shape(res_softfail$z, length(idx_softfail), num_ancs, "res_softfail$z")
            res_softfail$beta <- ensure_matrix_shape(res_softfail$beta, length(idx_softfail), num_ancs, "res_softfail$beta")
            res_softfail$se <- ensure_matrix_shape(res_softfail$se, length(idx_softfail), num_ancs, "res_softfail$se")
            pval_mat <- 2.0 * stats::pnorm(abs(z_mat), lower.tail = FALSE)
            if (cond_local) {
                laz_mat <- ensure_matrix_shape(res_softfail$laz, length(idx_softfail), num_ancs - 1L, "res_softfail$laz")
                res_softfail$laeff <- ensure_matrix_shape(res_softfail$laeff, length(idx_softfail), num_ancs - 1L, "res_softfail$laeff")
                res_softfail$lase <- ensure_matrix_shape(res_softfail$lase, length(idx_softfail), num_ancs - 1L, "res_softfail$lase")
                lapval_mat <- 2.0 * stats::pnorm(abs(laz_mat), lower.tail = FALSE)
            }
            wald_vec <- res_softfail$wald
            valid_wald <- is.finite(wald_vec) & (wald_vec >= 0)
            joint_pval_vec <- stats::pchisq(wald_vec, df = num_ancs, lower.tail = FALSE)
            joint_pval_vec[!valid_wald] <- NA_real_

            data.table::set(chunk_dt, i = idx_softfail, j = "joint_pval", value = joint_pval_vec)
            data.table::set(chunk_dt, i = idx_softfail, j = fallback_used_col_name, value = as.logical(res_softfail$used_qr_fallback))
            for (k in seq_len(num_ancs)) {
                data.table::set(chunk_dt, i = idx_softfail, j = beta_col_names[k], value = res_softfail$beta[, k])
                data.table::set(chunk_dt, i = idx_softfail, j = se_col_names[k], value = res_softfail$se[, k])
                data.table::set(chunk_dt, i = idx_softfail, j = pval_col_names[k], value = pval_mat[, k])
            }
            if (cond_local) {
                for (k in seq_len(num_ancs - 1L)) {
                    data.table::set(chunk_dt, i = idx_softfail, j = laeff_col_names[k], value = res_softfail$laeff[, k])
                    data.table::set(chunk_dt, i = idx_softfail, j = lase_col_names[k], value = res_softfail$lase[, k])
                    data.table::set(chunk_dt, i = idx_softfail, j = lapval_col_names[k], value = lapval_mat[, k])
                }
            }
        }

        if (length(idx_nosum) > 0L) {
            res_nosum <- if (method == "linear") {
                fill_chunk_without_sumstats_linear(
                    dos_list = dos_list,
                    hap_list = hap_list,
                    idx = idx_nosum,
                    precomp = precomp,
                    control = control,
                    has_covar = has_covar,
                    cond_local = cond_local,
                    num_ancs = num_ancs
                )
            } else {
                fill_chunk_without_sumstats_logistic(
                    dos_list = dos_list,
                    hap_list = hap_list,
                    idx = idx_nosum,
                    precomp = precomp,
                    control = control,
                    cond_local = cond_local,
                    num_ancs = num_ancs
                )
            }

            if (method == "linear") {
                t_mat <- ensure_matrix_shape(res_nosum$t, length(idx_nosum), num_ancs, "res_nosum$t")
                res_nosum$beta <- ensure_matrix_shape(res_nosum$beta, length(idx_nosum), num_ancs, "res_nosum$beta")
                res_nosum$se <- ensure_matrix_shape(res_nosum$se, length(idx_nosum), num_ancs, "res_nosum$se")
                df_resid_vec <- as.numeric(res_nosum$df_resid) # no need for conversion?
                valid_df <- is.finite(df_resid_vec) & (df_resid_vec > 0)
                pval_mat <- 2.0 * stats::pt(abs(t_mat), df = df_resid_vec, lower.tail = FALSE)
                pval_mat[!valid_df, ] <- NA_real_
                if (cond_local) {
                    # C++ returns local-ancestry t statistics as `lat`.
                    # Keep a fallback to `la_t` for backward compatibility with older builds.
                    la_t_src <- if (!is.null(res_nosum$lat)) res_nosum$lat else res_nosum$la_t
                    la_t_label <- if (!is.null(res_nosum$lat)) "res_nosum$lat" else "res_nosum$la_t"
                    la_t_mat <- ensure_matrix_shape(la_t_src, length(idx_nosum), num_ancs - 1L, la_t_label)
                    res_nosum$laeff <- ensure_matrix_shape(res_nosum$laeff, length(idx_nosum), num_ancs - 1L, "res_nosum$laeff")
                    res_nosum$lase <- ensure_matrix_shape(res_nosum$lase, length(idx_nosum), num_ancs - 1L, "res_nosum$lase")
                    lapval_mat <- 2.0 * stats::pt(abs(la_t_mat), df = df_resid_vec, lower.tail = FALSE)
                    lapval_mat[!valid_df, ] <- NA_real_
                }
            } else {
                z_mat <- ensure_matrix_shape(res_nosum$z, length(idx_nosum), num_ancs, "res_nosum$z")
                res_nosum$beta <- ensure_matrix_shape(res_nosum$beta, length(idx_nosum), num_ancs, "res_nosum$beta")
                res_nosum$se <- ensure_matrix_shape(res_nosum$se, length(idx_nosum), num_ancs, "res_nosum$se")
                pval_mat <- 2.0 * stats::pnorm(abs(z_mat), lower.tail = FALSE)
                if (cond_local) {
                    laz_mat <- ensure_matrix_shape(res_nosum$laz, length(idx_nosum), num_ancs - 1L, "res_nosum$laz")
                    res_nosum$laeff <- ensure_matrix_shape(res_nosum$laeff, length(idx_nosum), num_ancs - 1L, "res_nosum$laeff")
                    res_nosum$lase <- ensure_matrix_shape(res_nosum$lase, length(idx_nosum), num_ancs - 1L, "res_nosum$lase")
                    lapval_mat <- 2.0 * stats::pnorm(abs(laz_mat), lower.tail = FALSE)
                }                
            }

            wald_vec <- res_nosum$wald
            valid_wald <- is.finite(wald_vec) & (wald_vec >= 0)
            joint_pval_vec <- stats::pchisq(wald_vec, df = num_ancs, lower.tail = FALSE)
            joint_pval_vec[!valid_wald] <- NA_real_

            data.table::set(chunk_dt, i = idx_nosum, j = "joint_pval", value = joint_pval_vec)
            data.table::set(chunk_dt, i = idx_nosum, j = fallback_used_col_name, value = as.logical(res_nosum$used_qr_fallback))
            for (k in seq_len(num_ancs)) {
                data.table::set(chunk_dt, i = idx_nosum, j = beta_col_names[k], value = res_nosum$beta[, k])
                data.table::set(chunk_dt, i = idx_nosum, j = se_col_names[k], value = res_nosum$se[, k])
                data.table::set(chunk_dt, i = idx_nosum, j = pval_col_names[k], value = pval_mat[, k])
            }
            if (cond_local) {
                for (k in seq_len(num_ancs - 1L)) {
                    data.table::set(chunk_dt, i = idx_nosum, j = laeff_col_names[k], value = res_nosum$laeff[, k])
                    data.table::set(chunk_dt, i = idx_nosum, j = lase_col_names[k], value = res_nosum$lase[, k])
                    data.table::set(chunk_dt, i = idx_nosum, j = lapval_col_names[k], value = lapval_mat[, k])
                }
            }
        }

        data.table::fwrite(
            chunk_dt,
            file = output_path,
            sep = "\t",
            quote = FALSE,
            na = "NA",
            append = wrote_header,
            col.names = !wrote_header
        )
        wrote_header <- TRUE
        start <- end + 1L
    }
}