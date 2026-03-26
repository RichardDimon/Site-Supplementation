## Randomised supplementation allele and kinship curves 
## created by R. Dimon March 2026 with helper functions using GPT-5.3
## Focal target site + random supplementation from all available sites in the dataset
##
## Outputs:
##   1) allele summary table
##   2) allele replicate table
##   3) kinship summary table
##   4) kinship replicate table
##   5) combined 3-panel plot
##
## Notes:
## - This script assumes your dms object already contains: dms$gt and dms$meta$analyses
## - A full pairwise kinship matrix is calculated once at the start
## - At each step, kinship metrics (mean, min/max, proportion > threshold)
##   are calculated for the CURRENT pooled sample (target + added individuals)
## - Ribbons are min-max across replicates (but can be replaced for quantiles if needed)
## - The final plot shows mean kinship only
## - Proportion of related pairs is calculated and saved to csv outputs
## - Supplementation stops only when BOTH conditions are met:
##     1) min common allele capture across reps > 0.90
##     2) max mean kinship across reps < 0.0442, an unrelated kinship threshold based on Manichaikul et al. (2010)
## - There are 2 supplementation strategies: 
##      1) Sampling up to 1 individual per site (minimise total individuals to visit)
##      2) Sampling up to 5 individuals per site (minimise total sites to visit)



# 1. MAIN SETTINGS TO EDIT ==============================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# species working directory (containing 0_setup_variables.xlsx)
setwd("C:/Users/dimon/rrspecies/Acacia_decora/")

# focal target site (add exact site name below)
target_site <- "23"

# number of random replicate runs per sample size
max_steps <- 1000

# sample size range
# use "max" to run from all focal individuals up to all available individuals
N_t_vec <- "max"
# N_t_vec <- 6:60

# site-batch strategy: max number sampled per site before moving to another site
site_batch_size <- 5

# unrelated kinship threshold
kin_unrelated_thresh <- 0.0442 #value based on Manichaikul et al. (2010)

# minor allele frequency threshold for common vs rare
threshold_maf <- 0.05

# whether to include loci with missing data
IncludeNA <- "withNA"   # "withNA" or "noNA"




# 2. PACKAGES ===========================================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---


library(openxlsx)
library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(devtools)




# 3. LOAD SETUP VARIABLES ===============================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---


setup_variables <- read.xlsx("0_setup_variables.xlsx", colNames = TRUE)

maindir <- setup_variables[1, 2]
species <- setup_variables[2, 2]
dataset <- setup_variables[3, 2]
raw_meta_path <- setup_variables[4, 2]
species_col_name <- setup_variables[5, 2]
site_col_name <- setup_variables[6, 2]

locus_miss <- as.numeric(setup_variables[12, 2])
sample_miss <- as.numeric(setup_variables[13, 2])
maf_val <- as.numeric(setup_variables[14, 2])
clonal_threshold <- as.numeric(setup_variables[15, 2])
custom_meta <- setup_variables[16, 2]

RandRbase <- ""

setwd(maindir)



# 4. LOAD DATA ==========================================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---


dms <- readRDS(
  paste0(
    species, "/outputs_", site_col_name, "_", species_col_name, "/r_files/dms.RData"
  )
)

gt_sw_comp <- dms$gt
gt_ref <- 2 - gt_sw_comp   # reference allele count matrix



# 5. GET KINSHIP MATRIX ==================================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source_url("https://github.com/RichardDimon/R_Scripts/blob/main/individual_kinship_by_pop_sp.R?raw=TRUE")

kin <- individual_kinship_by_pop(
  dms,
  RandRbase,
  species,
  dataset,
  dms$meta$analyses[, species_col_name],
  maf = 0.05,
  mis = locus_miss,
  as_bigmat = TRUE
)

kin <- as.matrix(kin)
kin[is.na(kin)] <- 0

gt_names <- rownames(gt_sw_comp)

if (!all(gt_names %in% rownames(kin))) {
  missing_ids <- setdiff(gt_names, rownames(kin))
  stop(
    "Kinship matrix is missing samples found in gt_sw_comp. Example: ",
    paste(head(missing_ids, 5), collapse = ", ")
  )
}

K_full <- kin[gt_names, gt_names, drop = FALSE]

cat("\nK_full overall summary:\n")
print(summary(as.numeric(K_full[upper.tri(K_full)])))
cat("K_full min/max:", min(K_full, na.rm = TRUE), max(K_full, na.rm = TRUE), "\n")



# 6. ALLELE FUNCTIONS ===================================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

if (IncludeNA == "withNA") {
  source("https://raw.githubusercontent.com/RichardDimon/R_Scripts/main/get_minor_allele_frequenciesNARM.r")
  cat("Running with loci containing missing data. Total loci =", ncol(gt_sw_comp), "\n")
}

if (IncludeNA == "noNA") {
  gt_sw_comp <- gt_sw_comp[, which(colSums(gt_sw_comp) > 0 & colSums(gt_sw_comp) < nrow(gt_sw_comp) * 2), drop = FALSE]
  gt_ref <- 2 - gt_sw_comp
  source("https://raw.githubusercontent.com/RichardDimon/R_Scripts/main/get_minor_allele_frequencies.r")
  cat("Running with loci with no missing data. Total loci =", ncol(gt_sw_comp), "\n")
}

sw_maf <- get_minor_allele_frequencies(gt_sw_comp)
i_sw_common <- which(sw_maf > threshold_maf)
i_sw_rare <- which(sw_maf <= threshold_maf & sw_maf > 0)

cat("Common loci:", length(i_sw_common), "\n")
cat("Rare loci:", length(i_sw_rare), "\n")


# 7. ALIGN META TO GENOTYPE ORDER =======================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

meta_df <- as.data.frame(dms$meta$analyses, stringsAsFactors = FALSE)
meta_df$sample <- trimws(as.character(meta_df$sample))
meta_df[[site_col_name]] <- trimws(as.character(meta_df[[site_col_name]]))

site_vec <- meta_df[[site_col_name]][match(rownames(gt_sw_comp), meta_df$sample)]

if (anyNA(site_vec)) {
  warning(sum(is.na(site_vec)), " samples in gt_sw_comp did not match meta data.")
}

if (!(target_site %in% site_vec)) {
  stop("Target site not found in aligned metadata: ", target_site)
}


# 8. HELPER FUNCTIONS ===================================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# summarise a replicate x Nt matrix
summarise_matrix <- function(mat, metric_name, phase_name, strategy_name) {
  if (is.null(mat) || ncol(mat) == 0) return(NULL)
  
  safe_stat <- function(x, fun) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    fun(x)
  }
  
  out <- data.frame(
    N_t = as.numeric(colnames(mat)),
    metric = metric_name,
    mean = apply(mat, 2, safe_stat, fun = mean),
    min = apply(mat, 2, safe_stat, fun = min),
    max = apply(mat, 2, safe_stat, fun = max),
    stringsAsFactors = FALSE
  )
  
  out$phase <- phase_name
  out$strategy <- strategy_name
  out
}


# convert matrix to long format for export
matrix_to_long <- function(mat, metric_name, phase_name, strategy_name) {
  if (is.null(mat) || ncol(mat) == 0) return(NULL)
  
  out <- reshape2::melt(mat)
  colnames(out) <- c("replicate", "N_t_col", "value")
  out$N_t <- as.numeric(colnames(mat)[out$N_t_col])
  out$metric <- metric_name
  out$phase <- phase_name
  out$strategy <- strategy_name
  out$N_t_col <- NULL
  out
}


# pairwise kinship summary for the selected pooled sample
calc_kinship_stats <- function(sel_idx, Kmat, threshold) {
  if (is.null(Kmat) || length(sel_idx) < 2) {
    return(c(prop_above = NA_real_, mean_coeff = NA_real_, min_coeff = NA_real_, max_coeff = NA_real_))
  }
  
  kvals <- Kmat[sel_idx, sel_idx, drop = FALSE]
  kvals <- kvals[upper.tri(kvals)]
  kvals <- kvals[!is.na(kvals)]
  
  if (length(kvals) == 0) {
    return(c(prop_above = NA_real_, mean_coeff = NA_real_, min_coeff = NA_real_, max_coeff = NA_real_))
  }
  
  c(
    prop_above = mean(kvals > threshold),
    mean_coeff = mean(kvals),
    min_coeff = min(kvals),
    max_coeff = max(kvals)
  )
}


# allele capture for a selected set of individuals
calc_allele_capture <- function(sel_logical, gt_alt, gt_ref, i_common, i_rare) {
  alt_counts <- colSums(gt_alt[sel_logical, , drop = FALSE], na.rm = TRUE)
  ref_counts <- colSums(gt_ref[sel_logical, , drop = FALSE], na.rm = TRUE)
  minor_counts <- pmin(alt_counts, ref_counts)
  captured <- minor_counts > 0
  
  common_prop <- if (length(i_common) > 0) sum(captured[i_common]) / length(i_common) else NA_real_
  rare_prop <- if (length(i_rare) > 0) sum(captured[i_rare]) / length(i_rare) else NA_real_
  
  c(common_prop = common_prop, rare_prop = rare_prop)
}


# strategy 1: add one individual per site in a progressive way
sample_progressive_sites <- function(meta_sites, n_needed) {
  out <- rep(0, length(meta_sites))
  site_levels <- unique(meta_sites)
  
  while (sum(out) < n_needed && any(out == 0)) {
    for (s in site_levels) {
      idx <- which(meta_sites == s & out == 0)
      if (length(idx) > 0 && sum(out) < n_needed) {
        out[sample(idx, 1)] <- 1
      }
    }
  }
  
  out
}


# strategy 2: take several from one site, then move on
sample_site_batches <- function(meta_sites, n_needed, batch_size = 5) {
  out <- rep(0, length(meta_sites))
  available_idx <- seq_along(meta_sites)
  still_needed <- n_needed
  
  while (still_needed > 0 && length(available_idx) > 0) {
    available_sites <- unique(meta_sites[available_idx])
    chosen_site <- sample(available_sites, 1)
    
    site_idx <- available_idx[meta_sites[available_idx] == chosen_site]
    n_take <- min(batch_size, length(site_idx), still_needed)
    
    take_idx <- sample(site_idx, n_take, replace = FALSE)
    out[take_idx] <- 1
    
    still_needed <- still_needed - n_take
    available_idx <- which(out == 0)
  }
  
  out
}


# run one whole strategy across a vector of total sample sizes
run_strategy <- function(
    strategy_name,
    N_t_values,
    forced_idx,
    available_idx,
    site_vec,
    gt_sw_comp,
    gt_ref,
    i_sw_common,
    i_sw_rare,
    K_full,
    max_steps,
    kin_unrelated_thresh,
    sampler,
    batch_size = NULL
) {
  out_common <- matrix(NA_real_, nrow = max_steps, ncol = length(N_t_values))
  out_rare <- matrix(NA_real_, nrow = max_steps, ncol = length(N_t_values))
  out_kin_prop <- matrix(NA_real_, nrow = max_steps, ncol = length(N_t_values))
  out_kin_mean <- matrix(NA_real_, nrow = max_steps, ncol = length(N_t_values))
  
  colnames(out_common) <- N_t_values
  colnames(out_rare) <- N_t_values
  colnames(out_kin_prop) <- N_t_values
  colnames(out_kin_mean) <- N_t_values
  
  for (i in seq_along(N_t_values)) {
    iNt <- N_t_values[i]
    cat("\n[", strategy_name, "] N_t = ", iNt, " (", i, "/", length(N_t_values), ")\n", sep = "")
    
    pb <- txtProgressBar(min = 0, max = max_steps, style = 3)
    
    for (j in seq_len(max_steps)) {
      selected <- rep(FALSE, nrow(gt_sw_comp))
      selected[forced_idx] <- TRUE
      
      n_remaining <- iNt - sum(selected)
      
      if (n_remaining > 0) {
        avail_sites <- site_vec[available_idx]
        
        add_vec <- if (is.null(batch_size)) {
          sampler(avail_sites, n_remaining)
        } else {
          sampler(avail_sites, n_remaining, batch_size = batch_size)
        }
        
        selected[available_idx] <- as.logical(add_vec)
      }
      
      allele_stats <- calc_allele_capture(
        sel_logical = selected,
        gt_alt = gt_sw_comp,
        gt_ref = gt_ref,
        i_common = i_sw_common,
        i_rare = i_sw_rare
      )
      
      kin_stats <- calc_kinship_stats(
        sel_idx = which(selected),
        Kmat = K_full,
        threshold = kin_unrelated_thresh
      )
      
      out_common[j, i] <- allele_stats["common_prop"]
      out_rare[j, i] <- allele_stats["rare_prop"]
      out_kin_prop[j, i] <- kin_stats["prop_above"]
      out_kin_mean[j, i] <- kin_stats["mean_coeff"]
      
      if (j %% 50 == 0 || j == max_steps) setTxtProgressBar(pb, j)
    }
    
    close(pb)
    
    # stop once BOTH conditions are met:
    # 1) all replicates exceed 90% common allele capture
    # 2) all replicates are below the unrelated threshold for mean kinship
    cond_common <- min(out_common[, i], na.rm = TRUE) > 0.90
    cond_kin <- max(out_kin_mean[, i], na.rm = TRUE) < kin_unrelated_thresh
    
    if (cond_common && cond_kin) {
      message(
        "Stopping early for ", strategy_name, " at N_t = ", iNt,
        " because both criteria were met: ",
        "min common allele capture > 0.90 and ",
        "max mean kinship < ", kin_unrelated_thresh, "."
      )
      
      keep <- seq_len(i)
      out_common <- out_common[, keep, drop = FALSE]
      out_rare <- out_rare[, keep, drop = FALSE]
      out_kin_prop <- out_kin_prop[, keep, drop = FALSE]
      out_kin_mean <- out_kin_mean[, keep, drop = FALSE]
      break
    }
  }
  
  list(
    strategy = strategy_name,
    common = out_common,
    rare = out_rare,
    kin_prop = out_kin_prop,
    kin_mean = out_kin_mean
  )
}


# 9. MAIN ANALYSIS ======================================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

run_randomised_supplementation <- function(
    target_site,
    site_vec,
    gt_sw_comp,
    gt_ref,
    i_sw_common,
    i_sw_rare,
    K_full,
    N_t_vec,
    max_steps,
    kin_unrelated_thresh,
    site_batch_size,
    species,
    site_col_name,
    species_col_name
) {
  forced_idx <- which(site_vec == target_site)
  
  if (length(forced_idx) == 0) {
    stop("No samples found for target site: ", target_site)
  }
  
  all_idx <- seq_len(nrow(gt_sw_comp))
  available_idx <- setdiff(all_idx, forced_idx)
  
  n_forced <- length(forced_idx)
  start_n <- n_forced + 1
  max_n <- n_forced + length(available_idx)
  
  if (start_n > max_n) {
    stop("No additional individuals available beyond the focal site.")
  }
  
  if (length(N_t_vec) == 1 && N_t_vec[1] == "max") {
    N_t_run <- start_n:max_n
  } else {
    N_t_run <- sort(unique(as.integer(N_t_vec)))
    N_t_run <- N_t_run[N_t_run >= start_n & N_t_run <= max_n]
  }
  
  if (length(N_t_run) == 0) {
    stop("No valid values remain in N_t_vec after filtering.")
  }
  
  cat("\n--------------------------------------------------\n")
  cat("Target site:", target_site, "\n")
  cat("Focal individuals forced in:", n_forced, "\n")
  cat("N_t range:", min(N_t_run), "to", max(N_t_run), "\n")
  cat("Replicates per N_t:", max_steps, "\n")
  cat("--------------------------------------------------\n")
  
  # ---- target site alone ----
  N_t_focal <- seq_len(n_forced)
  
  focal_common <- matrix(NA_real_, nrow = max_steps, ncol = length(N_t_focal))
  focal_rare <- matrix(NA_real_, nrow = max_steps, ncol = length(N_t_focal))
  focal_kin_prop <- matrix(NA_real_, nrow = max_steps, ncol = length(N_t_focal))
  focal_kin_mean <- matrix(NA_real_, nrow = max_steps, ncol = length(N_t_focal))
  
  colnames(focal_common) <- N_t_focal
  colnames(focal_rare) <- N_t_focal
  colnames(focal_kin_prop) <- N_t_focal
  colnames(focal_kin_mean) <- N_t_focal
  
  cat("\nRunning target-site accumulation...\n")
  
  for (i in seq_along(N_t_focal)) {
    iNt <- N_t_focal[i]
    cat("\n[Target site] N_t = ", iNt, " (", i, "/", length(N_t_focal), ")\n", sep = "")
    
    pb <- txtProgressBar(min = 0, max = max_steps, style = 3)
    
    for (j in seq_len(max_steps)) {
      sel_idx <- sample(forced_idx, size = iNt, replace = FALSE)
      selected <- rep(FALSE, nrow(gt_sw_comp))
      selected[sel_idx] <- TRUE
      
      allele_stats <- calc_allele_capture(
        sel_logical = selected,
        gt_alt = gt_sw_comp,
        gt_ref = gt_ref,
        i_common = i_sw_common,
        i_rare = i_sw_rare
      )
      
      kin_stats <- calc_kinship_stats(
        sel_idx = sel_idx,
        Kmat = K_full,
        threshold = kin_unrelated_thresh
      )
      
      focal_common[j, i] <- allele_stats["common_prop"]
      focal_rare[j, i] <- allele_stats["rare_prop"]
      focal_kin_prop[j, i] <- kin_stats["prop_above"]
      focal_kin_mean[j, i] <- kin_stats["mean_coeff"]
      
      if (j %% 50 == 0 || j == max_steps) setTxtProgressBar(pb, j)
    }
    
    close(pb)
  }
  
  # ---- additional sampling strategies ----
  cat("\nRunning supplementation strategies...\n")
  
  res_progressive <- run_strategy(
    strategy_name = "1 indv. per site",
    N_t_values = N_t_run,
    forced_idx = forced_idx,
    available_idx = available_idx,
    site_vec = site_vec,
    gt_sw_comp = gt_sw_comp,
    gt_ref = gt_ref,
    i_sw_common = i_sw_common,
    i_sw_rare = i_sw_rare,
    K_full = K_full,
    max_steps = max_steps,
    kin_unrelated_thresh = kin_unrelated_thresh,
    sampler = sample_progressive_sites,
    batch_size = NULL
  )
  
  res_sitebatch <- run_strategy(
    strategy_name = "5 indv. per site",
    N_t_values = N_t_run,
    forced_idx = forced_idx,
    available_idx = available_idx,
    site_vec = site_vec,
    gt_sw_comp = gt_sw_comp,
    gt_ref = gt_ref,
    i_sw_common = i_sw_common,
    i_sw_rare = i_sw_rare,
    K_full = K_full,
    max_steps = max_steps,
    kin_unrelated_thresh = kin_unrelated_thresh,
    sampler = sample_site_batches,
    batch_size = site_batch_size
  )
  
  # ---- summaries ----
  curve_summary <- bind_rows(
    summarise_matrix(focal_common, "Common Alleles", "Target site", "Target site"),
    summarise_matrix(focal_rare, "Rare Alleles", "Target site", "Target site"),
    summarise_matrix(res_progressive$common, "Common Alleles", "Additional sampling", res_progressive$strategy),
    summarise_matrix(res_progressive$rare, "Rare Alleles", "Additional sampling", res_progressive$strategy),
    summarise_matrix(res_sitebatch$common, "Common Alleles", "Additional sampling", res_sitebatch$strategy),
    summarise_matrix(res_sitebatch$rare, "Rare Alleles", "Additional sampling", res_sitebatch$strategy),
    summarise_matrix(focal_kin_prop, paste0("Prop. kinship > ", kin_unrelated_thresh), "Target site", "Target site"),
    summarise_matrix(focal_kin_mean, "Mean kinship coefficient", "Target site", "Target site"),
    summarise_matrix(res_progressive$kin_prop, paste0("Prop. kinship > ", kin_unrelated_thresh), "Additional sampling", res_progressive$strategy),
    summarise_matrix(res_progressive$kin_mean, "Mean kinship coefficient", "Additional sampling", res_progressive$strategy),
    summarise_matrix(res_sitebatch$kin_prop, paste0("Prop. kinship > ", kin_unrelated_thresh), "Additional sampling", res_sitebatch$strategy),
    summarise_matrix(res_sitebatch$kin_mean, "Mean kinship coefficient", "Additional sampling", res_sitebatch$strategy)
  )
  
  curve_summary$phase <- factor(curve_summary$phase, levels = c("Target site", "Additional sampling"))
  curve_summary$strategy <- factor(
    curve_summary$strategy,
    levels = c("Target site", "5 indv. per site", "1 indv. per site")
  )
  
  curve_summary$metric <- factor(
    curve_summary$metric,
    levels = c(
      "Common Alleles",
      "Rare Alleles",
      "Mean kinship coefficient",
      paste0("Prop. kinship > ", kin_unrelated_thresh)
    )
  )
  
  replicate_long <- bind_rows(
    matrix_to_long(focal_common, "Common Alleles", "Target site", "Target site"),
    matrix_to_long(focal_rare, "Rare Alleles", "Target site", "Target site"),
    matrix_to_long(focal_kin_prop, paste0("Prop. kinship > ", kin_unrelated_thresh), "Target site", "Target site"),
    matrix_to_long(focal_kin_mean, "Mean kinship coefficient", "Target site", "Target site"),
    matrix_to_long(res_progressive$common, "Common Alleles", "Additional sampling", res_progressive$strategy),
    matrix_to_long(res_progressive$rare, "Rare Alleles", "Additional sampling", res_progressive$strategy),
    matrix_to_long(res_progressive$kin_prop, paste0("Prop. kinship > ", kin_unrelated_thresh), "Additional sampling", res_progressive$strategy),
    matrix_to_long(res_progressive$kin_mean, "Mean kinship coefficient", "Additional sampling", res_progressive$strategy),
    matrix_to_long(res_sitebatch$common, "Common Alleles", "Additional sampling", res_sitebatch$strategy),
    matrix_to_long(res_sitebatch$rare, "Rare Alleles", "Additional sampling", res_sitebatch$strategy),
    matrix_to_long(res_sitebatch$kin_prop, paste0("Prop. kinship > ", kin_unrelated_thresh), "Additional sampling", res_sitebatch$strategy),
    matrix_to_long(res_sitebatch$kin_mean, "Mean kinship coefficient", "Additional sampling", res_sitebatch$strategy)
  )
  
  replicate_long$phase <- factor(replicate_long$phase, levels = c("Target site", "Additional sampling"))
  replicate_long$strategy <- factor(
    replicate_long$strategy,
    levels = c("Target site", "5 indv. per site", "1 indv. per site")
  )
  replicate_long$metric <- factor(
    replicate_long$metric,
    levels = c(
      "Common Alleles",
      "Rare Alleles",
      "Mean kinship coefficient",
      paste0("Prop. kinship > ", kin_unrelated_thresh)
    )
  )
  
  # ---- output folder ----
  out_dir <- paste0(
    species, "/outputs_", site_col_name, "_", species_col_name,
    "/Supplementation/RandSites_AllSites/Site_", target_site
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_prefix <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_\\-]", "_", target_site), "_allsites"))
  
  
  
  
  # ---- plot ----
  
  plot_dat <- subset(
    curve_summary,
    metric %in% c("Common Alleles", "Rare Alleles", "Mean kinship coefficient")
  )
  
  plot_dat$metric <- factor(
    plot_dat$metric,
    levels = c("Common Alleles", "Rare Alleles", "Mean kinship coefficient")
  )
  
  p_summary <- ggplot(
    plot_dat,
    aes(x = N_t, y = mean, colour = strategy, fill = strategy)
  ) +
    geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1) +
    
    # target line for common allele capture
    geom_hline(
      data = subset(plot_dat, metric == "Common Alleles"),
      aes(yintercept = 0.90),
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    
    # target line for mean kinship
    geom_hline(
      data = subset(plot_dat, metric == "Mean kinship coefficient"),
      aes(yintercept = kin_unrelated_thresh),
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    
    facet_grid(metric ~ ., scales = "free_y") +
    scale_x_continuous(
      breaks = seq(
        floor(min(curve_summary$N_t, na.rm = TRUE)),
        ceiling(max(curve_summary$N_t, na.rm = TRUE)),
        by = 1
      )
    ) +
    scale_colour_manual(values = c(
      "Target site" = "black",
      "5 indv. per site" = "orange2",
      "1 indv. per site" = "skyblue2"
    )) +
    scale_fill_manual(values = c(
      "Target site" = "black",
      "5 indv. per site" = "orange2",
      "1 indv. per site" = "skyblue2"
    )) +
    labs(
      title = paste0("Randomised supplementation curves for site '", target_site, "'"),
      subtitle = "Final kinship panel shows mean kinship coefficient only",
      x = "Total number of individuals",
      y = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # ---- save outputs ----
  write.csv(
    subset(curve_summary, metric %in% c("Common Alleles", "Rare Alleles")),
    paste0(out_prefix, "_allele_curve_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    subset(replicate_long, metric %in% c("Common Alleles", "Rare Alleles")),
    paste0(out_prefix, "_allele_curve_replicates.csv"),
    row.names = FALSE
  )
  
  write.csv(
    subset(curve_summary, metric %in% c("Mean kinship coefficient", paste0("Prop. kinship > ", kin_unrelated_thresh))),
    paste0(out_prefix, "_kinship_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    subset(replicate_long, metric %in% c("Mean kinship coefficient", paste0("Prop. kinship > ", kin_unrelated_thresh))),
    paste0(out_prefix, "_kinship_replicates.csv"),
    row.names = FALSE
  )
  
  ggsave(
    filename = paste0(out_prefix, "_summary_curves.png"),
    plot = p_summary,
    width = 8,
    height = 10,
    dpi = 300
  )
  
  invisible(list(
    target_site = target_site,
    forced_idx = forced_idx,
    N_t_run = N_t_run,
    summary = curve_summary,
    replicates = replicate_long,
    plot = p_summary
  ))
}



# 10. RUN ===============================================
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

res_allsites <- run_randomised_supplementation(
  target_site = target_site,
  site_vec = site_vec,
  gt_sw_comp = gt_sw_comp,
  gt_ref = gt_ref,
  i_sw_common = i_sw_common,
  i_sw_rare = i_sw_rare,
  K_full = K_full,
  N_t_vec = N_t_vec,
  max_steps = max_steps,
  kin_unrelated_thresh = kin_unrelated_thresh,
  site_batch_size = site_batch_size,
  species = species,
  site_col_name = site_col_name,
  species_col_name = species_col_name
)
