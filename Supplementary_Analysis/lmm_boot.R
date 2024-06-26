# bootstrap 80% samples to validate LMMs
library(caret)
library(r2glmm)
library(lmerTest)
library(tidyverse)
library(doParallel)

setwd("H:/ABCD/Release4.0/Package_1194636/results/peer_environments/")
load("results.RData")

# boot function ----------------------------------------------------------------
boot_lmm <- function(boot_index, var) {
  
  abcd_year2_behavior <- abcd_year2_behavior[boot_behavior_mat[, boot_index], ]
  abcd_year2_smri <- abcd_year2_smri[boot_smri_mat[, boot_index], ]
  abcd_year2_rsfmri <- abcd_year2_rsfmri[boot_rsfmri_mat[, boot_index], ]
  
  # behavior
  lmm_behavior <- function(x, df) {
    cat("\n", x, "   Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ df[[var]] + sex + interview_age + race_ethnicity +
                  income_parent + edu_parent + (1|site_id_l/rel_family_id), data = df)
    fit_summary <- summary(fit)
    # 3 df; 4 t-value; 5 p-value
    t_value <- fit_summary$coefficients[2, 4]
    p_value <- fit_summary$coefficients[2, 5]
    effect_size <- r2beta(fit, method = "nsj")
    effect_size <- effect_size[rownames(effect_size) == "2", "Rsq"]
    return(c(t_value, p_value, effect_size))
  }
  
  # rs-fMRI
  lmm_rsfmri <- function(x, df) {
    cat("\n", x, "Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ df[[var]] + sex + interview_age + race_ethnicity +
                  income_parent + edu_parent + rsfmri_c_ngd_meanmotion +
                  (1|mri_info_deviceserialnumber/rel_family_id), data = df)
    fit_summary <- summary(fit)
    # 3 df; 4 t-value; 5 p-value
    t_value <- fit_summary$coefficients[2, 4]
    p_value <- fit_summary$coefficients[2, 5]
    effect_size <- r2beta(fit, method = "nsj")
    effect_size <- effect_size[rownames(effect_size) == "2", "Rsq"]
    return(c(t_value, p_value, effect_size))
  }
  
  # sMRI
  lmm_smri <- function(x, df) {
    cat("\n", x, "Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ df[[var]] + sex + interview_age + race_ethnicity + 
                  income_parent + edu_parent + smri_vol_scs_intracranialv + 
                  (1|mri_info_deviceserialnumber/rel_family_id), data = df)
    fit_summary <- summary(fit)
    # 3 df; 4 t-value; 5 p-value
    t_value <- fit_summary$coefficients[2, 4]
    p_value <- fit_summary$coefficients[2, 5]
    effect_size <- r2beta(fit, method = "nsj")
    effect_size <- effect_size[rownames(effect_size) == "2", "Rsq"]
    return(c(t_value, p_value, effect_size))
  }
  
  result_behavior <- sapply(behavior_year2, lmm_behavior, abcd_year2_behavior)
  result_networks <- sapply(fmri_networks, lmm_rsfmri, abcd_year2_rsfmri)
  result_network_subcortical <- sapply(fmri_network_subcortical, lmm_rsfmri, 
                                       abcd_year2_rsfmri)
  result_thick <- sapply(smri_thick, lmm_smri, abcd_year2_smri)
  result_area <- sapply(smri_area, lmm_smri, abcd_year2_smri)
  result_vol_cortical <- sapply(smri_vol_cortical, lmm_smri, abcd_year2_smri)
  result_vol_subcortical <- sapply(smri_vol_subcortical, lmm_smri, 
                                   abcd_year2_smri)
  result_smri <- sapply(smri_total, lmm_smri, abcd_year2_smri)
  
  # output
  t_behavior <- result_behavior[1, ]
  t_volume <- result_vol_cortical[1, ]
  t_subvolume <- result_vol_subcortical[1, ]
  t_area <- result_area[1, ]
  t_thickness <- result_thick[1, ]
  t_networks <- result_networks[1, ]
  t_subnetworks <- result_network_subcortical[1, ]
  t_smri <- result_smri[1, ]
  
  p_behavior <- result_behavior[2, ]
  p_volume <- result_vol_cortical[2, ]
  p_subvolume <- result_vol_subcortical[2, ]
  p_area <- result_area[2, ]
  p_thickness <- result_thick[2, ]
  p_networks <- result_networks[2, ]
  p_subnetworks <- result_network_subcortical[2, ]
  p_smri <- result_smri[2, ]
  
  R2_behavior <- result_behavior[3, ]
  R2_volume <- result_vol_cortical[3, ]
  R2_subvolume <- result_vol_subcortical[3, ]
  R2_area <- result_area[3, ]
  R2_thickness <- result_thick[3, ]
  R2_networks <- result_networks[3, ]
  R2_subnetworks <- result_network_subcortical[3, ]
  R2_smri <- result_smri[3, ]
  
  results <- c(t_behavior, t_volume, t_subvolume, t_area, t_thickness,
    t_networks, t_subnetworks, t_smri, R2_behavior, R2_volume, R2_subvolume, 
    R2_area, R2_thickness, R2_networks, R2_subnetworks, R2_smri, p_behavior, 
    p_volume, p_subvolume, p_area, p_thickness, p_networks, p_subnetworks, p_smri)
  
  all_vars <- c(behavior_year2, smri_vol_cortical, smri_vol_subcortical,
    smri_area, smri_thick, fmri_networks, fmri_network_subcortical, smri_total)
  names(results) <- c(paste0("t_", all_vars), paste0("R2_", all_vars), 
                      paste0("p_", all_vars))
  
  return(results)
}

# boot -------------------------------------------------------------------------
nboot <- 500

# boot matrix, 80% samples
set.seed(123)
boot_behavior_mat <- createDataPartition(seq(nrow(abcd_year2_behavior)), 
  nboot, 0.80, list = FALSE)
boot_smri_mat <- createDataPartition(seq(nrow(abcd_year2_smri)), 
  nboot, 0.80, list = FALSE)
boot_rsfmri_mat <- createDataPartition(seq(nrow(abcd_year2_rsfmri)), 
  nboot, 0.80, list = FALSE)

# set do parallel
cl <- makeCluster(14, type = 'PSOCK')
clusterExport(cl, c("abcd_year2_behavior", "abcd_year2_smri", 
  "abcd_year2_rsfmri", "boot_behavior_mat", "boot_smri_mat", 
  "boot_rsfmri_mat"))
registerDoParallel(cl)

# boot PFI
boot_PFI <- foreach(boot_seq = seq(nboot), 
                    .packages = c("tidyverse", "lmerTest", "r2glmm"), 
                    .combine = rbind) %dopar% 
  boot_lmm(boot_seq, "pbp_ss_prosocial_peers")
# boot DFI
boot_DFI <- foreach(boot_seq = seq(nboot), 
                    .packages = c("tidyverse", "lmerTest", "r2glmm"), 
                    .combine = rbind) %dopar% 
  boot_lmm(boot_seq, "pbp_ss_rule_break")

# mean and SD of significant regions -------------------------------------------
significant_boot_t <- function(data, data_fdr, data_boot, peer) {
  profiles <- function(var1, var2, data_boot) {
    index <- which(colnames(data_boot) %in% paste0(var2, var1))
    df <- data_boot[, index]
    df <- as.data.frame(df)
    return(df)
  }
  
  data_fdr[["sMRI"]] <- NULL
  data_fdr <- data_fdr[c("behavior", "vol_cortical", "vol_subcortical", "thick", 
    "area", "networks", "network_sub")]
  
  sigificants <- unlist(sapply(data_fdr, colnames))
  
  # primary results
  primary_t <- sapply(data_fdr, function(x) x <- unlist(x[1, ])) %>%
    unlist()
  primary_R2 <- sapply(data_fdr, function(x) x <- unlist(x[3, ])) %>%
    unlist()
  
  # boot t
  boot_sigificants_t <- profiles(sigificants, "t_", data_boot)
  boot_sigificants_t_mean <- colMeans(boot_sigificants_t)
  boot_sigificants_t_sd <- apply(boot_sigificants_t, 2, sd)
  
  # boot R2
  boot_sigificants_R2 <- profiles(sigificants, "R2_", data_boot)
  boot_sigificants_R2_mean <- colMeans(boot_sigificants_R2)
  boot_sigificants_R2_sd <- apply(boot_sigificants_R2, 2, sd)
  
  # boot sum 
  significant_boot_sum <- function(data, data_fdr, data_boot) {
    sigificants <- sapply(data_fdr, colnames) %>%
      unlist()
    sigificants_p <- paste0("p_", sigificants)
    
    sigficant_boot_times <- function(sets, data_boot, data) {
      vars <- colnames(data[[sets]])
      
      index_p <- which(colnames(data_boot) %in% paste0("p_", vars))
      
      df_p <- data_boot[, index_p]
      
      df_p_fdr <- apply(df_p, 1, p.adjust, method = "fdr")
      
      get_times <- function(x) {
        times <- ifelse(x < 0.05, 1, 0)
        times <- sum(times)
        return(times)
      }
      
      significant_times <- apply(df_p_fdr, 1, get_times)
      
      return(significant_times)
    }
    
    sigficant_times <- lapply(names(data)[1:7], sigficant_boot_times, data_boot, data) %>%
      unlist()
    sigficant_times <- sigficant_times[sigificants_p]
    sigficant_times <- sigficant_times / nboot
    return(sigficant_times)
  }
  t_boot_sum <- significant_boot_sum(data, data_fdr, data_boot)
  
  # summary results
  df <- data.frame(
    "Variables" = sigificants, 
    "t_primary" = primary_t, 
    "t_sampling_meansd" = paste0(round(boot_sigificants_t_mean, 2), 
      " (", round(boot_sigificants_t_sd, 2), ")"), 
    "R2_primary" = sprintf("%.1e", boot_sigificants_R2_mean),
    "R2_sampling_meansd" = paste0(sprintf("%.1e", boot_sigificants_R2_mean),
      " (", sprintf("%.1e", boot_sigificants_R2_sd), ")"),
    "t_sampling_sum" = paste0(t_boot_sum * 100, "%"))
  write_csv(df, paste0("supplementary_analysis/", peer, "_boot_summary.csv"))
  
  return(df)
}

# PFI
PFI_df_t <- significant_boot_t(PFI, PFI1, boot_PFI, "PFI")
# DFI
DFI1$vol_cortical <- as.matrix(DFI1$vol_cortical)
colnames(DFI1$vol_cortical) <- "smri_vol_cdk_locclh"
DFI_df_t <- significant_boot_t(DFI, DFI1, boot_DFI, "DFI")
