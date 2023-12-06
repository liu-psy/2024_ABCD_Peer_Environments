# Supplementary analyses
# addition confounds: the number of male and female friends,
# addition confounds: family conflict, neighborhood safety, 
# addition confounds: school protective and risk factors.
library(r2glmm)
library(lavaan)
library(lmerTest)
library(tidyverse)

# set working directory
setwd("H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/supplementary_analysis")
# load("supplementary_results.RData")

# load data
abcd <- read_csv("H:/ABCD/Relsease4.0/Package_1194636/abcd.csv")
base_information <- read_csv("H:/ABCD/Relsease4.0/Package_1194636/base_information.csv")
abcd_year2 <- filter(abcd, eventname == "2_year_follow_up_y_arm_1")
abcd_year3 <-  filter(abcd, eventname == "3_year_follow_up_y_arm_1")

# select variables -------------------------------------------------------------
# behavior: abcd_mhp02 abcd_mhy02 abcd_sscey01 abcd_sscep01
behavior_year2 <- select(abcd_year2, nihtbx_picvocab_uncorrected:peq_ss_relational_aggs,
  -pbp_ss_prosocial_peers, -pbp_ss_rule_break)
behavior_year3 <- select(abcd_year3, nihtbx_picvocab_uncorrected:peq_ss_relational_aggs,
  -pbp_ss_prosocial_peers, -pbp_ss_rule_break)

behavior_na <- data.frame(
  "Year2" = apply(behavior_year2, 2, function(x) sum(is.na(x))),
  "Year3" = apply(behavior_year3, 2, function(x) sum(is.na(x))))
# exclude variables with too much NA value
behavior_year2 <- names(behavior_year2[, behavior_na$Year2 <= 2329])
behavior_year2 <- behavior_year2[-which(behavior_year2 == "ple_y_ss_affected_bad_mean")]
behavior_year3 <- names(behavior_year3[, behavior_na$Year3 <= 887])
behavior_year3 <- behavior_year3[behavior_year3 %in% behavior_year2]
length(behavior_year3)

# basic variables
basic <- c("subjectkey", "interview_age", "sex", "race_ethnicity",
  "rel_family_id", "income_parent", "edu_parent")
peers <- c("pbp_ss_prosocial_peers", "pbp_ss_rule_break")

# additional confounding variables
confounds_year2 <- c("resiliency5a_y", "resiliency6a_y", "fes_p_ss_fc",
  "nsc_p_ss_mean_3_items", "srpf_y_ss_ses", "srpf_y_ss_iiss", "srpf_y_ss_dfs")
confounds_year3 <- c("resiliency5a_y", "resiliency6a_y", "fes_p_ss_fc",
  "srpf_y_ss_ses", "srpf_y_ss_iiss", "srpf_y_ss_dfs")

# MRI data
smri_thick <- names(select(abcd_year2, contains("_thick_cdk")))[-c(69:71)]
smri_area <- names(select(abcd_year2, contains("_area_cdk")))[-c(69:71)]
smri_vol_cortical <- names(select(abcd_year2, contains("_vol_cdk")))[-c(69:71)]
smri_vol_subcortical <- names(select(abcd_year2, contains(paste0("vol_scs_", 
  c("tp", "caudate", "putamen", "pallidum", "hpus", "amygdala", "aa", "vedc")))))
smri_total <- c("smri_vol_cdk_total", "smri_thick_cdk_mean", "smri_area_cdk_total")
smri <- c(smri_thick, smri_area, smri_vol_cortical, smri_vol_subcortical, 
  smri_total)

# RSFCs
fmri_networks <- names(select(abcd_year2, rsfmri_c_ngd_ad_ngd_ad:rsfmri_c_ngd_vta_ngd_vs))
fmri_network_subcortical <- names(select(abcd_year2, rsfmri_cor_ngd_au_scs_crcx:rsfmri_cor_ngd_vs_scs_bs))
fc_network <- c(fmri_networks, fmri_network_subcortical)

# exclude missing value and MRI QC
abcd_year2_behavior <- select(abcd_year2, all_of(basic), all_of(peers),
    all_of(confounds_year2), all_of(behavior_year2), site_id_l) %>%
  filter(complete.cases(.), !duplicated(.))

abcd_year2_smri <- select(abcd_year2, all_of(basic), all_of(peers), all_of(smri),
    all_of(confounds_year2), mri_info_deviceserialnumber,
    smri_vol_scs_intracranialv, imgincl_t1w_include) %>%
  filter(imgincl_t1w_include == 1, complete.cases(.), !duplicated(.))

abcd_year2_rsfmri <- select(abcd_year2, all_of(basic), all_of(peers),
    all_of(confounds_year2), all_of(fc_network), mri_info_deviceserialnumber,
    rsfmri_c_ngd_meanmotion, imgincl_rsfmri_include) %>%
  filter(imgincl_rsfmri_include == 1, complete.cases(.), !duplicated(.))

abcd_year3_behavior <- select(abcd_year3, all_of(basic), all_of(peers),
    all_of(confounds_year3), all_of(behavior_year3), site_id_l) %>%
  filter(complete.cases(.), !duplicated(.))

# longitudinal data
abcd_longitudinal_year2 <- abcd_year2_behavior[abcd_year2_behavior$subjectkey %in% abcd_year3_behavior$subjectkey, ]
abcd_longitudinal_year3 <- abcd_year3_behavior[abcd_year3_behavior$subjectkey %in% abcd_year2_behavior$subjectkey, ]

# 1. Association Analysis ------------------------------------------------------
lmm_vars <- function(var) {
  # behavior
  lmm_behavior <- function(x, df) {
    cat("\n", x, "   Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ df[[var]] + sex + interview_age + race_ethnicity +
      income_parent + edu_parent + resiliency5a_y + resiliency6a_y +
      fes_p_ss_fc + nsc_p_ss_mean_3_items + srpf_y_ss_ses + srpf_y_ss_iiss +
      srpf_y_ss_dfs + (1|site_id_l/rel_family_id), data = df)
    fit_summary <- summary(fit)
    # 3 df; 4 t-value; 5 p-value
    t_value <- fit_summary$coefficients[2, 4]
    p_value <- fit_summary$coefficients[2, 5]
    effect_size <- r2beta(fit, method = "nsj")$Rsq[2]
    return(list(t_value, p_value, effect_size))
  }

  # rs-fMRI
  lmm_rsfmri <- function(x, df) {
    cat("\n", x, "Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ df[[var]] + sex + interview_age + race_ethnicity +
      income_parent + edu_parent + resiliency5a_y + resiliency6a_y +
      fes_p_ss_fc + nsc_p_ss_mean_3_items + srpf_y_ss_ses + srpf_y_ss_iiss +
      srpf_y_ss_dfs + rsfmri_c_ngd_meanmotion +
      (1|mri_info_deviceserialnumber/rel_family_id), data = df)
    fit_summary <- summary(fit)
    # 3 df; 4 t-value; 5 p-value
    t_value <- fit_summary$coefficients[2, 4]
    p_value <- fit_summary$coefficients[2, 5]
    effect_size <- r2beta(fit, method = "nsj")$Rsq[2]
    return(list(t_value, p_value, effect_size))
  }

  # sMRI
  lmm_smri <- function(x, df) {
    cat("\n", x, "Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ df[[var]] + sex + interview_age + race_ethnicity +
      income_parent + edu_parent + resiliency5a_y + resiliency6a_y +
      fes_p_ss_fc + nsc_p_ss_mean_3_items + srpf_y_ss_ses +
      srpf_y_ss_iiss + srpf_y_ss_dfs + smri_vol_scs_intracranialv +
      (1|mri_info_deviceserialnumber/rel_family_id), data = df)
    fit_summary <- summary(fit)
    # 3 df; 4 t-value; 5 p-value
    t_value <- fit_summary$coefficients[2, 4]
    p_value <- fit_summary$coefficients[2, 5]
    effect_size <- r2beta(fit, method = "nsj")$Rsq[2]
    return(list(t_value, p_value, effect_size))
  }

  result_behavior <- sapply(behavior_year2, lmm_behavior, abcd_year2_behavior)
  result_networks <- sapply(fmri_networks, lmm_rsfmri, abcd_year2_rsfmri)
  result_network_subcortical <- sapply(fmri_network_subcortical, lmm_rsfmri, abcd_year2_rsfmri)
  result_thick <- sapply(smri_thick, lmm_smri, abcd_year2_smri)
  result_area <- sapply(smri_area, lmm_smri, abcd_year2_smri)
  result_vol_cortical <- sapply(smri_vol_cortical, lmm_smri, abcd_year2_smri)
  result_vol_subcortical <- sapply(smri_vol_subcortical, lmm_smri, abcd_year2_smri)
  result_smri <- sapply(smri_total, lmm_smri, abcd_year2_smri)

  # output
  results <- list(result_behavior, result_networks, result_network_subcortical,
    result_thick, result_area, result_vol_cortical, result_vol_subcortical,
    result_smri)
  names(results) <- c("behavior", "networks",  "network_sub", "thick", "area",
                      "vol_cortical", "vol_subcortical", "sMRI")

  # add rownames
  modify <- function(x) {
    rownames(x) <- c("t", "p", "effect_size")
    x[1, ] <- round(unlist(x[1, ]), 2)
    x[3, ] <- round(unlist(x[3, ]), 3)
    return(x)
  }
  results <- lapply(results, modify)
  
  return(results)
}
p_adjust <- function(x, p_correction) {
  x[2, ] <- p.adjust(x[2, ], method = p_correction)
  x <- x[, x[2, ] < 0.05]
  return(x)
}

# PFI
PFI <- lmm_vars("pbp_ss_prosocial_peers")
PFI1 <- lapply(PFI, p_adjust, "fdr")

# DFI
DFI <- lmm_vars("pbp_ss_rule_break")
DFI1 <- lapply(DFI, p_adjust, "fdr")

common_behaviors <- colnames(PFI1$behavior)[colnames(PFI1$behavior) %in% colnames(DFI1$behavior)]

# sub-cortical volumes
PFI_sub <- data.frame(
  "vars" = colnames(PFI1$vol_subcortical),
  "t" = round(unlist(PFI1$vol_subcortical[1, ]), 2),
  "R2" = round(unlist(PFI1$vol_subcortical[3, ]), 3), 
  "p" = sprintf("%.1e", unlist(PFI1$vol_subcortical[2, ])),
  "peer" = "PFI"
)
write_csv(PFI_sub, "association_sub.csv")

# 2. Neurotransmitters ---------------------------------------------------------
# brain t map
vol_PFI <- unlist(PFI$vol_cortical[1, ])
vol_DFI <- unlist(DFI$vol_cortical[1, ])
thick_PFI <- unlist(PFI$thick[1, ])
thick_DFI <- unlist(DFI$thick[1, ])
area_PFI <- unlist(PFI$area[1, ])
area_DFI <- unlist(DFI$area[1, ])

write.table(vol_PFI, "tmap_vol_PFI.txt", col.names = FALSE, row.names = FALSE)
write.table(vol_DFI, "tmap_vol_DFI.txt", col.names = FALSE, row.names = FALSE)
write.table(thick_PFI, "tmap_thick_PFI.txt", col.names = FALSE, row.names = FALSE)
write.table(thick_DFI, "tmap_thick_DFI.txt", col.names = FALSE, row.names = FALSE)
write.table(area_PFI, "tmap_area_PFI.txt", col.names = FALSE, row.names = FALSE)
write.table(area_DFI, "tmap_area_DFI.txt", col.names = FALSE, row.names = FALSE)

# load results (results from JuSpace) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Neurotransmitters <- read.csv("Resh.csv") %>%
  rename("Spearman" = Mean.Fisher.s.z..Spearman.rho., "PET" = PET.Map) %>%
  select(File, PET, Spearman, p_exact) %>%
  mutate(PET = str_replace(PET, "SERT", "5HTT")) %>%
  mutate(File = factor(File, levels = unique(File)))

# FDR corrections
p_fdr <- Neurotransmitters %>%
  group_by(File) %>%
  summarise_at(vars(p_exact), p.adjust, method = "fdr")
Neurotransmitters$p_fdr <- p_fdr$p_exact
Neurotransmitters$significants <- ifelse(Neurotransmitters$p_fdr < 0.01, "**",
  ifelse(Neurotransmitters$p_fdr < 0.05, "*", " "))

# volume
Neurotransmitters_vol <- Neurotransmitters[str_detect(Neurotransmitters$File, "vol"), ]
# area
Neurotransmitters_area <- Neurotransmitters[str_detect(Neurotransmitters$File, "area"), ]
# thick
Neurotransmitters_thick <- Neurotransmitters[str_detect(Neurotransmitters$File, "thick"), ]

# 3. Mediation Analysis --------------------------------------------------------
# cortical and subcortical volumes
mediation_volume <- select(abcd_year2, all_of(basic), all_of(confounds_year2),
    mri_info_deviceserialnumber, smri_vol_scs_intracranialv,  all_of(peers),
    all_of(behavior_year2), smri_vol_cdk_total, colnames(PFI1$vol_cortical),
    colnames(PFI1$vol_subcortical), smri_vol_cdk_locclh, imgincl_t1w_include) %>%
  filter(imgincl_t1w_include == 1, complete.cases(.)) %>%
  select(-imgincl_t1w_include) %>%
  mutate(sex = factor(sex, levels = c("F", "M"), labels = c(1, 0))) %>%
  mutate(mri_info_deviceserialnumber = factor(mri_info_deviceserialnumber,
    levels = unique(mri_info_deviceserialnumber),
    labels = seq(unique(mri_info_deviceserialnumber))))
write_csv(mediation_volume, "mediation_volume.csv")

# RSFC
mediation_rsfc <- select(abcd_year2, all_of(basic), all_of(confounds_year2),
    mri_info_deviceserialnumber, rsfmri_c_ngd_meanmotion, all_of(peers),
    all_of(behavior_year2), colnames(PFI1$networks), colnames(DFI1$networks),
    colnames(DFI1$network_sub), imgincl_rsfmri_include) %>%
  filter(imgincl_rsfmri_include == 1, complete.cases(.)) %>%
  select(-imgincl_rsfmri_include) %>%
  mutate(sex = factor(sex, levels = c("F", "M"), labels = c(1, 0))) %>%
  mutate(mri_info_deviceserialnumber = factor(mri_info_deviceserialnumber,
    levels = unique(mri_info_deviceserialnumber),
    labels = seq(unique(mri_info_deviceserialnumber))))
write_csv(mediation_rsfc, "mediation_rsfc.csv")

# mediation analysis in Matlab  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mediation_fdr <- function(strs, pattern, peers, replace_index) {
  
  # read mediation results
  read_csvs <- function(files) {
    # set row names
    paths_p <- paste0(c("PathA", "PathB", "PathCC", "PathC", "PathAB"), "_p")
    paths_t <- paste0(c("PathA", "PathB", "PathCC'", "PathC", "PathAB"), "_beta")
    CI <- paste0(rep(c("A", "B", "CC", "C", "AB"), 2), rep(1:2, each = 5))
    paths <- c(paths_p, paths_t, CI)
    
    # read files
    csv_file <- paste0(files, ".csv")
    result <- read.csv(csv_file, header = FALSE)
    
    colnames(result) <- paths
    rownames(result) <- paste(files, "-", behavior_year2)
    return(result)
  }
  
  sig_vars <- str_split(strs, pattern, simplify = TRUE)[, 2]
  sig_vars <- str_replace(sig_vars, replace_index, "")
  sig_vars <- paste0(peers, sig_vars)
  
  mediation_results <- lapply(sig_vars, read_csvs)
  df <- data.frame()
  for (i in seq(mediation_results)) {
    df <- rbind(df, mediation_results[[i]])
  }
  # FDR corrections
  df[1:5] <- apply(df[1:5], 2, p.adjust, method = "fdr")
  df <- filter(df, PathAB_p < 0.05, PathA_p < 0.05, PathB_p < 0.05,
               PathC_p < 0.05, PathCC_p < 0.05)
  return(df)
}

# mediation - PFI 7 regional volume (22 significant results)
PFI_sig_volumes <- c(colnames(PFI1$vol_cortical), colnames(PFI1$vol_subcortical)) %>%
  str_replace("scs", "cdk")
PFI_volume_all <- mediation_fdr(PFI_sig_volumes, "smri_vol_", 
  "PFI_", "cdk_")

# mediation - DFI 1 regional volume (no significant results)
DFI_volume_all <- mediation_fdr("smri_vol_cdk_locclh", "smri_vol_", 
  "DFI_", "cdk_")

# mediation - PFI 2 RSFCs (2 significant results)
PFI_fc <- mediation_fdr(colnames(PFI1$networks), "rsfmri_c_ngd_", 
  "PFI_", "ngd_")

# mediation - DFI 6 RSFCs (21 significant results)
DFI_fc <- mediation_fdr(colnames(DFI1$networks), "rsfmri_c_ngd_", 
  "DFI_", "ngd_")

# mediation - DFI 22 cortico-subcortical RSFCs (67 significant results)
DFI_sub <- mediation_fdr(colnames(DFI1$network_sub), "rsfmri_cor_ngd_", 
  "DFI_", "scs_")

# 4. Cross Lag Panel Model -----------------------------------------------------
# time: year2 - year3
# time: year2 - year3
clpms <- function(x) {

  # regressed out covariates (2YFU)
  residual_year2 <- function(y, data) {
    fit <- lmer(data[[y]] ~ sex + interview_age + race_ethnicity +
      income_parent + edu_parent + resiliency5a_y + resiliency6a_y +
      fes_p_ss_fc + nsc_p_ss_mean_3_items + srpf_y_ss_ses + srpf_y_ss_iiss +
      srpf_y_ss_dfs + (1|site_id_l/rel_family_id), data = data)
    return(residuals(fit))
  }
  # regressed out covariates (3YFU)
  residual_year3 <- function(y, data) {
    fit <- lmer(data[[y]] ~ sex + interview_age + race_ethnicity +
      income_parent + edu_parent + resiliency5a_y + resiliency6a_y +
      fes_p_ss_fc  + srpf_y_ss_ses + srpf_y_ss_iiss + srpf_y_ss_dfs +
      (1|site_id_l/rel_family_id), data = data)
    return(residuals(fit))
  }

  residual_year2_PFI <- residual_year2("pbp_ss_prosocial_peers", abcd_year2_behavior)
  residual_year3_PFI <- residual_year3("pbp_ss_prosocial_peers", abcd_year3_behavior)
  residual_year2_DFI <- residual_year2("pbp_ss_rule_break", abcd_year2_behavior)
  residual_year3_DFI <- residual_year3("pbp_ss_rule_break", abcd_year3_behavior)

  residual_year2_x <- residual_year2(x, abcd_year2_behavior)
  residual_year3_x <- residual_year3(x, abcd_year3_behavior)

  residualed_year2 <- data.frame(
    "subjectkey" = abcd_year2_behavior$subjectkey,
    "PFI_year2" = residual_year2_PFI,
    "DFI_year2" = residual_year2_DFI,
    "x_year2" = residual_year2_x
  )
  residualed_year3 <- data.frame(
    "subjectkey" = abcd_year3_behavior$subjectkey,
    "PFI_year3" = residual_year3_PFI,
    "DFI_year3" = residual_year3_DFI,
    "x_year3" = residual_year3_x
  )

  clpm_data <- inner_join(residualed_year2, residualed_year3, by = "subjectkey") %>%
    select(-subjectkey) %>%
    apply(2, scale)

  clpm_PFI <- '
    # Estimate the lagged effects between the observed variables.
    PFI_year3 + x_year3 ~ PFI_year2 + x_year2

    # Estimate the covariance between the observed variables at the first wave.
    PFI_year2 ~~ x_year2 # Covariance

    # Estimate the covariances between the residuals of the observed variables.
    PFI_year3 ~~ x_year3

    # Estimate the (residual) variance of the observed variables.
    PFI_year2 ~~ PFI_year2 # Variances
    x_year2 ~~ x_year2
    PFI_year3 ~~ PFI_year3 # Residual variances
    x_year3 ~~ x_year3
  '
  clpm_DFI <- '
    DFI_year3 + x_year3 ~ DFI_year2 + x_year2

    DFI_year2 ~~ x_year2 # Covariance

    DFI_year3 ~~ x_year3

    DFI_year2 ~~ DFI_year2
    x_year2 ~~ x_year2
    DFI_year3 ~~ DFI_year3
    x_year3 ~~ x_year3
  '

  CLPM_PFI <- lavaan(clpm_PFI, data = clpm_data, meanstructure = TRUE, int.ov.free = TRUE)
  CLPM_DFI <- lavaan(clpm_DFI, data = clpm_data, meanstructure = TRUE, int.ov.free = TRUE)
  result_PFI <- summary(CLPM_PFI, ci = TRUE)$pe
  result_DFI <- summary(CLPM_DFI, ci = TRUE)$pe

  result <- rbind(result_PFI[1:6, ], result_DFI[1:6, ])
  return(result)
}

clpm_all <- lapply(behavior_year3, clpms)
names(clpm_all) <- behavior_year3
clpm_all

# FDR corrections
time1_time2_PFI <- vector()
crosslag_x1_PFI <- vector()
crosslag_PFI_x1 <- vector()
time1_time2_x1 <- vector()
time1_x1_PFI <- vector()
time2_x1_PFI <- vector()

time1_time2_DFI <- vector()
crosslag_x2_DFI <- vector()
crosslag_DFI_x2 <- vector()
time1_time2_x2 <- vector()
time1_x2_DFI <- vector()
time2_x2_DFI <- vector()

# get p-values
for(i in seq(behavior_year3)) {
  time1_time2_PFI[i] <- clpm_all[[i]][1, 8]
  crosslag_x1_PFI[i] <- clpm_all[[i]][2, 8]
  crosslag_PFI_x1[i] <- clpm_all[[i]][3, 8]
  time1_time2_x1[i] <- clpm_all[[i]][4, 8]
  time1_x1_PFI[i] <- clpm_all[[i]][5, 8]
  time2_x1_PFI[i] <- clpm_all[[i]][6, 8]

  time1_time2_DFI[i] <- clpm_all[[i]][7, 8]
  crosslag_x2_DFI[i] <- clpm_all[[i]][8, 8]
  crosslag_DFI_x2[i] <- clpm_all[[i]][9, 8]
  time1_time2_x2[i] <- clpm_all[[i]][10, 8]
  time1_x2_DFI[i] <- clpm_all[[i]][11, 8]
  time2_x2_DFI[i] <- clpm_all[[i]][12, 8]
}

cross_lags <- data.frame(
  # PFI
  time1_time2_PFI,
  crosslag_x1_PFI,
  crosslag_PFI_x1,
  time1_time2_x1,
  time1_x1_PFI,
  time2_x1_PFI,
  # DFI
  time1_time2_DFI,
  crosslag_x2_DFI,
  crosslag_DFI_x2,
  time1_time2_x2,
  time1_x2_DFI,
  time2_x2_DFI
)

# FDR correction
cross_lags_fdr <- apply(cross_lags, 2, p.adjust, method = "fdr") %>%
  data.frame()
rownames(cross_lags_fdr) <- behavior_year3
cross_lags_fdr

# get betas of cross-lagged effect of peer environments
cross_lag_data <- function(sig_vars, peer) {
  beta_vector <- vector()

  if (peer == "PFI") {
    for(i in seq(sig_vars)) {
      beta_vector[i] <- clpm_all[[sig_vars[i]]][3, 5]
    }
  }
  if (peer == "DFI") {
    for(i in seq(sig_vars)) {
      beta_vector[i] <- clpm_all[[sig_vars[i]]][9, 5]
    }
  }

  df <- data.frame(
    "variables" = sig_vars,
    "beta" = round(beta_vector, 3)
  )
  return(df)
}

# PFI
x_PFI <- behavior_year3[cross_lags_fdr$crosslag_x1_PFI < 0.05]
PFI_x <- behavior_year3[cross_lags_fdr$crosslag_PFI_x1 < 0.05]
common_PFI <- x_PFI[x_PFI %in% PFI_x]
PFI_beta <- cross_lag_data(PFI_x, "PFI")

# DFI
x_DFI <- behavior_year3[cross_lags_fdr$crosslag_x2_DFI < 0.05]
DFI_x <- behavior_year3[cross_lags_fdr$crosslag_DFI_x2 < 0.05]
common_DFI <- x_DFI[x_DFI %in% DFI_x]
DFI_beta <- cross_lag_data(DFI_x, "DFI")

# no common variables
PFI_x[PFI_x %in% DFI_x]

# output -----------------------------------------------------------------------
# demographic
basic_stat <- function(data) {
  N <- nrow(data)
  female_percent <- (sum(data$sex %in% c("F", 1)) / nrow(data) * 100) %>% round(2)
  mean_age <- (mean(data$interview_age) / 12) %>% round(2)
  sd_age <- (sd(data$interview_age) / 12) %>% round(2)
  income_parent <- mean(data$income_parent) %>% round(2)
  sd_income_parent <- sd(data$income_parent) %>% round(2)
  edu_parent <- mean(data$edu_parent) %>% round(2)
  sd_edu_parent <- sd(data$edu_parent) %>% round(2)
  race <- unlist(table(data$race_ethnicity))
  return(c(N, female_percent, mean_age, sd_age, income_parent, sd_income_parent,
    edu_parent, sd_edu_parent, race))
}

demographic_table <- rbind(basic_stat(abcd_year2_behavior),
    basic_stat(abcd_year2_smri), basic_stat(abcd_year2_rsfmri),
    basic_stat(mediation_volume), basic_stat(mediation_rsfc),
    basic_stat(abcd_longitudinal_year2),basic_stat(abcd_longitudinal_year3)) %>%
  as.data.frame()
rownames(demographic_table) <- c("behavior", "smri", "rsfmri", 
  "mediation_volume", "mediation_rsfc", "longitudinal_year2",
  "longitudinal_year3")
write_csv(demographic_table, "demographic_table.csv")

# results
PFI_result <- list(PFI, PFI1, Neurotransmitters, PFI_volume_all, PFI_fc, 
  PFI_beta)
DFI_result <- list(DFI, DFI1, Neurotransmitters, DFI_volume_all, DFI_fc, 
  DFI_sub, DFI_beta)

# save working space
save.image("supplementary_results.RData")