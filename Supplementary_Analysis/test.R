# plot effect sizes
library(lmerTest)
library(r2glmm)
library(patchwork)
library(tidyverse)

# set working directory
setwd("H:/ABCD/Release4.0/Package_1194636/results/peer_environments")
load("CV_results.RData")

lmm_vars <- function(var) {
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
    effect_size_peer <- effect_size[rownames(effect_size) == "2", "Rsq"]
    effect_size_sex <- effect_size[rownames(effect_size) == "3", "Rsq"]
    effect_size_age <- effect_size[rownames(effect_size) == "4", "Rsq"]
    effect_size_race <- effect_size[rownames(effect_size) == "5", "Rsq"]
    effect_size_income <- effect_size[rownames(effect_size) == "6", "Rsq"]
    effect_size_education <- effect_size[rownames(effect_size) == "7", "Rsq"]
    return(list(t_value, p_value, effect_size_peer, effect_size_sex, effect_size_age,
                effect_size_race, effect_size_income, effect_size_education))
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
    effect_size_peer <- effect_size[rownames(effect_size) == "2", "Rsq"]
    effect_size_sex <- effect_size[rownames(effect_size) == "3", "Rsq"]
    effect_size_age <- effect_size[rownames(effect_size) == "4", "Rsq"]
    effect_size_race <- effect_size[rownames(effect_size) == "5", "Rsq"]
    effect_size_income <- effect_size[rownames(effect_size) == "6", "Rsq"]
    effect_size_education <- effect_size[rownames(effect_size) == "7", "Rsq"]
    return(list(t_value, p_value, effect_size_peer, effect_size_sex, effect_size_age,
                effect_size_race, effect_size_income, effect_size_education))
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
    effect_size_peer <- effect_size[rownames(effect_size) == "2", "Rsq"]
    effect_size_sex <- effect_size[rownames(effect_size) == "3", "Rsq"]
    effect_size_age <- effect_size[rownames(effect_size) == "4", "Rsq"]
    effect_size_race <- effect_size[rownames(effect_size) == "5", "Rsq"]
    effect_size_income <- effect_size[rownames(effect_size) == "6", "Rsq"]
    effect_size_education <- effect_size[rownames(effect_size) == "7", "Rsq"]
    return(list(t_value, p_value, effect_size_peer, effect_size_sex, effect_size_age,
                effect_size_race, effect_size_income, effect_size_education))
  }
  
  result_behavior <- sapply(behavior_year2, lmm_behavior, abcd_year2_behavior)
  result_vol_cortical <- sapply(smri_vol_cortical, lmm_smri, abcd_year2_smri)
  result_area <- sapply(smri_area, lmm_smri, abcd_year2_smri)
  result_thick <- sapply(smri_thick, lmm_smri, abcd_year2_smri)
  result_vol_subcortical <- sapply(smri_vol_subcortical, lmm_smri, 
                                   abcd_year2_smri)
  result_networks <- sapply(fmri_networks, lmm_rsfmri, abcd_year2_rsfmri)
  result_network_subcortical <- sapply(fmri_network_subcortical, lmm_rsfmri, 
                                       abcd_year2_rsfmri)
  result_smri <- sapply(smri_total, lmm_smri, abcd_year2_smri)
  
  # output
  results <- list(result_behavior, result_vol_cortical, result_area, 
                  result_thick, result_vol_subcortical, result_networks, 
                  result_network_subcortical, result_smri)
  names(results) <- c("behavior", "vol_cortical", "area", "thick", 
                      "vol_subcortical", "networks",  "network_sub", "sMRI")
  # add rownames
  modify <- function(x) {
    rownames(x) <- c("t", "p", "effect_peer", "effect_sex", "effect_age",
                     "effect_race", "effect_income", "effect_education")
    # x[1, ] <- round(unlist(x[1, ]), 2)
    # x[3, ] <- round(unlist(x[3, ]), 4)
    return(x)
  }
  results <- lapply(results, modify)
  
  return(results)
}

p_adjust <- function(x, p_correction) {
  x[2, ] <- p.adjust(x[2, ], method = p_correction)
  # x <- x[, x[2, ] < 0.05]
  return(x)
}

# PFI
PFI <- lmm_vars("pbp_ss_prosocial_peers")
PFI1 <- lapply(PFI, p_adjust, "fdr")
PFI1$sMRI <- NULL
PFI1 <- do.call(cbind, PFI1)
PFI1 <- PFI1[, PFI1[2, ] < 0.05]
# DFI
DFI <- lmm_vars("pbp_ss_rule_break")
DFI1 <- lapply(DFI, p_adjust, "fdr")
DFI1$sMRI <- NULL
DFI1 <- do.call(cbind, DFI1)
DFI1 <- DFI1[, DFI1[2, ] < 0.05]

# create data.frame for plot ---------------------------------------------------
effect_df <- function(data, peer, cv_r2) {
  sig_vars <- colnames(data)
  
  effect_peer <- data[3, ] %>% unlist
  effect_sex <- data[4, ] %>% unlist
  effect_age <- data[5, ] %>% unlist
  # effect_race <- data[6, ] %>% unlist
  effect_income <- data[7, ] %>% unlist
  effect_education <- data[8, ] %>% unlist
  
  # add parital R2 (cross-validation)
  effect_r2 <- cv_r2[match(colnames(data), names(cv_r2))]

  # effect_catgory <- rep(catgory, time = effect_catgory)
  effect_all <- c(effect_peer, effect_r2, effect_sex, effect_age, 
                  effect_income, effect_education)
  effect_sig <- rep(sig_vars, time = 6)
  effect_label <- rep(c(peer, "Cross-Validation", "Sex", "Age", 
                        "Parent_Income", "Parent_Education"), each = length(sig_vars))
  
  # data frame
  df <- data.frame("Effect" = effect_all, "Variables" = effect_sig, 
                   "Label" = effect_label) %>%
    mutate(Variables = factor(Variables, levels = sig_vars)) %>%
    mutate(Label = factor(Label, levels = c(peer, "Cross-Validation", "Sex", 
                                            "Age", "Parent_Income", 
                                            "Parent_Education"))) %>%
    mutate(X = as.numeric(Variables))
  

  return(df)
}

# significant results ----------------------------------------------------------
sig_PFI <- c(
  "Neurocognition - Picture Vocabulary",
  "CBCL - Anxious/Depressed problems",
  "CBCL - Withdrawn/Depressed problems",
  "CBCL - Somatic complaint problems",
  "CBCL - Social problems",
  "CBCL - Thought problems",
  "CBCL - Attention problems",
  "CBCL - Rule-Breaking problems",
  "CBCL - Aggressive problems",
  "CBCL - Internalizing problems",
  "CBCL - Externalizing problems",
  "CBCL - Total problems",
  "CBCL - Depression",
  "CBCL - Anxiety",
  "CBCL - Somatic",
  "CBCL - ADHD",
  "CBCL - Oppositional defiant problems",
  "CBCL - Conduct problems",
  "CBCL - Sluggish cognitive tempo",
  "CBCL - OCD",
  "CBCL - Stress problems",
  "Life Events (Parent) - Total events",
  "Life Events (Parent) - Bad events",
  "Life Events (Parent) - Bad affections",
  "Life Events (Parent) - Total affections",
  "Life Events - Total events",
  "Life Events - Bad affections",
  "Prodromal Psychosis - Total score",
  "Prodromal Psychosis - Distress score",
  "Impulsivity - Negative urgency",
  "Impulsivity - Lack of planning",
  "Impulsivity - Sensation seeking",
  "Impulsivity - Positive urgency",
  "Impulsivity - Lack of perseverance",
  "Inhibition and Reward-seeking - BIS sum score",
  "Inhibition and Reward-seeking - BAS reward responsiveness",
  "Inhibition and Reward-seeking - BAS drive",
  "Inhibition and Reward-seeking - BAS fun seeking",
  "Inhibition and Reward-seeking - BIS sum score (modified)",
  "Inhibition and Reward-seeking - BAS reward responsiveness (modified)",
  "Inhibition and Reward-seeking - BAS drive (modified)",
  "Adverse Peer Experience - Relational Victimization",
  "Adverse Peer Experience - Reputational Aggression",
  "Adverse Peer Experience - Reputational Victimization",
  "Adverse Peer Experience - Overt Aggression",
  "Adverse Peer Experience - Overt Victimization",
  "Adverse Peer Experience - Relational Aggression",
  "Inferior Parietal - Left (Volume)",
  "Paracentral - Left (Volume)",
  "Precentral – Left  (Volume)",
  "Superior Frontal - Left (Volume)",
  "Postcentral - Right (Volume)",
  "Insula - Right (Volume)",
  "Inferior Parietal - Left (Area)",
  "Paracentral - Left (Area)",
  "Insula - Left (Area)",
  "Postcentral - Right (Area)",
  "Posterior Cingulate - Right (Area)",
  "Frontal Pole - Right (Area)",
  "Insula - Right (Area)",
  "Putamen - Left (Subcortical Volume)",
  "Putamen - Right (Subcortical Volume)",
  "Pallidum - Right (Subcortical Volume)",
  "NAc - Right (Subcortical Volume)",
  "AN_CPN (RSFCs)",
  "AN_SN (RSFCs)",
  "CON_DAN (RSFCs)",
  "CPN_SHN (RSFCs)",
  "CPN_SMN (RSFCs)",
  "CPN_VN (RSFCs)",
  "DMN_FPN (RSFCs)",
  "DMN_RTN (RSFCs)",
  "DMN_Vtdc (Cortico-subcortical RSFCs)",
  "RTN_Amg (Cortico-subcortical RSFCs)",
  "RTN_Vtdc (Cortico-subcortical RSFCs)",
  "SMN_Cde (Cortico-subcortical RSFCs)"
)

sig_DFI <- c(
  "Neurocognition - Picture Vocabulary",
  "Neurocognition - Picture Sequence Memory",
  "Neurocognition - Oral reading recognition",
  "CBCL - Social problems",
  "CBCL - Thought problems",
  "CBCL - Attention problems",
  "CBCL - Rule-Breaking problems",
  "CBCL - Aggressive problems",
  "CBCL - Externalizing problems",
  "CBCL - Total problems",
  "CBCL - Depression",
  "CBCL - ADHD",
  "CBCL - Oppositional defiant problems",
  "CBCL - Conduct problems",
  "CBCL - Stress problems",
  "Subsyndromal Mania (Parent) - Total score",
  "Life Events (Parent) - Bad events",
  "Life Events (Parent) - Bad affections",
  "Life Events (Parent) - Mean affections",
  "Life Events - Total events",
  "Life Events - Bad events",
  "Life Events - Total affections",
  "Life Events - Good affections",
  "Life Events - Bad affections",
  "Prodromal Psychosis - Total score",
  "Prodromal Psychosis - Distress score",
  "Impulsivity - Negative urgency",
  "Impulsivity - Lack of planning",
  "Impulsivity - Sensation seeking",
  "Impulsivity - Positive urgency",
  "Impulsivity - Lack of perseverance",
  "Inhibition and Reward-seeking - BAS drive",
  "Inhibition and Reward-seeking - BAS fun seeking",
  "Inhibition and Reward-seeking - BAS drive (modified)",
  "Adverse Peer Experience - Relational Victimization",
  "Adverse Peer Experience - Reputational Aggression",
  "Adverse Peer Experience - Reputational Victimization",
  "Adverse Peer Experience - Overt Aggression",
  "Adverse Peer Experience - Overt Victimization",
  "Adverse Peer Experience - Relational Aggression",
  "Lateral Occipital - Left (Volume)",
  "Lateral Occipital - Left (Thickness)",
  "Lateral Occipital - Right (Thickness)",
  "DMN_DMN (RSFCs)",
  "DAN_DAN (RSFCs)",
  "SHN_SHN (RSFCs)",
  "CON_DMN (RSFCs)",
  "CON_VN (RSFCs)",
  "DMN_DAN (RSFCs)",
  "DAN_VAN (RSFCs)",
  "SHN_SMN (RSFCs)",
  "AN_Cde (Cortico-subcortical RSFCs)",
  "AN_Hip (Cortico-subcortical RSFCs)",
  "AN_NAc (Cortico-subcortical RSFCs)",
  "CON_Tha (Cortico-subcortical RSFCs)",
  "CON_Pt (Cortico-subcortical RSFCs)",
  "CON_Hip Cortico-subcortical (RSFCs)",
  "CON_Amg (Cortico-subcortical RSFCs)",
  "CPN_Crcx (Cortico-subcortical RSFCs)",
  "CPN_Cde (Cortico-subcortical RSFCs)",
  "CPN_Vtdc (Cortico-subcortical RSFCs)",
  "DMN_Pl (Cortico-subcortical RSFCs)",
  "DMN_Amg (Cortico-subcortical RSFCs)",
  "DMN_NAc (Cortico-subcortical RSFCs)",
  "DAN_Hip (Cortico-subcortical RSFCs)",
  "DAN_Amg (Cortico-subcortical RSFCs)",
  "DAN_Vtdc (Cortico-subcortical RSFCs)",
  "FPN_Crcx (Cortico-subcortical RSFCs)",
  "FPN_Tha (Cortico-subcortical RSFCs)",
  "FPN_Amg (Cortico-subcortical RSFCs)",
  "RTN_NAc (Cortico-subcortical RSFCs)",
  "SHN_Crcx (Cortico-subcortical RSFCs)",
  "SHN_Cde (Cortico-subcortical RSFCs)",
  "SHN_Pt (Cortico-subcortical RSFCs)",
  "SHN_Pl (Cortico-subcortical RSFCs)",
  "SHN_NAc (Cortico-subcortical RSFCs)",
  "SMN_Pl (Cortico-subcortical RSFCs)",
  "SMN_Hip (Cortico-subcortical RSFCs)",
  "SMN_Amg (Cortico-subcortical RSFCs)",
  "SN_Hip (Cortico-subcortical RSFCs)",
  "SN_NAc (Cortico-subcortical RSFCs)",
  "SN_Vtdc (Cortico-subcortical RSFCs)",
  "VAN_Cde (Cortico-subcortical RSFCs)",
  "CON_BS (Cortico-subcortical RSFCs)",
  "DAN_BS (Cortico-subcortical RSFCs)",
  "RTN_BS (Cortico-subcortical RSFCs)",
  "SHN_BS (Cortico-subcortical RSFCs)"
)

# ------------------------------------------------------------------------------
df_PFI <- effect_df(PFI1, "PFI", PFI_r2_test_sig) %>%
  mutate(Catgory = rep(rep(c("Behaviors", "Brain Structures", "RSFCs"), 
                           times = c(47, 17, 12)), time = 6)) %>%
  mutate(Catgory = factor(Catgory, levels = c("Behaviors", "Brain Structures", "RSFCs")))

df_DFI <- effect_df(DFI1, "DFI", DFI_r2_test_sig) %>%
  mutate(Catgory = rep(rep(c("Behaviors", "Brain Structures", "RSFCs"), 
                           times = c(40, 3, 44)), time = 6)) %>%
  mutate(Catgory = factor(Catgory, levels = c("Behaviors", "Brain Structures", "RSFCs")))

# 5 folds partial R2

# PFI_r2_df <- df_PFI[1:length(PFI_r2_test_sig), ]
# PFI_r2_df <- mutate(PFI_r2_df, Label = "Cross_Valiadtion") %>%
#   mutate(Effect = PFI_r2_test_sig)
# 
# DFI_r2_df <- df_DFI[1:length(DFI_r2_test_sig), ]
# DFI_r2_df <- mutate(DFI_r2_df, Label = "Cross_Valiadtion") %>%
#   mutate(Effect = DFI_r2_test_sig)
# 
# # combine data
# df_PFI <- rbind(df_PFI, PFI_r2_df)
# df_DFI <- rbind(df_DFI, DFI_r2_df)

# plot"Semi-Partial R^2"
effect_plot <- function(data, colors, sigvars) {
  
  plot_function <- function(df, colors, sigvars, show_legend) {
    p <- ggplot(df, aes(X, Effect, color = Label, shape = Label)) +
      geom_point(size = 2.5, alpha = 0.7, show.legend = show_legend) +
      geom_line(size = 1, alpha = 0.7, show.legend = show_legend) +
      scale_x_continuous(breaks = unique(df$X), labels = sigvars) +
      scale_color_manual(values = c(colors, "#E7298A", "#FEE08B", "#BF812D", 
                                    "#542788", "#FB6A4A", "#A50F15")) +
      labs(x = NULL, y = "Effect size") +
      # facet_grid(~Catgory, scales = "free", space = "free") +
      # facet_wrap(vars(Catgory), scales = "free", nrow = 3) +
      theme_classic() + 
      theme(
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15, face = "bold"),
        plot.tag = element_text(size = 30, face = "bold"),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.position = "top",
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_blank(),
        legend.box.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        panel.spacing = unit(1.5, "lines")
      )
    return(p)
  }
  # behavior data
  behaviors_data <- filter(data, Catgory == "Behaviors")
  behaivors_index <- which(data$Catgory == "Behaviors" & data$Label == behaviors_data$Label[1])
  # brain structures
  sMRI_data <- filter(data, Catgory == "Brain Structures")
  sMRI_index <- which(data$Catgory == "Brain Structures" & data$Label == behaviors_data$Label[1])
  # RSFCs
  fc_data <- filter(data, Catgory == "RSFCs")
  fc_index <- which(data$Catgory == "RSFCs" & data$Label == behaviors_data$Label[1])
  
  p1 <- plot_function(behaviors_data, colors, sigvars[behaivors_index], TRUE)
  p2 <- plot_function(sMRI_data, colors, sigvars[sMRI_index], FALSE)
  p3 <- plot_function(fc_data, colors, sigvars[fc_index], FALSE)
  p1 + p2 + p3 + 
    plot_layout(nrow = 3)
}
effect_plot(df_PFI, "#3B7346", sig_PFI)
ggsave("Effect_PFI.svg", width = 10, height = 16)
effect_plot(df_DFI, "#4286f4", sig_DFI)
ggsave("Effect_DFI.svg", width = 10, height = 16)
