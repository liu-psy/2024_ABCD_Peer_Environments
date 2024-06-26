# 5 folds cross-validation
library(caret)
library(lmerTest)
library(patchwork)
library(tidyverse)

setwd("H:/ABCD/Release4.0/Package_1194636/results/peer_environments")
# load("results.RData")

# reduced formulas (without PFI/DFI)
reduce_behavior <- "~ sex + interview_age + race_ethnicity + income_parent + 
  edu_parent + rel_family_id + (1|site_id_l)"
reduce_fc <- "~ sex + interview_age + race_ethnicity + income_parent + 
  edu_parent + rsfmri_c_ngd_meanmotion + rel_family_id + 
  (1|mri_info_deviceserialnumber)"
reduce_sMRI <- "~ sex + interview_age + race_ethnicity + income_parent + 
  edu_parent + smri_vol_scs_intracranialv + rel_family_id + 
  (1|mri_info_deviceserialnumber)"

reduce_behavior <- paste(behavior_year2, reduce_behavior)
reduce_fc <- paste(fc_network, reduce_fc)
reduce_sMRI <- paste(c(smri_vol_cortical, smri_area, smri_thick, smri_vol_subcortical), 
                     reduce_sMRI)
reduce_all <- c(reduce_behavior, reduce_sMRI, reduce_fc)

# 5 folds data
set.seed(123)
behavior_5folds <- createFolds(1:nrow(abcd_year2_behavior), 5, FALSE)
smri_5folds <- createFolds(1:nrow(abcd_year2_smri), 5, FALSE)
rsfmri_5folds <- createFolds(1:nrow(abcd_year2_rsfmri), 5, FALSE)

# get partial R2
get_r2 <- function(formulas, test_fold_index) {
  # full formula of PFI and DFI
  PFI_full <- str_replace(formulas, "~", paste("~", "pbp_ss_prosocial_peers", "+"))
  DFI_full <- str_replace(formulas, "~", paste("~", "pbp_ss_rule_break", "+"))
  
  # train set
  behavior_train <- abcd_year2_behavior[behavior_5folds != test_fold_index, ]
  rsfmri_train <- abcd_year2_rsfmri[rsfmri_5folds != test_fold_index, ]
  smri_train <- abcd_year2_smri[smri_5folds != test_fold_index, ]
  
  # test set
  behavior_test <- abcd_year2_behavior[behavior_5folds == test_fold_index, ]
  rsfmri_test <- abcd_year2_rsfmri[rsfmri_5folds == test_fold_index, ]
  smri_test <- abcd_year2_smri[smri_5folds == test_fold_index, ]
  
  # response variable
  y <- str_split(formulas, " ~ ", simplify = TRUE)[1]
  
  df_train <- ifelse(str_detect(formulas, "site_id_l"), "behavior_train", ifelse(
    str_detect(formulas, "rsfmri_c_ngd_meanmotion"), "rsfmri_train", "smri_train"
  ))
  
  df_test <- ifelse(str_detect(formulas, "site_id_l"), "behavior_test", ifelse(
    str_detect(formulas, "rsfmri_c_ngd_meanmotion"), "rsfmri_test", "smri_test"
  ))
  
  cat("Fold", test_fold_index, ": ", y, 
      df_train, " = ", nrow(get(df_train)), ";",
      df_test, " = ", nrow(get(df_test)), 
      "\n")
  # reduced model (in train set)
  reduce_train <- lmer(as.formula(formulas), data = get(df_train))
  # full model (in train set)
  PFI_full_train <- lmer(as.formula(PFI_full), data = get(df_train))
  DFI_full_train <- lmer(as.formula(DFI_full), data = get(df_train))
  
  # RSS (in train set)
  reduce_rss_train <- sum(summary(reduce_train)$residuals^2)
  PFI_full_rss_train <- sum(summary(PFI_full_train)$residuals^2)
  DFI_full_rss_train <- sum(summary(DFI_full_train)$residuals^2)
  
  # RSS (in test set)
  get_rss_test <- function(models, test_data) {
    predict_value <- predict(models, newdata = get(test_data))
    test_residuals <- get(test_data)[[y]] - predict_value
    rss_test <- sum(test_residuals^2)
    return(rss_test)
  }
  reduce_rss_test <- get_rss_test(reduce_train, df_test)
  PFI_full_rss_test <- get_rss_test(PFI_full_train, df_test)
  DFI_full_rss_test <- get_rss_test(DFI_full_train, df_test)
  
  # partial R2 (PFI)
  PFI_train_r2 <- (reduce_rss_train - PFI_full_rss_train) / reduce_rss_train
  PFI_test_r2 <- (reduce_rss_test - PFI_full_rss_test) / reduce_rss_test
  # partial R2 (DFI)
  DFI_train_r2 <- (reduce_rss_train - DFI_full_rss_train) / reduce_rss_train
  DFI_test_r2 <- (reduce_rss_test - DFI_full_rss_test) / reduce_rss_test
  
  return(c(reduce_rss_train, PFI_full_rss_train, DFI_full_rss_train, 
           reduce_rss_test, PFI_full_rss_test, DFI_full_rss_test, PFI_train_r2,
           PFI_test_r2, DFI_train_r2, DFI_test_r2))
}
CV_fold1 <- sapply(reduce_all, get_r2, 1)
CV_fold2 <- sapply(reduce_all, get_r2, 2)
CV_fold3 <- sapply(reduce_all, get_r2, 3)
CV_fold4 <- sapply(reduce_all, get_r2, 4)
CV_fold5 <- sapply(reduce_all, get_r2, 5)

# get mean partial R2 in test set
mean_r2 <- function(row_index) {
  df <- rbind(CV_fold1[row_index,], CV_fold2[row_index,], CV_fold3[row_index,],
              CV_fold4[row_index,], CV_fold5[row_index,])
  colnames(df) <- c(behavior_year2, smri_vol_cortical, smri_area, 
                    smri_thick, smri_vol_subcortical, fc_network)
  mean_r2 <- colMeans(df)
  return(mean_r2)
}
PFI_r2_test <- mean_r2(8)
DFI_r2_test <- mean_r2(10) 

# get raw effect sizes in significant results ----------------------------------
p_adjust <- function(x, p_correction) {
  x[2, ] <- p.adjust(x[2, ], method = p_correction)
  # x <- x[, x[2, ] < 0.05]
  return(x)
}

# PFI (significant results)
PFI1 <- lapply(PFI, p_adjust, "fdr")
PFI1$sMRI <- NULL
PFI1 <- do.call(cbind, PFI1)
PFI1 <- PFI1[, PFI1[2, ] < 0.05]
# DFI (significant results)
DFI1 <- lapply(DFI, p_adjust, "fdr")
DFI1$sMRI <- NULL
DFI1 <- do.call(cbind, DFI1)
DFI1 <- DFI1[, DFI1[2, ] < 0.05]

# significant variables
PFI1_sig <- colnames(PFI1)
DFI1_sig <- colnames(DFI1)

PFI_r2_test_sig <- PFI_r2_test[names(PFI_r2_test) %in% PFI1_sig]
DFI_r2_test_sig <- DFI_r2_test[names(DFI_r2_test) %in% DFI1_sig]

# create data.frame for plot
effect_df <- function(data, peer, cv_r2) {
  sig_vars <- colnames(data)
  
  effect_peer <- data[3, ] %>% unlist

  # add parital R2 (cross-validation)
  effect_r2 <- cv_r2[match(colnames(data), names(cv_r2))]
  
  # effect_catgory <- rep(catgory, time = effect_catgory)
  effect_all <- c(effect_peer, effect_r2)
  effect_sig <- rep(sig_vars, time = 2)
  effect_label <- rep(c(peer, "Cross-Validation"), each = length(sig_vars))
  
  # data frame
  df <- data.frame("Effect" = effect_all, "Variables" = effect_sig, 
                   "Label" = effect_label) %>%
    mutate(Variables = factor(Variables, levels = sig_vars)) %>%
    mutate(Label = factor(Label, levels = c(peer, "Cross-Validation"))) %>%
    mutate(X = as.numeric(Variables))
  
  return(df)
}

df_PFI <- effect_df(PFI1, "PFI", PFI_r2_test_sig) %>%
  mutate(Catgory = rep(rep(c("Behaviors", "Brain Structures", "RSFCs"), 
                           times = c(47, 17, 12)), time = 2)) %>%
  mutate(Catgory = factor(Catgory, levels = c("Behaviors", "Brain Structures", "RSFCs")))

df_DFI <- effect_df(DFI1, "DFI", DFI_r2_test_sig) %>%
  mutate(Catgory = rep(rep(c("Behaviors", "Brain Structures", "RSFCs"), 
                           times = c(40, 3, 44)), time = 2)) %>%
  mutate(Catgory = factor(Catgory, levels = c("Behaviors", "Brain Structures", "RSFCs")))

# correlations between raw effect sizes and CV-based effect sizes
cor.test(df_PFI$Effect[1:ncol(PFI1)], df_PFI$Effect[(1 + ncol(PFI1)):(ncol(PFI1) * 2)])
cor.test(df_DFI$Effect[1:ncol(DFI1)], df_DFI$Effect[(1 + ncol(DFI1)):(ncol(DFI1) * 2)])

# significant results ----------------------------------------------------------
# PFI
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
# DFI
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

# sFigure 12-13 ----------------------------------------------------------------
effect_plot <- function(data, colors, sigvars) {
  
  plot_function <- function(df, colors, sigvars, show_legend) {
    p <- ggplot(df, aes(X, Effect, color = Label, shape = Label)) +
      geom_point(size = 2.5, alpha = 0.7, show.legend = show_legend) +
      geom_line(size = 1, alpha = 0.7, show.legend = show_legend) +
      scale_x_continuous(breaks = unique(df$X), labels = sigvars) +
      scale_color_manual(values = c(colors, "#E7298A")) +
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
  behaivors_index <- which(data$Catgory == "Behaviors" & data$Label != "Cross-Validation")
  # brain strcutre
  sMRI_data <- filter(data, Catgory == "Brain Structures")
  sMRI_index <- which(data$Catgory == "Brain Structures" & data$Label != "Cross-Validation")
  # RSFCs
  fc_data <- filter(data, Catgory == "RSFCs")
  fc_index <- which(data$Catgory == "RSFCs" & data$Label != "Cross-Validation")
  
  
  p1 <- plot_function(behaviors_data, colors, sigvars[behaivors_index], TRUE)
  p2 <- plot_function(sMRI_data, colors, sigvars[sMRI_index], FALSE)
  p3 <- plot_function(fc_data, colors, sigvars[fc_index], FALSE)
  p1 + p2 + p3 + 
    plot_layout(nrow = 3)
}
effect_plot(df_PFI, "#3B7346", sig_PFI)
ggsave("sssssEffect_PFI.svg", width = 10, height = 12)
effect_plot(df_DFI, "#4286f4", sig_DFI)
ggsave("ssssssEffect_DFI.svg", width = 10, height = 12)