# PFI and DFI all included in LMMs
library(circlize)
library(ComplexHeatmap)
library(ggseg)
library(patchwork)
library(r2glmm)
library(lavaan)
library(lmerTest)
library(tidyverse)

# set working directory
setwd("H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/supplementary_analysis/")

# load data
abcd <- read_csv("H:/ABCD/Relsease4.0/Package_1194636/abcd.csv")
base_information <- read_csv("H:/ABCD/Relsease4.0/Package_1194636/base_information.csv")
abcd_year2 <- filter(abcd, eventname == "2_year_follow_up_y_arm_1")
abcd_year3 <- filter(abcd, eventname == "3_year_follow_up_y_arm_1")

# Select variables -------------------------------------------------------------
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
           "income_parent", "edu_parent", "rel_family_id")
# peer environments
peers <- c("pbp_ss_prosocial_peers", "pbp_ss_rule_break")

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
fmri_network_subcortical <- names(select(abcd_year2, contains("rsfmri_cor_ngd")))
fc_network <- c(fmri_networks, fmri_network_subcortical)

# exclude missing value and MRI QC
abcd_year2_behavior <- select(abcd_year2, all_of(basic), all_of(peers),
    all_of(behavior_year2), site_id_l) %>%
  filter(complete.cases(.), !duplicated(.))

abcd_year2_smri <- select(abcd_year2, all_of(basic), all_of(peers), all_of(smri),
    mri_info_deviceserialnumber, smri_vol_scs_intracranialv, imgincl_t1w_include) %>%
  filter(imgincl_t1w_include == 1, complete.cases(.), !duplicated(.))

abcd_year2_rsfmri <- select(abcd_year2, all_of(basic), all_of(peers),
    all_of(fc_network), mri_info_deviceserialnumber, rsfmri_c_ngd_meanmotion,
    imgincl_rsfmri_include) %>%
  filter(imgincl_rsfmri_include == 1, complete.cases(.), !duplicated(.))

abcd_year3_behavior <- select(abcd_year3, all_of(basic), all_of(peers),
    all_of(behavior_year3), site_id_l) %>%
  filter(complete.cases(.), !duplicated(.))

# longitudinal data
abcd_longitudinal_year2 <- abcd_year2_behavior[abcd_year2_behavior$subjectkey %in% abcd_year3_behavior$subjectkey, ]
abcd_longitudinal_year3 <- abcd_year3_behavior[abcd_year3_behavior$subjectkey %in% abcd_year2_behavior$subjectkey, ]

# 1. Association Analysis ------------------------------------------------------
# PFI and DFI in the models
lmm_vars <- function() {
  # behavior
  lmm_behavior <- function(x, df) {
    cat("\n", x, "   Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ pbp_ss_prosocial_peers + pbp_ss_rule_break + sex + 
                  interview_age + race_ethnicity + income_parent + edu_parent + 
                  (1|site_id_l/rel_family_id), data = df)
    fit_summary <- summary(fit)
    t_PFI <- fit_summary$coefficients[2, 4]
    p_PFI <- fit_summary$coefficients[2, 5]
    effect_PFI <- r2beta(fit, method = "nsj")$Rsq[2]
    
    t_DFI <- fit_summary$coefficients[3, 4]
    p_DFI <- fit_summary$coefficients[3, 5]
    effect_DFI <- r2beta(fit, method = "nsj")$Rsq[3]
    
    return(c(t_PFI, p_PFI, effect_PFI, t_DFI, p_DFI, effect_DFI))
  }
  
  # rs-fMRI
  lmm_rsfmri <- function(x, df) {
    cat("\n", x, "Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ pbp_ss_prosocial_peers + pbp_ss_rule_break + sex + 
                  interview_age + race_ethnicity + income_parent + edu_parent + 
                  rsfmri_c_ngd_meanmotion +
                  (1|mri_info_deviceserialnumber/rel_family_id), data = df)
    fit_summary <- summary(fit)
    # 3 df; 4 t-value; 5 p-value
    t_PFI <- fit_summary$coefficients[2, 4]
    p_PFI <- fit_summary$coefficients[2, 5]
    effect_PFI <- r2beta(fit, method = "nsj")$Rsq[2]
    
    t_DFI <- fit_summary$coefficients[3, 4]
    p_DFI <- fit_summary$coefficients[3, 5]
    effect_DFI <- r2beta(fit, method = "nsj")$Rsq[3]
    
    return(c(t_PFI, p_PFI, effect_PFI, t_DFI, p_DFI, effect_DFI))
  }
  
  # sMRI
  lmm_smri <- function(x, df) {
    cat("\n", x, "Particpant:", nrow(df))
    fit <- lmer(df[[x]] ~ pbp_ss_prosocial_peers + pbp_ss_rule_break + sex + 
                  interview_age + race_ethnicity + income_parent + edu_parent + 
                  smri_vol_scs_intracranialv + 
                  (1|mri_info_deviceserialnumber/rel_family_id), data = df)
    fit_summary <- summary(fit)
    # 3 df; 4 t-value; 5 p-value
    t_PFI <- fit_summary$coefficients[2, 4]
    p_PFI <- fit_summary$coefficients[2, 5]
    effect_PFI <- r2beta(fit, method = "nsj")$Rsq[2]
    
    t_DFI <- fit_summary$coefficients[3, 4]
    p_DFI <- fit_summary$coefficients[3, 5]
    effect_DFI <- r2beta(fit, method = "nsj")$Rsq[3]
    
    return(c(t_PFI, p_PFI, effect_PFI, t_DFI, p_DFI, effect_DFI))
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
  results <- list(result_behavior, result_networks, result_network_subcortical,
                  result_thick, result_area, result_vol_cortical, result_vol_subcortical,
                  result_smri)
  names(results) <- c("behavior", "networks",  "network_sub", "thick", "area",
                      "vol_cortical", "vol_subcortical", "sMRI")
  
  # add rownames
  modify <- function(x) {
    rownames(x) <- c("t_PFI", "p_PFI", "effect_size_PFI", "t_DFI", "p_DFI", 
                     "effect_size_DFI")
    x[1, ] <- round(unlist(x[1, ]), 2)
    x[3, ] <- round(unlist(x[3, ]), 4)
    x[4, ] <- round(unlist(x[4, ]), 2)
    x[6, ] <- round(unlist(x[6, ]), 4)
    return(x)
  }
  results <- lapply(results, modify)
  
  return(results)
}

both_results <- lmm_vars()

# select
select_row <- function(df, row_index) {
  df <- df[row_index, ]
  return(df)
}
PFI <- lapply(both_results, select_row, 1:3)
DFI <- lapply(both_results, select_row, 4:6)

# FDR corrections
p_adjust <- function(x, p_correction) {
  x[2, ] <- p.adjust(x[2, ], method = p_correction)
  x <- x[, x[2, ] < 0.05]
  return(x)
}
PFI1 <- lapply(PFI, p_adjust, "fdr")
DFI1 <- lapply(DFI, p_adjust, "fdr")

# sFiugre 12 -------------------------------------------------------------------
# peer environments
peers <- c("PFI", "DFI")

# new behavior names
behavior_names <- c(
  "Neurocognition - Picture Vocabulary",
  "Neurocognition - Picture Sequence Memory",
  "Neurocognition - Oral reading recognition",
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
  "Subsyndromal Mania (Parent) - Sum score",
  "Life Events (Parent) - Total events",
  "Life Events (Parent) - Good events",
  "Life Events (Parent) - Bad events",
  "Life Events (Parent) - Good affections",
  "Life Events (Parent) - Bad affections",
  "Life Events (Parent) - Mean affections",
  "Life Events (Parent) - Total affections",
  "Life Events - Total events",
  "Life Events - Good events",
  "Life Events - Bad events",
  "Life Events - Good affections",
  "Life Events - Bad affections",
  "Life Events - Mean affections",
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
  "Adverse Peer Experiences - Relational Victimization",
  "Adverse Peer Experiences - Reputational Aggression",
  "Adverse Peer Experiences - Reputational Victimization",
  "Adverse Peer Experiences - Overt Aggression",
  "Adverse Peer Experiences - Overt Victimization",
  "Adverse Peer Experiences - Relational Aggression"
)

mat_aossciations <- function(x, peer_type) {
  behavior_mat <- data.frame(
    "t" = unlist(x[["behavior"]][1, ]),
    "p" = unlist(x[["behavior"]][2, ])
  )
  behavior_mat$group <- c(
    rep("Neurocognition", 3),
    rep("CBCL", 20),
    "Subsyndromal\nMania (Parent)",
    rep("Life Events\n(Parent)", 7),
    rep("Life Events", 6),
    rep("Prodromal\nPsychosis", 2),
    rep("Impulsivity", 5),
    rep("Inhibition and\nReward-seeking", 7),
    rep("Adverse Peer\nExperiences", 6)
  )
  behavior_mat$group <- factor(
    behavior_mat$group,
    levels = c("CBCL",  "Life Events", "Life Events\n(Parent)",
               "Prodromal\nPsychosis", "Subsyndromal\nMania (Parent)", 
               "Inhibition and\nReward-seeking", "Impulsivity", 
               "Adverse Peer\nExperiences", "Neurocognition")
  )
  
  behavior_mat$pfdr <- p.adjust(behavior_mat$p, method = "fdr")
  behavior_mat$sig <- ifelse(behavior_mat$pfdr < 0.001, "***", ifelse(
    behavior_mat$pfdr < 0.01, "**", ifelse(
      behavior_mat$pfdr < 0.05, "*", ""
    )
  )
  )
  
  behavior_mat <- arrange(behavior_mat, group)
  behavior_mat$peers <- peer_type
  return(behavior_mat)
}

behavior_PFI <- mat_aossciations(PFI, "PFI")
behavior_DFI <- mat_aossciations(DFI, "DFI")
behavior_friends <- rbind(behavior_PFI, behavior_DFI)

behavior_friends$seq <- rep(seq(behavior_year2), time = length(peers))
behavior_friends$peers <- factor(behavior_friends$peers, levels = peers)

ggplot(behavior_friends, aes(group, t, color = peers)) +
  geom_point(size = 8, position = position_jitterdodge(), shape = 17,
             alpha = 0.8) +
  annotate("text", x = 9.2, y = 4, size = 7,
           label = expression(paste(bold(bolditalic("p")[fdr]), " ", bold("< 0.05")))) +
  scale_y_continuous(breaks = seq(-25, 20, 5)) +
  scale_color_manual(values = c("#3B7346", "#4286f4")) +
  labs(x = NULL, y = expression(paste(bolditalic("t"), " ", bold(value)))) +
  geom_hline(yintercept = 0, linewidth = 1) +
  geom_hline(yintercept = 2.27, linewidth = 1, linetype = 2) +
  geom_hline(yintercept = -2.27, linewidth = 1, linetype = 2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, hjust = 1, angle = 15, face = "bold"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.tag = element_text(size = 30, face = "bold"),
    legend.key.height = unit(1, 'cm'),
    legend.key.width = unit(1, 'cm'),
    legend.position = "top",
    legend.text = element_text(size = 35),
    legend.title = element_blank(),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(1.5, "lines")
  )  +
  guides(color = guide_legend(ncol = 2))
ggsave("sFigure12.svg", width = 20, height = 12, bg = "transparent")

# sFigure 13 -------------------------------------------------------------------
data(dk)
dk$data$region

extract_region <- function(data, modality) {
  ABCD_left_hemi <- str_split(smri_thick, "_cdk_", simplify = TRUE)[1:34, 2]
  ABCD_hemi <- str_extract(ABCD_left_hemi, "..$")
  ABCD_regions <- str_split(ABCD_left_hemi, ABCD_hemi, simplify = TRUE)[, 1]
  
  dk_regions <- c("bankssts", "caudal anterior cingulate",
                  "caudal middle frontal", "cuneus", "entorhinal", "entorhinal",
                  "inferior parietal", "inferior temporal", "isthmus cingulate",
                  "lateral occipital", "lateral orbitofrontal", "lingual",
                  "medial orbitofrontal", "middle temporal", "parahippocampal",
                  "paracentral", "pars opercularis", "pars orbitalis",
                  "pars triangularis", "pericalcarine", "postcentral",
                  "posterior cingulate", "precentral", "precuneus",
                  "rostral anterior cingulate", "rostral middle frontal", "superior frontal",
                  "superior parietal", "superior temporal", "supramarginal",
                  "frontal pole", "temporal pole", "transverse temporal", "insula")
  replace_table <- data.frame(ABCD_regions, dk_regions)
  
  ROI <- str_split(colnames(data[[modality]]), pattern = "_cdk_", simplify = TRUE)[, 2]
  hemi <- str_extract(ROI, "..$")
  regions <- str_split(ROI, hemi, simplify = TRUE)[, 1]
  
  brain_region <- tibble(
    "hemi" = ifelse(hemi == "lh", hemi <- "left", hemi <- "right"),
    "region" = regions,
    "t" = unlist(data[[modality]][1, ], use.names = FALSE)
  )
  
  for (i in 1:nrow(brain_region)) {
    brain_region$region[i] <- dk_regions[which(brain_region$region[i] == ABCD_regions)]
  }
  
  return(brain_region)
}

# PFI (volume)
PFI_vol <- extract_region(PFI1, "vol_cortical")
p1 <- ggplot(PFI_vol) +
  geom_brain(atlas = dk, aes(fill = t), size = 0.5, show.legend = FALSE,
             position = position_brain(. ~ side + hemi)) +
  scale_fill_gradient2(limits = c(-4.5, 4.5), low = "#9ECAE1", high = "#EF6548",
                       mid = "white", midpoint = 0, na.value = "white", space = "Lab", name = "t") +
  labs(tag = "a", title = "PFI (volume)") +
  theme_void() +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5, face = "bold"),
    plot.tag = element_text(size = 40, hjust = 0.5, face = "bold")
  )
# PFI (area)
PFI_area <- extract_region(PFI1, "area")
p2 <- ggplot(PFI_area) +
  geom_brain(atlas = dk, aes(fill = t), size = 0.5, show.legend = FALSE,
             position = position_brain(. ~ side + hemi)) +
  scale_fill_gradient2(limits = c(-4.5, 4.5), low = "#9ECAE1", high = "#EF6548",
                       mid = "white", midpoint = 0, na.value = "white", space = "Lab", name = "t") +
  labs(tag = "b", title = "PFI (area)") +
  theme_void() +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5, face = "bold"),
    plot.tag = element_text(size = 40, hjust = 0.5, face = "bold")
  )

# DFI (volume)
DFI_vol <- data.frame("hemi" = "left", region = "lateral occipital", "t" = -3.45)
p3 <- ggplot(DFI_vol) +
  geom_brain(atlas = dk, aes(fill = t), size = 0.5, show.legend = TRUE,
             position = position_brain(. ~ side + hemi)) +
  scale_fill_gradient2(limits = c(-4, 4), low = "#9ECAE1", high = "#EF6548",
                       mid = "white", midpoint = 0, na.value = "white", space = "Lab", name = "t") +
  labs(tag = "c", title = "DFI (volume)") +
  theme_void() +
  theme(plot.title = element_text(size = 26, hjust = 0.5)) + 
  theme_void() +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5, face = "bold"),
    plot.tag = element_text(size = 40, hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30),
    legend.key.height = unit(1.5, 'cm'),
    legend.key.width = unit(2, 'cm'),
  ) +
  guides(
    fill = guide_colorbar(
      title = expression(paste(bolditalic("t"), " ", bold(value))),
      title.position = "left",
      title.hjust	= 0.5,
      label.position = "bottom",
      ticks = FALSE
    )
  )

p1 + p2 + p3 + 
  plot_layout(nrow = 3)
ggsave("sFigure13.svg", width = 20, height = 15, bg = "transparent")

# sFigure 14 -------------------------------------------------------------------
network_edge <- function(data) {
  network <- c("ad", "cgc", "ca", "dt", "dla", "fo", "rspltp", "smh", "smm",
               "sa", "vta", "vs")
  network_label <- c("AN", "CON", "CPN", "DMN", "DAN", "FPN", "RTN", "SHN",
                     "SMN", "SN", "VAN", "VN")
  
  network_df <- data.frame(
    "from" = network,
    "to" = network
  )
  
  network_from <- str_split(colnames(data), "_ngd_", simplify = TRUE)[, 2]
  network_to <- str_split(colnames(data), "_ngd_", simplify = TRUE)[, 3]
  edges <- data.frame("from" = network_from, "to" = network_to)
  network_df <- rbind(edges, network_df)
  network_df$from <- factor(network_df$from, levels = network, labels = network_label)
  network_df$to <- factor(network_df$to, levels = network, labels = network_label)
  network_df$t_value <- c(unlist(data[1, ]), rep(NA, length(network_label)))
  
  return(network_df)
}
PFI_edges <- network_edge(PFI1$networks)
DFI_edges <- network_edge(DFI1$networks)

cortico_subcortical <- function(data) {
  network_subcortical <- c("au", "cerc", "copa", "df", "dsa", "fopa", "rst",
                           "smh", "smm", "sa", "vta", "vs", "crcx", "thp", "cde", "pt", "pl", "hp",
                           "ag", "aa", "vtdc", "bs")
  
  network_subcortical_label <- c("AN", "CON", "CPN", "DMN", "DAN", "FPN", "RTN",
                                 "SHN", "SMN", "SN", "VAN", "VN", "Crcx",  "Tha", "Cde", "Pt", "Pl", "Hip", 
                                 "Amg", "NAc", "Vtdc", "BS")
  
  template <- c("au", "cerc", "copa", "df", "dsa", "fopa", "rst", "smh", "smm", 
                "sa", "vta", "vs", "crcx", "thp", "cde", "pt", "pl", "hp", "ag", "aa", 
                "vtdc", "bs")
  
  network_df <- data.frame(
    "from" = template,
    "to" = template
  )
  
  network_from <- str_split(colnames(data), "_scs_", simplify = TRUE)
  network_from <- str_split(network_from[, 1], "_ngd_", simplify = TRUE)[, 2]
  network_to <- str_split(colnames(data), "_scs_", simplify = TRUE)[, 2]
  edges <- data.frame("from" = network_from, "to" = network_to)
  network_df <- rbind(edges, network_df)
  network_df$from <- factor(network_df$from, levels = network_subcortical, 
                            labels = network_subcortical_label)
  network_df$to <- factor(network_df$to, levels = network_subcortical, 
                          labels = network_subcortical_label)
  network_df$t_value <- c(unlist(data[1, ]), rep(NA, length(template)))
  
  return(network_df)
}

# DFI (cortico-subcortical RSFCs)
cs_DFI <- cortico_subcortical(DFI1$network_sub)

all_tvalue <- c(PFI_edges$t_value, DFI_edges$t_value, cs_DFI$t_value)
col_fun <- colorRamp2(
  c(min(all_tvalue, na.rm = TRUE), 0, max(all_tvalue, na.rm = TRUE)),
  c("#377EB8", "white", "#E41A1C")
)

# layouts
svg("sFigure14.svg", width = 16, height = 8)
lays <- layout(mat = matrix(c(1,1,2,2,3,3), nrow = 1, byrow = TRUE))
layout.show(lays)

# PFI
par(cex = 1.8)
circos.par(gap.after = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), start.degree = 44)
chordDiagram(
  PFI_edges[, c(1,2)],
  grid.col = "grey",
  annotationTrack = c("name", "grid"),
  annotationTrackHeight = c(0.01, 0.02),
  transparency = 0.7,
  col = col_fun(PFI_edges$t_value),
  order = c("AN", "SMN", "SHN", "VN", "DMN", "CON", "DAN", "CPN", "VAN", "FPN",
            "RTN", "SN"),
  scale = TRUE
)
title(main = "PFI", cex.main = 1.3, line = -2.5, font.main = 1)
circos.clear()

# DFI
circos.par(gap.after = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), start.degree = 44)
chordDiagram(
  DFI_edges[, c(1,2)],
  grid.col = "grey",
  annotationTrack = c("name", "grid"),
  annotationTrackHeight = c(0.01, 0.02),
  transparency = 0.7,
  col = col_fun(DFI_edges$t_value),
  order = c("AN", "SMN", "SHN", "VN", "DMN", "CON", "DAN", "CPN", "VAN", "FPN",
            "RTN", "SN"),
  scale = TRUE
)
title(main = "DFI", cex.main = 1.3, line = -2.5, font.main = 1)
circos.clear()

# DFI (cortico-subcortical RSFCs)
par(cex = 1.8)
circos.par(gap.after = c(rep(2, 11), 30, rep(2, 9), 30), start.degree = 263)
chordDiagram(
  cs_DFI[, c(1, 2)],
  grid.col = "grey",
  annotationTrack = c("name", "grid"),
  annotationTrackHeight = c(0.01, 0.02),
  transparency = 0.6,
  col = col_fun(cs_DFI$t_value),
  order = c("CON", "CPN", "DMN", "RTN", "VAN", "SN", "DAN", "FPN", "SHN",
            "SMN", "AN", "VN", "Crcx",  "Tha", "Hip", "Amg", "Pl", "Pt", "Cde", "NAc",
            "Vtdc", "BS"),
  scale = TRUE
)
title("DFI \n (cortico-subcortical RSFCs)", cex.main = 1.3, line = -3.5,
      font.main = 1)
circos.clear()

# add legend (run 2 times)
lgd_links <- Legend(
  at = c(-5.5, 0, 5.5),
  col_fun = col_fun,
  legend_width = unit(3, "in"), #
  grid_height = unit(10, "mm"),
  labels_gp = gpar(cex = 2, fontface = "bold"),
  legend_gp = gpar(alpha = 0.2, cex = 1),
  direction = "horizontal",
  title = expression(paste(italic("t"), " value")),
  title_position = "topcenter",
  title_gp = gpar(fontsize = 30, fontface = "bold"),
  title_gap = unit(2, "mm"),
  border = NA) %>%
  draw(x = unit(203.3, "mm"), y = unit(25, "mm"))
dev.off()