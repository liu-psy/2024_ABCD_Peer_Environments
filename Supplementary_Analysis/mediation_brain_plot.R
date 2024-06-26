# Theme: visualization mediation analysis
# Mediator = Brain variables
# sFigure 13-14
library(ggpubr)
library(reshape2)
library(tidyverse)

# load results from
setwd("H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/supplementary_analysis")
load("mediations_brain/mediation_brain_plot.RData")

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

mediation_fdr <- function(strs, pattern, peers, replace_index) {
  
  # read mediation results
  read_csvs <- function(files) {
    # set row names
    paths_p <- paste0(c("PathA", "PathB", "PathCC", "PathC", "PathAB"), "_p")
    paths_t <- paste0(c("PathA", "PathB", "PathCC", "PathC", "PathAB"), "_beta")
    CI <- paste0(rep(c("A", "B", "CC", "C", "AB"), 2), rep(1:2, each = 5))
    paths <- c(paths_p, paths_t, CI)
    
    # read files
    csv_file <- paste0("mediations_brain/", files, ".csv")
    result <- read.csv(csv_file, header = FALSE)
    
    colnames(result) <- paths
    rownames(result) <- paste0(files, "--", behavior_names)
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
  
  df$Mediator <- str_split(rownames(df), "_", simplify = TRUE)[, 1]
  df$Associations <- str_split(rownames(df), paste0(df$Mediator, "_"), 
                               simplify = TRUE)[, 2]
  df <- select(df, Mediator, Associations, everything())
  return(df)
}
# mediation - PFI 10 regional volume (13 significant results)
PFI_sig_volumes <- c(colnames(PFI1$vol_cortical), colnames(PFI1$vol_subcortical)) %>%
  str_replace("scs", "cdk")
PFI_volume_brain <- mediation_fdr(PFI_sig_volumes, "smri_vol_", "PFI_", "cdk_")

# mediation - DFI 8 RSFCs (20 significant results)
DFI_fc_brain <- mediation_fdr(colnames(DFI1$networks), "rsfmri_c_ngd_", 
                              "DFI_", "ngd_")

# mediation - DFI 36 cortico-subcortical RSFCs (119 significant results)
DFI_sub_brain <- mediation_fdr(colnames(DFI1$network_sub), "rsfmri_cor_ngd_", 
                         "DFI_", "scs_")


mediation_all_brain <- rbind(PFI_volume_brain, DFI_fc_brain, DFI_sub_brain) %>%
  mutate(
    independent = str_split(Associations, "--", simplify = TRUE)[, 1],
    dependent =  str_split(Associations, "--", simplify = TRUE)[, 2],
    proportion = paste0(round(abs(PathAB_beta/PathC_beta), 4) * 100, "%"),
    PathA_beta = round(PathA_beta, 3),
    PathB_beta = round(PathB_beta, 3),
    PathAB_beta = round(PathAB_beta, 4),
    PathC_beta = round(PathC_beta, 3),
    PathCC_beta = round(PathCC_beta, 3),
    PathAB_CI1 = round(AB1, 4),
    PathAB_CI2 = round(AB2, 4),
    PathAB_p = round(PathAB_p, 4)
  ) %>%
  select(Mediator, independent, dependent, PathA_beta, PathB_beta, PathAB_beta, 
         PathC_beta, PathCC_beta, proportion, PathAB_CI1, PathAB_CI2, PathAB_p)
# Brain as mediator
names(mediation_all_brain)[1:2] <- c("independent", "Mediator")  

# modify variables
mediation_vars_df <- data.frame("vars" = unique(mediation_all_brain$Mediator),
  "new_vars" = c("Inferior Parietal (L)", "Prencentral (L)",
    "Superior Frontal (L)", "Insula (R)", "Putamen (L)",
    "Putamen (R)", "Pallidum (R)", "DMN-DMN", "SHN-SHN", "CON-DMN",
    "DMN-DAN", "SHN-SHN", "AN-Cde", "CON-Thp", "CON-Pt", "CON-Hip", "CON-Amg", 
    "CPN-Crcx", "CPN-Cde", "CPN-Vtdc", "DMN-NAc", "DAN-Hip", "DAN-Amg", 
    "FPN-Crcx", "FPN-Thp", "FPN-Amg", "SHN-Crcx", "SHN-Cde", "SHN-Pt", "SHN-Pl", 
    "SHN-NAc", "SMN-Amg", "SN-Vtdc", "VAN-Cde", "CON-BS", "DAN-BS", "RTN-BS", 
    "SHN-BS"))
for(i in 1:nrow(mediation_vars_df)) {
  mediation_all_brain$Mediator <- str_replace(mediation_all_brain$Mediator, 
    mediation_vars_df[i, 1], mediation_vars_df[i, 2])
}
# output results
write_csv(mediation_all_brain, "mediation_all_brain.csv")

# 118 common variables: 13 for PFI, 105 for DFI
common_vars <- rownames(mediation_all_brain)[rownames(mediation_all_brain) %in% rownames(mediation_all)]

mediation_all_brain_common <- mediation_all_brain[rownames(mediation_all_brain) %in% common_vars, ]
mediation_all_common <- mediation_all[rownames(mediation_all) %in% common_vars, ]

# t.test
brain_proportion <- str_split(mediation_all_brain_common$proportion, "%", simplify = TRUE)[, 1] %>%
  as.numeric()
peer_proportion <- str_split(mediation_all_common$proportion, "%", simplify = TRUE)[, 1] %>%
  as.numeric()
t.test(brain_proportion, peer_proportion, paired = TRUE)

# plot -------------------------------------------------------------------------
plot_mat <- function(data, index, group) {
  df <- str_split(rownames(data), index, simplify = TRUE)[, 2] %>%
    str_split("--", simplify = TRUE) %>%
    as.data.frame()
  df$p <- data$PathAB_p
  df$peers <- group
  return(df)
}

modify_behavior_name <- function(data) {
  data$inventory <- str_split(data$V2, " - ", simplify = TRUE)[, 1]
  data$V2 <- str_split(data$V2, " - ", simplify = TRUE)[, 2]
  data$inventory <- factor(
    data$inventory,
    levels = c("CBCL", "Life Events", "Life Events (Parent)",
               "Adverse Peer Experiences", "Prodromal Psychosis",
               "Subsyndromal Mania (Parent)", "Inhibition and Reward-seeking",
               "Impulsivity", "Neurocognition"))
  data <- arrange(data, desc(inventory))
  return(data)
}

# volume - PFI
PFI_volume_all_plot <- plot_mat(PFI_volume_brain, "PFI_", "PFI")

# RSFC - DFI
DFI_fc_plot <- plot_mat(DFI_fc_brain, "DFI_", "DFI")

DFI_fc_plot$V1 <- str_replace(DFI_fc_plot$V1, "_", "-") %>%
  str_replace("dt", "DMN") %>%
  str_replace("dt", "DMN") %>%
  str_replace("cgc", "CON") %>%
  str_replace("ca", "CPN") %>%
  str_replace("dla", "DAN") %>%
  str_replace("vs", "VN") %>%
  str_replace("smm", "SMN") %>%
  str_replace("smh", "SHN") %>%
  str_replace("smh", "SHN")
DFI_fc_plot <- arrange(DFI_fc_plot, V1)

# combine volume and RSFC
mediation_plot <- rbind(DFI_fc_plot, PFI_volume_all_plot) %>%
  modify_behavior_name()

# modify variables
mediation_plot$V1 <- factor(mediation_plot$V1,
  levels = c(unique(DFI_fc_plot$V1), "ifpllh", "paracnlh", "precnlh",
    "sufrlh", "postcnrh", "insularh", "putamenlh", "putamenrh", "pallidumrh",
    "aar"),
  labels = c(unique(DFI_fc_plot$V1), "Inferior Parietal (L)",
    "Paracentral (L)", "Precentral (L)", "Superior Frontal (L)", 
    "Postcentral (R)", "Insula (R)", "Putamen (L)", "Putamen (R)", 
    "Pallidum (R)", "Accumbens (R)"))

mediation_plot$inventory <- mediation_plot$inventory %>%
  str_replace("Prodromal Psychosis", "Prodromal\nPsychosis") %>%
  str_replace("Inhibition and Reward-seeking", "Inhibition and\nReward-seeking") %>%
  str_replace("Adverse Peer Experiences", "Adverse Peer\nExperiences") %>%
  str_replace("Life Events ", "Life Events\n") %>%
  factor(levels = c("CBCL",  "Life Events", "Life Events\n(Parent)",
    "Adverse Peer\nExperiences", "Prodromal\nPsychosis",
    "Inhibition and\nReward-seeking", "Impulsivity", "Neurocognition"))

mediation_plot$V2 <- factor(mediation_plot$V2,
  levels = unique(arrange(mediation_plot, desc(inventory), desc(V2))$V2))

# sFigure 13c ------------------------------------------------------------------
ggplot(mediation_plot, aes(V2, V1, color = peers)) +
  geom_point(shape = 17, alpha = 0.7, size = 8) +
  labs(x = NULL, y = NULL) +
  facet_grid(~inventory, scales = "free", switch = "x") +
  theme_bw() +
  theme_cleveland() +
  scale_color_manual(name = "Independent Variables", values = c("#2166AC", "#7FC97F")) +
  guides(
    shape = guide_legend("Independent Variables"),
    color = guide_legend(title.hjust = 0.5, title.position = "top",
                         override.aes = list(size = 10),
                         title.theme = element_text(size = 25, face = "bold"),
                         byrow = TRUE, reverse = TRUE)) +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.key.spacing.y = unit(0.5, "cm"),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25, face = "bold"),
    legend.key.height = unit(1, 'cm'),
    legend.key.width = unit(1, 'cm'),
    legend.position = "top",
    legend.box.background = element_blank(),
    legend.background = element_rect(fill='transparent'),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    strip.placement = "outside",
    strip.background  = element_rect(fill='transparent', color=NA),
    strip.text = element_text(size = 20, face = "bold"),
    panel.spacing = unit(2, "lines"))
ggsave("sFigure13c.svg", width = 24, height = 13, bg = "transparent")

# sFigure14 b, cortico-subcoritcal RSFCs ---------------------------------------
DFI_sub_plot <- plot_mat(DFI_sub_brain, "DFI_", "DFI")
DFI_sub_plot$V1 <- str_replace(DFI_sub_plot$V1, "_", "-") %>%
  str_replace("df", "DMN") %>%
  str_replace("cerc", "CON") %>%
  str_replace("fopa", "FPN") %>%
  str_replace("copa", "CPN") %>%
  str_replace("vta", "VAN") %>%
  str_replace("dsa", "DAN") %>%
  str_replace("sa", "SN") %>%
  str_replace("rst", "RTN") %>%
  str_replace("vs", "VN") %>%
  str_replace("au", "AN") %>%
  str_replace("smh", "SHN") %>%
  str_replace("smm", "SMN") %>%
  str_replace("crcx", "Cerebellum Cortex") %>%
  str_replace("thp", "Thalamus") %>%
  str_replace("cde", "Caudate") %>%
  str_replace("pt", "Putamen") %>%
  str_replace("pl", "Pallidum") %>%
  str_replace("hp", "Hippocampus") %>%
  str_replace("ag", "Amygdala") %>%
  str_replace("aa", "Accumbens") %>%
  str_replace("vtdc", "Ventral Diencephalon") %>%
  str_replace("bs", "Brain Stem")
DFI_sub_plot$V1 <- paste0(str_split(DFI_sub_plot$V1, "-", simplify = TRUE)[, 2],
                          "-", str_split(DFI_sub_plot$V1, "-", simplify = TRUE)[, 1])

DFI_sub_plot <- modify_behavior_name(DFI_sub_plot)

DFI_sub_plot$inventory <- DFI_sub_plot$inventory %>%
  str_replace("Prodromal Psychosis", "Prodromal\nPsychosis") %>%
  str_replace("Inhibition and Reward-seeking", "Inhibition and\nReward-seeking") %>%
  str_replace("Adverse Peer Experiences", "Adverse Peer\nExperiences") %>%
  str_replace("Life Events ", "Life Events\n") %>%
  factor(levels = c("CBCL",  "Life Events", "Life Events\n(Parent)",
                    "Adverse Peer\nExperiences", "Prodromal\nPsychosis",
                    "Inhibition and\nReward-seeking", "Impulsivity", "Neurocognition"))


DFI_sub_plot$V2 <- factor(
  DFI_sub_plot$V2,
  levels = unique(arrange(DFI_sub_plot, desc(inventory), desc(V2))$V2)
)

ggplot(DFI_sub_plot, aes(V2, V1)) +
  geom_point(shape = 17, alpha = 0.8, size = 6, color = "#2166AC") +
  geom_hline(aes(yintercept = V2), color = "grey") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL) +
  facet_grid(~inventory, scales = "free", switch = "x") +
  theme_bw() +
  theme_cleveland() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    strip.placement = "outside",
    strip.background  = element_rect(fill='transparent', color=NA),
    strip.text = element_text(size = 15, face = "bold"),
    panel.spacing = unit(2, "lines")
  )
ggsave("sFigure14b.svg", width = 24, height = 13, bg = "transparent")
