# Theme: visualization results
library(brainconn)
library(circlize)
library(ComplexHeatmap)
library(ggseg)
library(ggradar)
library(ggpubr)
library(ggwordcloud)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(tidyverse)

# load results from analysis.R 
setwd("H:/ABCD/Release4.0/Package_1194636/results/peer_environments")
load("results.RData")

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
# Figure 2 Association analysis (behavior) -------------------------------------
# Figure 2a
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

p1 <- ggplot(behavior_friends, aes(group, t, color = peers)) +
  geom_point(size = 6, position = position_jitterdodge(), shape = 17,
    alpha = 0.8) +
  annotate("text", x = 9.3, y = 7, size = 5,
           label = expression(paste(bold(bolditalic("p")[fdr]), " ", bold("< 0.05")))) +
  scale_y_continuous(breaks = seq(-25, 20, 5)) +
  scale_color_manual(values = c("#3B7346", "#4286f4")) +
  labs(x = NULL, y = expression(paste(bolditalic("t"), " ", bold(value))),
       tag = "a") +
  geom_hline(yintercept = 0, linewidth = 1) +
  geom_hline(yintercept = 2.27, linewidth = 1, linetype = 2) +
  geom_hline(yintercept = -2.27, linewidth = 1, linetype = 2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, hjust = 1, angle = 15, face = "bold"),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.tag = element_text(size = 30, face = "bold"),
    legend.key.height = unit(1, 'cm'),
    legend.key.width = unit(1, 'cm'),
    legend.position = "top",
    legend.text = element_text(size = 25),
    legend.title = element_blank(),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(1.5, "lines")
  )  +
  guides(color = guide_legend(ncol = 2))

# Figure 2b (CBCL, PLE, PPS)
common_psychopathology <- str_detect(common_behaviors, "cbcl|ple|^pps")
common_psychopathology <- common_behaviors[common_psychopathology]

# common label
PFI_behaviors <- behavior_names[colnames(PFI$behavior) %in% colnames(PFI1$behavior)]
DFI_behaviors <- behavior_names[colnames(DFI$behavior) %in% colnames(DFI1$behavior)]
common_behaviors_label <- PFI_behaviors[DFI_behaviors %in% PFI_behaviors]
common_psychopathology_label <- str_detect(common_behaviors_label, "CBCL|Life|^Prodromal")
common_psychopathology_label <- common_behaviors_label[common_psychopathology_label]
common_psychopathology_label <- str_split(common_psychopathology_label, " - ", simplify = TRUE)[, 2]

common_p <- function(data) {
  common_data <- data[rownames(data) %in% common_psychopathology, ]
  common_data$seq <- seq(common_psychopathology)
  common_data$vars <- common_psychopathology_label
  return(common_data)
}
behavior_PFI_common <- common_p(behavior_PFI)
behavior_DFI_common <- common_p(behavior_DFI)
behavior_friends_common <- rbind(behavior_PFI_common, behavior_DFI_common)

p2 <- ggplot(behavior_friends_common, aes(vars, t, fill = peers)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) +
  scale_y_continuous(breaks = seq(-25, 20, 5)) +
  scale_fill_manual(values = c("#4286f4", "#3B7346")) +
  labs(x = NULL, y = expression(paste(bolditalic("t"), " ", bold(value))),
       tag = "b") +
  geom_hline(yintercept = 0, linewidth = 1) +
  facet_grid(~group, scales = "free_x", space = "free_x", switch = "x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, hjust = 1, angle = 35, face = "bold"),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, face = "plain"),
    plot.tag = element_text(size = 30, face = "bold"),
    legend.key.height = unit(1, 'cm'),
    legend.key.width = unit(1, 'cm'),
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1.5, "lines")
  )  +
  guides(fill = guide_legend(ncol = 2, reverse = TRUE))

p1 + p2 + plot_layout(nrow = 2)
ggsave("figures/Figure2.pdf", width = 16, height = 12, bg = "transparent")

# Figure 3 Association analysis (RSFCs) ----------------------------------------
# cortical RSFCs
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

# sub-cortical ROIs
# 1. cerebellum_cortex(crcx)                 2. thalamus_proper(thp)
# 3. caudate(cde)                            4. putamen(pt)
# 5. pallidum(pl)                            6. brain-stem(bs)
# 7. hippocampus(hp)                         8. amygdala(ag)
# 9. accumbens-area(aa)                      10.ventraldc(vtdc)
# networks
# 1. auditory network(au)                    2. cingulo-opercular network(cerc)
# 3. cingulo-parietal network(copa)          4. default network(df)
# 5. dorsal attention network(dsa)           6. fronto-parietal network(fopa)
# 7. retrosplenial temporal network(rst)     8. sensorimotor hand network(smh)
# 9. sensorimotor mouth network(smm)         10.salience network(sa)
# 11.ventral attention network(vta)          12.visual network(vs)
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
svg("figures/Figure3_RSFCs.svg", width = 16, height = 8)
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

# Figure 4 Neurotransmitters ---------------------------------------------------
Neurotransmitters <- read.csv("pet/Resh.csv") %>%
  rename("Spearman" = Mean.Fisher.s.z..Spearman.rho., "PET" = PET.Map) %>%
  select(File, PET, Spearman, p_exact) %>%
  mutate(PET = str_replace(PET, "SERT", "5HTT")) %>%
  mutate(PET = str_replace(PET, "FDOPA", "F-DOPA")) %>%
  mutate(File = factor(File, levels = unique(File)))

# FDR corrections
p_fdr <- Neurotransmitters %>%
  group_by(File) %>%
  summarise_at(vars(p_exact), p.adjust, method = "fdr")
Neurotransmitters$p_fdr <- p_fdr$p_exact
Neurotransmitters$significants <- ifelse(Neurotransmitters$p_fdr < 0.01, "**",
  ifelse(Neurotransmitters$p_fdr < 0.05, "* ", " "))
Neurotransmitters <- Neurotransmitters %>%
  mutate("PET" = paste0(Neurotransmitters$PET, Neurotransmitters$significants)) %>%
  select(File, PET, Spearman)

# select data by modality
Neurotransmitters_vol <- filter(Neurotransmitters, str_detect(File, "vol")) %>%
  mutate(File = str_replace(File, "vol_PFI.nii", "PFI")) %>%
  mutate(File = str_replace(File, "vol_DFI.nii", "DFI"))
Neurotransmitters_area <- filter(Neurotransmitters, str_detect(File, "area")) %>%
  mutate(File = str_replace(File, "area_PFI.nii", "PFI")) %>%
  mutate(File = str_replace(File, "area_DFI.nii", "DFI"))
Neurotransmitters_thick <- filter(Neurotransmitters, str_detect(File, "thick")) %>%
  mutate(File = str_replace(File, "thick_PFI.nii", "PFI")) %>%
  mutate(File = str_replace(File, "thick_DFI.nii", "DFI"))

# plot function
neurotranmitter_radar <- function(x, peer, color, title, tag) {
  data <- filter(x, File == peer) %>%
    pivot_wider(names_from = PET, values_from = Spearman)
  tags <- paste0(peer, "\n(", tag, ")")
  
  # plot radar plot
  ggradar(data, values.radar = c("-0.6", "0", "0.6"), grid.min = -0.6, 
    grid.mid = 0, grid.max = 0.6, plot.legend = FALSE, grid.label.size = 11, 
    axis.label.size = 6.5, axis.label.offset = 1.18,
    background.circle.colour = "white", gridline.mid.colour = "grey",
    group.colours = color, fill = TRUE, fill.alpha = 0.4) +
    labs(title = title, tag = tags) +
    theme(
      title = element_text(size = 30, face = "bold"),
      plot.tag = element_text(size = 20),
      plot.tag.position = "bottom"
    )
}

p1 <- neurotranmitter_radar(Neurotransmitters_vol, "DFI", "#9ECAE1", "a", 
  "Volume")
p2 <- neurotranmitter_radar(Neurotransmitters_thick, "DFI", "#9ECAE1", "b", 
  "Thicknes")
p3 <- neurotranmitter_radar(Neurotransmitters_area, "DFI", "#9ECAE1", "c", 
  "Area")
p1 + p2 + p3 + 
  plot_layout(ncol = 3)

ggsave("figures/Figure4.pdf", width = 18, height = 6)

# Figure 5 Mediation analysis --------------------------------------------------
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
    "frontal pole", "temporal pole", "transverse temporal", "insula"
  )
  replace_table <- data.frame(ABCD_regions, dk_regions)

  ROI <- str_split(colnames(data[[modality]]), pattern = "_cdk_",
    simplify = TRUE)[, 2]
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

# Figure 5a (one PFI brain)
PFI_vol <- extract_region(PFI1, "vol_cortical")
PFI_vol$t <- 4.5
ggplot(PFI_vol) +
  geom_brain(atlas = dk, aes(fill = t), size = 7, hemi = "left",
    side = "lateral", show.legend = FALSE) +
  scale_fill_gradient2(limits = c(-4.5, 4.5), low = "#9ECAE1", high = "#EF6548",
    mid = "white", midpoint = 0, na.value = "white", space = "Lab", name = "t") +
  theme_void() +
  theme(
    plot.tag = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold")
  )
ggsave("figures/Figure5a.svg", width = 15, height = 8)

# Figure 5a (one DFI brain)
DFI_vol <- data.frame("hemi" = "left", "region" = "lateral occipital", "t" = -4.5)
ggplot(DFI_vol) +
  geom_brain(atlas = dk, aes(fill = t), size = 7, hemi = "left",
    side = "lateral", show.legend = FALSE) +
  scale_fill_gradient2(limits = c(-4.5, 4.5),  low = "#9ECAE1", high = "#EF6548",
    mid = "white", midpoint = 0, na.value = "white", space = "Lab", name = "t") +
  theme_void() +
  theme(
    plot.tag = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold")
  )
ggsave("figures/Figure5b.svg", width = 15, height = 8)

# Figure 5c
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
PFI_volume_all_plot <- plot_mat(PFI_volume_all, "PFI_", "PFI")

# RSFC - DFI
DFI_fc_plot <- plot_mat(DFI_fc, "DFI_", "DFI")

DFI_fc_plot$V1 <- str_replace(DFI_fc_plot$V1, "_", "-") %>%
  str_replace("dt", "DMN") %>%
  str_replace("dt", "DMN") %>%
  str_replace("cgc", "CON") %>%
  str_replace("ca", "CPN") %>%
  str_replace("dla", "DAN") %>%
  str_replace("vs", "VN") %>%
  str_replace("smm", "SMN") %>%
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

ggplot(mediation_plot, aes(V2, V1, color = peers)) +
  geom_point(shape = 17, alpha = 0.7, size = 8) +
  labs(x = NULL, y = NULL) +
  facet_grid(~inventory, scales = "free", switch = "x") +
  theme_bw() +
  theme_cleveland() +
  scale_color_manual(name = "Mediator", values = c("#2166AC", "#7FC97F")) +
  guides(
    shape = guide_legend("Mediator"),
    color = guide_legend(title.hjust = 0.5, title.position = "top",
                         override.aes = list(size = 10),
                         title.theme = element_text(size = 25, face = "bold"),
                         byrow = TRUE, reverse = TRUE)) +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.spacing.y = unit(1, "cm"),
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
ggsave("figures/Figure5c.svg", width = 24, height = 13, bg = "transparent")

# Figure 6 Longitudinal analysis -----------------------------------------------
# PFI
df1 <- data.frame(
  "vars" = c("Withdrawal depression", "Internalizing", "Depression", 
    "Sluggish cognitive tempo", "Prodromal psychosis (total)",
    "Prodromal psychosis (distress)"),
  "Beta" = clpm_PFI_fdr$est
)
df1 <- df1[order(blcs_PFI_fdr$est), ]

ggplot(df1, aes(label = vars, size = Beta, color = Beta)) +
  geom_text_wordcloud(shape = "circle", show.legend = TRUE) +
  scale_size_area(max_size = 15) +
  scale_color_gradient(low = "#00441B", high = "#A6DBA0",
    breaks = seq(-0.06, -0.04, 0.01),
    label = c("-0.060", "-0.050", "-0.040")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    legend.title = element_text(size = 40, face = "bold"),
    legend.key.height = unit(1.5, 'cm'),
    legend.key.width = unit(2.5, 'cm')
  ) +
  guides(
    size = FALSE,
    color = guide_colorbar(
      title = expression(bolditalic(beta)),
      title.position = "left",
      label.position = "bottom",
      ticks = TRUE,
      reverse = TRUE,
      ticks.linewidth = 1.5
    )
  )
ggsave("figures/Figure6a.svg", width = 15, height = 8, bg = "transparent")

# DFI
df2 <- data.frame(
  "vars" = c("Social problem", "Rulebreak", "Aggressive", "Externalizing", 
    "Total problem", "ADHD", "Oppositional defiant", "Conduct problem", 
    "Total bad affections", "Total bad life events", "Prodromal psychosis (total)", 
    "Prodromal psychosis (distress)", "Reputation aggression", 
    "Reputation victimization", "Overt aggression", "Overt victimization", 
    "Relational aggression"),
  "Beta" = clpm_DFI_fdr$est
)
df2 <- df2[order(clpm_DFI_fdr$est, decreasing = TRUE), ]

ggplot(df2, aes(label = vars, size = Beta, color = Beta)) +
  geom_text_wordcloud(shape = "circle", show.legend = TRUE) +
  scale_size_area(max_size = 20) +
  scale_color_gradient(low = "#4286f4", high = "#373B44",
    breaks = seq(0.035, 0.10, 0.015)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    legend.title = element_text(size = 40, face = "bold"),
    legend.key.height = unit(1.5, 'cm'),
    legend.key.width = unit(2.5, 'cm')
  ) +
  guides(
    size = FALSE,
    color = guide_colorbar(
      title = expression(bolditalic(beta)),
      title.position = "left",
      label.position = "bottom",
      ticks = TRUE,
      ticks.linewidth = 1.5
    )
  )
ggsave("figures/Figure6b.svg", width = 15, height = 8, bg = "transparent")

# Supplementary Figure 2 (PFI-Area, DFI-Thickness) -----------------------------
brain_plot <- function(data, tags, title) {
  ggplot(data) +
    geom_brain(atlas = dk, aes(fill = t), size = 0.3, show.legend = TRUE,
      position = position_brain(. ~ side + hemi)) +
    scale_fill_gradient2(low = "#9ECAE1", high = "#EF6548", mid = "white", 
      midpoint = 0, na.value = "white", space = "Lab", name = "t") +
    labs(tag = tags, title = title) +
    theme_void() +
    theme(
      plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
      plot.tag = element_text(size = 30, hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold"),
      legend.key.height = unit(1.2, 'cm'),
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
}

# PFI (areas)
p1 <- extract_region(PFI1, "area") %>%
  brain_plot("a", "PFI (Area)")
# DFI (thicknes)
p2 <- extract_region(DFI1, "thick") %>%
  brain_plot("b", "DFI (Thickness)")

p1 + p2 + plot_layout(ncol = 1)
ggsave("figures/sFigure2_brain.svg", width = 16, height = 9)

# Supplementary Figure 3 (Neurotransmitters) -----------------------------------
p1 <- neurotranmitter_radar(Neurotransmitters_vol, "PFI", "#7FC97F", "a", 
  "Volume")
p2 <- neurotranmitter_radar(Neurotransmitters_thick, "PFI", "#7FC97F", "b", 
  "Thickness")
p3 <- neurotranmitter_radar(Neurotransmitters_area, "PFI", "#7FC97F", "c", 
  "Area")
p1 + p2 + p3 + 
  plot_layout(ncol = 3)
ggsave("figures/sFigure3_Neurotrasmitter.pdf", width = 18, height = 6)

# Supplementary Figure 4 (coritco-subcortical RSFCs, mediation)-----------------
# sFigure 4a (brain connectivity)
set.seed(12)
df1 <- matrix(rnorm(68^2), 68, 68)
df1 <- ifelse(df1 > 3, 1, 0)
df1[lower.tri(df1)] <- t(df1)[lower.tri(df1)]
brainconn(conmat = df1, atlas = "dk68", view = "left",
          node.size = 35, edge.alpha = 0.7, edge.width = 15, edge.color = "#2470a0",
          show.legend = FALSE) +
  scale_color_brewer(palette = "RdYlGn") +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
# ggsave("figures/sFigure4a.svg", width = 15, height = 8, bg = "transparent")

# sFigure 4b (brain connectivity)
set.seed(12)
df1 <- matrix(rnorm(68^2), 68, 68)
df1 <- ifelse(df1 > 2.85, 1, 0)
df1[lower.tri(df1)] <- t(df1)[lower.tri(df1)]
brainconn(conmat = df1, atlas = "dk68", view = "left",
          node.size = 35, edge.alpha = 0.7, edge.width = 15, edge.color = "#2470a0",
          show.legend = FALSE) +
  scale_color_brewer(palette = "RdYlGn") +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
# ggsave("figures/sFigure4b.svg", width = 15, height = 8, bg = "transparent")

# sFigure 4c
PFI_sub_plot <- plot_mat(PFI_sub, "PFI_", "PFI")
DFI_sub_plot <- plot_mat(DFI_sub, "DFI_", "DFI")
mediation_sub_plot <- rbind(PFI_sub_plot, DFI_sub_plot)

mediation_sub_plot$V1 <- str_replace(mediation_sub_plot$V1, "_", "-") %>%
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

mediation_sub_plot <- modify_behavior_name(mediation_sub_plot)
mediation_sub_plot$V2 <- factor(
  mediation_sub_plot$V2,
  levels = unique(arrange(mediation_sub_plot, desc(inventory), desc(V2))$V2)
)

mediation_sub_plot$inventory <- mediation_sub_plot$inventory %>%
  str_replace("Prodromal Psychosis", "Prodromal\nPsychosis") %>%
  str_replace("Inhibition and Reward-seeking", "Inhibition and\nReward-seeking") %>%
  str_replace("Adverse Peer Experiences", "Adverse Peer\nExperiences") %>%
  str_replace("Life Events ", "Life Events\n") %>%
  factor(levels = c("CBCL",  "Life Events", "Life Events\n(Parent)",
    "Adverse Peer\nExperiences", "Prodromal\nPsychosis",
    "Inhibition and\nReward-seeking", "Impulsivity", "Neurocognition"))

ggplot(mediation_sub_plot, aes(V2, V1, color = peers)) +
  geom_point(shape = 17, alpha = 0.8, size = 7) +
  geom_hline(aes(yintercept = V2), color = "grey") +
  scale_y_discrete(limits=rev) +
  labs(x = NULL, y = NULL) +
  scale_color_manual(name = "Mediator", values = c("#2166AC", "#7FC97F")) +
  facet_grid(~inventory, scales = "free", switch = "x") +
  guides(
    shape = guide_legend("Mediator"),
    color = guide_legend(title.hjust = 0.5, title.position = "top",
                         override.aes = list(size = 10),
                         title.theme = element_text(size = 25, face = "bold"),
                         byrow = TRUE, reverse = TRUE)) +
  theme_bw() +
  theme_cleveland() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.spacing.y = unit(1, "cm"),
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
ggsave("figures/sFigure4c.svg", width = 24, height = 13, bg = "transparent")

# Others -----------------------------------------------------------------------
# one random brain (neurotransmitters brain, Figure 1)
dk_atlas <- ggseg(atlas = dk, fill = "red", colour = "black", size = 1,
  hemi = "left", view = "lateral", show.legend = FALSE) + theme_void()

brains <- data.frame("hemi" = dk_atlas$data$hemi, "side" = dk_atlas$data$side,
  "region" = dk_atlas$data$region)
brains <- brains[!duplicated(brains), ]
brains$fills <- rnorm(nrow(brains))

ggplot(brains) +
  geom_brain(atlas = dk, aes(fill = fills), size = 1, hemi = "left",
    side = "lateral", show.legend = FALSE) +
  scale_fill_gradient2(limits = c(-2, 2), low = "#9ECAE1", high = "#EF6548",
    mid = "white", midpoint = 0, na.value = "white", space = "Lab") +
  theme_void()
ggsave("figures/random_brain.png", width = 2.17, height = 1.46, bg = "transparent")