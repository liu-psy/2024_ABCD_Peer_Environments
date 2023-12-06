# Supplementary analyses
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
setwd("H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/supplementary_analysis")
load("supplementary_results.RData")

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
  "Trauma Events (Parent) - Total events",
  "Trauma Events (Parent) - Good events",
  "Trauma Events (Parent) - Bad events",
  "Trauma Events (Parent) - Good affections",
  "Trauma Events (Parent) - Bad affections",
  "Trauma Events (Parent) - Mean affections",
  "Trauma Events (Parent) - Total affections",
  "Trauma Events - Total events",
  "Trauma Events - Good events",
  "Trauma Events - Bad events",
  "Trauma Events - Good affections",
  "Trauma Events - Bad affections",
  "Trauma Events - Mean affections",
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

# Association analysis (behavior, sFigure5) ------------------------------------
mat_aossciations <- function(x, peer_type) {
  behavior_mat <- data.frame(
    "t" = unlist(x[["behavior"]][1, ]),
    "p" = unlist(x[["behavior"]][2, ])
  )
  behavior_mat$group <- c(
    rep("Neurocognition", 3),
    rep("CBCL", 20),
    "Subsyndromal\nMania (Parent)",
    rep("Trauma Events\n(Parent)", 7),
    rep("Trauma Events", 6),
    rep("Prodromal\nPsychosis", 2),
    rep("Impulsivity", 5),
    rep("Inhibition and\nReward-seeking", 7),
    rep("Adverse Peer\nExperiences", 6)
  )
  behavior_mat$group <- factor(
    behavior_mat$group,
    levels = c("CBCL",  "Trauma Events", "Trauma Events\n(Parent)",
               "Adverse Peer\nExperiences", "Prodromal\nPsychosis",
               "Subsyndromal\nMania (Parent)", "Inhibition and\nReward-seeking",
               "Impulsivity", "Neurocognition")
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
  geom_point(size = 8, position = position_jitterdodge(), shape = 17,
             alpha = 0.8) +
  scale_y_continuous(breaks = seq(-25, 20, 5)) +
  scale_color_manual(values = c("#3B7346", "#4286f4")) +
  labs(x = NULL, y = expression(paste(bolditalic("t"), " ", bold(value)))) +
  geom_hline(yintercept = 0, linewidth = 1) +
  geom_hline(yintercept = 2.20, linewidth = 1, linetype = 2) +
  geom_hline(yintercept = -2.20, linewidth = 1, linetype = 2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, hjust = 1, angle = 15, face = "bold"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 25, face = "bold"),
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

ggsave("sFigure5.svg", width = 20, height = 12, bg = "transparent")

# Association analysis (sMRI, sFigure6) ----------------------------------------
# atlas: Desikan-Killia"y atlas (aparc)
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
DFI_vol <- data.frame("hemi" = "left", region = "lateral occipital", "t" = -3.50)
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
      title = expression(paste(italic("t"), " value")),
      title.position = "left",
      label.position = "bottom",
      ticks = FALSE
    )
  )

p1 + p2 + p3 + 
  plot_layout(nrow = 3)
ggsave("sFigure6.tiff", width = 20, height = 15, bg = "transparent")

# Association analysis (RSFCs, sFigure7) ---------------------------------------
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

# FDR
cs_DFI <- cortico_subcortical(DFI1$network_sub)

all_tvalue <- c(PFI_edges$t_value, DFI_edges$t_value, cs_DFI$t_value)
col_fun <- colorRamp2(
  c(min(all_tvalue, na.rm = TRUE), 0, max(all_tvalue, na.rm = TRUE)),
  c("#377EB8", "white", "#E41A1C")
)

# layouts
svg("sFigure7.svg", width = 16, height = 8)
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
  order = c("CON", "CPN", "DMN",  "RTN", "VAN", "SN", "DAN", "FPN", "SHN",
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

# Neurotransmitters (sFigure8) -------------------------------------------------
Neurotransmitters <- read.csv("Resh.csv") %>%
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

p1 <- neurotranmitter_radar(Neurotransmitters_vol, "PFI", "#7FC97F", "a", 
                            "Volume")
p2 <- neurotranmitter_radar(Neurotransmitters_vol, "DFI", "#9ECAE1", "d", 
                            "Volume")
p3 <- neurotranmitter_radar(Neurotransmitters_area, "PFI", "#7FC97F", "b", 
                            "Area")
p4 <- neurotranmitter_radar(Neurotransmitters_area, "DFI", "#9ECAE1", "e", 
                            "Area")
p5 <- neurotranmitter_radar(Neurotransmitters_thick, "PFI", "#7FC97F", "c", 
                            "Thicknes")
p6 <- neurotranmitter_radar(Neurotransmitters_thick, "DFI", "#9ECAE1", "f", 
                            "Thicknes")
p1 + p3 + p5 + p2 + p4 + p6 + 
  plot_layout(ncol = 3)
ggsave("sFigure8_Neurotrasmitter.svg", width = 18, height = 12)

# Mediation analysis (sFigure9) ------------------------------------------------
# one PFI brain (sFigure9a)
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
ggsave("sFigure9a.svg", width = 15, height = 8)

# one DFI brain (sFigure9b)
DFI_vol$t <- -4.5
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
ggsave("sFigure9b.svg", width = 15, height = 8)

# sFigure7c
plot_mat <- function(data, index, group) {
  df <- str_split(rownames(data), index, simplify = TRUE)[, 2] %>%
    str_split(" - ", simplify = TRUE) %>%
    as.data.frame()
  df$p <- data$PathAB_p
  df$peers <- group
  return(df)
}

modify_behavior_name <- function(data) {
  for(i in 1:nrow(data)) {
    data$V2[i] <- behavior_names[which(data$V2[i] == behavior_year2)]
  }
  data$inventory <- str_split(data$V2, " - ", simplify = TRUE)[, 1]
  data$V2 <- str_split(data$V2, " - ", simplify = TRUE)[, 2]
  data$inventory <- factor(
    data$inventory,
    levels = c("CBCL", "Trauma Events", "Trauma Events (Parent)",
      "Adverse Peer Experiences", "Prodromal Psychosis",
      "Subsyndromal Mania (Parent)", "Inhibition and Reward-seeking",
      "Impulsivity", "Neurocognition"))
  data <- arrange(data, desc(inventory))
  return(data)
}

# volume - PFI
PFI_volume_all_plot <- plot_mat(PFI_volume_all, "PFI_", "PFI")

# RSFC - PFI
PFI_fc_plot <- plot_mat(PFI_fc, "PFI_", "PFI")

# RSFC - DFI
DFI_fc_plot <- plot_mat(DFI_fc, "DFI_", "DFI")

# combine PFI and DFI (RSFC)
mediation_fc_plot <- rbind(PFI_fc_plot, DFI_fc_plot)
mediation_fc_plot$V1 <- str_replace(mediation_fc_plot$V1, "_", "-") %>%
  str_replace("dt", "dmn") %>%
  str_replace("dt", "dmn") %>%
  str_replace("cgc", "con") %>%
  str_replace("ca", "cpn") %>%
  str_replace("dla", "dan") %>%
  str_replace("vs", "VN") %>%
  str_replace("smm", "SMN") %>%
  str_replace("smh", "SHN") %>%
  toupper()
mediation_fc_plot <- arrange(mediation_fc_plot, V1)

# combine volume and RSFC
mediation_plot <- rbind(mediation_fc_plot, PFI_volume_all_plot) %>%
  modify_behavior_name()

# modify variables
mediation_plot$V1 <- factor(mediation_plot$V1,
  levels = c(unique(mediation_fc_plot$V1), "ifpllh", "precnlh", "insularh", 
    "putamenlh", "putamenrh", "pallidumrh", "aar"),
  labels = c(unique(mediation_fc_plot$V1), "Inferior Parietal (L)",
    "Precentral (L)", "Insula (R)", "Putamen (L)", "Putamen (R)",
    "Pallidum (R)", "Accumbens (R)"))
mediation_plot$V2 <- factor(mediation_plot$V2,
  levels = unique(arrange(mediation_plot, desc(inventory), desc(V2))$V2))

ggplot(mediation_plot, aes(V2, V1, color = peers)) +
  geom_point(shape = 17, alpha = 0.7, size = 8) +
  labs(x = NULL, y = NULL) +
  facet_grid(~inventory, scales = "free") +
  theme_bw() +
  theme_cleveland() +
  scale_color_manual(name = "Mediator", values = c("#2166AC", "#7FC97F")) +
  guides(
    shape = guide_legend("Mediator"),
    color = guide_legend(title.hjust = 0.5, title.position = "top",
    override.aes = list(size = 10),
    title.theme = element_text(size = 25, face = "bold"),
    byrow = TRUE)) +
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
    legend.position = "right",
    legend.box.background = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    strip.text = element_blank(),
    panel.spacing = unit(2, "lines"))
ggsave("sFigure9c.svg", width = 30, height = 12, bg = "transparent")

# Mediation analysis (cortico-subcortical RSFCs, sFigure10b) -------------------
DFI_sub_plot <- plot_mat(DFI_sub, "DFI_", "DFI")
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
DFI_sub_plot$V2 <- factor(
  DFI_sub_plot$V2,
  levels = unique(arrange(DFI_sub_plot, desc(inventory), desc(V2))$V2)
)

ggplot(DFI_sub_plot, aes(V2, V1)) +
  geom_point(shape = 17, alpha = 0.8, size = 8, color = "#2166AC") +
  geom_hline(aes(yintercept = V2), color = "grey") +
  scale_y_discrete(limits=rev) +
  labs(x = NULL, y = NULL) +
  facet_grid(~inventory, scales = "free", space = "free") +
  theme_bw() +
  theme_cleveland() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 21, face = "bold"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    strip.text = element_blank(),
    panel.spacing = unit(2, "lines")
  )
ggsave("sFigure10b.svg", width = 24, height = 12, bg = "transparent")

# Longitudinal analysis (sFigure11) --------------------------------------------
# PFI (sFigure11a)
df1 <- data.frame(
  "vars" = c("Withdrawal depression", "Depression", "Sluggish cognitive tempo",
    "Prodromal psychosis (total)", "Prodromal psychosis (distress)"),
  "Beta" = PFI_beta$beta)
df1 <- df1[order(df1$Beta), ]

ggplot(df1, aes(label = vars, size = Beta, color = Beta)) +
  geom_text_wordcloud(shape = "circle", show.legend = TRUE) +
  scale_size_area(max_size = 15) +
  scale_color_gradient(low = "#00441B", high = "#A6DBA0",
                       breaks = seq(-0.06, -0.04, 0.005)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    legend.title = element_text(size = 30, face = "bold"),
    legend.key.height = unit(1.5, 'cm'),
    legend.key.width = unit(2.5, 'cm')
  ) +
  guides(
    size = FALSE,
    color = guide_colorbar(
      title = "Beta",
      title.position = "left",
      label.position = "bottom",
      ticks = TRUE,
      reverse = TRUE,
      ticks.linewidth = 1.5
    )
  )
ggsave("sFigure11a.svg", width = 15, height = 8, bg = "transparent")

# DFI (sFigure9b)
df2 <- data.frame(
  "vars" = c("Rulebreak", "Aggressive", "Externalizing", "Oppositional Defiant",
    "Conduct", "Prodromal psychosis (total)", "Prodromal psychosis (distress)",
    "Reputation aggression", "Reputation victimization",
    "Overt aggression", "Overt victimization", "Relational aggression"),
  "Beta" = DFI_beta$beta
)
df2 <- df2[order(df2$Beta, decreasing = TRUE), ]

ggplot(df2, aes(label = vars, size = Beta, color = Beta)) +
  geom_text_wordcloud(shape = "circle", show.legend = TRUE) +
  scale_size_area(max_size = 20) +
  scale_color_gradient(low = "#4286f4", high = "#373B44",
                       breaks = seq(0.035, 0.10, 0.015)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    legend.title = element_text(size = 30, face = "bold"),
    legend.key.height = unit(1.5, 'cm'),
    legend.key.width = unit(2.5, 'cm')
  ) +
  guides(
    size = FALSE,
    color = guide_colorbar(
      title = "Beta",
      title.position = "left",
      label.position = "bottom",
      ticks = TRUE,
      ticks.linewidth = 1.5
    )
  )
ggsave("sFigure11b.svg", width = 15, height = 8, bg = "transparent")