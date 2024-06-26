# The validation of neurotransmitter analysis
library(ggradar)
library(ggpubr)
library(patchwork)
library(reshape2)
library(tidyverse)

setwd("H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/supplementary")

# load data
load_results <- function(files) {
  Neurotransmitters <- read.csv(files) %>%
    rename("Spearman" = Mean.Fisher.s.z..Spearman.rho., "PET" = PET.Map) %>%
    select(File, PET, Spearman, p_exact) %>%
    mutate(PET = str_replace(PET, "SERT", "5HTT")) %>%
    mutate(File = factor(File, levels = unique(File)))
  
  # FDR corrections
  p_fdr <- Neurotransmitters %>%
    group_by(File) %>%
    summarise_at(vars(p_exact), p.adjust, method = "fdr")
  
  p_fdr <- Neurotransmitters %>%
    group_by(File) %>%
    summarise_at(vars(p_exact), p.adjust, method = "fdr")
  Neurotransmitters$p_fdr <- p_fdr$p_exact
  Neurotransmitters$significants <- ifelse(Neurotransmitters$p_fdr < 0.01, "**",
                                           ifelse(Neurotransmitters$p_fdr < 0.05, "* ", " "))
  Neurotransmitters <- Neurotransmitters %>%
    mutate("PET_fdr" = paste0(Neurotransmitters$PET, Neurotransmitters$significants)) %>%
    select(File, PET, Spearman, p_fdr, PET_fdr)
  
  return(Neurotransmitters)
}
# primary results
primary_neuro <- load_results("pet/Resh.csv")
# validation results
validate_neuro <- load_results("pet/Resh_validate.csv")

# eight new maps
new_maps <- c("5HT1a", "5HT1b", "5HT2a", "5HTT", "D2", "GABAa", "VAChT", "mGluR5")
primary_neuro8 <- filter(primary_neuro, PET %in% new_maps)
validate_neuro8 <- filter(validate_neuro, PET %in% new_maps)

cor.test(primary_neuro8$Spearman, validate_neuro8$Spearman, method = "pearson")
t.test(primary_neuro8$Spearman, validate_neuro8$Spearman)

# plot -------------------------------------------------------------------------
neuro_plot <- function(files) {
  Neurotransmitters <- read.csv(files) %>%
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
  
  # plot
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
  p2 <- neurotranmitter_radar(Neurotransmitters_thick, "PFI", "#7FC97F", "b", 
                              "Thicknes")
  p3 <- neurotranmitter_radar(Neurotransmitters_area, "PFI", "#7FC97F", "c", 
                              "Area")
  p4 <- neurotranmitter_radar(Neurotransmitters_vol, "DFI", "#9ECAE1", "d", 
                              "Volume")
  p5 <- neurotranmitter_radar(Neurotransmitters_thick, "DFI", "#9ECAE1", "e", 
                              "Thickness")
  p6 <- neurotranmitter_radar(Neurotransmitters_area, "DFI", "#9ECAE1", "f", 
                              "Area")
  p1 + p2 + p3 + p4 + p5 + p6
}
neuro_plot("pet/Resh.csv")
ggsave("pet/neuro_all.svg", width = 18, height = 12)
neuro_plot("pet/Resh_validate.csv")
ggsave("sFigure11.svg", width = 18, height = 12)
