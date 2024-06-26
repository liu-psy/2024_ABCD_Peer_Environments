# The Bivariate Latent Change Score were used to validate the results of CLPM

library(lavaan)
library(lmerTest)
library(ggwordcloud)
library(tidyverse)

setwd("H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/supplementary_analysis")
load("H:/ABCD/Relsease4.0/Package_1194636/results/peer_environments/results.RData")

# Bivariate Latent Change Score ------------------------------------------------
# time: year2 - year3
blcs <- function(x) {
  
  # regressed out covariates
  residual_value <- function(y, data) {
    fit <- lmer(data[[y]] ~ sex + interview_age + race_ethnicity +
                  income_parent + edu_parent + (1|site_id_l/rel_family_id), data = data)
    return(residuals(fit))
  }
  
  residual_year2_PFI <- residual_value("pbp_ss_prosocial_peers", abcd_year2_behavior)
  residual_year3_PFI <- residual_value("pbp_ss_prosocial_peers", abcd_year3_behavior)
  residual_year2_DFI <- residual_value("pbp_ss_rule_break", abcd_year2_behavior)
  residual_year3_DFI <- residual_value("pbp_ss_rule_break", abcd_year3_behavior)
  
  residual_year2_x <- residual_value(x, abcd_year2_behavior)
  residual_year3_x <- residual_value(x, abcd_year3_behavior)
  
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
  
  blcs_data <- inner_join(residualed_year2, residualed_year3, by = "subjectkey") %>%
    select(-subjectkey)
  
  # the Bivariate Latent Change Score model
  BLCS_PFI <- '
    PFI_year3 ~ 1*PFI_year2     
    dPFI =~ 1*PFI_year3    
    dPFI ~ 1             
    PFI_year2 ~ 1         
    PFI_year3 ~ 0*1       

    x_year3 ~ 1*x_year2    
    dx =~ 1*x_year3     
    x_year3 ~ 0*1         
    x_year3 ~~ 0*x_year3  
    
    dPFI ~~ dPFI       
    PFI_year2 ~~ PFI_year2   
    PFI_year3 ~~ 0*PFI_year3    
    
    dx ~ 1           
    x_year2 ~ 1           
    dx ~~ dx        
    x_year2 ~~ x_year2      
    
    dx ~ PFI_year2 + x_year2   
    dPFI ~ x_year2 + PFI_year2   
    
    PFI_year2 ~~ x_year2    
    dPFI ~~ dx          
  '
  
  BLCS_DFI <- '
    DFI_year3 ~ 1*DFI_year2     
    dDFI =~ 1*DFI_year3    
    dDFI ~ 1             
    DFI_year2 ~ 1         
    DFI_year3 ~ 0*1       
    
    x_year3 ~ 1*x_year2    
    dx =~ 1*x_year3     
    x_year3 ~ 0*1         
    x_year3 ~~ 0*x_year3  
    
    dDFI ~~ dDFI       
    DFI_year2 ~~ DFI_year2   
    DFI_year3 ~~ 0*DFI_year3    
    
    dx ~ 1           
    x_year2 ~ 1           
    dx ~~ dx        
    x_year2 ~~ x_year2      
    
    dx ~ DFI_year2 + x_year2   
    dDFI ~ x_year2 + DFI_year2   
    
    DFI_year2 ~~ x_year2    
    dDFI ~~ dx          
  '
  
  
  fit_PFI <- lavaan(BLCS_PFI, data = blcs_data, estimator = "mlr", fixed.x = FALSE)
  fit_DFI <- lavaan(BLCS_DFI, data = blcs_data, estimator = "mlr", fixed.x = FALSE)
  
  standardized_PFI <- standardizedsolution(fit_PFI, level = 0.95)
  standardized_DFI <- standardizedsolution(fit_DFI, level = 0.95)
  
  result <- rbind(standardized_PFI[17:20, ], standardized_DFI[17:20, ])
  return(result)
}
blcs_all <- lapply(behavior_year3, blcs)
names(blcs_all) <- behavior_year3
blcs_all

blcs_results <- function(index, data) {
  df <- data.frame()
  
  for(i in seq(data)) {
    df <- rbind(df, unlist(data[[i]][index, ]))
  }
  colnames(df) <- colnames(blcs_all$cbcl_scr_syn_anxdep_r)
  
  
  # FDR corrections
  df$pfdr <- p.adjust(df$pvalue, method = "fdr")
  
  df <- select(df, est.std, ci.lower, ci.upper, se, pvalue, pfdr) %>%
    summarise_all(as.numeric)
  rownames(df) <- behavior_year3
  
  df[, 1:4] <- apply(df[, 1:4], 2, round, 3)
  return(df)
}
blcs_PFI <- blcs_results(1, blcs_all)
blcs_DFI <- blcs_results(5, blcs_all)

blcs_PFI_fdr <- filter(blcs_PFI, pfdr < 0.05)
blcs_PFI_fdr$Peer <- "PFI"
blcs_DFI_fdr <- filter(blcs_DFI, pfdr < 0.05)
blcs_DFI_fdr$Peer <- "DFI"

blcs_all_fdr <- rbind(blcs_PFI_fdr, blcs_DFI_fdr) %>%
  mutate("Variables" = rownames(.)) %>%
  select(Peer, Variables, everything()) %>%
  mutate(pvalue = sprintf("%.1e", pvalue)) %>%
  mutate(pfdr = sprintf("%.1e", pfdr))
write_csv(blcs_all_fdr, "blcs_all_fdr.csv")

# plot -------------------------------------------------------------------------
# PFI
df1 <- data.frame(
  "vars" = c("Withdrawal depression", "Internalizing", "Depression", 
             "Sluggish cognitive tempo", "Prodromal psychosis (total)",
             "Prodromal psychosis (distress)"),
  "Beta" = blcs_PFI_fdr$est
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
ggsave("sFigure15a.svg", width = 15, height = 8, bg = "transparent")

# DFI
df2 <- data.frame(
  "vars" = c("Rulebreak", "Aggressive", "Externalizing", "ADHD", 
    "Oppositional defiant", "Conduct problem", "Total bad life events", 
    "Prodromal psychosis (total)", "Prodromal psychosis (distress)", 
    "Reputation aggression", "Reputation victimization", "Overt aggression", 
    "Overt victimization", "Relational aggression"),
  "Beta" = blcs_DFI_fdr$est
)
df2 <- df2[order(blcs_DFI_fdr$est, decreasing = TRUE), ]

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
ggsave("sFigure15b.svg", width = 15, height = 8, bg = "transparent")
