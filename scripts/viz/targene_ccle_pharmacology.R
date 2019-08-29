#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggpmisc)})

option_list <- list(optparse::make_option(c("-c", "--classifier"),
                                          type = "character",
                                          help = "Location of classifier"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

classifier <- opt$classifier

pharm_file <- file.path(classifier,"tables", "pharmacology_predictions_targene_ccle.tsv")
pharm_full_df <- readr::read_tsv(pharm_file)
head(pharm_full_df)

comp <- as.matrix(pharm_full_df$Compound)
comp <- unique(comp)

setwd(file.path(classifier,"figures","cell_line"))

plot_drug <- function(pharm_df, compound, tissues = NULL, include_braf = FALSE, 
                      facet_tissue = TRUE, se = FALSE) {
  # Output scatter plots with correlations, visualizing drug activity
  # compared to Ras classifier Score
  #
  # Arguments:
  # pharm_df - dataframe of compound activity by cell line with Ras/BRAF status
  # compound - a specific compound to visualize
  # tissues - a list of tissues to consider plotting specifically in facets
  # include_braf - boolean to include BRAF in considering mutation status
  # facet_tissue - boolean of tissues to determine to plot in facet_wrap
  # se - boolean to plot standard error intervals in geom_smooth
  #
  # Output:
  # Scatter plot with correlation information

  pharm_subset_df <- pharm_df[pharm_df$Compound == compound, ]
  if (!is.null(tissues)) {
    pharm_subset_df <- pharm_subset_df %>%
      dplyr::filter(tissue %in% focus_tissues)
  }
  if (include_braf) {
    pharm_subset_df$targene_status[pharm_subset_df$BRAF_MUT == 1] <- 1
    legend_label <- "Targene/BRAF Status"
  } else {
    legend_label <- "Targene Status"
  }
  
  if (compound == "AZD6244") {
    compound <- "Selumetinib"
  }

formula <- y ~ x
     p <- ggplot(pharm_subset_df, aes(x = weight, y = ActArea,
                                   color = as.factor(targene_status),
                                   fill = as.factor(targene_status))) +
    geom_point(alpha = 0.5, size = 2) +
    scale_x_continuous(breaks = c(0, 0.5, 1),
                       limits = c(-0.1, 1.1)) +
    geom_smooth(method = "lm", se = se) +
    geom_segment(aes(x = 0.5, y = -0.1, xend = 0.5, yend = 6),
                 linetype = "dashed", color = "grey") +
    scale_fill_manual(values = c("#377eb8", "#ff7f00"),
                      name = legend_label,
                      breaks = c(0, 1),
                      labels = c("Wild-Type", "Mutant")) +
    scale_color_manual(values = c("#377eb8", "#ff7f00"),
                       name = legend_label,
                       breaks = c(0, 1),
                       labels = c("Wild-Type", "Mutant")) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 label.x.npc = 0.17, label.y.npc = 0.92,
                 formula = formula,
                 parse = TRUE, size = 4, na.rm = TRUE,
                 rr.digits = 1) +
    stat_fit_glance(method = "lm", geom = "text",
                    label.x.npc = 0.8, label.y.npc = 0.97,
                    method.args = list(formula = formula), size = 4,
                    aes(label = paste("P = ",
                                      signif(..p.value.., digits = 1),
                                      sep = ""))) +
    xlab("Targene Classifier Score") +
    ylab("Activity Area") +
    ggtitle(compound, subtitle = "CCLE Response") + 
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  if (facet_tissue) {
    p <- p + facet_wrap("tissue")
  }
  
  return(p)
}

# drug responses with targene mutational effect

focus_tissues <- c("CENTRAL NERVOUS SYSTEM", "SKIN", "BREAST", "HAEMATOPOIETIC AND LYMPHOID TISSUE", "LARGE INTESTINE",
                   "LUNG", "OVARY", "PANCREAS", "LIVER")

pdf("./targene_all_drug_response.pdf")

sel_fig <- plot_drug(pharm_full_df, "AZD6244", facet_tissue = FALSE, se = TRUE)
sel_fig

p_fig_AEW541 <- plot_drug(pharm_full_df, "AEW541", facet_tissue = FALSE, se = TRUE)
p_fig_AEW541

p_fig_Nilotinib <- plot_drug(pharm_full_df, "Nilotinib", facet_tissue = FALSE, se = TRUE)
p_fig_Nilotinib

p_fig_17 <- plot_drug(pharm_full_df, "17-AAG", facet_tissue = FALSE, se = TRUE)
p_fig_17

p_fig_PHA665752 <- plot_drug(pharm_full_df, "PHA-665752", facet_tissue = FALSE, se = TRUE)
p_fig_PHA665752

p_fig_Lapatinib <- plot_drug(pharm_full_df, "Lapatinib", facet_tissue = FALSE, se = TRUE)
p_fig_Lapatinib

p_fig_Nutlin <- plot_drug(pharm_full_df, "Nutlin-3", facet_tissue = FALSE, se = TRUE)
p_fig_Nutlin

p_fig_AZD0530 <- plot_drug(pharm_full_df, "AZD0530", facet_tissue = FALSE, se = TRUE)
p_fig_AZD0530

p_fig_PF2341066 <- plot_drug(pharm_full_df, "PF2341066", facet_tissue = FALSE, se = TRUE)
p_fig_PF2341066

p_fig_L685458 <- plot_drug(pharm_full_df, "L-685458", facet_tissue = FALSE, se = TRUE)
p_fig_L685458

p_fig_ZD6474 <- plot_drug(pharm_full_df, "ZD-6474", facet_tissue = FALSE, se = TRUE)
p_fig_ZD6474

p_fig_Panobinostat <- plot_drug(pharm_full_df, "Panobinostat", facet_tissue = FALSE, se = TRUE)
p_fig_Panobinostat

p_fig_Sorafenib <- plot_drug(pharm_full_df, "Sorafenib", facet_tissue = FALSE, se = TRUE)
p_fig_Sorafenib

p_fig_Irinotecan <- plot_drug(pharm_full_df, "Irinotecan", facet_tissue = FALSE, se = TRUE)
p_fig_Irinotecan

p_fig_Topotecan <- plot_drug(pharm_full_df, "Topotecan", facet_tissue = FALSE, se = TRUE)
p_fig_Topotecan

p_fig_LBW242 <- plot_drug(pharm_full_df, "LBW242", facet_tissue = FALSE, se = TRUE)
p_fig_LBW242

p_fig_PD0325901 <- plot_drug(pharm_full_df, "PD-0325901", facet_tissue = FALSE, se = TRUE)
p_fig_PD0325901

p_fig_PD0332991 <- plot_drug(pharm_full_df, "PD-0332991", facet_tissue = FALSE, se = TRUE)
p_fig_PD0332991

p_fig_Paclitaxel <- plot_drug(pharm_full_df, "Paclitaxel", facet_tissue = FALSE, se = TRUE)
p_fig_Paclitaxel

p_fig_PLX4720 <- plot_drug(pharm_full_df, "PLX4720", facet_tissue = FALSE, se = TRUE)
p_fig_PLX4720

p_fig_RAF265 <- plot_drug(pharm_full_df, "RAF265", facet_tissue = FALSE, se = TRUE)
p_fig_RAF265

p_fig_TAE684 <- plot_drug(pharm_full_df, "TAE684", facet_tissue = FALSE, se = TRUE)
p_fig_TAE684

p_fig_TKI258 <- plot_drug(pharm_full_df, "TKI258", facet_tissue = FALSE, se = TRUE)
p_fig_TKI258

dev.off()

# drug responses with targene + braf mutational effect
focus_tissues <- c("CENTRAL NERVOUS SYSTEM", "SKIN", "BREAST", "HAEMATOPOIETIC AND LYMPHOID TISSUE", "LARGE INTESTINE",
                   "LUNG", "OVARY", "PANCREAS", "LIVER")
pdf("./targene_braf_all_drug_response.pdf")

sel_fig <- plot_drug(pharm_full_df, "AZD6244", facet_tissue = FALSE, se = TRUE)
sel_fig

p_fig_AEW541 <- plot_drug(pharm_full_df, "AEW541", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_AEW541

p_fig_Nilotinib <- plot_drug(pharm_full_df, "Nilotinib", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_Nilotinib

p_fig_17 <- plot_drug(pharm_full_df, "17-AAG", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_17

p_fig_PHA665752 <- plot_drug(pharm_full_df, "PHA-665752", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_PHA665752

p_fig_Lapatinib <- plot_drug(pharm_full_df, "Lapatinib", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_Lapatinib

p_fig_Nutlin <- plot_drug(pharm_full_df, "Nutlin-3", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_Nutlin

p_fig_AZD0530 <- plot_drug(pharm_full_df, "AZD0530", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_AZD0530

p_fig_PF2341066 <- plot_drug(pharm_full_df, "PF2341066", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_PF2341066

p_fig_L685458 <- plot_drug(pharm_full_df, "L-685458", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_L685458

p_fig_ZD6474 <- plot_drug(pharm_full_df, "ZD-6474", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_ZD6474

p_fig_Panobinostat <- plot_drug(pharm_full_df, "Panobinostat", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_Panobinostat

p_fig_Sorafenib <- plot_drug(pharm_full_df, "Sorafenib", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_Sorafenib

p_fig_Irinotecan <- plot_drug(pharm_full_df, "Irinotecan", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_Irinotecan

p_fig_Topotecan <- plot_drug(pharm_full_df, "Topotecan", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_Topotecan

p_fig_LBW242 <- plot_drug(pharm_full_df, "LBW242", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_LBW242

p_fig_PD0325901 <- plot_drug(pharm_full_df, "PD-0325901", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_PD0325901

p_fig_PD0332991 <- plot_drug(pharm_full_df, "PD-0332991", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_PD0332991

p_fig_Paclitaxel <- plot_drug(pharm_full_df, "Paclitaxel", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_Paclitaxel

p_fig_PLX4720 <- plot_drug(pharm_full_df, "PLX4720", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_PLX4720

p_fig_RAF265 <- plot_drug(pharm_full_df, "RAF265", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_RAF265

p_fig_TAE684 <- plot_drug(pharm_full_df, "TAE684", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_TAE684

p_fig_TKI258 <- plot_drug(pharm_full_df, "TKI258", facet_tissue = FALSE, include_braf = TRUE, se = TRUE)
p_fig_TKI258

dev.off()

# Tissue responses targene mutational effect
focus_tissues <- c("CENTRAL NERVOUS SYSTEM", "SKIN", "BREAST", "HAEMATOPOIETIC AND LYMPHOID TISSUE", "LARGE INTESTINE",
                   "LUNG", "OVARY", "PANCREAS", "LIVER")

pdf("./targene_selected_tissue_response.pdf")
sel_fig <- plot_drug(pharm_full_df, "AZD6244",   tissues = focus_tissues)
sel_fig

p_fig_AEW541 <- plot_drug(pharm_full_df, "AEW541", tissues = focus_tissues)
p_fig_AEW541

p_fig_Nilotinib <- plot_drug(pharm_full_df, "Nilotinib",   tissues = focus_tissues)
p_fig_Nilotinib

p_fig_17 <- plot_drug(pharm_full_df, "17-AAG",   tissues = focus_tissues)
p_fig_17

p_fig_PHA665752 <- plot_drug(pharm_full_df, "PHA-665752",   tissues = focus_tissues)
p_fig_PHA665752

p_fig_Lapatinib <- plot_drug(pharm_full_df, "Lapatinib",   tissues = focus_tissues)
p_fig_Lapatinib

p_fig_Nutlin <- plot_drug(pharm_full_df, "Nutlin-3",   tissues = focus_tissues)
p_fig_Nutlin

p_fig_AZD0530 <- plot_drug(pharm_full_df, "AZD0530",   tissues = focus_tissues)
p_fig_AZD0530

p_fig_PF2341066 <- plot_drug(pharm_full_df, "PF2341066",   tissues = focus_tissues)
p_fig_PF2341066

p_fig_L685458 <- plot_drug(pharm_full_df, "L-685458",   tissues = focus_tissues)
p_fig_L685458

p_fig_ZD6474 <- plot_drug(pharm_full_df, "ZD-6474",   tissues = focus_tissues)
p_fig_ZD6474

p_fig_Panobinostat <- plot_drug(pharm_full_df, "Panobinostat",   tissues = focus_tissues)
p_fig_Panobinostat

p_fig_Sorafenib <- plot_drug(pharm_full_df, "Sorafenib",   tissues = focus_tissues)
p_fig_Sorafenib

p_fig_Irinotecan <- plot_drug(pharm_full_df, "Irinotecan",   tissues = focus_tissues)
p_fig_Irinotecan

p_fig_Topotecan <- plot_drug(pharm_full_df, "Topotecan",   tissues = focus_tissues)
p_fig_Topotecan

p_fig_LBW242 <- plot_drug(pharm_full_df, "LBW242",   tissues = focus_tissues)
p_fig_LBW242

p_fig_PD0325901 <- plot_drug(pharm_full_df, "PD-0325901",   tissues = focus_tissues)
p_fig_PD0325901

p_fig_PD0332991 <- plot_drug(pharm_full_df, "PD-0332991",   tissues = focus_tissues)
p_fig_PD0332991

p_fig_Paclitaxel <- plot_drug(pharm_full_df, "Paclitaxel",   tissues = focus_tissues)
p_fig_Paclitaxel

p_fig_PLX4720 <- plot_drug(pharm_full_df, "PLX4720",   tissues = focus_tissues)
p_fig_PLX4720

p_fig_RAF265 <- plot_drug(pharm_full_df, "RAF265",   tissues = focus_tissues)
p_fig_RAF265

p_fig_TAE684 <- plot_drug(pharm_full_df, "TAE684",   tissues = focus_tissues)
p_fig_TAE684

p_fig_TKI258 <- plot_drug(pharm_full_df, "TKI258",   tissues = focus_tissues)
p_fig_TKI258

dev.off()

# Tissue responses with targene + braf mutational effect
focus_tissues <- c("CENTRAL NERVOUS SYSTEM", "SKIN", "BREAST", "HAEMATOPOIETIC AND LYMPHOID TISSUE", "LARGE INTESTINE",
                   "LUNG", "OVARY", "PANCREAS", "LIVER")

pdf("./targene_braf_selected_tissue_response.pdf")

sel_fig <- plot_drug(pharm_full_df, "AZD6244",   tissues = focus_tissues)
sel_fig

p_fig_AEW541 <- plot_drug(pharm_full_df, "AEW541", include_braf = TRUE, tissues = focus_tissues)
p_fig_AEW541

p_fig_Nilotinib <- plot_drug(pharm_full_df, "Nilotinib",   include_braf = TRUE, tissues = focus_tissues)
p_fig_Nilotinib

p_fig_17 <- plot_drug(pharm_full_df, "17-AAG",   include_braf = TRUE, tissues = focus_tissues)
p_fig_17

p_fig_PHA665752 <- plot_drug(pharm_full_df, "PHA-665752",   include_braf = TRUE, tissues = focus_tissues)
p_fig_PHA665752

p_fig_Lapatinib <- plot_drug(pharm_full_df, "Lapatinib",   include_braf = TRUE, tissues = focus_tissues)
p_fig_Lapatinib

p_fig_Nutlin <- plot_drug(pharm_full_df, "Nutlin-3",   include_braf = TRUE, tissues = focus_tissues)
p_fig_Nutlin

p_fig_AZD0530 <- plot_drug(pharm_full_df, "AZD0530",   include_braf = TRUE, tissues = focus_tissues)
p_fig_AZD0530

p_fig_PF2341066 <- plot_drug(pharm_full_df, "PF2341066",   include_braf = TRUE, tissues = focus_tissues)
p_fig_PF2341066

p_fig_L685458 <- plot_drug(pharm_full_df, "L-685458",   include_braf = TRUE, tissues = focus_tissues)
p_fig_L685458

p_fig_ZD6474 <- plot_drug(pharm_full_df, "ZD-6474",   include_braf = TRUE, tissues = focus_tissues)
p_fig_ZD6474

p_fig_Panobinostat <- plot_drug(pharm_full_df, "Panobinostat",   include_braf = TRUE, tissues = focus_tissues)
p_fig_Panobinostat

p_fig_Sorafenib <- plot_drug(pharm_full_df, "Sorafenib",   include_braf = TRUE, tissues = focus_tissues)
p_fig_Sorafenib

p_fig_Irinotecan <- plot_drug(pharm_full_df, "Irinotecan",   include_braf = TRUE, tissues = focus_tissues)
p_fig_Irinotecan

p_fig_Topotecan <- plot_drug(pharm_full_df, "Topotecan",   include_braf = TRUE, tissues = focus_tissues)
p_fig_Topotecan

p_fig_LBW242 <- plot_drug(pharm_full_df, "LBW242",   include_braf = TRUE, tissues = focus_tissues)
p_fig_LBW242

p_fig_PD0325901 <- plot_drug(pharm_full_df, "PD-0325901",   include_braf = TRUE, tissues = focus_tissues)
p_fig_PD0325901

p_fig_PD0332991 <- plot_drug(pharm_full_df, "PD-0332991",   include_braf = TRUE, tissues = focus_tissues)
p_fig_PD0332991

p_fig_Paclitaxel <- plot_drug(pharm_full_df, "Paclitaxel",   include_braf = TRUE, tissues = focus_tissues)
p_fig_Paclitaxel

p_fig_PLX4720 <- plot_drug(pharm_full_df, "PLX4720",   include_braf = TRUE, tissues = focus_tissues)
p_fig_PLX4720

p_fig_RAF265 <- plot_drug(pharm_full_df, "RAF265",   include_braf = TRUE, tissues = focus_tissues)
p_fig_RAF265

p_fig_TAE684 <- plot_drug(pharm_full_df, "TAE684",   include_braf = TRUE, tissues = focus_tissues)
p_fig_TAE684

p_fig_TKI258 <- plot_drug(pharm_full_df, "TKI258",   include_braf = TRUE, tissues = focus_tissues)
p_fig_TKI258

dev.off()
