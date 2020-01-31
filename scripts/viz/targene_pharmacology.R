#!/usr/bin/env Rscript
# Pancancer_Aberrant_Pathway_Activity_Analysis scripts/viz/targene_pharmacology.R
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggpmisc)})

option_list <- list(optparse::make_option(c("-c", "--classifier"),
                                          type = "character",
                                          help = "Location of classifier"),
                    optparse::make_option(c("-p", "--compound"),
                                          type = "character",
                                          help = "list of compounds"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

classifier <- opt$classifier

gdsc1_pharm_file <- file.path(classifier,"tables", "gdsc1_targene_pharmacology_predictions.tsv")
gdsc1_pharm_full_df <- readr::read_tsv(gdsc1_pharm_file)

gdsc2_pharm_file <- file.path(classifier,"tables", "gdsc2_targene_pharmacology_predictions.tsv")
gdsc2_pharm_full_df <- readr::read_tsv(gdsc2_pharm_file)

gdsc1_ccle_file <- file.path(classifier,"tables", "gdsc1_ccle_targene_pharmacology_predictions.tsv")
gdsc1_ccle_full_df <- readr::read_tsv(gdsc1_ccle_file)

gdsc2_ccle_file <- file.path(classifier,"tables", "gdsc2_ccle_targene_pharmacology_predictions.tsv")
gdsc2_ccle_full_df <- readr::read_tsv(gdsc2_ccle_file)

comp_list <- as.matrix(read.table(opt$compound, sep = '\t',header = F))

setwd(file.path(classifier,"figures","cell_line"))

plot_drug <- function(pharm_df, compound, tissues = NULL, include_braf = FALSE, 
                      facet_tissue = TRUE, se = FALSE) {
  # Output scatter plots with correlations, visualizing drug activity
  # compared to Targene classifier Score
  #
  # Arguments:
  # pharm_df - dataframe of compound activity by cell line with targene status
  # compound - a specific compound to visualize
  # tissues - a list of tissues to consider plotting specifically in facets
  # include_braf - boolean to include BRAF in considering mutation status
  # facet_tissue - boolean of tissues to determine to plot in facet_wrap
  # se - boolean to plot standard error intervals in geom_smooth
  #
  # Output:
  # for stat_poly_eq(aes(label = paste(..rr.label..)))
  # Scatter plot with correlation information

  pharm_subset_df <- pharm_df[(pharm_df$Compound == compound),]
  if (!is.null(tissues)) {
    pharm_subset_df <- pharm_subset_df %>%
      dplyr::filter(tissue %in% focus_tissues)
  }
   else {
    legend_label <- "Targene Status"
  }
 
formula <- y ~ x
     p <- ggplot(pharm_subset_df, aes(x = weight, y = LN_IC50,
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
    stat_poly_eq(aes(label = NA),
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
    ylab("LN_IC50") +
    ggtitle(compound, subtitle = "GDSC Response") + 
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

# GDSC1 drug response for specific compounds 
focus_tissues <- c("CENTRAL NERVOUS SYSTEM", "SKIN", "BREAST", "HAEMATOPOIETIC AND LYMPHOID TISSUE", "LARGE INTESTINE",
                   "LUNG", "OVARY", "PANCREAS", "LIVER")

pdf("./GDSC1_targene_all_drug_response.pdf")

#comp_list <- as.matrix(read.csv("/data/vijay/git/pancancer/sanger_GDSC_data/final_data_for_testing/gdsc_ccle_common/compounds.csv", sep = ',',header = F))
#comp_list

for (i in comp_list)
{
    i <- plot_drug(gdsc1_pharm_full_df, i, facet_tissue = FALSE, se = TRUE)
    print(i)
}

dev.off()

# GDSC2 drug response for specific compounds 
focus_tissues <- c("CENTRAL NERVOUS SYSTEM", "SKIN", "BREAST", "HAEMATOPOIETIC AND LYMPHOID TISSUE", "LARGE INTESTINE",
                   "LUNG", "OVARY", "PANCREAS", "LIVER")

pdf("./GDSC2_targene_all_drug_response.pdf")

for (i in comp_list)
{
    i <- plot_drug(gdsc2_pharm_full_df, i, facet_tissue = FALSE, se = TRUE)
    print(i)
}

dev.off()

# GDSC1_ccle_drug response for specific compounds 
focus_tissues <- c("CENTRAL NERVOUS SYSTEM", "SKIN", "BREAST", "HAEMATOPOIETIC AND LYMPHOID TISSUE", "LARGE INTESTINE",
                   "LUNG", "OVARY", "PANCREAS", "LIVER")

pdf("./GDSC1_ccle_targene_all_drug_response.pdf")

for (i in comp_list)
{
    i <- plot_drug(gdsc1_ccle_full_df, i, facet_tissue = FALSE, se = TRUE)
    print(i)
}

dev.off()

# GDSC2_ccle_drug response for specific compounds 
focus_tissues <- c("CENTRAL NERVOUS SYSTEM", "SKIN", "BREAST", "HAEMATOPOIETIC AND LYMPHOID TISSUE", "LARGE INTESTINE",
                   "LUNG", "OVARY", "PANCREAS", "LIVER")

pdf("./GDSC2_ccle_targene_all_drug_response.pdf")

for (i in comp_list)
{
    i <- plot_drug(gdsc2_ccle_full_df, i, facet_tissue = FALSE, se = TRUE)
    print(i)
}

dev.off()
