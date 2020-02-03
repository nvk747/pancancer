#!/usr/bin/env Rscript

# Pancancer_Aberrant_Pathway_Activity_Analysis
# PanCancer Classifier
# scripts/compare_within_models.R
#
# Plots pancancer classifier performance compared to within cancer
#
# Usage: Run by assigning where the within classifier summary is and where the
#        Pan Cancer classifier summary is
#
#     Rscript scripts/compare_within_models.R 
#
#     with the required flags:
#         --pancan_summary      Directory of where classifier summary is
#         --within_dir          Directory of where within cancer-type data is
#
#     and optional flags:
#         --alt_gene            Directory of classifier summary for alt gene
#
# Output:
# Bar Plots for each comparison
# script scripts/compare_within_models.R --pancan_summary $gene_dir \
#        --within_dir $gene_dir'/within_disease' \
#        --alt_gene 'classifiers/TP53_ALT'

library(ggplot2)
library(dplyr)

# This function returns the absolute path of the script file called by Rscript
# Follows symlinks via normalizePath
get_script_path <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    cmdArgsTrailing <- commandArgs(trailingOnly = TRUE)
    cmdArgs <- cmdArgs[seq.int(from=1, length.out=length(cmdArgs) - length(cmdArgsTrailing))]
    res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmdArgs)
    res <- tail(res[res != ""], 1)
    if (length(res) > 0)
        return (normalizePath(res))
    NULL
}

source(file.path(dirname(get_script_path()),"util", "pancancer_util.R"))

option_list <- list(optparse::make_option(c("-p", "--pancan_summary"),
                                          type = "character",
                                          help = "location of pancan summary"),
                    optparse::make_option(c("-a", "--alt_gene"),
                                          type = "character",
                                          help = "location of alt gene summary",
                                          default = NULL),
                    make_papaa_version_option())

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

do_papaa_version_option(opt)

pan_summary_dir <- opt$pancan_summary
alt_gene_dir <- opt$alt_gene

# Process PanCancer Classifier and summary files
pan_summary <- file.path(pan_summary_dir, "classifier_summary.txt")
pancan_list <- parse_summary(pan_summary)
pancan_auroc_df <- process_classifier_summary(summary_list = pancan_list,
                                              model_type = "Pan",
                                              perf_type = "AUROC")
pancan_aupr_df <- process_classifier_summary(summary_list = pancan_list,
                                             model_type = "Pan",
                                             perf_type = "AUPR")
# Process Within Cancer Results
within_dir = file.path(opt$pancan_summary,"within_disease")
within_disease_files <- list.files(within_dir,
                                   pattern = "classifier_summary.txt",
                                   full.names = TRUE, recursive = TRUE)

within_disease_auroc <- data.frame()
within_disease_aupr <- data.frame()

for (file in within_disease_files) {
  disease_summary <- parse_summary(file)
  dis_auroc_df <- process_classifier_summary(summary_list = disease_summary,
                                             model_type = "Pan_within",
                                             perf_type = "AUROC")
  within_disease_auroc <- rbind(within_disease_auroc, dis_auroc_df)
  
  dis_aupr_df <- process_classifier_summary(summary_list = disease_summary,
                                            model_type = "Pan_within",
                                            perf_type = "AUPR")
  within_disease_aupr <- rbind(within_disease_aupr, dis_aupr_df)
}

# Plot a comparison of within disease classifiers to pan cancer classifier
auroc_plot <- plyr::rbind.fill(within_disease_auroc, pancan_auroc_df)
auroc_plot$AUROC <- as.numeric(paste(auroc_plot$Performance_type))
auroc_plot$Model <- factor(auroc_plot$Model, levels = c('Pan_within', 'Pan'))

aupr_plot <- plyr::rbind.fill(within_disease_aupr, pancan_aupr_df)
aupr_plot$AUPR <- as.numeric(paste(aupr_plot$Performance_type))
aupr_plot$Model <- factor(aupr_plot$Model, levels = c('Pan_within', 'Pan'))

base_comparison_theme <- theme_bw() + within_theme +
  theme(legend.position = c(1.07, 0.65),
        legend.background = element_rect(fill = alpha("white", 0)),
        plot.margin = unit(c(0.2, 1.5, 0, 0.1), "cm"))

# Create directory structure for figures
dir.create(file.path(pan_summary_dir, "figures"))

# AUROC comparison figure
options(repr.plot.res= 300)
auroc_comparison_plot <- ggplot(auroc_plot, aes(x = Disease, y = AUROC,
                                                fill = Model)) +
  base_comparison_theme +
  geom_bar(position = position_dodge2(preserve = "single"), stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("brown", "gold")) +
  geom_hline(yintercept = 0.5, linetype = "longdash", size = 0.4,
             color = "black") +
  ylab("CV AUROC") +
  coord_cartesian(ylim = c(0.4, 1)) +
  scale_y_continuous(breaks = seq(0.4, 1, 0.1))

auroc_comp_fig <- file.path(pan_summary_dir, "figures", "auroc_comparison.pdf")
pdf(auroc_comp_fig, width = 4.2, height = 1.5)
auroc_comparison_plot
dev.off()

# AUPR comparison figure
options(repr.plot.res= 300)
aupr_comparison_plot <- ggplot(aupr_plot, aes(x = Disease, y = AUPR,
                                              fill = Model)) +
  base_comparison_theme +
  geom_bar(position = position_dodge2(preserve = "single"), stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("brown", "gold")) +
  ylab("CV AUPR") +
  coord_cartesian(ylim = c(0.4, 1)) +
  scale_y_continuous(breaks = seq(0.4, 1, 0.1))

aupr_comp_fig <- file.path(pan_summary_dir, "figures", "aupr_comparison.pdf")
pdf(aupr_comp_fig, width = 4.2, height = 1.5)
aupr_comparison_plot
dev.off()

# processing altgene classifier summary and altgene within disease summary

if (!is.null(alt_gene_dir)) {
  # Compile alternative gene model performance dataframe
  pan_altgene_summary <- file.path(alt_gene_dir, "classifier_summary.txt")
  pancan_altgene_list <- parse_summary(pan_altgene_summary)
  
  pancan_altgene_auroc_df <- process_classifier_summary(summary_list = pancan_altgene_list,
                                                model_type = "altgene",
                                                perf_type = "AUROC")
  pancan_altgene_aupr_df <- process_classifier_summary(summary_list = pancan_altgene_list,
                                               model_type = "altgene",
                                               perf_type = "AUPR")
  
  # Process Within Cancer Results
  within_altgene_folder <- file.path(alt_gene_dir, "within_disease")
  within_altgene_disease_files <- list.files(within_altgene_folder,
                                    pattern = "classifier_summary.txt",
                                    full.names = TRUE, recursive = TRUE)
  
  within_disease_altgene_auroc <- data.frame()
  within_disease_altgene_aupr <- data.frame()

  for (file in within_altgene_disease_files) {
    dis_alt_summary <- parse_summary(file)
    dis_auroc_df <- process_classifier_summary(summary_list = dis_alt_summary,
                                            model_type = "alt_within",
                                            perf_type = "AUROC")
    within_disease_altgene_auroc <- rbind(within_disease_altgene_auroc, dis_auroc_df)
    
    dis_aupr_df <- process_classifier_summary(summary_list = dis_alt_summary,
                                              model_type = "alt_within",
                                              perf_type = "AUPR")
    within_disease_altgene_aupr <- rbind(within_disease_altgene_aupr, dis_aupr_df)
  }
  
  pancan_auroc_df <- plyr::rbind.fill(within_disease_altgene_auroc, pancan_altgene_auroc_df)
  pancan_aupr_df <- plyr::rbind.fill(within_disease_altgene_aupr, pancan_altgene_aupr_df)
  
  # processing pancan_alt_summary file
  # Determine alternative classifier prediction performance on alternative gene
  
  #classifier_gene <- paste(pancan_list[["Genes"]], collapse = "_")
  
  Pan_alt_auroc_df <- process_classifier_summary(summary_list = pancan_list,
                                             model_type = "Pan_alt",
                                             gene_type = "Alt gene performance",
                                             gene_class = "Alternative Genes",
                                             perf_type = "AUROC")
  Pan_alt_aupr_df <- process_classifier_summary(summary_list = pancan_list,
                                            model_type = "Pan_alt",
                                            gene_type = "Alt gene performance",
                                            gene_class = "Alternative Genes",
                                            perf_type = "AUPR")

  plot_auroc_alt <- rbind(Pan_alt_auroc_df, pancan_auroc_df)
  plot_aupr_alt <- rbind(Pan_alt_aupr_df, pancan_aupr_df)
  
  options(repr.plot.res= 300)
  # subjecting to specific cancer types
  # plot_auroc_alt$Model <- dplyr::recode(plot_auroc_alt$Model)
  # plot_aupr_alt$Model <- dplyr::recode(plot_aupr_alt$Model)
  
  plot_auroc_alt$Model <- factor(plot_auroc_alt$Model,
                                   levels = c('alt_within', 'altgene', 'Pan_alt'))
  plot_aupr_alt$Model <- factor(plot_aupr_alt$Model,
                                  levels = c('alt_within', 'altgene', 'Pan_alt'))
  # Subset to specific cancer_types
  # use_dis <- c("GBM", "LGG", "PCPG", "COAD", "OV", "UCEC")
  # plot_auroc_alt <- plot_auroc_alt %>% dplyr::filter(Disease %in% use_dis)
  # plot_aupr_alt <- plot_aupr_alt %>% dplyr::filter(Disease %in% use_dis)
  
  plot_auroc_alt$AUROC <- as.numeric(paste(plot_auroc_alt$Performance_type))
  plot_aupr_alt$AUPR <- as.numeric(paste(plot_aupr_alt$Performance_type))
  
  alt_auroc_plot <- ggplot(plot_auroc_alt, aes(x = Disease, y = AUROC,
                                               fill = Model)) +
    geom_bar(position = position_dodge2(preserve = "single"), stat = "identity", width = 0.8) +
    theme_bw() + within_theme +
    theme(legend.position = c(1.1, 0.65),
          legend.background = element_rect(fill = alpha("white", 0)),
          plot.margin = unit(c(0.2, 1.5, 0, 0.1), "cm")) +
    scale_fill_manual(values = c("brown", "gold", "cyan3")) +
    geom_hline(yintercept = 0.5, linetype = "longdash", size = 0.4,
               color = "black") +
    ylab("CV AUROC") +
    coord_cartesian(ylim = c(0.4, 1)) +
    scale_y_continuous(breaks = seq(0.4, 1, 0.1))
  
  alt_auroc_figure <- file.path(pan_summary_dir, "figures",
                               "alt_gene_auroc_comparison.pdf")
  pdf(alt_auroc_figure, width = 4.2, height = 1.5)
  print(alt_auroc_plot)
  dev.off()
  
  alt_aupr_plot <- ggplot(plot_aupr_alt, aes(x = Disease, y = AUPR,
                                             fill = Model)) +
    geom_bar(position = position_dodge2(preserve = "single"), stat = "identity", width = 0.8) +
    theme_bw() + within_theme +
    theme(legend.position = c(1.1, 0.65),
          legend.background = element_rect(fill = alpha("white", 0)),
          plot.margin = unit(c(0.2, 1.5, 0, 0.1), "cm")) +
    scale_fill_manual(values = c("brown", "gold", "cyan3")) +
    ylab("CV AUPR") +
    coord_cartesian(ylim = c(0.4, 1)) +
    scale_y_continuous(breaks = seq(0.4, 1, 0.1))
  
  alt_aupr_figure <- file.path(pan_summary_dir, "figures",
                               "alt_gene_aupr_comparison.pdf")
  pdf(alt_aupr_figure, width = 4.2 , height = 1.5)
  print(alt_aupr_plot)
  dev.off()
}
