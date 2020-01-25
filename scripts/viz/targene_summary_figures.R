# Gregory Way 2017
# PanCancer Classifier
# scripts/viz/ras_summary_figures.R
#
# Visualize summary for Ras Classifier Scores
#
# Usage: Run in command line
#
#     Rscript --vanilla scripts/viz/ras_summary_figures.R
#
# Output:
# Several figures to summarize Ras findings

library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(readr)
library(cowplot)
library(gridExtra)
library(Hmisc)
library(optparse)

# parse options
option_list = list(
  make_option(
    c("--classifier_folder"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Classifier base folder"
  )
)

opt <-parse_args(OptionParser(option_list = option_list))

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
print(get_script_path())  
print(dirname(get_script_path())) 
path = (dirname(get_script_path())) 
sp  <-  strsplit(path,'/viz') 
sp  
source(file.path(sp,"util","pancancer_util.R"))
set.seed(123)

results_folder <- opt$classifier_folder
results <- parse_summary(file.path(results_folder, "classifier_summary.txt"))
dir.create("figures")

# 1) Heatmap of the distribution of aberrant events across tumors
#heatmap_plot_file <- file.path(results_folder, "figures", "ras_heatmap.pdf")
#ras_heatmap_file <- file.path(results_folder, "figures", "all_ras_heatmap.pdf")

heatmap_plot_file <- file.path(results_folder, "figures", "targene_heatmap.pdf")
gene_heatmap_file <- file.path(results_folder, "figures", "all_targene_heatmap.pdf")

heat_file <- file.path(results_folder, "summary_counts.csv")
heat_df <- readr::read_csv(heat_file)

# processing summary counts file 

all_prop = heat_df[,grepl('_prop',colnames(heat_df))]
drop_total = all_prop[,!grepl('total',colnames(all_prop))]
prop_loss = drop_total[,!grepl('_gain',colnames(drop_total))]
prop = prop_loss[,!grepl('_loss',colnames(prop_loss))]
loss = drop_total[,grepl('_loss',colnames(drop_total))]
gain = drop_total[,grepl('_gain',colnames(drop_total))]

path_gain = as.data.frame(rowSums(gain))
path_prop = as.data.frame(rowSums(prop))
path_loss = as.data.frame(rowSums(loss))
heat_comb_df <- as.matrix(cbind(path_gain,path_loss,path_prop),index = "DISEASE")
colnames(heat_comb_df) <- c("Gain","Loss","Mutation")
rownames(heat_comb_df) <- heat_df$DISEASE

# All diseases that are used in building the classifier
targene_dis <- results[["Diseases"]]

# Build a vector for heatmap labels
classifier <- c()
for (disease in rownames(heat_comb_df)) {
  if (disease %in% targene_dis) {
    classifier <- c(classifier, "Training")
  } else {
    classifier <- c(classifier, "Dropped")
  }
}

classifier <- data.frame(classifier)
rownames(classifier) <- rownames(heat_comb_df)
classifier$classifier <- factor(classifier$classifier,
                                levels = c("Training", "Dropped"))
prop_matrix <- heat_comb_df[names(sort(heat_comb_df[,3], decreasing = TRUE)), ]

# Plot and save heatmap
options(repr.plot.width=8, repr.plot.height=2, repr.plot.res = 300)
pheatmap(t(prop_matrix * 100), scale = "none", cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.0f", fontsize_number = 8,
  
         number_color = "black", annotation_col = classifier,
         annotation_names_col = FALSE, legend = FALSE,
         filename = heatmap_plot_file,
         width = 8, height = 2)

# Plot heatmap without collapsing Ras genes
heat_targene_df <- as.matrix(cbind(gain,loss,prop))
rownames(heat_targene_df) <- heat_df$DISEASE

test_targene_df <- heat_targene_df[names(sort(heat_targene_df[,1], decreasing = TRUE)),]

# options(repr.plot.width=8, repr.plot.height=12, repr.plot.res = 300)
options(repr.plot.res = 300)
# Plot and save heatmap
h <- heat_df[grepl('_prop',colnames(heat_df))]
h = (ncol(h))
pheatmap(t(test_targene_df * 100), scale = "none", cluster_rows = FALSE,
         cluster_cols = FALSE, sort = test_targene_df[,1],
         display_numbers = TRUE, number_format = "%.0f", fontsize_number = 8,
         number_color = "black", annotation_col = classifier,
         annotation_names_col = FALSE, legend = FALSE,
         filename = gene_heatmap_file,
         width = 8, height = h * 0.5)

# 2) Coefficients contributing to the model
coef_plot_file <- file.path(results_folder, "figures", "targene_coef_plot.pdf")
coef_df <- results[["Coefficients"]]

#removing log10_mut from the coef_df_plot
coef_df <- coef_df[!coef_df$feature == 'log10_mut', -1]
coef_df <- coef_df[order(coef_df$weight, decreasing = FALSE), ]
coef_df$rank <- 1:nrow(coef_df)

# color_logic <- (coef_df$weight > 0.05 | coef_df$weight < -0.06) |
#   (coef_df$feature == 'log10_mut')
color_logic <- (coef_df$weight > 0.05 | coef_df$weight < -0.05) 

options(repr.plot.width=6, repr.plot.height=5, repr.plot.res = 300)

ggplot(coef_df, aes(x = 1:nrow(coef_df), y = weight)) +
  geom_point(color = ifelse(color_logic, 'red', 'lightgrey'),
             size = 0.01, alpha = 0.7) +
  ylab('Classifier Score') +
  xlab('Rank') +
  scale_y_continuous(breaks = seq(-0.25, 0.25, 0.05)) +
  scale_x_continuous(breaks = seq(0, 8000, 2000)) +
  geom_segment(aes(x = 0, y = 0, yend = 0, xend = nrow(coef_df)),
               colour = "navy", linetype = "dashed", size = 0.2) +
# geom_text_repel(data = subset(coef_df,
#                                (weight > 0.05 | weight < -0.05) |
#                                  coef_df$feature == 'log10_mut'),   
  geom_text_repel(data = subset(coef_df,
                                (weight > 0.05 | weight < -0.05)),
                  arrow = arrow(length = unit(0.01, 'npc')),
                  segment.size = 0.3,
                  segment.alpha = 0.3,
                  box.padding = 0.29,
                  point.padding = 0.3,
                  size = 1.3,
                  fontface = 'italic',
                  max.iter = 3e3,
                  force = 1,
                  direction = 'both',
                  xlim = c(0, 8000),
                  aes(x = rank, y = weight, label = feature)) +
  base_theme +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks = element_line(),
        axis.text = element_text(size = rel(0.55)),
        axis.title = element_text(size = rel(0.65)),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "cm"),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 3, r = 0, b = 0, l = 0)))
ggplot2::ggsave(coef_plot_file, dpi = 600, width = 1.7, height = 1.6)

 # 3) Plot distributions of predictions according to variant classification
var_gain_plot_file <- file.path(results_folder, "figures", "variant_gain_fill_map.pdf")
var_loss_plot_file <- file.path(results_folder, "figures", "variant_loss_fill_map.pdf")
mut_df <- readr::read_tsv(file.path(results_folder, "tables",
                                    "mutation_classification_scores.tsv"))

consider_mutations <- c("3'UTR", "5'UTR", "Intron", "Frame_Shift_Del",
                        "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                        "Missense_Mutation", "Nonsense_Mutation",
                        "Nonstop_Mutation", "RNA", "Splice_Site")


a = mut_df[,grepl("_gain",colnames(mut_df))]
a = a[,!grepl("copy_gain",colnames(a))]
b = mut_df[,grepl("_loss",colnames(mut_df))]
b = b[,!grepl("copy_loss",colnames(b))]

l1 = length(colnames(a))
l2 = length(colnames(b))

if (l1 >= 1){
    targene_gain <- apply(a,1,max)
    mut_df <- mut_df %>% mutate(targene_gain)
}

# if (length(colnames(b))> 1){
#    targene_loss = as.matrix(max(b[,grepl("_loss",colnames(b))]))
#    mut_df <- mut_df %>% mutate(targene_loss)
# }

if (l2 >= 1){
    targene_loss <- apply(b,1,max)
    mut_df <- mut_df %>% mutate(targene_loss)
}

silent_df <- mut_df %>% filter(Variant_Classification == "Silent") %>%
  filter(total_status == 0)
delet_df <- mut_df %>% filter(Variant_Classification %in% consider_mutations)

mut_filtered_df <- dplyr::bind_rows(delet_df, silent_df)

# Separate classes of mutations to summarize for abc_gain

a = mut_df[,grepl("targene_gain",colnames(mut_df))]
if (length(colnames(a))== 1) {
   copy_num_df <- mut_df %>% filter(targene_gain == 1) %>%
  # filter(TP53 == 0) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Loss")
  missense_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Missense_Mutation") %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Missense")
nonsense_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Nonsense_Mutation") %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Nonsense")
indel_df <- mut_filtered_df %>% filter(Variant_Classification %in%
                                       c("Frame_Shift_Del", "Frame_Shift_Ins",
                                         "In_Frame_Del", "In_Frame_Ins")) %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Indel")
utr_df <- mut_filtered_df %>%
  filter(Variant_Classification %in% c("3'UTR", "5'UTR", "Intron")) %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "UTR")
silent_df <- silent_df %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Silent")
splice_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Splice_Site") %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Splice")
wt_df <- mut_df %>% subset(total_status == 0) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "WT")
hyper_df <- mut_df %>%
  filter(hypermutated == 1) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Hyper")

final_gain_df <- dplyr::bind_rows(list(missense_df, nonsense_df, indel_df, utr_df,
                                  splice_df, silent_df, copy_num_df, wt_df,
                                  hyper_df))

colnames(final_gain_df) <- c("ID", "Gene", "Disease", "Weight", "HGVSc", "HGVSp",
                        "Class")
} else {
    print("no OG present in the dataset")
}

# Separate classes of mutations to summarize for abc_loss

b = mut_df[,grepl("targene_loss",colnames(mut_df))]
if (length(colnames(b))== 1) {
    copy_num_df <- mut_df %>% filter(targene_loss == 1) %>%
  # filter(TP53 == 0) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Loss")
missense_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Missense_Mutation") %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Missense")
nonsense_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Nonsense_Mutation") %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Nonsense")
indel_df <- mut_filtered_df %>% filter(Variant_Classification %in%
                                       c("Frame_Shift_Del", "Frame_Shift_Ins",
                                         "In_Frame_Del", "In_Frame_Ins")) %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Indel")
utr_df <- mut_filtered_df %>%
  filter(Variant_Classification %in% c("3'UTR", "5'UTR", "Intron")) %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "UTR")
silent_df <- silent_df %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Silent")
splice_df <- mut_filtered_df %>%
  filter(Variant_Classification == "Splice_Site") %>%
  filter(!(Variant_Classification %in%
             c(missense_df$Variant_Classification,
               nonsense_df$Variant_Classification))) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Splice")
wt_df <- mut_df %>% subset(total_status == 0) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "WT")
hyper_df <- mut_df %>%
  filter(hypermutated == 1) %>%
  select(Variant_Classification, Hugo_Symbol, DISEASE, weight, HGVSc, HGVSp) %>%
  mutate(classification = "Hyper")

final_loss_df <- dplyr::bind_rows(list(missense_df, nonsense_df, indel_df, utr_df,
                                  splice_df, silent_df, copy_num_df, wt_df,
                                  hyper_df))

colnames(final_loss_df) <- c("ID", "Gene", "Disease", "Weight", "HGVSc", "HGVSp",
                        "Class")

} else {
    print("no TSG are present in dataset")
}

options(repr.plot.width=4.5, repr.plot.height=4, repr.plot.res = 600)
# Plot summary distribution of variant classes prediction scores for gain
a = mut_df[,grepl("targene_gain",colnames(mut_df))]
if (length(colnames(a))== 1) {
  ggplot(final_gain_df, aes(Weight, ..count.., fill = Class)) +
  geom_density(position = "fill", size = 0.1) +
  geom_segment(aes(x = 0.5, y = 0, yend = 1, xend = 0.5), colour = "black",
               linetype = "dashed", size = 0.4) +
  labs(list(x = "Probability", y = "Proportion")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + base_theme +
  theme(legend.position = c(1.0, 0.65),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5),
        plot.margin = unit(c(0.2, 1.5, 0, 0.1),"cm"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8))
  ggsave(var_gain_plot_file, width = 4.5, height = 4, dpi = 600)
  dev.off()
} else {
    print("no OG variant plot")
}

#Plot summary distribution of variant classes prediction scores for loss
b = mut_df[,grepl("targene_loss",colnames(mut_df))]
if (length(colnames(b))== 1) {
   options(repr.plot.width=4, repr.plot.height=3.8, repr.plot.res = 300)
ggplot(final_loss_df, aes(Weight, ..count.., fill = Class)) +
  geom_density(position = "fill", size = 0.1) +
  geom_segment(aes(x = 0.5, y = 0, yend = 1, xend = 0.5), colour = "black",
               linetype = "dashed", size = 0.4) +
  labs(list(x = "Probability", y = "Proportion")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + base_theme +
  theme(legend.position = c(1.0, 0.65),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5),
        plot.margin = unit(c(0.2, 1.5, 0, 0.1),"cm"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8))
ggsave(var_loss_plot_file, width = 4.5, height = 4, dpi = 600)
dev.off() 
} else {
    print("no TSG variant plot")
}

# 4) Show mutation frequencies and scores
mut_weight_df <- mut_filtered_df %>% filter(!is.na(weight))
mut_weight_df <- mut_weight_df[mut_weight_df$hypermutated != 1, ]

aa_df <- mut_weight_df %>%
  group_by(HGVSp, Variant_Classification, Hugo_Symbol) %>%
  summarise(Mean = mean(weight, na.rm = TRUE),
            SD = sd(weight, na.rm = TRUE),
            count = n(),
            low_CI = get_boot(weight),
            high_CI = get_boot(weight, low = FALSE))

nuc_df <- mut_weight_df %>%
  group_by(HGVSc, Variant_Classification, Hugo_Symbol) %>%
  summarise(Mean = mean(weight),
            SD = sd(weight, na.rm = TRUE),
            count = n(),
            low_CI = get_boot(weight),
            high_CI = get_boot(weight, low = FALSE))

aa_df <- aa_df[order(aa_df$count, decreasing = TRUE),]
nuc_df <- nuc_df[order(nuc_df$count, decreasing = TRUE),]
write.table(aa_df, file = file.path(results_folder, "tables",
                                    "amino_acid_mutation_scores.tsv"),
            sep = "\t", row.names = FALSE)
write.table(nuc_df, file = file.path(results_folder, "tables",
                                     "nucleotide_mutation_scores.tsv"),
            sep = "\t", row.names = FALSE)

# Plot summary distribution of PTEN variants R130X and R233X prediction scores using OG based classifer

a = mut_df[,grepl("targene_gain",colnames(mut_df))]
if (length(colnames(a))== 1) {
 pten_df <- final_gain_df[complete.cases(final_gain_df), ]
 pten_R130G_df <- pten_df[pten_df$HGVSp == "p.Arg130Gly", ]
 pten_R130E_df <- pten_df[pten_df$HGVSp == "p.Arg130Gln", ]
 pten_R130L_df <- pten_df[pten_df$HGVSp == "p.Arg130Leu", ]
 pten_R233T_df <- pten_df[pten_df$HGVSp == "p.Arg233Ter", ]
 pten_germ_df <- dplyr::bind_rows(pten_R130G_df,pten_R130L_df,pten_R130E_df,pten_R233T_df)
 print(pten_germ_df$Disease, max = 200)
 print(nrow(pten_germ_df))
 pten_germ_df$Disease <- dplyr::recode(pten_germ_df$Disease,
                                 "BLCA" = "Other", "COAD" = "Other",
                                 "CESC" = "Other", "COAD" = "Other",
                                 "ESCA" = "Other", "LUAD" = "Other",
                                 "HNSC" = "Other", "STAD" = "Other", 
                                 "SARC" = "Other", "SKCM" = "Other")

 pten_germ_plot_file <- file.path(results_folder, "figures",
                            "PTEN_R130X_R233X_gain_distribution.pdf")

 pten_germ <- ggplot(pten_germ_df, aes(Weight, fill = Disease)) + 
         geom_density(alpha = 0.4) + theme_bw() +
         ylab("Density") + xlab("Classifier Score")
 
 ggsave(pten_germ_plot_file,pten_germ,width = 4, height = 3, dpi= 300)
 
}

# Plot summary distribution of PTEN variants R130X and R233X prediction scores using TSG based classifer

b = mut_df[,grepl("targene_loss",colnames(mut_df))]
if (length(colnames(b))== 1) {
 pten_df <- final_loss_df[complete.cases(final_loss_df), ]
 pten_R130G_df <- pten_df[pten_df$HGVSp == "p.Arg130Gly", ]
 pten_R130E_df <- pten_df[pten_df$HGVSp == "p.Arg130Gln", ]
 pten_R130L_df <- pten_df[pten_df$HGVSp == "p.Arg130Leu", ]
 pten_R233T_df <- pten_df[pten_df$HGVSp == "p.Arg233Ter", ]
 pten_germ_df <- dplyr::bind_rows(pten_R130G_df,pten_R130L_df,pten_R130E_df,pten_R233T_df)
 
 pten_germ_df$Disease <- dplyr::recode(pten_germ_df$Disease,
                                 "BLCA" = "Other", "COAD" = "Other",
                                 "CESC" = "Other", "COAD" = "Other",
                                 "ESCA" = "Other", "LUAD" = "Other",
                                 "HNSC" = "Other", "STAD" = "Other", 
                                 "SARC" = "Other", "SKCM" = "Other")
 pten_germ_plot_file <- file.path(results_folder, "figures",
                            "PTEN_R130X_R233X_loss_distribution.pdf")

 pten_germ <- ggplot(pten_germ_df, aes(Weight, fill = Disease)) +
  geom_density(alpha = 0.4) + theme_bw() +
  ylab("Density") + xlab("Classifier Score")

 ggsave(pten_germ_plot_file,pten_germ,width = 4, height = 3, dpi= 300)

}

# 5) Targene Summary Counts Distribution
targene_pathway_count_file <- file.path(results_folder, "tables",
                            "path_events_per_sample.tsv")
targene_pathway_summary_count_df <- readr::read_tsv(targene_pathway_count_file,
                                        col_types = cols(.default = "c",
                                                         "weight" = "d",
                                                         "total_status" = "c"))
targene_pathway_summary_count_df$copy_count <- factor(targene_pathway_summary_count_df$copy_count,
                                          levels = c("0", "1", "2", "3","4",
                                                     "5", "6", "7", "8", "9",
                                                     "10"))
targene_pathway_summary_count_df$copy_count <-
  dplyr::recode(targene_pathway_summary_count_df$copy_count, "6" = ">6", "7" = ">6",
                "8" = ">6", "9" = ">6", "10" = ">6")

# Get summary statistics for each comparison
mut_targene_pathway_prop <- targene_pathway_summary_count_df %>% group_by(mutation_count) %>%
  dplyr::summarize(mean_targene_pathway = round(mean(as.numeric(total_status)), 2))
cop_targene_pathway_prop <- targene_pathway_summary_count_df %>% group_by(copy_count) %>%
  dplyr::summarize(mean_targene_pathway = round(mean(as.numeric(total_status)), 2))

mut_targene_pathway_count <- targene_pathway_summary_count_df %>% group_by(mutation_count) %>% tally()
cop_targene_pathway_count <- targene_pathway_summary_count_df %>% group_by(copy_count) %>% tally()

# Combine to get summary tables
mut_sum <- dplyr::inner_join(mut_targene_pathway_count, mut_targene_pathway_prop, by = "mutation_count")
cop_sum <- dplyr::inner_join(cop_targene_pathway_count, cop_targene_pathway_prop, by = "copy_count")

med_weight <- median(targene_pathway_summary_count_df$weight)

#options(repr.plot.width=6.2, repr.plot.height=8.6, repr.plot.res = 300)
classifier_count_theme <- base_theme +
  theme(legend.title = element_text(size = rel(1.7)),
        legend.text = element_text(size = rel(0.9)),
        legend.key.size = unit(1, "cm"),
        legend.position = c(0.98, 0.7),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks = element_line(),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        plot.margin = unit(c(0.2, 3.6, 0.2, 0.2), "cm"))
mut <- ggplot(targene_pathway_summary_count_df, aes(x = mutation_count, y = weight)) +
  geom_boxplot(aes(fill = total_status)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(name = "Targene_Status", values = c("#3B9AB2", "#F2300F"),
                    labels = c("0" = "Wild-Type", "1" = "Activated")) +
  geom_text(data = mut_sum, aes(x = mutation_count, y = 1.08,
                                label = paste0(n, "\n", mean_targene_pathway))) +
  classifier_count_theme +
  labs(list(x = "Other Targene_pathway Mutations",
            y = "Targene_pathway Classifier Score"))

cop <- ggplot(targene_pathway_summary_count_df, aes(x = copy_count, y = weight)) +
  geom_boxplot(aes(fill = total_status)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(name = "Targene_Status", values = c("#3B9AB2", "#F2300F"),
                    labels = c("0" = "Wild-Type", "1" = "Activated")) +
  geom_text(data = cop_sum, aes(x = copy_count, y = 1.08,
                                label = paste0(n, "\n", mean_targene_pathway))) +
  classifier_count_theme +
  labs(list(x = "Other Targene_pathway Copy Number Events",
            y = "Targene_pathway Classifier Score"))

targene_pathway_counts_fig <- file.path(results_folder, "figures", "targene_pathway_events_counts.pdf")
pdf(targene_pathway_counts_fig, width = 6.5, height = 9)
plot_grid(mut, cop, align = "v", nrow = 2)
dev.off()

# 6) Performance Metrics Distribution across pathway members
options(repr.plot.width=6.2, repr.plot.height=8.6, repr.plot.res = 300)
perf_metric_file <- file.path(results_folder, "tables",
                              "all_gene_metric_ranks.tsv")
metric_ranks <- readr::read_tsv(perf_metric_file,
                                col_types = cols(.default = "c",
                                                 "AUROC" = "d",
                                                 "AUPRC" = "d",
                                                 "AUROC Rank" = "i",
                                                 "AUPRC Rank" = "i"))

# colnames(metric_ranks) <- c("Gene", "AUPR", "AUPR Rank", "ras", "AUROC",
#                             "AUROC Rank")

aupr_violin <- ggplot(metric_ranks, aes(y = AUPRC, x = paste(targene),
                                        fill = paste(targene))) +
  geom_violin() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("AUPR") +
  scale_x_discrete(labels = c("0" = "Other\nGenes",
                              "1" = "targene Pathway\nGenes")) +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks = element_line(),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)))

auroc_violin <- ggplot(metric_ranks, aes(y = AUROC, x = paste(targene),
                                         fill = paste(targene))) +
  geom_violin() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  xlab("") +
  scale_x_discrete(labels = c("0" = "Other", "1" = "targene Pathway Genes"))

aupr_plot <- ggplot(metric_ranks, aes(x = `AUPRC Rank`, y = AUPRC)) +
  geom_point(color = "darkgrey") +
  geom_point(data = metric_ranks[metric_ranks$targene == 1, ], color = "red") +
  xlab("AUPRC Rank") +
  ylab("AUPRC") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks = element_line(),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)))

auroc_plot <- ggplot(metric_ranks, aes(x = `AUROC Rank`, y = AUROC)) +
  geom_point(color = "darkgrey") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_point(data = metric_ranks[metric_ranks$targene == 1, ], color = "red")

# Get the top genes by both metrics
options(repr.plot.width=11.5, repr.plot.height=7.5, repr.plot.res = 300)
top_aupr_genes <- metric_ranks[order(metric_ranks$`AUPRC Rank`), 1:2]
top_aupr_genes <- top_aupr_genes %>% mutate(AUPRC = round(AUPRC, 2))

top_aupr_table_grob <- tableGrob(top_aupr_genes[1:15, ])
aupr_plot <- aupr_plot +
  annotation_custom(top_aupr_table_grob, xmin = 10000,
                    xmax = 15000, ymin = 0.1, ymax = 0.45)

top_auroc_genes <- metric_ranks[order(metric_ranks$`AUROC Rank`), c(1, 5)]
top_auroc_table_grob <- tableGrob(top_auroc_genes[1:10, ])
auroc_plot <- auroc_plot +
  annotation_custom(top_auroc_table_grob, xmin = 10000,
                    xmax = 15000, ymin = 0.6, ymax = 0.95)

aupr_distribution_fig <- file.path(results_folder, "figures",
                                   "aupr_distribution.pdf")

pdf(aupr_distribution_fig, width = 11.5, height = 7.5)
plot_grid(aupr_plot, aupr_violin, align = "h", ncol = 2)
dev.off()

auroc_distribution_fig <- file.path(results_folder, "figures",
                                    "auroc_distribution.pdf")

pdf(auroc_distribution_fig, width = 11, height = 7.5)
plot_grid(auroc_plot, auroc_violin, align = "h", ncol = 2)
dev.off()

# T-Test for AUPR between Ras pathway genes and Other genes
targene_pathway_genes_aupr <- metric_ranks %>%
  dplyr::filter(targene == 1) %>%
  dplyr::select(AUPRC)

other_genes_aupr <- metric_ranks %>%
  dplyr::filter(targene == 0) %>%
  dplyr::select(AUPRC)

t_test_file <- file.path(results_folder, "tables",
                         "targene_pathway_variant_AUPR_ttest.txt")
sink(t_test_file)
t.test(targene_pathway_genes_aupr$AUPRC, other_genes_aupr$AUPRC, alternative = "greater")
sink()
