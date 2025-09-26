#!/usr/bin/env Rscript

# NLRdynamics: Dynamics of NLRs vs normal genes in the Podospora species complex
#############################################################################
# S. Lorena Ament-Velasquez
# 2024-06-06
#############################################################################
# =======================================
# Load the necessary libraries
# ============================
library(ggplot2)
library(tidyverse)
library(plyr)
library(cowplot)
library(scales) # for percent_format()
# library(patchwork)
# ============================

# ============================
# Read file names
# ============================
### Snakemake
# Read data
geneevo_wd <- read.csv(snakemake@input$geneevo, sep = ";") # The manually curated table of gene statuses
dists <- read.table(snakemake@input$dists, header = TRUE) # Population genetics statistics

## Output
# Main figures
Dynamics_plots <- snakemake@output$dynamics
curation_plot <- snakemake@output$curation

# ### Local
# geneevo_wd <- read.csv("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/data/EvolutionAllGenes.csv", sep = ";")
# dists <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Orthologs_stats_anserina.txt", header = TRUE)
# 
# # Output
# Dynamics_plot <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/05_NLRdynamics/results/NLRvsRandomDynamics.pdf"

# ============================
# About the gene models used for the NLRs and randomly selected genes
# ============================
## ------ Distinguishing between usual NLRs and those with High Internal Conservation
nlrpalette =  c("LIC NLR" = "cadetblue3", "HIC NLR" = "#1F78B4", "Random" = "gray30")

## I decided to rename the types of NLRs
dists <- dists %>%
  mutate(Type = dplyr::recode(Type,
                              "NLR" = "LIC NLR",
                              "NLR_HIC" = "HIC NLR",
                              "Random" = "Random"))

# Add the HIC category to the geneevo dataframe from dists
geneevo_wd <- merge(geneevo_wd, dists %>% mutate(Gene_ID = Ortholog) %>% dplyr::select(Gene_ID, Type), by = "Gene_ID")

# Make a counts dataframe
genetype_counts <- geneevo_wd %>% dplyr::select(c(Final_model, Gene_type)) %>% group_by(Final_model, Gene_type) %>% plyr::count() 

# The Podan2 gene models are generally more reliable for normal genes than for NLRs.
curationhist <- ggplot(genetype_counts, aes(x=Gene_type, y=freq, fill=Final_model)) +
  geom_bar(stat="identity", position = "fill") + 
  scale_fill_brewer(palette = "Set2") +
  theme_bw() + ylab("Proportion (%)") + xlab("Gene sample")

ggsave(file = curation_plot, plot = curationhist, width = 5, height = 5)

# Turn the counts into a contigency table
Gene_type_ct <- xtabs(freq ~ Gene_type + Final_model, data = genetype_counts)
# Perform the chi-squared test
chi_squared_result <- chisq.test(Gene_type_ct)

Gene_type_ct <- xtabs(freq ~ Gene_type + Final_model, data = genetype_counts)
chi_squared_result <- chisq.test(Gene_type_ct)
# https://stats.stackexchange.com/questions/81483/warning-in-r-chi-squared-approximation-may-be-incorrect
chi_squared_result <- chisq.test(Gene_type_ct, simulate.p.value = TRUE)

# A Fisher test is more appropriate for small sample sizes!
# https://www.youtube.com/watch?v=p2c4C2D7ZZc&ab_channel=GraphPadSoftware
fisher.test(Gene_type_ct)

# ============================
# Now let's compare their dynamics
# ============================
# Make it long-format
geneevo_lg <- gather(geneevo_wd, Species, Status, anserina:pseudocomata, factor_key=TRUE)

## The name "Present" can be confusing so I decided to rename that category to "Conserved"
geneevo_lg <- geneevo_lg %>%
  mutate(Status = ifelse(Status == "Present", "Conserved", Status))

dynamics_counts <- geneevo_lg %>% select(!c(Final_model, Gene_ID, Chromosome)) %>% group_by(Type, Species) %>% plyr::count() 
# if I don't remove Chromosome it creates "duplicated"-ish lines for each category and that makes a plotting artifact like
# https://stackoverflow.com/questions/35630360/lines-overlaid-on-ggplot2-bar-plot-when-saved-with-ggsave

# Do a group-wise transform(), splitting on "Gene_type"
dynamics_counts_rel <- plyr::ddply(dynamics_counts, .(Type, Species), transform,
                                   percentf = freq / sum(freq) * 100)

# Reorder factor levels for plotting
dynamics_counts_rel$Status <- factor(dynamics_counts_rel$Status, levels = c("Conserved", "Poly", "Pseudo", "Moved", "Absent"))

# Change the names just for pretty plotting
dynamics_counts_rel <- dynamics_counts_rel %>% mutate(Gene_type2 = recode(Type, "LIC NLR" = "LIC NLR (n = 32)", "Random" = "Random (n = 100)", "HIC NLR" = "HIC NLR (n = 22)"))
dynamics_counts_rel <- dynamics_counts_rel %>% mutate(Species2 = recode(Species, "anserina" = "P. anserina", "pauciseta" = "P. pauciseta", "comata" = "P. comata", "bellae.mahoneyii" = "P. bellae-mahoneyi", "pseudoanserina" = "P. pseudoanserina", "pseudopauciseta" = "P. pseudopauciseta", "pseudocomata" = "P. pseudocomata"))

dynaplot <- ggplot(dynamics_counts_rel, aes(x=Species2, y=percentf, fill=Status)) +
    geom_bar(stat="identity") + 
    theme_bw() + ylab("Proportion of genes") +
    theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1, face = "italic"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background =element_rect(fill="white"), 
          strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = c("Conserved" = "#44AA99", 
                                 "Poly" = "#DDCC77", 
                                 "Pseudo" = "#CC6677",
                                 "Moved" = "#AA4499",
                                 "Absent" = "#332288")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    facet_grid(Gene_type2 ~ .) +
    xlab("Species")

# Turn the counts into a contigency table (but just NLRs vs Random)
## anserina
dynamics_ct <- xtabs(freq ~ Gene_type + Status, data = dynamics_counts_rel %>% dplyr::filter(Species == "anserina"))
fisher.test(dynamics_ct) # 1.623e-07

## the other species
fisher.test(xtabs(freq ~ Gene_type + Status, data = dynamics_counts_rel %>% dplyr::filter(Species == "pauciseta")))
fisher.test(xtabs(freq ~ Gene_type + Status, data = dynamics_counts_rel %>% dplyr::filter(Species == "comata")))
fisher.test(xtabs(freq ~ Gene_type + Status, data = dynamics_counts_rel %>% dplyr::filter(Species == "bellae.mahoneyii")))
fisher.test(xtabs(freq ~ Gene_type + Status, data = dynamics_counts_rel %>% dplyr::filter(Species == "pseudoanserina")))
fisher.test(xtabs(freq ~ Gene_type + Status, data = dynamics_counts_rel %>% dplyr::filter(Species == "pseudopauciseta")))
fisher.test(xtabs(freq ~ Gene_type + Status, data = dynamics_counts_rel %>% dplyr::filter(Species == "pseudocomata")))

# If I perform 7 tests, a Bonferroni correction of the p-value, alpha = 0.05 would be
0.05/7 # 0.007142857
# So thy are all significant!

# ============================
# How often a gene is actually fine in all species?
# ============================
# Apply the function to each row of the dataframe, excluding irrelevant columns
Conserved_counts <- apply(geneevo_wd %>% select(c(anserina, pauciseta, comata, bellae.mahoneyii, pseudoanserina, pseudopauciseta, pseudocomata)), 1, FUN=function(x) sum(x=='Present' | x=='Poly') )

# Create a new dataframe with the results
geneevo_wd_count <- mutate(geneevo_wd, Conserved_count = Conserved_counts, ) %>% 
  mutate(Type2 = recode(Type, "LIC NLR" = "LIC NLR (n = 32)", "HIC NLR" = "HIC NLR (n = 22)", "Random" = "Random (n = 100)"))

# Fix the levels for plotting
geneevo_wd_count$Type <- factor(geneevo_wd_count$Type, levels = c("Random", "LIC NLR", "HIC NLR"))

Conserved_hist <- ggplot(geneevo_wd_count, aes(x = Conserved_count, fill = Type))+ 
    geom_histogram(aes(y=0.5*..density..), position='dodge', binwidth=0.5) + # Make it relative to both groups # https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # legend.position="none",
          strip.background =element_rect(fill="white"), 
          strip.text = element_text(face = "bold")) +
    xlab("Number of species with functional (Conserved+Poly) orthologs") + #  
    ylab("Proportion of genes") +
    scale_fill_manual("Type", values = nlrpalette) +
    # scale_fill_brewer('Gene type', palette = "Paired", direction = -1) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_x_continuous(breaks = seq(min(geneevo_wd_count$Conserved_count), max(geneevo_wd_count$Conserved_count), by = 1)) # To have all values of the x axis shown

## Put them together
dynaConservedplot <- plot_grid(dynaplot, Conserved_hist, nrow=2, rel_heights = c(1.5, 0.5)) # , labels=c('a', 'b')
ggsave(file = Dynamics_plots, plot = dynaConservedplot, width = 5, height = 7)


