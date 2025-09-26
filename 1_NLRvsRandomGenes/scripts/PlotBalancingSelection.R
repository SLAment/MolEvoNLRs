#!/usr/bin/env Rscript

### PlotBalancingSelection: Looking for balancing selection evidence
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/03/04
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr)
library(tidyr)
library(ggplot2)
library("dunn.test")
library(ggpubr)
library(ggrepel)
library(cowplot)

# library(patchwork)
# library(reshape2)

# ============================

# ============================
# Reading the data
# ============================
## Snakemake
checksraw <- read.table(snakemake@input$checks, header = TRUE)
dists <- read.table(snakemake@input$dists, header = TRUE)
dNdSdists <- read.table(snakemake@input$dNdS, header = TRUE)
tajima_window <- read.table(snakemake@input$winTajima, header = TRUE)

geneevo_wd <- read.table(snakemake@input$geneevo, header = TRUE, sep = ";")

balselfile <- snakemake@output$plot
aligtajimafile <- snakemake@output$altajima
wintajimafile <- snakemake@output$wintajima
# 
# ## Local
# checksraw <- read.table("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Monophyly_check.txt", header = TRUE)
# dists <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Orthologs_stats_anserina.txt", header = TRUE)
# dNdSdists <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Orthologs_stats_spp.txt", header = TRUE)
# tajima_window <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/FocalGenes_WageningenTajimaMin5Kb.txt") # Produced by NLRvsRandomGenes.smk
# geneevo_wd <- read.csv("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/data/EvolutionAllGenes.csv", sep = ";")

# ## Output
# balselfile <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/results/Fig5_BalSelection.png"
# aligtajimafile <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/results/FigS_TajimaAlignDistribution.png"
# wintajimafile <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/results/FigS_TajimaWinDistribution.png"

# ============================
# Reading the data
# ============================
# Add the level of conservation in P. anserina in particular
dists <- merge(dists, geneevo_wd %>% mutate(Ortholog = Gene_ID, Status = anserina) %>% dplyr::select(c(Ortholog, Status, ATPase)), by = c("Ortholog"))

# Make a distinction with the pseudogenes and the poly genes 
dists <- dists %>% mutate(Conservation = case_when(Status == "Present" ~ Status, TRUE ~ "Degraded"))
dists$Conservation <- factor(dists$Conservation, levels = c("Present", "Degraded")) # for plotting

# Fix the levels for plotting
dists$Type <- factor(dists$Type, levels = c("Random", "NLR", "NLR_HIC"))

# Add the Tajima values calculated in windows
names(tajima_window) <- c("Ortholog", "TajimaWin")
dists <- merge(dists, tajima_window, by = "Ortholog")

# Fuse them with the monophyly check
checkdf <- merge(checksraw, dists, by = c("Ortholog"))

## I decided to rename the types of NLRs
dists <- dists %>%
  mutate(Type = recode(Type,
                       "NLR" = "LIC NLR",
                       "NLR_HIC" = "HIC NLR",
                       "Random" = "Random"))

checkdf <- checkdf %>%
  mutate(Type = recode(Type,
                       "NLR" = "LIC NLR",
                       "NLR_HIC" = "HIC NLR",
                       "Random" = "Random"))

checkcounts <- checkdf %>% dplyr::select(c(Type, Species, Check)) %>% group_by(Type, Species, Check) %>% plyr::count() 
# Fix the levels for plotting
checkcounts$Type <- factor(checkcounts$Type, levels = c("Random", "LIC NLR", "HIC NLR"))

# Make a subset to highlight where the het genes are
hetgenes <- c("het-e", "het-d", "het-r", "PaPlp1")
highlighted_genes <- dists %>% filter(Ortholog %in% hetgenes) %>% mutate(label_expr = paste0("italic('", Ortholog, "')"))

## ------ Distinguishing between Low internal Conservation NLRs and those with High Internal Conservation
nlrpalette =  c("LIC NLR" = "cadetblue3", "HIC NLR" = "#1F78B4", "Random" = "gray30")

# ============================
# Tajima's D (main figure)
# ============================
# Prepare significance values based on the 2 pairwise comparisons
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf)/2, symbols = c("****", "***", "**", "*", "ns"))
my_comparisons <- list( c("LIC NLR", "HIC NLR"), c("LIC NLR", "Random"), c("HIC NLR", "Random") )

dists %>% dplyr::summarize(kruskal_p = kruskal.test(TajimasD ~ Type)$p.value)
kruskal.test(TajimasD ~ Type, data = dists)

# If Kruskal-Wallis is significant, perform Dunn's test for pairwise comparisons
dunn.test(list( na.omit(dists$TajimasD[dists$Type == "Random"]), na.omit(dists$TajimasD[dists$Type == "LIC NLR"]), na.omit(dists$TajimasD[dists$Type == "HIC NLR"]) ), method = "bonferroni")

( tajimaplot <- ggplot(dists, aes(x = Type, y = TajimasD, group = Type, shape = ATPase)) +
    geom_hline(yintercept = 0, linetype='dashed', colour = "gray") +
    # facet_grid(. ~ Conservation) +
    theme_classic() + 
    geom_violin(aes(colour = Type)) + 
    geom_jitter(aes(colour = Type), width = 0.1, alpha = 0.6, size = 2) + 
    # guides(shape = guide_legend(title = NULL)) +  # Removes the Type legend (since fill represents Type)
    ggpubr::stat_compare_means(comparisons = my_comparisons, 
                       label.y = c(3, 3.5, 4.2),
                       label = "p.signif",
                       symnum.args = symnum.args) + # Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent
    ggrepel::geom_text_repel(data = highlighted_genes,
                             aes(label = label_expr),
                             parse = TRUE,
                             size = 3,
                             box.padding = 0.1,
                             point.padding = 0,
                             segment.color = 'grey50') +
    theme(axis.title.x = element_blank()) +
    ylab(expression(paste("Tajima's ", italic("D")))) + 
    coord_cartesian(ylim = c(-2.5, 4.7)) +
    scale_colour_manual(values = nlrpalette) )

## Save an alternative version separating them by conservation 

( tajimaplotconserv <- ggplot(dists, aes(x = Type, y = TajimasD, group = Type, shape = ATPase)) +
    geom_hline(yintercept = 0, linetype='dashed', colour = "gray") +
    facet_grid(. ~ Conservation) +
    theme_classic() + 
    geom_violin(aes(colour = Type)) + 
    geom_jitter(aes(colour = Type), width = 0.1, alpha = 0.6, size = 2) + 
    # guides(shape = guide_legend(title = NULL)) +  # Removes the Type legend (since fill represents Type)
    ggpubr::stat_compare_means(comparisons = my_comparisons, 
                       label.y = c(3, 3.5, 4.2),
                       label = "p.signif",
                       symnum.args = symnum.args) + # Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent
    ggrepel::geom_text_repel(data = highlighted_genes,
                             aes(label = label_expr),
                             parse = TRUE,
                             size = 3,
                             box.padding = 0.1,
                             point.padding = 0,
                             segment.color = 'grey50') +
    theme(axis.title.x = element_blank()) +
    ylab(expression(paste("Tajima's ", italic("D"), " (alignment based)"))) + 
    coord_cartesian(ylim = c(-2.5, 4.7)) +
    scale_colour_manual(values = nlrpalette) )

## Save the suppl. figure
ggsave(file = aligtajimafile, plot = tajimaplotconserv, width = 7.5, height = 3)

# ============================
# Tajima's D windows version (supplementary figure)
# ============================
dists %>% dplyr::summarize(kruskal_p = kruskal.test(TajimaWin ~ Type)$p.value)
kruskal.test(TajimaWin ~ Type, data = dists)

# If Kruskal-Wallis is significant, perform Dunn's test for pairwise comparisons
dunn.test(list( na.omit(dists$TajimaWin[dists$Type == "Random"]), na.omit(dists$TajimaWin[dists$Type == "LIC NLR"]), na.omit(dists$TajimaWin[dists$Type == "HIC NLR"]) ), method = "bonferroni")

( tajimawinplot <- ggplot(dists, aes(x = Type, y = TajimaWin, group = Type, shape = ATPase)) +
    geom_hline(yintercept = 0, linetype='dashed', colour = "gray") +
    facet_grid(. ~ Conservation) +
    theme_classic() + 
    geom_violin(aes(colour = Type)) + 
    geom_jitter(aes(colour = Type), width = 0.1, alpha = 0.6, size = 2) + 
    # guides(shape = guide_legend(title = NULL)) +  # Removes the Type legend (since fill represents Type)
    ggpubr::stat_compare_means(comparisons = my_comparisons, 
                       # label.y = c(3, 3.5, 4.2),
                       label = "p.signif",
                       symnum.args = symnum.args) + # Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent
    ggrepel::geom_text_repel(data = highlighted_genes,
                             aes(label = label_expr),
                             parse = TRUE,
                             size = 3,
                             box.padding = 0.1,
                             point.padding = 0,
                             segment.color = 'grey50') +
    theme(axis.title.x = element_blank()) +
    ylab(expression(paste("Tajima's ", italic("D"), " (windows based)"))) + 
    coord_cartesian(ylim = c(-3, 7)) +
    scale_colour_manual(values = nlrpalette) )

## Save the suppl. figure
ggsave(file = wintajimafile, plot = tajimawinplot, width = 7.5, height = 3)

# ============================
# Checking for monophyly
# ============================
## Comparing different species
# ggplot(checkcounts, aes(x=Type, y=freq, fill=Check)) +
#   facet_grid(.~Species) +
#   geom_bar(stat="identity", position = "fill") + # https://ggplot2-book.org/layers.html
#   theme_classic() +
#   scale_fill_brewer(palette = "Oranges")  # RdPu

## Barplot of anserina only
checkcounts_ans <- checkcounts %>% filter(Species == "anserina")

anserbarplot <- ggplot(checkcounts_ans, aes(x=Type, y=freq, fill=Check)) +
  geom_bar(stat="identity", position = "fill") + 
  scale_fill_brewer(palette = "Reds") + # RdPu
  theme_classic() + 
  ylab("Proportion (%)") +
  coord_cartesian(ylim = c(0, 1.15)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, 
                     # label.y = c(0.57, 0.63, 0.68),
                     label = "p.signif",
                     symnum.args = symnum.args) + # Add significance levels
  theme(legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.95), # Moves legend to bottom-right (deprecated)
        # legend.position.inside = c(0.5, 0.95), # Moves legend to bottom-right
        # legend.position = "top",
        legend.background = element_blank(),
        axis.title.x = element_blank())

## Prepare data for Fisher's test
table_counts <- xtabs(freq ~ Type + Check, data = checkcounts_ans)

# Fisher’s Exact Test
fisher.test(table_counts)

# Significant, so post hoc test:
p1 <- fisher.test(table_counts[c("Random", "LIC NLR"), c("Monophyletic", "Not_monophyletic")])$p.value
p2 <- fisher.test(table_counts[c("Random", "HIC NLR"), c("Monophyletic", "Not_monophyletic") ])$p.value
p3 <- fisher.test(table_counts[c("LIC NLR", "HIC NLR"), c("Monophyletic", "Not_monophyletic") ])$p.value

# Store p-values in a vector
p_values <- c(p1, p2, p3)

# Adjust p-values using Bonferroni, Holm, or FDR
p.adjust(p_values, method = "bonferroni")  # Strict
p.adjust(p_values, method = "holm")        # Less strict
p.adjust(p_values, method = "fdr")          # Controls false discovery rate

# They are all significant. Rather than putting it in the figure I'll have it in the caption

# Put together
balselplot <- cowplot::plot_grid(tajimaplot, anserbarplot, nrow=1, labels=c('a', 'b'), rel_widths = c(2, 1.5), align = "h", axis = "b") # 

## Save the figure
ggsave(file = balselfile, plot = balselplot, width = 7.5, height = 3)

