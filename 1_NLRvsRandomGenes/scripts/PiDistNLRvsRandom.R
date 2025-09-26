#!/usr/bin/env Rscript

### PiDistNLRvsRandom: Distribution of average pairwise nucleotide identity 
#############################################################################
# https://robjhyndman.com/hyndsight/transformations/
# https://www.dataduets.com/2023/12/how-to-and-not-to-log-transform-zero.html
# https://library.virginia.edu/data/articles/understanding-empirical-cumulative-distribution-functions
# https://www.healthknowledge.org.uk/public-health-textbook/research-methods/1b-statistical-methods/parametric-nonparametric-tests
# https://asaip.psu.edu/articles/beware-the-kolmogorov-smirnov-test/
# https://www.geeksforgeeks.org/how-to-conduct-an-anderson-darling-test-in-r/
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2024/10/17 - 2025/03/07
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(patchwork) # for insets
library(Matching) # for ks.boot
library(ggpubr)
library("dunn.test")
library(gamlss)
# library(rcompanion)
library(ggtext)
# ============================

# ============================
# Reading the data
# ============================
## Snakemake Input
dists <- read.table(snakemake@input$dists, header = TRUE) # Population genetics statistics
phylodists <- read.table(snakemake@input$phylo, header = TRUE) # Phylogenetic distances between genes
dNdSdists_raw <- read.table(snakemake@input$dNdS, header = TRUE) # Raw dN/dS ratios based on counts (Nei and Gojobori 1986)
TEcov <- read.table(snakemake@input$TEcov, header = TRUE) # Proportion covered by transposons in 2kb flanks
geneevo_wd <- read.csv(snakemake@input$geneevo, sep = ";") # The manually curated table of gene statuses
kaksdists <- read.table(snakemake@input$kaks, header = TRUE) # dN/dS as calculated by kakscalculator
muts <- read.table(snakemake@input$muts, header = TRUE) # Mutations based on gene alignments

## Output
# Main figures
figM_dist <- snakemake@output$main
# Supplementary figures
figS_RIPmuts <- snakemake@output$RIP
figS_SupDivStats <- snakemake@output$divstats
figS_dndS <- snakemake@output$dNdS
# figS_piNpiS <- snakemake@output$piSpiN

# ## Local
# dists <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Orthologs_stats_anserina.txt", header = TRUE)
# phylodists <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Orthologs_Distances.txt", header = TRUE)
# dNdSdists_raw <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Orthologs_stats_spp.txt", header = TRUE)
# TEcov <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/FlanksTEcontent.tab", header = TRUE)
# geneevo_wd <- read.csv("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/data/EvolutionAllGenes.csv", sep = ";")
# kaksdists <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Orthologs_kakscalculator.txt", header = TRUE)
# muts <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/reports/Orthologs_mutationprofile.txt", header = TRUE)
# 
# ## Output
# # Main figures
# figM_dist <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/results/Fig4_PopPhyloStats"
# 
# # Supplementary figures
# figS_RIPmuts <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/results/FigS1_RIPmuts"
# figS_SupDivStats <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/results/FigS2_SupDivStats"
# figS_dndS <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/results/FigS3_dNdSdegraded.png"
# figS_piNpiS <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/13_NLRvsRandomGenes/results/FigS_piSpiN.png"

# ============================

# ============================
# Preparing the data frame
# ============================
# Add the TE coverage
dists <- merge(dists, TEcov, by = c("Ortholog", "Type")) 

# Make geneevo_wd from wide to long format
geneevo_lg <- gather(geneevo_wd, Species, StatusSpp, anserina:pseudocomata, factor_key=TRUE) %>% 
  mutate(Ortholog = Gene_ID) %>% dplyr::select(-c(Gene_ID, Chromosome, Gene_type))

# Add the level of conservation in P. anserina in particular
dists <- merge(dists, geneevo_wd %>% mutate(Ortholog = Gene_ID, Status = anserina) %>% dplyr::select(c(Ortholog, Status, ATPase)), by = c("Ortholog"))

## I decided to rename the types of NLRs
dists <- dists %>%
  mutate(Type = dplyr::recode(Type,
                       "NLR" = "LIC NLR",
                       "NLR_HIC" = "HIC NLR",
                       "Random" = "Random"))

muts <- muts %>% 
  mutate(Type = dplyr::recode(Type,
                       "NLR" = "LIC NLR",
                       "NLR_HIC" = "HIC NLR",
                       "Random" = "Random"))

dNdSdists_raw <- dNdSdists_raw %>% 
  mutate(Type = dplyr::recode(Type,
                       "NLR" = "LIC NLR",
                       "NLR_HIC" = "HIC NLR",
                       "Random" = "Random"))

# Fix the levels for plotting
dists$Type <- factor(dists$Type, levels = c("Random", "LIC NLR", "HIC NLR"))
muts$Type <- factor(muts$Type, levels = c("Random", "LIC NLR", "HIC NLR"))

# Make a distinction with the pseudogenes and the poly genes in P. anserina in particular
dists <- dists %>% mutate(Conservation = case_when(Status == "Present" ~ "Conserved", TRUE ~ "Degraded"))
dists$Conservation <- factor(dists$Conservation, levels = c("Conserved", "Degraded")) # for plotting

# Make a version with only conserved genes
dists_conserved <- dists %>% filter(Conservation == "Conserved")

# For dNdS estimates, also keep track of the conservation status of each species
geneevo <- merge(geneevo_lg, dists %>% dplyr::select(Ortholog, Type) )
geneevo <- geneevo %>% mutate(Conservation = case_when(StatusSpp == "Present" ~ "Conserved", TRUE ~ "Degraded"))
geneevo$Conservation <- factor(geneevo$Conservation, levels = c("Conserved", "Degraded")) # for plotting

# Fix the species names
geneevo$Species <- as.character(geneevo$Species) # it doesn't work if I let it have Species as factor
geneevo$Species[geneevo$Species == "bellae.mahoneyii"] <- "bellae-mahoneyi"

## ------ Distinguishing between usual NLRs and those with High Internal Conservation
nlrpalette =  c("LIC NLR" = "cadetblue3", "HIC NLR" = "#1F78B4", "Random" = "gray30")

# Prepare significance values based on the 2 pairwise comparisons
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf)/2, symbols = c("****", "***", "**", "*", "ns"))
my_comparisons <- list( c("LIC NLR", "HIC NLR"), c("LIC NLR", "Random"), c("HIC NLR", "Random") )

# ==============================
# Transposon coverage
# ==============================
TEviolins_HIC <- ggplot(dists, aes(y = TEcov*100, x = Type)) + 
  theme_classic() + 
  geom_violin(aes(colour = Type)) + 
  geom_jitter(aes(colour = Type, shape = ATPase),width = 0.1, alpha = 0.6) + 
  theme(legend.position = "none", axis.title.x = element_blank(), plot.background = element_blank() ) +
  scale_x_discrete(labels = c("Random" = "Rand.", "LIC NLR" = "LIC", "HIC NLR" = "HIC")) +
  ylab("TE%") + 
  scale_colour_manual(values = nlrpalette) +
  stat_compare_means(comparisons = my_comparisons, 
                     p.adjust.method = "bonferroni", 
                     method = "wilcox.test",
                     label = "p.adj", 
                     label.y = c(63, 70, 83),
                     symnum.args = symnum.args) + # Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent
  coord_cartesian(ylim = c(0, 100))

TE_pecdf_HIC <- ggplot(dists, aes(x = TEcov*100, color = Type)) +
  stat_ecdf(geom = "step", linewidth = 1, alpha = 0.8) +
  labs(x = "TE% in 2 Kbp flanks",
       y = "Cumulative\nProbability") +
  theme_classic() + 
  theme(legend.position="none") + 
  scale_colour_manual(values = nlrpalette)

TEviolins_ecd_HIC <- TE_pecdf_HIC + inset_element(TEviolins_HIC, left = 0.35, bottom = 0.02, right = 0.98, top = 0.65)

### --- TEs statistics ---
## A Kruskal-Wallis test plus Dunn's test for pairwise comparisons

# Since Kruskal-Wallis test is significant, do a Dunn's test for pairwise comparisons
dists %>% dplyr::summarize(kruskal_p = kruskal.test(TEcov ~ Type)$p.value) # 1.237971e-08
dunn.test(list(dists$TEcov[dists$Type == "LIC NLR"], dists$TEcov[dists$Type == "HIC NLR"], dists$TEcov[dists$Type == "Random"]), method = "bonferroni")

# Only conserved genes
dists_conserved %>% dplyr::summarize(kruskal_p = kruskal.test(TEcov ~ Type)$p.value) # 2.791133e-07
dunn.test(list(dists_conserved$TEcov[dists_conserved$Type == "LIC NLR"], dists_conserved$TEcov[dists_conserved$Type == "HIC NLR"], dists_conserved$TEcov[dists_conserved$Type == "Random"]), method = "bonferroni")

# Scheirer-Ray-Hare test, a non-parameter alternative to multi-factorial ANOVA analyses
# rcompanion::scheirerRayHare(TEcov ~ Type + Conservation, data = dists) # Only Type has an effect
# rcompanion::scheirerRayHare(TEcov ~ Type + Conservation, data = dists)$p.value

# ==============================
# Average GC content
# ==============================
# Kruskal-Wallis test p value
dists %>% dplyr::summarize(kruskal_p = kruskal.test(avgGC ~ Type)$p.value) # 7.052045e-15

# If Kruskal-Wallis is significant, perform Dunn's test for pairwise comparisons
dunn.test(list(dists$avgGC[dists$Type == "LIC NLR"], dists$avgGC[dists$Type == "HIC NLR"], dists$avgGC[dists$Type == "Random"]), method = "bonferroni")

# What if I exclude genes that are not conserved?
dists_conserved %>% dplyr::summarize(kruskal_p = kruskal.test(avgGC ~ Type)$p.value) # 2.452514e-10
dunn.test(list(dists_conserved$avgGC[dists_conserved$Type == "LIC NLR"], dists_conserved$avgGC[dists_conserved$Type == "HIC NLR"], dists_conserved$avgGC[dists_conserved$Type == "Random"]), method = "bonferroni")

# Scheirer-Ray-Hare test, a non-parameter alternative to multi-factorial ANOVA analyses
# rcompanion::scheirerRayHare(avgGC ~ Type + Conservation, data = dists) # Both have an effect
# rcompanion::scheirerRayHare(avgGC ~ Type + Conservation, data = dists)$p.value

GCviolins_HIC <- ggplot(dists, aes(y = avgGC*100, x = Type)) + 
  theme_classic() + 
  geom_hline(yintercept = 50, linetype='dashed', color = "gray") +
  geom_violin(aes(colour = Type), show.legend = FALSE) + 
  geom_jitter(aes(colour = Type, shape = ATPase), width = 0.1, alpha = 0.6) +
  theme(axis.title.x = element_blank(),
        legend.background = element_blank(),
        # legend.direction = "horizontal",        # Make the legend horizontal
        legend.spacing.x = unit(0.01, "cm"),      # Adjust vertical spacing between legend items
        legend.position = c(0.15, 0.2)) +  # Adjust position inside the plot (bottom-left corner)
  guides(colour = "none") + # Ensure Type legend is removed from all layers
  labs(shape = NULL) + # Remove title of the ATPase legend
  ylab("Average GC\ncontent (%)") +
  xlab("Type of gene") + 
  scale_colour_manual(values = nlrpalette) +
  coord_cartesian(ylim = c(30, 75)) +
  stat_compare_means(comparisons = my_comparisons, 
                     p.adjust.method = "bonferroni", 
                     method = "wilcox.test",
                     label = "p.adj", 
                     label.y = c(57, 63, 68),
                     symnum.args = symnum.args) # Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent

# ============================
# What are the most abundant mutation types?
# ============================
## CpA:TpG is RIP
# C:G to T:A

mutation_order <- c("A>G", "G>A", "C>T", "T>C",  # Transitions
                    "A>C", "A>T", "C>A", "C>G",  # Transversions
                    "G>C", "G>T", "T>A", "T>G")

# Prepare mutation labels with coloring (I want to paint RIP-like mutations specifically)
mutation_labels <- setNames(mutation_order, mutation_order)
mutation_labels[c("G>A", "C>T")] <- c(
  "<span style='color:red;'>G&gt;A</span>",
  "<span style='color:red;'>C&gt;T</span>"
)

# Prepare summarized data
sum_muts <- muts %>%
  count(Mutation, Type, Status, name = "Count") %>%
  group_by(Status, Type) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  mutate(Mutation = factor(Mutation, levels = mutation_order),
         Status = recode(Status, 
                         "Fixed" = "Fixed\nsites", 
                         "Polymorphic" = "Polymorphic\nsites"),
         Status = factor(Status, levels = c("Polymorphic\nsites", "Fixed\nsites")) ) # Reorder

RIPmuts <- ggplot(sum_muts, aes(x = Mutation, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Status, ncol = 1) +  # still need facets to split panels; ncol = 1 keeps separate panels for Fixed/Polymorphic, but no strip shown.
  scale_x_discrete(labels = mutation_labels) +
  scale_fill_manual(values = nlrpalette) +
  theme_classic() +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank()
  ) +
  labs(x = "Mutation Type", y = "Mutations", fill = "Gene Type") +
  # Add Status as internal text inside each panel:
  geom_text(data = sum_muts %>% group_by(Status) %>% summarise(),
            # aes(x = 10, y = 0.25, label = Status),
            aes(x = 10, y = 350, label = Status),
            inherit.aes = FALSE, fontface = "bold")

RIPmutsProps <- ggplot(sum_muts, aes(x = Mutation, y = Proportion, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Status, ncol = 1) +  # still need facets to split panels; ncol = 1 keeps separate panels for Fixed/Polymorphic, but no strip shown.
  scale_x_discrete(labels = mutation_labels) +
  scale_fill_manual(values = nlrpalette) +
  theme_classic() +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    strip.text = element_blank(),  # Remove facet labels
    legend.position="none",
    strip.background = element_blank()
  ) +
  labs(x = "Mutation Type", y = "Proportion of\nmutations", fill = "Gene Type") +
  # Add Status as internal text inside each panel:
  geom_text(data = sum_muts %>% group_by(Status) %>% summarise(),
            aes(x = 10, y = 0.25, label = Status),
            inherit.aes = FALSE, fontface = "bold")

# Make a variable for RIP-like mutations
muts <- muts %>%
  mutate(
    RIP = case_when(
      Mutation == "C>T" & Right_base %in% c("A", "T") ~ TRUE,
      Mutation == "G>A" & Left_base  %in% c("A", "T") ~ TRUE,
      TRUE ~ FALSE
    )
  )

# Fix the lavels again
muts$Status <- factor(muts$Status, levels = c("Polymorphic", "Fixed"))

reallyRIP <- muts %>%
  count(RIP, Type, Status, name = "Count") %>%
  mutate(RIP = factor(RIP, levels = c(FALSE, TRUE), labels = c("Not RIP", "RIP"))) %>%
  ggplot(aes(x = RIP, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Status~.) +
  labs(x = "Presence of RIP mutation", y = "Number of\nMutations", fill = "Gene Type") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = nlrpalette)

### Put together a RIP supplementary figure
RIPsup <- plot_grid(GCviolins_HIC, RIPmutsProps, reallyRIP,  labels = c('a', 'b', "c"), ncol = 1, align = c("h"), axis = "tblr", rel_widths = c(1, 0.5))

# --- Are they statistically different? ---
# Totals
muts %>%
  count(Type, Status, name = "Count")

sum_muts %>%
  filter(Mutation %in% c("G>A", "C>T"), Status == "Polymorphic\nsites") 

## Logistic regression: make a binary variable of RIP vs not RIP
# Define the reference factors for the model 
muts <- muts %>%
  mutate(
    Type = factor(Type, levels = c("Random", "LIC NLR", "HIC NLR")), # Set the factors such that Random is the reference group
    Status = factor(Status, levels = c("Polymorphic", "Fixed")) ) # Fixed will be compared to Polymorphic

# Model
glm_model <- glm(RIP ~ Type * Status, data = muts, family = "binomial") %>% summary()

exp(coef(glm_model))  # odds ratios
muts %>% filter(RIP == TRUE) %>% 
  count(Mutation, Type, Status, name = "Count")

muts %>% filter(RIP == FALSE) %>% 
  count(Type, Status, name = "Count")

# ==============================
# About diversity within P. anserina
# ==============================

piviolins_HIC <- ggplot(dists, aes(y = piS, x = Type, colour = Type)) + 
  geom_violin(aes(colour = Type)) + 
  # scale_y_continuous(trans = "log1p") +  # Log-transform the y-axis to spread out small values
  geom_jitter(width = 0.1, alpha = 0.5, aes(shape = ATPase)) +
  theme_classic() + 
  theme(legend.position = "none", axis.title.x = element_blank(), plot.background = element_blank()) +
  # ylab(expression(atop("Average pairwise-nucleotide", "differences (" * pi * ")"))) +
  ylab(expression(pi[S])) + 
  scale_x_discrete(labels = c("Random" = "Rand.", "LIC NLR" = "LIC", "HIC NLR" = "HIC")) +
  scale_colour_manual(values = nlrpalette) +
  stat_compare_means(comparisons = my_comparisons, 
                     p.adjust.method = "bonferroni", 
                     method = "wilcox.test",
                     label = "p.adj", 
                     label.y = c(0.11, 0.13, 0.15),
                     symnum.args = symnum.args) + # Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent
  coord_cartesian(ylim = c(0, 0.17))

# The differences are hard to see, so I could add an empirical Cumulative Distribution functions (eCDFs) Plot
pi_pecdf_HIC <- ggplot(dists, aes(x = piS, color = Type)) +
  stat_ecdf(geom = "step", linewidth = 1, alpha = 0.8) +
  # scale_x_continuous(trans = "log1p") +
  labs(x = expression(pi[S]),
       y = "Cumulative\nProbability") +
  theme_classic() + 
  # scale_colour_brewer('Gene type', palette = "Paired", direction = -1) +
  theme(legend.position="none") + 
  # geom_vline(xintercept = 0.000492, linetype = "dashed", color = "gray") +
  scale_colour_manual(values = nlrpalette)

## *** Figure 4
piviolins_ecd_HIC <- pi_pecdf_HIC + inset_element(piviolins_HIC, left = 0.3, bottom = 0.02, right = 0.98, top = 0.75)

### --- piS statistics ---
# Kruskal-Wallis test p value
dists %>% dplyr::summarize(kruskal_p = kruskal.test(piS ~ Type)$p.value) # 2.189273e-10

# If Kruskal-Wallis is significant, perform Dunn's test for pairwise comparisons
(dt <- dunn.test(list(dists$piS[dists$Type == "LIC NLR"], dists$piS[dists$Type == "HIC NLR"], dists$piS[dists$Type == "Random"]), method = "bonferroni"))

# Kruskal-Wallis test but on conserved genes
dists_conserved %>% dplyr::summarize(kruskal_p = kruskal.test(piS ~ Type)$p.value) # 3.37117e-07
dunn.test(list(dists_conserved$piS[dists_conserved$Type == "LIC NLR"], dists_conserved$piS[dists_conserved$Type == "HIC NLR"], dists_conserved$piS[dists_conserved$Type == "Random"]), method = "bonferroni")

# ==============================
# Comparing ratios of synonymous diversity vs non-synonymous
# ==============================

## Supplementary figure for Present vs degenerate piS 
piSviolinsDegrade <- ggplot(dists, aes(x = Conservation, y = log(piS), colour = Type, shape = ATPase, group = Conservation)) +
  facet_grid(. ~ Type) +
  theme_classic() + 
  geom_violin(alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  labs(y = expression("log("*pi[S]*")")) +
  scale_colour_manual(values = nlrpalette) +
  scale_fill_manual(values = nlrpalette) +
  # guides(fill = "none",  # Removes the Type legend (since fill represents Type)
  #        colour = "none",  # Removes the Type legend (since colour represents Type),
  # shape = guide_legend(title = NULL)) +  # Removes title for Check legend
  theme(strip.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        # legend.position="none",
        # legend.position = c(0.15, 0.7), # Moves legend to bottom-right
        # plot.background = element_blank()
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background.x = element_blank(),
        legend.background = element_blank() ) +
  stat_compare_means(comparisons = list(c("Conserved", "Degraded")), 
                     p.adjust.method = "bonferroni", 
                     method = "wilcox.test",
                     label = "p.adj", 
                     symnum.args = symnum.args) # Add significance levels 

# Scheirer-Ray-Hare test, a non-parameter alternative to multi-factorial ANOVA analyses
# rcompanion::scheirerRayHare(piS ~ Type + Conservation, data = dists)
# rcompanion::scheirerRayHare(piS ~ Type + Conservation, data = dists)$p.value

wilcox.test(piS ~ Conservation, data = dists, subset = (Type == "HIC NLR")) # p-value = 0.1624
wilcox.test(piS ~ Conservation, data = dists, subset = (Type == "LIC NLR")) # p-value =  0.6755 but cannot compute exact p-value with ties
wilcox.test(piS ~ Conservation, data = dists, subset = (Type == "Random")) # p-value = 0.8847

### --- The piN part

piviolins_HIC_piN <- ggplot(dists, aes(y = piN, x = Type, colour = Type)) + 
  geom_violin(aes(colour = Type)) + 
  # scale_y_continuous(trans = "log1p") +  # Log-transform the y-axis to spread out small values
  geom_jitter(width = 0.1, alpha = 0.5, aes(shape = ATPase)) +
  theme_classic() + 
  theme(legend.position = "none", axis.title.x = element_blank(), plot.background = element_blank()) +
  # ylab(expression(atop("Average pairwise-nucleotide", "differences (" * pi * ")"))) +
  ylab(expression(pi[N])) + 
  scale_colour_manual(values = nlrpalette) +
  stat_compare_means(comparisons = my_comparisons, 
                     p.adjust.method = "bonferroni", 
                     method = "wilcox.test",
                     label = "p.adj", 
                     # label.y = c(0.57, 0.63, 0.68),
                     symnum.args = symnum.args) + # Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent
  coord_cartesian(ylim = c(0, 0.09))

# The differences are hard to see, so I could add an empirical Cumulative Distribution functions (eCDFs) Plot
pi_pecdf_HIC_piN <- ggplot(dists, aes(x = piN, color = Type)) +
  stat_ecdf(geom = "step", linewidth = 1, alpha = 0.8) +
  # scale_x_continuous(trans = "log1p") +
  labs(x = expression(pi[N]),
       y = "Cumulative\nProbability") +
  theme_classic() + 
  # scale_colour_brewer('Gene type', palette = "Paired", direction = -1) +
  theme(legend.position="none") + 
  # geom_vline(xintercept = 0.000492, linetype = "dashed", color = "gray") +
  scale_colour_manual(values = nlrpalette)

## *** Figure S1
piviolins_ecd_HIC_piN <- pi_pecdf_HIC_piN + inset_element(piviolins_HIC_piN, left = 0.3, bottom = 0.02, right = 0.98, top = 0.75)

## I decided the log of the ratio was better for visualization because log makes 
## the changes in magnitude in both directions (>1 and <1) symmetrical around 0.
(piNpiSviolins <- ggplot(dists, aes(x = Conservation, y = log(piNpiS), colour = Type, shape = ATPase, group = Conservation)) +
  geom_hline(yintercept = log(1), linetype='dashed', colour = "gray") +
  facet_grid(. ~ Type) +
  theme_classic() + 
  geom_violin(alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  labs(y = expression("log("*pi[N]/pi[S]*")")) +
  scale_colour_manual(values = nlrpalette) +
  scale_fill_manual(values = nlrpalette) +
  theme(strip.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        # legend.position="none",
        # legend.position = c(0.15, 0.7), # Moves legend to bottom-right
        # plot.background = element_blank()
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background.x = element_blank(),
        legend.background = element_blank() ) +
    guides(colour = "none") + # Ensure Type legend is removed from all layers
  stat_compare_means(comparisons = list(c("Conserved", "Degraded")), 
                     p.adjust.method = "bonferroni", 
                     method = "wilcox.test",
                     label = "p.adj", 
                     symnum.args = symnum.args) +
  coord_cartesian(ylim = c(-2.5, 2)) )


dists %>% dplyr::summarize(kruskal_p = kruskal.test(piNpiS ~ Conservation)$p.value)

kruskal.test(piNpiS ~ Type, data = dists)

# Scheirer-Ray-Hare test, a non-parameter alternative to multi-factorial ANOVA analyses
# rcompanion::scheirerRayHare(piNpiS ~ Type + Conservation, data = dists)
# rcompanion::scheirerRayHare(piNpiS ~ Type + Conservation, data = dists)$p.value

shapiro.test(dists$piNpiS)  # The data is not normally distributed 

# If Kruskal-Wallis is significant, perform Dunn's test for pairwise comparisons
dists_clean <- dists %>% dplyr::select(c(piNpiS, Type)) %>% na.omit()
dunn.test(list(dists_clean$piNpiS[dists_clean$Type == "Random"], dists_clean$piNpiS[dists_clean$Type == "LIC NLR"], dists_clean$piNpiS[dists_clean$Type == "HIC NLR"]), method = "bonferroni")

# ============================
# dN/dS 
# ============================
# Add the gene status data
dNdSdists <- merge(dNdSdists_raw, geneevo, by = c("Ortholog", "Type", "Species")) 

# Add the alternative dN/dS (Ka/Ks) estimation
dNdSdists <- merge(dNdSdists, dplyr::select(kaksdists, Ortholog, Species, Ka, Ks, KaKs, p_val, Divergence_Time), by = c("Ortholog", "Species"))


# Recode the strain names with the species
dNdSdists$Species2 <- recode_factor(dNdSdists$Species,
                                    "pauciseta" = "P. pauciseta",
                                    "comata" = "P. comata",
                                    "bellae-mahoneyi" = "P. bellae-mahoneyi",
                                    "pseudoanserina" = "P. pseudoanserina",
                                    "pseudopauciseta" = "P. pseudopauciseta",
                                    "pseudocomata" = "P. pseudocomata")

# Fix the levels for plotting
dNdSdists$Type <- factor(dNdSdists$Type, levels = c("Random", "LIC NLR", "HIC NLR"))

# Scheirer-Ray-Hare test, a non-parameter alternative to multi-factorial ANOVA analyses
# rcompanion::scheirerRayHare(dNdS ~ Species + Conservation, data = dNdSdists) # Species are not different
# rcompanion::scheirerRayHare(dNdS ~ Type + Conservation, data = dNdSdists) 

# And the ka/ks version? --> very similar
# rcompanion::scheirerRayHare(KaKs ~ Species + Conservation, data = dNdSdists) # Species are not different
# rcompanion::scheirerRayHare(KaKs ~ Type + Conservation, data = dNdSdists) # used in the paper (the interaction is significant!)

## --- Ka/Ks from KaKs_calculator MA model vs egglib (NG model)  ---
# There are couple of cases where KaKs_calculator gives KaKs = 50, because there 
# are no synonymous sites and they get corrected as tiiiny quantities. In the egglib
# version they are just NA, which I think it's fair given the few sites we have.

ggplot(dNdSdists %>% filter(KaKs < 10), aes(dNdS, KaKs)) + 
  geom_abline(colour = "gray") +
  geom_point() +
  geom_point(data = dNdSdists %>% filter(KaKs < 10, p_val < 0.05), colour = "red", alpha = 0.3) +
  ylab("Ka/Ks from KaKs_calculator MA model") +
  xlab("dN/dS from egglib NG model") 
  # geom_smooth(method='lm', formula= y~x) +
  # stat_cor(method = "pearson")

ggplot(dNdSdists %>% filter(KaKs < 10, Divergence_Time < 0.5), aes(KaKs, Divergence_Time)) +
  geom_point() +
    facet_grid(Conservation ~ .)

# The conservation state does have an influence on dN/dS. Hence, for the main figure
# I'll only consider present genes

(dNdSPlot <- ggplot(dNdSdists %>% filter(Conservation == "Conserved", KaKs < 40), aes(x = Type, y = KaKs, colour = Type)) +
    geom_hline(yintercept=1, colour = "gray", linetype='dashed') +
    geom_violin() +
    geom_point(alpha = 0.5) +
    facet_grid(.~Species2) +
    theme_classic() +
    theme(axis.text.x=element_blank(), 
          axis.title.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.position="none",
          strip.text = element_text(face = "bold.italic"),
          panel.border = element_rect(colour = "white", fill = NA),
          axis.ticks.x=element_blank()) +
    xlab("Gene type") +
    scale_colour_manual(values = nlrpalette) + 
    stat_compare_means(comparisons = my_comparisons, 
                       p.adjust.method = "bonferroni", 
                       method = "wilcox.test",
                       label = "p.adj",
                       label.y = c(1.6, 1.8, 2), # Adjust the position of the significance marks
                       symnum.args = symnum.args) +
    scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2)) +
    ylab(expression(d[N]/d[S])) +
    coord_cartesian(ylim = c(0, 2.2)) )

## --- dNdS version ---
ggplot(dNdSdists %>% filter(Conservation == "Conserved"), aes(x = Type, y = dNdS, colour = Type)) +
  geom_hline(yintercept=1, colour = "gray", linetype='dashed') +
  geom_violin() +
  geom_point(alpha = 0.5) +
  facet_grid(.~Species2) +
  theme_classic() +
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold.italic"),
        panel.border = element_rect(colour = "white", fill = NA),
        axis.ticks.x=element_blank()) +
  xlab("Gene type") +
  scale_colour_manual(values = nlrpalette) + 
  stat_compare_means(comparisons = my_comparisons, 
                     p.adjust.method = "bonferroni", 
                     method = "wilcox.test",
                     label = "p.adj",
                     label.y = c(1.6, 1.8, 2), # Adjust the position of the significance marks
                     symnum.args = symnum.args) +
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2)) +
  ylab(expression(d[N]/d[S])) +
  coord_cartesian(ylim = c(0, 2.2))
# -----

# dNdSdists %>% group_by(Type, Conservation, Species) %>% summarize(meandNdS = mean(KaKs, na.rm = TRUE))

## Are the distance distributions different within species?
krusPhyloDist <- dNdSdists %>% filter(Conservation == "Conserved") %>% dplyr::group_by(Species) %>% dplyr::summarize(kruskal_p = kruskal.test(KaKs ~ Type)$p.value)

# If Kruskal-Wallis is significant, perform Dunn's test for pairwise comparisons
dunntest_results_spp <- data.frame()
epitatenames <- c("pauciseta", "comata", "bellae-mahoneyi", "pseudoanserina", "pseudopauciseta", "pseudocomata")
for (species in epitatenames) {
  cat(species)
  pvals <- dunn.test(dNdSdists$KaKs[dNdSdists$Species == species],
                     dNdSdists$Type[dNdSdists$Species == species],
                     method = "bonferroni", list = TRUE)$P
  # method = "bh", list = TRUE)$P # similar results
  dunntest_results_spp <- rbind(dunntest_results_spp, data.frame(Species = species, Comparison = c("NLR-NLR_HIC", "NLR-Random", "NLR_HIC-Random"), pval = pvals))
}

# Only pauciseta is different from the plotting version...
dunntest_results_spp %>% mutate(Significant = pval < 1/18)

## For supplementary materials
(dNdSPlot_sup <- ggplot(dNdSdists %>% filter(Conservation != "Conserved"), aes(x = Type, y = KaKs, colour = Type)) +
    geom_hline(yintercept=1, colour = "gray", linetype='dashed') +
    geom_violin() +
    geom_point(alpha = 0.5) +
    facet_grid(.~Species2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.position="none",
          strip.text = element_text(face = "bold.italic"),
          panel.border = element_rect(colour = "white", fill = NA)) +
    xlab("Gene type") +
    scale_colour_manual(values = nlrpalette) + 
    stat_compare_means(comparisons = my_comparisons, 
                       p.adjust.method = "bonferroni", 
                       method = "wilcox.test",
                       label = "p.adj", 
                       label.y = c(1.6, 1.8, 2), # Adjust the position of the significance marks
                       symnum.args = symnum.args) +
    scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2)) +
    ylab(expression(d[N]/d[S])) +
    coord_cartesian(ylim = c(0, 2.2)) )

# Nothing is significant
# I'm guessing stat_compare_means() is comparing all the means and then 
# correcting for all 6 species comparisons.

### Save plot
ggsave(file = figS_dndS, plot = dNdSPlot_sup, width = 10, height = 3)

# ============================
# Phylogenetic distances
# ============================
# Rename to match the species
speciesnames <- c("P. pauciseta", "P. comata", "P. bellae-mahoneyi", "P. pseudoanserina", "P. pseudopauciseta", "P. pseudocomata")

# Make it a long format
phylodists_lf <- tidyr::gather(phylodists, Species, Distance, "pauciseta":"pseudocomata", factor_key=TRUE)

# Fix bellae-mahoneyi
phylodists_lf$Species <- as.character(phylodists_lf$Species) # it doesn't work if I let it have Species as factor
phylodists_lf$Species[phylodists_lf$Species == "bellae.mahoneyi"] <- "bellae-mahoneyi"

# Add the gene status data
phylodists_lf <- merge(phylodists_lf, geneevo, by = c("Ortholog", "Species"))

# Recode the strain names with the species
phylodists_lf$Species2 <- recode_factor(phylodists_lf$Species,
                                        "pauciseta" = "P. pauciseta",
                                        "comata" = "P. comata",
                                        "bellae-mahoneyi" = "P. bellae-mahoneyi",
                                        "pseudoanserina" = "P. pseudoanserina",
                                        "pseudopauciseta" = "P. pseudopauciseta",
                                        "pseudocomata" = "P. pseudocomata")


# Fix the levels for plotting
phylodists_lf$Type <- factor(phylodists_lf$Type, levels = c("Random", "LIC NLR", "HIC NLR"))

## Are the distance distributions different within species?
kruskal.test(Distance ~ Type, data = phylodists_lf) # but here it's putting all spp together
krusPhyloDist <- phylodists_lf %>% dplyr::group_by(Species) %>% dplyr::summarize(kruskal_p = kruskal.test(Distance ~ Type)$p.value)
# The Kruskal-Wallis test is used for comparing ordinal or non-Normal variables for more than two groups, and is a generalisation of the Mann-Whitney U test.

# Scheirer-Ray-Hare test, a non-parameter alternative to multi-factorial ANOVA analyses
# rcompanion::scheirerRayHare(Distance ~ Type + Species, data = phylodists_lf) # oh there is also an effect of species, maybe because pauci is so introgressed?
# rcompanion::scheirerRayHare(Distance ~ Type + Species, data = phylodists_lf %>% filter(Species2 != "P. pauciseta")) # it certainly has an effect!

# Make df to put a label the horizontal line
dat_text <- data.frame(
  label = "P. anserina \u03C0", # pi in unicode
  Type = "",
  Species2   = "P. pseudocomata",
  Distance   = "" )

# Fix the levels for plotting, again!
dat_text$Species2 <- factor(dat_text$Species2, levels = speciesnames)

( phylodistsPlot <- ggplot(phylodists_lf, aes(x = Type, y = log(Distance), colour = Type)) +
    geom_hline(yintercept = log(0.000497), linetype='dashed', colour = "gray") +
    geom_violin() +
    geom_point(alpha = 0.5) +
    facet_grid(.~Species2) +
    theme_classic() +
    theme(axis.text.x=element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.position="none",
          strip.text = element_blank(),
          panel.border = element_rect(colour = "white", fill = NA),
          axis.ticks.x=element_blank()) +
    xlab("Gene type") +
    scale_colour_manual(values = nlrpalette) + 
    stat_compare_means(comparisons = my_comparisons, 
                       p.adjust.method = "bonferroni", 
                       method = "wilcox.test",
                       label = "p.adj",
                       symnum.args = symnum.args) + # Add significance levels
    coord_cartesian(ylim = c(-13, 5)) +
    geom_text(data = dat_text, mapping = aes(x = 0.25, y = log(0.00005), label = label, parse = T), 
              fontface = "italic",
              hjust = -0.1, vjust = -1, size = 3) )

## *** Figure 4
# piNpiSfigures <- plot_grid(piviolins_ecd_HIC, beziplot, GCviolins_HIC, TEviolins_ecd_HIC, labels = c('a', 'b', 'c', "d"), nrow = 2, align = c("hv"), axis = "tblr")
piNpiSfigures <- plot_grid(TEviolins_ecd_HIC, RIPmuts, piviolins_ecd_HIC, piNpiSviolins, 
                           labels = c('a', 'b', 'c', "d"), 
                           label_size = 17,
                           nrow = 2, align = c("h"), 
                           axis = "tblr") 

## Put all the plots together
(distfigures <- plot_grid(piNpiSfigures, dNdSPlot, phylodistsPlot, 
                          labels = c('', 'e', 'f'), 
                          label_size = 17,
                          nrow = 3, 
                          rel_heights = c(4,1.5, 1.5) ))

## Finally save the full figure
ggsave(file = figM_dist, plot = distfigures, width = 9, height = 11) # png
# ggsave(file = figM_dist, plot = distfigures, width = 9, height = 11, device = cairo_pdf) # I need cairo to get the greek letter pi printed

## Save supplementary figures
ggsave(file = figS_RIPmuts, plot = RIPsup, width = 4, height = 10)

ggsave(file = figS_SupDivStats, 
       plot = plot_grid(piviolins_ecd_HIC_piN, piSviolinsDegrade, labels = c('a', 'b'), nrow = 1, axis = "b", align = "h"), 
       width = 9, height = 4)



