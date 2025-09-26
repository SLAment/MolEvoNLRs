#!/usr/bin/env Rscript

### PaIntrogressionProportions
#############################################################################
# What proportion of windows are introgressed from other species into P. anserina strains?
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2024/04/12
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(cowplot)
library(patchwork)
# ============================

# ============================
# Input files
# ============================
## Input
# Only chromosome 4
DxyWins <- read.csv(snakemake@input$dxy)
PopData <- read.table(snakemake@input$indvs, header=FALSE)
# Output
diverfullfile <- snakemake@output$plot

### Local
# setwd("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Manuscripts/11_MolEvoINWDs/GitHub/2_GenomeIdentity2024")
# DxyWins <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/06_GenomeIdentity2024/results/Dxy_chr4_2kb.csv.gz")
# PopData <- read.table("data/IndividualsIDsSpp.txt", header=FALSE)
## Output
# diverfullfile <- "figures/DxySppComplex_chr4.pdf"

PI <- 0.000497 # pi of P. anserina from Ament-Velásquez et al. (2022) NEE
minlen = 1000
# ============================
# Reading the data
# ============================
names(PopData) <- c("Strain", "Species")
# Metadata
PopDataStrain2 <- PopData %>% filter(Species != "anserina")
names(PopDataStrain2) <- c("Strain2", "Species")

anserinastrains <- PopData %>% filter(Species == "anserina") %>% .$Strain
# ============================
# Analysis
# ============================
thistitle = "Chromosome 4"; thischr = "chromosome_4"; legend = "legend"

# Add species info
pwdivergence <- merge(DxyWins, PopDataStrain2)
# Mark windows with fewer sites than the minimum as missing dagta
pwdivergence$dxy[which(pwdivergence$winlen <= minlen )] <- NA

# I calculated pauciseta vs comata too, but let's remove it for clarity
pwdivergence <- pwdivergence %>% filter(!(Strain1 %in% c("CBS237.71m", "CBS333.63p", "CBS451.62p") & Strain2 %in% c("PcTdp", "PcWa131m", "PcWa132p", "PcWa133m", "PcWa139m")))

# https://rgbacolorpicker.com/rgba-to-hex
farg <- c("pseudoanserina" = "#cc79a7", 
             "pseudopauciseta" = "#56b4e9",
             "pseudocomata" = "#e69f00",
             "bellae-mahoneyi" = "#f0e442",
             "comata" = "#d45500",
             "pauciseta" = "#009e73")

# Give an order to the Species levels and re-order the dataframe to follow that order
pwdivergence <- pwdivergence %>% mutate(Species = factor(Species, levels = c("pauciseta", "comata", "bellae-mahoneyi", "pseudoanserina", "pseudopauciseta", "pseudocomata")) ) %>% arrange(Species)
# Use the order in the dataframe to give levels to Strain2
pwdivergence <- pwdivergence %>% mutate(Strain2 = factor(Strain2, levels = unique(Strain2)))

# What is the distribution of divergence?
dxydistributions <- ggplot(pwdivergence %>% 
         filter(Strain1 %in% c("PaSp", "PaYp"), dxy < 0.1), # Some representative strains
       aes(x = dxy, fill = Species)) +
  geom_histogram() + 
  facet_grid(Strain2~Strain1) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust=1),
        # legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  # xlab(expression("Window pairwise-nucleotide differences ("*D['XY']*")") ) +
  xlab("Window pairwise-nucleotide differences") +
  ylab("Count") +
  scale_fill_manual(values= farg)

## --- What is the distribution of median divergence in the P. anserina collection? ----
medianDxy <- pwdivergence %>% group_by(Strain1, Strain2, Species) %>% summarise(median_dxy = median(dxy, na.rm = TRUE), n = n())

## Reorder to follow the species order that I like
# https://r-graph-gallery.com/267-reorder-a-variable-in-ggplot2.html
# Arrange following the factor order of Species
medianDxy <- medianDxy %>% 
  mutate(Species = factor(Species, levels=c("pauciseta", "comata", "pseudocomata", "pseudoanserina", "pseudopauciseta", "bellae-mahoneyi"))) %>% 
  arrange(Species)

# Now that it's re-arranged, get the desired order of the strains
SppOrder <- unique(medianDxy$Strain2)
# and fix the levels of the Strain2 column
medianDxy <- medianDxy %>% mutate(Strain2=factor(Strain2, levels=SppOrder))  # This trick update the factor levels

medianDxyPlot <- ggplot(medianDxy, aes(x = Strain2, y = median_dxy, colour = Species, fill = Species)) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") +
  # ylab( expression(D['XY']*" between "~italic("P. anserina")~" and other species")) +
  ylab("Median pairwise-nucleotide differences") +
  geom_jitter(width = 0.25, alpha = 0.5) +
  scale_color_manual(values= farg, guide=legend) +
  ylim(0, 0.017) +
  geom_hline(yintercept=PI, linetype="dashed", color = "gray50")
  # geom_violin()

## What proportion of windows are below the P.anserina pi?
IntroMedianDxy_raw <- pwdivergence %>% 
  filter(dxy <= PI) %>% 
  group_by(Strain1, Strain2, Species) %>% summarise(median_dxy_intro = median(dxy, na.rm = TRUE), n_intro = n())

## Many strains got no windows with introgression so I have to put them back
# All windows
allMedianDxy <- pwdivergence %>%
  group_by(Strain1, Strain2, Species) %>% summarise(median_dxy = median(dxy, na.rm = TRUE), n = n())
# Find the difference
missingstrains <- dplyr::setdiff(allMedianDxy[,c(1:3)], IntroMedianDxy_raw[,c(1:3)]) %>% 
  mutate(median_dxy_intro = NA, n_intro = 0)
# Incorporate them
IntroMedianDxy_raw <- rbind(IntroMedianDxy_raw, missingstrains)

# Add the total number of windows and calculate the fraction introgressed
IntroMedianDxy <- merge(IntroMedianDxy_raw, medianDxy) %>% mutate(PerIntrogression = (n_intro/n)*100 )

# Rearrange the order for plotting
IntroMedianDxy <- IntroMedianDxy %>% mutate(Strain2=factor(Strain2, levels=SppOrder))  # This trick update the factor levels

IdenticalDxyPlot <- ggplot(IntroMedianDxy %>% filter(Strain1 != "PcTdp") , aes(x = Strain2, y = PerIntrogression, colour = Species)) + 
  ylab("Proportion of identical 2kb windows (%)") + xlab("Strain") +
  geom_jitter(alpha = 0.5, width = 0.1, height = 0) +
  scale_color_manual(values= farg, guide=legend) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1), legend.position="none")
  #ylim(0, 10)

IntroMedianDxy %>% group_by(Species) %>% 
  summarise(max_introgression = max(PerIntrogression), min_introgressed = min(PerIntrogression))
IntroMedianDxy %>% filter(PerIntrogression > 5)

diverfullplot <- cowplot::plot_grid(dxydistributions, (medianDxyPlot / IdenticalDxyPlot), ncol=2, rel_widths = c(1,0.7), labels = c("a", "b"))

ggsave(plot = diverfullplot, 
       filename = diverfullfile, 
       width = 10, height = 8)

# P. pauciseta is likely the sister species of P. anserina, which could lead to some windows
# being identical because of very strong purifying selection or just by chance. However, the
# variance in the proportion of identical windows is quite large as compared to the other species.
# The proportion is also different for the different pauciseta and comata strains, which implies that
# those strains themselves have different degrees of introgression compared to their intraspp peers.
# In other words, PcTdp, PaWa132p and PaWa133m have more introgressed tracks than the other comata strains.

