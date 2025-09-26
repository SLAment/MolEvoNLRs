#!/usr/bin/env Rscript

# Exploring the correlation of divergence between Podospora species
#############################################################################
# Version for 2Kb non-overlapping windows
# =======================================
# S. Lorena Ament-Velasquez
# 2024-04-19
#############################################################################


# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(tidyr)

# ============================
# Input files
# ============================
# Input
DxyWinsChr1 <- read.csv(snakemake@input$chr1)
DxyWinsChr2 <- read.csv(snakemake@input$chr2)
DxyWinsChr3 <- read.csv(snakemake@input$chr3)
DxyWinsChr4 <- read.csv(snakemake@input$chr4)
DxyWinsChr5 <- read.csv(snakemake@input$chr5)
DxyWinsChr6 <- read.csv(snakemake@input$chr6)
DxyWinsChr7 <- read.csv(snakemake@input$chr7)

PopData <- read.table(snakemake@input$PopData, header=FALSE)

## Output
DxyPauciComata_chrs <- snakemake@output$corr
DxyPauciComata_genome <- snakemake@output$genome

# # Input
# DxyWinsChr1 <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/results/Dxy_chr1_2kb.csv.gz")
# DxyWinsChr2 <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/results/Dxy_chr2_2kb.csv.gz")
# DxyWinsChr3 <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/results/Dxy_chr3_2kb.csv.gz")
# DxyWinsChr4 <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/results/Dxy_chr4_2kb.csv.gz")
# DxyWinsChr5 <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/results/Dxy_chr5_2kb.csv.gz")
# DxyWinsChr6 <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/results/Dxy_chr6_2kb.csv.gz")
# DxyWinsChr7 <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/results/Dxy_chr7_2kb.csv.gz")

# PopData <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/06_GenomeIdentity2024/data/IndividualsIDsSpp.txt", header=FALSE)
# 
# ## Output
# DxyPauciComata_genome <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/figures/DxyPauciComata_genome.png"
# DxyPauciComata_chrs <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/07a_3SppIntrogression/figures/DxyPauciComata_chrs.png"

PI <- 0.000497
minlen = 1000
# ============================
# Reading the data
# ============================
DxyWins <- rbind(DxyWinsChr1, DxyWinsChr2, DxyWinsChr3, DxyWinsChr4, DxyWinsChr5, DxyWinsChr6, DxyWinsChr7)

names(PopData) <- c("Strain", "Species")
# Metadata
PopDataStrain2 <- PopData %>% filter(Species != "anserina") # remove the anserina strains
names(PopDataStrain2) <- c("Strain2", "Species")

# ============================
# Analysis
# ============================
# Add species info
pwdivergence <- merge(DxyWins, PopDataStrain2)
pwdivergence$dxy <- as.numeric(pwdivergence$dxy) # Make dxy numeric

# Remove windows with too much missing data
pwdivergence$dxy[which(pwdivergence$winlen <= minlen )] <- NA

# Select representative strain comparisons, trying to maximize introgression levels
meanDxy_rep <- pwdivergence %>% 
  filter(Strain2 %in% c("CBS237.71m", "PcTdp"), Strain1 %in% c("PaYp", "CBS455.64m") ) %>% 
  group_by(position, chr, Strain2, Strain1) %>% summarise(mean_dxy = mean(dxy, na.rm = TRUE), n = n())

# Break it into two variables for plotting
meanDxy_rep_corr_pa <- meanDxy_rep %>% filter(Strain2 == "CBS237.71m") %>% mutate(mean_dxy_pauci = mean_dxy) %>% ungroup() %>% select(!c(mean_dxy, Strain2))
meanDxy_rep_corr_co <- meanDxy_rep %>% filter(Strain2 == "PcTdp") %>% mutate(mean_dxy_comata = mean_dxy) %>% ungroup() %>% select(!c(mean_dxy, Strain2))
meanDxy_rep_corr <- merge(meanDxy_rep_corr_pa, meanDxy_rep_corr_co)

# Assign a color
meanDxy_rep_corr <- meanDxy_rep_corr %>% mutate(Dxy_bias = case_when(
  (mean_dxy_pauci <= PI & mean_dxy_comata <= PI) ~ "both",
  mean_dxy_pauci <= PI ~ "pauciseta",
  mean_dxy_comata <= PI ~ "comata",
  TRUE ~ "none")) %>% tidyr::drop_na()

#### By chromosome with R-square plotted
meanDxy_rep_corr_example <- meanDxy_rep_corr %>% filter(Strain1 == "PaYp")
cor.test(meanDxy_rep_corr_example$mean_dxy_comata, meanDxy_rep_corr_example$mean_dxy_pauci, method = "pearson")
cor.test(meanDxy_rep_corr_example$mean_dxy_comata, meanDxy_rep_corr_example$mean_dxy_pauci, method = "spearman")

linearRegre <- lm(mean_dxy_comata ~ mean_dxy_pauci*Strain1, data = meanDxy_rep_corr)
summary(linearRegre)


# Calculate R-squared values for each combination of Strain1 and chr
r_squared_data <- meanDxy_rep_corr %>%
  group_by(Strain1, chr) %>%
  summarize(
    r_squared = summary(lm(mean_dxy_comata ~ mean_dxy_pauci))$r.squared,
  )

# Convert r_squared_data to a dataframe for easier use in ggplot
r_squared_data <- r_squared_data %>%
  ungroup() %>%
  mutate(label = paste("R² =", round(r_squared, 3)))

# Fix the levels for plotting
meanDxy_rep_corr$Dxy_bias <- factor(meanDxy_rep_corr$Dxy_bias, levels = c("pauciseta", "comata", "both", "none"))

# Create the plot
dxycorr_chrs <- ggplot(meanDxy_rep_corr, aes(mean_dxy_pauci, mean_dxy_comata)) + 
  geom_point(alpha = 0.5, aes(colour = Dxy_bias, shape = Dxy_bias)) +
  geom_smooth(method='lm', color = "gray", na.rm = TRUE, data = meanDxy_rep_corr) +
  facet_grid(chr ~ Strain1) + 
  theme_bw() +
  xlab(expression("Mean divergence with"~italic("P. pauciseta")~"(CBS237.71m)")) +
  ylab(expression("Mean divergence with"~italic("P. comata")~"(PcTdp)")) +
  theme(
    strip.background = element_rect(fill="white"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    # legend.position="bottom"
  ) +
  guides(colour = guide_legend(title="Divergence bias"), shape = guide_legend(title="Divergence bias")) +
  scale_color_manual(values=c("pauciseta" = "deepskyblue3", "comata" = "red", "none" = "gray57", "both" = "gold")) +
  geom_text(data = r_squared_data, aes(x = Inf, y = 0.06, label = paste("R² = ", round(r_squared, 3))), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE)

ggsave(plot = dxycorr_chrs, 
       filename = DxyPauciComata_chrs, 
       width = 6, height = 10)

### -------
#### The whole genome

# Calculate R-squared values for each combination of Strain1 and full genome
r_squared_data <- meanDxy_rep_corr %>%
  group_by(Strain1) %>%
  summarize(
    lm_model = list(lm(mean_dxy_comata ~ mean_dxy_pauci)),
    r_squared = summary(lm(mean_dxy_comata ~ mean_dxy_pauci))$r.squared
  )

# Convert r_squared_data to a dataframe for easier use in ggplot
r_squared_data <- r_squared_data %>%
  ungroup() %>%
  mutate(label = paste("R² =", round(r_squared, 3)))


( dxycorr_genome <- ggplot(meanDxy_rep_corr, aes(mean_dxy_pauci, mean_dxy_comata, colour = Dxy_bias)) + 
    geom_point(alpha = 0.3) +
    facet_grid(Strain1~.) + theme_bw() +
    geom_smooth(method='lm', color = "gray48", na.rm = TRUE) +
    xlab(expression("Mean Divergence with"~italic("P. pauciseta")~"(CBS237.71m)")) +
    ylab(expression("Mean Divergence with"~italic("P. comata")~"(PcTdp)")) +
    theme(strip.background = element_rect(fill="white"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    scale_color_manual("Divergence bias", values=c("pauciseta" = "deepskyblue3", "comata" = "red", "none" = "gray87")) +
  geom_text(data = r_squared_data, aes(x = Inf, y = 0.06, label = paste("R² = ", round(r_squared, 3))), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE))

ggsave(plot = dxycorr_genome, 
       filename = DxyPauciComata_genome, 
       width = 5, height = 5)

## Some numbers
# How many windows are there per strain?
meanDxy_rep_corr %>% group_by(Strain1) %>% summarise(n = n()) # 14814 windows of 2kb, so 29.628 Mbp, a big loss for the total of 36 Mbp
nwinds_chrs <- meanDxy_rep_corr %>% group_by(Strain1, chr) %>% summarise(n_chr = n())  # number of windows per chr

# How many windows per chr and per strain?
nwinds_chrs_lowdiv <- meanDxy_rep_corr %>% filter(Dxy_bias == "pauciseta") %>% group_by(Strain1, chr) %>% summarise(n = n())
meanDxy_rep_corr %>% filter(Dxy_bias == "pauciseta") %>% group_by(Strain1) %>% summarise(n = n(), perdiv = n/14814 * 100)
meanDxy_rep_corr %>% filter(Dxy_bias == "comata") %>% group_by(Strain1) %>% summarise(n = n(), perdiv = n/14814 * 100)
merge(nwinds_chrs_lowdiv, nwinds_chrs) %>% mutate(perdiv = (n/n_chr * 100))





