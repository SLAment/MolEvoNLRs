#!/usr/bin/env Rscript

### PaIntrogressionExplore: Introgression tracks between the three *Podospora* species
#############################################################################
#
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2024/04/10
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(PopGenome)
library(dplyr, warn.conflicts = FALSE)
# library(ggplot2)
# ============================

# ============================
# Reading the data
# ============================

print("*** Reading data ***")
# Get the vcf directory
vcffolder <- dirname(snakemake@input$vcf) # the location of the vcf

# Window sizes after removing the missing data
winlens <- read.csv(snakemake@input$winlens)

# Variables
WIDTH <- snakemake@params$width
JUMP <- snakemake@params$jump
STRAIN1 <- snakemake@params$strain1
STRAIN2 <- snakemake@params$strain2
CHR <- snakemake@params$chr

# Output
outputname <- snakemake@output$df
# -----

###### ---- LOCAL ------
# Get the vcf directory
# vcffolder <- dirname("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/06_GenomeIdentity2024/popgenome/chr4/PodoPop-vs-Podan2-snps-NoTEs-gatkPASS-miss1-withNA-chr4.vcf")

# # Window sizes after removing the missing data
# winlens <- read.csv("/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/06_GenomeIdentity2024/filter/winlens_chr4_PaSp-CBS237.71m.txt")

# # Variables
# WIDTH <- 10000
# JUMP <- 10000
# STRAIN1 <- "PaSp"
# STRAIN2 <- "CBS237.71m"
# CHR <- "chromosome_4"

# # Output
# outputname <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/06_GenomeIdentity2024/reports/Dxy_chr4_1kb.csv"

# ============================
# Analysis
# ============================
cat("*** Reading VCF file ***\n")
GENOME.class <- readData(vcffolder, format="VCF", include.unknown = TRUE)

cat("*** Setting populations ***\n")
poplist <- list(STRAIN1,STRAIN2)
GENOME.class <- set.populations(GENOME.class, poplist, diploid = FALSE) # the global poplist

cat("*** Making sliding windows ***\n")
GENOME.class.slide <- sliding.window.transform(GENOME.class, width = WIDTH, jump = JUMP, type=2) # 1 scan only biallelic positions (SNPs), 2 scan the genome

# How many windows?
GENOME.class.slide@region.names %>% length() # 379

# The slot GENOME.class@region.names will store the genomic regions of each window
# as a character string. To convert those strings into a genomic numeric
# position we can apply the following
genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
  split <- strsplit(x," ")[[1]][c(1,3)]
  val <- mean(as.numeric(split))
  return(val)
})

cat("*** Dxy ***\n")
GENOME.class.slide <- diversity.stats.between(GENOME.class.slide)
dxypop <- data.frame(GENOME.class.slide@nuc.diversity.between/winlens$winlen) # Between pop

# Make a dataframe 
pwdivergence <- data.frame(Strain1 = poplist[[1]], Strain2 = poplist[[2]], 
                           position = genome.pos, 
                           winlen = winlens$winlen,
                           chr = CHR,
                           dxy = dxypop$pop1.pop2, 
                           row.names = NULL) %>%
  mutate(start = position - ((WIDTH/2)-0.5), end = position + ((WIDTH/2)-0.5))

# Write it into a file
write.table(pwdivergence, sep=",", file = outputname, row.names = FALSE, col.names = FALSE)
