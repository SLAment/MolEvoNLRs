#!/usr/bin/env Rscript

# PerSample_vcfR_plotter: Plot the coverage distribution of the *P. anserina* samples
#############################################################################

# A modified version of DiversityStats_vcfR_plotter.R to handle a single given
# strain to reduce memory problems. The goal here is the table, not the
# plot.

# =======================================
# Sandra Lorena Ament Velasquez
# 2024-04-09
#############################################################################

# ============================
# Check input
# ============================
vcffile <- snakemake@input$vcf
tablecov <- snakemake@output$table

thissample <- snakemake@params$sample
lquantile = snakemake@params$lowquantile
uquantile = snakemake@params$uppquantile

# vcffile <- "/Users/lorena/Dropbox/VRwork/Analyses/VRpipelines/06_GenomeIdentity2024/temp/PodoPop-vs-Podan2-snps_PaSp.vcf.gz"
# thissample <- "PaSp"
# lquantile <- 0.25
# uquantile <- 0.995
# ============================
# Load the necessary libraries
# ============================
library(vcfR)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
# library(ggplot2)
# library(Hmisc)
# ============================
# Reading the data
# ============================
vcf <- read.vcfR(vcffile, verbose = FALSE)

# ============================
# Prepare the functions
# ============================
# vignette('sequence_coverage')
# It requires 
coolviolins <- function(dp){
  if( require(reshape2) & require(ggplot2) ){
    dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
    dpf <- dpf[ dpf$Depth > 0,]
    p <- ggplot(dpf, aes(x=Sample, y=Depth)) + geom_violin(fill="#C0C0C0", adjust=1.0,
                                                           scale = "count", trim=TRUE)
    p <- p + theme_bw()
    p <- p + theme(axis.title.x = element_blank(), 
                   axis.text.x = element_text(angle = 60, hjust = 1))
    p <- p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")
    p <- p + scale_y_continuous(trans=scales::log2_trans(), breaks=c(1, 10, 100, 800))
    p
    return(p)
  } else {
    message("The packages reshape2 and ggplot2 are required for this example but do not appear
          to be installed. Please use install.packages(c('reshape2', 'ggplot2', 'scales')) if you would
          like to install them.")
  }
}

vcfsubsetbyindv <- function(vcf, indv = c("PaSp")){
  indvcf <- vcf
  indvcf@gt <- vcf@gt[, c("FORMAT", indv)] # Get only the relevant genotypes using the range of the chromosome
  return(indvcf)
}

# ============================
# Analysis
# ============================

## Print violin plots of coverage
dp <- extract.gt(vcf, element="DP", as.numeric = TRUE) # Make data frame of depth of coverage
# ggsave(plotcov, plot = coolviolins(dp), width = 30, height = 10, units = "cm")

# Print a summary table of the coverage distribution per sample
# I chose 0.985 out of trial and error. I tested 0.9-0.99
dpsum <- data.frame(sample = colnames(dp), median = apply(dp, 2, median, na.rm = TRUE), mean = apply(dp, 2, mean, na.rm = TRUE), lowquantile = apply(dp, 2, quantile, probs = lquantile, na.rm = TRUE), upquantile = apply(dp, 2, quantile, probs = uquantile, na.rm = TRUE))

write.table(dpsum, file = tablecov, sep = "\t", row.names=FALSE, quote=FALSE, col.names = FALSE)


