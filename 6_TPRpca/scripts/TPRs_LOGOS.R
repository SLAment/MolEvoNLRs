#!/usr/bin/env Rscript

### LOGOS of the WD40 repeats in the HNWD family
#############################################################################

# =======================================
# Sandra Lorena Ament Velasquez
# 2024-06-18
# =======================================

# ============================
# Load the necessary libraries
# ============================
library(ggplot2)
library(dplyr)
library(ggseqlogo)

# ============================
# Inputs and outputs
# ============================
# Snakemake Input
Pa_2_10340 <- snakemake@input$Pa_2_10340
Pa_3_8170 <- snakemake@input$Pa_3_8170
Pa_6_7950 <- snakemake@input$Pa_6_7950
Pa_6_8850 <- snakemake@input$Pa_6_8850
Pa_7_3550 <- snakemake@input$Pa_7_3550
PaFnt1 <- snakemake@input$PaFnt1
PaPnt1 <- snakemake@input$PaPnt1

# Snakemake Output
eachgene <- snakemake@output$eachgene

# # Local Input
# hete <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_het-e_aa.fa"
# hetd <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_het-d_aa.fa"
# hetr <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_het-r_aa.fa"
# hnwd1 <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_hnwd1_aa.fa"
# hnwd3 <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_hnwd3_aa.fa"
# nwd1 <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_nwd1_aa.fa"
# nwd2 <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_nwd2_aa.fa"
# nwd3 <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_nwd3_aa.fa"
# nwd5 <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_nwd5_aa.fa"
# nwd6 <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_nwd6_aa.fa"
# nwdp2 <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/alignments/PerSpecies/Protein/WDrepeats_10-11-12-14-30-32-39_anserina_nwdp-2_aa.fa"

# Local output
# eachgene <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/01_NLRpca/results/WD40_LOGOS.pdf"
# ============================
# Functions 
# ============================
## Read the fasta into a vector character for ggseqlogo
fasta2vector <- function(fasta_file){
  fasta_contents <- readLines(fasta_file)
  # Initialize variables
  sequences <- c()
  current_sequence <- ""
  
  # Parse the FASTA file contents
  for (line in fasta_contents) {
    if (startsWith(line, ">")) {
      # If a sequence was being read, add it to the list
      if (current_sequence != "") {
        sequences <- c(sequences, current_sequence)
        current_sequence <- ""
      }
    } else {
      # Concatenate lines to form the sequence
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  # Add the last sequence read to the list
  if (current_sequence != "") {
    sequences <- c(sequences, current_sequence)
  }
  return(sequences)
}


# # Print the character vector to verify
# hetd <- fasta2vector(fasta_file)
# 
# ggseqlogo( hetd, method = 'bits' ) # default
# ggseqlogo( hetd, method = 'prob' )
# ggseqlogo( hetd, col_scheme='chemistry' ) # For amino acids you can pick chemistry, hydrophobicity, clustalx, taylor
# 
# ggseqlogo( hetd, col_scheme='chemistry' ) + theme(legend.position="none") # Remove legend

# ============================
# Make Logos 
# ============================
# Make a list with all the genes
allseqs <- list(Pa_2_10340 = fasta2vector(Pa_2_10340), 
                Pa_3_8170 = fasta2vector(Pa_3_8170),
                "Pa_6_7950+Pa_6_7955" = fasta2vector(Pa_6_7950),
                "Pa_6_8850+Pa_6_8860" = fasta2vector(Pa_6_8850), 
                Pa_7_3550 = fasta2vector(Pa_7_3550), 
                PaFnt1 = fasta2vector(PaFnt1),
                PaPnt1 = fasta2vector(PaPnt1)
                )

anserina <- as.vector(unlist(allseqs))
eachgene_plot <- ggseqlogo(c(allseqs, list(All = anserina)), ncol = 1) + 
  # theme(strip.text = element_text(size = 25, face = "italic"))
  theme(strip.text = element_text(size = 25), 
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

ggsave(plot = eachgene_plot, 
       filename = eachgene, 
       width = 12, height = 15)
