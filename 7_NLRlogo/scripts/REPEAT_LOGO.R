#!/usr/bin/env Rscript

### REPEAT_LOGO.R - Produce a logo from input alignment of repeats
#############################################################################
# =======================================
# Sandra Lorena Ament Velasquez
# 2024-08-13
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
repeats <- snakemake@input$fasta
genename <- snakemake@params$gene

# Snakemake Output
logoplot <- snakemake@output$logo

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

# ============================
# Make Logos 
# ============================

# Print the character vector to verify
heat <- fasta2vector(repeats)

# To give it a title, make it a list
allseqs <- list()
allseqs[[genename]]  <- fasta2vector(repeats) # So that I can use the variable to assign a name
# allseqs <- list(genename = fasta2vector(repeats))

heatplot <- ggseqlogo(allseqs) + 
  theme(strip.text = element_text(size = 25), 
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

ggsave(plot = heatplot, 
       filename = logoplot, 
       width = 11, height = 3)


# ggseqlogo( heat, method = 'bits' ) # default
# ggseqlogo( heat, method = 'prob' )
# ggseqlogo( heat, col_scheme='chemistry' ) # For amino acids you can pick chemistry, hydrophobicity, clustalx, taylor
# ggseqlogo( heat, col_scheme='chemistry' ) + theme(legend.position="none") # Remove legend