#!/usr/bin/env Rscript

### PCA of the TPR repeats in P. anserina NB-ARC NLRs
#############################################################################

# =======================================
# Sandra Lorena Ament Velasquez
# 2025-03-27
# =======================================
# https://arftrhmn.net/quick-pca-analysis-from-sequence-alignment-data-in-r/
# https://stats.stackexchange.com/questions/222/what-are-principal-component-scores
# https://thomasadventure.blog/posts/turning-your-ggplot2-code-into-a-function/
# https://builtin.com/data-science/step-step-explanation-principal-component-analysis

# ============================
# Load the necessary libraries
# ============================
library(adegenet) # It loads the package ade4 too
library(ggplot2)
library(dplyr)
library(ape)
library(patchwork)
library(cowplot)

# ============================
# Read input data
# ============================
# meta file that contains sequence features that we are interested to visualize 
# in the PCA. The meta file must have a column containing isolate ids that 
# matched fasta id as well.

## Input
meta <- read.table(snakemake@input$meta, sep='\t', header = T)
dnaseqs <- ape::read.dna(snakemake@input$fasta, format = "fasta")
# meta <- read.table("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/09_TPRpca/alignments/REPEATrepeats_10-14-15-16-17-19-26-27-33_metadata.txt", sep='\t', header = T)
# dnaseqs <- ape::read.dna("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/3_NLRdescription/09_TPRpca/alignments/PerSpecies/REPEATrepeats_10-14-15-16-17-19-26-27-33_anserina.fa", format = "fasta")

## Output
paperfigfile <- snakemake@output$paper

# ============================
# Define some functions
# ============================
## General function to make a PCA dataframe for plotting
makePCAfromAlignment <- function(dnaseqs, nf = 10){
  # Make genind object from the DNAbin object
  multisites <- adegenet::DNAbin2genind(dnaseqs, polyThres=0)
  
  # Deal with missing data before making PCA
  # sum(is.na(multisites$tab)) 
  X <- adegenet::tab(multisites, freq = TRUE, NA.method = "mean") # replace NA by the mean allele frequencies
  
  # Calculate the actual Principal Component Analysis
  pca <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = nf)
  
  # Make a dataframe from the results
  pca.dataset <- data.frame(Sequence = row.names(pca$li),
                            as.data.frame(pca$li)) # the row coordinates i.e. the principal components
  
  # Re-assign the name of the sequence based on the Sequence variable
  pca.dataset <- base::merge(meta, pca.dataset, by = "Sequence")
  
  return(list(pca.dataset, pca))
}

## Get variance explained by a PC
getvarPCA <- function(pca, axis = 1){
  Eigenvalues <- pca$eig # Note to self: weights = loadings = eigenvalues
  Variance <- Eigenvalues / sum(Eigenvalues) 
  VarAxis <- (100 * signif(Variance[axis], 4)) %>% round(digits=2)
  return(VarAxis)
}

## Make a theme for the PCA
clean_theme <- function() {
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
}

colorgenes = c("Pa_2_10340" = "#f867cdff", 
               "Pa_3_8170" = "#f8766cff", 
               "Pa_6_7950+Pa_6_7955" = "#00bd5cff" , 
               "Pa_6_8850+Pa_6_8860" = "#00a6ffff", 
               "Pa_7_3550" = "#de950fff", 
               "PaFnt1" = "#b88dffff", 
               "PaPnt1" = "#ffd42aff")

shapesgenes = c("Pa_2_10340" = 15, 
                "Pa_3_8170" = 16, 
                "Pa_6_7950+Pa_6_7955" = 17 , 
                "Pa_6_8850+Pa_6_8860" = 18, 
                "Pa_7_3550" = 20, 
                "PaFnt1" = 21, 
                "PaPnt1" = 24)


## Plot different combinations of PCs
# The Switch allows you to turn on and off the shape of the points:
# https://stackoverflow.com/questions/22915337/if-else-condition-in-ggplot-to-add-an-extra-layer
plot4pcas <- function(dnaseqs, classi = Gene, Switch = TRUE){
  pcadf <- makePCAfromAlignment(dnaseqs)[[1]]
  pcaraw <- makePCAfromAlignment(dnaseqs)[[2]]
  
  # # Make a list of shapes (more general but stochastic depending on the dataset)
  # nuclassi <- nrow(unique(pcadf[ deparse(substitute(classi)) ]))
  # # shapevals <- seq(1,nuclassi, 1) # i don't like the first shapes
  # shapevals <- c(seq(15, 15+nuclassi, 1), seq(1,nuclassi, 1))
  
  # Aesthetics
  alphaval <- 0.5
  
  p1 <- ggplot(pcadf, aes(Axis1, Axis2, colour= {{classi}})) +  # The !! is necessary to interpret the variable (unquote)
    geom_point({if(Switch)aes(shape={{classi}})}, size=3, alpha=alphaval) +
    scale_shape_manual(values=shapesgenes) +
    scale_color_manual(values = colorgenes) +
    xlab(paste0("PC1 (", getvarPCA(pcaraw, 1), " %)")) +
    ylab(paste0("PC2 (", getvarPCA(pcaraw, 2), " %)")) +
    theme_bw() +
    clean_theme() +
    theme(legend.position="none")
  
  p2 <- ggplot(pcadf, aes(Axis1, Axis3, colour= {{classi}})) + 
    geom_point({if(Switch)aes(shape={{classi}})}, size=3, alpha=alphaval) +
    scale_shape_manual(values=shapesgenes) +
    scale_color_manual(values = colorgenes) +
    xlab(paste0("PC1 (", getvarPCA(pcaraw, 1), " %)")) +
    ylab(paste0("PC3 (", getvarPCA(pcaraw, 3), " %)")) +
    theme_bw() +
    clean_theme() +
    theme(legend.position="none")
  
  p3 <- ggplot(pcadf, aes(Axis2, Axis3, colour= {{classi}})) + 
    geom_point({if(Switch)aes(shape={{classi}})}, size=3, alpha=alphaval) +
    scale_shape_manual(values=shapesgenes) +
    scale_color_manual(values = colorgenes) +
    xlab(paste0("PC2 (", getvarPCA(pcaraw, 2), " %)")) +
    ylab(paste0("PC3 (", getvarPCA(pcaraw, 3), " %)")) +
    theme_bw() +
    clean_theme() 
  # theme(legend.position="none")
  
  p4 <- ggplot(pcadf, aes(Axis3, Axis4, colour= {{classi}})) + 
    geom_point({if(Switch)aes(shape={{classi}})}, size=3, alpha=alphaval) +
    scale_shape_manual(values=shapesgenes) +
    scale_color_manual(values = colorgenes) +
    xlab(paste0("PC3 (", getvarPCA(pcaraw, 3), " %)")) +
    ylab(paste0("PC4 (", getvarPCA(pcaraw, 4), " %)")) +
    theme_bw() +
    clean_theme() +
    theme(legend.position="none")
  
  return((p1 + p2 + p3))
  # return((p1 + p2) / ( p3 + p4))
}

# ---- Make different PCA combinations ---
allTPRs <- plot4pcas(dnaseqs)

# What's up with PaPnt1 and Pa_7_3550?
suspish <- meta %>% 
  filter(Gene %in% c("PaPnt1", "Pa_7_3550")) %>% 
  filter(Species %in% c("anserina"))
dnaseqs_suspish <- dnaseqs[which(labels(dnaseqs) %in% suspish$Sequence),]

suspishTPRs <- plot4pcas(dnaseqs_suspish)

## Put them together for the paper
paperfig <- plot_grid(allTPRs, suspishTPRs, nrow=2, labels=c('a', 'b'), align = "hv")
ggsave(file = paperfigfile, plot = paperfig, width = 9.5, height = 4.5)

