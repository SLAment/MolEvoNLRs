#!/usr/bin/env Rscript

### EnrichOddRatios_alla: Display all the results of an enrichment analysis
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2024/07/29
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr)
library(ggplot2)
library(cowplot)
# library(forcats) # for fct_rev()
# ============================

# ============================
# Reading the data
# ============================
## Snakemake
# Read enrichment file: it gives some warnings
enrichment_both <- read.table(snakemake@input$both, head = TRUE, sep = "\t", quote = "")
enrichment_pau <- read.table(snakemake@input$pau, head = TRUE, sep = "\t", quote = "")
enrichment_com <- read.table(snakemake@input$com, head = TRUE, sep = "\t", quote = "")
enrichment_taj <- read.table(snakemake@input$taj, head = TRUE, sep = "\t", quote = "")
enrichment_miss <- read.table(snakemake@input$miss, head = TRUE, sep = "\t", quote = "")
outplot <- snakemake@output$plot
TOP <- snakemake@params$top

# ============================
# Functions
# ============================
enrichplotfun <- function(enrichment, legend = FALSE, ptitle = "Introgression regions", top = TOP, ymax = 75){
  # Remove the cases where there are less than 4 occurrences in the whole genome (to avoid Odd_ratios == Inf)
  enrichment_small <- enrichment %>% mutate(total_Count = Subset_Count + Rest_Count) %>% filter(total_Count > 3) %>% head(n = top) %>% mutate(DescPFAM = paste0(Description, "\n(", PFAM, ") " )) 
  
  # Set some boundaries for the plot
  minratio <- min(enrichment_small$Odds_ratio, 0)
  # maxratio <- max(enrichment_small$Odds_ratio) + 10 # the 10 extra is to allow the P-value to be fully printed within the plot
  maxratio <- ymax
  
  enrichplot <- ggplot(enrichment_small, aes(x = Odds_ratio, y = reorder(DescPFAM, Odds_ratio), fill = Type)) +
    # ggplot(enrichment_small, aes(x = Odds_ratio, y = fct_rev(reorder(DescPFAM, P_Value_Adj)), fill = Type)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = ifelse(P_Value_Adj < 0.05, 
                                 paste0("bolditalic('p=", format(round(P_Value_Adj, 4), nsmall = 4), "')"), # format() ensures that they are displayed with three decimal places.
                                 paste0("italic('p=", format(round(P_Value_Adj, 4), nsmall = 4), "')"))), 
              parse = TRUE, hjust = -0.1, vjust = 0.3, color = "black", size = 3.5) +
    theme_minimal() +
    labs(x = "Odds ratio", y = "") +
    # scale_fill_manual(values = c("gray36","gray")) +
    scale_fill_brewer('', palette = "Paired", direction = -1) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      plot.title.position = 'plot', #  use the entire plot as a reference for the centering
      plot.title = element_text(hjust = 0.5)
    ) + 
    ggtitle(ptitle) +
    xlim(minratio, maxratio)
  
  if (legend) {
    enrichplot <- enrichplot + theme(legend.position = c(0.75, 0.55))
  } else {
    enrichplot <- enrichplot + theme(legend.position = "none")
  }
  
  # format the labels using paste0 to include HTML tags for bold text where needed, 
  # and geom_text will need parse=TRUE to correctly interpret the HTML tags.
  
  return(enrichplot)
}
# ============================
# Plot
# ============================

allplots <- plot_grid(enrichplotfun(enrichment_both, ptitle = "Introgression regions", legend = TRUE), 
                      enrichplotfun(enrichment_pau, ptitle = expression("Introgressed from"~italic("P. pauciseta"))),
                      enrichplotfun(enrichment_taj, ptitle = expression("Tajima's D">=2~"regions in the Wageningen Collection"), ymax = 100),
                      enrichplotfun(enrichment_com, ptitle = expression("Introgressed from"~italic("P. comata"))),
                      enrichplotfun(enrichment_miss , ptitle = "Regions of missing data"),
                      ncol = 2, align = "hv")

ggsave(plot = allplots, 
       filename = outplot, 
       width = 15, height = 12)


