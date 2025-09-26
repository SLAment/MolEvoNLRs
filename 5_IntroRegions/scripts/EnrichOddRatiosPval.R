#!/usr/bin/env Rscript

### EnrichOddRatiosPval: Display all the results of an enrichment analysis
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
enrichment <- read.table(snakemake@input$bed, head = TRUE, sep = "\t", quote = "")
outplot <- snakemake@output$plot

# ============================
# Plot
# ============================
# Odds ratio calculation example: For the HET domain -- (25/(756-25))/(97/(11265-97))

enrichment_df <- enrichment %>%
  mutate(Gene_number = Subset_Count + Rest_Count) %>% filter(P_Value < 0.05)

# Make some plot limits
maxOdds <- enrichment_df %>% filter(Odds_ratio != Inf) %>% .$Odds_ratio %>% max()
maxPval <- enrichment_df %>% mutate(logPval = -log(P_Value_Adj)) %>% .$logPval %>% max()

# Fix names to make them shorter
enrichment_df$Name <- recode(enrichment_df$Description, 
       'Heterokaryon incompatibility protein (HET)' = "HET", 
       "NACHT domain" = "NACHT",
       "NB-ARC domain" = "NB-ARC",
       "Tetratricopeptide repeat" = "TPR",
       "Fungal N-terminal domain of STAND proteins" = "HeLo-like")

nlrpalette =  c("NLR-associated" = "#1F78B4", "Other" = "gray50")

# Create a new data frame with the top 3 points per Type
top_labels <- enrichment_df %>% filter(P_Value_Adj < 1, Type == "NLR-associated")

(oddspval <- ggplot(enrichment_df %>% filter(P_Value < 0.05, Odds_ratio != Inf),
                    aes(x = Odds_ratio, y = -log(P_Value_Adj), size = Gene_number, colour = Type)) +
    geom_point(alpha = 0.5) + theme_bw() +
    geom_hline(yintercept= -log(0.05), linetype='dashed', color = "gray", size = 0.3) +
    geom_text(data = top_labels, 
              aes(label = Name), 
              show.legend = FALSE, # Don't put a letter a on top of the colored legend
              size = 3, vjust = -0.6, hjust = 0.7) +  # Adjust text position
    scale_colour_manual(values = nlrpalette) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          # legend.position = c(0.6, 0.8), # Adjust the x and y position of the legends inside the plot
          # legend.direction = "horizontal", # Arrange legends side by side
          # legend.box = "vertical" # Make sure both legends share the same horizontal space
          ) +
    coord_cartesian(xlim = c(0, maxOdds + maxOdds*0.05), ylim = c(0, maxPval + maxPval*0.05) ) +
    xlab("Odds ratio") + ylab("-log(adjusted p-values)") +
    guides(size = guide_legend(title="Gene number")) ) 

ggsave(plot = oddspval, 
       filename = outplot, 
       width = 4.5, height = 2.8)
       # width = 5.5, height = 3)
