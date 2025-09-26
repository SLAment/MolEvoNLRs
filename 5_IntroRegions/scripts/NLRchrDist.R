#!/usr/bin/env Rscript

# NLRchrDist: Visualizing NLR-associated domains in Podan2
#############################################################################
# Version for 10Kb non-overlapping windows
# =======================================
# S. Lorena Ament-Velasquez
# 2024-08-01
#############################################################################

# ============================
# Load the necessary libraries
# ============================
library(ggplot2)
library(dplyr)
library(scales)
library(ggnewscale)
# ============================
# Input files
# ============================
# Snakemake
introbed <- read.table(snakemake@input$intro, header = FALSE, stringsAsFactors = FALSE) %>% mutate(Region = "Introgressed")
missbed <- read.table(snakemake@input$miss, header = FALSE, stringsAsFactors = FALSE) %>% mutate(Region = "Missing")
gene_data <- read.table(snakemake@input$genes, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
genes_coords <- read.table(snakemake@input$coords, header = TRUE, sep = "\t")

## Output
domainplotfile <- snakemake@output$domain
NLRplotfile <- snakemake@output$NLRs

# ============================
# Process 
# ============================

bed_data <- rbind(introbed, missbed)
colnames(bed_data) <- c("chromosome", "start", "end", "Region")

colnames(gene_data) <- c("chromosome", "start", "end", "Domain", "PFAM", "geneID")

chromosome_lengths <- data.frame(chromosome = c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7"),
                                 length = c(8813524, 5165621, 4137471, 3808395, 4734292, 4264132, 4087160))

centromeres <- data.frame(chromosome = c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7"),
                          start = c(4479205, 236468.5, 675115.5, 1236388.5, 2062808.5, 238150, 3562141.5)) %>% 
  mutate(end = start, Domain = "Centromere", PFAM = "Centromere", geneID = "Centromere", ATPase = "Centromere", Type = "Centromere")

# Define the desired order for the Domain levels
domain_order <- c("HET", "NACHT", "NB-ARC", "ANK", "TPR", "WD40", "Centromere")

## Add a chr number column for convenience
bed_data <- bed_data %>% mutate(chrnum = as.numeric(factor(chromosome)))
chromosome_lengths <- chromosome_lengths %>% mutate(chrnum = as.numeric(factor(chromosome)))
centromeres <- centromeres %>% mutate(chrnum = as.numeric(factor(chromosome)))
genes_coords <- genes_coords %>% mutate(chrnum = as.numeric(factor(chromosome)))

# ============================
# Location of NLR-associated domains
# ============================

# Define the width of the rectangles representing the chromosomes
rect_width <- 0.3

(domainplot <- ggplot() +
    
    # Initialize the base plot with gray rectangles representing the chromosomes
    geom_rect(data = chromosome_lengths, 
              aes(xmin = as.numeric(factor(chromosome)) - rect_width / 2, 
                  xmax = as.numeric(factor(chromosome)) + rect_width / 2, 
                  ymin = 0, 
                  ymax = length), 
              fill = "gray95", color = "gray60", linewidth = 0.2, linetype = "solid") +
    theme_classic() +
    labs(x = "Chromosome", y = "Position") +
    coord_flip() + # Make the chrs horizontal
    scale_x_continuous(breaks = 1:7, labels = paste0("Chr", 1:7)) + # Ensure all chromosome numbers are displayed on the y-axis
    scale_y_continuous(expand = c(0.01, 0), labels = ~ format(.x, scientific = FALSE)) + # Ensure y-axis starts at 0 and ends at max value and remove scientific notation
    
    # Add the introgressed and missing regions
    geom_segment(data = bed_data, 
                 aes(x = as.numeric(factor(chromosome)), 
                     xend = as.numeric(factor(chromosome)), 
                     y = start, yend = end, color = Region), linewidth = 5) +
    scale_color_manual(values= c("Introgressed" = "black", "Missing" = "gray60")) +
    ggnewscale::new_scale("color") + # Geoms added to a plot after this function will use a new scale definition.
    
    # Add the NLR domains
    geom_point(data = gene_data %>% filter(Domain %in% c("HET")), 
               aes(x = as.numeric(factor(chromosome)) + -rect_width+0.06, 
                   y = (start + end) / 2, color = Domain, shape = Domain), 
               size = 1.2, alpha = 0.5) + #, width = 0.03
    geom_point(data = gene_data %>% filter(Domain %in% c("NACHT", "NB-ARC")), 
               aes(x = as.numeric(factor(chromosome)) + -rect_width-0.05, 
                   y = (start + end) / 2, color = Domain, shape = Domain), 
               size = 1.2, alpha = 0.7) +
    # geom_point(data = gene_data %>% filter(Domain %in% c("TPR", "ANK", "WD40")), 
    #            aes(x = as.numeric(factor(chromosome)) + -rect_width-0.18, 
    #                y = (start + end) / 2, color = Domain, shape = Domain), 
    #            size = 1.2, alpha = 0.7) +
    
    # Add the centromeres
    geom_point(data = centromeres,
    aes(x = as.numeric(factor(chromosome)),
        y = start, color = Domain, shape = Domain),
    size = 2, alpha = 1) +
   
    ## Set the colors (the breaks are necessary to force the order!)
    scale_shape_manual("Features", breaks = domain_order, values = c("Centromere" = 19, "ANK" = 17, "HET" = 8, "NACHT" = 15, "NB-ARC" = 24, "TPR" = 1, "WD40" = 23 )) +
    scale_color_manual("Features", breaks = domain_order, values = c("Centromere" = "black", "ANK" = "#E6AB02", "HET" = "#D95F02", "NACHT" = "#E7298A", "NB-ARC" = "#66A61E", "TPR" = "#7570B3", "WD40" = "#1B9E77" )) +
    # Position the legend inside the plot
    labs(x = NULL, y = "Position (bp)") +
    guides(color = guide_legend(order = 0), shape = guide_legend(order = 0)) +
    theme(legend.position = c(0.8, 0.58), # Adjust the position as needed
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.5)))  )# Add a background to the legend


ggsave(plot = domainplot, 
       filename = domainplotfile, 
       width = 8, height = 4)

# ============================
### Location of the studied NLR
# ============================

NLRplot <- (ggplot() +
              # Initialize the base plot with gray rectangles representing the chromosomes
              geom_rect(data = chromosome_lengths, 
                        aes(xmin = chrnum - rect_width / 2, 
                            xmax = chrnum + rect_width / 2, 
                            ymin = 0, 
                            ymax = length), 
                        fill = "gray95", color = "gray60", linewidth = 0.2, linetype = "solid") +
              theme_classic() +
              labs(x = "Chromosome", y = "Position") +
              coord_flip() + # Make the chrs horizontal
              scale_x_continuous(breaks = 1:7, labels = paste0("Chr", 1:7)) + # Ensure all chromosome numbers are displayed on the y-axis
              scale_y_continuous(expand = c(0.01, 0), labels = ~ format(.x, scientific = FALSE)) + # Ensure y-axis starts at 0 and ends at max value and remove scientific notation
              
              # Add the introgressed and missing regions
              geom_segment(data = bed_data, 
                           aes(x = chrnum, 
                               xend = chrnum, 
                               y = start, yend = end, color = Region), linewidth = 5) +
              scale_color_manual(values= c("Introgressed" = "black", "Missing" = "gray60")) +
              ggnewscale::new_scale("color") + # Geoms added to a plot after this function will use a new scale definition.
              
              # Add the gene locations based on their ATPase/lack of
              geom_point(data = genes_coords %>% filter(ATPase %in% c("NACHT")),
                         aes(x = chrnum + -rect_width+0.06, 
                             y = (start + end) / 2, color = Type, shape = ATPase), 
                         size = 2, alpha = 0.8) + #, width = 0.03
              geom_point(data = genes_coords %>% filter(ATPase %in% c("NB-ARC")),
                         aes(x = chrnum + -rect_width-0.09, 
                             y = (start + end) / 2, color = Type, shape = ATPase), 
                         size = 2, alpha = 0.8) + #, width = 0.03
              # geom_point(data = genes_coords %>% filter(ATPase %in% c("Random")),
              #            aes(x = chrnum + -rect_width-0.18, 
              #                y = (start + end) / 2, color = Type, shape = ATPase), 
              #            size = 1.2, alpha = 0.5) + #, width = 0.03
              
              # Add the centromeres
              geom_point(data = centromeres,
                         aes(x = chrnum,
                             y = start, color = Type, shape = Type),
                         size = 2, alpha = 1) +
              
              # ## Set the colors (the breaks are necessary to force the order!)
              scale_color_manual("Type", values = c("Centromere" = "indianred4", "NLR" = "cadetblue3", "NLR_HIC" = "#1F78B4", "Random" = "gray30")) +
              
              # Position the legend inside the plot
              labs(x = NULL, y = "Position (bp)") +
              guides(color = guide_legend(order = 1, title = "Features"), shape = guide_legend(order = 1)) +
              theme(legend.position = c(0.8, 0.58), # Adjust the position as needed
                    axis.line.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.background = element_rect(fill = alpha('white', 0.5))) )  # Add a background to the legend

ggsave(plot = NLRplot, 
       filename = NLRplotfile, 
       width = 8, height = 4.5)

### Some numbers
bed_data <- bed_data %>% mutate(size = end - start)
totalPodan2size <- chromosome_lengths$length %>% sum

bed_data %>% group_by(Region) %>% dplyr::summarise(total = sum(size)) %>% mutate(perc = (total*100)/totalPodan2size)
# 1 Introgressed 2056000  5.87
# 2 Missing      2927000  8.36

# How many of each group?
gene_data %>% filter(Domain %in% c("HET", "TPR", "ANK", "WD40", "NACHT", "NB-ARC")) %>% group_by(Domain) %>% summarise(n = n())

# How many NLRs?
classifiedNLRs <- gene_data[grep('NLR', gene_data$geneID),]
# classifiedNLRs %>% group_by(geneID) %>% dplyr::select(-PFAM) %>% distinct() %>% View

gene_data %>% group_by(geneID) %>% 
  filter(Domain %in% c("TPR", "ANK", "WD40", "NACHT", "NB-ARC", "HET")) %>% 
  dplyr::select(Domain, geneID) %>% 
  distinct() %>% 
  group_by(geneID) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) %>% # To exclude non-NLRs that only have one of the domains, but it's not ideal
  nrow()

# gene_data %>% filter(geneID == "NLR002")


### If not by domain, is the proportion of NLRs different in introgressed regions?
# There are 11655 protein-coding genes in Podan2-nice.3.02. From those 571
# are located in missing data regions that are not overlapping with
# introgression. So effectively we can see 11655 - 571
# = 11084 genes. I annotated 54 NLRs, but from those 11 are located in missing data
# regions (and not introgressed). So effectively we can only see 54-11 = 43
# NLRs. From those, 13 are in introgressed regions and 9 are in high (>=2) Tajima's
# D regions

# 11084 genes in total - 43 NLRs = 11041 From those, 756 fell into introgressed regions
# and 203 fell into high Tajima's D regions

total <- 11655
totalNLRs <- 54
totalIntro <- 756
totalMissing <- 571
introNLRs <- 13
missingNLRs <- 11

restNLRs <- totalNLRs-introNLRs-missingNLRs

introOthers <- totalIntro - introNLRs
missingOthers <- totalMissing - missingNLRs
restOthers <- total - totalNLRs - introOthers - missingOthers 

NLRsIntrogressed <- matrix(c(introNLRs, restNLRs, introOthers, restOthers),
                           nrow = 2,
                           dimnames = list(c("Introgressed", "Rest"), c("NLR", "Other")) )

cat("Are NLRs overrepresented in the introgressed set?")
fisher.test(NLRsIntrogressed, alternative = "two.sided") # p-value = 3.359e-06

totalTajima <- 1188
TajimaNLRs <- 9

restNLRs <- totalNLRs-TajimaNLRs-missingNLRs

TajimaOthers <- totalTajima - TajimaNLRs
restOthers <- total - totalNLRs - missingOthers - TajimaOthers

NLRsTajima <- matrix(c(TajimaNLRs, restNLRs, TajimaOthers, restOthers),
                     nrow = 2,
                     dimnames = list(c("Introgressed", "Rest"), c("NLR", "Other")) )
cat("Are NLRs overrepresented in the high Tajima's D set?")
fisher.test(NLRsTajima, alternative = "two.sided") # p-value = 0.04339
