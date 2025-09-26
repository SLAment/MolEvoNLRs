#!/usr/bin/env Rscript

### PaNLRchrPopStats2Spp: Exploring introgression in the *P. anserina* species complex
#############################################################################
# https://stackoverflow.com/questions/14563989/force-r-to-stop-plotting-abbreviated-axis-labels-scientific-notation-e-g-1e
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2024/10/17
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(cowplot) # For plot_grid
library(RColorBrewer) # For brewer.pal
library(scales) 
# ============================

# ============================
# Input files
# ============================
# Input
## Pauciseta
DxyWinsChr1_pau <- read.csv(snakemake@input$chr1sp1)
DxyWinsChr2_pau <- read.csv(snakemake@input$chr2sp1)
DxyWinsChr3_pau <- read.csv(snakemake@input$chr3sp1)
DxyWinsChr4_pau <- read.csv(snakemake@input$chr4sp1)
DxyWinsChr5_pau <- read.csv(snakemake@input$chr5sp1)
DxyWinsChr6_pau <- read.csv(snakemake@input$chr6sp1)
DxyWinsChr7_pau <- read.csv(snakemake@input$chr7sp1)

## Comata
DxyWinsChr1_com <- read.csv(snakemake@input$chr1sp2)
DxyWinsChr2_com <- read.csv(snakemake@input$chr2sp2)
DxyWinsChr3_com <- read.csv(snakemake@input$chr3sp2)
DxyWinsChr4_com <- read.csv(snakemake@input$chr4sp2)
DxyWinsChr5_com <- read.csv(snakemake@input$chr5sp2)
DxyWinsChr6_com <- read.csv(snakemake@input$chr6sp2)
DxyWinsChr7_com <- read.csv(snakemake@input$chr7sp2)

# Population genetic statistics
PopDataChr1 <- read.csv(snakemake@input$pop1, header = TRUE, sep=",")
PopDataChr2 <- read.csv(snakemake@input$pop2, header = TRUE, sep=",")
PopDataChr3 <- read.csv(snakemake@input$pop3, header = TRUE, sep=",")
PopDataChr4 <- read.csv(snakemake@input$pop4, header = TRUE, sep=",")
PopDataChr5 <- read.csv(snakemake@input$pop5, header = TRUE, sep=",")
PopDataChr6 <- read.csv(snakemake@input$pop6, header = TRUE, sep=",")
PopDataChr7 <- read.csv(snakemake@input$pop7, header = TRUE, sep=",")

genesall <- read.table(snakemake@input$nlrmeta, header = TRUE)
PopData <- read.table(snakemake@input$indvs, header=FALSE)

# ============================
# Reading the data
# ============================

# Put together
DxyWins_pau <- rbind(DxyWinsChr1_pau, DxyWinsChr2_pau, DxyWinsChr3_pau, DxyWinsChr4_pau, DxyWinsChr5_pau, DxyWinsChr6_pau, DxyWinsChr7_pau)
DxyWins_com <- rbind(DxyWinsChr1_com, DxyWinsChr2_com, DxyWinsChr3_com, DxyWinsChr4_com, DxyWinsChr5_com, DxyWinsChr6_com, DxyWinsChr7_com)
DxyWins_pau$Species <- "pauciseta"
DxyWins_com$Species <- "comata"
DxyWins <- rbind(DxyWins_pau, DxyWins_com)

# Process the strain data frame
names(PopData) <- c("Strain", "Species")
# Metadata
PopDataStrain2 <- PopData %>% filter(Species != "anserina")
names(PopDataStrain2) <- c("Strain2", "Species")

# The coordinates of the genes, used by the plotting functions
genes <- genesall %>% filter(infig == 1)

# Make a dataframe with the lengths of the Podan2 chromosomes
chrlines <- data.frame(pos =c(0), endchr = c(8813524, 5165621, 4137471, 3808395, 4734292, 4264132, 4087160) , chr = c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7"), Strain = character(1))

# Upper limit for plotting divergence
maxPi <- 0.015 

# ============================
# Preparing the data
# ============================
# Add species info
pwdivergence <- merge(DxyWins, PopDataStrain2) # This takes a while
pwdivergence$position <- as.numeric(pwdivergence$position)
pwdivergence$dxy <- as.numeric(pwdivergence$dxy)

# https://rgbacolorpicker.com/rgba-to-hex
farg <- c("pseudoanserina" = "#cc79a7", 
          "pseudopauciseta" = "#56b4e9",
          "pseudocomata" = "#e69f00",
          "bellae-mahoneyi" = "#f0e442",
          "comata" = "#d45500",
          "pauciseta" = "#009e73")

# I need a lot of colors to cover all the strains
nb.cols <- length(unique(pwdivergence$Strain1))
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)

minlen = 5000 # As in Ament-Velásquez et al. (2022) NEE. But it is too stringent to see many het and NLR genes.

# ============================
# Functions
# ============================
# To print the scale with more digits
scaleFUN <- function(x) sprintf("%.3f", x)

# Amazing function to italize some labels in legend by user mrFlick
# https://stackoverflow.com/questions/44641870/ggplot2-formatting-legend-categories
# But I couldn't get it to interact with the shapes :/
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Centromere", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

# dxydata = pwdivergence
# pitajdata = PopDataChr3
# thischr = "chromosome_3"
# minlen = 5000
# mincoord = 3500000 
# maxcoord = 3537000
# genes = genes
# annotations = annotations_PaFnt1
# pilegend = TRUE
# maxDxy = 0.044

## Big function to plot the population genomics statistics and the divergence along the genome
popgenANDdiverg <- function(dxydata, pitajdata, thischr = "chromosome_4", minlen = 5000, mincoord = 1, maxcoord = NA, genes = genes, annotations = NA, pilegend = TRUE, pilegendpos = 0.5, pilegendposy = 0.7, rmAllyTitles = FALSE, maxDxy = 0.044, maxpi = maxPi, legenddirection = "vertical"){
  # preparing variables
  if (is.na(maxcoord)) {
    maxcoord <- max(dxydata$position)
  }
  
  chosenstrains <- c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m")
  
  if (is.character(annotations)) {
    annotations <- data.frame(
      Strain2 = chosenstrains,
      label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
      x = mincoord + 20000,  # Position for the text on x-axis
      y = 0.038         # Position for the text on y-axis
    )
  } 
  
  if (pilegend) {
    pileg <- guide_legend(order = 2)
  } else {
    pileg <- "none"
  }
  
  dxydatachr <- dxydata %>% filter(chr == thischr)
  dxydatachr$dxy[which(dxydatachr$winlen <= minlen )] <- NA
  
  geneschr <- genes[order(genes$Locus),] %>% filter(chr == thischr) %>% mutate(value = 0, dxy = 0) %>% filter(Start >= mincoord, End <= maxcoord)
  
  # Chromosome number for plotting
  chrnumb <- substr(thischr, nchar(thischr), nchar(thischr))
  
  # -------
  ## Plots
  # -------
  ### Dxy
  Dxyplot <- ggplot(dxydatachr %>% filter(Strain2 %in% chosenstrains),
         aes(x = position, y = dxy, colour = Strain1)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_grid(Strain2 ~ .) +
    scale_colour_manual(values = mycolors, guide = "none") +
    scale_x_continuous(expand = c(0.02, 0.02), labels = function(x) format(x, scientific = FALSE)) + # Remove the white space on the sides
    scale_y_continuous(labels = scales::label_number(accuracy = 0.001), limits = c(0, maxDxy)) +  # 3 significant digits
    theme_bw() +
    theme(axis.title= element_text(size = rel(1.2)),
          # axis.text.x = element_text(size = rel(0.9)),
          # axis.text.y = element_text(size = rel(1)),
          axis.text = element_text(size = rel(1)),
          strip.background = element_rect(fill="white"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          legend.key=element_blank(), # remove background of the key in the legend
          legend.title=element_blank(),
          legend.text=element_text(size=rel(0.9)),
          legend.spacing.y = unit(0.005, 'cm'), # Put the items of the legend closer to each other
          strip.text = element_text(size = rel(0.8)), # Adjust facet title size
          legend.background=element_blank()) +
    # ylab(expression(paste(D[XY]~"between ", italic("P. anserina"), " and other species"))) +
    ylab(expression(paste("Divergence between ", italic("P. anserina"), " and other species"))) +
    xlab(paste("Position (bp) in Chromosome", chrnumb)) +
    geom_point(data = geneschr, aes(x = position, shape = Locus, fill = Locus), colour = "black", size = 3) +
    scale_shape_manual(values = geneschr$shapes) +
    coord_cartesian(xlim = c(mincoord, maxcoord)) +
    geom_text(data = annotations, aes(x = x, y = y, label = label), 
              color = "black", size = 4, hjust = 0, fontface = "italic")
  
  ### Pi and Watersons theta
  PopDataFiltered <- pitajdata %>% filter(Pop == "All") %>% filter(position >= mincoord, position <= maxcoord) # "All" is all of Wageningen, as opposed to the groups defined by het-v
  PopDataFiltered$value[which(PopDataFiltered$winlens <= minlen )] <- NA 
  
  piplot <- ggplot(PopDataFiltered %>% filter(variable != "Tajima"), aes(x = position, y = value, colour = variable)) +
    theme_bw() + # Remove background color
    labs(y = expression(pi*" or "*theta["W"])) +
    scale_x_continuous(expand = c(0.02, 0.02), labels = function(x) format(x, scientific = FALSE)) + # Remove the white space on the sides
    theme(plot.title = element_text(size = rel(1.3), face="bold", family = "Helvetica", hjust = 0.5), 
          # plot.margin = margin(t = 0, r = 10, b = 0, l = 3, unit = "pt"),
          axis.title.y= element_text(size = rel(1.2)),
          axis.title.x=element_blank(), 
          axis.text = element_text(size = rel(1)),
          panel.grid.major = element_blank(), # remove grid lines
          panel.grid.minor = element_blank(), # remove grid lines
          legend.key=element_blank(),
          legend.title=element_blank(),
          legend.key.height=unit(0.8,"line"), # Make the items *within* a legend closer to each other
          legend.text=element_text(size=rel(0.9)),
          legend.position=c(pilegendpos, pilegendposy),
          legend.spacing.y = unit(0.005, 'cm'), # Put the items of the legend closer to each other
          legend.spacing.x = unit(0.005, 'cm'), # Put the labels of the legend closer to the symbols
          legend.box = "horizontal",
          legend.direction = legenddirection, # OJO
          legend.background=element_blank()) +
    scale_color_manual(values= c("firebrick1", "blue4"), labels=c(expression(pi["  "]), expression(theta["W"]))) + # The space in pi is to center it better...
    ylim(0, maxpi) +
    geom_line(alpha = 0.7) + 
    geom_point(data = geneschr, aes(x = position, shape = Locus, fill = Locus), colour = "black", size = 3) + 
    scale_shape_manual(values = geneschr$shapes) + # Give it specific shapes (but re-order the dataframe so it maches)
    guides(colour = pileg,  # Line legend
           fill = guide_legend(order = 1),    # Point legend moved next to the line legend
           shape = guide_legend(order = 1))   # Point legend moved next to the line legend
  
  tajimaplot <- ggplot(PopDataFiltered %>% filter(variable == "Tajima"), aes(x = position, y = value, colour = variable, fill = variable)) +
    geom_line(alpha = 0.4) + 
    theme_bw() + # Remove background color
    labs(y = "Tajima's D") + 
    scale_x_continuous(expand = c(0.02, 0.02), labels = function(x) format(x, scientific = FALSE)) + # Remove the white space on the sides, and force full numbers instead of scientific notation
    theme(plot.title = element_text(size = rel(1.3), face="bold", family = "Helvetica", hjust = 0.5), 
          # plot.margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"),
          axis.title.y= element_text(size = rel(1.2),),
          axis.title.x=element_blank(),
          axis.text = element_text(size = rel(1)),
          panel.grid.major = element_blank(), # remove grid lines
          panel.grid.minor = element_blank(), # remove grid lines
          legend.position="none") +
    scale_color_manual(values="darkorange3", labels= "Tajima's D", guide= 'none') +
    geom_area(alpha = 0.5, fill = "darkorange3") +
    scale_y_continuous(labels=scaleFUN, limits = c(-3, 4)) +
    geom_point(data = geneschr, aes(x = position, shape = Locus, fill = Locus), colour = "black", size = 3,
               show.legend = "none") + guides(fill="none") + # The guides remove part of the legend
    scale_shape_manual(values= geneschr$shapes) # Give it specific shapes (but re-order the dataframe so it maches)
  
  if (rmAllyTitles) {
    Dxyplot <- Dxyplot + theme(axis.title.y=element_blank())
    piplot <- piplot + theme(axis.title.y=element_blank())
    tajimaplot <- tajimaplot + theme(axis.title.y=element_blank())
  }
  
  together <- plot_grid(piplot, tajimaplot, Dxyplot, ncol = 1, align = "v", rel_heights = c(1,1,3), axis = "lr") # the align part ensures the x axis matches
  return(together)
}

# ============================
# Main figure: het-e and hnwd3
# ============================

annotations_chr4L <- data.frame(
  Strain2 = c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m"),
  label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
  x = 10000,  # Position for the text on x-axis
  y = 0.037         # Position for the text on y-axis
)

cool_chr4L <- popgenANDdiverg(pwdivergence, PopDataChr4, 
                thischr = "chromosome_4", 
                mincoord = 1, 
                maxcoord = 800000, 
                genes = genes, 
                annotations = annotations_chr4L)

annotations_hnwd3 <- data.frame(
  Strain2 = c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m"),
  label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
  x = 652000,  # Position for the text on x-axis
  y = 0.037         # Position for the text on y-axis
)

cool_hnwd3 <- popgenANDdiverg(pwdivergence, PopDataChr7, 
                thischr = "chromosome_7", 
                mincoord = 650000, 
                maxcoord = 800000, 
                genes = genes, 
                annotations = annotations_hnwd3,
                pilegend = FALSE,
                pilegendpos = 0.5,
                rmAllyTitles = TRUE)

mainfig <- plot_grid(cool_chr4L, cool_hnwd3, labels = c('a', 'b'), label_x = c(0, -0.02), nrow = 1, rel_widths = c(2,1), label_size = 15)

ggsave(plot = mainfig,
       filename = snakemake@output$main,
       width = 10, height = 8)

# ============================
# Supp figure 7: nwd5 and nwd6
# ============================

annotations_nwd5 <- data.frame(
  Strain2 = c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m"),
  label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
  x = 602000,  # Position for the text on x-axis
  y = 0.037         # Position for the text on y-axis
)

cool_nwd5 <- popgenANDdiverg(pwdivergence, PopDataChr5, 
                              thischr = "chromosome_5", 
                              mincoord = 600000, 
                              maxcoord = 700000, 
                              genes = genes, 
                              pilegendpos = 0.5,
                              pilegendposy = 0.8,
                              annotations = annotations_nwd5)

annotations_nwd6 <- data.frame(
  Strain2 = c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m"),
  label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
  x = 652000,  # Position for the text on x-axis
  y = 0.037         # Position for the text on y-axis
)

cool_nwd6 <- popgenANDdiverg(pwdivergence, PopDataChr5, 
                             thischr = "chromosome_5", 
                             mincoord = 3000000, 
                             maxcoord = 3500000, 
                             genes = genes, 
                             annotations = annotations_nwd6,
                             pilegendpos = 0.5,
                             pilegendposy = 0.75,
                             pilegend = FALSE,
                             legenddirection = "horizontal",
                             rmAllyTitles = TRUE)

supfig7 <- plot_grid(cool_nwd5, cool_nwd6, labels = c('a', 'b'), label_x = c(0, -0.02), nrow = 1, rel_widths = c(1,2), label_size = 15)

ggsave(plot = supfig7,
       filename = snakemake@output$figs7,
       width = 10, height = 8)

# ============================
# Supp figure 8: Good cases
# ============================
annotations_chr3R1 <- data.frame(
  Strain2 = c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m"),
  label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
  x = 2820000,  # Position for the text on x-axis
  y = 0.037         # Position for the text on y-axis
)

# het-b
hetB <- popgenANDdiverg(pwdivergence, PopDataChr3, 
                        thischr = "chromosome_3", 
                        mincoord = 920000, 
                        maxcoord = 1000000, 
                        genes = genes, 
                        pilegendpos = 0.5,
                        pilegendposy = 0.6,
                        minlen = 3000,
                        annotations = annotations_chr3R1)
# This one is the best to show-case a gene with trans-species polymorphism but 
# with alleles that are still different within class (so no introgression)

cool_PNDUDP <- popgenANDdiverg(pwdivergence, PopDataChr3, 
                               thischr = "chromosome_3", 
                               mincoord = 2800000, 
                               maxcoord = 3100000, 
                               genes = genes, 
                               pilegendpos = 0.5,
                               pilegendposy = 0.6,
                               annotations = annotations_chr3R1)

supfig8 <- plot_grid(hetB, cool_PNDUDP, labels = c('a', 'b'), label_x = c(0, -0.02), nrow = 1, rel_widths = c(1,1), label_size = 15)

ggsave(plot = supfig8,
       filename = snakemake@output$figs8,
       width = 8, height = 8)

# ============================
# Supp figure 8: het-z
# ============================

annotations_hetz <- data.frame(
  Strain2 = c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m"),
  label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
  x = 3930000,  # Position for the text on x-axis
  y = 0.037         # Position for the text on y-axis
)


(hetz <- popgenANDdiverg(pwdivergence, PopDataChr1,
                         thischr = "chromosome_1", 
                         mincoord = 3925000, 
                         maxcoord = 4000000, 
                         genes = genes, 
                         maxDxy = 0.05, 
                         maxpi = 0.027,
                         pilegendpos = 0.35,
                         annotations = annotations_hetz))

ggsave(plot = hetz,
       filename = snakemake@output$hetz,
       width = 4, height = 8)

# ============================
# Pa_2_8180
# ============================

annotations_Pa_2_8180 <- data.frame(
  Strain2 = c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m"),
  label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
  x = 602000,  # Position for the text on x-axis
  y = 0.037         # Position for the text on y-axis
)

cool_Pa_2_8180 <- popgenANDdiverg(pwdivergence, PopDataChr2,
                             thischr = "chromosome_2", 
                             mincoord = 3080000, 
                             maxcoord = 3300000, 
                             genes = genes,
                             annotations = annotations_Pa_2_8180)

ggsave(plot = cool_Pa_2_8180,
       filename = snakemake@output$pnpudp,
       width = 6, height = 8)

# ============================
# PaFnt1 and HELLF
# ============================

annotations_chr3R2 <- data.frame(
  Strain2 = c("CBS237.71m", "CBS333.63p", "PcTdp", "PcWa139m"),
  label = c("P. pauciseta", "P. pauciseta", "P. comata", "P. comata"),
  x = 3500000,  # Position for the text on x-axis
  y = 0.037         # Position for the text on y-axis
)

cool_PaFnt1 <- popgenANDdiverg(pwdivergence, PopDataChr3, 
                                thischr = "chromosome_3", 
                                mincoord = 3450000, 
                                maxcoord = 3600000, 
                                genes = genes, 
                                pilegendpos = 0.5,
                                pilegendposy = 0.6,
                                annotations = annotations_chr3R2,
                                pilegend = FALSE, 
                                rmAllyTitles = TRUE)

# This one is tricky because there is a lot of missing data (e.g. absent in the 
# Argentinian paucis). But from the sequence alignment it looks like it is not 
# introgressed. Israel is very different.  
# However there are multiple, somewhat different alleles within P. anserina. 

ggsave(plot = cool_PaFnt1,
       filename = snakemake@output$Fnt1,
       width = 4, height = 8)

