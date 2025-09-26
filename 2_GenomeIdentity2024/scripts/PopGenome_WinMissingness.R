#!/usr/bin/env Rscript

### PopGenome_WinMissingness
#############################################################################
# After modifying a vcf file with `badsites2vcf.py`, the sites with missing data 
# (as opposed to monomorphic sites) should be explicit. Hence, we can calculate
# the actual number of sites in a window of a given size.
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2024/04/10
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(PopGenome)
library(dplyr, warn.conflicts = FALSE)
# ============================

# ============================
# Reading the data
# ============================

print("*** Reading data ***")
# Read entire vcf file and gff file
vcffolder <- dirname(snakemake@input$vcf) # the location of the vcf

# Bed with sites to be excluded
badsites <- read.table(snakemake@input$bed, header = FALSE, sep="\t")

# Windows parameters
width <- snakemake@params$width
jump <- snakemake@params$jump
chr <- snakemake@params$chr

# Output
outputname <- snakemake@output$df
# -----

# ============================
# Define functions
# ============================

#### Function to calculate the amount of missing data per window
## The idea is to reduce the size of each window with the overlapping regions in badsites
winmissing <- function(GENOME.class.slide, badsites){ 
  # Get coordinates of the windows for this GENOME.class.slide object
  coords <- c()
  for (win in GENOME.class.slide@region.names) {
    start <- strsplit(win," ")[[1]][c(1)] %>% as.numeric()
    end <- strsplit(win," ")[[1]][c(3)] %>% as.numeric()
    coords <- rbind(coords, cbind(start, end))
  }
  coords <- data.frame(coords)
  
  winlens <- c()
  for (win in seq(1, nrow(coords))){
    start1 <- coords[win,1] %>% as.numeric()
    end1 <- coords[win,2] %>% as.numeric()
    
    fulllenwin <- end1 - start1 + 1 # Plus one because the range includes the ends
    for (bad in seq(1, nrow(badsites))) {
      start2 <- badsites[bad,2] %>% as.numeric()
      end2 <- badsites[bad,3] %>% as.numeric()
      
      # Is there overlap between ranges?
      if ((start1 <= end2) & (end1 >= start2)) { # Yes
        # How much?
        if (start2 <= start1 & end2 <= end1) { # The bad range starts before the window and ends within
          fulllenwin <- fulllenwin - (end2 - start1 + 1)
        } else if (start2 < start1 & end2 > end1) { # The bad range is larger and includes the window
          fulllenwin <- 0
          break
        } else if (start2 > start1 & end2 < end1) { # The bad range is fully contained in the window
          fulllenwin <- fulllenwin - (end2 - start2 + 1)
        } else if (start2 > start1 & end2 > end1) { # The bad range starts within the window but ends after
          fulllenwin <- fulllenwin - (end1 - start2 + 1)
        }
      } else if (start2 > end1){ # Are the bad ranges after the window? then there is no point in continue
        break
      }
    }
    winlens <- c(winlens, fulllenwin)
  }
  return(winlens)
}

# ============================

GENOME.class <- readData(vcffolder, format="VCF", include.unknown = TRUE)

cat("*** Make windows ***\n")
GENOME.class.slide <- sliding.window.transform(GENOME.class, width = width, jump = jump, type=2) # 1 scan only biallelic positions (SNPs), 2 scan the genome

# # How many windows?
# GENOME.class.slide@region.names %>% length() # 379

cat("*** Calculate missing data in each window ***\n")
winlens <- winmissing(GENOME.class.slide, badsites)

# The slot GENOME.class@region.names will store the genomic regions of each window
# as a character string. To convert those strings into a genomic numeric
# position we can apply the following
genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
  split <- strsplit(x," ")[[1]][c(1,3)]
  val <- mean(as.numeric(split))
  return(val)
})

# Put it in a dataframe
missingness <- data.frame(chr = chr, pos = genome.pos, winlen = winlens, row.names = NULL)

# Write it into a file
write.csv(missingness, file = outputname, row.names = FALSE)
