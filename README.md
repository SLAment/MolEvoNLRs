# Molecular evolution of NLRs in the *Podospora anserina* species complex

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

Here you'll find the code associated with the manuscript:

> Ament-Velásquez & Saupe (2025) "NOD-like receptor genes evolve under diversity-enhancing mechanisms in a fungal species complex", biorxiv, [https://doi.org/10.1101/2025.09.29.679196](https://doi.org/10.1101/2025.09.29.679196).

-----------

Here you'll find all the code necessary to replicate the figures in the paper (except for the phylogenies and a few I put together in Inkscape). They are all [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines that depend on [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environments.

The standard with snakemake is that the pipeline is a file called `Snakefile`. However, I actually name them something else, like `mypipeline.smk` simply because it helps me to keep track of what pipeline is doing what other than based on it's path. Feel free to rename the files when you are working.

The directories in the repository are ordered to reflect the order of analyses in the paper. But more specifically, the figures were made in the following pipelines:

* Fig 3 - NLRvsRandomGenes.smk
* Fig 4 - NLRvsRandomGenes.smk
* Fig 5 - NLRvsRandomGenes.smk
* Fig 6 - 2SppIntrogression.smk
* Fig 7 - IntroRegions.smk
* Fig 8 - REPEATexplorerTPR.smk
* Fig 9 - smartHIC.smk
* Fig S1 - NLRvsRandomGenes.smk
* Fig S2 - NLRvsRandomGenes.smk
* Fig S3 - NLRvsRandomGenes.smk
* Fig S5 - GenomeIdentity2024.smk
* Fig S6 - 3SppIntrogression.smk
* Fig S7 - 2SppIntrogression.smk
* Fig S8 - 2SppIntrogression.smk
* Fig S9 - 2SppIntrogression.smk (partially, the other part was with IQ-TREE outside a pipeline) 
* Fig S10 - IntroRegions.smk
* Fig S11 - REPEATexplorerTPR.smk
* Fig S12 - REPEATexplorer.smk
* Fig S13 - smartHIC.smk
* Fig S14 - smartHIC.smk

----

Disclaimer: These scripts and files are provided "as is" and without any express or implied warranties, including, without limitation, the implied warranties of merchantability and fitness for a particular purpose.
