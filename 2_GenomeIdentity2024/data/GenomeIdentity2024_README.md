# GenomeIdentity2024: A pipeline to compare identity along the genome against a reference

The aim of the pipeline is to find introgression (or ancestral polymorphism) tracks, but updated to 2024. By calling SNPs against the reference and then plotting the identity patterns, tracks of identity should (hopefully) pop-up.

    $ cd /proj/sllstore2017101/b2015200/VRwork/VRpipelines/06_GenomeIdentity2024

In all these pipelines, I didn't put the plotting part because I wasn't sure what would be the final thing, but also because it takes forever to run them... but they all have a folder called `scripts` with R scripts for plotting stuff.

## Building the environment


This pipeline relies on the VariantCalling2024 [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) environment. But first, I will update my base conda:

    $ conda update -n base conda

Normally installing a ton of things in a single environment can be problematic, but I'm trying to avoid installing the same things many times in the cluster. So here goes the main global environment for most of the tasks (rules) in the pipeline:

    $ mamba create -n VariantCalling2024 -c bioconda snakemake-minimal=8.10.4 bwa=0.7.17 samtools=1.19.2 picard=2.19.0 gatk4=4.5.0.0 repeatmasker=4.1.5 bedtools=2.31.1 vcftools=0.1.16 bcftools=1.10 htslib=1.19.1 bedops=2.4.41

Because there are so many programs, it takes a whileee. Notice that I couldn't use a more updated bcftools version because I kept having issues with 

    bcftools: error while loading shared libraries: libgsl.so.25: cannot open shared object file: No such file or directory

In order to use profiles in the slurm cluster Uppmax, I also installed:

    $ mamba install snakemake-executor-plugin-cluster-generic=1.0.7

    $ conda activate VariantCalling2024

For the R scripts I used an independent little environment:

    $ cat envs/popgenome.yaml
```yaml
channels:
  - bioconda
  - conda-forge
  - r
dependencies:
  - r-ggplot2=3.3.0
  - r-dplyr=0.8.5
  - r-popgenome=2.6.1
  - r-vcfr=1.10.0
  - r-reshape2=1.4.4
  - r-hmisc=4.2_0
```

## Configuration file

This pipeline depends on a given configuration file including the samples, the path to the data, and the reference. 


## Run pipeline in Uppmax

Get into the folder:

    $ cd /proj/sllstore2017101/b2015200/VRwork/VRpipelines/06_GenomeIdentity2024

To get an idea of how the pipeline looks like we can make a rulegraph:

    $ snakemake --profile profile --rulegraph | dot -Tpng > rulegraph.png

To check that the files for the pipeline are in order:

    $ snakemake --profile profile -pn

There are many ways of running the pipeline. In this case I'm using the profile file defined above:

    $ screen -R popgen
    $ mamba activate VariantCalling2024
    $ snakemake --profile profile &> snakemake.log &
    [1] 38448

Notice that sometimes GATK freaks out and crashes if more than one thread is used. So I ended up running the whole thing with one thread.
