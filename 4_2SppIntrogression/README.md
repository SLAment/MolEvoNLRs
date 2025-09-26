# 2SppIntrogression: A pipeline to compare pair-wise identity of P. anserina/pauciseta/comata

The aim of the pipeline is to find introgression (or ancestral polymorphism) tracks. It uses the output of the `GenomeIdentity2024.smk` pipeline and re-runs some of the steps but only with the two species to reduce the amount of mising data. Unlike `3SppIntrogression.smk`, *P. pauciseta* and *P. comata* are treated separately. Putting them together has the advantage that it allows direct comparisons of genome proportion under introgression. But the cost is the ability to detect introgression of NLRs, which are presence-absence and might be in one species but not the other. Another difference is that here I will use OVERLAPPING windows to make them equivalent to the [NEE paper](https://doi.org/10.1038/s41559-022-01734-x).

## Building the environment

This pipeline relies on the VariantCalling2024 [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) environment I made for the `GenomeIdentity2024.smk` pipeline, which was made like so:

    $ mamba create -n VariantCalling2024 -c bioconda snakemake-minimal=8.10.4 bwa=0.7.17 samtools=1.19.2 picard=2.19.0 gatk4=4.5.0.0 repeatmasker=4.1.5 bedtools=2.31.1 vcftools=0.1.16 bcftools=1.10 htslib=1.19.1 bedops=2.4.41

Because there are so many programs, it takes a whileee. Notice that I couldn't use a more updated bcftools version because I kept having issues with 

    bcftools: error while loading shared libraries: libgsl.so.25: cannot open shared object file: No such file or directory

In order to use profiles in the slurm cluster Uppmax, I also installed:

    $ mamba install snakemake-executor-plugin-cluster-generic=1.0.7

    $ conda activate VariantCalling2024

For the R scripts I used an independent little environment, but it relies on the ones used by `3SppIntrogression.smk` so it expects this path (but you can change it in the Snakemake file):

    $ cat ../3_3SppIntrogression/envs/popgenome.yaml
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

    % cat ../3_3SppIntrogression/envs/simpleplot.yaml
```yaml
channels:
  - bioconda
  - conda-forge
  - r
dependencies:
  - r-dplyr=0.8.3
  - r-tidyverse=1.2.1 # It comes with ggplot 3.0.0
  - r-cowplot=1.0.0
  - r-patchwork
```

## The profile

For this pipeline I use a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), which has the information necessary to run the pipeline in a SLURM server. It depends on a file called `config.v8+.yaml` in the `profile` folder.

The config.yaml file contains:

```yaml
snakefile: 2SppIntrogression.smk
executor: cluster-generic

cluster-generic-submit-cmd:
  mkdir -p logs/{rule} &&
  sbatch
    --account={resources.account}
    --partition={resources.partition}
    --cpus-per-task={resources.threads}
    --mem={resources.mem_mb}
    --job-name={rule}
    --error=logs/{rule}/{rule}-{wildcards}-%j.err
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
    --time={resources.time}
    --parsable
default-resources:
  - account="XXXXXXXXX"
  - partition="core"
  - time="1:00:00"
  - threads=1
  - mem_mb=2GB

restart-times: 0
max-jobs-per-second: 10
max-status-checks-per-second: 1
jobs: 100
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
use-conda: True
cluster-cancel: scancel # To automatically cancel all running jobs when you cancel the main Snakemake process 
cluster-cancel-nargs: 50
```

Where XXXXXXXXX is your cluster account. Replace that! Be also mindful of the name of the partitions used in your cluster and the memory given per thread.

You can run the pipeline without the profile (using `--snakefile 2SppIntrogression.smk`), but keep in mind that I designed the rules to call for resources from it. But this can be easily modified if you are a bit familiar with Snakemake.

I also used a configuration file, but that is called internally in the pipeline. It's in the `config` folder.

## Configuration file

This pipeline depends on a given configuration file pointing to the results from the `GenomeIdentity2024.smk` pipeline.

```yaml
# File with a list of samples to analyze
SampleIDs_file: "../2_GenomeIdentity2024/data/IndividualsIDsSpp.txt"

## Species being compared to P. anserina
comparedspecies: ["pauciseta", "comata"]

# VCF file before filtering for missing data
VCFfile: "/proj/sllstore2017101/b2015200/VRwork/VRpipelines/06_GenomeIdentity2024/vcfs/PodoPop-vs-Podan2-snps-NoTEs-gatkPASS.vcf.gz"
PathToCovBED: "/proj/sllstore2017101/b2015200/VRwork/VRpipelines/06_GenomeIdentity2024/genomecov"
TEsites: "/proj/sllstore2017101/b2015200/VRwork/VRpipelines/06_GenomeIdentity2024/data/cov/TEsites.bed"

# Output name
OUTPUTbasename: "PodoPop-vs-Podan2-snps-NoTEs-gatkPASS"

# Variables
WINSIZE: 10000 # Size of overlapping windows to calculate dxy, like in the het-v NEE paper
JUMP: 1000
MINWINSIZE: 5000
PI: 0.00049

# Scripts
PopGenome_WinMissingness: "../2_GenomeIdentity2024/scripts/PopGenome_WinMissingness.R"
PopGenome_PairwiseDivergence: "../2_GenomeIdentity2024/scripts/PopGenome_PairwiseDivergence.R"
```

## Run pipeline in a SLURM cluster

Get into the folder where the pipeline is. 

To get an idea of how the pipeline looks like we can make a rulegraph:

    $ snakemake --profile profile --rulegraph | dot -Tpng > rulegraph.png

To check that the files for the pipeline are in order:

    $ snakemake --profile profile -pn

There are many ways of running the pipeline. In this case I'm using the profile file defined above:

    $ screen -R popgen2spp
    $ mamba activate VariantCalling2024
    $ snakemake --profile profile &> snakemake.log &
    [1] 1532

Naturally it can also be run locally, but the rules were designed to call for resources. Use with care or modify the rules.

    $ snakemake --snakefile 2SppIntrogression.smk --rulegraph | dot -Tpng > rulegraph.png
    % snakemake --snakefile 2SppIntrogression.smk --use-conda -j8

## Results

The folder `figures` contains the main Figure 6, and the supplementary figures S7, S8, and part of S9. The `results` folder contains files used by the pipeline `IntroRegions.smk` and the files used by the script `PaNLRchrPopStats2Spp.R` to make the plots.


