# GenomeIdentity2024: A pipeline to compare identity along the genome against a reference

The aim of the pipeline is to find introgression (or ancestral polymorphism) tracks, but updated to 2024. By calling SNPs against the reference and then plotting the identity patterns, tracks of identity should (hopefully) pop-up.

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

## The profile

For this pipeline I use a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), which has the information necessary to run the pipeline in a SLURM server. It depends on a file called `config.v8+.yaml` in the `profile` folder.

The config.yaml file contains:

```yaml
snakefile: GenomeIdentity2024.smk
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

You can run the pipeline without the profile (using `--snakefile GenomeIdentity2024.smk`), but keep in mind that I designed the rules to call for resources from it. But this can be easily modified if you are a bit familiar with Snakemake.

I also used a configuration file, but that is called internally in the pipeline. It's in the `config` folder.

## Configuration file

All the input files should be specified in a configuration file in yaml format, called `config/config.yaml` that looks like so:

```yaml
# File with a list of samples to analyze
SampleIDs_file: "data/IndividualsIDsSpp.txt"

# Path to Illumina samples from Ament-Velásquez et al. (2022) NEE
IlluminaPath: "/proj/sllstore2017101/b2015200/SnakePipelines/0_CleanIllumina/postQC-data"

# The reference genome
refgenome: "path/to/Podan2_AssemblyScaffoldsmt.fa"
refname: "Podan2"

# Name of the mitochondrial scaffold to exclude it from the final vcf file
mitoctg: "PaMt_NC_001329.3"

# TE library
TElib: "path/to/PodoTE-1.00.lib"

# Variables
WINSIZE: 2000 # Size of non-overlapping windows to calculate dxy

# Scripts
vcfR_plotter: "scripts/PerSample_vcfR_plotter.R"
PopGenome_WinMissingness: "scripts/PopGenome_WinMissingness.R"
PopGenome_PairwiseDivergence: "scripts/PopGenome_PairwiseDivergence.R"
PaIntrogressionTracks: "scripts/PaIntrogressionTracks2024.R"
PaIntrogressionProportions: "scripts/PaIntrogressionProportions.R"
```

The Illumina data can be retrieved from NCBI, using the accession numbers available in the Supplementary Table 1 of [Ament-Velásquez et al. (2022) NEE](https://doi.org/10.1038/s41559-022-01734-x). 

The reference assembly is called Podan2. It corresponds to the second iteration of the reference genome of the strain S+ ([Espagne et al. 2008](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-5-r77)), produced by JGI's Mycocosm. The sequence seems to be identical to the latest version Podan3, but the newest annotation in JGI discarded the original gene codes (e.g., Pa_Z_XXXX, where Z is the chromosome and XXXX is the gene number) that the community has been using since 2008. It's available [here](https://github.com/johannessonlab/HetVPaper/blob/master/SNPpop/data/Podan2/Podan2_AssemblyScaffoldsmt.fa) and in the legacy websites of JGI. There is a [newer version of the S+ reference genome](https://doi.org/10.1186/s12864-022-09085-4), but using this version allows full compatibility with out previous work.

The TE library `PodoTE-1.00.lib` is available [here](https://github.com/johannessonlab/SpokBlockPaper/blob/master/Annotation/data/PodoTE-1.00.lib).


## Run pipeline in Uppmax

Get into the folder where the pipeline is. 

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

To run locally do:

    $ snakemake --snakefile GenomeIdentity2024.smk --rulegraph | dot -Tpng > rulegraph.png

    $ snakemake --snakefile GenomeIdentity2024.smk --use-conda -j8

But notice the pipeline is designed for a cluster, adjust as necessary.

## Results

The pipeline will produce vcf files with various degrees of filtering. The vcf `vcfs/PodoPop-vs-Podan2-snps-NoTEs-gatkPASS-miss1.vcf.gz` contains the variants after removing missing data (i.e. no missing data in some but not other individuals is allowed). I create a non-standard VCF file with this one as input. Basically I re-entered all sites that were marked as missing data as `.:0,0:0:.:0,0`. Hence the missing data is explicitly coded in the VCF file and all other sites that are not in the VCF are invariant compared to the reference. I make this file so that I can calculate correctly the nucleotide diversity (pi), or in this case Dxy, while taking into account missing data and without having to handle huge files with *all sites* explicitly. The resulting file is `vcfs/PodoPop-vs-Podan2-snps-NoTEs-gatkPASS-miss1-withNA.vcf.gz`, which is later split by chromosome so that PopGenome can handle it in R.

This file is later used by subsequent pipelines. Here it is also used to calculate divergence between all the strains from other *Podospora* species vs *P. anserina* but only in chromosome 4, since it was computationally demanding.


