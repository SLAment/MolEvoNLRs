# IntroRegions: Characterizing introgression regions between *P. anserina* and *P. comata* or *P. pauciseta*

The goal of this pipeline is to determine if the NLRs are enriched in the introgressed regions between *P. anserina* and *P. comata* or *P. pauciseta*. First the introgressed regions as well as regions with missing data per species comparisons are obtained from the pipeline `2SppIntrogression.smk`. I also use similar windows base on the Tajima's D calculation in Wageningen collection from Ament-Velásquez et al. (2022) NEE, in bed format. Then, using the manually curated annotation `Podan2.nice-3.02.gff3`, I obtain what genes overlap with said regions. I use InterProScan to re-annotate the proteome and I use the result to make an enrichment analysis of the introgressed regions.

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline.

# Environment

I ran this pipeline in a SLURM cluster. I created most of the environment with the [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) implementation of [conda](https://docs.conda.io/en/latest/).

To create the environment I mostly use Conda, although the installation of InterProScan is independent of that environment.

	% mamba create -n annotation -c bioconda snakemake-minimal=8.30.0 snakemake-executor-plugin-cluster-generic=1.0.9 biopython=1.85 gffutils=0.13 scipy=1.14.0 pandas=2.2.2 statsmodels=0.14.2 bedtools=2.31.1

And I used the InterProScan already installed in the cluster:

	$ module load interproscan/5.75-106.0

Unfortunately, the strict channel priority of conda is too strict to install the R packages used by the plotting rules, which depend on `envs/plot_scales.yaml`. Hence, I have to do:

	$ conda config --set channel_priority true

## The profile

For this pipeline I use a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), which has the information necessary to run the pipeline in a SLURM server. It depends on a file called `config.v8+.yaml` in the `profile` folder.

The config.yaml file contains:

```yaml
snakefile: IntroRegions.smk
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

You can run the pipeline without the profile (using `--snakefile IntroRegions.smk`), but keep in mind that I designed the rules to call for resources from it. But this can be easily modified if you are a bit familiar with Snakemake.

I also used a configuration file, but that is called internally in the pipeline. It's in the `config` folder.

## Configuration file

All the input files should be specified in a configuration file in yaml format, called `config/config.yaml` that looks like so:

```yaml
## Input files
# BED files from 2SppIntrogression.smk
INTROBED_both: "data/Dxy_10kb_introgression-all.bed"
INTROBED_pau: "data/Dxy_pauciseta_10kb_introgression.bed"
INTROBED_com: "data/Dxy_comata_10kb_introgression.bed"
INTROBED_miss: "data/Dxy_10kb_missingness-all.bed"
INTROBED_miss_pau: "data/Dxy_pauciseta_10kb_missingness.bed" # Missing data when calling variants in P. anserina and P. pauciseta
INTROBED_miss_com: "data/Dxy_comata_10kb_missingness.bed" # Missing data when calling variants in P. anserina and P. comata
INTROBED_taj: "data/High_Tajima_Wa.bed" # From the NEE paper, 10 kb windows with 1 kb steps (from which at least 5kb have data) and Tajima's D >= 2

FULLTajima: "data/WageningenTajimaMin5Kb.bed"

ANNOTATION: "path/to/Podan2.nice-3.02.gff3"
REFGENOME: "path/to/Podan2_AssemblyScaffoldsmt.fa"

PFAMref: "../HICproteins/data/PFAM_reference.txt"

GenesNLRsCoords: "data/GenesNLRsCoords.txt"

# Scripts
EnrichOddRatios: "scripts/EnrichOddRatios_alla.R"
EnrichOddRatiosPval: "scripts/EnrichOddRatiosPval.R"
NLRchrDist: "scripts/NLRchrDist.R"
```

The reference assembly is called Podan2. It corresponds to the second iteration of the reference genome of the strain S+ ([Espagne et al. 2008](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-5-r77)), produced by JGI's Mycocosm. The sequence seems to be identical to the latest version Podan3, but the newest annotation in JGI discarded the original gene codes (e.g., Pa_Z_XXXX, where Z is the chromosome and XXXX is the gene number) that the community has been using since 2008. It's available [here](https://github.com/johannessonlab/HetVPaper/blob/master/SNPpop/data/Podan2/Podan2_AssemblyScaffoldsmt.fa) and in the legacy websites of JGI. There is a [newer version of the S+ reference genome](https://doi.org/10.1186/s12864-022-09085-4), but using this version allows full compatibility with out previous work.

I updated my latest annotation of Podan2 with the manually curated models, resulting in the gff3 file `Podan2.nice-3.02.gff3`.


## Run pipeline

Get into the folder:

    $ cd IntroRegions

To get an idea of how the pipeline looks like we can make a rulegraph:

    $ snakemake --snakefile IntroRegions.smk --rulegraph | dot -Tpng > rulegraph.png

To check that the files for the pipeline are in order:

    $ snakemake --snakefile IntroRegions.smk -pn

There are many ways of running the pipeline. In this case I'm using the profile file defined above:

    $ screen -R interpro
    $ mamba activate annotation
    $ module load interproscan/5.75-106.0
    $ snakemake --snakefile IntroRegions.smk &> snakemake.log &
    [1] 3080255

If for some reason the profile is not working, one can do (but the InterProScan rule will not run as a job).

	$ snakemake --snakefile IntroRegions.smk -j4 --use-conda

## Results

The pipeline will produce some tables and figures.

- `enrichmnent_allplots.pdf` - a plot of the top 10 domains in the Enrichment analyses for various types of genome regions (e.g. the areas introgressed in either *P. comata* or *P. pauciseta*, or both). I didn't include this figure in the paper in the end.

- `enrichmnent_OddsVsPval_both.pdf` - Final enrichment analysis figure, part of the main Figure 7, for both *P. comata* or *P. pauciseta* introgressed regions.

- `enrichmnent_OddsVsPval_miss.pdf` - Final enrichment analysis figure, for regions marked as missing data in both species. Not used in the paper.

- `enrichmnent_both.txt` - Enrichment values of all domains in the `Podan2.nice-3.02.gff3` annotation ordered by p-value, for both *P. comata* or *P. pauciseta* introgressed regions.

- `enrichmnent_com.txt` - Enrichment values of all domains in the `Podan2.nice-3.02.gff3` annotation ordered by p-value, for *P. comata* introgressed regions.

- `enrichmnent_pau.txt` - Enrichment values of all domains in the `Podan2.nice-3.02.gff3` annotation ordered by p-value, for *P. pauciseta* introgressed regions.

- `enrichmnent_miss.txt` - Enrichment values of all domains in the `Podan2.nice-3.02.gff3` annotation ordered by p-value, for regions marked as missing data in both species.

- `enrichmnent_taj.txt` - Enrichment values of all domains in the `Podan2.nice-3.02.gff3` annotation ordered by p-value, for regions with Tajima's D > 2 in the window-based analysis of [Ament-Velásquez et al. (2022) NEE](https://doi.org/10.1038/s41559-022-01734-x). Not used in the paper.

- `Genes_NLRdomains_coords.bed` - Coordinates of genes with NLR-associated domains.

- `NLRdomains_chrs.png` - Figure S10, chromosomal distribution of some NLR-associated domains (as annotated by InterProScan) in relation to introgressed regions from *P. pauciseta* or *P. comata*.

- `NLRgenes_chrs.pdf` - Part of main Figure 7, chromosomal distribution of P. anserina NLRs in relation to introgressed regions from *P. pauciseta* or *P. comata*.

I put together the final Figure 7 in Inkscape.

