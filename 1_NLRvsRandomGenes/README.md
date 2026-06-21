# NLRvsRandomGenes - A pipeline to compare NLRs from randomly selected genes in the *Podospora anserina* species complex

I manually produced alignments of a set of NLRs in *P. anserina* and related species, as well as randomly selected genes. In this pipeline I try to get some numbers to compare their dynamics, under the expectation that NLRs evolve faster and under a birth-death dynamic.

The pipeline is not computationally demanding, so I ran it locally on a MacOS.

## Input files

All the input files should be specified in a configuration file in yaml format, called `config/config.yaml`. The involved files are as follows:

- There are three files with a list of genes IDs (in a single column), each for the three kinds of genes: HIC NLRs, LIC NLRs and random genes. these are all available in the `data` folder.

- A csv file with the different classifications of conservation state per gene and per species, called `EvolutionAllGenes.csv` also in `data`.

- A txt file derived from `../0_InterProPodo/reports/GenePFAMCounts.txt` but with the species names added. It contains the PFAM counts of NACHT and NB-ARC genes for the strains with chromosome-level assemblies.

- The path of a folder contaning three other folders, each with the alignments in fasta format of each NLR and random gene: 
	- Folder `nlrs` contains the LIC NLRs.
	- Folder `nlrsHIC` contains the HIC NLRs.
	- folder `random` contains the randomly selected genes.

Each gene has a full alignment including introns, and another with only the coding sequences (CDS). The HIC NLRs also include alignments with the HIC repeats removed except for the first one (`noHICrepts`) or the full sensor domain removed till the stop codon (`_noReptDomain`).

- A bed file of the regions that match 5 Kbp windows (overlapping, steps of 1 Kbp) with Tajima's D >= 2 in the Wageningen population of *Podospora anserina*, from [Ament-Velásquez et al. (2022) NEE](https://doi.org/10.1038/s41559-022-01734-x). This file is to calculate enrichment of particular domains, but it wasn't used in the paper.

- The manually curated annotation file `data/Podan2.nice-3.02.gff3`. The gene IDs in the list of genes correspond to this file.

- The assembly Podan2, corresponding to the second iteration of the reference genome of the strain S+ ([Espagne et al. 2008](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-5-r77)), produced by JGI's Mycocosm. The sequence seems to be identical to the latest version Podan3, but the newest annotation in JGI discarded the original gene codes (e.g., Pa_Z_XXXX, where Z is the chromosome and XXXX is the gene number) that the community has been using since 2008. It's available [here](https://github.com/johannessonlab/HetVPaper/blob/master/SNPpop/data/Podan2/Podan2_AssemblyScaffoldsmt.fa) and in the legacy websites of JGI. There is a [newer version of the S+ reference genome](https://doi.org/10.1186/s12864-022-09085-4), but using this version allows full compatibility with our previous work.

- The annotation (in gff3 format) of the transposable elements (from [Ament-Velásquez et al. 2024](https://doi.org/10.1093/gbe/evae034)), produced with RepeatMasker using the library PodoTE-1.00 ([Vogan et al. 2021](https://genome.cshlp.org/content/31/5/789)), available [here](https://datadryad.org/dataset/doi:10.5061/dryad.1vhhmgr0j). The file is called `Podan2.repeatmasker.gff3`.

- The size of the flanks of each gene to calculate the TE coverage (`winsize`). In the paper I used `winsize: 2000`.

- Several ploting scripts in the `scripts` folder: 
	* `NLRdynamics.R`
	* `PiDistNLRvsRandom.R`
	* `PlotBalancingSelection.R`


The pipeline will automatically download my script `fasta2axt.py` to format the alignments into the input format (axn) for kakscalculator.

This is how the `config.yaml` file should look like:

```yaml
# Input metadata
randomfile: "data/RandomGenesList.txt"
HICfile: "data/HICGenesList.txt"
NLRfile: "data/LICGenesList.txt"
geneevo: "data/EvolutionAllGenes.csv"
expansions: "data/GeneFamilyExpansions.txt" # PFAM counts of the all Podospora species

path2alignments: "data"

# Tajima calculated in 10Kb windows with less than 50% missing data
FULLTajima: "../IntroRegions/data/WageningenTajimaMin5Kb.bed"

ANNOTATION: "../IntroRegions/data/Podan2.nice-3.02.gff3"
refgenome: "path/to/Podan2_AssemblyScaffoldsmt.fa"
TEs: "path/to/Podan2.repeatmasker.gff3"

winsize: 2000

PlotBalancingSelection: "scripts/PlotBalancingSelection.R"
PiDistNLRvsRandom: "scripts/PiDistNLRvsRandom.R"
NLRdynamics: "scripts/NLRdynamics.R"
```
**IMPORTANT**

If you want to add a new strain or species, then the actual code in the snakefile `NLRvsRandomGenes.smk` has to be modified to handle the taxonomy and names correctly. That is defined at the top of the snakefile as a python dictionary that looks like so:

```python
straindic = {"anserina": ["PaSp", "Podan2", "PaWa21m", "PaWa28m", "PaWa46p", "PaWa53m", "PaWa58m", "PaWa63p", "PaWa87p", "PaWa100p", "PaWa137m", "PaZp", "PaTgp", "PaYp"], 
			 "pauciseta": ["CBS237.71m", "CBS333.63p", "CBS451.62p"], 
			 "comata": ["PcTdp", "PODCO", "PcWa131m", "PcWa132p", "PcWa133m", "PcWa139m"], 
			 "bellae-mahoneyi": ["CBS112042p"], 
			 "pseudoanserina": ["CBS124.78p", "CBS253.71p"],
			 "pseudopauciseta": ["CBS411.78m"], 
			 "pseudocomata": ["CBS415.72m", "IMI230595m"]}
```

The order of the strain matters: the first strain will be given priority to use it as a reference to make trees and calculate dN/dS. For example, whenever "PaSp" (i.e., the strain S+) is present in an alignment, that one will be used as a *P. anserina* reference, otherwise "PaWa21m" will be used (Podan2 is the same strain as PaSp but the first one comes from JGI and the second from a SPAdes assembly of S+ illumina data).

Say that tomorrow a new strain of *P. pseudopauciseta* is found and sequenced. Then we would add it as `"pseudopauciseta": ["CBS411.78m", "new_strain_ID"]` in the dictionary above (and in all the input alignments in the `data` folder, of course).

## Environment

### With conda

At the time of working some conda packages (egglib, really) were not available on the channels for M1 architectures (the Macbook chip) so I had to use `CONDA_SUBDIR=osx-64` before calling mamba. But this won't be necessary for you if you are working on Linux!

Because I'm a rebel, I had to remove the strict channel priority for this to work.

	$ conda config --set channel_priority true

Make the environment:

	$ CONDA_SUBDIR=osx-64 mamba create -n nlrstats -c bioconda snakemake-minimal=7.32.4 biopython=1.85 bedtools=2.31.1 pandas=2.3.2 iqtree=3.0.1 conda-forge::ete3=3.1.3 gffutils=0.13 wget=1.21.4 egglib=3.1.0 kakscalculator2=2.0.1

It will take a good while to resolve the environment. Be patient.

The R scripts used within the pipeline have their own mini environment. However, when running conda, it will try to use normal mamba, without `CONDA_SUBDIR=osx-64`. So it will throw an error about the environment folder already existing. This is a [bug](https://github.com/mamba-org/mamba/issues/2736) associated to the old snakemake version. Using newer Snakemake versions at the time of writing was not an option because of conflicts with the other packages. But there is a workaround below.

(After you're done running the pipeline, put the channel priority back as strict.)

	$ conda config --set channel_priority strict

### With pixi

An easier alternative is to use [pixi](https://pixi.prefix.dev/latest/workspace/environment/), which is like an upgraded conda. Once you have pixi installed in your environment you can go to the `envs` folder where you will find the `pixi.lock` and `pixi.toml` files, which you need to reproduce the environment. Just do:

	$ cd envs
	$ pixi install
	$ cd ..

The environment should be active!

## Running the pipeline locally

Go to working directory and activate the environment.

	$ mamba activate nlrstats

First, to get an idea of how the pipeline looks like we can make a rulegraph:

	$ snakemake --snakefile NLRvsRandomGenes.smk --rulegraph | dot -Tpng > rulegraph.png

![rulegraph](rulegraph.png "rulegraph")

To check that the files for the pipeline are in order:

	$ snakemake --snakefile NLRvsRandomGenes.smk -pn

Let's run it for real:

	$ snakemake --snakefile NLRvsRandomGenes.smk --use-conda -pj8

This is when the error happens. It will look like this:

	Building DAG of jobs...
	Your conda installation is not configured to use strict channel priorities. This is however crucial for having robust and correct environments (for details, see https://conda-forge.org/docs/user/tipsandtricks.html). Please consider to configure strict priorities by executing 'conda config --set channel_priority strict'.
	Creating conda environment envs/plot.yaml...
	Downloading and installing remote packages.
	CreateCondaEnvironmentException:
	Could not create conda environment from /Users/lorena/Library/CloudStorage/Dropbox/VRwork/Manuscripts/11_MolEvoINWDs/GitHub/NLRvsRandomGenes/envs/plot.yaml:
	Command:
	mamba env create --quiet --file "yourpath/GitHub/NLRvsRandomGenes/.snakemake/conda/1e2bf99ec593c0b1bf8fe0b9f3de0692_.yaml" --prefix "yourpath/GitHub/NLRvsRandomGenes/.snakemake/conda/1e2bf99ec593c0b1bf8fe0b9f3de0692_"
	Output:
	error    libmamba Non-conda folder exists at prefix
	critical libmamba Aborting.

Simply re-run the command yourself with the ID given by Snakemake in the error and `CONDA_SUBDIR=osx-64` before it:

	$ CONDA_SUBDIR=osx-64 mamba env create --quiet --file "yourpath/GitHub/NLRvsRandomGenes/.snakemake/conda/1e2bf99ec593c0b1bf8fe0b9f3de0692_.yaml" --prefix "yourpath/GitHub/NLRvsRandomGenes/.snakemake/conda/1e2bf99ec593c0b1bf8fe0b9f3de0692_"
	Confirm changes: [Y/n] y

That will take some time... but now you should be able to run it for real:

	$ snakemake --snakefile NLRvsRandomGenes.smk --keep-going --use-conda -pj8


## Results

The pipeline will produce a number of figures in the folder `results`.

- Figure 3: Partially, it has almost all the elements of the figure but I added the phylogeny in Inkscape
- Figure 5: Tajima's D and test of monophyly
- Figure S1: About mutations in the 3 gene classes and the RIP patterns
- Figure S2: Nucleotide diversity in non-synonymous sites and also Nucleotide diversity in synonymous sites for genes classified as conserved vs degraded.
- Figure S3: dN/dS ratios between Podospora anserina and the other Podospora species, but only for degraded genes.

The pipeline also produces additional figures that didn't make it to the paper in the folder `figures`.

