# Phylogenies in supplementary figures

The alignments and maximum likelihood phylogenies of *het-s*, *nwd2*, and *PaPlp1* are provided here. Notice that these alignments are derived from the [alignments in the `data` directory](https://github.com/SLAment/MolEvoNLRs/tree/main/1_NLRvsRandomGenes/data), and the trees are produced directly by the pipeline too. However, you can re-make them independently. For example with het-s:

Make an environment with the chosen IQ-TREE version, in this case 2.3.6 (but there are much newer versions now!):

	$ mamba create -n iqtree236 -c bioconda iqtree=2.3.6
	$ mamba activate iqtree236

Run IQ-TREE to get 100 standard bootstraps with automatic selection of evolutionary model:

	$ iqtree -s het-s.fa -m MFP -seed 1234 -b 100 -nt 4 -pre IQTREE/het-s_GBE/het-s

