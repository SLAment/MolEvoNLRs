# Phylogeny of *het-s* and *nwd2*

The alignments and maximum likelihood phylogenies of *het-s* and *nwd2* are provided here. The alignment of *nwd2* is found [here]() as part of the dataset used for other figures in the paper.

To make a tree, for example with het-s:

Make an environment with the chosen IQ-TREE version, in this case 2.3.6 (but there are much newer versions now!):

	$ mamba create -n iqtree236 -c bioconda iqtree=2.3.6
	$ mamba activate iqtree236

Run IQ-TREE to get 100 standard bootstraps with automatic selection of evolutionary model:

	$ iqtree -s het-s.fa -m MFP -seed 1234 -b 100 -nt 4 -pre IQTREE/het-s_GBE/het-s

