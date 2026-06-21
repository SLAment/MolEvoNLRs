# Phylogenies in figures

The alignments and maximum likelihood phylogenies of *het-s*, *nwd2*, and *PaPlp1* are provided here. Notice that these alignments are derived from the [alignments in the `data` directory](https://github.com/SLAment/MolEvoNLRs/tree/main/1_NLRvsRandomGenes/data), and the trees are produced directly by the pipeline too. However, you can re-make them independently. For example with het-s:

Make an environment with the chosen IQ-TREE version, in this case 2.3.6 (but there are much newer versions now!):

	$ mamba create -n iqtree236 -c bioconda iqtree=2.3.6
	$ mamba activate iqtree236

Run IQ-TREE to get 100 standard bootstraps with automatic selection of evolutionary model:

	$ iqtree -s het-s.fa -m MFP -seed 1234 -b 100 -nt 4 -pre IQTREE/het-s_GBE/het-s

The alignments in nexus format and their phylogenies in newick format used for Figures 1 and 2 are also included here.

- Figure 1: alignment `NWDfam_ref_aa_corrected_extended_clean_GxxxxGKT_hetEextend_subset_outgroup.nxs` and tree `NACHT_NWDs.tre`
- Figure 2: alignment `NB-ARCfam_ref_aa.nxs` and tree `NB-ARC_NLRs.tre`.

You can open nexus files, manipulate and transform them into fasta file using for example [SeaView](https://doua.prabi.fr/software/seaview). The advantage of the nexus files is that you can define subsets within the alignment. In this case, the subsets that say `core` where used to make the trees.