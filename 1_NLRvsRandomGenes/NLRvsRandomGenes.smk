# -*- snakemake -*-

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools
import egglib
from Bio import Phylo
import os
import statistics
from ete3 import Tree

import gffutils # CONDA_SUBDIR=osx-64 mamba install -c bioconda gffutils=0.13 

### NLRvsRandomGenes
#############################################################################
# A pipeline to compare NLRs from randomly selected genes in the Podospora anserina species complex
#############################################################################

# https://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#trees

# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2024/10/16
# ++++++++++++++++++++++++++++++++++++++++++++++

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"

randomfile = config["randomfile"]
HICfile = config["HICfile"]
NLRfile = config["NLRfile"]
geneevo = config["geneevo"]

path2alignments = config["path2alignments"]

FULLTajima = config["FULLTajima"]

ANNOTATION = config["ANNOTATION"]
refgenome = config["refgenome"]
TEs = config["TEs"]

winsize = config["winsize"]

# Scripts
PiDistNLRvsRandom = config["PiDistNLRvsRandom"]
PlotBalancingSelection = config["PlotBalancingSelection"]
NLRdynamics = config["NLRdynamics"]

# -------------------------------------------------

straindic = {"anserina": ["PaSp", "Podan2", "PaWa21m", "PaWa28m", "PaWa46p", "PaWa53m", "PaWa58m", "PaWa63p", "PaWa87p", "PaWa100p", "PaWa137m", "PaZp", "PaTgp", "PaYp"], 
			 "pauciseta": ["CBS237.71m", "CBS333.63p", "CBS451.62p"], 
			 "comata": ["PcTdp", "PODCO", "PcWa131m", "PcWa132p", "PcWa133m", "PcWa139m"], 
			 "bellae-mahoneyi": ["CBS112042p"], 
			 "pseudoanserina": ["CBS124.78p", "CBS253.71p"],
			 "pseudopauciseta": ["CBS411.78m"], 
			 "pseudocomata": ["CBS415.72m"]}

anserina_strains = straindic["anserina"]
sppstrains = list(straindic.keys())
sppstrains_sans_anserina = sppstrains
sppstrains_sans_anserina.remove("anserina")

# -------------------------------------------------

# Prepare the lists
random_selection = [line.rstrip("\n") for line in open(randomfile, 'r')]
random_selection.sort()

HIC_selection = [line.rstrip("\n") for line in open(HICfile, 'r')]
NLR_selection = [line.rstrip("\n") for line in open(NLRfile, 'r')]

focalgenes = random_selection + HIC_selection + NLR_selection

# ----------
rule all:
	input:
		"results/Fig3_NLRvsRandomDynamics_raw.pdf",
		"results/Fig4_PopPhyloStats.png", # This one requires a change in the R code to plot pdf with cairo device
		"results/FigS1_RIPmuts.png",
		"results/FigS2_SupDivStats.png",
		"results/FigS3_dNdSdegraded.png",
		"results/Fig5_BalSelection.png", # Figure 5 showing Tajima and monophyly tests

# --- Get population statistics of P. anserina ---

def calculate_gc_content(seq):
	"""
	Calculate GC content for a single sequence, ignoring gaps ('-').
	
	:param seq: The sequence (Bio.Seq object or string) to analyze.
	:return: GC content as a fraction (float).
	"""
	seq_upper = seq.upper()
	gc_count = 0
	valid_bases = 0
	
	for base in seq_upper:
		if base != '-':
			valid_bases += 1
			if base in ['G', 'C']:
				gc_count += 1

	if valid_bases > 0:
		gc_content = gc_count / valid_bases
	else:
		gc_content = 0
	
	return gc_content

def GC_content_alignment(alignment):
	# Calculate GC content for each sequence in the alignment
	gc_contents = []
	for record in alignment:
		gc_content = calculate_gc_content(record.seq)
		gc_contents.append(gc_content)

	# Calculate average GC content
	if gc_contents:
		avg_gc_content = sum(gc_contents) / len(gc_contents)
	else:
		avg_gc_content = 0

	return avg_gc_content

def get_input_fasta_cds(wildcards):
	if wildcards.gene in random_selection:
		return f"{path2alignments}/random/{wildcards.gene}_CDS.fa"
	elif wildcards.gene in HIC_selection:
		return f"{path2alignments}/nlrsHIC/{wildcards.gene}_noHICrepts_CDS.fa"
		# return f"{path2alignments}/nlrsHIC/{wildcards.gene}_noHICrepts_CDS_noReptDomain.fa"
	elif wildcards.gene in NLR_selection:
		return f"{path2alignments}/nlrs/{wildcards.gene}_CDS.fa"


def calc_stats(alignment, propmissing = 0, skipstop = False):
	""" Modified from https://github.com/topfm/selection-stats/blob/main/piNpiS.py """

	## use EggLib3 to handle the new alignment
	aln = egglib.io.from_fasta(alignment, egglib.alphabets.DNA)
	aln_codons = egglib.tools.to_codons(aln) # This method overwrites the current instance. After the conversion, there will one site per codon.

	aln_div = egglib.stats.CodingDiversity(aln_codons, multiple_alleles=True, multiple_hits=True, max_missing= propmissing, skipstop = skipstop)
	# multiple_alleles		include coding sites with more than two alleles (regardless of whether mutations hit the same position within the triplet). If there are more than two different codons, they must either encode for the same amino acid or all encode for different amino acids (otherwise the site is rejected). 
	# multiple_hits		 include coding sites for which more than one of the three positions has changed (regardless of the number of alleles). If there are more than two alleles, the option multiple_alleles is also considered.
	# max_missing		maximum relative proportion of missing data (per codon site) to allow (including stop codons if skipstop if true). By default, all codon sites with any missing data are excluded.
	# skipstop		if True, stop codons are treated as missing data and skipped (default). If so, potential mutations to stop codons are not taken into account when estimating the number of non-synonymous sites. Warning (this may be counter-intuitive): it actually assumes that stop codons are not biologically plausible and considers them as missing data. On the other hand, if skipstop is False, it considers stop codons as if they were valid amino acids. This option has no effect if raise_stop is True.
	
	# Notice: Missing data include all ambiguity characters and **alignment gaps**. Sites not passing this threshold are ignored for all analyses!

	# Make a new ComputeStats object and get stats
	statDict = {}
	cs = egglib.stats.ComputeStats(multi_hits=True)
	# print(cs.list_stats())
	cs.add_stats('Pi','lseff', 'nseff', 'D')
	# nseff: Average number of effectively analyzed sequences
	# lseff: Number of effectively analyzed sites
	# D: Tajima's D

	statsall = cs.process_align(aln, max_missing = propmissing)
	statsS = cs.process_sites(aln_div.sites_S)
	statsNS = cs.process_sites(aln_div.sites_NS)

	statDict['nseff'] = statsall['nseff']
	statDict['lseff'] = statsall['lseff']
	
	if statsall['D'] == None:
		statDict['D'] = 'NA' # This happens if there is no variation in the alignment
	else:
		statDict['D'] = f"{statsall['D']:.6f}"

	try:
		statDict['pi'] = float(statsall['Pi'])/statsall['lseff']
	##sometimes Pi = None because there isn't enough information because the alignment or codon alignment didn't pass the max_missing threshhold.
	##This is when a gene has really bad alignment. Stats aren't able to be calculated, so we just get "None."
	except TypeError:	
		statDict['pi'] = None

	# Non-synonymous sites
	statDict['num_sites_NS'] = round(aln_div.num_sites_NS) # make it an integer
	try:
		statDict['piN'] = float(statsNS['Pi'])/aln_div.num_sites_NS
		statDict['sites_NS'] = round(statsNS['Pi'])  # make it an integer
	except TypeError: # apparently I also get None if there is 0 variation in the alignment
		statDict['piN'] = 0 
		statDict['sites_NS'] = 0 
	
	# Synonymous sites
	statDict['num_sites_S'] = round(aln_div.num_sites_S)  # make it an integer
	try:
		statDict['piS'] = float(statsS['Pi'])/aln_div.num_sites_S
		statDict['sites_S'] = round(statsS['Pi'])  # make it an integer
	except TypeError:
		statDict['piS'] = 0
		statDict['sites_S'] = 0 
	##making sure that 0/0 = 0 and 0/1= 0 and 0/None = None, in terms of piNpiS.

	try: 
		statDict['piNpiS'] = statDict['piN']/statDict['piS']
	except ZeroDivisionError:
		if statDict['piN'] == 0:
			statDict['piNpiS'] = 0
		else:
			statDict['piNpiS'] = None

	return statDict


noPagenes_chrs = {'Pa_0_80': 6, 'nwd5': 5, 'nwd2': 3, 'nwd3': 5, 
				  'nwd1': 3, 'het-e': 4, 'hnwd3': 7, 'nwd6': 5,
				  'nwdp-2': 3, 'het-d': 2, 'hnwd1': 4, 'het-r': 2,
				  'nwd4': 5, "PaPnt1": 5, "PaFnt1": 3, "PaPlp1": 1}

rule subset_aln_anserina:
	""" Get only the anserina strains into a fasta file """
	input:
		fasta = get_input_fasta_cds,
		listar = randomfile, # to monitor that it's the latest version
		listah = HICfile,
	output:
		fasta = "alignments/anserina/{gene}_CDS_anserina.fa"
	run:
		# Load alignment from a fasta file
		alignment = AlignIO.read(input.fasta, "fasta")

		# Get just the anserina sequences
		anserina_seqs = []
		for strain in anserina_strains:
			for seq in alignment:
				if strain in seq.id: 
					anserina_seqs.append(seq.id)
					break

		# Write a new alignment fasta file with just Podospora anserina strains
		new_alignment = MultipleSeqAlignment([seq for seq in alignment if seq.id in anserina_seqs])

		with open(output.fasta, "w") as handle: 
			SeqIO.write(new_alignment, handle, "fasta")


rule popstats_anserina:
	""" Get pop statistics from the anserina strains """
	## Modified from https://github.com/topfm/selection-stats/blob/main/piNpiS.py
	input:
		fasta = "alignments/anserina/{gene}_CDS_anserina.fa",
		listar = randomfile, # to monitor that it's the latest version
		listah = HICfile,
	output:
		stats = "temp/stats/{gene}_stats_anserina.txt"
	run:
		# Load alignment from a fasta file
		alignment = AlignIO.read(input.fasta, "fasta")

		# Infer chromosome
		if wildcards.gene in noPagenes_chrs.keys():
			chrom = noPagenes_chrs[wildcards.gene]
		elif 'Pa_' in wildcards.gene:
			chrom = wildcards.gene.split('_')[1]
		else:
			chrom = "NA"

		# Type of gene
		if wildcards.gene in random_selection:
			genetype = "Random"
		elif wildcards.gene in HIC_selection:
			genetype = "NLR_HIC"
		elif wildcards.gene in NLR_selection:
			genetype = "NLR"

		# Calculate the average pairwise identity and average GC content
		average_gc_content = GC_content_alignment(alignment)

		## use EggLib3 to handle the alignment
		alignDict = calc_stats(input.fasta, propmissing = 1, skipstop = False)

		# Sometimes piS is 0
		if alignDict['piNpiS'] == None:
			piNpiS = "NA"
		else:
			piNpiS = f"{alignDict['piNpiS']:.6f}"

		# Make the report
		ofile = open(output.stats, 'w')
		ofile.write(f"{wildcards.gene}\t{genetype}\t{chrom}\t{average_gc_content:.6f}\t{alignDict['D']}\t{alignDict['pi']:.8f}\t{alignDict['piN']:.6f}\t{alignDict['piS']:.6f}\t{piNpiS}\t{alignDict['sites_S']}\t{alignDict['num_sites_S']}\t{alignDict['sites_NS']}\t{alignDict['num_sites_NS']}\t{len(alignment)}\t{alignDict['nseff']:.6f}\t{alignDict['lseff']}\n")

rule get_all_genes:
	input:
		reports = expand("temp/stats/{gene}_stats_anserina.txt", gene = focalgenes),
		listar = randomfile, # to monitor that it's the latest version
		listah = HICfile, 
	output:
		"reports/Orthologs_stats_anserina.txt"
	shell:
		"printf 'Ortholog\tType\tChr\tavgGC\tTajimasD\tPi_anserina\tpiN\tpiS\tpiNpiS\tsites_S\tnum_sites_S\tsites_NS\tnum_sites_NS\tN_seqs\tnseff\tlseff\n' > {output}; "
		"cat {input.reports} >> {output}"
	# nseff: Average number of effectively analyzed sequences
	# lseff: Number of effectively analyzed sites


# --- Plot NLRs vs random genes dynamics ---

rule Plot_NLRdynamics:
	input:
		dists = "reports/Orthologs_stats_anserina.txt", # Population genetics statistics
		geneevo = geneevo, # proportion covered by transposons in 2kb flanks
	output:
		dynamics = "results/Fig3_NLRvsRandomDynamics_raw.pdf",
		curation = "figures/Curation_proportions.pdf"
	conda: 
		"envs/plot.yaml"
	script:
		NLRdynamics


# --- Get TE% in reference genome ---

rule get_gene_intervals:
	""" Get the gene coordinates and infer the flanking coords based on that """ 
	input:
		refgenome = refgenome,
		genes = ANNOTATION,
	output:
		bed = "reports/FlanksGenesCoords.bed"
	params:
		winsize = winsize
	run:
		ofile = open(output.bed, 'w')
		for line in open(input.genes, 'r'):
			if '\tgene\t' in line:
				for gene in focalgenes:
					if f"={gene};" in line:
						tabs = line.rstrip("\n").split("\t")
						contig = tabs[0]
						start = int(tabs[3])
						end = int(tabs[4])
						ofile.write( f"{contig}\t{start - params.winsize}\t{start - 1}\t{gene}\n" )
						ofile.write( f"{contig}\t{end + 1}\t{end + params.winsize}\t{gene}\n" )
						break

rule intersect_bed:
	input:
		bed = "reports/FlanksGenesCoords.bed",
		TEs = TEs,
	output:
		bed = "reports/FlanksTEintersect.bed"
	shell:
		"bedtools intersect -a {input.bed} -b {input.TEs} >> {output.bed}"

rule calc_TEcoverage:
	input:
		bed = "reports/FlanksTEintersect.bed"
	output:
		tab = "reports/FlanksTEcontent.tab"
	params:
		winsize = winsize
	run:
		genesdic = {}
		for line in open(input.bed, 'r'):
			contig, start, end, gene = line.rstrip("\n").split("\t")
			area = int(end) - int(start)
			if gene in genesdic.keys():
				genesdic[gene] += area
			else:
				genesdic[gene] = area

		ofile = open(output.tab, 'w')
		ofile.write("Ortholog\tType\tTEcov\n")
		for gene in focalgenes:
			if gene not in genesdic.keys():
				TEcov = 0
			else:
				TEcov = genesdic[gene]/(params.winsize*2)

			if gene in random_selection:
				genetype = "Random"
			elif gene in HIC_selection:
				genetype = "NLR_HIC"
			elif gene in NLR_selection:
				genetype = "NLR"

			ofile.write(f"{gene}\t{genetype}\t{TEcov:.6f}\n")

# --- Get phylogenetic distances ---

rule subset_aln_spp:
	""" Get only the one strain per species into a fasta file """
	## Modified from https://github.com/topfm/selection-stats/blob/main/piNpiS.py
	input:
		fasta = get_input_fasta_cds,
		listar = randomfile, # to monitor that it's the latest version
		listah = HICfile,
	output:
		fasta = "alignments/sppcomplex/{gene}_CDS_spp.fa"
	log:
		file = "logs/subset_aln_spp/{gene}.log"
	run:
		records_dict = SeqIO.to_dict(SeqIO.parse(open(input.fasta, 'r'), "fasta"))

		nomsdic = {}
		for record in records_dict.keys():
			if '__' in record:
				Strain = record.split('__')[1].split('_')[0]
				nomsdic[Strain] = record
			else:
				Strain = record.split('_')[0]
				nomsdic[Strain] = record

		# print(nomsdic)

		anserina = ""
		comata = ""
		pauciseta = ""
		pseudoanserina = ""

		finalrepresentatives = {}
		for species in straindic.keys():
			for strain in straindic[species]:
				if strain in nomsdic.keys():
					finalrepresentatives[strain] = species
					# finalrepresentatives[species] = strain
					break

		# Keep track of who was chosen somewhere
		sr = str(finalrepresentatives)
		logf = open(log.file, 'w')
		logf.write(sr)

		output_handle = open(output.fasta, "w")
		for strain in nomsdic.keys():
			if strain in finalrepresentatives.keys():
				seqname = nomsdic[strain]
				record = records_dict[seqname]
				record.id = finalrepresentatives[strain] # change the name of the strain to the species
				record.description = ""
				SeqIO.write(record, output_handle, "fasta")

rule IQTreePerGene:
	""" Run IQTree for each gene """
	input:
		aln = "alignments/sppcomplex/{gene}_CDS_spp.fa"
	output:
		tree = "iqtree/{gene}/{gene}_CDS_spp.treefile"
	threads: 2
	params:
		bootstraps = 1000, # UFBoot
	run:
		# Check that there are at least 3 species in there
		records_dict = SeqIO.to_dict(SeqIO.parse(input.aln, "fasta"))

		if len(records_dict) > 3: # you can't make a tree with two taxa...
			cmd = f"iqtree -s {input.aln} -m MFP -seed 1234 -bb {params.bootstraps} -nt {threads} -bnni -pre 'iqtree/{wildcards.gene}/{wildcards.gene}_CDS_spp' --keep-ident -redo"
			shell(cmd)
		elif len(records_dict) > 2: # no bootstraps 
			cmd = f"iqtree -s {input.aln} -m MFP -seed 1234 -nt {threads} -bnni -pre 'iqtree/{wildcards.gene}/{wildcards.gene}_CDS_spp' --keep-ident -redo"
			shell(cmd)
		else:
			ofile = open(output.tree, "w") # make a dummy empty file

rule get_sppdistances:
	""" what is the distance between P. anserina and the other species for each gene """
	input:
		nw = "iqtree/{gene}/{gene}_CDS_spp.treefile"
	output:
		txt = "temp/treedistances/{gene}_dist.txt"
	run:
		sppstrains = list(straindic.keys())

		# Check if the file is empty
		if os.stat(input.nw).st_size == 0:
			geneline = wildcards.gene + "\tNA"*(len(sppstrains) - 1) + '\n'
		else:
			tree = Phylo.read(input.nw, "newick")

			# Find all terminal nodes (leaves)
			terminal_nodes = tree.get_terminals()

			# Calculate distances from the reference leaf to each other terminal node
			distances = {}
			for leaf in terminal_nodes:
				if leaf.name != "anserina":  # Skip the reference leaf itself
					# Calculate the distance from the reference to the current leaf
					distance = tree.distance("anserina", leaf.name)
					distances[leaf.name] = distance


			geneline = f"{wildcards.gene}"
			for strain in sppstrains_sans_anserina: 
				if strain in distances.keys():
					geneline += f"\t{distances[strain]}"
				else:
					geneline += "\tNA"
			geneline += '\n'	

		ofile = open(output.txt, "w")
		ofile.write(geneline)	

rule get_all_distances:
	input:
		reports = expand("temp/treedistances/{gene}_dist.txt", gene = focalgenes),
		listar = randomfile, # to monitor that it's the latest version
		listah = HICfile, 
	output:
		"reports/Orthologs_Distances.txt"
	params:
		header = "Ortholog\t" + '\t'.join(sppstrains_sans_anserina)
	shell:
		"echo '{params.header}' >> {output}; "
		"cat {input.reports} >> {output}"

## ------ get crude dN/dS -------

def calc_stats_spp(alignment, propmissing = 0, skipstop = False):
	""" Modified from https://github.com/topfm/selection-stats/blob/main/piNpiS.py """

	## use EggLib3 to handle the new alignment
	aln = egglib.io.from_fasta(alignment, egglib.alphabets.DNA)

	# Find the reference
	# anserinaseq = aln.find(anserinaref)

	sample_names = [seq.name for seq in aln] # a list of all the names in the alignment
	DicOfdicst = {}

	for seq in aln:
		if seq.name != "anserina":
			# Subset the alignment to have only those two samples
			subset_aln = aln.subset([sample_names.index(name) for name in [seq.name, "anserina"]])

			aln_codons = egglib.tools.to_codons(subset_aln) # This method overwrites the current instance. After the conversion, there will one site per codon.

			aln_div = egglib.stats.CodingDiversity(aln_codons, multiple_alleles=True, multiple_hits=True, max_missing= propmissing, skipstop = skipstop)
			# multiple_alleles		include coding sites with more than two alleles (regardless of whether mutations hit the same position within the triplet). If there are more than two different codons, they must either encode for the same amino acid or all encode for different amino acids (otherwise the site is rejected). 
			# multiple_hits		 include coding sites for which more than one of the three positions has changed (regardless of the number of alleles). If there are more than two alleles, the option multiple_alleles is also considered.
			# max_missing		maximum relative proportion of missing data (per codon site) to allow (including stop codons if skipstop if true). By default, all codon sites with any missing data are excluded.
			# skipstop		if True, stop codons are treated as missing data and skipped (default). If so, potential mutations to stop codons are not taken into account when estimating the number of non-synonymous sites. Warning (this may be counter-intuitive): it actually assumes that stop codons are not biologically plausible and considers them as missing data. On the other hand, if skipstop is False, it considers stop codons as if they were valid amino acids. This option has no effect if raise_stop is True.
			
			# Notice: Missing data include all ambiguity characters and **alignment gaps**. Sites not passing this threshold are ignored for all analyses!

			# Make a new ComputeStats object and get stats
			statDict = {}
			cs = egglib.stats.ComputeStats(multi_hits=True)
			cs.add_stats('Pi','lseff', 'nseff')
			# nseff: Average number of effectively analyzed sequences
			# lseff: Number of effectively analyzed sites

			statsall = cs.process_align(subset_aln, max_missing= propmissing)
			statsS = cs.process_sites(aln_div.sites_S)
			statsNS = cs.process_sites(aln_div.sites_NS)

			statDict['nseff'] = statsall['nseff']
			statDict['lseff'] = statsall['lseff']

			try:
				statDict['pi'] = float(statsall['Pi'])/statsall['lseff']
			##sometimes Pi = None because there isn't enough information because the alignment or codon alignment didn't pass the max_missing threshhold.
			##This is when a gene has really bad alignment. Stats aren't able to be calculated, so we just get "None."
			except TypeError:	
				statDict['pi'] = None

			try:
				statDict['piN'] = float(statsNS['Pi'])/aln_div.num_sites_NS
			except TypeError:
				statDict['piN'] = 0

			try:
				statDict['piS'] = float(statsS['Pi'])/aln_div.num_sites_S
			except TypeError:
				statDict['piS'] = 0
			##making sure that 0/0 = 0 and 0/1= 0 and 0/None = None, in terms of piNpiS.

			try: 
				statDict['piNpiS'] = statDict['piN']/statDict['piS']
			except ZeroDivisionError:
				if statDict['piN'] == 0:
					statDict['piNpiS'] = 0
				else:
					statDict['piNpiS'] = None

			DicOfdicst.update({seq.name: statDict})

	return DicOfdicst

rule popstats_spp:
	""" Calculate a rough dN/dS with the piN/piS code """
	## Modified from https://github.com/topfm/selection-stats/blob/main/piNpiS.py
	input:
		fasta = "alignments/sppcomplex/{gene}_CDS_spp.fa",
		listar = randomfile, # to monitor that it's the latest version
		listah = HICfile,
	output:
		stats = temp("temp/stats/{gene}_stats_spp.txt")
	run:
		# Infer chromosome
		if wildcards.gene in noPagenes_chrs.keys():
			chrom = noPagenes_chrs[wildcards.gene]
		elif 'Pa_' in wildcards.gene:
			chrom = wildcards.gene.split('_')[1]
		else:
			chrom = "NA"

		# Type of gene
		if wildcards.gene in random_selection:
			genetype = "Random"
		elif wildcards.gene in HIC_selection:
			genetype = "NLR_HIC"
		elif wildcards.gene in NLR_selection:
			genetype = "NLR"

		DictSpp = calc_stats_spp(input.fasta, propmissing = 1, skipstop = False)

		# Make the report
		ofile = open(output.stats, 'w')

		for sp in sppstrains_sans_anserina: 
			if sp in DictSpp.keys():
				if DictSpp[sp]['piNpiS'] == None: # Sometimes piS is 0
					piNpiS = "NA"
				else:
					piNpiS = f"{DictSpp[sp]['piNpiS']:.6f}"				

				ofile.write(f"{wildcards.gene}\t{genetype}\t{chrom}\t{sp}\t{DictSpp[sp]['pi']:.6f}\t{DictSpp[sp]['piN']:.6f}\t{DictSpp[sp]['piS']:.6f}\t{piNpiS}\n")
			else:
				ofile.write(f"{wildcards.gene}\t{genetype}\t{chrom}\t{sp}\tNA\tNA\tNA\tNA\n")


rule get_all_popstats_spp:
	input:
		reports = expand("temp/stats/{gene}_stats_spp.txt", gene = focalgenes),
		listar = randomfile, # to monitor that it's the latest version
		listah = HICfile, 
	output:
		"reports/Orthologs_stats_spp.txt"
	shell:
		"printf 'Ortholog\tType\tChr\tSpecies\tD_spp\tdN\tdS\tdNdS\n' > {output}; "
		"cat {input.reports} >> {output}"
	# nseff: Average number of effectively analyzed sequences
	# lseff: Number of effectively analyzed sites

rule get_fasta2axt:
	""" Get script from my GitHub """
	output:
		"scripts/fasta2axt.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/refs/heads/master/Miscellaneous/fasta2axt.py"

rule fasta2axt:
	input:
		script = "scripts/fasta2axt.py",
		fasta = "alignments/sppcomplex/{gene}_CDS_spp.fa",
	output:
		axt = "kakscalculator/axts/{gene}_CDS_spp.axt"
	shell:
		"python {input.script} {input.fasta} --focus anserina > {output.axt}"

rule kakscalculator2:
	input:
		axt = "kakscalculator/axts/{gene}_CDS_spp.axt"
	output:
		kaks = "kakscalculator/runs/{gene}_CDS_spp.axt.kaks"
	shell:
		"KaKs_Calculator -i {input.axt} -o {output.kaks} -m MA"


rule gather_kaks:
	input:
		expand("kakscalculator/runs/{gene}_CDS_spp.axt.kaks", gene = focalgenes)
	output:
		txt = "reports/Orthologs_kakscalculator.txt"
	run:
		with open(output.txt, 'w') as ofile:
			# Start with the general header
			ofile.write("Ortholog\tSpecies\tMethod\tKa\tKs\tKaKs\tp_val\tLength\tS_Sites\tN_Sites\tFold_Sites\tSubstitutions\tS_Substitutions\tN_Substitutions\tFold_S_Substitutions\tFold_N_Substitutions\tDivergence_Time\tSubstitution_Rate_Ratio\tGC\tML_Score\tAICc\tAkaike_Weight\tModel\n")

			for ortho in focalgenes:
				with open(f"kakscalculator/runs/{ortho}_CDS_spp.axt.kaks") as kaks:
					for line in kaks:
						if 'Sequence\tMethod' not in line: #avoid the header
							tabs = line.rstrip("\n").split("\t")
							siblingsp = tabs[0].replace("anserina_vs_", "")
							newline = ortho + '\t' + siblingsp + '\t' + "\t".join(tabs[1:]) + "\n"
							ofile.write(newline)

### ---- Is there a mutational bias in the variable sites? ----

transtrans_dic = {
	"AA>CA": "Transversion",
	"AC>CC": "Transversion",
	"AG>CG": "Transversion",
	"AT>CT": "Transversion",

	"AA>TA": "Transversion",
	"AC>TC": "Transversion",
	"AG>TG": "Transversion",
	"AT>TT": "Transversion",

	"GA>CA": "Transversion",
	"GC>CC": "Transversion",
	"GG>CG": "Transversion",
	"GT>CT": "Transversion",

	"GA>TA": "Transversion",
	"GC>TC": "Transversion",
	"GG>TG": "Transversion",
	"GT>TT": "Transversion",

	"CA>AA": "Transversion",
	"CC>AC": "Transversion",
	"CG>AG": "Transversion",
	"CT>AT": "Transversion",

	"CA>GA": "Transversion",
	"CC>GC": "Transversion",
	"CG>GG": "Transversion",
	"CT>GT": "Transversion",

	"TA>AA": "Transversion",
	"TC>AC": "Transversion",
	"TG>AG": "Transversion",
	"TT>AT": "Transversion",

	"TA>GA": "Transversion",
	"TC>GC": "Transversion",
	"TG>GG": "Transversion",
	"TT>GT": "Transversion",

	"AA>GA": "Transition",
	"AC>GC": "Transition",
	"AG>GG": "Transition",
	"AT>GT": "Transition",

	"GA>AA": "Transition",
	"GC>AC": "Transition",
	"GG>AG": "Transition",
	"GT>AT": "Transition",

	"CA>TA": "Transition",
	"CC>TC": "Transition",
	"CG>TG": "Transition",
	"CT>TT": "Transition",

	"TA>CA": "Transition",
	"TC>CC": "Transition",
	"TG>CG": "Transition",
	"TT>CT": "Transition",
}

rule count_variable_sites:
	input:
		anserina = "alignments/anserina/{gene}_CDS_anserina.fa",
		otherspp = "alignments/sppcomplex/{gene}_CDS_spp.fa"
	output:
		muts = "temp/{gene}_mutations.txt"
	run:
		# What type of gene is this?
		if wildcards.gene in random_selection:
			genetype = 'Random'
		elif wildcards.gene in HIC_selection:
			genetype = 'NLR_HIC'
		elif wildcards.gene in NLR_selection:
			genetype = 'NLR'

		anserina_aln = AlignIO.read(input.anserina, "fasta")
		otherspp_aln_raw = AlignIO.read(input.otherspp, "fasta")
		otherspp_aln = MultipleSeqAlignment([record for record in otherspp_aln_raw if record.id != 'anserina']) # remove anserina

		alignment_length = anserina_aln.get_alignment_length() # both the anserina and the otherspp should be the same length

		with open(output.muts, "w") as ofile: 
			for i in range(alignment_length):
				column_anserina = [record.seq[i].upper() for record in anserina_aln]
				unique_bases_anserina = set(column_anserina)

				column_otherspp = [record.seq[i].upper() for record in otherspp_aln]
				unique_bases_otherspp = set(column_otherspp)

				# First the easiest case to define the ancestral state
				if len(unique_bases_otherspp) == 1: # if there are more than one bases in the other spp it becomes harder to assign the ancestral state 
					ancestral_state = list(unique_bases_otherspp)[0]
					new_mutations = set(unique_bases_anserina) - set(unique_bases_otherspp)

					# is it polymorphic?
					if len(unique_bases_anserina) > 1:
						Status = "Polymorphic"
					elif len(unique_bases_anserina) == 1:
						Status = "Fixed"

					if len(new_mutations) > 0: 
						# Extract context from outgroups (ancestral)
						if i == 0: 
							context_left = []
						else:
							context_left = list(set([record.seq[i - 1].upper() for record in otherspp_aln]))

						if i == alignment_length - 1: 	
							context_right = []
						else:
							context_right = list(set([record.seq[i + 1].upper() for record in otherspp_aln]))

						if len(context_left) == 1 and len(context_right) == 1: # It could happen that it's polymorphic too ... but the easy case first: monomorphic
							left_base = context_left[0]
							right_base = context_right[0]

							for mutation in new_mutations:
								if mutation != '-' and ancestral_state != '-' and right_base != '-':
									nucpairs = ancestral_state + '>' + mutation

									context = f'{left_base}[{ancestral_state}>{mutation}]{right_base}'

									ContextMutation = f"{ancestral_state}{right_base}>{mutation}{right_base}"
									transtrans = transtrans_dic[ContextMutation]

									line = f'{wildcards.gene}\t{genetype}\t{i + 1}\t{mutation}\t{nucpairs}\t{left_base}\t{right_base}\t{context}\t{mutation}p{right_base}\t{ContextMutation}\t{transtrans}\t{Status}\n'
									ofile.write(line)
						else:
							# Skip cases with ambiguous context (could also handle with gaps if desired)
							continue

rule merge_variable_sites:
	input:
		expand("temp/{gene}_mutations.txt", gene = focalgenes)
	output:
		muts = "reports/Orthologs_mutationprofile.txt"
	shell:
		"printf 'Ortholog\tType\tPosition\tNew_base\tMutation\tLeft_base\tRight_base\tContext\tDinucleotide\tContextMutation\tMutType\tStatus\n' > {output}; "
		"cat {input} >> {output}"

### ---------- Putting it all together -----------

rule Plot_PiDistNLRvsRandom:
	input:
		dists = "reports/Orthologs_stats_anserina.txt", # Population genetics statistics
		phylo = "reports/Orthologs_Distances.txt", # Phylogenetic distances between genes
		dNdS = "reports/Orthologs_stats_spp.txt", # Raw dN/dS ratios based on counts (Nei and Gojobori 1986)
		TEcov = "reports/FlanksTEcontent.tab", # proportion covered by transposons in 2kb flanks
		geneevo = geneevo, # proportion covered by transposons in 2kb flanks
		kaksdists = "reports/Orthologs_kakscalculator.txt", # dN/dS as calculated by kakscalculator
		muts = "reports/Orthologs_mutationprofile.txt" # Mutations based on gene alignments
	output:
		main = "results/Fig4_PopPhyloStats.png",
		RIP = "results/FigS1_RIPmuts.png",
		divstats = "results/FigS2_SupDivStats.png",
		dNdS = "results/FigS3_dNdSdegraded.png",
		# piSpiN = "results/FigS_piSpiN.png",
	conda: 
		"envs/plot.yaml"
	script:
		PiDistNLRvsRandom


### ---------- Quantifying trans-species polymorphism/incomplete lineage sorting -----------

rule clean_alignment:
	""" My baseline alignments contain a reference sequence so I need to clean it """
	input:
		fasta = get_input_fasta_cds,
		listar = randomfile, # to monitor that it's the latest version
		listah = HICfile,
	output:
		fasta = "alignments/all/{gene}_CDS_all.fa"
	run:
		# Load alignment from a fasta file
		alignment = AlignIO.read(input.fasta, "fasta")

		# All strains
		allstrains = sum(straindic.values(), []) # sum(iterable, start) adds all elements in iterable to start. straindic.values() gives a list of lists so I need to merge them

		# Get just the relevant sequences
		relevant_seqs = []
		for strain in allstrains:
			for seq in alignment:
				if strain in seq.id: 
					relevant_seqs.append(seq.id)
					break

		# Write a new alignment fasta file with just Podospora anserina strains
		new_alignment = MultipleSeqAlignment([seq for seq in alignment if seq.id in relevant_seqs])

		with open(output.fasta, "w") as handle: 
			SeqIO.write(new_alignment, handle, "fasta")


rule IQTreePerGeneAll:
	""" Run IQTree for each gene """
	input:
		aln = "alignments/all/{gene}_CDS_all.fa"
	output:
		tree = "iqtree_all/{gene}/{gene}_CDS_all.treefile"
	threads: 2
	params:
		bootstraps = 1000, # UFBoot
	run:
		# Check that there are at least 3 species in there
		records_dict = SeqIO.to_dict(SeqIO.parse(input.aln, "fasta"))

		if len(records_dict) > 3: # you can't make a tree with two taxa...
			cmd = f"iqtree -s {input.aln} -m MFP -seed 1234 -bb {params.bootstraps} -nt {threads} -bnni -pre 'iqtree_all/{wildcards.gene}/{wildcards.gene}_CDS_all' --keep-ident -redo"
			shell(cmd)
		elif len(records_dict) > 2: # no bootstraps 
			cmd = f"iqtree -s {input.aln} -m MFP -seed 1234 -nt {threads} -bnni -pre 'iqtree_all/{wildcards.gene}/{wildcards.gene}_CDS_all' --keep-ident -redo"
			shell(cmd)
		else:
			ofile = open(output.tree, "w") # make a dummy empty file


rule renametrees:
	""" Rename the leaves in the tree to match the strains """
	input:
		tree = "iqtree_all/{gene}/{gene}_CDS_all.treefile"
	output:
		tree = "iqtree_all/{gene}/{gene}_CDS_all_nom.tre"
	run:
		tree = Tree(input.tree, format = 0)

		# All strains
		allstrains = sum(straindic.values(), []) # sum(iterable, start) adds all elements in iterable to start. straindic.values() gives a list of lists so I need to merge them

		# Loop through the leaves and change the names
		for node in tree.iter_leaves():
			for sample in allstrains:
				if sample in node.name:
					node.name = sample
					break

		# Print tree in a new file
		namedtree = open(output.tree, 'w')
		namedtree.write(tree.write() + "\n")

rule monophyly_check:
	""" Check if the species are monophyletic in the tree """
	input:
		tree = "iqtree_all/{gene}/{gene}_CDS_all_nom.tre"
	output:
		check = temp("temp/monophyly/{gene}_check.txt")
	run:
		tree = Tree(input.tree)

		# Only check monophyly for species with more than one strain
		relevantspp = [species for species in straindic.keys() if len(straindic[species]) > 1]

		# Start the output file
		ofile = open(output.check, "w")

		# Set an arbitrary outgroup
		rooted = False
		for loner in ["CBS112042p", "CBS411.78m", "CBS415.72m"]:
			print(f"Attempting to root with {loner}")
			if tree.search_nodes(name=loner):
				tree.set_outgroup(loner)
				print("...Success")
				rooted = True
				break
			else:
				print("...Fail") # The loner species are not in the tree

		if not rooted: # We can try to root with P. anserina itself, which should always be present
			print(f"Attempting to root with P. anserina")
			anserinastrains = straindic["anserina"]
			allothers = []

			for node in tree.iter_leaves():
				if node.name not in anserinastrains: 
					allothers.append(node.name)

			if len(allothers) > 0: # there is something else apart from P. anserina in the tree
				ancestor = tree.get_common_ancestor(allothers) # and by exclusion, the anserinas might become monophyletic
				tree.set_outgroup(ancestor) # set is root

				# Check if P. anserina really is monophyletic
				is_monophyletic, clade_type, node = tree.check_monophyly(anserinastrains, target_attr="name", ignore_missing = True)

				if is_monophyletic:
					rooted = True
					print("...Success")
				else:
					print("...Fail") # P. anserina is not monophyletic in this tree
			# If there is only anserina, leave the tree as is because it will be "monophyletic" (the IQ-Tree trees are rooted newicks)

		Monophyletic = True
		for species in relevantspp:
			# Make a list of strains for this species
			strains = [] 
			sppstrains = straindic[species] # all possible names, including redundant (e.g., PaSp and Podan2)
			
			# Check if there are any missing strains
			missing = [name for name in sppstrains if not tree.search_nodes(name=name)]
			# print("Missing nodes:", species, missing) # debugging

			if sppstrains == missing: # The gene is absent in this species!
				status = "Absent"
			elif len(sppstrains) - 1 == len(missing): # There is a single strain from the species there
				status = "NA"
			else:

				for node in tree.iter_leaves():
					if node.name in sppstrains:
						strains.append(node.name)

				# Check if this species is monophyletic
				is_monophyletic, clade_type, node = tree.check_monophyly(strains, target_attr="name")

				if is_monophyletic:
					status = "Monophyletic"
				else:
					status = "Not_monophyletic"
					Monophyletic = False

			oline = f"{wildcards.gene}\t{species}\t{status}\n"
			ofile.write(oline)

		if Monophyletic:
			ofile.write(f"{wildcards.gene}\tAll\tMonophyletic\n")
		else:
			ofile.write(f"{wildcards.gene}\tAll\tNot_monophyletic\n")

rule get_monophyly_checks:
	input:
		expand("temp/monophyly/{gene}_check.txt", gene = focalgenes)
	output:
		"reports/Monophyly_check.txt"
	shell:
		"printf 'Ortholog\tSpecies\tCheck\n' > {output}; "
		"cat {input} >> {output}"

# http://daler.github.io/gffutils/database-ids.html
id_spec={"gene": ["ID", "Name"], 
	"mRNA": ["ID", "transcript_id"], 
	"repeat": ["repeat", "simple_repeat", "similarity"]
	} 

## Check Tajima's D calculated from windows in the NEE paper
rule get_Tajima_windows:
	""" Find genes in the given BED file """
	input:
		gff = ANNOTATION,
		bed = FULLTajima,
	output:
		txt = "reports/FocalGenes_WageningenTajimaMin5Kb.txt"
	run:
		tabs = [line.rstrip("\n").split("\t") for line in open(input.bed) if '#' not in line]

		db = gffutils.create_db(data = input.gff, 
			dbfn = ':memory:',
			force = True, # force=True overwrite any existing databases.
			id_spec = id_spec, 
			merge_strategy = "create_unique")

		## Get iterator of focal genes
		genes = []
		for gene in db.features_of_type("gene"):
			for thing in gene.attributes:
				if gene.attributes[thing][0] in focalgenes: # loop through the attributes, to see if you find the Name
					gene.name = thing
					genes.append(gene)
					continue # move to the next gene

		# Make a dictionary of the contig ranges in the BED file
		beddic = {}
		for line in open(input.bed):
			if 'Contig' not in line:
				contig, start, end, tajima = line.rstrip("\n").split("\t")
				if contig not in beddic.keys():
					beddic[contig] = [[start, end, tajima]]
				else:
					beddic[contig].append([start, end, tajima])
			
		ofile = open(output.txt, 'w')

		# Find the windows overlapping with the genes
		for gene in set(genes):
			genename = gene['Name'][0]
			genectg = gene.seqid
			genestart = int(gene.start)
			geneend = int(gene.end)

			geneTajimas = []
			# print(genename, genectg, genestart, geneend)
			for ranges in beddic[genectg]:
				if genestart >= float(ranges[0]) and geneend <= float(ranges[1]):
					geneTajimas.append(float(ranges[2]))

			if len(geneTajimas) == 0:
				ofile.write(f"{genename}\tNA\n")
			else:
				# ofile.write('\t'.join([str(num) for num in geneTajimas]) + '\n') # for debugging
				ofile.write(f"{genename}\t{sum(geneTajimas)/len(geneTajimas)}\n") # This dilutes the signal too much

rule PlotBalancingSelection:
	input:
		checks = "reports/Monophyly_check.txt",
		dists = "reports/Orthologs_stats_anserina.txt",
		dNdS = "reports/Orthologs_stats_spp.txt",
		winTajima = "reports/FocalGenes_WageningenTajimaMin5Kb.txt",
		geneevo = geneevo
	output:
		plot = "results/Fig5_BalSelection.png",
		altajima = "figures/FigS_TajimaAlignDistribution.png",
		wintajima = "figures/FigS_TajimaWinDistribution.png",
	conda: 
		"envs/plot.yaml"
	script:
		PlotBalancingSelection



# ### ---------- Original code, abandoned ----------

# def calculate_pairwise_differences(alignment):
# 	num_sequences = len(alignment)
# 	pairwise_identities = []

# 	# Generate all pairs of sequences
# 	for seq1, seq2 in itertools.combinations(alignment, 2):
# 		identical_positions = 0
# 		total_positions = 0

# 		# Convert both sequences to uppercase
# 		seq1_upper = seq1.seq.upper()
# 		seq2_upper = seq2.seq.upper()

# 		# Compare the two sequences
# 		for base1, base2 in zip(seq1_upper, seq2_upper):
# 			# Count positions that are not gaps in both sequences
# 			if (base1 != '-' and base2 != '-') or (base1 != 'N' or base2 != 'N'):
# 				total_positions += 1
# 				if base1 == base2:
# 					identical_positions += 1
		
# 		if total_positions > 0:
# 			identity = identical_positions / total_positions
# 			# print(seq1.id, seq2.id, identical_positions, total_positions, identity, len(seq1_upper)) # Debugging
# 			pairwise_identities.append(identity)

# 	# Calculate average pairwise identity
# 	if pairwise_identities:
# 		avg_identity = sum(pairwise_identities) / len(pairwise_identities)
# 	else:
# 		avg_identity = 0

# 	return 1 - avg_identity

# def get_input_fasta(wildcards):
# 	if wildcards.gene in random_selection:
# 		return f"data/random/{wildcards.gene}.fa"
# 	elif wildcards.gene in HIC_selection:
# 		return f"data/nlrsHIC/{wildcards.gene}_noHICrepts.fa"
# 	elif wildcards.gene in NLR_selection:
# 		return f"data/nlrs/{wildcards.gene}.fa"


# rule anserina_pi:
# 	""" Calculate average pairwise-identity within P. anserina """
# 	input:
# 		fasta = get_input_fasta,
# 		listar = randomfile, # to monitor that it's the latest version
# 		listah = HICfile, 
# 	output:
# 		txt = "temp/{gene}_anserina_pi.txt"
# 		# txt = temp("temp/{gene}_anserina_pi.txt")
# 	run:
# 		# Load alignment from a fasta file
# 		alignment = AlignIO.read(input.fasta, "fasta")

# 		# Get just the anserina sequences
# 		anserina_seqs = []
# 		for strain in anserina_strains:
# 			for seq in alignment:
# 				if strain in seq.id: 
# 					anserina_seqs.append(seq.id)
# 					break

# 		# Infer chromosome
# 		if wildcards.gene in noPagenes_chrs.keys():
# 			chrom = noPagenes_chrs[wildcards.gene]
# 		elif 'Pa_' in wildcards.gene:
# 			chrom = wildcards.gene.split('_')[1]
# 		else:
# 			chrom = "NA"

# 		new_alignment = MultipleSeqAlignment([seq for seq in alignment if seq.id in anserina_seqs])

# 		# Calculate the average pairwise identity
# 		average_differences = calculate_pairwise_differences(new_alignment)

# 		# Calculate the average pairwise identity and average GC content
# 		average_gc_content = GC_content_alignment(new_alignment)

# 		# Type of gene
# 		if wildcards.gene in random_selection:
# 			genetype = "Random"
# 		elif wildcards.gene in HIC_selection:
# 			genetype = "NLR_HIC"
# 		elif wildcards.gene in NLR_selection:
# 			genetype = "NLR"

# 		# Make the report
# 		ofile = open(output.txt, 'w')
# 		ofile.write(f"{wildcards.gene}\t{genetype}\t{average_differences:.6f}\t{average_gc_content:.6f}\t{len(new_alignment)}\t{len(alignment)}\t{len(alignment[0])}\t{chrom}\n")



