# -*- snakemake -*-

### smartHIC: Quantifying HIC repeats across the Tree of Life
#############################################################################
# Proteins with WD40, ANK, and TPR repeatd were extracted from the SMART
# database. The idea is to look for instances of HIC with T-reks and report
# some numbers.
#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# 2025-06-09
# +++++++++++++++++++++++++++++++++++++++++++++++++

import re
import time
from Bio import Entrez

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"
profilefile: "profile/config.v8+.yaml"

# Input parameters
TAXA = config["TAXA"]
path2files = config["path2files"]
PFAMref = config["PFAMref"]

# REPEATs = ["WD40"]
REPEATs = config["repeat_types"]

PSIMmin = config["PSIMmin"]

# Scripts
HICinSMART = config["HICinSMART"]

# -------------------------------------------------


# ----------
# Rules not submitted to a job
localrules: merge_tables, rename_sequences, looking_for_NLRs, merge_annotations, filterTreks, plot

# ----------


rule all:
	input:
		"results/Barplot_HICproteins_propotions.png",
		"results/AllvsHICproteins_propotions.png",
		"results/AllvsHICproteins_totals.png"

taxodic = {
	"Actinobacteria": "Bacteria",
	"Archaea": "Archaea",
	"Ascomycetes": "Fungi",
	"Bacteroidota": "Bacteria",
	"Basidiomycetes": "Fungi",
	"Chloroflexi": "Bacteria", 
	"Cyanobacteria": "Bacteria",
	"Metazoa": "Animals",
	"OtherEukaryotes": "Other Eukaryotes",
	"OtherFungi": "Fungi",
	"Planctomycetes": "Bacteria",
	"Proteobacteria": "Bacteria",
	"Verrucomicrobia": "Bacteria",
	"Viridiplantae": "Viridiplantae"
}

rule filterTreks:
	""" Parse the T-reks output into a managable table """
	input:
		treks = path2files + "{repeat}/Treks_{taxon}.txt",
		fasta = path2files + "{repeat}/diamond_deepclust_{taxon}.fa",
	output:
		treks = "processing/{repeat}/Filtered_Treks_{taxon}.txt"
	params:
		minlen = 21,
		maxlen = 100
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		# mem_mb = "16GB",
		time = "3:00:00",
	run:
		sppname = re.compile("([a-zA-Z0-9:\\-_()\\[\\],\\. ]+)\\[([a-zA-Z0-9\\.\\-\\/:=_#()'+ ]+)\\]")

		# Compare to the fasta file
		allseqs = {}
		with open(input.fasta, 'r') as fasta:
			for line in fasta:
				if '>' in line:
					try:
						db, seqID, fulldescription = line.rstrip("\n").split("|")
					except:
						tabs = line.rstrip("\n").split("|")
						db = tabs[0]
						seqID = tabs[1]
						fulldescription = '_'.join(tabs[2:])

					matchy = sppname.search(fulldescription)
					if matchy:
						annotation = matchy.group(1)
						taxonomy = matchy.group(2)

					allseqs[seqID] = (taxonomy, annotation)

		# Filter the T-reks output
		treksdic = {}
		with open(input.treks, 'r') as treks:
			for line in treks: 
				if '>' in line: # Assuming Uniprot database
					try:
						db, seqID, fulldescription = line.rstrip("\n").split("|")
					except:
						tabs = line.rstrip("\n").split("|")
						db = tabs[0]
						seqID = tabs[1]
						fulldescription = '_'.join(tabs[2:])

					## The fasta files are the reduced set, so ignore this sequence if it's not in the reduced set
					if seqID in allseqs.keys():
						# Get the description from the fasta because T-reks tends to truncate long headers
						taxonomy, annotation = allseqs[seqID]
					else:
						taxonomy = []
						annotation = []
					
					treksdic[seqID] = (taxonomy, annotation, [])
					# treksdic[seqID] = ([], [], [])

				elif 'Length:' in line: #the repeats start
					tabs = line.rstrip("\n").split(" ")
					lenrepeat = int(tabs[1])
					nb = int(tabs[5])
					psim = tabs[13].lstrip('Psim:')
					replen = int(tabs[15].lstrip('Length:'))

					if lenrepeat >= params.minlen and lenrepeat <= params.maxlen:
						treksdic[seqID][2].append([lenrepeat, nb, psim, replen])

		# Make dictionary with sequences NOT in the T-reks file but present in the fasta file
		noHICseqs = {}
		for seqID in allseqs.keys():
			if seqID not in treksdic.keys():
				noHICseqs[seqID] = allseqs[seqID]

		# differences = set(allseqs.keys()) ^ set(noHICseqs.keys())
		# # print(differences)
		# if len(differences) > 0:
		# 	print(f"{len(differences)} missing seqs in fasta file for repeat {wildcards.repeat} and taxon {wildcards.taxon}")
		# 	sys.exit(1)

		# Print the table
		with open(output.treks, 'w') as ofile: 
			for seqID in treksdic.keys(): # Print the prots with HIC repeats
				if seqID in allseqs.keys(): # if not it didn't make it to the set with reduced redundancy
					taxonomy, annotation, predictions = treksdic[seqID]
					maxregionlen = 0
					maxindex = 0

					genus = taxonomy.lstrip('uncultured ').split(" ")[0].lstrip("'")
					if genus == 'fungal': genus = 'NA'

					phylum = "NA"

					if len(predictions) != 0:
						for repindex in range(0, len(predictions)):
							if predictions[repindex][3] > maxregionlen:
								maxregionlen = predictions[repindex][3]
								maxindex = repindex
						
						winner = predictions[maxindex]
						line = f'{wildcards.repeat}\t{taxodic[wildcards.taxon]}\t{wildcards.taxon}\t{phylum}\t{seqID}\t{genus}\t{taxonomy}\t{winner[0]}\t{winner[1]}\t{winner[2]}\t{winner[3]}\n'
						ofile.write(line)

			for seqID in noHICseqs.keys(): # Print the rest
				taxonomy, annotation = noHICseqs[seqID]
				genus = taxonomy.lstrip('uncultured ').split(" ")[0].lstrip("'")
				if genus == 'fungal': genus = 'NA'

				phylum = "NA"

				line = f'{wildcards.repeat}\t{taxodic[wildcards.taxon]}\t{wildcards.taxon}\t{phylum}\t{seqID}\t{genus}\t{taxonomy}\t0\t0\tNA\tNA\n'
				ofile.write(line)			

rule merge_tables:
	input:
		expand("processing/{repeat}/Filtered_Treks_{taxon}.txt", taxon = TAXA, repeat = REPEATs)
	output:
		"reports/HICinSMART.txt"
	shell:
		"cat {input} | sed 's/#/no/' > {output}"


# ---- Domain annotation ----

rule rename_sequences:
	""" Rename the sequences to have only the sequence ID """
	input:
		fasta = path2files + "{repeat}/diamond_deepclust_{taxon}.fa",
	output:
		fasta = "temp/{repeat}/diamond_deepclust_{taxon}_renamed.fa"
	run:
		with open(input.fasta, 'r') as fasta, open(output.fasta, 'w') as ofile:
			for line in fasta:
				if '>' in line:
					tabs = line.split("|")
					seqID = tabs[1]
					ofile.write(f'>{seqID}\n')
				else:
					ofile.write(line)

rule InterProScan:
	""" Annotate the whole genome """
	input:
		fasta = "temp/{repeat}/diamond_deepclust_{taxon}_renamed.fa"
	output:
		"InterProScan/{repeat}/diamond_deepclust_{taxon}_IPS.tsv"
	threads: 4
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "16GB",
		time = "08:00:00",
	log: "InterProScan/{repeat}/diamond_deepclust_{taxon}_IPS.log"
	shell:
		"interproscan.sh -i {input} --output-file-base InterProScan/{wildcards.repeat}/diamond_deepclust_{wildcards.taxon}_IPS --cpu {threads} -appl Pfam > {log}"

# -appl Pfam 		restrict the analysis to just PFAMs
# -goterms,--goterms										Optional, switch on lookup of corresponding Gene Ontology
#														   annotation (IMPLIES -iprlookup option)
# -iprlookup,--iprlookup								    Also include lookup of corresponding InterPro annotation in
#														   the TSV and GFF3 output formats.
# -pa,--pathways										    Optional, switch on lookup of corresponding Pathway
#														    annotation (IMPLIES -iprlookup option)

rule looking_for_NLRs:
	input:
		ips = "InterProScan/{repeat}/diamond_deepclust_{taxon}_IPS.tsv",
		pfamref = PFAMref
	output:
		tsv = "processing/InterProScan/{repeat}_{taxon}_annotation.tsv"
	run:
		# Read the PFAM reference
		pfamsref = [line.rstrip("\n").split("\t") for line in open(input.pfamref) if 'Pfam' not in line]

		ATPasedomains = ['NACHT', 'NB-ARC', 'ATPase']
		pfamsrefdic = {}
		atpaserefdic = {}
		pfamsreptsdic = {}
		for ref in pfamsref:
			if ref[0] in ATPasedomains: 
				atpaserefdic[ref[2]] = ref[0]
			else: # non-atpase to try to find Nterm domains later
				pfamsrefdic[ref[2]] = ref[0]

			if ref[0] in REPEATs:
				if ref[0] != wildcards.repeat: # To check if there are more repeat types other than the expected one
					pfamsreptsdic[ref[2]] = ref[0]

		# Read the InterProScan annotation
		tabs = [line.rstrip("\n").split("\t") for line in open(input.ips) if '#' not in line]

		annotated_genes = {}
		for tab in tabs:
			seqID = tab[0]
			pfamcode = tab[4]
			description = tab[5]
			interprocode = tab[10]
			interprodescription = tab[11]

			ATPase = 'None'
			if pfamcode in atpaserefdic.keys():
				ATPase = atpaserefdic[pfamcode]

			# Look for N-term domains
			nlrdomain = 'None'
			for pfam in pfamsrefdic.keys():
				if pfam == pfamcode: 
					nlrdomain = pfamsrefdic[pfam]
					break

			# Look for repeats out of place
			ExtraRepeats = 'None'
			if pfamcode in pfamsreptsdic.keys():
				ExtraRepeats = pfamsreptsdic[pfamcode]

			# Is it an NLR sensu lato?
			if seqID in annotated_genes.keys():
				annotated_genes[seqID][0].append(pfamcode)
				annotated_genes[seqID][1].append(ATPase)
				annotated_genes[seqID][2].append(nlrdomain)
				annotated_genes[seqID][3].append(description)
				annotated_genes[seqID][4].append(ExtraRepeats)
			else:
				annotated_genes[seqID] = ([pfamcode], [ATPase], [nlrdomain], [description], [ExtraRepeats])

		with open(output.tsv, 'w') as ofile:
			for gene in annotated_genes.keys():

				pfams = list(set(annotated_genes[gene][0]))
				ATPases = list(set(annotated_genes[gene][1]))
				nlrdomains = list(set(annotated_genes[gene][2]))
				descriptions = list(set(annotated_genes[gene][3]))
				ExtraRepeats = list(set(annotated_genes[gene][4]))

				# Make string of pfams
				pfam_str = ','.join(pfams)

				# Deal with the None clases
				ATPases_str = ','.join(ATPases)
				ATPases_str = ATPases_str.replace('None,', '').replace(',None', '')

				# Some cases are ambiguous
				if ATPases_str == 'NB-ARC,NACHT' or ATPases_str == 'NACHT,NB-ARC' or ATPases_str == 'ATPase,NACHT' or ATPases_str == 'NACHT,ATPase':
					ATPases_str = 'ATPase'

				nlrdomains_str = ','.join(nlrdomains)
				nlrdomains_str = nlrdomains_str.replace('None,', '').replace(',None', '')

				descriptions_str = ','.join(descriptions)

				# Deal with the None clases
				ExtraRepeats_str = ','.join(ExtraRepeats)
				ExtraRepeats_str = ExtraRepeats_str.replace('None,', '').replace(',None', '')

				# Write a simple file
				ofile.write(f'{gene}\t{wildcards.taxon}\t{wildcards.repeat}\t{ATPases_str}\t{nlrdomains_str}\t{descriptions_str}\t{pfam_str}\t{ExtraRepeats_str}\n')

rule merge_annotations:
	input:
		expand("processing/InterProScan/{repeat}_{taxon}_annotation.tsv", taxon = TAXA, repeat = REPEATs)
	output:
		"reports/SMART_IPSannotation.txt"
	shell:
		"printf 'SeqID\tGroup\tRepeat\tATPase\tNLRdomains\tDescriptions\tPfams\tExtraRepeats\n' > {output}; "
		"cat {input} >> {output}"


rule plot:
	""" Plot the abundances of HIC with and without ATPases """
	input:
		annotation  = "reports/SMART_IPSannotation.txt",
		treks = "reports/HICinSMART.txt",
	output:
		hicprop = "results/Barplot_HICproteins_propotions.png",
		allvshic = "results/AllvsHICproteins_propotions.png",
		totals = "results/AllvsHICproteins_totals.png"
	params:
		psimmin = PSIMmin
	conda:
		"envs/plot.yaml"
	script:
		HICinSMART



