# -*- snakemake -*-

### IntroRegions
#############################################################################
# Characterizing introgression regions between Podospora anserina and P. comata or P. pauciseta
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2024/06/27
# ++++++++++++++++++++++++++++++++++++++++++++++

# -------------------------------------------------
## Variables set from the configuration file
configfile: "config/config.yaml"
profilefile: "/shared/projects/podosporagenomics/pipelines/Goby/profile/config.v8+.yaml"
# -------------------------------------------------
# Input files
INTROBED_both = config["INTROBED_both"]
INTROBED_pau = config["INTROBED_pau"]
INTROBED_com = config["INTROBED_com"]
INTROBED_miss = config["INTROBED_miss"]
INTROBED_taj = config["INTROBED_taj"]

INTROBED_missp = config["INTROBED_miss_pau"]
INTROBED_missc = config["INTROBED_miss_com"]

FULLTajima = config["FULLTajima"]

ANNOTATION = config["ANNOTATION"]
REFGENOME = config["REFGENOME"]
PFAMref = config["PFAMref"]

GenesNLRsCoords = config["GenesNLRsCoords"]

# Scripts
EnrichOddRatios = config["EnrichOddRatios"]
EnrichOddRatiosPval = config["EnrichOddRatiosPval"]
NLRchrDist = config["NLRchrDist"]
# -------------------------------------------------

import gffutils
from collections import Counter
from scipy.stats import fisher_exact
import pandas as pd
from statsmodels.stats.multitest import multipletests

# ----------
# Rules not submitted to a job
localrules: fix_missingness, 
			plot_NLRchrdist, 
			bed_of_domains, 
			plot_enrichment_alla, 
			enrichmnent_by_domain, 
			get_data, 
			get_data_taj, 
			get_gff_introgressed, 
			get_gffutils2fasta, 
			get_protseqs, 
			get_protseqs_genome, 
			process_pfams, 
			get_Tajima_NWDgenes
# ----------

rule all:
	input:
		"results/NLRdomains_chrs.png",
		"results/NLRgenes_chrs.pdf",
		"regions/NWDrelatedTajima.txt",
		"results/enrichmnent_allplots.pdf",
		expand("results/enrichmnent_OddsVsPval_{type}.pdf", type = ["both", "miss"])

rule get_data:
	""" This is just for convenience with the rule expansion """
	input:
		both = INTROBED_both,
		pau = INTROBED_pau,
		com = INTROBED_com,
		miss = INTROBED_miss,
	output:
		both = "baseline/targetregions_both.bed",
		pau = "baseline/targetregions_pau.bed",
		com = "baseline/targetregions_com.bed",
		miss = "baseline/targetregions_miss_raw.bed",
	shell:
		"ln -sf ../{input.both} {output.both}; " 
		"ln -sf ../{input.pau} {output.pau}; " 
		"ln -sf ../{input.com} {output.com}; " 
		"ln -sf ../{input.miss} {output.miss}; " 

# The naming is missleading, it's not about introgression or 
# Dxy but I just want it to go through the pipeline (it has nothing to do with the target regions)
rule get_data_taj:
	""" This is just for convenience with the rule expansion """
	input:
		bed = INTROBED_taj, # 10 kb windows with 1 kb steps (from which at least 5kb have data) and Tajima's D >= 2
	output:
		bed = "baseline/targetregions_taj.bed",
	shell:
		"ln -sf ../{input.bed} {output.bed}; " # it didn't like making a symlink to the same folder...

rule true_missing:
	""" Find the windows that are missing in both pauci and comata """
	# Windows that are not missing in one of the species implies that the other is missing because it's most likely not introgressed (otherwise it wouldn't be missing)
	input:
		pau = INTROBED_missp,
		com = INTROBED_missc
	output:
		"baseline/missedInBothSpp.bed"
	shell:
		"bedtools intersect -a {input.pau} -b {input.com} > {output}"

rule fix_missingness:
	input:
		# miss = "baseline/targetregions_miss_raw.bed",
		miss = "baseline/missedInBothSpp.bed",
		intro = "baseline/targetregions_both.bed"
	output:
		"baseline/targetregions_miss.bed"
	shell:
		"bedtools subtract -a {input.miss} -b {input.intro} > {output}"

	# The bed file represent a combination of 10kb windows that were missing
	# in either comata or pauciseta compared to anserina. However, the
	# introgression windows also represent a combination. Hence, only those
	# windows that were not detected as introgressed are truly missing from
	# the dataset.

# ----
def remove_overlap(ranges):
	""" Simplify a list of ranges; I got it from https://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap """
	result = []
	current_start = -1
	current_stop = -1 

	for start, stop in sorted(ranges):
		if start > current_stop:
			# this segment starts after the last segment stops
			# just add a new segment
			result.append( (start, stop) )
			current_start, current_stop = start, stop
		else:
			# current_start already guaranteed to be lower
			current_stop = max(current_stop, stop)
			# segments overlap, replace
			result[-1] = (current_start, current_stop) # SLAV: I modified this to update the stop too.
	return(result)
# ----

# http://daler.github.io/gffutils/database-ids.html
id_spec={"gene": ["ID", "Name"], 
	"mRNA": ["ID", "transcript_id"], 
	"repeat": ["repeat", "simple_repeat", "similarity"]
	} 


rule get_gff_introgressed:
	""" Find genes in the given BED file """
	input:
		gff = ANNOTATION,
		bed = "baseline/targetregions_{clase}.bed",
	output:
		gff = "regions/Podan2_targetregions_{clase}.gff",
		genes = "regions/Podan2_targetregions_{clase}_genesIDs.txt"
	run:
		# Here I had to make decisions. If a gene is partially outside the
		# missing data then it means that the rest of the gene can still be
		# classified as either introgressed, under a Tajima peak or normal. I chose to
		# give priority to that little bit of remaining information outside
		# of the missing-data window, rather than call it all missing data.

		tabs = [line.rstrip("\n").split("\t") for line in open(input.bed) if '#' not in line]

		# Read the whole annotation
		db = gffutils.create_db(data = input.gff, 
			dbfn = ':memory:',
			force = True, # force=True overwrite any existing databases.
			id_spec = id_spec, 
			merge_strategy = "create_unique")

		## Get iterator of all genes
		genes = [gene for gene in db.features_of_type("gene")]

		targetregions_genes =[]

		for gene in genes:
			geneID = gene['ID'][0]
			feature_contig = gene.seqid
				
			for tab in tabs: # Loop through the regions and see if the gene falls there
				contig, bed_start, bed_end  = tab[:]

				if feature_contig == contig: #the right contig	
					# We want from the start to the stop codon, including the introns, but excluding the UTRs
					allchildren = [child for child in db.children(gene, featuretype='CDS', order_by='start')]

					if len(allchildren) >= 1: # Some genes might not have CDS (e.g. tRNAs)
						feature_start = allchildren[0].start # Start of gene excluding UTRs
						feature_end = allchildren[len(allchildren) - 1].end # End of gene excluding UTRs

						if wildcards.clase != "miss": # Other categories can touch partially a gene and keep it
							if (feature_end < int(bed_start)) or (feature_start > int(bed_end)): # not overlapping features
								pass
							else: # there is some overlap!
								targetregions_genes.append(geneID)
								continue # Go to the next gene
						elif (feature_start > int(bed_start)) and (feature_end < int(bed_end)): # missing data only keeps genes totally inside
							targetregions_genes.append(geneID)
							continue # Go to the next gene

		targetregions_genes = list(set(targetregions_genes))

		# Make a list of genes
		with open(output.genes, 'w') as geneslist:
			for ID in targetregions_genes:
				geneslist.write(ID + '\n')

		# Start an output gff file
		with open(output.gff, 'w') as ofile:
			ofile.write("##gff-version 3\n")

			for gene in genes:
				if gene.id in targetregions_genes:
					ofile.write(str(gene) + "\n")
					for child in list(db.children(gene)):
						ofile.write(str(child) + "\n")
					ofile.write("\n")

rule get_gffutils2fasta:
	output:
		"scripts/gffutils2fasta.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/master/GenomeAnnotation/gffutils2fasta.py"

rule get_protseqs:
	input:
		gff = "regions/Podan2_targetregions_{clase}.gff",
		fas = REFGENOME,
		gffutils2fasta = "scripts/gffutils2fasta.py"
	output:
		fas = "regions/Podan2_targetregions_{clase}.fas"
	shell:
		"python {input.gffutils2fasta} {input.fas} {input.gff} -t CDS --proteinon --join --onlyids --output {output.fas}"

rule get_protseqs_genome:
	input:
		gff = ANNOTATION,
		fas = REFGENOME,
		gffutils2fasta = "scripts/gffutils2fasta.py"
	output:
		fas = "annotation/Full_proteome.fas"
	shell:
		"python {input.gffutils2fasta} {input.fas} {input.gff} -t CDS --proteinon --join --onlyids --output {output.fas}; "
		"sed -i 's/*//g' {output.fas}"

rule InterProScan: # Modify for the IFB Core Cluster (Montpellier)
	""" Annotate the proteome """
	input:
		"annotation/Full_proteome.fas"
	output:
		"InterProScan/Full_proteome.tsv"
	threads: 4
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "16GB",
		time = "08:00:00",
	log: "InterProScan/Full_proteome_IPS.log"
	shell:
		# "module load interproscan/5.75-106.0; "
		"interproscan.sh -i {input} --output-file-base InterProScan/Full_proteome --cpu {threads} -appl Pfam > {log}"

# -appl Pfam 		restrict the analysis to just PFAMs
# -goterms,--goterms                                        Optional, switch on lookup of corresponding Gene Ontology
#                                                           annotation (IMPLIES -iprlookup option)
# -iprlookup,--iprlookup                                    Also include lookup of corresponding InterPro annotation in
#                                                           the TSV and GFF3 output formats.
# -pa,--pathways                                            Optional, switch on lookup of corresponding Pathway
#                                                            annotation (IMPLIES -iprlookup option)

rule process_pfams:
	input:
		gff = ANNOTATION,
		ips = "InterProScan/Full_proteome.tsv",
		pfamref = PFAMref
	output:
		ips = "InterProScan/Full_proteome_pretty.tsv"
	run:
		# Make a database to retrive the Pa_X_XXXX codes 
		db = gffutils.create_db(data = input.gff, 
			dbfn = ':memory:',
			force = True, # force=True overwrite any existing databases.
			id_spec = id_spec, 
			merge_strategy = "create_unique")

		## Get iterator of all genes
		# genes = [gene for gene in db.features_of_type("gene")]
		genes = []
		for gene in db.features_of_type("gene"):
			allchildren = [child for child in db.children(gene, featuretype='CDS', order_by='start')]
			if len(allchildren) >= 1: # Make sure you get only genes with CDS (protein-coding)
				genes.append(gene)


		# Try to find the Pa_X_XXXX codes for each gene (gene names of P. anserina reference)
		allgenes = {}
		for gene in genes:
			geneid = gene.id
			pacode = 'NA'
			alias = 'NA'

			for attribute in ['Note', 'Name', 'Alias']:
				if attribute in gene.attributes:
					for item in gene.attributes[attribute]:
						if 'Pa_' in item:
							pacode = item
							break
						elif attribute == 'Name':
							if item != geneid: alias = item
						elif attribute == 'Alias' and 'jgi' not in item:
							if item != geneid: alias = item

			allgenes[geneid] = [pacode, alias]
		
		# Read the PFAM reference
		pfamsref = [line.rstrip("\n").split("\t") for line in open(input.pfamref) if 'Pfam' not in line]

		pfamsrefdic = {}
		for ref in pfamsref:
			pfamsrefdic[ref[2]] = ref[0]

		# Read the InterProScan annotation
		tabs = [line.rstrip("\n").split("\t") for line in open(input.ips) if '#' not in line]

		# Now let's look at the annotation
		annotated_tabs = {}
		for tab in tabs:
			geneid = tab[0]
			pfamcode = tab[4]
			description = tab[5]
			NLRbool = "Other" # For plotting

			## Keep it agnostic
			clase = pfamcode

			# In InterProScan 5.75-106.0, the PFAM PF24883 (NPHP3_N, Nephrocystin 3, N-terminal) was 
			# added and took over everything previously called NACHT. So I will re-name it as NACHT for consistency
			if pfamcode == 'PF24883':
				tab[4] = 'PF05729'
				tab[5] = "NACHT domain"

			## Or classify the domains as NLR-related
			if pfamcode in pfamsrefdic.keys():
				NLRbool = "NLR-associated"
			# 	clase = pfamsrefdic[pfamcode]
			elif 'ANK' in description or 'Ankyrin' in description:
				NLRbool = "NLR-associated"
			# 	clase = "ANK"
			elif 'HEAT' in description:
				NLRbool = "NLR-associated"
			# 	clase = "HEAT"
			elif 'Heterokaryon incompatibility' in description:
				NLRbool = "NLR-associated"
			# 	clase = "HET"
			elif 'NACHT_N' in description:
				NLRbool = "NLR-associated"
			# 	clase = "NAD"
			elif 'NACHT' in description: 
				NLRbool = "NLR-associated"
			# 	clase = "NACHT"
			elif 'NB-ARC' in description:
				NLRbool = "NLR-associated"
			# 	clase = "NB-ARC"
			elif 'TPR' in description or 'Tetratricopeptide repeat' in description:
				NLRbool = "NLR-associated"
				# 	clase = "TPR"
			elif 'PNP_UDP' in description:
				NLRbool = "NLR-associated"
			# 	clase = "PNP_UDP"
			elif 'Abhydrolase' in description or 'Alpha/beta hydrolase' in description or 'AB_hydrolase' in description:
				NLRbool = "NLR-associated"
			# 	clase = "SesB"
			# elif 'WD40' in description or 'WD domain' in description: # There are so many other proteins with WD40 that are not NLRs
			# 	NLRbool = "NLR-associated"
			# 	clase = "WD40"
			# elif 'ATPase' in description or 'NTPase' in description:
			# 	if 'AAA' in description or 'P-' in description:
			# 		clase = "NLR"
			# 		clase = "ATPase"

			# # In case I want to use the NLR as a category rather than the domains
			# if NLRbool == "NLR-associated":
			# 	clase = "NLR"

			newline = '\t'.join([geneid] + allgenes[geneid] + [clase] + [NLRbool] + tab[1:]) + '\n'

			# There might be more than one PFAM per gene
			if geneid not in annotated_tabs.keys():
				annotated_tabs[geneid] = [newline]
			else:
				annotated_tabs[geneid].append(newline)

		ofile = open(output.ips, 'w')
		for gene in genes:
			geneid = gene.id
			if geneid in annotated_tabs.keys():
				for pfam in annotated_tabs[geneid]:
					ofile.write(pfam)
			else:
				newline = '\t'.join([geneid] + allgenes[geneid] + ["NA"]*(len(tabs[0]) + 1) ) + '\n'
				ofile.write(newline)


rule enrichmnent_by_domain:
	input:
		ips = "InterProScan/Full_proteome_pretty.tsv",
		subset = "regions/Podan2_targetregions_{clase}_genesIDs.txt",
		missing = "regions/Podan2_targetregions_miss_genesIDs.txt",
	output:
		table = "results/enrichmnent_{clase}.txt"
	run:
		genes_miss = [line.rstrip("\n") for line in open(input.missing)] 
		
		# If it's in the introgressed region, then genes *partially*
		# overlapping with missing data are still counted as introgressed.
		# The missing regions are already defined to be mutually exclusive with the introgressed regions. 
		# But that's not the case for the Tajima regions
		if input.subset == input.missing:
			genes_subset = [line.rstrip("\n") for line in open(input.subset)]
			tabs = [line.rstrip("\n").split("\t") for line in open(input.ips)]
		else:
			genes_subset = [line.rstrip("\n") for line in open(input.subset) if line.rstrip("\n").split("\t")[0] not in genes_miss]
			tabs = [line.rstrip("\n").split("\t") for line in open(input.ips) if line.rstrip("\n").split("\t")[0] not in genes_miss] # I do have to remove all the genes that fall into missing data for the whole proteome

		genedic = {} # All genes in the proteome
		pfam_meanings = {} # What is the description of each PFAM code?
		pfam_NLRbool = {} # PFAM classified as NLR-associated or Other
		for tab in tabs:
			geneid = tab[0]
			pfam = tab[3]

			if geneid not in genedic.keys():
				genedic[geneid] = [pfam]
			else:
				genedic[geneid].append(pfam)

			if pfam not in pfam_meanings.keys():
				pfam_meanings[pfam] = tab[9]
				pfam_NLRbool[pfam] = tab[4]

		# ----
		quegenes = list(genedic.keys())
		print(len(quegenes))
		print(len(set(quegenes)))
		# ----

		# Get PFAM codes for subset 
		pfam_subset = []
		for gene in genes_subset:
			pfams = [pfam for pfam in set(genedic[gene]) if pfam != 'NA'] # because I'm using set() that means that a given PFAM can only occur once in a gene, and it truly represents how many genes have that PFAM
			pfam_subset.extend(pfams)

		# Get PFAM codes for the rest of the proteome
		genes_rest = set(genedic.keys()) - set(genes_subset)
		pfam_rest = []
		for gene in genes_rest:
			pfams = [pfam for pfam in set(genedic[gene]) if pfam != 'NA']
			pfam_rest.extend(pfams)		

		# Count frequencies of PFAMs
		pfam_subset_counts = Counter(pfam_subset)
		pfam_rest_counts = Counter(pfam_rest)
		
		# Create a DataFrame to hold results
		results = []

		# Total counts of PFAMs
		total_subset = len(pfam_subset)
		total_rest = len(pfam_rest)

		# Total number of genes compared
		total_genes_subset = len(genes_subset)
		total_genes_rest = len(genes_rest)

		# Perform Fisher's exact test for each PFAM code
		pfam_codes = set(pfam_subset_counts.keys()).union(set(pfam_rest_counts.keys()))

		for pfam in pfam_codes:
			subset_count = pfam_subset_counts.get(pfam, 0) # with the get() method, if the key is not found in the dictionary, it returns a default value 0
			rest_count = pfam_rest_counts.get(pfam, 0)

			## Version using the number of genes as a base. Maybe less biased? because some genes will have no annotation
			contingency_table = [
				[subset_count, total_genes_subset - subset_count],
				[rest_count, total_genes_rest - rest_count]
			]

			odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater') # greater: the odds ratio of the underlying population is greater than one
	
			# Get a description of what the PFAM means
			if pfam not in pfam_meanings.keys():
				descrip = pfam
				NLRbool = 'Other'
			else:
				descrip = pfam_meanings[pfam]
				NLRbool = pfam_NLRbool[pfam]

			results.append((pfam, descrip, NLRbool, subset_count, total_subset, total_genes_subset, rest_count, total_rest, total_genes_rest, "{:.4f}".format(odds_ratio), p_value)) # Gene version

		# Convert results to DataFrame
		results_df = pd.DataFrame(results, columns=['PFAM', 'Description', 'Type', 'Subset_Count', "total_PFAMs_subset", "total_genes_subset", 'Rest_Count', "total_PFAMs_rest", "total_genes_rest", "Odds_ratio", 'P_Value'])

		# Adjust p-values for multiple testing
		results_df['P_Value_Adj'] = multipletests(results_df['P_Value'], method='fdr_bh')[1] # The fdr_bh method specifically refers to the Benjamini-Hochberg procedure

		# Sort them by the p-value
		results_df_sorted = results_df.sort_values(by=['P_Value_Adj', 'P_Value'])

		# Write the DataFrame to a tab-separated file
		results_df_sorted.to_csv(output.table, sep='\t', index=False)

		# Filter for significant PFAM codes
		significant_pfam = results_df_sorted[results_df_sorted['P_Value_Adj'] < 0.05]

		# Print results to screen
		print(output.table)
		print(significant_pfam)


rule plot_enrichment_alla:
	""" Plot the top enrichment categories based on the Fisher test """
	# Notice that I excluded domains that have less than 3 occurances in the whole genome
	input:
		both = "results/enrichmnent_both.txt",
		pau = "results/enrichmnent_pau.txt",
		com = "results/enrichmnent_com.txt",
		miss = "results/enrichmnent_miss.txt",
		taj = "results/enrichmnent_taj.txt",
	output:
		plot = "results/enrichmnent_allplots.pdf"
	params:
		top = 10 
	conda: 
		"envs/plot_scales.yaml"
	script:
		EnrichOddRatios

# A plotting alternative
rule plot_EnrichOddRatiosPval:
	input:
		bed = "results/enrichmnent_{type}.txt",
	output:
		plot = "results/enrichmnent_OddsVsPval_{type}.pdf"
	conda: 
		"envs/plot_scales.yaml"
	script:
		EnrichOddRatiosPval


rule bed_of_domains:
	""" Make a bed file for the introgressed domains """
	input:
		ips = "InterProScan/Full_proteome_pretty.tsv",
		gff = ANNOTATION,
		pfamref = PFAMref
	output:
		bed = "results/Genes_NLRdomains_coords.bed"
	run:
		# domains = [line.rstrip("\n").split("\t") for line in open(input.domains)]
		tabs = [line.rstrip("\n").split("\t") for line in open(input.ips)]
		gff = [line.rstrip("\n").split("\t") for line in open(input.gff) if '\tgene\t' in line]

		# # Get domains associated with NLRs as defined above
		# NLRdomains = [dom[0] for dom in domains if dom[2] == "NLR-associated"]
		# Read the PFAM reference
		pfamsref = [line.rstrip("\n").split("\t") for line in open(input.pfamref) if 'Pfam' not in line]

		pfamsrefdic = {}
		for ref in pfamsref:
			pfamsrefdic[ref[2]] = ref[0]

		domaindic = {}
		for tab in tabs:
			geneid = tab[0]
			genepfam = tab[3]
			genedomaindescription = tab[9]

			if genepfam in pfamsrefdic.keys():
				if geneid not in domaindic:
					domaindic[geneid] = [genepfam]
				else:
					domaindic[geneid].append(genepfam)

		ofile = open(output.bed, 'w')

		for gene in gff:
			contig = gene[0]
			startgene = gene[3]
			endgene = gene[4]
			geneid = gene[8].split(";")[0].lstrip("ID=")

			if geneid in domaindic.keys():
				for pfam in set(domaindic[geneid]):
					newline = f"{contig}\t{startgene}\t{endgene}\t{pfamsrefdic[pfam]}\t{pfam}\t{geneid}\n"
					ofile.write(newline)
		## We want
		# chromosome	startgene	endgene	Description	PFAM	geneid

## -------
# Plot chromosomal distribution of NLRs
## -------


rule plot_NLRchrdist:
	input:
		intro = "baseline/targetregions_both.bed",
		miss = "baseline/targetregions_miss.bed",
		genes = "results/Genes_NLRdomains_coords.bed",
		coords = GenesNLRsCoords
	output:
		domain = "results/NLRdomains_chrs.png",
		NLRs = "results/NLRgenes_chrs.pdf",
	conda: 
		"envs/plot_scales.yaml"
	script:
		NLRchrDist

## -------
# Get a report of Tajima for the NWD-related genes
## -------

rule get_Tajima_NWDgenes:
	""" Find genes in the given BED file """
	input:
		gff = ANNOTATION,
		bed = FULLTajima,
	output:
		txt = "regions/NWDrelatedTajima.txt"
	run:
		tabs = [line.rstrip("\n").split("\t") for line in open(input.bed) if '#' not in line]

		db = gffutils.create_db(data = input.gff, 
			dbfn = ':memory:',
			force = True, # force=True overwrite any existing databases.
			id_spec = id_spec, 
			merge_strategy = "create_unique")

		## Get iterator of NLR genes
		genes = [gene for gene in db.features_of_type("gene") if 'NLR' in gene.id]

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
		for gene in genes:
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
				# ofile.write(f"{genename}\t{max(geneTajimas)}\n")
				# ofile.write('\t'.join([str(num) for num in geneTajimas]) + '\n') # for debugging
				ofile.write(f"{genename}\t{sum(geneTajimas)/len(geneTajimas)}\n") # This dilutes the signal too much


