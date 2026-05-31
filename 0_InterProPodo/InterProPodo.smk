# -*- snakemake -*-

### InterProPodo
#############################################################################
# Running InterProScan on a given proteome
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/12/12
# ++++++++++++++++++++++++++++++++++++++++++++++

# -------------------------------------------------
## Variables set from the configuration file
configfile: "config/config.yaml"
# -------------------------------------------------
# Path to protein-coding genes in GFF3 format (nice-3.00 style)
path2annotations = config["path2annotations"]
# Path to genome assemblies that matches annotation in fasta format
path2genomes = config["path2genomes"]
# Sample IDs
SAMPLES = config["SAMPLES"]
# Full path to interproscan.sh
interproscan = config["interproscan"]
# -------------------------------------------------

import gffutils

# ----------
# Rules not submitted to a job
localrules: get_gffutils2fasta,
			gffutils2fasta,
			get_reference_PFAMs,
			process_pfams_nice,
			describe_pfams_nice,
			merge_to_report,
# ----------


rule all:
	input:
		"reports/GenePFAMCounts.txt"

rule get_gffutils2fasta:
	output:
		"scripts/gffutils2fasta.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/master/GenomeAnnotation/gffutils2fasta.py"


def get_annotation(wildcards):
	if wildcards.sample == "Podan2":
		return f"{path2annotations}/Podan2.nice-3.02.gff3" # The updated annotation 
	elif wildcards.sample == "Podan2016":
		return f"{path2annotations}/genome_annotation_PODANS_v2016.gff" # From here: https://github.com/Podospora-anserina/transcript_annotation_2022
	else:
		return f"{path2annotations}/{wildcards.sample}.nice-3.00.gff3" # As in Ament-Velasquez et al. (2024)

def get_assembly(wildcards):
	if wildcards.sample in ["CBS415.72m", "CBS411.78m", "CBS112042p"]:
		return f"{path2genomes}/{wildcards.sample}.nice.fa" # As in Ament-Velasquez et al. (2024)
		# return f"{path2genomes}/{wildcards.sample}.nice_mt.fa" # in case these are available in the path, the ones submitted in NCBI with a corrected mitochondrion
	elif wildcards.sample in ["Podan2"]:
		return f"{path2genomes}/Podan2_AssemblyScaffoldsmt.fa" # If you got it from JGI
	elif wildcards.sample in ["PODCO"]:
		return f"{path2genomes}/PODCO_genomic.fas" # If you got it from NCBI
	elif wildcards.sample == "Podan2016":
		return f"{path2genomes}/genomePodoMatPlus.fasta" # From here: https://github.com/Podospora-anserina/transcript_annotation_2022
	else:
		return f"{path2genomes}/{wildcards.sample}.nice.fa" # As in Ament-Velasquez et al. (2024)

rule gffutils2fasta:
	input:
		gff = get_annotation,
		fas = get_assembly,
		gffutils2fasta = "scripts/gffutils2fasta.py"
	output:
		fas = "annotation/{sample}.fa"
	shell:
		"python {input.gffutils2fasta} {input.fas} {input.gff} -t CDS --proteinon --join --onlyids --output {output.fas}; "
		"sed -i 's/*//g' {output.fas}" # InterProScan doesn't like the stop codons

rule InterProScan:
	""" Annotate the proteome """
	input:
		"annotation/{sample}.fa"
	output:
		"InterProScan/{sample}/{sample}.tsv"
	threads: 4
	resources:
		cpus_per_task = lambda wildcards, threads: threads,
		runtime = "8h",
	log: "InterProScan/{sample}_IPS.log"
	shell:
		"{interproscan} -i {input} --output-file-base InterProScan/{wildcards.sample}/{wildcards.sample} --cpu {threads} -appl Pfam &> {log}"
		# -appl Pfam 		restrict the analysis to just PFAMs
		# -goterms,--goterms                                        Optional, switch on lookup of corresponding Gene Ontology
		#                                                           annotation (IMPLIES -iprlookup option)
		# -iprlookup,--iprlookup                                    Also include lookup of corresponding InterPro annotation in
		#                                                           the TSV and GFF3 output formats.
		# -pa,--pathways                                            Optional, switch on lookup of corresponding Pathway
		#                                                            annotation (IMPLIES -iprlookup option)

rule get_reference_PFAMs:
	output:
		temp("data/PFAM_reference.txt")
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/MolEvoNLRs/refs/heads/main/8_HICproteins/data/PFAM_reference.txt"

# http://daler.github.io/gffutils/database-ids.html
id_spec={"gene": ["ID", "Name"], "mRNA": ["ID"]} 

# TE-like keywords
TEkeywods = ["ransposase", "DDE superfamily", "Reverse transcriptase", "LTR", "Pol polyprotein", "virus", "RNase H", "ntegrase", "transposable"]
suspishkeywords = ["nuclease"]

rule process_pfams_nice:
	input:
		gff = get_annotation,
		ips = "InterProScan/{sample}/{sample}.tsv",
		pfamref = "data/PFAM_reference.txt"
	output:
		ips = "results/{sample}_pretty.tsv",
		gff = "results/{sample}_colored.gff",
	run:
		# Make a database of the annotation
		db = gffutils.create_db(data = input.gff, 
			dbfn = ':memory:',
			force = True, # force=True overwrite any existing databases.
			id_spec = id_spec, 
			merge_strategy = "create_unique")

		## Get iterator of all genes
		# genes = [gene for gene in db.features_of_type("gene")] # I know there are no genes without CDS because it comes from Augustus
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

			## Specific to "nice-3.00" annotation
			if '_' in geneid:
				geneidcode = geneid.split("_")[1]

				if len(geneidcode) == 6: # A Pa_X_XXXX homolog was identified when I named the genes
					chromosome = geneidcode[0]
					genenum = geneidcode[1:].lstrip('0')
					pacode = f"Pa_{chromosome}_{genenum}"
			else: # e.g., Sp001
				pacode = gene['Name'][0]

			allgenes[geneid] = [pacode]
		
		# Read the PFAM reference
		pfamsref = [line.rstrip("\n").split("\t") for line in open(input.pfamref) if 'Pfam' not in line]

		pfamsrefdic = {}
		for ref in pfamsref:
			pfamsrefdic[ref[2]] = ref[0]

		# Read the InterProScan annotation
		tabs = [line.rstrip("\n").split("\t") for line in open(input.ips) if '#' not in line]
		lentab = len(tabs[0])

		# Now let's look at the annotation
		annotated_tabs = {}
		pfamdic = {}
		TEGenes = []
		SuspishGenes = []
		NLRGenes = []

		HETcount = []

		for tab in tabs:
			geneid = tab[0]
			pfamcode = tab[4]
			description = tab[5]
			NLRbool = "Other" # For plotting

			# In InterProScan 5.75-106.0, the PFAM PF24883 (NPHP3_N, Nephrocystin 3, N-terminal) was 
			# added and took over everything previously called NACHT. So I will re-name it as NACHT for consistency
			if pfamcode == 'PF24883':
				tab[4] = 'PF05729'
				tab[5] = "NACHT domain"

			if pfamcode == 'PF06985':
				HETcount.append(geneid)

			## Or classify the domains as NLR-related
			if pfamcode in pfamsrefdic.keys():
				NLRbool = "NLR-associated"
			elif 'ANK' in description or 'Ankyrin' in description:
				NLRbool = "NLR-associated"
			elif 'HEAT' in description:
				NLRbool = "NLR-associated"
			elif 'Heterokaryon incompatibility' in description:
				NLRbool = "NLR-associated"
			elif 'NACHT_N' in description:
				NLRbool = "NLR-associated"
			elif 'NACHT' in description: 
				NLRbool = "NLR-associated"
			elif 'NB-ARC' in description:
				NLRbool = "NLR-associated"
			elif 'TPR' in description or 'Tetratricopeptide repeat' in description:
				NLRbool = "NLR-associated"
			elif 'PNP_UDP' in description:
				NLRbool = "NLR-associated"
			elif 'Abhydrolase' in description or 'Alpha/beta hydrolase' in description or 'AB_hydrolase' in description:
				NLRbool = "NLR-associated"

			# Make a dictionary of pfams and their descriptions
			if pfamcode not in pfamdic.keys():
				pfamdic[pfamcode] = description

			if NLRbool != "Other": NLRGenes.append(geneid)

			## Is the name suspish?
			repeatpfam = "None"

			# I'm not sure they are TE, but they could be
			for keyword in suspishkeywords:
				if keyword in description:
					repeatpfam = "Suspicious"
					break
			
			# Look very much like a TE (it will rewrite suspish if necessary)
			for keyword in TEkeywods:
				if keyword in description:
					repeatpfam = "TE"
					break

			if repeatpfam == "Suspicious": SuspishGenes.append(geneid)
			if repeatpfam == "TE": TEGenes.append(geneid)
			# print(geneid, allgenes[geneid], NLRbool, repeatpfam)
			newline = '\t'.join([geneid] + allgenes[geneid] + [NLRbool] + [repeatpfam] + tab[1:]) + '\n'

			# There might be more than one PFAM per gene
			if geneid not in annotated_tabs.keys():
				annotated_tabs[geneid] = [newline]
			else:
				annotated_tabs[geneid].append(newline)

		# What genes are annotated but normal otherwise?
		normalgenes = list(set(annotated_tabs.keys()) - set(TEGenes) - set(SuspishGenes) - set(NLRGenes))

		# Report some numbers
		print(f"Number of TE genes: {len(set(TEGenes))}")
		print(f"Number of suspish (TE-like) genes: {len(set(SuspishGenes))}")
		print(f"Number of NLR-like genes: {len(set(NLRGenes))}")
		print(f"\tNumber of HET-domain genes: {len(set(HETcount))}")
		print(f"Other annotated genes: {len(normalgenes)}")

		# Make the output files
		with open(output.ips, 'w') as ofile, open(output.gff, 'w') as ogff: 
			ogff.write("##gff-version 3\n")

			for gene in genes:
				# Write the tsv
				geneid = gene.id
				if geneid in annotated_tabs.keys():
					genedescriptions = [] # Make a list to store the PFAMs of this gene
					for pfamline in annotated_tabs[geneid]:
						ofile.write(pfamline)
						pfam = pfamline.rstrip("\n").split("\t")[7]
						genedescriptions.append(pfamdic[pfam])
				else:
					newline = '\t'.join([geneid] + allgenes[geneid] + ["NA"]*(lentab + 3) ) + '\n'
					ofile.write(newline)

				# Write the gff3
				colour = "#000000" # black
				if geneid in TEGenes: 
					colour = "#C82909"
				elif geneid in SuspishGenes:
					colour = "#FAA18F"
				elif geneid in NLRGenes:
					colour = "#279583"
				elif geneid in normalgenes: # annotated genes that are not suspish or NLRs
					colour = "#688678"

				gene['color'] = colour
				
				# Write the gene body
				if geneid in annotated_tabs.keys(): 
					descriptionline = ','.join(set(genedescriptions))
					gene['description'] = descriptionline
				# Add the Pa_ name if available
				if allgenes[geneid] != ['NA']: 
					gene['Alias'] = allgenes[geneid]

				# Write the gene in the new gff file
				ogff.write(str(gene) + "\n")

				# Write the children
				for child in list(db.children(gene)):
					if child.featuretype == "mRNA":
						child['color'] = colour
						if geneid in annotated_tabs.keys():
							child['description'] = descriptionline

					ogff.write(str(child) + "\n")

rule describe_pfams_nice:
	input:
		# gff = "results/{sample}_colored.gff",
		ips = "results/{sample}_pretty.tsv",
		pfamref = "data/PFAM_reference.txt"
	output:
		report = temp("temp/{sample}_GenePFAMCounts.txt")
	threads: 1
	run:
		# Read the PFAM reference
		pfamsref = [line.rstrip("\n").split("\t") for line in open(input.pfamref) if 'Pfam' not in line]
		pfamsrefdic = {}
		for ref in pfamsref:
			pfamsrefdic[ref[2]] = ref[0]

		one_to_one_obvious_orthologs = 0
		pfamdic = {}
		genedic = {}

		# Read the InterProScan annotation
		with open(input.ips, 'r') as ips:
			for line in ips:
				tabs = line.rstrip("\n").split("\t")
				geneid, Pacode, NLRbool = tabs[0:3]
				pfamcode, description = tabs[7:9]

				# Make a dictionary of pfams and their descriptions
				if pfamcode not in pfamdic.keys():
					pfamdic[pfamcode] = description

				if geneid in genedic.keys():
					genedic[geneid][1].append(pfamcode)
				else:
					genedic[geneid] = [Pacode, [pfamcode]]

		# Interesting metrics to consider
		NACHTcount = 0
		NACHTcountOrtho = 0
		HETcount = 0
		HETcountOrtho = 0
		NBARCcount = 0
		NBARCcountOrtho = 0

		for gene in genedic.keys():
			Pacode = genedic[gene][0]
			if 'Pa_' in Pacode:
				one_to_one_obvious_orthologs += 1

			geneNLRdomains = []
			for pfam in set(genedic[gene][1]):
				if pfam in pfamsrefdic.keys():
					geneNLRdomains.append(pfamsrefdic[pfam])

			if 'NACHT' in geneNLRdomains: 
				NACHTcount+= 1
				if 'Pa_' in Pacode: NACHTcountOrtho +=1
			if 'HET' in geneNLRdomains: 
				HETcount+= 1
				if 'Pa_' in Pacode: HETcountOrtho +=1
			if 'NB-ARC' in geneNLRdomains: 
				NBARCcount+= 1
				if 'Pa_' in Pacode: NBARCcountOrtho +=1

		with open(output.report, 'w') as ofile:
			newline = f'{wildcards.sample}\t{len(genedic.keys())}\t{one_to_one_obvious_orthologs}\t{HETcount}\t{HETcountOrtho}\t{NACHTcount}\t{NACHTcountOrtho}\t{NBARCcount}\t{NBARCcountOrtho}\n'
			ofile.write(newline)


rule merge_to_report:
	input:
		expand("temp/{sample}_GenePFAMCounts.txt", sample = SAMPLES)
	output:
		"reports/GenePFAMCounts.txt"
	shell:
		"printf 'Strain\tNo_genes\tObvious_Orthologs\tHET\tHET_ortholog\tNACHT\tNACHT_ortholog\tNB-ARC\tNB-ARC_ortholog\n' > {output}; "
		"cat {input} >> {output}"

