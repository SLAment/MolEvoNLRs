# -*- snakemake -*-

### 3SppIntrogression: A pipeline to compare pair-wise identity of P. anserina/pauciseta/comata
#############################################################################
# The aim of the pipeline is to find introgression (or ancestral polymorphism)
# tracks. By calling SNPs against the reference and then plotting the identity
# patterns, tracks of identity should (hopefully) pop-up.

### Sources:
# https://software.broadinstitute.org/gatk/documentation/quickstart.php
# First steps: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
# Germline short variant discovery (SNPs + Indels): https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
# HaplotypeCaller: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
# GenomicsDBImport: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php
# GenotypeGVCFs: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php
# GATK hard filtering: https://software.broadinstitute.org/gatk/documentation/article?id=23216#2
# VCF INFO annotations: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_annotator_ReadPosRankSumTest.php
# http://corearray.sourceforge.net/tutorials/SNPRelate/


#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# 2024/04/03
# ---------------

import gzip
import csv

# -------------------------------------------------
## Variables set from the configuration file
configfile: "config/config.yaml"
# -------------------------------------------------
samples_file = config["SampleIDs_file"]

# VCF file before filtering for missing data
VCFfile = config["VCFfile"]
PathToCovBED = config["PathToCovBED"]
TEsites = config["TEsites"]

# Output name
OUTPUTbasename = config["OUTPUTbasename"]

# Variables
WINSIZE = config["WINSIZE"]

# Scripts
PopGenome_WinMissingness = config["PopGenome_WinMissingness"]
PopGenome_PairwiseDivergence = config["PopGenome_PairwiseDivergence"]
PaDivergenceCorr2kb = config["PaDivergenceCorr2kb"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: removeoverlap_all_spp, cat_beds_all, cat_beds, reduce_bed, bedOfIntrogression, removeoverlap, extractchrcov, get_badsites2vcf, compressvcf, cat_PopGenome_PairwiseDivergence, compress_PairwiseDivergence
# ----------

WINSIZEkb = int(round(WINSIZE/1000))
COMPAREDSPECIES = ["pauciseta", "comata"]

SppStrainDic = {}
for line in open(samples_file, 'r'):
	tab = line.rstrip("\n").split('\t')
	if tab[1] in SppStrainDic.keys():
		SppStrainDic[tab[1]].append(tab[0])
	else:
		SppStrainDic[tab[1]] = [tab[0]]

## -------------------------

rule all:
	input:
		"figures/FigS6_DxyPauciComata_chrs.png",
		expand("results/Dxy_{WINSIZEkb}kb_{tipo}-all.bed", tipo = ["introgression", "missingness"], WINSIZEkb = WINSIZEkb )

### --- Produced a new bed file of missing data ---

rule BEDOPS:
	""" Take all the BED files and overlap them into one single BED file"""
	# https://bedops.readthedocs.io/en/latest/content/overview.html#overview
	# https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#merge-m-merge
	input:
		TEsites,
		expand(PathToCovBED + "/{sample}_cov_filtered.bed", sample = SppStrainDic["anserina"] + SppStrainDic["pauciseta"] + SppStrainDic["comata"])
	output:
		temp("filter/excludedsites-raw.bed")
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "1:00:00",
	shell:
		"bedops --merge {input} > {output}"

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

rule removeoverlap:
	""" Fuse overlapping ranges in a BED file """
	# Basically the same as the script totalcovergff.py but for bed files
	input:
		"filter/excludedsites-raw.bed"
	output:
		"filter/excludedsites.bed"
	run:
		# Make a dictionary of ranges
		beddic = {} # key: chromosome, value: (start, end)

		with open(input[0], 'r') as file: 
			for line in file:
				if '#' in line:
					pass
				elif line not in ['\n', '\r\n']: # Ignore empty lines
					cols = line.rstrip("\n").split("\t")

					contig = cols[0]
					start = int(cols[1])
					end = int(cols[2])

					if contig in list(beddic.keys()):
						beddic[contig].append([start, end])
					else: # contig is new
						beddic[contig] = [[start, end]]

			# Reduce the overlaps
			for ctg in beddic.keys():
				beddic[ctg] = remove_overlap(beddic[ctg])

			# Print them in a new file
			with open(output[0], 'w') as result:
				result.write('#Contig\tStart\tEnd\n') # header
				for ctg in beddic.keys():
					for interval in beddic[ctg]:
						result.write('{}\t{}\t{}\n'.format(ctg, interval[0], interval[1]))

rule extractchrcov:
	""" Extract one individual chromosome from the bed file """
	input:
		"filter/excludedsites.bed"
	output:
		"filter/excludedsites-chr{nchr}.bed"
	shell:
		"""
		head -n1 {input} > {output}
		grep 'chromosome_{wildcards.nchr}' {input} >> {output}
		""" 

### ----

rule subset_strains:
	input:
		vcf = VCFfile
	output:
		vcf = f"vcfs/{OUTPUTbasename}-3Spp.vcf.gz"
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	run:
		samplelist = ""
		for sample in SppStrainDic["anserina"] + SppStrainDic["pauciseta"] + SppStrainDic["comata"]:
			samplelist += f"{sample},"
		# for sample in otherSpp:
		# 	samplelist += f"^{sample}," # Didn't work :/

		bcfcmd = f"bcftools view -Oz -s '{samplelist.rstrip(',')}' {input.vcf} > {output.vcf}"
		shell(bcfcmd)

# -O, --output-type b|u|z|v[0-9]
# Output compressed BCF (b), uncompressed BCF (u), compressed VCF
# (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools
# subcommands to speed up performance by removing unnecessary
# compression/decompression and VCF←→BCF conversion.   The compression level
# of the compressed formats (b and z) can be set by by appending a number
# between 0-9.

rule snpsvcfnomiss:
	""" Filter out sites with missing data """
	input:
		vcf = f"vcfs/{OUTPUTbasename}-3Spp.vcf.gz"
	output:
		vcf = f"vcfs/{OUTPUTbasename}-3Spp-miss1.vcf.gz"
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	shell:
		"vcftools --gzvcf {input.vcf} --max-missing 1 --recode --recode-INFO-all --stdout | bgzip > {output.vcf}; "
		"tabix -p vcf {output.vcf}"

rule get_badsites2vcf:
	output:
		temp("scripts/badsites2vcf.py")
	shell:
		"wget -O {output} https://raw.githubusercontent.com/johannessonlab/HetVPaper/master/DiversityStats/scripts/badsites2vcf.py"

rule inputmissingsites:
	""" Input the sites that should be excluded as explicit missing data sites """
	input:
		vcf = f"vcfs/{OUTPUTbasename}-3Spp-miss1.vcf.gz",
		bed = "filter/excludedsites.bed",
		badsites2vcf = "scripts/badsites2vcf.py"
	output:
		f"vcfs/{OUTPUTbasename}-3Spp-miss1-withNA.vcf"
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	shell:
		"python {input.badsites2vcf} {input.vcf} {input.bed} --bed > {output}"

rule compressvcf:
	""" Compress and index vcf file """
	input:
		f"vcfs/{OUTPUTbasename}-3Spp-miss1-withNA.vcf"
	output:
		f"vcfs/{OUTPUTbasename}-3Spp-miss1-withNA.vcf.gz"
	shell:
		"bgzip {input}"

# ------- Prepare the files for PopGenome -------

rule extractchrvcf:
	""" Extract one individual chromosome from the vcf file """
	input:
		vcf = f"vcfs/{OUTPUTbasename}-3Spp-miss1-withNA.vcf.gz",
	output:
		"popgenome/chr{nchr}/" + OUTPUTbasename + "-3Spp-miss1-withNA-chr{nchr}.vcf"
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "15:00",
	shell:
		"vcftools --gzvcf {input.vcf} --recode --recode-INFO-all --chr chromosome_{wildcards.nchr} --stdout > {output}" # Notice it has to be uncompressed for PopGenome

rule PopGenome_WinMissingness:
	""" This is a demanding step that should only be calculated once per chromosome """
	input:
		vcf = "popgenome/chr{nchr}/" + OUTPUTbasename + "-3Spp-miss1-withNA-chr{nchr}.vcf",
		bed = "filter/excludedsites-chr{nchr}.bed",
	output:
		df = "filter/winlens_chr{nchr}_{WINSIZEkb}kb.txt"
	threads: 3 # some chrs need more memory
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "2-00:00:00",
	params:
		chr = "chr{nchr}",
		width = WINSIZE,
		jump = WINSIZE, # Non-overlapping windows
	conda: 
		"envs/popgenome.yaml"
	script:
		PopGenome_WinMissingness

# ------- Calculate divergence -------

rule PopGenome_PairwiseDivergence:
	""" This is a demanding step that should only be calculated once per chromosome """
	input:
		vcf = "popgenome/chr{nchr}/" + OUTPUTbasename + "-3Spp-miss1-withNA-chr{nchr}.vcf",
		winlens = "filter/winlens_chr{nchr}_{WINSIZEkb}kb.txt",
	output:
		df = "reports/Dxy_chr{nchr}_{WINSIZEkb}kb-{strain1}-{strain2}.csv"
	threads: 2 # some chrs need more memory
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	params:
		width = WINSIZE,
		jump = WINSIZE, # Non-overlapping windows
		strain1 = "{strain1}",
		strain2 = "{strain2}",
		chr = "chromosome_{nchr}"
	conda: 
		"envs/popgenome.yaml"
	script:
		PopGenome_PairwiseDivergence


rule cat_PopGenome_PairwiseDivergence:
	input:
		expand("reports/Dxy_chr{{nchr}}_{WINSIZEkb}kb-{strain1}-{strain2}.csv", WINSIZEkb = WINSIZEkb, strain1 = SppStrainDic["anserina"], strain2 = SppStrainDic["pauciseta"] + SppStrainDic["comata"]),
		expand("reports/Dxy_chr{{nchr}}_{WINSIZEkb}kb-{strain1}-{strain2}.csv", WINSIZEkb = WINSIZEkb, strain1 = SppStrainDic["pauciseta"], strain2 = SppStrainDic["comata"])
	output:
		"results/Dxy_chr{nchr}_{WINSIZEkb}kb.csv"
	shell:
		"""
		echo '"Strain1","Strain2","position","winlen","chr","dxy","start","end"' > {output}
		cat {input} >> {output}
		"""

rule compress_PairwiseDivergence:
	input:
		"results/Dxy_chr{nchr}_{WINSIZEkb}kb.csv"
	output:
		"results/Dxy_chr{nchr}_{WINSIZEkb}kb.csv.gz"
	shell:
		"bgzip {input}"

## ---- Plot the correlation between the divergence of P. anserina vs either P. comata or P. pauciseta ----

rule Plot_PaDivergenceCorr2kb:
	input:
		chr1 = f"results/Dxy_chr1_{WINSIZEkb}kb.csv.gz",
		chr2 = f"results/Dxy_chr2_{WINSIZEkb}kb.csv.gz",
		chr3 = f"results/Dxy_chr3_{WINSIZEkb}kb.csv.gz",
		chr4 = f"results/Dxy_chr4_{WINSIZEkb}kb.csv.gz",
		chr5 = f"results/Dxy_chr5_{WINSIZEkb}kb.csv.gz",
		chr6 = f"results/Dxy_chr6_{WINSIZEkb}kb.csv.gz",
		chr7 = f"results/Dxy_chr7_{WINSIZEkb}kb.csv.gz",
		PopData = samples_file,
	output:
		corr = "figures/FigS6_DxyPauciComata_chrs.png",
		genome = "figures/DxyPauciComata_genome.png",
	conda: 
		"envs/simpleplot.yaml"
	script:
		PaDivergenceCorr2kb

# ----
def whatSpecies(strain):
	if strain in SppStrainDic["anserina"]:
		return "anserina"
	elif strain in SppStrainDic["pauciseta"]:
		return "pauciseta"
	elif strain in SppStrainDic["comata"]:
		return "comata"
	else:
		return "other"
# ----

rule bedOfIntrogression:
	""" Find the windows involved with potential introgression """
	input:
		csv = "results/Dxy_chr{nchr}_{WINSIZEkb}kb.csv.gz"
	output:
		rawbed = temp("filter/Dxy_chr{nchr}_{WINSIZEkb}kb_introgression_{species}.bedraw"),
		missing = temp("filter/Dxy_chr{nchr}_{WINSIZEkb}kb_missingness_{species}.bedraw"),
		# bed = temp("filter/Dxy_chr{nchr}_{WINSIZEkb}kb_introgression_{species}.bed"),
	params:
		pi = 0.0004,
		minwin = 500
	run:
		# Start the files
		rawbedfile = open(output.rawbed, 'w')
		missfile = open(output.missing, 'w')

		# Open the bgzip file and read its contents
		with gzip.open(input.csv, mode='rt') as file:
			reader = csv.reader(file)
			for row in reader:
				if 'Strain1' not in row:
					# print(row)
					winlen, chrom, dxy, start, end = row[3:8]
					outline = f"{chrom}\t{int(float(start))}\t{int(float(end))}\t{dxy}\t{whatSpecies(row[1])}\n"

					if dxy != 'NA': # Not total missing data
						if float(dxy) <= params.pi and float(winlen) >= params.minwin:
							rawbedfile.write(outline)
						elif float(winlen) < params.minwin: # too much missing data so basically NA too
							missfile.write(outline)
					else:
						missfile.write(outline)

		rawbedfile.close()
		missfile.close()

rule reduce_bed:
	input:
		rawbed = "filter/Dxy_chr{nchr}_{WINSIZEkb}kb_{tipo}_{species}.bedraw"
	output:
		bed = temp("filter/Dxy_chr{nchr}_{WINSIZEkb}kb_{tipo}_{species}.bed")
	run:
		## Fuse windows like in `removeoverlap`
		# Make a dictionary of ranges
		beddic = {} # key: chromosome, value: (start, end)

		with open(input.rawbed, 'r') as file: 
			for line in file:
				cols = line.rstrip("\n").split("\t")

				contig = cols[0]
				start = int(cols[1]) - 1 # The positions are originally designed to avoid overlap but in this case I want them to overlap, hence the -1
				end = int(cols[2])
				species = cols[4]

				if species == wildcards.species:
					if contig in list(beddic.keys()):
						beddic[contig].append([start, end])
					else: # contig is new
						beddic[contig] = [[start, end]]

			# Reduce the overlaps
			for ctg in beddic.keys():
				beddic[ctg] = remove_overlap(beddic[ctg])

			# Print them in a new file
			with open(output.bed, 'w') as result:
				result.write('#Contig\tStart\tEnd\n') # header
				for ctg in beddic.keys():
					for interval in beddic[ctg]:
						result.write('{}\t{}\t{}\n'.format(ctg, interval[0], interval[1]))

rule cat_beds:
	input:
		expand("filter/Dxy_chr{nchr}_{{WINSIZEkb}}kb_{{tipo}}_{{species}}.bed", nchr = list(range(1, 8)))
	output:
		"results/Dxy_{WINSIZEkb}kb_{tipo}_{species}.bed"
	shell:
		"cat {input} > {output}"

rule cat_beds_all:
	input:
		expand("results/Dxy_{{WINSIZEkb}}kb_{{tipo}}_{species}.bed", species = COMPAREDSPECIES)
	output:
		"filter/Dxy_{WINSIZEkb}kb_{tipo}-raw.bed"
	shell:
		"cat {input} > {output}"

rule removeoverlap_all_spp:
	# This one is redundant right now because by design the missing data is across all species
	""" Fuse overlapping ranges in a BED file """
	# Basically the same as the script totalcovergff.py but for bed files
	input:
		"filter/Dxy_{WINSIZEkb}kb_{tipo}-raw.bed"
	output:
		"results/Dxy_{WINSIZEkb}kb_{tipo}-all.bed"
	run:
		# Make a dictionary of ranges
		beddic = {} # key: chromosome, value: (start, end)

		with open(input[0], 'r') as file: 
			for line in file:
				if '#' in line:
					pass
				elif line not in ['\n', '\r\n']: # Ignore empty lines
					cols = line.rstrip("\n").split("\t")

					contig = cols[0]
					start = int(cols[1])
					end = int(cols[2])

					if contig in list(beddic.keys()):
						beddic[contig].append([start, end])
					else: # contig is new
						beddic[contig] = [[start, end]]

			# Reduce the overlaps
			for ctg in beddic.keys():
				beddic[ctg] = remove_overlap(beddic[ctg])

			# Print them in a new file
			with open(output[0], 'w') as result:
				result.write('#Contig\tStart\tEnd\n') # header
				for ctg in beddic.keys():
					for interval in beddic[ctg]:
						result.write('{}\t{}\t{}\n'.format(ctg, interval[0], interval[1]))

