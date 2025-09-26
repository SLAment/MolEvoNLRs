# -*- snakemake -*-

### 2SppIntrogression: A pipeline to compare pair-wise identity of P. anserina/pauciseta/comata
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

## Species being compared to P. anserina
COMPAREDSPECIES = config["comparedspecies"]

# VCF file before filtering for missing data
VCFfile = config["VCFfile"]
PathToCovBED = config["PathToCovBED"]
TEsites = config["TEsites"]

# Name of the output vcf files after filtering
OUTPUTbasename = config["OUTPUTbasename"]

# Path to the diversity statistics calculated in Ament-Velásquez et al. (2022) NEE
path2divstats = config["path2divstats"]

# Coordinates of NLRgenes
NLRmetadata = config["NLRmetadata"]

# Variables
WINSIZE = config["WINSIZE"]
JUMP = config["JUMP"]
MINWINSIZE = config["MINWINSIZE"]
PI = config["PI"]

# Scripts
PopGenome_WinMissingness = config["PopGenome_WinMissingness"]
PopGenome_PairwiseDivergence = config["PopGenome_PairwiseDivergence"]
PaNLRchrPopStats2Spp = config["PaNLRchrPopStats2Spp"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: reduce2nonoverlap, removeoverlap_all_spp, cat_beds_all, cat_beds, reduce_bed, bedOfIntrogression, removeoverlap, extractchrcov, get_badsites2vcf, compressvcf, cat_PopGenome_PairwiseDivergence, compress_PairwiseDivergence
# ----------

WINSIZEkb = int(round(WINSIZE/1000))
CHRNUMS = list(range(1, 8))

SppStrainDic = {}
for line in open(samples_file, 'r'):
	tab = line.rstrip("\n").split('\t')
	if tab[1] in SppStrainDic.keys():
		SppStrainDic[tab[1]].append(tab[0])
	else:
		SppStrainDic[tab[1]] = [tab[0]]

## -------------------------

# I need this to avoid ambiguities with names
wildcard_constraints:
    species ="\\w+" # equivalent to {sample,\w+} - limits to alphabet letters

rule all:
	input:
		"figures/Fig6_Stats_NLRs.png",
		"figures/FigS7_Stats_NLRs_nwd5.png",
		"figures/FigS8_Stats_NLRs_hetb.png",
		"figures/FigS9_Stats_NLR_hetz.png",

		# Input for IntroRegions.smk
		expand("results/Dxy_{WINSIZEkb}kb_{tipo}-all.bed", tipo = ["introgression", "missingness"], WINSIZEkb = WINSIZEkb )
		
		# expand("results/Dxy_chr{nchr}_{WINSIZEkb}kb_nonoverlap.csv.gz", nchr = CHRNUMS, WINSIZEkb = WINSIZEkb),


### --- Produced a new bed file of missing data ---

def get_samples_bed(wildcards):
	files = [TEsites]
	for sample in SppStrainDic["anserina"] + SppStrainDic[wildcards.species]:
		newfile = f"{PathToCovBED}/{sample}_cov_filtered.bed"
		files.append(newfile)
	return(files)

rule BEDOPS:
	""" Take all the BED files and overlap them into one single BED file"""
	# https://bedops.readthedocs.io/en/latest/content/overview.html#overview
	# https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#merge-m-merge
	input:
		get_samples_bed
	output:
		temp("filter/excludedsites-raw_{species}.bed")
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
		"filter/excludedsites-raw_{species}.bed"
	output:
		"filter/excludedsites_{species}.bed"
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
		"filter/excludedsites_{species}.bed"
	output:
		"filter/excludedsites_{species}-chr{nchr}.bed"
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
		vcf = "vcfs/{OUTPUTbasename}_{species}.vcf.gz"
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	run:
		samplelist = ""
		for sample in SppStrainDic["anserina"] + SppStrainDic[wildcards.species]:
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
		vcf = "vcfs/{OUTPUTbasename}_{species}.vcf.gz"
	output:
		vcf = "vcfs/{OUTPUTbasename}_{species}-miss1.vcf.gz"
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
		vcf = "vcfs/{OUTPUTbasename}_{species}-miss1.vcf.gz",
		bed = "filter/excludedsites_{species}.bed",
		badsites2vcf = "scripts/badsites2vcf.py"
	output:
		"vcfs/{OUTPUTbasename}_{species}-miss1-withNA.vcf"
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
		"vcfs/{OUTPUTbasename}_{species}-miss1-withNA.vcf"
	output:
		"vcfs/{OUTPUTbasename}_{species}-miss1-withNA.vcf.gz"
	shell:
		"bgzip {input}"

# ------- Prepare the files for PopGenome -------

rule extractchrvcf:
	""" Extract one individual chromosome from the vcf file """
	input:
		vcf = "vcfs/{OUTPUTbasename}_{species}-miss1-withNA.vcf.gz",
	output:
		"popgenome/chr{nchr}/{species}/{OUTPUTbasename}_{species}-miss1-withNA-chr{nchr}.vcf"
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
		vcf = "popgenome/chr{nchr}/{species}/" + OUTPUTbasename + "_{species}-miss1-withNA-chr{nchr}.vcf",
		bed = "filter/excludedsites_{species}-chr{nchr}.bed",
	output:
		df = "filter/winlens_chr{nchr}_{species}_{WINSIZEkb}kb.txt"
	threads: 12 # some chrs need more memory (most where fine with 6 threads)
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "2-00:00:00",
	params:
		chr = "chr{nchr}",
		width = WINSIZE,
		jump = JUMP, # overlapping windows
	conda: 
		"../07a_3SppIntrogression/envs/popgenome.yaml"
	script:
		PopGenome_WinMissingness

# ------- Calculate divergence -------

rule PopGenome_PairwiseDivergence:
	input:
		vcf = "popgenome/chr{nchr}/{species}/" + OUTPUTbasename + "_{species}-miss1-withNA-chr{nchr}.vcf",
		winlens = "filter/winlens_chr{nchr}_{species}_{WINSIZEkb}kb.txt"
	output:
		df = "reports/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb-{strain1}-{strain2}.csv"
	threads: 2 # some chrs need more memory
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	params:
		width = WINSIZE,
		jump = JUMP, # overlapping windows
		strain1 = "{strain1}",
		strain2 = "{strain2}",
		chr = "chromosome_{nchr}"
	conda: 
		"../07a_3SppIntrogression/envs/popgenome.yaml"
	script:
		PopGenome_PairwiseDivergence

def get_samples_Pairwise(wildcards):
	files = []
	for strain1 in SppStrainDic["anserina"]:
		for strain2 in SppStrainDic[wildcards.species]:
			newfile = f"reports/Dxy_chr{wildcards.nchr}_{wildcards.species}_{wildcards.WINSIZEkb}kb-{strain1}-{strain2}.csv"
			files.append(newfile)
	return(files)

rule cat_PopGenome_PairwiseDivergence:
	input:
		get_samples_Pairwise
	output:
		"results/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb.csv"
	shell:
		"""
		echo '"Strain1","Strain2","position","winlen","chr","dxy","start","end"' > {output}
		cat {input} >> {output}
		"""

rule compress_PairwiseDivergence:
	input:
		"results/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb.csv"
	output:
		"results/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb.csv.gz"
	shell:
		"bgzip {input}"

# ---- 

rule reduce2nonoverlap:
	""" Simplify the file to get non-overlapping windows for both species """
	input:
		expand("results/Dxy_chr{{nchr}}_{species}_{{WINSIZEkb}}kb.csv.gz", species = COMPAREDSPECIES)
	output:
		csv = "results/Dxy_chr{nchr}_{WINSIZEkb}kb_nonoverlap.csv.gz"
	run:
		# Open the output file in gzip format
		with gzip.open(output.csv, 'wt', newline='') as output_file:
			nonoverlaps = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

			# Pre-calculate the inverse of WINSIZE for efficiency (multiplication is less computationally expensive than division)
			inv_WINSIZE = 1 / WINSIZE

			for inputcsv in input:
				# Open the bgzip file and read its contents
				with gzip.open(inputcsv, mode='rt') as file:
					reader = csv.reader(file)
					for row in reader:
						if 'Strain1' in row: # Header
							nonoverlaps.writerow(row)
						else:
							endwin = float(row[7])
							if (endwin * inv_WINSIZE).is_integer(): # The window is not overlapping
								nonoverlaps.writerow(row)


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
		csv = "results/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb.csv.gz"
	output:
		rawbed = temp("filter/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb_introgression.bedraw"),
		missing = temp("filter/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb_missingness.bedraw"),
	params:
		pi = PI,
		minwin = float(MINWINSIZE)
	run:
		# Start the files
		rawbedfile = open(output.rawbed, 'w')
		missfile = open(output.missing, 'w')

		# Open the bgzip file and read its contents
		with gzip.open(input.csv, mode='rt') as file:
			reader = csv.reader(file)
			for row in reader:
				if 'Strain1' not in row: # Avoide the header
					# print(row)
					winlen, chrom, dxy, start, end = row[3:8]
					outline = f"{chrom}\t{int(float(start))}\t{int(float(end))}\t{dxy}\t{wildcards.species}\n"

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
	""" Fuse windows like in `removeoverlap`"""
	input:
		rawbed = "filter/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb_{tipo}.bedraw"
	output:
		bed = "filter/Dxy_chr{nchr}_{species}_{WINSIZEkb}kb_{tipo}.bed"
	run:
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
		expand("filter/Dxy_chr{nchr}_{{species}}_{{WINSIZEkb}}kb_{{tipo}}.bed", nchr = CHRNUMS)
	output:
		"results/Dxy_{species}_{WINSIZEkb}kb_{tipo}.bed"
	shell:
		"cat {input} > {output}"

rule cat_beds_all:
	input:
		expand("results/Dxy_{species}_{{WINSIZEkb}}kb_{{tipo}}.bed", species = COMPAREDSPECIES)
	output:
		temp("filter/Dxy_{WINSIZEkb}kb_{tipo}-raw.bed")
	shell:
		"cat {input} > {output}"

rule removeoverlap_all_spp:
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

# ---- Plot the divergence along chromosomes for chosen regions ----

rule Plot_PaNLRchrPopStats2Spp:
	input:
		chr1sp1 = f"results/Dxy_chr1_{COMPAREDSPECIES[0]}_{WINSIZEkb}kb.csv.gz",
		chr2sp1 = f"results/Dxy_chr2_{COMPAREDSPECIES[0]}_{WINSIZEkb}kb.csv.gz",
		chr3sp1 = f"results/Dxy_chr3_{COMPAREDSPECIES[0]}_{WINSIZEkb}kb.csv.gz",
		chr4sp1 = f"results/Dxy_chr4_{COMPAREDSPECIES[0]}_{WINSIZEkb}kb.csv.gz",
		chr5sp1 = f"results/Dxy_chr5_{COMPAREDSPECIES[0]}_{WINSIZEkb}kb.csv.gz",
		chr6sp1 = f"results/Dxy_chr6_{COMPAREDSPECIES[0]}_{WINSIZEkb}kb.csv.gz",
		chr7sp1 = f"results/Dxy_chr7_{COMPAREDSPECIES[0]}_{WINSIZEkb}kb.csv.gz",
		chr1sp2 = f"results/Dxy_chr1_{COMPAREDSPECIES[1]}_{WINSIZEkb}kb.csv.gz",
		chr2sp2 = f"results/Dxy_chr2_{COMPAREDSPECIES[1]}_{WINSIZEkb}kb.csv.gz",
		chr3sp2 = f"results/Dxy_chr3_{COMPAREDSPECIES[1]}_{WINSIZEkb}kb.csv.gz",
		chr4sp2 = f"results/Dxy_chr4_{COMPAREDSPECIES[1]}_{WINSIZEkb}kb.csv.gz",
		chr5sp2 = f"results/Dxy_chr5_{COMPAREDSPECIES[1]}_{WINSIZEkb}kb.csv.gz",
		chr6sp2 = f"results/Dxy_chr6_{COMPAREDSPECIES[1]}_{WINSIZEkb}kb.csv.gz",
		chr7sp2 = f"results/Dxy_chr7_{COMPAREDSPECIES[1]}_{WINSIZEkb}kb.csv.gz",
		pop1 = path2divstats + "/PodoPop-snps-withNA-chr1-stats.csv",
		pop2 = path2divstats + "/PodoPop-snps-withNA-chr2-stats.csv",
		pop3 = path2divstats + "/PodoPop-snps-withNA-chr3-stats.csv",
		pop4 = path2divstats + "/PodoPop-snps-withNA-chr4-stats.csv",
		pop5 = path2divstats + "/PodoPop-snps-withNA-chr5-stats.csv",
		pop6 = path2divstats + "/PodoPop-snps-withNA-chr6-stats.csv",
		pop7 = path2divstats + "/PodoPop-snps-withNA-chr7-stats.csv",
		nlrmeta = NLRmetadata,
		indvs = samples_file
	output:
		main = "figures/Fig6_Stats_NLRs.png",
		figs7 = "figures/FigS7_Stats_NLRs_nwd5.png",
		figs8 = "figures/FigS8_Stats_NLRs_hetb.png",
		hetz = "figures/FigS9_Stats_NLR_hetz.png",
		pnpudp = "figures/Stats_NLRs_Pa_2_8180.png",
		Fnt1 = "figures/Stats_NLRs_Fnt1.png"
	conda: 
		"envs/plot.yaml"
	script:
		PaNLRchrPopStats2Spp	


