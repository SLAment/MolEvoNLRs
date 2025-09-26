# -*- snakemake -*-

### GenomeIdentity2024: A pipeline to compare identity along the genome against a reference, 2024 version
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

# -------------------------------------------------
## Variables set from the configuration file
configfile: "config/config.yaml"
# -------------------------------------------------
samples_file = config["SampleIDs_file"]

# Illumina reads path:
IlluminaPath = config["IlluminaPath"]

# The reference genome
refgenome = config["refgenome"]
refname = config["refname"]
mitoctg = config["mitoctg"]

# TE library
TElib = config["TElib"]

# Variables
WINSIZE = config["WINSIZE"]

# Scripts
vcfR_plotter = config["vcfR_plotter"]
PopGenome_WinMissingness = config["PopGenome_WinMissingness"]
PopGenome_PairwiseDivergence = config["PopGenome_PairwiseDivergence"]
PaIntrogressionTracks = config["PaIntrogressionTracks"]
PaIntrogressionProportions = config["PaIntrogressionProportions"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: rawdata_illumina, referencegenome, indexbwa, bedtoolsTEs, bgzip_tabix, VariantFiltration, getPASSsites, bgzip_tabix2, makeintervals, filtergenomecov, totalcovergff, removeoverlap, extractchrcov, get_totalcovergff, CatGetCoverage, get_badsites2vcf, cat_PopGenome_PairwiseDivergence, compress_PairwiseDivergence, put_together_PairwiseDivergence
# ----------

SppStrainDic = {}
for line in open(samples_file, 'r'):
	tab = line.rstrip("\n").split('\t')
	if tab[1] in SppStrainDic.keys():
		SppStrainDic[tab[1]].append(tab[0])
	else:
		SppStrainDic[tab[1]] = [tab[0]]

otherSpp = []
samples = []
for key in SppStrainDic.keys():
	samples.extend(SppStrainDic[key])
	if key != "anserina":
		otherSpp.extend(SppStrainDic[key])

# samples = [line.rstrip("\n").split('\t')[0] for line in open(samples_file, 'r')] 			# Read tab file into a list

WINSIZEkb = int(round(WINSIZE/1000))

rule all:
	input:
		"figures/FigS5_DxySppComplex_chr4.png"
		

# ------- PREPARE ALL DATA --------
rule rawdata_illumina:
	""" Prepare a folder with the Illumina data files """
	output:
		read1 = "data/Illumina/{sample}_postQC.1.fq.gz",
		read2 = "data/Illumina/{sample}_postQC.2.fq.gz"	
	params:
		illuminapath = IlluminaPath
	shell:
		"""
		ln -sf {params.illuminapath}/{wildcards.sample}_postQC.1.fq.gz data/Illumina
		ln -sf {params.illuminapath}/{wildcards.sample}_postQC.2.fq.gz data/Illumina
		"""

rule referencegenome:
	""" Prepare a link of the reference genome """
	output:
		"data/reference/" + refname + ".fa"
	params:
		refgenome = refgenome
	shell:
		"""
		ln -s {params.refgenome} {output}
		"""

# ---------------------------------

rule indexbwa:
	""" Index genome with BWA """
	input:
		genome = "data/reference/" + refname + ".fa"
	output:
		index = "data/reference/" + refname + ".fa.bwt"
	shell:
		"""
		bwa index {input.genome}
		"""

rule bwa_mem:
	""" Map Illumina reads with BWA """
	input:
		genome = "data/reference/" + refname + ".fa",
		index = "data/reference/" + refname + ".fa.bwt",
		read1 = IlluminaPath + "/{sample}_postQC.1.fq.gz",
		read2 = IlluminaPath + "/{sample}_postQC.2.fq.gz",
	output:
		bwaoutput = temp("mapping/{sample}/{sample}-to-" + refname + ".sorted.bam"),
	log:
		"logs/bwa_mem/{sample}.log"
	threads: 10,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "3:30:00",
	params:
		refbase = refname,
		rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina",
	shell:
		"""
		(bwa mem {input.genome} {input.read1} {input.read2} -t {threads} -R '{params.rg}' -M | samtools view -Su - | samtools sort -l 5 -O bam -T {wildcards.sample}'-to-'{params.refbase} -@ {threads} > {output.bwaoutput}) 2> {log}
		# -l 5 following Doug
		"""

rule markduplicates:
	""" Mark duplicates in BAM """
	# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
	input:
		bwaoutput = "mapping/{sample}/{sample}-to-" + refname + ".sorted.bam"
	output:
		mdoutput = "mapping/{sample}/{sample}-to-" + refname + ".sorted.debup.bam",
		mdmetrics = "mapping/{sample}/{sample}-to-" + refname + ".sorted.metrics.txt"
	threads: 10,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "1:30:00",
	params:
		unpatterned = ""
	run:
		# Some samples are older and come from unpatterned HiSeq platforms
		if wildcards.sample in ["PaWa1p", "PaWa1m", "PaWa53p", "PaWa53m", "PaWa87p", "PaWa87m", "PaWa28p", "PaWa28m", "PaWa85p", "PaWa85m", "PaWa21p", "PaWa21m", "PaWa27p", "PaWa27m", "PaWa46p", "PaWa46m", "PaYp", "PaYm", "PaWa47p", "PaWa47m", "PaWa100p", "PaWa100m", "PaWa58p", "PaWa58m", "PaZp", "PaWa3p", "PaWa3m", "PaWa32p", "PaWa32m", "PaWa63p", "PaWa63m", "Pa170m", "Pa130p", "Pa200p", "Pa180p"]:
			shell(f"picard MarkDuplicates I={input.bwaoutput} O={output.mdoutput} M={output.mdmetrics} ASSUME_SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR='temp' REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500")
		else: # The rest are HiSeq X (patterned)
			shell(f"picard MarkDuplicates I={input.bwaoutput} O={output.mdoutput} M={output.mdmetrics} ASSUME_SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR='temp' REMOVE_DUPLICATES=true")

		# """
		# if 
		# picard MarkDuplicates -I {input.bwaoutput} -O {output.mdoutput} -M {output.mdmetrics} --ASSUME_SORT_ORDER coordinate --CREATE_INDEX true --TMP_DIR "temp" --REMOVE_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

		# """	
		# # VALIDATION_STRINGENCY=ValidationStringency
		# #                               Validation stringency for all SAM files read by this program.  Setting stringency to
		# #                               SILENT can improve performance when processing a BAM file in which variable-length data
		# #                               (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This
		# #                               option can be set to 'null' to clear the default value. Possible values: STRICT,
		# #                               LENIENT, SILENT
		# # CREATE_INDEX=Boolean          Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value:
		# #                               false. This option can be set to 'null' to clear the default value. Possible values:
		# #                               true, false
		# # TMP_DIR (File)  Default value: null. This option may be specified 0 or more times.
		# # --REMOVE_DUPLICATES:Boolean   If true do not write duplicates to the output file instead of writing them with
        # # 		                      appropriate flags set.  Default value: false. Possible values: {true, false}
        # # --OPTICAL_DUPLICATE_PIXEL_DISTANCE:Integer
        # #                               The maximum offset between two duplicate clusters in order to consider them optical
        # #                               duplicates. The default is appropriate for unpatterned versions of the Illumina platform.
        # #                               For the patterned flowcell models, 2500 is more appropriate. For other platforms and
        # #                               models, users should experiment to find what works best.  Default value: 100.
        # #	--ASSUME_SORT_ORDER,-ASO:SortOrder
        # #                               If not null, assume that the input file has this order even if the header says otherwise.
        # #                               Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate,
        # #                               unknown}  Cannot be used in conjunction with argument(s) ASSUME_SORTED (AS)

		# Patterned platforms: the NovaSeq 6000 System, the NovaSeq 5000 System, the HiSeq X System, the HiSeq 4000 System, and the HiSeq 3000 System
		# Most of the pop Podospora data is HiSeq X. However the very first strains were sequenced with an unkown earlier version of HiSeq. It might have been anything between HiSeq 2000 and 3000/4000

# -------------- GATK4 -----------------
rule indexsanddict:
	""" Index reference for GATK """ 
	input:
		genome = "data/reference/" + refname + ".fa",
	output:
		indexsamtools = "data/reference/" + refname + ".fa.fai",
		diction = "data/reference/" + refname + ".dict"
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "1:00:00",
	shell:
		"""
		# Make a reference index
		samtools faidx {input.genome}

		# Make a reference dictionary
		picard CreateSequenceDictionary R={input.genome} O={output.diction}
		"""	

rule HaplotypeCaller:
	""" Produce a GVCF file from BAM - haploid """
	# https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller
	input:
		bam = "mapping/{sample}/{sample}-to-" + refname + ".sorted.debup.bam",
		ref = "data/reference/" + refname + ".fa",
		indexsamtools = "data/reference/" + refname + ".fa.fai",
		diction = "data/reference/" + refname + ".dict"
	output:
		gvcf = "gvcfs/{sample}.g.vcf",
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "1-10:00:00",
	params:
		JavaMem = lambda wildcards, threads: int(threads * 6.8),
		ploidy = 1
	shell:
		"""
		gatk --java-options "-Xmx{params.JavaMem}G" HaplotypeCaller \\
		-I {input.bam} -R {input.ref} \\
		-O {output.gvcf} \\
		-ploidy {params.ploidy} \\
		-ERC GVCF

		# -nct {threads} \\
		# --bam-output,-bamout:String   File to which assembled haplotypes should be written  Default value: null.
		# --annotateNDA 	Annotate number of alleles observed
		# --useNewAFCalculator	Use new AF model instead of the so-called exact model
		# --emitRefConfidence GVCF 	Mode for emitting reference confidence scores
		"""

rule makeintervals:
	""" Make an interval file for GenomicsDBImport """
	input:
		"data/reference/" + refname + ".fa.fai"
	output:
		"data/reference/" + refname + ".intervals"
	shell:
		"cat {input} | cut -f1 > {output}"

rule GenomicsDBImport:
	""" Import single-sample GVCFs into GenomicsDB before joint genotyping """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.1/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php
	input:
		intervals = "data/reference/" + refname + ".intervals",
		gvcfs = expand("gvcfs/{sample}.g.vcf", sample = samples)
	output:
		directory("genomicsdb")
	threads: 2,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "2:00:00",
	params:
		JavaMem = lambda wildcards, threads: int(threads * 6.8),
		ploidy = 1
	run:
		# Create a string in the format --variant path/to/gvcf/sample1 --variant path/to/gvcf/sample2 etc...
		variantlist = ""
		for sample in input.gvcfs:
			variantlist += "--variant " + sample + " "
		
		shell("mkdir -p tmp")
		gatkcommand = f'gatk GenomicsDBImport --java-options "-Xmx{params.JavaMem}G" --tmp-dir tmp --genomicsdb-workspace-path {output[0]} -L {input.intervals} {variantlist}'
		shell(gatkcommand) # Notice GATK will create the output directory


rule GenotypeGVCFs:
	""" Perform joint genotyping """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php
	input:
		my_database = "genomicsdb",
		ref = "data/reference/" + refname + ".fa",
	output:
		rawvcf = temp("vcfs/PodoPop-vs-" + refname + "raw_align.vcf.gz")
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "1-00:00:00",
	params:
		JavaMem = lambda wildcards, threads: int(threads * 6.8),
		ploidy = 1
	shell:
		"""
		gatk --java-options "-Xmx{params.JavaMem}G" GenotypeGVCFs \\
		-R {input.ref} \\
		-V gendb://{input.my_database} \\
		-O {output.rawvcf} \\
		-ploidy {params.ploidy} \\
		--create-output-variant-index 

		# --create-output-variant-index	If true, create a VCF index when writing a coordinate-sorted VCF file.
		# --use-new-qual-calculator / -new-qual 	Use the new AF model instead of the so-called exact model. Default: true
		# By default, GATK HaplotypeCaller and GenotypeGVCFs do not emit variants with QUAL < 10, controlled with -stand-call-conf
		"""

rule LeftAlign: # The same indel can often be placed at multiple positions and still represent the same haplotype.
	""" Left align INDELs and trim alleles """
	# https://gatk.broadinstitute.org/hc/en-us/articles/21905044109467-LeftAlignAndTrimVariants
	input:
		ref = "data/reference/" + refname + ".fa",
		vcf = "vcfs/PodoPop-vs-" + refname + "raw_align.vcf.gz"
	output:
		vcf = "vcfs/PodoPop-vs-" + refname + ".vcf.gz"
	threads: 2,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "10:00:00",
	shell:
		"gatk LeftAlignAndTrimVariants -R {input.ref} -V {input.vcf} -O {output.vcf}" 
		# --dont-trim-alleles, -no-trim	false	Do not Trim alleles to remove bases common to all of them

# ------- RepeatMasking -------

rule repeatmasker:
	""" Use RepeatMasker to find regions that should be filtered out """
	# https://www.animalgenome.org/bioinfo/resources/manuals/RepeatMasker.html
	input:
		"data/reference/" + refname + ".fa"
	output:
		"RepeatMasker/" + refname + ".fa.out.gff"
	threads: 16,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "2:00:00",
	params:
		TElib = TElib
	shell:
		""" 
		RepeatMasker -pa {threads} -a -xsmall -gccalc -gff -excln -lib {params.TElib} -dir RepeatMasker {input}
		"""
	# -gccalc        Program calculates the GC content even for batch files/small 
	#                   sequences.
	# -a      shows the alignments in a .align output file; -ali(gnments) also works
	# -xsmall returns repetitive regions in lowercase (rest capitals) rather than masked
	# -gff    creates an additional General Feature Finding format output
	# -excln	The percentages displayed in the .tbl file are calculated using a total sequence length excluding runs of 25 Ns or more.

# ------- Filtering -------

rule snpsvcf:
	""" Filter resulting vcf file to get only SNPs """
	input:
		rawvcf = "vcfs/PodoPop-vs-" + refname + ".vcf.gz"
	output:
		snpsvcf = temp("vcfs/PodoPop-vs-" + refname + "-snps-raw.vcf.gz"),
	threads: 1, # It needs a lot of memory
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "20:00", # The file can be pretty big
	params:
		mitoctg = mitoctg
	shell:
		"vcftools --gzvcf {input.rawvcf} --remove-indels --recode --recode-INFO-all --not-chr {params.mitoctg} --stdout | bgzip > {output.snpsvcf}" # Notice I filtered out the mitochondria too

rule fix_missingdata:
	""" GATK changed the formatting of missing data and now it's not standard, ugh. Fix it. """
	# https://github.com/broadinstitute/gatk/issues/8328
	# https://twitter.com/ksamuk/status/1754698564884996590
	input:
		"vcfs/PodoPop-vs-" + refname + "-snps-raw.vcf.gz",
	output:
		"vcfs/PodoPop-vs-" + refname + "-snps.vcf.gz",
	threads: 1, # It needs a lot of memory
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "10:00",
	shell:
		"bcftools +setGT {input} -- -t q -n . -e 'FMT/DP>=1' > {output}"
	# -t, --target-gt <type>      Genotypes to change, see above
	#	q    .. select genotypes using -i/-e options
	# -n, --new-gt <type>         Genotypes to set, see above
	# 	.    .. partially or completely missing
	# -e, --exclude <expr>        Exclude a genotype if true (requires -t q)


### ---- Use coverage to filter bad sites


# rule get_vcfR_plotter:
# 	output:
# 		temp("scripts/DiversityStats_vcfR_plotter_raw.R")
# 	shell:
# 		"wget -O {output} https://raw.githubusercontent.com/johannessonlab/HetVPaper/master/DiversityStats/scripts/DiversityStats_vcfR_plotter.R"

# rule fix_vcfR_plotter:
# 	input:
# 		"scripts/DiversityStats_vcfR_plotter_raw.R"
# 	output:
# 		temp("scripts/DiversityStats_vcfR_plotter.R")	
# 	shell:
# 		"sed 's/input\\[\\[1\\]\\]/input$vcf/' {input} > {output}"

rule break_vcf:
	""" The vcf file is too heavy for my little R script so let's break it by sample first """
	input:
		vcf = "vcfs/PodoPop-vs-" + refname + "-snps.vcf.gz",
	output:
		vcf = temp("temp/PodoPop-vs-" + refname + "-snps_{sample}.vcf.gz")
	shell:
		"bcftools view -s {wildcards.sample} {input} | bgzip > {output}"

rule GetCoverage_persample:
	""" Get coverage distribution values per sample """
	input:
		vcf = "temp/PodoPop-vs-" + refname + "-snps_{sample}.vcf.gz",
	output:
		table = "processing/cov/CoverageSummary_{sample}.tab",
	threads: 2, # It needs a lot of memory
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "1:00:00",
	params:
		lowquantile = 0.25,
		uppquantile = 0.995,
		sample = "{sample}"
	conda: 
		"envs/popgenome.yaml"
	script:
		vcfR_plotter

rule CatGetCoverage:
	input:
		expand("processing/cov/CoverageSummary_{sample}.tab", sample = samples)
	output:
		"data/cov/CoverageSummary.tab"
	shell:
		"printf 'sample\\tmedian\\tmean\\tlowquantile\\tupquantile\\n' > {output}; "
		"cat {input} >> {output}"

# rule PlotCoverage:
# 	""" Plot coverage distribution for all samples and report """
# 	input:  # In this order!
# 		vcfR_plotter = "scripts/DiversityStats_vcfR_plotter.R",
# 		vcf = "temp/PodoPop-vs-" + refname + "-snps_{sample}.vcf.gz",
# 	output:
# 		"figures/Coverage/Coverage_{sample}.pdf",
# 		"processing/cov/CoverageSummary_{sample}.tab",
# 	threads: 8, # It needs a lot of memory
# 	resources:
# 		threads = lambda wildcards, threads: threads,
# 		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
# 		time = "4:00:00",
# 	params:
# 		lowquantile = 0.25,
# 		uppquantile = 0.995,
# 	conda: 
# 		"envs/popgenome.yaml"
# 	script:
# 		"scripts/DiversityStats_vcfR_plotter.R" # I can only make it work if I give the name of the input file directly


rule genomecov:
	""" genomecov makes BED with intervals of coverage """ 
	input:
		"mapping/{sample}/{sample}-to-" + refname + ".sorted.debup.bam"
	output:
		"genomecov/{sample}_cov.bed"
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "15:00",
	shell:
		"bedtools genomecov -ibam {input} -bga > {output}"
		# bga - report coverage in intervals, including sites with 0 coverage

rule filtergenomecov:
	""" Filter the coverage sites based on the quantile of each sample's distribution """
	input:
		covsum = "data/cov/CoverageSummary.tab",
		bed = "genomecov/{sample}_cov.bed",
	output:
		"genomecov/{sample}_cov_filtered.bed"
	shell:
		"""
		lowquantile=$(grep {wildcards.sample} {input.covsum} | cut -f4)
		uppquantile=$(grep {wildcards.sample} {input.covsum} | cut -f5)

		awk -v low="$lowquantile" -v upp="$uppquantile" '$4 < low || $4 > upp' {input.bed} > {output}
		"""

rule get_totalcovergff:
	output:
		temp("scripts/totalcovergff.py")
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/f6923f552d643414e20602e6d3aeadd7840b982c/GenomeAnnotation/totalcovergff.py"

rule totalcovergff:
	""" Make bed file of sites overlapping transposable elements """
	input:
		totalcovergff = "scripts/totalcovergff.py",
		gff = "RepeatMasker/" + refname + ".fa.out.gff"
	output:
		"data/cov/TEsites.bed"
	shell:
		"python {input.totalcovergff} {input.gff} > {output}"

rule BEDOPS:
	""" Take all the BED files and overlap them into one single BED file"""
	# https://bedops.readthedocs.io/en/latest/content/overview.html#overview
	# https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#merge-m-merge
	input:
		"data/cov/TEsites.bed",
		expand("genomecov/{sample}_cov_filtered.bed", sample = samples)
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

rule bedtoolsTEs:
	""" Filter vcf with the repeats from RepeatMasker """
	input:
		vcf = "vcfs/PodoPop-vs-" + refname + "-snps.vcf.gz",
		gfffile = "RepeatMasker/" + refname + ".fa.out.gff"
	output:
		filteredvcf = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs.vcf",
	shell:
		"""
		# Get the header
		bcftools view -h {input.vcf} > {output.filteredvcf}
		
		# Filter out the repeats
		bedtools intersect -a {input.vcf} -b {input.gfffile} -v >> {output.filteredvcf}
		"""

rule bgzip_tabix:
	""" Compress and index vcf file """
	input:
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs.vcf"
	output:
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs.vcf.gz"
	shell:
		"bgzip {input}; "
		"tabix -p vcf {output}"


rule VariantFiltration:
	""" Use GATK's VariantFiltration to mark out sites by INFO annotations """
	input:
		vcf = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs.vcf.gz",
	output:
		filteredvcf = temp("vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatk.vcf.gz"),
	shell:
		"""
		gatk VariantFiltration \\
		-V {input.vcf} \\
		-O {output.filteredvcf} \\
		--filter-expression "QD < 2.0" --filter-name "QD2" \\
		--filter-expression "FS > 60.0" --filter-name "FS60" \\
		--filter-expression "MQ < 40.0" --filter-name "MQ40" \\
		--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \\
		--filter-expression "SOR > 3.0" --filter-name "SOR3" \\
		--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
		
		##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
		##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
		##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
		##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
		##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
		# RankSum annotations can only be calculated for REF/ALT heterozygous sites and therefore will be absent from records that do not present read counts towards heterozygous genotypes.
		##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
		##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
		##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">	
		# In practice, we only filter out low negative values when evaluating
		# variant quality because the idea is to filter out variants for which the
		# quality of the data supporting the alternate allele is comparatively low.
		# The reverse case, where it is the quality of data supporting the reference
		# allele that is lower (resulting in positive ranksum scores), is not really
		# informative for filtering variants.
		##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
		"""

rule getPASSsites:
	""" Remove sites rejected in VariantFiltration """
	input:
		vcf = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatk.vcf.gz",
	output:
		vcf = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS.vcf"
	shell:
		"""
		# Get the header
		bcftools view -h {input.vcf} > {output.vcf}
		
		# Filter out the repeats
		bcftools view -H {input.vcf} | grep 'PASS' >> {output.vcf}
		"""

rule bgzip_tabix2:
	""" Compress and index vcf file """
	input:
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS.vcf"
	output:
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS.vcf.gz"
	shell:
		"bgzip {input}; "
		"tabix -p vcf {output}"

rule snpsvcfnomiss:
	""" Filter out sites with missing data """
	input:
		vcf = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS.vcf.gz"
	output:
		filteredvcf = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1.vcf.gz"
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	shell:
		"vcftools --gzvcf {input.vcf} --max-missing 1 --recode --recode-INFO-all --stdout | bgzip > {output.filteredvcf}; "
		"tabix -p vcf {output.filteredvcf}"

rule get_badsites2vcf:
	output:
		temp("scripts/badsites2vcf.py")
	shell:
		"wget -O {output} https://raw.githubusercontent.com/johannessonlab/HetVPaper/master/DiversityStats/scripts/badsites2vcf.py"

rule inputmissingsites:
	""" Input the sites that should be excluded as explicit missing data sites """
	input:
		vcf = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1.vcf.gz",
		bed = "filter/excludedsites.bed",
		badsites2vcf = "scripts/badsites2vcf.py"
	output:
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA.vcf"
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
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA.vcf"
	output:
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA.vcf.gz"
	shell:
		"bgzip {input}"

rule tabix:
	""" Sanity check, it shouldn't complain """
	input:
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA.vcf.gz"
	output:
		"vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA.vcf.gz.tbi"
	shell:
		"tabix -p vcf {input}" # This is a sanity check, it shouldn't complain	

# ------- Prepare the files for PopGenome -------

rule extractchrvcf:
	""" Extract one individual chromosome from the vcf file """
	input:
		vcf = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA.vcf.gz",
		tbi = "vcfs/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA.vcf.gz.tbi" # just to force the previous rule to produce it
	output:
		"popgenome/chr{nchr}/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA-chr{nchr}.vcf"
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
		vcf = "popgenome/chr{nchr}/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA-chr{nchr}.vcf",
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
	input:
		vcf = "popgenome/chr{nchr}/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA-chr{nchr}.vcf",
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
		expand("reports/Dxy_chr{{nchr}}_{WINSIZEkb}kb-{strain1}-{strain2}.csv", WINSIZEkb = WINSIZEkb, strain1 = SppStrainDic["anserina"], strain2 = otherSpp),
		expand("reports/Dxy_chr{{nchr}}_{WINSIZEkb}kb-{strain1}-{strain2}.csv", WINSIZEkb = WINSIZEkb, strain1 = SppStrainDic["pauciseta"], strain2 = SppStrainDic["comata"]),
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


# ------- Plot divergence for chromosome 4 -------

rule Plot_PaIntrogressionProportions:
	""" Plot the proportion of introgressed regions between P. anserina and other Podospora species """
	input:
		dxy = expand("results/Dxy_chr{nchr}_{WINSIZEkb}kb.csv.gz", nchr = 4, WINSIZEkb = WINSIZEkb),
		indvs = samples_file
	output:
		plot = "figures/FigS5_DxySppComplex_chr4.png"
	conda: 
		"envs/simpleplot.yaml"
	script:
		PaIntrogressionProportions

# ------- Put together the output of the other chromosomes, but I didn't use it in the end -------

# rule put_together_PairwiseDivergence:
# 	input:
# 		expand("results/Dxy_chr{nchr}_{{WINSIZEkb}}kb.csv", nchr = list(range(1,8)))
# 	output:
# 		"results/Dxy_{WINSIZEkb}kb.csv"
# 	shell:
# 		"""
# 		echo '"Strain1","Strain2","position","winlen","chr","dxy","start","end"' > {output}
# 		cat {input} | grep -v 'Strain1' >> {output}
# 		"""

# rule final_compress:
# 	input:
# 		"results/Dxy_{WINSIZEkb}kb.csv"
# 	output:
# 		"results/Dxy_{WINSIZEkb}kb.csv.gz"
# 	shell:
# 		"bgzip {input}"	


# rule CalcDivergence:
# 	""" Calculate diversity statistics """
# 	input:  # In this order!
# 		vcf = "popgenome/chr{nchr}/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-withNA-chr{nchr}.vcf",
# 		bed = "filter/excludedsites-chr{nchr}.bed",
# 	output:
# 		pi = "results/stats/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-chr{nchr}-pi.csv",
# 		dxy = "results/stats/PodoPop-vs-" + refname + "-snps-NoTEs-gatkPASS-miss1-chr{nchr}-dxy.csv",
# 	params:
# 		time = "1:30:00",
# 		threads = 3 #5, # I need to give these many in case memory is needed
# 	conda: 
# 		"envs/statsplot.yaml"
# 	script:
# 		PaIntrogressionTracks
