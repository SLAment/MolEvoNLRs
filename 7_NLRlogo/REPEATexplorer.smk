### REPEATexplorer: Exploring the diversity of domain repeats in the Podospora anserina species complex
#############################################################################
#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# 2024-08-14
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 3

# from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import re

## For the Psim calculation
# # Nucleotide motif version
# from Bio import AlignIO, SeqIO
# from Bio.Align.Applications import ClustalwCommandline
# from Bio.Seq import Seq
# from Bio.motifs import create

# Depricated but functional version
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment

# I get a lot of warnings because I use depricated code
import warnings
from Bio import BiopythonDeprecationWarning
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"

# Input files
basealignment = config["basealignment"]

# Scripts
descriptor = config["descriptor"]
REPEAT_LOGO = config["REPEAT_LOGO"]

## Variables
GENE = config["GENE"]

# Definition of the repeats
REGEX = config["REGEX"]
REPLEN = config["REPLEN"]
REPMINLEN = config["REPMINLEN"]
REPMAXLEN = config["REPMAXLEN"]

# The amino acid positions used to classify the repeats (influenced by the REGEX definition of the repeat)
AMINOS = config["AMINOS"]

# -------------------------------------------------


rule all:
	input:
		f"results/{GENE}_LOGO.pdf",
		expand("reports/Psim_{aminoacids}_anserina.txt", aminoacids = AMINOS)


rule descriptor:
	""" Characterize the WD40 repeats """
	# The alignment is special: only the WD40 repeats, no more, no less, and aligned, but without totally empty columns
	input:
		al = basealignment,
		descriptor = descriptor
	output:
		report = "reports/REPEATreport_{aminoacids}.txt",
		fasta = "alignments/REPEATrepeats_{aminoacids}.fa",
		meta = "alignments/REPEATrepeats_{aminoacids}_metadata.txt",
	shell:
		"""
		aminos=$(echo {wildcards.aminoacids} | sed 's/-/,/g'); 
		printf "{input.descriptor} {input.al} -o 'alignments/REPEATrepeats_{wildcards.aminoacids}' -a $aminos --regex '{REGEX}' --replen {REPLEN} --repminlen {REPMINLEN} --repmaxlen {REPMAXLEN} \n" 

		python {input.descriptor} {input.al} -o 'alignments/REPEATrepeats_{wildcards.aminoacids}' -a $aminos --regex '{REGEX}' --replen {REPLEN} --repminlen {REPMINLEN} --repmaxlen {REPMAXLEN}  > {output.report}
		"""
	

# Function to translate nucleotide sequence considering gaps and frameshifts
# Thanks ChatGPT
# See also https://chatgpt.com/c/eba0aff2-cb16-4047-b1e1-08cbef511c39
def custom_translate(nuc_seq): 
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    protein_seq = []

    # Translate each codon
    for i in range(0, len(nuc_seq), 3):
        codon = str(nuc_seq[i:i+3])
        if len(codon) < 3 or '-' in codon:  # Incomplete codon or indel
            protein_seq.append('X')  # 'X' for unknown amino acid
        else:
            aa = standard_table.forward_table.get(codon, 'X')  # Translate or 'X' if invalid codon
            protein_seq.append(aa)

    return ''.join(protein_seq)


rule anserina_alignment:
	""" Get an anserina-only alignment """
	input:
		fasta = "alignments/REPEATrepeats_{aminoacids}.fa",
		meta = "alignments/REPEATrepeats_{aminoacids}_metadata.txt",
	output:
		fasta = "alignments/PerSpecies/REPEATrepeats_{aminoacids}_anserina.fa",
		report = "reports/REPEATrepeats_{aminoacids}_counts.txt"
	run:
		tabs = [line.rstrip("\n").split("\t") for line in open(input.meta, 'r')] 			# Read tab file into a list

		# Filter for the sequences that are part of the paper (the backcrosses only count for het-e and het-d)
		accepted_strains = ["Podan2", "PaYp", "PaZp.R10", "PaTgp", "PaWa21m", "PaWa28m", "PaWa46p", "PaWa53m", "PaWa58m", "PaWa63p", "PaWa87p", "PaWa100p", "PaWa137m"] # , "CmEmm.R10", "CoEcp.R10", "CoEfp.R10", "ChEhDap.R10", "CaDam.R10", "CsDfp.R10"

		ans_tab = []
		genedic = {}
		for tab in tabs[1:]: # ignore the header (first line)
			SequenceName = tab[0]
			Strain = SequenceName.split('__')[1].split('_')[0]
			Species = tab[3]
			Gene = tab[4]

			if Strain in accepted_strains:
				ans_tab.append(SequenceName)

				if Gene in genedic.keys():
					genedic[Gene].append(SequenceName)
				else:
					genedic[Gene] = [SequenceName]

		# This stores in memory
		records_dict = SeqIO.to_dict(SeqIO.parse(input.fasta, "fasta"))

		# All sequences
		total_anserina = 0
		total_anserina_filtered = 0
		outputfasta = open(output.fasta, "w")
		for seq in ans_tab:
			total_anserina += 1
			protein_seq = custom_translate(records_dict[seq].seq)

			# Accept all repeats
			SeqIO.write(records_dict[seq], outputfasta, "fasta")
			total_anserina_filtered +=1

			# if 'X' not in protein_seq: # Keep only repeats that are intact
			# 	SeqIO.write(records_dict[seq], outputfasta, "fasta")
			# 	total_anserina_filtered +=1

		# print a little report
		reportfile = open(output.report, 'w')
		reportfile.write(f"Total number of processed repeats: {len(records_dict.keys())}\n")
		reportfile.write(f"Total number of anserina repeats: {total_anserina}\n")
		reportfile.write(f"Total number of anserina filtered repeats: {total_anserina_filtered}\n")


rule nc2aa:
	""" I want to have the repeats in amino acids as well """
	input:
		fasta = "alignments/PerSpecies/REPEATrepeats_{aminoacids}_anserina.fa"
	output:
		fasta = "alignments/PerSpecies/Protein/REPEATrepeats_{aminoacids}_anserina_aa.fa"
	run:
		# Initialize a list to hold the protein sequences
		protein_records = []

		ofile = open(output.fasta, 'w')
		# Read and translate the nucleotide sequences
		for record in SeqIO.parse(input.fasta, "fasta"):
			# Translate the nucleotide sequence to a protein sequence with custom function to tolerate frame-shifts
			protein_seq = custom_translate(record.seq)
			protseq = ">" + record.id + "\n" + protein_seq + "\n" # Using f-strings makes weird errors with appending extra spaces 
			ofile.write(protseq)

rule aa2aln:
	""" align the repeats """
	# Careful!! The names of the sequences get chopped and lost, but it's not important for the logo
	input:
		fasta = "alignments/PerSpecies/Protein/REPEATrepeats_{aminoacids}_anserina_aa.fa"
	output:
		fasta = "alignments/PerSpecies/Protein/REPEATrepeats_{aminoacids}_anserina_aa_al.fa",
		aln = temp("alignments/PerSpecies/Protein/REPEATrepeats_{aminoacids}_anserina_aa.aln"),
		dnd = temp("alignments/PerSpecies/Protein/REPEATrepeats_{aminoacids}_anserina_aa.dnd"), # remove the guide tree that ClustalW uses to guide the alignment process
	run:
		# Perform Multiple Sequence Alignment (you need ClustalW installed)
		clustalw_cline = ClustalwCommandline("clustalw2", infile=input.fasta)
		stdout, stderr = clustalw_cline()

		alignment = AlignIO.read(output.aln, "clustal") # Careful!! The names of the sequences get chopped and lost!
		# alnname = os.path.splitext(input.fasta)[0] + ".aln"

		# Read the alignment and make it a fasta file 
		alignment = AlignIO.read(output.aln, "clustal")
		AlignIO.write(alignment, output.fasta, "fasta")

rule makeLOGOS:
	input:
		fasta = expand("alignments/PerSpecies/Protein/REPEATrepeats_{aminoacids}_anserina_aa_al.fa", aminoacids = AMINOS),
	output:
		logo = f"results/{GENE}_LOGO.pdf"
	params:
		gene = GENE
	script:
		REPEAT_LOGO

### ---- Calculate the Psim metric of Jorda & Kajava (2009) using in R-REKS  https://academic.oup.com/bioinformatics/article/25/20/2632/193638 ----


# Calculate Hamming distance for each sequence
def hamming_distance(seq1, seq2):
	assert len(seq1) == len(seq2)
	return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def Psim(alignment):
	# ---- Depricated version, but it works with this biopython version
	# Generate consensus sequence
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = summary_align.gap_consensus(threshold=0) # the ties sites become X

	# # --- Updated version, but only works with nucleotides
	# # For some reason I can't import Bio.Align from Bio, so I have to read the alignment differently
	# # Read the sequences from a FASTA file
	# alignment = []
	# with open(output.aln) as handle:
	# 	for record in SeqIO.parse(handle, "clustal"):
	# 		print(record)
	# 		alignment.append(record.seq)

	# # Create a Motif object from the sequences
	# motif = create(alignment)

	# # Get the consensus sequence
	# consensus = motif.consensus
	# print(f"Consensus sequence of {wildcards.gene}:", consensus)
	# ---- 
	
	D = [hamming_distance(str(consensus), str(aln_seq.seq)) for aln_seq in alignment]

	# Calculate P_sim
	m = len(alignment)
	l = len(consensus)
	N = m * l

	P_sim = (N - sum(D)) / N

	return(P_sim)

rule calculate_psim:
	# https://biopython.org/docs/1.75/api/Bio.Align.html
	input:
		fasta = "alignments/PerSpecies/Protein/REPEATrepeats_{aminoacids}_anserina_aa_al.fa",
	output:
		report = "reports/Psim_{aminoacids}_anserina.txt"
	run:
		alignment = AlignIO.read(input.fasta, "fasta")

		P_sim = Psim(alignment)

		ofile = open(output.report, 'w')
		ofile.write(f"Gene\tSpecies\tPsim\tNo_repeats\n")
		ofile.write(f"{GENE}\tanserina\t{P_sim}\t{len(alignment)}\n")
	
		# ---- Re-calculate without the flanking repeat
		# A key assumption here is that the repeats in the alignment follow
		# the linear arrangement of the repeats in each allele

		straindic = {}
		for record in alignment:
			strain = record.id.split('__')[1].split('_')[0]
			
			if strain not in straindic.keys():
				straindic[strain] = [record]
			else:
				straindic[strain].append(record)

		# Make a new subalignment without it
		sansfirstrep = []
		sanslastrep = []
		sansfirstandlastrep = []
		sansfirstSeclastrep = []
		sansmisc = []
		for strain in straindic.keys():
			strainrepeats = straindic[strain]
			sanslastrep.extend( strainrepeats[:len(strainrepeats) - 1] )
			sansfirstrep.extend( strainrepeats[1:len(strainrepeats)] )
			sansfirstandlastrep.extend( strainrepeats[1:len(strainrepeats)-1] )
			sansfirstSeclastrep.extend( strainrepeats[2:len(strainrepeats)-1] )
			sansmisc.extend( strainrepeats[3:len(strainrepeats)-2] )

		sansFirst_alignment = MultipleSeqAlignment(sansfirstrep)
		ofile.write(f"{GENE}_sansFirstRepeat\tanserina\t{Psim(sansFirst_alignment)}\t{len(sansFirst_alignment)}\n")

		sansLastr_alignment = MultipleSeqAlignment(sanslastrep)
		ofile.write(f"{GENE}_sansLastRepeat\tanserina\t{Psim(sansLastr_alignment)}\t{len(sansLastr_alignment)}\n")

		sansFirstLastr_alignment = MultipleSeqAlignment(sansfirstandlastrep)
		ofile.write(f"{GENE}_sansFirstLastRepeat\tanserina\t{Psim(sansFirstLastr_alignment)}\t{len(sansFirstLastr_alignment)}\n")

		sansFirstSecLastr_alignment = MultipleSeqAlignment(sansfirstSeclastrep)
		ofile.write(f"{GENE}_sansFirstSecLastRepeat\tanserina\t{Psim(sansFirstSecLastr_alignment)}\t{len(sansFirstSecLastr_alignment)}\n")

		# sansmisc_alignment = MultipleSeqAlignment(sansmisc)
		# ofile.write(f"{GENE}_sansmisc\tanserina\t{Psim(sansmisc_alignment)}\t{len(sansmisc_alignment)}\n")


