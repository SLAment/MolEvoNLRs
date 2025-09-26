#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function

# ================== REPEATescriptor =================
# Script to characterize the repeat domain of the NLR genes in the Podospora species complex. 
# Inspired on HetRdescriptor.py

# The input is a special alignment of nucleotides but aligned to preserve the
# coding frame of the repeats.

# ==================================================
# Sandra Lorena Ament-Velasquez
# 2023/02/22
# +++++++++++++++++++++++++++++++++++++++++++++++++
# 2024/03/19 - v1.32, fixed an issue with the repeat naming for unknown repeats ("?")
# 2024/03/19 - v1.31, finally removed code from original Chevanne et al. 2010 classification
# 2024/02/28 - v1.3, added recognition of nwd1
# ------------------------------------------------------
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

import sys # For reading the input
import re
import os # For the input name
from collections import defaultdict  # This tells python that the dictionary contains a list so you can freely append things to the value of each key
import argparse # For the fancy options
# ------------------------------------------------------

version = 2.1
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Characterize the REPEAT diversity of protein alignment *", epilog="Warning: the script doesn't work properly if the input sequences have tracks of IUPAC symbols at the beginning and end of the sequences, or '?' symbols anywhere.") # Create the object using class argparse

# Mandatory options
parser.add_argument('fasta', help="Nucleotide alignment of an NLR, but with gaps introduced to maintain frame in case of introns or indels")
# Extra options
parser.add_argument("--output", "-o", help="Base of output name (nothing by default, so reporting to working directory)", default="alignment")
parser.add_argument("--aminoacids", "-a", help="String of comma-separated amino acid residues to classify the REPEAT repeats, base 1. Default: 5,6,7,9,25,27", default = "5,6,7,9,25,27", type=str)
parser.add_argument("--regex", "-r", help="Regex used to find the repeats. Default: regex for ANK repeats of Pa_2_8180.", default = r"([NGRWDSYLQFI]{2,5})(TAL)(Q|E)([ATVG]{3})([\w]{2})(G|S)(H|Y)(L)(E|K|R)([VI]{2})([EK]{2})([L]{2})([VAT]{2})([GATDVNI]{4})([\w]{3,6})", type=str)
parser.add_argument("--replen", "-l", help="Expected typical repeat length (for classification, not repeat detection). Default: 33", default = 33, type=int)
parser.add_argument("--repminlen", "-m", help="Expected minimum repeat length (for classification, not repeat detection). Default: replen", type=int)
parser.add_argument("--repmaxlen", "-M", help="Expected maximum repeat length (for classification, not repeat detection). Default: replen", type=int)

# Input from console
try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	fastaopen = open(args.fasta, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------

# ---------------------------------
# Get sequences in file
# ---------------------------------
# Function to translate nucleotide sequence considering gaps and frameshifts
# ThREPEATs ChatGPT
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

seqdic = {}

# Older version
for seq_record in SeqIO.parse(fastaopen, "fasta"):
	geneseq = seq_record

	# Remove the spaces I have for me to see the WDs better, and the actual
	# frame shifts make them ambiguous bases so that the translation works 
	# Indels are also removed
	ncseq = str(geneseq.seq).upper().strip("-").replace("---","")

	# Make a new record
	geneseq.seq = Seq(ncseq)

	# Finally translate
	aaseq = custom_translate(geneseq.seq)

	# Add to the working dictionary
	seqdic[seq_record.id] = (ncseq, aaseq)

# for seq in seqdic.keys(): print(seq, seqdic[seq])
# print(seqdic)

# ---------------------------------
# Functions
# ---------------------------------
seqidregex = re.compile(r"([h]?nwd[a-zA-Z0-9\-]{0,5}|het-.{1,2})(_{1,2})([a-zA-Z0-9\.]+)_(.*)")

def simplifynames(name):
	namematch = seqidregex.search(name.replace("New|", ""))
	if namematch:
		newname = f"{namematch.group(1).rstrip('_')}__{namematch.group(3)}"
	else:
		newname = name
	return(newname)

# Variables to keep track of numbers and names
REPEATid_dic = {}

REPEATcount4names = 0
REPEATcount4names_Pnt1 = 0
REPEATcount4names_73550 = 0

def assignREPEATid(REPEATrepeat, baseseqid):
	global REPEATid_dic
	global REPEATcount4names
	global REPEATcount4names_Pnt1
	global REPEATcount4names_73550

	aminoacids = args.aminoacids.split(",")
	
	# Substract the chosen amino acids from the REPEAT repeat in a specific string
	specificstring = ''
	for amino in args.aminoacids.split(","):
		# print(int(amino), file=sys.stderr)
		specificstring += REPEATrepeat[int(amino) - 1]

	if specificstring not in REPEATid_dic.keys():
		if 'PaPnt1' in baseseqid:
			REPEATcount4names_Pnt1 += 1
			REPEATid_dic[specificstring] = "p" + str(REPEATcount4names_Pnt1)
		elif 'Pa_7_3550' in baseseqid:
			REPEATcount4names_73550 += 1
			REPEATid_dic[specificstring] = "q" + str(REPEATcount4names_73550)
		else:
			REPEATcount4names += 1
			REPEATid_dic[specificstring] = "a" + str(REPEATcount4names)

	return(REPEATid_dic[specificstring])

def assigngene(baseseqid):
	if "__" in baseseqid:
		gene, baseseqid = baseseqid.split("__")
	else:
		gene = "?"

	# if "Pa_" in gene:
	# 	allele = "NA"
	# if "Pa_6_7270" in baseseqid: # Here the gene and allele is the same, Pa_6_7270
	# 	allele = "Pa_6_7270"
	# 	gene = allele
	# else:
	# 	allele = baseseqid
		# gene = "?"
	
	# return([gene, allele])
	return gene

def assignspp(strainID):
	species = "?"
	if "Pa" in strainID:
		species = "anserina"
	elif "Pc" in strainID:
		species = "comata"
	elif strainID in ["CmEmm.R10", "CoEcp.R10", "CoEcm.R10", "CoEfp.R10", "ChEhDap.R10", "ChEhDa.R10", "CaDam.R10", "CsDfp.R10", "CsDfm.R10"]:
		species = "anserina"
	elif "CBS333.63" in strainID or "CBS451.62" in strainID or "CBS237.71" in strainID:
		species = "pauciseta"
	elif "CBS411.78" in strainID:
		species = "pseudoanserina"
	elif "CBS415.72" in strainID:
		species = "pseudocomata"
	elif "CBS112042" in strainID:
		species = "bellae-mahoneyi"
	elif "CBS124.78" in strainID or "CBS253.71" in strainID:
		species = "pseudopauciseta"
	return(species)


# ---------------------------------
# Find the repeats
# ---------------------------------

# Output files
outputfasta = open( f"{args.output}.fa", "w")
outputmetadata = open(f"{args.output}_metadata.txt", "w")

outputmetadata.write("Sequence\tRepeatType\tStrain\tSpecies\tGene\tSeqID\tPosition\n") # Print a header for the metadata

# REGEX and expected length to find REPEAT repeats
StrictREPEATmotif = re.compile(r'{}'.format(args.regex))
EXPlenREPEAT = args.replen

# The range of minimum an maximum sizes impacts the classification
if args.repminlen:
	minlenREPEAT = args.repminlen
else:
	minlenREPEAT = EXPlenREPEAT

if args.repmaxlen: # The range() function I use use this as a non-inclusive limit, so add 1
	maxlenREPEAT = args.repmaxlen + 1
else:
	maxlenREPEAT = EXPlenREPEAT + 1


# Regex for species and strain classification for Podospora spp
strainregex = re.compile(r"(Pa|Pc|CBS)([a-zA-Z0-9\.]+)")
strainR10regex = re.compile(r"([a-zA-Z]*)(p|m)(.R10)")

# Keep track of the origin of everything
generalcount = 0

# Header of results
sys.stdout.write(f'Seqid\tRPdomainLen_aa\tExpected_no_repeats\tCount_full_repeats\tClassification\tCterm\n') 

for seqid in seqdic.keys():
	baseseqid = simplifynames(seqid)
	REPEATseq = seqdic[seqid][1]

	# Loop through the expected length of the REPEAT C-term:
	chippingREPEATseq = REPEATseq
	chippingREPEATseqNC = seqdic[seqid][0]

	countREPEAT = 0
	REPEATstring = ''
	badREPEATstring = ''
	RPdomainStart = 0 
	for n in range(0, len(REPEATseq)): # Loop the whole sequence (a bit slow if there are sections that are not repeats)
		REPEATmatch = StrictREPEATmotif.search(chippingREPEATseq) # Look for a specific shape of repeat as defined above
		if REPEATmatch:
			currentRP = REPEATmatch.group(0)
			rpstart, rpend = REPEATmatch.span()

			if 'X' not in currentRP:
				if rpstart > 0: # something is wrong with that repeat
					if len(currentRP) in list(range(minlenREPEAT, maxlenREPEAT)): # Tolerate indels to consider it a complete repeat
						rpid = assignREPEATid(currentRP, seqid)
						if REPEATstring == '': # This is the first encountered repeat
							REPEATstring += rpid + ','
							RPdomainStart = rpstart # There is seq before the repeat domain in this alignment
						else:
							REPEATstring += "?," + rpid + ','
					else: # Probably the repeat is broken
						REPEATstring += "?,"
						rpid = "?"
				else: # The Repeat is well-behaved
					rpid = assignREPEATid(currentRP, seqid)
					REPEATstring += rpid + ','
			else: # There is a stop codon so we can't tedermine what repeat this is
				REPEATstring += "!,"
				rpid = "!"

			# Retrieve sequence into the fasta file
			generalcount += 1 # To make every single repeat unique
			currentRPnc = chippingREPEATseqNC[rpstart * 3 : rpend * 3] # KEY assumption: the base nucleotide alignment is made in such a way that it conserves the triplets
			outputfasta.write(f">{rpid}_{baseseqid}_{generalcount}\n")
			outputfasta.write(f"{currentRPnc}\n")

			# Get the metadata
			strainID = "" # Re-start so it doesn't carry from previous loop
			species = ""

			## Infer the strain and the species
			if "Podan" in baseseqid or "CAL30215.1" in baseseqid:
				strainID = "PaSp"
				species = "anserina"
			elif "PODCO" in baseseqid:
				strainID = "PcTdp"
				species = "comata"
			elif "FJ897789" in baseseqid:
				strainID = "PaA"
				species = "anserina"
			else:
				if "__" in baseseqid: # it follows the proper format
					gene, actualbaseseqid = baseseqid.split("__")
					strainmatch = strainregex.search(actualbaseseqid)
					strainR10match = strainR10regex.search(actualbaseseqid)
				else:
					strainmatch = strainregex.search(baseseqid)
					strainR10match = strainR10regex.search(baseseqid)

				if strainR10match:
					strainID = strainR10match.group(1) + strainR10match.group(2)
					species = "anserina"
				elif strainmatch:
					strainID = strainmatch.group(1) + strainmatch.group(2)
					species = assignspp(strainID)
				else:
					species = "Podospora"

			outputmetadata.write(f"{rpid}_{baseseqid}_{generalcount}\t{rpid}\t{strainID}\t{species}\t{assigngene(baseseqid)}\t{generalcount}\t{n+1}\n")
			# outputmetadata.write(f"{rpid}_{baseseqid}_{generalcount}\t{rpid}\t{strainID}\t{species}\t{assigngene(baseseqid)[0]}\t{assigngene(baseseqid)[1]}\t{generalcount}\t{n+1}\n")

			# Move to the next section of the sequence
			chippingREPEATseq = chippingREPEATseq[rpend:] # Move to the next one
			chippingREPEATseqNC = chippingREPEATseqNC[rpend * 3:] # Move to the next one

			countREPEAT += 1 # Count this repeat as whole, for reporting
			currentRP = "" # Re-start so it doesn't carry from previous loop
			strainID = "" # Re-start so it doesn't carry from previous loop
		# else: # Hopefully just the Cterm
		# 	print(chippingREPEATseq)


	# Put something if the REPEATstring is empty
	if REPEATstring == '': REPEATstring = '.'

	# Report
	lenREPEATdomain =  len(REPEATseq) - RPdomainStart - len(chippingREPEATseq)
	resultline = seqid + '\t' + str(lenREPEATdomain) + '\t' + format(lenREPEATdomain/EXPlenREPEAT, '0.2f') + '\t' + str(countREPEAT) + '\t' + REPEATstring.rstrip(',') + '\t' + chippingREPEATseq + '\n'
	
	sys.stdout.write(resultline) 

outputfasta.close()
outputmetadata.close()

# Reporting
print("-----------------------------------")
print(f"Amino acids used for classification: {args.aminoacids}")
print("Repeat number of Unknown:", REPEATcount4names)
print("Total", len(list(REPEATid_dic.values())), '\n')

print(REPEATid_dic)
