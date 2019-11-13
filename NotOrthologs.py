#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(description='This script takes two fasta files, one that is all possible sequences for a species and one that contains the ')
parser.add_argument("-or","--orthologues",
					type=argparse.FileType('r'),
					help="Fasta with sequences from the first individual (should be the species in the first column of the .tsv file)")
parser.add_argument("-af","--allfasta",
					type=argparse.FileType('r'),
					help="Fasta with sequences from the second individual (should be the species in the second column of the .tsv file)")
parser.add_argument("-ot","--outfasta",
					type=str,
					help="output file (e.g. output.fasta and you have to include the .fasta) ")
args=parser.parse_args()

orthos = list(SeqIO.parse(args.orthologues,"fasta"))

all_fasta = list(SeqIO.parse(args.allfasta,"fasta"))

singles = []
for seq1 in all_fasta:
	x = 0
	for seq2 in orthos:
		if seq2.id == seq1.id:
			x = 1
			break
	if x == 0:
		singles.append(seq1)
	
write_fasta1 = args.outfasta
handle=open(write_fasta1, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(singles)
handle.close()


