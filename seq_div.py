#!/usr/bin/env python

import argparse
import sys, os, shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from sets import Set
import glob
import csv


parser = argparse.ArgumentParser(description='This script takes two fasta files and an orthologues .tsv file from orthofinder and will write a fasta for each species that contains only the single orthologues (this includes the orthologues identified in orthologous groups as opposed to just the single copy orthologies)')
parser.add_argument("-f1","--fasta1",
					type=argparse.FileType('r'),
					help="Fasta with sequences from the first individual (should be the species in the first column of the .tsv file)")
parser.add_argument("-f2","--fasta2",
					type=argparse.FileType('r'),
					help="Fasta with sequences from the second individual (should be the species in the second column of the .tsv file)")
args=parser.parse_args()


fasta1 = list(SeqIO.parse(args.fasta1,"fasta"))
#fasta1 = list(SeqIO.parse("../../../../../Bnube-CLUSTERED_CURATED.fasta","fasta"))

fasta2 = list(SeqIO.parse(args.fasta2,"fasta"))
#fasta2 = list(SeqIO.parse("../../../../../Bnigr-CLUSTERED_CURATED.fasta","fasta"))

if len(fasta1) != len(fasta2):
	print "Error: Sequence files are different lengths. Can't be orthologous databases"
	quit()

command = "mkdir pairwise_fastas"
subprocess.call(command,shell=True)

for seq1, seq2 in zip(fasta1, fasta2):
	seqs = [seq1, seq2]
	write_fasta = "pairwise_fastas/" + seq1.id + ".fasta"
	handle=open(write_fasta, "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(seqs)
	handle.close()


command = "mkdir alignments"
subprocess.call(command,shell=True)

cwd = os.getcwd()
files=sorted(glob.glob(cwd+'/pairwise_fastas/*'))

for file in files:
	print(file)
	name = file.split('/')
	name = name[-1]
	name = name.split('.')
	print(name)
	name = name[0] 
	command = "prank -d=" + file + " -o=alignments/" + name + "_codon-align -codon -F"
	subprocess.call(command,shell=True)


cwd = os.getcwd()
files=sorted(glob.glob(cwd+'/alignments/*'))


output = []
header = ['otholog1','ortholog2','#_conserved_sites','percent_identity']
output.append(header)
for file in files:
	alignment = AlignIO.read(file,"fasta")
	x = 0
	for site1, site2 in zip(alignment[0].seq, alignment[1].seq):
		if site1 == site2:
			x = x + 1
	prop = float(x) / float(len(alignment[0]))
	row = [alignment[0].id, alignment[1].id,x, prop]
	output.append(row)
	

with open('Distances.csv', mode='w') as csv_file:
	csv_writer = csv.writer(csv_file, delimiter=',')
	
	for line in output:
		csv_writer.writerow(line)		
	
