#!/usr/bin/env python

import argparse
import sys, os, shutil
import subprocess
import csv
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
import time

# Command line options
parser = argparse.ArgumentParser(description='This script takes two fasta files and an orthologues .tsv file from orthofinder and will write a fasta for each species that contains only the single orthologues (this includes the orthologues identified in orthologous groups as opposed to just the single copy orthologies)')
parser.add_argument("-f1","--fasta1",
					type=argparse.FileType('r'),
					help="Fasta with sequences from the first individual (should be the species in the first column of the .tsv file)")
parser.add_argument("-f2","--fasta2",
					type=argparse.FileType('r'),
					help="Fasta with sequences from the second individual (should be the species in the second column of the .tsv file)")
parser.add_argument("-t","--tsv",
					type=str,
					help="Orhtologues .tsv file")
parser.add_argument("-nd","--namedel",
					type=str,
                    default=" ",
                    help="if there is part of the name that needs removed from the end of each sequence (e.g. _AA)")
args=parser.parse_args()


fasta1 = list(SeqIO.parse(args.fasta1,"fasta"))
#fasta1 = list(SeqIO.parse("../../../../../Bnube-CLUSTERED_CURATED.fasta","fasta"))

fasta2 = list(SeqIO.parse(args.fasta2,"fasta"))
#fasta2 = list(SeqIO.parse("../../../../../Bnigr-CLUSTERED_CURATED.fasta","fasta"))

## Longest substring problem, taken from stack exchange user thefourtheye on page https://stackoverflow.com/questions/18715688/find-common-substring-between-two-strings
def longestSubstringFinder(string1, string2):
    answer = ""
    len1, len2 = len(string1), len(string2)
    for i in range(len1):
        match = ""
        for j in range(len2):
            if (i + j < len1 and string1[i + j] == string2[j]):
                match += string2[j]
            else:
                if (len(match) > len(answer)): answer = match
                match = ""
    return answer


all_orthologs = []
with open(args.tsv) as OF:
#with open('Bothriechis_nubestris__v__Bothriechis_nigroviridis.tsv') as OF:
	reader = csv.reader(OF, delimiter='\t')
	for row in reader:
		all_orthologs.append(row)

name1=all_orthologs[0][1]
name2=all_orthologs[0][2]

all_orthologs = all_orthologs[1:]

## Identify 
bad_lines = []
for ortho1 in all_orthologs:
	for ortho2 in all_orthologs:
		if ortho1[1] == ortho2[1] and ortho1 != ortho2:
			bad_lines.append(ortho1)
			bad_lines.append(ortho2)
		elif ortho1[2] == ortho2[2] and ortho1 != ortho2:
			bad_lines.append(ortho1)
			bad_lines.append(ortho2)

all_orthologs = [x for x in all_orthologs if x not in bad_lines]  

## Identify and eliminate rows that have multiple matches (based on commas)
strict=[]
for ortho in all_orthologs:
	if ',' not in ortho[1] and ',' not in ortho[2]:
		strict.append(ortho)


## Remove Amino Acid designators from names
for row in strict:
	#row[1] = row[1].split("_AA")
	row[1] = row[1].split(args.namedel)
	row[1] = row[1][0]
	
for row in strict:
	#row[2] = row[2].split("_AA")
	row[2] = row[2].split(args.namedel)
	row[2] = row[2][0]




def MakeConcensusAlignment(pair):
	write_tmpfasta = "tmp.fasta"
	handle=open(write_tmpfasta, "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(pair)
	handle1.close
	
	while os.path.exists('tmp.fasta') == False:
		time.sleep(1)	
	command = "mafft " + write_tmpfasta +" > tmp_2.fasta"
	print(command)
	subprocess.call(command,shell=True)
		
	return()
		

## Make list of sequences for each of the fastas. Also write fasta of clustered sequences
out_fasta1 = []
out_fasta2 = []
y=0
tmp_list = []
for ortho in strict:
	pair = []
	for seq in fasta1:
		if ortho[1] in seq.id:
			out_fasta1.append(seq)
			pair.append(seq)
	for seq in fasta2:
		if ortho[2] in seq.id:
			out_fasta2.append(seq)
			pair.append(seq)
	write_tmpfasta = "tmp.fasta" + str(y)
	tmp_list.append(write_tmpfasta)
	handle1=open(write_tmpfasta, "w")
	writer1 = FastaIO.FastaWriter(handle1, wrap=None)
	writer1.write_file(pair)
	handle1.close()
	y = y + 1
	
seqids = []
x=1
consensuses = []
for tmp in tmp_list:
	command = "mafft --auto " + tmp +" > tmp_2.fasta"
	print(command)
	subprocess.call(command,shell=True)
	command = "rm " + tmp
	subprocess.call(command,shell=True)
	alignment = AlignIO.read("tmp_2.fasta","fasta")
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = summary_align.dumb_consensus(ambiguous='n')
	seqid = longestSubstringFinder(alignment[0].id,alignment[1].id)
	for seqid1 in seqids:
		if seqid == seqid1:
				seqid = seqid + str(x)
				x = x + 1
	seqids.append(seqid)			
	thing = name1 + "_" + name2 + "_" + seqid
	record = SeqRecord(consensus, id = thing, name = thing, description = thing)
	consensuses.append(record)

write_fasta = name1 + "_" + name2 + "_CLUSTERED_ORTHO.fasta"
SeqIO.write(consensuses, write_fasta, "fasta-2line")

	
	
	#subprocess.call(command,shell=True)
	#command = "cat tmp_cluster.fasta >> " + name1 + "_" + name2 + "_CLUSTERED_ORTHO.fasta"
	#subprocess.call(command,shell=True)
	
#command = "rm tmp.fasta tmp_cluster.fasta tmp_cluster.fasta.clstr" 
#subprocess.call(command,shell=True)

			
## Write orthofastas	
write_fasta1 = name1 + "-ortho-only.fasta"
#write_fasta1 = "Bnube-ortho-only.fasta"
handle=open(write_fasta1, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(out_fasta1)
handle.close()
	
write_fasta2 = name2 + "-ortho-only.fasta"
#write_fasta2 = "Bnigr-ortho-only.fasta"
handle=open(write_fasta2, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(out_fasta2)
handle.close()
