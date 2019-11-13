#!/usr/bin/env python

import argparse
import sys, os, shutil
import subprocess
import csv
import numpy as np
from Bio import SeqIO


# Command line options
parser = argparse.ArgumentParser(description='This script takes two fasta files and an orthologues .tsv file from orthofinder and will write a fasta for each species that contains only the single orthologues (this includes the orthologues identified in orthologous groups as opposed to just the single copy orthologies)')
parser.add_argument("-r1","--RSEM1",
					type=str,
					help="Fasta with sequences from the first individual (should be the species in the first column of the .tsv file)")
parser.add_argument("-r2","--RSEM2",
					type=str,
					help="Fasta with sequences from the second individual (should be the species in the second column of the .tsv file)")
parser.add_argument("-t","--tsv",
					type=str,
					help="Orhtologues .tsv file")
parser.add_argument("-nd","--namedel",
					type=str,
                    default=" ",
                    help="if there is part of the name that needs removed from the end of each sequence (e.g. _AA)")
args=parser.parse_args()


all_orthologs = []
with open(args.tsv) as OF:
#with open('Bothriechis_nubestris__v__Bothriechis_nigroviridis.tsv') as OF:
	reader = csv.reader(OF, delimiter='\t')
	for row in reader:
		all_orthologs.append(row)

name1=all_orthologs[0][1]
name2=all_orthologs[0][2]

all_orthologs = all_orthologs[1:]


taxon1 = []
with open(args.RSEM1) as OF:
#with open('Bnubestris_RSEM.csv') as OF:
	reader = csv.reader(OF, delimiter=',')
	for row in reader:
		taxon1.append(row)





taxon2 = []
with open(args.RSEM2) as OF:
#with open('Bnigroviridis_RSEM.csv') as OF:
	reader = csv.reader(OF, delimiter=',')
	for row in reader:
		taxon2.append(row)



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


## Identify orthologues in sheet1
orthos_sect = []
header = []
header.extend(taxon1[0])
header.insert(1,taxon2[0][0])
header.extend(taxon2[0][3:len(row)])
orthos_sect.append(header)
for ortho in strict:
	for row in taxon1:
		if ortho[1] in row[0] or ortho[2] in row[0]:
			print("got1")
			ortho_row = []
			ortho_row.extend(row)
			len(ortho_row)
			for row_2 in taxon2:
				if ortho[1] in row_2[0] or ortho[2] in row_2[0]:
					ortho_row.insert(1,row_2[0])
					ortho_row.extend(row_2[3:len(row)])					
	orthos_sect.append(ortho_row)
					
print(len(orthos_sect))

print(len(strict))
					
## remove orthologs			
for ortho in strict:
	for row in taxon1:
		if ortho[1] in row[0] or ortho[2] in row[0]:
			del taxon1[taxon1.index(row)]
			

for ortho in strict:	
	for row in taxon2:
		if ortho[1] in row[0] or ortho[2] in row[0]:
			del taxon2[taxon2.index(row)]


diff_1 = (len(orthos_sect[0])-len(taxon1[0]) - 1)
for row in taxon1:
	row.insert(1,'Paralogue')
	row.extend('0'* (diff_1))


diff_2 = len(orthos_sect[0])-diff_1-4
for row in taxon2:
	row.insert(0,'Paralogue')
	for i in range(4,(4+diff_2)):
		row.insert(i,'0')

final_set = []
for row in orthos_sect:
	final_set.append(row)
	
for row in taxon1:
	final_set.append(row)

for row in taxon2:
	final_set.append(row)


with open('Combined_RSEM.csv', 'w') as csv_file:
	csv_writer = csv.writer(csv_file, delimiter = ',')
	for row in final_set:
		csv_writer.writerow(row)
		
csv_file.close()
	

