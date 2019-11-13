

import argparse
from Bio import AlignIO

parser = argparse.ArgumentParser(description='This script takes two fasta files and an orthologues .tsv file from orthofinder and will write a fasta for each species that contains only the single orthologues (this includes the orthologues identified in orthologous groups as opposed to just the single copy orthologies)')
parser.add_argument("-f","--fasta",
					type=str,
					help="Fasta alignment to trim")
args=parser.parse_args()

name = args.fasta
name = name.split('.')
name = name[0]

alignment = AlignIO.read(args.fasta,"fasta")
for seq in alignment:
	seq.seq = seq.seq[0:-3]

outfile = name + "_no_stop.fasta"
handle = open(outfile, "w")
AlignIO.write(alignment, outfile, "fasta")
handle.close()

