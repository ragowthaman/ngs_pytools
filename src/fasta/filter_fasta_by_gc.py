#!/Users/gramasamy/Software/virtualenvs/ngs_pytools/bin/python
'''
    filter fasta file by GC content of the sequences
    input:
        fasta file to be filtered
        min read length (inclusive); max read length (inclusive)
    author: ragowthaman@gmail.com
'''

import os
import datetime
import re
from Bio import SeqIO
from Bio.SeqUtils import GC


from datetime import date
today = date.today()

####################################################################
#
# Functions
#
####################################################################
def argparse():
    '''parses the command line arguments'''
    import argparse
    parser = argparse.ArgumentParser(description='filter fasta file by GC content of the reads')
    parser.add_argument('--fasta', required=True, help='input Fasta file', type=str)
    parser.add_argument('--minimum', required=True, help='minimum GC content to be included', type=int)
    parser.add_argument('--maximum', required=True, help='maximum GC content to be included', type=int)
    parser.add_argument('--outfile', required=False, help='name of the outfile', type=str)
    parser.add_argument('--debug', required=False, default=False, type=bool)
    parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    return args

####################################################################
#
# main
#
####################################################################
args = argparse()

fasta = args.fasta
minimum = args.minimum
maximum = args.maximum

# log lists
seq_lt_minimum = []
seq_mt_maximum = []
seq_in_range = []

seq_records_in_range = []
# loop thru fasta
for seq_record in SeqIO.parse(fasta, "fasta"):
    seq = seq_record.seq
    seqId = seq_record.id
    seqDescription = seq_record.description
    seqGC = GC(seq)

    if minimum <= seqGC <= maximum:
        seq_in_range.append(seqId)
        seq_records_in_range.append(seq_record)
    elif seqGC < minimum:
        seq_lt_minimum.append(seqId)
    elif seqGC > maximum:
        seq_mt_maximum.append(seqId)


# wirte output
if not args.outfile:
    output_filename = args.fasta + '_filtered_gc'
else:
    output_filename = args.outfile

output_FH = open(output_filename, "a")
SeqIO.write(seq_records_in_range, output_FH, "fasta")
