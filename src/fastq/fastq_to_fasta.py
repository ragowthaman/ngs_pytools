#!/Users/gramasamy/Software/virtualenvs/ngs_pytools/bin/python
'''
    convert fastq to fasta format
    input:
        fastq file to be converted
    author: ragowthaman@gmail.com
'''

import os
import datetime
import re
from Bio import SeqIO

from datetime import date
today = date.today()

####################################################################
#
# Functions
#
####################################################################
def argparse():
    '''convert fastq to fasta format'''
    import argparse
    parser = argparse.ArgumentParser(description='filter fastq file by length of the reads')
    parser.add_argument('--fastq', required=True, help='input Fastq file', type=str)
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

fastq = args.fastq

outformat = 'fasta'
output_filename = fastq.rstrip('.fastq') + '.fasta'

# log lists
filtered_seq_records = []

# read in
seq_records = SeqIO.parse(fastq, "fastq")
SeqIO.write(seq_records, output_filename, "fasta")



# # loop thru fasta; rename it
# for seq_record in SeqIO.parse(fastq, "fastq"):
#     seq = seq_record.seq
#     seqId = seq_record.id
#     seqDescription = seq_record.description
#
#             filtered_seq_records.append(seq_record)
#     if args.idfile:
#         if seqId in id_dict:
#             filtered_seq_records.append(seq_record)
#
#
#
# # wirte output
# if not args.outfile:
#     output_filename = args.fastq + '_filtered_id'
# else:
#     output_filename = args.outfile
#
# output_FH = open(output_filename, "a")
# SeqIO.write(filtered_seq_records, output_FH, outformat)
