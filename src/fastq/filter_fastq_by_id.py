#!/Users/gramasamy/Software/virtualenvs/ngs_pytools/bin/python
'''
    filter fastq file by id of the sequences
    input:
        fasta file to be filtered
        file of ids or single id to be pulled out
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
    '''parses the command line arguments'''
    import argparse
    parser = argparse.ArgumentParser(description='filter fastq file by length of the reads')
    parser.add_argument('--fastq', required=True, help='input Fastq file', type=str)
    parser.add_argument('--id', required=False, help='single id to be filtered', type=str)
    parser.add_argument('--idfile', required=False, help='file of ids to be filtered', type=str)
    parser.add_argument('--outfile', required=False, help='name of the outfile', type=str)
    parser.add_argument('--outformat', required=False, help='format for the output eg., fasta, fastq, gb', type=str)
    parser.add_argument('--debug', required=False, default=False, type=bool)
    parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    return args


def hasify_one_column_file(file):
    '''
    hasify a file with only one column; one word per row; newline should be the only space character
    each row will become keyword holding a value of int(1)
    '''

    dict = {}
    with open(file, 'r') as f:
        for line in f.readlines():
            dict[line] = 1

    return dict


####################################################################
#
# main
#
####################################################################
args = argparse()

fastq = args.fastq
if args.id:
    id = args.id
if args.idfile:
    idfile = args.idfile

outformat = args.outformat if args.outformat else 'fasta'

# log lists
filtered_seq_records = []

if args.idfile:
    id_dict = hasify_one_column_file(args.idfile)

# loop thru fasta; rename it
for seq_record in SeqIO.parse(fastq, "fastq"):
    seq = seq_record.seq
    seqId = seq_record.id
    seqDescription = seq_record.description
    seqLen = len(seq)

    if args.id:
        if args.id in seqId:
            filtered_seq_records.append(seq_record)

    if args.idfile:
        if seqId in id_dict:
            filtered_seq_records.append(seq_record)



# wirte output
if not args.outfile:
    output_filename = args.fastq + '_filtered_id'
else:
    output_filename = args.outfile

output_FH = open(output_filename, "a")
SeqIO.write(filtered_seq_records, output_FH, outformat)
