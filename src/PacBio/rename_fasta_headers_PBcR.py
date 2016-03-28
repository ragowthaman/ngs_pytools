#!/Users/gramasamy/Software/virtualenvs/ngs_pytools/bin/python
'''
    renames fasta headers per tab file from PBcR read correction
    PBcR pipeline renames the input pacbio read ids. This script takes
    the corrected fasta and renames using the information from log file

    input:
        fasta file to be renamed;
        tab file with two column of names. Column1: From Name; Column2: To Name; Column4,5 are appeneded at the end of "original" or "old" readid.sequences
    author: ragowthaman@gmail.com
'''

import os
import datetime
import re
from Bio import SeqIO
import collections

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
    parser = argparse.ArgumentParser(
        description='renames fasta headers per names in the tab file')
    parser.add_argument('--fasta', required=True, help='input Fasta file', type=str)
    parser.add_argument('--names', required=True, help='tab file with names <FROM> and <TO>', type=str)
    parser.add_argument('--reverse_name_order', required=False, default=False, type=bool, help='set this if first column is <TO> and second column is <FROM>')
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
names = args.names

# hasify names
names_dict = dict()
with open(names, mode='r') as f:
    for line in f.readlines():
        cols = line.split("\t")
        name_from = cols[0]
        name_to = cols[1]
        subread_start = cols[3]
        subread_end = cols[4]

        if args.reverse_name_order is True:
            name_from = cols[1]
            name_to = cols[0]
        name_from_suffix = name_from.split("_")
        name_to = name_to + '#' + subread_start + '_' + subread_end
        names_dict[name_from] = name_to

print names_dict


# log lists
fasta_notin_rename = []
renamed = []

renamed_seq_records = []
# loop thru fasta; rename it
for seq_record in SeqIO.parse(fasta, "fasta"):
    seq = seq_record.seq
    seqId = seq_record.id
    seqDescription = seq_record.description
    # print seqId
    seqId = re.sub(r'/\d+_\d+', '', seqId)
    # print seqId

    if seqId in names_dict:
        print seqId
        print "\tyes"
        new_seqId = names_dict[seqId]
        print "\t"+new_seqId
        seq_record.id = new_seqId
        renamed_seq_records.append(seq_record)
        renamed.append(seqId)
    else:
        fasta_notin_rename.append(seqId)

output_filename = args.fasta + '.renamed'
renamed_seq_record_filehandle = open(output_filename, "a")
SeqIO.write(renamed_seq_records, renamed_seq_record_filehandle, "fasta")
