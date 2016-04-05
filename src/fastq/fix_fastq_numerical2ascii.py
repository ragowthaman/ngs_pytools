#!/Users/gramasamy/Software/virtualenvs/ngs_pytools/bin/python
'''
    convert fastq quality line with numerical values into ascii values (assuming its Phred quality)
    input:
        fastq file to be fixed
    author: ragowthaman@gmail.com
'''

import os
import datetime
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


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
    parser = argparse.ArgumentParser(description='convert fastq quality line with numerical values into ascii values (assuming its Phred quality)')
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
line_count = 0

with open(fastq, 'r') as f:
    for line in f:
        line_count += 1
        if line_count == 1:
            id = line.lstrip('@')
        if line_count == 2:
            seq = line.rstrip("\n")
        if line_count == 4:
            line_count = 0
            line = line.rstrip('\n').rstrip(" ").lstrip(" ")
            numerical_quality_values = []
            for item in line.split(" "):
                numerical_quality_values.append(float(item))

            new_record = SeqRecord(Seq(seq, generic_dna), id=id, description="")
            new_record.letter_annotations["phred_quality"] = numerical_quality_values
            print new_record.format("fastq")