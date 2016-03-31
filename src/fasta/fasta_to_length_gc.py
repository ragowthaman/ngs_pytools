#!/Users/gramasamy/Software/virtualenvs/ngs_pytools/bin/python
'''
    find length and gc for the entries in fasta file
    input:
        fasta file
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
    parser.add_argument('--outfile', required=False, help='name of the outfile', type=str)
    parser.add_argument('--noplots', required=False, help='set not to make plots', type=bool)
    parser.add_argument('--plotly_username', required=False, help='Username for plotly', type=str)
    parser.add_argument('--plotly_passwd', required=False, help='password for plotly', type=str)
    parser.add_argument('--debug', required=False, default=False, type=bool)
    parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    return args

def make_histogram(list, title, filename):
    '''make histogram and save it locally'''
    import plotly.plotly as py
    import plotly.graph_objs as go
    import numpy as np

    py.sign_in(plotly_username, plotly_passwd)
    data = [ go.Histogram( x=list) ]
    layout = go.Layout(
            title=title,
            width=800, height=640,
            yaxis=dict(title='Count')
        )

    fig = go.Figure(data=data, layout=layout)
    py.image.save_as(fig, filename=filename)

####################################################################
#
# main
#
####################################################################
args = argparse()
plotly_username = args.plotly_username
plotly_passwd = args.plotly_passwd

fasta = args.fasta
lengths = []
gcs = []

infile_prefix = fasta.rstrip('.fasta').rstrip('.fsa')
outfile = infile_prefix + '_length_gc.tab'
if args.outfile:
    outfile = args.outfile
outfileH = open(outfile, 'w')

# loop thru fasta
for seq_record in SeqIO.parse(fasta, "fasta"):
    seq = seq_record.seq
    seqId = seq_record.id
    seqDescription = seq_record.description

    seqGC = GC(seq)
    gcs.append(seqGC)

    seqLen = len(seq)
    lengths.append(seqLen)


    outfileH.write("%s\t%s\t%s\n" % (seqId, seqLen, seqGC))

# make plots
if not args.noplots:
    plot_outfile_name = infile_prefix + '_length.png'
    make_histogram(lengths, "Fasta seq Lengths", plot_outfile_name)

    plot_outfile_name = infile_prefix + '_gc.png'
    make_histogram(gcs, "Fasta seq GC distribution", plot_outfile_name)
