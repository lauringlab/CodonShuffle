#!/usr/bin/python2.7

from __future__ import division
import os
import numpy
import scipy
import matplotlib
import pandas
import statsmodels
import patsy
import sys
import argparse
import matplotlib.pyplot as plt

from random import shuffle,random,randint,choice
from collections import Counter
from os import system
from pandas import DataFrame
from pandas import *
from Bio import SeqIO
#from ggplot import *
from subprocess import call
from Bio.SeqRecord import SeqRecord




#Input Data
parser = argparse.ArgumentParser(description='RNA sliding window.')
parser.add_argument('-i', nargs='?', help='Input Filename', required=True, dest="filename")
parser.add_argument('-s', nargs='*', help='Number of the sequence in the fasta file', dest="seq_number")
args = parser.parse_args()

inseq_file=open(args.filename,'rU')
wanted_id=args.seq_number

seq_records=[]
wanted_seq=[]


print wanted_id

seq_file = open('seq_file.fas', 'w') 
for seq_record in SeqIO.parse(inseq_file, "fasta"):
    seq_records.append(seq_record)
    if seq_record.id in wanted_id:
        SeqIO.write([seq_record], seq_file, "fasta")
seq_file.close()

#Split the sequence 
cmd = "pyfasta split -k 100 -o 80 seq_file.fas -n 1"
os.system(cmd)

filename = 'seq_file.split.100mer.80overlap.fas'
foldname = 'seq_file.fold'

i = open(filename, "r")
o = open(foldname, "w")

call(["./lib/ViennaRNA-2.1.9/Progs/RNAfold", "--noPS"], stdin=i, stdout=o) #Insert path before RNAfold in this line (MFOLD)

i.close
o.close

