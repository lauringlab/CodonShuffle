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
import pyfasta
import re

from random import shuffle,random,randint,choice
from collections import Counter
from os import system
from pandas import DataFrame
from pandas import *
from Bio import SeqIO
from ggplot import *
from subprocess import call
from Bio.SeqRecord import SeqRecord
from matplotlib.patches import Rectangle



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
#fold_mfe = 'seq_fold_mfe.txt'

#Run RNAfold
i = open(filename, "r")
o = open(foldname, "w")
call(["./lib/ViennaRNA-2.1.9/Progs/RNAfold", "--noPS"], stdin=i, stdout=o) #Insert path before RNAfold in this line (MFOLD)
i.close
o.close

#Read fold file
fold=open(foldname,"rU").readlines()
#fold_mfe = open(fold_mfe, "w")
fold_table = pandas.DataFrame()
matches = []
pattern = re.compile(r'.*\(\s{0,2}(\-?\d{1,3}\.\d{1,3})')
for line in fold:
    matches += pattern.findall(line)
    #    fold_mfe.write(str(matches))
#fold_mfe.close
fold_table.insert(0, 'mfe', matches)
fold_table['mfe'] = fold_table['mfe'].astype(float)
fold_mean = fold_table['mfe'].mean()
fold_std = fold_table['mfe'].std()



matches_seq_pos = []
pattern_seq_pos = re.compile(r'^\>.*_(\d+)$')
for lines in fold:
    matches_seq_pos += pattern_seq_pos.findall(lines)

fold_table.insert(0, 'Position', matches_seq_pos)
fold_table['Position'] = fold_table['Position'].astype(float)






#Make a graph 
#Ggplot
fold_graphname = os.path.splitext(args.filename)[0]+'_'+str(wanted_id[0])+'_fold.pdf'
fold_graph = ggplot(aes(x='Position', y='mfe'), fold_table)+geom_bar(stat='identity')+labs("Nucleotide position","Minimum free energy")+ggtitle('RNA sliding window \n Mean='+str(round(fold_mean, 2))+'\n Standard Deviation ='+str(round(fold_std, 2)))+theme_seaborn()
ggsave(fold_graph, fold_graphname)    





# #Matplotlib
# extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.bar(fold_table['Position'], fold_table['mfe'])
# ax.set_title('RNA sliding window')
# ax.set_xlabel('Nucleotide position')
# ax.set_ylabel('Minimum free energy')
# ax.legend([extra, extra], ('Mean='+str(fold_mean), 'Standard Deviation='+str(fold_std)), loc=4)
# #plt.show()
# #ax.text(right, top, 'right bottom', horizontalalignment='right', verticalalignment='bottom')
# fig.savefig('seq_file_fold.pdf')
# plt.close(fig)
# plt.show()











