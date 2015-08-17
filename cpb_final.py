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
# import rpy2.robjects as robjects
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
# pandas2ri.activate()
from Bio import SeqIO
from ggplot import *
from subprocess import call
from Bio.SeqRecord import SeqRecord

#Constant values

#nt list
nts = ['A','C','G','T']

#Translation Table
tt = {"TTT":"F|Phe","TTC":"F|Phe","TTA":"L|Leu","TTG":"L|Leu","TCT":"S|Ser","TCC":"S|Ser","TCA":"S|Ser","TCG":"S|Ser", "TAT":"Y|Tyr","TAC":"Y|Tyr","TAA":"*|Stp","TAG":"*|Stp","TGT":"C|Cys","TGC":"C|Cys","TGA":"*|Stp","TGG":"W|Trp", "CTT":"L|Leu","CTC":"L|Leu","CTA":"L|Leu","CTG":"L|Leu","CCT":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro","CAT":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGT":"R|Arg","CGC":"R|Arg","CGA":"R|Arg","CGG":"R|Arg", "ATT":"I|Ile","ATC":"I|Ile","ATA":"I|Ile","ATG":"M|Met","ACT":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr", "AAT":"N|Asn","AAC":"N|Asn","AAA":"K|Lys","AAG":"K|Lys","AGT":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg","GTT":"V|Val","GTC":"V|Val","GTA":"V|Val","GTG":"V|Val","GCT":"A|Ala","GCC":"A|Ala","GCA":"A|Ala","GCG":"A|Ala", "GAT":"D|Asp","GAC":"D|Asp","GAA":"E|Glu","GAG":"E|Glu","GGT":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}
#Following functions or their combinations produce randomized or scrambled nucleotide sequence from input sequence.
#Amino-acid sequence of derived sequence is identical to the input sequence, but nucleotide composition (GC-, nucleotide, or dinucleotide content) may differ slightly for randomized sequences.

cpbname = 'seq.cpb'
cpb_file=open(cpbname,'w')

#Human CPS from Coleman et al 2008, pmid: 18583614
cps_human = pandas.read_csv("Coleman_CPS.csv", sep=';')
#Delete column 
cps_human = cps_human.drop('Aapair', 1)
cps_human = cps_human.drop('Expected', 1)
cps_human = cps_human.drop('Observed', 1)
cps_human = cps_human.drop('Observed/Expected', 1)
cps_human = cps_human.sort(['CodonPair'], ascending=[True])


for nuc_rec in SeqIO.parse("dn23.fas", "fasta"):
    dicodons = [str(nuc_rec.seq[i:i+6]) for i in range(0,len(nuc_rec.seq)-3, 3)]


    dicodon_counts = Counter(dicodons)

#    cps = []
    for cp in dicodons:
        cod1 = cp[:3]
        cod2 = cp[3:]
        if cod1 in ['TAG', 'TGA', 'TAA'] or cod2 in ['TAG', 'TGA', 'TAA']:
            continue #skip stop codons
        aa1 = tt[cod1].split('|')[0]
        aa2 = tt[cod2].split('|')[0]
        ap = aa1+aa2
        #score = np.log(dicodons.count(cp) / ( ((codons.count(cod1) * codons.count(cod2)) / (aas.count(aa1) * aas.count(aa2))) * diaas.count(ap)))
#        score = np.log(dicodon_counts[cp] / ( ((codon_counts[cod1] * codon_counts[cod2]) / (aa_counts[aa1] * aa_counts[aa2])) * diaa_counts[ap]))
#        cps.append(score)
        dicodon_df = pandas.DataFrame.from_dict(dicodon_counts, orient='index').reset_index()
        dicodon_df = dicodon_df.sort(['index'], ascending=[True])
        dicodon_df.columns = ['CodonPair', 'Obs']
        cps_tb_final = pandas.merge(cps_human, dicodon_df, on='CodonPair', how='inner')
        cps_tb_final['CPS'] = cps_tb_final['CPS'].replace({',':'.'}, regex=True)
        cps_tb_final['CPS'] = cps_tb_final['CPS'].astype(float)
        cps_tb_final['result'] = cps_tb_final.CPS * cps_tb_final.Obs
    cpb = sum(cps_tb_final['result'])/sum(cps_tb_final['Obs'])
    print str(cpb)
    cpb_file.write(str(cpb)+"\n")
cpb_file.close()

