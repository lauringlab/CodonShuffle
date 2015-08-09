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

# utils = importr("utils")
# plyr = importr("plyr")
# seqinr = importr("seqinr")

#Constant values

#nt list
nts = ['A','C','G','T']

#Translation Table
tt = {"TTT":"F|Phe","TTC":"F|Phe","TTA":"L|Leu","TTG":"L|Leu","TCT":"S|Ser","TCC":"S|Ser","TCA":"S|Ser","TCG":"S|Ser", "TAT":"Y|Tyr","TAC":"Y|Tyr","TAA":"*|Stp","TAG":"*|Stp","TGT":"C|Cys","TGC":"C|Cys","TGA":"*|Stp","TGG":"W|Trp", "CTT":"L|Leu","CTC":"L|Leu","CTA":"L|Leu","CTG":"L|Leu","CCT":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro","CAT":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGT":"R|Arg","CGC":"R|Arg","CGA":"R|Arg","CGG":"R|Arg", "ATT":"I|Ile","ATC":"I|Ile","ATA":"I|Ile","ATG":"M|Met","ACT":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr", "AAT":"N|Asn","AAC":"N|Asn","AAA":"K|Lys","AAG":"K|Lys","AGT":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg","GTT":"V|Val","GTC":"V|Val","GTA":"V|Val","GTG":"V|Val","GCT":"A|Ala","GCC":"A|Ala","GCA":"A|Ala","GCG":"A|Ala", "GAT":"D|Asp","GAC":"D|Asp","GAA":"E|Glu","GAG":"E|Glu","GGT":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}
#Following functions or their combinations produce randomized or scrambled nucleotide sequence from input sequence.
#Amino-acid sequence of derived sequence is identical to the input sequence, but nucleotide composition (GC-, nucleotide, or dinucleotide content) may differ slightly for randomized sequences.

def gc3(seq):#this function creates sequence with GC-content close to the GC-content of the input sequence, but counts of nucleotides may differ from input sequence.
    gc=at=0
    for num in xrange(2,len(seq),3):#first calculating A+T and G+C of the input sequence in third codon position
        if seq[num]=='A'or seq[num]=='T':
            at+=1
        elif seq[num]=='G'or seq[num]=='C':
            gc+=1
    at=at/(len(seq)/3.)
    gc=gc/(len(seq)/3.)
    
    seq1=[]
    for num in xrange(2,len(seq),3):#list "seq1" will contain the first two nt of codon, third codon position will containe flags for subsequent randomization. Flags ('_Y_','_R_','_H_',or '_N_') correspond to IUPAC single-letter code, Y-Pyrimindine(C or T), R-Purine(A or G), H-Not G(A or C or T), N-any.
        seq1+=seq[num-2:num],
        if (seq[num]=='T'or seq[num]=='C')and(seq[num-2:num]=='TT'or seq[num-2:num]=='TA'or seq[num-2:num]=='TG'or seq[num-2:num]=='CA'or seq[num-2:num]=='AA'or seq[num-2:num]=='AG'or seq[num-2:num]=='GA'):
            seq1+='_Y_',
        elif (seq[num]=='A'or seq[num]=='G')and(seq[num-2:num]=='TT'or seq[num-2:num]=='CA'or seq[num-2:num]=='AA'or seq[num-2:num]=='AG'or seq[num-2:num]=='GA'):
            seq1+='_R_',
        elif seq[num-2:num+1]=='ATT'or seq[num-2:num+1]=='ATC'or seq[num-2:num+1]=='ATA':
            seq1+='_H_',
        elif (seq[num]=='A'or seq[num]=='G'or seq[num]=='T'or seq[num]=='C')and(seq[num-2:num]=='TC'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CG'or seq[num-2:num]=='AC'or seq[num-2:num]=='GT'or seq[num-2:num]=='GC'or seq[num-2:num]=='GG'):
            seq1+='_N_',
        else: seq1+=seq[num],
    seq2=''#"seq2" will contain the derived sequence, approproate nucleotide is chosen for flags in "seq1", according to GC-content
    for i in seq1:
        if i == '_Y_':
            x=random()
            if x<=gc:
                seq2+='C'
            elif gc<x<=gc+at:
                seq2+='T'
            else: seq2+=choice('TC')
        elif i == '_R_':
            x=random()
            if x<=gc:
                seq2+='G'
            elif gc<x<=gc+at:
                seq2+='A'
            else: seq2+=choice('AG')
        elif i == '_H_':
            x=random()
            if x<=gc:
                seq2+='C'
            elif gc<x<=gc+at:
                seq2+=choice('AT')
            else: seq2+=choice('ATC')
        elif i == '_N_':
            x=random()
            if x<=gc:
                seq2+=choice('GC')
            elif gc<x<=gc+at:
                seq2+=choice('AT')
            else: seq2+=choice('AGTC')
        else: seq2+=i
    seq=seq2
    return seq   


def third_simple(seq):#this function creates scrambled sequence with the numbers of each nucleotide identical to the input sequence.
    Y=[]
    seq1=[]
    for num in xrange(2,len(seq),3):
        if (seq[num]=='T' or seq[num]=='C')and(seq[num-2:num]=='TT'or seq[num-2:num]=='TC'or seq[num-2:num]=='TA'or seq[num-2:num]=='TG'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CA'or seq[num-2:num]=='CG'or seq[num-2:num]=='AT'or seq[num-2:num]=='AC'or seq[num-2:num]=='AA'or seq[num-2:num]=='AG'or seq[num-2:num]=='GU'or seq[num-2:num]=='GC'or seq[num-2:num]=='GA'or seq[num-2:num]=='GG'):
            Y+=seq[num],
            seq1+=seq[num-2:num],'_Y_',
        else:seq1+=seq[num-2:num+1],
    #now "seq1" contains flag '_Y_' in the third position of all codons, where C->T or T->C shuffling preserves the aminoacid sequence (i.e. PHE, SER etc.).
    #C and T from the original sequence in this case would be extracted into list "Y"
    shuffle(Y)#shuffling of list "Y". For example, before shuffling "Y" is ['C','T','C']; after - ['T','C','C']or['C','C','T']or['C','T','C']
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='_Y_':seq2+=Y.pop(0)#now elements of "Y" are inserted back into the sequence instead of '_Y_', but in a different order compared to the input sequence
        else:seq2+=seq1[i]
    seq=seq2

    R=[]#similar to the previous step, but A and G are shuffled
    seq1=[]
    for num in xrange(2,len(seq),3):
        if (seq[num]=='A' or seq[num]=='G')and(seq[num-2:num]=='TT'or seq[num-2:num]=='TC'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CA'or seq[num-2:num]=='CG'or seq[num-2:num]=='AC'or seq[num-2:num]=='AA'or seq[num-2:num]=='AG'or seq[num-2:num]=='GT'or seq[num-2:num]=='GC'or seq[num-2:num]=='GA'or seq[num-2:num]=='GG'):
            R+=seq[num],
            seq1+=seq[num-2:num],'_R_',
        else:seq1+=seq[num-2:num+1],
    shuffle(R)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='_R_':seq2+=R.pop(0)
        else:seq2+=seq1[i]
    seq=seq2
    

    H=[]#similar to the previous step, but A,C, and T are shuffled. Affected aminoacids are ILE (three codons), four-codon and four-codon portion of six-codon aminoacids.
    seq1=[]
    for num in xrange(2,len(seq),3):
        if (seq[num]=='A'or seq[num]=='C'or seq[num]=='T')and(seq[num-2:num]=='TC'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CG'or seq[num-2:num]=='AT'or seq[num-2:num]=='AC'or seq[num-2:num]=='GT'or seq[num-2:num]=='GC'or seq[num-2:num]=='GG'):
            H+=seq[num],
            seq1+=seq[num-2:num],'_H_',
        else:seq1+=seq[num-2:num+1],
    shuffle(H)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='_H_':seq2+=H.pop(0)
        else:seq2+=seq1[i]
    seq=seq2

    N=[]#Shuffling of all four nucleotides, where possible. Affected aminoacids are four-codons and four-codon portion of six-codon aminoacids.
    seq1=[]
    for num in xrange(2,len(seq),3):
        if (seq[num]=='A'or seq[num]=='C'or seq[num]=='T'or seq[num]=='G')and(seq[num-2:num]=='TC'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CG'or seq[num-2:num]=='AC'or seq[num-2:num]=='GT'or seq[num-2:num]=='GC'or seq[num-2:num]=='GG'):
            N+=seq[num],
            seq1+=seq[num-2:num],'_N_',
        else:seq1+=seq[num-2:num+1],
    shuffle(N)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='_N_':seq2+=N.pop(0)
        else:seq2+=seq1[i]
    seq=seq2        
    return seq

def dn23(seq):#this function creates a randomized sequence, with dinucleotide frequences in codon position 2-3 close to those of the input sequence.
    aa=ag=ac=at=ga=gg=gc=gt=ca=cg=cc=ct=ta=tg=tc=tt=0
    for num in xrange(2,len(seq),3):#first calculating dinucleotide frequences in codon position 2-3 
        if seq[num-1]=='A':
            if seq[num]=='A':
                aa+=1
            elif seq[num]=='G':
                ag+=1
            elif seq[num]=='C':
                ac+=1
            elif seq[num]=='T':
                at+=1
        elif seq[num-1]=='G':
            if seq[num]=='A':
                ga+=1
            elif seq[num]=='G':
                gg+=1
            elif seq[num]=='C':
                gc+=1
            elif seq[num]=='T':
                gt+=1
        elif seq[num-1]=='C':
            if seq[num]=='A':
                ca+=1
            elif seq[num]=='G':
                cg+=1
            elif seq[num]=='C':
                cc+=1
            elif seq[num]=='T':
                ct+=1
        elif seq[num-1]=='T':
            if seq[num]=='A':
                ta+=1
            elif seq[num]=='G':
                tg+=1
            elif seq[num]=='C':
                tc+=1
            elif seq[num]=='T':
                tt+=1
    aa,ag,ac,at,ga,gg,gc,gt,ca,cg,cc,ct,ta,tg,tc,tt=aa/(len(seq)/3.),ag/(len(seq)/3.),ac/(len(seq)/3.),at/(len(seq)/3.),ga/(len(seq)/3.),gg/(len(seq)/3.),gc/(len(seq)/3.),gt/(len(seq)/3.),ca/(len(seq)/3.),cg/(len(seq)/3.),cc/(len(seq)/3.),ct/(len(seq)/3.),ta/(len(seq)/3.),tg/(len(seq)/3.),tc/(len(seq)/3.),tt/(len(seq)/3.)
    seq2=''
    for num in xrange(2,len(seq),3):#now each codon is replaced with a synonimous codon according to the dinucleotide frequences in codon position 2-3 
        seq2+=seq[num-2:num]
        if seq[num-1]=='A'and seq[num-2:num+1]!='TAA'and seq[num-2:num+1]!='TAG':
            if seq[num]=='T'or seq[num]=='C':
                space=at+ac
                AT,AC=at/space,ac/space
                x=random()
                if x<=AT:
                    seq2+='T'
                elif AT<x<=AT+AC:
                    seq2+='C'
            elif seq[num]=='A'or seq[num]=='G':
                space=aa+ag
                AA,AG=aa/space,ag/space
                x=random()
                if x<=AA:
                    seq2+='A'
                elif AA<x<=AA+AG:
                    seq2+='G'
            else:seq2+=seq[num]
        elif seq[num-1]=='G'and seq[num-2:num+1]!='TGA'and seq[num-2:num+1]!='TGG':
            if (seq[num-2]=='T'or seq[num-2]=='A')and(seq[num]=='C'or seq[num]=='T'):
                space = gt+gc
                GT,GC=gt/space,gc/space
                x=random()
                if x<=GT:
                    seq2+='T'
                elif GT<x<=GT+GC:
                    seq2+='C'
            elif seq[num-2:num+1]=='AGA'or seq[num-2:num+1]=='AGG':
                space=ga+gg
                GA,GG=ga/space,gg/space
                x=random()
                if x<=GA:
                    seq2+='A'
                elif GA<x<=GA+GG:
                    seq2+='G'
            elif seq[num-2]=='C'or seq[num-2]=='G':
                space=ga+gg+gc+gt
                GA,GG,GC,GT=ga/space,gg/space,gc/space,gt/space
                x=random()
                if x<=GA:seq2+='A'
                elif GA<x<=GA+GG:seq2+='G'
                elif GA+GG<x<=GA+GG+GC:seq2+='C'
                elif GA+GG+GC<x<=GA+GG+GC+GT:seq2+='T'
            else:seq2+=seq[num]
        elif seq[num-1]=='C':
            space=ca+cg+cc+ct
            CA,CG,CC,CT=ca/space,cg/space,cc/space,ct/space
            x=random()
            if x<=CA:seq2+='A'
            elif CA<x<=CA+CG:seq2+='G'
            elif CA+CG<x<=CA+CG+CC:seq2+='C'
            elif CA+CG+CC<x<=CA+CG+CC+CT:seq2+='T'
        elif seq[num-1]=='T':
            if seq[num-2:num+1]=='TTT'or seq[num-2:num+1]=='TTC':
                space = tt+tc
                TT,TC=tt/space,tc/space
                x=random()
                if x<=TT:
                    seq2+='T'
                elif TT<x<=TT+TC:
                    seq2+='C'
            elif seq[num-2:num+1]=='TTA'or seq[num-2:num+1]=='TTG':
                space = ta+tg
                TA,TG=ta/space,tg/space
                x=random()
                if x<=TA:
                    seq2+='A'
                elif TA<x<=TA+TG:
                    seq2+='G'
            elif seq[num-2:num+1]=='ATT'or seq[num-2:num+1]=='ATC'or seq[num-2:num+1]=='ATA':
                space=tt+tc+ta
                TT,TC,TA=tt/space,tc/space,ta/space
                x=random()
                if x<=TA:seq2+='A'
                elif TA<x<=TA+TC:seq2+='C'
                elif TA+TC<x<=TA+TC+TT:seq2+='T'
            elif seq[num-2]=='C'or seq[num-2]=='G':
                space=ta+tg+tc+tt
                TA,TG,TC,TT=ta/space,tg/space,tc/space,tt/space
                x=random()
                if x<=TA:seq2+='A'
                elif TA<x<=TA+TG:seq2+='G'
                elif TA+TG<x<=TA+TG+TC:seq2+='C'
                elif TA+TG+TC<x<=TA+TG+TC+TT:seq2+='T'
            else:seq2+=seq[num]
        else:seq2+=seq[num]
    seq=seq2
    return seq


def third(seq):#this function creates sequence with dinucleotide content in codon positions 2-3 and 3-1 identical to that of the input sequence.
    #logic of this function similar to "third_simple"
    seq1=[]
    TNT,TNC,TNA,TNG,GNT,GNA,GNC,GNG,CNG,CNA,CNT,CNC=[],[],[],[],[],[],[],[],[],[],[],[]#these lists will contain nucleotides from third codon position, with given [-1]and[+1]nucleotides to preserve dinucleotide content in conserved positions
    #four-codon and four-codon portion of six-codon aminoacids are affected
    for num in xrange(2,len(seq)-3,3):
        seq1+=seq[num-2:num]        
        if seq[num]=='T' or seq[num]=='C' or seq[num]=='A' or seq[num]=='G':
            if seq[num-2:num]=='CT' or seq[num-2:num]=='GT':#LEU4 or VAL
                if seq[num+1]=='T':
                    seq1+='TNT',
                    TNT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='TNC',
                    TNC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='TNA',
                    TNA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='TNG',
                    TNG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='C':#SER4 or PRO or THR or ALA
                if seq[num+1]=='T':
                    seq1+='CNT',
                    CNT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='CNC',
                    CNC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='CNA',
                    CNA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='CNG',
                    CNG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-2:num]=='CG' or seq[num-2:num]=='GG':#ARG4 or GLY
                if seq[num+1]=='T':
                    seq1+='GNT',
                    GNT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='GNC',
                    GNC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='GNA',
                    GNA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='GNG',
                    GNG+=seq[num],
                else:
                    seq1+=seq[num]
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-3:]
    
    shuffle(TNT),shuffle(TNC),shuffle(TNA),shuffle(TNG),shuffle(GNG),shuffle(GNA),shuffle(GNT),shuffle(GNC),shuffle(CNT),shuffle(CNC),shuffle(CNA),shuffle(CNG)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='TNT':seq2+=TNT.pop(0)
        elif seq1[i]=='TNC':seq2+=TNC.pop(0)
        elif seq1[i]=='TNG':seq2+=TNG.pop(0)
        elif seq1[i]=='TNA':seq2+=TNA.pop(0)
        elif seq1[i]=='GNT':seq2+=GNT.pop(0)
        elif seq1[i]=='GNA':seq2+=GNA.pop(0)
        elif seq1[i]=='GNC':seq2+=GNC.pop(0)
        elif seq1[i]=='GNG':seq2+=GNG.pop(0)
        elif seq1[i]=='CNT':seq2+=CNT.pop(0)
        elif seq1[i]=='CNC':seq2+=CNC.pop(0)
        elif seq1[i]=='CNG':seq2+=CNG.pop(0)
        elif seq1[i]=='CNA':seq2+=CNA.pop(0)
        else:seq2+=seq1[i]
    seq=seq2
    seq1=[]
    THT,THC,THA,THG,GHT,GHA,GHC,GHG,CHG,CHA,CHT,CHC=[],[],[],[],[],[],[],[],[],[],[],[]
    for num in xrange(2,len(seq)-3,3):
        seq1+=seq[num-2:num]
        if seq[num]=='T' or seq[num]=='C' or seq[num]=='A':
            if seq[num-2:num]=='CT' or seq[num-2:num]=='GT' or seq[num-2:num]=='AT':#ILE3 or LEU4 or VAL
                if seq[num+1]=='T':
                    seq1+='THT',
                    THT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='THC',
                    THC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='THA',
                    THA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='THG',
                    THG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='C':#SER4 or PRO or THR or ALA
                if seq[num+1]=='T':
                    seq1+='CHT',
                    CHT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='CHC',
                    CHC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='CHA',
                    CHA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='CHG',
                    CHG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-2:num]=='CG' or seq[num-2:num]=='GG':#ARG4 or GLY
                if seq[num+1]=='T':
                    seq1+='GHT',
                    GHT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='GHC',
                    GHC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='GHA',
                    GHA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='GHG',
                    GHG+=seq[num],
                else:
                    seq1+=seq[num]
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-3:]
    
    shuffle(THT),shuffle(THC),shuffle(THA),shuffle(THG),shuffle(GHG),shuffle(GHA),shuffle(GHT),shuffle(GHC),shuffle(CHT),shuffle(CHC),shuffle(CHA),shuffle(CHG)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='THT':seq2+=THT.pop(0)
        elif seq1[i]=='THC':seq2+=THC.pop(0)
        elif seq1[i]=='THA':seq2+=THA.pop(0)
        elif seq1[i]=='THG':seq2+=THG.pop(0)
        elif seq1[i]=='GHT':seq2+=GHT.pop(0)
        elif seq1[i]=='GHA':seq2+=GHA.pop(0)
        elif seq1[i]=='GHC':seq2+=GHC.pop(0)
        elif seq1[i]=='GHG':seq2+=GHG.pop(0)
        elif seq1[i]=='CHT':seq2+=CHT.pop(0)
        elif seq1[i]=='CHC':seq2+=CHC.pop(0)
        elif seq1[i]=='CHG':seq2+=CHG.pop(0)
        elif seq1[i]=='CHA':seq2+=CHA.pop(0)
        else:seq2+=seq1[i]
    seq=seq2
    seq1=[]
    TRT,TRC,TRA,TRG,ART,ARC,ARG,ARA,GRT,GRA,GRC,GRG,CRG,CRA,CRT,CRC=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    for num in xrange(2,len(seq)-3,3):
        seq1+=seq[num-2:num]
        if seq[num]=='A' or seq[num]=='G':
            if seq[num-1]=='T' and seq[num-2:num]!='AT':#not MET
                if seq[num+1]=='T':
                    seq1+='TRT',
                    TRT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='TRC',
                    TRC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='TRA',
                    TRA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='TRG',
                    TRG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='C':
                if seq[num+1]=='T':
                    seq1+='CRT',
                    CRT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='CRC',
                    CRC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='CRA',
                    CRA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='CRG',
                    CRG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='A' and seq[num-2:num]!='TA':#not Amber, not Ochre
                if seq[num+1]=='T':
                    seq1+='ART',
                    ART+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='ARC',
                    ARC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='ARA',
                    ARA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='ARG',
                    ARG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='G' and seq[num-2:num]!='TG':#not TRP, not Opal
                if seq[num+1]=='T':
                    seq1+='GRT',
                    GRT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='GRC',
                    GRC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='GRA',
                    GRA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='GRG',
                    GRG+=seq[num],
                else:
                    seq1+=seq[num]
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-3:]
    
    shuffle(TRT),shuffle(TRC),shuffle(TRA),shuffle(TRG),shuffle(GRG),shuffle(GRA),shuffle(GRT),shuffle(GRC),shuffle(ARG),shuffle(ARC),shuffle(ART),shuffle(ARA),shuffle(CRT),shuffle(CRC),shuffle(CRA),shuffle(CRG)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='TRT':seq2+=TRT.pop(0)
        elif seq1[i]=='TRC':seq2+=TRC.pop(0)
        elif seq1[i]=='TRA':seq2+=TRA.pop(0)
        elif seq1[i]=='TRG':seq2+=TRG.pop(0)
        elif seq1[i]=='ART':seq2+=ART.pop(0)
        elif seq1[i]=='ARC':seq2+=ARC.pop(0)
        elif seq1[i]=='ARG':seq2+=ARG.pop(0)
        elif seq1[i]=='ARA':seq2+=ARA.pop(0)
        elif seq1[i]=='GRT':seq2+=GRT.pop(0)
        elif seq1[i]=='GRA':seq2+=GRA.pop(0)
        elif seq1[i]=='GRC':seq2+=GRC.pop(0)
        elif seq1[i]=='GRG':seq2+=GRG.pop(0)
        elif seq1[i]=='CRT':seq2+=CRT.pop(0)
        elif seq1[i]=='CRC':seq2+=CRC.pop(0)
        elif seq1[i]=='CRG':seq2+=CRG.pop(0)
        elif seq1[i]=='CRA':seq2+=CRA.pop(0)
        else:seq2+=seq1[i]
    seq=seq2

    seq1=[]
    TYT,TYC,TYA,TYG,AYT,AYC,AYG,AYA,GYT,GYA,GYC,GYG,CYG,CYA,CYT,CYC=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    for num in xrange(2,len(seq)-3,3):
        seq1+=seq[num-2:num]
        if seq[num]=='T' or seq[num]=='C':
            if seq[num-1]=='T':
                if seq[num+1]=='T':
                    seq1+='TYT',
                    TYT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='TYC',
                    TYC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='TYA',
                    TYA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='TYG',
                    TYG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='C':
                if seq[num+1]=='T':
                    seq1+='CYT',
                    CYT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='CYC',
                    CYC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='CYA',
                    CYA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='CYG',
                    CYG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='A':
                if seq[num+1]=='T':
                    seq1+='AYT',
                    AYT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='AYC',
                    AYC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='AYA',
                    AYA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='AYG',
                    AYG+=seq[num],
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='G':
                if seq[num+1]=='T':
                    seq1+='GYT',
                    GYT+=seq[num],
                elif seq[num+1]=='C':
                    seq1+='GYC',
                    GYC+=seq[num],
                elif seq[num+1]=='A':
                    seq1+='GYA',
                    GYA+=seq[num],
                elif seq[num+1]=='G':
                    seq1+='GYG',
                    GYG+=seq[num],
                else:
                    seq1+=seq[num]
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-3:]
    
    shuffle(TYT),shuffle(TYC),shuffle(TYA),shuffle(TYG),shuffle(GYG),shuffle(GYA),shuffle(GYT),shuffle(GYC),shuffle(AYG),shuffle(AYC),shuffle(AYT),shuffle(AYA),shuffle(CYT),shuffle(CYC),shuffle(CYA),shuffle(CYG)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='TYT':seq2+=TYT.pop(0)
        elif seq1[i]=='TYC':seq2+=TYC.pop(0)
        elif seq1[i]=='TYA':seq2+=TYA.pop(0)
        elif seq1[i]=='TYG':seq2+=TYG.pop(0)
        elif seq1[i]=='AYT':seq2+=AYT.pop(0)
        elif seq1[i]=='AYC':seq2+=AYC.pop(0)
        elif seq1[i]=='AYG':seq2+=AYG.pop(0)
        elif seq1[i]=='AYA':seq2+=AYA.pop(0)
        elif seq1[i]=='GYT':seq2+=GYT.pop(0)
        elif seq1[i]=='GYA':seq2+=GYA.pop(0)
        elif seq1[i]=='GYC':seq2+=GYC.pop(0)
        elif seq1[i]=='GYG':seq2+=GYG.pop(0)
        elif seq1[i]=='CYT':seq2+=CYT.pop(0)
        elif seq1[i]=='CYC':seq2+=CYC.pop(0)
        elif seq1[i]=='CYG':seq2+=CYG.pop(0)
        elif seq1[i]=='CYA':seq2+=CYA.pop(0)
        else:seq2+=seq1[i]
    return seq2


def exchange6deg(seq):#this function shuffles the first nucleotide for two six-codon aminoacids (LEU and ARG) with the third codon position,
    #preserving overall dinucleotide content, but not position-specific dinucleotide content. For SER (TCN+AGY)shuffling is more complicated, see below.
    #after this shuffling, the 'third'-function should be used as above

    #LEU
    seq1=[]
    #LEU has TTR and CTN codons
    #first convertion of CTY to CTR is required to make them compatible with TTR codons and improve the shuffling efficiency;
    #CTY LEU codons might be re-introduced later upon the third-position shuffling 
    CTYT=[]#these lists will contain third nucleotides of CTY-codons 
    CTYA=[]
    CTTG=[]
    CTCG=[]#T/C separation for shuffling of ARG, (A->C)
    CTYC=[]
    
    TAT,TAA,TAG,TAC,TGT,TGA,TGG,TGC=[],[],[],[],[],[],[],[]#these lists will contain the third R-nucleotide of ILE and VAL to exchange with the first position of Leu
    TAG_arg=[]#only 'A' to 'C'. Arginine has AGY and CGN codons, first position can be shuffled. 
    
    seq1+=seq[:2]
    for num in xrange(2,len(seq)-2,1):
        #seq1+=seq[num-2:num]
        if num%3==2 and seq[num-2:num]=='CT' and (seq[num]=='C' or seq[num]=='T'):
            if seq[num+1]=='A':
                CTYA+=seq[num],
                seq1+='CTYA',
            elif seq[num+1]=='G':
                if seq[num]=='C':#T/C separation because of ARG, (A->C)
                    CTCG+=seq[num],
                    seq1+='CTYG',
                elif seq[num]=='T':
                    CTTG+=seq[num],
                    seq1+='CTYG',
            elif seq[num+1]=='C':
                CTYC+=seq[num],
                seq1+='CTYC',
            elif seq[num+1]=='T':
                CTYT+=seq[num],
                seq1+='CTYT',
            else:
                seq1+=seq[num],
        elif num%3==2 and (seq[num-2:num+1]=='GTA'or seq[num-2:num+1]=='ATA'):
            if seq[num+1]=='A':
                TAA+=seq[num],
                seq1+='TAA',
            elif seq[num+1]=='G':
                TAG+=seq[num],
                seq1+='TAG',
            elif seq[num+1]=='C':
                TAC+=seq[num],
                seq1+='TAC',
            elif seq[num+1]=='T':
                TAT+=seq[num],
                seq1+='TAT',
            else:
                seq1+=seq[num],
        elif num%3==2 and seq[num-2:num+1]=='GTG':
            if seq[num+1]=='A':
                TGA+=seq[num],
                seq1+='TGA',
            elif seq[num+1]=='G':
                TGG+=seq[num],
                seq1+='TGG',
            elif seq[num+1]=='C':
                TGC+=seq[num],
                seq1+='TGC',
            elif seq[num+1]=='T':
                TGT+=seq[num],
                seq1+='TGT',
            else:
                seq1+=seq[num],
        elif num%3==0 and (seq[num:num+3]=='AGA' or seq[num:num+3]=='AGG') and seq[num-1]=='T'and seq[num-2:num]!='TT'and seq[num-3:num]!='CTT':
            TAG_arg+=seq[num],
            seq1+='TAG_arg',
        else:
            seq1+=seq[num],
    seq1+=seq[-2:]
    
    #now replacing the third position Y with R in LEU 
    CTAG=[]
    change_num=min(len(CTCG),len(TAG_arg))
    CTAG,TAG_arg[:change_num],CTCG=TAG_arg[:change_num],CTCG[:change_num],CTCG[-1*(len(CTCG)-change_num):]
    CTYG=CTCG+CTTG
    shuffle(CTYG)
    change_num=min(len(CTYG),len(TAG))                                    #first TAN,
    CTYG[:change_num],TAG[:change_num]=TAG[:change_num],CTYG[:change_num] #then TGN,
    CTYG.reverse()                                                        #important,
    change_num=min(len(CTYG),len(TGG))                                    #Arginine,
    CTYG[:change_num],TGG[:change_num]=TGG[:change_num],CTYG[:change_num] #first TAN,
    CTYG+=CTAG                                                            #then TGN,
                                                                          #important,
    change_num=min(len(CTYT),len(TAT))                                    #first TAN,
    CTYT[:change_num],TAT[:change_num]=TAT[:change_num],CTYT[:change_num] #then TGN,
    CTYT.reverse()                                                        #important,
    change_num=min(len(CTYT),len(TGT))                                    #first TAN,
    CTYT[:change_num],TGT[:change_num]=TGT[:change_num],CTYT[:change_num] #then TGN,
                                                                          #important,
                                                                          #first TAN,
    change_num=min(len(CTYA),len(TAA))                                    #then TGN,
    CTYA[:change_num],TAA[:change_num]=TAA[:change_num],CTYA[:change_num] #important,
    CTYA.reverse()                                                        #Arginine,
    change_num=min(len(CTYA),len(TGA))                                    #then TGN,
    CTYA[:change_num],TGA[:change_num]=TGA[:change_num],CTYA[:change_num] #important,
                                                                          #first TAN,
                                                                          #then TGN,
    change_num=min(len(CTYC),len(TAC))                                    #important,
    CTYC[:change_num],TAC[:change_num]=TAC[:change_num],CTYC[:change_num] #first TAN,
    CTYC.reverse()                                                        #important,
    change_num=min(len(CTYC),len(TGC))                                    #Arginine,
    CTYC[:change_num],TGC[:change_num]=TGC[:change_num],CTYC[:change_num] #important
    
    shuffle(CTYT),shuffle(CTYA),shuffle(CTYG),shuffle(CTYC),shuffle(TAT),shuffle(TAA),shuffle(TAG),shuffle(TAC),shuffle(TGT),shuffle(TGA),shuffle(TGG),shuffle(TGC),shuffle(TAG_arg)
    
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='CTYT':seq2+=CTYT.pop(0)
        elif seq1[i]=='CTYC':seq2+=CTYC.pop(0)
        elif seq1[i]=='CTYA':seq2+=CTYA.pop(0)
        elif seq1[i]=='CTYG':seq2+=CTYG.pop(0)
        elif seq1[i]=='TAG':seq2+=TAG.pop(0)
        elif seq1[i]=='TAA':seq2+=TAA.pop(0)
        elif seq1[i]=='TAC':seq2+=TAC.pop(0)
        elif seq1[i]=='TAT':seq2+=TAT.pop(0)
        elif seq1[i]=='TGG':seq2+=TGG.pop(0)
        elif seq1[i]=='TGA':seq2+=TGA.pop(0)
        elif seq1[i]=='TGC':seq2+=TGC.pop(0)
        elif seq1[i]=='TGT':seq2+=TGT.pop(0)
        elif seq1[i]=='TAG_arg':seq2+=TAG_arg.pop(0)
        else:seq2+=seq1[i]
    seq=seq2#convertion of CTY to CTR is finished

    #now shuffling the first nucleotide of TTR and CTR LEU codons, and the third Y of other codons
    seq1=[]
    TYT,AYT,GYT,CYT=[],[],[],[]
    seq1+=seq[:2]
    
    for num in xrange(2,len(seq)-2,1):    
        if ((seq[num]=='T' or seq[num]=='C') and seq[num+1]=='T' and num%3==2 and seq[num+1:num+4]!='TTA' and seq[num+1:num+4]!='TTG' and seq[num+1:num+4]!='CTA' and seq[num+1:num+4]!='CTG') or (num%3==0 and (seq[num:num+3]=='TTA' or seq[num:num+3]=='TTG' or seq[num:num+3]=='CTA' or seq[num:num+3]=='CTG')):
            if seq[num-1]=='T':
                seq1+='TYT',
                TYT+=seq[num],
            elif seq[num-1]=='C':
                seq1+='CYT',
                CYT+=seq[num],
            elif seq[num-1]=='A':
                seq1+='AYT',
                AYT+=seq[num],
            elif seq[num-1]=='G':
                seq1+='GYT',
                GYT+=seq[num],
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-2:]
    
    shuffle(TYT),shuffle(GYT),shuffle(AYT),shuffle(CYT)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='TYT':seq2+=TYT.pop(0)
        elif seq1[i]=='AYT':seq2+=AYT.pop(0)
        elif seq1[i]=='GYT':seq2+=GYT.pop(0)
        elif seq1[i]=='CYT':seq2+=CYT.pop(0)
        else:seq2+=seq1[i]
    seq=seq2#shuffling is finished


    #SER
    seq1=[]
    
    TCRC,TCRA,TCRG,TCRT=[],[],[],[]#SER has TCN and AGY codons, first TCR will be converted to TCY
    CYC,CYA,CYG,CYT=[],[],[],[]
    
    for num in xrange(2,len(seq)-2,3):
        seq1+=seq[num-2:num]
        if seq[num-2:num+1]=='TCA' or seq[num-2:num+1]=='TCG':
            if seq[num+1]=='C':
                TCRC+=seq[num],
                seq1+='TCRC',
            elif seq[num+1]=='G':
                TCRG+=seq[num],
                seq1+='TCRG',
            elif seq[num+1]=='A':
                TCRA+=seq[num],
                seq1+='TCRA',
            elif seq[num+1]=='T':
                TCRT+=seq[num],
                seq1+='TCRT',
            else:
                seq1+=seq[num],
        elif (seq[num-2:num]=='CC'or seq[num-2:num]=='AC' or seq[num-2:num]=='GC')and(seq[num]=='T'or seq[num]=='C'):
            if seq[num+1]=='C':
                CYC+=seq[num],
                seq1+='CYC',
            elif seq[num+1]=='G':
                CYG+=seq[num],
                seq1+='CYG',
            elif seq[num+1]=='A':
                CYA+=seq[num],
                seq1+='CYA',
            elif seq[num+1]=='T':
                CYT+=seq[num],
                seq1+='CYT',
            else:
                seq1+=seq[num],
        else:
            seq1+=seq[num],
    seq1+=seq[-3:]
    
    change_num=min(len(TCRC),len(CYC))
    TCRC[:change_num],CYC[:change_num]=CYC[:change_num],TCRC[:change_num]
    
    change_num=min(len(TCRA),len(CYA))
    TCRA[:change_num],CYA[:change_num]=CYA[:change_num],TCRA[:change_num]
    
    change_num=min(len(TCRG),len(CYG))
    TCRG[:change_num],CYG[:change_num]=CYG[:change_num],TCRG[:change_num]
    
    change_num=min(len(TCRT),len(CYT))
    TCRT[:change_num],CYT[:change_num]=CYT[:change_num],TCRT[:change_num]
    
    shuffle(TCRC),shuffle(TCRA),shuffle(TCRG),shuffle(TCRT),shuffle(CYC),shuffle(CYA),shuffle(CYG),shuffle(CYT)
    
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='TCRC':seq2+=TCRC.pop(0)
        elif seq1[i]=='TCRA':seq2+=TCRA.pop(0)
        elif seq1[i]=='TCRG':seq2+=TCRG.pop(0)
        elif seq1[i]=='TCRT':seq2+=TCRT.pop(0)
        elif seq1[i]=='CYT':seq2+=CYT.pop(0)
        elif seq1[i]=='CYA':seq2+=CYA.pop(0)
        elif seq1[i]=='CYG':seq2+=CYG.pop(0)
        elif seq1[i]=='CYC':seq2+=CYC.pop(0)
        else:seq2+=seq1[i]
    
    seq=seq2#convertion is finished
    
    
    #ATCY,AAGY conversion to BTCY,BAGY
    seq1=[]
    AAA_R,GAA_N,GAA_R,CAA_N,TAA_R,TAA_N,TAA_H,AAT_R,GAT_N,GAT_R,CAT_N,TAT_R,TAT_N,TAT_H=[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    AGA_R,GBA_N,GGA_R,TGA_R,TBA_N,TYA_H,AGT_R,GBT_N,GGT_R,CBT_N,TGT_R,TBT_N,TYT_H,CBA_N=[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    
    for num in xrange(2,len(seq)-2,3):
        seq1+=seq[num-2:num]
        if seq[num:num+4]=='AAGC'or seq[num:num+4]=='AAGT':
            if seq[num-2:num]=='CA'or seq[num-2:num]=='AA'or seq[num-2:num]=='GA':
                AAA_R+=seq[num],
                seq1+='AAA_R',
            elif seq[num-2:num]=='CG'or seq[num-2:num]=='GG':
                GAA_N+=seq[num],
                seq1+='GAA_N',
            elif seq[num-2:num]=='AG':
                GAA_R+=seq[num],
                seq1+='GAA_R',
            elif seq[num-2:num]=='CC'or seq[num-2:num]=='AC'or seq[num-2:num]=='GC':
                CAA_N+=seq[num],
                seq1+='CAA_N',
            elif seq[num-2:num]=='TT':
                TAA_R+=seq[num],
                seq1+='TAA_R',
            elif seq[num-2:num]=='CT'or seq[num-2:num]=='GT':
                TAA_N+=seq[num],
                seq1+='TAA_N',
            elif seq[num-2:num]=='AT':
                TAA_H+=seq[num],
                seq1+='TAA_H',
            else:
                seq1+=seq[num]
        elif seq[num:num+4]=='ATCC'or seq[num:num+4]=='ATCT':
            if seq[num-2:num]=='CA'or seq[num-2:num]=='AA'or seq[num-2:num]=='GA':
                AAT_R+=seq[num],
                seq1+='AAT_R',
            elif seq[num-2:num]=='CG'or seq[num-2:num]=='GG':
                GAT_N+=seq[num],
                seq1+='GAT_N',
            elif seq[num-2:num]=='AG':
                GAT_R+=seq[num],
                seq1+='GAT_R',
            elif seq[num-2:num]=='CC'or seq[num-2:num]=='AC'or seq[num-2:num]=='GC':
                CAT_N+=seq[num],
                seq1+='CAT_N',
            elif seq[num-2:num]=='TT':
                TAT_R+=seq[num],
                seq1+='TAT_R',
            elif seq[num-2:num]=='CT'or seq[num-2:num]=='GT':
                TAT_N+=seq[num],
                seq1+='TAT_N',
            elif seq[num-2:num]=='AT':
                TAT_H+=seq[num],
                seq1+='TAT_H',
            else:
                seq1+=seq[num]
        elif seq[num+1:num+4]!='AGC'and seq[num+1:num+4]!='AGT'and seq[num+1]=='A'and seq[num]!='A':
            if seq[num-2:num+1]=='CAG'or seq[num-2:num+1]=='AAG'or seq[num-2:num+1]=='GAG':
                AGA_R+=seq[num],
                seq1+='AGA_R',
            elif seq[num-2:num]=='CG'or seq[num-2:num]=='GG':
                GBA_N+=seq[num],
                seq1+='GBA_N',
            elif seq[num-2:num+1]=='AGG':
                GGA_R+=seq[num],
                seq1+='GGA_R',
            elif seq[num-2:num]=='CC'or seq[num-2:num]=='AC'or seq[num-2:num]=='GC':
                CBA_N+=seq[num],
                seq1+='CBA_N',
            elif seq[num-2:num+1]=='TTG':
                TGA_R+=seq[num],
                seq1+='TGA_R',
            elif seq[num-2:num]=='CT'or seq[num-2:num]=='GT':
                TBA_N+=seq[num],
                seq1+='TBA_N',
            elif seq[num-2:num]=='AT'and seq[num]!='G':
                TYA_H+=seq[num],
                seq1+='TYA_H',
            else:
                seq1+=seq[num]
        elif seq[num+1:num+4]!='TCC'and seq[num+1:num+4]!='TCT'and seq[num+1]=='T'and seq[num]!='A':
            if seq[num-2:num+1]=='CAG'or seq[num-2:num+1]=='AAG'or seq[num-2:num+1]=='GAG':
                AGT_R+=seq[num],
                seq1+='AGT_R',
            elif seq[num-2:num]=='CG'or seq[num-2:num]=='GG':
                GBT_N+=seq[num],
                seq1+='GBT_N',
            elif seq[num-2:num+1]=='AGG':
                GGT_R+=seq[num],
                seq1+='GGT_R',
            elif seq[num-2:num]=='CC'or seq[num-2:num]=='AC'or seq[num-2:num]=='GC':
                CBT_N+=seq[num],
                seq1+='CBT_N',
            elif seq[num-2:num+1]=='TTG':
                TGT_R+=seq[num],
                seq1+='TGT_R',
            elif seq[num-2:num]=='CT'or seq[num-2:num]=='GT':
                TBT_N+=seq[num],
                seq1+='TBT_N',
            elif seq[num-2:num]=='AT'and seq[num]!='G':
                TYT_H+=seq[num],
                seq1+='TYT_H',
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-3:]
    
    
    change_num=min(len(AAA_R),len(AGA_R))
    AAA_R[:change_num],AGA_R[:change_num]=AGA_R[:change_num],AAA_R[:change_num]
    
    change_num=min(len(GAA_R),len(GGA_R))
    GAA_R[:change_num],GGA_R[:change_num]=GGA_R[:change_num],GAA_R[:change_num]
    
    change_num=min(len(TAA_R),len(TGA_R))
    TAA_R[:change_num],TGA_R[:change_num]=TGA_R[:change_num],TAA_R[:change_num]
    
    change_num=min(len(AAT_R),len(AGT_R))
    AAT_R[:change_num],AGT_R[:change_num]=AGT_R[:change_num],AAT_R[:change_num]
    
    change_num=min(len(GAT_R),len(GGT_R))
    GAT_R[:change_num],GGT_R[:change_num]=GGT_R[:change_num],GAT_R[:change_num]
    
    change_num=min(len(TAT_R),len(TGT_R))
    TAT_R[:change_num],TGT_R[:change_num]=TGT_R[:change_num],TAT_R[:change_num]
    
    change_num=min(len(TAT_H),len(TYT_H))
    TAT_H[:change_num],TYT_H[:change_num]=TYT_H[:change_num],TAT_H[:change_num]
    
    change_num=min(len(TAA_H),len(TYA_H))
    TAA_H[:change_num],TYA_H[:change_num]=TYA_H[:change_num],TAA_H[:change_num]
    
    change_num=min(len(CAA_N),len(CBA_N))
    CAA_N[:change_num],CBA_N[:change_num]=CBA_N[:change_num],CAA_N[:change_num]
    
    change_num=min(len(GAA_N),len(GBA_N))
    GAA_N[:change_num],GBA_N[:change_num]=GBA_N[:change_num],GAA_N[:change_num]
    
    change_num=min(len(TAA_N),len(TBA_N))
    TAA_N[:change_num],TBA_N[:change_num]=TBA_N[:change_num],TAA_N[:change_num]
    
    change_num=min(len(GAT_N),len(GBT_N))
    GAT_N[:change_num],GBT_N[:change_num]=GBT_N[:change_num],GAT_N[:change_num]
    
    change_num=min(len(CAT_N),len(CBT_N))
    CAT_N[:change_num],CBT_N[:change_num]=CBT_N[:change_num],CAT_N[:change_num]
    
    change_num=min(len(TAT_N),len(TBT_N))
    TAT_N[:change_num],TBT_N[:change_num]=TBT_N[:change_num],TAT_N[:change_num]
    
    shuffle(AAA_R),shuffle(GAA_N),shuffle(GAA_R),shuffle(CAA_N),shuffle(TAA_R),shuffle(TAA_N),shuffle(TAA_H),shuffle(AAT_R),shuffle(GAT_N),shuffle(GAT_R),shuffle(CAT_N),shuffle(TAT_R),shuffle(TAT_N),shuffle(TAT_H),shuffle(AGA_R),shuffle(GBA_N),shuffle(GGA_R),shuffle(TGA_R),shuffle(TBA_N),shuffle(TYA_H),shuffle(AGT_R),shuffle(GBT_N),shuffle(GGT_R),shuffle(CBT_N),shuffle(TGT_R),shuffle(TBT_N),shuffle(TYT_H),shuffle(CBA_N)
    
    
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='AAA_R':seq2+=AAA_R.pop(0)
        elif seq1[i]=='GAA_N':seq2+=GAA_N.pop(0)
        elif seq1[i]=='GAA_R':seq2+=GAA_R.pop(0)    
        elif seq1[i]=='CAA_N':seq2+=CAA_N.pop(0)
        elif seq1[i]=='TAA_R':seq2+=TAA_R.pop(0)
        elif seq1[i]=='TAA_N':seq2+=TAA_N.pop(0)
        elif seq1[i]=='TAA_H':seq2+=TAA_H.pop(0)
        elif seq1[i]=='AAT_R':seq2+=AAT_R.pop(0)
        elif seq1[i]=='GAT_N':seq2+=GAT_N.pop(0)
        elif seq1[i]=='GAT_R':seq2+=GAT_R.pop(0)
        elif seq1[i]=='CAT_N':seq2+=CAT_N.pop(0)
        elif seq1[i]=='TAT_R':seq2+=TAT_R.pop(0)
        elif seq1[i]=='TAT_N':seq2+=TAT_N.pop(0)
        elif seq1[i]=='TAT_H':seq2+=TAT_H.pop(0)
        elif seq1[i]=='AGA_R':seq2+=AGA_R.pop(0)
        elif seq1[i]=='GBA_N':seq2+=GBA_N.pop(0)
        elif seq1[i]=='GGA_R':seq2+=GGA_R.pop(0)
        elif seq1[i]=='TGA_R':seq2+=TGA_R.pop(0)
        elif seq1[i]=='TBA_N':seq2+=TBA_N.pop(0)
        elif seq1[i]=='TYA_H':seq2+=TYA_H.pop(0)
        elif seq1[i]=='AGT_R':seq2+=AGT_R.pop(0)
        elif seq1[i]=='GBT_N':seq2+=GBT_N.pop(0)
        elif seq1[i]=='GGT_R':seq2+=GGT_R.pop(0)
        elif seq1[i]=='CBT_N':seq2+=CBT_N.pop(0)
        elif seq1[i]=='TGT_R':seq2+=TGT_R.pop(0)
        elif seq1[i]=='TBT_N':seq2+=TBT_N.pop(0)
        elif seq1[i]=='TYT_H':seq2+=TYT_H.pop(0)
        elif seq1[i]=='CBA_N':seq2+=CBA_N.pop(0)
        else:seq2+=seq1[i]
    seq=seq2#ATCY,AAGY conversion to BTCY,BAGY is finished

    
    seq1=[]#tAGT SER codon could be converted to tTCT, preserving overall dinucleotide content, if in other positions in sequence
    #tTg will be replaced with tAG, and tCt will be replaced with tGt, and vice versa.
    #Using the same logic, SER codons ending with Y could be replaced, according to the third nucleotide of previous codon.
    
    GTG,GAG=[],[]#G_SER
    CTG,CAG=[],[]#C_SER
    TTG,TAG=[],[]#T_SER
    
    TCC,TGC=[],[]#SER_C
    TCT,TGT=[],[]#SER_T
    
    GAGC,CAGC,TAGC,GAGT,CAGT,TAGT=[],[],[],[],[],[]
    GTCC,CTCC,TTCC,GTCT,CTCT,TTCT=[],[],[],[],[],[]
    
    for num in xrange(2,len(seq)-3,3):
        if seq[num-2:num]=='AG'and(seq[num]=='C'or seq[num]=='T')and(seq[num-3]=='G'or seq[num-3]=='C'or seq[num-3]=='T'):
            if seq[num-3]=='G':
                if seq[num]=='C':
                    seq1+='GAGC','C'
                    GAGC+=seq[num-2:num],
                elif seq[num]=='T':
                    seq1+='GAGT','T'
                    GAGT+=seq[num-2:num],
            elif seq[num-3]=='C':
                if seq[num]=='C':
                    seq1+='CAGC','C'
                    CAGC+=seq[num-2:num],
                elif seq[num]=='T':
                    seq1+='CAGT','T'
                    CAGT+=seq[num-2:num],
            elif seq[num-3]=='T':
                if seq[num]=='C':
                    seq1+='TAGC','C'
                    TAGC+=seq[num-2:num],
                elif seq[num]=='T':
                    seq1+='TAGT','T'
                    TAGT+=seq[num-2:num],
        elif seq[num-2:num]=='TC'and(seq[num]=='C'or seq[num]=='T')and(seq[num-3]=='G'or seq[num-3]=='C'or seq[num-3]=='T'):
            if seq[num-3]=='G':
                if seq[num]=='C':
                    seq1+='GTCC','C'
                    GTCC+=seq[num-2:num],
                elif seq[num]=='T':
                    seq1+='GTCT','T'
                    GTCT+=seq[num-2:num],
            elif seq[num-3]=='C':
                if seq[num]=='C':
                    seq1+='CTCC','C'
                    CTCC+=seq[num-2:num],
                elif seq[num]=='T':
                    seq1+='CTCT','T'
                    CTCT+=seq[num-2:num],
            elif seq[num-3]=='T':
                if seq[num]=='C':
                    seq1+='TTCC','C'
                    TTCC+=seq[num-2:num],
                elif seq[num]=='T':
                    seq1+='TTCT','T'
                    TTCT+=seq[num-2:num],
        elif seq[num-1:num+2]=='GTG'and(seq[num-2:num]=='CG'or seq[num-2:num]=='GG'):
            seq1+=seq[num-2:num],'GTG',
            GTG+=seq[num],
        elif seq[num-1:num+2]=='GAG'and(seq[num-2:num]=='CG'or seq[num-2:num]=='GG'):
            seq1+=seq[num-2:num],'GAG',
            GAG+=seq[num],
        elif seq[num-1:num+2]=='CTG'and(seq[num-2:num]=='AC'or seq[num-2:num]=='GC'or seq[num-2:num]=='CC'):
            seq1+=seq[num-2:num],'CTG',
            CTG+=seq[num],
        elif seq[num-1:num+2]=='CAG'and(seq[num-2:num]=='AC'or seq[num-2:num]=='GC'or seq[num-2:num]=='CC'):
            seq1+=seq[num-2:num],'CAG',
            CAG+=seq[num],
        elif seq[num-1:num+2]=='TTG'and(seq[num-2:num]=='AT'or seq[num-2:num]=='GT'or seq[num-2:num]=='CT'):
            seq1+=seq[num-2:num],'TTG',
            TTG+=seq[num],
        elif seq[num-1:num+2]=='TAG'and(seq[num-2:num]=='AT'or seq[num-2:num]=='GT'or seq[num-2:num]=='CT'):
            seq1+=seq[num-2:num],'TAG',
            TAG+=seq[num],
        elif seq[num-1:num+2]=='TCC'and(seq[num-2:num]=='GT'or seq[num-2:num]=='CT'):
            seq1+=seq[num-2:num],'TCC',
            TCC+=seq[num],
        elif seq[num-1:num+2]=='TGC'and(seq[num-2:num]=='GT'or seq[num-2:num]=='CT'):
            seq1+=seq[num-2:num],'TGC',
            TGC+=seq[num],
        elif seq[num-1:num+2]=='TCT'and(seq[num-2:num]=='GT'or seq[num-2:num]=='CT')and seq[num+1:num+4]!='TCC'and seq[num+1:num+4]!='TCT':
            seq1+=seq[num-2:num],'TCT',
            TCT+=seq[num],
        elif seq[num-1:num+2]=='TGT'and(seq[num-2:num]=='GT'or seq[num-2:num]=='CT')and seq[num+1:num+4]!='TCC'and seq[num+1:num+4]!='TCT':
            seq1+=seq[num-2:num],'TGT',
            TGT+=seq[num],
        else:
            seq1+=seq[num-2:num+1]
    seq1+=seq[-3:]
    
    gAGc=min(len(GAGC),len(GTG),len(TCC))
    gTCc=min(len(GTCC),len(GAG),len(TGC))
    
    gAGt=min(len(GAGT),len(GTG),len(TCT))
    gTCt=min(len(GTCT),len(GAG),len(TGT))
    
    cAGc=min(len(CAGC),len(CTG),len(TCC))
    cTCc=min(len(CTCC),len(CAG),len(TGC))
    
    cAGt=min(len(CAGT),len(CTG),len(TCT))
    cTCt=min(len(CTCT),len(CAG),len(TGT))
    
    tAGc=min(len(TAGC),len(TTG),len(TCC))
    tTCc=min(len(TTCC),len(TAG),len(TGC))
    
    tAGt=min(len(TAGC),len(TTG),len(TCT))
    tTCt=min(len(TTCC),len(TAG),len(TGT))
    
    if gAGc+gAGt>len(GTG):
        gAGc,gAGt=round(gAGc*(len(GTG)/float(gAGc+gAGt))),round(gAGt*(len(GTG)/float(gAGc+gAGt)))
        gAGc,gAGt=int(gAGc),int(gAGt)
    if gTCc+gTCt>len(GAG):
        gTCc,gTCt=round(gTCc*(len(GAG)/float(gTCc+gTCt))),round(gTCt*(len(GAG)/float(gTCc+gTCt)))
        gTCc,gTCt=int(gTCc),int(gTCt)
    
    if cAGc+cAGt>len(CTG):
        cAGc,cAGt=round(cAGc*(len(CTG)/float(cAGc+cAGt))),round(cAGt*(len(CTG)/float(cAGc+cAGt)))
        cAGc,cAGt=int(cAGc),int(cAGt)
    if cTCc+cTCt>len(CAG):
        cTCc,cTCt=round(cTCc*(len(CAG)/float(cTCc+cTCt))),round(cTCt*(len(CAG)/float(cTCc+cTCt)))
        cTCc,cTCt=int(cTCc),int(cTCt)
    
    if tAGc+tAGt>len(TTG):
        tAGc,tAGt=round(tAGc*(len(TTG)/float(tAGc+tAGt))),round(tAGt*(len(TTG)/float(tAGc+tAGt)))
        tAGc,tAGt=int(tAGc),int(tAGt)
    if tTCc+tTCt>len(TAG):
        tTCc,tTCt=round(tTCc*(len(TAG)/float(tTCc+tTCt))),round(tTCt*(len(TAG)/float(tTCc+tTCt)))
        tTCc,tTCt=int(tTCc),int(tTCt)
    
    
    
    if gAGc+cAGc+tAGc>len(TCC):
        gAGc,cAGc,tAGc=round(gAGc*(len(TCC)/float(gAGc+cAGc+tAGc))),round(cAGc*(len(TCC)/float(gAGc+cAGc+tAGc))),round(tAGc*(len(TCC)/float(gAGc+cAGc+tAGc)))
        gAGc,cAGc,tAGc=int(gAGc),int(cAGc),int(tAGc)
    
    if gAGt+cAGt+tAGt>len(TCT):
        gAGt,cAGt,tAGt=round(gAGt*(len(TCT)/float(gAGt+cAGt+tAGt))),round(cAGt*(len(TCT)/float(gAGt+cAGt+tAGt))),round(tAGt*(len(TCT)/float(gAGt+cAGt+tAGt)))
        gAGt,cAGt,tAGt=int(gAGt),int(cAGt),int(tAGt)
    
    if gTCc+cTCc+tTCc>len(TGC):
        gTCc,cTCc,tTCc=round(gTCc*(len(TGC)/float(gTCc+cTCc+tTCc))),round(cTCc*(len(TGC)/float(gTCc+cTCc+tTCc))),round(tTCc*(len(TGC)/float(gTCc+cTCc+tTCc)))
        gTCc,cTCc,tTCc=int(gTCc),int(cTCc),int(tTCc)
    
    if gTCt+cTCt+tTCt>len(TGT):
        gTCt,cTCt,tTCt=round(gTCt*(len(TGT)/float(gTCt+cTCt+tTCt))),round(cTCt*(len(TGT)/float(gTCt+cTCt+tTCt))),round(tTCt*(len(TGT)/float(gTCt+cTCt+tTCt)))
        gTCt,cTCt,tTCt=int(gTCt),int(cTCt),int(tTCt)
    
    #GAGC,GTG,TCC
    change_num=randint(0,gAGc)
    for i in xrange(change_num):
        GAGC[i]='TC'
        del GTG[0];GTG.append('A')
        del TCC[0];TCC.append('G')
    
    #GTCC,GAG,TGC
    change_num=randint(0,gTCc)
    for i in xrange(change_num):
        GTCC[i]='AG'
        del GAG[0];GAG.append('T')
        del TGC[0];TGC.append('C')
    
    #GAGT,GTG,TCT
    change_num=randint(0,gAGt)
    for i in xrange(change_num):
        GAGT[i]='TC'
        del GTG[0];GTG.append('A')
        del TCT[0];TCT.append('G')
    
    #GTCT,GAG,TGT
    change_num=randint(0,gTCt)
    for i in xrange(change_num):
        GTCT[i]='AG'
        del GAG[0];GAG.append('T')
        del TGT[0];TGT.append('C')
    
    #CAGC,CTG,TCC
    change_num=randint(0,cAGc)
    for i in xrange(change_num):
        CAGC[i]='TC'
        del CTG[0];CTG.append('A')
        del TCC[0];TCC.append('G')
    
    #CTCC,CAG,TGC
    change_num=randint(0,cTCc)
    for i in xrange(change_num):
        CTCC[i]='AG'
        del CAG[0];CAG.append('T')
        del TGC[0];TGC.append('C')
    
    #CAGT,CTG,TCT
    change_num=randint(0,cAGt)
    for i in xrange(change_num):
        CAGT[i]='TC'
        del CTG[0];CTG.append('A')
        del TCT[0];TCT.append('G')
    
    #CTCT,CAG,TGT
    change_num=randint(0,cTCt)
    for i in xrange(change_num):
        CTCT[i]='AG'
        del CAG[0];CAG.append('T')
        del TGT[0];TGT.append('C')
    
    #TAGC,TTG,TCC
    change_num=randint(0,tAGc)
    for i in xrange(change_num):
        TAGC[i]='TC'
        del TTG[0];TTG.append('A')
        del TCC[0];TCC.append('G')
    
    #TTCC,TAG,TGC
    change_num=randint(0,tTCc)
    for i in xrange(change_num):
        TTCC[i]='AG'
        del TAG[0];TAG.append('T')
        del TGC[0];TGC.append('C')
    
    #TAGC,TTG,TCT
    change_num=randint(0,tAGt)
    for i in xrange(change_num):
        TAGC[i]='TC'
        del TTG[0];TTG.append('A')
        del TCT[0];TCT.append('G')
    
    #tTCt=min(len(TTCC),len(TAG),len(TGT))
    change_num=randint(0,tTCt)
    for i in xrange(change_num):
        TTCC[i]='AG'
        del TAG[0];TAG.append('T')
        del TGT[0];TGT.append('C')
    
    shuffle(GTG),shuffle(GAG),shuffle(CTG),shuffle(CAG),shuffle(TTG),shuffle(TAG),shuffle(TCC),shuffle(TGC),shuffle(TCT),shuffle(TGT)
    shuffle(GAGC),shuffle(CAGC),shuffle(TAGC),shuffle(GAGT),shuffle(CAGT),shuffle(TAGT)
    shuffle(GTCC),shuffle(CTCC),shuffle(TTCC),shuffle(GTCT),shuffle(CTCT),shuffle(TTCT)
    
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='GTG':seq2+=GTG.pop(0)
        elif seq1[i]=='GAG':seq2+=GAG.pop(0)
        elif seq1[i]=='CTG':seq2+=CTG.pop(0)
        elif seq1[i]=='CAG':seq2+=CAG.pop(0)
        elif seq1[i]=='TTG':seq2+=TTG.pop(0)
        elif seq1[i]=='TAG':seq2+=TAG.pop(0)
        elif seq1[i]=='TCC':seq2+=TCC.pop(0)
        elif seq1[i]=='TGC':seq2+=TGC.pop(0)
        elif seq1[i]=='TCT':seq2+=TCT.pop(0)
        elif seq1[i]=='TGT':seq2+=TGT.pop(0)
        elif seq1[i]=='GAGC':seq2+=GAGC.pop(0)
        elif seq1[i]=='CAGC':seq2+=CAGC.pop(0)    
        elif seq1[i]=='TAGC':seq2+=TAGC.pop(0)
        elif seq1[i]=='GAGT':seq2+=GAGT.pop(0)
        elif seq1[i]=='CAGT':seq2+=CAGT.pop(0)
        elif seq1[i]=='TAGT':seq2+=TAGT.pop(0)
        elif seq1[i]=='GTCC':seq2+=GTCC.pop(0)
        elif seq1[i]=='CTCC':seq2+=CTCC.pop(0)
        elif seq1[i]=='TTCC':seq2+=TTCC.pop(0)
        elif seq1[i]=='GTCT':seq2+=GTCT.pop(0)
        elif seq1[i]=='CTCT':seq2+=CTCT.pop(0)
        elif seq1[i]=='TTCT':seq2+=TTCT.pop(0)
        else:seq2+=seq1[i]
    seq=seq2


    #ARG
    #shuffling of the first codon position of aginine codons requires convertion of CGY to CGR.
    #then the first codon positions of AGR and CGR are shuffled with the third codon position of other codons
    #aminoacid sequence and overall dinucleotide content are maintained

    seq1=[]
    
    GYA,GYT,GYC,GYG=[],[],[],[]
    
    GRA=[]#only gly
    GRT,GRC,GRG=[],[],[]
    
    for num in xrange(2,len(seq)-3,3):
        seq1+=seq[num-2:num]
        if seq[num-2:num]=='CG' and(seq[num]=='C'or seq[num]=='T'):
            if seq[num+1]=='A':
                GYA+=seq[num],
                seq1+='GYA',
            elif seq[num+1]=='T':
                GYT+=seq[num],
                seq1+='GYT',
            elif seq[num+1]=='C':
                GYC+=seq[num],
                seq1+='GYC',
            elif seq[num+1]=='G':
                GYG+=seq[num],
                seq1+='GYG',
            else:
                seq1+=seq[num],
        elif seq[num-2:num]=='GG' and(seq[num]=='A'or seq[num]=='G'):#only gly
            if seq[num+1]=='A':
                GRA+=seq[num],
                seq1+='GRA',        
            elif seq[num+1]=='T':
                GRT+=seq[num],
                seq1+='GRT',
            elif seq[num+1]=='C':
                GRC+=seq[num],
                seq1+='GRC',
            elif seq[num+1]=='G':
                GRG+=seq[num],
                seq1+='GRG',
            else:
                seq1+=seq[num],
        else:
            seq1+=seq[num],
    seq1+=seq[-3:]
    
    change_num=min(len(GYA),len(GRA))
    GYA[:change_num],GRA[:change_num]=GRA[:change_num],GYA[:change_num]
    
    change_num=min(len(GYT),len(GRT))
    GYT[:change_num],GRT[:change_num]=GRT[:change_num],GYT[:change_num]
    
    change_num=min(len(GYC),len(GRC))
    GYC[:change_num],GRC[:change_num]=GRC[:change_num],GYC[:change_num]
    
    change_num=min(len(GYG),len(GRG))
    GYG[:change_num],GRG[:change_num]=GRG[:change_num],GYG[:change_num]
    
    shuffle(GYA),shuffle(GYT),shuffle(GYC),shuffle(GYG),shuffle(GRA),shuffle(GRG),shuffle(GRC),shuffle(GRT)
    
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='GYA':seq2+=GYA.pop(0)
        elif seq1[i]=='GRA':seq2+=GRA.pop(0)
        elif seq1[i]=='GYT':seq2+=GYT.pop(0)
        elif seq1[i]=='GRT':seq2+=GRT.pop(0)
        elif seq1[i]=='GYC':seq2+=GYC.pop(0)
        elif seq1[i]=='GRC':seq2+=GRC.pop(0)
        elif seq1[i]=='GYG':seq2+=GYG.pop(0)
        elif seq1[i]=='GRG':seq2+=GRG.pop(0)
        else:seq2+=seq1[i]
    
    seq=seq2
    seq1=[]
    
    CGTG,CGCG,GTG,GCG,GAG=[],[],[],[],[]
    
    seq1+=seq[:2]
    for num in xrange(2,len(seq)-2,1):
        if seq[num-2:num+2]=='CGCG' and num%3==2:
            CGCG+=seq[num],
            seq1+='CGYG',
        elif seq[num-2:num+2]=='CGTG' and num%3==2:
            CGTG+=seq[num],
            seq1+='CGYG',
        elif seq[num-1:num+2]=='GTG' and seq[num-2]!='C' and num%3==2:
            GTG+=seq[num],
            seq1+='GTG',
        elif seq[num-1:num+2]=='GCG' and seq[num-2]!='C' and num%3==2:
            GCG+=seq[num],
            seq1+='GCG',
        elif (seq[num-1:num+3]=='GAGA' or seq[num-1:num+3]=='GAGG')and num%3==0:
            GAG+=seq[num],
            seq1+='GAG',
        else:
            seq1+=seq[num],
    seq1+=seq[-2:]
    
    CGYG=[]
    change_num=min(len(CGTG),len(GCG))
    CGYG[:change_num],GCG[:change_num]=GCG[:change_num],CGTG[:change_num]
    CGTG_num=change_num
    
    CGYG+=CGCG
    change_num=min(len(CGYG),len(GAG))
    CGYG[:change_num],GAG[:change_num]=GAG[:change_num],CGYG[:change_num]
    CGYG+=CGTG[CGTG_num:]
    
    shuffle(CGYG),shuffle(GTG),shuffle(GCG),shuffle(GAG)
    
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='CGYG':seq2+=CGYG.pop(0)
        elif seq1[i]=='GTG':seq2+=GTG.pop(0)
        elif seq1[i]=='GCG':seq2+=GCG.pop(0)
        elif seq1[i]=='GAG':seq2+=GAG.pop(0)
        else:seq2+=seq1[i]
    seq=seq2
    
    seq1=[]
    TMG,AMG,GMG,CMG=[],[],[],[]
    seq1+=seq[:2],
    for num in xrange(2,len(seq)-2,1):
        if (num%3==0 and(seq[num:num+3]=='CGA' or seq[num:num+3]=='CGG' or seq[num:num+3]=='AGA' or seq[num:num+3]=='AGG'))or(num%3==2 and (seq[num:num+2]=='CG' or seq[num:num+2]=='AG')and((seq[num-1]=='T'and seq[num-2]!='T')or seq[num-1]=='C'or seq[num-2:num]=='GG')):
            if seq[num-1]=='T':
                seq1+='TMG',
                TMG+=seq[num],
            elif seq[num-1]=='C':
                seq1+='CMG',
                CMG+=seq[num],
            elif seq[num-1]=='A':
                seq1+='AMG',
                AMG+=seq[num],
            elif seq[num-1]=='G':
                seq1+='GMG',
                GMG+=seq[num],
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-2:]
    shuffle(TMG),shuffle(GMG),shuffle(AMG),shuffle(CMG)
    seq2=''
    for i in xrange(len(seq1)):
        if seq1[i]=='TMG':seq2+=TMG.pop(0)
        elif seq1[i]=='AMG':seq2+=AMG.pop(0)
        elif seq1[i]=='GMG':seq2+=GMG.pop(0)
        elif seq1[i]=='CMG':seq2+=CMG.pop(0)
        else:seq2+=seq1[i]
    return seq2

def get_difference(seq1,seq2):
    assert len(seq1) == len(seq2)
    return sum(seq1 != seq2 for seq1, seq2 in zip(seq1,seq2))

def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq = nuc_record.seq.translate(to_stop=True), \
                     id = "trans_" + nuc_record.id, \
                     description = "translation of CDS, using default table")


# main script body that calls the above subroutines

#Shuffle Script
# -*- coding: cp1251 -*-

#Input Data
parser = argparse.ArgumentParser(description='CodonShuffle.')
parser.add_argument('-i', nargs='?', help='Input Filename', required=True, dest="input_file_name")
parser.add_argument('-s', choices=['dn23', 'dn31', 'dn231', 'n3'], nargs='?', help='Type of shuffle', default="dn23", dest="random_type")
parser.add_argument('-r', nargs='?', help='Number of replications (int)', default='1000', dest="reps", type=int)
parser.add_argument('-m', choices=['CAI', 'CPB', 'DN', 'ENC', 'VFOLD', 'UFOLD', 'all'], nargs='*', help='Control Features [select one, multiple, or all]', default='all', dest="modules")
parser.add_argument('-g', dest="graphics", help='Generate Feature Plots', action="store_true")
args = parser.parse_args()

types_of_rnd=args.random_type

infile=open(args.input_file_name,'r')
names_list=[]
data={}
for line in infile:
    if line[0]=='>':
        strain=line
        data[strain]=''
        names_list+=strain,
    else:
        for liter in line:
            if liter!='\n' and liter!='-' and liter!='~':
                data[strain]+=liter
infile.close()


out_names=[]
for strain in names_list:
    seq_name=''
    for liter in strain:
        if liter =='\\' or liter =='/' or liter ==' ' or liter =='-' or liter ==',' or liter=='|': #replacing non-DOS chars in sequence names
            seq_name+='_'
        else:
            seq_name+=liter
    inseq_file=open(seq_name[1:-1]+'.fas','w')
    inseq_file.write(strain+data[strain]+'\n')
    inseq_file.close()
#     bat_enc='chips  -seqall '+seq_name[1:-1]+'.fas -nosum -outfile '+seq_name[1:-1]+'.enc -auto\n'
#     system(bat_enc)

    outfile=open(seq_name[1:-1]+'_'+args.random_type+'.fas','w') #Create the file with wild type in the first position
    outfile.write(strain+data[strain]+'\n')
    outfile.close()


    outfile=open(seq_name[1:-1]+'_'+args.random_type+'.fas','a') #Append permuted sequence
    for i in xrange(args.reps):
        outseq=data[strain]
        if args.random_type =='gc3':
            outseq=gc3(outseq)
        elif args.random_type == 'dn23':
            outseq=dn23(outseq)
        elif args.random_type=='dn31':
            outseq=third(outseq)
        elif args.random_type=='dn231':
            outseq=third(exchange6deg(outseq))
        elif args.random_type=='3n':
            outseq=third_simple(outseq)
        outfile.write('>replicate_'+str(i+1)+'\n'+outseq+'\n')
    outfile.close()
#     bat_enc='chips  -seqall '+seq_name[1:-1]+'_'+args.random_type+'.fas -nosum -outfile '+seq_name[1:-1]+'_'+args.random_type+'.enc -auto\n'
#     system(bat_enc)

first_tot_out_string=`args.reps`+'_replicas\nstrain\t\tstrain_Nc\t'
#if 'gc3'in types_of_rnd:first_tot_out_string+='\tmean_Nc_gc3\tsd\t'
if '3n'in types_of_rnd:first_tot_out_string+='\tmean_Nc_3n\tsd\t'
if 'all'in types_of_rnd:first_tot_out_string+='\tmean_Nc_all\tsd\t'
if 'dn23'in types_of_rnd:first_tot_out_string+='\tmean_Nc_dn23\tsd\t'
if 'dn31'in types_of_rnd:first_tot_out_string+='\tmean_Nc_dn31\tsd\t'
if 'dn231'in types_of_rnd:first_tot_out_string+='\tmean_Nc_dn231\tsd\t'

# tot_out=open(args.input_file_name[:-4]+'.out','w') #output file writing
# tot_out.write(first_tot_out_string+'\n')
# for strain in names_list:
#     seq_name=''
#     for liter in strain:
#         if liter =='\\' or liter =='/' or liter ==' ' or liter =='-' or liter ==',' or liter=='|':
#             seq_name+='_'
#         else:
#             seq_name+=liter
#     enc_in=open(seq_name[1:-1]+'.enc','r')
#     inseq_enc=0.
#     for i in enc_in:
#         if strain[1:-1] in i:
#             inseq_enc=float(i.split(' Nc = ')[1])
#     enc_in.close()
# 
#     out_string=strain[1:-1]+'\t\t'+`inseq_enc`+'\t\t'
#     enc_in=open(seq_name[1:-1]+'_'+args.random_type+'.enc','r')
#     enc=[]
#     for i in enc_in:
#         enc+=float(i.split(' Nc = ')[1]),
#     enc_in.close()
# 
#     mean_enc=sum(enc)/len(enc)
#     dd=0
#     for i in enc:d=(i-mean_enc)**2;dd+=d
#     sd_enc=(dd/(len(enc)-1))**(0.5)
#     out_string+=`mean_enc`+'\t'+`sd_enc`+'\t\t'
#     out_string+='\n'
#     tot_out.write(out_string)
# tot_out.close()

#filename.input(first_tot_out_string+'\n', inplace=1) #Add the wild type sequence in the variable

filename = seq_name[1:-1]+'_'+args.random_type+'.fas'

final_table=pandas.DataFrame()

#Calculate Sequence distance
seq_records=[]

inseq_file=open(filename,'rU')
outnt_file=open(filename+'.hamming', 'w')

for seq_record in SeqIO.parse(inseq_file, "fasta"):
    seq_records.append(seq_record.seq)

n = len(seq_records)
least_squares = pandas.DataFrame(numpy.zeros(n))
 
my_array = np.zeros((n,n))
for i in range(0,n):
    for j in range(i+1,n):
        difference = get_difference(seq_records[i], seq_records[j])
        outnt_file.write(str(difference)+'\n')
        my_array[i, j] = difference
        my_array[j, i] = difference
#nuc_dist = my_array.iloc[:,0]
outnt_file.close()


if args.graphics:
    hamming_graphname = filename+'.hamming.pdf'
    hamming_table = pandas.read_csv(filename+'.hamming', sep='\t', names=['distance'])
    hamming_graph = ggplot(hamming_table, aes('distance'))+geom_density() +labs("Hamming distance","Frequency")+ geom_vline(xintercept = [hamming_table['distance'].iloc[0]] , colour="red", linetype = "dashed") +ggtitle(seq_name[1:-1]+'_'+args.random_type+' (Hamming)') #Get the name of the script
    ggsave(hamming_graph, hamming_graphname)

if 'CAI' in args.modules or 'all' in args.modules:
    #Run CAI and read result table
    cainame = seq_name[1:-1]+'_'+args.random_type+'.cai'
    call(["./lib/EMBOSS-6.6.0/emboss/cai", "-seqall="+filename, "-outfile="+cainame, "-cfile=Eyeast_cai.cut"]) #Insert path before cai in this line (CAI)
    u_cols = ['a', 'sequence', 'b', 'cai']
    cai_table = pandas.read_csv(cainame, sep=' ', names=u_cols)
    cai_table = cai_table.drop('a', 1)
    cai_table = cai_table.drop('b', 1)
    cai_table_z =  ((cai_table['cai'] - cai_table['cai'].mean()) / cai_table['cai'].std())
    cai_table_z_ls = (cai_table_z-cai_table_z[0])**2

    least_squares = least_squares.add(cai_table_z_ls, axis=0)
#    final_table=final_table.append(cai_table['cai'])
    final_table.insert(0, "Cai", cai_table['cai'])
    
    if args.graphics:
        cai_graphname = cainame +'.pdf'
        u_cols = ['a', 'sequence', 'b', 'cai']
        cai_table = pandas.read_csv(cainame, sep=' ', names=u_cols)
        cai_table = cai_table.drop('a', 1)
        cai_table = cai_table.drop('b', 1)

        cai_graph = ggplot(cai_table, aes('cai'))+geom_density() +labs("CAI","Frequency")+ geom_vline(xintercept = [cai_table['cai'].iloc[0]] , colour="red", linetype = "dashed") +ggtitle(seq_name[1:-1]+'_'+args.random_type+' (CAI)') #Get the name of the script
        ggsave(cai_graph, cai_graphname)

if 'ENC' in args.modules or 'all' in args.modules:
    #  Run ENC and result table
    call(["./lib/codonW/codonw", filename, "-enc", "-nomenu", "-nowarn", "-silent"]) #Insert path before codonw in this line (ENC)
    enc_filename = seq_name[1:-1]+'_'+args.random_type+'.out'
    enc_table = pandas.read_csv(enc_filename, sep='\t')
    enc_table = enc_table.drop('Unnamed: 2',1)
    enc_table_z =  ((enc_table['Nc'] - enc_table['Nc'].mean()) / enc_table['Nc'].std())
    enc_table_z_ls = (enc_table_z-enc_table_z[0])**2
    
    least_squares = least_squares.add(enc_table_z_ls, axis=0)
#    final_table=final_table.append(enc_table['Nc'])
    final_table.insert(0, "ENC", enc_table['Nc'])
    
    if args.graphics:
        enc_filename = seq_name[1:-1]+'_'+args.random_type+'.out'
        enc_graphname = enc_filename+'.enc.pdf'
        enc_table = pandas.read_csv(enc_filename, sep='\t')
        enc_table = enc_table.drop('Unnamed: 2',1)
        enc_graph = ggplot(enc_table, aes('Nc'))+geom_density() +labs("ENC","Frequency")+ geom_vline(xintercept = [enc_table['Nc'].iloc[0]] , colour="red", linetype = "dashed") +ggtitle(seq_name[1:-1]+'_'+args.random_type+' (ENC)') #Get the name of the script
        ggsave(enc_graph, enc_graphname)


if 'VFOLD' in args.modules or 'all' in args.modules:
#Read FOLD (minimum free energy) table
    foldname = seq_name[1:-1]+'_'+args.random_type+'.fold'
    mfename = seq_name[1:-1]+'_'+args.random_type+'.mfe'
    i = open(filename, "r")
    o = open(foldname, "w")
    call(["./lib/ViennaRNA-2.1.9/Progs/RNAfold", "--noPS"], stdin=i, stdout=o) #Insert path before RNAfold in this line (MFOLD)
    i.close
    o.close
#    os.system("cat Poliovirus_1_Mahoney_P1_dn23.fold | sed 'N;N;s/\\n/ /g' | cut -f 4 -d ' ' | tr -d '()' > " + mfename)

    fold_tb = open(foldname, "r").read().split()
    fold_file=open(filename+'fold_table_mfe.txt', 'w')
    for i in range(3, len(fold_tb)-3, 4):
        fold_mfe = fold_tb[3]+'\n'+fold_tb[i+4]
        fold_file.write(str(fold_mfe)+'\n')
    fold_file.close()
    
    fold_table = pandas.read_csv(filename+'fold_table_mfe.txt', names=['mfe'])
    fold_table['mfe'] = fold_table['mfe'].map(lambda x: x.lstrip('(').rstrip(')'))
    fold_table['mfe']=fold_table.apply(lambda row: float(row['mfe']), axis=1)
    fold_table_z = ((fold_table['mfe'] - fold_table['mfe'].mean()) / fold_table['mfe'].std())
    fold_table_z_ls = (fold_table_z-fold_table_z[0])**2

    least_squares = least_squares.add(fold_table_z_ls, axis=0)
#    final_table=final_table.append(fold_table['mfe'])
    final_table.insert(0, "VFOLD (mfe)", fold_table['mfe'])
        
    if args.graphics:
        fold_graphname = seq_name[1:-1]+'_'+args.random_type+'.fold.pdf'
#        fold_table = pandas.read_csv(mfename, names=['mfe'])
        fold_graph = ggplot(fold_table, aes('mfe'))+geom_density() +labs("MFE","Frequency")+ geom_vline(xintercept = [fold_table['mfe'].iloc[0]] , colour="red", linetype = "dashed") +ggtitle(seq_name[1:-1]+'_'+args.random_type+' (FOLD)') #Get the name of the script
        ggsave(fold_graph, fold_graphname)


if 'UFOLD' in args.modules:
#Read FOLD (minimum free energy) table
    foldname = seq_name[1:-1]+'_'+args.random_type+'.fold'
    mfename = seq_name[1:-1]+'_'+args.random_type+'.fasta.dG'
    i = open(filename, "r")
    o = open(foldname, "w")
    call(["hybrid-ss", "-E", "--output="+mfename, filename]) #Insert path before hybrid-ss in this line (UNAFOLD)
    i.close
    o.close
#    os.system("cat Poliovirus_1_Mahoney_P1_dn23.fold | sed 'N;N;s/\\n/ /g' | cut -f 4 -d ' ' | tr -d '()' > " + mfename)
    
    ufold_table = pandas.read_csv(mfename, sep='	')
    ufold_table_z =  ((fold_table['-RT ln Z'] - fold_table['-RT ln Z'].mean()) / fold_table['-RT ln Z'].std())
    ufold_table_z_ls = (fold_table_z-fold_table_z[0])**2

    least_squares = least_squares.add(fold_table_z_ls, axis=0)
#    final_table=final_table.append(ufold_table['-RT ln Z'])
    final_table.insert(0, "UFOLD (mfe)", ufold_table['-RT ln Z'])
        
    if args.graphics:
        fold_graphname = seq_name[1:-1]+'_'+args.random_type+'.fold.pdf'
        fold_table = pandas.read_csv(mfename, mfename, sep='	')
        fold_graph = ggplot(fold_table, aes('-RT ln Z'))+geom_density() +labs("MFE","Frequency")+ geom_vline(xintercept = [fold_table['-RT ln Z'].iloc[0]] , colour="red", linetype = "dashed") +ggtitle(seq_name[1:-1]+'_'+args.random_type+' (FOLD)') #Get the name of the script
        ggsave(fold_graph, fold_graphname)


if 'DN' in args.modules or 'all' in args.modules:
    dnname = seq_name[1:-1]+'_'+args.random_type+'.dn'
    dn_file=open(dnname, 'w')
       
    dn_file.write("id");
    for nt1 in nts:
        for nt2 in nts:
            dinut = nt1+nt2
            dn_file.write("\t"+dinut)                
    dn_file.write("\n") 

    for nuc_rec in SeqIO.parse(filename, "fasta"):
        nucs = [str(nuc_rec.seq[i:i+1]) for i in range(0,len(nuc_rec.seq),1)]
        dinucs = [str(nuc_rec.seq[i:i+2]) for i in range(0,len(nuc_rec.seq)-1,1)]  
        
        nuc_counts = Counter(nucs)
        dinuc_counts = Counter(dinucs)
        seq_len = len(nuc_rec.seq)

        dn_file.write(nuc_rec.id)
        for nt1 in nts:
            for nt2 in nts:
                dinut = nt1+nt2
                if (dinut in dinuc_counts):
                    freq = (dinuc_counts[dinut] / (seq_len - 1)) / ( (nuc_counts[nt1] / seq_len) * (nuc_counts[nt2] / seq_len))
                    #print(nt1 + " " + nt2 + " " + str(freq) + " " + str(dinuc_counts[dinut]) + " " + str(nuc_counts[nt1]) + " " + str(nuc_counts[nt2]) + "\n")
                    dn_file.write("\t"+str(freq))
                else:
                    dn_file.write("\t0")
        dn_file.write("\n")
    dn_file.close()

    dnlsname = seq_name[1:- 1]+'_'+args.random_type+'.dnls'
#    dn_least_file=open(dnlsname, 'w')
    dn_table = pandas.read_csv(dnname, sep='	')
    dn_table_least = np.sqrt((dn_table['AA']-dn_table.iloc[0,1])**2+(dn_table['AC']-dn_table.iloc[0,2])**2+(dn_table['AG']-dn_table.iloc[0,3])**2+(dn_table['AT']-dn_table.iloc[0,4])**2+(dn_table['CA']-dn_table.iloc[0,5])**2+(dn_table['CC']-dn_table.iloc[0,6])**2+(dn_table['CG']-dn_table.iloc[0,7])**2+(dn_table['CT']-dn_table.iloc[0,8])**2+(dn_table['GA']-dn_table.iloc[0,9])**2+(dn_table['GC']-dn_table.iloc[0,10])**2+(dn_table['GG']-dn_table.iloc[0,11])**2+(dn_table['GT']-dn_table.iloc[0,12])**2+(dn_table['TA']-dn_table.iloc[0,13])**2+(dn_table['TC']-dn_table.iloc[0,14])**2+(dn_table['TG']-dn_table.iloc[0,15])**2+(dn_table['TT']-dn_table.iloc[0,16])**2)
    dn_table_least.to_csv(dnlsname, sep="\t")
#    dn_least_file.write(dn_table_least)
#    dn_least_file.close()
    dn_table_least_z =  ((dn_table_least - dn_table_least.mean()) / dn_table_least.std())
    dn_table_least_z_ls = (dn_table_least_z-dn_table_least_z[0])**2

    
    least_squares = least_squares.add(dn_table_least_z_ls, axis=0)
    u_cols = ['Replication', 'DN_least_square']
    dn_table_ls = pandas.read_csv(dnlsname, sep='	', names=u_cols)
#    final_table=final_table.append(dn_table_ls['DN_least_square'])
    final_table.insert(0, "DN_least_square", dn_table_ls['DN_least_square'])
    final_table.insert(1, "DN (TT)", dn_table['TT'])
    final_table.insert(1, "DN (TG)", dn_table['TG'])
    final_table.insert(1, "DN (TC)", dn_table['TC'])
    final_table.insert(1, "DN (TA)", dn_table['TA'])
    final_table.insert(1, "DN (GT)", dn_table['GT'])
    final_table.insert(1, "DN (GG)", dn_table['GG'])
    final_table.insert(1, "DN (GC)", dn_table['GC'])
    final_table.insert(1, "DN (GA)", dn_table['GA'])
    final_table.insert(1, "DN (CT)", dn_table['CT'])
    final_table.insert(1, "DN (CG)", dn_table['CG'])
    final_table.insert(1, "DN (CC)", dn_table['CC'])
    final_table.insert(1, "DN (CA)", dn_table['CA'])
    final_table.insert(1, "DN (AT)", dn_table['AT'])
    final_table.insert(1, "DN (AG)", dn_table['AG'])
    final_table.insert(1, "DN (AC)", dn_table['AC'])
    final_table.insert(1, "DN (AA)", dn_table['AA'])
    

    if args.graphics:
        dnls_graphname = dnlsname + '.pdf'

        #--bug in python ggplot for this, use rpy2 instead--
        dn_table_least = pandas.read_csv(dnlsname, sep='	', names=['Rep', 'DN']) 
        dn_graph = ggplot(dn_table_least,aes('DN')) + geom_density() + xlab("dinucleotide") + ylab('Dinucleotide frequency')+ geom_vline(xintercept = [dn_table_least['DN'].iloc[0]] , colour="red", linetype = "dashed") +ggtitle(seq_name[1:-1]+'_'+args.random_type+' (DN)')
        ggsave(dn_graph, dnls_graphname)
        
        dn_table = pandas.read_csv(dnname, sep="	") 
        fig, ax = plt.subplots()
        dn_table.boxplot(return_type='axes')
        ax.set_title(seq_name[1:- 1]+'_'+args.random_type+' (Dinucleotide Frequency)')
        ax.set_xlabel('Dinucleotide')
        ax.set_ylabel('Dinuc obser/Dinuc expec')
        fig.savefig(seq_name[1:-1]+'_'+args.random_type+'_dn.pdf')
        plt.close(fig)

       
#         r = robjects.r
#         r.library("ggplot2")
#         r.library("reshape2")
#         robjects.r('dn_table=read.csv("Poliovirus_1_Mahoney_P1_dn23.dn", sep="\t", header=T)')
#         r.pdf(dn_graphname)
#         robjects.r('p<-ggplot(melt(dn_table, "id"), aes(variable, value)) + geom_boxplot() + xlab("dinucleotide") + ylab("dinucleotide weight")')
#         robjects.r('print(p)')
#         robjects.r['dev.off']()
 
if 'CPB' in args.modules or 'all' in args.modules:
    #CPB, from Coleman et al 2008, pmid: 18583614

    #DNA to protein to CPB analysis
    #proteins = list((make_protein_record(nuc_rec) for nuc_rec in \
    #                 SeqIO.parse(filename, "fasta")))
    #SeqIO.write(proteins, filename+"_prot", "fasta")
    
    cpbname = seq_name[1:-1]+'_'+args.random_type+'.cpb'
    cpb_file=open(cpbname,'w')
    
    for nuc_rec in SeqIO.parse(filename, "fasta"):
        prot_rec = make_protein_record(nuc_rec)
        codons = [str(nuc_rec.seq[i:i+3]) for i in range(0,len(nuc_rec.seq)-3,3)]
        dicodons = [str(nuc_rec.seq[i:i+6]) for i in range(0,len(nuc_rec.seq)-6,6)]
        aas = [str(prot_rec.seq[i:i+1]) for i in range(0,len(prot_rec.seq),1)]
        diaas = [str(prot_rec.seq[i:i+2]) for i in range(0,len(prot_rec.seq)-1,1)]

        codon_counts = Counter(codons)
        dicodon_counts = Counter(dicodons)
        aa_counts = Counter(aas)
        diaa_counts = Counter(diaas)
        
        cps = []
        for cp in dicodons:
            cod1 = cp[:3]
            cod2 = cp[3:]
            if cod1 in ['TAG', 'TGA', 'TAA'] or cod2 in ['TAG', 'TGA', 'TAA']:
                continue #skip stop codons
            aa1 = tt[cod1].split('|')[0]
            aa2 = tt[cod2].split('|')[0]
            ap = aa1+aa2
            #score = np.log(dicodons.count(cp) / ( ((codons.count(cod1) * codons.count(cod2)) / (aas.count(aa1) * aas.count(aa2))) * diaas.count(ap)))
            score = np.log(dicodon_counts[cp] / ( ((codon_counts[cod1] * codon_counts[cod2]) / (aa_counts[aa1] * aa_counts[aa2])) * diaa_counts[ap]))
            cps.append(score)
        cpb = sum(cps)/len(cps)
        print str(cpb)
        cpb_file.write(str(cpb)+"\n")
    cpb_file.close()
    
    u_cols = ['cpb']
    cpb_table = pandas.read_csv(cpbname, sep=' ', names=u_cols)
    cpb_table_z =  ((cpb_table['cpb'] - cpb_table['cpb'].mean()) / cpb_table['cpb'].std())
    cpb_table_z_ls = (cpb_table_z-cpb_table_z[0])**2
    
    least_squares = least_squares.add(cpb_table_z_ls, axis=0)
#    final_table=final_table.append(cpb_table['cpb'])
    final_table.insert(0, "CPB", cpb_table['cpb'])
        
    if args.graphics:
        cpb_graphname = cpbname +'.pdf'
        u_cols = ['cpb']
        cpb_table = pandas.read_csv(cpbname, sep=' ', names=u_cols)
        cpb_graph = ggplot(cpb_table, aes('cpb'))+geom_density() +labs("CPB","Frequency")+ geom_vline(xintercept = [cpb_table['cpb'].iloc[0]] , colour="red", linetype = "dashed") +ggtitle(seq_name[1:-1]+'_'+args.random_type+' (CPB)') #Get the name of the script
        ggsave(cpb_graph, cpb_graphname)


#Calculation least square

least_squares = np.sqrt(least_squares)
least_squares.columns = ['distance']

least_table_name = seq_name[1:-1]+'_'+args.random_type+'_least_square.txt'
least_squares.to_csv(least_table_name, sep="\t")
#final_table=final_table.append(least_squares['distance'])
final_table.insert(0, "Distance(ls)", least_squares['distance'])
#least_table_file=open(least_table_name,'w')
#least_table_file.write(least_squares)
#least_table_file.close()
#least_table = sqrt()


# Create final graph and table 

nuc_distance_name=filename+'_distance_table.txt'
nuc_distance_file=open(nuc_distance_name,'w')
nuc_distance_table = np.zeros((n,n))
for j in range(1,n):
    difference = get_difference(seq_records[0], seq_records[j])
    nuc_distance_file.write(str(difference)+'\n')
    my_array[0, j] = difference
nuc_distance_file.close()

col_name = ['Nucleotide_difference']
new_nuc_distance_table = pandas.read_csv(nuc_distance_name, sep=' ', names=col_name)



new_nuc_distance_table.loc[-1]=[0]  
new_nuc_distance_table.index = new_nuc_distance_table.index + 1  
new_nuc_distance_table = new_nuc_distance_table.sort()
new_table=pandas.DataFrame()
new_table.insert(0, "Distance", least_squares['distance'])
new_table.insert(1, "Nucleotide_difference", new_nuc_distance_table['Nucleotide_difference'])

new_table_name = seq_name[1:-1]+'_'+args.random_type+'new_table_final_graph.txt'
new_table.to_csv(new_table_name, sep="\t")
  
#final_table=final_table.append(new_nuc_distance_table['Nucleotide_difference'])
final_table.insert(1, "Nucleotide_difference", new_nuc_distance_table['Nucleotide_difference'])

final_tb_name = seq_name[1:-1]+'_'+args.random_type+'_final_table.csv'
final_table.to_csv(final_tb_name, sep='\t')

if args.graphics:
    final_graphname = filename +'final_graph.pdf'
    final_graph = ggplot(new_table, aes('Distance', 'Nucleotide_difference'))+geom_point()+labs("Least Square Distance","Hamming Distance (nt)")+ggtitle(seq_name[1:-1]+'_'+args.random_type)
    ggsave(final_graph, final_graphname)

























