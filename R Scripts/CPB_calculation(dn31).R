library("plyr")
library("seqinr")

#Not counting stop codons


#Nucleotide
#dn231 <- read.fasta(file = "3dn_6deg/P1_3dn_6deg_with_original_seq.fas")                     
#n3 <- read.fasta(file = "3N/P1_3n_with_original_seq.fas")                     
#dn23 <- read.fasta(file = "dN23/P1_dn23_with_original_seq.fas")                     
dn31 <- read.fasta(file = "dN31/P1_dn31_with_original_seq.fas")                     

#AA
#dn231.aa <- read.fasta(file = "3dn_6deg/P1_3dn_6deg_with_original_seq_aa")                     
#n3.aa <- read.fasta(file = "3N/P1_3n_with_original_seq_aa.fas")                     
#dn23.aa <- read.fasta(file = "dN23/P1_dn23_with_original_seq_aa.fas")                     
dn31.aa <- read.fasta(file = "dN31/P1_dn31_with_original_seq_aa.fas")

#Count trinucleotide frequency in the sequence
#lapply(n3, function(x) count(x,3))
#Generate a Table with trinucleotide frequency 
#dn231.trinuc<-ldply(dn231, function(x) count(x,3)) 
#n3.trinuc<-ldply(n3, function(x) count(x,3)) 
#dn23.trinuc<-ldply(dn23, function(x) count(x,3)) 
dn31.trinuc<-ldply(dn31, function(x) count(x,3)) 

#Count frequency of codon pair(two trinucleotide) 
#dn231.codon.pair<-ldply(dn231, function(x) count(x,6)) 
#n3.codon.pair<-ldply(n3, function(x) count(x,6)) 
#dn23.codon.pair<-ldply(dn23, function(x) count(x,6)) 
dn31.codon.pair<-ldply(dn31, function(x) count(x,6))

#Generate a Table with AA frequency 
#dn231.count.aa<-ldply(dn231.aa, function(x) count(x,1, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 
#n3.count.aa<-ldply(n3.aa, function(x) count(x,1, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 
#dn23.count.aa<-ldply(dn23.aa, function(x) count(x,1, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 
dn31.count.aa<-ldply(dn31.aa, function(x) count(x,1, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 

#Generate a Table with diaa frequency 
#dn231.count.diaa<-ldply(dn231.aa, function(x) count(x,2, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 
#n3.count.diaa<-ldply(n3.aa, function(x) count(x,2)) 
#dn23.count.diaa<-ldply(dn23.aa, function(x) count(x,2)) 
dn31.count.diaa<-ldply(dn31.aa, function(x) count(x,2)) 

#CPS expected
#((F(A)xF(B) expected frequency)/F(X)xF(Y))xF(XY)
cps.aaa.aaa<-((dn31.trinuc$aaa*dn31.trinuc$aaa)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aaa.aac<-((dn31.trinuc$aaa*dn31.trinuc$aac)/(dn31.count.aa$k*dn31.count.aa$n))/dn31.count.diaa$kn
cps.aaa.aag<-((dn31.trinuc$aaa*dn31.trinuc$aag)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aaa.aat<-((dn31.trinuc$aaa*dn31.trinuc$aat)/(dn31.count.aa$k*dn31.count.aa$n))/dn31.count.diaa$kn

cps.aaa.aca<-((dn31.trinuc$aaa*dn31.trinuc$aca)/(dn31.count.aa$k*dn31.count.aa$t))/dn31.count.diaa$kt
cps.aaa.acc<-((dn31.trinuc$aaa*dn31.trinuc$acc)/(dn31.count.aa$k*dn31.count.aa$t))/dn31.count.diaa$kt
cps.aaa.acg<-((dn31.trinuc$aaa*dn31.trinuc$acg)/(dn31.count.aa$k*dn31.count.aa$t))/dn31.count.diaa$kt
cps.aaa.act<-((dn31.trinuc$aaa*dn31.trinuc$act)/(dn31.count.aa$k*dn31.count.aa$t))/dn31.count.diaa$kt

cps.aaa.aga<-((dn31.trinuc$aaa*dn31.trinuc$aga)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aaa.agc<-((dn31.trinuc$aaa*dn31.trinuc$agc)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks
cps.aaa.agg<-((dn31.trinuc$aaa*dn31.trinuc$agg)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aaa.agt<-((dn31.trinuc$aaa*dn31.trinuc$agt)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks

cps.aaa.ata<-((dn31.trinuc$aaa*dn31.trinuc$ata)/(dn31.count.aa$k*dn31.count.aa$i))/dn31.count.diaa$ki
cps.aaa.atc<-((dn31.trinuc$aaa*dn31.trinuc$atc)/(dn31.count.aa$k*dn31.count.aa$i))/dn31.count.diaa$ki
cps.aaa.atg<-((dn31.trinuc$aaa*dn31.trinuc$atg)/(dn31.count.aa$k*dn31.count.aa$m))/dn31.count.diaa$km
cps.aaa.att<-((dn31.trinuc$aaa*dn31.trinuc$att)/(dn31.count.aa$k*dn31.count.aa$i))/dn31.count.diaa$ki

cps.aaa.caa<-((dn31.trinuc$aaa*dn31.trinuc$caa)/(dn31.count.aa$k*dn31.count.aa$q))/dn31.count.diaa$kq
cps.aaa.cac<-((dn31.trinuc$aaa*dn31.trinuc$cac)/(dn31.count.aa$k*dn31.count.aa$h))/dn31.count.diaa$kh
cps.aaa.cag<-((dn31.trinuc$aaa*dn31.trinuc$cag)/(dn31.count.aa$k*dn31.count.aa$q))/dn31.count.diaa$kq
cps.aaa.cat<-((dn31.trinuc$aaa*dn31.trinuc$cat)/(dn31.count.aa$k*dn31.count.aa$h))/dn31.count.diaa$kh

cps.aaa.cca<-((dn31.trinuc$aaa*dn31.trinuc$cca)/(dn31.count.aa$k*dn31.count.aa$p))/dn31.count.diaa$kp
cps.aaa.ccc<-((dn31.trinuc$aaa*dn31.trinuc$ccc)/(dn31.count.aa$k*dn31.count.aa$p))/dn31.count.diaa$kp
cps.aaa.ccg<-((dn31.trinuc$aaa*dn31.trinuc$ccg)/(dn31.count.aa$k*dn31.count.aa$p))/dn31.count.diaa$kp
cps.aaa.cct<-((dn31.trinuc$aaa*dn31.trinuc$cct)/(dn31.count.aa$k*dn31.count.aa$p))/dn31.count.diaa$kp

cps.aaa.cga<-((dn31.trinuc$aaa*dn31.trinuc$cga)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aaa.cgc<-((dn31.trinuc$aaa*dn31.trinuc$cgc)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aaa.cgg<-((dn31.trinuc$aaa*dn31.trinuc$cgg)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aaa.cgt<-((dn31.trinuc$aaa*dn31.trinuc$cgt)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr

cps.aaa.cta<-((dn31.trinuc$aaa*dn31.trinuc$cta)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aaa.ctc<-((dn31.trinuc$aaa*dn31.trinuc$ctc)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aaa.ctg<-((dn31.trinuc$aaa*dn31.trinuc$ctg)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aaa.ctt<-((dn31.trinuc$aaa*dn31.trinuc$ctt)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl

cps.aaa.gaa<-((dn31.trinuc$aaa*dn31.trinuc$gaa)/(dn31.count.aa$k*dn31.count.aa$e))/dn31.count.diaa$ke
cps.aaa.gac<-((dn31.trinuc$aaa*dn31.trinuc$gac)/(dn31.count.aa$k*dn31.count.aa$d))/dn31.count.diaa$kd
cps.aaa.gag<-((dn31.trinuc$aaa*dn31.trinuc$gag)/(dn31.count.aa$k*dn31.count.aa$e))/dn31.count.diaa$ke
cps.aaa.gat<-((dn31.trinuc$aaa*dn31.trinuc$gat)/(dn31.count.aa$k*dn31.count.aa$d))/dn31.count.diaa$kd

cps.aaa.gca<-((dn31.trinuc$aaa*dn31.trinuc$gca)/(dn31.count.aa$k*dn31.count.aa$a))/dn31.count.diaa$ka
cps.aaa.gcc<-((dn31.trinuc$aaa*dn31.trinuc$gcc)/(dn31.count.aa$k*dn31.count.aa$a))/dn31.count.diaa$ka
cps.aaa.gcg<-((dn31.trinuc$aaa*dn31.trinuc$gcg)/(dn31.count.aa$k*dn31.count.aa$a))/dn31.count.diaa$ka
cps.aaa.gct<-((dn31.trinuc$aaa*dn31.trinuc$gct)/(dn31.count.aa$k*dn31.count.aa$a))/dn31.count.diaa$ka

cps.aaa.gga<-((dn31.trinuc$aaa*dn31.trinuc$gga)/(dn31.count.aa$k*dn31.count.aa$g))/dn31.count.diaa$kg
cps.aaa.ggc<-((dn31.trinuc$aaa*dn31.trinuc$ggc)/(dn31.count.aa$k*dn31.count.aa$g))/dn31.count.diaa$kg
cps.aaa.ggg<-((dn31.trinuc$aaa*dn31.trinuc$ggg)/(dn31.count.aa$k*dn31.count.aa$g))/dn31.count.diaa$kg
cps.aaa.ggt<-((dn31.trinuc$aaa*dn31.trinuc$ggt)/(dn31.count.aa$k*dn31.count.aa$g))/dn31.count.diaa$kg

cps.aaa.gta<-((dn31.trinuc$aaa*dn31.trinuc$gta)/(dn31.count.aa$k*dn31.count.aa$v))/dn31.count.diaa$kv
cps.aaa.gtc<-((dn31.trinuc$aaa*dn31.trinuc$gtc)/(dn31.count.aa$k*dn31.count.aa$v))/dn31.count.diaa$kv
cps.aaa.gtg<-((dn31.trinuc$aaa*dn31.trinuc$gtg)/(dn31.count.aa$k*dn31.count.aa$v))/dn31.count.diaa$kv
cps.aaa.gtt<-((dn31.trinuc$aaa*dn31.trinuc$gtt)/(dn31.count.aa$k*dn31.count.aa$v))/dn31.count.diaa$kv

#stop codon
#cps.aaa.taa<-((dn31.trinuc$aaa*dn31.trinuc$taa)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aaa.tac<-((dn31.trinuc$aaa*dn31.trinuc$tac)/(dn31.count.aa$k*dn31.count.aa$y))/dn31.count.diaa$ky
#stop codon
#cps.aaa.tag<-((dn31.trinuc$aaa*dn31.trinuc$tag)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aaa.tat<-((dn31.trinuc$aaa*dn31.trinuc$tat)/(dn31.count.aa$k*dn31.count.aa$y))/dn31.count.diaa$ky

cps.aaa.tca<-((dn31.trinuc$aaa*dn31.trinuc$tca)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks
cps.aaa.tcc<-((dn31.trinuc$aaa*dn31.trinuc$tcc)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks
cps.aaa.tcg<-((dn31.trinuc$aaa*dn31.trinuc$tcg)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks
cps.aaa.tct<-((dn31.trinuc$aaa*dn31.trinuc$tct)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks

#Stop codon
#cps.aaa.tga<-((dn31.trinuc$aaa*dn31.trinuc$tga)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aaa.tgc<-((dn31.trinuc$aaa*dn31.trinuc$tgc)/(dn31.count.aa$k*dn31.count.aa$c))/dn31.count.diaa$kc
cps.aaa.tgg<-((dn31.trinuc$aaa*dn31.trinuc$tgg)/(dn31.count.aa$k*dn31.count.aa$w))/dn31.count.diaa$kw
cps.aaa.tgt<-((dn31.trinuc$aaa*dn31.trinuc$tgt)/(dn31.count.aa$k*dn31.count.aa$c))/dn31.count.diaa$kc

cps.aaa.tta<-((dn31.trinuc$aaa*dn31.trinuc$tta)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aaa.ttc<-((dn31.trinuc$aaa*dn31.trinuc$ttc)/(dn31.count.aa$k*dn31.count.aa$f))/dn31.count.diaa$kf
cps.aaa.ttg<-((dn31.trinuc$aaa*dn31.trinuc$ttg)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aaa.ttt<-((dn31.trinuc$aaa*dn31.trinuc$ttt)/(dn31.count.aa$k*dn31.count.aa$f))/dn31.count.diaa$kf




cps.aac.aaa<-((dn31.trinuc$aac*dn31.trinuc$aaa)/(dn31.count.aa$n*dn31.count.aa$k))/dn31.count.diaa$nk
cps.aac.aac<-((dn31.trinuc$aac*dn31.trinuc$aac)/(dn31.count.aa$n*dn31.count.aa$n))/dn31.count.diaa$nn
cps.aac.aag<-((dn31.trinuc$aac*dn31.trinuc$aag)/(dn31.count.aa$n*dn31.count.aa$k))/dn31.count.diaa$nk
cps.aac.aat<-((dn31.trinuc$aac*dn31.trinuc$aat)/(dn31.count.aa$n*dn31.count.aa$n))/dn31.count.diaa$nn

cps.aac.aca<-((dn31.trinuc$aac*dn31.trinuc$aca)/(dn31.count.aa$n*dn31.count.aa$t))/dn31.count.diaa$nt
cps.aac.acc<-((dn31.trinuc$aac*dn31.trinuc$acc)/(dn31.count.aa$n*dn31.count.aa$t))/dn31.count.diaa$nt
cps.aac.acg<-((dn31.trinuc$aac*dn31.trinuc$acg)/(dn31.count.aa$n*dn31.count.aa$t))/dn31.count.diaa$nt
cps.aac.act<-((dn31.trinuc$aac*dn31.trinuc$act)/(dn31.count.aa$n*dn31.count.aa$t))/dn31.count.diaa$nt

cps.aac.aga<-((dn31.trinuc$aac*dn31.trinuc$aga)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aac.agc<-((dn31.trinuc$aac*dn31.trinuc$agc)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns
cps.aac.agg<-((dn31.trinuc$aac*dn31.trinuc$agg)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aac.agt<-((dn31.trinuc$aac*dn31.trinuc$agt)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns

cps.aac.ata<-((dn31.trinuc$aac*dn31.trinuc$ata)/(dn31.count.aa$n*dn31.count.aa$i))/dn31.count.diaa$ni
cps.aac.atc<-((dn31.trinuc$aac*dn31.trinuc$atc)/(dn31.count.aa$n*dn31.count.aa$i))/dn31.count.diaa$ni
cps.aac.atg<-((dn31.trinuc$aac*dn31.trinuc$atg)/(dn31.count.aa$n*dn31.count.aa$m))/dn31.count.diaa$nm
cps.aac.att<-((dn31.trinuc$aac*dn31.trinuc$att)/(dn31.count.aa$n*dn31.count.aa$i))/dn31.count.diaa$ni

cps.aac.caa<-((dn31.trinuc$aac*dn31.trinuc$caa)/(dn31.count.aa$n*dn31.count.aa$q))/dn31.count.diaa$nq
cps.aac.cac<-((dn31.trinuc$aac*dn31.trinuc$cac)/(dn31.count.aa$n*dn31.count.aa$h))/dn31.count.diaa$nh
cps.aac.cag<-((dn31.trinuc$aac*dn31.trinuc$cag)/(dn31.count.aa$n*dn31.count.aa$q))/dn31.count.diaa$nq
cps.aac.cat<-((dn31.trinuc$aac*dn31.trinuc$cat)/(dn31.count.aa$n*dn31.count.aa$h))/dn31.count.diaa$nh

cps.aac.cca<-((dn31.trinuc$aac*dn31.trinuc$cca)/(dn31.count.aa$n*dn31.count.aa$p))/dn31.count.diaa$np
cps.aac.ccc<-((dn31.trinuc$aac*dn31.trinuc$ccc)/(dn31.count.aa$n*dn31.count.aa$p))/dn31.count.diaa$np
cps.aac.ccg<-((dn31.trinuc$aac*dn31.trinuc$ccg)/(dn31.count.aa$n*dn31.count.aa$p))/dn31.count.diaa$np
cps.aac.cct<-((dn31.trinuc$aac*dn31.trinuc$cct)/(dn31.count.aa$n*dn31.count.aa$p))/dn31.count.diaa$np

cps.aac.cga<-((dn31.trinuc$aac*dn31.trinuc$cga)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aac.cgc<-((dn31.trinuc$aac*dn31.trinuc$cgc)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aac.cgg<-((dn31.trinuc$aac*dn31.trinuc$cgg)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aac.cgt<-((dn31.trinuc$aac*dn31.trinuc$cgt)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr

cps.aac.cta<-((dn31.trinuc$aac*dn31.trinuc$cta)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aac.ctc<-((dn31.trinuc$aac*dn31.trinuc$ctc)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aac.ctg<-((dn31.trinuc$aac*dn31.trinuc$ctg)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aac.ctt<-((dn31.trinuc$aac*dn31.trinuc$ctt)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl

cps.aac.gaa<-((dn31.trinuc$aac*dn31.trinuc$gaa)/(dn31.count.aa$n*dn31.count.aa$e))/dn31.count.diaa$ne
cps.aac.gac<-((dn31.trinuc$aac*dn31.trinuc$gac)/(dn31.count.aa$n*dn31.count.aa$d))/dn31.count.diaa$nd
cps.aac.gag<-((dn31.trinuc$aac*dn31.trinuc$gag)/(dn31.count.aa$n*dn31.count.aa$e))/dn31.count.diaa$ne
cps.aac.gat<-((dn31.trinuc$aac*dn31.trinuc$gat)/(dn31.count.aa$n*dn31.count.aa$d))/dn31.count.diaa$nd

cps.aac.gca<-((dn31.trinuc$aac*dn31.trinuc$gca)/(dn31.count.aa$n*dn31.count.aa$a))/dn31.count.diaa$na
cps.aac.gcc<-((dn31.trinuc$aac*dn31.trinuc$gcc)/(dn31.count.aa$n*dn31.count.aa$a))/dn31.count.diaa$na
cps.aac.gcg<-((dn31.trinuc$aac*dn31.trinuc$gcg)/(dn31.count.aa$n*dn31.count.aa$a))/dn31.count.diaa$na
cps.aac.gct<-((dn31.trinuc$aac*dn31.trinuc$gct)/(dn31.count.aa$n*dn31.count.aa$a))/dn31.count.diaa$na

cps.aac.gga<-((dn31.trinuc$aac*dn31.trinuc$gga)/(dn31.count.aa$n*dn31.count.aa$g))/dn31.count.diaa$ng
cps.aac.ggc<-((dn31.trinuc$aac*dn31.trinuc$ggc)/(dn31.count.aa$n*dn31.count.aa$g))/dn31.count.diaa$ng
cps.aac.ggg<-((dn31.trinuc$aac*dn31.trinuc$ggg)/(dn31.count.aa$n*dn31.count.aa$g))/dn31.count.diaa$ng
cps.aac.ggt<-((dn31.trinuc$aac*dn31.trinuc$ggt)/(dn31.count.aa$n*dn31.count.aa$g))/dn31.count.diaa$ng

cps.aac.gta<-((dn31.trinuc$aac*dn31.trinuc$gta)/(dn31.count.aa$n*dn31.count.aa$v))/dn31.count.diaa$nv
cps.aac.gtc<-((dn31.trinuc$aac*dn31.trinuc$gtc)/(dn31.count.aa$n*dn31.count.aa$v))/dn31.count.diaa$nv
cps.aac.gtg<-((dn31.trinuc$aac*dn31.trinuc$gtg)/(dn31.count.aa$n*dn31.count.aa$v))/dn31.count.diaa$nv
cps.aac.gtt<-((dn31.trinuc$aac*dn31.trinuc$gtt)/(dn31.count.aa$n*dn31.count.aa$v))/dn31.count.diaa$nv

#Stop codon
#cps.aac.taa<-((dn31.trinuc$aac*dn31.trinuc$taa)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aac.tac<-((dn31.trinuc$aac*dn31.trinuc$tac)/(dn31.count.aa$n*dn31.count.aa$y))/dn31.count.diaa$ny
#Stop codon
#cps.aac.tag<-((dn31.trinuc$aac*dn31.trinuc$tag)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aac.tat<-((dn31.trinuc$aac*dn31.trinuc$tat)/(dn31.count.aa$n*dn31.count.aa$y))/dn31.count.diaa$ny

cps.aac.tca<-((dn31.trinuc$aac*dn31.trinuc$tca)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns
cps.aac.tcc<-((dn31.trinuc$aac*dn31.trinuc$tcc)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns
cps.aac.tcg<-((dn31.trinuc$aac*dn31.trinuc$tcg)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns
cps.aac.tct<-((dn31.trinuc$aac*dn31.trinuc$tct)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns

#Stop codon
#cps.aac.tga<-((dn31.trinuc$aac*dn31.trinuc$tga)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aac.tgc<-((dn31.trinuc$aac*dn31.trinuc$tgc)/(dn31.count.aa$n*dn31.count.aa$c))/dn31.count.diaa$nc
cps.aac.tgg<-((dn31.trinuc$aac*dn31.trinuc$tgg)/(dn31.count.aa$n*dn31.count.aa$w))/dn31.count.diaa$nw
cps.aac.tgt<-((dn31.trinuc$aac*dn31.trinuc$tgt)/(dn31.count.aa$n*dn31.count.aa$c))/dn31.count.diaa$nc

cps.aac.tta<-((dn31.trinuc$aac*dn31.trinuc$tta)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aac.ttc<-((dn31.trinuc$aac*dn31.trinuc$ttc)/(dn31.count.aa$n*dn31.count.aa$f))/dn31.count.diaa$nf
cps.aac.ttg<-((dn31.trinuc$aac*dn31.trinuc$ttg)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aac.ttt<-((dn31.trinuc$aac*dn31.trinuc$ttt)/(dn31.count.aa$n*dn31.count.aa$f))/dn31.count.diaa$nf





cps.aag.aaa<-((dn31.trinuc$aag*dn31.trinuc$aaa)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aag.aac<-((dn31.trinuc$aag*dn31.trinuc$aac)/(dn31.count.aa$k*dn31.count.aa$n))/dn31.count.diaa$kn
cps.aag.aag<-((dn31.trinuc$aag*dn31.trinuc$aag)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aag.aat<-((dn31.trinuc$aag*dn31.trinuc$aat)/(dn31.count.aa$k*dn31.count.aa$n))/dn31.count.diaa$kn

cps.aag.aca<-((dn31.trinuc$aag*dn31.trinuc$aca)/(dn31.count.aa$k*dn31.count.aa$t))/dn31.count.diaa$kt
cps.aag.acc<-((dn31.trinuc$aag*dn31.trinuc$acc)/(dn31.count.aa$k*dn31.count.aa$t))/dn31.count.diaa$kt
cps.aag.acg<-((dn31.trinuc$aag*dn31.trinuc$acg)/(dn31.count.aa$k*dn31.count.aa$t))/dn31.count.diaa$kt
cps.aag.act<-((dn31.trinuc$aag*dn31.trinuc$act)/(dn31.count.aa$k*dn31.count.aa$t))/dn31.count.diaa$kt

cps.aag.aga<-((dn31.trinuc$aag*dn31.trinuc$aga)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aag.agc<-((dn31.trinuc$aag*dn31.trinuc$agc)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks
cps.aag.agg<-((dn31.trinuc$aag*dn31.trinuc$agg)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aag.agt<-((dn31.trinuc$aag*dn31.trinuc$agt)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks

cps.aag.ata<-((dn31.trinuc$aag*dn31.trinuc$ata)/(dn31.count.aa$k*dn31.count.aa$i))/dn31.count.diaa$ki
cps.aag.atc<-((dn31.trinuc$aag*dn31.trinuc$atc)/(dn31.count.aa$k*dn31.count.aa$i))/dn31.count.diaa$ki
cps.aag.atg<-((dn31.trinuc$aag*dn31.trinuc$atg)/(dn31.count.aa$k*dn31.count.aa$m))/dn31.count.diaa$km
cps.aag.att<-((dn31.trinuc$aag*dn31.trinuc$att)/(dn31.count.aa$k*dn31.count.aa$i))/dn31.count.diaa$ki

cps.aag.caa<-((dn31.trinuc$aag*dn31.trinuc$caa)/(dn31.count.aa$k*dn31.count.aa$q))/dn31.count.diaa$kq
cps.aag.cac<-((dn31.trinuc$aag*dn31.trinuc$cac)/(dn31.count.aa$k*dn31.count.aa$h))/dn31.count.diaa$kh
cps.aag.cag<-((dn31.trinuc$aag*dn31.trinuc$cag)/(dn31.count.aa$k*dn31.count.aa$q))/dn31.count.diaa$kq
cps.aag.cat<-((dn31.trinuc$aag*dn31.trinuc$cat)/(dn31.count.aa$k*dn31.count.aa$h))/dn31.count.diaa$kh

cps.aag.cca<-((dn31.trinuc$aag*dn31.trinuc$cca)/(dn31.count.aa$k*dn31.count.aa$p))/dn31.count.diaa$kp
cps.aag.ccc<-((dn31.trinuc$aag*dn31.trinuc$ccc)/(dn31.count.aa$k*dn31.count.aa$p))/dn31.count.diaa$kp
cps.aag.ccg<-((dn31.trinuc$aag*dn31.trinuc$ccg)/(dn31.count.aa$k*dn31.count.aa$p))/dn31.count.diaa$kp
cps.aag.cct<-((dn31.trinuc$aag*dn31.trinuc$cct)/(dn31.count.aa$k*dn31.count.aa$p))/dn31.count.diaa$kp

cps.aag.cga<-((dn31.trinuc$aag*dn31.trinuc$cga)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aag.cgc<-((dn31.trinuc$aag*dn31.trinuc$cgc)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aag.cgg<-((dn31.trinuc$aag*dn31.trinuc$cgg)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr
cps.aag.cgt<-((dn31.trinuc$aag*dn31.trinuc$cgt)/(dn31.count.aa$k*dn31.count.aa$r))/dn31.count.diaa$kr

cps.aag.cta<-((dn31.trinuc$aag*dn31.trinuc$cta)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aag.ctc<-((dn31.trinuc$aag*dn31.trinuc$ctc)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aag.ctg<-((dn31.trinuc$aag*dn31.trinuc$ctg)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aag.ctt<-((dn31.trinuc$aag*dn31.trinuc$ctt)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl

cps.aag.gaa<-((dn31.trinuc$aag*dn31.trinuc$gaa)/(dn31.count.aa$k*dn31.count.aa$e))/dn31.count.diaa$ke
cps.aag.gac<-((dn31.trinuc$aag*dn31.trinuc$gac)/(dn31.count.aa$k*dn31.count.aa$d))/dn31.count.diaa$kd
cps.aag.gag<-((dn31.trinuc$aag*dn31.trinuc$gag)/(dn31.count.aa$k*dn31.count.aa$e))/dn31.count.diaa$ke
cps.aag.gat<-((dn31.trinuc$aag*dn31.trinuc$gat)/(dn31.count.aa$k*dn31.count.aa$d))/dn31.count.diaa$kd

cps.aag.gca<-((dn31.trinuc$aag*dn31.trinuc$gca)/(dn31.count.aa$k*dn31.count.aa$a))/dn31.count.diaa$ka
cps.aag.gcc<-((dn31.trinuc$aag*dn31.trinuc$gcc)/(dn31.count.aa$k*dn31.count.aa$a))/dn31.count.diaa$ka
cps.aag.gcg<-((dn31.trinuc$aag*dn31.trinuc$gcg)/(dn31.count.aa$k*dn31.count.aa$a))/dn31.count.diaa$ka
cps.aag.gct<-((dn31.trinuc$aag*dn31.trinuc$gct)/(dn31.count.aa$k*dn31.count.aa$a))/dn31.count.diaa$ka

cps.aag.gga<-((dn31.trinuc$aag*dn31.trinuc$gga)/(dn31.count.aa$k*dn31.count.aa$g))/dn31.count.diaa$kg
cps.aag.ggc<-((dn31.trinuc$aag*dn31.trinuc$ggc)/(dn31.count.aa$k*dn31.count.aa$g))/dn31.count.diaa$kg
cps.aag.ggg<-((dn31.trinuc$aag*dn31.trinuc$ggg)/(dn31.count.aa$k*dn31.count.aa$g))/dn31.count.diaa$kg
cps.aag.ggt<-((dn31.trinuc$aag*dn31.trinuc$ggt)/(dn31.count.aa$k*dn31.count.aa$g))/dn31.count.diaa$kg

cps.aag.gta<-((dn31.trinuc$aag*dn31.trinuc$gta)/(dn31.count.aa$k*dn31.count.aa$v))/dn31.count.diaa$kv
cps.aag.gtc<-((dn31.trinuc$aag*dn31.trinuc$gtc)/(dn31.count.aa$k*dn31.count.aa$v))/dn31.count.diaa$kv
cps.aag.gtg<-((dn31.trinuc$aag*dn31.trinuc$gtg)/(dn31.count.aa$k*dn31.count.aa$v))/dn31.count.diaa$kv
cps.aag.gtt<-((dn31.trinuc$aag*dn31.trinuc$gtt)/(dn31.count.aa$k*dn31.count.aa$v))/dn31.count.diaa$kv

#Stop codon
#cps.aag.taa<-((dn31.trinuc$aag*dn31.trinuc$taa)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aag.tac<-((dn31.trinuc$aag*dn31.trinuc$tac)/(dn31.count.aa$k*dn31.count.aa$y))/dn31.count.diaa$ky
#Stop codon
#cps.aag.tag<-((dn31.trinuc$aag*dn31.trinuc$tag)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aag.tat<-((dn31.trinuc$aag*dn31.trinuc$tat)/(dn31.count.aa$k*dn31.count.aa$y))/dn31.count.diaa$ky

cps.aag.tca<-((dn31.trinuc$aag*dn31.trinuc$tca)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks
cps.aag.tcc<-((dn31.trinuc$aag*dn31.trinuc$tcc)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks
cps.aag.tcg<-((dn31.trinuc$aag*dn31.trinuc$tcg)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks
cps.aag.tct<-((dn31.trinuc$aag*dn31.trinuc$tct)/(dn31.count.aa$k*dn31.count.aa$s))/dn31.count.diaa$ks

#Stop codon
#cps.aag.tga<-((dn31.trinuc$aag*dn31.trinuc$tga)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aag.tgc<-((dn31.trinuc$aag*dn31.trinuc$tgc)/(dn31.count.aa$k*dn31.count.aa$c))/dn31.count.diaa$kc
cps.aag.tgg<-((dn31.trinuc$aag*dn31.trinuc$tgg)/(dn31.count.aa$k*dn31.count.aa$w))/dn31.count.diaa$kw
cps.aag.tgt<-((dn31.trinuc$aag*dn31.trinuc$tgt)/(dn31.count.aa$k*dn31.count.aa$c))/dn31.count.diaa$kc

cps.aag.tta<-((dn31.trinuc$aag*dn31.trinuc$tta)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aag.ttc<-((dn31.trinuc$aag*dn31.trinuc$ttc)/(dn31.count.aa$k*dn31.count.aa$f))/dn31.count.diaa$kf
cps.aag.ttg<-((dn31.trinuc$aag*dn31.trinuc$ttg)/(dn31.count.aa$k*dn31.count.aa$l))/dn31.count.diaa$kl
cps.aag.ttt<-((dn31.trinuc$aag*dn31.trinuc$ttt)/(dn31.count.aa$k*dn31.count.aa$f))/dn31.count.diaa$kf








cps.aat.aaa<-((dn31.trinuc$aat*dn31.trinuc$aaa)/(dn31.count.aa$n*dn31.count.aa$k))/dn31.count.diaa$nk
cps.aat.aac<-((dn31.trinuc$aat*dn31.trinuc$aac)/(dn31.count.aa$n*dn31.count.aa$n))/dn31.count.diaa$nn
cps.aat.aag<-((dn31.trinuc$aat*dn31.trinuc$aag)/(dn31.count.aa$n*dn31.count.aa$k))/dn31.count.diaa$nk
cps.aat.aat<-((dn31.trinuc$aat*dn31.trinuc$aat)/(dn31.count.aa$n*dn31.count.aa$n))/dn31.count.diaa$nn

cps.aat.aca<-((dn31.trinuc$aat*dn31.trinuc$aca)/(dn31.count.aa$n*dn31.count.aa$t))/dn31.count.diaa$nt
cps.aat.acc<-((dn31.trinuc$aat*dn31.trinuc$acc)/(dn31.count.aa$n*dn31.count.aa$t))/dn31.count.diaa$nt
cps.aat.acg<-((dn31.trinuc$aat*dn31.trinuc$acg)/(dn31.count.aa$n*dn31.count.aa$t))/dn31.count.diaa$nt
cps.aat.act<-((dn31.trinuc$aat*dn31.trinuc$act)/(dn31.count.aa$n*dn31.count.aa$t))/dn31.count.diaa$nt

cps.aat.aga<-((dn31.trinuc$aat*dn31.trinuc$aga)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aat.agc<-((dn31.trinuc$aat*dn31.trinuc$agc)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns
cps.aat.agg<-((dn31.trinuc$aat*dn31.trinuc$agg)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aat.agt<-((dn31.trinuc$aat*dn31.trinuc$agt)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns

cps.aat.ata<-((dn31.trinuc$aat*dn31.trinuc$ata)/(dn31.count.aa$n*dn31.count.aa$i))/dn31.count.diaa$ni
cps.aat.atc<-((dn31.trinuc$aat*dn31.trinuc$atc)/(dn31.count.aa$n*dn31.count.aa$i))/dn31.count.diaa$ni
cps.aat.atg<-((dn31.trinuc$aat*dn31.trinuc$atg)/(dn31.count.aa$n*dn31.count.aa$m))/dn31.count.diaa$nm
cps.aat.att<-((dn31.trinuc$aat*dn31.trinuc$att)/(dn31.count.aa$n*dn31.count.aa$i))/dn31.count.diaa$ni

cps.aat.caa<-((dn31.trinuc$aat*dn31.trinuc$caa)/(dn31.count.aa$n*dn31.count.aa$q))/dn31.count.diaa$nq
cps.aat.cac<-((dn31.trinuc$aat*dn31.trinuc$cac)/(dn31.count.aa$n*dn31.count.aa$h))/dn31.count.diaa$nh
cps.aat.cag<-((dn31.trinuc$aat*dn31.trinuc$cag)/(dn31.count.aa$n*dn31.count.aa$q))/dn31.count.diaa$nq
cps.aat.cat<-((dn31.trinuc$aat*dn31.trinuc$cat)/(dn31.count.aa$n*dn31.count.aa$h))/dn31.count.diaa$nh

cps.aat.cca<-((dn31.trinuc$aat*dn31.trinuc$cca)/(dn31.count.aa$n*dn31.count.aa$p))/dn31.count.diaa$np
cps.aat.ccc<-((dn31.trinuc$aat*dn31.trinuc$ccc)/(dn31.count.aa$n*dn31.count.aa$p))/dn31.count.diaa$np
cps.aat.ccg<-((dn31.trinuc$aat*dn31.trinuc$ccg)/(dn31.count.aa$n*dn31.count.aa$p))/dn31.count.diaa$np
cps.aat.cct<-((dn31.trinuc$aat*dn31.trinuc$cct)/(dn31.count.aa$n*dn31.count.aa$p))/dn31.count.diaa$np

cps.aat.cga<-((dn31.trinuc$aat*dn31.trinuc$cga)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aat.cgc<-((dn31.trinuc$aat*dn31.trinuc$cgc)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aat.cgg<-((dn31.trinuc$aat*dn31.trinuc$cgg)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr
cps.aat.cgt<-((dn31.trinuc$aat*dn31.trinuc$cgt)/(dn31.count.aa$n*dn31.count.aa$r))/dn31.count.diaa$nr

cps.aat.cta<-((dn31.trinuc$aat*dn31.trinuc$cta)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aat.ctc<-((dn31.trinuc$aat*dn31.trinuc$ctc)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aat.ctg<-((dn31.trinuc$aat*dn31.trinuc$ctg)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aat.ctt<-((dn31.trinuc$aat*dn31.trinuc$ctt)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl

cps.aat.gaa<-((dn31.trinuc$aat*dn31.trinuc$gaa)/(dn31.count.aa$n*dn31.count.aa$e))/dn31.count.diaa$ne
cps.aat.gac<-((dn31.trinuc$aat*dn31.trinuc$gac)/(dn31.count.aa$n*dn31.count.aa$d))/dn31.count.diaa$nd
cps.aat.gag<-((dn31.trinuc$aat*dn31.trinuc$gag)/(dn31.count.aa$n*dn31.count.aa$e))/dn31.count.diaa$ne
cps.aat.gat<-((dn31.trinuc$aat*dn31.trinuc$gat)/(dn31.count.aa$n*dn31.count.aa$d))/dn31.count.diaa$nd

cps.aat.gca<-((dn31.trinuc$aat*dn31.trinuc$gca)/(dn31.count.aa$n*dn31.count.aa$a))/dn31.count.diaa$na
cps.aat.gcc<-((dn31.trinuc$aat*dn31.trinuc$gcc)/(dn31.count.aa$n*dn31.count.aa$a))/dn31.count.diaa$na
cps.aat.gcg<-((dn31.trinuc$aat*dn31.trinuc$gcg)/(dn31.count.aa$n*dn31.count.aa$a))/dn31.count.diaa$na
cps.aat.gct<-((dn31.trinuc$aat*dn31.trinuc$gct)/(dn31.count.aa$n*dn31.count.aa$a))/dn31.count.diaa$na

cps.aat.gga<-((dn31.trinuc$aat*dn31.trinuc$gga)/(dn31.count.aa$n*dn31.count.aa$g))/dn31.count.diaa$ng
cps.aat.ggc<-((dn31.trinuc$aat*dn31.trinuc$ggc)/(dn31.count.aa$n*dn31.count.aa$g))/dn31.count.diaa$ng
cps.aat.ggg<-((dn31.trinuc$aat*dn31.trinuc$ggg)/(dn31.count.aa$n*dn31.count.aa$g))/dn31.count.diaa$ng
cps.aat.ggt<-((dn31.trinuc$aat*dn31.trinuc$ggt)/(dn31.count.aa$n*dn31.count.aa$g))/dn31.count.diaa$ng

cps.aat.gta<-((dn31.trinuc$aat*dn31.trinuc$gta)/(dn31.count.aa$n*dn31.count.aa$v))/dn31.count.diaa$nv
cps.aat.gtc<-((dn31.trinuc$aat*dn31.trinuc$gtc)/(dn31.count.aa$n*dn31.count.aa$v))/dn31.count.diaa$nv
cps.aat.gtg<-((dn31.trinuc$aat*dn31.trinuc$gtg)/(dn31.count.aa$n*dn31.count.aa$v))/dn31.count.diaa$nv
cps.aat.gtt<-((dn31.trinuc$aat*dn31.trinuc$gtt)/(dn31.count.aa$n*dn31.count.aa$v))/dn31.count.diaa$nv

#Stop codon
#cps.aat.taa<-((dn31.trinuc$aat*dn31.trinuc$taa)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aat.tac<-((dn31.trinuc$aat*dn31.trinuc$tac)/(dn31.count.aa$n*dn31.count.aa$y))/dn31.count.diaa$ny
#Stop codon
#cps.aat.tag<-((dn31.trinuc$aat*dn31.trinuc$tag)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aat.tat<-((dn31.trinuc$aat*dn31.trinuc$tat)/(dn31.count.aa$n*dn31.count.aa$y))/dn31.count.diaa$ny

cps.aat.tca<-((dn31.trinuc$aat*dn31.trinuc$tca)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns
cps.aat.tcc<-((dn31.trinuc$aat*dn31.trinuc$tcc)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns
cps.aat.tcg<-((dn31.trinuc$aat*dn31.trinuc$tcg)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns
cps.aat.tct<-((dn31.trinuc$aat*dn31.trinuc$tct)/(dn31.count.aa$n*dn31.count.aa$s))/dn31.count.diaa$ns

#Stop codon
#cps.aat.tga<-((dn31.trinuc$aat*dn31.trinuc$tga)/(dn31.count.aa$k*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aat.tgc<-((dn31.trinuc$aat*dn31.trinuc$tgc)/(dn31.count.aa$n*dn31.count.aa$c))/dn31.count.diaa$nc
cps.aat.tgg<-((dn31.trinuc$aat*dn31.trinuc$tgg)/(dn31.count.aa$n*dn31.count.aa$w))/dn31.count.diaa$nw
cps.aat.tgt<-((dn31.trinuc$aat*dn31.trinuc$tgt)/(dn31.count.aa$n*dn31.count.aa$c))/dn31.count.diaa$nc

cps.aat.tta<-((dn31.trinuc$aat*dn31.trinuc$tta)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aat.ttc<-((dn31.trinuc$aat*dn31.trinuc$ttc)/(dn31.count.aa$n*dn31.count.aa$f))/dn31.count.diaa$nf
cps.aat.ttg<-((dn31.trinuc$aat*dn31.trinuc$ttg)/(dn31.count.aa$n*dn31.count.aa$l))/dn31.count.diaa$nl
cps.aat.ttt<-((dn31.trinuc$aat*dn31.trinuc$ttt)/(dn31.count.aa$n*dn31.count.aa$f))/dn31.count.diaa$nf






cps.aca.aaa<-((dn31.trinuc$aca*dn31.trinuc$aaa)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$tk
cps.aca.aac<-((dn31.trinuc$aca*dn31.trinuc$aac)/(dn31.count.aa$t*dn31.count.aa$n))/dn31.count.diaa$tn
cps.aca.aag<-((dn31.trinuc$aca*dn31.trinuc$aag)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$tk
cps.aca.aat<-((dn31.trinuc$aca*dn31.trinuc$aat)/(dn31.count.aa$t*dn31.count.aa$n))/dn31.count.diaa$tn

cps.aca.aca<-((dn31.trinuc$aca*dn31.trinuc$aca)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.aca.acc<-((dn31.trinuc$aca*dn31.trinuc$acc)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.aca.acg<-((dn31.trinuc$aca*dn31.trinuc$acg)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.aca.act<-((dn31.trinuc$aca*dn31.trinuc$act)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt

cps.aca.aga<-((dn31.trinuc$aca*dn31.trinuc$aga)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.aca.agc<-((dn31.trinuc$aca*dn31.trinuc$agc)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.aca.agg<-((dn31.trinuc$aca*dn31.trinuc$agg)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.aca.agt<-((dn31.trinuc$aca*dn31.trinuc$agt)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts

cps.aca.ata<-((dn31.trinuc$aca*dn31.trinuc$ata)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti
cps.aca.atc<-((dn31.trinuc$aca*dn31.trinuc$atc)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti
cps.aca.atg<-((dn31.trinuc$aca*dn31.trinuc$atg)/(dn31.count.aa$t*dn31.count.aa$m))/dn31.count.diaa$tm
cps.aca.att<-((dn31.trinuc$aca*dn31.trinuc$att)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti

cps.aca.caa<-((dn31.trinuc$aca*dn31.trinuc$caa)/(dn31.count.aa$t*dn31.count.aa$q))/dn31.count.diaa$tq
cps.aca.cac<-((dn31.trinuc$aca*dn31.trinuc$cac)/(dn31.count.aa$t*dn31.count.aa$h))/dn31.count.diaa$th
cps.aca.cag<-((dn31.trinuc$aca*dn31.trinuc$cag)/(dn31.count.aa$t*dn31.count.aa$q))/dn31.count.diaa$tq
cps.aca.cat<-((dn31.trinuc$aca*dn31.trinuc$cat)/(dn31.count.aa$t*dn31.count.aa$h))/dn31.count.diaa$th

cps.aca.cca<-((dn31.trinuc$aca*dn31.trinuc$cca)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.aca.ccc<-((dn31.trinuc$aca*dn31.trinuc$ccc)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.aca.ccg<-((dn31.trinuc$aca*dn31.trinuc$ccg)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.aca.cct<-((dn31.trinuc$aca*dn31.trinuc$cct)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp

cps.aca.cga<-((dn31.trinuc$aca*dn31.trinuc$cga)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.aca.cgc<-((dn31.trinuc$aca*dn31.trinuc$cgc)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.aca.cgg<-((dn31.trinuc$aca*dn31.trinuc$cgg)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.aca.cgt<-((dn31.trinuc$aca*dn31.trinuc$cgt)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr

cps.aca.cta<-((dn31.trinuc$aca*dn31.trinuc$cta)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.aca.ctc<-((dn31.trinuc$aca*dn31.trinuc$ctc)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.aca.ctg<-((dn31.trinuc$aca*dn31.trinuc$ctg)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.aca.ctt<-((dn31.trinuc$aca*dn31.trinuc$ctt)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl

cps.aca.gaa<-((dn31.trinuc$aca*dn31.trinuc$gaa)/(dn31.count.aa$t*dn31.count.aa$e))/dn31.count.diaa$te
cps.aca.gac<-((dn31.trinuc$aca*dn31.trinuc$gac)/(dn31.count.aa$t*dn31.count.aa$d))/dn31.count.diaa$td
cps.aca.gag<-((dn31.trinuc$aca*dn31.trinuc$gag)/(dn31.count.aa$t*dn31.count.aa$e))/dn31.count.diaa$te
cps.aca.gat<-((dn31.trinuc$aca*dn31.trinuc$gat)/(dn31.count.aa$t*dn31.count.aa$d))/dn31.count.diaa$td

cps.aca.gca<-((dn31.trinuc$aca*dn31.trinuc$gca)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.aca.gcc<-((dn31.trinuc$aca*dn31.trinuc$gcc)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.aca.gcg<-((dn31.trinuc$aca*dn31.trinuc$gcg)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.aca.gct<-((dn31.trinuc$aca*dn31.trinuc$gct)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta

cps.aca.gga<-((dn31.trinuc$aca*dn31.trinuc$gga)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.aca.ggc<-((dn31.trinuc$aca*dn31.trinuc$ggc)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.aca.ggg<-((dn31.trinuc$aca*dn31.trinuc$ggg)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.aca.ggt<-((dn31.trinuc$aca*dn31.trinuc$ggt)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg

cps.aca.gta<-((dn31.trinuc$aca*dn31.trinuc$gta)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.aca.gtc<-((dn31.trinuc$aca*dn31.trinuc$gtc)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.aca.gtg<-((dn31.trinuc$aca*dn31.trinuc$gtg)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.aca.gtt<-((dn31.trinuc$aca*dn31.trinuc$gtt)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv

#Stop codon
#cps.aca.taa<-((dn31.trinuc$aca*dn31.trinuc$taa)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aca.tac<-((dn31.trinuc$aca*dn31.trinuc$tac)/(dn31.count.aa$t*dn31.count.aa$y))/dn31.count.diaa$ty
#Stop codon
#cps.aca.tag<-((dn31.trinuc$aca*dn31.trinuc$tag)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aca.tat<-((dn31.trinuc$aca*dn31.trinuc$tat)/(dn31.count.aa$t*dn31.count.aa$y))/dn31.count.diaa$ty

cps.aca.tca<-((dn31.trinuc$aca*dn31.trinuc$tca)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.aca.tcc<-((dn31.trinuc$aca*dn31.trinuc$tcc)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.aca.tcg<-((dn31.trinuc$aca*dn31.trinuc$tcg)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.aca.tct<-((dn31.trinuc$aca*dn31.trinuc$tct)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts

#Stop codon
#cps.aca.tga<-((dn31.trinuc$aca*dn31.trinuc$tga)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aca.tgc<-((dn31.trinuc$aca*dn31.trinuc$tgc)/(dn31.count.aa$t*dn31.count.aa$c))/dn31.count.diaa$tc
cps.aca.tgg<-((dn31.trinuc$aca*dn31.trinuc$tgg)/(dn31.count.aa$t*dn31.count.aa$w))/dn31.count.diaa$tw
cps.aca.tgt<-((dn31.trinuc$aca*dn31.trinuc$tgt)/(dn31.count.aa$t*dn31.count.aa$c))/dn31.count.diaa$tc

cps.aca.tta<-((dn31.trinuc$aca*dn31.trinuc$tta)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.aca.ttc<-((dn31.trinuc$aca*dn31.trinuc$ttc)/(dn31.count.aa$t*dn31.count.aa$f))/dn31.count.diaa$tf
cps.aca.ttg<-((dn31.trinuc$aca*dn31.trinuc$ttg)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.aca.ttt<-((dn31.trinuc$aca*dn31.trinuc$ttt)/(dn31.count.aa$t*dn31.count.aa$f))/dn31.count.diaa$tf







cps.acc.aaa<-((dn31.trinuc$acc*dn31.trinuc$aaa)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$tk
cps.acc.aac<-((dn31.trinuc$acc*dn31.trinuc$aac)/(dn31.count.aa$t*dn31.count.aa$n))/dn31.count.diaa$tn
cps.acc.aag<-((dn31.trinuc$acc*dn31.trinuc$aag)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$tk
cps.acc.aat<-((dn31.trinuc$acc*dn31.trinuc$aat)/(dn31.count.aa$t*dn31.count.aa$n))/dn31.count.diaa$tn

cps.acc.aca<-((dn31.trinuc$acc*dn31.trinuc$aca)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.acc.acc<-((dn31.trinuc$acc*dn31.trinuc$acc)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.acc.acg<-((dn31.trinuc$acc*dn31.trinuc$acg)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.acc.act<-((dn31.trinuc$acc*dn31.trinuc$act)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt

cps.acc.aga<-((dn31.trinuc$acc*dn31.trinuc$aga)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acc.agc<-((dn31.trinuc$acc*dn31.trinuc$agc)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.acc.agg<-((dn31.trinuc$acc*dn31.trinuc$agg)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acc.agt<-((dn31.trinuc$acc*dn31.trinuc$agt)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts

cps.acc.ata<-((dn31.trinuc$acc*dn31.trinuc$ata)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti
cps.acc.atc<-((dn31.trinuc$acc*dn31.trinuc$atc)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti
cps.acc.atg<-((dn31.trinuc$acc*dn31.trinuc$atg)/(dn31.count.aa$t*dn31.count.aa$m))/dn31.count.diaa$tm
cps.acc.att<-((dn31.trinuc$acc*dn31.trinuc$att)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti

cps.acc.caa<-((dn31.trinuc$acc*dn31.trinuc$caa)/(dn31.count.aa$t*dn31.count.aa$q))/dn31.count.diaa$tq
cps.acc.cac<-((dn31.trinuc$acc*dn31.trinuc$cac)/(dn31.count.aa$t*dn31.count.aa$h))/dn31.count.diaa$th
cps.acc.cag<-((dn31.trinuc$acc*dn31.trinuc$cag)/(dn31.count.aa$t*dn31.count.aa$q))/dn31.count.diaa$tq
cps.acc.cat<-((dn31.trinuc$acc*dn31.trinuc$cat)/(dn31.count.aa$t*dn31.count.aa$h))/dn31.count.diaa$th

cps.acc.cca<-((dn31.trinuc$acc*dn31.trinuc$cca)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.acc.ccc<-((dn31.trinuc$acc*dn31.trinuc$ccc)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.acc.ccg<-((dn31.trinuc$acc*dn31.trinuc$ccg)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.acc.cct<-((dn31.trinuc$acc*dn31.trinuc$cct)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp

cps.acc.cga<-((dn31.trinuc$acc*dn31.trinuc$cga)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acc.cgc<-((dn31.trinuc$acc*dn31.trinuc$cgc)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acc.cgg<-((dn31.trinuc$acc*dn31.trinuc$cgg)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acc.cgt<-((dn31.trinuc$acc*dn31.trinuc$cgt)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr

cps.acc.cta<-((dn31.trinuc$acc*dn31.trinuc$cta)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acc.ctc<-((dn31.trinuc$acc*dn31.trinuc$ctc)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acc.ctg<-((dn31.trinuc$acc*dn31.trinuc$ctg)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acc.ctt<-((dn31.trinuc$acc*dn31.trinuc$ctt)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl

cps.acc.gaa<-((dn31.trinuc$acc*dn31.trinuc$gaa)/(dn31.count.aa$t*dn31.count.aa$e))/dn31.count.diaa$te
cps.acc.gac<-((dn31.trinuc$acc*dn31.trinuc$gac)/(dn31.count.aa$t*dn31.count.aa$d))/dn31.count.diaa$td
cps.acc.gag<-((dn31.trinuc$acc*dn31.trinuc$gag)/(dn31.count.aa$t*dn31.count.aa$e))/dn31.count.diaa$te
cps.acc.gat<-((dn31.trinuc$acc*dn31.trinuc$gat)/(dn31.count.aa$t*dn31.count.aa$d))/dn31.count.diaa$td

cps.acc.gca<-((dn31.trinuc$acc*dn31.trinuc$gca)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.acc.gcc<-((dn31.trinuc$acc*dn31.trinuc$gcc)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.acc.gcg<-((dn31.trinuc$acc*dn31.trinuc$gcg)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.acc.gct<-((dn31.trinuc$acc*dn31.trinuc$gct)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta

cps.acc.gga<-((dn31.trinuc$acc*dn31.trinuc$gga)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.acc.ggc<-((dn31.trinuc$acc*dn31.trinuc$ggc)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.acc.ggg<-((dn31.trinuc$acc*dn31.trinuc$ggg)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.acc.ggt<-((dn31.trinuc$acc*dn31.trinuc$ggt)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg

cps.acc.gta<-((dn31.trinuc$acc*dn31.trinuc$gta)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.acc.gtc<-((dn31.trinuc$acc*dn31.trinuc$gtc)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.acc.gtg<-((dn31.trinuc$acc*dn31.trinuc$gtg)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.acc.gtt<-((dn31.trinuc$acc*dn31.trinuc$gtt)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv

#Stop codon
#cps.acc.taa<-((dn31.trinuc$acc*dn31.trinuc$taa)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.acc.tac<-((dn31.trinuc$acc*dn31.trinuc$tac)/(dn31.count.aa$t*dn31.count.aa$y))/dn31.count.diaa$ty
#stop codon
#cps.acc.tag<-((dn31.trinuc$acc*dn31.trinuc$tag)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.acc.tat<-((dn31.trinuc$acc*dn31.trinuc$tat)/(dn31.count.aa$t*dn31.count.aa$y))/dn31.count.diaa$ty

cps.acc.tca<-((dn31.trinuc$acc*dn31.trinuc$tca)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.acc.tcc<-((dn31.trinuc$acc*dn31.trinuc$tcc)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.acc.tcg<-((dn31.trinuc$acc*dn31.trinuc$tcg)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.acc.tct<-((dn31.trinuc$acc*dn31.trinuc$tct)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts

#Stop codon
#cps.acc.tga<-((dn31.trinuc$acc*dn31.trinuc$tga)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.acc.tgc<-((dn31.trinuc$acc*dn31.trinuc$tgc)/(dn31.count.aa$t*dn31.count.aa$c))/dn31.count.diaa$tc
cps.acc.tgg<-((dn31.trinuc$acc*dn31.trinuc$tgg)/(dn31.count.aa$t*dn31.count.aa$w))/dn31.count.diaa$tw
cps.acc.tgt<-((dn31.trinuc$acc*dn31.trinuc$tgt)/(dn31.count.aa$t*dn31.count.aa$c))/dn31.count.diaa$tc

cps.acc.tta<-((dn31.trinuc$acc*dn31.trinuc$tta)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acc.ttc<-((dn31.trinuc$acc*dn31.trinuc$ttc)/(dn31.count.aa$t*dn31.count.aa$f))/dn31.count.diaa$tf
cps.acc.ttg<-((dn31.trinuc$acc*dn31.trinuc$ttg)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acc.ttt<-((dn31.trinuc$acc*dn31.trinuc$ttt)/(dn31.count.aa$t*dn31.count.aa$f))/dn31.count.diaa$tf







cps.acg.aaa<-((dn31.trinuc$acg*dn31.trinuc$aaa)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$tk
cps.acg.aac<-((dn31.trinuc$acg*dn31.trinuc$aac)/(dn31.count.aa$t*dn31.count.aa$n))/dn31.count.diaa$tn
cps.acg.aag<-((dn31.trinuc$acg*dn31.trinuc$aag)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$tk
cps.acg.aat<-((dn31.trinuc$acg*dn31.trinuc$aat)/(dn31.count.aa$t*dn31.count.aa$n))/dn31.count.diaa$tn

cps.acg.aca<-((dn31.trinuc$acg*dn31.trinuc$aca)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.acg.acc<-((dn31.trinuc$acg*dn31.trinuc$acc)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.acg.acg<-((dn31.trinuc$acg*dn31.trinuc$acg)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.acg.act<-((dn31.trinuc$acg*dn31.trinuc$act)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt

cps.acg.aga<-((dn31.trinuc$acg*dn31.trinuc$aga)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acg.agc<-((dn31.trinuc$acg*dn31.trinuc$agc)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.acg.agg<-((dn31.trinuc$acg*dn31.trinuc$agg)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acg.agt<-((dn31.trinuc$acg*dn31.trinuc$agt)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts

cps.acg.ata<-((dn31.trinuc$acg*dn31.trinuc$ata)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti
cps.acg.atc<-((dn31.trinuc$acg*dn31.trinuc$atc)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti
cps.acg.atg<-((dn31.trinuc$acg*dn31.trinuc$atg)/(dn31.count.aa$t*dn31.count.aa$m))/dn31.count.diaa$tm
cps.acg.att<-((dn31.trinuc$acg*dn31.trinuc$att)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti

cps.acg.caa<-((dn31.trinuc$acg*dn31.trinuc$caa)/(dn31.count.aa$t*dn31.count.aa$q))/dn31.count.diaa$tq
cps.acg.cac<-((dn31.trinuc$acg*dn31.trinuc$cac)/(dn31.count.aa$t*dn31.count.aa$h))/dn31.count.diaa$th
cps.acg.cag<-((dn31.trinuc$acg*dn31.trinuc$cag)/(dn31.count.aa$t*dn31.count.aa$q))/dn31.count.diaa$tq
cps.acg.cat<-((dn31.trinuc$acg*dn31.trinuc$cat)/(dn31.count.aa$t*dn31.count.aa$h))/dn31.count.diaa$th

cps.acg.cca<-((dn31.trinuc$acg*dn31.trinuc$cca)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.acg.ccc<-((dn31.trinuc$acg*dn31.trinuc$ccc)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.acg.ccg<-((dn31.trinuc$acg*dn31.trinuc$ccg)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.acg.cct<-((dn31.trinuc$acg*dn31.trinuc$cct)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp

cps.acg.cga<-((dn31.trinuc$acg*dn31.trinuc$cga)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acg.cgc<-((dn31.trinuc$acg*dn31.trinuc$cgc)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acg.cgg<-((dn31.trinuc$acg*dn31.trinuc$cgg)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.acg.cgt<-((dn31.trinuc$acg*dn31.trinuc$cgt)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr

cps.acg.cta<-((dn31.trinuc$acg*dn31.trinuc$cta)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acg.ctc<-((dn31.trinuc$acg*dn31.trinuc$ctc)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acg.ctg<-((dn31.trinuc$acg*dn31.trinuc$ctg)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acg.ctt<-((dn31.trinuc$acg*dn31.trinuc$ctt)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl

cps.acg.gaa<-((dn31.trinuc$acg*dn31.trinuc$gaa)/(dn31.count.aa$t*dn31.count.aa$e))/dn31.count.diaa$te
cps.acg.gac<-((dn31.trinuc$acg*dn31.trinuc$gac)/(dn31.count.aa$t*dn31.count.aa$d))/dn31.count.diaa$td
cps.acg.gag<-((dn31.trinuc$acg*dn31.trinuc$gag)/(dn31.count.aa$t*dn31.count.aa$e))/dn31.count.diaa$te
cps.acg.gat<-((dn31.trinuc$acg*dn31.trinuc$gat)/(dn31.count.aa$t*dn31.count.aa$d))/dn31.count.diaa$td

cps.acg.gca<-((dn31.trinuc$acg*dn31.trinuc$gca)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.acg.gcc<-((dn31.trinuc$acg*dn31.trinuc$gcc)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.acg.gcg<-((dn31.trinuc$acg*dn31.trinuc$gcg)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.acg.gct<-((dn31.trinuc$acg*dn31.trinuc$gct)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta

cps.acg.gga<-((dn31.trinuc$acg*dn31.trinuc$gga)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.acg.ggc<-((dn31.trinuc$acg*dn31.trinuc$ggc)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.acg.ggg<-((dn31.trinuc$acg*dn31.trinuc$ggg)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.acg.ggt<-((dn31.trinuc$acg*dn31.trinuc$ggt)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg

cps.acg.gta<-((dn31.trinuc$acg*dn31.trinuc$gta)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.acg.gtc<-((dn31.trinuc$acg*dn31.trinuc$gtc)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.acg.gtg<-((dn31.trinuc$acg*dn31.trinuc$gtg)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.acg.gtt<-((dn31.trinuc$acg*dn31.trinuc$gtt)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv

#Stop codon
#cps.acg.taa<-((dn31.trinuc$acg*dn31.trinuc$taa)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.acg.tac<-((dn31.trinuc$acg*dn31.trinuc$tac)/(dn31.count.aa$t*dn31.count.aa$y))/dn31.count.diaa$ty
#Stop codon
#cps.acg.tag<-((dn31.trinuc$acg*dn31.trinuc$tag)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.acg.tat<-((dn31.trinuc$acg*dn31.trinuc$tat)/(dn31.count.aa$t*dn31.count.aa$y))/dn31.count.diaa$ty

cps.acg.tca<-((dn31.trinuc$acg*dn31.trinuc$tca)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.acg.tcc<-((dn31.trinuc$acg*dn31.trinuc$tcc)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.acg.tcg<-((dn31.trinuc$acg*dn31.trinuc$tcg)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.acg.tct<-((dn31.trinuc$acg*dn31.trinuc$tct)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts

#Stop codon
#cps.acg.tga<-((dn31.trinuc$acg*dn31.trinuc$tga)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.acg.tgc<-((dn31.trinuc$acg*dn31.trinuc$tgc)/(dn31.count.aa$t*dn31.count.aa$c))/dn31.count.diaa$tc
cps.acg.tgg<-((dn31.trinuc$acg*dn31.trinuc$tgg)/(dn31.count.aa$t*dn31.count.aa$w))/dn31.count.diaa$tw
cps.acg.tgt<-((dn31.trinuc$acg*dn31.trinuc$tgt)/(dn31.count.aa$t*dn31.count.aa$c))/dn31.count.diaa$tc

cps.acg.tta<-((dn31.trinuc$acg*dn31.trinuc$tta)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acg.ttc<-((dn31.trinuc$acg*dn31.trinuc$ttc)/(dn31.count.aa$t*dn31.count.aa$f))/dn31.count.diaa$tf
cps.acg.ttg<-((dn31.trinuc$acg*dn31.trinuc$ttg)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.acg.ttt<-((dn31.trinuc$acg*dn31.trinuc$ttt)/(dn31.count.aa$t*dn31.count.aa$f))/dn31.count.diaa$tf







cps.act.aaa<-((dn31.trinuc$act*dn31.trinuc$aaa)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$tk
cps.act.aac<-((dn31.trinuc$act*dn31.trinuc$aac)/(dn31.count.aa$t*dn31.count.aa$n))/dn31.count.diaa$tn
cps.act.aag<-((dn31.trinuc$act*dn31.trinuc$aag)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$tk
cps.act.aat<-((dn31.trinuc$act*dn31.trinuc$aat)/(dn31.count.aa$t*dn31.count.aa$n))/dn31.count.diaa$tn

cps.act.aca<-((dn31.trinuc$act*dn31.trinuc$aca)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.act.acc<-((dn31.trinuc$act*dn31.trinuc$acc)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.act.acg<-((dn31.trinuc$act*dn31.trinuc$acg)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt
cps.act.act<-((dn31.trinuc$act*dn31.trinuc$act)/(dn31.count.aa$t*dn31.count.aa$t))/dn31.count.diaa$tt

cps.act.aga<-((dn31.trinuc$act*dn31.trinuc$aga)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.act.agc<-((dn31.trinuc$act*dn31.trinuc$agc)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.act.agg<-((dn31.trinuc$act*dn31.trinuc$agg)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.act.agt<-((dn31.trinuc$act*dn31.trinuc$agt)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts

cps.act.ata<-((dn31.trinuc$act*dn31.trinuc$ata)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti
cps.act.atc<-((dn31.trinuc$act*dn31.trinuc$atc)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti
cps.act.atg<-((dn31.trinuc$act*dn31.trinuc$atg)/(dn31.count.aa$t*dn31.count.aa$m))/dn31.count.diaa$tw
cps.act.att<-((dn31.trinuc$act*dn31.trinuc$att)/(dn31.count.aa$t*dn31.count.aa$i))/dn31.count.diaa$ti

cps.act.caa<-((dn31.trinuc$act*dn31.trinuc$caa)/(dn31.count.aa$t*dn31.count.aa$q))/dn31.count.diaa$tq
cps.act.cac<-((dn31.trinuc$act*dn31.trinuc$cac)/(dn31.count.aa$t*dn31.count.aa$h))/dn31.count.diaa$th
cps.act.cag<-((dn31.trinuc$act*dn31.trinuc$cag)/(dn31.count.aa$t*dn31.count.aa$q))/dn31.count.diaa$tq
cps.act.cat<-((dn31.trinuc$act*dn31.trinuc$cat)/(dn31.count.aa$t*dn31.count.aa$h))/dn31.count.diaa$th

cps.act.cca<-((dn31.trinuc$act*dn31.trinuc$cca)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.act.ccc<-((dn31.trinuc$act*dn31.trinuc$ccc)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.act.ccg<-((dn31.trinuc$act*dn31.trinuc$ccg)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp
cps.act.cct<-((dn31.trinuc$act*dn31.trinuc$cct)/(dn31.count.aa$t*dn31.count.aa$p))/dn31.count.diaa$tp

cps.act.cga<-((dn31.trinuc$act*dn31.trinuc$cga)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.act.cgc<-((dn31.trinuc$act*dn31.trinuc$cgc)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.act.cgg<-((dn31.trinuc$act*dn31.trinuc$cgg)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr
cps.act.cgt<-((dn31.trinuc$act*dn31.trinuc$cgt)/(dn31.count.aa$t*dn31.count.aa$r))/dn31.count.diaa$tr

cps.act.cta<-((dn31.trinuc$act*dn31.trinuc$cta)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.act.ctc<-((dn31.trinuc$act*dn31.trinuc$ctc)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.act.ctg<-((dn31.trinuc$act*dn31.trinuc$ctg)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl
cps.act.ctt<-((dn31.trinuc$act*dn31.trinuc$ctt)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$tl

cps.act.gaa<-((dn31.trinuc$act*dn31.trinuc$gaa)/(dn31.count.aa$t*dn31.count.aa$e))/dn31.count.diaa$te
cps.act.gac<-((dn31.trinuc$act*dn31.trinuc$gac)/(dn31.count.aa$t*dn31.count.aa$d))/dn31.count.diaa$td
cps.act.gag<-((dn31.trinuc$act*dn31.trinuc$gag)/(dn31.count.aa$t*dn31.count.aa$e))/dn31.count.diaa$te
cps.act.gat<-((dn31.trinuc$act*dn31.trinuc$gat)/(dn31.count.aa$t*dn31.count.aa$d))/dn31.count.diaa$td

cps.act.gca<-((dn31.trinuc$act*dn31.trinuc$gca)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.act.gcc<-((dn31.trinuc$act*dn31.trinuc$gcc)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.act.gcg<-((dn31.trinuc$act*dn31.trinuc$gcg)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta
cps.act.gct<-((dn31.trinuc$act*dn31.trinuc$gct)/(dn31.count.aa$t*dn31.count.aa$a))/dn31.count.diaa$ta

cps.act.gga<-((dn31.trinuc$act*dn31.trinuc$gga)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.act.ggc<-((dn31.trinuc$act*dn31.trinuc$ggc)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.act.ggg<-((dn31.trinuc$act*dn31.trinuc$ggg)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg
cps.act.ggt<-((dn31.trinuc$act*dn31.trinuc$ggt)/(dn31.count.aa$t*dn31.count.aa$g))/dn31.count.diaa$tg

cps.act.gta<-((dn31.trinuc$act*dn31.trinuc$gta)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.act.gtc<-((dn31.trinuc$act*dn31.trinuc$gtc)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.act.gtg<-((dn31.trinuc$act*dn31.trinuc$gtg)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv
cps.act.gtt<-((dn31.trinuc$act*dn31.trinuc$gtt)/(dn31.count.aa$t*dn31.count.aa$v))/dn31.count.diaa$tv

#Stop codon
#cps.act.taa<-((dn31.trinuc$act*dn31.trinuc$taa)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.act.tac<-((dn31.trinuc$act*dn31.trinuc$tac)/(dn31.count.aa$t*dn31.count.aa$y))/dn31.count.diaa$ty
#Stop codon
#cps.act.tag<-((dn31.trinuc$act*dn31.trinuc$tag)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.act.tat<-((dn31.trinuc$act*dn31.trinuc$tat)/(dn31.count.aa$t*dn31.count.aa$y))/dn31.count.diaa$ty

cps.act.tca<-((dn31.trinuc$act*dn31.trinuc$tca)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.act.tcc<-((dn31.trinuc$act*dn31.trinuc$tcc)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.act.tcg<-((dn31.trinuc$act*dn31.trinuc$tcg)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts
cps.act.tct<-((dn31.trinuc$act*dn31.trinuc$tct)/(dn31.count.aa$t*dn31.count.aa$s))/dn31.count.diaa$ts

#Stop codon
#cps.act.tga<-((dn31.trinuc$act*dn31.trinuc$tga)/(dn31.count.aa$t*dn31.count.aa$k))/dn31.count.diaa$kk
cps.act.tgc<-((dn31.trinuc$act*dn31.trinuc$tgc)/(dn31.count.aa$t*dn31.count.aa$c))/dn31.count.diaa$ts
cps.act.tgg<-((dn31.trinuc$act*dn31.trinuc$tgg)/(dn31.count.aa$t*dn31.count.aa$w))/dn31.count.diaa$ts
cps.act.tgt<-((dn31.trinuc$act*dn31.trinuc$tgt)/(dn31.count.aa$t*dn31.count.aa$c))/dn31.count.diaa$ts

cps.act.tta<-((dn31.trinuc$act*dn31.trinuc$tta)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$ts
cps.act.ttc<-((dn31.trinuc$act*dn31.trinuc$ttc)/(dn31.count.aa$t*dn31.count.aa$f))/dn31.count.diaa$ts
cps.act.ttg<-((dn31.trinuc$act*dn31.trinuc$ttg)/(dn31.count.aa$t*dn31.count.aa$l))/dn31.count.diaa$ts
cps.act.ttt<-((dn31.trinuc$act*dn31.trinuc$ttt)/(dn31.count.aa$t*dn31.count.aa$f))/dn31.count.diaa$ts













cps.aga.aaa<-((dn31.trinuc$aga*dn31.trinuc$aaa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.aga.aac<-((dn31.trinuc$aga*dn31.trinuc$aac)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn
cps.aga.aag<-((dn31.trinuc$aga*dn31.trinuc$aag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.aga.aat<-((dn31.trinuc$aga*dn31.trinuc$aat)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn

cps.aga.aca<-((dn31.trinuc$aga*dn31.trinuc$aca)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.aga.acc<-((dn31.trinuc$aga*dn31.trinuc$acc)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.aga.acg<-((dn31.trinuc$aga*dn31.trinuc$acg)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.aga.act<-((dn31.trinuc$aga*dn31.trinuc$act)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt

cps.aga.aga<-((dn31.trinuc$aga*dn31.trinuc$aga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.aga.agc<-((dn31.trinuc$aga*dn31.trinuc$agc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.aga.agg<-((dn31.trinuc$aga*dn31.trinuc$agg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.aga.agt<-((dn31.trinuc$aga*dn31.trinuc$agt)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

cps.aga.ata<-((dn31.trinuc$aga*dn31.trinuc$ata)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.aga.atc<-((dn31.trinuc$aga*dn31.trinuc$atc)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.aga.atg<-((dn31.trinuc$aga*dn31.trinuc$atg)/(dn31.count.aa$r*dn31.count.aa$m))/dn31.count.diaa$rm
cps.aga.att<-((dn31.trinuc$aga*dn31.trinuc$att)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri

cps.aga.caa<-((dn31.trinuc$aga*dn31.trinuc$caa)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.aga.cac<-((dn31.trinuc$aga*dn31.trinuc$cac)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh
cps.aga.cag<-((dn31.trinuc$aga*dn31.trinuc$cag)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.aga.cat<-((dn31.trinuc$aga*dn31.trinuc$cat)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh

cps.aga.cca<-((dn31.trinuc$aga*dn31.trinuc$cca)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.aga.ccc<-((dn31.trinuc$aga*dn31.trinuc$ccc)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.aga.ccg<-((dn31.trinuc$aga*dn31.trinuc$ccg)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.aga.cct<-((dn31.trinuc$aga*dn31.trinuc$cct)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp

cps.aga.cga<-((dn31.trinuc$aga*dn31.trinuc$cga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.aga.cgc<-((dn31.trinuc$aga*dn31.trinuc$cgc)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.aga.cgg<-((dn31.trinuc$aga*dn31.trinuc$cgg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.aga.cgt<-((dn31.trinuc$aga*dn31.trinuc$cgt)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr

cps.aga.cta<-((dn31.trinuc$aga*dn31.trinuc$cta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.aga.ctc<-((dn31.trinuc$aga*dn31.trinuc$ctc)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.aga.ctg<-((dn31.trinuc$aga*dn31.trinuc$ctg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.aga.ctt<-((dn31.trinuc$aga*dn31.trinuc$ctt)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl

cps.aga.gaa<-((dn31.trinuc$aga*dn31.trinuc$gaa)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.aga.gac<-((dn31.trinuc$aga*dn31.trinuc$gac)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd
cps.aga.gag<-((dn31.trinuc$aga*dn31.trinuc$gag)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.aga.gat<-((dn31.trinuc$aga*dn31.trinuc$gat)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd

cps.aga.gca<-((dn31.trinuc$aga*dn31.trinuc$gca)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.aga.gcc<-((dn31.trinuc$aga*dn31.trinuc$gcc)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.aga.gcg<-((dn31.trinuc$aga*dn31.trinuc$gcg)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.aga.gct<-((dn31.trinuc$aga*dn31.trinuc$gct)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra

cps.aga.gga<-((dn31.trinuc$aga*dn31.trinuc$gga)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.aga.ggc<-((dn31.trinuc$aga*dn31.trinuc$ggc)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.aga.ggg<-((dn31.trinuc$aga*dn31.trinuc$ggg)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.aga.ggt<-((dn31.trinuc$aga*dn31.trinuc$ggt)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg

cps.aga.gta<-((dn31.trinuc$aga*dn31.trinuc$gta)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.aga.gtc<-((dn31.trinuc$aga*dn31.trinuc$gtc)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.aga.gtg<-((dn31.trinuc$aga*dn31.trinuc$gtg)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.aga.gtt<-((dn31.trinuc$aga*dn31.trinuc$gtt)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv

#Stop codon
#cps.aga.taa<-((dn31.trinuc$aga*dn31.trinuc$taa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aga.tac<-((dn31.trinuc$aga*dn31.trinuc$tac)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry
#Stop codon
#cps.aga.tag<-((dn31.trinuc$aga*dn31.trinuc$tag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aga.tat<-((dn31.trinuc$aga*dn31.trinuc$tat)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry

cps.aga.tca<-((dn31.trinuc$aga*dn31.trinuc$tca)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.aga.tcc<-((dn31.trinuc$aga*dn31.trinuc$tcc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.aga.tcg<-((dn31.trinuc$aga*dn31.trinuc$tcg)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.aga.tct<-((dn31.trinuc$aga*dn31.trinuc$tct)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

#Stop codon
#cps.aga.tga<-((dn31.trinuc$aga*dn31.trinuc$tga)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.aga.tgc<-((dn31.trinuc$aga*dn31.trinuc$tgc)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc
cps.aga.tgg<-((dn31.trinuc$aga*dn31.trinuc$tgg)/(dn31.count.aa$r*dn31.count.aa$w))/dn31.count.diaa$rw
cps.aga.tgt<-((dn31.trinuc$aga*dn31.trinuc$tgt)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc

cps.aga.tta<-((dn31.trinuc$aga*dn31.trinuc$tta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.aga.ttc<-((dn31.trinuc$aga*dn31.trinuc$ttc)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf
cps.aga.ttg<-((dn31.trinuc$aga*dn31.trinuc$ttg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.aga.ttt<-((dn31.trinuc$aga*dn31.trinuc$ttt)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf







cps.agc.aaa<-((dn31.trinuc$agc*dn31.trinuc$aaa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.agc.aac<-((dn31.trinuc$acg*dn31.trinuc$aac)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn
cps.agc.aag<-((dn31.trinuc$agc*dn31.trinuc$aag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.agc.aat<-((dn31.trinuc$agc*dn31.trinuc$aat)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn

cps.agc.aca<-((dn31.trinuc$agc*dn31.trinuc$aca)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.agc.acc<-((dn31.trinuc$agc*dn31.trinuc$acc)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.agc.acg<-((dn31.trinuc$agc*dn31.trinuc$acg)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.agc.act<-((dn31.trinuc$agc*dn31.trinuc$act)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st

cps.agc.aga<-((dn31.trinuc$agc*dn31.trinuc$aga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agc.agc<-((dn31.trinuc$agc*dn31.trinuc$agc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.agc.agg<-((dn31.trinuc$agc*dn31.trinuc$agg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agc.agt<-((dn31.trinuc$agc*dn31.trinuc$agt)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

cps.agc.ata<-((dn31.trinuc$agc*dn31.trinuc$ata)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.agc.atc<-((dn31.trinuc$agc*dn31.trinuc$atc)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.agc.atg<-((dn31.trinuc$agc*dn31.trinuc$atg)/(dn31.count.aa$s*dn31.count.aa$m))/dn31.count.diaa$sm
cps.agc.att<-((dn31.trinuc$agc*dn31.trinuc$att)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si

cps.agc.caa<-((dn31.trinuc$agc*dn31.trinuc$caa)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.agc.cac<-((dn31.trinuc$agc*dn31.trinuc$cac)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh
cps.agc.cag<-((dn31.trinuc$agc*dn31.trinuc$cag)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.agc.cat<-((dn31.trinuc$agc*dn31.trinuc$cat)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh

cps.agc.cca<-((dn31.trinuc$agc*dn31.trinuc$cca)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.agc.ccc<-((dn31.trinuc$agc*dn31.trinuc$ccc)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.agc.ccg<-((dn31.trinuc$agc*dn31.trinuc$ccg)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.agc.cct<-((dn31.trinuc$agc*dn31.trinuc$cct)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp

cps.agc.cga<-((dn31.trinuc$agc*dn31.trinuc$cga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agc.cgc<-((dn31.trinuc$agc*dn31.trinuc$cgc)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agc.cgg<-((dn31.trinuc$agc*dn31.trinuc$cgg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agc.cgt<-((dn31.trinuc$agc*dn31.trinuc$cgt)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr

cps.agc.cta<-((dn31.trinuc$agc*dn31.trinuc$cta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agc.ctc<-((dn31.trinuc$agc*dn31.trinuc$ctc)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agc.ctg<-((dn31.trinuc$agc*dn31.trinuc$ctg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agc.ctt<-((dn31.trinuc$agc*dn31.trinuc$ctt)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl

cps.agc.gaa<-((dn31.trinuc$agc*dn31.trinuc$gaa)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.agc.gac<-((dn31.trinuc$agc*dn31.trinuc$gac)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd
cps.agc.gag<-((dn31.trinuc$agc*dn31.trinuc$gag)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.agc.gat<-((dn31.trinuc$agc*dn31.trinuc$gat)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd

cps.agc.gca<-((dn31.trinuc$agc*dn31.trinuc$gca)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.agc.gcc<-((dn31.trinuc$agc*dn31.trinuc$gcc)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.agc.gcg<-((dn31.trinuc$agc*dn31.trinuc$gcg)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.agc.gct<-((dn31.trinuc$agc*dn31.trinuc$gct)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa

cps.agc.gga<-((dn31.trinuc$agc*dn31.trinuc$gga)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.agc.ggc<-((dn31.trinuc$agc*dn31.trinuc$ggc)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.agc.ggg<-((dn31.trinuc$agc*dn31.trinuc$ggg)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.agc.ggt<-((dn31.trinuc$agc*dn31.trinuc$ggt)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg

cps.agc.gta<-((dn31.trinuc$agc*dn31.trinuc$gta)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.agc.gtc<-((dn31.trinuc$agc*dn31.trinuc$gtc)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.agc.gtg<-((dn31.trinuc$agc*dn31.trinuc$gtg)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.agc.gtt<-((dn31.trinuc$agc*dn31.trinuc$gtt)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv

#Stop codon
#cps.agc.taa<-((dn31.trinuc$agc*dn31.trinuc$taa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agc.tac<-((dn31.trinuc$agc*dn31.trinuc$tac)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy
#Stop codon
#cps.agc.tag<-((dn31.trinuc$agc*dn31.trinuc$tag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agc.tat<-((dn31.trinuc$agc*dn31.trinuc$tat)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy

cps.agc.tca<-((dn31.trinuc$agc*dn31.trinuc$tca)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.agc.tcc<-((dn31.trinuc$agc*dn31.trinuc$tcc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.agc.tcg<-((dn31.trinuc$agc*dn31.trinuc$tcg)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.agc.tct<-((dn31.trinuc$agc*dn31.trinuc$tct)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

#Stop codon
#cps.agc.tga<-((dn31.trinuc$agc*dn31.trinuc$tga)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agc.tgc<-((dn31.trinuc$agc*dn31.trinuc$tgc)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc
cps.agc.tgg<-((dn31.trinuc$agc*dn31.trinuc$tgg)/(dn31.count.aa$s*dn31.count.aa$w))/dn31.count.diaa$sw
cps.agc.tgt<-((dn31.trinuc$acg*dn31.trinuc$tgt)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc

cps.agc.tta<-((dn31.trinuc$agc*dn31.trinuc$tta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agc.ttc<-((dn31.trinuc$agc*dn31.trinuc$ttc)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf
cps.agc.ttg<-((dn31.trinuc$agc*dn31.trinuc$ttg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agc.ttt<-((dn31.trinuc$agc*dn31.trinuc$ttt)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf







cps.agg.aaa<-((dn31.trinuc$agg*dn31.trinuc$aaa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.agg.aac<-((dn31.trinuc$agg*dn31.trinuc$aac)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn
cps.agg.aag<-((dn31.trinuc$agg*dn31.trinuc$aag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.agg.aat<-((dn31.trinuc$agg*dn31.trinuc$aat)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn

cps.agg.aca<-((dn31.trinuc$agg*dn31.trinuc$aca)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.agg.acc<-((dn31.trinuc$agg*dn31.trinuc$acc)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.agg.acg<-((dn31.trinuc$agg*dn31.trinuc$acg)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.agg.act<-((dn31.trinuc$agg*dn31.trinuc$act)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt

cps.agg.aga<-((dn31.trinuc$agg*dn31.trinuc$aga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.agg.agc<-((dn31.trinuc$agg*dn31.trinuc$agc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.agg.agg<-((dn31.trinuc$agg*dn31.trinuc$agg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.agg.agt<-((dn31.trinuc$agg*dn31.trinuc$agt)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

cps.agg.ata<-((dn31.trinuc$agg*dn31.trinuc$ata)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.agg.atc<-((dn31.trinuc$agg*dn31.trinuc$atc)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.agg.atg<-((dn31.trinuc$agg*dn31.trinuc$atg)/(dn31.count.aa$r*dn31.count.aa$m))/dn31.count.diaa$rm
cps.agg.att<-((dn31.trinuc$agg*dn31.trinuc$att)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri

cps.agg.caa<-((dn31.trinuc$agg*dn31.trinuc$caa)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.agg.cac<-((dn31.trinuc$agg*dn31.trinuc$cac)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh
cps.agg.cag<-((dn31.trinuc$agg*dn31.trinuc$cag)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.agg.cat<-((dn31.trinuc$agg*dn31.trinuc$cat)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh

cps.agg.cca<-((dn31.trinuc$agg*dn31.trinuc$cca)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.agg.ccc<-((dn31.trinuc$agg*dn31.trinuc$ccc)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.agg.ccg<-((dn31.trinuc$agg*dn31.trinuc$ccg)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.agg.cct<-((dn31.trinuc$agg*dn31.trinuc$cct)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp

cps.agg.cga<-((dn31.trinuc$agg*dn31.trinuc$cga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.agg.cgc<-((dn31.trinuc$agg*dn31.trinuc$cgc)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.agg.cgg<-((dn31.trinuc$agg*dn31.trinuc$cgg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.agg.cgt<-((dn31.trinuc$agg*dn31.trinuc$cgt)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr

cps.agg.cta<-((dn31.trinuc$agg*dn31.trinuc$cta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.agg.ctc<-((dn31.trinuc$agg*dn31.trinuc$ctc)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.agg.ctg<-((dn31.trinuc$agg*dn31.trinuc$ctg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.agg.ctt<-((dn31.trinuc$agg*dn31.trinuc$ctt)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl

cps.agg.gaa<-((dn31.trinuc$agg*dn31.trinuc$gaa)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.agg.gac<-((dn31.trinuc$agg*dn31.trinuc$gac)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd
cps.agg.gag<-((dn31.trinuc$agg*dn31.trinuc$gag)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.agg.gat<-((dn31.trinuc$agg*dn31.trinuc$gat)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd

cps.agg.gca<-((dn31.trinuc$agg*dn31.trinuc$gca)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.agg.gcc<-((dn31.trinuc$agg*dn31.trinuc$gcc)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.agg.gcg<-((dn31.trinuc$agg*dn31.trinuc$gcg)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.agg.gct<-((dn31.trinuc$agg*dn31.trinuc$gct)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra

cps.agg.gga<-((dn31.trinuc$agg*dn31.trinuc$gga)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.agg.ggc<-((dn31.trinuc$agg*dn31.trinuc$ggc)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.agg.ggg<-((dn31.trinuc$agg*dn31.trinuc$ggg)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.agg.ggt<-((dn31.trinuc$agg*dn31.trinuc$ggt)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg

cps.agg.gta<-((dn31.trinuc$agg*dn31.trinuc$gta)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.agg.gtc<-((dn31.trinuc$agg*dn31.trinuc$gtc)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.agg.gtg<-((dn31.trinuc$agg*dn31.trinuc$gtg)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.agg.gtt<-((dn31.trinuc$agg*dn31.trinuc$gtt)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv

#Stop codon
#cps.agg.taa<-((dn31.trinuc$agg*dn31.trinuc$taa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agg.tac<-((dn31.trinuc$agg*dn31.trinuc$tac)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry
#Stop codon
#cps.agg.tag<-((dn31.trinuc$agg*dn31.trinuc$tag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agg.tat<-((dn31.trinuc$agg*dn31.trinuc$tat)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry

cps.agg.tca<-((dn31.trinuc$agg*dn31.trinuc$tca)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.agg.tcc<-((dn31.trinuc$agg*dn31.trinuc$tcc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.agg.tcg<-((dn31.trinuc$agg*dn31.trinuc$tcg)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.agg.tct<-((dn31.trinuc$agg*dn31.trinuc$tct)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

#Stop codon
#cps.agg.tga<-((dn31.trinuc$agg*dn31.trinuc$tga)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agg.tgc<-((dn31.trinuc$agg*dn31.trinuc$tgc)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc
cps.agg.tgg<-((dn31.trinuc$agg*dn31.trinuc$tgg)/(dn31.count.aa$r*dn31.count.aa$w))/dn31.count.diaa$rw
cps.agg.tgt<-((dn31.trinuc$agg*dn31.trinuc$tgt)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc

cps.agg.tta<-((dn31.trinuc$agg*dn31.trinuc$tta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.agg.ttc<-((dn31.trinuc$agg*dn31.trinuc$ttc)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf
cps.agg.ttg<-((dn31.trinuc$agg*dn31.trinuc$ttg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.agg.ttt<-((dn31.trinuc$agg*dn31.trinuc$ttt)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf







cps.agt.aaa<-((dn31.trinuc$agt*dn31.trinuc$aaa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.agt.aac<-((dn31.trinuc$agt*dn31.trinuc$aac)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn
cps.agt.aag<-((dn31.trinuc$agt*dn31.trinuc$aag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.agt.aat<-((dn31.trinuc$agt*dn31.trinuc$aat)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn

cps.agt.aca<-((dn31.trinuc$agt*dn31.trinuc$aca)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.agt.acc<-((dn31.trinuc$agt*dn31.trinuc$acc)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.agt.acg<-((dn31.trinuc$agt*dn31.trinuc$acg)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.agt.act<-((dn31.trinuc$agt*dn31.trinuc$act)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st

cps.agt.aga<-((dn31.trinuc$agt*dn31.trinuc$aga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agt.agc<-((dn31.trinuc$agt*dn31.trinuc$agc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.agt.agg<-((dn31.trinuc$agt*dn31.trinuc$agg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agt.agt<-((dn31.trinuc$agt*dn31.trinuc$agt)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

cps.agt.ata<-((dn31.trinuc$agt*dn31.trinuc$ata)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.agt.atc<-((dn31.trinuc$agt*dn31.trinuc$atc)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.agt.atg<-((dn31.trinuc$agt*dn31.trinuc$atg)/(dn31.count.aa$s*dn31.count.aa$m))/dn31.count.diaa$sm
cps.agt.att<-((dn31.trinuc$agt*dn31.trinuc$att)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si

cps.agt.caa<-((dn31.trinuc$agt*dn31.trinuc$caa)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.agt.cac<-((dn31.trinuc$agt*dn31.trinuc$cac)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh
cps.agt.cag<-((dn31.trinuc$agt*dn31.trinuc$cag)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.agt.cat<-((dn31.trinuc$agt*dn31.trinuc$cat)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh

cps.agt.cca<-((dn31.trinuc$agt*dn31.trinuc$cca)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.agt.ccc<-((dn31.trinuc$agt*dn31.trinuc$ccc)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.agt.ccg<-((dn31.trinuc$agt*dn31.trinuc$ccg)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.agt.cct<-((dn31.trinuc$agt*dn31.trinuc$cct)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp

cps.agt.cga<-((dn31.trinuc$agt*dn31.trinuc$cga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agt.cgc<-((dn31.trinuc$agt*dn31.trinuc$cgc)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agt.cgg<-((dn31.trinuc$agt*dn31.trinuc$cgg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.agt.cgt<-((dn31.trinuc$agt*dn31.trinuc$cgt)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr

cps.agt.cta<-((dn31.trinuc$agt*dn31.trinuc$cta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agt.ctc<-((dn31.trinuc$agt*dn31.trinuc$ctc)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agt.ctg<-((dn31.trinuc$agt*dn31.trinuc$ctg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agt.ctt<-((dn31.trinuc$agt*dn31.trinuc$ctt)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl

cps.agt.gaa<-((dn31.trinuc$agt*dn31.trinuc$gaa)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.agt.gac<-((dn31.trinuc$agt*dn31.trinuc$gac)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd
cps.agt.gag<-((dn31.trinuc$agt*dn31.trinuc$gag)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.agt.gat<-((dn31.trinuc$agt*dn31.trinuc$gat)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd

cps.agt.gca<-((dn31.trinuc$agt*dn31.trinuc$gca)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.agt.gcc<-((dn31.trinuc$agt*dn31.trinuc$gcc)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.agt.gcg<-((dn31.trinuc$agt*dn31.trinuc$gcg)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.agt.gct<-((dn31.trinuc$agt*dn31.trinuc$gct)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa

cps.agt.gga<-((dn31.trinuc$agt*dn31.trinuc$gga)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.agt.ggc<-((dn31.trinuc$agt*dn31.trinuc$ggc)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.agt.ggg<-((dn31.trinuc$agt*dn31.trinuc$ggg)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.agt.ggt<-((dn31.trinuc$agt*dn31.trinuc$ggt)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg

cps.agt.gta<-((dn31.trinuc$agt*dn31.trinuc$gta)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.agt.gtc<-((dn31.trinuc$agt*dn31.trinuc$gtc)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.agt.gtg<-((dn31.trinuc$agt*dn31.trinuc$gtg)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.agt.gtt<-((dn31.trinuc$agt*dn31.trinuc$gtt)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv

#Stop codon
#cps.agt.taa<-((dn31.trinuc$agt*dn31.trinuc$taa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agt.tac<-((dn31.trinuc$agt*dn31.trinuc$tac)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy
#Stop codon
#cps.agt.tag<-((dn31.trinuc$agt*dn31.trinuc$tag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agt.tat<-((dn31.trinuc$agt*dn31.trinuc$tat)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy

cps.agt.tca<-((dn31.trinuc$agt*dn31.trinuc$tca)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.agt.tcc<-((dn31.trinuc$agt*dn31.trinuc$tcc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.agt.tcg<-((dn31.trinuc$agt*dn31.trinuc$tcg)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.agt.tct<-((dn31.trinuc$agt*dn31.trinuc$tct)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

#Stop codon
#cps.agt.tga<-((dn31.trinuc$agt*dn31.trinuc$tga)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$kk
cps.agt.tgc<-((dn31.trinuc$agt*dn31.trinuc$tgc)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc
cps.agt.tgg<-((dn31.trinuc$agt*dn31.trinuc$tgg)/(dn31.count.aa$s*dn31.count.aa$w))/dn31.count.diaa$sw
cps.agt.tgt<-((dn31.trinuc$agt*dn31.trinuc$tgt)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc

cps.agt.tta<-((dn31.trinuc$agt*dn31.trinuc$tta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agt.ttc<-((dn31.trinuc$agt*dn31.trinuc$ttc)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf
cps.agt.ttg<-((dn31.trinuc$agt*dn31.trinuc$ttg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.agt.ttt<-((dn31.trinuc$agt*dn31.trinuc$ttt)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf













cps.ata.aaa<-((dn31.trinuc$ata*dn31.trinuc$aaa)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.ata.aac<-((dn31.trinuc$ata*dn31.trinuc$aac)/(dn31.count.aa$i*dn31.count.aa$n))/dn31.count.diaa$"in"
cps.ata.aag<-((dn31.trinuc$ata*dn31.trinuc$aag)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.ata.aat<-((dn31.trinuc$ata*dn31.trinuc$aat)/(dn31.count.aa$i*dn31.count.aa$n))/dn31.count.diaa$"in"

cps.ata.aca<-((dn31.trinuc$ata*dn31.trinuc$aca)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.ata.acc<-((dn31.trinuc$ata*dn31.trinuc$acc)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.ata.acg<-((dn31.trinuc$ata*dn31.trinuc$acg)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.ata.act<-((dn31.trinuc$ata*dn31.trinuc$act)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it

cps.ata.aga<-((dn31.trinuc$ata*dn31.trinuc$aga)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.ata.agc<-((dn31.trinuc$ata*dn31.trinuc$agc)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.ata.agg<-((dn31.trinuc$ata*dn31.trinuc$agg)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.ata.agt<-((dn31.trinuc$ata*dn31.trinuc$agt)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is

cps.ata.ata<-((dn31.trinuc$ata*dn31.trinuc$ata)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii
cps.ata.atc<-((dn31.trinuc$ata*dn31.trinuc$atc)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii
cps.ata.atg<-((dn31.trinuc$ata*dn31.trinuc$atg)/(dn31.count.aa$i*dn31.count.aa$m))/dn31.count.diaa$im
cps.ata.att<-((dn31.trinuc$ata*dn31.trinuc$att)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii

cps.ata.caa<-((dn31.trinuc$ata*dn31.trinuc$caa)/(dn31.count.aa$i*dn31.count.aa$q))/dn31.count.diaa$iq
cps.ata.cac<-((dn31.trinuc$ata*dn31.trinuc$cac)/(dn31.count.aa$i*dn31.count.aa$h))/dn31.count.diaa$ih
cps.ata.cag<-((dn31.trinuc$ata*dn31.trinuc$cag)/(dn31.count.aa$i*dn31.count.aa$q))/dn31.count.diaa$iq
cps.ata.cat<-((dn31.trinuc$ata*dn31.trinuc$cat)/(dn31.count.aa$i*dn31.count.aa$h))/dn31.count.diaa$ih

cps.ata.cca<-((dn31.trinuc$ata*dn31.trinuc$cca)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.ata.ccc<-((dn31.trinuc$ata*dn31.trinuc$ccc)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.ata.ccg<-((dn31.trinuc$ata*dn31.trinuc$ccg)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.ata.cct<-((dn31.trinuc$ata*dn31.trinuc$cct)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip

cps.ata.cga<-((dn31.trinuc$ata*dn31.trinuc$cga)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.ata.cgc<-((dn31.trinuc$ata*dn31.trinuc$cgc)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.ata.cgg<-((dn31.trinuc$ata*dn31.trinuc$cgg)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.ata.cgt<-((dn31.trinuc$ata*dn31.trinuc$cgt)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir

cps.ata.cta<-((dn31.trinuc$ata*dn31.trinuc$cta)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.ata.ctc<-((dn31.trinuc$ata*dn31.trinuc$ctc)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.ata.ctg<-((dn31.trinuc$ata*dn31.trinuc$ctg)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.ata.ctt<-((dn31.trinuc$ata*dn31.trinuc$ctt)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il

cps.ata.gaa<-((dn31.trinuc$ata*dn31.trinuc$gaa)/(dn31.count.aa$i*dn31.count.aa$e))/dn31.count.diaa$ie
cps.ata.gac<-((dn31.trinuc$ata*dn31.trinuc$gac)/(dn31.count.aa$i*dn31.count.aa$d))/dn31.count.diaa$id
cps.ata.gag<-((dn31.trinuc$ata*dn31.trinuc$gag)/(dn31.count.aa$i*dn31.count.aa$e))/dn31.count.diaa$ie
cps.ata.gat<-((dn31.trinuc$ata*dn31.trinuc$gat)/(dn31.count.aa$i*dn31.count.aa$d))/dn31.count.diaa$id

cps.ata.gca<-((dn31.trinuc$ata*dn31.trinuc$gca)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.ata.gcc<-((dn31.trinuc$ata*dn31.trinuc$gcc)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.ata.gcg<-((dn31.trinuc$ata*dn31.trinuc$gcg)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.ata.gct<-((dn31.trinuc$ata*dn31.trinuc$gct)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia

cps.ata.gga<-((dn31.trinuc$ata*dn31.trinuc$gga)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.ata.ggc<-((dn31.trinuc$ata*dn31.trinuc$ggc)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.ata.ggg<-((dn31.trinuc$ata*dn31.trinuc$ggg)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.ata.ggt<-((dn31.trinuc$ata*dn31.trinuc$ggt)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig

cps.ata.gta<-((dn31.trinuc$ata*dn31.trinuc$gta)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.ata.gtc<-((dn31.trinuc$ata*dn31.trinuc$gtc)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.ata.gtg<-((dn31.trinuc$ata*dn31.trinuc$gtg)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.ata.gtt<-((dn31.trinuc$ata*dn31.trinuc$gtt)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv

#Stop codon
#cps.ata.taa<-((dn31.trinuc$ata*dn31.trinuc$taa)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.ata.tac<-((dn31.trinuc$ata*dn31.trinuc$tac)/(dn31.count.aa$i*dn31.count.aa$y))/dn31.count.diaa$iy
#Stop codon
#cps.ata.tag<-((dn31.trinuc$ata*dn31.trinuc$tag)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.ata.tat<-((dn31.trinuc$ata*dn31.trinuc$tat)/(dn31.count.aa$i*dn31.count.aa$y))/dn31.count.diaa$iy

cps.ata.tca<-((dn31.trinuc$ata*dn31.trinuc$tca)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.ata.tcc<-((dn31.trinuc$ata*dn31.trinuc$tcc)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.ata.tcg<-((dn31.trinuc$ata*dn31.trinuc$tcg)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.ata.tct<-((dn31.trinuc$ata*dn31.trinuc$tct)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is

#Stop codon
#cps.ata.tga<-((dn31.trinuc$ata*dn31.trinuc$tga)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.ata.tgc<-((dn31.trinuc$ata*dn31.trinuc$tgc)/(dn31.count.aa$i*dn31.count.aa$c))/dn31.count.diaa$ic
cps.ata.tgg<-((dn31.trinuc$ata*dn31.trinuc$tgg)/(dn31.count.aa$i*dn31.count.aa$w))/dn31.count.diaa$iw
cps.ata.tgt<-((dn31.trinuc$ata*dn31.trinuc$tgt)/(dn31.count.aa$i*dn31.count.aa$c))/dn31.count.diaa$ic

cps.ata.tta<-((dn31.trinuc$ata*dn31.trinuc$tta)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.ata.ttc<-((dn31.trinuc$ata*dn31.trinuc$ttc)/(dn31.count.aa$i*dn31.count.aa$f))/dn31.count.diaa$"if"
cps.ata.ttg<-((dn31.trinuc$ata*dn31.trinuc$ttg)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.ata.ttt<-((dn31.trinuc$ata*dn31.trinuc$ttt)/(dn31.count.aa$i*dn31.count.aa$f))/dn31.count.diaa$"if"










cps.atc.aaa<-((dn31.trinuc$atc*dn31.trinuc$aaa)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.atc.aac<-((dn31.trinuc$atc*dn31.trinuc$aac)/(dn31.count.aa$i*dn31.count.aa$n))/dn31.count.diaa$"in"
cps.atc.aag<-((dn31.trinuc$atc*dn31.trinuc$aag)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.atc.aat<-((dn31.trinuc$atc*dn31.trinuc$aat)/(dn31.count.aa$i*dn31.count.aa$n))/dn31.count.diaa$"in"

cps.atc.aca<-((dn31.trinuc$atc*dn31.trinuc$aca)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.atc.acc<-((dn31.trinuc$atc*dn31.trinuc$acc)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.atc.acg<-((dn31.trinuc$atc*dn31.trinuc$acg)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.atc.act<-((dn31.trinuc$atc*dn31.trinuc$act)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it

cps.atc.aga<-((dn31.trinuc$atc*dn31.trinuc$aga)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.atc.agc<-((dn31.trinuc$atc*dn31.trinuc$agc)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.atc.agg<-((dn31.trinuc$atc*dn31.trinuc$agg)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.atc.agt<-((dn31.trinuc$atc*dn31.trinuc$agt)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is

cps.atc.ata<-((dn31.trinuc$atc*dn31.trinuc$ata)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii
cps.atc.atc<-((dn31.trinuc$atc*dn31.trinuc$atc)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii
cps.atc.atg<-((dn31.trinuc$atc*dn31.trinuc$atg)/(dn31.count.aa$i*dn31.count.aa$m))/dn31.count.diaa$im
cps.atc.att<-((dn31.trinuc$atc*dn31.trinuc$att)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii

cps.atc.caa<-((dn31.trinuc$atc*dn31.trinuc$caa)/(dn31.count.aa$i*dn31.count.aa$q))/dn31.count.diaa$iq
cps.atc.cac<-((dn31.trinuc$atc*dn31.trinuc$cac)/(dn31.count.aa$i*dn31.count.aa$h))/dn31.count.diaa$ih
cps.atc.cag<-((dn31.trinuc$atc*dn31.trinuc$cag)/(dn31.count.aa$i*dn31.count.aa$q))/dn31.count.diaa$iq
cps.atc.cat<-((dn31.trinuc$atc*dn31.trinuc$cat)/(dn31.count.aa$i*dn31.count.aa$h))/dn31.count.diaa$ih

cps.atc.cca<-((dn31.trinuc$atc*dn31.trinuc$cca)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.atc.ccc<-((dn31.trinuc$atc*dn31.trinuc$ccc)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.atc.ccg<-((dn31.trinuc$atc*dn31.trinuc$ccg)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.atc.cct<-((dn31.trinuc$atc*dn31.trinuc$cct)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip

cps.atc.cga<-((dn31.trinuc$atc*dn31.trinuc$cga)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.atc.cgc<-((dn31.trinuc$atc*dn31.trinuc$cgc)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.atc.cgg<-((dn31.trinuc$atc*dn31.trinuc$cgg)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.atc.cgt<-((dn31.trinuc$atc*dn31.trinuc$cgt)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir

cps.atc.cta<-((dn31.trinuc$atc*dn31.trinuc$cta)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.atc.ctc<-((dn31.trinuc$atc*dn31.trinuc$ctc)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.atc.ctg<-((dn31.trinuc$atc*dn31.trinuc$ctg)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.atc.ctt<-((dn31.trinuc$atc*dn31.trinuc$ctt)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il

cps.atc.gaa<-((dn31.trinuc$atc*dn31.trinuc$gaa)/(dn31.count.aa$i*dn31.count.aa$e))/dn31.count.diaa$ie
cps.atc.gac<-((dn31.trinuc$atc*dn31.trinuc$gac)/(dn31.count.aa$i*dn31.count.aa$d))/dn31.count.diaa$id
cps.atc.gag<-((dn31.trinuc$atc*dn31.trinuc$gag)/(dn31.count.aa$i*dn31.count.aa$e))/dn31.count.diaa$ie
cps.atc.gat<-((dn31.trinuc$atc*dn31.trinuc$gat)/(dn31.count.aa$i*dn31.count.aa$d))/dn31.count.diaa$id

cps.atc.gca<-((dn31.trinuc$atc*dn31.trinuc$gca)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.atc.gcc<-((dn31.trinuc$atc*dn31.trinuc$gcc)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.atc.gcg<-((dn31.trinuc$atc*dn31.trinuc$gcg)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.atc.gct<-((dn31.trinuc$atc*dn31.trinuc$gct)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia

cps.atc.gga<-((dn31.trinuc$atc*dn31.trinuc$gga)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.atc.ggc<-((dn31.trinuc$atc*dn31.trinuc$ggc)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.atc.ggg<-((dn31.trinuc$atc*dn31.trinuc$ggg)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.atc.ggt<-((dn31.trinuc$atc*dn31.trinuc$ggt)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig

cps.atc.gta<-((dn31.trinuc$atc*dn31.trinuc$gta)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.atc.gtc<-((dn31.trinuc$atc*dn31.trinuc$gtc)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.atc.gtg<-((dn31.trinuc$atc*dn31.trinuc$gtg)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.atc.gtt<-((dn31.trinuc$atc*dn31.trinuc$gtt)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv

#Stop codon
#cps.atc.taa<-((dn31.trinuc$atc*dn31.trinuc$taa)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.atc.tac<-((dn31.trinuc$atc*dn31.trinuc$tac)/(dn31.count.aa$i*dn31.count.aa$y))/dn31.count.diaa$iy
#Stop codon
#cps.atc.tag<-((dn31.trinuc$atc*dn31.trinuc$tag)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.atc.tat<-((dn31.trinuc$atc*dn31.trinuc$tat)/(dn31.count.aa$i*dn31.count.aa$y))/dn31.count.diaa$iy

cps.atc.tca<-((dn31.trinuc$atc*dn31.trinuc$tca)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.atc.tcc<-((dn31.trinuc$atc*dn31.trinuc$tcc)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.atc.tcg<-((dn31.trinuc$atc*dn31.trinuc$tcg)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.atc.tct<-((dn31.trinuc$atc*dn31.trinuc$tct)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is

#Stop codon
#cps.atc.tga<-((dn31.trinuc$atc*dn31.trinuc$tga)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.atc.tgc<-((dn31.trinuc$atc*dn31.trinuc$tgc)/(dn31.count.aa$i*dn31.count.aa$c))/dn31.count.diaa$ic
cps.atc.tgg<-((dn31.trinuc$atc*dn31.trinuc$tgg)/(dn31.count.aa$i*dn31.count.aa$w))/dn31.count.diaa$iw
cps.atc.tgt<-((dn31.trinuc$atc*dn31.trinuc$tgt)/(dn31.count.aa$i*dn31.count.aa$c))/dn31.count.diaa$ic

cps.atc.tta<-((dn31.trinuc$atc*dn31.trinuc$tta)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.atc.ttc<-((dn31.trinuc$atc*dn31.trinuc$ttc)/(dn31.count.aa$i*dn31.count.aa$f))/dn31.count.diaa$"if"
cps.atc.ttg<-((dn31.trinuc$atc*dn31.trinuc$ttg)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.atc.ttt<-((dn31.trinuc$atc*dn31.trinuc$ttt)/(dn31.count.aa$i*dn31.count.aa$f))/dn31.count.diaa$"if"










cps.atg.aaa<-((dn31.trinuc$atg*dn31.trinuc$aaa)/(dn31.count.aa$m*dn31.count.aa$k))/dn31.count.diaa$mk
cps.atg.aac<-((dn31.trinuc$atg*dn31.trinuc$aac)/(dn31.count.aa$m*dn31.count.aa$n))/dn31.count.diaa$mn
cps.atg.aag<-((dn31.trinuc$atg*dn31.trinuc$aag)/(dn31.count.aa$m*dn31.count.aa$k))/dn31.count.diaa$mk
cps.atg.aat<-((dn31.trinuc$atg*dn31.trinuc$aat)/(dn31.count.aa$m*dn31.count.aa$n))/dn31.count.diaa$mn

cps.atg.aca<-((dn31.trinuc$atg*dn31.trinuc$aca)/(dn31.count.aa$m*dn31.count.aa$t))/dn31.count.diaa$mt
cps.atg.acc<-((dn31.trinuc$atg*dn31.trinuc$acc)/(dn31.count.aa$m*dn31.count.aa$t))/dn31.count.diaa$mt
cps.atg.acg<-((dn31.trinuc$atg*dn31.trinuc$acg)/(dn31.count.aa$m*dn31.count.aa$t))/dn31.count.diaa$mt
cps.atg.act<-((dn31.trinuc$atg*dn31.trinuc$act)/(dn31.count.aa$m*dn31.count.aa$t))/dn31.count.diaa$mt

cps.atg.aga<-((dn31.trinuc$atg*dn31.trinuc$aga)/(dn31.count.aa$m*dn31.count.aa$r))/dn31.count.diaa$mr
cps.atg.agc<-((dn31.trinuc$atg*dn31.trinuc$agc)/(dn31.count.aa$m*dn31.count.aa$s))/dn31.count.diaa$ms
cps.atg.agg<-((dn31.trinuc$atg*dn31.trinuc$agg)/(dn31.count.aa$m*dn31.count.aa$r))/dn31.count.diaa$mr
cps.atg.agt<-((dn31.trinuc$atg*dn31.trinuc$agt)/(dn31.count.aa$m*dn31.count.aa$s))/dn31.count.diaa$ms

cps.atg.ata<-((dn31.trinuc$atg*dn31.trinuc$ata)/(dn31.count.aa$m*dn31.count.aa$i))/dn31.count.diaa$mi
cps.atg.atc<-((dn31.trinuc$atg*dn31.trinuc$atc)/(dn31.count.aa$m*dn31.count.aa$i))/dn31.count.diaa$mi
cps.atg.atg<-((dn31.trinuc$atg*dn31.trinuc$atg)/(dn31.count.aa$m*dn31.count.aa$m))/dn31.count.diaa$mm
cps.atg.att<-((dn31.trinuc$atg*dn31.trinuc$att)/(dn31.count.aa$m*dn31.count.aa$i))/dn31.count.diaa$mi

cps.atg.caa<-((dn31.trinuc$atg*dn31.trinuc$caa)/(dn31.count.aa$m*dn31.count.aa$q))/dn31.count.diaa$mq
cps.atg.cac<-((dn31.trinuc$atg*dn31.trinuc$cac)/(dn31.count.aa$m*dn31.count.aa$h))/dn31.count.diaa$mh
cps.atg.cag<-((dn31.trinuc$atg*dn31.trinuc$cag)/(dn31.count.aa$m*dn31.count.aa$q))/dn31.count.diaa$mq
cps.atg.cat<-((dn31.trinuc$atg*dn31.trinuc$cat)/(dn31.count.aa$m*dn31.count.aa$h))/dn31.count.diaa$mh

cps.atg.cca<-((dn31.trinuc$atg*dn31.trinuc$cca)/(dn31.count.aa$m*dn31.count.aa$p))/dn31.count.diaa$mp
cps.atg.ccc<-((dn31.trinuc$atg*dn31.trinuc$ccc)/(dn31.count.aa$m*dn31.count.aa$p))/dn31.count.diaa$mp
cps.atg.ccg<-((dn31.trinuc$atg*dn31.trinuc$ccg)/(dn31.count.aa$m*dn31.count.aa$p))/dn31.count.diaa$mp
cps.atg.cct<-((dn31.trinuc$atg*dn31.trinuc$cct)/(dn31.count.aa$m*dn31.count.aa$p))/dn31.count.diaa$mp

cps.atg.cga<-((dn31.trinuc$atg*dn31.trinuc$cga)/(dn31.count.aa$m*dn31.count.aa$r))/dn31.count.diaa$mr
cps.atg.cgc<-((dn31.trinuc$atg*dn31.trinuc$cgc)/(dn31.count.aa$m*dn31.count.aa$r))/dn31.count.diaa$mr
cps.atg.cgg<-((dn31.trinuc$atg*dn31.trinuc$cgg)/(dn31.count.aa$m*dn31.count.aa$r))/dn31.count.diaa$mr
cps.atg.cgt<-((dn31.trinuc$atg*dn31.trinuc$cgt)/(dn31.count.aa$m*dn31.count.aa$r))/dn31.count.diaa$mr

cps.atg.cta<-((dn31.trinuc$atg*dn31.trinuc$cta)/(dn31.count.aa$m*dn31.count.aa$l))/dn31.count.diaa$ml
cps.atg.ctc<-((dn31.trinuc$atg*dn31.trinuc$ctc)/(dn31.count.aa$m*dn31.count.aa$l))/dn31.count.diaa$ml
cps.atg.ctg<-((dn31.trinuc$atg*dn31.trinuc$ctg)/(dn31.count.aa$m*dn31.count.aa$l))/dn31.count.diaa$ml
cps.atg.ctt<-((dn31.trinuc$atg*dn31.trinuc$ctt)/(dn31.count.aa$m*dn31.count.aa$l))/dn31.count.diaa$ml

cps.atg.gaa<-((dn31.trinuc$atg*dn31.trinuc$gaa)/(dn31.count.aa$m*dn31.count.aa$e))/dn31.count.diaa$me
cps.atg.gac<-((dn31.trinuc$atg*dn31.trinuc$gac)/(dn31.count.aa$m*dn31.count.aa$d))/dn31.count.diaa$md
cps.atg.gag<-((dn31.trinuc$atg*dn31.trinuc$gag)/(dn31.count.aa$m*dn31.count.aa$e))/dn31.count.diaa$me
cps.atg.gat<-((dn31.trinuc$atg*dn31.trinuc$gat)/(dn31.count.aa$m*dn31.count.aa$d))/dn31.count.diaa$md

cps.atg.gca<-((dn31.trinuc$atg*dn31.trinuc$gca)/(dn31.count.aa$m*dn31.count.aa$a))/dn31.count.diaa$ma
cps.atg.gcc<-((dn31.trinuc$atg*dn31.trinuc$gcc)/(dn31.count.aa$m*dn31.count.aa$a))/dn31.count.diaa$ma
cps.atg.gcg<-((dn31.trinuc$atg*dn31.trinuc$gcg)/(dn31.count.aa$m*dn31.count.aa$a))/dn31.count.diaa$ma
cps.atg.gct<-((dn31.trinuc$atg*dn31.trinuc$gct)/(dn31.count.aa$m*dn31.count.aa$a))/dn31.count.diaa$ma

cps.atg.gga<-((dn31.trinuc$atg*dn31.trinuc$gga)/(dn31.count.aa$m*dn31.count.aa$g))/dn31.count.diaa$mg
cps.atg.ggc<-((dn31.trinuc$atg*dn31.trinuc$ggc)/(dn31.count.aa$m*dn31.count.aa$g))/dn31.count.diaa$mg
cps.atg.ggg<-((dn31.trinuc$atg*dn31.trinuc$ggg)/(dn31.count.aa$m*dn31.count.aa$g))/dn31.count.diaa$mg
cps.atg.ggt<-((dn31.trinuc$atg*dn31.trinuc$ggt)/(dn31.count.aa$m*dn31.count.aa$g))/dn31.count.diaa$mg

cps.atg.gta<-((dn31.trinuc$atg*dn31.trinuc$gta)/(dn31.count.aa$m*dn31.count.aa$v))/dn31.count.diaa$mv
cps.atg.gtc<-((dn31.trinuc$atg*dn31.trinuc$gtc)/(dn31.count.aa$m*dn31.count.aa$v))/dn31.count.diaa$mv
cps.atg.gtg<-((dn31.trinuc$atg*dn31.trinuc$gtg)/(dn31.count.aa$m*dn31.count.aa$v))/dn31.count.diaa$mv
cps.atg.gtt<-((dn31.trinuc$atg*dn31.trinuc$gtt)/(dn31.count.aa$m*dn31.count.aa$v))/dn31.count.diaa$mv

#Stop codon
#cps.atg.taa<-((dn31.trinuc$atg*dn31.trinuc$taa)/(dn31.count.aa$m*dn31.count.aa$k))/dn31.count.diaa$mk
cps.atg.tac<-((dn31.trinuc$atg*dn31.trinuc$tac)/(dn31.count.aa$m*dn31.count.aa$y))/dn31.count.diaa$my
#Stop codon
#cps.atg.tag<-((dn31.trinuc$atg*dn31.trinuc$tag)/(dn31.count.aa$m*dn31.count.aa$k))/dn31.count.diaa$mk
cps.atg.tat<-((dn31.trinuc$atg*dn31.trinuc$tat)/(dn31.count.aa$m*dn31.count.aa$y))/dn31.count.diaa$my

cps.atg.tca<-((dn31.trinuc$atg*dn31.trinuc$tca)/(dn31.count.aa$m*dn31.count.aa$s))/dn31.count.diaa$ms
cps.atg.tcc<-((dn31.trinuc$atg*dn31.trinuc$tcc)/(dn31.count.aa$m*dn31.count.aa$s))/dn31.count.diaa$ms
cps.atg.tcg<-((dn31.trinuc$atg*dn31.trinuc$tcg)/(dn31.count.aa$m*dn31.count.aa$s))/dn31.count.diaa$ms
cps.atg.tct<-((dn31.trinuc$atg*dn31.trinuc$tct)/(dn31.count.aa$m*dn31.count.aa$s))/dn31.count.diaa$ms

#Stop codon
#cps.atg.tga<-((dn31.trinuc$atg*dn31.trinuc$tga)/(dn31.count.aa$m*dn31.count.aa$k))/dn31.count.diaa$mk
cps.atg.tgc<-((dn31.trinuc$atg*dn31.trinuc$tgc)/(dn31.count.aa$m*dn31.count.aa$c))/dn31.count.diaa$mc
cps.atg.tgg<-((dn31.trinuc$atg*dn31.trinuc$tgg)/(dn31.count.aa$m*dn31.count.aa$w))/dn31.count.diaa$mw
cps.atg.tgt<-((dn31.trinuc$atg*dn31.trinuc$tgt)/(dn31.count.aa$m*dn31.count.aa$c))/dn31.count.diaa$mc

cps.atg.tta<-((dn31.trinuc$atg*dn31.trinuc$tta)/(dn31.count.aa$m*dn31.count.aa$l))/dn31.count.diaa$ml
cps.atg.ttc<-((dn31.trinuc$atg*dn31.trinuc$ttc)/(dn31.count.aa$m*dn31.count.aa$f))/dn31.count.diaa$mf
cps.atg.ttg<-((dn31.trinuc$atg*dn31.trinuc$ttg)/(dn31.count.aa$m*dn31.count.aa$l))/dn31.count.diaa$ml
cps.atg.ttt<-((dn31.trinuc$atg*dn31.trinuc$ttt)/(dn31.count.aa$m*dn31.count.aa$f))/dn31.count.diaa$mf








cps.att.aaa<-((dn31.trinuc$att*dn31.trinuc$aaa)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.att.aac<-((dn31.trinuc$att*dn31.trinuc$aac)/(dn31.count.aa$i*dn31.count.aa$n))/dn31.count.diaa$"in"
cps.att.aag<-((dn31.trinuc$att*dn31.trinuc$aag)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.att.aat<-((dn31.trinuc$att*dn31.trinuc$aat)/(dn31.count.aa$i*dn31.count.aa$n))/dn31.count.diaa$"in"

cps.att.aca<-((dn31.trinuc$att*dn31.trinuc$aca)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.att.acc<-((dn31.trinuc$att*dn31.trinuc$acc)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.att.acg<-((dn31.trinuc$att*dn31.trinuc$acg)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it
cps.att.act<-((dn31.trinuc$att*dn31.trinuc$act)/(dn31.count.aa$i*dn31.count.aa$t))/dn31.count.diaa$it

cps.att.aga<-((dn31.trinuc$att*dn31.trinuc$aga)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.att.agc<-((dn31.trinuc$att*dn31.trinuc$agc)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.att.agg<-((dn31.trinuc$att*dn31.trinuc$agg)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.att.agt<-((dn31.trinuc$att*dn31.trinuc$agt)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is

cps.att.ata<-((dn31.trinuc$att*dn31.trinuc$ata)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii
cps.att.atc<-((dn31.trinuc$att*dn31.trinuc$atc)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii
cps.att.atg<-((dn31.trinuc$att*dn31.trinuc$atg)/(dn31.count.aa$i*dn31.count.aa$m))/dn31.count.diaa$im
cps.att.att<-((dn31.trinuc$att*dn31.trinuc$att)/(dn31.count.aa$i*dn31.count.aa$i))/dn31.count.diaa$ii

cps.att.caa<-((dn31.trinuc$att*dn31.trinuc$caa)/(dn31.count.aa$i*dn31.count.aa$q))/dn31.count.diaa$iq
cps.att.cac<-((dn31.trinuc$att*dn31.trinuc$cac)/(dn31.count.aa$i*dn31.count.aa$h))/dn31.count.diaa$ih
cps.att.cag<-((dn31.trinuc$att*dn31.trinuc$cag)/(dn31.count.aa$i*dn31.count.aa$q))/dn31.count.diaa$iq
cps.att.cat<-((dn31.trinuc$att*dn31.trinuc$cat)/(dn31.count.aa$i*dn31.count.aa$h))/dn31.count.diaa$ih

cps.att.cca<-((dn31.trinuc$att*dn31.trinuc$cca)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.att.ccc<-((dn31.trinuc$att*dn31.trinuc$ccc)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.att.ccg<-((dn31.trinuc$att*dn31.trinuc$ccg)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip
cps.att.cct<-((dn31.trinuc$att*dn31.trinuc$cct)/(dn31.count.aa$i*dn31.count.aa$p))/dn31.count.diaa$ip

cps.att.cga<-((dn31.trinuc$att*dn31.trinuc$cga)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.att.cgc<-((dn31.trinuc$att*dn31.trinuc$cgc)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.att.cgg<-((dn31.trinuc$att*dn31.trinuc$cgg)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir
cps.att.cgt<-((dn31.trinuc$att*dn31.trinuc$cgt)/(dn31.count.aa$i*dn31.count.aa$r))/dn31.count.diaa$ir

cps.att.cta<-((dn31.trinuc$att*dn31.trinuc$cta)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.att.ctc<-((dn31.trinuc$att*dn31.trinuc$ctc)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.att.ctg<-((dn31.trinuc$att*dn31.trinuc$ctg)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.att.ctt<-((dn31.trinuc$att*dn31.trinuc$ctt)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il

cps.att.gaa<-((dn31.trinuc$att*dn31.trinuc$gaa)/(dn31.count.aa$i*dn31.count.aa$e))/dn31.count.diaa$ie
cps.att.gac<-((dn31.trinuc$att*dn31.trinuc$gac)/(dn31.count.aa$i*dn31.count.aa$d))/dn31.count.diaa$id
cps.att.gag<-((dn31.trinuc$att*dn31.trinuc$gag)/(dn31.count.aa$i*dn31.count.aa$e))/dn31.count.diaa$ie
cps.att.gat<-((dn31.trinuc$att*dn31.trinuc$gat)/(dn31.count.aa$i*dn31.count.aa$d))/dn31.count.diaa$id

cps.att.gca<-((dn31.trinuc$att*dn31.trinuc$gca)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.att.gcc<-((dn31.trinuc$att*dn31.trinuc$gcc)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.att.gcg<-((dn31.trinuc$att*dn31.trinuc$gcg)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia
cps.att.gct<-((dn31.trinuc$att*dn31.trinuc$gct)/(dn31.count.aa$i*dn31.count.aa$a))/dn31.count.diaa$ia

cps.att.gga<-((dn31.trinuc$att*dn31.trinuc$gga)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.att.ggc<-((dn31.trinuc$att*dn31.trinuc$ggc)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.att.ggg<-((dn31.trinuc$att*dn31.trinuc$ggg)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig
cps.att.ggt<-((dn31.trinuc$att*dn31.trinuc$ggt)/(dn31.count.aa$i*dn31.count.aa$g))/dn31.count.diaa$ig

cps.att.gta<-((dn31.trinuc$att*dn31.trinuc$gta)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.att.gtc<-((dn31.trinuc$att*dn31.trinuc$gtc)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.att.gtg<-((dn31.trinuc$att*dn31.trinuc$gtg)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv
cps.att.gtt<-((dn31.trinuc$att*dn31.trinuc$gtt)/(dn31.count.aa$i*dn31.count.aa$v))/dn31.count.diaa$iv

#Stop codon
#cps.att.taa<-((dn31.trinuc$att*dn31.trinuc$taa)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.att.tac<-((dn31.trinuc$att*dn31.trinuc$tac)/(dn31.count.aa$i*dn31.count.aa$y))/dn31.count.diaa$iy
#Stop codon
#cps.att.tag<-((dn31.trinuc$att*dn31.trinuc$tag)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.att.tat<-((dn31.trinuc$att*dn31.trinuc$tat)/(dn31.count.aa$i*dn31.count.aa$y))/dn31.count.diaa$iy

cps.att.tca<-((dn31.trinuc$att*dn31.trinuc$tca)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.att.tcc<-((dn31.trinuc$att*dn31.trinuc$tcc)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.att.tcg<-((dn31.trinuc$att*dn31.trinuc$tcg)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is
cps.att.tct<-((dn31.trinuc$att*dn31.trinuc$tct)/(dn31.count.aa$i*dn31.count.aa$s))/dn31.count.diaa$is

#Stop codon
#cps.att.tga<-((dn31.trinuc$att*dn31.trinuc$tga)/(dn31.count.aa$i*dn31.count.aa$k))/dn31.count.diaa$ik
cps.att.tgc<-((dn31.trinuc$att*dn31.trinuc$tgc)/(dn31.count.aa$i*dn31.count.aa$c))/dn31.count.diaa$ic
cps.att.tgg<-((dn31.trinuc$att*dn31.trinuc$tgg)/(dn31.count.aa$i*dn31.count.aa$w))/dn31.count.diaa$iw
cps.att.tgt<-((dn31.trinuc$att*dn31.trinuc$tgt)/(dn31.count.aa$i*dn31.count.aa$c))/dn31.count.diaa$ic

cps.att.tta<-((dn31.trinuc$att*dn31.trinuc$tta)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.att.ttc<-((dn31.trinuc$att*dn31.trinuc$ttc)/(dn31.count.aa$i*dn31.count.aa$f))/dn31.count.diaa$"if"
cps.att.ttg<-((dn31.trinuc$att*dn31.trinuc$ttg)/(dn31.count.aa$i*dn31.count.aa$l))/dn31.count.diaa$il
cps.att.ttt<-((dn31.trinuc$att*dn31.trinuc$ttt)/(dn31.count.aa$i*dn31.count.aa$f))/dn31.count.diaa$"if"








cps.caa.aaa<-((dn31.trinuc$caa*dn31.trinuc$aaa)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.caa.aac<-((dn31.trinuc$caa*dn31.trinuc$aac)/(dn31.count.aa$q*dn31.count.aa$n))/dn31.count.diaa$qn
cps.caa.aag<-((dn31.trinuc$caa*dn31.trinuc$aag)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.caa.aat<-((dn31.trinuc$caa*dn31.trinuc$aat)/(dn31.count.aa$q*dn31.count.aa$n))/dn31.count.diaa$qn

cps.caa.aca<-((dn31.trinuc$caa*dn31.trinuc$aca)/(dn31.count.aa$q*dn31.count.aa$t))/dn31.count.diaa$qt
cps.caa.acc<-((dn31.trinuc$caa*dn31.trinuc$acc)/(dn31.count.aa$q*dn31.count.aa$t))/dn31.count.diaa$qt
cps.caa.acg<-((dn31.trinuc$caa*dn31.trinuc$acg)/(dn31.count.aa$q*dn31.count.aa$t))/dn31.count.diaa$qt
cps.caa.act<-((dn31.trinuc$caa*dn31.trinuc$act)/(dn31.count.aa$q*dn31.count.aa$t))/dn31.count.diaa$qt

cps.caa.aga<-((dn31.trinuc$caa*dn31.trinuc$aga)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.caa.agc<-((dn31.trinuc$caa*dn31.trinuc$agc)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs
cps.caa.agg<-((dn31.trinuc$caa*dn31.trinuc$agg)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.caa.agt<-((dn31.trinuc$caa*dn31.trinuc$agt)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs

cps.caa.ata<-((dn31.trinuc$caa*dn31.trinuc$ata)/(dn31.count.aa$q*dn31.count.aa$i))/dn31.count.diaa$qi
cps.caa.atc<-((dn31.trinuc$caa*dn31.trinuc$atc)/(dn31.count.aa$q*dn31.count.aa$i))/dn31.count.diaa$qi
cps.caa.atg<-((dn31.trinuc$caa*dn31.trinuc$atg)/(dn31.count.aa$q*dn31.count.aa$m))/dn31.count.diaa$qm
cps.caa.att<-((dn31.trinuc$caa*dn31.trinuc$att)/(dn31.count.aa$q*dn31.count.aa$i))/dn31.count.diaa$qi

cps.caa.caa<-((dn31.trinuc$caa*dn31.trinuc$caa)/(dn31.count.aa$q*dn31.count.aa$q))/dn31.count.diaa$qq
cps.caa.cac<-((dn31.trinuc$caa*dn31.trinuc$cac)/(dn31.count.aa$q*dn31.count.aa$h))/dn31.count.diaa$qh
cps.caa.cag<-((dn31.trinuc$caa*dn31.trinuc$cag)/(dn31.count.aa$q*dn31.count.aa$q))/dn31.count.diaa$qq
cps.caa.cat<-((dn31.trinuc$caa*dn31.trinuc$cat)/(dn31.count.aa$q*dn31.count.aa$h))/dn31.count.diaa$qh

cps.caa.cca<-((dn31.trinuc$caa*dn31.trinuc$cca)/(dn31.count.aa$q*dn31.count.aa$p))/dn31.count.diaa$qp
cps.caa.ccc<-((dn31.trinuc$caa*dn31.trinuc$ccc)/(dn31.count.aa$q*dn31.count.aa$p))/dn31.count.diaa$qp
cps.caa.ccg<-((dn31.trinuc$caa*dn31.trinuc$ccg)/(dn31.count.aa$q*dn31.count.aa$p))/dn31.count.diaa$qp
cps.caa.cct<-((dn31.trinuc$caa*dn31.trinuc$cct)/(dn31.count.aa$q*dn31.count.aa$p))/dn31.count.diaa$qp

cps.caa.cga<-((dn31.trinuc$caa*dn31.trinuc$cga)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.caa.cgc<-((dn31.trinuc$caa*dn31.trinuc$cgc)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.caa.cgg<-((dn31.trinuc$caa*dn31.trinuc$cgg)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.caa.cgt<-((dn31.trinuc$caa*dn31.trinuc$cgt)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr

cps.caa.cta<-((dn31.trinuc$caa*dn31.trinuc$cta)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.caa.ctc<-((dn31.trinuc$caa*dn31.trinuc$ctc)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.caa.ctg<-((dn31.trinuc$caa*dn31.trinuc$ctg)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.caa.ctt<-((dn31.trinuc$caa*dn31.trinuc$ctt)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql

cps.caa.gaa<-((dn31.trinuc$caa*dn31.trinuc$gaa)/(dn31.count.aa$q*dn31.count.aa$e))/dn31.count.diaa$qe
cps.caa.gac<-((dn31.trinuc$caa*dn31.trinuc$gac)/(dn31.count.aa$q*dn31.count.aa$d))/dn31.count.diaa$qd
cps.caa.gag<-((dn31.trinuc$caa*dn31.trinuc$gag)/(dn31.count.aa$q*dn31.count.aa$e))/dn31.count.diaa$qe
cps.caa.gat<-((dn31.trinuc$caa*dn31.trinuc$gat)/(dn31.count.aa$q*dn31.count.aa$d))/dn31.count.diaa$qd

cps.caa.gca<-((dn31.trinuc$caa*dn31.trinuc$gca)/(dn31.count.aa$q*dn31.count.aa$a))/dn31.count.diaa$qa
cps.caa.gcc<-((dn31.trinuc$caa*dn31.trinuc$gcc)/(dn31.count.aa$q*dn31.count.aa$a))/dn31.count.diaa$qa
cps.caa.gcg<-((dn31.trinuc$caa*dn31.trinuc$gcg)/(dn31.count.aa$q*dn31.count.aa$a))/dn31.count.diaa$qa
cps.caa.gct<-((dn31.trinuc$caa*dn31.trinuc$gct)/(dn31.count.aa$q*dn31.count.aa$a))/dn31.count.diaa$qa

cps.caa.gga<-((dn31.trinuc$caa*dn31.trinuc$gga)/(dn31.count.aa$q*dn31.count.aa$g))/dn31.count.diaa$qg
cps.caa.ggc<-((dn31.trinuc$caa*dn31.trinuc$ggc)/(dn31.count.aa$q*dn31.count.aa$g))/dn31.count.diaa$qg
cps.caa.ggg<-((dn31.trinuc$caa*dn31.trinuc$ggg)/(dn31.count.aa$q*dn31.count.aa$g))/dn31.count.diaa$qg
cps.caa.ggt<-((dn31.trinuc$caa*dn31.trinuc$ggt)/(dn31.count.aa$q*dn31.count.aa$g))/dn31.count.diaa$qg

cps.caa.gta<-((dn31.trinuc$caa*dn31.trinuc$gta)/(dn31.count.aa$q*dn31.count.aa$v))/dn31.count.diaa$qv
cps.caa.gtc<-((dn31.trinuc$caa*dn31.trinuc$gtc)/(dn31.count.aa$q*dn31.count.aa$v))/dn31.count.diaa$qv
cps.caa.gtg<-((dn31.trinuc$caa*dn31.trinuc$gtg)/(dn31.count.aa$q*dn31.count.aa$v))/dn31.count.diaa$qv
cps.caa.gtt<-((dn31.trinuc$caa*dn31.trinuc$gtt)/(dn31.count.aa$q*dn31.count.aa$v))/dn31.count.diaa$qv

#Stop codon
#cps.caa.taa<-((dn31.trinuc$caa*dn31.trinuc$taa)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.caa.tac<-((dn31.trinuc$caa*dn31.trinuc$tac)/(dn31.count.aa$q*dn31.count.aa$y))/dn31.count.diaa$qy
#Stop codon
#cps.caa.tag<-((dn31.trinuc$caa*dn31.trinuc$tag)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.caa.tat<-((dn31.trinuc$caa*dn31.trinuc$tat)/(dn31.count.aa$q*dn31.count.aa$y))/dn31.count.diaa$qy

cps.caa.tca<-((dn31.trinuc$caa*dn31.trinuc$tca)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs
cps.caa.tcc<-((dn31.trinuc$caa*dn31.trinuc$tcc)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs
cps.caa.tcg<-((dn31.trinuc$caa*dn31.trinuc$tcg)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs
cps.caa.tct<-((dn31.trinuc$caa*dn31.trinuc$tct)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs

#Stop codon
#cps.caa.tga<-((dn31.trinuc$caa*dn31.trinuc$tga)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.caa.tgc<-((dn31.trinuc$caa*dn31.trinuc$tgc)/(dn31.count.aa$q*dn31.count.aa$c))/dn31.count.diaa$qc
cps.caa.tgg<-((dn31.trinuc$caa*dn31.trinuc$tgg)/(dn31.count.aa$q*dn31.count.aa$w))/dn31.count.diaa$qw
cps.caa.tgt<-((dn31.trinuc$caa*dn31.trinuc$tgt)/(dn31.count.aa$q*dn31.count.aa$c))/dn31.count.diaa$qc

cps.caa.tta<-((dn31.trinuc$caa*dn31.trinuc$tta)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.caa.ttc<-((dn31.trinuc$caa*dn31.trinuc$ttc)/(dn31.count.aa$q*dn31.count.aa$f))/dn31.count.diaa$qf
cps.caa.ttg<-((dn31.trinuc$caa*dn31.trinuc$ttg)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.caa.ttt<-((dn31.trinuc$caa*dn31.trinuc$ttt)/(dn31.count.aa$q*dn31.count.aa$f))/dn31.count.diaa$qf








cps.cac.aaa<-((dn31.trinuc$cac*dn31.trinuc$aaa)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cac.aac<-((dn31.trinuc$cac*dn31.trinuc$aac)/(dn31.count.aa$h*dn31.count.aa$n))/dn31.count.diaa$hn
cps.cac.aag<-((dn31.trinuc$cac*dn31.trinuc$aag)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cac.aat<-((dn31.trinuc$cac*dn31.trinuc$aat)/(dn31.count.aa$h*dn31.count.aa$n))/dn31.count.diaa$hn

cps.cac.aca<-((dn31.trinuc$cac*dn31.trinuc$aca)/(dn31.count.aa$h*dn31.count.aa$t))/dn31.count.diaa$ht
cps.cac.acc<-((dn31.trinuc$cac*dn31.trinuc$acc)/(dn31.count.aa$h*dn31.count.aa$t))/dn31.count.diaa$ht
cps.cac.acg<-((dn31.trinuc$cac*dn31.trinuc$acg)/(dn31.count.aa$h*dn31.count.aa$t))/dn31.count.diaa$ht
cps.cac.act<-((dn31.trinuc$cac*dn31.trinuc$act)/(dn31.count.aa$h*dn31.count.aa$t))/dn31.count.diaa$ht

cps.cac.aga<-((dn31.trinuc$cac*dn31.trinuc$aga)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cac.agc<-((dn31.trinuc$cac*dn31.trinuc$agc)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs
cps.cac.agg<-((dn31.trinuc$cac*dn31.trinuc$agg)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cac.agt<-((dn31.trinuc$cac*dn31.trinuc$agt)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs

cps.cac.ata<-((dn31.trinuc$cac*dn31.trinuc$ata)/(dn31.count.aa$h*dn31.count.aa$i))/dn31.count.diaa$hi
cps.cac.atc<-((dn31.trinuc$cac*dn31.trinuc$atc)/(dn31.count.aa$h*dn31.count.aa$i))/dn31.count.diaa$hi
cps.cac.atg<-((dn31.trinuc$cac*dn31.trinuc$atg)/(dn31.count.aa$h*dn31.count.aa$m))/dn31.count.diaa$hm
cps.cac.att<-((dn31.trinuc$cac*dn31.trinuc$att)/(dn31.count.aa$h*dn31.count.aa$i))/dn31.count.diaa$hi

cps.cac.caa<-((dn31.trinuc$cac*dn31.trinuc$caa)/(dn31.count.aa$h*dn31.count.aa$q))/dn31.count.diaa$hq
cps.cac.cac<-((dn31.trinuc$cac*dn31.trinuc$cac)/(dn31.count.aa$h*dn31.count.aa$h))/dn31.count.diaa$hh
cps.cac.cag<-((dn31.trinuc$cac*dn31.trinuc$cag)/(dn31.count.aa$h*dn31.count.aa$q))/dn31.count.diaa$hq
cps.cac.cat<-((dn31.trinuc$cac*dn31.trinuc$cat)/(dn31.count.aa$h*dn31.count.aa$h))/dn31.count.diaa$hh

cps.cac.cca<-((dn31.trinuc$cac*dn31.trinuc$cca)/(dn31.count.aa$h*dn31.count.aa$p))/dn31.count.diaa$hp
cps.cac.ccc<-((dn31.trinuc$cac*dn31.trinuc$ccc)/(dn31.count.aa$h*dn31.count.aa$p))/dn31.count.diaa$hp
cps.cac.ccg<-((dn31.trinuc$cac*dn31.trinuc$ccg)/(dn31.count.aa$h*dn31.count.aa$p))/dn31.count.diaa$hp
cps.cac.cct<-((dn31.trinuc$cac*dn31.trinuc$cct)/(dn31.count.aa$h*dn31.count.aa$p))/dn31.count.diaa$hp

cps.cac.cga<-((dn31.trinuc$cac*dn31.trinuc$cga)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cac.cgc<-((dn31.trinuc$cac*dn31.trinuc$cgc)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cac.cgg<-((dn31.trinuc$cac*dn31.trinuc$cgg)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cac.cgt<-((dn31.trinuc$cac*dn31.trinuc$cgt)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr

cps.cac.cta<-((dn31.trinuc$cac*dn31.trinuc$cta)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cac.ctc<-((dn31.trinuc$cac*dn31.trinuc$ctc)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cac.ctg<-((dn31.trinuc$cac*dn31.trinuc$ctg)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cac.ctt<-((dn31.trinuc$cac*dn31.trinuc$ctt)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl

cps.cac.gaa<-((dn31.trinuc$cac*dn31.trinuc$gaa)/(dn31.count.aa$h*dn31.count.aa$e))/dn31.count.diaa$he
cps.cac.gac<-((dn31.trinuc$cac*dn31.trinuc$gac)/(dn31.count.aa$h*dn31.count.aa$d))/dn31.count.diaa$hd
cps.cac.gag<-((dn31.trinuc$cac*dn31.trinuc$gag)/(dn31.count.aa$h*dn31.count.aa$e))/dn31.count.diaa$he
cps.cac.gat<-((dn31.trinuc$cac*dn31.trinuc$gat)/(dn31.count.aa$h*dn31.count.aa$d))/dn31.count.diaa$hd

cps.cac.gca<-((dn31.trinuc$cac*dn31.trinuc$gca)/(dn31.count.aa$h*dn31.count.aa$a))/dn31.count.diaa$ha
cps.cac.gcc<-((dn31.trinuc$cac*dn31.trinuc$gcc)/(dn31.count.aa$h*dn31.count.aa$a))/dn31.count.diaa$ha
cps.cac.gcg<-((dn31.trinuc$cac*dn31.trinuc$gcg)/(dn31.count.aa$h*dn31.count.aa$a))/dn31.count.diaa$ha
cps.cac.gct<-((dn31.trinuc$cac*dn31.trinuc$gct)/(dn31.count.aa$h*dn31.count.aa$a))/dn31.count.diaa$ha

cps.cac.gga<-((dn31.trinuc$cac*dn31.trinuc$gga)/(dn31.count.aa$h*dn31.count.aa$g))/dn31.count.diaa$hg
cps.cac.ggc<-((dn31.trinuc$cac*dn31.trinuc$ggc)/(dn31.count.aa$h*dn31.count.aa$g))/dn31.count.diaa$hg
cps.cac.ggg<-((dn31.trinuc$cac*dn31.trinuc$ggg)/(dn31.count.aa$h*dn31.count.aa$g))/dn31.count.diaa$hg
cps.cac.ggt<-((dn31.trinuc$cac*dn31.trinuc$ggt)/(dn31.count.aa$h*dn31.count.aa$g))/dn31.count.diaa$hg

cps.cac.gta<-((dn31.trinuc$cac*dn31.trinuc$gta)/(dn31.count.aa$h*dn31.count.aa$v))/dn31.count.diaa$hv
cps.cac.gtc<-((dn31.trinuc$cac*dn31.trinuc$gtc)/(dn31.count.aa$h*dn31.count.aa$v))/dn31.count.diaa$hv
cps.cac.gtg<-((dn31.trinuc$cac*dn31.trinuc$gtg)/(dn31.count.aa$h*dn31.count.aa$v))/dn31.count.diaa$hv
cps.cac.gtt<-((dn31.trinuc$cac*dn31.trinuc$gtt)/(dn31.count.aa$h*dn31.count.aa$v))/dn31.count.diaa$hv

#Stop codon
#cps.cac.taa<-((dn31.trinuc$cac*dn31.trinuc$taa)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cac.tac<-((dn31.trinuc$cac*dn31.trinuc$tac)/(dn31.count.aa$h*dn31.count.aa$y))/dn31.count.diaa$hy
#Stop codon
#cps.cac.tag<-((dn31.trinuc$cac*dn31.trinuc$tag)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cac.tat<-((dn31.trinuc$cac*dn31.trinuc$tat)/(dn31.count.aa$h*dn31.count.aa$y))/dn31.count.diaa$hy

cps.cac.tca<-((dn31.trinuc$cac*dn31.trinuc$tca)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs
cps.cac.tcc<-((dn31.trinuc$cac*dn31.trinuc$tcc)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs
cps.cac.tcg<-((dn31.trinuc$cac*dn31.trinuc$tcg)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs
cps.cac.tct<-((dn31.trinuc$cac*dn31.trinuc$tct)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs

#Stop codon
#cps.cac.tga<-((dn31.trinuc$cac*dn31.trinuc$tga)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cac.tgc<-((dn31.trinuc$cac*dn31.trinuc$tgc)/(dn31.count.aa$h*dn31.count.aa$c))/dn31.count.diaa$hc
cps.cac.tgg<-((dn31.trinuc$cac*dn31.trinuc$tgg)/(dn31.count.aa$h*dn31.count.aa$w))/dn31.count.diaa$hw
cps.cac.tgt<-((dn31.trinuc$cac*dn31.trinuc$tgt)/(dn31.count.aa$h*dn31.count.aa$c))/dn31.count.diaa$hc

cps.cac.tta<-((dn31.trinuc$cac*dn31.trinuc$tta)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cac.ttc<-((dn31.trinuc$cac*dn31.trinuc$ttc)/(dn31.count.aa$h*dn31.count.aa$f))/dn31.count.diaa$hf
cps.cac.ttg<-((dn31.trinuc$cac*dn31.trinuc$ttg)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cac.ttt<-((dn31.trinuc$cac*dn31.trinuc$ttt)/(dn31.count.aa$h*dn31.count.aa$f))/dn31.count.diaa$hf








cps.cag.aaa<-((dn31.trinuc$cag*dn31.trinuc$aaa)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.cag.aac<-((dn31.trinuc$cag*dn31.trinuc$aac)/(dn31.count.aa$q*dn31.count.aa$n))/dn31.count.diaa$qn
cps.cag.aag<-((dn31.trinuc$cag*dn31.trinuc$aag)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.cag.aat<-((dn31.trinuc$cag*dn31.trinuc$aat)/(dn31.count.aa$q*dn31.count.aa$n))/dn31.count.diaa$qn

cps.cag.aca<-((dn31.trinuc$cag*dn31.trinuc$aca)/(dn31.count.aa$q*dn31.count.aa$t))/dn31.count.diaa$qt
cps.cag.acc<-((dn31.trinuc$cag*dn31.trinuc$acc)/(dn31.count.aa$q*dn31.count.aa$t))/dn31.count.diaa$qt
cps.cag.acg<-((dn31.trinuc$cag*dn31.trinuc$acg)/(dn31.count.aa$q*dn31.count.aa$t))/dn31.count.diaa$qt
cps.cag.act<-((dn31.trinuc$cag*dn31.trinuc$act)/(dn31.count.aa$q*dn31.count.aa$t))/dn31.count.diaa$qt

cps.cag.aga<-((dn31.trinuc$cag*dn31.trinuc$aga)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.cag.agc<-((dn31.trinuc$cag*dn31.trinuc$agc)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs
cps.cag.agg<-((dn31.trinuc$cag*dn31.trinuc$agg)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.cag.agt<-((dn31.trinuc$cag*dn31.trinuc$agt)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs

cps.cag.ata<-((dn31.trinuc$cag*dn31.trinuc$ata)/(dn31.count.aa$q*dn31.count.aa$i))/dn31.count.diaa$qi
cps.cag.atc<-((dn31.trinuc$cag*dn31.trinuc$atc)/(dn31.count.aa$q*dn31.count.aa$i))/dn31.count.diaa$qi
cps.cag.atg<-((dn31.trinuc$cag*dn31.trinuc$atg)/(dn31.count.aa$q*dn31.count.aa$m))/dn31.count.diaa$qm
cps.cag.att<-((dn31.trinuc$cag*dn31.trinuc$att)/(dn31.count.aa$q*dn31.count.aa$i))/dn31.count.diaa$qi

cps.cag.caa<-((dn31.trinuc$cag*dn31.trinuc$caa)/(dn31.count.aa$q*dn31.count.aa$q))/dn31.count.diaa$qq
cps.cag.cac<-((dn31.trinuc$cag*dn31.trinuc$cac)/(dn31.count.aa$q*dn31.count.aa$h))/dn31.count.diaa$qh
cps.cag.cag<-((dn31.trinuc$cag*dn31.trinuc$cag)/(dn31.count.aa$q*dn31.count.aa$q))/dn31.count.diaa$qq
cps.cag.cat<-((dn31.trinuc$cag*dn31.trinuc$cat)/(dn31.count.aa$q*dn31.count.aa$h))/dn31.count.diaa$qh

cps.cag.cca<-((dn31.trinuc$cag*dn31.trinuc$cca)/(dn31.count.aa$q*dn31.count.aa$p))/dn31.count.diaa$qp
cps.cag.ccc<-((dn31.trinuc$cag*dn31.trinuc$ccc)/(dn31.count.aa$q*dn31.count.aa$p))/dn31.count.diaa$qp
cps.cag.ccg<-((dn31.trinuc$cag*dn31.trinuc$ccg)/(dn31.count.aa$q*dn31.count.aa$p))/dn31.count.diaa$qp
cps.cag.cct<-((dn31.trinuc$cag*dn31.trinuc$cct)/(dn31.count.aa$q*dn31.count.aa$p))/dn31.count.diaa$qp

cps.cag.cga<-((dn31.trinuc$cag*dn31.trinuc$cga)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.cag.cgc<-((dn31.trinuc$cag*dn31.trinuc$cgc)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.cag.cgg<-((dn31.trinuc$cag*dn31.trinuc$cgg)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr
cps.cag.cgt<-((dn31.trinuc$cag*dn31.trinuc$cgt)/(dn31.count.aa$q*dn31.count.aa$r))/dn31.count.diaa$qr

cps.cag.cta<-((dn31.trinuc$cag*dn31.trinuc$cta)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.cag.ctc<-((dn31.trinuc$cag*dn31.trinuc$ctc)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.cag.ctg<-((dn31.trinuc$cag*dn31.trinuc$ctg)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.cag.ctt<-((dn31.trinuc$cag*dn31.trinuc$ctt)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql

cps.cag.gaa<-((dn31.trinuc$cag*dn31.trinuc$gaa)/(dn31.count.aa$q*dn31.count.aa$e))/dn31.count.diaa$qe
cps.cag.gac<-((dn31.trinuc$cag*dn31.trinuc$gac)/(dn31.count.aa$q*dn31.count.aa$d))/dn31.count.diaa$qd
cps.cag.gag<-((dn31.trinuc$cag*dn31.trinuc$gag)/(dn31.count.aa$q*dn31.count.aa$e))/dn31.count.diaa$qe
cps.cag.gat<-((dn31.trinuc$cag*dn31.trinuc$gat)/(dn31.count.aa$q*dn31.count.aa$d))/dn31.count.diaa$qd

cps.cag.gca<-((dn31.trinuc$cag*dn31.trinuc$gca)/(dn31.count.aa$q*dn31.count.aa$a))/dn31.count.diaa$qa
cps.cag.gcc<-((dn31.trinuc$cag*dn31.trinuc$gcc)/(dn31.count.aa$q*dn31.count.aa$a))/dn31.count.diaa$qa
cps.cag.gcg<-((dn31.trinuc$cag*dn31.trinuc$gcg)/(dn31.count.aa$q*dn31.count.aa$a))/dn31.count.diaa$qa
cps.cag.gct<-((dn31.trinuc$cag*dn31.trinuc$gct)/(dn31.count.aa$q*dn31.count.aa$a))/dn31.count.diaa$qa

cps.cag.gga<-((dn31.trinuc$cag*dn31.trinuc$gga)/(dn31.count.aa$q*dn31.count.aa$g))/dn31.count.diaa$qg
cps.cag.ggc<-((dn31.trinuc$cag*dn31.trinuc$ggc)/(dn31.count.aa$q*dn31.count.aa$g))/dn31.count.diaa$qg
cps.cag.ggg<-((dn31.trinuc$cag*dn31.trinuc$ggg)/(dn31.count.aa$q*dn31.count.aa$g))/dn31.count.diaa$qg
cps.cag.ggt<-((dn31.trinuc$cag*dn31.trinuc$ggt)/(dn31.count.aa$q*dn31.count.aa$g))/dn31.count.diaa$qg

cps.cag.gta<-((dn31.trinuc$cag*dn31.trinuc$gta)/(dn31.count.aa$q*dn31.count.aa$v))/dn31.count.diaa$qv
cps.cag.gtc<-((dn31.trinuc$cag*dn31.trinuc$gtc)/(dn31.count.aa$q*dn31.count.aa$v))/dn31.count.diaa$qv
cps.cag.gtg<-((dn31.trinuc$cag*dn31.trinuc$gtg)/(dn31.count.aa$q*dn31.count.aa$v))/dn31.count.diaa$qv
cps.cag.gtt<-((dn31.trinuc$cag*dn31.trinuc$gtt)/(dn31.count.aa$q*dn31.count.aa$v))/dn31.count.diaa$qv

#Stop codon
#cps.cag.taa<-((dn31.trinuc$cag*dn31.trinuc$taa)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.cag.tac<-((dn31.trinuc$cag*dn31.trinuc$tac)/(dn31.count.aa$q*dn31.count.aa$y))/dn31.count.diaa$qy
#Stop codon
#cps.cag.tag<-((dn31.trinuc$cag*dn31.trinuc$tag)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.cag.tat<-((dn31.trinuc$cag*dn31.trinuc$tat)/(dn31.count.aa$q*dn31.count.aa$y))/dn31.count.diaa$qy

cps.cag.tca<-((dn31.trinuc$cag*dn31.trinuc$tca)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs
cps.cag.tcc<-((dn31.trinuc$cag*dn31.trinuc$tcc)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs
cps.cag.tcg<-((dn31.trinuc$cag*dn31.trinuc$tcg)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs
cps.cag.tct<-((dn31.trinuc$cag*dn31.trinuc$tct)/(dn31.count.aa$q*dn31.count.aa$s))/dn31.count.diaa$qs

#Stop codon
#cps.cag.tga<-((dn31.trinuc$cag*dn31.trinuc$tga)/(dn31.count.aa$q*dn31.count.aa$k))/dn31.count.diaa$qk
cps.cag.tgc<-((dn31.trinuc$cag*dn31.trinuc$tgc)/(dn31.count.aa$q*dn31.count.aa$c))/dn31.count.diaa$qc
cps.cag.tgg<-((dn31.trinuc$cag*dn31.trinuc$tgg)/(dn31.count.aa$q*dn31.count.aa$w))/dn31.count.diaa$qw
cps.cag.tgt<-((dn31.trinuc$cag*dn31.trinuc$tgt)/(dn31.count.aa$q*dn31.count.aa$c))/dn31.count.diaa$qc

cps.cag.tta<-((dn31.trinuc$cag*dn31.trinuc$tta)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.cag.ttc<-((dn31.trinuc$cag*dn31.trinuc$ttc)/(dn31.count.aa$q*dn31.count.aa$f))/dn31.count.diaa$qf
cps.cag.ttg<-((dn31.trinuc$cag*dn31.trinuc$ttg)/(dn31.count.aa$q*dn31.count.aa$l))/dn31.count.diaa$ql
cps.cag.ttt<-((dn31.trinuc$cag*dn31.trinuc$ttt)/(dn31.count.aa$q*dn31.count.aa$f))/dn31.count.diaa$qf








cps.cat.aaa<-((dn31.trinuc$cat*dn31.trinuc$aaa)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cat.aac<-((dn31.trinuc$cat*dn31.trinuc$aac)/(dn31.count.aa$h*dn31.count.aa$n))/dn31.count.diaa$hn
cps.cat.aag<-((dn31.trinuc$cat*dn31.trinuc$aag)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cat.aat<-((dn31.trinuc$cat*dn31.trinuc$aat)/(dn31.count.aa$h*dn31.count.aa$n))/dn31.count.diaa$hn

cps.cat.aca<-((dn31.trinuc$cat*dn31.trinuc$aca)/(dn31.count.aa$h*dn31.count.aa$t))/dn31.count.diaa$ht
cps.cat.acc<-((dn31.trinuc$cat*dn31.trinuc$acc)/(dn31.count.aa$h*dn31.count.aa$t))/dn31.count.diaa$ht
cps.cat.acg<-((dn31.trinuc$cat*dn31.trinuc$acg)/(dn31.count.aa$h*dn31.count.aa$t))/dn31.count.diaa$ht
cps.cat.act<-((dn31.trinuc$cat*dn31.trinuc$act)/(dn31.count.aa$h*dn31.count.aa$t))/dn31.count.diaa$ht

cps.cat.aga<-((dn31.trinuc$cat*dn31.trinuc$aga)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cat.agc<-((dn31.trinuc$cat*dn31.trinuc$agc)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs
cps.cat.agg<-((dn31.trinuc$cat*dn31.trinuc$agg)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cat.agt<-((dn31.trinuc$cat*dn31.trinuc$agt)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs

cps.cat.ata<-((dn31.trinuc$cat*dn31.trinuc$ata)/(dn31.count.aa$h*dn31.count.aa$i))/dn31.count.diaa$hi
cps.cat.atc<-((dn31.trinuc$cat*dn31.trinuc$atc)/(dn31.count.aa$h*dn31.count.aa$i))/dn31.count.diaa$hi
cps.cat.atg<-((dn31.trinuc$cat*dn31.trinuc$atg)/(dn31.count.aa$h*dn31.count.aa$m))/dn31.count.diaa$hm
cps.cat.att<-((dn31.trinuc$cat*dn31.trinuc$att)/(dn31.count.aa$h*dn31.count.aa$i))/dn31.count.diaa$hi

cps.cat.caa<-((dn31.trinuc$cat*dn31.trinuc$caa)/(dn31.count.aa$h*dn31.count.aa$q))/dn31.count.diaa$hq
cps.cat.cac<-((dn31.trinuc$cat*dn31.trinuc$cac)/(dn31.count.aa$h*dn31.count.aa$h))/dn31.count.diaa$hh
cps.cat.cag<-((dn31.trinuc$cat*dn31.trinuc$cag)/(dn31.count.aa$h*dn31.count.aa$q))/dn31.count.diaa$hq
cps.cat.cat<-((dn31.trinuc$cat*dn31.trinuc$cat)/(dn31.count.aa$h*dn31.count.aa$h))/dn31.count.diaa$hh

cps.cat.cca<-((dn31.trinuc$cat*dn31.trinuc$cca)/(dn31.count.aa$h*dn31.count.aa$p))/dn31.count.diaa$hp
cps.cat.ccc<-((dn31.trinuc$cat*dn31.trinuc$ccc)/(dn31.count.aa$h*dn31.count.aa$p))/dn31.count.diaa$hp
cps.cat.ccg<-((dn31.trinuc$cat*dn31.trinuc$ccg)/(dn31.count.aa$h*dn31.count.aa$p))/dn31.count.diaa$hp
cps.cat.cct<-((dn31.trinuc$cat*dn31.trinuc$cct)/(dn31.count.aa$h*dn31.count.aa$p))/dn31.count.diaa$hp

cps.cat.cga<-((dn31.trinuc$cat*dn31.trinuc$cga)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cat.cgc<-((dn31.trinuc$cat*dn31.trinuc$cgc)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cat.cgg<-((dn31.trinuc$cat*dn31.trinuc$cgg)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr
cps.cat.cgt<-((dn31.trinuc$cat*dn31.trinuc$cgt)/(dn31.count.aa$h*dn31.count.aa$r))/dn31.count.diaa$hr

cps.cat.cta<-((dn31.trinuc$cat*dn31.trinuc$cta)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cat.ctc<-((dn31.trinuc$cat*dn31.trinuc$ctc)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cat.ctg<-((dn31.trinuc$cat*dn31.trinuc$ctg)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cat.ctt<-((dn31.trinuc$cat*dn31.trinuc$ctt)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl

cps.cat.gaa<-((dn31.trinuc$cat*dn31.trinuc$gaa)/(dn31.count.aa$h*dn31.count.aa$e))/dn31.count.diaa$he
cps.cat.gac<-((dn31.trinuc$cat*dn31.trinuc$gac)/(dn31.count.aa$h*dn31.count.aa$d))/dn31.count.diaa$hd
cps.cat.gag<-((dn31.trinuc$cat*dn31.trinuc$gag)/(dn31.count.aa$h*dn31.count.aa$e))/dn31.count.diaa$he
cps.cat.gat<-((dn31.trinuc$cat*dn31.trinuc$gat)/(dn31.count.aa$h*dn31.count.aa$d))/dn31.count.diaa$hd

cps.cat.gca<-((dn31.trinuc$cat*dn31.trinuc$gca)/(dn31.count.aa$h*dn31.count.aa$a))/dn31.count.diaa$ha
cps.cat.gcc<-((dn31.trinuc$cat*dn31.trinuc$gcc)/(dn31.count.aa$h*dn31.count.aa$a))/dn31.count.diaa$ha
cps.cat.gcg<-((dn31.trinuc$cat*dn31.trinuc$gcg)/(dn31.count.aa$h*dn31.count.aa$a))/dn31.count.diaa$ha
cps.cat.gct<-((dn31.trinuc$cat*dn31.trinuc$gct)/(dn31.count.aa$h*dn31.count.aa$a))/dn31.count.diaa$ha

cps.cat.gga<-((dn31.trinuc$cat*dn31.trinuc$gga)/(dn31.count.aa$h*dn31.count.aa$g))/dn31.count.diaa$hg
cps.cat.ggc<-((dn31.trinuc$cat*dn31.trinuc$ggc)/(dn31.count.aa$h*dn31.count.aa$g))/dn31.count.diaa$hg
cps.cat.ggg<-((dn31.trinuc$cat*dn31.trinuc$ggg)/(dn31.count.aa$h*dn31.count.aa$g))/dn31.count.diaa$hg
cps.cat.ggt<-((dn31.trinuc$cat*dn31.trinuc$ggt)/(dn31.count.aa$h*dn31.count.aa$g))/dn31.count.diaa$hg

cps.cat.gta<-((dn31.trinuc$cat*dn31.trinuc$gta)/(dn31.count.aa$h*dn31.count.aa$v))/dn31.count.diaa$hv
cps.cat.gtc<-((dn31.trinuc$cat*dn31.trinuc$gtc)/(dn31.count.aa$h*dn31.count.aa$v))/dn31.count.diaa$hv
cps.cat.gtg<-((dn31.trinuc$cat*dn31.trinuc$gtg)/(dn31.count.aa$h*dn31.count.aa$v))/dn31.count.diaa$hv
cps.cat.gtt<-((dn31.trinuc$cat*dn31.trinuc$gtt)/(dn31.count.aa$h*dn31.count.aa$v))/dn31.count.diaa$hv

#Stop codon
#cps.cat.taa<-((dn31.trinuc$cat*dn31.trinuc$taa)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cat.tac<-((dn31.trinuc$cat*dn31.trinuc$tac)/(dn31.count.aa$h*dn31.count.aa$y))/dn31.count.diaa$hy
#Stop codon
#cps.cat.tag<-((dn31.trinuc$cat*dn31.trinuc$tag)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cat.tat<-((dn31.trinuc$cat*dn31.trinuc$tat)/(dn31.count.aa$h*dn31.count.aa$y))/dn31.count.diaa$hy

cps.cat.tca<-((dn31.trinuc$cat*dn31.trinuc$tca)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs
cps.cat.tcc<-((dn31.trinuc$cat*dn31.trinuc$tcc)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs
cps.cat.tcg<-((dn31.trinuc$cat*dn31.trinuc$tcg)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs
cps.cat.tct<-((dn31.trinuc$cat*dn31.trinuc$tct)/(dn31.count.aa$h*dn31.count.aa$s))/dn31.count.diaa$hs

#Stop codon
#cps.cat.tga<-((dn31.trinuc$cat*dn31.trinuc$tga)/(dn31.count.aa$h*dn31.count.aa$k))/dn31.count.diaa$hk
cps.cat.tgc<-((dn31.trinuc$cat*dn31.trinuc$tgc)/(dn31.count.aa$h*dn31.count.aa$c))/dn31.count.diaa$hc
cps.cat.tgg<-((dn31.trinuc$cat*dn31.trinuc$tgg)/(dn31.count.aa$h*dn31.count.aa$w))/dn31.count.diaa$hw
cps.cat.tgt<-((dn31.trinuc$cat*dn31.trinuc$tgt)/(dn31.count.aa$h*dn31.count.aa$c))/dn31.count.diaa$hc

cps.cat.tta<-((dn31.trinuc$cat*dn31.trinuc$tta)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cat.ttc<-((dn31.trinuc$cat*dn31.trinuc$ttc)/(dn31.count.aa$h*dn31.count.aa$f))/dn31.count.diaa$hf
cps.cat.ttg<-((dn31.trinuc$cat*dn31.trinuc$ttg)/(dn31.count.aa$h*dn31.count.aa$l))/dn31.count.diaa$hl
cps.cat.ttt<-((dn31.trinuc$cat*dn31.trinuc$ttt)/(dn31.count.aa$h*dn31.count.aa$f))/dn31.count.diaa$hf








cps.cca.aaa<-((dn31.trinuc$cca*dn31.trinuc$aaa)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cca.aac<-((dn31.trinuc$cca*dn31.trinuc$aac)/(dn31.count.aa$p*dn31.count.aa$n))/dn31.count.diaa$pn
cps.cca.aag<-((dn31.trinuc$cca*dn31.trinuc$aag)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cca.aat<-((dn31.trinuc$cca*dn31.trinuc$aat)/(dn31.count.aa$p*dn31.count.aa$n))/dn31.count.diaa$pn

cps.cca.aca<-((dn31.trinuc$cca*dn31.trinuc$aca)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.cca.acc<-((dn31.trinuc$cca*dn31.trinuc$acc)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.cca.acg<-((dn31.trinuc$cca*dn31.trinuc$acg)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.cca.act<-((dn31.trinuc$cca*dn31.trinuc$act)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt

cps.cca.aga<-((dn31.trinuc$cca*dn31.trinuc$aga)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cca.agc<-((dn31.trinuc$cca*dn31.trinuc$agc)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.cca.agg<-((dn31.trinuc$cca*dn31.trinuc$agg)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cca.agt<-((dn31.trinuc$cca*dn31.trinuc$agt)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps

cps.cca.ata<-((dn31.trinuc$cca*dn31.trinuc$ata)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi
cps.cca.atc<-((dn31.trinuc$cca*dn31.trinuc$atc)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi
cps.cca.atg<-((dn31.trinuc$cca*dn31.trinuc$atg)/(dn31.count.aa$p*dn31.count.aa$m))/dn31.count.diaa$pm
cps.cca.att<-((dn31.trinuc$cca*dn31.trinuc$att)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi

cps.cca.caa<-((dn31.trinuc$cca*dn31.trinuc$caa)/(dn31.count.aa$p*dn31.count.aa$q))/dn31.count.diaa$pq
cps.cca.cac<-((dn31.trinuc$cca*dn31.trinuc$cac)/(dn31.count.aa$p*dn31.count.aa$h))/dn31.count.diaa$ph
cps.cca.cag<-((dn31.trinuc$cca*dn31.trinuc$cag)/(dn31.count.aa$p*dn31.count.aa$q))/dn31.count.diaa$pq
cps.cca.cat<-((dn31.trinuc$cca*dn31.trinuc$cat)/(dn31.count.aa$p*dn31.count.aa$h))/dn31.count.diaa$ph

cps.cca.cca<-((dn31.trinuc$cca*dn31.trinuc$cca)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.cca.ccc<-((dn31.trinuc$cca*dn31.trinuc$ccc)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.cca.ccg<-((dn31.trinuc$cca*dn31.trinuc$ccg)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.cca.cct<-((dn31.trinuc$cca*dn31.trinuc$cct)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp

cps.cca.cga<-((dn31.trinuc$cca*dn31.trinuc$cga)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cca.cgc<-((dn31.trinuc$cca*dn31.trinuc$cgc)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cca.cgg<-((dn31.trinuc$cca*dn31.trinuc$cgg)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cca.cgt<-((dn31.trinuc$cca*dn31.trinuc$cgt)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr

cps.cca.cta<-((dn31.trinuc$cca*dn31.trinuc$cta)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cca.ctc<-((dn31.trinuc$cca*dn31.trinuc$ctc)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cca.ctg<-((dn31.trinuc$cca*dn31.trinuc$ctg)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cca.ctt<-((dn31.trinuc$cca*dn31.trinuc$ctt)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl

cps.cca.gaa<-((dn31.trinuc$cca*dn31.trinuc$gaa)/(dn31.count.aa$p*dn31.count.aa$e))/dn31.count.diaa$pe
cps.cca.gac<-((dn31.trinuc$cca*dn31.trinuc$gac)/(dn31.count.aa$p*dn31.count.aa$d))/dn31.count.diaa$pd
cps.cca.gag<-((dn31.trinuc$cca*dn31.trinuc$gag)/(dn31.count.aa$p*dn31.count.aa$e))/dn31.count.diaa$pe
cps.cca.gat<-((dn31.trinuc$cca*dn31.trinuc$gat)/(dn31.count.aa$p*dn31.count.aa$d))/dn31.count.diaa$pd

cps.cca.gca<-((dn31.trinuc$cca*dn31.trinuc$gca)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.cca.gcc<-((dn31.trinuc$cca*dn31.trinuc$gcc)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.cca.gcg<-((dn31.trinuc$cca*dn31.trinuc$gcg)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.cca.gct<-((dn31.trinuc$cca*dn31.trinuc$gct)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa

cps.cca.gga<-((dn31.trinuc$cca*dn31.trinuc$gga)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.cca.ggc<-((dn31.trinuc$cca*dn31.trinuc$ggc)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.cca.ggg<-((dn31.trinuc$cca*dn31.trinuc$ggg)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.cca.ggt<-((dn31.trinuc$cca*dn31.trinuc$ggt)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg

cps.cca.gta<-((dn31.trinuc$cca*dn31.trinuc$gta)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.cca.gtc<-((dn31.trinuc$cca*dn31.trinuc$gtc)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.cca.gtg<-((dn31.trinuc$cca*dn31.trinuc$gtg)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.cca.gtt<-((dn31.trinuc$cca*dn31.trinuc$gtt)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv

#Stop codon
#cps.cca.taa<-((dn31.trinuc$cca*dn31.trinuc$taa)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cca.tac<-((dn31.trinuc$cca*dn31.trinuc$tac)/(dn31.count.aa$p*dn31.count.aa$y))/dn31.count.diaa$py
#Stop codon
#cps.cca.tag<-((dn31.trinuc$cca*dn31.trinuc$tag)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cca.tat<-((dn31.trinuc$cca*dn31.trinuc$tat)/(dn31.count.aa$p*dn31.count.aa$y))/dn31.count.diaa$py

cps.cca.tca<-((dn31.trinuc$cca*dn31.trinuc$tca)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.cca.tcc<-((dn31.trinuc$cca*dn31.trinuc$tcc)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.cca.tcg<-((dn31.trinuc$cca*dn31.trinuc$tcg)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.cca.tct<-((dn31.trinuc$cca*dn31.trinuc$tct)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps

#Stop codon
#cps.cca.tga<-((dn31.trinuc$cca*dn31.trinuc$tga)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cca.tgc<-((dn31.trinuc$cca*dn31.trinuc$tgc)/(dn31.count.aa$p*dn31.count.aa$c))/dn31.count.diaa$pc
cps.cca.tgg<-((dn31.trinuc$cca*dn31.trinuc$tgg)/(dn31.count.aa$p*dn31.count.aa$w))/dn31.count.diaa$pw
cps.cca.tgt<-((dn31.trinuc$cca*dn31.trinuc$tgt)/(dn31.count.aa$p*dn31.count.aa$c))/dn31.count.diaa$pc

cps.cca.tta<-((dn31.trinuc$cca*dn31.trinuc$tta)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cca.ttc<-((dn31.trinuc$cca*dn31.trinuc$ttc)/(dn31.count.aa$p*dn31.count.aa$f))/dn31.count.diaa$pf
cps.cca.ttg<-((dn31.trinuc$cca*dn31.trinuc$ttg)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cca.ttt<-((dn31.trinuc$cca*dn31.trinuc$ttt)/(dn31.count.aa$p*dn31.count.aa$f))/dn31.count.diaa$pf








cps.ccc.aaa<-((dn31.trinuc$ccc*dn31.trinuc$aaa)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccc.aac<-((dn31.trinuc$ccc*dn31.trinuc$aac)/(dn31.count.aa$p*dn31.count.aa$n))/dn31.count.diaa$pn
cps.ccc.aag<-((dn31.trinuc$ccc*dn31.trinuc$aag)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccc.aat<-((dn31.trinuc$ccc*dn31.trinuc$aat)/(dn31.count.aa$p*dn31.count.aa$n))/dn31.count.diaa$pn

cps.ccc.aca<-((dn31.trinuc$ccc*dn31.trinuc$aca)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.ccc.acc<-((dn31.trinuc$ccc*dn31.trinuc$acc)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.ccc.acg<-((dn31.trinuc$ccc*dn31.trinuc$acg)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.ccc.act<-((dn31.trinuc$ccc*dn31.trinuc$act)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt

cps.ccc.aga<-((dn31.trinuc$ccc*dn31.trinuc$aga)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccc.agc<-((dn31.trinuc$ccc*dn31.trinuc$agc)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.ccc.agg<-((dn31.trinuc$ccc*dn31.trinuc$agg)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccc.agt<-((dn31.trinuc$ccc*dn31.trinuc$agt)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps

cps.ccc.ata<-((dn31.trinuc$ccc*dn31.trinuc$ata)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi
cps.ccc.atc<-((dn31.trinuc$ccc*dn31.trinuc$atc)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi
cps.ccc.atg<-((dn31.trinuc$ccc*dn31.trinuc$atg)/(dn31.count.aa$p*dn31.count.aa$m))/dn31.count.diaa$pm
cps.ccc.att<-((dn31.trinuc$ccc*dn31.trinuc$att)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi

cps.ccc.caa<-((dn31.trinuc$ccc*dn31.trinuc$caa)/(dn31.count.aa$p*dn31.count.aa$q))/dn31.count.diaa$pq
cps.ccc.cac<-((dn31.trinuc$ccc*dn31.trinuc$cac)/(dn31.count.aa$p*dn31.count.aa$h))/dn31.count.diaa$ph
cps.ccc.cag<-((dn31.trinuc$ccc*dn31.trinuc$cag)/(dn31.count.aa$p*dn31.count.aa$q))/dn31.count.diaa$pq
cps.ccc.cat<-((dn31.trinuc$ccc*dn31.trinuc$cat)/(dn31.count.aa$p*dn31.count.aa$h))/dn31.count.diaa$ph

cps.ccc.cca<-((dn31.trinuc$ccc*dn31.trinuc$cca)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.ccc.ccc<-((dn31.trinuc$ccc*dn31.trinuc$ccc)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.ccc.ccg<-((dn31.trinuc$ccc*dn31.trinuc$ccg)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.ccc.cct<-((dn31.trinuc$ccc*dn31.trinuc$cct)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp

cps.ccc.cga<-((dn31.trinuc$ccc*dn31.trinuc$cga)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccc.cgc<-((dn31.trinuc$ccc*dn31.trinuc$cgc)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccc.cgg<-((dn31.trinuc$ccc*dn31.trinuc$cgg)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccc.cgt<-((dn31.trinuc$ccc*dn31.trinuc$cgt)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr

cps.ccc.cta<-((dn31.trinuc$ccc*dn31.trinuc$cta)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccc.ctc<-((dn31.trinuc$ccc*dn31.trinuc$ctc)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccc.ctg<-((dn31.trinuc$ccc*dn31.trinuc$ctg)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccc.ctt<-((dn31.trinuc$ccc*dn31.trinuc$ctt)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl

cps.ccc.gaa<-((dn31.trinuc$ccc*dn31.trinuc$gaa)/(dn31.count.aa$p*dn31.count.aa$e))/dn31.count.diaa$pe
cps.ccc.gac<-((dn31.trinuc$ccc*dn31.trinuc$gac)/(dn31.count.aa$p*dn31.count.aa$d))/dn31.count.diaa$pd
cps.ccc.gag<-((dn31.trinuc$ccc*dn31.trinuc$gag)/(dn31.count.aa$p*dn31.count.aa$e))/dn31.count.diaa$pe
cps.ccc.gat<-((dn31.trinuc$ccc*dn31.trinuc$gat)/(dn31.count.aa$p*dn31.count.aa$d))/dn31.count.diaa$pd

cps.ccc.gca<-((dn31.trinuc$ccc*dn31.trinuc$gca)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.ccc.gcc<-((dn31.trinuc$ccc*dn31.trinuc$gcc)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.ccc.gcg<-((dn31.trinuc$ccc*dn31.trinuc$gcg)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.ccc.gct<-((dn31.trinuc$ccc*dn31.trinuc$gct)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa

cps.ccc.gga<-((dn31.trinuc$ccc*dn31.trinuc$gga)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.ccc.ggc<-((dn31.trinuc$ccc*dn31.trinuc$ggc)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.ccc.ggg<-((dn31.trinuc$ccc*dn31.trinuc$ggg)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.ccc.ggt<-((dn31.trinuc$ccc*dn31.trinuc$ggt)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg

cps.ccc.gta<-((dn31.trinuc$ccc*dn31.trinuc$gta)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.ccc.gtc<-((dn31.trinuc$ccc*dn31.trinuc$gtc)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.ccc.gtg<-((dn31.trinuc$ccc*dn31.trinuc$gtg)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.ccc.gtt<-((dn31.trinuc$ccc*dn31.trinuc$gtt)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv

#Stop codon
#cps.ccc.taa<-((dn31.trinuc$ccc*dn31.trinuc$taa)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccc.tac<-((dn31.trinuc$ccc*dn31.trinuc$tac)/(dn31.count.aa$p*dn31.count.aa$y))/dn31.count.diaa$py
#Stop codon
#cps.ccc.tag<-((dn31.trinuc$ccc*dn31.trinuc$tag)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccc.tat<-((dn31.trinuc$ccc*dn31.trinuc$tat)/(dn31.count.aa$p*dn31.count.aa$y))/dn31.count.diaa$py

cps.ccc.tca<-((dn31.trinuc$ccc*dn31.trinuc$tca)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.ccc.tcc<-((dn31.trinuc$ccc*dn31.trinuc$tcc)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.ccc.tcg<-((dn31.trinuc$ccc*dn31.trinuc$tcg)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.ccc.tct<-((dn31.trinuc$ccc*dn31.trinuc$tct)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps

#Stop codon
#cps.ccc.tga<-((dn31.trinuc$ccc*dn31.trinuc$tga)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccc.tgc<-((dn31.trinuc$ccc*dn31.trinuc$tgc)/(dn31.count.aa$p*dn31.count.aa$c))/dn31.count.diaa$pc
cps.ccc.tgg<-((dn31.trinuc$ccc*dn31.trinuc$tgg)/(dn31.count.aa$p*dn31.count.aa$w))/dn31.count.diaa$pw
cps.ccc.tgt<-((dn31.trinuc$ccc*dn31.trinuc$tgt)/(dn31.count.aa$p*dn31.count.aa$c))/dn31.count.diaa$pc

cps.ccc.tta<-((dn31.trinuc$ccc*dn31.trinuc$tta)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccc.ttc<-((dn31.trinuc$ccc*dn31.trinuc$ttc)/(dn31.count.aa$p*dn31.count.aa$f))/dn31.count.diaa$pf
cps.ccc.ttg<-((dn31.trinuc$ccc*dn31.trinuc$ttg)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccc.ttt<-((dn31.trinuc$ccc*dn31.trinuc$ttt)/(dn31.count.aa$p*dn31.count.aa$f))/dn31.count.diaa$pf








cps.ccg.aaa<-((dn31.trinuc$ccg*dn31.trinuc$aaa)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccg.aac<-((dn31.trinuc$ccg*dn31.trinuc$aac)/(dn31.count.aa$p*dn31.count.aa$n))/dn31.count.diaa$pn
cps.ccg.aag<-((dn31.trinuc$ccg*dn31.trinuc$aag)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccg.aat<-((dn31.trinuc$ccg*dn31.trinuc$aat)/(dn31.count.aa$p*dn31.count.aa$n))/dn31.count.diaa$pn

cps.ccg.aca<-((dn31.trinuc$ccg*dn31.trinuc$aca)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.ccg.acc<-((dn31.trinuc$ccg*dn31.trinuc$acc)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.ccg.acg<-((dn31.trinuc$ccg*dn31.trinuc$acg)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.ccg.act<-((dn31.trinuc$ccg*dn31.trinuc$act)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt

cps.ccg.aga<-((dn31.trinuc$ccg*dn31.trinuc$aga)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccg.agc<-((dn31.trinuc$ccg*dn31.trinuc$agc)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.ccg.agg<-((dn31.trinuc$ccg*dn31.trinuc$agg)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccg.agt<-((dn31.trinuc$ccg*dn31.trinuc$agt)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps

cps.ccg.ata<-((dn31.trinuc$ccg*dn31.trinuc$ata)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi
cps.ccg.atc<-((dn31.trinuc$ccg*dn31.trinuc$atc)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi
cps.ccg.atg<-((dn31.trinuc$ccg*dn31.trinuc$atg)/(dn31.count.aa$p*dn31.count.aa$m))/dn31.count.diaa$pm
cps.ccg.att<-((dn31.trinuc$ccg*dn31.trinuc$att)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi

cps.ccg.caa<-((dn31.trinuc$ccg*dn31.trinuc$caa)/(dn31.count.aa$p*dn31.count.aa$q))/dn31.count.diaa$pq
cps.ccg.cac<-((dn31.trinuc$ccg*dn31.trinuc$cac)/(dn31.count.aa$p*dn31.count.aa$h))/dn31.count.diaa$ph
cps.ccg.cag<-((dn31.trinuc$ccg*dn31.trinuc$cag)/(dn31.count.aa$p*dn31.count.aa$q))/dn31.count.diaa$pq
cps.ccg.cat<-((dn31.trinuc$ccg*dn31.trinuc$cat)/(dn31.count.aa$p*dn31.count.aa$h))/dn31.count.diaa$ph

cps.ccg.cca<-((dn31.trinuc$ccg*dn31.trinuc$cca)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.ccg.ccc<-((dn31.trinuc$ccg*dn31.trinuc$ccc)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.ccg.ccg<-((dn31.trinuc$ccg*dn31.trinuc$ccg)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.ccg.cct<-((dn31.trinuc$ccg*dn31.trinuc$cct)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp

cps.ccg.cga<-((dn31.trinuc$ccg*dn31.trinuc$cga)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccg.cgc<-((dn31.trinuc$ccg*dn31.trinuc$cgc)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccg.cgg<-((dn31.trinuc$ccg*dn31.trinuc$cgg)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.ccg.cgt<-((dn31.trinuc$ccg*dn31.trinuc$cgt)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr

cps.ccg.cta<-((dn31.trinuc$ccg*dn31.trinuc$cta)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccg.ctc<-((dn31.trinuc$ccg*dn31.trinuc$ctc)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccg.ctg<-((dn31.trinuc$ccg*dn31.trinuc$ctg)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccg.ctt<-((dn31.trinuc$ccg*dn31.trinuc$ctt)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl

cps.ccg.gaa<-((dn31.trinuc$ccg*dn31.trinuc$gaa)/(dn31.count.aa$p*dn31.count.aa$e))/dn31.count.diaa$pe
cps.ccg.gac<-((dn31.trinuc$ccg*dn31.trinuc$gac)/(dn31.count.aa$p*dn31.count.aa$d))/dn31.count.diaa$pd
cps.ccg.gag<-((dn31.trinuc$ccg*dn31.trinuc$gag)/(dn31.count.aa$p*dn31.count.aa$e))/dn31.count.diaa$pe
cps.ccg.gat<-((dn31.trinuc$ccg*dn31.trinuc$gat)/(dn31.count.aa$p*dn31.count.aa$d))/dn31.count.diaa$pd

cps.ccg.gca<-((dn31.trinuc$ccg*dn31.trinuc$gca)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.ccg.gcc<-((dn31.trinuc$ccg*dn31.trinuc$gcc)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.ccg.gcg<-((dn31.trinuc$ccg*dn31.trinuc$gcg)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.ccg.gct<-((dn31.trinuc$ccg*dn31.trinuc$gct)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa

cps.ccg.gga<-((dn31.trinuc$ccg*dn31.trinuc$gga)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.ccg.ggc<-((dn31.trinuc$ccg*dn31.trinuc$ggc)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.ccg.ggg<-((dn31.trinuc$ccg*dn31.trinuc$ggg)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.ccg.ggt<-((dn31.trinuc$ccg*dn31.trinuc$ggt)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg

cps.ccg.gta<-((dn31.trinuc$ccg*dn31.trinuc$gta)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.ccg.gtc<-((dn31.trinuc$ccg*dn31.trinuc$gtc)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.ccg.gtg<-((dn31.trinuc$ccg*dn31.trinuc$gtg)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.ccg.gtt<-((dn31.trinuc$ccg*dn31.trinuc$gtt)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv

#Stop codon
#cps.ccg.taa<-((dn31.trinuc$ccg*dn31.trinuc$taa)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccg.tac<-((dn31.trinuc$ccg*dn31.trinuc$tac)/(dn31.count.aa$p*dn31.count.aa$y))/dn31.count.diaa$py
#Stop codon
#cps.ccg.tag<-((dn31.trinuc$ccg*dn31.trinuc$tag)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccg.tat<-((dn31.trinuc$ccg*dn31.trinuc$tat)/(dn31.count.aa$p*dn31.count.aa$y))/dn31.count.diaa$py

cps.ccg.tca<-((dn31.trinuc$ccg*dn31.trinuc$tca)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.ccg.tcc<-((dn31.trinuc$ccg*dn31.trinuc$tcc)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.ccg.tcg<-((dn31.trinuc$ccg*dn31.trinuc$tcg)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.ccg.tct<-((dn31.trinuc$ccg*dn31.trinuc$tct)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps

#Stop codon
#cps.ccg.tga<-((dn31.trinuc$ccg*dn31.trinuc$tga)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.ccg.tgc<-((dn31.trinuc$ccg*dn31.trinuc$tgc)/(dn31.count.aa$p*dn31.count.aa$c))/dn31.count.diaa$pc
cps.ccg.tgg<-((dn31.trinuc$ccg*dn31.trinuc$tgg)/(dn31.count.aa$p*dn31.count.aa$w))/dn31.count.diaa$pw
cps.ccg.tgt<-((dn31.trinuc$ccg*dn31.trinuc$tgt)/(dn31.count.aa$p*dn31.count.aa$c))/dn31.count.diaa$pc

cps.ccg.tta<-((dn31.trinuc$ccg*dn31.trinuc$tta)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccg.ttc<-((dn31.trinuc$ccg*dn31.trinuc$ttc)/(dn31.count.aa$p*dn31.count.aa$f))/dn31.count.diaa$pf
cps.ccg.ttg<-((dn31.trinuc$ccg*dn31.trinuc$ttg)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.ccg.ttt<-((dn31.trinuc$ccg*dn31.trinuc$ttt)/(dn31.count.aa$p*dn31.count.aa$f))/dn31.count.diaa$pf








cps.cct.aaa<-((dn31.trinuc$cct*dn31.trinuc$aaa)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cct.aac<-((dn31.trinuc$cct*dn31.trinuc$aac)/(dn31.count.aa$p*dn31.count.aa$n))/dn31.count.diaa$pn
cps.cct.aag<-((dn31.trinuc$cct*dn31.trinuc$aag)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cct.aat<-((dn31.trinuc$cct*dn31.trinuc$aat)/(dn31.count.aa$p*dn31.count.aa$n))/dn31.count.diaa$pn

cps.cct.aca<-((dn31.trinuc$cct*dn31.trinuc$aca)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.cct.acc<-((dn31.trinuc$cct*dn31.trinuc$acc)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.cct.acg<-((dn31.trinuc$cct*dn31.trinuc$acg)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt
cps.cct.act<-((dn31.trinuc$cct*dn31.trinuc$act)/(dn31.count.aa$p*dn31.count.aa$t))/dn31.count.diaa$pt

cps.cct.aga<-((dn31.trinuc$cct*dn31.trinuc$aga)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cct.agc<-((dn31.trinuc$cct*dn31.trinuc$agc)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.cct.agg<-((dn31.trinuc$cct*dn31.trinuc$agg)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cct.agt<-((dn31.trinuc$cct*dn31.trinuc$agt)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps

cps.cct.ata<-((dn31.trinuc$cct*dn31.trinuc$ata)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi
cps.cct.atc<-((dn31.trinuc$cct*dn31.trinuc$atc)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi
cps.cct.atg<-((dn31.trinuc$cct*dn31.trinuc$atg)/(dn31.count.aa$p*dn31.count.aa$m))/dn31.count.diaa$pm
cps.cct.att<-((dn31.trinuc$cct*dn31.trinuc$att)/(dn31.count.aa$p*dn31.count.aa$i))/dn31.count.diaa$pi

cps.cct.caa<-((dn31.trinuc$cct*dn31.trinuc$caa)/(dn31.count.aa$p*dn31.count.aa$q))/dn31.count.diaa$pq
cps.cct.cac<-((dn31.trinuc$cct*dn31.trinuc$cac)/(dn31.count.aa$p*dn31.count.aa$h))/dn31.count.diaa$ph
cps.cct.cag<-((dn31.trinuc$cct*dn31.trinuc$cag)/(dn31.count.aa$p*dn31.count.aa$q))/dn31.count.diaa$pq
cps.cct.cat<-((dn31.trinuc$cct*dn31.trinuc$cat)/(dn31.count.aa$p*dn31.count.aa$h))/dn31.count.diaa$ph

cps.cct.cca<-((dn31.trinuc$cct*dn31.trinuc$cca)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.cct.ccc<-((dn31.trinuc$cct*dn31.trinuc$ccc)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.cct.ccg<-((dn31.trinuc$cct*dn31.trinuc$ccg)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp
cps.cct.cct<-((dn31.trinuc$cct*dn31.trinuc$cct)/(dn31.count.aa$p*dn31.count.aa$p))/dn31.count.diaa$pp

cps.cct.cga<-((dn31.trinuc$cct*dn31.trinuc$cga)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cct.cgc<-((dn31.trinuc$cct*dn31.trinuc$cgc)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cct.cgg<-((dn31.trinuc$cct*dn31.trinuc$cgg)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr
cps.cct.cgt<-((dn31.trinuc$cct*dn31.trinuc$cgt)/(dn31.count.aa$p*dn31.count.aa$r))/dn31.count.diaa$pr

cps.cct.cta<-((dn31.trinuc$cct*dn31.trinuc$cta)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cct.ctc<-((dn31.trinuc$cct*dn31.trinuc$ctc)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cct.ctg<-((dn31.trinuc$cct*dn31.trinuc$ctg)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cct.ctt<-((dn31.trinuc$cct*dn31.trinuc$ctt)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl

cps.cct.gaa<-((dn31.trinuc$cct*dn31.trinuc$gaa)/(dn31.count.aa$p*dn31.count.aa$e))/dn31.count.diaa$pe
cps.cct.gac<-((dn31.trinuc$cct*dn31.trinuc$gac)/(dn31.count.aa$p*dn31.count.aa$d))/dn31.count.diaa$pd
cps.cct.gag<-((dn31.trinuc$cct*dn31.trinuc$gag)/(dn31.count.aa$p*dn31.count.aa$e))/dn31.count.diaa$pe
cps.cct.gat<-((dn31.trinuc$cct*dn31.trinuc$gat)/(dn31.count.aa$p*dn31.count.aa$d))/dn31.count.diaa$pd

cps.cct.gca<-((dn31.trinuc$cct*dn31.trinuc$gca)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.cct.gcc<-((dn31.trinuc$cct*dn31.trinuc$gcc)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.cct.gcg<-((dn31.trinuc$cct*dn31.trinuc$gcg)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa
cps.cct.gct<-((dn31.trinuc$cct*dn31.trinuc$gct)/(dn31.count.aa$p*dn31.count.aa$a))/dn31.count.diaa$pa

cps.cct.gga<-((dn31.trinuc$cct*dn31.trinuc$gga)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.cct.ggc<-((dn31.trinuc$cct*dn31.trinuc$ggc)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.cct.ggg<-((dn31.trinuc$cct*dn31.trinuc$ggg)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg
cps.cct.ggt<-((dn31.trinuc$cct*dn31.trinuc$ggt)/(dn31.count.aa$p*dn31.count.aa$g))/dn31.count.diaa$pg

cps.cct.gta<-((dn31.trinuc$cct*dn31.trinuc$gta)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.cct.gtc<-((dn31.trinuc$cct*dn31.trinuc$gtc)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.cct.gtg<-((dn31.trinuc$cct*dn31.trinuc$gtg)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv
cps.cct.gtt<-((dn31.trinuc$cct*dn31.trinuc$gtt)/(dn31.count.aa$p*dn31.count.aa$v))/dn31.count.diaa$pv

#Stop codon
#cps.cct.taa<-((dn31.trinuc$cct*dn31.trinuc$taa)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cct.tac<-((dn31.trinuc$cct*dn31.trinuc$tac)/(dn31.count.aa$p*dn31.count.aa$y))/dn31.count.diaa$py
#Stop codon
#cps.cct.tag<-((dn31.trinuc$cct*dn31.trinuc$tag)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cct.tat<-((dn31.trinuc$cct*dn31.trinuc$tat)/(dn31.count.aa$p*dn31.count.aa$y))/dn31.count.diaa$py

cps.cct.tca<-((dn31.trinuc$cct*dn31.trinuc$tca)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.cct.tcc<-((dn31.trinuc$cct*dn31.trinuc$tcc)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.cct.tcg<-((dn31.trinuc$cct*dn31.trinuc$tcg)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps
cps.cct.tct<-((dn31.trinuc$cct*dn31.trinuc$tct)/(dn31.count.aa$p*dn31.count.aa$s))/dn31.count.diaa$ps

#Stop codon
#cps.cct.tga<-((dn31.trinuc$cct*dn31.trinuc$tga)/(dn31.count.aa$p*dn31.count.aa$k))/dn31.count.diaa$pk
cps.cct.tgc<-((dn31.trinuc$cct*dn31.trinuc$tgc)/(dn31.count.aa$p*dn31.count.aa$c))/dn31.count.diaa$pc
cps.cct.tgg<-((dn31.trinuc$cct*dn31.trinuc$tgg)/(dn31.count.aa$p*dn31.count.aa$w))/dn31.count.diaa$pw
cps.cct.tgt<-((dn31.trinuc$cct*dn31.trinuc$tgt)/(dn31.count.aa$p*dn31.count.aa$c))/dn31.count.diaa$pc

cps.cct.tta<-((dn31.trinuc$cct*dn31.trinuc$tta)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cct.ttc<-((dn31.trinuc$cct*dn31.trinuc$ttc)/(dn31.count.aa$p*dn31.count.aa$f))/dn31.count.diaa$pf
cps.cct.ttg<-((dn31.trinuc$cct*dn31.trinuc$ttg)/(dn31.count.aa$p*dn31.count.aa$l))/dn31.count.diaa$pl
cps.cct.ttt<-((dn31.trinuc$cct*dn31.trinuc$ttt)/(dn31.count.aa$p*dn31.count.aa$f))/dn31.count.diaa$pf








cps.cga.aaa<-((dn31.trinuc$cga*dn31.trinuc$aaa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.cga.aac<-((dn31.trinuc$cga*dn31.trinuc$aac)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn
cps.cga.aag<-((dn31.trinuc$cga*dn31.trinuc$aag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.cga.aat<-((dn31.trinuc$cga*dn31.trinuc$aat)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn

cps.cga.aca<-((dn31.trinuc$cga*dn31.trinuc$aca)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cga.acc<-((dn31.trinuc$cga*dn31.trinuc$acc)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cga.acg<-((dn31.trinuc$cga*dn31.trinuc$acg)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cga.act<-((dn31.trinuc$cga*dn31.trinuc$act)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt

cps.cga.aga<-((dn31.trinuc$cga*dn31.trinuc$aga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cga.agc<-((dn31.trinuc$cga*dn31.trinuc$agc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cga.agg<-((dn31.trinuc$cga*dn31.trinuc$agg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cga.agt<-((dn31.trinuc$cga*dn31.trinuc$agt)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

cps.cga.ata<-((dn31.trinuc$cga*dn31.trinuc$ata)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.cga.atc<-((dn31.trinuc$cga*dn31.trinuc$atc)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.cga.atg<-((dn31.trinuc$cga*dn31.trinuc$atg)/(dn31.count.aa$r*dn31.count.aa$m))/dn31.count.diaa$rm
cps.cga.att<-((dn31.trinuc$cga*dn31.trinuc$att)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri

cps.cga.caa<-((dn31.trinuc$cga*dn31.trinuc$caa)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.cga.cac<-((dn31.trinuc$cga*dn31.trinuc$cac)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh
cps.cga.cag<-((dn31.trinuc$cga*dn31.trinuc$cag)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.cga.cat<-((dn31.trinuc$cga*dn31.trinuc$cat)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh

cps.cga.cca<-((dn31.trinuc$cga*dn31.trinuc$cca)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cga.ccc<-((dn31.trinuc$cga*dn31.trinuc$ccc)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cga.ccg<-((dn31.trinuc$cga*dn31.trinuc$ccg)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cga.cct<-((dn31.trinuc$cga*dn31.trinuc$cct)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp

cps.cga.cga<-((dn31.trinuc$cga*dn31.trinuc$cga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cga.cgc<-((dn31.trinuc$cga*dn31.trinuc$cgc)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cga.cgg<-((dn31.trinuc$cga*dn31.trinuc$cgg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cga.cgt<-((dn31.trinuc$cga*dn31.trinuc$cgt)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr

cps.cga.cta<-((dn31.trinuc$cga*dn31.trinuc$cta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cga.ctc<-((dn31.trinuc$cga*dn31.trinuc$ctc)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cga.ctg<-((dn31.trinuc$cga*dn31.trinuc$ctg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cga.ctt<-((dn31.trinuc$cga*dn31.trinuc$ctt)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl

cps.cga.gaa<-((dn31.trinuc$cga*dn31.trinuc$gaa)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.cga.gac<-((dn31.trinuc$cga*dn31.trinuc$gac)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd
cps.cga.gag<-((dn31.trinuc$cga*dn31.trinuc$gag)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.cga.gat<-((dn31.trinuc$cga*dn31.trinuc$gat)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd

cps.cga.gca<-((dn31.trinuc$cga*dn31.trinuc$gca)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cga.gcc<-((dn31.trinuc$cga*dn31.trinuc$gcc)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cga.gcg<-((dn31.trinuc$cga*dn31.trinuc$gcg)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cga.gct<-((dn31.trinuc$cga*dn31.trinuc$gct)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra

cps.cga.gga<-((dn31.trinuc$cga*dn31.trinuc$gga)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cga.ggc<-((dn31.trinuc$cga*dn31.trinuc$ggc)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cga.ggg<-((dn31.trinuc$cga*dn31.trinuc$ggg)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cga.ggt<-((dn31.trinuc$cga*dn31.trinuc$ggt)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg

cps.cga.gta<-((dn31.trinuc$cga*dn31.trinuc$gta)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cga.gtc<-((dn31.trinuc$cga*dn31.trinuc$gtc)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cga.gtg<-((dn31.trinuc$cga*dn31.trinuc$gtg)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cga.gtt<-((dn31.trinuc$cga*dn31.trinuc$gtt)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv

#Stop codon
#cps.cga.taa<-((dn31.trinuc$cga*dn31.trinuc$taa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cga.tac<-((dn31.trinuc$cga*dn31.trinuc$tac)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry
#Stop codon
#cps.cga.tag<-((dn31.trinuc$cga*dn31.trinuc$tag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cga.tat<-((dn31.trinuc$cga*dn31.trinuc$tat)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry

cps.cga.tca<-((dn31.trinuc$cga*dn31.trinuc$tca)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cga.tcc<-((dn31.trinuc$cga*dn31.trinuc$tcc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cga.tcg<-((dn31.trinuc$cga*dn31.trinuc$tcg)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cga.tct<-((dn31.trinuc$cga*dn31.trinuc$tct)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

#Stop codon
#cps.cga.tga<-((dn31.trinuc$cga*dn31.trinuc$tga)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cga.tgc<-((dn31.trinuc$cga*dn31.trinuc$tgc)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc
cps.cga.tgg<-((dn31.trinuc$cga*dn31.trinuc$tgg)/(dn31.count.aa$r*dn31.count.aa$w))/dn31.count.diaa$rw
cps.cga.tgt<-((dn31.trinuc$cga*dn31.trinuc$tgt)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc

cps.cga.tta<-((dn31.trinuc$cga*dn31.trinuc$tta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cga.ttc<-((dn31.trinuc$cga*dn31.trinuc$ttc)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf
cps.cga.ttg<-((dn31.trinuc$cga*dn31.trinuc$ttg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cga.ttt<-((dn31.trinuc$cga*dn31.trinuc$ttt)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf








cps.cgc.aaa<-((dn31.trinuc$cgc*dn31.trinuc$aaa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.cgc.aac<-((dn31.trinuc$cgc*dn31.trinuc$aac)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn
cps.cgc.aag<-((dn31.trinuc$cgc*dn31.trinuc$aag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.cgc.aat<-((dn31.trinuc$cgc*dn31.trinuc$aat)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn

cps.cgc.aca<-((dn31.trinuc$cgc*dn31.trinuc$aca)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgc.acc<-((dn31.trinuc$cgc*dn31.trinuc$acc)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgc.acg<-((dn31.trinuc$cgc*dn31.trinuc$acg)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgc.act<-((dn31.trinuc$cgc*dn31.trinuc$act)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt

cps.cgc.aga<-((dn31.trinuc$cgc*dn31.trinuc$aga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgc.agc<-((dn31.trinuc$cgc*dn31.trinuc$agc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgc.agg<-((dn31.trinuc$cgc*dn31.trinuc$agg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgc.agt<-((dn31.trinuc$cgc*dn31.trinuc$agt)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

cps.cgc.ata<-((dn31.trinuc$cgc*dn31.trinuc$ata)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.cgc.atc<-((dn31.trinuc$cgc*dn31.trinuc$atc)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.cgc.atg<-((dn31.trinuc$cgc*dn31.trinuc$atg)/(dn31.count.aa$r*dn31.count.aa$m))/dn31.count.diaa$rm
cps.cgc.att<-((dn31.trinuc$cgc*dn31.trinuc$att)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri

cps.cgc.caa<-((dn31.trinuc$cgc*dn31.trinuc$caa)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.cgc.cac<-((dn31.trinuc$cgc*dn31.trinuc$cac)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh
cps.cgc.cag<-((dn31.trinuc$cgc*dn31.trinuc$cag)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.cgc.cat<-((dn31.trinuc$cgc*dn31.trinuc$cat)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh

cps.cgc.cca<-((dn31.trinuc$cgc*dn31.trinuc$cca)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgc.ccc<-((dn31.trinuc$cgc*dn31.trinuc$ccc)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgc.ccg<-((dn31.trinuc$cgc*dn31.trinuc$ccg)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgc.cct<-((dn31.trinuc$cgc*dn31.trinuc$cct)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp

cps.cgc.cga<-((dn31.trinuc$cgc*dn31.trinuc$cga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgc.cgc<-((dn31.trinuc$cgc*dn31.trinuc$cgc)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgc.cgg<-((dn31.trinuc$cgc*dn31.trinuc$cgg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgc.cgt<-((dn31.trinuc$cgc*dn31.trinuc$cgt)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr

cps.cgc.cta<-((dn31.trinuc$cgc*dn31.trinuc$cta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgc.ctc<-((dn31.trinuc$cgc*dn31.trinuc$ctc)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgc.ctg<-((dn31.trinuc$cgc*dn31.trinuc$ctg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgc.ctt<-((dn31.trinuc$cgc*dn31.trinuc$ctt)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl

cps.cgc.gaa<-((dn31.trinuc$cgc*dn31.trinuc$gaa)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.cgc.gac<-((dn31.trinuc$cgc*dn31.trinuc$gac)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd
cps.cgc.gag<-((dn31.trinuc$cgc*dn31.trinuc$gag)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.cgc.gat<-((dn31.trinuc$cgc*dn31.trinuc$gat)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd

cps.cgc.gca<-((dn31.trinuc$cgc*dn31.trinuc$gca)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgc.gcc<-((dn31.trinuc$cgc*dn31.trinuc$gcc)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgc.gcg<-((dn31.trinuc$cgc*dn31.trinuc$gcg)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgc.gct<-((dn31.trinuc$cgc*dn31.trinuc$gct)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra

cps.cgc.gga<-((dn31.trinuc$cgc*dn31.trinuc$gga)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgc.ggc<-((dn31.trinuc$cgc*dn31.trinuc$ggc)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgc.ggg<-((dn31.trinuc$cgc*dn31.trinuc$ggg)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgc.ggt<-((dn31.trinuc$cgc*dn31.trinuc$ggt)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg

cps.cgc.gta<-((dn31.trinuc$cgc*dn31.trinuc$gta)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgc.gtc<-((dn31.trinuc$cgc*dn31.trinuc$gtc)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgc.gtg<-((dn31.trinuc$cgc*dn31.trinuc$gtg)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgc.gtt<-((dn31.trinuc$cgc*dn31.trinuc$gtt)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv

#Stop codon
#cps.cgc.taa<-((dn31.trinuc$cgc*dn31.trinuc$taa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgc.tac<-((dn31.trinuc$cgc*dn31.trinuc$tac)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry
#Stop codon
#cps.cgc.tag<-((dn31.trinuc$cgc*dn31.trinuc$tag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgc.tat<-((dn31.trinuc$cgc*dn31.trinuc$tat)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry

cps.cgc.tca<-((dn31.trinuc$cgc*dn31.trinuc$tca)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgc.tcc<-((dn31.trinuc$cgc*dn31.trinuc$tcc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgc.tcg<-((dn31.trinuc$cgc*dn31.trinuc$tcg)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgc.tct<-((dn31.trinuc$cgc*dn31.trinuc$tct)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

#Stop codon
#cps.cgc.tga<-((dn31.trinuc$cgc*dn31.trinuc$tga)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgc.tgc<-((dn31.trinuc$cgc*dn31.trinuc$tgc)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc
cps.cgc.tgg<-((dn31.trinuc$cgc*dn31.trinuc$tgg)/(dn31.count.aa$r*dn31.count.aa$w))/dn31.count.diaa$rw
cps.cgc.tgt<-((dn31.trinuc$cgc*dn31.trinuc$tgt)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc

cps.cgc.tta<-((dn31.trinuc$cgc*dn31.trinuc$tta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgc.ttc<-((dn31.trinuc$cgc*dn31.trinuc$ttc)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf
cps.cgc.ttg<-((dn31.trinuc$cgc*dn31.trinuc$ttg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgc.ttt<-((dn31.trinuc$cgc*dn31.trinuc$ttt)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf








cps.cgg.aaa<-((dn31.trinuc$cgg*dn31.trinuc$aaa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.cgg.aac<-((dn31.trinuc$cgg*dn31.trinuc$aac)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn
cps.cgg.aag<-((dn31.trinuc$cgg*dn31.trinuc$aag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.cgg.aat<-((dn31.trinuc$cgg*dn31.trinuc$aat)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn

cps.cgg.aca<-((dn31.trinuc$cgg*dn31.trinuc$aca)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgg.acc<-((dn31.trinuc$cgg*dn31.trinuc$acc)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgg.acg<-((dn31.trinuc$cgg*dn31.trinuc$acg)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgg.act<-((dn31.trinuc$cgg*dn31.trinuc$act)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt

cps.cgg.aga<-((dn31.trinuc$cgg*dn31.trinuc$aga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgg.agc<-((dn31.trinuc$cgg*dn31.trinuc$agc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgg.agg<-((dn31.trinuc$cgg*dn31.trinuc$agg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgg.agt<-((dn31.trinuc$cgg*dn31.trinuc$agt)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

cps.cgg.ata<-((dn31.trinuc$cgg*dn31.trinuc$ata)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.cgg.atc<-((dn31.trinuc$cgg*dn31.trinuc$atc)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.cgg.atg<-((dn31.trinuc$cgg*dn31.trinuc$atg)/(dn31.count.aa$r*dn31.count.aa$m))/dn31.count.diaa$rm
cps.cgg.att<-((dn31.trinuc$cgg*dn31.trinuc$att)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri

cps.cgg.caa<-((dn31.trinuc$cgg*dn31.trinuc$caa)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.cgg.cac<-((dn31.trinuc$cgg*dn31.trinuc$cac)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh
cps.cgg.cag<-((dn31.trinuc$cgg*dn31.trinuc$cag)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.cgg.cat<-((dn31.trinuc$cgg*dn31.trinuc$cat)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh

cps.cgg.cca<-((dn31.trinuc$cgg*dn31.trinuc$cca)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgg.ccc<-((dn31.trinuc$cgg*dn31.trinuc$ccc)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgg.ccg<-((dn31.trinuc$cgg*dn31.trinuc$ccg)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgg.cct<-((dn31.trinuc$cgg*dn31.trinuc$cct)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp

cps.cgg.cga<-((dn31.trinuc$cgg*dn31.trinuc$cga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgg.cgc<-((dn31.trinuc$cgg*dn31.trinuc$cgc)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgg.cgg<-((dn31.trinuc$cgg*dn31.trinuc$cgg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgg.cgt<-((dn31.trinuc$cgg*dn31.trinuc$cgt)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr

cps.cgg.cta<-((dn31.trinuc$cgg*dn31.trinuc$cta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgg.ctc<-((dn31.trinuc$cgg*dn31.trinuc$ctc)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgg.ctg<-((dn31.trinuc$cgg*dn31.trinuc$ctg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgg.ctt<-((dn31.trinuc$cgg*dn31.trinuc$ctt)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl

cps.cgg.gaa<-((dn31.trinuc$cgg*dn31.trinuc$gaa)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.cgg.gac<-((dn31.trinuc$cgg*dn31.trinuc$gac)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd
cps.cgg.gag<-((dn31.trinuc$cgg*dn31.trinuc$gag)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.cgg.gat<-((dn31.trinuc$cgg*dn31.trinuc$gat)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd

cps.cgg.gca<-((dn31.trinuc$cgg*dn31.trinuc$gca)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgg.gcc<-((dn31.trinuc$cgg*dn31.trinuc$gcc)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgg.gcg<-((dn31.trinuc$cgg*dn31.trinuc$gcg)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgg.gct<-((dn31.trinuc$cgg*dn31.trinuc$gct)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra

cps.cgg.gga<-((dn31.trinuc$cgg*dn31.trinuc$gga)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgg.ggc<-((dn31.trinuc$cgg*dn31.trinuc$ggc)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgg.ggg<-((dn31.trinuc$cgg*dn31.trinuc$ggg)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgg.ggt<-((dn31.trinuc$cgg*dn31.trinuc$ggt)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg

cps.cgg.gta<-((dn31.trinuc$cgg*dn31.trinuc$gta)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgg.gtc<-((dn31.trinuc$cgg*dn31.trinuc$gtc)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgg.gtg<-((dn31.trinuc$cgg*dn31.trinuc$gtg)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgg.gtt<-((dn31.trinuc$cgg*dn31.trinuc$gtt)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv

#Stop codon
#cps.cgg.taa<-((dn31.trinuc$cgg*dn31.trinuc$taa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgg.tac<-((dn31.trinuc$cgg*dn31.trinuc$tac)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry
#Stop codon
#cps.cgg.tag<-((dn31.trinuc$cgg*dn31.trinuc$tag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgg.tat<-((dn31.trinuc$cgg*dn31.trinuc$tat)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry

cps.cgg.tca<-((dn31.trinuc$cgg*dn31.trinuc$tca)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgg.tcc<-((dn31.trinuc$cgg*dn31.trinuc$tcc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgg.tcg<-((dn31.trinuc$cgg*dn31.trinuc$tcg)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgg.tct<-((dn31.trinuc$cgg*dn31.trinuc$tct)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

#Stop codon
#cps.cgg.tga<-((dn31.trinuc$cgg*dn31.trinuc$tga)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgg.tgc<-((dn31.trinuc$cgg*dn31.trinuc$tgc)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc
cps.cgg.tgg<-((dn31.trinuc$cgg*dn31.trinuc$tgg)/(dn31.count.aa$r*dn31.count.aa$w))/dn31.count.diaa$rw
cps.cgg.tgt<-((dn31.trinuc$cgg*dn31.trinuc$tgt)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc

cps.cgg.tta<-((dn31.trinuc$cgg*dn31.trinuc$tta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgg.ttc<-((dn31.trinuc$cgg*dn31.trinuc$ttc)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf
cps.cgg.ttg<-((dn31.trinuc$cgg*dn31.trinuc$ttg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgg.ttt<-((dn31.trinuc$cgg*dn31.trinuc$ttt)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf








cps.cgt.aaa<-((dn31.trinuc$cgt*dn31.trinuc$aaa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.cgt.aac<-((dn31.trinuc$cgt*dn31.trinuc$aac)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn
cps.cgt.aag<-((dn31.trinuc$cgt*dn31.trinuc$aag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$rk
cps.cgt.aat<-((dn31.trinuc$cgt*dn31.trinuc$aat)/(dn31.count.aa$r*dn31.count.aa$n))/dn31.count.diaa$rn

cps.cgt.aca<-((dn31.trinuc$cgt*dn31.trinuc$aca)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgt.acc<-((dn31.trinuc$cgt*dn31.trinuc$acc)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgt.acg<-((dn31.trinuc$cgt*dn31.trinuc$acg)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt
cps.cgt.act<-((dn31.trinuc$cgt*dn31.trinuc$act)/(dn31.count.aa$r*dn31.count.aa$t))/dn31.count.diaa$rt

cps.cgt.aga<-((dn31.trinuc$cgt*dn31.trinuc$aga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgt.agc<-((dn31.trinuc$cgt*dn31.trinuc$agc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgt.agg<-((dn31.trinuc$cgt*dn31.trinuc$agg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgt.agt<-((dn31.trinuc$cgt*dn31.trinuc$agt)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

cps.cgt.ata<-((dn31.trinuc$cgt*dn31.trinuc$ata)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.cgt.atc<-((dn31.trinuc$cgt*dn31.trinuc$atc)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri
cps.cgt.atg<-((dn31.trinuc$cgt*dn31.trinuc$atg)/(dn31.count.aa$r*dn31.count.aa$m))/dn31.count.diaa$rm
cps.cgt.att<-((dn31.trinuc$cgt*dn31.trinuc$att)/(dn31.count.aa$r*dn31.count.aa$i))/dn31.count.diaa$ri

cps.cgt.caa<-((dn31.trinuc$cgt*dn31.trinuc$caa)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.cgt.cac<-((dn31.trinuc$cgt*dn31.trinuc$cac)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh
cps.cgt.cag<-((dn31.trinuc$cgt*dn31.trinuc$cag)/(dn31.count.aa$r*dn31.count.aa$q))/dn31.count.diaa$rq
cps.cgt.cat<-((dn31.trinuc$cgt*dn31.trinuc$cat)/(dn31.count.aa$r*dn31.count.aa$h))/dn31.count.diaa$rh

cps.cgt.cca<-((dn31.trinuc$cgt*dn31.trinuc$cca)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgt.ccc<-((dn31.trinuc$cgt*dn31.trinuc$ccc)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgt.ccg<-((dn31.trinuc$cgt*dn31.trinuc$ccg)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp
cps.cgt.cct<-((dn31.trinuc$cgt*dn31.trinuc$cct)/(dn31.count.aa$r*dn31.count.aa$p))/dn31.count.diaa$rp

cps.cgt.cga<-((dn31.trinuc$cgt*dn31.trinuc$cga)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgt.cgc<-((dn31.trinuc$cgt*dn31.trinuc$cgc)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgt.cgg<-((dn31.trinuc$cgt*dn31.trinuc$cgg)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr
cps.cgt.cgt<-((dn31.trinuc$cgt*dn31.trinuc$cgt)/(dn31.count.aa$r*dn31.count.aa$r))/dn31.count.diaa$rr

cps.cgt.cta<-((dn31.trinuc$cgt*dn31.trinuc$cta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgt.ctc<-((dn31.trinuc$cgt*dn31.trinuc$ctc)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgt.ctg<-((dn31.trinuc$cgt*dn31.trinuc$ctg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgt.ctt<-((dn31.trinuc$cgt*dn31.trinuc$ctt)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl

cps.cgt.gaa<-((dn31.trinuc$cgt*dn31.trinuc$gaa)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.cgt.gac<-((dn31.trinuc$cgt*dn31.trinuc$gac)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd
cps.cgt.gag<-((dn31.trinuc$cgt*dn31.trinuc$gag)/(dn31.count.aa$r*dn31.count.aa$e))/dn31.count.diaa$re
cps.cgt.gat<-((dn31.trinuc$cgt*dn31.trinuc$gat)/(dn31.count.aa$r*dn31.count.aa$d))/dn31.count.diaa$rd

cps.cgt.gca<-((dn31.trinuc$cgt*dn31.trinuc$gca)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgt.gcc<-((dn31.trinuc$cgt*dn31.trinuc$gcc)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgt.gcg<-((dn31.trinuc$cgt*dn31.trinuc$gcg)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra
cps.cgt.gct<-((dn31.trinuc$cgt*dn31.trinuc$gct)/(dn31.count.aa$r*dn31.count.aa$a))/dn31.count.diaa$ra

cps.cgt.gga<-((dn31.trinuc$cgt*dn31.trinuc$gga)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgt.ggc<-((dn31.trinuc$cgt*dn31.trinuc$ggc)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgt.ggg<-((dn31.trinuc$cgt*dn31.trinuc$ggg)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg
cps.cgt.ggt<-((dn31.trinuc$cgt*dn31.trinuc$ggt)/(dn31.count.aa$r*dn31.count.aa$g))/dn31.count.diaa$rg

cps.cgt.gta<-((dn31.trinuc$cgt*dn31.trinuc$gta)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgt.gtc<-((dn31.trinuc$cgt*dn31.trinuc$gtc)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgt.gtg<-((dn31.trinuc$cgt*dn31.trinuc$gtg)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv
cps.cgt.gtt<-((dn31.trinuc$cgt*dn31.trinuc$gtt)/(dn31.count.aa$r*dn31.count.aa$v))/dn31.count.diaa$rv

#Stop codon
#cps.cgt.taa<-((dn31.trinuc$cgt*dn31.trinuc$taa)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgt.tac<-((dn31.trinuc$cgt*dn31.trinuc$tac)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry
#Stop codon
#cps.cgt.tag<-((dn31.trinuc$cgt*dn31.trinuc$tag)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgt.tat<-((dn31.trinuc$cgt*dn31.trinuc$tat)/(dn31.count.aa$r*dn31.count.aa$y))/dn31.count.diaa$ry

cps.cgt.tca<-((dn31.trinuc$cgt*dn31.trinuc$tca)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgt.tcc<-((dn31.trinuc$cgt*dn31.trinuc$tcc)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgt.tcg<-((dn31.trinuc$cgt*dn31.trinuc$tcg)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs
cps.cgt.tct<-((dn31.trinuc$cgt*dn31.trinuc$tct)/(dn31.count.aa$r*dn31.count.aa$s))/dn31.count.diaa$rs

#Stop codon
#cps.cgt.tga<-((dn31.trinuc$cgt*dn31.trinuc$tga)/(dn31.count.aa$r*dn31.count.aa$k))/dn31.count.diaa$kk
cps.cgt.tgc<-((dn31.trinuc$cgt*dn31.trinuc$tgc)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc
cps.cgt.tgg<-((dn31.trinuc$cgt*dn31.trinuc$tgg)/(dn31.count.aa$r*dn31.count.aa$w))/dn31.count.diaa$rw
cps.cgt.tgt<-((dn31.trinuc$cgt*dn31.trinuc$tgt)/(dn31.count.aa$r*dn31.count.aa$c))/dn31.count.diaa$rc

cps.cgt.tta<-((dn31.trinuc$cgt*dn31.trinuc$tta)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgt.ttc<-((dn31.trinuc$cgt*dn31.trinuc$ttc)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf
cps.cgt.ttg<-((dn31.trinuc$cgt*dn31.trinuc$ttg)/(dn31.count.aa$r*dn31.count.aa$l))/dn31.count.diaa$rl
cps.cgt.ttt<-((dn31.trinuc$cgt*dn31.trinuc$ttt)/(dn31.count.aa$r*dn31.count.aa$f))/dn31.count.diaa$rf








cps.cta.aaa<-((dn31.trinuc$cta*dn31.trinuc$aaa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.cta.aac<-((dn31.trinuc$cta*dn31.trinuc$aac)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln
cps.cta.aag<-((dn31.trinuc$cta*dn31.trinuc$aag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.cta.aat<-((dn31.trinuc$cta*dn31.trinuc$aat)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln

cps.cta.aca<-((dn31.trinuc$cta*dn31.trinuc$aca)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.cta.acc<-((dn31.trinuc$cta*dn31.trinuc$acc)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.cta.acg<-((dn31.trinuc$cta*dn31.trinuc$acg)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.cta.act<-((dn31.trinuc$cta*dn31.trinuc$act)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt

cps.cta.aga<-((dn31.trinuc$cta*dn31.trinuc$aga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.cta.agc<-((dn31.trinuc$cta*dn31.trinuc$agc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.cta.agg<-((dn31.trinuc$cta*dn31.trinuc$agg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.cta.agt<-((dn31.trinuc$cta*dn31.trinuc$agt)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

cps.cta.ata<-((dn31.trinuc$cta*dn31.trinuc$ata)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.cta.atc<-((dn31.trinuc$cta*dn31.trinuc$atc)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.cta.atg<-((dn31.trinuc$cta*dn31.trinuc$atg)/(dn31.count.aa$l*dn31.count.aa$m))/dn31.count.diaa$lm
cps.cta.att<-((dn31.trinuc$cta*dn31.trinuc$att)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li

cps.cta.caa<-((dn31.trinuc$cta*dn31.trinuc$caa)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.cta.cac<-((dn31.trinuc$cta*dn31.trinuc$cac)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh
cps.cta.cag<-((dn31.trinuc$cta*dn31.trinuc$cag)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.cta.cat<-((dn31.trinuc$cta*dn31.trinuc$cat)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh

cps.cta.cca<-((dn31.trinuc$cta*dn31.trinuc$cca)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.cta.ccc<-((dn31.trinuc$cta*dn31.trinuc$ccc)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.cta.ccg<-((dn31.trinuc$cta*dn31.trinuc$ccg)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.cta.cct<-((dn31.trinuc$cta*dn31.trinuc$cct)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp

cps.cta.cga<-((dn31.trinuc$cta*dn31.trinuc$cga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.cta.cgc<-((dn31.trinuc$cta*dn31.trinuc$cgc)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.cta.cgg<-((dn31.trinuc$cta*dn31.trinuc$cgg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.cta.cgt<-((dn31.trinuc$cta*dn31.trinuc$cgt)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr

cps.cta.cta<-((dn31.trinuc$cta*dn31.trinuc$cta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.cta.ctc<-((dn31.trinuc$cta*dn31.trinuc$ctc)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.cta.ctg<-((dn31.trinuc$cta*dn31.trinuc$ctg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.cta.ctt<-((dn31.trinuc$cta*dn31.trinuc$ctt)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll

cps.cta.gaa<-((dn31.trinuc$cta*dn31.trinuc$gaa)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.cta.gac<-((dn31.trinuc$cta*dn31.trinuc$gac)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld
cps.cta.gag<-((dn31.trinuc$cta*dn31.trinuc$gag)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.cta.gat<-((dn31.trinuc$cta*dn31.trinuc$gat)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld

cps.cta.gca<-((dn31.trinuc$cta*dn31.trinuc$gca)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.cta.gcc<-((dn31.trinuc$cta*dn31.trinuc$gcc)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.cta.gcg<-((dn31.trinuc$cta*dn31.trinuc$gcg)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.cta.gct<-((dn31.trinuc$cta*dn31.trinuc$gct)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la

cps.cta.gga<-((dn31.trinuc$cta*dn31.trinuc$gga)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.cta.ggc<-((dn31.trinuc$cta*dn31.trinuc$ggc)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.cta.ggg<-((dn31.trinuc$cta*dn31.trinuc$ggg)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.cta.ggt<-((dn31.trinuc$cta*dn31.trinuc$ggt)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg

cps.cta.gta<-((dn31.trinuc$cta*dn31.trinuc$gta)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.cta.gtc<-((dn31.trinuc$cta*dn31.trinuc$gtc)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.cta.gtg<-((dn31.trinuc$cta*dn31.trinuc$gtg)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.cta.gtt<-((dn31.trinuc$cta*dn31.trinuc$gtt)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv

#Stop codon
#cps.cta.taa<-((dn31.trinuc$cta*dn31.trinuc$taa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.cta.tac<-((dn31.trinuc$cta*dn31.trinuc$tac)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly
#Stop codon
#cps.cta.tag<-((dn31.trinuc$cta*dn31.trinuc$tag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.cta.tat<-((dn31.trinuc$cta*dn31.trinuc$tat)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly

cps.cta.tca<-((dn31.trinuc$cta*dn31.trinuc$tca)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.cta.tcc<-((dn31.trinuc$cta*dn31.trinuc$tcc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.cta.tcg<-((dn31.trinuc$cta*dn31.trinuc$tcg)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.cta.tct<-((dn31.trinuc$cta*dn31.trinuc$tct)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

#Stop codon
#cps.cta.tga<-((dn31.trinuc$cta*dn31.trinuc$tga)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.cta.tgc<-((dn31.trinuc$cta*dn31.trinuc$tgc)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc
cps.cta.tgg<-((dn31.trinuc$cta*dn31.trinuc$tgg)/(dn31.count.aa$l*dn31.count.aa$w))/dn31.count.diaa$lw
cps.cta.tgt<-((dn31.trinuc$cta*dn31.trinuc$tgt)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc

cps.cta.tta<-((dn31.trinuc$cta*dn31.trinuc$tta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.cta.ttc<-((dn31.trinuc$cta*dn31.trinuc$ttc)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf
cps.cta.ttg<-((dn31.trinuc$cta*dn31.trinuc$ttg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.cta.ttt<-((dn31.trinuc$cta*dn31.trinuc$ttt)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf








cps.ctc.aaa<-((dn31.trinuc$ctc*dn31.trinuc$aaa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctc.aac<-((dn31.trinuc$ctc*dn31.trinuc$aac)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln
cps.ctc.aag<-((dn31.trinuc$ctc*dn31.trinuc$aag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctc.aat<-((dn31.trinuc$ctc*dn31.trinuc$aat)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln

cps.ctc.aca<-((dn31.trinuc$ctc*dn31.trinuc$aca)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctc.acc<-((dn31.trinuc$ctc*dn31.trinuc$acc)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctc.acg<-((dn31.trinuc$ctc*dn31.trinuc$acg)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctc.act<-((dn31.trinuc$ctc*dn31.trinuc$act)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt

cps.ctc.aga<-((dn31.trinuc$ctc*dn31.trinuc$aga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctc.agc<-((dn31.trinuc$ctc*dn31.trinuc$agc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctc.agg<-((dn31.trinuc$ctc*dn31.trinuc$agg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctc.agt<-((dn31.trinuc$ctc*dn31.trinuc$agt)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

cps.ctc.ata<-((dn31.trinuc$ctc*dn31.trinuc$ata)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.ctc.atc<-((dn31.trinuc$ctc*dn31.trinuc$atc)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.ctc.atg<-((dn31.trinuc$ctc*dn31.trinuc$atg)/(dn31.count.aa$l*dn31.count.aa$m))/dn31.count.diaa$lm
cps.ctc.att<-((dn31.trinuc$ctc*dn31.trinuc$att)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li

cps.ctc.caa<-((dn31.trinuc$ctc*dn31.trinuc$caa)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.ctc.cac<-((dn31.trinuc$ctc*dn31.trinuc$cac)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh
cps.ctc.cag<-((dn31.trinuc$ctc*dn31.trinuc$cag)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.ctc.cat<-((dn31.trinuc$ctc*dn31.trinuc$cat)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh

cps.ctc.cca<-((dn31.trinuc$ctc*dn31.trinuc$cca)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctc.ccc<-((dn31.trinuc$ctc*dn31.trinuc$ccc)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctc.ccg<-((dn31.trinuc$ctc*dn31.trinuc$ccg)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctc.cct<-((dn31.trinuc$ctc*dn31.trinuc$cct)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp

cps.ctc.cga<-((dn31.trinuc$ctc*dn31.trinuc$cga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctc.cgc<-((dn31.trinuc$ctc*dn31.trinuc$cgc)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctc.cgg<-((dn31.trinuc$ctc*dn31.trinuc$cgg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctc.cgt<-((dn31.trinuc$ctc*dn31.trinuc$cgt)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr

cps.ctc.cta<-((dn31.trinuc$ctc*dn31.trinuc$cta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctc.ctc<-((dn31.trinuc$ctc*dn31.trinuc$ctc)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctc.ctg<-((dn31.trinuc$ctc*dn31.trinuc$ctg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctc.ctt<-((dn31.trinuc$ctc*dn31.trinuc$ctt)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll

cps.ctc.gaa<-((dn31.trinuc$ctc*dn31.trinuc$gaa)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.ctc.gac<-((dn31.trinuc$ctc*dn31.trinuc$gac)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld
cps.ctc.gag<-((dn31.trinuc$ctc*dn31.trinuc$gag)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.ctc.gat<-((dn31.trinuc$ctc*dn31.trinuc$gat)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld

cps.ctc.gca<-((dn31.trinuc$ctc*dn31.trinuc$gca)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctc.gcc<-((dn31.trinuc$ctc*dn31.trinuc$gcc)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctc.gcg<-((dn31.trinuc$ctc*dn31.trinuc$gcg)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctc.gct<-((dn31.trinuc$ctc*dn31.trinuc$gct)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la

cps.ctc.gga<-((dn31.trinuc$ctc*dn31.trinuc$gga)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctc.ggc<-((dn31.trinuc$ctc*dn31.trinuc$ggc)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctc.ggg<-((dn31.trinuc$ctc*dn31.trinuc$ggg)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctc.ggt<-((dn31.trinuc$ctc*dn31.trinuc$ggt)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg

cps.ctc.gta<-((dn31.trinuc$ctc*dn31.trinuc$gta)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctc.gtc<-((dn31.trinuc$ctc*dn31.trinuc$gtc)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctc.gtg<-((dn31.trinuc$ctc*dn31.trinuc$gtg)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctc.gtt<-((dn31.trinuc$ctc*dn31.trinuc$gtt)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv

#Stop codon
#cps.ctc.taa<-((dn31.trinuc$ctc*dn31.trinuc$taa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctc.tac<-((dn31.trinuc$ctc*dn31.trinuc$tac)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly
#Stop codon
#cps.ctc.tag<-((dn31.trinuc$ctc*dn31.trinuc$tag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctc.tat<-((dn31.trinuc$ctc*dn31.trinuc$tat)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly

cps.ctc.tca<-((dn31.trinuc$ctc*dn31.trinuc$tca)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctc.tcc<-((dn31.trinuc$ctc*dn31.trinuc$tcc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctc.tcg<-((dn31.trinuc$ctc*dn31.trinuc$tcg)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctc.tct<-((dn31.trinuc$ctc*dn31.trinuc$tct)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

#Stop codon
#cps.ctc.tga<-((dn31.trinuc$ctc*dn31.trinuc$tga)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctc.tgc<-((dn31.trinuc$ctc*dn31.trinuc$tgc)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc
cps.ctc.tgg<-((dn31.trinuc$ctc*dn31.trinuc$tgg)/(dn31.count.aa$l*dn31.count.aa$w))/dn31.count.diaa$lw
cps.ctc.tgt<-((dn31.trinuc$ctc*dn31.trinuc$tgt)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc

cps.ctc.tta<-((dn31.trinuc$ctc*dn31.trinuc$tta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctc.ttc<-((dn31.trinuc$ctc*dn31.trinuc$ttc)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf
cps.ctc.ttg<-((dn31.trinuc$ctc*dn31.trinuc$ttg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctc.ttt<-((dn31.trinuc$ctc*dn31.trinuc$ttt)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf








cps.ctg.aaa<-((dn31.trinuc$ctg*dn31.trinuc$aaa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctg.aac<-((dn31.trinuc$ctg*dn31.trinuc$aac)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln
cps.ctg.aag<-((dn31.trinuc$ctg*dn31.trinuc$aag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctg.aat<-((dn31.trinuc$ctg*dn31.trinuc$aat)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln

cps.ctg.aca<-((dn31.trinuc$ctg*dn31.trinuc$aca)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctg.acc<-((dn31.trinuc$ctg*dn31.trinuc$acc)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctg.acg<-((dn31.trinuc$ctg*dn31.trinuc$acg)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctg.act<-((dn31.trinuc$ctg*dn31.trinuc$act)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt

cps.ctg.aga<-((dn31.trinuc$ctg*dn31.trinuc$aga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctg.agc<-((dn31.trinuc$ctg*dn31.trinuc$agc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctg.agg<-((dn31.trinuc$ctg*dn31.trinuc$agg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctg.agt<-((dn31.trinuc$ctg*dn31.trinuc$agt)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

cps.ctg.ata<-((dn31.trinuc$ctg*dn31.trinuc$ata)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.ctg.atc<-((dn31.trinuc$ctg*dn31.trinuc$atc)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.ctg.atg<-((dn31.trinuc$ctg*dn31.trinuc$atg)/(dn31.count.aa$l*dn31.count.aa$m))/dn31.count.diaa$lm
cps.ctg.att<-((dn31.trinuc$ctg*dn31.trinuc$att)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li

cps.ctg.caa<-((dn31.trinuc$ctg*dn31.trinuc$caa)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.ctg.cac<-((dn31.trinuc$ctg*dn31.trinuc$cac)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh
cps.ctg.cag<-((dn31.trinuc$ctg*dn31.trinuc$cag)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.ctg.cat<-((dn31.trinuc$ctg*dn31.trinuc$cat)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh

cps.ctg.cca<-((dn31.trinuc$ctg*dn31.trinuc$cca)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctg.ccc<-((dn31.trinuc$ctg*dn31.trinuc$ccc)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctg.ccg<-((dn31.trinuc$ctg*dn31.trinuc$ccg)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctg.cct<-((dn31.trinuc$ctg*dn31.trinuc$cct)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp

cps.ctg.cga<-((dn31.trinuc$ctg*dn31.trinuc$cga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctg.cgc<-((dn31.trinuc$ctg*dn31.trinuc$cgc)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctg.cgg<-((dn31.trinuc$ctg*dn31.trinuc$cgg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctg.cgt<-((dn31.trinuc$ctg*dn31.trinuc$cgt)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr

cps.ctg.cta<-((dn31.trinuc$ctg*dn31.trinuc$cta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctg.ctc<-((dn31.trinuc$ctg*dn31.trinuc$ctc)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctg.ctg<-((dn31.trinuc$ctg*dn31.trinuc$ctg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctg.ctt<-((dn31.trinuc$ctg*dn31.trinuc$ctt)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll

cps.ctg.gaa<-((dn31.trinuc$ctg*dn31.trinuc$gaa)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.ctg.gac<-((dn31.trinuc$ctg*dn31.trinuc$gac)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld
cps.ctg.gag<-((dn31.trinuc$ctg*dn31.trinuc$gag)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.ctg.gat<-((dn31.trinuc$ctg*dn31.trinuc$gat)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld

cps.ctg.gca<-((dn31.trinuc$ctg*dn31.trinuc$gca)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctg.gcc<-((dn31.trinuc$ctg*dn31.trinuc$gcc)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctg.gcg<-((dn31.trinuc$ctg*dn31.trinuc$gcg)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctg.gct<-((dn31.trinuc$ctg*dn31.trinuc$gct)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la

cps.ctg.gga<-((dn31.trinuc$ctg*dn31.trinuc$gga)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctg.ggc<-((dn31.trinuc$ctg*dn31.trinuc$ggc)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctg.ggg<-((dn31.trinuc$ctg*dn31.trinuc$ggg)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctg.ggt<-((dn31.trinuc$ctg*dn31.trinuc$ggt)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg

cps.ctg.gta<-((dn31.trinuc$ctg*dn31.trinuc$gta)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctg.gtc<-((dn31.trinuc$ctg*dn31.trinuc$gtc)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctg.gtg<-((dn31.trinuc$ctg*dn31.trinuc$gtg)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctg.gtt<-((dn31.trinuc$ctg*dn31.trinuc$gtt)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv

#Stop codon
#cps.ctg.taa<-((dn31.trinuc$ctg*dn31.trinuc$taa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctg.tac<-((dn31.trinuc$ctg*dn31.trinuc$tac)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly
#Stop codon
#cps.ctg.tag<-((dn31.trinuc$ctg*dn31.trinuc$tag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctg.tat<-((dn31.trinuc$ctg*dn31.trinuc$tat)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly

cps.ctg.tca<-((dn31.trinuc$ctg*dn31.trinuc$tca)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctg.tcc<-((dn31.trinuc$ctg*dn31.trinuc$tcc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctg.tcg<-((dn31.trinuc$ctg*dn31.trinuc$tcg)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctg.tct<-((dn31.trinuc$ctg*dn31.trinuc$tct)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

#Stop codon
#cps.ctg.tga<-((dn31.trinuc$ctg*dn31.trinuc$tga)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctg.tgc<-((dn31.trinuc$ctg*dn31.trinuc$tgc)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc
cps.ctg.tgg<-((dn31.trinuc$ctg*dn31.trinuc$tgg)/(dn31.count.aa$l*dn31.count.aa$w))/dn31.count.diaa$lw
cps.ctg.tgt<-((dn31.trinuc$ctg*dn31.trinuc$tgt)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc

cps.ctg.tta<-((dn31.trinuc$ctg*dn31.trinuc$tta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctg.ttc<-((dn31.trinuc$ctg*dn31.trinuc$ttc)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf
cps.ctg.ttg<-((dn31.trinuc$ctg*dn31.trinuc$ttg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctg.ttt<-((dn31.trinuc$ctg*dn31.trinuc$ttt)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf








cps.ctt.aaa<-((dn31.trinuc$ctt*dn31.trinuc$aaa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctt.aac<-((dn31.trinuc$ctt*dn31.trinuc$aac)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln
cps.ctt.aag<-((dn31.trinuc$ctt*dn31.trinuc$aag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctt.aat<-((dn31.trinuc$ctt*dn31.trinuc$aat)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln

cps.ctt.aca<-((dn31.trinuc$ctt*dn31.trinuc$aca)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctt.acc<-((dn31.trinuc$ctt*dn31.trinuc$acc)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctt.acg<-((dn31.trinuc$ctt*dn31.trinuc$acg)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ctt.act<-((dn31.trinuc$ctt*dn31.trinuc$act)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt

cps.ctt.aga<-((dn31.trinuc$ctt*dn31.trinuc$aga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctt.agc<-((dn31.trinuc$ctt*dn31.trinuc$agc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctt.agg<-((dn31.trinuc$ctt*dn31.trinuc$agg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctt.agt<-((dn31.trinuc$ctt*dn31.trinuc$agt)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

cps.ctt.ata<-((dn31.trinuc$ctt*dn31.trinuc$ata)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.ctt.atc<-((dn31.trinuc$ctt*dn31.trinuc$atc)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.ctt.atg<-((dn31.trinuc$ctt*dn31.trinuc$atg)/(dn31.count.aa$l*dn31.count.aa$m))/dn31.count.diaa$lm
cps.ctt.att<-((dn31.trinuc$ctt*dn31.trinuc$att)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li

cps.ctt.caa<-((dn31.trinuc$ctt*dn31.trinuc$caa)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.ctt.cac<-((dn31.trinuc$ctt*dn31.trinuc$cac)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh
cps.ctt.cag<-((dn31.trinuc$ctt*dn31.trinuc$cag)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.ctt.cat<-((dn31.trinuc$ctt*dn31.trinuc$cat)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh

cps.ctt.cca<-((dn31.trinuc$ctt*dn31.trinuc$cca)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctt.ccc<-((dn31.trinuc$ctt*dn31.trinuc$ccc)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctt.ccg<-((dn31.trinuc$ctt*dn31.trinuc$ccg)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ctt.cct<-((dn31.trinuc$ctt*dn31.trinuc$cct)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp

cps.ctt.cga<-((dn31.trinuc$ctt*dn31.trinuc$cga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctt.cgc<-((dn31.trinuc$ctt*dn31.trinuc$cgc)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctt.cgg<-((dn31.trinuc$ctt*dn31.trinuc$cgg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ctt.cgt<-((dn31.trinuc$ctt*dn31.trinuc$cgt)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr

cps.ctt.cta<-((dn31.trinuc$ctt*dn31.trinuc$cta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctt.ctc<-((dn31.trinuc$ctt*dn31.trinuc$ctc)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctt.ctg<-((dn31.trinuc$ctt*dn31.trinuc$ctg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctt.ctt<-((dn31.trinuc$ctt*dn31.trinuc$ctt)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll

cps.ctt.gaa<-((dn31.trinuc$ctt*dn31.trinuc$gaa)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.ctt.gac<-((dn31.trinuc$ctt*dn31.trinuc$gac)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld
cps.ctt.gag<-((dn31.trinuc$ctt*dn31.trinuc$gag)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.ctt.gat<-((dn31.trinuc$ctt*dn31.trinuc$gat)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld

cps.ctt.gca<-((dn31.trinuc$ctt*dn31.trinuc$gca)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctt.gcc<-((dn31.trinuc$ctt*dn31.trinuc$gcc)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctt.gcg<-((dn31.trinuc$ctt*dn31.trinuc$gcg)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ctt.gct<-((dn31.trinuc$ctt*dn31.trinuc$gct)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la

cps.ctt.gga<-((dn31.trinuc$ctt*dn31.trinuc$gga)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctt.ggc<-((dn31.trinuc$ctt*dn31.trinuc$ggc)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctt.ggg<-((dn31.trinuc$ctt*dn31.trinuc$ggg)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ctt.ggt<-((dn31.trinuc$ctt*dn31.trinuc$ggt)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg

cps.ctt.gta<-((dn31.trinuc$ctt*dn31.trinuc$gta)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctt.gtc<-((dn31.trinuc$ctt*dn31.trinuc$gtc)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctt.gtg<-((dn31.trinuc$ctt*dn31.trinuc$gtg)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ctt.gtt<-((dn31.trinuc$ctt*dn31.trinuc$gtt)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv

#Stop codon
#cps.ctt.taa<-((dn31.trinuc$ctt*dn31.trinuc$taa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctt.tac<-((dn31.trinuc$ctt*dn31.trinuc$tac)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly
#Stop codon
#cps.ctt.tag<-((dn31.trinuc$ctt*dn31.trinuc$tag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctt.tat<-((dn31.trinuc$ctt*dn31.trinuc$tat)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly

cps.ctt.tca<-((dn31.trinuc$ctt*dn31.trinuc$tca)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctt.tcc<-((dn31.trinuc$ctt*dn31.trinuc$tcc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctt.tcg<-((dn31.trinuc$ctt*dn31.trinuc$tcg)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ctt.tct<-((dn31.trinuc$ctt*dn31.trinuc$tct)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

#Stop codon
#cps.ctt.tga<-((dn31.trinuc$ctt*dn31.trinuc$tga)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ctt.tgc<-((dn31.trinuc$ctt*dn31.trinuc$tgc)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc
cps.ctt.tgg<-((dn31.trinuc$ctt*dn31.trinuc$tgg)/(dn31.count.aa$l*dn31.count.aa$w))/dn31.count.diaa$lw
cps.ctt.tgt<-((dn31.trinuc$ctt*dn31.trinuc$tgt)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc

cps.ctt.tta<-((dn31.trinuc$ctt*dn31.trinuc$tta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctt.ttc<-((dn31.trinuc$ctt*dn31.trinuc$ttc)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf
cps.ctt.ttg<-((dn31.trinuc$ctt*dn31.trinuc$ttg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ctt.ttt<-((dn31.trinuc$ctt*dn31.trinuc$ttt)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf








cps.gaa.aaa<-((dn31.trinuc$gaa*dn31.trinuc$aaa)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gaa.aac<-((dn31.trinuc$gaa*dn31.trinuc$aac)/(dn31.count.aa$e*dn31.count.aa$n))/dn31.count.diaa$en
cps.gaa.aag<-((dn31.trinuc$gaa*dn31.trinuc$aag)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gaa.aat<-((dn31.trinuc$gaa*dn31.trinuc$aat)/(dn31.count.aa$e*dn31.count.aa$n))/dn31.count.diaa$en

cps.gaa.aca<-((dn31.trinuc$gaa*dn31.trinuc$aca)/(dn31.count.aa$e*dn31.count.aa$t))/dn31.count.diaa$et
cps.gaa.acc<-((dn31.trinuc$gaa*dn31.trinuc$acc)/(dn31.count.aa$e*dn31.count.aa$t))/dn31.count.diaa$et
cps.gaa.acg<-((dn31.trinuc$gaa*dn31.trinuc$acg)/(dn31.count.aa$e*dn31.count.aa$t))/dn31.count.diaa$et
cps.gaa.act<-((dn31.trinuc$gaa*dn31.trinuc$act)/(dn31.count.aa$e*dn31.count.aa$t))/dn31.count.diaa$et

cps.gaa.aga<-((dn31.trinuc$gaa*dn31.trinuc$aga)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gaa.agc<-((dn31.trinuc$gaa*dn31.trinuc$agc)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es
cps.gaa.agg<-((dn31.trinuc$gaa*dn31.trinuc$agg)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gaa.agt<-((dn31.trinuc$gaa*dn31.trinuc$agt)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es

cps.gaa.ata<-((dn31.trinuc$gaa*dn31.trinuc$ata)/(dn31.count.aa$e*dn31.count.aa$i))/dn31.count.diaa$ei
cps.gaa.atc<-((dn31.trinuc$gaa*dn31.trinuc$atc)/(dn31.count.aa$e*dn31.count.aa$i))/dn31.count.diaa$ei
cps.gaa.atg<-((dn31.trinuc$gaa*dn31.trinuc$atg)/(dn31.count.aa$e*dn31.count.aa$m))/dn31.count.diaa$em
cps.gaa.att<-((dn31.trinuc$gaa*dn31.trinuc$att)/(dn31.count.aa$e*dn31.count.aa$i))/dn31.count.diaa$ei

cps.gaa.caa<-((dn31.trinuc$gaa*dn31.trinuc$caa)/(dn31.count.aa$e*dn31.count.aa$q))/dn31.count.diaa$eq
cps.gaa.cac<-((dn31.trinuc$gaa*dn31.trinuc$cac)/(dn31.count.aa$e*dn31.count.aa$h))/dn31.count.diaa$eh
cps.gaa.cag<-((dn31.trinuc$gaa*dn31.trinuc$cag)/(dn31.count.aa$e*dn31.count.aa$q))/dn31.count.diaa$eq
cps.gaa.cat<-((dn31.trinuc$gaa*dn31.trinuc$cat)/(dn31.count.aa$e*dn31.count.aa$h))/dn31.count.diaa$eh

cps.gaa.cca<-((dn31.trinuc$gaa*dn31.trinuc$cca)/(dn31.count.aa$e*dn31.count.aa$p))/dn31.count.diaa$ep
cps.gaa.ccc<-((dn31.trinuc$gaa*dn31.trinuc$ccc)/(dn31.count.aa$e*dn31.count.aa$p))/dn31.count.diaa$ep
cps.gaa.ccg<-((dn31.trinuc$gaa*dn31.trinuc$ccg)/(dn31.count.aa$e*dn31.count.aa$p))/dn31.count.diaa$ep
cps.gaa.cct<-((dn31.trinuc$gaa*dn31.trinuc$cct)/(dn31.count.aa$e*dn31.count.aa$p))/dn31.count.diaa$ep

cps.gaa.cga<-((dn31.trinuc$gaa*dn31.trinuc$cga)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gaa.cgc<-((dn31.trinuc$gaa*dn31.trinuc$cgc)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gaa.cgg<-((dn31.trinuc$gaa*dn31.trinuc$cgg)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gaa.cgt<-((dn31.trinuc$gaa*dn31.trinuc$cgt)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er

cps.gaa.cta<-((dn31.trinuc$gaa*dn31.trinuc$cta)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gaa.ctc<-((dn31.trinuc$gaa*dn31.trinuc$ctc)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gaa.ctg<-((dn31.trinuc$gaa*dn31.trinuc$ctg)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gaa.ctt<-((dn31.trinuc$gaa*dn31.trinuc$ctt)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el

cps.gaa.gaa<-((dn31.trinuc$gaa*dn31.trinuc$gaa)/(dn31.count.aa$e*dn31.count.aa$e))/dn31.count.diaa$ee
cps.gaa.gac<-((dn31.trinuc$gaa*dn31.trinuc$gac)/(dn31.count.aa$e*dn31.count.aa$d))/dn31.count.diaa$ed
cps.gaa.gag<-((dn31.trinuc$gaa*dn31.trinuc$gag)/(dn31.count.aa$e*dn31.count.aa$e))/dn31.count.diaa$ee
cps.gaa.gat<-((dn31.trinuc$gaa*dn31.trinuc$gat)/(dn31.count.aa$e*dn31.count.aa$d))/dn31.count.diaa$ed

cps.gaa.gca<-((dn31.trinuc$gaa*dn31.trinuc$gca)/(dn31.count.aa$e*dn31.count.aa$a))/dn31.count.diaa$ea
cps.gaa.gcc<-((dn31.trinuc$gaa*dn31.trinuc$gcc)/(dn31.count.aa$e*dn31.count.aa$a))/dn31.count.diaa$ea
cps.gaa.gcg<-((dn31.trinuc$gaa*dn31.trinuc$gcg)/(dn31.count.aa$e*dn31.count.aa$a))/dn31.count.diaa$ea
cps.gaa.gct<-((dn31.trinuc$gaa*dn31.trinuc$gct)/(dn31.count.aa$e*dn31.count.aa$a))/dn31.count.diaa$ea

cps.gaa.gga<-((dn31.trinuc$gaa*dn31.trinuc$gga)/(dn31.count.aa$e*dn31.count.aa$g))/dn31.count.diaa$eg
cps.gaa.ggc<-((dn31.trinuc$gaa*dn31.trinuc$ggc)/(dn31.count.aa$e*dn31.count.aa$g))/dn31.count.diaa$eg
cps.gaa.ggg<-((dn31.trinuc$gaa*dn31.trinuc$ggg)/(dn31.count.aa$e*dn31.count.aa$g))/dn31.count.diaa$eg
cps.gaa.ggt<-((dn31.trinuc$gaa*dn31.trinuc$ggt)/(dn31.count.aa$e*dn31.count.aa$g))/dn31.count.diaa$eg

cps.gaa.gta<-((dn31.trinuc$gaa*dn31.trinuc$gta)/(dn31.count.aa$e*dn31.count.aa$v))/dn31.count.diaa$ev
cps.gaa.gtc<-((dn31.trinuc$gaa*dn31.trinuc$gtc)/(dn31.count.aa$e*dn31.count.aa$v))/dn31.count.diaa$ev
cps.gaa.gtg<-((dn31.trinuc$gaa*dn31.trinuc$gtg)/(dn31.count.aa$e*dn31.count.aa$v))/dn31.count.diaa$ev
cps.gaa.gtt<-((dn31.trinuc$gaa*dn31.trinuc$gtt)/(dn31.count.aa$e*dn31.count.aa$v))/dn31.count.diaa$ev

#Stop codon
#cps.gaa.taa<-((dn31.trinuc$gaa*dn31.trinuc$taa)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gaa.tac<-((dn31.trinuc$gaa*dn31.trinuc$tac)/(dn31.count.aa$e*dn31.count.aa$y))/dn31.count.diaa$ey
#Stop codon
#cps.gaa.tag<-((dn31.trinuc$gaa*dn31.trinuc$tag)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gaa.tat<-((dn31.trinuc$gaa*dn31.trinuc$tat)/(dn31.count.aa$e*dn31.count.aa$y))/dn31.count.diaa$ey

cps.gaa.tca<-((dn31.trinuc$gaa*dn31.trinuc$tca)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es
cps.gaa.tcc<-((dn31.trinuc$gaa*dn31.trinuc$tcc)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es
cps.gaa.tcg<-((dn31.trinuc$gaa*dn31.trinuc$tcg)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es
cps.gaa.tct<-((dn31.trinuc$gaa*dn31.trinuc$tct)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es

#Stop codon
#cps.gaa.tga<-((dn31.trinuc$gaa*dn31.trinuc$tga)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gaa.tgc<-((dn31.trinuc$gaa*dn31.trinuc$tgc)/(dn31.count.aa$e*dn31.count.aa$c))/dn31.count.diaa$ec
cps.gaa.tgg<-((dn31.trinuc$gaa*dn31.trinuc$tgg)/(dn31.count.aa$e*dn31.count.aa$w))/dn31.count.diaa$ew
cps.gaa.tgt<-((dn31.trinuc$gaa*dn31.trinuc$tgt)/(dn31.count.aa$e*dn31.count.aa$c))/dn31.count.diaa$ec

cps.gaa.tta<-((dn31.trinuc$gaa*dn31.trinuc$tta)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gaa.ttc<-((dn31.trinuc$gaa*dn31.trinuc$ttc)/(dn31.count.aa$e*dn31.count.aa$f))/dn31.count.diaa$ef
cps.gaa.ttg<-((dn31.trinuc$gaa*dn31.trinuc$ttg)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gaa.ttt<-((dn31.trinuc$gaa*dn31.trinuc$ttt)/(dn31.count.aa$e*dn31.count.aa$f))/dn31.count.diaa$ef








cps.gac.aaa<-((dn31.trinuc$gac*dn31.trinuc$aaa)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gac.aac<-((dn31.trinuc$gac*dn31.trinuc$aac)/(dn31.count.aa$d*dn31.count.aa$n))/dn31.count.diaa$dn
cps.gac.aag<-((dn31.trinuc$gac*dn31.trinuc$aag)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gac.aat<-((dn31.trinuc$gac*dn31.trinuc$aat)/(dn31.count.aa$d*dn31.count.aa$n))/dn31.count.diaa$dn

cps.gac.aca<-((dn31.trinuc$gac*dn31.trinuc$aca)/(dn31.count.aa$d*dn31.count.aa$t))/dn31.count.diaa$dt
cps.gac.acc<-((dn31.trinuc$gac*dn31.trinuc$acc)/(dn31.count.aa$d*dn31.count.aa$t))/dn31.count.diaa$dt
cps.gac.acg<-((dn31.trinuc$gac*dn31.trinuc$acg)/(dn31.count.aa$d*dn31.count.aa$t))/dn31.count.diaa$dt
cps.gac.act<-((dn31.trinuc$gac*dn31.trinuc$act)/(dn31.count.aa$d*dn31.count.aa$t))/dn31.count.diaa$dt

cps.gac.aga<-((dn31.trinuc$gac*dn31.trinuc$aga)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gac.agc<-((dn31.trinuc$gac*dn31.trinuc$agc)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds
cps.gac.agg<-((dn31.trinuc$gac*dn31.trinuc$agg)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gac.agt<-((dn31.trinuc$gac*dn31.trinuc$agt)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds

cps.gac.ata<-((dn31.trinuc$gac*dn31.trinuc$ata)/(dn31.count.aa$d*dn31.count.aa$i))/dn31.count.diaa$di
cps.gac.atc<-((dn31.trinuc$gac*dn31.trinuc$atc)/(dn31.count.aa$d*dn31.count.aa$i))/dn31.count.diaa$di
cps.gac.atg<-((dn31.trinuc$gac*dn31.trinuc$atg)/(dn31.count.aa$d*dn31.count.aa$m))/dn31.count.diaa$dm
cps.gac.att<-((dn31.trinuc$gac*dn31.trinuc$att)/(dn31.count.aa$d*dn31.count.aa$i))/dn31.count.diaa$di

cps.gac.caa<-((dn31.trinuc$gac*dn31.trinuc$caa)/(dn31.count.aa$d*dn31.count.aa$q))/dn31.count.diaa$dq
cps.gac.cac<-((dn31.trinuc$gac*dn31.trinuc$cac)/(dn31.count.aa$d*dn31.count.aa$h))/dn31.count.diaa$dh
cps.gac.cag<-((dn31.trinuc$gac*dn31.trinuc$cag)/(dn31.count.aa$d*dn31.count.aa$q))/dn31.count.diaa$dq
cps.gac.cat<-((dn31.trinuc$gac*dn31.trinuc$cat)/(dn31.count.aa$d*dn31.count.aa$h))/dn31.count.diaa$dh

cps.gac.cca<-((dn31.trinuc$gac*dn31.trinuc$cca)/(dn31.count.aa$d*dn31.count.aa$p))/dn31.count.diaa$dp
cps.gac.ccc<-((dn31.trinuc$gac*dn31.trinuc$ccc)/(dn31.count.aa$d*dn31.count.aa$p))/dn31.count.diaa$dp
cps.gac.ccg<-((dn31.trinuc$gac*dn31.trinuc$ccg)/(dn31.count.aa$d*dn31.count.aa$p))/dn31.count.diaa$dp
cps.gac.cct<-((dn31.trinuc$gac*dn31.trinuc$cct)/(dn31.count.aa$d*dn31.count.aa$p))/dn31.count.diaa$dp

cps.gac.cga<-((dn31.trinuc$gac*dn31.trinuc$cga)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gac.cgc<-((dn31.trinuc$gac*dn31.trinuc$cgc)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gac.cgg<-((dn31.trinuc$gac*dn31.trinuc$cgg)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gac.cgt<-((dn31.trinuc$gac*dn31.trinuc$cgt)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr

cps.gac.cta<-((dn31.trinuc$gac*dn31.trinuc$cta)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gac.ctc<-((dn31.trinuc$gac*dn31.trinuc$ctc)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gac.ctg<-((dn31.trinuc$gac*dn31.trinuc$ctg)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gac.ctt<-((dn31.trinuc$gac*dn31.trinuc$ctt)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl

cps.gac.gaa<-((dn31.trinuc$gac*dn31.trinuc$gaa)/(dn31.count.aa$d*dn31.count.aa$e))/dn31.count.diaa$de
cps.gac.gac<-((dn31.trinuc$gac*dn31.trinuc$gac)/(dn31.count.aa$d*dn31.count.aa$d))/dn31.count.diaa$dd
cps.gac.gag<-((dn31.trinuc$gac*dn31.trinuc$gag)/(dn31.count.aa$d*dn31.count.aa$e))/dn31.count.diaa$de
cps.gac.gat<-((dn31.trinuc$gac*dn31.trinuc$gat)/(dn31.count.aa$d*dn31.count.aa$d))/dn31.count.diaa$dd

cps.gac.gca<-((dn31.trinuc$gac*dn31.trinuc$gca)/(dn31.count.aa$d*dn31.count.aa$a))/dn31.count.diaa$da
cps.gac.gcc<-((dn31.trinuc$gac*dn31.trinuc$gcc)/(dn31.count.aa$d*dn31.count.aa$a))/dn31.count.diaa$da
cps.gac.gcg<-((dn31.trinuc$gac*dn31.trinuc$gcg)/(dn31.count.aa$d*dn31.count.aa$a))/dn31.count.diaa$da
cps.gac.gct<-((dn31.trinuc$gac*dn31.trinuc$gct)/(dn31.count.aa$d*dn31.count.aa$a))/dn31.count.diaa$da

cps.gac.gga<-((dn31.trinuc$gac*dn31.trinuc$gga)/(dn31.count.aa$d*dn31.count.aa$g))/dn31.count.diaa$dg
cps.gac.ggc<-((dn31.trinuc$gac*dn31.trinuc$ggc)/(dn31.count.aa$d*dn31.count.aa$g))/dn31.count.diaa$dg
cps.gac.ggg<-((dn31.trinuc$gac*dn31.trinuc$ggg)/(dn31.count.aa$d*dn31.count.aa$g))/dn31.count.diaa$dg
cps.gac.ggt<-((dn31.trinuc$gac*dn31.trinuc$ggt)/(dn31.count.aa$d*dn31.count.aa$g))/dn31.count.diaa$dg

cps.gac.gta<-((dn31.trinuc$gac*dn31.trinuc$gta)/(dn31.count.aa$d*dn31.count.aa$v))/dn31.count.diaa$dv
cps.gac.gtc<-((dn31.trinuc$gac*dn31.trinuc$gtc)/(dn31.count.aa$d*dn31.count.aa$v))/dn31.count.diaa$dv
cps.gac.gtg<-((dn31.trinuc$gac*dn31.trinuc$gtg)/(dn31.count.aa$d*dn31.count.aa$v))/dn31.count.diaa$dv
cps.gac.gtt<-((dn31.trinuc$gac*dn31.trinuc$gtt)/(dn31.count.aa$d*dn31.count.aa$v))/dn31.count.diaa$dv

#Stop codon
#cps.gac.taa<-((dn31.trinuc$gac*dn31.trinuc$taa)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gac.tac<-((dn31.trinuc$gac*dn31.trinuc$tac)/(dn31.count.aa$d*dn31.count.aa$y))/dn31.count.diaa$dy
#Stop codon
#cps.gac.tag<-((dn31.trinuc$gac*dn31.trinuc$tag)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gac.tat<-((dn31.trinuc$gac*dn31.trinuc$tat)/(dn31.count.aa$d*dn31.count.aa$y))/dn31.count.diaa$dy

cps.gac.tca<-((dn31.trinuc$gac*dn31.trinuc$tca)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds
cps.gac.tcc<-((dn31.trinuc$gac*dn31.trinuc$tcc)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds
cps.gac.tcg<-((dn31.trinuc$gac*dn31.trinuc$tcg)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds
cps.gac.tct<-((dn31.trinuc$gac*dn31.trinuc$tct)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds

#Stop codon
#cps.gac.tga<-((dn31.trinuc$gac*dn31.trinuc$tga)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gac.tgc<-((dn31.trinuc$gac*dn31.trinuc$tgc)/(dn31.count.aa$d*dn31.count.aa$c))/dn31.count.diaa$dc
cps.gac.tgg<-((dn31.trinuc$gac*dn31.trinuc$tgg)/(dn31.count.aa$d*dn31.count.aa$w))/dn31.count.diaa$dw
cps.gac.tgt<-((dn31.trinuc$gac*dn31.trinuc$tgt)/(dn31.count.aa$d*dn31.count.aa$c))/dn31.count.diaa$dc

cps.gac.tta<-((dn31.trinuc$gac*dn31.trinuc$tta)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gac.ttc<-((dn31.trinuc$gac*dn31.trinuc$ttc)/(dn31.count.aa$d*dn31.count.aa$f))/dn31.count.diaa$df
cps.gac.ttg<-((dn31.trinuc$gac*dn31.trinuc$ttg)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gac.ttt<-((dn31.trinuc$gac*dn31.trinuc$ttt)/(dn31.count.aa$d*dn31.count.aa$f))/dn31.count.diaa$df








cps.gag.aaa<-((dn31.trinuc$gag*dn31.trinuc$aaa)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gag.aac<-((dn31.trinuc$gag*dn31.trinuc$aac)/(dn31.count.aa$e*dn31.count.aa$n))/dn31.count.diaa$en
cps.gag.aag<-((dn31.trinuc$gag*dn31.trinuc$aag)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gag.aat<-((dn31.trinuc$gag*dn31.trinuc$aat)/(dn31.count.aa$e*dn31.count.aa$n))/dn31.count.diaa$en

cps.gag.aca<-((dn31.trinuc$gag*dn31.trinuc$aca)/(dn31.count.aa$e*dn31.count.aa$t))/dn31.count.diaa$et
cps.gag.acc<-((dn31.trinuc$gag*dn31.trinuc$acc)/(dn31.count.aa$e*dn31.count.aa$t))/dn31.count.diaa$et
cps.gag.acg<-((dn31.trinuc$gag*dn31.trinuc$acg)/(dn31.count.aa$e*dn31.count.aa$t))/dn31.count.diaa$et
cps.gag.act<-((dn31.trinuc$gag*dn31.trinuc$act)/(dn31.count.aa$e*dn31.count.aa$t))/dn31.count.diaa$et

cps.gag.aga<-((dn31.trinuc$gag*dn31.trinuc$aga)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gag.agc<-((dn31.trinuc$gag*dn31.trinuc$agc)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es
cps.gag.agg<-((dn31.trinuc$gag*dn31.trinuc$agg)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gag.agt<-((dn31.trinuc$gag*dn31.trinuc$agt)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es

cps.gag.ata<-((dn31.trinuc$gag*dn31.trinuc$ata)/(dn31.count.aa$e*dn31.count.aa$i))/dn31.count.diaa$ei
cps.gag.atc<-((dn31.trinuc$gag*dn31.trinuc$atc)/(dn31.count.aa$e*dn31.count.aa$i))/dn31.count.diaa$ei
cps.gag.atg<-((dn31.trinuc$gag*dn31.trinuc$atg)/(dn31.count.aa$e*dn31.count.aa$m))/dn31.count.diaa$em
cps.gag.att<-((dn31.trinuc$gag*dn31.trinuc$att)/(dn31.count.aa$e*dn31.count.aa$i))/dn31.count.diaa$ei

cps.gag.caa<-((dn31.trinuc$gag*dn31.trinuc$caa)/(dn31.count.aa$e*dn31.count.aa$q))/dn31.count.diaa$eq
cps.gag.cac<-((dn31.trinuc$gag*dn31.trinuc$cac)/(dn31.count.aa$e*dn31.count.aa$h))/dn31.count.diaa$eh
cps.gag.cag<-((dn31.trinuc$gag*dn31.trinuc$cag)/(dn31.count.aa$e*dn31.count.aa$q))/dn31.count.diaa$eq
cps.gag.cat<-((dn31.trinuc$gag*dn31.trinuc$cat)/(dn31.count.aa$e*dn31.count.aa$h))/dn31.count.diaa$eh

cps.gag.cca<-((dn31.trinuc$gag*dn31.trinuc$cca)/(dn31.count.aa$e*dn31.count.aa$p))/dn31.count.diaa$ep
cps.gag.ccc<-((dn31.trinuc$gag*dn31.trinuc$ccc)/(dn31.count.aa$e*dn31.count.aa$p))/dn31.count.diaa$ep
cps.gag.ccg<-((dn31.trinuc$gag*dn31.trinuc$ccg)/(dn31.count.aa$e*dn31.count.aa$p))/dn31.count.diaa$ep
cps.gag.cct<-((dn31.trinuc$gag*dn31.trinuc$cct)/(dn31.count.aa$e*dn31.count.aa$p))/dn31.count.diaa$ep

cps.gag.cga<-((dn31.trinuc$gag*dn31.trinuc$cga)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gag.cgc<-((dn31.trinuc$gag*dn31.trinuc$cgc)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gag.cgg<-((dn31.trinuc$gag*dn31.trinuc$cgg)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er
cps.gag.cgt<-((dn31.trinuc$gag*dn31.trinuc$cgt)/(dn31.count.aa$e*dn31.count.aa$r))/dn31.count.diaa$er

cps.gag.cta<-((dn31.trinuc$gag*dn31.trinuc$cta)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gag.ctc<-((dn31.trinuc$gag*dn31.trinuc$ctc)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gag.ctg<-((dn31.trinuc$gag*dn31.trinuc$ctg)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gag.ctt<-((dn31.trinuc$gag*dn31.trinuc$ctt)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el

cps.gag.gaa<-((dn31.trinuc$gag*dn31.trinuc$gaa)/(dn31.count.aa$e*dn31.count.aa$e))/dn31.count.diaa$ee
cps.gag.gac<-((dn31.trinuc$gag*dn31.trinuc$gac)/(dn31.count.aa$e*dn31.count.aa$d))/dn31.count.diaa$ed
cps.gag.gag<-((dn31.trinuc$gag*dn31.trinuc$gag)/(dn31.count.aa$e*dn31.count.aa$e))/dn31.count.diaa$ee
cps.gag.gat<-((dn31.trinuc$gag*dn31.trinuc$gat)/(dn31.count.aa$e*dn31.count.aa$d))/dn31.count.diaa$ed

cps.gag.gca<-((dn31.trinuc$gag*dn31.trinuc$gca)/(dn31.count.aa$e*dn31.count.aa$a))/dn31.count.diaa$ea
cps.gag.gcc<-((dn31.trinuc$gag*dn31.trinuc$gcc)/(dn31.count.aa$e*dn31.count.aa$a))/dn31.count.diaa$ea
cps.gag.gcg<-((dn31.trinuc$gag*dn31.trinuc$gcg)/(dn31.count.aa$e*dn31.count.aa$a))/dn31.count.diaa$ea
cps.gag.gct<-((dn31.trinuc$gag*dn31.trinuc$gct)/(dn31.count.aa$e*dn31.count.aa$a))/dn31.count.diaa$ea

cps.gag.gga<-((dn31.trinuc$gag*dn31.trinuc$gga)/(dn31.count.aa$e*dn31.count.aa$g))/dn31.count.diaa$eg
cps.gag.ggc<-((dn31.trinuc$gag*dn31.trinuc$ggc)/(dn31.count.aa$e*dn31.count.aa$g))/dn31.count.diaa$eg
cps.gag.ggg<-((dn31.trinuc$gag*dn31.trinuc$ggg)/(dn31.count.aa$e*dn31.count.aa$g))/dn31.count.diaa$eg
cps.gag.ggt<-((dn31.trinuc$gag*dn31.trinuc$ggt)/(dn31.count.aa$e*dn31.count.aa$g))/dn31.count.diaa$eg

cps.gag.gta<-((dn31.trinuc$gag*dn31.trinuc$gta)/(dn31.count.aa$e*dn31.count.aa$v))/dn31.count.diaa$ev
cps.gag.gtc<-((dn31.trinuc$gag*dn31.trinuc$gtc)/(dn31.count.aa$e*dn31.count.aa$v))/dn31.count.diaa$ev
cps.gag.gtg<-((dn31.trinuc$gag*dn31.trinuc$gtg)/(dn31.count.aa$e*dn31.count.aa$v))/dn31.count.diaa$ev
cps.gag.gtt<-((dn31.trinuc$gag*dn31.trinuc$gtt)/(dn31.count.aa$e*dn31.count.aa$v))/dn31.count.diaa$ev

#Stop codon
#cps.gag.taa<-((dn31.trinuc$gag*dn31.trinuc$taa)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gag.tac<-((dn31.trinuc$gag*dn31.trinuc$tac)/(dn31.count.aa$e*dn31.count.aa$y))/dn31.count.diaa$ey
#Stop codon
#cps.gag.tag<-((dn31.trinuc$gag*dn31.trinuc$tag)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gag.tat<-((dn31.trinuc$gag*dn31.trinuc$tat)/(dn31.count.aa$e*dn31.count.aa$y))/dn31.count.diaa$ey

cps.gag.tca<-((dn31.trinuc$gag*dn31.trinuc$tca)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es
cps.gag.tcc<-((dn31.trinuc$gag*dn31.trinuc$tcc)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es
cps.gag.tcg<-((dn31.trinuc$gag*dn31.trinuc$tcg)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es
cps.gag.tct<-((dn31.trinuc$gag*dn31.trinuc$tct)/(dn31.count.aa$e*dn31.count.aa$s))/dn31.count.diaa$es

#Stop codon
#cps.gag.tga<-((dn31.trinuc$gag*dn31.trinuc$tga)/(dn31.count.aa$e*dn31.count.aa$k))/dn31.count.diaa$ek
cps.gag.tgc<-((dn31.trinuc$gag*dn31.trinuc$tgc)/(dn31.count.aa$e*dn31.count.aa$c))/dn31.count.diaa$ec
cps.gag.tgg<-((dn31.trinuc$gag*dn31.trinuc$tgg)/(dn31.count.aa$e*dn31.count.aa$w))/dn31.count.diaa$ew
cps.gag.tgt<-((dn31.trinuc$gag*dn31.trinuc$tgt)/(dn31.count.aa$e*dn31.count.aa$c))/dn31.count.diaa$ec

cps.gag.tta<-((dn31.trinuc$gag*dn31.trinuc$tta)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gag.ttc<-((dn31.trinuc$gag*dn31.trinuc$ttc)/(dn31.count.aa$e*dn31.count.aa$f))/dn31.count.diaa$ef
cps.gag.ttg<-((dn31.trinuc$gag*dn31.trinuc$ttg)/(dn31.count.aa$e*dn31.count.aa$l))/dn31.count.diaa$el
cps.gag.ttt<-((dn31.trinuc$gag*dn31.trinuc$ttt)/(dn31.count.aa$e*dn31.count.aa$f))/dn31.count.diaa$ef








cps.gat.aaa<-((dn31.trinuc$gat*dn31.trinuc$aaa)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gat.aac<-((dn31.trinuc$gat*dn31.trinuc$aac)/(dn31.count.aa$d*dn31.count.aa$n))/dn31.count.diaa$dn
cps.gat.aag<-((dn31.trinuc$gat*dn31.trinuc$aag)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gat.aat<-((dn31.trinuc$gat*dn31.trinuc$aat)/(dn31.count.aa$d*dn31.count.aa$n))/dn31.count.diaa$dn

cps.gat.aca<-((dn31.trinuc$gat*dn31.trinuc$aca)/(dn31.count.aa$d*dn31.count.aa$t))/dn31.count.diaa$dt
cps.gat.acc<-((dn31.trinuc$gat*dn31.trinuc$acc)/(dn31.count.aa$d*dn31.count.aa$t))/dn31.count.diaa$dt
cps.gat.acg<-((dn31.trinuc$gat*dn31.trinuc$acg)/(dn31.count.aa$d*dn31.count.aa$t))/dn31.count.diaa$dt
cps.gat.act<-((dn31.trinuc$gat*dn31.trinuc$act)/(dn31.count.aa$d*dn31.count.aa$t))/dn31.count.diaa$dt

cps.gat.aga<-((dn31.trinuc$gat*dn31.trinuc$aga)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gat.agc<-((dn31.trinuc$gat*dn31.trinuc$agc)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds
cps.gat.agg<-((dn31.trinuc$gat*dn31.trinuc$agg)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gat.agt<-((dn31.trinuc$gat*dn31.trinuc$agt)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds

cps.gat.ata<-((dn31.trinuc$gat*dn31.trinuc$ata)/(dn31.count.aa$d*dn31.count.aa$i))/dn31.count.diaa$di
cps.gat.atc<-((dn31.trinuc$gat*dn31.trinuc$atc)/(dn31.count.aa$d*dn31.count.aa$i))/dn31.count.diaa$di
cps.gat.atg<-((dn31.trinuc$gat*dn31.trinuc$atg)/(dn31.count.aa$d*dn31.count.aa$m))/dn31.count.diaa$dm
cps.gat.att<-((dn31.trinuc$gat*dn31.trinuc$att)/(dn31.count.aa$d*dn31.count.aa$i))/dn31.count.diaa$di

cps.gat.caa<-((dn31.trinuc$gat*dn31.trinuc$caa)/(dn31.count.aa$d*dn31.count.aa$q))/dn31.count.diaa$dq
cps.gat.cac<-((dn31.trinuc$gat*dn31.trinuc$cac)/(dn31.count.aa$d*dn31.count.aa$h))/dn31.count.diaa$dh
cps.gat.cag<-((dn31.trinuc$gat*dn31.trinuc$cag)/(dn31.count.aa$d*dn31.count.aa$q))/dn31.count.diaa$dq
cps.gat.cat<-((dn31.trinuc$gat*dn31.trinuc$cat)/(dn31.count.aa$d*dn31.count.aa$h))/dn31.count.diaa$dh

cps.gat.cca<-((dn31.trinuc$gat*dn31.trinuc$cca)/(dn31.count.aa$d*dn31.count.aa$p))/dn31.count.diaa$dp
cps.gat.ccc<-((dn31.trinuc$gat*dn31.trinuc$ccc)/(dn31.count.aa$d*dn31.count.aa$p))/dn31.count.diaa$dp
cps.gat.ccg<-((dn31.trinuc$gat*dn31.trinuc$ccg)/(dn31.count.aa$d*dn31.count.aa$p))/dn31.count.diaa$dp
cps.gat.cct<-((dn31.trinuc$gat*dn31.trinuc$cct)/(dn31.count.aa$d*dn31.count.aa$p))/dn31.count.diaa$dp

cps.gat.cga<-((dn31.trinuc$gat*dn31.trinuc$cga)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gat.cgc<-((dn31.trinuc$gat*dn31.trinuc$cgc)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gat.cgg<-((dn31.trinuc$gat*dn31.trinuc$cgg)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr
cps.gat.cgt<-((dn31.trinuc$gat*dn31.trinuc$cgt)/(dn31.count.aa$d*dn31.count.aa$r))/dn31.count.diaa$dr

cps.gat.cta<-((dn31.trinuc$gat*dn31.trinuc$cta)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gat.ctc<-((dn31.trinuc$gat*dn31.trinuc$ctc)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gat.ctg<-((dn31.trinuc$gat*dn31.trinuc$ctg)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gat.ctt<-((dn31.trinuc$gat*dn31.trinuc$ctt)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl

cps.gat.gaa<-((dn31.trinuc$gat*dn31.trinuc$gaa)/(dn31.count.aa$d*dn31.count.aa$e))/dn31.count.diaa$de
cps.gat.gac<-((dn31.trinuc$gat*dn31.trinuc$gac)/(dn31.count.aa$d*dn31.count.aa$d))/dn31.count.diaa$dd
cps.gat.gag<-((dn31.trinuc$gat*dn31.trinuc$gag)/(dn31.count.aa$d*dn31.count.aa$e))/dn31.count.diaa$de
cps.gat.gat<-((dn31.trinuc$gat*dn31.trinuc$gat)/(dn31.count.aa$d*dn31.count.aa$d))/dn31.count.diaa$dd

cps.gat.gca<-((dn31.trinuc$gat*dn31.trinuc$gca)/(dn31.count.aa$d*dn31.count.aa$a))/dn31.count.diaa$da
cps.gat.gcc<-((dn31.trinuc$gat*dn31.trinuc$gcc)/(dn31.count.aa$d*dn31.count.aa$a))/dn31.count.diaa$da
cps.gat.gcg<-((dn31.trinuc$gat*dn31.trinuc$gcg)/(dn31.count.aa$d*dn31.count.aa$a))/dn31.count.diaa$da
cps.gat.gct<-((dn31.trinuc$gat*dn31.trinuc$gct)/(dn31.count.aa$d*dn31.count.aa$a))/dn31.count.diaa$da

cps.gat.gga<-((dn31.trinuc$gat*dn31.trinuc$gga)/(dn31.count.aa$d*dn31.count.aa$g))/dn31.count.diaa$dg
cps.gat.ggc<-((dn31.trinuc$gat*dn31.trinuc$ggc)/(dn31.count.aa$d*dn31.count.aa$g))/dn31.count.diaa$dg
cps.gat.ggg<-((dn31.trinuc$gat*dn31.trinuc$ggg)/(dn31.count.aa$d*dn31.count.aa$g))/dn31.count.diaa$dg
cps.gat.ggt<-((dn31.trinuc$gat*dn31.trinuc$ggt)/(dn31.count.aa$d*dn31.count.aa$g))/dn31.count.diaa$dg

cps.gat.gta<-((dn31.trinuc$gat*dn31.trinuc$gta)/(dn31.count.aa$d*dn31.count.aa$v))/dn31.count.diaa$dv
cps.gat.gtc<-((dn31.trinuc$gat*dn31.trinuc$gtc)/(dn31.count.aa$d*dn31.count.aa$v))/dn31.count.diaa$dv
cps.gat.gtg<-((dn31.trinuc$gat*dn31.trinuc$gtg)/(dn31.count.aa$d*dn31.count.aa$v))/dn31.count.diaa$dv
cps.gat.gtt<-((dn31.trinuc$gat*dn31.trinuc$gtt)/(dn31.count.aa$d*dn31.count.aa$v))/dn31.count.diaa$dv

#Stop codon
#cps.gat.taa<-((dn31.trinuc$gat*dn31.trinuc$taa)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gat.tac<-((dn31.trinuc$gat*dn31.trinuc$tac)/(dn31.count.aa$d*dn31.count.aa$y))/dn31.count.diaa$dy
#Stop codon
#cps.gat.tag<-((dn31.trinuc$gat*dn31.trinuc$tag)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gat.tat<-((dn31.trinuc$gat*dn31.trinuc$tat)/(dn31.count.aa$d*dn31.count.aa$y))/dn31.count.diaa$dy

cps.gat.tca<-((dn31.trinuc$gat*dn31.trinuc$tca)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds
cps.gat.tcc<-((dn31.trinuc$gat*dn31.trinuc$tcc)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds
cps.gat.tcg<-((dn31.trinuc$gat*dn31.trinuc$tcg)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds
cps.gat.tct<-((dn31.trinuc$gat*dn31.trinuc$tct)/(dn31.count.aa$d*dn31.count.aa$s))/dn31.count.diaa$ds

#Stop codon
#cps.gat.tga<-((dn31.trinuc$gat*dn31.trinuc$tga)/(dn31.count.aa$d*dn31.count.aa$k))/dn31.count.diaa$dk
cps.gat.tgc<-((dn31.trinuc$gat*dn31.trinuc$tgc)/(dn31.count.aa$d*dn31.count.aa$c))/dn31.count.diaa$dc
cps.gat.tgg<-((dn31.trinuc$gat*dn31.trinuc$tgg)/(dn31.count.aa$d*dn31.count.aa$w))/dn31.count.diaa$dw
cps.gat.tgt<-((dn31.trinuc$gat*dn31.trinuc$tgt)/(dn31.count.aa$d*dn31.count.aa$c))/dn31.count.diaa$dc

cps.gat.tta<-((dn31.trinuc$gat*dn31.trinuc$tta)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gat.ttc<-((dn31.trinuc$gat*dn31.trinuc$ttc)/(dn31.count.aa$d*dn31.count.aa$f))/dn31.count.diaa$df
cps.gat.ttg<-((dn31.trinuc$gat*dn31.trinuc$ttg)/(dn31.count.aa$d*dn31.count.aa$l))/dn31.count.diaa$dl
cps.gat.ttt<-((dn31.trinuc$gat*dn31.trinuc$ttt)/(dn31.count.aa$d*dn31.count.aa$f))/dn31.count.diaa$df








cps.gca.aaa<-((dn31.trinuc$gca*dn31.trinuc$aaa)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gca.aac<-((dn31.trinuc$gca*dn31.trinuc$aac)/(dn31.count.aa$a*dn31.count.aa$n))/dn31.count.diaa$an
cps.gca.aag<-((dn31.trinuc$gca*dn31.trinuc$aag)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gca.aat<-((dn31.trinuc$gca*dn31.trinuc$aat)/(dn31.count.aa$a*dn31.count.aa$n))/dn31.count.diaa$an

cps.gca.aca<-((dn31.trinuc$gca*dn31.trinuc$aca)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gca.acc<-((dn31.trinuc$gca*dn31.trinuc$acc)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gca.acg<-((dn31.trinuc$gca*dn31.trinuc$acg)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gca.act<-((dn31.trinuc$gca*dn31.trinuc$act)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at

cps.gca.aga<-((dn31.trinuc$gca*dn31.trinuc$aga)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gca.agc<-((dn31.trinuc$gca*dn31.trinuc$agc)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gca.agg<-((dn31.trinuc$gca*dn31.trinuc$agg)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gca.agt<-((dn31.trinuc$gca*dn31.trinuc$agt)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as

cps.gca.ata<-((dn31.trinuc$gca*dn31.trinuc$ata)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai
cps.gca.atc<-((dn31.trinuc$gca*dn31.trinuc$atc)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai
cps.gca.atg<-((dn31.trinuc$gca*dn31.trinuc$atg)/(dn31.count.aa$a*dn31.count.aa$m))/dn31.count.diaa$am
cps.gca.att<-((dn31.trinuc$gca*dn31.trinuc$att)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai

cps.gca.caa<-((dn31.trinuc$gca*dn31.trinuc$caa)/(dn31.count.aa$a*dn31.count.aa$q))/dn31.count.diaa$aq
cps.gca.cac<-((dn31.trinuc$gca*dn31.trinuc$cac)/(dn31.count.aa$a*dn31.count.aa$h))/dn31.count.diaa$ah
cps.gca.cag<-((dn31.trinuc$gca*dn31.trinuc$cag)/(dn31.count.aa$a*dn31.count.aa$q))/dn31.count.diaa$aq
cps.gca.cat<-((dn31.trinuc$gca*dn31.trinuc$cat)/(dn31.count.aa$a*dn31.count.aa$h))/dn31.count.diaa$ah

cps.gca.cca<-((dn31.trinuc$gca*dn31.trinuc$cca)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gca.ccc<-((dn31.trinuc$gca*dn31.trinuc$ccc)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gca.ccg<-((dn31.trinuc$gca*dn31.trinuc$ccg)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gca.cct<-((dn31.trinuc$gca*dn31.trinuc$cct)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap

cps.gca.cga<-((dn31.trinuc$gca*dn31.trinuc$cga)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gca.cgc<-((dn31.trinuc$gca*dn31.trinuc$cgc)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gca.cgg<-((dn31.trinuc$gca*dn31.trinuc$cgg)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gca.cgt<-((dn31.trinuc$gca*dn31.trinuc$cgt)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar

cps.gca.cta<-((dn31.trinuc$gca*dn31.trinuc$cta)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gca.ctc<-((dn31.trinuc$gca*dn31.trinuc$ctc)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gca.ctg<-((dn31.trinuc$gca*dn31.trinuc$ctg)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gca.ctt<-((dn31.trinuc$gca*dn31.trinuc$ctt)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al

cps.gca.gaa<-((dn31.trinuc$gca*dn31.trinuc$gaa)/(dn31.count.aa$a*dn31.count.aa$e))/dn31.count.diaa$ae
cps.gca.gac<-((dn31.trinuc$gca*dn31.trinuc$gac)/(dn31.count.aa$a*dn31.count.aa$d))/dn31.count.diaa$ad
cps.gca.gag<-((dn31.trinuc$gca*dn31.trinuc$gag)/(dn31.count.aa$a*dn31.count.aa$e))/dn31.count.diaa$ae
cps.gca.gat<-((dn31.trinuc$gca*dn31.trinuc$gat)/(dn31.count.aa$a*dn31.count.aa$d))/dn31.count.diaa$ad

cps.gca.gca<-((dn31.trinuc$gca*dn31.trinuc$gca)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gca.gcc<-((dn31.trinuc$gca*dn31.trinuc$gcc)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gca.gcg<-((dn31.trinuc$gca*dn31.trinuc$gcg)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gca.gct<-((dn31.trinuc$gca*dn31.trinuc$gct)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa

cps.gca.gga<-((dn31.trinuc$gca*dn31.trinuc$gga)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gca.ggc<-((dn31.trinuc$gca*dn31.trinuc$ggc)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gca.ggg<-((dn31.trinuc$gca*dn31.trinuc$ggg)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gca.ggt<-((dn31.trinuc$gca*dn31.trinuc$ggt)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag

cps.gca.gta<-((dn31.trinuc$gca*dn31.trinuc$gta)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gca.gtc<-((dn31.trinuc$gca*dn31.trinuc$gtc)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gca.gtg<-((dn31.trinuc$gca*dn31.trinuc$gtg)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gca.gtt<-((dn31.trinuc$gca*dn31.trinuc$gtt)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av

#Stop codon
#cps.gca.taa<-((dn31.trinuc$gca*dn31.trinuc$taa)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gca.tac<-((dn31.trinuc$gca*dn31.trinuc$tac)/(dn31.count.aa$a*dn31.count.aa$y))/dn31.count.diaa$ay
#Stop codon
#cps.gca.tag<-((dn31.trinuc$gca*dn31.trinuc$tag)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gca.tat<-((dn31.trinuc$gca*dn31.trinuc$tat)/(dn31.count.aa$a*dn31.count.aa$y))/dn31.count.diaa$ay

cps.gca.tca<-((dn31.trinuc$gca*dn31.trinuc$tca)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gca.tcc<-((dn31.trinuc$gca*dn31.trinuc$tcc)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gca.tcg<-((dn31.trinuc$gca*dn31.trinuc$tcg)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gca.tct<-((dn31.trinuc$gca*dn31.trinuc$tct)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as

#Stop codon
#cps.gca.tga<-((dn31.trinuc$gca*dn31.trinuc$tga)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gca.tgc<-((dn31.trinuc$gca*dn31.trinuc$tgc)/(dn31.count.aa$a*dn31.count.aa$c))/dn31.count.diaa$ac
cps.gca.tgg<-((dn31.trinuc$gca*dn31.trinuc$tgg)/(dn31.count.aa$a*dn31.count.aa$w))/dn31.count.diaa$aw
cps.gca.tgt<-((dn31.trinuc$gca*dn31.trinuc$tgt)/(dn31.count.aa$a*dn31.count.aa$c))/dn31.count.diaa$ac

cps.gca.tta<-((dn31.trinuc$gca*dn31.trinuc$tta)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gca.ttc<-((dn31.trinuc$gca*dn31.trinuc$ttc)/(dn31.count.aa$a*dn31.count.aa$f))/dn31.count.diaa$af
cps.gca.ttg<-((dn31.trinuc$gca*dn31.trinuc$ttg)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gca.ttt<-((dn31.trinuc$gca*dn31.trinuc$ttt)/(dn31.count.aa$a*dn31.count.aa$f))/dn31.count.diaa$af








cps.gcc.aaa<-((dn31.trinuc$gcc*dn31.trinuc$aaa)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcc.aac<-((dn31.trinuc$gcc*dn31.trinuc$aac)/(dn31.count.aa$a*dn31.count.aa$n))/dn31.count.diaa$an
cps.gcc.aag<-((dn31.trinuc$gcc*dn31.trinuc$aag)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcc.aat<-((dn31.trinuc$gcc*dn31.trinuc$aat)/(dn31.count.aa$a*dn31.count.aa$n))/dn31.count.diaa$an

cps.gcc.aca<-((dn31.trinuc$gcc*dn31.trinuc$aca)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gcc.acc<-((dn31.trinuc$gcc*dn31.trinuc$acc)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gcc.acg<-((dn31.trinuc$gcc*dn31.trinuc$acg)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gcc.act<-((dn31.trinuc$gcc*dn31.trinuc$act)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at

cps.gcc.aga<-((dn31.trinuc$gcc*dn31.trinuc$aga)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcc.agc<-((dn31.trinuc$gcc*dn31.trinuc$agc)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gcc.agg<-((dn31.trinuc$gcc*dn31.trinuc$agg)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcc.agt<-((dn31.trinuc$gcc*dn31.trinuc$agt)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as

cps.gcc.ata<-((dn31.trinuc$gcc*dn31.trinuc$ata)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai
cps.gcc.atc<-((dn31.trinuc$gcc*dn31.trinuc$atc)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai
cps.gcc.atg<-((dn31.trinuc$gcc*dn31.trinuc$atg)/(dn31.count.aa$a*dn31.count.aa$m))/dn31.count.diaa$am
cps.gcc.att<-((dn31.trinuc$gcc*dn31.trinuc$att)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai

cps.gcc.caa<-((dn31.trinuc$gcc*dn31.trinuc$caa)/(dn31.count.aa$a*dn31.count.aa$q))/dn31.count.diaa$aq
cps.gcc.cac<-((dn31.trinuc$gcc*dn31.trinuc$cac)/(dn31.count.aa$a*dn31.count.aa$h))/dn31.count.diaa$ah
cps.gcc.cag<-((dn31.trinuc$gcc*dn31.trinuc$cag)/(dn31.count.aa$a*dn31.count.aa$q))/dn31.count.diaa$aq
cps.gcc.cat<-((dn31.trinuc$gcc*dn31.trinuc$cat)/(dn31.count.aa$a*dn31.count.aa$h))/dn31.count.diaa$ah

cps.gcc.cca<-((dn31.trinuc$gcc*dn31.trinuc$cca)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gcc.ccc<-((dn31.trinuc$gcc*dn31.trinuc$ccc)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gcc.ccg<-((dn31.trinuc$gcc*dn31.trinuc$ccg)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gcc.cct<-((dn31.trinuc$gcc*dn31.trinuc$cct)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap

cps.gcc.cga<-((dn31.trinuc$gcc*dn31.trinuc$cga)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcc.cgc<-((dn31.trinuc$gcc*dn31.trinuc$cgc)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcc.cgg<-((dn31.trinuc$gcc*dn31.trinuc$cgg)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcc.cgt<-((dn31.trinuc$gcc*dn31.trinuc$cgt)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar

cps.gcc.cta<-((dn31.trinuc$gcc*dn31.trinuc$cta)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcc.ctc<-((dn31.trinuc$gcc*dn31.trinuc$ctc)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcc.ctg<-((dn31.trinuc$gcc*dn31.trinuc$ctg)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcc.ctt<-((dn31.trinuc$gcc*dn31.trinuc$ctt)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al

cps.gcc.gaa<-((dn31.trinuc$gcc*dn31.trinuc$gaa)/(dn31.count.aa$a*dn31.count.aa$e))/dn31.count.diaa$ae
cps.gcc.gac<-((dn31.trinuc$gcc*dn31.trinuc$gac)/(dn31.count.aa$a*dn31.count.aa$d))/dn31.count.diaa$ad
cps.gcc.gag<-((dn31.trinuc$gcc*dn31.trinuc$gag)/(dn31.count.aa$a*dn31.count.aa$e))/dn31.count.diaa$ae
cps.gcc.gat<-((dn31.trinuc$gcc*dn31.trinuc$gat)/(dn31.count.aa$a*dn31.count.aa$d))/dn31.count.diaa$ad

cps.gcc.gca<-((dn31.trinuc$gcc*dn31.trinuc$gca)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gcc.gcc<-((dn31.trinuc$gcc*dn31.trinuc$gcc)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gcc.gcg<-((dn31.trinuc$gcc*dn31.trinuc$gcg)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gcc.gct<-((dn31.trinuc$gcc*dn31.trinuc$gct)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa

cps.gcc.gga<-((dn31.trinuc$gcc*dn31.trinuc$gga)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gcc.ggc<-((dn31.trinuc$gcc*dn31.trinuc$ggc)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gcc.ggg<-((dn31.trinuc$gcc*dn31.trinuc$ggg)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gcc.ggt<-((dn31.trinuc$gcc*dn31.trinuc$ggt)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag

cps.gcc.gta<-((dn31.trinuc$gcc*dn31.trinuc$gta)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gcc.gtc<-((dn31.trinuc$gcc*dn31.trinuc$gtc)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gcc.gtg<-((dn31.trinuc$gcc*dn31.trinuc$gtg)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gcc.gtt<-((dn31.trinuc$gcc*dn31.trinuc$gtt)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av

#Stop codon
#cps.gcc.taa<-((dn31.trinuc$gcc*dn31.trinuc$taa)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcc.tac<-((dn31.trinuc$gcc*dn31.trinuc$tac)/(dn31.count.aa$a*dn31.count.aa$y))/dn31.count.diaa$ay
#Stop codon
#cps.gcc.tag<-((dn31.trinuc$gcc*dn31.trinuc$tag)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcc.tat<-((dn31.trinuc$gcc*dn31.trinuc$tat)/(dn31.count.aa$a*dn31.count.aa$y))/dn31.count.diaa$ay

cps.gcc.tca<-((dn31.trinuc$gcc*dn31.trinuc$tca)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gcc.tcc<-((dn31.trinuc$gcc*dn31.trinuc$tcc)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gcc.tcg<-((dn31.trinuc$gcc*dn31.trinuc$tcg)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gcc.tct<-((dn31.trinuc$gcc*dn31.trinuc$tct)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as

#Stop codon
#cps.gcc.tga<-((dn31.trinuc$gcc*dn31.trinuc$tga)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcc.tgc<-((dn31.trinuc$gcc*dn31.trinuc$tgc)/(dn31.count.aa$a*dn31.count.aa$c))/dn31.count.diaa$ac
cps.gcc.tgg<-((dn31.trinuc$gcc*dn31.trinuc$tgg)/(dn31.count.aa$a*dn31.count.aa$w))/dn31.count.diaa$aw
cps.gcc.tgt<-((dn31.trinuc$gcc*dn31.trinuc$tgt)/(dn31.count.aa$a*dn31.count.aa$c))/dn31.count.diaa$ac

cps.gcc.tta<-((dn31.trinuc$gcc*dn31.trinuc$tta)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcc.ttc<-((dn31.trinuc$gcc*dn31.trinuc$ttc)/(dn31.count.aa$a*dn31.count.aa$f))/dn31.count.diaa$af
cps.gcc.ttg<-((dn31.trinuc$gcc*dn31.trinuc$ttg)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcc.ttt<-((dn31.trinuc$gcc*dn31.trinuc$ttt)/(dn31.count.aa$a*dn31.count.aa$f))/dn31.count.diaa$af








cps.gcg.aaa<-((dn31.trinuc$gcg*dn31.trinuc$aaa)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcg.aac<-((dn31.trinuc$gcg*dn31.trinuc$aac)/(dn31.count.aa$a*dn31.count.aa$n))/dn31.count.diaa$an
cps.gcg.aag<-((dn31.trinuc$gcg*dn31.trinuc$aag)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcg.aat<-((dn31.trinuc$gcg*dn31.trinuc$aat)/(dn31.count.aa$a*dn31.count.aa$n))/dn31.count.diaa$an

cps.gcg.aca<-((dn31.trinuc$gcg*dn31.trinuc$aca)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gcg.acc<-((dn31.trinuc$gcg*dn31.trinuc$acc)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gcg.acg<-((dn31.trinuc$gcg*dn31.trinuc$acg)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gcg.act<-((dn31.trinuc$gcg*dn31.trinuc$act)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at

cps.gcg.aga<-((dn31.trinuc$gcg*dn31.trinuc$aga)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcg.agc<-((dn31.trinuc$gcg*dn31.trinuc$agc)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gcg.agg<-((dn31.trinuc$gcg*dn31.trinuc$agg)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcg.agt<-((dn31.trinuc$gcg*dn31.trinuc$agt)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as

cps.gcg.ata<-((dn31.trinuc$gcg*dn31.trinuc$ata)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai
cps.gcg.atc<-((dn31.trinuc$gcg*dn31.trinuc$atc)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai
cps.gcg.atg<-((dn31.trinuc$gcg*dn31.trinuc$atg)/(dn31.count.aa$a*dn31.count.aa$m))/dn31.count.diaa$am
cps.gcg.att<-((dn31.trinuc$gcg*dn31.trinuc$att)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai

cps.gcg.caa<-((dn31.trinuc$gcg*dn31.trinuc$caa)/(dn31.count.aa$a*dn31.count.aa$q))/dn31.count.diaa$aq
cps.gcg.cac<-((dn31.trinuc$gcg*dn31.trinuc$cac)/(dn31.count.aa$a*dn31.count.aa$h))/dn31.count.diaa$ah
cps.gcg.cag<-((dn31.trinuc$gcg*dn31.trinuc$cag)/(dn31.count.aa$a*dn31.count.aa$q))/dn31.count.diaa$aq
cps.gcg.cat<-((dn31.trinuc$gcg*dn31.trinuc$cat)/(dn31.count.aa$a*dn31.count.aa$h))/dn31.count.diaa$ah

cps.gcg.cca<-((dn31.trinuc$gcg*dn31.trinuc$cca)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gcg.ccc<-((dn31.trinuc$gcg*dn31.trinuc$ccc)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gcg.ccg<-((dn31.trinuc$gcg*dn31.trinuc$ccg)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gcg.cct<-((dn31.trinuc$gcg*dn31.trinuc$cct)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap

cps.gcg.cga<-((dn31.trinuc$gcg*dn31.trinuc$cga)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcg.cgc<-((dn31.trinuc$gcg*dn31.trinuc$cgc)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcg.cgg<-((dn31.trinuc$gcg*dn31.trinuc$cgg)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gcg.cgt<-((dn31.trinuc$gcg*dn31.trinuc$cgt)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar

cps.gcg.cta<-((dn31.trinuc$gcg*dn31.trinuc$cta)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcg.ctc<-((dn31.trinuc$gcg*dn31.trinuc$ctc)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcg.ctg<-((dn31.trinuc$gcg*dn31.trinuc$ctg)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcg.ctt<-((dn31.trinuc$gcg*dn31.trinuc$ctt)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al

cps.gcg.gaa<-((dn31.trinuc$gcg*dn31.trinuc$gaa)/(dn31.count.aa$a*dn31.count.aa$e))/dn31.count.diaa$ae
cps.gcg.gac<-((dn31.trinuc$gcg*dn31.trinuc$gac)/(dn31.count.aa$a*dn31.count.aa$d))/dn31.count.diaa$ad
cps.gcg.gag<-((dn31.trinuc$gcg*dn31.trinuc$gag)/(dn31.count.aa$a*dn31.count.aa$e))/dn31.count.diaa$ae
cps.gcg.gat<-((dn31.trinuc$gcg*dn31.trinuc$gat)/(dn31.count.aa$a*dn31.count.aa$d))/dn31.count.diaa$ad

cps.gcg.gca<-((dn31.trinuc$gcg*dn31.trinuc$gca)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gcg.gcc<-((dn31.trinuc$gcg*dn31.trinuc$gcc)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gcg.gcg<-((dn31.trinuc$gcg*dn31.trinuc$gcg)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gcg.gct<-((dn31.trinuc$gcg*dn31.trinuc$gct)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa

cps.gcg.gga<-((dn31.trinuc$gcg*dn31.trinuc$gga)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gcg.ggc<-((dn31.trinuc$gcg*dn31.trinuc$ggc)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gcg.ggg<-((dn31.trinuc$gcg*dn31.trinuc$ggg)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gcg.ggt<-((dn31.trinuc$gcg*dn31.trinuc$ggt)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag

cps.gcg.gta<-((dn31.trinuc$gcg*dn31.trinuc$gta)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gcg.gtc<-((dn31.trinuc$gcg*dn31.trinuc$gtc)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gcg.gtg<-((dn31.trinuc$gcg*dn31.trinuc$gtg)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gcg.gtt<-((dn31.trinuc$gcg*dn31.trinuc$gtt)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av

#Stop codon
#cps.gcg.taa<-((dn31.trinuc$gcg*dn31.trinuc$taa)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcg.tac<-((dn31.trinuc$gcg*dn31.trinuc$tac)/(dn31.count.aa$a*dn31.count.aa$y))/dn31.count.diaa$ay
#Stop codon
#cps.gcg.tag<-((dn31.trinuc$gcg*dn31.trinuc$tag)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcg.tat<-((dn31.trinuc$gcg*dn31.trinuc$tat)/(dn31.count.aa$a*dn31.count.aa$y))/dn31.count.diaa$ay

cps.gcg.tca<-((dn31.trinuc$gcg*dn31.trinuc$tca)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gcg.tcc<-((dn31.trinuc$gcg*dn31.trinuc$tcc)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gcg.tcg<-((dn31.trinuc$gcg*dn31.trinuc$tcg)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gcg.tct<-((dn31.trinuc$gcg*dn31.trinuc$tct)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as

#Stop codon
#cps.gcg.tga<-((dn31.trinuc$gcg*dn31.trinuc$tga)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gcg.tgc<-((dn31.trinuc$gcg*dn31.trinuc$tgc)/(dn31.count.aa$a*dn31.count.aa$c))/dn31.count.diaa$ac
cps.gcg.tgg<-((dn31.trinuc$gcg*dn31.trinuc$tgg)/(dn31.count.aa$a*dn31.count.aa$w))/dn31.count.diaa$aw
cps.gcg.tgt<-((dn31.trinuc$gcg*dn31.trinuc$tgt)/(dn31.count.aa$a*dn31.count.aa$c))/dn31.count.diaa$ac

cps.gcg.tta<-((dn31.trinuc$gcg*dn31.trinuc$tta)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcg.ttc<-((dn31.trinuc$gcg*dn31.trinuc$ttc)/(dn31.count.aa$a*dn31.count.aa$f))/dn31.count.diaa$af
cps.gcg.ttg<-((dn31.trinuc$gcg*dn31.trinuc$ttg)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gcg.ttt<-((dn31.trinuc$gcg*dn31.trinuc$ttt)/(dn31.count.aa$a*dn31.count.aa$f))/dn31.count.diaa$af








cps.gct.aaa<-((dn31.trinuc$gct*dn31.trinuc$aaa)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gct.aac<-((dn31.trinuc$gct*dn31.trinuc$aac)/(dn31.count.aa$a*dn31.count.aa$n))/dn31.count.diaa$an
cps.gct.aag<-((dn31.trinuc$gct*dn31.trinuc$aag)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gct.aat<-((dn31.trinuc$gct*dn31.trinuc$aat)/(dn31.count.aa$a*dn31.count.aa$n))/dn31.count.diaa$an

cps.gct.aca<-((dn31.trinuc$gct*dn31.trinuc$aca)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gct.acc<-((dn31.trinuc$gct*dn31.trinuc$acc)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gct.acg<-((dn31.trinuc$gct*dn31.trinuc$acg)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at
cps.gct.act<-((dn31.trinuc$gct*dn31.trinuc$act)/(dn31.count.aa$a*dn31.count.aa$t))/dn31.count.diaa$at

cps.gct.aga<-((dn31.trinuc$gct*dn31.trinuc$aga)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gct.agc<-((dn31.trinuc$gct*dn31.trinuc$agc)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gct.agg<-((dn31.trinuc$gct*dn31.trinuc$agg)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gct.agt<-((dn31.trinuc$gct*dn31.trinuc$agt)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as

cps.gct.ata<-((dn31.trinuc$gct*dn31.trinuc$ata)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai
cps.gct.atc<-((dn31.trinuc$gct*dn31.trinuc$atc)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai
cps.gct.atg<-((dn31.trinuc$gct*dn31.trinuc$atg)/(dn31.count.aa$a*dn31.count.aa$m))/dn31.count.diaa$am
cps.gct.att<-((dn31.trinuc$gct*dn31.trinuc$att)/(dn31.count.aa$a*dn31.count.aa$i))/dn31.count.diaa$ai

cps.gct.caa<-((dn31.trinuc$gct*dn31.trinuc$caa)/(dn31.count.aa$a*dn31.count.aa$q))/dn31.count.diaa$aq
cps.gct.cac<-((dn31.trinuc$gct*dn31.trinuc$cac)/(dn31.count.aa$a*dn31.count.aa$h))/dn31.count.diaa$ah
cps.gct.cag<-((dn31.trinuc$gct*dn31.trinuc$cag)/(dn31.count.aa$a*dn31.count.aa$q))/dn31.count.diaa$aq
cps.gct.cat<-((dn31.trinuc$gct*dn31.trinuc$cat)/(dn31.count.aa$a*dn31.count.aa$h))/dn31.count.diaa$ah

cps.gct.cca<-((dn31.trinuc$gct*dn31.trinuc$cca)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gct.ccc<-((dn31.trinuc$gct*dn31.trinuc$ccc)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gct.ccg<-((dn31.trinuc$gct*dn31.trinuc$ccg)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap
cps.gct.cct<-((dn31.trinuc$gct*dn31.trinuc$cct)/(dn31.count.aa$a*dn31.count.aa$p))/dn31.count.diaa$ap

cps.gct.cga<-((dn31.trinuc$gct*dn31.trinuc$cga)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gct.cgc<-((dn31.trinuc$gct*dn31.trinuc$cgc)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gct.cgg<-((dn31.trinuc$gct*dn31.trinuc$cgg)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar
cps.gct.cgt<-((dn31.trinuc$gct*dn31.trinuc$cgt)/(dn31.count.aa$a*dn31.count.aa$r))/dn31.count.diaa$ar

cps.gct.cta<-((dn31.trinuc$gct*dn31.trinuc$cta)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gct.ctc<-((dn31.trinuc$gct*dn31.trinuc$ctc)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gct.ctg<-((dn31.trinuc$gct*dn31.trinuc$ctg)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gct.ctt<-((dn31.trinuc$gct*dn31.trinuc$ctt)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al

cps.gct.gaa<-((dn31.trinuc$gct*dn31.trinuc$gaa)/(dn31.count.aa$a*dn31.count.aa$e))/dn31.count.diaa$ae
cps.gct.gac<-((dn31.trinuc$gct*dn31.trinuc$gac)/(dn31.count.aa$a*dn31.count.aa$d))/dn31.count.diaa$ad
cps.gct.gag<-((dn31.trinuc$gct*dn31.trinuc$gag)/(dn31.count.aa$a*dn31.count.aa$e))/dn31.count.diaa$ae
cps.gct.gat<-((dn31.trinuc$gct*dn31.trinuc$gat)/(dn31.count.aa$a*dn31.count.aa$d))/dn31.count.diaa$ad

cps.gct.gca<-((dn31.trinuc$gct*dn31.trinuc$gca)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gct.gcc<-((dn31.trinuc$gct*dn31.trinuc$gcc)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gct.gcg<-((dn31.trinuc$gct*dn31.trinuc$gcg)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa
cps.gct.gct<-((dn31.trinuc$gct*dn31.trinuc$gct)/(dn31.count.aa$a*dn31.count.aa$a))/dn31.count.diaa$aa

cps.gct.gga<-((dn31.trinuc$gct*dn31.trinuc$gga)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gct.ggc<-((dn31.trinuc$gct*dn31.trinuc$ggc)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gct.ggg<-((dn31.trinuc$gct*dn31.trinuc$ggg)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag
cps.gct.ggt<-((dn31.trinuc$gct*dn31.trinuc$ggt)/(dn31.count.aa$a*dn31.count.aa$g))/dn31.count.diaa$ag

cps.gct.gta<-((dn31.trinuc$gct*dn31.trinuc$gta)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gct.gtc<-((dn31.trinuc$gct*dn31.trinuc$gtc)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gct.gtg<-((dn31.trinuc$gct*dn31.trinuc$gtg)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av
cps.gct.gtt<-((dn31.trinuc$gct*dn31.trinuc$gtt)/(dn31.count.aa$a*dn31.count.aa$v))/dn31.count.diaa$av

#Stop codon
#cps.gct.taa<-((dn31.trinuc$gct*dn31.trinuc$taa)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gct.tac<-((dn31.trinuc$gct*dn31.trinuc$tac)/(dn31.count.aa$a*dn31.count.aa$y))/dn31.count.diaa$ay
#Stop codon
#cps.gct.tag<-((dn31.trinuc$gct*dn31.trinuc$tag)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gct.tat<-((dn31.trinuc$gct*dn31.trinuc$tat)/(dn31.count.aa$a*dn31.count.aa$y))/dn31.count.diaa$ay

cps.gct.tca<-((dn31.trinuc$gct*dn31.trinuc$tca)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gct.tcc<-((dn31.trinuc$gct*dn31.trinuc$tcc)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gct.tcg<-((dn31.trinuc$gct*dn31.trinuc$tcg)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as
cps.gct.tct<-((dn31.trinuc$gct*dn31.trinuc$tct)/(dn31.count.aa$a*dn31.count.aa$s))/dn31.count.diaa$as

#Stop codon
#cps.gct.tga<-((dn31.trinuc$gct*dn31.trinuc$tga)/(dn31.count.aa$a*dn31.count.aa$k))/dn31.count.diaa$ak
cps.gct.tgc<-((dn31.trinuc$gct*dn31.trinuc$tgc)/(dn31.count.aa$a*dn31.count.aa$c))/dn31.count.diaa$ac
cps.gct.tgg<-((dn31.trinuc$gct*dn31.trinuc$tgg)/(dn31.count.aa$a*dn31.count.aa$w))/dn31.count.diaa$aw
cps.gct.tgt<-((dn31.trinuc$gct*dn31.trinuc$tgt)/(dn31.count.aa$a*dn31.count.aa$c))/dn31.count.diaa$ac

cps.gct.tta<-((dn31.trinuc$gct*dn31.trinuc$tta)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gct.ttc<-((dn31.trinuc$gct*dn31.trinuc$ttc)/(dn31.count.aa$a*dn31.count.aa$f))/dn31.count.diaa$af
cps.gct.ttg<-((dn31.trinuc$gct*dn31.trinuc$ttg)/(dn31.count.aa$a*dn31.count.aa$l))/dn31.count.diaa$al
cps.gct.ttt<-((dn31.trinuc$gct*dn31.trinuc$ttt)/(dn31.count.aa$a*dn31.count.aa$f))/dn31.count.diaa$af








cps.gga.aaa<-((dn31.trinuc$gga*dn31.trinuc$aaa)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.gga.aac<-((dn31.trinuc$gga*dn31.trinuc$aac)/(dn31.count.aa$g*dn31.count.aa$n))/dn31.count.diaa$gn
cps.gga.aag<-((dn31.trinuc$gga*dn31.trinuc$aag)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.gga.aat<-((dn31.trinuc$gga*dn31.trinuc$aat)/(dn31.count.aa$g*dn31.count.aa$n))/dn31.count.diaa$gn

cps.gga.aca<-((dn31.trinuc$gga*dn31.trinuc$aca)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.gga.acc<-((dn31.trinuc$gga*dn31.trinuc$acc)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.gga.acg<-((dn31.trinuc$gga*dn31.trinuc$acg)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.gga.act<-((dn31.trinuc$gga*dn31.trinuc$act)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt

cps.gga.aga<-((dn31.trinuc$gga*dn31.trinuc$aga)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.gga.agc<-((dn31.trinuc$gga*dn31.trinuc$agc)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.gga.agg<-((dn31.trinuc$gga*dn31.trinuc$agg)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.gga.agt<-((dn31.trinuc$gga*dn31.trinuc$agt)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs

cps.gga.ata<-((dn31.trinuc$gga*dn31.trinuc$ata)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi
cps.gga.atc<-((dn31.trinuc$gga*dn31.trinuc$atc)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi
cps.gga.atg<-((dn31.trinuc$gga*dn31.trinuc$atg)/(dn31.count.aa$g*dn31.count.aa$m))/dn31.count.diaa$gm
cps.gga.att<-((dn31.trinuc$gga*dn31.trinuc$att)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi

cps.gga.caa<-((dn31.trinuc$gga*dn31.trinuc$caa)/(dn31.count.aa$g*dn31.count.aa$q))/dn31.count.diaa$gq
cps.gga.cac<-((dn31.trinuc$gga*dn31.trinuc$cac)/(dn31.count.aa$g*dn31.count.aa$h))/dn31.count.diaa$gh
cps.gga.cag<-((dn31.trinuc$gga*dn31.trinuc$cag)/(dn31.count.aa$g*dn31.count.aa$q))/dn31.count.diaa$gq
cps.gga.cat<-((dn31.trinuc$gga*dn31.trinuc$cat)/(dn31.count.aa$g*dn31.count.aa$h))/dn31.count.diaa$gh

cps.gga.cca<-((dn31.trinuc$gga*dn31.trinuc$cca)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.gga.ccc<-((dn31.trinuc$gga*dn31.trinuc$ccc)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.gga.ccg<-((dn31.trinuc$gga*dn31.trinuc$ccg)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.gga.cct<-((dn31.trinuc$gga*dn31.trinuc$cct)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp

cps.gga.cga<-((dn31.trinuc$gga*dn31.trinuc$cga)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.gga.cgc<-((dn31.trinuc$gga*dn31.trinuc$cgc)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.gga.cgg<-((dn31.trinuc$gga*dn31.trinuc$cgg)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.gga.cgt<-((dn31.trinuc$gga*dn31.trinuc$cgt)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr

cps.gga.cta<-((dn31.trinuc$gga*dn31.trinuc$cta)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.gga.ctc<-((dn31.trinuc$gga*dn31.trinuc$ctc)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.gga.ctg<-((dn31.trinuc$gga*dn31.trinuc$ctg)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.gga.ctt<-((dn31.trinuc$gga*dn31.trinuc$ctt)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl

cps.gga.gaa<-((dn31.trinuc$gga*dn31.trinuc$gaa)/(dn31.count.aa$g*dn31.count.aa$e))/dn31.count.diaa$ge
cps.gga.gac<-((dn31.trinuc$gga*dn31.trinuc$gac)/(dn31.count.aa$g*dn31.count.aa$d))/dn31.count.diaa$gd
cps.gga.gag<-((dn31.trinuc$gga*dn31.trinuc$gag)/(dn31.count.aa$g*dn31.count.aa$e))/dn31.count.diaa$ge
cps.gga.gat<-((dn31.trinuc$gga*dn31.trinuc$gat)/(dn31.count.aa$g*dn31.count.aa$d))/dn31.count.diaa$gd

cps.gga.gca<-((dn31.trinuc$gga*dn31.trinuc$gca)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.gga.gcc<-((dn31.trinuc$gga*dn31.trinuc$gcc)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.gga.gcg<-((dn31.trinuc$gga*dn31.trinuc$gcg)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.gga.gct<-((dn31.trinuc$gga*dn31.trinuc$gct)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga

cps.gga.gga<-((dn31.trinuc$gga*dn31.trinuc$gga)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.gga.ggc<-((dn31.trinuc$gga*dn31.trinuc$ggc)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.gga.ggg<-((dn31.trinuc$gga*dn31.trinuc$ggg)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.gga.ggt<-((dn31.trinuc$gga*dn31.trinuc$ggt)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg

cps.gga.gta<-((dn31.trinuc$gga*dn31.trinuc$gta)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.gga.gtc<-((dn31.trinuc$gga*dn31.trinuc$gtc)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.gga.gtg<-((dn31.trinuc$gga*dn31.trinuc$gtg)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.gga.gtt<-((dn31.trinuc$gga*dn31.trinuc$gtt)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv

#Stop codon
#cps.gga.taa<-((dn31.trinuc$gga*dn31.trinuc$taa)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.gga.tac<-((dn31.trinuc$gga*dn31.trinuc$tac)/(dn31.count.aa$g*dn31.count.aa$y))/dn31.count.diaa$gy
#Stop codon
#cps.gga.tag<-((dn31.trinuc$gga*dn31.trinuc$tag)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.gga.tat<-((dn31.trinuc$gga*dn31.trinuc$tat)/(dn31.count.aa$g*dn31.count.aa$y))/dn31.count.diaa$gy

cps.gga.tca<-((dn31.trinuc$gga*dn31.trinuc$tca)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.gga.tcc<-((dn31.trinuc$gga*dn31.trinuc$tcc)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.gga.tcg<-((dn31.trinuc$gga*dn31.trinuc$tcg)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.gga.tct<-((dn31.trinuc$gga*dn31.trinuc$tct)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs

#Stop codon
#cps.gga.tga<-((dn31.trinuc$gga*dn31.trinuc$tga)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.gga.tgc<-((dn31.trinuc$gga*dn31.trinuc$tgc)/(dn31.count.aa$g*dn31.count.aa$c))/dn31.count.diaa$gc
cps.gga.tgg<-((dn31.trinuc$gga*dn31.trinuc$tgg)/(dn31.count.aa$g*dn31.count.aa$w))/dn31.count.diaa$gw
cps.gga.tgt<-((dn31.trinuc$gga*dn31.trinuc$tgt)/(dn31.count.aa$g*dn31.count.aa$c))/dn31.count.diaa$gc

cps.gga.tta<-((dn31.trinuc$gga*dn31.trinuc$tta)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.gga.ttc<-((dn31.trinuc$gga*dn31.trinuc$ttc)/(dn31.count.aa$g*dn31.count.aa$f))/dn31.count.diaa$gf
cps.gga.ttg<-((dn31.trinuc$gga*dn31.trinuc$ttg)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.gga.ttt<-((dn31.trinuc$gga*dn31.trinuc$ttt)/(dn31.count.aa$g*dn31.count.aa$f))/dn31.count.diaa$gf








cps.ggc.aaa<-((dn31.trinuc$ggc*dn31.trinuc$aaa)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggc.aac<-((dn31.trinuc$ggc*dn31.trinuc$aac)/(dn31.count.aa$g*dn31.count.aa$n))/dn31.count.diaa$gn
cps.ggc.aag<-((dn31.trinuc$ggc*dn31.trinuc$aag)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggc.aat<-((dn31.trinuc$ggc*dn31.trinuc$aat)/(dn31.count.aa$g*dn31.count.aa$n))/dn31.count.diaa$gn

cps.ggc.aca<-((dn31.trinuc$ggc*dn31.trinuc$aca)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggc.acc<-((dn31.trinuc$ggc*dn31.trinuc$acc)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggc.acg<-((dn31.trinuc$ggc*dn31.trinuc$acg)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggc.act<-((dn31.trinuc$ggc*dn31.trinuc$act)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt

cps.ggc.aga<-((dn31.trinuc$ggc*dn31.trinuc$aga)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggc.agc<-((dn31.trinuc$ggc*dn31.trinuc$agc)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggc.agg<-((dn31.trinuc$ggc*dn31.trinuc$agg)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggc.agt<-((dn31.trinuc$ggc*dn31.trinuc$agt)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs

cps.ggc.ata<-((dn31.trinuc$ggc*dn31.trinuc$ata)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi
cps.ggc.atc<-((dn31.trinuc$ggc*dn31.trinuc$atc)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi
cps.ggc.atg<-((dn31.trinuc$ggc*dn31.trinuc$atg)/(dn31.count.aa$g*dn31.count.aa$m))/dn31.count.diaa$gm
cps.ggc.att<-((dn31.trinuc$ggc*dn31.trinuc$att)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi

cps.ggc.caa<-((dn31.trinuc$ggc*dn31.trinuc$caa)/(dn31.count.aa$g*dn31.count.aa$q))/dn31.count.diaa$gq
cps.ggc.cac<-((dn31.trinuc$ggc*dn31.trinuc$cac)/(dn31.count.aa$g*dn31.count.aa$h))/dn31.count.diaa$gh
cps.ggc.cag<-((dn31.trinuc$ggc*dn31.trinuc$cag)/(dn31.count.aa$g*dn31.count.aa$q))/dn31.count.diaa$gq
cps.ggc.cat<-((dn31.trinuc$ggc*dn31.trinuc$cat)/(dn31.count.aa$g*dn31.count.aa$h))/dn31.count.diaa$gh

cps.ggc.cca<-((dn31.trinuc$ggc*dn31.trinuc$cca)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggc.ccc<-((dn31.trinuc$ggc*dn31.trinuc$ccc)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggc.ccg<-((dn31.trinuc$ggc*dn31.trinuc$ccg)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggc.cct<-((dn31.trinuc$ggc*dn31.trinuc$cct)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp

cps.ggc.cga<-((dn31.trinuc$ggc*dn31.trinuc$cga)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggc.cgc<-((dn31.trinuc$ggc*dn31.trinuc$cgc)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggc.cgg<-((dn31.trinuc$ggc*dn31.trinuc$cgg)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggc.cgt<-((dn31.trinuc$ggc*dn31.trinuc$cgt)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr

cps.ggc.cta<-((dn31.trinuc$ggc*dn31.trinuc$cta)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggc.ctc<-((dn31.trinuc$ggc*dn31.trinuc$ctc)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggc.ctg<-((dn31.trinuc$ggc*dn31.trinuc$ctg)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggc.ctt<-((dn31.trinuc$ggc*dn31.trinuc$ctt)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl

cps.ggc.gaa<-((dn31.trinuc$ggc*dn31.trinuc$gaa)/(dn31.count.aa$g*dn31.count.aa$e))/dn31.count.diaa$ge
cps.ggc.gac<-((dn31.trinuc$ggc*dn31.trinuc$gac)/(dn31.count.aa$g*dn31.count.aa$d))/dn31.count.diaa$gd
cps.ggc.gag<-((dn31.trinuc$ggc*dn31.trinuc$gag)/(dn31.count.aa$g*dn31.count.aa$e))/dn31.count.diaa$ge
cps.ggc.gat<-((dn31.trinuc$ggc*dn31.trinuc$gat)/(dn31.count.aa$g*dn31.count.aa$d))/dn31.count.diaa$gd

cps.ggc.gca<-((dn31.trinuc$ggc*dn31.trinuc$gca)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggc.gcc<-((dn31.trinuc$ggc*dn31.trinuc$gcc)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggc.gcg<-((dn31.trinuc$ggc*dn31.trinuc$gcg)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggc.gct<-((dn31.trinuc$ggc*dn31.trinuc$gct)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga

cps.ggc.gga<-((dn31.trinuc$ggc*dn31.trinuc$gga)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggc.ggc<-((dn31.trinuc$ggc*dn31.trinuc$ggc)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggc.ggg<-((dn31.trinuc$ggc*dn31.trinuc$ggg)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggc.ggt<-((dn31.trinuc$ggc*dn31.trinuc$ggt)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg

cps.ggc.gta<-((dn31.trinuc$ggc*dn31.trinuc$gta)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggc.gtc<-((dn31.trinuc$ggc*dn31.trinuc$gtc)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggc.gtg<-((dn31.trinuc$ggc*dn31.trinuc$gtg)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggc.gtt<-((dn31.trinuc$ggc*dn31.trinuc$gtt)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv

#Stop codon
#cps.ggc.taa<-((dn31.trinuc$ggc*dn31.trinuc$taa)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggc.tac<-((dn31.trinuc$ggc*dn31.trinuc$tac)/(dn31.count.aa$g*dn31.count.aa$y))/dn31.count.diaa$gy
#Stop codon
#cps.ggc.tag<-((dn31.trinuc$ggc*dn31.trinuc$tag)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggc.tat<-((dn31.trinuc$ggc*dn31.trinuc$tat)/(dn31.count.aa$g*dn31.count.aa$y))/dn31.count.diaa$gy

cps.ggc.tca<-((dn31.trinuc$ggc*dn31.trinuc$tca)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggc.tcc<-((dn31.trinuc$ggc*dn31.trinuc$tcc)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggc.tcg<-((dn31.trinuc$ggc*dn31.trinuc$tcg)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggc.tct<-((dn31.trinuc$ggc*dn31.trinuc$tct)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs

#Stop codon
#cps.ggc.tga<-((dn31.trinuc$ggc*dn31.trinuc$tga)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggc.tgc<-((dn31.trinuc$ggc*dn31.trinuc$tgc)/(dn31.count.aa$g*dn31.count.aa$c))/dn31.count.diaa$gc
cps.ggc.tgg<-((dn31.trinuc$ggc*dn31.trinuc$tgg)/(dn31.count.aa$g*dn31.count.aa$w))/dn31.count.diaa$gw
cps.ggc.tgt<-((dn31.trinuc$ggc*dn31.trinuc$tgt)/(dn31.count.aa$g*dn31.count.aa$c))/dn31.count.diaa$gc

cps.ggc.tta<-((dn31.trinuc$ggc*dn31.trinuc$tta)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggc.ttc<-((dn31.trinuc$ggc*dn31.trinuc$ttc)/(dn31.count.aa$g*dn31.count.aa$f))/dn31.count.diaa$gf
cps.ggc.ttg<-((dn31.trinuc$ggc*dn31.trinuc$ttg)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggc.ttt<-((dn31.trinuc$ggc*dn31.trinuc$ttt)/(dn31.count.aa$g*dn31.count.aa$f))/dn31.count.diaa$gf








cps.ggg.aaa<-((dn31.trinuc$ggg*dn31.trinuc$aaa)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggg.aac<-((dn31.trinuc$ggg*dn31.trinuc$aac)/(dn31.count.aa$g*dn31.count.aa$n))/dn31.count.diaa$gn
cps.ggg.aag<-((dn31.trinuc$ggg*dn31.trinuc$aag)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggg.aat<-((dn31.trinuc$ggg*dn31.trinuc$aat)/(dn31.count.aa$g*dn31.count.aa$n))/dn31.count.diaa$gn

cps.ggg.aca<-((dn31.trinuc$ggg*dn31.trinuc$aca)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggg.acc<-((dn31.trinuc$ggg*dn31.trinuc$acc)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggg.acg<-((dn31.trinuc$ggg*dn31.trinuc$acg)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggg.act<-((dn31.trinuc$ggg*dn31.trinuc$act)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt

cps.ggg.aga<-((dn31.trinuc$ggg*dn31.trinuc$aga)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggg.agc<-((dn31.trinuc$ggg*dn31.trinuc$agc)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggg.agg<-((dn31.trinuc$ggg*dn31.trinuc$agg)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggg.agt<-((dn31.trinuc$ggg*dn31.trinuc$agt)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs

cps.ggg.ata<-((dn31.trinuc$ggg*dn31.trinuc$ata)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi
cps.ggg.atc<-((dn31.trinuc$ggg*dn31.trinuc$atc)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi
cps.ggg.atg<-((dn31.trinuc$ggg*dn31.trinuc$atg)/(dn31.count.aa$g*dn31.count.aa$m))/dn31.count.diaa$gm
cps.ggg.att<-((dn31.trinuc$ggg*dn31.trinuc$att)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi

cps.ggg.caa<-((dn31.trinuc$ggg*dn31.trinuc$caa)/(dn31.count.aa$g*dn31.count.aa$q))/dn31.count.diaa$gq
cps.ggg.cac<-((dn31.trinuc$ggg*dn31.trinuc$cac)/(dn31.count.aa$g*dn31.count.aa$h))/dn31.count.diaa$gh
cps.ggg.cag<-((dn31.trinuc$ggg*dn31.trinuc$cag)/(dn31.count.aa$g*dn31.count.aa$q))/dn31.count.diaa$gq
cps.ggg.cat<-((dn31.trinuc$ggg*dn31.trinuc$cat)/(dn31.count.aa$g*dn31.count.aa$h))/dn31.count.diaa$gh

cps.ggg.cca<-((dn31.trinuc$ggg*dn31.trinuc$cca)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggg.ccc<-((dn31.trinuc$ggg*dn31.trinuc$ccc)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggg.ccg<-((dn31.trinuc$ggg*dn31.trinuc$ccg)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggg.cct<-((dn31.trinuc$ggg*dn31.trinuc$cct)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp

cps.ggg.cga<-((dn31.trinuc$ggg*dn31.trinuc$cga)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggg.cgc<-((dn31.trinuc$ggg*dn31.trinuc$cgc)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggg.cgg<-((dn31.trinuc$ggg*dn31.trinuc$cgg)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggg.cgt<-((dn31.trinuc$ggg*dn31.trinuc$cgt)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr

cps.ggg.cta<-((dn31.trinuc$ggg*dn31.trinuc$cta)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggg.ctc<-((dn31.trinuc$ggg*dn31.trinuc$ctc)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggg.ctg<-((dn31.trinuc$ggg*dn31.trinuc$ctg)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggg.ctt<-((dn31.trinuc$ggg*dn31.trinuc$ctt)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl

cps.ggg.gaa<-((dn31.trinuc$ggg*dn31.trinuc$gaa)/(dn31.count.aa$g*dn31.count.aa$e))/dn31.count.diaa$ge
cps.ggg.gac<-((dn31.trinuc$ggg*dn31.trinuc$gac)/(dn31.count.aa$g*dn31.count.aa$d))/dn31.count.diaa$gd
cps.ggg.gag<-((dn31.trinuc$ggg*dn31.trinuc$gag)/(dn31.count.aa$g*dn31.count.aa$e))/dn31.count.diaa$ge
cps.ggg.gat<-((dn31.trinuc$ggg*dn31.trinuc$gat)/(dn31.count.aa$g*dn31.count.aa$d))/dn31.count.diaa$gd

cps.ggg.gca<-((dn31.trinuc$ggg*dn31.trinuc$gca)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggg.gcc<-((dn31.trinuc$ggg*dn31.trinuc$gcc)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggg.gcg<-((dn31.trinuc$ggg*dn31.trinuc$gcg)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggg.gct<-((dn31.trinuc$ggg*dn31.trinuc$gct)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga

cps.ggg.gga<-((dn31.trinuc$ggg*dn31.trinuc$gga)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggg.ggc<-((dn31.trinuc$ggg*dn31.trinuc$ggc)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggg.ggg<-((dn31.trinuc$ggg*dn31.trinuc$ggg)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggg.ggt<-((dn31.trinuc$ggg*dn31.trinuc$ggt)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg

cps.ggg.gta<-((dn31.trinuc$ggg*dn31.trinuc$gta)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggg.gtc<-((dn31.trinuc$ggg*dn31.trinuc$gtc)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggg.gtg<-((dn31.trinuc$ggg*dn31.trinuc$gtg)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggg.gtt<-((dn31.trinuc$ggg*dn31.trinuc$gtt)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv

#Stop codon
#cps.ggg.taa<-((dn31.trinuc$ggg*dn31.trinuc$taa)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggg.tac<-((dn31.trinuc$ggg*dn31.trinuc$tac)/(dn31.count.aa$g*dn31.count.aa$y))/dn31.count.diaa$gy
#Stop codon
#cps.ggg.tag<-((dn31.trinuc$ggg*dn31.trinuc$tag)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggg.tat<-((dn31.trinuc$ggg*dn31.trinuc$tat)/(dn31.count.aa$g*dn31.count.aa$y))/dn31.count.diaa$gy

cps.ggg.tca<-((dn31.trinuc$ggg*dn31.trinuc$tca)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggg.tcc<-((dn31.trinuc$ggg*dn31.trinuc$tcc)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggg.tcg<-((dn31.trinuc$ggg*dn31.trinuc$tcg)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggg.tct<-((dn31.trinuc$ggg*dn31.trinuc$tct)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs

#Stop codon
#cps.ggg.tga<-((dn31.trinuc$ggg*dn31.trinuc$tga)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggg.tgc<-((dn31.trinuc$ggg*dn31.trinuc$tgc)/(dn31.count.aa$g*dn31.count.aa$c))/dn31.count.diaa$gc
cps.ggg.tgg<-((dn31.trinuc$ggg*dn31.trinuc$tgg)/(dn31.count.aa$g*dn31.count.aa$w))/dn31.count.diaa$gw
cps.ggg.tgt<-((dn31.trinuc$ggg*dn31.trinuc$tgt)/(dn31.count.aa$g*dn31.count.aa$c))/dn31.count.diaa$gc

cps.ggg.tta<-((dn31.trinuc$ggg*dn31.trinuc$tta)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggg.ttc<-((dn31.trinuc$ggg*dn31.trinuc$ttc)/(dn31.count.aa$g*dn31.count.aa$f))/dn31.count.diaa$gf
cps.ggg.ttg<-((dn31.trinuc$ggg*dn31.trinuc$ttg)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggg.ttt<-((dn31.trinuc$ggg*dn31.trinuc$ttt)/(dn31.count.aa$g*dn31.count.aa$f))/dn31.count.diaa$gf








cps.ggt.aaa<-((dn31.trinuc$ggt*dn31.trinuc$aaa)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggt.aac<-((dn31.trinuc$ggt*dn31.trinuc$aac)/(dn31.count.aa$g*dn31.count.aa$n))/dn31.count.diaa$gn
cps.ggt.aag<-((dn31.trinuc$ggt*dn31.trinuc$aag)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggt.aat<-((dn31.trinuc$ggt*dn31.trinuc$aat)/(dn31.count.aa$g*dn31.count.aa$n))/dn31.count.diaa$gn

cps.ggt.aca<-((dn31.trinuc$ggt*dn31.trinuc$aca)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggt.acc<-((dn31.trinuc$ggt*dn31.trinuc$acc)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggt.acg<-((dn31.trinuc$ggt*dn31.trinuc$acg)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt
cps.ggt.act<-((dn31.trinuc$ggt*dn31.trinuc$act)/(dn31.count.aa$g*dn31.count.aa$t))/dn31.count.diaa$gt

cps.ggt.aga<-((dn31.trinuc$ggt*dn31.trinuc$aga)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggt.agc<-((dn31.trinuc$ggt*dn31.trinuc$agc)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggt.agg<-((dn31.trinuc$ggt*dn31.trinuc$agg)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggt.agt<-((dn31.trinuc$ggt*dn31.trinuc$agt)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs

cps.ggt.ata<-((dn31.trinuc$ggt*dn31.trinuc$ata)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi
cps.ggt.atc<-((dn31.trinuc$ggt*dn31.trinuc$atc)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi
cps.ggt.atg<-((dn31.trinuc$ggt*dn31.trinuc$atg)/(dn31.count.aa$g*dn31.count.aa$m))/dn31.count.diaa$gm
cps.ggt.att<-((dn31.trinuc$ggt*dn31.trinuc$att)/(dn31.count.aa$g*dn31.count.aa$i))/dn31.count.diaa$gi

cps.ggt.caa<-((dn31.trinuc$ggt*dn31.trinuc$caa)/(dn31.count.aa$g*dn31.count.aa$q))/dn31.count.diaa$gq
cps.ggt.cac<-((dn31.trinuc$ggt*dn31.trinuc$cac)/(dn31.count.aa$g*dn31.count.aa$h))/dn31.count.diaa$gh
cps.ggt.cag<-((dn31.trinuc$ggt*dn31.trinuc$cag)/(dn31.count.aa$g*dn31.count.aa$q))/dn31.count.diaa$gq
cps.ggt.cat<-((dn31.trinuc$ggt*dn31.trinuc$cat)/(dn31.count.aa$g*dn31.count.aa$h))/dn31.count.diaa$gh

cps.ggt.cca<-((dn31.trinuc$ggt*dn31.trinuc$cca)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggt.ccc<-((dn31.trinuc$ggt*dn31.trinuc$ccc)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggt.ccg<-((dn31.trinuc$ggt*dn31.trinuc$ccg)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp
cps.ggt.cct<-((dn31.trinuc$ggt*dn31.trinuc$cct)/(dn31.count.aa$g*dn31.count.aa$p))/dn31.count.diaa$gp

cps.ggt.cga<-((dn31.trinuc$ggt*dn31.trinuc$cga)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggt.cgc<-((dn31.trinuc$ggt*dn31.trinuc$cgc)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggt.cgg<-((dn31.trinuc$ggt*dn31.trinuc$cgg)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr
cps.ggt.cgt<-((dn31.trinuc$ggt*dn31.trinuc$cgt)/(dn31.count.aa$g*dn31.count.aa$r))/dn31.count.diaa$gr

cps.ggt.cta<-((dn31.trinuc$ggt*dn31.trinuc$cta)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggt.ctc<-((dn31.trinuc$ggt*dn31.trinuc$ctc)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggt.ctg<-((dn31.trinuc$ggt*dn31.trinuc$ctg)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggt.ctt<-((dn31.trinuc$ggt*dn31.trinuc$ctt)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl

cps.ggt.gaa<-((dn31.trinuc$ggt*dn31.trinuc$gaa)/(dn31.count.aa$g*dn31.count.aa$e))/dn31.count.diaa$ge
cps.ggt.gac<-((dn31.trinuc$ggt*dn31.trinuc$gac)/(dn31.count.aa$g*dn31.count.aa$d))/dn31.count.diaa$gd
cps.ggt.gag<-((dn31.trinuc$ggt*dn31.trinuc$gag)/(dn31.count.aa$g*dn31.count.aa$e))/dn31.count.diaa$ge
cps.ggt.gat<-((dn31.trinuc$ggt*dn31.trinuc$gat)/(dn31.count.aa$g*dn31.count.aa$d))/dn31.count.diaa$gd

cps.ggt.gca<-((dn31.trinuc$ggt*dn31.trinuc$gca)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggt.gcc<-((dn31.trinuc$ggt*dn31.trinuc$gcc)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggt.gcg<-((dn31.trinuc$ggt*dn31.trinuc$gcg)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga
cps.ggt.gct<-((dn31.trinuc$ggt*dn31.trinuc$gct)/(dn31.count.aa$g*dn31.count.aa$a))/dn31.count.diaa$ga

cps.ggt.gga<-((dn31.trinuc$ggt*dn31.trinuc$gga)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggt.ggc<-((dn31.trinuc$ggt*dn31.trinuc$ggc)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggt.ggg<-((dn31.trinuc$ggt*dn31.trinuc$ggg)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg
cps.ggt.ggt<-((dn31.trinuc$ggt*dn31.trinuc$ggt)/(dn31.count.aa$g*dn31.count.aa$g))/dn31.count.diaa$gg

cps.ggt.gta<-((dn31.trinuc$ggt*dn31.trinuc$gta)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggt.gtc<-((dn31.trinuc$ggt*dn31.trinuc$gtc)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggt.gtg<-((dn31.trinuc$ggt*dn31.trinuc$gtg)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv
cps.ggt.gtt<-((dn31.trinuc$ggt*dn31.trinuc$gtt)/(dn31.count.aa$g*dn31.count.aa$v))/dn31.count.diaa$gv

#Stop codon
#cps.ggt.taa<-((dn31.trinuc$ggt*dn31.trinuc$taa)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggt.tac<-((dn31.trinuc$ggt*dn31.trinuc$tac)/(dn31.count.aa$g*dn31.count.aa$y))/dn31.count.diaa$gy
#Stop codon
#cps.ggt.tag<-((dn31.trinuc$ggt*dn31.trinuc$tag)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggt.tat<-((dn31.trinuc$ggt*dn31.trinuc$tat)/(dn31.count.aa$g*dn31.count.aa$y))/dn31.count.diaa$gy

cps.ggt.tca<-((dn31.trinuc$ggt*dn31.trinuc$tca)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggt.tcc<-((dn31.trinuc$ggt*dn31.trinuc$tcc)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggt.tcg<-((dn31.trinuc$ggt*dn31.trinuc$tcg)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs
cps.ggt.tct<-((dn31.trinuc$ggt*dn31.trinuc$tct)/(dn31.count.aa$g*dn31.count.aa$s))/dn31.count.diaa$gs

#Stop codon
#cps.ggt.tga<-((dn31.trinuc$ggt*dn31.trinuc$tga)/(dn31.count.aa$g*dn31.count.aa$k))/dn31.count.diaa$gk
cps.ggt.tgc<-((dn31.trinuc$ggt*dn31.trinuc$tgc)/(dn31.count.aa$g*dn31.count.aa$c))/dn31.count.diaa$gc
cps.ggt.tgg<-((dn31.trinuc$ggt*dn31.trinuc$tgg)/(dn31.count.aa$g*dn31.count.aa$w))/dn31.count.diaa$gw
cps.ggt.tgt<-((dn31.trinuc$ggt*dn31.trinuc$tgt)/(dn31.count.aa$g*dn31.count.aa$c))/dn31.count.diaa$gc

cps.ggt.tta<-((dn31.trinuc$ggt*dn31.trinuc$tta)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggt.ttc<-((dn31.trinuc$ggt*dn31.trinuc$ttc)/(dn31.count.aa$g*dn31.count.aa$f))/dn31.count.diaa$gf
cps.ggt.ttg<-((dn31.trinuc$ggt*dn31.trinuc$ttg)/(dn31.count.aa$g*dn31.count.aa$l))/dn31.count.diaa$gl
cps.ggt.ttt<-((dn31.trinuc$ggt*dn31.trinuc$ttt)/(dn31.count.aa$g*dn31.count.aa$f))/dn31.count.diaa$gf








cps.gta.aaa<-((dn31.trinuc$gta*dn31.trinuc$aaa)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gta.aac<-((dn31.trinuc$gta*dn31.trinuc$aac)/(dn31.count.aa$v*dn31.count.aa$n))/dn31.count.diaa$vn
cps.gta.aag<-((dn31.trinuc$gta*dn31.trinuc$aag)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gta.aat<-((dn31.trinuc$gta*dn31.trinuc$aat)/(dn31.count.aa$v*dn31.count.aa$n))/dn31.count.diaa$vn

cps.gta.aca<-((dn31.trinuc$gta*dn31.trinuc$aca)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gta.acc<-((dn31.trinuc$gta*dn31.trinuc$acc)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gta.acg<-((dn31.trinuc$gta*dn31.trinuc$acg)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gta.act<-((dn31.trinuc$gta*dn31.trinuc$act)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt

cps.gta.aga<-((dn31.trinuc$gta*dn31.trinuc$aga)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gta.agc<-((dn31.trinuc$gta*dn31.trinuc$agc)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gta.agg<-((dn31.trinuc$gta*dn31.trinuc$agg)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gta.agt<-((dn31.trinuc$gta*dn31.trinuc$agt)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs

cps.gta.ata<-((dn31.trinuc$gta*dn31.trinuc$ata)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi
cps.gta.atc<-((dn31.trinuc$gta*dn31.trinuc$atc)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi
cps.gta.atg<-((dn31.trinuc$gta*dn31.trinuc$atg)/(dn31.count.aa$v*dn31.count.aa$m))/dn31.count.diaa$vm
cps.gta.att<-((dn31.trinuc$gta*dn31.trinuc$att)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi

cps.gta.caa<-((dn31.trinuc$gta*dn31.trinuc$caa)/(dn31.count.aa$v*dn31.count.aa$q))/dn31.count.diaa$vq
cps.gta.cac<-((dn31.trinuc$gta*dn31.trinuc$cac)/(dn31.count.aa$v*dn31.count.aa$h))/dn31.count.diaa$vh
cps.gta.cag<-((dn31.trinuc$gta*dn31.trinuc$cag)/(dn31.count.aa$v*dn31.count.aa$q))/dn31.count.diaa$vq
cps.gta.cat<-((dn31.trinuc$gta*dn31.trinuc$cat)/(dn31.count.aa$v*dn31.count.aa$h))/dn31.count.diaa$vh

cps.gta.cca<-((dn31.trinuc$gta*dn31.trinuc$cca)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gta.ccc<-((dn31.trinuc$gta*dn31.trinuc$ccc)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gta.ccg<-((dn31.trinuc$gta*dn31.trinuc$ccg)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gta.cct<-((dn31.trinuc$gta*dn31.trinuc$cct)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp

cps.gta.cga<-((dn31.trinuc$gta*dn31.trinuc$cga)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gta.cgc<-((dn31.trinuc$gta*dn31.trinuc$cgc)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gta.cgg<-((dn31.trinuc$gta*dn31.trinuc$cgg)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gta.cgt<-((dn31.trinuc$gta*dn31.trinuc$cgt)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr

cps.gta.cta<-((dn31.trinuc$gta*dn31.trinuc$cta)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gta.ctc<-((dn31.trinuc$gta*dn31.trinuc$ctc)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gta.ctg<-((dn31.trinuc$gta*dn31.trinuc$ctg)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gta.ctt<-((dn31.trinuc$gta*dn31.trinuc$ctt)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl

cps.gta.gaa<-((dn31.trinuc$gta*dn31.trinuc$gaa)/(dn31.count.aa$v*dn31.count.aa$e))/dn31.count.diaa$ve
cps.gta.gac<-((dn31.trinuc$gta*dn31.trinuc$gac)/(dn31.count.aa$v*dn31.count.aa$d))/dn31.count.diaa$vd
cps.gta.gag<-((dn31.trinuc$gta*dn31.trinuc$gag)/(dn31.count.aa$v*dn31.count.aa$e))/dn31.count.diaa$ve
cps.gta.gat<-((dn31.trinuc$gta*dn31.trinuc$gat)/(dn31.count.aa$v*dn31.count.aa$d))/dn31.count.diaa$vd

cps.gta.gca<-((dn31.trinuc$gta*dn31.trinuc$gca)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gta.gcc<-((dn31.trinuc$gta*dn31.trinuc$gcc)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gta.gcg<-((dn31.trinuc$gta*dn31.trinuc$gcg)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gta.gct<-((dn31.trinuc$gta*dn31.trinuc$gct)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va

cps.gta.gga<-((dn31.trinuc$gta*dn31.trinuc$gga)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gta.ggc<-((dn31.trinuc$gta*dn31.trinuc$ggc)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gta.ggg<-((dn31.trinuc$gta*dn31.trinuc$ggg)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gta.ggt<-((dn31.trinuc$gta*dn31.trinuc$ggt)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg

cps.gta.gta<-((dn31.trinuc$gta*dn31.trinuc$gta)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gta.gtc<-((dn31.trinuc$gta*dn31.trinuc$gtc)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gta.gtg<-((dn31.trinuc$gta*dn31.trinuc$gtg)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gta.gtt<-((dn31.trinuc$gta*dn31.trinuc$gtt)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv

#Stop codon
#cps.gta.taa<-((dn31.trinuc$gta*dn31.trinuc$taa)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gta.tac<-((dn31.trinuc$gta*dn31.trinuc$tac)/(dn31.count.aa$v*dn31.count.aa$y))/dn31.count.diaa$vy
#Stop codon
#cps.gta.tag<-((dn31.trinuc$gta*dn31.trinuc$tag)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gta.tat<-((dn31.trinuc$gta*dn31.trinuc$tat)/(dn31.count.aa$v*dn31.count.aa$y))/dn31.count.diaa$vy

cps.gta.tca<-((dn31.trinuc$gta*dn31.trinuc$tca)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gta.tcc<-((dn31.trinuc$gta*dn31.trinuc$tcc)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gta.tcg<-((dn31.trinuc$gta*dn31.trinuc$tcg)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gta.tct<-((dn31.trinuc$gta*dn31.trinuc$tct)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs

#Stop codon
#cps.gta.tga<-((dn31.trinuc$gta*dn31.trinuc$tga)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gta.tgc<-((dn31.trinuc$gta*dn31.trinuc$tgc)/(dn31.count.aa$v*dn31.count.aa$c))/dn31.count.diaa$vc
cps.gta.tgg<-((dn31.trinuc$gta*dn31.trinuc$tgg)/(dn31.count.aa$v*dn31.count.aa$w))/dn31.count.diaa$vw
cps.gta.tgt<-((dn31.trinuc$gta*dn31.trinuc$tgt)/(dn31.count.aa$v*dn31.count.aa$c))/dn31.count.diaa$vc

cps.gta.tta<-((dn31.trinuc$gta*dn31.trinuc$tta)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gta.ttc<-((dn31.trinuc$gta*dn31.trinuc$ttc)/(dn31.count.aa$v*dn31.count.aa$f))/dn31.count.diaa$vf
cps.gta.ttg<-((dn31.trinuc$gta*dn31.trinuc$ttg)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gta.ttt<-((dn31.trinuc$gta*dn31.trinuc$ttt)/(dn31.count.aa$v*dn31.count.aa$f))/dn31.count.diaa$vf








cps.gtc.aaa<-((dn31.trinuc$gtc*dn31.trinuc$aaa)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtc.aac<-((dn31.trinuc$gtc*dn31.trinuc$aac)/(dn31.count.aa$v*dn31.count.aa$n))/dn31.count.diaa$vn
cps.gtc.aag<-((dn31.trinuc$gtc*dn31.trinuc$aag)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtc.aat<-((dn31.trinuc$gtc*dn31.trinuc$aat)/(dn31.count.aa$v*dn31.count.aa$n))/dn31.count.diaa$vn

cps.gtc.aca<-((dn31.trinuc$gtc*dn31.trinuc$aca)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtc.acc<-((dn31.trinuc$gtc*dn31.trinuc$acc)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtc.acg<-((dn31.trinuc$gtc*dn31.trinuc$acg)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtc.act<-((dn31.trinuc$gtc*dn31.trinuc$act)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt

cps.gtc.aga<-((dn31.trinuc$gtc*dn31.trinuc$aga)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtc.agc<-((dn31.trinuc$gtc*dn31.trinuc$agc)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtc.agg<-((dn31.trinuc$gtc*dn31.trinuc$agg)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtc.agt<-((dn31.trinuc$gtc*dn31.trinuc$agt)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs

cps.gtc.ata<-((dn31.trinuc$gtc*dn31.trinuc$ata)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi
cps.gtc.atc<-((dn31.trinuc$gtc*dn31.trinuc$atc)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi
cps.gtc.atg<-((dn31.trinuc$gtc*dn31.trinuc$atg)/(dn31.count.aa$v*dn31.count.aa$m))/dn31.count.diaa$vm
cps.gtc.att<-((dn31.trinuc$gtc*dn31.trinuc$att)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi

cps.gtc.caa<-((dn31.trinuc$gtc*dn31.trinuc$caa)/(dn31.count.aa$v*dn31.count.aa$q))/dn31.count.diaa$vq
cps.gtc.cac<-((dn31.trinuc$gtc*dn31.trinuc$cac)/(dn31.count.aa$v*dn31.count.aa$h))/dn31.count.diaa$vh
cps.gtc.cag<-((dn31.trinuc$gtc*dn31.trinuc$cag)/(dn31.count.aa$v*dn31.count.aa$q))/dn31.count.diaa$vq
cps.gtc.cat<-((dn31.trinuc$gtc*dn31.trinuc$cat)/(dn31.count.aa$v*dn31.count.aa$h))/dn31.count.diaa$vh

cps.gtc.cca<-((dn31.trinuc$gtc*dn31.trinuc$cca)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtc.ccc<-((dn31.trinuc$gtc*dn31.trinuc$ccc)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtc.ccg<-((dn31.trinuc$gtc*dn31.trinuc$ccg)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtc.cct<-((dn31.trinuc$gtc*dn31.trinuc$cct)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp

cps.gtc.cga<-((dn31.trinuc$gtc*dn31.trinuc$cga)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtc.cgc<-((dn31.trinuc$gtc*dn31.trinuc$cgc)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtc.cgg<-((dn31.trinuc$gtc*dn31.trinuc$cgg)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtc.cgt<-((dn31.trinuc$gtc*dn31.trinuc$cgt)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr

cps.gtc.cta<-((dn31.trinuc$gtc*dn31.trinuc$cta)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtc.ctc<-((dn31.trinuc$gtc*dn31.trinuc$ctc)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtc.ctg<-((dn31.trinuc$gtc*dn31.trinuc$ctg)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtc.ctt<-((dn31.trinuc$gtc*dn31.trinuc$ctt)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl

cps.gtc.gaa<-((dn31.trinuc$gtc*dn31.trinuc$gaa)/(dn31.count.aa$v*dn31.count.aa$e))/dn31.count.diaa$ve
cps.gtc.gac<-((dn31.trinuc$gtc*dn31.trinuc$gac)/(dn31.count.aa$v*dn31.count.aa$d))/dn31.count.diaa$vd
cps.gtc.gag<-((dn31.trinuc$gtc*dn31.trinuc$gag)/(dn31.count.aa$v*dn31.count.aa$e))/dn31.count.diaa$ve
cps.gtc.gat<-((dn31.trinuc$gtc*dn31.trinuc$gat)/(dn31.count.aa$v*dn31.count.aa$d))/dn31.count.diaa$vd

cps.gtc.gca<-((dn31.trinuc$gtc*dn31.trinuc$gca)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtc.gcc<-((dn31.trinuc$gtc*dn31.trinuc$gcc)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtc.gcg<-((dn31.trinuc$gtc*dn31.trinuc$gcg)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtc.gct<-((dn31.trinuc$gtc*dn31.trinuc$gct)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va

cps.gtc.gga<-((dn31.trinuc$gtc*dn31.trinuc$gga)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtc.ggc<-((dn31.trinuc$gtc*dn31.trinuc$ggc)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtc.ggg<-((dn31.trinuc$gtc*dn31.trinuc$ggg)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtc.ggt<-((dn31.trinuc$gtc*dn31.trinuc$ggt)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg

cps.gtc.gta<-((dn31.trinuc$gtc*dn31.trinuc$gta)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtc.gtc<-((dn31.trinuc$gtc*dn31.trinuc$gtc)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtc.gtg<-((dn31.trinuc$gtc*dn31.trinuc$gtg)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtc.gtt<-((dn31.trinuc$gtc*dn31.trinuc$gtt)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv

#Stop codon
#cps.gtc.taa<-((dn31.trinuc$gtc*dn31.trinuc$taa)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtc.tac<-((dn31.trinuc$gtc*dn31.trinuc$tac)/(dn31.count.aa$v*dn31.count.aa$y))/dn31.count.diaa$vy
#Stop codon
#cps.gtc.tag<-((dn31.trinuc$gtc*dn31.trinuc$tag)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtc.tat<-((dn31.trinuc$gtc*dn31.trinuc$tat)/(dn31.count.aa$v*dn31.count.aa$y))/dn31.count.diaa$vy

cps.gtc.tca<-((dn31.trinuc$gtc*dn31.trinuc$tca)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtc.tcc<-((dn31.trinuc$gtc*dn31.trinuc$tcc)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtc.tcg<-((dn31.trinuc$gtc*dn31.trinuc$tcg)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtc.tct<-((dn31.trinuc$gtc*dn31.trinuc$tct)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs

#Stop codon
#cps.gtc.tga<-((dn31.trinuc$gtc*dn31.trinuc$tga)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtc.tgc<-((dn31.trinuc$gtc*dn31.trinuc$tgc)/(dn31.count.aa$v*dn31.count.aa$c))/dn31.count.diaa$vc
cps.gtc.tgg<-((dn31.trinuc$gtc*dn31.trinuc$tgg)/(dn31.count.aa$v*dn31.count.aa$w))/dn31.count.diaa$vw
cps.gtc.tgt<-((dn31.trinuc$gtc*dn31.trinuc$tgt)/(dn31.count.aa$v*dn31.count.aa$c))/dn31.count.diaa$vc

cps.gtc.tta<-((dn31.trinuc$gtc*dn31.trinuc$tta)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtc.ttc<-((dn31.trinuc$gtc*dn31.trinuc$ttc)/(dn31.count.aa$v*dn31.count.aa$f))/dn31.count.diaa$vf
cps.gtc.ttg<-((dn31.trinuc$gtc*dn31.trinuc$ttg)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtc.ttt<-((dn31.trinuc$gtc*dn31.trinuc$ttt)/(dn31.count.aa$v*dn31.count.aa$f))/dn31.count.diaa$vf








cps.gtg.aaa<-((dn31.trinuc$gtg*dn31.trinuc$aaa)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtg.aac<-((dn31.trinuc$gtg*dn31.trinuc$aac)/(dn31.count.aa$v*dn31.count.aa$n))/dn31.count.diaa$vn
cps.gtg.aag<-((dn31.trinuc$gtg*dn31.trinuc$aag)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtg.aat<-((dn31.trinuc$gtg*dn31.trinuc$aat)/(dn31.count.aa$v*dn31.count.aa$n))/dn31.count.diaa$vn

cps.gtg.aca<-((dn31.trinuc$gtg*dn31.trinuc$aca)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtg.acc<-((dn31.trinuc$gtg*dn31.trinuc$acc)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtg.acg<-((dn31.trinuc$gtg*dn31.trinuc$acg)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtg.act<-((dn31.trinuc$gtg*dn31.trinuc$act)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt

cps.gtg.aga<-((dn31.trinuc$gtg*dn31.trinuc$aga)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtg.agc<-((dn31.trinuc$gtg*dn31.trinuc$agc)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtg.agg<-((dn31.trinuc$gtg*dn31.trinuc$agg)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtg.agt<-((dn31.trinuc$gtg*dn31.trinuc$agt)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs

cps.gtg.ata<-((dn31.trinuc$gtg*dn31.trinuc$ata)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi
cps.gtg.atc<-((dn31.trinuc$gtg*dn31.trinuc$atc)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi
cps.gtg.atg<-((dn31.trinuc$gtg*dn31.trinuc$atg)/(dn31.count.aa$v*dn31.count.aa$m))/dn31.count.diaa$vm
cps.gtg.att<-((dn31.trinuc$gtg*dn31.trinuc$att)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi

cps.gtg.caa<-((dn31.trinuc$gtg*dn31.trinuc$caa)/(dn31.count.aa$v*dn31.count.aa$q))/dn31.count.diaa$vq
cps.gtg.cac<-((dn31.trinuc$gtg*dn31.trinuc$cac)/(dn31.count.aa$v*dn31.count.aa$h))/dn31.count.diaa$vh
cps.gtg.cag<-((dn31.trinuc$gtg*dn31.trinuc$cag)/(dn31.count.aa$v*dn31.count.aa$q))/dn31.count.diaa$vq
cps.gtg.cat<-((dn31.trinuc$gtg*dn31.trinuc$cat)/(dn31.count.aa$v*dn31.count.aa$h))/dn31.count.diaa$vh

cps.gtg.cca<-((dn31.trinuc$gtg*dn31.trinuc$cca)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtg.ccc<-((dn31.trinuc$gtg*dn31.trinuc$ccc)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtg.ccg<-((dn31.trinuc$gtg*dn31.trinuc$ccg)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtg.cct<-((dn31.trinuc$gtg*dn31.trinuc$cct)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp

cps.gtg.cga<-((dn31.trinuc$gtg*dn31.trinuc$cga)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtg.cgc<-((dn31.trinuc$gtg*dn31.trinuc$cgc)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtg.cgg<-((dn31.trinuc$gtg*dn31.trinuc$cgg)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtg.cgt<-((dn31.trinuc$gtg*dn31.trinuc$cgt)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr

cps.gtg.cta<-((dn31.trinuc$gtg*dn31.trinuc$cta)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtg.ctc<-((dn31.trinuc$gtg*dn31.trinuc$ctc)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtg.ctg<-((dn31.trinuc$gtg*dn31.trinuc$ctg)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtg.ctt<-((dn31.trinuc$gtg*dn31.trinuc$ctt)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl

cps.gtg.gaa<-((dn31.trinuc$gtg*dn31.trinuc$gaa)/(dn31.count.aa$v*dn31.count.aa$e))/dn31.count.diaa$ve
cps.gtg.gac<-((dn31.trinuc$gtg*dn31.trinuc$gac)/(dn31.count.aa$v*dn31.count.aa$d))/dn31.count.diaa$vd
cps.gtg.gag<-((dn31.trinuc$gtg*dn31.trinuc$gag)/(dn31.count.aa$v*dn31.count.aa$e))/dn31.count.diaa$ve
cps.gtg.gat<-((dn31.trinuc$gtg*dn31.trinuc$gat)/(dn31.count.aa$v*dn31.count.aa$d))/dn31.count.diaa$vd

cps.gtg.gca<-((dn31.trinuc$gtg*dn31.trinuc$gca)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtg.gcc<-((dn31.trinuc$gtg*dn31.trinuc$gcc)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtg.gcg<-((dn31.trinuc$gtg*dn31.trinuc$gcg)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtg.gct<-((dn31.trinuc$gtg*dn31.trinuc$gct)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va

cps.gtg.gga<-((dn31.trinuc$gtg*dn31.trinuc$gga)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtg.ggc<-((dn31.trinuc$gtg*dn31.trinuc$ggc)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtg.ggg<-((dn31.trinuc$gtg*dn31.trinuc$ggg)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtg.ggt<-((dn31.trinuc$gtg*dn31.trinuc$ggt)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg

cps.gtg.gta<-((dn31.trinuc$gtg*dn31.trinuc$gta)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtg.gtc<-((dn31.trinuc$gtg*dn31.trinuc$gtc)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtg.gtg<-((dn31.trinuc$gtg*dn31.trinuc$gtg)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtg.gtt<-((dn31.trinuc$gtg*dn31.trinuc$gtt)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv

#Stop codon
#cps.gtg.taa<-((dn31.trinuc$gtg*dn31.trinuc$taa)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtg.tac<-((dn31.trinuc$gtg*dn31.trinuc$tac)/(dn31.count.aa$v*dn31.count.aa$y))/dn31.count.diaa$vy
#Stop codon
#cps.gtg.tag<-((dn31.trinuc$gtg*dn31.trinuc$tag)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtg.tat<-((dn31.trinuc$gtg*dn31.trinuc$tat)/(dn31.count.aa$v*dn31.count.aa$y))/dn31.count.diaa$vy

cps.gtg.tca<-((dn31.trinuc$gtg*dn31.trinuc$tca)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtg.tcc<-((dn31.trinuc$gtg*dn31.trinuc$tcc)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtg.tcg<-((dn31.trinuc$gtg*dn31.trinuc$tcg)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtg.tct<-((dn31.trinuc$gtg*dn31.trinuc$tct)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs

#Stop codon
#cps.gtg.tga<-((dn31.trinuc$gtg*dn31.trinuc$tga)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtg.tgc<-((dn31.trinuc$gtg*dn31.trinuc$tgc)/(dn31.count.aa$v*dn31.count.aa$c))/dn31.count.diaa$vc
cps.gtg.tgg<-((dn31.trinuc$gtg*dn31.trinuc$tgg)/(dn31.count.aa$v*dn31.count.aa$w))/dn31.count.diaa$vw
cps.gtg.tgt<-((dn31.trinuc$gtg*dn31.trinuc$tgt)/(dn31.count.aa$v*dn31.count.aa$c))/dn31.count.diaa$vc

cps.gtg.tta<-((dn31.trinuc$gtg*dn31.trinuc$tta)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtg.ttc<-((dn31.trinuc$gtg*dn31.trinuc$ttc)/(dn31.count.aa$v*dn31.count.aa$f))/dn31.count.diaa$vf
cps.gtg.ttg<-((dn31.trinuc$gtg*dn31.trinuc$ttg)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtg.ttt<-((dn31.trinuc$gtg*dn31.trinuc$ttt)/(dn31.count.aa$v*dn31.count.aa$f))/dn31.count.diaa$vf








cps.gtt.aaa<-((dn31.trinuc$gtt*dn31.trinuc$aaa)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtt.aac<-((dn31.trinuc$gtt*dn31.trinuc$aac)/(dn31.count.aa$v*dn31.count.aa$n))/dn31.count.diaa$vn
cps.gtt.aag<-((dn31.trinuc$gtt*dn31.trinuc$aag)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtt.aat<-((dn31.trinuc$gtt*dn31.trinuc$aat)/(dn31.count.aa$v*dn31.count.aa$n))/dn31.count.diaa$vn

cps.gtt.aca<-((dn31.trinuc$gtt*dn31.trinuc$aca)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtt.acc<-((dn31.trinuc$gtt*dn31.trinuc$acc)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtt.acg<-((dn31.trinuc$gtt*dn31.trinuc$acg)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt
cps.gtt.act<-((dn31.trinuc$gtt*dn31.trinuc$act)/(dn31.count.aa$v*dn31.count.aa$t))/dn31.count.diaa$vt

cps.gtt.aga<-((dn31.trinuc$gtt*dn31.trinuc$aga)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtt.agc<-((dn31.trinuc$gtt*dn31.trinuc$agc)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtt.agg<-((dn31.trinuc$gtt*dn31.trinuc$agg)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtt.agt<-((dn31.trinuc$gtt*dn31.trinuc$agt)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs

cps.gtt.ata<-((dn31.trinuc$gtt*dn31.trinuc$ata)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi
cps.gtt.atc<-((dn31.trinuc$gtt*dn31.trinuc$atc)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi
cps.gtt.atg<-((dn31.trinuc$gtt*dn31.trinuc$atg)/(dn31.count.aa$v*dn31.count.aa$m))/dn31.count.diaa$vm
cps.gtt.att<-((dn31.trinuc$gtt*dn31.trinuc$att)/(dn31.count.aa$v*dn31.count.aa$i))/dn31.count.diaa$vi

cps.gtt.caa<-((dn31.trinuc$gtt*dn31.trinuc$caa)/(dn31.count.aa$v*dn31.count.aa$q))/dn31.count.diaa$vq
cps.gtt.cac<-((dn31.trinuc$gtt*dn31.trinuc$cac)/(dn31.count.aa$v*dn31.count.aa$h))/dn31.count.diaa$vh
cps.gtt.cag<-((dn31.trinuc$gtt*dn31.trinuc$cag)/(dn31.count.aa$v*dn31.count.aa$q))/dn31.count.diaa$vq
cps.gtt.cat<-((dn31.trinuc$gtt*dn31.trinuc$cat)/(dn31.count.aa$v*dn31.count.aa$h))/dn31.count.diaa$vh

cps.gtt.cca<-((dn31.trinuc$gtt*dn31.trinuc$cca)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtt.ccc<-((dn31.trinuc$gtt*dn31.trinuc$ccc)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtt.ccg<-((dn31.trinuc$gtt*dn31.trinuc$ccg)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp
cps.gtt.cct<-((dn31.trinuc$gtt*dn31.trinuc$cct)/(dn31.count.aa$v*dn31.count.aa$p))/dn31.count.diaa$vp

cps.gtt.cga<-((dn31.trinuc$gtt*dn31.trinuc$cga)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtt.cgc<-((dn31.trinuc$gtt*dn31.trinuc$cgc)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtt.cgg<-((dn31.trinuc$gtt*dn31.trinuc$cgg)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr
cps.gtt.cgt<-((dn31.trinuc$gtt*dn31.trinuc$cgt)/(dn31.count.aa$v*dn31.count.aa$r))/dn31.count.diaa$vr

cps.gtt.cta<-((dn31.trinuc$gtt*dn31.trinuc$cta)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtt.ctc<-((dn31.trinuc$gtt*dn31.trinuc$ctc)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtt.ctg<-((dn31.trinuc$gtt*dn31.trinuc$ctg)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtt.ctt<-((dn31.trinuc$gtt*dn31.trinuc$ctt)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl

cps.gtt.gaa<-((dn31.trinuc$gtt*dn31.trinuc$gaa)/(dn31.count.aa$v*dn31.count.aa$e))/dn31.count.diaa$ve
cps.gtt.gac<-((dn31.trinuc$gtt*dn31.trinuc$gac)/(dn31.count.aa$v*dn31.count.aa$d))/dn31.count.diaa$vd
cps.gtt.gag<-((dn31.trinuc$gtt*dn31.trinuc$gag)/(dn31.count.aa$v*dn31.count.aa$e))/dn31.count.diaa$ve
cps.gtt.gat<-((dn31.trinuc$gtt*dn31.trinuc$gat)/(dn31.count.aa$v*dn31.count.aa$d))/dn31.count.diaa$vd

cps.gtt.gca<-((dn31.trinuc$gtt*dn31.trinuc$gca)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtt.gcc<-((dn31.trinuc$gtt*dn31.trinuc$gcc)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtt.gcg<-((dn31.trinuc$gtt*dn31.trinuc$gcg)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va
cps.gtt.gct<-((dn31.trinuc$gtt*dn31.trinuc$gct)/(dn31.count.aa$v*dn31.count.aa$a))/dn31.count.diaa$va

cps.gtt.gga<-((dn31.trinuc$gtt*dn31.trinuc$gga)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtt.ggc<-((dn31.trinuc$gtt*dn31.trinuc$ggc)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtt.ggg<-((dn31.trinuc$gtt*dn31.trinuc$ggg)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg
cps.gtt.ggt<-((dn31.trinuc$gtt*dn31.trinuc$ggt)/(dn31.count.aa$v*dn31.count.aa$g))/dn31.count.diaa$vg

cps.gtt.gta<-((dn31.trinuc$gtt*dn31.trinuc$gta)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtt.gtc<-((dn31.trinuc$gtt*dn31.trinuc$gtc)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtt.gtg<-((dn31.trinuc$gtt*dn31.trinuc$gtg)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv
cps.gtt.gtt<-((dn31.trinuc$gtt*dn31.trinuc$gtt)/(dn31.count.aa$v*dn31.count.aa$v))/dn31.count.diaa$vv

#Stop codon
#cps.gtt.taa<-((dn31.trinuc$gtt*dn31.trinuc$taa)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtt.tac<-((dn31.trinuc$gtt*dn31.trinuc$tac)/(dn31.count.aa$v*dn31.count.aa$y))/dn31.count.diaa$vy
#Stop codon
#cps.gtt.tag<-((dn31.trinuc$gtt*dn31.trinuc$tag)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtt.tat<-((dn31.trinuc$gtt*dn31.trinuc$tat)/(dn31.count.aa$v*dn31.count.aa$y))/dn31.count.diaa$vy

cps.gtt.tca<-((dn31.trinuc$gtt*dn31.trinuc$tca)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtt.tcc<-((dn31.trinuc$gtt*dn31.trinuc$tcc)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtt.tcg<-((dn31.trinuc$gtt*dn31.trinuc$tcg)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs
cps.gtt.tct<-((dn31.trinuc$gtt*dn31.trinuc$tct)/(dn31.count.aa$v*dn31.count.aa$s))/dn31.count.diaa$vs

#Stop codon
#cps.gtt.tga<-((dn31.trinuc$gtt*dn31.trinuc$tga)/(dn31.count.aa$v*dn31.count.aa$k))/dn31.count.diaa$vk
cps.gtt.tgc<-((dn31.trinuc$gtt*dn31.trinuc$tgc)/(dn31.count.aa$v*dn31.count.aa$c))/dn31.count.diaa$vc
cps.gtt.tgg<-((dn31.trinuc$gtt*dn31.trinuc$tgg)/(dn31.count.aa$v*dn31.count.aa$w))/dn31.count.diaa$vw
cps.gtt.tgt<-((dn31.trinuc$gtt*dn31.trinuc$tgt)/(dn31.count.aa$v*dn31.count.aa$c))/dn31.count.diaa$vc

cps.gtt.tta<-((dn31.trinuc$gtt*dn31.trinuc$tta)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtt.ttc<-((dn31.trinuc$gtt*dn31.trinuc$ttc)/(dn31.count.aa$v*dn31.count.aa$f))/dn31.count.diaa$vf
cps.gtt.ttg<-((dn31.trinuc$gtt*dn31.trinuc$ttg)/(dn31.count.aa$v*dn31.count.aa$l))/dn31.count.diaa$vl
cps.gtt.ttt<-((dn31.trinuc$gtt*dn31.trinuc$ttt)/(dn31.count.aa$v*dn31.count.aa$f))/dn31.count.diaa$vf










#Stop codon
#cps.taa.aaa<-((dn31.trinuc$taa*dn31.trinuc$aaa)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.taa.aac<-((dn31.trinuc$taa*dn31.trinuc$aac)/(dn31.count.aa$y*dn31.count.aa$n))/dn31.count.diaa$yn
#cps.taa.aag<-((dn31.trinuc$taa*dn31.trinuc$aag)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.taa.aat<-((dn31.trinuc$taa*dn31.trinuc$aat)/(dn31.count.aa$y*dn31.count.aa$n))/dn31.count.diaa$yn

#cps.taa.aca<-((dn31.trinuc$taa*dn31.trinuc$aca)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
#cps.taa.acc<-((dn31.trinuc$taa*dn31.trinuc$acc)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
#cps.taa.acg<-((dn31.trinuc$taa*dn31.trinuc$acg)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
#cps.taa.act<-((dn31.trinuc$taa*dn31.trinuc$act)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt

#cps.taa.aga<-((dn31.trinuc$taa*dn31.trinuc$aga)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.taa.agc<-((dn31.trinuc$taa*dn31.trinuc$agc)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
#cps.taa.agg<-((dn31.trinuc$taa*dn31.trinuc$agg)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.taa.agt<-((dn31.trinuc$taa*dn31.trinuc$agt)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys

#cps.taa.ata<-((dn31.trinuc$taa*dn31.trinuc$ata)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi
#cps.taa.atc<-((dn31.trinuc$taa*dn31.trinuc$atc)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi
#cps.taa.atg<-((dn31.trinuc$taa*dn31.trinuc$atg)/(dn31.count.aa$y*dn31.count.aa$m))/dn31.count.diaa$ym
#cps.taa.att<-((dn31.trinuc$taa*dn31.trinuc$att)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi

#cps.taa.caa<-((dn31.trinuc$taa*dn31.trinuc$caa)/(dn31.count.aa$y*dn31.count.aa$q))/dn31.count.diaa$yq
#cps.taa.cac<-((dn31.trinuc$taa*dn31.trinuc$cac)/(dn31.count.aa$y*dn31.count.aa$h))/dn31.count.diaa$yh
#cps.taa.cag<-((dn31.trinuc$taa*dn31.trinuc$cag)/(dn31.count.aa$y*dn31.count.aa$q))/dn31.count.diaa$yq
#cps.taa.cat<-((dn31.trinuc$taa*dn31.trinuc$cat)/(dn31.count.aa$y*dn31.count.aa$h))/dn31.count.diaa$yh

#cps.taa.cca<-((dn31.trinuc$taa*dn31.trinuc$cca)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
#cps.taa.ccc<-((dn31.trinuc$taa*dn31.trinuc$ccc)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
#cps.taa.ccg<-((dn31.trinuc$taa*dn31.trinuc$ccg)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
#cps.taa.cct<-((dn31.trinuc$taa*dn31.trinuc$cct)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp

#cps.taa.cga<-((dn31.trinuc$taa*dn31.trinuc$cga)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.taa.cgc<-((dn31.trinuc$taa*dn31.trinuc$cgc)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.taa.cgg<-((dn31.trinuc$taa*dn31.trinuc$cgg)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.taa.cgt<-((dn31.trinuc$taa*dn31.trinuc$cgt)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr

#cps.taa.cta<-((dn31.trinuc$taa*dn31.trinuc$cta)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.taa.ctc<-((dn31.trinuc$taa*dn31.trinuc$ctc)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.taa.ctg<-((dn31.trinuc$taa*dn31.trinuc$ctg)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.taa.ctt<-((dn31.trinuc$taa*dn31.trinuc$ctt)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl

#cps.taa.gaa<-((dn31.trinuc$taa*dn31.trinuc$gaa)/(dn31.count.aa$y*dn31.count.aa$e))/dn31.count.diaa$ye
#cps.taa.gac<-((dn31.trinuc$taa*dn31.trinuc$gac)/(dn31.count.aa$y*dn31.count.aa$d))/dn31.count.diaa$yd
#cps.taa.gag<-((dn31.trinuc$taa*dn31.trinuc$gag)/(dn31.count.aa$y*dn31.count.aa$e))/dn31.count.diaa$ye
#cps.taa.gat<-((dn31.trinuc$taa*dn31.trinuc$gat)/(dn31.count.aa$y*dn31.count.aa$d))/dn31.count.diaa$yd

#cps.taa.gca<-((dn31.trinuc$taa*dn31.trinuc$gca)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
#cps.taa.gcc<-((dn31.trinuc$taa*dn31.trinuc$gcc)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
#cps.taa.gcg<-((dn31.trinuc$taa*dn31.trinuc$gcg)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
#cps.taa.gct<-((dn31.trinuc$taa*dn31.trinuc$gct)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya

#cps.taa.gga<-((dn31.trinuc$taa*dn31.trinuc$gga)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
#cps.taa.ggc<-((dn31.trinuc$taa*dn31.trinuc$ggc)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
#cps.taa.ggg<-((dn31.trinuc$taa*dn31.trinuc$ggg)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
#cps.taa.ggt<-((dn31.trinuc$taa*dn31.trinuc$ggt)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg

#cps.taa.gta<-((dn31.trinuc$taa*dn31.trinuc$gta)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
#cps.taa.gtc<-((dn31.trinuc$taa*dn31.trinuc$gtc)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
#cps.taa.gtg<-((dn31.trinuc$taa*dn31.trinuc$gtg)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
#cps.taa.gtt<-((dn31.trinuc$taa*dn31.trinuc$gtt)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv

#Stop codon
#cps.taa.taa<-((dn31.trinuc$taa*dn31.trinuc$taa)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.taa.tac<-((dn31.trinuc$taa*dn31.trinuc$tac)/(dn31.count.aa$y*dn31.count.aa$y))/dn31.count.diaa$yy
#Stop codon
#cps.taa.tag<-((dn31.trinuc$taa*dn31.trinuc$tag)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.taa.tat<-((dn31.trinuc$taa*dn31.trinuc$tat)/(dn31.count.aa$y*dn31.count.aa$y))/dn31.count.diaa$yy

#cps.taa.tca<-((dn31.trinuc$taa*dn31.trinuc$tca)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
#cps.taa.tcc<-((dn31.trinuc$taa*dn31.trinuc$tcc)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
#cps.taa.tcg<-((dn31.trinuc$taa*dn31.trinuc$tcg)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
#cps.taa.tct<-((dn31.trinuc$taa*dn31.trinuc$tct)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys

#Stop codon
#cps.taa.tga<-((dn31.trinuc$taa*dn31.trinuc$tga)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.taa.tgc<-((dn31.trinuc$taa*dn31.trinuc$tgc)/(dn31.count.aa$y*dn31.count.aa$c))/dn31.count.diaa$yc
#cps.taa.tgg<-((dn31.trinuc$taa*dn31.trinuc$tgg)/(dn31.count.aa$y*dn31.count.aa$w))/dn31.count.diaa$yw
#cps.taa.tgt<-((dn31.trinuc$taa*dn31.trinuc$tgt)/(dn31.count.aa$y*dn31.count.aa$c))/dn31.count.diaa$yc

#cps.taa.tta<-((dn31.trinuc$taa*dn31.trinuc$tta)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.taa.ttc<-((dn31.trinuc$taa*dn31.trinuc$ttc)/(dn31.count.aa$y*dn31.count.aa$f))/dn31.count.diaa$yf
#cps.taa.ttg<-((dn31.trinuc$taa*dn31.trinuc$ttg)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.taa.ttt<-((dn31.trinuc$taa*dn31.trinuc$ttt)/(dn31.count.aa$y*dn31.count.aa$f))/dn31.count.diaa$yf



















cps.tac.aaa<-((dn31.trinuc$tac*dn31.trinuc$aaa)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tac.aac<-((dn31.trinuc$tac*dn31.trinuc$aac)/(dn31.count.aa$y*dn31.count.aa$n))/dn31.count.diaa$yn
cps.tac.aag<-((dn31.trinuc$tac*dn31.trinuc$aag)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tac.aat<-((dn31.trinuc$tac*dn31.trinuc$aat)/(dn31.count.aa$y*dn31.count.aa$n))/dn31.count.diaa$yn

cps.tac.aca<-((dn31.trinuc$tac*dn31.trinuc$aca)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
cps.tac.acc<-((dn31.trinuc$tac*dn31.trinuc$acc)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
cps.tac.acg<-((dn31.trinuc$tac*dn31.trinuc$acg)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
cps.tac.act<-((dn31.trinuc$tac*dn31.trinuc$act)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt

cps.tac.aga<-((dn31.trinuc$tac*dn31.trinuc$aga)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tac.agc<-((dn31.trinuc$tac*dn31.trinuc$agc)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
cps.tac.agg<-((dn31.trinuc$tac*dn31.trinuc$agg)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tac.agt<-((dn31.trinuc$tac*dn31.trinuc$agt)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys

cps.tac.ata<-((dn31.trinuc$tac*dn31.trinuc$ata)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi
cps.tac.atc<-((dn31.trinuc$tac*dn31.trinuc$atc)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi
cps.tac.atg<-((dn31.trinuc$tac*dn31.trinuc$atg)/(dn31.count.aa$y*dn31.count.aa$m))/dn31.count.diaa$ym
cps.tac.att<-((dn31.trinuc$tac*dn31.trinuc$att)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi

cps.tac.caa<-((dn31.trinuc$tac*dn31.trinuc$caa)/(dn31.count.aa$y*dn31.count.aa$q))/dn31.count.diaa$yq
cps.tac.cac<-((dn31.trinuc$tac*dn31.trinuc$cac)/(dn31.count.aa$y*dn31.count.aa$h))/dn31.count.diaa$yh
cps.tac.cag<-((dn31.trinuc$tac*dn31.trinuc$cag)/(dn31.count.aa$y*dn31.count.aa$q))/dn31.count.diaa$yq
cps.tac.cat<-((dn31.trinuc$tac*dn31.trinuc$cat)/(dn31.count.aa$y*dn31.count.aa$h))/dn31.count.diaa$yh

cps.tac.cca<-((dn31.trinuc$tac*dn31.trinuc$cca)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
cps.tac.ccc<-((dn31.trinuc$tac*dn31.trinuc$ccc)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
cps.tac.ccg<-((dn31.trinuc$tac*dn31.trinuc$ccg)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
cps.tac.cct<-((dn31.trinuc$tac*dn31.trinuc$cct)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp

cps.tac.cga<-((dn31.trinuc$tac*dn31.trinuc$cga)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tac.cgc<-((dn31.trinuc$tac*dn31.trinuc$cgc)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tac.cgg<-((dn31.trinuc$tac*dn31.trinuc$cgg)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tac.cgt<-((dn31.trinuc$tac*dn31.trinuc$cgt)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr

cps.tac.cta<-((dn31.trinuc$tac*dn31.trinuc$cta)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tac.ctc<-((dn31.trinuc$tac*dn31.trinuc$ctc)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tac.ctg<-((dn31.trinuc$tac*dn31.trinuc$ctg)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tac.ctt<-((dn31.trinuc$tac*dn31.trinuc$ctt)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl

cps.tac.gaa<-((dn31.trinuc$tac*dn31.trinuc$gaa)/(dn31.count.aa$y*dn31.count.aa$e))/dn31.count.diaa$ye
cps.tac.gac<-((dn31.trinuc$tac*dn31.trinuc$gac)/(dn31.count.aa$y*dn31.count.aa$d))/dn31.count.diaa$yd
cps.tac.gag<-((dn31.trinuc$tac*dn31.trinuc$gag)/(dn31.count.aa$y*dn31.count.aa$e))/dn31.count.diaa$ye
cps.tac.gat<-((dn31.trinuc$tac*dn31.trinuc$gat)/(dn31.count.aa$y*dn31.count.aa$d))/dn31.count.diaa$yd

cps.tac.gca<-((dn31.trinuc$tac*dn31.trinuc$gca)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
cps.tac.gcc<-((dn31.trinuc$tac*dn31.trinuc$gcc)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
cps.tac.gcg<-((dn31.trinuc$tac*dn31.trinuc$gcg)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
cps.tac.gct<-((dn31.trinuc$tac*dn31.trinuc$gct)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya

cps.tac.gga<-((dn31.trinuc$tac*dn31.trinuc$gga)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
cps.tac.ggc<-((dn31.trinuc$tac*dn31.trinuc$ggc)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
cps.tac.ggg<-((dn31.trinuc$tac*dn31.trinuc$ggg)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
cps.tac.ggt<-((dn31.trinuc$tac*dn31.trinuc$ggt)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg

cps.tac.gta<-((dn31.trinuc$tac*dn31.trinuc$gta)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
cps.tac.gtc<-((dn31.trinuc$tac*dn31.trinuc$gtc)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
cps.tac.gtg<-((dn31.trinuc$tac*dn31.trinuc$gtg)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
cps.tac.gtt<-((dn31.trinuc$tac*dn31.trinuc$gtt)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv

#Stop codon
#cps.tac.taa<-((dn31.trinuc$tac*dn31.trinuc$taa)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tac.tac<-((dn31.trinuc$tac*dn31.trinuc$tac)/(dn31.count.aa$y*dn31.count.aa$y))/dn31.count.diaa$yy
#Stop codon
#cps.tac.tag<-((dn31.trinuc$tac*dn31.trinuc$tag)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tac.tat<-((dn31.trinuc$tac*dn31.trinuc$tat)/(dn31.count.aa$y*dn31.count.aa$y))/dn31.count.diaa$yy

cps.tac.tca<-((dn31.trinuc$tac*dn31.trinuc$tca)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
cps.tac.tcc<-((dn31.trinuc$tac*dn31.trinuc$tcc)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
cps.tac.tcg<-((dn31.trinuc$tac*dn31.trinuc$tcg)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
cps.tac.tct<-((dn31.trinuc$tac*dn31.trinuc$tct)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys

#Stop codon
#cps.tac.tga<-((dn31.trinuc$tac*dn31.trinuc$tga)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tac.tgc<-((dn31.trinuc$tac*dn31.trinuc$tgc)/(dn31.count.aa$y*dn31.count.aa$c))/dn31.count.diaa$yc
cps.tac.tgg<-((dn31.trinuc$tac*dn31.trinuc$tgg)/(dn31.count.aa$y*dn31.count.aa$w))/dn31.count.diaa$yw
cps.tac.tgt<-((dn31.trinuc$tac*dn31.trinuc$tgt)/(dn31.count.aa$y*dn31.count.aa$c))/dn31.count.diaa$yc

cps.tac.tta<-((dn31.trinuc$tac*dn31.trinuc$tta)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tac.ttc<-((dn31.trinuc$tac*dn31.trinuc$ttc)/(dn31.count.aa$y*dn31.count.aa$f))/dn31.count.diaa$yf
cps.tac.ttg<-((dn31.trinuc$tac*dn31.trinuc$ttg)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tac.ttt<-((dn31.trinuc$tac*dn31.trinuc$ttt)/(dn31.count.aa$y*dn31.count.aa$f))/dn31.count.diaa$yf








#Stop codon

#cps.tag.aaa<-((dn31.trinuc$tag*dn31.trinuc$aaa)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.tag.aac<-((dn31.trinuc$tag*dn31.trinuc$aac)/(dn31.count.aa$y*dn31.count.aa$n))/dn31.count.diaa$yn
#cps.tag.aag<-((dn31.trinuc$tag*dn31.trinuc$aag)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.tag.aat<-((dn31.trinuc$tag*dn31.trinuc$aat)/(dn31.count.aa$y*dn31.count.aa$n))/dn31.count.diaa$yn

#cps.tag.aca<-((dn31.trinuc$tag*dn31.trinuc$aca)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
#cps.tag.acc<-((dn31.trinuc$tag*dn31.trinuc$acc)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
#cps.tag.acg<-((dn31.trinuc$tag*dn31.trinuc$acg)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
#cps.tag.act<-((dn31.trinuc$tag*dn31.trinuc$act)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt

#cps.tag.aga<-((dn31.trinuc$tag*dn31.trinuc$aga)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.tag.agc<-((dn31.trinuc$tag*dn31.trinuc$agc)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
#cps.tag.agg<-((dn31.trinuc$tag*dn31.trinuc$agg)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.tag.agt<-((dn31.trinuc$tag*dn31.trinuc$agt)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys

#cps.tag.ata<-((dn31.trinuc$tag*dn31.trinuc$ata)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi
#cps.tag.atc<-((dn31.trinuc$tag*dn31.trinuc$atc)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi
#cps.tag.atg<-((dn31.trinuc$tag*dn31.trinuc$atg)/(dn31.count.aa$y*dn31.count.aa$m))/dn31.count.diaa$ym
#cps.tag.att<-((dn31.trinuc$tag*dn31.trinuc$att)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi

#cps.tag.caa<-((dn31.trinuc$tag*dn31.trinuc$caa)/(dn31.count.aa$y*dn31.count.aa$q))/dn31.count.diaa$yq
#cps.tag.cac<-((dn31.trinuc$tag*dn31.trinuc$cac)/(dn31.count.aa$y*dn31.count.aa$h))/dn31.count.diaa$yh
#cps.tag.cag<-((dn31.trinuc$tag*dn31.trinuc$cag)/(dn31.count.aa$y*dn31.count.aa$q))/dn31.count.diaa$yq
#cps.tag.cat<-((dn31.trinuc$tag*dn31.trinuc$cat)/(dn31.count.aa$y*dn31.count.aa$h))/dn31.count.diaa$yh

#cps.tag.cca<-((dn31.trinuc$tag*dn31.trinuc$cca)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
#cps.tag.ccc<-((dn31.trinuc$tag*dn31.trinuc$ccc)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
#cps.tag.ccg<-((dn31.trinuc$tag*dn31.trinuc$ccg)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
#cps.tag.cct<-((dn31.trinuc$tag*dn31.trinuc$cct)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp

#cps.tag.cga<-((dn31.trinuc$tag*dn31.trinuc$cga)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.tag.cgc<-((dn31.trinuc$tag*dn31.trinuc$cgc)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.tag.cgg<-((dn31.trinuc$tag*dn31.trinuc$cgg)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
#cps.tag.cgt<-((dn31.trinuc$tag*dn31.trinuc$cgt)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr

#cps.tag.cta<-((dn31.trinuc$tag*dn31.trinuc$cta)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.tag.ctc<-((dn31.trinuc$tag*dn31.trinuc$ctc)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.tag.ctg<-((dn31.trinuc$tag*dn31.trinuc$ctg)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.tag.ctt<-((dn31.trinuc$tag*dn31.trinuc$ctt)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl

#cps.tag.gaa<-((dn31.trinuc$tag*dn31.trinuc$gaa)/(dn31.count.aa$y*dn31.count.aa$e))/dn31.count.diaa$ye
#cps.tag.gac<-((dn31.trinuc$tag*dn31.trinuc$gac)/(dn31.count.aa$y*dn31.count.aa$d))/dn31.count.diaa$yd
#cps.tag.gag<-((dn31.trinuc$tag*dn31.trinuc$gag)/(dn31.count.aa$y*dn31.count.aa$e))/dn31.count.diaa$ye
#cps.tag.gat<-((dn31.trinuc$tag*dn31.trinuc$gat)/(dn31.count.aa$y*dn31.count.aa$d))/dn31.count.diaa$yd

#cps.tag.gca<-((dn31.trinuc$tag*dn31.trinuc$gca)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
#cps.tag.gcc<-((dn31.trinuc$tag*dn31.trinuc$gcc)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
#cps.tag.gcg<-((dn31.trinuc$tag*dn31.trinuc$gcg)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
#cps.tag.gct<-((dn31.trinuc$tag*dn31.trinuc$gct)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya

#cps.tag.gga<-((dn31.trinuc$tag*dn31.trinuc$gga)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
#cps.tag.ggc<-((dn31.trinuc$tag*dn31.trinuc$ggc)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
#cps.tag.ggg<-((dn31.trinuc$tag*dn31.trinuc$ggg)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
#cps.tag.ggt<-((dn31.trinuc$tag*dn31.trinuc$ggt)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg

#cps.tag.gta<-((dn31.trinuc$tag*dn31.trinuc$gta)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
#cps.tag.gtc<-((dn31.trinuc$tag*dn31.trinuc$gtc)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
#cps.tag.gtg<-((dn31.trinuc$tag*dn31.trinuc$gtg)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
#cps.tag.gtt<-((dn31.trinuc$tag*dn31.trinuc$gtt)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv

#Stop codon
#cps.tag.taa<-((dn31.trinuc$tag*dn31.trinuc$taa)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.tag.tac<-((dn31.trinuc$tag*dn31.trinuc$tac)/(dn31.count.aa$y*dn31.count.aa$y))/dn31.count.diaa$yy
#Stop codon
#cps.tag.tag<-((dn31.trinuc$tag*dn31.trinuc$tag)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.tag.tat<-((dn31.trinuc$tag*dn31.trinuc$tat)/(dn31.count.aa$y*dn31.count.aa$y))/dn31.count.diaa$yy

#cps.tag.tca<-((dn31.trinuc$tag*dn31.trinuc$tca)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
#cps.tag.tcc<-((dn31.trinuc$tag*dn31.trinuc$tcc)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
#cps.tag.tcg<-((dn31.trinuc$tag*dn31.trinuc$tcg)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
#cps.tag.tct<-((dn31.trinuc$tag*dn31.trinuc$tct)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys

#Stop codon
#cps.tag.tga<-((dn31.trinuc$tag*dn31.trinuc$tga)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
#cps.tag.tgc<-((dn31.trinuc$tag*dn31.trinuc$tgc)/(dn31.count.aa$y*dn31.count.aa$c))/dn31.count.diaa$yc
#cps.tag.tgg<-((dn31.trinuc$tag*dn31.trinuc$tgg)/(dn31.count.aa$y*dn31.count.aa$w))/dn31.count.diaa$yw
#cps.tag.tgt<-((dn31.trinuc$tag*dn31.trinuc$tgt)/(dn31.count.aa$y*dn31.count.aa$c))/dn31.count.diaa$yc

#cps.tag.tta<-((dn31.trinuc$tag*dn31.trinuc$tta)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.tag.ttc<-((dn31.trinuc$tag*dn31.trinuc$ttc)/(dn31.count.aa$y*dn31.count.aa$f))/dn31.count.diaa$yf
#cps.tag.ttg<-((dn31.trinuc$tag*dn31.trinuc$ttg)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
#cps.tag.ttt<-((dn31.trinuc$tag*dn31.trinuc$ttt)/(dn31.count.aa$y*dn31.count.aa$f))/dn31.count.diaa$yf

















cps.tat.aaa<-((dn31.trinuc$tat*dn31.trinuc$aaa)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tat.aac<-((dn31.trinuc$tat*dn31.trinuc$aac)/(dn31.count.aa$y*dn31.count.aa$n))/dn31.count.diaa$yn
cps.tat.aag<-((dn31.trinuc$tat*dn31.trinuc$aag)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tat.aat<-((dn31.trinuc$tat*dn31.trinuc$aat)/(dn31.count.aa$y*dn31.count.aa$n))/dn31.count.diaa$yn

cps.tat.aca<-((dn31.trinuc$tat*dn31.trinuc$aca)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
cps.tat.acc<-((dn31.trinuc$tat*dn31.trinuc$acc)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
cps.tat.acg<-((dn31.trinuc$tat*dn31.trinuc$acg)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt
cps.tat.act<-((dn31.trinuc$tat*dn31.trinuc$act)/(dn31.count.aa$y*dn31.count.aa$t))/dn31.count.diaa$yt

cps.tat.aga<-((dn31.trinuc$tat*dn31.trinuc$aga)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tat.agc<-((dn31.trinuc$tat*dn31.trinuc$agc)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
cps.tat.agg<-((dn31.trinuc$tat*dn31.trinuc$agg)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tat.agt<-((dn31.trinuc$tat*dn31.trinuc$agt)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys

cps.tat.ata<-((dn31.trinuc$tat*dn31.trinuc$ata)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi
cps.tat.atc<-((dn31.trinuc$tat*dn31.trinuc$atc)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi
cps.tat.atg<-((dn31.trinuc$tat*dn31.trinuc$atg)/(dn31.count.aa$y*dn31.count.aa$m))/dn31.count.diaa$ym
cps.tat.att<-((dn31.trinuc$tat*dn31.trinuc$att)/(dn31.count.aa$y*dn31.count.aa$i))/dn31.count.diaa$yi

cps.tat.caa<-((dn31.trinuc$tat*dn31.trinuc$caa)/(dn31.count.aa$y*dn31.count.aa$q))/dn31.count.diaa$yq
cps.tat.cac<-((dn31.trinuc$tat*dn31.trinuc$cac)/(dn31.count.aa$y*dn31.count.aa$h))/dn31.count.diaa$yh
cps.tat.cag<-((dn31.trinuc$tat*dn31.trinuc$cag)/(dn31.count.aa$y*dn31.count.aa$q))/dn31.count.diaa$yq
cps.tat.cat<-((dn31.trinuc$tat*dn31.trinuc$cat)/(dn31.count.aa$y*dn31.count.aa$h))/dn31.count.diaa$yh

cps.tat.cca<-((dn31.trinuc$tat*dn31.trinuc$cca)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
cps.tat.ccc<-((dn31.trinuc$tat*dn31.trinuc$ccc)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
cps.tat.ccg<-((dn31.trinuc$tat*dn31.trinuc$ccg)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp
cps.tat.cct<-((dn31.trinuc$tat*dn31.trinuc$cct)/(dn31.count.aa$y*dn31.count.aa$p))/dn31.count.diaa$yp

cps.tat.cga<-((dn31.trinuc$tat*dn31.trinuc$cga)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tat.cgc<-((dn31.trinuc$tat*dn31.trinuc$cgc)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tat.cgg<-((dn31.trinuc$tat*dn31.trinuc$cgg)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr
cps.tat.cgt<-((dn31.trinuc$tat*dn31.trinuc$cgt)/(dn31.count.aa$y*dn31.count.aa$r))/dn31.count.diaa$yr

cps.tat.cta<-((dn31.trinuc$tat*dn31.trinuc$cta)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tat.ctc<-((dn31.trinuc$tat*dn31.trinuc$ctc)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tat.ctg<-((dn31.trinuc$tat*dn31.trinuc$ctg)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tat.ctt<-((dn31.trinuc$tat*dn31.trinuc$ctt)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl

cps.tat.gaa<-((dn31.trinuc$tat*dn31.trinuc$gaa)/(dn31.count.aa$y*dn31.count.aa$e))/dn31.count.diaa$ye
cps.tat.gac<-((dn31.trinuc$tat*dn31.trinuc$gac)/(dn31.count.aa$y*dn31.count.aa$d))/dn31.count.diaa$yd
cps.tat.gag<-((dn31.trinuc$tat*dn31.trinuc$gag)/(dn31.count.aa$y*dn31.count.aa$e))/dn31.count.diaa$ye
cps.tat.gat<-((dn31.trinuc$tat*dn31.trinuc$gat)/(dn31.count.aa$y*dn31.count.aa$d))/dn31.count.diaa$yd

cps.tat.gca<-((dn31.trinuc$tat*dn31.trinuc$gca)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
cps.tat.gcc<-((dn31.trinuc$tat*dn31.trinuc$gcc)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
cps.tat.gcg<-((dn31.trinuc$tat*dn31.trinuc$gcg)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya
cps.tat.gct<-((dn31.trinuc$tat*dn31.trinuc$gct)/(dn31.count.aa$y*dn31.count.aa$a))/dn31.count.diaa$ya

cps.tat.gga<-((dn31.trinuc$tat*dn31.trinuc$gga)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
cps.tat.ggc<-((dn31.trinuc$tat*dn31.trinuc$ggc)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
cps.tat.ggg<-((dn31.trinuc$tat*dn31.trinuc$ggg)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg
cps.tat.ggt<-((dn31.trinuc$tat*dn31.trinuc$ggt)/(dn31.count.aa$y*dn31.count.aa$g))/dn31.count.diaa$yg

cps.tat.gta<-((dn31.trinuc$tat*dn31.trinuc$gta)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
cps.tat.gtc<-((dn31.trinuc$tat*dn31.trinuc$gtc)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
cps.tat.gtg<-((dn31.trinuc$tat*dn31.trinuc$gtg)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv
cps.tat.gtt<-((dn31.trinuc$tat*dn31.trinuc$gtt)/(dn31.count.aa$y*dn31.count.aa$v))/dn31.count.diaa$yv

#Stop codon
#cps.tat.taa<-((dn31.trinuc$tat*dn31.trinuc$taa)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tat.tac<-((dn31.trinuc$tat*dn31.trinuc$tac)/(dn31.count.aa$y*dn31.count.aa$y))/dn31.count.diaa$yy
#Stop codon
#cps.tat.tag<-((dn31.trinuc$tat*dn31.trinuc$tag)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tat.tat<-((dn31.trinuc$tat*dn31.trinuc$tat)/(dn31.count.aa$y*dn31.count.aa$y))/dn31.count.diaa$yy

cps.tat.tca<-((dn31.trinuc$tat*dn31.trinuc$tca)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
cps.tat.tcc<-((dn31.trinuc$tat*dn31.trinuc$tcc)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
cps.tat.tcg<-((dn31.trinuc$tat*dn31.trinuc$tcg)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys
cps.tat.tct<-((dn31.trinuc$tat*dn31.trinuc$tct)/(dn31.count.aa$y*dn31.count.aa$s))/dn31.count.diaa$ys

#Stop codon
#cps.tat.tga<-((dn31.trinuc$tat*dn31.trinuc$tga)/(dn31.count.aa$y*dn31.count.aa$k))/dn31.count.diaa$yk
cps.tat.tgc<-((dn31.trinuc$tat*dn31.trinuc$tgc)/(dn31.count.aa$y*dn31.count.aa$c))/dn31.count.diaa$yc
cps.tat.tgg<-((dn31.trinuc$tat*dn31.trinuc$tgg)/(dn31.count.aa$y*dn31.count.aa$w))/dn31.count.diaa$yw
cps.tat.tgt<-((dn31.trinuc$tat*dn31.trinuc$tgt)/(dn31.count.aa$y*dn31.count.aa$c))/dn31.count.diaa$yc

cps.tat.tta<-((dn31.trinuc$tat*dn31.trinuc$tta)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tat.ttc<-((dn31.trinuc$tat*dn31.trinuc$ttc)/(dn31.count.aa$y*dn31.count.aa$f))/dn31.count.diaa$yf
cps.tat.ttg<-((dn31.trinuc$tat*dn31.trinuc$ttg)/(dn31.count.aa$y*dn31.count.aa$l))/dn31.count.diaa$yl
cps.tat.ttt<-((dn31.trinuc$tat*dn31.trinuc$ttt)/(dn31.count.aa$y*dn31.count.aa$f))/dn31.count.diaa$yf







cps.tca.aaa<-((dn31.trinuc$tca*dn31.trinuc$aaa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tca.aac<-((dn31.trinuc$tca*dn31.trinuc$aac)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn
cps.tca.aag<-((dn31.trinuc$tca*dn31.trinuc$aag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tca.aat<-((dn31.trinuc$tca*dn31.trinuc$aat)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn

cps.tca.aca<-((dn31.trinuc$tca*dn31.trinuc$aca)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tca.acc<-((dn31.trinuc$tca*dn31.trinuc$acc)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tca.acg<-((dn31.trinuc$tca*dn31.trinuc$acg)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tca.act<-((dn31.trinuc$tca*dn31.trinuc$act)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st

cps.tca.aga<-((dn31.trinuc$tca*dn31.trinuc$aga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tca.agc<-((dn31.trinuc$tca*dn31.trinuc$agc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tca.agg<-((dn31.trinuc$tca*dn31.trinuc$agg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tca.agt<-((dn31.trinuc$tca*dn31.trinuc$agt)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

cps.tca.ata<-((dn31.trinuc$tca*dn31.trinuc$ata)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.tca.atc<-((dn31.trinuc$tca*dn31.trinuc$atc)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.tca.atg<-((dn31.trinuc$tca*dn31.trinuc$atg)/(dn31.count.aa$s*dn31.count.aa$m))/dn31.count.diaa$sm
cps.tca.att<-((dn31.trinuc$tca*dn31.trinuc$att)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si

cps.tca.caa<-((dn31.trinuc$tca*dn31.trinuc$caa)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.tca.cac<-((dn31.trinuc$tca*dn31.trinuc$cac)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh
cps.tca.cag<-((dn31.trinuc$tca*dn31.trinuc$cag)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.tca.cat<-((dn31.trinuc$tca*dn31.trinuc$cat)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh

cps.tca.cca<-((dn31.trinuc$tca*dn31.trinuc$cca)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tca.ccc<-((dn31.trinuc$tca*dn31.trinuc$ccc)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tca.ccg<-((dn31.trinuc$tca*dn31.trinuc$ccg)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tca.cct<-((dn31.trinuc$tca*dn31.trinuc$cct)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp

cps.tca.cga<-((dn31.trinuc$tca*dn31.trinuc$cga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tca.cgc<-((dn31.trinuc$tca*dn31.trinuc$cgc)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tca.cgg<-((dn31.trinuc$tca*dn31.trinuc$cgg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tca.cgt<-((dn31.trinuc$tca*dn31.trinuc$cgt)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr

cps.tca.cta<-((dn31.trinuc$tca*dn31.trinuc$cta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tca.ctc<-((dn31.trinuc$tca*dn31.trinuc$ctc)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tca.ctg<-((dn31.trinuc$tca*dn31.trinuc$ctg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tca.ctt<-((dn31.trinuc$tca*dn31.trinuc$ctt)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl

cps.tca.gaa<-((dn31.trinuc$tca*dn31.trinuc$gaa)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.tca.gac<-((dn31.trinuc$tca*dn31.trinuc$gac)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd
cps.tca.gag<-((dn31.trinuc$tca*dn31.trinuc$gag)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.tca.gat<-((dn31.trinuc$tca*dn31.trinuc$gat)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd

cps.tca.gca<-((dn31.trinuc$tca*dn31.trinuc$gca)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tca.gcc<-((dn31.trinuc$tca*dn31.trinuc$gcc)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tca.gcg<-((dn31.trinuc$tca*dn31.trinuc$gcg)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tca.gct<-((dn31.trinuc$tca*dn31.trinuc$gct)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa

cps.tca.gga<-((dn31.trinuc$tca*dn31.trinuc$gga)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tca.ggc<-((dn31.trinuc$tca*dn31.trinuc$ggc)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tca.ggg<-((dn31.trinuc$tca*dn31.trinuc$ggg)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tca.ggt<-((dn31.trinuc$tca*dn31.trinuc$ggt)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg

cps.tca.gta<-((dn31.trinuc$tca*dn31.trinuc$gta)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tca.gtc<-((dn31.trinuc$tca*dn31.trinuc$gtc)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tca.gtg<-((dn31.trinuc$tca*dn31.trinuc$gtg)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tca.gtt<-((dn31.trinuc$tca*dn31.trinuc$gtt)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv

#Stop codon
#cps.tca.taa<-((dn31.trinuc$tca*dn31.trinuc$taa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tca.tac<-((dn31.trinuc$tca*dn31.trinuc$tac)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy
#Stop codon
#cps.tca.tag<-((dn31.trinuc$tca*dn31.trinuc$tag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tca.tat<-((dn31.trinuc$tca*dn31.trinuc$tat)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy

cps.tca.tca<-((dn31.trinuc$tca*dn31.trinuc$tca)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tca.tcc<-((dn31.trinuc$tca*dn31.trinuc$tcc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tca.tcg<-((dn31.trinuc$tca*dn31.trinuc$tcg)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tca.tct<-((dn31.trinuc$tca*dn31.trinuc$tct)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

#Stop codon
#cps.tca.tga<-((dn31.trinuc$tca*dn31.trinuc$tga)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tca.tgc<-((dn31.trinuc$tca*dn31.trinuc$tgc)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc
cps.tca.tgg<-((dn31.trinuc$tca*dn31.trinuc$tgg)/(dn31.count.aa$s*dn31.count.aa$w))/dn31.count.diaa$sw
cps.tca.tgt<-((dn31.trinuc$tca*dn31.trinuc$tgt)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc

cps.tca.tta<-((dn31.trinuc$tca*dn31.trinuc$tta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tca.ttc<-((dn31.trinuc$tca*dn31.trinuc$ttc)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf
cps.tca.ttg<-((dn31.trinuc$tca*dn31.trinuc$ttg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tca.ttt<-((dn31.trinuc$tca*dn31.trinuc$ttt)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf








cps.tcc.aaa<-((dn31.trinuc$tcc*dn31.trinuc$aaa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcc.aac<-((dn31.trinuc$tcc*dn31.trinuc$aac)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn
cps.tcc.aag<-((dn31.trinuc$tcc*dn31.trinuc$aag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcc.aat<-((dn31.trinuc$tcc*dn31.trinuc$aat)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn

cps.tcc.aca<-((dn31.trinuc$tcc*dn31.trinuc$aca)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tcc.acc<-((dn31.trinuc$tcc*dn31.trinuc$acc)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tcc.acg<-((dn31.trinuc$tcc*dn31.trinuc$acg)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tcc.act<-((dn31.trinuc$tcc*dn31.trinuc$act)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st

cps.tcc.aga<-((dn31.trinuc$tcc*dn31.trinuc$aga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcc.agc<-((dn31.trinuc$tcc*dn31.trinuc$agc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tcc.agg<-((dn31.trinuc$tcc*dn31.trinuc$agg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcc.agt<-((dn31.trinuc$tcc*dn31.trinuc$agt)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

cps.tcc.ata<-((dn31.trinuc$tcc*dn31.trinuc$ata)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.tcc.atc<-((dn31.trinuc$tcc*dn31.trinuc$atc)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.tcc.atg<-((dn31.trinuc$tcc*dn31.trinuc$atg)/(dn31.count.aa$s*dn31.count.aa$m))/dn31.count.diaa$sm
cps.tcc.att<-((dn31.trinuc$tcc*dn31.trinuc$att)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si

cps.tcc.caa<-((dn31.trinuc$tcc*dn31.trinuc$caa)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.tcc.cac<-((dn31.trinuc$tcc*dn31.trinuc$cac)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh
cps.tcc.cag<-((dn31.trinuc$tcc*dn31.trinuc$cag)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.tcc.cat<-((dn31.trinuc$tcc*dn31.trinuc$cat)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh

cps.tcc.cca<-((dn31.trinuc$tcc*dn31.trinuc$cca)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tcc.ccc<-((dn31.trinuc$tcc*dn31.trinuc$ccc)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tcc.ccg<-((dn31.trinuc$tcc*dn31.trinuc$ccg)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tcc.cct<-((dn31.trinuc$tcc*dn31.trinuc$cct)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp

cps.tcc.cga<-((dn31.trinuc$tcc*dn31.trinuc$cga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcc.cgc<-((dn31.trinuc$tcc*dn31.trinuc$cgc)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcc.cgg<-((dn31.trinuc$tcc*dn31.trinuc$cgg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcc.cgt<-((dn31.trinuc$tcc*dn31.trinuc$cgt)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr

cps.tcc.cta<-((dn31.trinuc$tcc*dn31.trinuc$cta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcc.ctc<-((dn31.trinuc$tcc*dn31.trinuc$ctc)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcc.ctg<-((dn31.trinuc$tcc*dn31.trinuc$ctg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcc.ctt<-((dn31.trinuc$tcc*dn31.trinuc$ctt)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl

cps.tcc.gaa<-((dn31.trinuc$tcc*dn31.trinuc$gaa)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.tcc.gac<-((dn31.trinuc$tcc*dn31.trinuc$gac)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd
cps.tcc.gag<-((dn31.trinuc$tcc*dn31.trinuc$gag)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.tcc.gat<-((dn31.trinuc$tcc*dn31.trinuc$gat)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd

cps.tcc.gca<-((dn31.trinuc$tcc*dn31.trinuc$gca)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tcc.gcc<-((dn31.trinuc$tcc*dn31.trinuc$gcc)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tcc.gcg<-((dn31.trinuc$tcc*dn31.trinuc$gcg)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tcc.gct<-((dn31.trinuc$tcc*dn31.trinuc$gct)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa

cps.tcc.gga<-((dn31.trinuc$tcc*dn31.trinuc$gga)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tcc.ggc<-((dn31.trinuc$tcc*dn31.trinuc$ggc)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tcc.ggg<-((dn31.trinuc$tcc*dn31.trinuc$ggg)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tcc.ggt<-((dn31.trinuc$tcc*dn31.trinuc$ggt)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg

cps.tcc.gta<-((dn31.trinuc$tcc*dn31.trinuc$gta)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tcc.gtc<-((dn31.trinuc$tcc*dn31.trinuc$gtc)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tcc.gtg<-((dn31.trinuc$tcc*dn31.trinuc$gtg)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tcc.gtt<-((dn31.trinuc$tcc*dn31.trinuc$gtt)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv

#Stop codon
#cps.tcc.taa<-((dn31.trinuc$tcc*dn31.trinuc$taa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcc.tac<-((dn31.trinuc$tcc*dn31.trinuc$tac)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy
#Stop codon
#cps.tcc.tag<-((dn31.trinuc$tcc*dn31.trinuc$tag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcc.tat<-((dn31.trinuc$tcc*dn31.trinuc$tat)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy

cps.tcc.tca<-((dn31.trinuc$tcc*dn31.trinuc$tca)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tcc.tcc<-((dn31.trinuc$tcc*dn31.trinuc$tcc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tcc.tcg<-((dn31.trinuc$tcc*dn31.trinuc$tcg)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tcc.tct<-((dn31.trinuc$tcc*dn31.trinuc$tct)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

#Stop codon
#cps.tcc.tga<-((dn31.trinuc$tcc*dn31.trinuc$tga)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcc.tgc<-((dn31.trinuc$tcc*dn31.trinuc$tgc)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc
cps.tcc.tgg<-((dn31.trinuc$tcc*dn31.trinuc$tgg)/(dn31.count.aa$s*dn31.count.aa$w))/dn31.count.diaa$sw
cps.tcc.tgt<-((dn31.trinuc$tcc*dn31.trinuc$tgt)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc

cps.tcc.tta<-((dn31.trinuc$tcc*dn31.trinuc$tta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcc.ttc<-((dn31.trinuc$tcc*dn31.trinuc$ttc)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf
cps.tcc.ttg<-((dn31.trinuc$tcc*dn31.trinuc$ttg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcc.ttt<-((dn31.trinuc$tcc*dn31.trinuc$ttt)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf








cps.tcg.aaa<-((dn31.trinuc$tcg*dn31.trinuc$aaa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcg.aac<-((dn31.trinuc$tcg*dn31.trinuc$aac)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn
cps.tcg.aag<-((dn31.trinuc$tcg*dn31.trinuc$aag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcg.aat<-((dn31.trinuc$tcg*dn31.trinuc$aat)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn

cps.tcg.aca<-((dn31.trinuc$tcg*dn31.trinuc$aca)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tcg.acc<-((dn31.trinuc$tcg*dn31.trinuc$acc)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tcg.acg<-((dn31.trinuc$tcg*dn31.trinuc$acg)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tcg.act<-((dn31.trinuc$tcg*dn31.trinuc$act)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st

cps.tcg.aga<-((dn31.trinuc$tcg*dn31.trinuc$aga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcg.agc<-((dn31.trinuc$tcg*dn31.trinuc$agc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tcg.agg<-((dn31.trinuc$tcg*dn31.trinuc$agg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcg.agt<-((dn31.trinuc$tcg*dn31.trinuc$agt)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

cps.tcg.ata<-((dn31.trinuc$tcg*dn31.trinuc$ata)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.tcg.atc<-((dn31.trinuc$tcg*dn31.trinuc$atc)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.tcg.atg<-((dn31.trinuc$tcg*dn31.trinuc$atg)/(dn31.count.aa$s*dn31.count.aa$m))/dn31.count.diaa$sm
cps.tcg.att<-((dn31.trinuc$tcg*dn31.trinuc$att)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si

cps.tcg.caa<-((dn31.trinuc$tcg*dn31.trinuc$caa)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.tcg.cac<-((dn31.trinuc$tcg*dn31.trinuc$cac)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh
cps.tcg.cag<-((dn31.trinuc$tcg*dn31.trinuc$cag)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.tcg.cat<-((dn31.trinuc$tcg*dn31.trinuc$cat)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh

cps.tcg.cca<-((dn31.trinuc$tcg*dn31.trinuc$cca)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tcg.ccc<-((dn31.trinuc$tcg*dn31.trinuc$ccc)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tcg.ccg<-((dn31.trinuc$tcg*dn31.trinuc$ccg)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tcg.cct<-((dn31.trinuc$tcg*dn31.trinuc$cct)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp

cps.tcg.cga<-((dn31.trinuc$tcg*dn31.trinuc$cga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcg.cgc<-((dn31.trinuc$tcg*dn31.trinuc$cgc)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcg.cgg<-((dn31.trinuc$tcg*dn31.trinuc$cgg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tcg.cgt<-((dn31.trinuc$tcg*dn31.trinuc$cgt)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr

cps.tcg.cta<-((dn31.trinuc$tcg*dn31.trinuc$cta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcg.ctc<-((dn31.trinuc$tcg*dn31.trinuc$ctc)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcg.ctg<-((dn31.trinuc$tcg*dn31.trinuc$ctg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcg.ctt<-((dn31.trinuc$tcg*dn31.trinuc$ctt)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl

cps.tcg.gaa<-((dn31.trinuc$tcg*dn31.trinuc$gaa)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.tcg.gac<-((dn31.trinuc$tcg*dn31.trinuc$gac)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd
cps.tcg.gag<-((dn31.trinuc$tcg*dn31.trinuc$gag)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.tcg.gat<-((dn31.trinuc$tcg*dn31.trinuc$gat)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd

cps.tcg.gca<-((dn31.trinuc$tcg*dn31.trinuc$gca)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tcg.gcc<-((dn31.trinuc$tcg*dn31.trinuc$gcc)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tcg.gcg<-((dn31.trinuc$tcg*dn31.trinuc$gcg)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tcg.gct<-((dn31.trinuc$tcg*dn31.trinuc$gct)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa

cps.tcg.gga<-((dn31.trinuc$tcg*dn31.trinuc$gga)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tcg.ggc<-((dn31.trinuc$tcg*dn31.trinuc$ggc)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tcg.ggg<-((dn31.trinuc$tcg*dn31.trinuc$ggg)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tcg.ggt<-((dn31.trinuc$tcg*dn31.trinuc$ggt)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg

cps.tcg.gta<-((dn31.trinuc$tcg*dn31.trinuc$gta)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tcg.gtc<-((dn31.trinuc$tcg*dn31.trinuc$gtc)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tcg.gtg<-((dn31.trinuc$tcg*dn31.trinuc$gtg)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tcg.gtt<-((dn31.trinuc$tcg*dn31.trinuc$gtt)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv

#Stop codon
#cps.tcg.taa<-((dn31.trinuc$tcg*dn31.trinuc$taa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcg.tac<-((dn31.trinuc$tcg*dn31.trinuc$tac)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy
#Stop codon
#cps.tcg.tag<-((dn31.trinuc$tcg*dn31.trinuc$tag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcg.tat<-((dn31.trinuc$tcg*dn31.trinuc$tat)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy

cps.tcg.tca<-((dn31.trinuc$tcg*dn31.trinuc$tca)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tcg.tcc<-((dn31.trinuc$tcg*dn31.trinuc$tcc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tcg.tcg<-((dn31.trinuc$tcg*dn31.trinuc$tcg)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tcg.tct<-((dn31.trinuc$tcg*dn31.trinuc$tct)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

#Stop codon
#cps.tcg.tga<-((dn31.trinuc$tcg*dn31.trinuc$tga)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tcg.tgc<-((dn31.trinuc$tcg*dn31.trinuc$tgc)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc
cps.tcg.tgg<-((dn31.trinuc$tcg*dn31.trinuc$tgg)/(dn31.count.aa$s*dn31.count.aa$w))/dn31.count.diaa$sw
cps.tcg.tgt<-((dn31.trinuc$tcg*dn31.trinuc$tgt)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc

cps.tcg.tta<-((dn31.trinuc$tcg*dn31.trinuc$tta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcg.ttc<-((dn31.trinuc$tcg*dn31.trinuc$ttc)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf
cps.tcg.ttg<-((dn31.trinuc$tcg*dn31.trinuc$ttg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tcg.ttt<-((dn31.trinuc$tcg*dn31.trinuc$ttt)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf








cps.tct.aaa<-((dn31.trinuc$tct*dn31.trinuc$aaa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tct.aac<-((dn31.trinuc$tct*dn31.trinuc$aac)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn
cps.tct.aag<-((dn31.trinuc$tct*dn31.trinuc$aag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tct.aat<-((dn31.trinuc$tct*dn31.trinuc$aat)/(dn31.count.aa$s*dn31.count.aa$n))/dn31.count.diaa$sn

cps.tct.aca<-((dn31.trinuc$tct*dn31.trinuc$aca)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tct.acc<-((dn31.trinuc$tct*dn31.trinuc$acc)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tct.acg<-((dn31.trinuc$tct*dn31.trinuc$acg)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st
cps.tct.act<-((dn31.trinuc$tct*dn31.trinuc$act)/(dn31.count.aa$s*dn31.count.aa$t))/dn31.count.diaa$st

cps.tct.aga<-((dn31.trinuc$tct*dn31.trinuc$aga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tct.agc<-((dn31.trinuc$tct*dn31.trinuc$agc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tct.agg<-((dn31.trinuc$tct*dn31.trinuc$agg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tct.agt<-((dn31.trinuc$tct*dn31.trinuc$agt)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

cps.tct.ata<-((dn31.trinuc$tct*dn31.trinuc$ata)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.tct.atc<-((dn31.trinuc$tct*dn31.trinuc$atc)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si
cps.tct.atg<-((dn31.trinuc$tct*dn31.trinuc$atg)/(dn31.count.aa$s*dn31.count.aa$m))/dn31.count.diaa$sm
cps.tct.att<-((dn31.trinuc$tct*dn31.trinuc$att)/(dn31.count.aa$s*dn31.count.aa$i))/dn31.count.diaa$si

cps.tct.caa<-((dn31.trinuc$tct*dn31.trinuc$caa)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.tct.cac<-((dn31.trinuc$tct*dn31.trinuc$cac)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh
cps.tct.cag<-((dn31.trinuc$tct*dn31.trinuc$cag)/(dn31.count.aa$s*dn31.count.aa$q))/dn31.count.diaa$sq
cps.tct.cat<-((dn31.trinuc$tct*dn31.trinuc$cat)/(dn31.count.aa$s*dn31.count.aa$h))/dn31.count.diaa$sh

cps.tct.cca<-((dn31.trinuc$tct*dn31.trinuc$cca)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tct.ccc<-((dn31.trinuc$tct*dn31.trinuc$ccc)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tct.ccg<-((dn31.trinuc$tct*dn31.trinuc$ccg)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp
cps.tct.cct<-((dn31.trinuc$tct*dn31.trinuc$cct)/(dn31.count.aa$s*dn31.count.aa$p))/dn31.count.diaa$sp

cps.tct.cga<-((dn31.trinuc$tct*dn31.trinuc$cga)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tct.cgc<-((dn31.trinuc$tct*dn31.trinuc$cgc)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tct.cgg<-((dn31.trinuc$tct*dn31.trinuc$cgg)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr
cps.tct.cgt<-((dn31.trinuc$tct*dn31.trinuc$cgt)/(dn31.count.aa$s*dn31.count.aa$r))/dn31.count.diaa$sr

cps.tct.cta<-((dn31.trinuc$tct*dn31.trinuc$cta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tct.ctc<-((dn31.trinuc$tct*dn31.trinuc$ctc)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tct.ctg<-((dn31.trinuc$tct*dn31.trinuc$ctg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tct.ctt<-((dn31.trinuc$tct*dn31.trinuc$ctt)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl

cps.tct.gaa<-((dn31.trinuc$tct*dn31.trinuc$gaa)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.tct.gac<-((dn31.trinuc$tct*dn31.trinuc$gac)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd
cps.tct.gag<-((dn31.trinuc$tct*dn31.trinuc$gag)/(dn31.count.aa$s*dn31.count.aa$e))/dn31.count.diaa$se
cps.tct.gat<-((dn31.trinuc$tct*dn31.trinuc$gat)/(dn31.count.aa$s*dn31.count.aa$d))/dn31.count.diaa$sd

cps.tct.gca<-((dn31.trinuc$tct*dn31.trinuc$gca)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tct.gcc<-((dn31.trinuc$tct*dn31.trinuc$gcc)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tct.gcg<-((dn31.trinuc$tct*dn31.trinuc$gcg)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa
cps.tct.gct<-((dn31.trinuc$tct*dn31.trinuc$gct)/(dn31.count.aa$s*dn31.count.aa$a))/dn31.count.diaa$sa

cps.tct.gga<-((dn31.trinuc$tct*dn31.trinuc$gga)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tct.ggc<-((dn31.trinuc$tct*dn31.trinuc$ggc)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tct.ggg<-((dn31.trinuc$tct*dn31.trinuc$ggg)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg
cps.tct.ggt<-((dn31.trinuc$tct*dn31.trinuc$ggt)/(dn31.count.aa$s*dn31.count.aa$g))/dn31.count.diaa$sg

cps.tct.gta<-((dn31.trinuc$tct*dn31.trinuc$gta)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tct.gtc<-((dn31.trinuc$tct*dn31.trinuc$gtc)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tct.gtg<-((dn31.trinuc$tct*dn31.trinuc$gtg)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv
cps.tct.gtt<-((dn31.trinuc$tct*dn31.trinuc$gtt)/(dn31.count.aa$s*dn31.count.aa$v))/dn31.count.diaa$sv

#Stop codon
#cps.tct.taa<-((dn31.trinuc$tct*dn31.trinuc$taa)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tct.tac<-((dn31.trinuc$tct*dn31.trinuc$tac)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy
#Stop codon
#cps.tct.tag<-((dn31.trinuc$tct*dn31.trinuc$tag)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tct.tat<-((dn31.trinuc$tct*dn31.trinuc$tat)/(dn31.count.aa$s*dn31.count.aa$y))/dn31.count.diaa$sy

cps.tct.tca<-((dn31.trinuc$tct*dn31.trinuc$tca)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tct.tcc<-((dn31.trinuc$tct*dn31.trinuc$tcc)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tct.tcg<-((dn31.trinuc$tct*dn31.trinuc$tcg)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss
cps.tct.tct<-((dn31.trinuc$tct*dn31.trinuc$tct)/(dn31.count.aa$s*dn31.count.aa$s))/dn31.count.diaa$ss

#Stop codon
#cps.tct.tga<-((dn31.trinuc$tct*dn31.trinuc$tga)/(dn31.count.aa$s*dn31.count.aa$k))/dn31.count.diaa$sk
cps.tct.tgc<-((dn31.trinuc$tct*dn31.trinuc$tgc)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc
cps.tct.tgg<-((dn31.trinuc$tct*dn31.trinuc$tgg)/(dn31.count.aa$s*dn31.count.aa$w))/dn31.count.diaa$sw
cps.tct.tgt<-((dn31.trinuc$tct*dn31.trinuc$tgt)/(dn31.count.aa$s*dn31.count.aa$c))/dn31.count.diaa$sc

cps.tct.tta<-((dn31.trinuc$tct*dn31.trinuc$tta)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tct.ttc<-((dn31.trinuc$tct*dn31.trinuc$ttc)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf
cps.tct.ttg<-((dn31.trinuc$tct*dn31.trinuc$ttg)/(dn31.count.aa$s*dn31.count.aa$l))/dn31.count.diaa$sl
cps.tct.ttt<-((dn31.trinuc$tct*dn31.trinuc$ttt)/(dn31.count.aa$s*dn31.count.aa$f))/dn31.count.diaa$sf











#Stop codon
#cps.tga.aaa<-((dn31.trinuc$tga*dn31.trinuc$aaa)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
#cps.tga.aac<-((dn31.trinuc$tga*dn31.trinuc$aac)/(dn31.count.aa$c*dn31.count.aa$n))/dn31.count.diaa$cn
#cps.tga.aag<-((dn31.trinuc$tga*dn31.trinuc$aag)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
#cps.tga.aat<-((dn31.trinuc$tga*dn31.trinuc$aat)/(dn31.count.aa$c*dn31.count.aa$n))/dn31.count.diaa$cn

#cps.tga.aca<-((dn31.trinuc$tga*dn31.trinuc$aca)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
#cps.tga.acc<-((dn31.trinuc$tga*dn31.trinuc$acc)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
#cps.tga.acg<-((dn31.trinuc$tga*dn31.trinuc$acg)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
#cps.tga.act<-((dn31.trinuc$tga*dn31.trinuc$act)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct

#cps.tga.aga<-((dn31.trinuc$tga*dn31.trinuc$aga)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
#cps.tga.agc<-((dn31.trinuc$tga*dn31.trinuc$agc)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
#cps.tga.agg<-((dn31.trinuc$tga*dn31.trinuc$agg)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
#cps.tga.agt<-((dn31.trinuc$tga*dn31.trinuc$agt)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs

#cps.tga.ata<-((dn31.trinuc$tga*dn31.trinuc$ata)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci
#cps.tga.atc<-((dn31.trinuc$tga*dn31.trinuc$atc)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci
#cps.tga.atg<-((dn31.trinuc$tga*dn31.trinuc$atg)/(dn31.count.aa$c*dn31.count.aa$m))/dn31.count.diaa$cm
#cps.tga.att<-((dn31.trinuc$tga*dn31.trinuc$att)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci

#cps.tga.caa<-((dn31.trinuc$tga*dn31.trinuc$caa)/(dn31.count.aa$c*dn31.count.aa$q))/dn31.count.diaa$cq
#cps.tga.cac<-((dn31.trinuc$tga*dn31.trinuc$cac)/(dn31.count.aa$c*dn31.count.aa$h))/dn31.count.diaa$ch
#cps.tga.cag<-((dn31.trinuc$tga*dn31.trinuc$cag)/(dn31.count.aa$c*dn31.count.aa$q))/dn31.count.diaa$cq
#cps.tga.cat<-((dn31.trinuc$tga*dn31.trinuc$cat)/(dn31.count.aa$c*dn31.count.aa$h))/dn31.count.diaa$ch

#cps.tga.cca<-((dn31.trinuc$tga*dn31.trinuc$cca)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
#cps.tga.ccc<-((dn31.trinuc$tga*dn31.trinuc$ccc)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
#cps.tga.ccg<-((dn31.trinuc$tga*dn31.trinuc$ccg)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
#cps.tga.cct<-((dn31.trinuc$tga*dn31.trinuc$cct)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp

#cps.tga.cga<-((dn31.trinuc$tga*dn31.trinuc$cga)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
#cps.tga.cgc<-((dn31.trinuc$tga*dn31.trinuc$cgc)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
#cps.tga.cgg<-((dn31.trinuc$tga*dn31.trinuc$cgg)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
#cps.tga.cgt<-((dn31.trinuc$tga*dn31.trinuc$cgt)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr

#cps.tga.cta<-((dn31.trinuc$tga*dn31.trinuc$cta)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
#cps.tga.ctc<-((dn31.trinuc$tga*dn31.trinuc$ctc)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
#cps.tga.ctg<-((dn31.trinuc$tga*dn31.trinuc$ctg)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
#cps.tga.ctt<-((dn31.trinuc$tga*dn31.trinuc$ctt)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl

#cps.tga.gaa<-((dn31.trinuc$tga*dn31.trinuc$gaa)/(dn31.count.aa$c*dn31.count.aa$e))/dn31.count.diaa$ce
#cps.tga.gac<-((dn31.trinuc$tga*dn31.trinuc$gac)/(dn31.count.aa$c*dn31.count.aa$d))/dn31.count.diaa$cd
#cps.tga.gag<-((dn31.trinuc$tga*dn31.trinuc$gag)/(dn31.count.aa$c*dn31.count.aa$e))/dn31.count.diaa$ce
#cps.tga.gat<-((dn31.trinuc$tga*dn31.trinuc$gat)/(dn31.count.aa$c*dn31.count.aa$d))/dn31.count.diaa$cd

#cps.tga.gca<-((dn31.trinuc$tga*dn31.trinuc$gca)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
#cps.tga.gcc<-((dn31.trinuc$tga*dn31.trinuc$gcc)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
#cps.tga.gcg<-((dn31.trinuc$tga*dn31.trinuc$gcg)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
#cps.tga.gct<-((dn31.trinuc$tga*dn31.trinuc$gct)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca

#cps.tga.gga<-((dn31.trinuc$tga*dn31.trinuc$gga)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
#cps.tga.ggc<-((dn31.trinuc$tga*dn31.trinuc$ggc)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
#cps.tga.ggg<-((dn31.trinuc$tga*dn31.trinuc$ggg)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
#cps.tga.ggt<-((dn31.trinuc$tga*dn31.trinuc$ggt)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg

#cps.tga.gta<-((dn31.trinuc$tga*dn31.trinuc$gta)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
#cps.tga.gtc<-((dn31.trinuc$tga*dn31.trinuc$gtc)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
#cps.tga.gtg<-((dn31.trinuc$tga*dn31.trinuc$gtg)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
#cps.tga.gtt<-((dn31.trinuc$tga*dn31.trinuc$gtt)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv

#Stop codon
#cps.tga.taa<-((dn31.trinuc$tga*dn31.trinuc$taa)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
#cps.tga.tac<-((dn31.trinuc$tga*dn31.trinuc$tac)/(dn31.count.aa$c*dn31.count.aa$y))/dn31.count.diaa$cy
#Stop codon
#cps.tga.tag<-((dn31.trinuc$tga*dn31.trinuc$tag)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
#cps.tga.tat<-((dn31.trinuc$tga*dn31.trinuc$tat)/(dn31.count.aa$c*dn31.count.aa$y))/dn31.count.diaa$cy

#cps.tga.tca<-((dn31.trinuc$tga*dn31.trinuc$tca)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
#cps.tga.tcc<-((dn31.trinuc$tga*dn31.trinuc$tcc)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
#cps.tga.tcg<-((dn31.trinuc$tga*dn31.trinuc$tcg)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
#cps.tga.tct<-((dn31.trinuc$tga*dn31.trinuc$tct)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs

#Stop codon
#cps.tga.tga<-((dn31.trinuc$tga*dn31.trinuc$tga)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
#cps.tga.tgc<-((dn31.trinuc$tga*dn31.trinuc$tgc)/(dn31.count.aa$c*dn31.count.aa$c))/dn31.count.diaa$cc
#cps.tga.tgg<-((dn31.trinuc$tga*dn31.trinuc$tgg)/(dn31.count.aa$c*dn31.count.aa$w))/dn31.count.diaa$cw
#cps.tga.tgt<-((dn31.trinuc$tga*dn31.trinuc$tgt)/(dn31.count.aa$c*dn31.count.aa$c))/dn31.count.diaa$cc

#cps.tga.tta<-((dn31.trinuc$tga*dn31.trinuc$tta)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
#cps.tga.ttc<-((dn31.trinuc$tga*dn31.trinuc$ttc)/(dn31.count.aa$c*dn31.count.aa$f))/dn31.count.diaa$cf
#cps.tga.ttg<-((dn31.trinuc$tga*dn31.trinuc$ttg)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
#cps.tga.ttt<-((dn31.trinuc$tga*dn31.trinuc$ttt)/(dn31.count.aa$c*dn31.count.aa$f))/dn31.count.diaa$cf












cps.tgc.aaa<-((dn31.trinuc$tgc*dn31.trinuc$aaa)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgc.aac<-((dn31.trinuc$tgc*dn31.trinuc$aac)/(dn31.count.aa$c*dn31.count.aa$n))/dn31.count.diaa$cn
cps.tgc.aag<-((dn31.trinuc$tgc*dn31.trinuc$aag)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgc.aat<-((dn31.trinuc$tgc*dn31.trinuc$aat)/(dn31.count.aa$c*dn31.count.aa$n))/dn31.count.diaa$cn

cps.tgc.aca<-((dn31.trinuc$tgc*dn31.trinuc$aca)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
cps.tgc.acc<-((dn31.trinuc$tgc*dn31.trinuc$acc)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
cps.tgc.acg<-((dn31.trinuc$tgc*dn31.trinuc$acg)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
cps.tgc.act<-((dn31.trinuc$tgc*dn31.trinuc$act)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct

cps.tgc.aga<-((dn31.trinuc$tgc*dn31.trinuc$aga)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgc.agc<-((dn31.trinuc$tgc*dn31.trinuc$agc)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
cps.tgc.agg<-((dn31.trinuc$tgc*dn31.trinuc$agg)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgc.agt<-((dn31.trinuc$tgc*dn31.trinuc$agt)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs

cps.tgc.ata<-((dn31.trinuc$tgc*dn31.trinuc$ata)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci
cps.tgc.atc<-((dn31.trinuc$tgc*dn31.trinuc$atc)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci
cps.tgc.atg<-((dn31.trinuc$tgc*dn31.trinuc$atg)/(dn31.count.aa$c*dn31.count.aa$m))/dn31.count.diaa$cm
cps.tgc.att<-((dn31.trinuc$tgc*dn31.trinuc$att)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci

cps.tgc.caa<-((dn31.trinuc$tgc*dn31.trinuc$caa)/(dn31.count.aa$c*dn31.count.aa$q))/dn31.count.diaa$cq
cps.tgc.cac<-((dn31.trinuc$tgc*dn31.trinuc$cac)/(dn31.count.aa$c*dn31.count.aa$h))/dn31.count.diaa$ch
cps.tgc.cag<-((dn31.trinuc$tgc*dn31.trinuc$cag)/(dn31.count.aa$c*dn31.count.aa$q))/dn31.count.diaa$cq
cps.tgc.cat<-((dn31.trinuc$tgc*dn31.trinuc$cat)/(dn31.count.aa$c*dn31.count.aa$h))/dn31.count.diaa$ch

cps.tgc.cca<-((dn31.trinuc$tgc*dn31.trinuc$cca)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
cps.tgc.ccc<-((dn31.trinuc$tgc*dn31.trinuc$ccc)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
cps.tgc.ccg<-((dn31.trinuc$tgc*dn31.trinuc$ccg)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
cps.tgc.cct<-((dn31.trinuc$tgc*dn31.trinuc$cct)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp

cps.tgc.cga<-((dn31.trinuc$tgc*dn31.trinuc$cga)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgc.cgc<-((dn31.trinuc$tgc*dn31.trinuc$cgc)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgc.cgg<-((dn31.trinuc$tgc*dn31.trinuc$cgg)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgc.cgt<-((dn31.trinuc$tgc*dn31.trinuc$cgt)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr

cps.tgc.cta<-((dn31.trinuc$tgc*dn31.trinuc$cta)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgc.ctc<-((dn31.trinuc$tgc*dn31.trinuc$ctc)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgc.ctg<-((dn31.trinuc$tgc*dn31.trinuc$ctg)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgc.ctt<-((dn31.trinuc$tgc*dn31.trinuc$ctt)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl

cps.tgc.gaa<-((dn31.trinuc$tgc*dn31.trinuc$gaa)/(dn31.count.aa$c*dn31.count.aa$e))/dn31.count.diaa$ce
cps.tgc.gac<-((dn31.trinuc$tgc*dn31.trinuc$gac)/(dn31.count.aa$c*dn31.count.aa$d))/dn31.count.diaa$cd
cps.tgc.gag<-((dn31.trinuc$tgc*dn31.trinuc$gag)/(dn31.count.aa$c*dn31.count.aa$e))/dn31.count.diaa$ce
cps.tgc.gat<-((dn31.trinuc$tgc*dn31.trinuc$gat)/(dn31.count.aa$c*dn31.count.aa$d))/dn31.count.diaa$cd

cps.tgc.gca<-((dn31.trinuc$tgc*dn31.trinuc$gca)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
cps.tgc.gcc<-((dn31.trinuc$tgc*dn31.trinuc$gcc)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
cps.tgc.gcg<-((dn31.trinuc$tgc*dn31.trinuc$gcg)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
cps.tgc.gct<-((dn31.trinuc$tgc*dn31.trinuc$gct)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca

cps.tgc.gga<-((dn31.trinuc$tgc*dn31.trinuc$gga)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
cps.tgc.ggc<-((dn31.trinuc$tgc*dn31.trinuc$ggc)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
cps.tgc.ggg<-((dn31.trinuc$tgc*dn31.trinuc$ggg)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
cps.tgc.ggt<-((dn31.trinuc$tgc*dn31.trinuc$ggt)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg

cps.tgc.gta<-((dn31.trinuc$tgc*dn31.trinuc$gta)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
cps.tgc.gtc<-((dn31.trinuc$tgc*dn31.trinuc$gtc)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
cps.tgc.gtg<-((dn31.trinuc$tgc*dn31.trinuc$gtg)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
cps.tgc.gtt<-((dn31.trinuc$tgc*dn31.trinuc$gtt)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv

#Stop codon
#cps.tgc.taa<-((dn31.trinuc$tgc*dn31.trinuc$taa)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgc.tac<-((dn31.trinuc$tgc*dn31.trinuc$tac)/(dn31.count.aa$c*dn31.count.aa$y))/dn31.count.diaa$cy
#Stop codon
#cps.tgc.tag<-((dn31.trinuc$tgc*dn31.trinuc$tag)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgc.tat<-((dn31.trinuc$tgc*dn31.trinuc$tat)/(dn31.count.aa$c*dn31.count.aa$y))/dn31.count.diaa$cy

cps.tgc.tca<-((dn31.trinuc$tgc*dn31.trinuc$tca)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
cps.tgc.tcc<-((dn31.trinuc$tgc*dn31.trinuc$tcc)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
cps.tgc.tcg<-((dn31.trinuc$tgc*dn31.trinuc$tcg)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
cps.tgc.tct<-((dn31.trinuc$tgc*dn31.trinuc$tct)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs

#Stop codon
#cps.tgc.tga<-((dn31.trinuc$tgc*dn31.trinuc$tga)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgc.tgc<-((dn31.trinuc$tgc*dn31.trinuc$tgc)/(dn31.count.aa$c*dn31.count.aa$c))/dn31.count.diaa$cc
cps.tgc.tgg<-((dn31.trinuc$tgc*dn31.trinuc$tgg)/(dn31.count.aa$c*dn31.count.aa$w))/dn31.count.diaa$cw
cps.tgc.tgt<-((dn31.trinuc$tgc*dn31.trinuc$tgt)/(dn31.count.aa$c*dn31.count.aa$c))/dn31.count.diaa$cc

cps.tgc.tta<-((dn31.trinuc$tgc*dn31.trinuc$tta)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgc.ttc<-((dn31.trinuc$tgc*dn31.trinuc$ttc)/(dn31.count.aa$c*dn31.count.aa$f))/dn31.count.diaa$cf
cps.tgc.ttg<-((dn31.trinuc$tgc*dn31.trinuc$ttg)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgc.ttt<-((dn31.trinuc$tgc*dn31.trinuc$ttt)/(dn31.count.aa$c*dn31.count.aa$f))/dn31.count.diaa$cf












cps.tgg.aaa<-((dn31.trinuc$tgg*dn31.trinuc$aaa)/(dn31.count.aa$w*dn31.count.aa$k))/dn31.count.diaa$wk
cps.tgg.aac<-((dn31.trinuc$tgg*dn31.trinuc$aac)/(dn31.count.aa$w*dn31.count.aa$n))/dn31.count.diaa$wn
cps.tgg.aag<-((dn31.trinuc$tgg*dn31.trinuc$aag)/(dn31.count.aa$w*dn31.count.aa$k))/dn31.count.diaa$wk
cps.tgg.aat<-((dn31.trinuc$tgg*dn31.trinuc$aat)/(dn31.count.aa$w*dn31.count.aa$n))/dn31.count.diaa$wn

cps.tgg.aca<-((dn31.trinuc$tgg*dn31.trinuc$aca)/(dn31.count.aa$w*dn31.count.aa$t))/dn31.count.diaa$wt
cps.tgg.acc<-((dn31.trinuc$tgg*dn31.trinuc$acc)/(dn31.count.aa$w*dn31.count.aa$t))/dn31.count.diaa$wt
cps.tgg.acg<-((dn31.trinuc$tgg*dn31.trinuc$acg)/(dn31.count.aa$w*dn31.count.aa$t))/dn31.count.diaa$wt
cps.tgg.act<-((dn31.trinuc$tgg*dn31.trinuc$act)/(dn31.count.aa$w*dn31.count.aa$t))/dn31.count.diaa$wt

cps.tgg.aga<-((dn31.trinuc$tgg*dn31.trinuc$aga)/(dn31.count.aa$w*dn31.count.aa$r))/dn31.count.diaa$wr
cps.tgg.agc<-((dn31.trinuc$tgg*dn31.trinuc$agc)/(dn31.count.aa$w*dn31.count.aa$s))/dn31.count.diaa$ws
cps.tgg.agg<-((dn31.trinuc$tgg*dn31.trinuc$agg)/(dn31.count.aa$w*dn31.count.aa$r))/dn31.count.diaa$wr
cps.tgg.agt<-((dn31.trinuc$tgg*dn31.trinuc$agt)/(dn31.count.aa$w*dn31.count.aa$s))/dn31.count.diaa$ws

cps.tgg.ata<-((dn31.trinuc$tgg*dn31.trinuc$ata)/(dn31.count.aa$w*dn31.count.aa$i))/dn31.count.diaa$wi
cps.tgg.atc<-((dn31.trinuc$tgg*dn31.trinuc$atc)/(dn31.count.aa$w*dn31.count.aa$i))/dn31.count.diaa$wi
cps.tgg.atg<-((dn31.trinuc$tgg*dn31.trinuc$atg)/(dn31.count.aa$w*dn31.count.aa$m))/dn31.count.diaa$wm
cps.tgg.att<-((dn31.trinuc$tgg*dn31.trinuc$att)/(dn31.count.aa$w*dn31.count.aa$i))/dn31.count.diaa$wi

cps.tgg.caa<-((dn31.trinuc$tgg*dn31.trinuc$caa)/(dn31.count.aa$w*dn31.count.aa$q))/dn31.count.diaa$wq
cps.tgg.cac<-((dn31.trinuc$tgg*dn31.trinuc$cac)/(dn31.count.aa$w*dn31.count.aa$h))/dn31.count.diaa$wh
cps.tgg.cag<-((dn31.trinuc$tgg*dn31.trinuc$cag)/(dn31.count.aa$w*dn31.count.aa$q))/dn31.count.diaa$wq
cps.tgg.cat<-((dn31.trinuc$tgg*dn31.trinuc$cat)/(dn31.count.aa$w*dn31.count.aa$h))/dn31.count.diaa$wh

cps.tgg.cca<-((dn31.trinuc$tgg*dn31.trinuc$cca)/(dn31.count.aa$w*dn31.count.aa$p))/dn31.count.diaa$wp
cps.tgg.ccc<-((dn31.trinuc$tgg*dn31.trinuc$ccc)/(dn31.count.aa$w*dn31.count.aa$p))/dn31.count.diaa$wp
cps.tgg.ccg<-((dn31.trinuc$tgg*dn31.trinuc$ccg)/(dn31.count.aa$w*dn31.count.aa$p))/dn31.count.diaa$wp
cps.tgg.cct<-((dn31.trinuc$tgg*dn31.trinuc$cct)/(dn31.count.aa$w*dn31.count.aa$p))/dn31.count.diaa$wp

cps.tgg.cga<-((dn31.trinuc$tgg*dn31.trinuc$cga)/(dn31.count.aa$w*dn31.count.aa$r))/dn31.count.diaa$wr
cps.tgg.cgc<-((dn31.trinuc$tgg*dn31.trinuc$cgc)/(dn31.count.aa$w*dn31.count.aa$r))/dn31.count.diaa$wr
cps.tgg.cgg<-((dn31.trinuc$tgg*dn31.trinuc$cgg)/(dn31.count.aa$w*dn31.count.aa$r))/dn31.count.diaa$wr
cps.tgg.cgt<-((dn31.trinuc$tgg*dn31.trinuc$cgt)/(dn31.count.aa$w*dn31.count.aa$r))/dn31.count.diaa$wr

cps.tgg.cta<-((dn31.trinuc$tgg*dn31.trinuc$cta)/(dn31.count.aa$w*dn31.count.aa$l))/dn31.count.diaa$wl
cps.tgg.ctc<-((dn31.trinuc$tgg*dn31.trinuc$ctc)/(dn31.count.aa$w*dn31.count.aa$l))/dn31.count.diaa$wl
cps.tgg.ctg<-((dn31.trinuc$tgg*dn31.trinuc$ctg)/(dn31.count.aa$w*dn31.count.aa$l))/dn31.count.diaa$wl
cps.tgg.ctt<-((dn31.trinuc$tgg*dn31.trinuc$ctt)/(dn31.count.aa$w*dn31.count.aa$l))/dn31.count.diaa$wl

cps.tgg.gaa<-((dn31.trinuc$tgg*dn31.trinuc$gaa)/(dn31.count.aa$w*dn31.count.aa$e))/dn31.count.diaa$we
cps.tgg.gac<-((dn31.trinuc$tgg*dn31.trinuc$gac)/(dn31.count.aa$w*dn31.count.aa$d))/dn31.count.diaa$wd
cps.tgg.gag<-((dn31.trinuc$tgg*dn31.trinuc$gag)/(dn31.count.aa$w*dn31.count.aa$e))/dn31.count.diaa$we
cps.tgg.gat<-((dn31.trinuc$tgg*dn31.trinuc$gat)/(dn31.count.aa$w*dn31.count.aa$d))/dn31.count.diaa$wd

cps.tgg.gca<-((dn31.trinuc$tgg*dn31.trinuc$gca)/(dn31.count.aa$w*dn31.count.aa$a))/dn31.count.diaa$wa
cps.tgg.gcc<-((dn31.trinuc$tgg*dn31.trinuc$gcc)/(dn31.count.aa$w*dn31.count.aa$a))/dn31.count.diaa$wa
cps.tgg.gcg<-((dn31.trinuc$tgg*dn31.trinuc$gcg)/(dn31.count.aa$w*dn31.count.aa$a))/dn31.count.diaa$wa
cps.tgg.gct<-((dn31.trinuc$tgg*dn31.trinuc$gct)/(dn31.count.aa$w*dn31.count.aa$a))/dn31.count.diaa$wa

cps.tgg.gga<-((dn31.trinuc$tgg*dn31.trinuc$gga)/(dn31.count.aa$w*dn31.count.aa$g))/dn31.count.diaa$wg
cps.tgg.ggc<-((dn31.trinuc$tgg*dn31.trinuc$ggc)/(dn31.count.aa$w*dn31.count.aa$g))/dn31.count.diaa$wg
cps.tgg.ggg<-((dn31.trinuc$tgg*dn31.trinuc$ggg)/(dn31.count.aa$w*dn31.count.aa$g))/dn31.count.diaa$wg
cps.tgg.ggt<-((dn31.trinuc$tgg*dn31.trinuc$ggt)/(dn31.count.aa$w*dn31.count.aa$g))/dn31.count.diaa$wg

cps.tgg.gta<-((dn31.trinuc$tgg*dn31.trinuc$gta)/(dn31.count.aa$w*dn31.count.aa$v))/dn31.count.diaa$wv
cps.tgg.gtc<-((dn31.trinuc$tgg*dn31.trinuc$gtc)/(dn31.count.aa$w*dn31.count.aa$v))/dn31.count.diaa$wv
cps.tgg.gtg<-((dn31.trinuc$tgg*dn31.trinuc$gtg)/(dn31.count.aa$w*dn31.count.aa$v))/dn31.count.diaa$wv
cps.tgg.gtt<-((dn31.trinuc$tgg*dn31.trinuc$gtt)/(dn31.count.aa$w*dn31.count.aa$v))/dn31.count.diaa$wv

#Stop codon
#cps.tgg.taa<-((dn31.trinuc$tgg*dn31.trinuc$taa)/(dn31.count.aa$w*dn31.count.aa$k))/dn31.count.diaa$wk
cps.tgg.tac<-((dn31.trinuc$tgg*dn31.trinuc$tac)/(dn31.count.aa$w*dn31.count.aa$y))/dn31.count.diaa$wy
#Stop codon
#cps.tgg.tag<-((dn31.trinuc$tgg*dn31.trinuc$tag)/(dn31.count.aa$w*dn31.count.aa$k))/dn31.count.diaa$wk
cps.tgg.tat<-((dn31.trinuc$tgg*dn31.trinuc$tat)/(dn31.count.aa$w*dn31.count.aa$y))/dn31.count.diaa$wy

cps.tgg.tca<-((dn31.trinuc$tgg*dn31.trinuc$tca)/(dn31.count.aa$w*dn31.count.aa$s))/dn31.count.diaa$ws
cps.tgg.tcc<-((dn31.trinuc$tgg*dn31.trinuc$tcc)/(dn31.count.aa$w*dn31.count.aa$s))/dn31.count.diaa$ws
cps.tgg.tcg<-((dn31.trinuc$tgg*dn31.trinuc$tcg)/(dn31.count.aa$w*dn31.count.aa$s))/dn31.count.diaa$ws
cps.tgg.tct<-((dn31.trinuc$tgg*dn31.trinuc$tct)/(dn31.count.aa$w*dn31.count.aa$s))/dn31.count.diaa$ws

#Stop codon
#cps.tgg.tga<-((dn31.trinuc$tgg*dn31.trinuc$tga)/(dn31.count.aa$w*dn31.count.aa$k))/dn31.count.diaa$wk
cps.tgg.tgc<-((dn31.trinuc$tgg*dn31.trinuc$tgc)/(dn31.count.aa$w*dn31.count.aa$c))/dn31.count.diaa$wc
cps.tgg.tgg<-((dn31.trinuc$tgg*dn31.trinuc$tgg)/(dn31.count.aa$w*dn31.count.aa$w))/dn31.count.diaa$ww
cps.tgg.tgt<-((dn31.trinuc$tgg*dn31.trinuc$tgt)/(dn31.count.aa$w*dn31.count.aa$c))/dn31.count.diaa$wc

cps.tgg.tta<-((dn31.trinuc$tgg*dn31.trinuc$tta)/(dn31.count.aa$w*dn31.count.aa$l))/dn31.count.diaa$wl
cps.tgg.ttc<-((dn31.trinuc$tgg*dn31.trinuc$ttc)/(dn31.count.aa$w*dn31.count.aa$f))/dn31.count.diaa$wf
cps.tgg.ttg<-((dn31.trinuc$tgg*dn31.trinuc$ttg)/(dn31.count.aa$w*dn31.count.aa$l))/dn31.count.diaa$wl
cps.tgg.ttt<-((dn31.trinuc$tgg*dn31.trinuc$ttt)/(dn31.count.aa$w*dn31.count.aa$f))/dn31.count.diaa$wf








cps.tgt.aaa<-((dn31.trinuc$tgt*dn31.trinuc$aaa)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgt.aac<-((dn31.trinuc$tgt*dn31.trinuc$aac)/(dn31.count.aa$c*dn31.count.aa$n))/dn31.count.diaa$cn
cps.tgt.aag<-((dn31.trinuc$tgt*dn31.trinuc$aag)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgt.aat<-((dn31.trinuc$tgt*dn31.trinuc$aat)/(dn31.count.aa$c*dn31.count.aa$n))/dn31.count.diaa$cn

cps.tgt.aca<-((dn31.trinuc$tgt*dn31.trinuc$aca)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
cps.tgt.acc<-((dn31.trinuc$tgt*dn31.trinuc$acc)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
cps.tgt.acg<-((dn31.trinuc$tgt*dn31.trinuc$acg)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct
cps.tgt.act<-((dn31.trinuc$tgt*dn31.trinuc$act)/(dn31.count.aa$c*dn31.count.aa$t))/dn31.count.diaa$ct

cps.tgt.aga<-((dn31.trinuc$tgt*dn31.trinuc$aga)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgt.agc<-((dn31.trinuc$tgt*dn31.trinuc$agc)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
cps.tgt.agg<-((dn31.trinuc$tgt*dn31.trinuc$agg)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgt.agt<-((dn31.trinuc$tgt*dn31.trinuc$agt)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs

cps.tgt.ata<-((dn31.trinuc$tgt*dn31.trinuc$ata)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci
cps.tgt.atc<-((dn31.trinuc$tgt*dn31.trinuc$atc)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci
cps.tgt.atg<-((dn31.trinuc$tgt*dn31.trinuc$atg)/(dn31.count.aa$c*dn31.count.aa$m))/dn31.count.diaa$cm
cps.tgt.att<-((dn31.trinuc$tgt*dn31.trinuc$att)/(dn31.count.aa$c*dn31.count.aa$i))/dn31.count.diaa$ci

cps.tgt.caa<-((dn31.trinuc$tgt*dn31.trinuc$caa)/(dn31.count.aa$c*dn31.count.aa$q))/dn31.count.diaa$cq
cps.tgt.cac<-((dn31.trinuc$tgt*dn31.trinuc$cac)/(dn31.count.aa$c*dn31.count.aa$h))/dn31.count.diaa$ch
cps.tgt.cag<-((dn31.trinuc$tgt*dn31.trinuc$cag)/(dn31.count.aa$c*dn31.count.aa$q))/dn31.count.diaa$cq
cps.tgt.cat<-((dn31.trinuc$tgt*dn31.trinuc$cat)/(dn31.count.aa$c*dn31.count.aa$h))/dn31.count.diaa$ch

cps.tgt.cca<-((dn31.trinuc$tgt*dn31.trinuc$cca)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
cps.tgt.ccc<-((dn31.trinuc$tgt*dn31.trinuc$ccc)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
cps.tgt.ccg<-((dn31.trinuc$tgt*dn31.trinuc$ccg)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp
cps.tgt.cct<-((dn31.trinuc$tgt*dn31.trinuc$cct)/(dn31.count.aa$c*dn31.count.aa$p))/dn31.count.diaa$cp

cps.tgt.cga<-((dn31.trinuc$tgt*dn31.trinuc$cga)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgt.cgc<-((dn31.trinuc$tgt*dn31.trinuc$cgc)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgt.cgg<-((dn31.trinuc$tgt*dn31.trinuc$cgg)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr
cps.tgt.cgt<-((dn31.trinuc$tgt*dn31.trinuc$cgt)/(dn31.count.aa$c*dn31.count.aa$r))/dn31.count.diaa$cr

cps.tgt.cta<-((dn31.trinuc$tgt*dn31.trinuc$cta)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgt.ctc<-((dn31.trinuc$tgt*dn31.trinuc$ctc)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgt.ctg<-((dn31.trinuc$tgt*dn31.trinuc$ctg)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgt.ctt<-((dn31.trinuc$tgt*dn31.trinuc$ctt)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl

cps.tgt.gaa<-((dn31.trinuc$tgt*dn31.trinuc$gaa)/(dn31.count.aa$c*dn31.count.aa$e))/dn31.count.diaa$ce
cps.tgt.gac<-((dn31.trinuc$tgt*dn31.trinuc$gac)/(dn31.count.aa$c*dn31.count.aa$d))/dn31.count.diaa$cd
cps.tgt.gag<-((dn31.trinuc$tgt*dn31.trinuc$gag)/(dn31.count.aa$c*dn31.count.aa$e))/dn31.count.diaa$ce
cps.tgt.gat<-((dn31.trinuc$tgt*dn31.trinuc$gat)/(dn31.count.aa$c*dn31.count.aa$d))/dn31.count.diaa$cd

cps.tgt.gca<-((dn31.trinuc$tgt*dn31.trinuc$gca)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
cps.tgt.gcc<-((dn31.trinuc$tgt*dn31.trinuc$gcc)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
cps.tgt.gcg<-((dn31.trinuc$tgt*dn31.trinuc$gcg)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca
cps.tgt.gct<-((dn31.trinuc$tgt*dn31.trinuc$gct)/(dn31.count.aa$c*dn31.count.aa$a))/dn31.count.diaa$ca

cps.tgt.gga<-((dn31.trinuc$tgt*dn31.trinuc$gga)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
cps.tgt.ggc<-((dn31.trinuc$tgt*dn31.trinuc$ggc)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
cps.tgt.ggg<-((dn31.trinuc$tgt*dn31.trinuc$ggg)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg
cps.tgt.ggt<-((dn31.trinuc$tgt*dn31.trinuc$ggt)/(dn31.count.aa$c*dn31.count.aa$g))/dn31.count.diaa$cg

cps.tgt.gta<-((dn31.trinuc$tgt*dn31.trinuc$gta)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
cps.tgt.gtc<-((dn31.trinuc$tgt*dn31.trinuc$gtc)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
cps.tgt.gtg<-((dn31.trinuc$tgt*dn31.trinuc$gtg)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv
cps.tgt.gtt<-((dn31.trinuc$tgt*dn31.trinuc$gtt)/(dn31.count.aa$c*dn31.count.aa$v))/dn31.count.diaa$cv

#Stop codon
#cps.tgt.taa<-((dn31.trinuc$tgt*dn31.trinuc$taa)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgt.tac<-((dn31.trinuc$tgt*dn31.trinuc$tac)/(dn31.count.aa$c*dn31.count.aa$y))/dn31.count.diaa$cy
#Stop codon
#cps.tgt.tag<-((dn31.trinuc$tgt*dn31.trinuc$tag)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgt.tat<-((dn31.trinuc$tgt*dn31.trinuc$tat)/(dn31.count.aa$c*dn31.count.aa$y))/dn31.count.diaa$cy

cps.tgt.tca<-((dn31.trinuc$tgt*dn31.trinuc$tca)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
cps.tgt.tcc<-((dn31.trinuc$tgt*dn31.trinuc$tcc)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
cps.tgt.tcg<-((dn31.trinuc$tgt*dn31.trinuc$tcg)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs
cps.tgt.tct<-((dn31.trinuc$tgt*dn31.trinuc$tct)/(dn31.count.aa$c*dn31.count.aa$s))/dn31.count.diaa$cs

#Stop codon
#cps.tgt.tga<-((dn31.trinuc$tgt*dn31.trinuc$tga)/(dn31.count.aa$c*dn31.count.aa$k))/dn31.count.diaa$ck
cps.tgt.tgc<-((dn31.trinuc$tgt*dn31.trinuc$tgc)/(dn31.count.aa$c*dn31.count.aa$c))/dn31.count.diaa$cc
cps.tgt.tgg<-((dn31.trinuc$tgt*dn31.trinuc$tgg)/(dn31.count.aa$c*dn31.count.aa$w))/dn31.count.diaa$cw
cps.tgt.tgt<-((dn31.trinuc$tgt*dn31.trinuc$tgt)/(dn31.count.aa$c*dn31.count.aa$c))/dn31.count.diaa$cc

cps.tgt.tta<-((dn31.trinuc$tgt*dn31.trinuc$tta)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgt.ttc<-((dn31.trinuc$tgt*dn31.trinuc$ttc)/(dn31.count.aa$c*dn31.count.aa$f))/dn31.count.diaa$cf
cps.tgt.ttg<-((dn31.trinuc$tgt*dn31.trinuc$ttg)/(dn31.count.aa$c*dn31.count.aa$l))/dn31.count.diaa$cl
cps.tgt.ttt<-((dn31.trinuc$tgt*dn31.trinuc$ttt)/(dn31.count.aa$c*dn31.count.aa$f))/dn31.count.diaa$cf








cps.tta.aaa<-((dn31.trinuc$tta*dn31.trinuc$aaa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.tta.aac<-((dn31.trinuc$tta*dn31.trinuc$aac)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln
cps.tta.aag<-((dn31.trinuc$tta*dn31.trinuc$aag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.tta.aat<-((dn31.trinuc$tta*dn31.trinuc$aat)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln

cps.tta.aca<-((dn31.trinuc$tta*dn31.trinuc$aca)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.tta.acc<-((dn31.trinuc$tta*dn31.trinuc$acc)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.tta.acg<-((dn31.trinuc$tta*dn31.trinuc$acg)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.tta.act<-((dn31.trinuc$tta*dn31.trinuc$act)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt

cps.tta.aga<-((dn31.trinuc$tta*dn31.trinuc$aga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.tta.agc<-((dn31.trinuc$tta*dn31.trinuc$agc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.tta.agg<-((dn31.trinuc$tta*dn31.trinuc$agg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.tta.agt<-((dn31.trinuc$tta*dn31.trinuc$agt)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

cps.tta.ata<-((dn31.trinuc$tta*dn31.trinuc$ata)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.tta.atc<-((dn31.trinuc$tta*dn31.trinuc$atc)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.tta.atg<-((dn31.trinuc$tta*dn31.trinuc$atg)/(dn31.count.aa$l*dn31.count.aa$m))/dn31.count.diaa$lm
cps.tta.att<-((dn31.trinuc$tta*dn31.trinuc$att)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li

cps.tta.caa<-((dn31.trinuc$tta*dn31.trinuc$caa)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.tta.cac<-((dn31.trinuc$tta*dn31.trinuc$cac)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh
cps.tta.cag<-((dn31.trinuc$tta*dn31.trinuc$cag)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.tta.cat<-((dn31.trinuc$tta*dn31.trinuc$cat)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh

cps.tta.cca<-((dn31.trinuc$tta*dn31.trinuc$cca)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.tta.ccc<-((dn31.trinuc$tta*dn31.trinuc$ccc)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.tta.ccg<-((dn31.trinuc$tta*dn31.trinuc$ccg)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.tta.cct<-((dn31.trinuc$tta*dn31.trinuc$cct)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp

cps.tta.cga<-((dn31.trinuc$tta*dn31.trinuc$cga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.tta.cgc<-((dn31.trinuc$tta*dn31.trinuc$cgc)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.tta.cgg<-((dn31.trinuc$tta*dn31.trinuc$cgg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.tta.cgt<-((dn31.trinuc$tta*dn31.trinuc$cgt)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr

cps.tta.cta<-((dn31.trinuc$tta*dn31.trinuc$cta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.tta.ctc<-((dn31.trinuc$tta*dn31.trinuc$ctc)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.tta.ctg<-((dn31.trinuc$tta*dn31.trinuc$ctg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.tta.ctt<-((dn31.trinuc$tta*dn31.trinuc$ctt)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll

cps.tta.gaa<-((dn31.trinuc$tta*dn31.trinuc$gaa)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.tta.gac<-((dn31.trinuc$tta*dn31.trinuc$gac)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld
cps.tta.gag<-((dn31.trinuc$tta*dn31.trinuc$gag)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.tta.gat<-((dn31.trinuc$tta*dn31.trinuc$gat)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld

cps.tta.gca<-((dn31.trinuc$tta*dn31.trinuc$gca)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.tta.gcc<-((dn31.trinuc$tta*dn31.trinuc$gcc)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.tta.gcg<-((dn31.trinuc$tta*dn31.trinuc$gcg)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.tta.gct<-((dn31.trinuc$tta*dn31.trinuc$gct)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la

cps.tta.gga<-((dn31.trinuc$tta*dn31.trinuc$gga)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.tta.ggc<-((dn31.trinuc$tta*dn31.trinuc$ggc)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.tta.ggg<-((dn31.trinuc$tta*dn31.trinuc$ggg)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.tta.ggt<-((dn31.trinuc$tta*dn31.trinuc$ggt)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg

cps.tta.gta<-((dn31.trinuc$tta*dn31.trinuc$gta)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.tta.gtc<-((dn31.trinuc$tta*dn31.trinuc$gtc)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.tta.gtg<-((dn31.trinuc$tta*dn31.trinuc$gtg)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.tta.gtt<-((dn31.trinuc$tta*dn31.trinuc$gtt)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv

#Stop codon
#cps.tta.taa<-((dn31.trinuc$tta*dn31.trinuc$taa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.tta.tac<-((dn31.trinuc$tta*dn31.trinuc$tac)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly
#Stop codon
#cps.tta.tag<-((dn31.trinuc$tta*dn31.trinuc$tag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.tta.tat<-((dn31.trinuc$tta*dn31.trinuc$tat)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly

cps.tta.tca<-((dn31.trinuc$tta*dn31.trinuc$tca)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.tta.tcc<-((dn31.trinuc$tta*dn31.trinuc$tcc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.tta.tcg<-((dn31.trinuc$tta*dn31.trinuc$tcg)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.tta.tct<-((dn31.trinuc$tta*dn31.trinuc$tct)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

#Stop codon
#cps.tta.tga<-((dn31.trinuc$tta*dn31.trinuc$tga)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.tta.tgc<-((dn31.trinuc$tta*dn31.trinuc$tgc)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc
cps.tta.tgg<-((dn31.trinuc$tta*dn31.trinuc$tgg)/(dn31.count.aa$l*dn31.count.aa$w))/dn31.count.diaa$lw
cps.tta.tgt<-((dn31.trinuc$tta*dn31.trinuc$tgt)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc

cps.tta.tta<-((dn31.trinuc$tta*dn31.trinuc$tta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.tta.ttc<-((dn31.trinuc$tta*dn31.trinuc$ttc)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf
cps.tta.ttg<-((dn31.trinuc$tta*dn31.trinuc$ttg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.tta.ttt<-((dn31.trinuc$tta*dn31.trinuc$ttt)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf








cps.ttc.aaa<-((dn31.trinuc$ttc*dn31.trinuc$aaa)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttc.aac<-((dn31.trinuc$ttc*dn31.trinuc$aac)/(dn31.count.aa$f*dn31.count.aa$n))/dn31.count.diaa$fn
cps.ttc.aag<-((dn31.trinuc$ttc*dn31.trinuc$aag)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttc.aat<-((dn31.trinuc$ttc*dn31.trinuc$aat)/(dn31.count.aa$f*dn31.count.aa$n))/dn31.count.diaa$fn

cps.ttc.aca<-((dn31.trinuc$ttc*dn31.trinuc$aca)/(dn31.count.aa$f*dn31.count.aa$t))/dn31.count.diaa$ft
cps.ttc.acc<-((dn31.trinuc$ttc*dn31.trinuc$acc)/(dn31.count.aa$f*dn31.count.aa$t))/dn31.count.diaa$ft
cps.ttc.acg<-((dn31.trinuc$ttc*dn31.trinuc$acg)/(dn31.count.aa$f*dn31.count.aa$t))/dn31.count.diaa$ft
cps.ttc.act<-((dn31.trinuc$ttc*dn31.trinuc$act)/(dn31.count.aa$f*dn31.count.aa$t))/dn31.count.diaa$ft

cps.ttc.aga<-((dn31.trinuc$ttc*dn31.trinuc$aga)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttc.agc<-((dn31.trinuc$ttc*dn31.trinuc$agc)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs
cps.ttc.agg<-((dn31.trinuc$ttc*dn31.trinuc$agg)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttc.agt<-((dn31.trinuc$ttc*dn31.trinuc$agt)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs

cps.ttc.ata<-((dn31.trinuc$ttc*dn31.trinuc$ata)/(dn31.count.aa$f*dn31.count.aa$i))/dn31.count.diaa$fi
cps.ttc.atc<-((dn31.trinuc$ttc*dn31.trinuc$atc)/(dn31.count.aa$f*dn31.count.aa$i))/dn31.count.diaa$fi
cps.ttc.atg<-((dn31.trinuc$ttc*dn31.trinuc$atg)/(dn31.count.aa$f*dn31.count.aa$m))/dn31.count.diaa$fm
cps.ttc.att<-((dn31.trinuc$ttc*dn31.trinuc$att)/(dn31.count.aa$f*dn31.count.aa$i))/dn31.count.diaa$fi

cps.ttc.caa<-((dn31.trinuc$ttc*dn31.trinuc$caa)/(dn31.count.aa$f*dn31.count.aa$q))/dn31.count.diaa$fq
cps.ttc.cac<-((dn31.trinuc$ttc*dn31.trinuc$cac)/(dn31.count.aa$f*dn31.count.aa$h))/dn31.count.diaa$fh
cps.ttc.cag<-((dn31.trinuc$ttc*dn31.trinuc$cag)/(dn31.count.aa$f*dn31.count.aa$q))/dn31.count.diaa$fq
cps.ttc.cat<-((dn31.trinuc$ttc*dn31.trinuc$cat)/(dn31.count.aa$f*dn31.count.aa$h))/dn31.count.diaa$fh

cps.ttc.cca<-((dn31.trinuc$ttc*dn31.trinuc$cca)/(dn31.count.aa$f*dn31.count.aa$p))/dn31.count.diaa$fp
cps.ttc.ccc<-((dn31.trinuc$ttc*dn31.trinuc$ccc)/(dn31.count.aa$f*dn31.count.aa$p))/dn31.count.diaa$fp
cps.ttc.ccg<-((dn31.trinuc$ttc*dn31.trinuc$ccg)/(dn31.count.aa$f*dn31.count.aa$p))/dn31.count.diaa$fp
cps.ttc.cct<-((dn31.trinuc$ttc*dn31.trinuc$cct)/(dn31.count.aa$f*dn31.count.aa$p))/dn31.count.diaa$fp

cps.ttc.cga<-((dn31.trinuc$ttc*dn31.trinuc$cga)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttc.cgc<-((dn31.trinuc$ttc*dn31.trinuc$cgc)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttc.cgg<-((dn31.trinuc$ttc*dn31.trinuc$cgg)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttc.cgt<-((dn31.trinuc$ttc*dn31.trinuc$cgt)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr

cps.ttc.cta<-((dn31.trinuc$ttc*dn31.trinuc$cta)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttc.ctc<-((dn31.trinuc$ttc*dn31.trinuc$ctc)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttc.ctg<-((dn31.trinuc$ttc*dn31.trinuc$ctg)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttc.ctt<-((dn31.trinuc$ttc*dn31.trinuc$ctt)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl

cps.ttc.gaa<-((dn31.trinuc$ttc*dn31.trinuc$gaa)/(dn31.count.aa$f*dn31.count.aa$e))/dn31.count.diaa$fe
cps.ttc.gac<-((dn31.trinuc$ttc*dn31.trinuc$gac)/(dn31.count.aa$f*dn31.count.aa$d))/dn31.count.diaa$fd
cps.ttc.gag<-((dn31.trinuc$ttc*dn31.trinuc$gag)/(dn31.count.aa$f*dn31.count.aa$e))/dn31.count.diaa$fe
cps.ttc.gat<-((dn31.trinuc$ttc*dn31.trinuc$gat)/(dn31.count.aa$f*dn31.count.aa$d))/dn31.count.diaa$fd

cps.ttc.gca<-((dn31.trinuc$ttc*dn31.trinuc$gca)/(dn31.count.aa$f*dn31.count.aa$a))/dn31.count.diaa$fa
cps.ttc.gcc<-((dn31.trinuc$ttc*dn31.trinuc$gcc)/(dn31.count.aa$f*dn31.count.aa$a))/dn31.count.diaa$fa
cps.ttc.gcg<-((dn31.trinuc$ttc*dn31.trinuc$gcg)/(dn31.count.aa$f*dn31.count.aa$a))/dn31.count.diaa$fa
cps.ttc.gct<-((dn31.trinuc$ttc*dn31.trinuc$gct)/(dn31.count.aa$f*dn31.count.aa$a))/dn31.count.diaa$fa

cps.ttc.gga<-((dn31.trinuc$ttc*dn31.trinuc$gga)/(dn31.count.aa$f*dn31.count.aa$g))/dn31.count.diaa$fg
cps.ttc.ggc<-((dn31.trinuc$ttc*dn31.trinuc$ggc)/(dn31.count.aa$f*dn31.count.aa$g))/dn31.count.diaa$fg
cps.ttc.ggg<-((dn31.trinuc$ttc*dn31.trinuc$ggg)/(dn31.count.aa$f*dn31.count.aa$g))/dn31.count.diaa$fg
cps.ttc.ggt<-((dn31.trinuc$ttc*dn31.trinuc$ggt)/(dn31.count.aa$f*dn31.count.aa$g))/dn31.count.diaa$fg

cps.ttc.gta<-((dn31.trinuc$ttc*dn31.trinuc$gta)/(dn31.count.aa$f*dn31.count.aa$v))/dn31.count.diaa$fv
cps.ttc.gtc<-((dn31.trinuc$ttc*dn31.trinuc$gtc)/(dn31.count.aa$f*dn31.count.aa$v))/dn31.count.diaa$fv
cps.ttc.gtg<-((dn31.trinuc$ttc*dn31.trinuc$gtg)/(dn31.count.aa$f*dn31.count.aa$v))/dn31.count.diaa$fv
cps.ttc.gtt<-((dn31.trinuc$ttc*dn31.trinuc$gtt)/(dn31.count.aa$f*dn31.count.aa$v))/dn31.count.diaa$fv

#Stop codon
#cps.ttc.taa<-((dn31.trinuc$ttc*dn31.trinuc$taa)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttc.tac<-((dn31.trinuc$ttc*dn31.trinuc$tac)/(dn31.count.aa$f*dn31.count.aa$y))/dn31.count.diaa$fy
#Stop codon
#cps.ttc.tag<-((dn31.trinuc$ttc*dn31.trinuc$tag)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttc.tat<-((dn31.trinuc$ttc*dn31.trinuc$tat)/(dn31.count.aa$f*dn31.count.aa$y))/dn31.count.diaa$fy

cps.ttc.tca<-((dn31.trinuc$ttc*dn31.trinuc$tca)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs
cps.ttc.tcc<-((dn31.trinuc$ttc*dn31.trinuc$tcc)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs
cps.ttc.tcg<-((dn31.trinuc$ttc*dn31.trinuc$tcg)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs
cps.ttc.tct<-((dn31.trinuc$ttc*dn31.trinuc$tct)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs

#Stop codon
#cps.ttc.tga<-((dn31.trinuc$ttc*dn31.trinuc$tga)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttc.tgc<-((dn31.trinuc$ttc*dn31.trinuc$tgc)/(dn31.count.aa$f*dn31.count.aa$c))/dn31.count.diaa$fc
cps.ttc.tgg<-((dn31.trinuc$ttc*dn31.trinuc$tgg)/(dn31.count.aa$f*dn31.count.aa$w))/dn31.count.diaa$fw
cps.ttc.tgt<-((dn31.trinuc$ttc*dn31.trinuc$tgt)/(dn31.count.aa$f*dn31.count.aa$c))/dn31.count.diaa$fc

cps.ttc.tta<-((dn31.trinuc$ttc*dn31.trinuc$tta)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttc.ttc<-((dn31.trinuc$ttc*dn31.trinuc$ttc)/(dn31.count.aa$f*dn31.count.aa$f))/dn31.count.diaa$ff
cps.ttc.ttg<-((dn31.trinuc$ttc*dn31.trinuc$ttg)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttc.ttt<-((dn31.trinuc$ttc*dn31.trinuc$ttt)/(dn31.count.aa$f*dn31.count.aa$f))/dn31.count.diaa$ff








cps.ttg.aaa<-((dn31.trinuc$ttg*dn31.trinuc$aaa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ttg.aac<-((dn31.trinuc$ttg*dn31.trinuc$aac)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln
cps.ttg.aag<-((dn31.trinuc$ttg*dn31.trinuc$aag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ttg.aat<-((dn31.trinuc$ttg*dn31.trinuc$aat)/(dn31.count.aa$l*dn31.count.aa$n))/dn31.count.diaa$ln

cps.ttg.aca<-((dn31.trinuc$ttg*dn31.trinuc$aca)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ttg.acc<-((dn31.trinuc$ttg*dn31.trinuc$acc)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ttg.acg<-((dn31.trinuc$ttg*dn31.trinuc$acg)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt
cps.ttg.act<-((dn31.trinuc$ttg*dn31.trinuc$act)/(dn31.count.aa$l*dn31.count.aa$t))/dn31.count.diaa$lt

cps.ttg.aga<-((dn31.trinuc$ttg*dn31.trinuc$aga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ttg.agc<-((dn31.trinuc$ttg*dn31.trinuc$agc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ttg.agg<-((dn31.trinuc$ttg*dn31.trinuc$agg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ttg.agt<-((dn31.trinuc$ttg*dn31.trinuc$agt)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

cps.ttg.ata<-((dn31.trinuc$ttg*dn31.trinuc$ata)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.ttg.atc<-((dn31.trinuc$ttg*dn31.trinuc$atc)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li
cps.ttg.atg<-((dn31.trinuc$ttg*dn31.trinuc$atg)/(dn31.count.aa$l*dn31.count.aa$m))/dn31.count.diaa$lm
cps.ttg.att<-((dn31.trinuc$ttg*dn31.trinuc$att)/(dn31.count.aa$l*dn31.count.aa$i))/dn31.count.diaa$li

cps.ttg.caa<-((dn31.trinuc$ttg*dn31.trinuc$caa)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.ttg.cac<-((dn31.trinuc$ttg*dn31.trinuc$cac)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh
cps.ttg.cag<-((dn31.trinuc$ttg*dn31.trinuc$cag)/(dn31.count.aa$l*dn31.count.aa$q))/dn31.count.diaa$lq
cps.ttg.cat<-((dn31.trinuc$ttg*dn31.trinuc$cat)/(dn31.count.aa$l*dn31.count.aa$h))/dn31.count.diaa$lh

cps.ttg.cca<-((dn31.trinuc$ttg*dn31.trinuc$cca)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ttg.ccc<-((dn31.trinuc$ttg*dn31.trinuc$ccc)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ttg.ccg<-((dn31.trinuc$ttg*dn31.trinuc$ccg)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp
cps.ttg.cct<-((dn31.trinuc$ttg*dn31.trinuc$cct)/(dn31.count.aa$l*dn31.count.aa$p))/dn31.count.diaa$lp

cps.ttg.cga<-((dn31.trinuc$ttg*dn31.trinuc$cga)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ttg.cgc<-((dn31.trinuc$ttg*dn31.trinuc$cgc)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ttg.cgg<-((dn31.trinuc$ttg*dn31.trinuc$cgg)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr
cps.ttg.cgt<-((dn31.trinuc$ttg*dn31.trinuc$cgt)/(dn31.count.aa$l*dn31.count.aa$r))/dn31.count.diaa$lr

cps.ttg.cta<-((dn31.trinuc$ttg*dn31.trinuc$cta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ttg.ctc<-((dn31.trinuc$ttg*dn31.trinuc$ctc)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ttg.ctg<-((dn31.trinuc$ttg*dn31.trinuc$ctg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ttg.ctt<-((dn31.trinuc$ttg*dn31.trinuc$ctt)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll

cps.ttg.gaa<-((dn31.trinuc$ttg*dn31.trinuc$gaa)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.ttg.gac<-((dn31.trinuc$ttg*dn31.trinuc$gac)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld
cps.ttg.gag<-((dn31.trinuc$ttg*dn31.trinuc$gag)/(dn31.count.aa$l*dn31.count.aa$e))/dn31.count.diaa$le
cps.ttg.gat<-((dn31.trinuc$ttg*dn31.trinuc$gat)/(dn31.count.aa$l*dn31.count.aa$d))/dn31.count.diaa$ld

cps.ttg.gca<-((dn31.trinuc$ttg*dn31.trinuc$gca)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ttg.gcc<-((dn31.trinuc$ttg*dn31.trinuc$gcc)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ttg.gcg<-((dn31.trinuc$ttg*dn31.trinuc$gcg)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la
cps.ttg.gct<-((dn31.trinuc$ttg*dn31.trinuc$gct)/(dn31.count.aa$l*dn31.count.aa$a))/dn31.count.diaa$la

cps.ttg.gga<-((dn31.trinuc$ttg*dn31.trinuc$gga)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ttg.ggc<-((dn31.trinuc$ttg*dn31.trinuc$ggc)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ttg.ggg<-((dn31.trinuc$ttg*dn31.trinuc$ggg)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg
cps.ttg.ggt<-((dn31.trinuc$ttg*dn31.trinuc$ggt)/(dn31.count.aa$l*dn31.count.aa$g))/dn31.count.diaa$lg

cps.ttg.gta<-((dn31.trinuc$ttg*dn31.trinuc$gta)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ttg.gtc<-((dn31.trinuc$ttg*dn31.trinuc$gtc)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ttg.gtg<-((dn31.trinuc$ttg*dn31.trinuc$gtg)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv
cps.ttg.gtt<-((dn31.trinuc$ttg*dn31.trinuc$gtt)/(dn31.count.aa$l*dn31.count.aa$v))/dn31.count.diaa$lv

#Stop codon
#cps.ttg.taa<-((dn31.trinuc$ttg*dn31.trinuc$taa)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ttg.tac<-((dn31.trinuc$ttg*dn31.trinuc$tac)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly
#Stop codon
#cps.ttg.tag<-((dn31.trinuc$ttg*dn31.trinuc$tag)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ttg.tat<-((dn31.trinuc$ttg*dn31.trinuc$tat)/(dn31.count.aa$l*dn31.count.aa$y))/dn31.count.diaa$ly

cps.ttg.tca<-((dn31.trinuc$ttg*dn31.trinuc$tca)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ttg.tcc<-((dn31.trinuc$ttg*dn31.trinuc$tcc)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ttg.tcg<-((dn31.trinuc$ttg*dn31.trinuc$tcg)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls
cps.ttg.tct<-((dn31.trinuc$ttg*dn31.trinuc$tct)/(dn31.count.aa$l*dn31.count.aa$s))/dn31.count.diaa$ls

#Stop codon
#cps.ttg.tga<-((dn31.trinuc$ttg*dn31.trinuc$tga)/(dn31.count.aa$l*dn31.count.aa$k))/dn31.count.diaa$lk
cps.ttg.tgc<-((dn31.trinuc$ttg*dn31.trinuc$tgc)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc
cps.ttg.tgg<-((dn31.trinuc$ttg*dn31.trinuc$tgg)/(dn31.count.aa$l*dn31.count.aa$w))/dn31.count.diaa$lw
cps.ttg.tgt<-((dn31.trinuc$ttg*dn31.trinuc$tgt)/(dn31.count.aa$l*dn31.count.aa$c))/dn31.count.diaa$lc

cps.ttg.tta<-((dn31.trinuc$ttg*dn31.trinuc$tta)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ttg.ttc<-((dn31.trinuc$ttg*dn31.trinuc$ttc)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf
cps.ttg.ttg<-((dn31.trinuc$ttg*dn31.trinuc$ttg)/(dn31.count.aa$l*dn31.count.aa$l))/dn31.count.diaa$ll
cps.ttg.ttt<-((dn31.trinuc$ttg*dn31.trinuc$ttt)/(dn31.count.aa$l*dn31.count.aa$f))/dn31.count.diaa$lf








cps.ttt.aaa<-((dn31.trinuc$ttt*dn31.trinuc$aaa)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttt.aac<-((dn31.trinuc$ttt*dn31.trinuc$aac)/(dn31.count.aa$f*dn31.count.aa$n))/dn31.count.diaa$fn
cps.ttt.aag<-((dn31.trinuc$ttt*dn31.trinuc$aag)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttt.aat<-((dn31.trinuc$ttt*dn31.trinuc$aat)/(dn31.count.aa$f*dn31.count.aa$n))/dn31.count.diaa$fn

cps.ttt.aca<-((dn31.trinuc$ttt*dn31.trinuc$aca)/(dn31.count.aa$f*dn31.count.aa$t))/dn31.count.diaa$ft
cps.ttt.acc<-((dn31.trinuc$ttt*dn31.trinuc$acc)/(dn31.count.aa$f*dn31.count.aa$t))/dn31.count.diaa$ft
cps.ttt.acg<-((dn31.trinuc$ttt*dn31.trinuc$acg)/(dn31.count.aa$f*dn31.count.aa$t))/dn31.count.diaa$ft
cps.ttt.act<-((dn31.trinuc$ttt*dn31.trinuc$act)/(dn31.count.aa$f*dn31.count.aa$t))/dn31.count.diaa$ft

cps.ttt.aga<-((dn31.trinuc$ttt*dn31.trinuc$aga)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttt.agc<-((dn31.trinuc$ttt*dn31.trinuc$agc)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs
cps.ttt.agg<-((dn31.trinuc$ttt*dn31.trinuc$agg)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttt.agt<-((dn31.trinuc$ttt*dn31.trinuc$agt)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs

cps.ttt.ata<-((dn31.trinuc$ttt*dn31.trinuc$ata)/(dn31.count.aa$f*dn31.count.aa$i))/dn31.count.diaa$fi
cps.ttt.atc<-((dn31.trinuc$ttt*dn31.trinuc$atc)/(dn31.count.aa$f*dn31.count.aa$i))/dn31.count.diaa$fi
cps.ttt.atg<-((dn31.trinuc$ttt*dn31.trinuc$atg)/(dn31.count.aa$f*dn31.count.aa$m))/dn31.count.diaa$fm
cps.ttt.att<-((dn31.trinuc$ttt*dn31.trinuc$att)/(dn31.count.aa$f*dn31.count.aa$i))/dn31.count.diaa$fi

cps.ttt.caa<-((dn31.trinuc$ttt*dn31.trinuc$caa)/(dn31.count.aa$f*dn31.count.aa$q))/dn31.count.diaa$fq
cps.ttt.cac<-((dn31.trinuc$ttt*dn31.trinuc$cac)/(dn31.count.aa$f*dn31.count.aa$h))/dn31.count.diaa$fh
cps.ttt.cag<-((dn31.trinuc$ttt*dn31.trinuc$cag)/(dn31.count.aa$f*dn31.count.aa$q))/dn31.count.diaa$fq
cps.ttt.cat<-((dn31.trinuc$ttt*dn31.trinuc$cat)/(dn31.count.aa$f*dn31.count.aa$h))/dn31.count.diaa$fh

cps.ttt.cca<-((dn31.trinuc$ttt*dn31.trinuc$cca)/(dn31.count.aa$f*dn31.count.aa$p))/dn31.count.diaa$fp
cps.ttt.ccc<-((dn31.trinuc$ttt*dn31.trinuc$ccc)/(dn31.count.aa$f*dn31.count.aa$p))/dn31.count.diaa$fp
cps.ttt.ccg<-((dn31.trinuc$ttt*dn31.trinuc$ccg)/(dn31.count.aa$f*dn31.count.aa$p))/dn31.count.diaa$fp
cps.ttt.cct<-((dn31.trinuc$ttt*dn31.trinuc$cct)/(dn31.count.aa$f*dn31.count.aa$p))/dn31.count.diaa$fp

cps.ttt.cga<-((dn31.trinuc$ttt*dn31.trinuc$cga)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttt.cgc<-((dn31.trinuc$ttt*dn31.trinuc$cgc)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttt.cgg<-((dn31.trinuc$ttt*dn31.trinuc$cgg)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr
cps.ttt.cgt<-((dn31.trinuc$ttt*dn31.trinuc$cgt)/(dn31.count.aa$f*dn31.count.aa$r))/dn31.count.diaa$fr

cps.ttt.cta<-((dn31.trinuc$ttt*dn31.trinuc$cta)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttt.ctc<-((dn31.trinuc$ttt*dn31.trinuc$ctc)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttt.ctg<-((dn31.trinuc$ttt*dn31.trinuc$ctg)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttt.ctt<-((dn31.trinuc$ttt*dn31.trinuc$ctt)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl

cps.ttt.gaa<-((dn31.trinuc$ttt*dn31.trinuc$gaa)/(dn31.count.aa$f*dn31.count.aa$e))/dn31.count.diaa$fe
cps.ttt.gac<-((dn31.trinuc$ttt*dn31.trinuc$gac)/(dn31.count.aa$f*dn31.count.aa$d))/dn31.count.diaa$fd
cps.ttt.gag<-((dn31.trinuc$ttt*dn31.trinuc$gag)/(dn31.count.aa$f*dn31.count.aa$e))/dn31.count.diaa$fe
cps.ttt.gat<-((dn31.trinuc$ttt*dn31.trinuc$gat)/(dn31.count.aa$f*dn31.count.aa$d))/dn31.count.diaa$fd

cps.ttt.gca<-((dn31.trinuc$ttt*dn31.trinuc$gca)/(dn31.count.aa$f*dn31.count.aa$a))/dn31.count.diaa$fa
cps.ttt.gcc<-((dn31.trinuc$ttt*dn31.trinuc$gcc)/(dn31.count.aa$f*dn31.count.aa$a))/dn31.count.diaa$fa
cps.ttt.gcg<-((dn31.trinuc$ttt*dn31.trinuc$gcg)/(dn31.count.aa$f*dn31.count.aa$a))/dn31.count.diaa$fa
cps.ttt.gct<-((dn31.trinuc$ttt*dn31.trinuc$gct)/(dn31.count.aa$f*dn31.count.aa$a))/dn31.count.diaa$fa

cps.ttt.gga<-((dn31.trinuc$ttt*dn31.trinuc$gga)/(dn31.count.aa$f*dn31.count.aa$g))/dn31.count.diaa$fg
cps.ttt.ggc<-((dn31.trinuc$ttt*dn31.trinuc$ggc)/(dn31.count.aa$f*dn31.count.aa$g))/dn31.count.diaa$fg
cps.ttt.ggg<-((dn31.trinuc$ttt*dn31.trinuc$ggg)/(dn31.count.aa$f*dn31.count.aa$g))/dn31.count.diaa$fg
cps.ttt.ggt<-((dn31.trinuc$ttt*dn31.trinuc$ggt)/(dn31.count.aa$f*dn31.count.aa$g))/dn31.count.diaa$fg

cps.ttt.gta<-((dn31.trinuc$ttt*dn31.trinuc$gta)/(dn31.count.aa$f*dn31.count.aa$v))/dn31.count.diaa$fv
cps.ttt.gtc<-((dn31.trinuc$ttt*dn31.trinuc$gtc)/(dn31.count.aa$f*dn31.count.aa$v))/dn31.count.diaa$fv
cps.ttt.gtg<-((dn31.trinuc$ttt*dn31.trinuc$gtg)/(dn31.count.aa$f*dn31.count.aa$v))/dn31.count.diaa$fv
cps.ttt.gtt<-((dn31.trinuc$ttt*dn31.trinuc$gtt)/(dn31.count.aa$f*dn31.count.aa$v))/dn31.count.diaa$fv

#Stop codon
#cps.ttt.taa<-((dn31.trinuc$ttt*dn31.trinuc$taa)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttt.tac<-((dn31.trinuc$ttt*dn31.trinuc$tac)/(dn31.count.aa$f*dn31.count.aa$y))/dn31.count.diaa$fy
#Stop codon
#cps.ttt.tag<-((dn31.trinuc$ttt*dn31.trinuc$tag)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttt.tat<-((dn31.trinuc$ttt*dn31.trinuc$tat)/(dn31.count.aa$f*dn31.count.aa$y))/dn31.count.diaa$fy

cps.ttt.tca<-((dn31.trinuc$ttt*dn31.trinuc$tca)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs
cps.ttt.tcc<-((dn31.trinuc$ttt*dn31.trinuc$tcc)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs
cps.ttt.tcg<-((dn31.trinuc$ttt*dn31.trinuc$tcg)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs
cps.ttt.tct<-((dn31.trinuc$ttt*dn31.trinuc$tct)/(dn31.count.aa$f*dn31.count.aa$s))/dn31.count.diaa$fs

#Stop codon
#cps.ttt.tga<-((dn31.trinuc$ttt*dn31.trinuc$tga)/(dn31.count.aa$f*dn31.count.aa$k))/dn31.count.diaa$fk
cps.ttt.tgc<-((dn31.trinuc$ttt*dn31.trinuc$tgc)/(dn31.count.aa$f*dn31.count.aa$c))/dn31.count.diaa$fc
cps.ttt.tgg<-((dn31.trinuc$ttt*dn31.trinuc$tgg)/(dn31.count.aa$f*dn31.count.aa$w))/dn31.count.diaa$fw
cps.ttt.tgt<-((dn31.trinuc$ttt*dn31.trinuc$tgt)/(dn31.count.aa$f*dn31.count.aa$c))/dn31.count.diaa$fc

cps.ttt.tta<-((dn31.trinuc$ttt*dn31.trinuc$tta)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttt.ttc<-((dn31.trinuc$ttt*dn31.trinuc$ttc)/(dn31.count.aa$f*dn31.count.aa$f))/dn31.count.diaa$ff
cps.ttt.ttg<-((dn31.trinuc$ttt*dn31.trinuc$ttg)/(dn31.count.aa$f*dn31.count.aa$l))/dn31.count.diaa$fl
cps.ttt.ttt<-((dn31.trinuc$ttt*dn31.trinuc$ttt)/(dn31.count.aa$f*dn31.count.aa$f))/dn31.count.diaa$ff


























































#CPS and CPS observed
#cps.observed/cps.expected
cps.obs.exp.aaaaaa<-log(dn31.codon.pair$aaaaaa/cps.aaa.aaa)
cps.obs.exp.aaaaac<-log(dn31.codon.pair$aaaaac/cps.aaa.aac)
cps.obs.exp.aaaaag<-log(dn31.codon.pair$aaaaag/cps.aaa.aag)
cps.obs.exp.aaaaat<-log(dn31.codon.pair$aaaaat/cps.aaa.aat)
cps.obs.exp.aaaaca<-log(dn31.codon.pair$aaaaca/cps.aaa.aca)
cps.obs.exp.aaaacc<-log(dn31.codon.pair$aaaacc/cps.aaa.acc)
cps.obs.exp.aaaacg<-log(dn31.codon.pair$aaaacg/cps.aaa.acg)
cps.obs.exp.aaaact<-log(dn31.codon.pair$aaaact/cps.aaa.act)
cps.obs.exp.aaaaga<-log(dn31.codon.pair$aaaaga/cps.aaa.aga)
cps.obs.exp.aaaagc<-log(dn31.codon.pair$aaaagc/cps.aaa.agc)
cps.obs.exp.aaaagg<-log(dn31.codon.pair$aaaagg/cps.aaa.agg)
cps.obs.exp.aaaagt<-log(dn31.codon.pair$aaaagt/cps.aaa.agt)
cps.obs.exp.aaaata<-log(dn31.codon.pair$aaaata/cps.aaa.ata)
cps.obs.exp.aaaatc<-log(dn31.codon.pair$aaaatc/cps.aaa.atc)
cps.obs.exp.aaaatg<-log(dn31.codon.pair$aaaatg/cps.aaa.atg)
cps.obs.exp.aaaatt<-log(dn31.codon.pair$aaaatt/cps.aaa.att)

cps.obs.exp.aaacaa<-log(dn31.codon.pair$aaacaa/cps.aaa.caa)
cps.obs.exp.aaacac<-log(dn31.codon.pair$aaacac/cps.aaa.cac)
cps.obs.exp.aaacag<-log(dn31.codon.pair$aaacag/cps.aaa.cag)
cps.obs.exp.aaacat<-log(dn31.codon.pair$aaacat/cps.aaa.cat)
cps.obs.exp.aaacca<-log(dn31.codon.pair$aaacca/cps.aaa.cca)
cps.obs.exp.aaaccc<-log(dn31.codon.pair$aaaccc/cps.aaa.ccc)
cps.obs.exp.aaaccg<-log(dn31.codon.pair$aaaccg/cps.aaa.ccg)
cps.obs.exp.aaacct<-log(dn31.codon.pair$aaacct/cps.aaa.cct)
cps.obs.exp.aaacga<-log(dn31.codon.pair$aaacga/cps.aaa.cga)
cps.obs.exp.aaacgc<-log(dn31.codon.pair$aaacgc/cps.aaa.cgc)
cps.obs.exp.aaacgg<-log(dn31.codon.pair$aaacgg/cps.aaa.cgg)
cps.obs.exp.aaacgt<-log(dn31.codon.pair$aaacgt/cps.aaa.cgt)
cps.obs.exp.aaacta<-log(dn31.codon.pair$aaacta/cps.aaa.cta)
cps.obs.exp.aaactc<-log(dn31.codon.pair$aaactc/cps.aaa.ctc)
cps.obs.exp.aaactg<-log(dn31.codon.pair$aaactg/cps.aaa.ctg)
cps.obs.exp.aaactt<-log(dn31.codon.pair$aaactt/cps.aaa.ctt)

cps.obs.exp.aaagaa<-log(dn31.codon.pair$aaagaa/cps.aaa.gaa)
cps.obs.exp.aaagac<-log(dn31.codon.pair$aaagac/cps.aaa.gac)
cps.obs.exp.aaagag<-log(dn31.codon.pair$aaagag/cps.aaa.gag)
cps.obs.exp.aaagat<-log(dn31.codon.pair$aaagat/cps.aaa.gat)
cps.obs.exp.aaagca<-log(dn31.codon.pair$aaagca/cps.aaa.gca)
cps.obs.exp.aaagcc<-log(dn31.codon.pair$aaagcc/cps.aaa.gcc)
cps.obs.exp.aaagcg<-log(dn31.codon.pair$aaagcg/cps.aaa.gcg)
cps.obs.exp.aaagct<-log(dn31.codon.pair$aaagct/cps.aaa.gct)
cps.obs.exp.aaagga<-log(dn31.codon.pair$aaagga/cps.aaa.gga)
cps.obs.exp.aaaggc<-log(dn31.codon.pair$aaaggc/cps.aaa.ggc)
cps.obs.exp.aaaggg<-log(dn31.codon.pair$aaaggg/cps.aaa.ggg)
cps.obs.exp.aaaggt<-log(dn31.codon.pair$aaaggt/cps.aaa.ggt)
cps.obs.exp.aaagta<-log(dn31.codon.pair$aaagta/cps.aaa.gta)
cps.obs.exp.aaagtc<-log(dn31.codon.pair$aaagtc/cps.aaa.gtc)
cps.obs.exp.aaagtg<-log(dn31.codon.pair$aaagtg/cps.aaa.gtg)
cps.obs.exp.aaagtt<-log(dn31.codon.pair$aaagtt/cps.aaa.gtt)

#cps.obs.exp.aaataa<-log(dn31.codon.pair$aaataa/cps.aaa.taa)
cps.obs.exp.aaatac<-log(dn31.codon.pair$aaatac/cps.aaa.tac)
#cps.obs.exp.aaatag<-log(dn31.codon.pair$aaatag/cps.aaa.tag)
cps.obs.exp.aaatat<-log(dn31.codon.pair$aaatat/cps.aaa.tat)
cps.obs.exp.aaatca<-log(dn31.codon.pair$aaatca/cps.aaa.tca)
cps.obs.exp.aaatcc<-log(dn31.codon.pair$aaatcc/cps.aaa.tcc)
cps.obs.exp.aaatcg<-log(dn31.codon.pair$aaatcg/cps.aaa.tcg)
cps.obs.exp.aaatct<-log(dn31.codon.pair$aaatct/cps.aaa.tct)
#cps.obs.exp.aaatga<-log(dn31.codon.pair$aaatga/cps.aaa.tga)
cps.obs.exp.aaatgc<-log(dn31.codon.pair$aaatgc/cps.aaa.tgc)
cps.obs.exp.aaatgg<-log(dn31.codon.pair$aaatgg/cps.aaa.tgg)
cps.obs.exp.aaatgt<-log(dn31.codon.pair$aaatgt/cps.aaa.tgt)
cps.obs.exp.aaatta<-log(dn31.codon.pair$aaatta/cps.aaa.tta)
cps.obs.exp.aaattc<-log(dn31.codon.pair$aaattc/cps.aaa.ttc)
cps.obs.exp.aaattg<-log(dn31.codon.pair$aaattg/cps.aaa.ttg)
cps.obs.exp.aaattt<-log(dn31.codon.pair$aaattt/cps.aaa.ttt)




cps.obs.exp.aacaaa<-log(dn31.codon.pair$aacaaa/cps.aac.aaa)
cps.obs.exp.aacaac<-log(dn31.codon.pair$aacaac/cps.aac.aac)
cps.obs.exp.aacaag<-log(dn31.codon.pair$aacaag/cps.aac.aag)
cps.obs.exp.aacaat<-log(dn31.codon.pair$aacaat/cps.aac.aat)
cps.obs.exp.aacaca<-log(dn31.codon.pair$aacaca/cps.aac.aca)
cps.obs.exp.aacacc<-log(dn31.codon.pair$aacacc/cps.aac.acc)
cps.obs.exp.aacacg<-log(dn31.codon.pair$aacacg/cps.aac.acg)
cps.obs.exp.aacact<-log(dn31.codon.pair$aacact/cps.aac.act)
cps.obs.exp.aacaga<-log(dn31.codon.pair$aacaga/cps.aac.aga)
cps.obs.exp.aacagc<-log(dn31.codon.pair$aacagc/cps.aac.agc)
cps.obs.exp.aacagg<-log(dn31.codon.pair$aacagg/cps.aac.agg)
cps.obs.exp.aacagt<-log(dn31.codon.pair$aacagt/cps.aac.agt)
cps.obs.exp.aacata<-log(dn31.codon.pair$aacata/cps.aac.ata)
cps.obs.exp.aacatc<-log(dn31.codon.pair$aacatc/cps.aac.atc)
cps.obs.exp.aacatg<-log(dn31.codon.pair$aacatg/cps.aac.atg)
cps.obs.exp.aacatt<-log(dn31.codon.pair$aacatt/cps.aac.att)

cps.obs.exp.aaccaa<-log(dn31.codon.pair$aaccaa/cps.aac.caa)
cps.obs.exp.aaccac<-log(dn31.codon.pair$aaccac/cps.aac.cac)
cps.obs.exp.aaccag<-log(dn31.codon.pair$aaccag/cps.aac.cag)
cps.obs.exp.aaccat<-log(dn31.codon.pair$aaccat/cps.aac.cat)
cps.obs.exp.aaccca<-log(dn31.codon.pair$aaccca/cps.aac.cca)
cps.obs.exp.aacccc<-log(dn31.codon.pair$aacccc/cps.aac.ccc)
cps.obs.exp.aacccg<-log(dn31.codon.pair$aacccg/cps.aac.ccg)
cps.obs.exp.aaccct<-log(dn31.codon.pair$aaccct/cps.aac.cct)
cps.obs.exp.aaccga<-log(dn31.codon.pair$aaccga/cps.aac.cga)
cps.obs.exp.aaccgc<-log(dn31.codon.pair$aaccgc/cps.aac.cgc)
cps.obs.exp.aaccgg<-log(dn31.codon.pair$aaccgg/cps.aac.cgg)
cps.obs.exp.aaccgt<-log(dn31.codon.pair$aaccgt/cps.aac.cgt)
cps.obs.exp.aaccta<-log(dn31.codon.pair$aaccta/cps.aac.cta)
cps.obs.exp.aacctc<-log(dn31.codon.pair$aacctc/cps.aac.ctc)
cps.obs.exp.aacctg<-log(dn31.codon.pair$aacctg/cps.aac.ctg)
cps.obs.exp.aacctt<-log(dn31.codon.pair$aacctt/cps.aac.ctt)

cps.obs.exp.aacgaa<-log(dn31.codon.pair$aacgaa/cps.aac.gaa)
cps.obs.exp.aacgac<-log(dn31.codon.pair$aacgac/cps.aac.gac)
cps.obs.exp.aacgag<-log(dn31.codon.pair$aacgag/cps.aac.gag)
cps.obs.exp.aacgat<-log(dn31.codon.pair$aacgat/cps.aac.gat)
cps.obs.exp.aacgca<-log(dn31.codon.pair$aacgca/cps.aac.gca)
cps.obs.exp.aacgcc<-log(dn31.codon.pair$aacgcc/cps.aac.gcc)
cps.obs.exp.aacgcg<-log(dn31.codon.pair$aacgcg/cps.aac.gcg)
cps.obs.exp.aacgct<-log(dn31.codon.pair$aacgct/cps.aac.gct)
cps.obs.exp.aacgga<-log(dn31.codon.pair$aacgga/cps.aac.gga)
cps.obs.exp.aacggc<-log(dn31.codon.pair$aacggc/cps.aac.ggc)
cps.obs.exp.aacggg<-log(dn31.codon.pair$aacggg/cps.aac.ggg)
cps.obs.exp.aacggt<-log(dn31.codon.pair$aacggt/cps.aac.ggt)
cps.obs.exp.aacgta<-log(dn31.codon.pair$aacgta/cps.aac.gta)
cps.obs.exp.aacgtc<-log(dn31.codon.pair$aacgtc/cps.aac.gtc)
cps.obs.exp.aacgtg<-log(dn31.codon.pair$aacgtg/cps.aac.gtg)
cps.obs.exp.aacgtt<-log(dn31.codon.pair$aacgtt/cps.aac.gtt)

#cps.obs.exp.aactaa<-log(dn31.codon.pair$aactaa/cps.aac.taa)
cps.obs.exp.aactac<-log(dn31.codon.pair$aactac/cps.aac.tac)
#cps.obs.exp.aactag<-log(dn31.codon.pair$aactag/cps.aac.tag)
cps.obs.exp.aactat<-log(dn31.codon.pair$aactat/cps.aac.tat)
cps.obs.exp.aactca<-log(dn31.codon.pair$aactca/cps.aac.tca)
cps.obs.exp.aactcc<-log(dn31.codon.pair$aactcc/cps.aac.tcc)
cps.obs.exp.aactcg<-log(dn31.codon.pair$aactcg/cps.aac.tcg)
cps.obs.exp.aactct<-log(dn31.codon.pair$aactct/cps.aac.tct)
#cps.obs.exp.aactga<-log(dn31.codon.pair$aactga/cps.aac.tga)
cps.obs.exp.aactgc<-log(dn31.codon.pair$aactgc/cps.aac.tgc)
cps.obs.exp.aactgg<-log(dn31.codon.pair$aactgg/cps.aac.tgg)
cps.obs.exp.aactgt<-log(dn31.codon.pair$aactgt/cps.aac.tgt)
cps.obs.exp.aactta<-log(dn31.codon.pair$aactta/cps.aac.tta)
cps.obs.exp.aacttc<-log(dn31.codon.pair$aacttc/cps.aac.ttc)
cps.obs.exp.aacttg<-log(dn31.codon.pair$aacttg/cps.aac.ttg)
cps.obs.exp.aacttt<-log(dn31.codon.pair$aacttt/cps.aac.ttt)








cps.obs.exp.aagaaa<-log(dn31.codon.pair$aagaaa/cps.aag.aaa)
cps.obs.exp.aagaac<-log(dn31.codon.pair$aagaac/cps.aag.aac)
cps.obs.exp.aagaag<-log(dn31.codon.pair$aagaag/cps.aag.aag)
cps.obs.exp.aagaat<-log(dn31.codon.pair$aagaat/cps.aag.aat)
cps.obs.exp.aagaca<-log(dn31.codon.pair$aagaca/cps.aag.aca)
cps.obs.exp.aagacc<-log(dn31.codon.pair$aagacc/cps.aag.acc)
cps.obs.exp.aagacg<-log(dn31.codon.pair$aagacg/cps.aag.acg)
cps.obs.exp.aagact<-log(dn31.codon.pair$aagact/cps.aag.act)
cps.obs.exp.aagaga<-log(dn31.codon.pair$aagaga/cps.aag.aga)
cps.obs.exp.aagagc<-log(dn31.codon.pair$aagagc/cps.aag.agc)
cps.obs.exp.aagagg<-log(dn31.codon.pair$aagagg/cps.aag.agg)
cps.obs.exp.aagagt<-log(dn31.codon.pair$aagagt/cps.aag.agt)
cps.obs.exp.aagata<-log(dn31.codon.pair$aagata/cps.aag.ata)
cps.obs.exp.aagatc<-log(dn31.codon.pair$aagatc/cps.aag.atc)
cps.obs.exp.aagatg<-log(dn31.codon.pair$aagatg/cps.aag.atg)
cps.obs.exp.aagatt<-log(dn31.codon.pair$aagatt/cps.aag.att)

cps.obs.exp.aagcaa<-log(dn31.codon.pair$aagcaa/cps.aag.caa)
cps.obs.exp.aagcac<-log(dn31.codon.pair$aagcac/cps.aag.cac)
cps.obs.exp.aagcag<-log(dn31.codon.pair$aagcag/cps.aag.cag)
cps.obs.exp.aagcat<-log(dn31.codon.pair$aagcat/cps.aag.cat)
cps.obs.exp.aagcca<-log(dn31.codon.pair$aagcca/cps.aag.cca)
cps.obs.exp.aagccc<-log(dn31.codon.pair$aagccc/cps.aag.ccc)
cps.obs.exp.aagccg<-log(dn31.codon.pair$aagccg/cps.aag.ccg)
cps.obs.exp.aagcct<-log(dn31.codon.pair$aagcct/cps.aag.cct)
cps.obs.exp.aagcga<-log(dn31.codon.pair$aagcga/cps.aag.cga)
cps.obs.exp.aagcgc<-log(dn31.codon.pair$aagcgc/cps.aag.cgc)
cps.obs.exp.aagcgg<-log(dn31.codon.pair$aagcgg/cps.aag.cgg)
cps.obs.exp.aagcgt<-log(dn31.codon.pair$aagcgt/cps.aag.cgt)
cps.obs.exp.aagcta<-log(dn31.codon.pair$aagcta/cps.aag.cta)
cps.obs.exp.aagctc<-log(dn31.codon.pair$aagctc/cps.aag.ctc)
cps.obs.exp.aagctg<-log(dn31.codon.pair$aagctg/cps.aag.ctg)
cps.obs.exp.aagctt<-log(dn31.codon.pair$aagctt/cps.aag.ctt)

cps.obs.exp.aaggaa<-log(dn31.codon.pair$aaggaa/cps.aag.gaa)
cps.obs.exp.aaggac<-log(dn31.codon.pair$aaggac/cps.aag.gac)
cps.obs.exp.aaggag<-log(dn31.codon.pair$aaggag/cps.aag.gag)
cps.obs.exp.aaggat<-log(dn31.codon.pair$aaggat/cps.aag.gat)
cps.obs.exp.aaggca<-log(dn31.codon.pair$aaggca/cps.aag.gca)
cps.obs.exp.aaggcc<-log(dn31.codon.pair$aaggcc/cps.aag.gcc)
cps.obs.exp.aaggcg<-log(dn31.codon.pair$aaggcg/cps.aag.gcg)
cps.obs.exp.aaggct<-log(dn31.codon.pair$aaggct/cps.aag.gct)
cps.obs.exp.aaggga<-log(dn31.codon.pair$aaggga/cps.aag.gga)
cps.obs.exp.aagggc<-log(dn31.codon.pair$aagggc/cps.aag.ggc)
cps.obs.exp.aagggg<-log(dn31.codon.pair$aagggg/cps.aag.ggg)
cps.obs.exp.aagggt<-log(dn31.codon.pair$aagggt/cps.aag.ggt)
cps.obs.exp.aaggta<-log(dn31.codon.pair$aaggta/cps.aag.gta)
cps.obs.exp.aaggtc<-log(dn31.codon.pair$aaggtc/cps.aag.gtc)
cps.obs.exp.aaggtg<-log(dn31.codon.pair$aaggtg/cps.aag.gtg)
cps.obs.exp.aaggtt<-log(dn31.codon.pair$aaggtt/cps.aag.gtt)

#cps.obs.exp.aagtaa<-log(dn31.codon.pair$aagtaa/cps.aag.taa)
cps.obs.exp.aagtac<-log(dn31.codon.pair$aagtac/cps.aag.tac)
#cps.obs.exp.aagtag<-log(dn31.codon.pair$aagtag/cps.aag.tag)
cps.obs.exp.aagtat<-log(dn31.codon.pair$aagtat/cps.aag.tat)
cps.obs.exp.aagtca<-log(dn31.codon.pair$aagtca/cps.aag.tca)
cps.obs.exp.aagtcc<-log(dn31.codon.pair$aagtcc/cps.aag.tcc)
cps.obs.exp.aagtcg<-log(dn31.codon.pair$aagtcg/cps.aag.tcg)
cps.obs.exp.aagtct<-log(dn31.codon.pair$aagtct/cps.aag.tct)
#cps.obs.exp.aagtga<-log(dn31.codon.pair$aagtga/cps.aag.tga)
cps.obs.exp.aagtgc<-log(dn31.codon.pair$aagtgc/cps.aag.tgc)
cps.obs.exp.aagtgg<-log(dn31.codon.pair$aagtgg/cps.aag.tgg)
cps.obs.exp.aagtgt<-log(dn31.codon.pair$aagtgt/cps.aag.tgt)
cps.obs.exp.aagtta<-log(dn31.codon.pair$aagtta/cps.aag.tta)
cps.obs.exp.aagttc<-log(dn31.codon.pair$aagttc/cps.aag.ttc)
cps.obs.exp.aagttg<-log(dn31.codon.pair$aagttg/cps.aag.ttg)
cps.obs.exp.aagttt<-log(dn31.codon.pair$aagttt/cps.aag.ttt)







cps.obs.exp.aataaa<-log(dn31.codon.pair$aataaa/cps.aat.aaa)
cps.obs.exp.aataac<-log(dn31.codon.pair$aataac/cps.aat.aac)
cps.obs.exp.aataag<-log(dn31.codon.pair$aataag/cps.aat.aag)
cps.obs.exp.aataat<-log(dn31.codon.pair$aataat/cps.aat.aat)
cps.obs.exp.aataca<-log(dn31.codon.pair$aataca/cps.aat.aca)
cps.obs.exp.aatacc<-log(dn31.codon.pair$aatacc/cps.aat.acc)
cps.obs.exp.aatacg<-log(dn31.codon.pair$aatacg/cps.aat.acg)
cps.obs.exp.aatact<-log(dn31.codon.pair$aatact/cps.aat.act)
cps.obs.exp.aataga<-log(dn31.codon.pair$aataga/cps.aat.aga)
cps.obs.exp.aatagc<-log(dn31.codon.pair$aatagc/cps.aat.agc)
cps.obs.exp.aatagg<-log(dn31.codon.pair$aatagg/cps.aat.agg)
cps.obs.exp.aatagt<-log(dn31.codon.pair$aatagt/cps.aat.agt)
cps.obs.exp.aatata<-log(dn31.codon.pair$aatata/cps.aat.ata)
cps.obs.exp.aatatc<-log(dn31.codon.pair$aatatc/cps.aat.atc)
cps.obs.exp.aatatg<-log(dn31.codon.pair$aatatg/cps.aat.atg)
cps.obs.exp.aatatt<-log(dn31.codon.pair$aatatt/cps.aat.att)

cps.obs.exp.aatcaa<-log(dn31.codon.pair$aatcaa/cps.aat.caa)
cps.obs.exp.aatcac<-log(dn31.codon.pair$aatcac/cps.aat.cac)
cps.obs.exp.aatcag<-log(dn31.codon.pair$aatcag/cps.aat.cag)
cps.obs.exp.aatcat<-log(dn31.codon.pair$aatcat/cps.aat.cat)
cps.obs.exp.aatcca<-log(dn31.codon.pair$aatcca/cps.aat.cca)
cps.obs.exp.aatccc<-log(dn31.codon.pair$aatccc/cps.aat.ccc)
cps.obs.exp.aatccg<-log(dn31.codon.pair$aatccg/cps.aat.ccg)
cps.obs.exp.aatcct<-log(dn31.codon.pair$aatcct/cps.aat.cct)
cps.obs.exp.aatcga<-log(dn31.codon.pair$aatcga/cps.aat.cga)
cps.obs.exp.aatcgc<-log(dn31.codon.pair$aatcgc/cps.aat.cgc)
cps.obs.exp.aatcgg<-log(dn31.codon.pair$aatcgg/cps.aat.cgg)
cps.obs.exp.aatcgt<-log(dn31.codon.pair$aatcgt/cps.aat.cgt)
cps.obs.exp.aatcta<-log(dn31.codon.pair$aatcta/cps.aat.cta)
cps.obs.exp.aatctc<-log(dn31.codon.pair$aatctc/cps.aat.ctc)
cps.obs.exp.aatctg<-log(dn31.codon.pair$aatctg/cps.aat.ctg)
cps.obs.exp.aatctt<-log(dn31.codon.pair$aatctt/cps.aat.ctt)

cps.obs.exp.aatgaa<-log(dn31.codon.pair$aatgaa/cps.aat.gaa)
cps.obs.exp.aatgac<-log(dn31.codon.pair$aatgac/cps.aat.gac)
cps.obs.exp.aatgag<-log(dn31.codon.pair$aatgag/cps.aat.gag)
cps.obs.exp.aatgat<-log(dn31.codon.pair$aatgat/cps.aat.gat)
cps.obs.exp.aatgca<-log(dn31.codon.pair$aatgca/cps.aat.gca)
cps.obs.exp.aatgcc<-log(dn31.codon.pair$aatgcc/cps.aat.gcc)
cps.obs.exp.aatgcg<-log(dn31.codon.pair$aatgcg/cps.aat.gcg)
cps.obs.exp.aatgct<-log(dn31.codon.pair$aatgct/cps.aat.gct)
cps.obs.exp.aatgga<-log(dn31.codon.pair$aatgga/cps.aat.gga)
cps.obs.exp.aatggc<-log(dn31.codon.pair$aatggc/cps.aat.ggc)
cps.obs.exp.aatggg<-log(dn31.codon.pair$aatggg/cps.aat.ggg)
cps.obs.exp.aatggt<-log(dn31.codon.pair$aatggt/cps.aat.ggt)
cps.obs.exp.aatgta<-log(dn31.codon.pair$aatgta/cps.aat.gta)
cps.obs.exp.aatgtc<-log(dn31.codon.pair$aatgtc/cps.aat.gtc)
cps.obs.exp.aatgtg<-log(dn31.codon.pair$aatgtg/cps.aat.gtg)
cps.obs.exp.aatgtt<-log(dn31.codon.pair$aatgtt/cps.aat.gtt)

#cps.obs.exp.aattaa<-log(dn31.codon.pair$aattaa/cps.aat.taa)
cps.obs.exp.aattac<-log(dn31.codon.pair$aattac/cps.aat.tac)
#cps.obs.exp.aattag<-log(dn31.codon.pair$aattag/cps.aat.tag)
cps.obs.exp.aattat<-log(dn31.codon.pair$aattat/cps.aat.tat)
cps.obs.exp.aattca<-log(dn31.codon.pair$aattca/cps.aat.tca)
cps.obs.exp.aattcc<-log(dn31.codon.pair$aattcc/cps.aat.tcc)
cps.obs.exp.aattcg<-log(dn31.codon.pair$aattcg/cps.aat.tcg)
cps.obs.exp.aattct<-log(dn31.codon.pair$aattct/cps.aat.tct)
#cps.obs.exp.aattga<-log(dn31.codon.pair$aattga/cps.aat.tga)
cps.obs.exp.aattgc<-log(dn31.codon.pair$aattgc/cps.aat.tgc)
cps.obs.exp.aattgg<-log(dn31.codon.pair$aattgg/cps.aat.tgg)
cps.obs.exp.aattgt<-log(dn31.codon.pair$aattgt/cps.aat.tgt)
cps.obs.exp.aattta<-log(dn31.codon.pair$aattta/cps.aat.tta)
cps.obs.exp.aatttc<-log(dn31.codon.pair$aatttc/cps.aat.ttc)
cps.obs.exp.aatttg<-log(dn31.codon.pair$aatttg/cps.aat.ttg)
cps.obs.exp.aatttt<-log(dn31.codon.pair$aatttt/cps.aat.ttt)

















cps.obs.exp.acaaaa<-log(dn31.codon.pair$acaaaa/cps.aca.aaa)
cps.obs.exp.acaaac<-log(dn31.codon.pair$acaaac/cps.aca.aac)
cps.obs.exp.acaaag<-log(dn31.codon.pair$acaaag/cps.aca.aag)
cps.obs.exp.acaaat<-log(dn31.codon.pair$acaaat/cps.aca.aat)
cps.obs.exp.acaaca<-log(dn31.codon.pair$acaaca/cps.aca.aca)
cps.obs.exp.acaacc<-log(dn31.codon.pair$acaacc/cps.aca.acc)
cps.obs.exp.acaacg<-log(dn31.codon.pair$acaacg/cps.aca.acg)
cps.obs.exp.acaact<-log(dn31.codon.pair$acaact/cps.aca.act)
cps.obs.exp.acaaga<-log(dn31.codon.pair$acaaga/cps.aca.aga)
cps.obs.exp.acaagc<-log(dn31.codon.pair$acaagc/cps.aca.agc)
cps.obs.exp.acaagg<-log(dn31.codon.pair$acaagg/cps.aca.agg)
cps.obs.exp.acaagt<-log(dn31.codon.pair$acaagt/cps.aca.agt)
cps.obs.exp.acaata<-log(dn31.codon.pair$acaata/cps.aca.ata)
cps.obs.exp.acaatc<-log(dn31.codon.pair$acaatc/cps.aca.atc)
cps.obs.exp.acaatg<-log(dn31.codon.pair$acaatg/cps.aca.atg)
cps.obs.exp.acaatt<-log(dn31.codon.pair$acaatt/cps.aca.att)

cps.obs.exp.acacaa<-log(dn31.codon.pair$acacaa/cps.aca.caa)
cps.obs.exp.acacac<-log(dn31.codon.pair$acacac/cps.aca.cac)
cps.obs.exp.acacag<-log(dn31.codon.pair$acacag/cps.aca.cag)
cps.obs.exp.acacat<-log(dn31.codon.pair$acacat/cps.aca.cat)
cps.obs.exp.acacca<-log(dn31.codon.pair$acacca/cps.aca.cca)
cps.obs.exp.acaccc<-log(dn31.codon.pair$acaccc/cps.aca.ccc)
cps.obs.exp.acaccg<-log(dn31.codon.pair$acaccg/cps.aca.ccg)
cps.obs.exp.acacct<-log(dn31.codon.pair$acacct/cps.aca.cct)
cps.obs.exp.acacga<-log(dn31.codon.pair$acacga/cps.aca.cga)
cps.obs.exp.acacgc<-log(dn31.codon.pair$acacgc/cps.aca.cgc)
cps.obs.exp.acacgg<-log(dn31.codon.pair$acacgg/cps.aca.cgg)
cps.obs.exp.acacgt<-log(dn31.codon.pair$acacgt/cps.aca.cgt)
cps.obs.exp.acacta<-log(dn31.codon.pair$acacta/cps.aca.cta)
cps.obs.exp.acactc<-log(dn31.codon.pair$acactc/cps.aca.ctc)
cps.obs.exp.acactg<-log(dn31.codon.pair$acactg/cps.aca.ctg)
cps.obs.exp.acactt<-log(dn31.codon.pair$acactt/cps.aca.ctt)

cps.obs.exp.acagaa<-log(dn31.codon.pair$acagaa/cps.aca.gaa)
cps.obs.exp.acagac<-log(dn31.codon.pair$acagac/cps.aca.gac)
cps.obs.exp.acagag<-log(dn31.codon.pair$acagag/cps.aca.gag)
cps.obs.exp.acagat<-log(dn31.codon.pair$acagat/cps.aca.gat)
cps.obs.exp.acagca<-log(dn31.codon.pair$acagca/cps.aca.gca)
cps.obs.exp.acagcc<-log(dn31.codon.pair$acagcc/cps.aca.gcc)
cps.obs.exp.acagcg<-log(dn31.codon.pair$acagcg/cps.aca.gcg)
cps.obs.exp.acagct<-log(dn31.codon.pair$acagct/cps.aca.gct)
cps.obs.exp.acagga<-log(dn31.codon.pair$acagga/cps.aca.gga)
cps.obs.exp.acaggc<-log(dn31.codon.pair$acaggc/cps.aca.ggc)
cps.obs.exp.acaggg<-log(dn31.codon.pair$acaggg/cps.aca.ggg)
cps.obs.exp.acaggt<-log(dn31.codon.pair$acaggt/cps.aca.ggt)
cps.obs.exp.acagta<-log(dn31.codon.pair$acagta/cps.aca.gta)
cps.obs.exp.acagtc<-log(dn31.codon.pair$acagtc/cps.aca.gtc)
cps.obs.exp.acagtg<-log(dn31.codon.pair$acagtg/cps.aca.gtg)
cps.obs.exp.acagtt<-log(dn31.codon.pair$acagtt/cps.aca.gtt)

#cps.obs.exp.acataa<-log(dn31.codon.pair$acataa/cps.aca.taa)
cps.obs.exp.acatac<-log(dn31.codon.pair$acatac/cps.aca.tac)
#cps.obs.exp.acatag<-log(dn31.codon.pair$acatag/cps.aca.tag)
cps.obs.exp.acatat<-log(dn31.codon.pair$acatat/cps.aca.tat)
cps.obs.exp.acatca<-log(dn31.codon.pair$acatca/cps.aca.tca)
cps.obs.exp.acatcc<-log(dn31.codon.pair$acatcc/cps.aca.tcc)
cps.obs.exp.acatcg<-log(dn31.codon.pair$acatcg/cps.aca.tcg)
cps.obs.exp.acatct<-log(dn31.codon.pair$acatct/cps.aca.tct)
#cps.obs.exp.acatga<-log(dn31.codon.pair$acatga/cps.aca.tga)
cps.obs.exp.acatgc<-log(dn31.codon.pair$acatgc/cps.aca.tgc)
cps.obs.exp.acatgg<-log(dn31.codon.pair$acatgg/cps.aca.tgg)
cps.obs.exp.acatgt<-log(dn31.codon.pair$acatgt/cps.aca.tgt)
cps.obs.exp.acatta<-log(dn31.codon.pair$acatta/cps.aca.tta)
cps.obs.exp.acattc<-log(dn31.codon.pair$acattc/cps.aca.ttc)
cps.obs.exp.acattg<-log(dn31.codon.pair$acattg/cps.aca.ttg)
cps.obs.exp.acattt<-log(dn31.codon.pair$acattt/cps.aca.ttt)









cps.obs.exp.accaaa<-log(dn31.codon.pair$accaaa/cps.acc.aaa)
cps.obs.exp.accaac<-log(dn31.codon.pair$accaac/cps.acc.aac)
cps.obs.exp.accaag<-log(dn31.codon.pair$accaag/cps.acc.aag)
cps.obs.exp.accaat<-log(dn31.codon.pair$accaat/cps.acc.aat)
cps.obs.exp.accaca<-log(dn31.codon.pair$accaca/cps.acc.aca)
cps.obs.exp.accacc<-log(dn31.codon.pair$accacc/cps.acc.acc)
cps.obs.exp.accacg<-log(dn31.codon.pair$accacg/cps.acc.acg)
cps.obs.exp.accact<-log(dn31.codon.pair$accact/cps.acc.act)
cps.obs.exp.accaga<-log(dn31.codon.pair$accaga/cps.acc.aga)
cps.obs.exp.accagc<-log(dn31.codon.pair$accagc/cps.acc.agc)
cps.obs.exp.accagg<-log(dn31.codon.pair$accagg/cps.acc.agg)
cps.obs.exp.accagt<-log(dn31.codon.pair$accagt/cps.acc.agt)
cps.obs.exp.accata<-log(dn31.codon.pair$accata/cps.acc.ata)
cps.obs.exp.accatc<-log(dn31.codon.pair$accatc/cps.acc.atc)
cps.obs.exp.accatg<-log(dn31.codon.pair$accatg/cps.acc.atg)
cps.obs.exp.accatt<-log(dn31.codon.pair$accatt/cps.acc.att)

cps.obs.exp.acccaa<-log(dn31.codon.pair$acccaa/cps.acc.caa)
cps.obs.exp.acccac<-log(dn31.codon.pair$acccac/cps.acc.cac)
cps.obs.exp.acccag<-log(dn31.codon.pair$acccag/cps.acc.cag)
cps.obs.exp.acccat<-log(dn31.codon.pair$acccat/cps.acc.cat)
cps.obs.exp.acccca<-log(dn31.codon.pair$acccca/cps.acc.cca)
cps.obs.exp.accccc<-log(dn31.codon.pair$accccc/cps.acc.ccc)
cps.obs.exp.accccg<-log(dn31.codon.pair$accccg/cps.acc.ccg)
cps.obs.exp.acccct<-log(dn31.codon.pair$acccct/cps.acc.cct)
cps.obs.exp.acccga<-log(dn31.codon.pair$acccga/cps.acc.cga)
cps.obs.exp.acccgc<-log(dn31.codon.pair$acccgc/cps.acc.cgc)
cps.obs.exp.acccgg<-log(dn31.codon.pair$acccgg/cps.acc.cgg)
cps.obs.exp.acccgt<-log(dn31.codon.pair$acccgt/cps.acc.cgt)
cps.obs.exp.acccta<-log(dn31.codon.pair$acccta/cps.acc.cta)
cps.obs.exp.accctc<-log(dn31.codon.pair$accctc/cps.acc.ctc)
cps.obs.exp.accctg<-log(dn31.codon.pair$accctg/cps.acc.ctg)
cps.obs.exp.accctt<-log(dn31.codon.pair$accctt/cps.acc.ctt)

cps.obs.exp.accgaa<-log(dn31.codon.pair$accgaa/cps.acc.gaa)
cps.obs.exp.accgac<-log(dn31.codon.pair$accgac/cps.acc.gac)
cps.obs.exp.accgag<-log(dn31.codon.pair$accgag/cps.acc.gag)
cps.obs.exp.accgat<-log(dn31.codon.pair$accgat/cps.acc.gat)
cps.obs.exp.accgca<-log(dn31.codon.pair$accgca/cps.acc.gca)
cps.obs.exp.accgcc<-log(dn31.codon.pair$accgcc/cps.acc.gcc)
cps.obs.exp.accgcg<-log(dn31.codon.pair$accgcg/cps.acc.gcg)
cps.obs.exp.accgct<-log(dn31.codon.pair$accgct/cps.acc.gct)
cps.obs.exp.accgga<-log(dn31.codon.pair$accgga/cps.acc.gga)
cps.obs.exp.accggc<-log(dn31.codon.pair$accggc/cps.acc.ggc)
cps.obs.exp.accggg<-log(dn31.codon.pair$accggg/cps.acc.ggg)
cps.obs.exp.accggt<-log(dn31.codon.pair$accggt/cps.acc.ggt)
cps.obs.exp.accgta<-log(dn31.codon.pair$accgta/cps.acc.gta)
cps.obs.exp.accgtc<-log(dn31.codon.pair$accgtc/cps.acc.gtc)
cps.obs.exp.accgtg<-log(dn31.codon.pair$accgtg/cps.acc.gtg)
cps.obs.exp.accgtt<-log(dn31.codon.pair$accgtt/cps.acc.gtt)

#cps.obs.exp.acctaa<-log(dn31.codon.pair$acctaa/cps.acc.taa)
cps.obs.exp.acctac<-log(dn31.codon.pair$acctac/cps.acc.tac)
#cps.obs.exp.acctag<-log(dn31.codon.pair$acctag/cps.acc.tag)
cps.obs.exp.acctat<-log(dn31.codon.pair$acctat/cps.acc.tat)
cps.obs.exp.acctca<-log(dn31.codon.pair$acctca/cps.acc.tca)
cps.obs.exp.acctcc<-log(dn31.codon.pair$acctcc/cps.acc.tcc)
cps.obs.exp.acctcg<-log(dn31.codon.pair$acctcg/cps.acc.tcg)
cps.obs.exp.acctct<-log(dn31.codon.pair$acctct/cps.acc.tct)
#cps.obs.exp.acctga<-log(dn31.codon.pair$acctga/cps.acc.tga)
cps.obs.exp.acctgc<-log(dn31.codon.pair$acctgc/cps.acc.tgc)
cps.obs.exp.acctgg<-log(dn31.codon.pair$acctgg/cps.acc.tgg)
cps.obs.exp.acctgt<-log(dn31.codon.pair$acctgt/cps.acc.tgt)
cps.obs.exp.acctta<-log(dn31.codon.pair$acctta/cps.acc.tta)
cps.obs.exp.accttc<-log(dn31.codon.pair$accttc/cps.acc.ttc)
cps.obs.exp.accttg<-log(dn31.codon.pair$accttg/cps.acc.ttg)
cps.obs.exp.accttt<-log(dn31.codon.pair$accttt/cps.acc.ttt)








cps.obs.exp.acgaaa<-log(dn31.codon.pair$acgaaa/cps.acg.aaa)
cps.obs.exp.acgaac<-log(dn31.codon.pair$acgaac/cps.acg.aac)
cps.obs.exp.acgaag<-log(dn31.codon.pair$acgaag/cps.acg.aag)
cps.obs.exp.acgaat<-log(dn31.codon.pair$acgaat/cps.acg.aat)
cps.obs.exp.acgaca<-log(dn31.codon.pair$acgaca/cps.acg.aca)
cps.obs.exp.acgacc<-log(dn31.codon.pair$acgacc/cps.acg.acc)
cps.obs.exp.acgacg<-log(dn31.codon.pair$acgacg/cps.acg.acg)
cps.obs.exp.acgact<-log(dn31.codon.pair$acgact/cps.acg.act)
cps.obs.exp.acgaga<-log(dn31.codon.pair$acgaga/cps.acg.aga)
cps.obs.exp.acgagc<-log(dn31.codon.pair$acgagc/cps.acg.agc)
cps.obs.exp.acgagg<-log(dn31.codon.pair$acgagg/cps.acg.agg)
cps.obs.exp.acgagt<-log(dn31.codon.pair$acgagt/cps.acg.agt)
cps.obs.exp.acgata<-log(dn31.codon.pair$acgata/cps.acg.ata)
cps.obs.exp.acgatc<-log(dn31.codon.pair$acgatc/cps.acg.atc)
cps.obs.exp.acgatg<-log(dn31.codon.pair$acgatg/cps.acg.atg)
cps.obs.exp.acgatt<-log(dn31.codon.pair$acgatt/cps.acg.att)

cps.obs.exp.acgcaa<-log(dn31.codon.pair$acgcaa/cps.acg.caa)
cps.obs.exp.acgcac<-log(dn31.codon.pair$acgcac/cps.acg.cac)
cps.obs.exp.acgcag<-log(dn31.codon.pair$acgcag/cps.acg.cag)
cps.obs.exp.acgcat<-log(dn31.codon.pair$acgcat/cps.acg.cat)
cps.obs.exp.acgcca<-log(dn31.codon.pair$acgcca/cps.acg.cca)
cps.obs.exp.acgccc<-log(dn31.codon.pair$acgccc/cps.acg.ccc)
cps.obs.exp.acgccg<-log(dn31.codon.pair$acgccg/cps.acg.ccg)
cps.obs.exp.acgcct<-log(dn31.codon.pair$acgcct/cps.acg.cct)
cps.obs.exp.acgcga<-log(dn31.codon.pair$acgcga/cps.acg.cga)
cps.obs.exp.acgcgc<-log(dn31.codon.pair$acgcgc/cps.acg.cgc)
cps.obs.exp.acgcgg<-log(dn31.codon.pair$acgcgg/cps.acg.cgg)
cps.obs.exp.acgcgt<-log(dn31.codon.pair$acgcgt/cps.acg.cgt)
cps.obs.exp.acgcta<-log(dn31.codon.pair$acgcta/cps.acg.cta)
cps.obs.exp.acgctc<-log(dn31.codon.pair$acgctc/cps.acg.ctc)
cps.obs.exp.acgctg<-log(dn31.codon.pair$acgctg/cps.acg.ctg)
cps.obs.exp.acgctt<-log(dn31.codon.pair$acgctt/cps.acg.ctt)

cps.obs.exp.acggaa<-log(dn31.codon.pair$acggaa/cps.acg.gaa)
cps.obs.exp.acggac<-log(dn31.codon.pair$acggac/cps.acg.gac)
cps.obs.exp.acggag<-log(dn31.codon.pair$acggag/cps.acg.gag)
cps.obs.exp.acggat<-log(dn31.codon.pair$acggat/cps.acg.gat)
cps.obs.exp.acggca<-log(dn31.codon.pair$acggca/cps.acg.gca)
cps.obs.exp.acggcc<-log(dn31.codon.pair$acggcc/cps.acg.gcc)
cps.obs.exp.acggcg<-log(dn31.codon.pair$acggcg/cps.acg.gcg)
cps.obs.exp.acggct<-log(dn31.codon.pair$acggct/cps.acg.gct)
cps.obs.exp.acggga<-log(dn31.codon.pair$acggga/cps.acg.gga)
cps.obs.exp.acgggc<-log(dn31.codon.pair$acgggc/cps.acg.ggc)
cps.obs.exp.acgggg<-log(dn31.codon.pair$acgggg/cps.acg.ggg)
cps.obs.exp.acgggt<-log(dn31.codon.pair$acgggt/cps.acg.ggt)
cps.obs.exp.acggta<-log(dn31.codon.pair$acggta/cps.acg.gta)
cps.obs.exp.acggtc<-log(dn31.codon.pair$acggtc/cps.acg.gtc)
cps.obs.exp.acggtg<-log(dn31.codon.pair$acggtg/cps.acg.gtg)
cps.obs.exp.acggtt<-log(dn31.codon.pair$acggtt/cps.acg.gtt)

#cps.obs.exp.acgtaa<-log(dn31.codon.pair$acgtaa/cps.acg.taa)
cps.obs.exp.acgtac<-log(dn31.codon.pair$acgtac/cps.acg.tac)
#cps.obs.exp.acgtag<-log(dn31.codon.pair$acgtag/cps.acg.tag)
cps.obs.exp.acgtat<-log(dn31.codon.pair$acgtat/cps.acg.tat)
cps.obs.exp.acgtca<-log(dn31.codon.pair$acgtca/cps.acg.tca)
cps.obs.exp.acgtcc<-log(dn31.codon.pair$acgtcc/cps.acg.tcc)
cps.obs.exp.acgtcg<-log(dn31.codon.pair$acgtcg/cps.acg.tcg)
cps.obs.exp.acgtct<-log(dn31.codon.pair$acgtct/cps.acg.tct)
#cps.obs.exp.acgtga<-log(dn31.codon.pair$acgtga/cps.acg.tga)
cps.obs.exp.acgtgc<-log(dn31.codon.pair$acgtgc/cps.acg.tgc)
cps.obs.exp.acgtgg<-log(dn31.codon.pair$acgtgg/cps.acg.tgg)
cps.obs.exp.acgtgt<-log(dn31.codon.pair$acgtgt/cps.acg.tgt)
cps.obs.exp.acgtta<-log(dn31.codon.pair$acgtta/cps.acg.tta)
cps.obs.exp.acgttc<-log(dn31.codon.pair$acgttc/cps.acg.ttc)
cps.obs.exp.acgttg<-log(dn31.codon.pair$acgttg/cps.acg.ttg)
cps.obs.exp.acgttt<-log(dn31.codon.pair$acgttt/cps.acg.ttt)








cps.obs.exp.actaaa<-log(dn31.codon.pair$actaaa/cps.act.aaa)
cps.obs.exp.actaac<-log(dn31.codon.pair$actaac/cps.act.aac)
cps.obs.exp.actaag<-log(dn31.codon.pair$actaag/cps.act.aag)
cps.obs.exp.actaat<-log(dn31.codon.pair$actaat/cps.act.aat)
cps.obs.exp.actaca<-log(dn31.codon.pair$actaca/cps.act.aca)
cps.obs.exp.actacc<-log(dn31.codon.pair$actacc/cps.act.acc)
cps.obs.exp.actacg<-log(dn31.codon.pair$actacg/cps.act.acg)
cps.obs.exp.actact<-log(dn31.codon.pair$actact/cps.act.act)
cps.obs.exp.actaga<-log(dn31.codon.pair$actaga/cps.act.aga)
cps.obs.exp.actagc<-log(dn31.codon.pair$actagc/cps.act.agc)
cps.obs.exp.actagg<-log(dn31.codon.pair$actagg/cps.act.agg)
cps.obs.exp.actagt<-log(dn31.codon.pair$actagt/cps.act.agt)
cps.obs.exp.actata<-log(dn31.codon.pair$actata/cps.act.ata)
cps.obs.exp.actatc<-log(dn31.codon.pair$actatc/cps.act.atc)
cps.obs.exp.actatg<-log(dn31.codon.pair$actatg/cps.act.atg)
cps.obs.exp.actatt<-log(dn31.codon.pair$actatt/cps.act.att)

cps.obs.exp.actcaa<-log(dn31.codon.pair$actcaa/cps.act.caa)
cps.obs.exp.actcac<-log(dn31.codon.pair$actcac/cps.act.cac)
cps.obs.exp.actcag<-log(dn31.codon.pair$actcag/cps.act.cag)
cps.obs.exp.actcat<-log(dn31.codon.pair$actcat/cps.act.cat)
cps.obs.exp.actcca<-log(dn31.codon.pair$actcca/cps.act.cca)
cps.obs.exp.actccc<-log(dn31.codon.pair$actccc/cps.act.ccc)
cps.obs.exp.actccg<-log(dn31.codon.pair$actccg/cps.act.ccg)
cps.obs.exp.actcct<-log(dn31.codon.pair$actcct/cps.act.cct)
cps.obs.exp.actcga<-log(dn31.codon.pair$actcga/cps.act.cga)
cps.obs.exp.actcgc<-log(dn31.codon.pair$actcgc/cps.act.cgc)
cps.obs.exp.actcgg<-log(dn31.codon.pair$actcgg/cps.act.cgg)
cps.obs.exp.actcgt<-log(dn31.codon.pair$actcgt/cps.act.cgt)
cps.obs.exp.actcta<-log(dn31.codon.pair$actcta/cps.act.cta)
cps.obs.exp.actctc<-log(dn31.codon.pair$actctc/cps.act.ctc)
cps.obs.exp.actctg<-log(dn31.codon.pair$actctg/cps.act.ctg)
cps.obs.exp.actctt<-log(dn31.codon.pair$actctt/cps.act.ctt)

cps.obs.exp.actgaa<-log(dn31.codon.pair$actgaa/cps.act.gaa)
cps.obs.exp.actgac<-log(dn31.codon.pair$actgac/cps.act.gac)
cps.obs.exp.actgag<-log(dn31.codon.pair$actgag/cps.act.gag)
cps.obs.exp.actgat<-log(dn31.codon.pair$actgat/cps.act.gat)
cps.obs.exp.actgca<-log(dn31.codon.pair$actgca/cps.act.gca)
cps.obs.exp.actgcc<-log(dn31.codon.pair$actgcc/cps.act.gcc)
cps.obs.exp.actgcg<-log(dn31.codon.pair$actgcg/cps.act.gcg)
cps.obs.exp.actgct<-log(dn31.codon.pair$actgct/cps.act.gct)
cps.obs.exp.actgga<-log(dn31.codon.pair$actgga/cps.act.gga)
cps.obs.exp.actggc<-log(dn31.codon.pair$actggc/cps.act.ggc)
cps.obs.exp.actggg<-log(dn31.codon.pair$actggg/cps.act.ggg)
cps.obs.exp.actggt<-log(dn31.codon.pair$actggt/cps.act.ggt)
cps.obs.exp.actgta<-log(dn31.codon.pair$actgta/cps.act.gta)
cps.obs.exp.actgtc<-log(dn31.codon.pair$actgtc/cps.act.gtc)
cps.obs.exp.actgtg<-log(dn31.codon.pair$actgtg/cps.act.gtg)
cps.obs.exp.actgtt<-log(dn31.codon.pair$actgtt/cps.act.gtt)

#cps.obs.exp.acttaa<-log(dn31.codon.pair$acttaa/cps.act.taa)
cps.obs.exp.acttac<-log(dn31.codon.pair$acttac/cps.act.tac)
#cps.obs.exp.acttag<-log(dn31.codon.pair$acttag/cps.act.tag)
cps.obs.exp.acttat<-log(dn31.codon.pair$acttat/cps.act.tat)
cps.obs.exp.acttca<-log(dn31.codon.pair$acttca/cps.act.tca)
cps.obs.exp.acttcc<-log(dn31.codon.pair$acttcc/cps.act.tcc)
cps.obs.exp.acttcg<-log(dn31.codon.pair$acttcg/cps.act.tcg)
cps.obs.exp.acttct<-log(dn31.codon.pair$acttct/cps.act.tct)
#cps.obs.exp.acttga<-log(dn31.codon.pair$acttga/cps.act.tga)
cps.obs.exp.acttgc<-log(dn31.codon.pair$acttgc/cps.act.tgc)
cps.obs.exp.acttgg<-log(dn31.codon.pair$acttgg/cps.act.tgg)
cps.obs.exp.acttgt<-log(dn31.codon.pair$acttgt/cps.act.tgt)
cps.obs.exp.acttta<-log(dn31.codon.pair$acttta/cps.act.tta)
cps.obs.exp.actttc<-log(dn31.codon.pair$actttc/cps.act.ttc)
cps.obs.exp.actttg<-log(dn31.codon.pair$actttg/cps.act.ttg)
cps.obs.exp.actttt<-log(dn31.codon.pair$actttt/cps.act.ttt)



















cps.obs.exp.agaaaa<-log(dn31.codon.pair$agaaaa/cps.aga.aaa)
cps.obs.exp.agaaac<-log(dn31.codon.pair$agaaac/cps.aga.aac)
cps.obs.exp.agaaag<-log(dn31.codon.pair$agaaag/cps.aga.aag)
cps.obs.exp.agaaat<-log(dn31.codon.pair$agaaat/cps.aga.aat)
cps.obs.exp.agaaca<-log(dn31.codon.pair$agaaca/cps.aga.aca)
cps.obs.exp.agaacc<-log(dn31.codon.pair$agaacc/cps.aga.acc)
cps.obs.exp.agaacg<-log(dn31.codon.pair$agaacg/cps.aga.acg)
cps.obs.exp.agaact<-log(dn31.codon.pair$agaact/cps.aga.act)
cps.obs.exp.agaaga<-log(dn31.codon.pair$agaaga/cps.aga.aga)
cps.obs.exp.agaagc<-log(dn31.codon.pair$agaagc/cps.aga.agc)
cps.obs.exp.agaagg<-log(dn31.codon.pair$agaagg/cps.aga.agg)
cps.obs.exp.agaagt<-log(dn31.codon.pair$agaagt/cps.aga.agt)
cps.obs.exp.agaata<-log(dn31.codon.pair$agaata/cps.aga.ata)
cps.obs.exp.agaatc<-log(dn31.codon.pair$agaatc/cps.aga.atc)
cps.obs.exp.agaatg<-log(dn31.codon.pair$agaatg/cps.aga.atg)
cps.obs.exp.agaatt<-log(dn31.codon.pair$agaatt/cps.aga.att)

cps.obs.exp.agacaa<-log(dn31.codon.pair$agacaa/cps.aga.caa)
cps.obs.exp.agacac<-log(dn31.codon.pair$agacac/cps.aga.cac)
cps.obs.exp.agacag<-log(dn31.codon.pair$agacag/cps.aga.cag)
cps.obs.exp.agacat<-log(dn31.codon.pair$agacat/cps.aga.cat)
cps.obs.exp.agacca<-log(dn31.codon.pair$agacca/cps.aga.cca)
cps.obs.exp.agaccc<-log(dn31.codon.pair$agaccc/cps.aga.ccc)
cps.obs.exp.agaccg<-log(dn31.codon.pair$agaccg/cps.aga.ccg)
cps.obs.exp.agacct<-log(dn31.codon.pair$agacct/cps.aga.cct)
cps.obs.exp.agacga<-log(dn31.codon.pair$agacga/cps.aga.cga)
cps.obs.exp.agacgc<-log(dn31.codon.pair$agacgc/cps.aga.cgc)
cps.obs.exp.agacgg<-log(dn31.codon.pair$agacgg/cps.aga.cgg)
cps.obs.exp.agacgt<-log(dn31.codon.pair$agacgt/cps.aga.cgt)
cps.obs.exp.agacta<-log(dn31.codon.pair$agacta/cps.aga.cta)
cps.obs.exp.agactc<-log(dn31.codon.pair$agactc/cps.aga.ctc)
cps.obs.exp.agactg<-log(dn31.codon.pair$agactg/cps.aga.ctg)
cps.obs.exp.agactt<-log(dn31.codon.pair$agactt/cps.aga.ctt)

cps.obs.exp.agagaa<-log(dn31.codon.pair$agagaa/cps.aga.gaa)
cps.obs.exp.agagac<-log(dn31.codon.pair$agagac/cps.aga.gac)
cps.obs.exp.agagag<-log(dn31.codon.pair$agagag/cps.aga.gag)
cps.obs.exp.agagat<-log(dn31.codon.pair$agagat/cps.aga.gat)
cps.obs.exp.agagca<-log(dn31.codon.pair$agagca/cps.aga.gca)
cps.obs.exp.agagcc<-log(dn31.codon.pair$agagcc/cps.aga.gcc)
cps.obs.exp.agagcg<-log(dn31.codon.pair$agagcg/cps.aga.gcg)
cps.obs.exp.agagct<-log(dn31.codon.pair$agagct/cps.aga.gct)
cps.obs.exp.agagga<-log(dn31.codon.pair$agagga/cps.aga.gga)
cps.obs.exp.agaggc<-log(dn31.codon.pair$agaggc/cps.aga.ggc)
cps.obs.exp.agaggg<-log(dn31.codon.pair$agaggg/cps.aga.ggg)
cps.obs.exp.agaggt<-log(dn31.codon.pair$agaggt/cps.aga.ggt)
cps.obs.exp.agagta<-log(dn31.codon.pair$agagta/cps.aga.gta)
cps.obs.exp.agagtc<-log(dn31.codon.pair$agagtc/cps.aga.gtc)
cps.obs.exp.agagtg<-log(dn31.codon.pair$agagtg/cps.aga.gtg)
cps.obs.exp.agagtt<-log(dn31.codon.pair$agagtt/cps.aga.gtt)

#cps.obs.exp.agataa<-log(dn31.codon.pair$agataa/cps.aga.taa)
cps.obs.exp.agatac<-log(dn31.codon.pair$agatac/cps.aga.tac)
#cps.obs.exp.agatag<-log(dn31.codon.pair$agatag/cps.aga.tag)
cps.obs.exp.agatat<-log(dn31.codon.pair$agatat/cps.aga.tat)
cps.obs.exp.agatca<-log(dn31.codon.pair$agatca/cps.aga.tca)
cps.obs.exp.agatcc<-log(dn31.codon.pair$agatcc/cps.aga.tcc)
cps.obs.exp.agatcg<-log(dn31.codon.pair$agatcg/cps.aga.tcg)
cps.obs.exp.agatct<-log(dn31.codon.pair$agatct/cps.aga.tct)
#cps.obs.exp.agatga<-log(dn31.codon.pair$agatga/cps.aga.tga)
cps.obs.exp.agatgc<-log(dn31.codon.pair$agatgc/cps.aga.tgc)
cps.obs.exp.agatgg<-log(dn31.codon.pair$agatgg/cps.aga.tgg)
cps.obs.exp.agatgt<-log(dn31.codon.pair$agatgt/cps.aga.tgt)
cps.obs.exp.agatta<-log(dn31.codon.pair$agatta/cps.aga.tta)
cps.obs.exp.agattc<-log(dn31.codon.pair$agattc/cps.aga.ttc)
cps.obs.exp.agattg<-log(dn31.codon.pair$agattg/cps.aga.ttg)
cps.obs.exp.agattt<-log(dn31.codon.pair$agattt/cps.aga.ttt)







cps.obs.exp.agcaaa<-log(dn31.codon.pair$agcaaa/cps.agc.aaa)
cps.obs.exp.agcaac<-log(dn31.codon.pair$agcaac/cps.agc.aac)
cps.obs.exp.agcaag<-log(dn31.codon.pair$agcaag/cps.agc.aag)
cps.obs.exp.agcaat<-log(dn31.codon.pair$agcaat/cps.agc.aat)
cps.obs.exp.agcaca<-log(dn31.codon.pair$agcaca/cps.agc.aca)
cps.obs.exp.agcacc<-log(dn31.codon.pair$agcacc/cps.agc.acc)
cps.obs.exp.agcacg<-log(dn31.codon.pair$agcacg/cps.agc.acg)
cps.obs.exp.agcact<-log(dn31.codon.pair$agcact/cps.agc.act)
cps.obs.exp.agcaga<-log(dn31.codon.pair$agcaga/cps.agc.aga)
cps.obs.exp.agcagc<-log(dn31.codon.pair$agcagc/cps.agc.agc)
cps.obs.exp.agcagg<-log(dn31.codon.pair$agcagg/cps.agc.agg)
cps.obs.exp.agcagt<-log(dn31.codon.pair$agcagt/cps.agc.agt)
cps.obs.exp.agcata<-log(dn31.codon.pair$agcata/cps.agc.ata)
cps.obs.exp.agcatc<-log(dn31.codon.pair$agcatc/cps.agc.atc)
cps.obs.exp.agcatg<-log(dn31.codon.pair$agcatg/cps.agc.atg)
cps.obs.exp.agcatt<-log(dn31.codon.pair$agcatt/cps.agc.att)

cps.obs.exp.agccaa<-log(dn31.codon.pair$agccaa/cps.agc.caa)
cps.obs.exp.agccac<-log(dn31.codon.pair$agccac/cps.agc.cac)
cps.obs.exp.agccag<-log(dn31.codon.pair$agccag/cps.agc.cag)
cps.obs.exp.agccat<-log(dn31.codon.pair$agccat/cps.agc.cat)
cps.obs.exp.agccca<-log(dn31.codon.pair$agccca/cps.agc.cca)
cps.obs.exp.agcccc<-log(dn31.codon.pair$agcccc/cps.agc.ccc)
cps.obs.exp.agcccg<-log(dn31.codon.pair$agcccg/cps.agc.ccg)
cps.obs.exp.agccct<-log(dn31.codon.pair$agccct/cps.agc.cct)
cps.obs.exp.agccga<-log(dn31.codon.pair$agccga/cps.agc.cga)
cps.obs.exp.agccgc<-log(dn31.codon.pair$agccgc/cps.agc.cgc)
cps.obs.exp.agccgg<-log(dn31.codon.pair$agccgg/cps.agc.cgg)
cps.obs.exp.agccgt<-log(dn31.codon.pair$agccgt/cps.agc.cgt)
cps.obs.exp.agccta<-log(dn31.codon.pair$agccta/cps.agc.cta)
cps.obs.exp.agcctc<-log(dn31.codon.pair$agcctc/cps.agc.ctc)
cps.obs.exp.agcctg<-log(dn31.codon.pair$agcctg/cps.agc.ctg)
cps.obs.exp.agcctt<-log(dn31.codon.pair$agcctt/cps.agc.ctt)

cps.obs.exp.agcgaa<-log(dn31.codon.pair$agcgaa/cps.agc.gaa)
cps.obs.exp.agcgac<-log(dn31.codon.pair$agcgac/cps.agc.gac)
cps.obs.exp.agcgag<-log(dn31.codon.pair$agcgag/cps.agc.gag)
cps.obs.exp.agcgat<-log(dn31.codon.pair$agcgat/cps.agc.gat)
cps.obs.exp.agcgca<-log(dn31.codon.pair$agcgca/cps.agc.gca)
cps.obs.exp.agcgcc<-log(dn31.codon.pair$agcgcc/cps.agc.gcc)
cps.obs.exp.agcgcg<-log(dn31.codon.pair$agcgcg/cps.agc.gcg)
cps.obs.exp.agcgct<-log(dn31.codon.pair$agcgct/cps.agc.gct)
cps.obs.exp.agcgga<-log(dn31.codon.pair$agcgga/cps.agc.gga)
cps.obs.exp.agcggc<-log(dn31.codon.pair$agcggc/cps.agc.ggc)
cps.obs.exp.agcggg<-log(dn31.codon.pair$agcggg/cps.agc.ggg)
cps.obs.exp.agcggt<-log(dn31.codon.pair$agcggt/cps.agc.ggt)
cps.obs.exp.agcgta<-log(dn31.codon.pair$agcgta/cps.agc.gta)
cps.obs.exp.agcgtc<-log(dn31.codon.pair$agcgtc/cps.agc.gtc)
cps.obs.exp.agcgtg<-log(dn31.codon.pair$agcgtg/cps.agc.gtg)
cps.obs.exp.agcgtt<-log(dn31.codon.pair$agcgtt/cps.agc.gtt)

#cps.obs.exp.agctaa<-log(dn31.codon.pair$agctaa/cps.agc.taa)
cps.obs.exp.agctac<-log(dn31.codon.pair$agctac/cps.agc.tac)
#cps.obs.exp.agctag<-log(dn31.codon.pair$agctag/cps.agc.tag)
cps.obs.exp.agctat<-log(dn31.codon.pair$agctat/cps.agc.tat)
cps.obs.exp.agctca<-log(dn31.codon.pair$agctca/cps.agc.tca)
cps.obs.exp.agctcc<-log(dn31.codon.pair$agctcc/cps.agc.tcc)
cps.obs.exp.agctcg<-log(dn31.codon.pair$agctcg/cps.agc.tcg)
cps.obs.exp.agctct<-log(dn31.codon.pair$agctct/cps.agc.tct)
#cps.obs.exp.agctga<-log(dn31.codon.pair$agctga/cps.agc.tga)
cps.obs.exp.agctgc<-log(dn31.codon.pair$agctgc/cps.agc.tgc)
cps.obs.exp.agctgg<-log(dn31.codon.pair$agctgg/cps.agc.tgg)
cps.obs.exp.agctgt<-log(dn31.codon.pair$agctgt/cps.agc.tgt)
cps.obs.exp.agctta<-log(dn31.codon.pair$agctta/cps.agc.tta)
cps.obs.exp.agcttc<-log(dn31.codon.pair$agcttc/cps.agc.ttc)
cps.obs.exp.agcttg<-log(dn31.codon.pair$agcttg/cps.agc.ttg)
cps.obs.exp.agcttt<-log(dn31.codon.pair$agcttt/cps.agc.ttt)









cps.obs.exp.aggaaa<-log(dn31.codon.pair$aggaaa/cps.agg.aaa)
cps.obs.exp.aggaac<-log(dn31.codon.pair$aggaac/cps.agg.aac)
cps.obs.exp.aggaag<-log(dn31.codon.pair$aggaag/cps.agg.aag)
cps.obs.exp.aggaat<-log(dn31.codon.pair$aggaat/cps.agg.aat)
cps.obs.exp.aggaca<-log(dn31.codon.pair$aggaca/cps.agg.aca)
cps.obs.exp.aggacc<-log(dn31.codon.pair$aggacc/cps.agg.acc)
cps.obs.exp.aggacg<-log(dn31.codon.pair$aggacg/cps.agg.acg)
cps.obs.exp.aggact<-log(dn31.codon.pair$aggact/cps.agg.act)
cps.obs.exp.aggaga<-log(dn31.codon.pair$aggaga/cps.agg.aga)
cps.obs.exp.aggagc<-log(dn31.codon.pair$aggagc/cps.agg.agc)
cps.obs.exp.aggagg<-log(dn31.codon.pair$aggagg/cps.agg.agg)
cps.obs.exp.aggagt<-log(dn31.codon.pair$aggagt/cps.agg.agt)
cps.obs.exp.aggata<-log(dn31.codon.pair$aggata/cps.agg.ata)
cps.obs.exp.aggatc<-log(dn31.codon.pair$aggatc/cps.agg.atc)
cps.obs.exp.aggatg<-log(dn31.codon.pair$aggatg/cps.agg.atg)
cps.obs.exp.aggatt<-log(dn31.codon.pair$aggatt/cps.agg.att)

cps.obs.exp.aggcaa<-log(dn31.codon.pair$aggcaa/cps.agg.caa)
cps.obs.exp.aggcac<-log(dn31.codon.pair$aggcac/cps.agg.cac)
cps.obs.exp.aggcag<-log(dn31.codon.pair$aggcag/cps.agg.cag)
cps.obs.exp.aggcat<-log(dn31.codon.pair$aggcat/cps.agg.cat)
cps.obs.exp.aggcca<-log(dn31.codon.pair$aggcca/cps.agg.cca)
cps.obs.exp.aggccc<-log(dn31.codon.pair$aggccc/cps.agg.ccc)
cps.obs.exp.aggccg<-log(dn31.codon.pair$aggccg/cps.agg.ccg)
cps.obs.exp.aggcct<-log(dn31.codon.pair$aggcct/cps.agg.cct)
cps.obs.exp.aggcga<-log(dn31.codon.pair$aggcga/cps.agg.cga)
cps.obs.exp.aggcgc<-log(dn31.codon.pair$aggcgc/cps.agg.cgc)
cps.obs.exp.aggcgg<-log(dn31.codon.pair$aggcgg/cps.agg.cgg)
cps.obs.exp.aggcgt<-log(dn31.codon.pair$aggcgt/cps.agg.cgt)
cps.obs.exp.aggcta<-log(dn31.codon.pair$aggcta/cps.agg.cta)
cps.obs.exp.aggctc<-log(dn31.codon.pair$aggctc/cps.agg.ctc)
cps.obs.exp.aggctg<-log(dn31.codon.pair$aggctg/cps.agg.ctg)
cps.obs.exp.aggctt<-log(dn31.codon.pair$aggctt/cps.agg.ctt)

cps.obs.exp.agggaa<-log(dn31.codon.pair$agggaa/cps.agg.gaa)
cps.obs.exp.agggac<-log(dn31.codon.pair$agggac/cps.agg.gac)
cps.obs.exp.agggag<-log(dn31.codon.pair$agggag/cps.agg.gag)
cps.obs.exp.agggat<-log(dn31.codon.pair$agggat/cps.agg.gat)
cps.obs.exp.agggca<-log(dn31.codon.pair$agggca/cps.agg.gca)
cps.obs.exp.agggcc<-log(dn31.codon.pair$agggcc/cps.agg.gcc)
cps.obs.exp.agggcg<-log(dn31.codon.pair$agggcg/cps.agg.gcg)
cps.obs.exp.agggct<-log(dn31.codon.pair$agggct/cps.agg.gct)
cps.obs.exp.agggga<-log(dn31.codon.pair$agggga/cps.agg.gga)
cps.obs.exp.aggggc<-log(dn31.codon.pair$aggggc/cps.agg.ggc)
cps.obs.exp.aggggg<-log(dn31.codon.pair$aggggg/cps.agg.ggg)
cps.obs.exp.aggggt<-log(dn31.codon.pair$aggggt/cps.agg.ggt)
cps.obs.exp.agggta<-log(dn31.codon.pair$agggta/cps.agg.gta)
cps.obs.exp.agggtc<-log(dn31.codon.pair$agggtc/cps.agg.gtc)
cps.obs.exp.agggtg<-log(dn31.codon.pair$agggtg/cps.agg.gtg)
cps.obs.exp.agggtt<-log(dn31.codon.pair$agggtt/cps.agg.gtt)

#cps.obs.exp.aggtaa<-log(dn31.codon.pair$aggtaa/cps.agg.taa)
cps.obs.exp.aggtac<-log(dn31.codon.pair$aggtac/cps.agg.tac)
#cps.obs.exp.aggtag<-log(dn31.codon.pair$aggtag/cps.agg.tag)
cps.obs.exp.aggtat<-log(dn31.codon.pair$aggtat/cps.agg.tat)
cps.obs.exp.aggtca<-log(dn31.codon.pair$aggtca/cps.agg.tca)
cps.obs.exp.aggtcc<-log(dn31.codon.pair$aggtcc/cps.agg.tcc)
cps.obs.exp.aggtcg<-log(dn31.codon.pair$aggtcg/cps.agg.tcg)
cps.obs.exp.aggtct<-log(dn31.codon.pair$aggtct/cps.agg.tct)
#cps.obs.exp.aggtga<-log(dn31.codon.pair$aggtga/cps.agg.tga)
cps.obs.exp.aggtgc<-log(dn31.codon.pair$aggtgc/cps.agg.tgc)
cps.obs.exp.aggtgg<-log(dn31.codon.pair$aggtgg/cps.agg.tgg)
cps.obs.exp.aggtgt<-log(dn31.codon.pair$aggtgt/cps.agg.tgt)
cps.obs.exp.aggtta<-log(dn31.codon.pair$aggtta/cps.agg.tta)
cps.obs.exp.aggttc<-log(dn31.codon.pair$aggttc/cps.agg.ttc)
cps.obs.exp.aggttg<-log(dn31.codon.pair$aggttg/cps.agg.ttg)
cps.obs.exp.aggttt<-log(dn31.codon.pair$aggttt/cps.agg.ttt)







cps.obs.exp.agtaaa<-log(dn31.codon.pair$agtaaa/cps.agt.aaa)
cps.obs.exp.agtaac<-log(dn31.codon.pair$agtaac/cps.agt.aac)
cps.obs.exp.agtaag<-log(dn31.codon.pair$agtaag/cps.agt.aag)
cps.obs.exp.agtaat<-log(dn31.codon.pair$agtaat/cps.agt.aat)
cps.obs.exp.agtaca<-log(dn31.codon.pair$agtaca/cps.agt.aca)
cps.obs.exp.agtacc<-log(dn31.codon.pair$agtacc/cps.agt.acc)
cps.obs.exp.agtacg<-log(dn31.codon.pair$agtacg/cps.agt.acg)
cps.obs.exp.agtact<-log(dn31.codon.pair$agtact/cps.agt.act)
cps.obs.exp.agtaga<-log(dn31.codon.pair$agtaga/cps.agt.aga)
cps.obs.exp.agtagc<-log(dn31.codon.pair$agtagc/cps.agt.agc)
cps.obs.exp.agtagg<-log(dn31.codon.pair$agtagg/cps.agt.agg)
cps.obs.exp.agtagt<-log(dn31.codon.pair$agtagt/cps.agt.agt)
cps.obs.exp.agtata<-log(dn31.codon.pair$agtata/cps.agt.ata)
cps.obs.exp.agtatc<-log(dn31.codon.pair$agtatc/cps.agt.atc)
cps.obs.exp.agtatg<-log(dn31.codon.pair$agtatg/cps.agt.atg)
cps.obs.exp.agtatt<-log(dn31.codon.pair$agtatt/cps.agt.att)

cps.obs.exp.agtcaa<-log(dn31.codon.pair$agtcaa/cps.agt.caa)
cps.obs.exp.agtcac<-log(dn31.codon.pair$agtcac/cps.agt.cac)
cps.obs.exp.agtcag<-log(dn31.codon.pair$agtcag/cps.agt.cag)
cps.obs.exp.agtcat<-log(dn31.codon.pair$agtcat/cps.agt.cat)
cps.obs.exp.agtcca<-log(dn31.codon.pair$agtcca/cps.agt.cca)
cps.obs.exp.agtccc<-log(dn31.codon.pair$agtccc/cps.agt.ccc)
cps.obs.exp.agtccg<-log(dn31.codon.pair$agtccg/cps.agt.ccg)
cps.obs.exp.agtcct<-log(dn31.codon.pair$agtcct/cps.agt.cct)
cps.obs.exp.agtcga<-log(dn31.codon.pair$agtcga/cps.agt.cga)
cps.obs.exp.agtcgc<-log(dn31.codon.pair$agtcgc/cps.agt.cgc)
cps.obs.exp.agtcgg<-log(dn31.codon.pair$agtcgg/cps.agt.cgg)
cps.obs.exp.agtcgt<-log(dn31.codon.pair$agtcgt/cps.agt.cgt)
cps.obs.exp.agtcta<-log(dn31.codon.pair$agtcta/cps.agt.cta)
cps.obs.exp.agtctc<-log(dn31.codon.pair$agtctc/cps.agt.ctc)
cps.obs.exp.agtctg<-log(dn31.codon.pair$agtctg/cps.agt.ctg)
cps.obs.exp.agtctt<-log(dn31.codon.pair$agtctt/cps.agt.ctt)

cps.obs.exp.agtgaa<-log(dn31.codon.pair$agtgaa/cps.agt.gaa)
cps.obs.exp.agtgac<-log(dn31.codon.pair$agtgac/cps.agt.gac)
cps.obs.exp.agtgag<-log(dn31.codon.pair$agtgag/cps.agt.gag)
cps.obs.exp.agtgat<-log(dn31.codon.pair$agtgat/cps.agt.gat)
cps.obs.exp.agtgca<-log(dn31.codon.pair$agtgca/cps.agt.gca)
cps.obs.exp.agtgcc<-log(dn31.codon.pair$agtgcc/cps.agt.gcc)
cps.obs.exp.agtgcg<-log(dn31.codon.pair$agtgcg/cps.agt.gcg)
cps.obs.exp.agtgct<-log(dn31.codon.pair$agtgct/cps.agt.gct)
cps.obs.exp.agtgga<-log(dn31.codon.pair$agtgga/cps.agt.gga)
cps.obs.exp.agtggc<-log(dn31.codon.pair$agtggc/cps.agt.ggc)
cps.obs.exp.agtggg<-log(dn31.codon.pair$agtggg/cps.agt.ggg)
cps.obs.exp.agtggt<-log(dn31.codon.pair$agtggt/cps.agt.ggt)
cps.obs.exp.agtgta<-log(dn31.codon.pair$agtgta/cps.agt.gta)
cps.obs.exp.agtgtc<-log(dn31.codon.pair$agtgtc/cps.agt.gtc)
cps.obs.exp.agtgtg<-log(dn31.codon.pair$agtgtg/cps.agt.gtg)
cps.obs.exp.agtgtt<-log(dn31.codon.pair$agtgtt/cps.agt.gtt)

#cps.obs.exp.agttaa<-log(dn31.codon.pair$agttaa/cps.agt.taa)
cps.obs.exp.agttac<-log(dn31.codon.pair$agttac/cps.agt.tac)
#cps.obs.exp.agttag<-log(dn31.codon.pair$agttag/cps.agt.tag)
cps.obs.exp.agttat<-log(dn31.codon.pair$agttat/cps.agt.tat)
cps.obs.exp.agttca<-log(dn31.codon.pair$agttca/cps.agt.tca)
cps.obs.exp.agttcc<-log(dn31.codon.pair$agttcc/cps.agt.tcc)
cps.obs.exp.agttcg<-log(dn31.codon.pair$agttcg/cps.agt.tcg)
cps.obs.exp.agttct<-log(dn31.codon.pair$agttct/cps.agt.tct)
#cps.obs.exp.agttga<-log(dn31.codon.pair$agttga/cps.agt.tga)
cps.obs.exp.agttgc<-log(dn31.codon.pair$agttgc/cps.agt.tgc)
cps.obs.exp.agttgg<-log(dn31.codon.pair$agttgg/cps.agt.tgg)
cps.obs.exp.agttgt<-log(dn31.codon.pair$agttgt/cps.agt.tgt)
cps.obs.exp.agttta<-log(dn31.codon.pair$agttta/cps.agt.tta)
cps.obs.exp.agtttc<-log(dn31.codon.pair$agtttc/cps.agt.ttc)
cps.obs.exp.agtttg<-log(dn31.codon.pair$agtttg/cps.agt.ttg)
cps.obs.exp.agtttt<-log(dn31.codon.pair$agtttt/cps.agt.ttt)




















cps.obs.exp.ataaaa<-log(dn31.codon.pair$ataaaa/cps.ata.aaa)
cps.obs.exp.ataaac<-log(dn31.codon.pair$ataaac/cps.ata.aac)
cps.obs.exp.ataaag<-log(dn31.codon.pair$ataaag/cps.ata.aag)
cps.obs.exp.ataaat<-log(dn31.codon.pair$ataaat/cps.ata.aat)
cps.obs.exp.ataaca<-log(dn31.codon.pair$ataaca/cps.ata.aca)
cps.obs.exp.ataacc<-log(dn31.codon.pair$ataacc/cps.ata.acc)
cps.obs.exp.ataacg<-log(dn31.codon.pair$ataacg/cps.ata.acg)
cps.obs.exp.ataact<-log(dn31.codon.pair$ataact/cps.ata.act)
cps.obs.exp.ataaga<-log(dn31.codon.pair$ataaga/cps.ata.aga)
cps.obs.exp.ataagc<-log(dn31.codon.pair$ataagc/cps.ata.agc)
cps.obs.exp.ataagg<-log(dn31.codon.pair$ataagg/cps.ata.agg)
cps.obs.exp.ataagt<-log(dn31.codon.pair$ataagt/cps.ata.agt)
cps.obs.exp.ataata<-log(dn31.codon.pair$ataata/cps.ata.ata)
cps.obs.exp.ataatc<-log(dn31.codon.pair$ataatc/cps.ata.atc)
cps.obs.exp.ataatg<-log(dn31.codon.pair$ataatg/cps.ata.atg)
cps.obs.exp.ataatt<-log(dn31.codon.pair$ataatt/cps.ata.att)

cps.obs.exp.atacaa<-log(dn31.codon.pair$atacaa/cps.ata.caa)
cps.obs.exp.atacac<-log(dn31.codon.pair$atacac/cps.ata.cac)
cps.obs.exp.atacag<-log(dn31.codon.pair$atacag/cps.ata.cag)
cps.obs.exp.atacat<-log(dn31.codon.pair$atacat/cps.ata.cat)
cps.obs.exp.atacca<-log(dn31.codon.pair$atacca/cps.ata.cca)
cps.obs.exp.ataccc<-log(dn31.codon.pair$ataccc/cps.ata.ccc)
cps.obs.exp.ataccg<-log(dn31.codon.pair$ataccg/cps.ata.ccg)
cps.obs.exp.atacct<-log(dn31.codon.pair$atacct/cps.ata.cct)
cps.obs.exp.atacga<-log(dn31.codon.pair$atacga/cps.ata.cga)
cps.obs.exp.atacgc<-log(dn31.codon.pair$atacgc/cps.ata.cgc)
cps.obs.exp.atacgg<-log(dn31.codon.pair$atacgg/cps.ata.cgg)
cps.obs.exp.atacgt<-log(dn31.codon.pair$atacgt/cps.ata.cgt)
cps.obs.exp.atacta<-log(dn31.codon.pair$atacta/cps.ata.cta)
cps.obs.exp.atactc<-log(dn31.codon.pair$atactc/cps.ata.ctc)
cps.obs.exp.atactg<-log(dn31.codon.pair$atactg/cps.ata.ctg)
cps.obs.exp.atactt<-log(dn31.codon.pair$atactt/cps.ata.ctt)

cps.obs.exp.atagaa<-log(dn31.codon.pair$atagaa/cps.ata.gaa)
cps.obs.exp.atagac<-log(dn31.codon.pair$atagac/cps.ata.gac)
cps.obs.exp.atagag<-log(dn31.codon.pair$atagag/cps.ata.gag)
cps.obs.exp.atagat<-log(dn31.codon.pair$atagat/cps.ata.gat)
cps.obs.exp.atagca<-log(dn31.codon.pair$atagca/cps.ata.gca)
cps.obs.exp.atagcc<-log(dn31.codon.pair$atagcc/cps.ata.gcc)
cps.obs.exp.atagcg<-log(dn31.codon.pair$atagcg/cps.ata.gcg)
cps.obs.exp.atagct<-log(dn31.codon.pair$atagct/cps.ata.gct)
cps.obs.exp.atagga<-log(dn31.codon.pair$atagga/cps.ata.gga)
cps.obs.exp.ataggc<-log(dn31.codon.pair$ataggc/cps.ata.ggc)
cps.obs.exp.ataggg<-log(dn31.codon.pair$ataggg/cps.ata.ggg)
cps.obs.exp.ataggt<-log(dn31.codon.pair$ataggt/cps.ata.ggt)
cps.obs.exp.atagta<-log(dn31.codon.pair$atagta/cps.ata.gta)
cps.obs.exp.atagtc<-log(dn31.codon.pair$atagtc/cps.ata.gtc)
cps.obs.exp.atagtg<-log(dn31.codon.pair$atagtg/cps.ata.gtg)
cps.obs.exp.atagtt<-log(dn31.codon.pair$atagtt/cps.ata.gtt)

#cps.obs.exp.atataa<-log(dn31.codon.pair$atataa/cps.ata.taa)
cps.obs.exp.atatac<-log(dn31.codon.pair$atatac/cps.ata.tac)
#cps.obs.exp.atatag<-log(dn31.codon.pair$atatag/cps.ata.tag)
cps.obs.exp.atatat<-log(dn31.codon.pair$atatat/cps.ata.tat)
cps.obs.exp.atatca<-log(dn31.codon.pair$atatca/cps.ata.tca)
cps.obs.exp.atatcc<-log(dn31.codon.pair$atatcc/cps.ata.tcc)
cps.obs.exp.atatcg<-log(dn31.codon.pair$atatcg/cps.ata.tcg)
cps.obs.exp.atatct<-log(dn31.codon.pair$atatct/cps.ata.tct)
#cps.obs.exp.atatga<-log(dn31.codon.pair$atatga/cps.ata.tga)
cps.obs.exp.atatgc<-log(dn31.codon.pair$atatgc/cps.ata.tgc)
cps.obs.exp.atatgg<-log(dn31.codon.pair$atatgg/cps.ata.tgg)
cps.obs.exp.atatgt<-log(dn31.codon.pair$atatgt/cps.ata.tgt)
cps.obs.exp.atatta<-log(dn31.codon.pair$atatta/cps.ata.tta)
cps.obs.exp.atattc<-log(dn31.codon.pair$atattc/cps.ata.ttc)
cps.obs.exp.atattg<-log(dn31.codon.pair$atattg/cps.ata.ttg)
cps.obs.exp.atattt<-log(dn31.codon.pair$atattt/cps.ata.ttt)







cps.obs.exp.atcaaa<-log(dn31.codon.pair$atcaaa/cps.atc.aaa)
cps.obs.exp.atcaac<-log(dn31.codon.pair$atcaac/cps.atc.aac)
cps.obs.exp.atcaag<-log(dn31.codon.pair$atcaag/cps.atc.aag)
cps.obs.exp.atcaat<-log(dn31.codon.pair$atcaat/cps.atc.aat)
cps.obs.exp.atcaca<-log(dn31.codon.pair$atcaca/cps.atc.aca)
cps.obs.exp.atcacc<-log(dn31.codon.pair$atcacc/cps.atc.acc)
cps.obs.exp.atcacg<-log(dn31.codon.pair$atcacg/cps.atc.acg)
cps.obs.exp.atcact<-log(dn31.codon.pair$atcact/cps.atc.act)
cps.obs.exp.atcaga<-log(dn31.codon.pair$atcaga/cps.atc.aga)
cps.obs.exp.atcagc<-log(dn31.codon.pair$atcagc/cps.atc.agc)
cps.obs.exp.atcagg<-log(dn31.codon.pair$atcagg/cps.atc.agg)
cps.obs.exp.atcagt<-log(dn31.codon.pair$atcagt/cps.atc.agt)
cps.obs.exp.atcata<-log(dn31.codon.pair$atcata/cps.atc.ata)
cps.obs.exp.atcatc<-log(dn31.codon.pair$atcatc/cps.atc.atc)
cps.obs.exp.atcatg<-log(dn31.codon.pair$atcatg/cps.atc.atg)
cps.obs.exp.atcatt<-log(dn31.codon.pair$atcatt/cps.atc.att)

cps.obs.exp.atccaa<-log(dn31.codon.pair$atccaa/cps.atc.caa)
cps.obs.exp.atccac<-log(dn31.codon.pair$atccac/cps.atc.cac)
cps.obs.exp.atccag<-log(dn31.codon.pair$atccag/cps.atc.cag)
cps.obs.exp.atccat<-log(dn31.codon.pair$atccat/cps.atc.cat)
cps.obs.exp.atccca<-log(dn31.codon.pair$atccca/cps.atc.cca)
cps.obs.exp.atcccc<-log(dn31.codon.pair$atcccc/cps.atc.ccc)
cps.obs.exp.atcccg<-log(dn31.codon.pair$atcccg/cps.atc.ccg)
cps.obs.exp.atccct<-log(dn31.codon.pair$atccct/cps.atc.cct)
cps.obs.exp.atccga<-log(dn31.codon.pair$atccga/cps.atc.cga)
cps.obs.exp.atccgc<-log(dn31.codon.pair$atccgc/cps.atc.cgc)
cps.obs.exp.atccgg<-log(dn31.codon.pair$atccgg/cps.atc.cgg)
cps.obs.exp.atccgt<-log(dn31.codon.pair$atccgt/cps.atc.cgt)
cps.obs.exp.atccta<-log(dn31.codon.pair$atccta/cps.atc.cta)
cps.obs.exp.atcctc<-log(dn31.codon.pair$atcctc/cps.atc.ctc)
cps.obs.exp.atcctg<-log(dn31.codon.pair$atcctg/cps.atc.ctg)
cps.obs.exp.atcctt<-log(dn31.codon.pair$atcctt/cps.atc.ctt)

cps.obs.exp.atcgaa<-log(dn31.codon.pair$atcgaa/cps.atc.gaa)
cps.obs.exp.atcgac<-log(dn31.codon.pair$atcgac/cps.atc.gac)
cps.obs.exp.atcgag<-log(dn31.codon.pair$atcgag/cps.atc.gag)
cps.obs.exp.atcgat<-log(dn31.codon.pair$atcgat/cps.atc.gat)
cps.obs.exp.atcgca<-log(dn31.codon.pair$atcgca/cps.atc.gca)
cps.obs.exp.atcgcc<-log(dn31.codon.pair$atcgcc/cps.atc.gcc)
cps.obs.exp.atcgcg<-log(dn31.codon.pair$atcgcg/cps.atc.gcg)
cps.obs.exp.atcgct<-log(dn31.codon.pair$atcgct/cps.atc.gct)
cps.obs.exp.atcgga<-log(dn31.codon.pair$atcgga/cps.atc.gga)
cps.obs.exp.atcggc<-log(dn31.codon.pair$atcggc/cps.atc.ggc)
cps.obs.exp.atcggg<-log(dn31.codon.pair$atcggg/cps.atc.ggg)
cps.obs.exp.atcggt<-log(dn31.codon.pair$atcggt/cps.atc.ggt)
cps.obs.exp.atcgta<-log(dn31.codon.pair$atcgta/cps.atc.gta)
cps.obs.exp.atcgtc<-log(dn31.codon.pair$atcgtc/cps.atc.gtc)
cps.obs.exp.atcgtg<-log(dn31.codon.pair$atcgtg/cps.atc.gtg)
cps.obs.exp.atcgtt<-log(dn31.codon.pair$atcgtt/cps.atc.gtt)

#cps.obs.exp.atctaa<-log(dn31.codon.pair$atctaa/cps.atc.taa)
cps.obs.exp.atctac<-log(dn31.codon.pair$atctac/cps.atc.tac)
#cps.obs.exp.atctag<-log(dn31.codon.pair$atctag/cps.atc.tag)
cps.obs.exp.atctat<-log(dn31.codon.pair$atctat/cps.atc.tat)
cps.obs.exp.atctca<-log(dn31.codon.pair$atctca/cps.atc.tca)
cps.obs.exp.atctcc<-log(dn31.codon.pair$atctcc/cps.atc.tcc)
cps.obs.exp.atctcg<-log(dn31.codon.pair$atctcg/cps.atc.tcg)
cps.obs.exp.atctct<-log(dn31.codon.pair$atctct/cps.atc.tct)
#cps.obs.exp.atctga<-log(dn31.codon.pair$atctga/cps.atc.tga)
cps.obs.exp.atctgc<-log(dn31.codon.pair$atctgc/cps.atc.tgc)
cps.obs.exp.atctgg<-log(dn31.codon.pair$atctgg/cps.atc.tgg)
cps.obs.exp.atctgt<-log(dn31.codon.pair$atctgt/cps.atc.tgt)
cps.obs.exp.atctta<-log(dn31.codon.pair$atctta/cps.atc.tta)
cps.obs.exp.atcttc<-log(dn31.codon.pair$atcttc/cps.atc.ttc)
cps.obs.exp.atcttg<-log(dn31.codon.pair$atcttg/cps.atc.ttg)
cps.obs.exp.atcttt<-log(dn31.codon.pair$atcttt/cps.atc.ttt)








cps.obs.exp.atgaaa<-log(dn31.codon.pair$atgaaa/cps.atg.aaa)
cps.obs.exp.atgaac<-log(dn31.codon.pair$atgaac/cps.atg.aac)
cps.obs.exp.atgaag<-log(dn31.codon.pair$atgaag/cps.atg.aag)
cps.obs.exp.atgaat<-log(dn31.codon.pair$atgaat/cps.atg.aat)
cps.obs.exp.atgaca<-log(dn31.codon.pair$atgaca/cps.atg.aca)
cps.obs.exp.atgacc<-log(dn31.codon.pair$atgacc/cps.atg.acc)
cps.obs.exp.atgacg<-log(dn31.codon.pair$atgacg/cps.atg.acg)
cps.obs.exp.atgact<-log(dn31.codon.pair$atgact/cps.atg.act)
cps.obs.exp.atgaga<-log(dn31.codon.pair$atgaga/cps.atg.aga)
cps.obs.exp.atgagc<-log(dn31.codon.pair$atgagc/cps.atg.agc)
cps.obs.exp.atgagg<-log(dn31.codon.pair$atgagg/cps.atg.agg)
cps.obs.exp.atgagt<-log(dn31.codon.pair$atgagt/cps.atg.agt)
cps.obs.exp.atgata<-log(dn31.codon.pair$atgata/cps.atg.ata)
cps.obs.exp.atgatc<-log(dn31.codon.pair$atgatc/cps.atg.atc)
cps.obs.exp.atgatg<-log(dn31.codon.pair$atgatg/cps.atg.atg)
cps.obs.exp.atgatt<-log(dn31.codon.pair$atgatt/cps.atg.att)

cps.obs.exp.atgcaa<-log(dn31.codon.pair$atgcaa/cps.atg.caa)
cps.obs.exp.atgcac<-log(dn31.codon.pair$atgcac/cps.atg.cac)
cps.obs.exp.atgcag<-log(dn31.codon.pair$atgcag/cps.atg.cag)
cps.obs.exp.atgcat<-log(dn31.codon.pair$atgcat/cps.atg.cat)
cps.obs.exp.atgcca<-log(dn31.codon.pair$atgcca/cps.atg.cca)
cps.obs.exp.atgccc<-log(dn31.codon.pair$atgccc/cps.atg.ccc)
cps.obs.exp.atgccg<-log(dn31.codon.pair$atgccg/cps.atg.ccg)
cps.obs.exp.atgcct<-log(dn31.codon.pair$atgcct/cps.atg.cct)
cps.obs.exp.atgcga<-log(dn31.codon.pair$atgcga/cps.atg.cga)
cps.obs.exp.atgcgc<-log(dn31.codon.pair$atgcgc/cps.atg.cgc)
cps.obs.exp.atgcgg<-log(dn31.codon.pair$atgcgg/cps.atg.cgg)
cps.obs.exp.atgcgt<-log(dn31.codon.pair$atgcgt/cps.atg.cgt)
cps.obs.exp.atgcta<-log(dn31.codon.pair$atgcta/cps.atg.cta)
cps.obs.exp.atgctc<-log(dn31.codon.pair$atgctc/cps.atg.ctc)
cps.obs.exp.atgctg<-log(dn31.codon.pair$atgctg/cps.atg.ctg)
cps.obs.exp.atgctt<-log(dn31.codon.pair$atgctt/cps.atg.ctt)

cps.obs.exp.atggaa<-log(dn31.codon.pair$atggaa/cps.atg.gaa)
cps.obs.exp.atggac<-log(dn31.codon.pair$atggac/cps.atg.gac)
cps.obs.exp.atggag<-log(dn31.codon.pair$atggag/cps.atg.gag)
cps.obs.exp.atggat<-log(dn31.codon.pair$atggat/cps.atg.gat)
cps.obs.exp.atggca<-log(dn31.codon.pair$atggca/cps.atg.gca)
cps.obs.exp.atggcc<-log(dn31.codon.pair$atggcc/cps.atg.gcc)
cps.obs.exp.atggcg<-log(dn31.codon.pair$atggcg/cps.atg.gcg)
cps.obs.exp.atggct<-log(dn31.codon.pair$atggct/cps.atg.gct)
cps.obs.exp.atggga<-log(dn31.codon.pair$atggga/cps.atg.gga)
cps.obs.exp.atgggc<-log(dn31.codon.pair$atgggc/cps.atg.ggc)
cps.obs.exp.atgggg<-log(dn31.codon.pair$atgggg/cps.atg.ggg)
cps.obs.exp.atgggt<-log(dn31.codon.pair$atgggt/cps.atg.ggt)
cps.obs.exp.atggta<-log(dn31.codon.pair$atggta/cps.atg.gta)
cps.obs.exp.atggtc<-log(dn31.codon.pair$atggtc/cps.atg.gtc)
cps.obs.exp.atggtg<-log(dn31.codon.pair$atggtg/cps.atg.gtg)
cps.obs.exp.atggtt<-log(dn31.codon.pair$atggtt/cps.atg.gtt)

#cps.obs.exp.atgtaa<-log(dn31.codon.pair$atgtaa/cps.atg.taa)
cps.obs.exp.atgtac<-log(dn31.codon.pair$atgtac/cps.atg.tac)
#cps.obs.exp.atgtag<-log(dn31.codon.pair$atgtag/cps.atg.tag)
cps.obs.exp.atgtat<-log(dn31.codon.pair$atgtat/cps.atg.tat)
cps.obs.exp.atgtca<-log(dn31.codon.pair$atgtca/cps.atg.tca)
cps.obs.exp.atgtcc<-log(dn31.codon.pair$atgtcc/cps.atg.tcc)
cps.obs.exp.atgtcg<-log(dn31.codon.pair$atgtcg/cps.atg.tcg)
cps.obs.exp.atgtct<-log(dn31.codon.pair$atgtct/cps.atg.tct)
#cps.obs.exp.atgtga<-log(dn31.codon.pair$atgtga/cps.atg.tga)
cps.obs.exp.atgtgc<-log(dn31.codon.pair$atgtgc/cps.atg.tgc)
cps.obs.exp.atgtgg<-log(dn31.codon.pair$atgtgg/cps.atg.tgg)
cps.obs.exp.atgtgt<-log(dn31.codon.pair$atgtgt/cps.atg.tgt)
cps.obs.exp.atgtta<-log(dn31.codon.pair$atgtta/cps.atg.tta)
cps.obs.exp.atgttc<-log(dn31.codon.pair$atgttc/cps.atg.ttc)
cps.obs.exp.atgttg<-log(dn31.codon.pair$atgttg/cps.atg.ttg)
cps.obs.exp.atgttt<-log(dn31.codon.pair$atgttt/cps.atg.ttt)







cps.obs.exp.attaaa<-log(dn31.codon.pair$attaaa/cps.att.aaa)
cps.obs.exp.attaac<-log(dn31.codon.pair$attaac/cps.att.aac)
cps.obs.exp.attaag<-log(dn31.codon.pair$attaag/cps.att.aag)
cps.obs.exp.attaat<-log(dn31.codon.pair$attaat/cps.att.aat)
cps.obs.exp.attaca<-log(dn31.codon.pair$attaca/cps.att.aca)
cps.obs.exp.attacc<-log(dn31.codon.pair$attacc/cps.att.acc)
cps.obs.exp.attacg<-log(dn31.codon.pair$attacg/cps.att.acg)
cps.obs.exp.attact<-log(dn31.codon.pair$attact/cps.att.act)
cps.obs.exp.attaga<-log(dn31.codon.pair$attaga/cps.att.aga)
cps.obs.exp.attagc<-log(dn31.codon.pair$attagc/cps.att.agc)
cps.obs.exp.attagg<-log(dn31.codon.pair$attagg/cps.att.agg)
cps.obs.exp.attagt<-log(dn31.codon.pair$attagt/cps.att.agt)
cps.obs.exp.attata<-log(dn31.codon.pair$attata/cps.att.ata)
cps.obs.exp.attatc<-log(dn31.codon.pair$attatc/cps.att.atc)
cps.obs.exp.attatg<-log(dn31.codon.pair$attatg/cps.att.atg)
cps.obs.exp.attatt<-log(dn31.codon.pair$attatt/cps.att.att)

cps.obs.exp.attcaa<-log(dn31.codon.pair$attcaa/cps.att.caa)
cps.obs.exp.attcac<-log(dn31.codon.pair$attcac/cps.att.cac)
cps.obs.exp.attcag<-log(dn31.codon.pair$attcag/cps.att.cag)
cps.obs.exp.attcat<-log(dn31.codon.pair$attcat/cps.att.cat)
cps.obs.exp.attcca<-log(dn31.codon.pair$attcca/cps.att.cca)
cps.obs.exp.attccc<-log(dn31.codon.pair$attccc/cps.att.ccc)
cps.obs.exp.attccg<-log(dn31.codon.pair$attccg/cps.att.ccg)
cps.obs.exp.attcct<-log(dn31.codon.pair$attcct/cps.att.cct)
cps.obs.exp.attcga<-log(dn31.codon.pair$attcga/cps.att.cga)
cps.obs.exp.attcgc<-log(dn31.codon.pair$attcgc/cps.att.cgc)
cps.obs.exp.attcgg<-log(dn31.codon.pair$attcgg/cps.att.cgg)
cps.obs.exp.attcgt<-log(dn31.codon.pair$attcgt/cps.att.cgt)
cps.obs.exp.attcta<-log(dn31.codon.pair$attcta/cps.att.cta)
cps.obs.exp.attctc<-log(dn31.codon.pair$attctc/cps.att.ctc)
cps.obs.exp.attctg<-log(dn31.codon.pair$attctg/cps.att.ctg)
cps.obs.exp.attctt<-log(dn31.codon.pair$attctt/cps.att.ctt)

cps.obs.exp.attgaa<-log(dn31.codon.pair$attgaa/cps.att.gaa)
cps.obs.exp.attgac<-log(dn31.codon.pair$attgac/cps.att.gac)
cps.obs.exp.attgag<-log(dn31.codon.pair$attgag/cps.att.gag)
cps.obs.exp.attgat<-log(dn31.codon.pair$attgat/cps.att.gat)
cps.obs.exp.attgca<-log(dn31.codon.pair$attgca/cps.att.gca)
cps.obs.exp.attgcc<-log(dn31.codon.pair$attgcc/cps.att.gcc)
cps.obs.exp.attgcg<-log(dn31.codon.pair$attgcg/cps.att.gcg)
cps.obs.exp.attgct<-log(dn31.codon.pair$attgct/cps.att.gct)
cps.obs.exp.attgga<-log(dn31.codon.pair$attgga/cps.att.gga)
cps.obs.exp.attggc<-log(dn31.codon.pair$attggc/cps.att.ggc)
cps.obs.exp.attggg<-log(dn31.codon.pair$attggg/cps.att.ggg)
cps.obs.exp.attggt<-log(dn31.codon.pair$attggt/cps.att.ggt)
cps.obs.exp.attgta<-log(dn31.codon.pair$attgta/cps.att.gta)
cps.obs.exp.attgtc<-log(dn31.codon.pair$attgtc/cps.att.gtc)
cps.obs.exp.attgtg<-log(dn31.codon.pair$attgtg/cps.att.gtg)
cps.obs.exp.attgtt<-log(dn31.codon.pair$attgtt/cps.att.gtt)

#cps.obs.exp.atttaa<-log(dn31.codon.pair$atttaa/cps.att.taa)
cps.obs.exp.atttac<-log(dn31.codon.pair$atttac/cps.att.tac)
#cps.obs.exp.atttag<-log(dn31.codon.pair$atttag/cps.att.tag)
cps.obs.exp.atttat<-log(dn31.codon.pair$atttat/cps.att.tat)
cps.obs.exp.atttca<-log(dn31.codon.pair$atttca/cps.att.tca)
cps.obs.exp.atttcc<-log(dn31.codon.pair$atttcc/cps.att.tcc)
cps.obs.exp.atttcg<-log(dn31.codon.pair$atttcg/cps.att.tcg)
cps.obs.exp.atttct<-log(dn31.codon.pair$atttct/cps.att.tct)
#cps.obs.exp.atttga<-log(dn31.codon.pair$atttga/cps.att.tga)
cps.obs.exp.atttgc<-log(dn31.codon.pair$atttgc/cps.att.tgc)
cps.obs.exp.atttgg<-log(dn31.codon.pair$atttgg/cps.att.tgg)
cps.obs.exp.atttgt<-log(dn31.codon.pair$atttgt/cps.att.tgt)
cps.obs.exp.atttta<-log(dn31.codon.pair$atttta/cps.att.tta)
cps.obs.exp.attttc<-log(dn31.codon.pair$attttc/cps.att.ttc)
cps.obs.exp.attttg<-log(dn31.codon.pair$attttg/cps.att.ttg)
cps.obs.exp.attttt<-log(dn31.codon.pair$attttt/cps.att.ttt)



















cps.obs.exp.caaaaa<-log(dn31.codon.pair$caaaaa/cps.caa.aaa)
cps.obs.exp.caaaac<-log(dn31.codon.pair$caaaac/cps.caa.aac)
cps.obs.exp.caaaag<-log(dn31.codon.pair$caaaag/cps.caa.aag)
cps.obs.exp.caaaat<-log(dn31.codon.pair$caaaat/cps.caa.aat)
cps.obs.exp.caaaca<-log(dn31.codon.pair$caaaca/cps.caa.aca)
cps.obs.exp.caaacc<-log(dn31.codon.pair$caaacc/cps.caa.acc)
cps.obs.exp.caaacg<-log(dn31.codon.pair$caaacg/cps.caa.acg)
cps.obs.exp.caaact<-log(dn31.codon.pair$caaact/cps.caa.act)
cps.obs.exp.caaaga<-log(dn31.codon.pair$caaaga/cps.caa.aga)
cps.obs.exp.caaagc<-log(dn31.codon.pair$caaagc/cps.caa.agc)
cps.obs.exp.caaagg<-log(dn31.codon.pair$caaagg/cps.caa.agg)
cps.obs.exp.caaagt<-log(dn31.codon.pair$caaagt/cps.caa.agt)
cps.obs.exp.caaata<-log(dn31.codon.pair$caaata/cps.caa.ata)
cps.obs.exp.caaatc<-log(dn31.codon.pair$caaatc/cps.caa.atc)
cps.obs.exp.caaatg<-log(dn31.codon.pair$caaatg/cps.caa.atg)
cps.obs.exp.caaatt<-log(dn31.codon.pair$caaatt/cps.caa.att)

cps.obs.exp.caacaa<-log(dn31.codon.pair$caacaa/cps.caa.caa)
cps.obs.exp.caacac<-log(dn31.codon.pair$caacac/cps.caa.cac)
cps.obs.exp.caacag<-log(dn31.codon.pair$caacag/cps.caa.cag)
cps.obs.exp.caacat<-log(dn31.codon.pair$caacat/cps.caa.cat)
cps.obs.exp.caacca<-log(dn31.codon.pair$caacca/cps.caa.cca)
cps.obs.exp.caaccc<-log(dn31.codon.pair$caaccc/cps.caa.ccc)
cps.obs.exp.caaccg<-log(dn31.codon.pair$caaccg/cps.caa.ccg)
cps.obs.exp.caacct<-log(dn31.codon.pair$caacct/cps.caa.cct)
cps.obs.exp.caacga<-log(dn31.codon.pair$caacga/cps.caa.cga)
cps.obs.exp.caacgc<-log(dn31.codon.pair$caacgc/cps.caa.cgc)
cps.obs.exp.caacgg<-log(dn31.codon.pair$caacgg/cps.caa.cgg)
cps.obs.exp.caacgt<-log(dn31.codon.pair$caacgt/cps.caa.cgt)
cps.obs.exp.caacta<-log(dn31.codon.pair$caacta/cps.caa.cta)
cps.obs.exp.caactc<-log(dn31.codon.pair$caactc/cps.caa.ctc)
cps.obs.exp.caactg<-log(dn31.codon.pair$caactg/cps.caa.ctg)
cps.obs.exp.caactt<-log(dn31.codon.pair$caactt/cps.caa.ctt)

cps.obs.exp.caagaa<-log(dn31.codon.pair$caagaa/cps.caa.gaa)
cps.obs.exp.caagac<-log(dn31.codon.pair$caagac/cps.caa.gac)
cps.obs.exp.caagag<-log(dn31.codon.pair$caagag/cps.caa.gag)
cps.obs.exp.caagat<-log(dn31.codon.pair$caagat/cps.caa.gat)
cps.obs.exp.caagca<-log(dn31.codon.pair$caagca/cps.caa.gca)
cps.obs.exp.caagcc<-log(dn31.codon.pair$caagcc/cps.caa.gcc)
cps.obs.exp.caagcg<-log(dn31.codon.pair$caagcg/cps.caa.gcg)
cps.obs.exp.caagct<-log(dn31.codon.pair$caagct/cps.caa.gct)
cps.obs.exp.caagga<-log(dn31.codon.pair$caagga/cps.caa.gga)
cps.obs.exp.caaggc<-log(dn31.codon.pair$caaggc/cps.caa.ggc)
cps.obs.exp.caaggg<-log(dn31.codon.pair$caaggg/cps.caa.ggg)
cps.obs.exp.caaggt<-log(dn31.codon.pair$caaggt/cps.caa.ggt)
cps.obs.exp.caagta<-log(dn31.codon.pair$caagta/cps.caa.gta)
cps.obs.exp.caagtc<-log(dn31.codon.pair$caagtc/cps.caa.gtc)
cps.obs.exp.caagtg<-log(dn31.codon.pair$caagtg/cps.caa.gtg)
cps.obs.exp.caagtt<-log(dn31.codon.pair$caagtt/cps.caa.gtt)

#cps.obs.exp.caataa<-log(dn31.codon.pair$caataa/cps.caa.taa)
cps.obs.exp.caatac<-log(dn31.codon.pair$caatac/cps.caa.tac)
#cps.obs.exp.caatag<-log(dn31.codon.pair$caatag/cps.caa.tag)
cps.obs.exp.caatat<-log(dn31.codon.pair$caatat/cps.caa.tat)
cps.obs.exp.caatca<-log(dn31.codon.pair$caatca/cps.caa.tca)
cps.obs.exp.caatcc<-log(dn31.codon.pair$caatcc/cps.caa.tcc)
cps.obs.exp.caatcg<-log(dn31.codon.pair$caatcg/cps.caa.tcg)
cps.obs.exp.caatct<-log(dn31.codon.pair$caatct/cps.caa.tct)
#cps.obs.exp.caatga<-log(dn31.codon.pair$caatga/cps.caa.tga)
cps.obs.exp.caatgc<-log(dn31.codon.pair$caatgc/cps.caa.tgc)
cps.obs.exp.caatgg<-log(dn31.codon.pair$caatgg/cps.caa.tgg)
cps.obs.exp.caatgt<-log(dn31.codon.pair$caatgt/cps.caa.tgt)
cps.obs.exp.caatta<-log(dn31.codon.pair$caatta/cps.caa.tta)
cps.obs.exp.caattc<-log(dn31.codon.pair$caattc/cps.caa.ttc)
cps.obs.exp.caattg<-log(dn31.codon.pair$caattg/cps.caa.ttg)
cps.obs.exp.caattt<-log(dn31.codon.pair$caattt/cps.caa.ttt)







cps.obs.exp.cacaaa<-log(dn31.codon.pair$cacaaa/cps.cac.aaa)
cps.obs.exp.cacaac<-log(dn31.codon.pair$cacaac/cps.cac.aac)
cps.obs.exp.cacaag<-log(dn31.codon.pair$cacaag/cps.cac.aag)
cps.obs.exp.cacaat<-log(dn31.codon.pair$cacaat/cps.cac.aat)
cps.obs.exp.cacaca<-log(dn31.codon.pair$cacaca/cps.cac.aca)
cps.obs.exp.cacacc<-log(dn31.codon.pair$cacacc/cps.cac.acc)
cps.obs.exp.cacacg<-log(dn31.codon.pair$cacacg/cps.cac.acg)
cps.obs.exp.cacact<-log(dn31.codon.pair$cacact/cps.cac.act)
cps.obs.exp.cacaga<-log(dn31.codon.pair$cacaga/cps.cac.aga)
cps.obs.exp.cacagc<-log(dn31.codon.pair$cacagc/cps.cac.agc)
cps.obs.exp.cacagg<-log(dn31.codon.pair$cacagg/cps.cac.agg)
cps.obs.exp.cacagt<-log(dn31.codon.pair$cacagt/cps.cac.agt)
cps.obs.exp.cacata<-log(dn31.codon.pair$cacata/cps.cac.ata)
cps.obs.exp.cacatc<-log(dn31.codon.pair$cacatc/cps.cac.atc)
cps.obs.exp.cacatg<-log(dn31.codon.pair$cacatg/cps.cac.atg)
cps.obs.exp.cacatt<-log(dn31.codon.pair$cacatt/cps.cac.att)

cps.obs.exp.caccaa<-log(dn31.codon.pair$caccaa/cps.cac.caa)
cps.obs.exp.caccac<-log(dn31.codon.pair$caccac/cps.cac.cac)
cps.obs.exp.caccag<-log(dn31.codon.pair$caccag/cps.cac.cag)
cps.obs.exp.caccat<-log(dn31.codon.pair$caccat/cps.cac.cat)
cps.obs.exp.caccca<-log(dn31.codon.pair$caccca/cps.cac.cca)
cps.obs.exp.cacccc<-log(dn31.codon.pair$cacccc/cps.cac.ccc)
cps.obs.exp.cacccg<-log(dn31.codon.pair$cacccg/cps.cac.ccg)
cps.obs.exp.caccct<-log(dn31.codon.pair$caccct/cps.cac.cct)
cps.obs.exp.caccga<-log(dn31.codon.pair$caccga/cps.cac.cga)
cps.obs.exp.caccgc<-log(dn31.codon.pair$caccgc/cps.cac.cgc)
cps.obs.exp.caccgg<-log(dn31.codon.pair$caccgg/cps.cac.cgg)
cps.obs.exp.caccgt<-log(dn31.codon.pair$caccgt/cps.cac.cgt)
cps.obs.exp.caccta<-log(dn31.codon.pair$caccta/cps.cac.cta)
cps.obs.exp.cacctc<-log(dn31.codon.pair$cacctc/cps.cac.ctc)
cps.obs.exp.cacctg<-log(dn31.codon.pair$cacctg/cps.cac.ctg)
cps.obs.exp.cacctt<-log(dn31.codon.pair$cacctt/cps.cac.ctt)

cps.obs.exp.cacgaa<-log(dn31.codon.pair$cacgaa/cps.cac.gaa)
cps.obs.exp.cacgac<-log(dn31.codon.pair$cacgac/cps.cac.gac)
cps.obs.exp.cacgag<-log(dn31.codon.pair$cacgag/cps.cac.gag)
cps.obs.exp.cacgat<-log(dn31.codon.pair$cacgat/cps.cac.gat)
cps.obs.exp.cacgca<-log(dn31.codon.pair$cacgca/cps.cac.gca)
cps.obs.exp.cacgcc<-log(dn31.codon.pair$cacgcc/cps.cac.gcc)
cps.obs.exp.cacgcg<-log(dn31.codon.pair$cacgcg/cps.cac.gcg)
cps.obs.exp.cacgct<-log(dn31.codon.pair$cacgct/cps.cac.gct)
cps.obs.exp.cacgga<-log(dn31.codon.pair$cacgga/cps.cac.gga)
cps.obs.exp.cacggc<-log(dn31.codon.pair$cacggc/cps.cac.ggc)
cps.obs.exp.cacggg<-log(dn31.codon.pair$cacggg/cps.cac.ggg)
cps.obs.exp.cacggt<-log(dn31.codon.pair$cacggt/cps.cac.ggt)
cps.obs.exp.cacgta<-log(dn31.codon.pair$cacgta/cps.cac.gta)
cps.obs.exp.cacgtc<-log(dn31.codon.pair$cacgtc/cps.cac.gtc)
cps.obs.exp.cacgtg<-log(dn31.codon.pair$cacgtg/cps.cac.gtg)
cps.obs.exp.cacgtt<-log(dn31.codon.pair$cacgtt/cps.cac.gtt)

#cps.obs.exp.cactaa<-log(dn31.codon.pair$cactaa/cps.cac.taa)
cps.obs.exp.cactac<-log(dn31.codon.pair$cactac/cps.cac.tac)
#cps.obs.exp.cactag<-log(dn31.codon.pair$cactag/cps.cac.tag)
cps.obs.exp.cactat<-log(dn31.codon.pair$cactat/cps.cac.tat)
cps.obs.exp.cactca<-log(dn31.codon.pair$cactca/cps.cac.tca)
cps.obs.exp.cactcc<-log(dn31.codon.pair$cactcc/cps.cac.tcc)
cps.obs.exp.cactcg<-log(dn31.codon.pair$cactcg/cps.cac.tcg)
cps.obs.exp.cactct<-log(dn31.codon.pair$cactct/cps.cac.tct)
#cps.obs.exp.cactga<-log(dn31.codon.pair$cactga/cps.cac.tga)
cps.obs.exp.cactgc<-log(dn31.codon.pair$cactgc/cps.cac.tgc)
cps.obs.exp.cactgg<-log(dn31.codon.pair$cactgg/cps.cac.tgg)
cps.obs.exp.cactgt<-log(dn31.codon.pair$cactgt/cps.cac.tgt)
cps.obs.exp.cactta<-log(dn31.codon.pair$cactta/cps.cac.tta)
cps.obs.exp.cacttc<-log(dn31.codon.pair$cacttc/cps.cac.ttc)
cps.obs.exp.cacttg<-log(dn31.codon.pair$cacttg/cps.cac.ttg)
cps.obs.exp.cacttt<-log(dn31.codon.pair$cacttt/cps.cac.ttt)









cps.obs.exp.cagaaa<-log(dn31.codon.pair$cagaaa/cps.cag.aaa)
cps.obs.exp.cagaac<-log(dn31.codon.pair$cagaac/cps.cag.aac)
cps.obs.exp.cagaag<-log(dn31.codon.pair$cagaag/cps.cag.aag)
cps.obs.exp.cagaat<-log(dn31.codon.pair$cagaat/cps.cag.aat)
cps.obs.exp.cagaca<-log(dn31.codon.pair$cagaca/cps.cag.aca)
cps.obs.exp.cagacc<-log(dn31.codon.pair$cagacc/cps.cag.acc)
cps.obs.exp.cagacg<-log(dn31.codon.pair$cagacg/cps.cag.acg)
cps.obs.exp.cagact<-log(dn31.codon.pair$cagact/cps.cag.act)
cps.obs.exp.cagaga<-log(dn31.codon.pair$cagaga/cps.cag.aga)
cps.obs.exp.cagagc<-log(dn31.codon.pair$cagagc/cps.cag.agc)
cps.obs.exp.cagagg<-log(dn31.codon.pair$cagagg/cps.cag.agg)
cps.obs.exp.cagagt<-log(dn31.codon.pair$cagagt/cps.cag.agt)
cps.obs.exp.cagata<-log(dn31.codon.pair$cagata/cps.cag.ata)
cps.obs.exp.cagatc<-log(dn31.codon.pair$cagatc/cps.cag.atc)
cps.obs.exp.cagatg<-log(dn31.codon.pair$cagatg/cps.cag.atg)
cps.obs.exp.cagatt<-log(dn31.codon.pair$cagatt/cps.cag.att)

cps.obs.exp.cagcaa<-log(dn31.codon.pair$cagcaa/cps.cag.caa)
cps.obs.exp.cagcac<-log(dn31.codon.pair$cagcac/cps.cag.cac)
cps.obs.exp.cagcag<-log(dn31.codon.pair$cagcag/cps.cag.cag)
cps.obs.exp.cagcat<-log(dn31.codon.pair$cagcat/cps.cag.cat)
cps.obs.exp.cagcca<-log(dn31.codon.pair$cagcca/cps.cag.cca)
cps.obs.exp.cagccc<-log(dn31.codon.pair$cagccc/cps.cag.ccc)
cps.obs.exp.cagccg<-log(dn31.codon.pair$cagccg/cps.cag.ccg)
cps.obs.exp.cagcct<-log(dn31.codon.pair$cagcct/cps.cag.cct)
cps.obs.exp.cagcga<-log(dn31.codon.pair$cagcga/cps.cag.cga)
cps.obs.exp.cagcgc<-log(dn31.codon.pair$cagcgc/cps.cag.cgc)
cps.obs.exp.cagcgg<-log(dn31.codon.pair$cagcgg/cps.cag.cgg)
cps.obs.exp.cagcgt<-log(dn31.codon.pair$cagcgt/cps.cag.cgt)
cps.obs.exp.cagcta<-log(dn31.codon.pair$cagcta/cps.cag.cta)
cps.obs.exp.cagctc<-log(dn31.codon.pair$cagctc/cps.cag.ctc)
cps.obs.exp.cagctg<-log(dn31.codon.pair$cagctg/cps.cag.ctg)
cps.obs.exp.cagctt<-log(dn31.codon.pair$cagctt/cps.cag.ctt)

cps.obs.exp.caggaa<-log(dn31.codon.pair$caggaa/cps.cag.gaa)
cps.obs.exp.caggac<-log(dn31.codon.pair$caggac/cps.cag.gac)
cps.obs.exp.caggag<-log(dn31.codon.pair$caggag/cps.cag.gag)
cps.obs.exp.caggat<-log(dn31.codon.pair$caggat/cps.cag.gat)
cps.obs.exp.caggca<-log(dn31.codon.pair$caggca/cps.cag.gca)
cps.obs.exp.caggcc<-log(dn31.codon.pair$caggcc/cps.cag.gcc)
cps.obs.exp.caggcg<-log(dn31.codon.pair$caggcg/cps.cag.gcg)
cps.obs.exp.caggct<-log(dn31.codon.pair$caggct/cps.cag.gct)
cps.obs.exp.caggga<-log(dn31.codon.pair$caggga/cps.cag.gga)
cps.obs.exp.cagggc<-log(dn31.codon.pair$cagggc/cps.cag.ggc)
cps.obs.exp.cagggg<-log(dn31.codon.pair$cagggg/cps.cag.ggg)
cps.obs.exp.cagggt<-log(dn31.codon.pair$cagggt/cps.cag.ggt)
cps.obs.exp.caggta<-log(dn31.codon.pair$caggta/cps.cag.gta)
cps.obs.exp.caggtc<-log(dn31.codon.pair$caggtc/cps.cag.gtc)
cps.obs.exp.caggtg<-log(dn31.codon.pair$caggtg/cps.cag.gtg)
cps.obs.exp.caggtt<-log(dn31.codon.pair$caggtt/cps.cag.gtt)

#cps.obs.exp.cagtaa<-log(dn31.codon.pair$cagtaa/cps.cag.taa)
cps.obs.exp.cagtac<-log(dn31.codon.pair$cagtac/cps.cag.tac)
#cps.obs.exp.cagtag<-log(dn31.codon.pair$cagtag/cps.cag.tag)
cps.obs.exp.cagtat<-log(dn31.codon.pair$cagtat/cps.cag.tat)
cps.obs.exp.cagtca<-log(dn31.codon.pair$cagtca/cps.cag.tca)
cps.obs.exp.cagtcc<-log(dn31.codon.pair$cagtcc/cps.cag.tcc)
cps.obs.exp.cagtcg<-log(dn31.codon.pair$cagtcg/cps.cag.tcg)
cps.obs.exp.cagtct<-log(dn31.codon.pair$cagtct/cps.cag.tct)
#cps.obs.exp.cagtga<-log(dn31.codon.pair$cagtga/cps.cag.tga)
cps.obs.exp.cagtgc<-log(dn31.codon.pair$cagtgc/cps.cag.tgc)
cps.obs.exp.cagtgg<-log(dn31.codon.pair$cagtgg/cps.cag.tgg)
cps.obs.exp.cagtgt<-log(dn31.codon.pair$cagtgt/cps.cag.tgt)
cps.obs.exp.cagtta<-log(dn31.codon.pair$cagtta/cps.cag.tta)
cps.obs.exp.cagttc<-log(dn31.codon.pair$cagttc/cps.cag.ttc)
cps.obs.exp.cagttg<-log(dn31.codon.pair$cagttg/cps.cag.ttg)
cps.obs.exp.cagttt<-log(dn31.codon.pair$cagttt/cps.cag.ttt)








cps.obs.exp.cataaa<-log(dn31.codon.pair$cataaa/cps.cat.aaa)
cps.obs.exp.cataac<-log(dn31.codon.pair$cataac/cps.cat.aac)
cps.obs.exp.cataag<-log(dn31.codon.pair$cataag/cps.cat.aag)
cps.obs.exp.cataat<-log(dn31.codon.pair$cataat/cps.cat.aat)
cps.obs.exp.cataca<-log(dn31.codon.pair$cataca/cps.cat.aca)
cps.obs.exp.catacc<-log(dn31.codon.pair$catacc/cps.cat.acc)
cps.obs.exp.catacg<-log(dn31.codon.pair$catacg/cps.cat.acg)
cps.obs.exp.catact<-log(dn31.codon.pair$catact/cps.cat.act)
cps.obs.exp.cataga<-log(dn31.codon.pair$cataga/cps.cat.aga)
cps.obs.exp.catagc<-log(dn31.codon.pair$catagc/cps.cat.agc)
cps.obs.exp.catagg<-log(dn31.codon.pair$catagg/cps.cat.agg)
cps.obs.exp.catagt<-log(dn31.codon.pair$catagt/cps.cat.agt)
cps.obs.exp.catata<-log(dn31.codon.pair$catata/cps.cat.ata)
cps.obs.exp.catatc<-log(dn31.codon.pair$catatc/cps.cat.atc)
cps.obs.exp.catatg<-log(dn31.codon.pair$catatg/cps.cat.atg)
cps.obs.exp.catatt<-log(dn31.codon.pair$catatt/cps.cat.att)

cps.obs.exp.catcaa<-log(dn31.codon.pair$catcaa/cps.cat.caa)
cps.obs.exp.catcac<-log(dn31.codon.pair$catcac/cps.cat.cac)
cps.obs.exp.catcag<-log(dn31.codon.pair$catcag/cps.cat.cag)
cps.obs.exp.catcat<-log(dn31.codon.pair$catcat/cps.cat.cat)
cps.obs.exp.catcca<-log(dn31.codon.pair$catcca/cps.cat.cca)
cps.obs.exp.catccc<-log(dn31.codon.pair$catccc/cps.cat.ccc)
cps.obs.exp.catccg<-log(dn31.codon.pair$catccg/cps.cat.ccg)
cps.obs.exp.catcct<-log(dn31.codon.pair$catcct/cps.cat.cct)
cps.obs.exp.catcga<-log(dn31.codon.pair$catcga/cps.cat.cga)
cps.obs.exp.catcgc<-log(dn31.codon.pair$catcgc/cps.cat.cgc)
cps.obs.exp.catcgg<-log(dn31.codon.pair$catcgg/cps.cat.cgg)
cps.obs.exp.catcgt<-log(dn31.codon.pair$catcgt/cps.cat.cgt)
cps.obs.exp.catcta<-log(dn31.codon.pair$catcta/cps.cat.cta)
cps.obs.exp.catctc<-log(dn31.codon.pair$catctc/cps.cat.ctc)
cps.obs.exp.catctg<-log(dn31.codon.pair$catctg/cps.cat.ctg)
cps.obs.exp.catctt<-log(dn31.codon.pair$catctt/cps.cat.ctt)

cps.obs.exp.catgaa<-log(dn31.codon.pair$catgaa/cps.cat.gaa)
cps.obs.exp.catgac<-log(dn31.codon.pair$catgac/cps.cat.gac)
cps.obs.exp.catgag<-log(dn31.codon.pair$catgag/cps.cat.gag)
cps.obs.exp.catgat<-log(dn31.codon.pair$catgat/cps.cat.gat)
cps.obs.exp.catgca<-log(dn31.codon.pair$catgca/cps.cat.gca)
cps.obs.exp.catgcc<-log(dn31.codon.pair$catgcc/cps.cat.gcc)
cps.obs.exp.catgcg<-log(dn31.codon.pair$catgcg/cps.cat.gcg)
cps.obs.exp.catgct<-log(dn31.codon.pair$catgct/cps.cat.gct)
cps.obs.exp.catgga<-log(dn31.codon.pair$catgga/cps.cat.gga)
cps.obs.exp.catggc<-log(dn31.codon.pair$catggc/cps.cat.ggc)
cps.obs.exp.catggg<-log(dn31.codon.pair$catggg/cps.cat.ggg)
cps.obs.exp.catggt<-log(dn31.codon.pair$catggt/cps.cat.ggt)
cps.obs.exp.catgta<-log(dn31.codon.pair$catgta/cps.cat.gta)
cps.obs.exp.catgtc<-log(dn31.codon.pair$catgtc/cps.cat.gtc)
cps.obs.exp.catgtg<-log(dn31.codon.pair$catgtg/cps.cat.gtg)
cps.obs.exp.catgtt<-log(dn31.codon.pair$catgtt/cps.cat.gtt)

#cps.obs.exp.cattaa<-log(dn31.codon.pair$cattaa/cps.cat.taa)
cps.obs.exp.cattac<-log(dn31.codon.pair$cattac/cps.cat.tac)
#cps.obs.exp.cattag<-log(dn31.codon.pair$cattag/cps.cat.tag)
cps.obs.exp.cattat<-log(dn31.codon.pair$cattat/cps.cat.tat)
cps.obs.exp.cattca<-log(dn31.codon.pair$cattca/cps.cat.tca)
cps.obs.exp.cattcc<-log(dn31.codon.pair$cattcc/cps.cat.tcc)
cps.obs.exp.cattcg<-log(dn31.codon.pair$cattcg/cps.cat.tcg)
cps.obs.exp.cattct<-log(dn31.codon.pair$cattct/cps.cat.tct)
#cps.obs.exp.cattga<-log(dn31.codon.pair$cattga/cps.cat.tga)
cps.obs.exp.cattgc<-log(dn31.codon.pair$cattgc/cps.cat.tgc)
cps.obs.exp.cattgg<-log(dn31.codon.pair$cattgg/cps.cat.tgg)
cps.obs.exp.cattgt<-log(dn31.codon.pair$cattgt/cps.cat.tgt)
cps.obs.exp.cattta<-log(dn31.codon.pair$cattta/cps.cat.tta)
cps.obs.exp.catttc<-log(dn31.codon.pair$catttc/cps.cat.ttc)
cps.obs.exp.catttg<-log(dn31.codon.pair$catttg/cps.cat.ttg)
cps.obs.exp.catttt<-log(dn31.codon.pair$catttt/cps.cat.ttt)

















cps.obs.exp.ccaaaa<-log(dn31.codon.pair$ccaaaa/cps.cca.aaa)
cps.obs.exp.ccaaac<-log(dn31.codon.pair$ccaaac/cps.cca.aac)
cps.obs.exp.ccaaag<-log(dn31.codon.pair$ccaaag/cps.cca.aag)
cps.obs.exp.ccaaat<-log(dn31.codon.pair$ccaaat/cps.cca.aat)
cps.obs.exp.ccaaca<-log(dn31.codon.pair$ccaaca/cps.cca.aca)
cps.obs.exp.ccaacc<-log(dn31.codon.pair$ccaacc/cps.cca.acc)
cps.obs.exp.ccaacg<-log(dn31.codon.pair$ccaacg/cps.cca.acg)
cps.obs.exp.ccaact<-log(dn31.codon.pair$ccaact/cps.cca.act)
cps.obs.exp.ccaaga<-log(dn31.codon.pair$ccaaga/cps.cca.aga)
cps.obs.exp.ccaagc<-log(dn31.codon.pair$ccaagc/cps.cca.agc)
cps.obs.exp.ccaagg<-log(dn31.codon.pair$ccaagg/cps.cca.agg)
cps.obs.exp.ccaagt<-log(dn31.codon.pair$ccaagt/cps.cca.agt)
cps.obs.exp.ccaata<-log(dn31.codon.pair$ccaata/cps.cca.ata)
cps.obs.exp.ccaatc<-log(dn31.codon.pair$ccaatc/cps.cca.atc)
cps.obs.exp.ccaatg<-log(dn31.codon.pair$ccaatg/cps.cca.atg)
cps.obs.exp.ccaatt<-log(dn31.codon.pair$ccaatt/cps.cca.att)

cps.obs.exp.ccacaa<-log(dn31.codon.pair$ccacaa/cps.cca.caa)
cps.obs.exp.ccacac<-log(dn31.codon.pair$ccacac/cps.cca.cac)
cps.obs.exp.ccacag<-log(dn31.codon.pair$ccacag/cps.cca.cag)
cps.obs.exp.ccacat<-log(dn31.codon.pair$ccacat/cps.cca.cat)
cps.obs.exp.ccacca<-log(dn31.codon.pair$ccacca/cps.cca.cca)
cps.obs.exp.ccaccc<-log(dn31.codon.pair$ccaccc/cps.cca.ccc)
cps.obs.exp.ccaccg<-log(dn31.codon.pair$ccaccg/cps.cca.ccg)
cps.obs.exp.ccacct<-log(dn31.codon.pair$ccacct/cps.cca.cct)
cps.obs.exp.ccacga<-log(dn31.codon.pair$ccacga/cps.cca.cga)
cps.obs.exp.ccacgc<-log(dn31.codon.pair$ccacgc/cps.cca.cgc)
cps.obs.exp.ccacgg<-log(dn31.codon.pair$ccacgg/cps.cca.cgg)
cps.obs.exp.ccacgt<-log(dn31.codon.pair$ccacgt/cps.cca.cgt)
cps.obs.exp.ccacta<-log(dn31.codon.pair$ccacta/cps.cca.cta)
cps.obs.exp.ccactc<-log(dn31.codon.pair$ccactc/cps.cca.ctc)
cps.obs.exp.ccactg<-log(dn31.codon.pair$ccactg/cps.cca.ctg)
cps.obs.exp.ccactt<-log(dn31.codon.pair$ccactt/cps.cca.ctt)

cps.obs.exp.ccagaa<-log(dn31.codon.pair$ccagaa/cps.cca.gaa)
cps.obs.exp.ccagac<-log(dn31.codon.pair$ccagac/cps.cca.gac)
cps.obs.exp.ccagag<-log(dn31.codon.pair$ccagag/cps.cca.gag)
cps.obs.exp.ccagat<-log(dn31.codon.pair$ccagat/cps.cca.gat)
cps.obs.exp.ccagca<-log(dn31.codon.pair$ccagca/cps.cca.gca)
cps.obs.exp.ccagcc<-log(dn31.codon.pair$ccagcc/cps.cca.gcc)
cps.obs.exp.ccagcg<-log(dn31.codon.pair$ccagcg/cps.cca.gcg)
cps.obs.exp.ccagct<-log(dn31.codon.pair$ccagct/cps.cca.gct)
cps.obs.exp.ccagga<-log(dn31.codon.pair$ccagga/cps.cca.gga)
cps.obs.exp.ccaggc<-log(dn31.codon.pair$ccaggc/cps.cca.ggc)
cps.obs.exp.ccaggg<-log(dn31.codon.pair$ccaggg/cps.cca.ggg)
cps.obs.exp.ccaggt<-log(dn31.codon.pair$ccaggt/cps.cca.ggt)
cps.obs.exp.ccagta<-log(dn31.codon.pair$ccagta/cps.cca.gta)
cps.obs.exp.ccagtc<-log(dn31.codon.pair$ccagtc/cps.cca.gtc)
cps.obs.exp.ccagtg<-log(dn31.codon.pair$ccagtg/cps.cca.gtg)
cps.obs.exp.ccagtt<-log(dn31.codon.pair$ccagtt/cps.cca.gtt)

#cps.obs.exp.ccataa<-log(dn31.codon.pair$ccataa/cps.cca.taa)
cps.obs.exp.ccatac<-log(dn31.codon.pair$ccatac/cps.cca.tac)
#cps.obs.exp.ccatag<-log(dn31.codon.pair$ccatag/cps.cca.tag)
cps.obs.exp.ccatat<-log(dn31.codon.pair$ccatat/cps.cca.tat)
cps.obs.exp.ccatca<-log(dn31.codon.pair$ccatca/cps.cca.tca)
cps.obs.exp.ccatcc<-log(dn31.codon.pair$ccatcc/cps.cca.tcc)
cps.obs.exp.ccatcg<-log(dn31.codon.pair$ccatcg/cps.cca.tcg)
cps.obs.exp.ccatct<-log(dn31.codon.pair$ccatct/cps.cca.tct)
#cps.obs.exp.ccatga<-log(dn31.codon.pair$ccatga/cps.cca.tga)
cps.obs.exp.ccatgc<-log(dn31.codon.pair$ccatgc/cps.cca.tgc)
cps.obs.exp.ccatgg<-log(dn31.codon.pair$ccatgg/cps.cca.tgg)
cps.obs.exp.ccatgt<-log(dn31.codon.pair$ccatgt/cps.cca.tgt)
cps.obs.exp.ccatta<-log(dn31.codon.pair$ccatta/cps.cca.tta)
cps.obs.exp.ccattc<-log(dn31.codon.pair$ccattc/cps.cca.ttc)
cps.obs.exp.ccattg<-log(dn31.codon.pair$ccattg/cps.cca.ttg)
cps.obs.exp.ccattt<-log(dn31.codon.pair$ccattt/cps.cca.ttt)







cps.obs.exp.cccaaa<-log(dn31.codon.pair$cccaaa/cps.ccc.aaa)
cps.obs.exp.cccaac<-log(dn31.codon.pair$cccaac/cps.ccc.aac)
cps.obs.exp.cccaag<-log(dn31.codon.pair$cccaag/cps.ccc.aag)
cps.obs.exp.cccaat<-log(dn31.codon.pair$cccaat/cps.ccc.aat)
cps.obs.exp.cccaca<-log(dn31.codon.pair$cccaca/cps.ccc.aca)
cps.obs.exp.cccacc<-log(dn31.codon.pair$cccacc/cps.ccc.acc)
cps.obs.exp.cccacg<-log(dn31.codon.pair$cccacg/cps.ccc.acg)
cps.obs.exp.cccact<-log(dn31.codon.pair$cccact/cps.ccc.act)
cps.obs.exp.cccaga<-log(dn31.codon.pair$cccaga/cps.ccc.aga)
cps.obs.exp.cccagc<-log(dn31.codon.pair$cccagc/cps.ccc.agc)
cps.obs.exp.cccagg<-log(dn31.codon.pair$cccagg/cps.ccc.agg)
cps.obs.exp.cccagt<-log(dn31.codon.pair$cccagt/cps.ccc.agt)
cps.obs.exp.cccata<-log(dn31.codon.pair$cccata/cps.ccc.ata)
cps.obs.exp.cccatc<-log(dn31.codon.pair$cccatc/cps.ccc.atc)
cps.obs.exp.cccatg<-log(dn31.codon.pair$cccatg/cps.ccc.atg)
cps.obs.exp.cccatt<-log(dn31.codon.pair$cccatt/cps.ccc.att)

cps.obs.exp.ccccaa<-log(dn31.codon.pair$ccccaa/cps.ccc.caa)
cps.obs.exp.ccccac<-log(dn31.codon.pair$ccccac/cps.ccc.cac)
cps.obs.exp.ccccag<-log(dn31.codon.pair$ccccag/cps.ccc.cag)
cps.obs.exp.ccccat<-log(dn31.codon.pair$ccccat/cps.ccc.cat)
cps.obs.exp.ccccca<-log(dn31.codon.pair$ccccca/cps.ccc.cca)
cps.obs.exp.cccccc<-log(dn31.codon.pair$cccccc/cps.ccc.ccc)
cps.obs.exp.cccccg<-log(dn31.codon.pair$cccccg/cps.ccc.ccg)
cps.obs.exp.ccccct<-log(dn31.codon.pair$ccccct/cps.ccc.cct)
cps.obs.exp.ccccga<-log(dn31.codon.pair$ccccga/cps.ccc.cga)
cps.obs.exp.ccccgc<-log(dn31.codon.pair$ccccgc/cps.ccc.cgc)
cps.obs.exp.ccccgg<-log(dn31.codon.pair$ccccgg/cps.ccc.cgg)
cps.obs.exp.ccccgt<-log(dn31.codon.pair$ccccgt/cps.ccc.cgt)
cps.obs.exp.ccccta<-log(dn31.codon.pair$ccccta/cps.ccc.cta)
cps.obs.exp.cccctc<-log(dn31.codon.pair$cccctc/cps.ccc.ctc)
cps.obs.exp.cccctg<-log(dn31.codon.pair$cccctg/cps.ccc.ctg)
cps.obs.exp.cccctt<-log(dn31.codon.pair$cccctt/cps.ccc.ctt)

cps.obs.exp.cccgaa<-log(dn31.codon.pair$cccgaa/cps.ccc.gaa)
cps.obs.exp.cccgac<-log(dn31.codon.pair$cccgac/cps.ccc.gac)
cps.obs.exp.cccgag<-log(dn31.codon.pair$cccgag/cps.ccc.gag)
cps.obs.exp.cccgat<-log(dn31.codon.pair$cccgat/cps.ccc.gat)
cps.obs.exp.cccgca<-log(dn31.codon.pair$cccgca/cps.ccc.gca)
cps.obs.exp.cccgcc<-log(dn31.codon.pair$cccgcc/cps.ccc.gcc)
cps.obs.exp.cccgcg<-log(dn31.codon.pair$cccgcg/cps.ccc.gcg)
cps.obs.exp.cccgct<-log(dn31.codon.pair$cccgct/cps.ccc.gct)
cps.obs.exp.cccgga<-log(dn31.codon.pair$cccgga/cps.ccc.gga)
cps.obs.exp.cccggc<-log(dn31.codon.pair$cccggc/cps.ccc.ggc)
cps.obs.exp.cccggg<-log(dn31.codon.pair$cccggg/cps.ccc.ggg)
cps.obs.exp.cccggt<-log(dn31.codon.pair$cccggt/cps.ccc.ggt)
cps.obs.exp.cccgta<-log(dn31.codon.pair$cccgta/cps.ccc.gta)
cps.obs.exp.cccgtc<-log(dn31.codon.pair$cccgtc/cps.ccc.gtc)
cps.obs.exp.cccgtg<-log(dn31.codon.pair$cccgtg/cps.ccc.gtg)
cps.obs.exp.cccgtt<-log(dn31.codon.pair$cccgtt/cps.ccc.gtt)

#cps.obs.exp.ccctaa<-log(dn31.codon.pair$ccctaa/cps.ccc.taa)
cps.obs.exp.ccctac<-log(dn31.codon.pair$ccctac/cps.ccc.tac)
#cps.obs.exp.ccctag<-log(dn31.codon.pair$ccctag/cps.ccc.tag)
cps.obs.exp.ccctat<-log(dn31.codon.pair$ccctat/cps.ccc.tat)
cps.obs.exp.ccctca<-log(dn31.codon.pair$ccctca/cps.ccc.tca)
cps.obs.exp.ccctcc<-log(dn31.codon.pair$ccctcc/cps.ccc.tcc)
cps.obs.exp.ccctcg<-log(dn31.codon.pair$ccctcg/cps.ccc.tcg)
cps.obs.exp.ccctct<-log(dn31.codon.pair$ccctct/cps.ccc.tct)
#cps.obs.exp.ccctga<-log(dn31.codon.pair$ccctga/cps.ccc.tga)
cps.obs.exp.ccctgc<-log(dn31.codon.pair$ccctgc/cps.ccc.tgc)
cps.obs.exp.ccctgg<-log(dn31.codon.pair$ccctgg/cps.ccc.tgg)
cps.obs.exp.ccctgt<-log(dn31.codon.pair$ccctgt/cps.ccc.tgt)
cps.obs.exp.ccctta<-log(dn31.codon.pair$ccctta/cps.ccc.tta)
cps.obs.exp.cccttc<-log(dn31.codon.pair$cccttc/cps.ccc.ttc)
cps.obs.exp.cccttg<-log(dn31.codon.pair$cccttg/cps.ccc.ttg)
cps.obs.exp.cccttt<-log(dn31.codon.pair$cccttt/cps.ccc.ttt)









cps.obs.exp.ccgaaa<-log(dn31.codon.pair$ccgaaa/cps.ccg.aaa)
cps.obs.exp.ccgaac<-log(dn31.codon.pair$ccgaac/cps.ccg.aac)
cps.obs.exp.ccgaag<-log(dn31.codon.pair$ccgaag/cps.ccg.aag)
cps.obs.exp.ccgaat<-log(dn31.codon.pair$ccgaat/cps.ccg.aat)
cps.obs.exp.ccgaca<-log(dn31.codon.pair$ccgaca/cps.ccg.aca)
cps.obs.exp.ccgacc<-log(dn31.codon.pair$ccgacc/cps.ccg.acc)
cps.obs.exp.ccgacg<-log(dn31.codon.pair$ccgacg/cps.ccg.acg)
cps.obs.exp.ccgact<-log(dn31.codon.pair$ccgact/cps.ccg.act)
cps.obs.exp.ccgaga<-log(dn31.codon.pair$ccgaga/cps.ccg.aga)
cps.obs.exp.ccgagc<-log(dn31.codon.pair$ccgagc/cps.ccg.agc)
cps.obs.exp.ccgagg<-log(dn31.codon.pair$ccgagg/cps.ccg.agg)
cps.obs.exp.ccgagt<-log(dn31.codon.pair$ccgagt/cps.ccg.agt)
cps.obs.exp.ccgata<-log(dn31.codon.pair$ccgata/cps.ccg.ata)
cps.obs.exp.ccgatc<-log(dn31.codon.pair$ccgatc/cps.ccg.atc)
cps.obs.exp.ccgatg<-log(dn31.codon.pair$ccgatg/cps.ccg.atg)
cps.obs.exp.ccgatt<-log(dn31.codon.pair$ccgatt/cps.ccg.att)

cps.obs.exp.ccgcaa<-log(dn31.codon.pair$ccgcaa/cps.ccg.caa)
cps.obs.exp.ccgcac<-log(dn31.codon.pair$ccgcac/cps.ccg.cac)
cps.obs.exp.ccgcag<-log(dn31.codon.pair$ccgcag/cps.ccg.cag)
cps.obs.exp.ccgcat<-log(dn31.codon.pair$ccgcat/cps.ccg.cat)
cps.obs.exp.ccgcca<-log(dn31.codon.pair$ccgcca/cps.ccg.cca)
cps.obs.exp.ccgccc<-log(dn31.codon.pair$ccgccc/cps.ccg.ccc)
cps.obs.exp.ccgccg<-log(dn31.codon.pair$ccgccg/cps.ccg.ccg)
cps.obs.exp.ccgcct<-log(dn31.codon.pair$ccgcct/cps.ccg.cct)
cps.obs.exp.ccgcga<-log(dn31.codon.pair$ccgcga/cps.ccg.cga)
cps.obs.exp.ccgcgc<-log(dn31.codon.pair$ccgcgc/cps.ccg.cgc)
cps.obs.exp.ccgcgg<-log(dn31.codon.pair$ccgcgg/cps.ccg.cgg)
cps.obs.exp.ccgcgt<-log(dn31.codon.pair$ccgcgt/cps.ccg.cgt)
cps.obs.exp.ccgcta<-log(dn31.codon.pair$ccgcta/cps.ccg.cta)
cps.obs.exp.ccgctc<-log(dn31.codon.pair$ccgctc/cps.ccg.ctc)
cps.obs.exp.ccgctg<-log(dn31.codon.pair$ccgctg/cps.ccg.ctg)
cps.obs.exp.ccgctt<-log(dn31.codon.pair$ccgctt/cps.ccg.ctt)

cps.obs.exp.ccggaa<-log(dn31.codon.pair$ccggaa/cps.ccg.gaa)
cps.obs.exp.ccggac<-log(dn31.codon.pair$ccggac/cps.ccg.gac)
cps.obs.exp.ccggag<-log(dn31.codon.pair$ccggag/cps.ccg.gag)
cps.obs.exp.ccggat<-log(dn31.codon.pair$ccggat/cps.ccg.gat)
cps.obs.exp.ccggca<-log(dn31.codon.pair$ccggca/cps.ccg.gca)
cps.obs.exp.ccggcc<-log(dn31.codon.pair$ccggcc/cps.ccg.gcc)
cps.obs.exp.ccggcg<-log(dn31.codon.pair$ccggcg/cps.ccg.gcg)
cps.obs.exp.ccggct<-log(dn31.codon.pair$ccggct/cps.ccg.gct)
cps.obs.exp.ccggga<-log(dn31.codon.pair$ccggga/cps.ccg.gga)
cps.obs.exp.ccgggc<-log(dn31.codon.pair$ccgggc/cps.ccg.ggc)
cps.obs.exp.ccgggg<-log(dn31.codon.pair$ccgggg/cps.ccg.ggg)
cps.obs.exp.ccgggt<-log(dn31.codon.pair$ccgggt/cps.ccg.ggt)
cps.obs.exp.ccggta<-log(dn31.codon.pair$ccggta/cps.ccg.gta)
cps.obs.exp.ccggtc<-log(dn31.codon.pair$ccggtc/cps.ccg.gtc)
cps.obs.exp.ccggtg<-log(dn31.codon.pair$ccggtg/cps.ccg.gtg)
cps.obs.exp.ccggtt<-log(dn31.codon.pair$ccggtt/cps.ccg.gtt)

#cps.obs.exp.ccgtaa<-log(dn31.codon.pair$ccgtaa/cps.ccg.taa)
cps.obs.exp.ccgtac<-log(dn31.codon.pair$ccgtac/cps.ccg.tac)
#cps.obs.exp.ccgtag<-log(dn31.codon.pair$ccgtag/cps.ccg.tag)
cps.obs.exp.ccgtat<-log(dn31.codon.pair$ccgtat/cps.ccg.tat)
cps.obs.exp.ccgtca<-log(dn31.codon.pair$ccgtca/cps.ccg.tca)
cps.obs.exp.ccgtcc<-log(dn31.codon.pair$ccgtcc/cps.ccg.tcc)
cps.obs.exp.ccgtcg<-log(dn31.codon.pair$ccgtcg/cps.ccg.tcg)
cps.obs.exp.ccgtct<-log(dn31.codon.pair$ccgtct/cps.ccg.tct)
#cps.obs.exp.ccgtga<-log(dn31.codon.pair$ccgtga/cps.ccg.tga)
cps.obs.exp.ccgtgc<-log(dn31.codon.pair$ccgtgc/cps.ccg.tgc)
cps.obs.exp.ccgtgg<-log(dn31.codon.pair$ccgtgg/cps.ccg.tgg)
cps.obs.exp.ccgtgt<-log(dn31.codon.pair$ccgtgt/cps.ccg.tgt)
cps.obs.exp.ccgtta<-log(dn31.codon.pair$ccgtta/cps.ccg.tta)
cps.obs.exp.ccgttc<-log(dn31.codon.pair$ccgttc/cps.ccg.ttc)
cps.obs.exp.ccgttg<-log(dn31.codon.pair$ccgttg/cps.ccg.ttg)
cps.obs.exp.ccgttt<-log(dn31.codon.pair$ccgttt/cps.ccg.ttt)







cps.obs.exp.cctaaa<-log(dn31.codon.pair$cctaaa/cps.cct.aaa)
cps.obs.exp.cctaac<-log(dn31.codon.pair$cctaac/cps.cct.aac)
cps.obs.exp.cctaag<-log(dn31.codon.pair$cctaag/cps.cct.aag)
cps.obs.exp.cctaat<-log(dn31.codon.pair$cctaat/cps.cct.aat)
cps.obs.exp.cctaca<-log(dn31.codon.pair$cctaca/cps.cct.aca)
cps.obs.exp.cctacc<-log(dn31.codon.pair$cctacc/cps.cct.acc)
cps.obs.exp.cctacg<-log(dn31.codon.pair$cctacg/cps.cct.acg)
cps.obs.exp.cctact<-log(dn31.codon.pair$cctact/cps.cct.act)
cps.obs.exp.cctaga<-log(dn31.codon.pair$cctaga/cps.cct.aga)
cps.obs.exp.cctagc<-log(dn31.codon.pair$cctagc/cps.cct.agc)
cps.obs.exp.cctagg<-log(dn31.codon.pair$cctagg/cps.cct.agg)
cps.obs.exp.cctagt<-log(dn31.codon.pair$cctagt/cps.cct.agt)
cps.obs.exp.cctata<-log(dn31.codon.pair$cctata/cps.cct.ata)
cps.obs.exp.cctatc<-log(dn31.codon.pair$cctatc/cps.cct.atc)
cps.obs.exp.cctatg<-log(dn31.codon.pair$cctatg/cps.cct.atg)
cps.obs.exp.cctatt<-log(dn31.codon.pair$cctatt/cps.cct.att)

cps.obs.exp.cctcaa<-log(dn31.codon.pair$cctcaa/cps.cct.caa)
cps.obs.exp.cctcac<-log(dn31.codon.pair$cctcac/cps.cct.cac)
cps.obs.exp.cctcag<-log(dn31.codon.pair$cctcag/cps.cct.cag)
cps.obs.exp.cctcat<-log(dn31.codon.pair$cctcat/cps.cct.cat)
cps.obs.exp.cctcca<-log(dn31.codon.pair$cctcca/cps.cct.cca)
cps.obs.exp.cctccc<-log(dn31.codon.pair$cctccc/cps.cct.ccc)
cps.obs.exp.cctccg<-log(dn31.codon.pair$cctccg/cps.cct.ccg)
cps.obs.exp.cctcct<-log(dn31.codon.pair$cctcct/cps.cct.cct)
cps.obs.exp.cctcga<-log(dn31.codon.pair$cctcga/cps.cct.cga)
cps.obs.exp.cctcgc<-log(dn31.codon.pair$cctcgc/cps.cct.cgc)
cps.obs.exp.cctcgg<-log(dn31.codon.pair$cctcgg/cps.cct.cgg)
cps.obs.exp.cctcgt<-log(dn31.codon.pair$cctcgt/cps.cct.cgt)
cps.obs.exp.cctcta<-log(dn31.codon.pair$cctcta/cps.cct.cta)
cps.obs.exp.cctctc<-log(dn31.codon.pair$cctctc/cps.cct.ctc)
cps.obs.exp.cctctg<-log(dn31.codon.pair$cctctg/cps.cct.ctg)
cps.obs.exp.cctctt<-log(dn31.codon.pair$cctctt/cps.cct.ctt)

cps.obs.exp.cctgaa<-log(dn31.codon.pair$cctgaa/cps.cct.gaa)
cps.obs.exp.cctgac<-log(dn31.codon.pair$cctgac/cps.cct.gac)
cps.obs.exp.cctgag<-log(dn31.codon.pair$cctgag/cps.cct.gag)
cps.obs.exp.cctgat<-log(dn31.codon.pair$cctgat/cps.cct.gat)
cps.obs.exp.cctgca<-log(dn31.codon.pair$cctgca/cps.cct.gca)
cps.obs.exp.cctgcc<-log(dn31.codon.pair$cctgcc/cps.cct.gcc)
cps.obs.exp.cctgcg<-log(dn31.codon.pair$cctgcg/cps.cct.gcg)
cps.obs.exp.cctgct<-log(dn31.codon.pair$cctgct/cps.cct.gct)
cps.obs.exp.cctgga<-log(dn31.codon.pair$cctgga/cps.cct.gga)
cps.obs.exp.cctggc<-log(dn31.codon.pair$cctggc/cps.cct.ggc)
cps.obs.exp.cctggg<-log(dn31.codon.pair$cctggg/cps.cct.ggg)
cps.obs.exp.cctggt<-log(dn31.codon.pair$cctggt/cps.cct.ggt)
cps.obs.exp.cctgta<-log(dn31.codon.pair$cctgta/cps.cct.gta)
cps.obs.exp.cctgtc<-log(dn31.codon.pair$cctgtc/cps.cct.gtc)
cps.obs.exp.cctgtg<-log(dn31.codon.pair$cctgtg/cps.cct.gtg)
cps.obs.exp.cctgtt<-log(dn31.codon.pair$cctgtt/cps.cct.gtt)

#cps.obs.exp.ccttaa<-log(dn31.codon.pair$ccttaa/cps.cct.taa)
cps.obs.exp.ccttac<-log(dn31.codon.pair$ccttac/cps.cct.tac)
#cps.obs.exp.ccttag<-log(dn31.codon.pair$ccttag/cps.cct.tag)
cps.obs.exp.ccttat<-log(dn31.codon.pair$ccttat/cps.cct.tat)
cps.obs.exp.ccttca<-log(dn31.codon.pair$ccttca/cps.cct.tca)
cps.obs.exp.ccttcc<-log(dn31.codon.pair$ccttcc/cps.cct.tcc)
cps.obs.exp.ccttcg<-log(dn31.codon.pair$ccttcg/cps.cct.tcg)
cps.obs.exp.ccttct<-log(dn31.codon.pair$ccttct/cps.cct.tct)
#cps.obs.exp.ccttga<-log(dn31.codon.pair$ccttga/cps.cct.tga)
cps.obs.exp.ccttgc<-log(dn31.codon.pair$ccttgc/cps.cct.tgc)
cps.obs.exp.ccttgg<-log(dn31.codon.pair$ccttgg/cps.cct.tgg)
cps.obs.exp.ccttgt<-log(dn31.codon.pair$ccttgt/cps.cct.tgt)
cps.obs.exp.ccttta<-log(dn31.codon.pair$ccttta/cps.cct.tta)
cps.obs.exp.cctttc<-log(dn31.codon.pair$cctttc/cps.cct.ttc)
cps.obs.exp.cctttg<-log(dn31.codon.pair$cctttg/cps.cct.ttg)
cps.obs.exp.cctttt<-log(dn31.codon.pair$cctttt/cps.cct.ttt)



















cps.obs.exp.cgaaaa<-log(dn31.codon.pair$cgaaaa/cps.cga.aaa)
cps.obs.exp.cgaaac<-log(dn31.codon.pair$cgaaac/cps.cga.aac)
cps.obs.exp.cgaaag<-log(dn31.codon.pair$cgaaag/cps.cga.aag)
cps.obs.exp.cgaaat<-log(dn31.codon.pair$cgaaat/cps.cga.aat)
cps.obs.exp.cgaaca<-log(dn31.codon.pair$cgaaca/cps.cga.aca)
cps.obs.exp.cgaacc<-log(dn31.codon.pair$cgaacc/cps.cga.acc)
cps.obs.exp.cgaacg<-log(dn31.codon.pair$cgaacg/cps.cga.acg)
cps.obs.exp.cgaact<-log(dn31.codon.pair$cgaact/cps.cga.act)
cps.obs.exp.cgaaga<-log(dn31.codon.pair$cgaaga/cps.cga.aga)
cps.obs.exp.cgaagc<-log(dn31.codon.pair$cgaagc/cps.cga.agc)
cps.obs.exp.cgaagg<-log(dn31.codon.pair$cgaagg/cps.cga.agg)
cps.obs.exp.cgaagt<-log(dn31.codon.pair$cgaagt/cps.cga.agt)
cps.obs.exp.cgaata<-log(dn31.codon.pair$cgaata/cps.cga.ata)
cps.obs.exp.cgaatc<-log(dn31.codon.pair$cgaatc/cps.cga.atc)
cps.obs.exp.cgaatg<-log(dn31.codon.pair$cgaatg/cps.cga.atg)
cps.obs.exp.cgaatt<-log(dn31.codon.pair$cgaatt/cps.cga.att)

cps.obs.exp.cgacaa<-log(dn31.codon.pair$cgacaa/cps.cga.caa)
cps.obs.exp.cgacac<-log(dn31.codon.pair$cgacac/cps.cga.cac)
cps.obs.exp.cgacag<-log(dn31.codon.pair$cgacag/cps.cga.cag)
cps.obs.exp.cgacat<-log(dn31.codon.pair$cgacat/cps.cga.cat)
cps.obs.exp.cgacca<-log(dn31.codon.pair$cgacca/cps.cga.cca)
cps.obs.exp.cgaccc<-log(dn31.codon.pair$cgaccc/cps.cga.ccc)
cps.obs.exp.cgaccg<-log(dn31.codon.pair$cgaccg/cps.cga.ccg)
cps.obs.exp.cgacct<-log(dn31.codon.pair$cgacct/cps.cga.cct)
cps.obs.exp.cgacga<-log(dn31.codon.pair$cgacga/cps.cga.cga)
cps.obs.exp.cgacgc<-log(dn31.codon.pair$cgacgc/cps.cga.cgc)
cps.obs.exp.cgacgg<-log(dn31.codon.pair$cgacgg/cps.cga.cgg)
cps.obs.exp.cgacgt<-log(dn31.codon.pair$cgacgt/cps.cga.cgt)
cps.obs.exp.cgacta<-log(dn31.codon.pair$cgacta/cps.cga.cta)
cps.obs.exp.cgactc<-log(dn31.codon.pair$cgactc/cps.cga.ctc)
cps.obs.exp.cgactg<-log(dn31.codon.pair$cgactg/cps.cga.ctg)
cps.obs.exp.cgactt<-log(dn31.codon.pair$cgactt/cps.cga.ctt)

cps.obs.exp.cgagaa<-log(dn31.codon.pair$cgagaa/cps.cga.gaa)
cps.obs.exp.cgagac<-log(dn31.codon.pair$cgagac/cps.cga.gac)
cps.obs.exp.cgagag<-log(dn31.codon.pair$cgagag/cps.cga.gag)
cps.obs.exp.cgagat<-log(dn31.codon.pair$cgagat/cps.cga.gat)
cps.obs.exp.cgagca<-log(dn31.codon.pair$cgagca/cps.cga.gca)
cps.obs.exp.cgagcc<-log(dn31.codon.pair$cgagcc/cps.cga.gcc)
cps.obs.exp.cgagcg<-log(dn31.codon.pair$cgagcg/cps.cga.gcg)
cps.obs.exp.cgagct<-log(dn31.codon.pair$cgagct/cps.cga.gct)
cps.obs.exp.cgagga<-log(dn31.codon.pair$cgagga/cps.cga.gga)
cps.obs.exp.cgaggc<-log(dn31.codon.pair$cgaggc/cps.cga.ggc)
cps.obs.exp.cgaggg<-log(dn31.codon.pair$cgaggg/cps.cga.ggg)
cps.obs.exp.cgaggt<-log(dn31.codon.pair$cgaggt/cps.cga.ggt)
cps.obs.exp.cgagta<-log(dn31.codon.pair$cgagta/cps.cga.gta)
cps.obs.exp.cgagtc<-log(dn31.codon.pair$cgagtc/cps.cga.gtc)
cps.obs.exp.cgagtg<-log(dn31.codon.pair$cgagtg/cps.cga.gtg)
cps.obs.exp.cgagtt<-log(dn31.codon.pair$cgagtt/cps.cga.gtt)

#cps.obs.exp.cgataa<-log(dn31.codon.pair$cgataa/cps.cga.taa)
cps.obs.exp.cgatac<-log(dn31.codon.pair$cgatac/cps.cga.tac)
#cps.obs.exp.cgatag<-log(dn31.codon.pair$cgatag/cps.cga.tag)
cps.obs.exp.cgatat<-log(dn31.codon.pair$cgatat/cps.cga.tat)
cps.obs.exp.cgatca<-log(dn31.codon.pair$cgatca/cps.cga.tca)
cps.obs.exp.cgatcc<-log(dn31.codon.pair$cgatcc/cps.cga.tcc)
cps.obs.exp.cgatcg<-log(dn31.codon.pair$cgatcg/cps.cga.tcg)
cps.obs.exp.cgatct<-log(dn31.codon.pair$cgatct/cps.cga.tct)
#cps.obs.exp.cgatga<-log(dn31.codon.pair$cgatga/cps.cga.tga)
cps.obs.exp.cgatgc<-log(dn31.codon.pair$cgatgc/cps.cga.tgc)
cps.obs.exp.cgatgg<-log(dn31.codon.pair$cgatgg/cps.cga.tgg)
cps.obs.exp.cgatgt<-log(dn31.codon.pair$cgatgt/cps.cga.tgt)
cps.obs.exp.cgatta<-log(dn31.codon.pair$cgatta/cps.cga.tta)
cps.obs.exp.cgattc<-log(dn31.codon.pair$cgattc/cps.cga.ttc)
cps.obs.exp.cgattg<-log(dn31.codon.pair$cgattg/cps.cga.ttg)
cps.obs.exp.cgattt<-log(dn31.codon.pair$cgattt/cps.cga.ttt)








cps.obs.exp.cgcaaa<-log(dn31.codon.pair$cgcaaa/cps.cgc.aaa)
cps.obs.exp.cgcaac<-log(dn31.codon.pair$cgcaac/cps.cgc.aac)
cps.obs.exp.cgcaag<-log(dn31.codon.pair$cgcaag/cps.cgc.aag)
cps.obs.exp.cgcaat<-log(dn31.codon.pair$cgcaat/cps.cgc.aat)
cps.obs.exp.cgcaca<-log(dn31.codon.pair$cgcaca/cps.cgc.aca)
cps.obs.exp.cgcacc<-log(dn31.codon.pair$cgcacc/cps.cgc.acc)
cps.obs.exp.cgcacg<-log(dn31.codon.pair$cgcacg/cps.cgc.acg)
cps.obs.exp.cgcact<-log(dn31.codon.pair$cgcact/cps.cgc.act)
cps.obs.exp.cgcaga<-log(dn31.codon.pair$cgcaga/cps.cgc.aga)
cps.obs.exp.cgcagc<-log(dn31.codon.pair$cgcagc/cps.cgc.agc)
cps.obs.exp.cgcagg<-log(dn31.codon.pair$cgcagg/cps.cgc.agg)
cps.obs.exp.cgcagt<-log(dn31.codon.pair$cgcagt/cps.cgc.agt)
cps.obs.exp.cgcata<-log(dn31.codon.pair$cgcata/cps.cgc.ata)
cps.obs.exp.cgcatc<-log(dn31.codon.pair$cgcatc/cps.cgc.atc)
cps.obs.exp.cgcatg<-log(dn31.codon.pair$cgcatg/cps.cgc.atg)
cps.obs.exp.cgcatt<-log(dn31.codon.pair$cgcatt/cps.cgc.att)

cps.obs.exp.cgccaa<-log(dn31.codon.pair$cgccaa/cps.cgc.caa)
cps.obs.exp.cgccac<-log(dn31.codon.pair$cgccac/cps.cgc.cac)
cps.obs.exp.cgccag<-log(dn31.codon.pair$cgccag/cps.cgc.cag)
cps.obs.exp.cgccat<-log(dn31.codon.pair$cgccat/cps.cgc.cat)
cps.obs.exp.cgccca<-log(dn31.codon.pair$cgccca/cps.cgc.cca)
cps.obs.exp.cgcccc<-log(dn31.codon.pair$cgcccc/cps.cgc.ccc)
cps.obs.exp.cgcccg<-log(dn31.codon.pair$cgcccg/cps.cgc.ccg)
cps.obs.exp.cgccct<-log(dn31.codon.pair$cgccct/cps.cgc.cct)
cps.obs.exp.cgccga<-log(dn31.codon.pair$cgccga/cps.cgc.cga)
cps.obs.exp.cgccgc<-log(dn31.codon.pair$cgccgc/cps.cgc.cgc)
cps.obs.exp.cgccgg<-log(dn31.codon.pair$cgccgg/cps.cgc.cgg)
cps.obs.exp.cgccgt<-log(dn31.codon.pair$cgccgt/cps.cgc.cgt)
cps.obs.exp.cgccta<-log(dn31.codon.pair$cgccta/cps.cgc.cta)
cps.obs.exp.cgcctc<-log(dn31.codon.pair$cgcctc/cps.cgc.ctc)
cps.obs.exp.cgcctg<-log(dn31.codon.pair$cgcctg/cps.cgc.ctg)
cps.obs.exp.cgcctt<-log(dn31.codon.pair$cgcctt/cps.cgc.ctt)

cps.obs.exp.cgcgaa<-log(dn31.codon.pair$cgcgaa/cps.cgc.gaa)
cps.obs.exp.cgcgac<-log(dn31.codon.pair$cgcgac/cps.cgc.gac)
cps.obs.exp.cgcgag<-log(dn31.codon.pair$cgcgag/cps.cgc.gag)
cps.obs.exp.cgcgat<-log(dn31.codon.pair$cgcgat/cps.cgc.gat)
cps.obs.exp.cgcgca<-log(dn31.codon.pair$cgcgca/cps.cgc.gca)
cps.obs.exp.cgcgcc<-log(dn31.codon.pair$cgcgcc/cps.cgc.gcc)
cps.obs.exp.cgcgcg<-log(dn31.codon.pair$cgcgcg/cps.cgc.gcg)
cps.obs.exp.cgcgct<-log(dn31.codon.pair$cgcgct/cps.cgc.gct)
cps.obs.exp.cgcgga<-log(dn31.codon.pair$cgcgga/cps.cgc.gga)
cps.obs.exp.cgcggc<-log(dn31.codon.pair$cgcggc/cps.cgc.ggc)
cps.obs.exp.cgcggg<-log(dn31.codon.pair$cgcggg/cps.cgc.ggg)
cps.obs.exp.cgcggt<-log(dn31.codon.pair$cgcggt/cps.cgc.ggt)
cps.obs.exp.cgcgta<-log(dn31.codon.pair$cgcgta/cps.cgc.gta)
cps.obs.exp.cgcgtc<-log(dn31.codon.pair$cgcgtc/cps.cgc.gtc)
cps.obs.exp.cgcgtg<-log(dn31.codon.pair$cgcgtg/cps.cgc.gtg)
cps.obs.exp.cgcgtt<-log(dn31.codon.pair$cgcgtt/cps.cgc.gtt)

#cps.obs.exp.cgctaa<-log(dn31.codon.pair$cgctaa/cps.cgc.taa)
cps.obs.exp.cgctac<-log(dn31.codon.pair$cgctac/cps.cgc.tac)
#cps.obs.exp.cgctag<-log(dn31.codon.pair$cgctag/cps.cgc.tag)
cps.obs.exp.cgctat<-log(dn31.codon.pair$cgctat/cps.cgc.tat)
cps.obs.exp.cgctca<-log(dn31.codon.pair$cgctca/cps.cgc.tca)
cps.obs.exp.cgctcc<-log(dn31.codon.pair$cgctcc/cps.cgc.tcc)
cps.obs.exp.cgctcg<-log(dn31.codon.pair$cgctcg/cps.cgc.tcg)
cps.obs.exp.cgctct<-log(dn31.codon.pair$cgctct/cps.cgc.tct)
#cps.obs.exp.cgctga<-log(dn31.codon.pair$cgctga/cps.cgc.tga)
cps.obs.exp.cgctgc<-log(dn31.codon.pair$cgctgc/cps.cgc.tgc)
cps.obs.exp.cgctgg<-log(dn31.codon.pair$cgctgg/cps.cgc.tgg)
cps.obs.exp.cgctgt<-log(dn31.codon.pair$cgctgt/cps.cgc.tgt)
cps.obs.exp.cgctta<-log(dn31.codon.pair$cgctta/cps.cgc.tta)
cps.obs.exp.cgcttc<-log(dn31.codon.pair$cgcttc/cps.cgc.ttc)
cps.obs.exp.cgcttg<-log(dn31.codon.pair$cgcttg/cps.cgc.ttg)
cps.obs.exp.cgcttt<-log(dn31.codon.pair$cgcttt/cps.cgc.ttt)










cps.obs.exp.cggaaa<-log(dn31.codon.pair$cggaaa/cps.cgg.aaa)
cps.obs.exp.cggaac<-log(dn31.codon.pair$cggaac/cps.cgg.aac)
cps.obs.exp.cggaag<-log(dn31.codon.pair$cggaag/cps.cgg.aag)
cps.obs.exp.cggaat<-log(dn31.codon.pair$cggaat/cps.cgg.aat)
cps.obs.exp.cggaca<-log(dn31.codon.pair$cggaca/cps.cgg.aca)
cps.obs.exp.cggacc<-log(dn31.codon.pair$cggacc/cps.cgg.acc)
cps.obs.exp.cggacg<-log(dn31.codon.pair$cggacg/cps.cgg.acg)
cps.obs.exp.cggact<-log(dn31.codon.pair$cggact/cps.cgg.act)
cps.obs.exp.cggaga<-log(dn31.codon.pair$cggaga/cps.cgg.aga)
cps.obs.exp.cggagc<-log(dn31.codon.pair$cggagc/cps.cgg.agc)
cps.obs.exp.cggagg<-log(dn31.codon.pair$cggagg/cps.cgg.agg)
cps.obs.exp.cggagt<-log(dn31.codon.pair$cggagt/cps.cgg.agt)
cps.obs.exp.cggata<-log(dn31.codon.pair$cggata/cps.cgg.ata)
cps.obs.exp.cggatc<-log(dn31.codon.pair$cggatc/cps.cgg.atc)
cps.obs.exp.cggatg<-log(dn31.codon.pair$cggatg/cps.cgg.atg)
cps.obs.exp.cggatt<-log(dn31.codon.pair$cggatt/cps.cgg.att)

cps.obs.exp.cggcaa<-log(dn31.codon.pair$cggcaa/cps.cgg.caa)
cps.obs.exp.cggcac<-log(dn31.codon.pair$cggcac/cps.cgg.cac)
cps.obs.exp.cggcag<-log(dn31.codon.pair$cggcag/cps.cgg.cag)
cps.obs.exp.cggcat<-log(dn31.codon.pair$cggcat/cps.cgg.cat)
cps.obs.exp.cggcca<-log(dn31.codon.pair$cggcca/cps.cgg.cca)
cps.obs.exp.cggccc<-log(dn31.codon.pair$cggccc/cps.cgg.ccc)
cps.obs.exp.cggccg<-log(dn31.codon.pair$cggccg/cps.cgg.ccg)
cps.obs.exp.cggcct<-log(dn31.codon.pair$cggcct/cps.cgg.cct)
cps.obs.exp.cggcga<-log(dn31.codon.pair$cggcga/cps.cgg.cga)
cps.obs.exp.cggcgc<-log(dn31.codon.pair$cggcgc/cps.cgg.cgc)
cps.obs.exp.cggcgg<-log(dn31.codon.pair$cggcgg/cps.cgg.cgg)
cps.obs.exp.cggcgt<-log(dn31.codon.pair$cggcgt/cps.cgg.cgt)
cps.obs.exp.cggcta<-log(dn31.codon.pair$cggcta/cps.cgg.cta)
cps.obs.exp.cggctc<-log(dn31.codon.pair$cggctc/cps.cgg.ctc)
cps.obs.exp.cggctg<-log(dn31.codon.pair$cggctg/cps.cgg.ctg)
cps.obs.exp.cggctt<-log(dn31.codon.pair$cggctt/cps.cgg.ctt)

cps.obs.exp.cgggaa<-log(dn31.codon.pair$cgggaa/cps.cgg.gaa)
cps.obs.exp.cgggac<-log(dn31.codon.pair$cgggac/cps.cgg.gac)
cps.obs.exp.cgggag<-log(dn31.codon.pair$cgggag/cps.cgg.gag)
cps.obs.exp.cgggat<-log(dn31.codon.pair$cgggat/cps.cgg.gat)
cps.obs.exp.cgggca<-log(dn31.codon.pair$cgggca/cps.cgg.gca)
cps.obs.exp.cgggcc<-log(dn31.codon.pair$cgggcc/cps.cgg.gcc)
cps.obs.exp.cgggcg<-log(dn31.codon.pair$cgggcg/cps.cgg.gcg)
cps.obs.exp.cgggct<-log(dn31.codon.pair$cgggct/cps.cgg.gct)
cps.obs.exp.cgggga<-log(dn31.codon.pair$cgggga/cps.cgg.gga)
cps.obs.exp.cggggc<-log(dn31.codon.pair$cggggc/cps.cgg.ggc)
cps.obs.exp.cggggg<-log(dn31.codon.pair$cggggg/cps.cgg.ggg)
cps.obs.exp.cggggt<-log(dn31.codon.pair$cggggt/cps.cgg.ggt)
cps.obs.exp.cgggta<-log(dn31.codon.pair$cgggta/cps.cgg.gta)
cps.obs.exp.cgggtc<-log(dn31.codon.pair$cgggtc/cps.cgg.gtc)
cps.obs.exp.cgggtg<-log(dn31.codon.pair$cgggtg/cps.cgg.gtg)
cps.obs.exp.cgggtt<-log(dn31.codon.pair$cgggtt/cps.cgg.gtt)

#cps.obs.exp.cggtaa<-log(dn31.codon.pair$cggtaa/cps.cgg.taa)
cps.obs.exp.cggtac<-log(dn31.codon.pair$cggtac/cps.cgg.tac)
#cps.obs.exp.cggtag<-log(dn31.codon.pair$cggtag/cps.cgg.tag)
cps.obs.exp.cggtat<-log(dn31.codon.pair$cggtat/cps.cgg.tat)
cps.obs.exp.cggtca<-log(dn31.codon.pair$cggtca/cps.cgg.tca)
cps.obs.exp.cggtcc<-log(dn31.codon.pair$cggtcc/cps.cgg.tcc)
cps.obs.exp.cggtcg<-log(dn31.codon.pair$cggtcg/cps.cgg.tcg)
cps.obs.exp.cggtct<-log(dn31.codon.pair$cggtct/cps.cgg.tct)
#cps.obs.exp.cggtga<-log(dn31.codon.pair$cggtga/cps.cgg.tga)
cps.obs.exp.cggtgc<-log(dn31.codon.pair$cggtgc/cps.cgg.tgc)
cps.obs.exp.cggtgg<-log(dn31.codon.pair$cggtgg/cps.cgg.tgg)
cps.obs.exp.cggtgt<-log(dn31.codon.pair$cggtgt/cps.cgg.tgt)
cps.obs.exp.cggtta<-log(dn31.codon.pair$cggtta/cps.cgg.tta)
cps.obs.exp.cggttc<-log(dn31.codon.pair$cggttc/cps.cgg.ttc)
cps.obs.exp.cggttg<-log(dn31.codon.pair$cggttg/cps.cgg.ttg)
cps.obs.exp.cggttt<-log(dn31.codon.pair$cggttt/cps.cgg.ttt)








cps.obs.exp.cgtaaa<-log(dn31.codon.pair$cgtaaa/cps.cgt.aaa)
cps.obs.exp.cgtaac<-log(dn31.codon.pair$cgtaac/cps.cgt.aac)
cps.obs.exp.cgtaag<-log(dn31.codon.pair$cgtaag/cps.cgt.aag)
cps.obs.exp.cgtaat<-log(dn31.codon.pair$cgtaat/cps.cgt.aat)
cps.obs.exp.cgtaca<-log(dn31.codon.pair$cgtaca/cps.cgt.aca)
cps.obs.exp.cgtacc<-log(dn31.codon.pair$cgtacc/cps.cgt.acc)
cps.obs.exp.cgtacg<-log(dn31.codon.pair$cgtacg/cps.cgt.acg)
cps.obs.exp.cgtact<-log(dn31.codon.pair$cgtact/cps.cgt.act)
cps.obs.exp.cgtaga<-log(dn31.codon.pair$cgtaga/cps.cgt.aga)
cps.obs.exp.cgtagc<-log(dn31.codon.pair$cgtagc/cps.cgt.agc)
cps.obs.exp.cgtagg<-log(dn31.codon.pair$cgtagg/cps.cgt.agg)
cps.obs.exp.cgtagt<-log(dn31.codon.pair$cgtagt/cps.cgt.agt)
cps.obs.exp.cgtata<-log(dn31.codon.pair$cgtata/cps.cgt.ata)
cps.obs.exp.cgtatc<-log(dn31.codon.pair$cgtatc/cps.cgt.atc)
cps.obs.exp.cgtatg<-log(dn31.codon.pair$cgtatg/cps.cgt.atg)
cps.obs.exp.cgtatt<-log(dn31.codon.pair$cgtatt/cps.cgt.att)

cps.obs.exp.cgtcaa<-log(dn31.codon.pair$cgtcaa/cps.cgt.caa)
cps.obs.exp.cgtcac<-log(dn31.codon.pair$cgtcac/cps.cgt.cac)
cps.obs.exp.cgtcag<-log(dn31.codon.pair$cgtcag/cps.cgt.cag)
cps.obs.exp.cgtcat<-log(dn31.codon.pair$cgtcat/cps.cgt.cat)
cps.obs.exp.cgtcca<-log(dn31.codon.pair$cgtcca/cps.cgt.cca)
cps.obs.exp.cgtccc<-log(dn31.codon.pair$cgtccc/cps.cgt.ccc)
cps.obs.exp.cgtccg<-log(dn31.codon.pair$cgtccg/cps.cgt.ccg)
cps.obs.exp.cgtcct<-log(dn31.codon.pair$cgtcct/cps.cgt.cct)
cps.obs.exp.cgtcga<-log(dn31.codon.pair$cgtcga/cps.cgt.cga)
cps.obs.exp.cgtcgc<-log(dn31.codon.pair$cgtcgc/cps.cgt.cgc)
cps.obs.exp.cgtcgg<-log(dn31.codon.pair$cgtcgg/cps.cgt.cgg)
cps.obs.exp.cgtcgt<-log(dn31.codon.pair$cgtcgt/cps.cgt.cgt)
cps.obs.exp.cgtcta<-log(dn31.codon.pair$cgtcta/cps.cgt.cta)
cps.obs.exp.cgtctc<-log(dn31.codon.pair$cgtctc/cps.cgt.ctc)
cps.obs.exp.cgtctg<-log(dn31.codon.pair$cgtctg/cps.cgt.ctg)
cps.obs.exp.cgtctt<-log(dn31.codon.pair$cgtctt/cps.cgt.ctt)

cps.obs.exp.cgtgaa<-log(dn31.codon.pair$cgtgaa/cps.cgt.gaa)
cps.obs.exp.cgtgac<-log(dn31.codon.pair$cgtgac/cps.cgt.gac)
cps.obs.exp.cgtgag<-log(dn31.codon.pair$cgtgag/cps.cgt.gag)
cps.obs.exp.cgtgat<-log(dn31.codon.pair$cgtgat/cps.cgt.gat)
cps.obs.exp.cgtgca<-log(dn31.codon.pair$cgtgca/cps.cgt.gca)
cps.obs.exp.cgtgcc<-log(dn31.codon.pair$cgtgcc/cps.cgt.gcc)
cps.obs.exp.cgtgcg<-log(dn31.codon.pair$cgtgcg/cps.cgt.gcg)
cps.obs.exp.cgtgct<-log(dn31.codon.pair$cgtgct/cps.cgt.gct)
cps.obs.exp.cgtgga<-log(dn31.codon.pair$cgtgga/cps.cgt.gga)
cps.obs.exp.cgtggc<-log(dn31.codon.pair$cgtggc/cps.cgt.ggc)
cps.obs.exp.cgtggg<-log(dn31.codon.pair$cgtggg/cps.cgt.ggg)
cps.obs.exp.cgtggt<-log(dn31.codon.pair$cgtggt/cps.cgt.ggt)
cps.obs.exp.cgtgta<-log(dn31.codon.pair$cgtgta/cps.cgt.gta)
cps.obs.exp.cgtgtc<-log(dn31.codon.pair$cgtgtc/cps.cgt.gtc)
cps.obs.exp.cgtgtg<-log(dn31.codon.pair$cgtgtg/cps.cgt.gtg)
cps.obs.exp.cgtgtt<-log(dn31.codon.pair$cgtgtt/cps.cgt.gtt)

#cps.obs.exp.cgttaa<-log(dn31.codon.pair$cgttaa/cps.cgt.taa)
cps.obs.exp.cgttac<-log(dn31.codon.pair$cgttac/cps.cgt.tac)
#cps.obs.exp.cgttag<-log(dn31.codon.pair$cgttag/cps.cgt.tag)
cps.obs.exp.cgttat<-log(dn31.codon.pair$cgttat/cps.cgt.tat)
cps.obs.exp.cgttca<-log(dn31.codon.pair$cgttca/cps.cgt.tca)
cps.obs.exp.cgttcc<-log(dn31.codon.pair$cgttcc/cps.cgt.tcc)
cps.obs.exp.cgttcg<-log(dn31.codon.pair$cgttcg/cps.cgt.tcg)
cps.obs.exp.cgttct<-log(dn31.codon.pair$cgttct/cps.cgt.tct)
#cps.obs.exp.cgttga<-log(dn31.codon.pair$cgttga/cps.cgt.tga)
cps.obs.exp.cgttgc<-log(dn31.codon.pair$cgttgc/cps.cgt.tgc)
cps.obs.exp.cgttgg<-log(dn31.codon.pair$cgttgg/cps.cgt.tgg)
cps.obs.exp.cgttgt<-log(dn31.codon.pair$cgttgt/cps.cgt.tgt)
cps.obs.exp.cgttta<-log(dn31.codon.pair$cgttta/cps.cgt.tta)
cps.obs.exp.cgtttc<-log(dn31.codon.pair$cgtttc/cps.cgt.ttc)
cps.obs.exp.cgtttg<-log(dn31.codon.pair$cgtttg/cps.cgt.ttg)
cps.obs.exp.cgtttt<-log(dn31.codon.pair$cgtttt/cps.cgt.ttt)




















cps.obs.exp.ctaaaa<-log(dn31.codon.pair$ctaaaa/cps.cta.aaa)
cps.obs.exp.ctaaac<-log(dn31.codon.pair$ctaaac/cps.cta.aac)
cps.obs.exp.ctaaag<-log(dn31.codon.pair$ctaaag/cps.cta.aag)
cps.obs.exp.ctaaat<-log(dn31.codon.pair$ctaaat/cps.cta.aat)
cps.obs.exp.ctaaca<-log(dn31.codon.pair$ctaaca/cps.cta.aca)
cps.obs.exp.ctaacc<-log(dn31.codon.pair$ctaacc/cps.cta.acc)
cps.obs.exp.ctaacg<-log(dn31.codon.pair$ctaacg/cps.cta.acg)
cps.obs.exp.ctaact<-log(dn31.codon.pair$ctaact/cps.cta.act)
cps.obs.exp.ctaaga<-log(dn31.codon.pair$ctaaga/cps.cta.aga)
cps.obs.exp.ctaagc<-log(dn31.codon.pair$ctaagc/cps.cta.agc)
cps.obs.exp.ctaagg<-log(dn31.codon.pair$ctaagg/cps.cta.agg)
cps.obs.exp.ctaagt<-log(dn31.codon.pair$ctaagt/cps.cta.agt)
cps.obs.exp.ctaata<-log(dn31.codon.pair$ctaata/cps.cta.ata)
cps.obs.exp.ctaatc<-log(dn31.codon.pair$ctaatc/cps.cta.atc)
cps.obs.exp.ctaatg<-log(dn31.codon.pair$ctaatg/cps.cta.atg)
cps.obs.exp.ctaatt<-log(dn31.codon.pair$ctaatt/cps.cta.att)

cps.obs.exp.ctacaa<-log(dn31.codon.pair$ctacaa/cps.cta.caa)
cps.obs.exp.ctacac<-log(dn31.codon.pair$ctacac/cps.cta.cac)
cps.obs.exp.ctacag<-log(dn31.codon.pair$ctacag/cps.cta.cag)
cps.obs.exp.ctacat<-log(dn31.codon.pair$ctacat/cps.cta.cat)
cps.obs.exp.ctacca<-log(dn31.codon.pair$ctacca/cps.cta.cca)
cps.obs.exp.ctaccc<-log(dn31.codon.pair$ctaccc/cps.cta.ccc)
cps.obs.exp.ctaccg<-log(dn31.codon.pair$ctaccg/cps.cta.ccg)
cps.obs.exp.ctacct<-log(dn31.codon.pair$ctacct/cps.cta.cct)
cps.obs.exp.ctacga<-log(dn31.codon.pair$ctacga/cps.cta.cga)
cps.obs.exp.ctacgc<-log(dn31.codon.pair$ctacgc/cps.cta.cgc)
cps.obs.exp.ctacgg<-log(dn31.codon.pair$ctacgg/cps.cta.cgg)
cps.obs.exp.ctacgt<-log(dn31.codon.pair$ctacgt/cps.cta.cgt)
cps.obs.exp.ctacta<-log(dn31.codon.pair$ctacta/cps.cta.cta)
cps.obs.exp.ctactc<-log(dn31.codon.pair$ctactc/cps.cta.ctc)
cps.obs.exp.ctactg<-log(dn31.codon.pair$ctactg/cps.cta.ctg)
cps.obs.exp.ctactt<-log(dn31.codon.pair$ctactt/cps.cta.ctt)

cps.obs.exp.ctagaa<-log(dn31.codon.pair$ctagaa/cps.cta.gaa)
cps.obs.exp.ctagac<-log(dn31.codon.pair$ctagac/cps.cta.gac)
cps.obs.exp.ctagag<-log(dn31.codon.pair$ctagag/cps.cta.gag)
cps.obs.exp.ctagat<-log(dn31.codon.pair$ctagat/cps.cta.gat)
cps.obs.exp.ctagca<-log(dn31.codon.pair$ctagca/cps.cta.gca)
cps.obs.exp.ctagcc<-log(dn31.codon.pair$ctagcc/cps.cta.gcc)
cps.obs.exp.ctagcg<-log(dn31.codon.pair$ctagcg/cps.cta.gcg)
cps.obs.exp.ctagct<-log(dn31.codon.pair$ctagct/cps.cta.gct)
cps.obs.exp.ctagga<-log(dn31.codon.pair$ctagga/cps.cta.gga)
cps.obs.exp.ctaggc<-log(dn31.codon.pair$ctaggc/cps.cta.ggc)
cps.obs.exp.ctaggg<-log(dn31.codon.pair$ctaggg/cps.cta.ggg)
cps.obs.exp.ctaggt<-log(dn31.codon.pair$ctaggt/cps.cta.ggt)
cps.obs.exp.ctagta<-log(dn31.codon.pair$ctagta/cps.cta.gta)
cps.obs.exp.ctagtc<-log(dn31.codon.pair$ctagtc/cps.cta.gtc)
cps.obs.exp.ctagtg<-log(dn31.codon.pair$ctagtg/cps.cta.gtg)
cps.obs.exp.ctagtt<-log(dn31.codon.pair$ctagtt/cps.cta.gtt)

#cps.obs.exp.ctataa<-log(dn31.codon.pair$ctataa/cps.cta.taa)
cps.obs.exp.ctatac<-log(dn31.codon.pair$ctatac/cps.cta.tac)
#cps.obs.exp.ctatag<-log(dn31.codon.pair$ctatag/cps.cta.tag)
cps.obs.exp.ctatat<-log(dn31.codon.pair$ctatat/cps.cta.tat)
cps.obs.exp.ctatca<-log(dn31.codon.pair$ctatca/cps.cta.tca)
cps.obs.exp.ctatcc<-log(dn31.codon.pair$ctatcc/cps.cta.tcc)
cps.obs.exp.ctatcg<-log(dn31.codon.pair$ctatcg/cps.cta.tcg)
cps.obs.exp.ctatct<-log(dn31.codon.pair$ctatct/cps.cta.tct)
#cps.obs.exp.ctatga<-log(dn31.codon.pair$ctatga/cps.cta.tga)
cps.obs.exp.ctatgc<-log(dn31.codon.pair$ctatgc/cps.cta.tgc)
cps.obs.exp.ctatgg<-log(dn31.codon.pair$ctatgg/cps.cta.tgg)
cps.obs.exp.ctatgt<-log(dn31.codon.pair$ctatgt/cps.cta.tgt)
cps.obs.exp.ctatta<-log(dn31.codon.pair$ctatta/cps.cta.tta)
cps.obs.exp.ctattc<-log(dn31.codon.pair$ctattc/cps.cta.ttc)
cps.obs.exp.ctattg<-log(dn31.codon.pair$ctattg/cps.cta.ttg)
cps.obs.exp.ctattt<-log(dn31.codon.pair$ctattt/cps.cta.ttt)









cps.obs.exp.ctcaaa<-log(dn31.codon.pair$ctcaaa/cps.ctc.aaa)
cps.obs.exp.ctcaac<-log(dn31.codon.pair$ctcaac/cps.ctc.aac)
cps.obs.exp.ctcaag<-log(dn31.codon.pair$ctcaag/cps.ctc.aag)
cps.obs.exp.ctcaat<-log(dn31.codon.pair$ctcaat/cps.ctc.aat)
cps.obs.exp.ctcaca<-log(dn31.codon.pair$ctcaca/cps.ctc.aca)
cps.obs.exp.ctcacc<-log(dn31.codon.pair$ctcacc/cps.ctc.acc)
cps.obs.exp.ctcacg<-log(dn31.codon.pair$ctcacg/cps.ctc.acg)
cps.obs.exp.ctcact<-log(dn31.codon.pair$ctcact/cps.ctc.act)
cps.obs.exp.ctcaga<-log(dn31.codon.pair$ctcaga/cps.ctc.aga)
cps.obs.exp.ctcagc<-log(dn31.codon.pair$ctcagc/cps.ctc.agc)
cps.obs.exp.ctcagg<-log(dn31.codon.pair$ctcagg/cps.ctc.agg)
cps.obs.exp.ctcagt<-log(dn31.codon.pair$ctcagt/cps.ctc.agt)
cps.obs.exp.ctcata<-log(dn31.codon.pair$ctcata/cps.ctc.ata)
cps.obs.exp.ctcatc<-log(dn31.codon.pair$ctcatc/cps.ctc.atc)
cps.obs.exp.ctcatg<-log(dn31.codon.pair$ctcatg/cps.ctc.atg)
cps.obs.exp.ctcatt<-log(dn31.codon.pair$ctcatt/cps.ctc.att)

cps.obs.exp.ctccaa<-log(dn31.codon.pair$ctccaa/cps.ctc.caa)
cps.obs.exp.ctccac<-log(dn31.codon.pair$ctccac/cps.ctc.cac)
cps.obs.exp.ctccag<-log(dn31.codon.pair$ctccag/cps.ctc.cag)
cps.obs.exp.ctccat<-log(dn31.codon.pair$ctccat/cps.ctc.cat)
cps.obs.exp.ctccca<-log(dn31.codon.pair$ctccca/cps.ctc.cca)
cps.obs.exp.ctcccc<-log(dn31.codon.pair$ctcccc/cps.ctc.ccc)
cps.obs.exp.ctcccg<-log(dn31.codon.pair$ctcccg/cps.ctc.ccg)
cps.obs.exp.ctccct<-log(dn31.codon.pair$ctccct/cps.ctc.cct)
cps.obs.exp.ctccga<-log(dn31.codon.pair$ctccga/cps.ctc.cga)
cps.obs.exp.ctccgc<-log(dn31.codon.pair$ctccgc/cps.ctc.cgc)
cps.obs.exp.ctccgg<-log(dn31.codon.pair$ctccgg/cps.ctc.cgg)
cps.obs.exp.ctccgt<-log(dn31.codon.pair$ctccgt/cps.ctc.cgt)
cps.obs.exp.ctccta<-log(dn31.codon.pair$ctccta/cps.ctc.cta)
cps.obs.exp.ctcctc<-log(dn31.codon.pair$ctcctc/cps.ctc.ctc)
cps.obs.exp.ctcctg<-log(dn31.codon.pair$ctcctg/cps.ctc.ctg)
cps.obs.exp.ctcctt<-log(dn31.codon.pair$ctcctt/cps.ctc.ctt)

cps.obs.exp.ctcgaa<-log(dn31.codon.pair$ctcgaa/cps.ctc.gaa)
cps.obs.exp.ctcgac<-log(dn31.codon.pair$ctcgac/cps.ctc.gac)
cps.obs.exp.ctcgag<-log(dn31.codon.pair$ctcgag/cps.ctc.gag)
cps.obs.exp.ctcgat<-log(dn31.codon.pair$ctcgat/cps.ctc.gat)
cps.obs.exp.ctcgca<-log(dn31.codon.pair$ctcgca/cps.ctc.gca)
cps.obs.exp.ctcgcc<-log(dn31.codon.pair$ctcgcc/cps.ctc.gcc)
cps.obs.exp.ctcgcg<-log(dn31.codon.pair$ctcgcg/cps.ctc.gcg)
cps.obs.exp.ctcgct<-log(dn31.codon.pair$ctcgct/cps.ctc.gct)
cps.obs.exp.ctcgga<-log(dn31.codon.pair$ctcgga/cps.ctc.gga)
cps.obs.exp.ctcggc<-log(dn31.codon.pair$ctcggc/cps.ctc.ggc)
cps.obs.exp.ctcggg<-log(dn31.codon.pair$ctcggg/cps.ctc.ggg)
cps.obs.exp.ctcggt<-log(dn31.codon.pair$ctcggt/cps.ctc.ggt)
cps.obs.exp.ctcgta<-log(dn31.codon.pair$ctcgta/cps.ctc.gta)
cps.obs.exp.ctcgtc<-log(dn31.codon.pair$ctcgtc/cps.ctc.gtc)
cps.obs.exp.ctcgtg<-log(dn31.codon.pair$ctcgtg/cps.ctc.gtg)
cps.obs.exp.ctcgtt<-log(dn31.codon.pair$ctcgtt/cps.ctc.gtt)

#cps.obs.exp.ctctaa<-log(dn31.codon.pair$ctctaa/cps.ctc.taa)
cps.obs.exp.ctctac<-log(dn31.codon.pair$ctctac/cps.ctc.tac)
#cps.obs.exp.ctctag<-log(dn31.codon.pair$ctctag/cps.ctc.tag)
cps.obs.exp.ctctat<-log(dn31.codon.pair$ctctat/cps.ctc.tat)
cps.obs.exp.ctctca<-log(dn31.codon.pair$ctctca/cps.ctc.tca)
cps.obs.exp.ctctcc<-log(dn31.codon.pair$ctctcc/cps.ctc.tcc)
cps.obs.exp.ctctcg<-log(dn31.codon.pair$ctctcg/cps.ctc.tcg)
cps.obs.exp.ctctct<-log(dn31.codon.pair$ctctct/cps.ctc.tct)
#cps.obs.exp.ctctga<-log(dn31.codon.pair$ctctga/cps.ctc.tga)
cps.obs.exp.ctctgc<-log(dn31.codon.pair$ctctgc/cps.ctc.tgc)
cps.obs.exp.ctctgg<-log(dn31.codon.pair$ctctgg/cps.ctc.tgg)
cps.obs.exp.ctctgt<-log(dn31.codon.pair$ctctgt/cps.ctc.tgt)
cps.obs.exp.ctctta<-log(dn31.codon.pair$ctctta/cps.ctc.tta)
cps.obs.exp.ctcttc<-log(dn31.codon.pair$ctcttc/cps.ctc.ttc)
cps.obs.exp.ctcttg<-log(dn31.codon.pair$ctcttg/cps.ctc.ttg)
cps.obs.exp.ctcttt<-log(dn31.codon.pair$ctcttt/cps.ctc.ttt)









cps.obs.exp.ctgaaa<-log(dn31.codon.pair$ctgaaa/cps.ctg.aaa)
cps.obs.exp.ctgaac<-log(dn31.codon.pair$ctgaac/cps.ctg.aac)
cps.obs.exp.ctgaag<-log(dn31.codon.pair$ctgaag/cps.ctg.aag)
cps.obs.exp.ctgaat<-log(dn31.codon.pair$ctgaat/cps.ctg.aat)
cps.obs.exp.ctgaca<-log(dn31.codon.pair$ctgaca/cps.ctg.aca)
cps.obs.exp.ctgacc<-log(dn31.codon.pair$ctgacc/cps.ctg.acc)
cps.obs.exp.ctgacg<-log(dn31.codon.pair$ctgacg/cps.ctg.acg)
cps.obs.exp.ctgact<-log(dn31.codon.pair$ctgact/cps.ctg.act)
cps.obs.exp.ctgaga<-log(dn31.codon.pair$ctgaga/cps.ctg.aga)
cps.obs.exp.ctgagc<-log(dn31.codon.pair$ctgagc/cps.ctg.agc)
cps.obs.exp.ctgagg<-log(dn31.codon.pair$ctgagg/cps.ctg.agg)
cps.obs.exp.ctgagt<-log(dn31.codon.pair$ctgagt/cps.ctg.agt)
cps.obs.exp.ctgata<-log(dn31.codon.pair$ctgata/cps.ctg.ata)
cps.obs.exp.ctgatc<-log(dn31.codon.pair$ctgatc/cps.ctg.atc)
cps.obs.exp.ctgatg<-log(dn31.codon.pair$ctgatg/cps.ctg.atg)
cps.obs.exp.ctgatt<-log(dn31.codon.pair$ctgatt/cps.ctg.att)

cps.obs.exp.ctgcaa<-log(dn31.codon.pair$ctgcaa/cps.ctg.caa)
cps.obs.exp.ctgcac<-log(dn31.codon.pair$ctgcac/cps.ctg.cac)
cps.obs.exp.ctgcag<-log(dn31.codon.pair$ctgcag/cps.ctg.cag)
cps.obs.exp.ctgcat<-log(dn31.codon.pair$ctgcat/cps.ctg.cat)
cps.obs.exp.ctgcca<-log(dn31.codon.pair$ctgcca/cps.ctg.cca)
cps.obs.exp.ctgccc<-log(dn31.codon.pair$ctgccc/cps.ctg.ccc)
cps.obs.exp.ctgccg<-log(dn31.codon.pair$ctgccg/cps.ctg.ccg)
cps.obs.exp.ctgcct<-log(dn31.codon.pair$ctgcct/cps.ctg.cct)
cps.obs.exp.ctgcga<-log(dn31.codon.pair$ctgcga/cps.ctg.cga)
cps.obs.exp.ctgcgc<-log(dn31.codon.pair$ctgcgc/cps.ctg.cgc)
cps.obs.exp.ctgcgg<-log(dn31.codon.pair$ctgcgg/cps.ctg.cgg)
cps.obs.exp.ctgcgt<-log(dn31.codon.pair$ctgcgt/cps.ctg.cgt)
cps.obs.exp.ctgcta<-log(dn31.codon.pair$ctgcta/cps.ctg.cta)
cps.obs.exp.ctgctc<-log(dn31.codon.pair$ctgctc/cps.ctg.ctc)
cps.obs.exp.ctgctg<-log(dn31.codon.pair$ctgctg/cps.ctg.ctg)
cps.obs.exp.ctgctt<-log(dn31.codon.pair$ctgctt/cps.ctg.ctt)

cps.obs.exp.ctggaa<-log(dn31.codon.pair$ctggaa/cps.ctg.gaa)
cps.obs.exp.ctggac<-log(dn31.codon.pair$ctggac/cps.ctg.gac)
cps.obs.exp.ctggag<-log(dn31.codon.pair$ctggag/cps.ctg.gag)
cps.obs.exp.ctggat<-log(dn31.codon.pair$ctggat/cps.ctg.gat)
cps.obs.exp.ctggca<-log(dn31.codon.pair$ctggca/cps.ctg.gca)
cps.obs.exp.ctggcc<-log(dn31.codon.pair$ctggcc/cps.ctg.gcc)
cps.obs.exp.ctggcg<-log(dn31.codon.pair$ctggcg/cps.ctg.gcg)
cps.obs.exp.ctggct<-log(dn31.codon.pair$ctggct/cps.ctg.gct)
cps.obs.exp.ctggga<-log(dn31.codon.pair$ctggga/cps.ctg.gga)
cps.obs.exp.ctgggc<-log(dn31.codon.pair$ctgggc/cps.ctg.ggc)
cps.obs.exp.ctgggg<-log(dn31.codon.pair$ctgggg/cps.ctg.ggg)
cps.obs.exp.ctgggt<-log(dn31.codon.pair$ctgggt/cps.ctg.ggt)
cps.obs.exp.ctggta<-log(dn31.codon.pair$ctggta/cps.ctg.gta)
cps.obs.exp.ctggtc<-log(dn31.codon.pair$ctggtc/cps.ctg.gtc)
cps.obs.exp.ctggtg<-log(dn31.codon.pair$ctggtg/cps.ctg.gtg)
cps.obs.exp.ctggtt<-log(dn31.codon.pair$ctggtt/cps.ctg.gtt)

#cps.obs.exp.ctgtaa<-log(dn31.codon.pair$ctgtaa/cps.ctg.taa)
cps.obs.exp.ctgtac<-log(dn31.codon.pair$ctgtac/cps.ctg.tac)
#cps.obs.exp.ctgtag<-log(dn31.codon.pair$ctgtag/cps.ctg.tag)
cps.obs.exp.ctgtat<-log(dn31.codon.pair$ctgtat/cps.ctg.tat)
cps.obs.exp.ctgtca<-log(dn31.codon.pair$ctgtca/cps.ctg.tca)
cps.obs.exp.ctgtcc<-log(dn31.codon.pair$ctgtcc/cps.ctg.tcc)
cps.obs.exp.ctgtcg<-log(dn31.codon.pair$ctgtcg/cps.ctg.tcg)
cps.obs.exp.ctgtct<-log(dn31.codon.pair$ctgtct/cps.ctg.tct)
#cps.obs.exp.ctgtga<-log(dn31.codon.pair$ctgtga/cps.ctg.tga)
cps.obs.exp.ctgtgc<-log(dn31.codon.pair$ctgtgc/cps.ctg.tgc)
cps.obs.exp.ctgtgg<-log(dn31.codon.pair$ctgtgg/cps.ctg.tgg)
cps.obs.exp.ctgtgt<-log(dn31.codon.pair$ctgtgt/cps.ctg.tgt)
cps.obs.exp.ctgtta<-log(dn31.codon.pair$ctgtta/cps.ctg.tta)
cps.obs.exp.ctgttc<-log(dn31.codon.pair$ctgttc/cps.ctg.ttc)
cps.obs.exp.ctgttg<-log(dn31.codon.pair$ctgttg/cps.ctg.ttg)
cps.obs.exp.ctgttt<-log(dn31.codon.pair$ctgttt/cps.ctg.ttt)








cps.obs.exp.cttaaa<-log(dn31.codon.pair$cttaaa/cps.ctt.aaa)
cps.obs.exp.cttaac<-log(dn31.codon.pair$cttaac/cps.ctt.aac)
cps.obs.exp.cttaag<-log(dn31.codon.pair$cttaag/cps.ctt.aag)
cps.obs.exp.cttaat<-log(dn31.codon.pair$cttaat/cps.ctt.aat)
cps.obs.exp.cttaca<-log(dn31.codon.pair$cttaca/cps.ctt.aca)
cps.obs.exp.cttacc<-log(dn31.codon.pair$cttacc/cps.ctt.acc)
cps.obs.exp.cttacg<-log(dn31.codon.pair$cttacg/cps.ctt.acg)
cps.obs.exp.cttact<-log(dn31.codon.pair$cttact/cps.ctt.act)
cps.obs.exp.cttaga<-log(dn31.codon.pair$cttaga/cps.ctt.aga)
cps.obs.exp.cttagc<-log(dn31.codon.pair$cttagc/cps.ctt.agc)
cps.obs.exp.cttagg<-log(dn31.codon.pair$cttagg/cps.ctt.agg)
cps.obs.exp.cttagt<-log(dn31.codon.pair$cttagt/cps.ctt.agt)
cps.obs.exp.cttata<-log(dn31.codon.pair$cttata/cps.ctt.ata)
cps.obs.exp.cttatc<-log(dn31.codon.pair$cttatc/cps.ctt.atc)
cps.obs.exp.cttatg<-log(dn31.codon.pair$cttatg/cps.ctt.atg)
cps.obs.exp.cttatt<-log(dn31.codon.pair$cttatt/cps.ctt.att)

cps.obs.exp.cttcaa<-log(dn31.codon.pair$cttcaa/cps.ctt.caa)
cps.obs.exp.cttcac<-log(dn31.codon.pair$cttcac/cps.ctt.cac)
cps.obs.exp.cttcag<-log(dn31.codon.pair$cttcag/cps.ctt.cag)
cps.obs.exp.cttcat<-log(dn31.codon.pair$cttcat/cps.ctt.cat)
cps.obs.exp.cttcca<-log(dn31.codon.pair$cttcca/cps.ctt.cca)
cps.obs.exp.cttccc<-log(dn31.codon.pair$cttccc/cps.ctt.ccc)
cps.obs.exp.cttccg<-log(dn31.codon.pair$cttccg/cps.ctt.ccg)
cps.obs.exp.cttcct<-log(dn31.codon.pair$cttcct/cps.ctt.cct)
cps.obs.exp.cttcga<-log(dn31.codon.pair$cttcga/cps.ctt.cga)
cps.obs.exp.cttcgc<-log(dn31.codon.pair$cttcgc/cps.ctt.cgc)
cps.obs.exp.cttcgg<-log(dn31.codon.pair$cttcgg/cps.ctt.cgg)
cps.obs.exp.cttcgt<-log(dn31.codon.pair$cttcgt/cps.ctt.cgt)
cps.obs.exp.cttcta<-log(dn31.codon.pair$cttcta/cps.ctt.cta)
cps.obs.exp.cttctc<-log(dn31.codon.pair$cttctc/cps.ctt.ctc)
cps.obs.exp.cttctg<-log(dn31.codon.pair$cttctg/cps.ctt.ctg)
cps.obs.exp.cttctt<-log(dn31.codon.pair$cttctt/cps.ctt.ctt)

cps.obs.exp.cttgaa<-log(dn31.codon.pair$cttgaa/cps.ctt.gaa)
cps.obs.exp.cttgac<-log(dn31.codon.pair$cttgac/cps.ctt.gac)
cps.obs.exp.cttgag<-log(dn31.codon.pair$cttgag/cps.ctt.gag)
cps.obs.exp.cttgat<-log(dn31.codon.pair$cttgat/cps.ctt.gat)
cps.obs.exp.cttgca<-log(dn31.codon.pair$cttgca/cps.ctt.gca)
cps.obs.exp.cttgcc<-log(dn31.codon.pair$cttgcc/cps.ctt.gcc)
cps.obs.exp.cttgcg<-log(dn31.codon.pair$cttgcg/cps.ctt.gcg)
cps.obs.exp.cttgct<-log(dn31.codon.pair$cttgct/cps.ctt.gct)
cps.obs.exp.cttgga<-log(dn31.codon.pair$cttgga/cps.ctt.gga)
cps.obs.exp.cttggc<-log(dn31.codon.pair$cttggc/cps.ctt.ggc)
cps.obs.exp.cttggg<-log(dn31.codon.pair$cttggg/cps.ctt.ggg)
cps.obs.exp.cttggt<-log(dn31.codon.pair$cttggt/cps.ctt.ggt)
cps.obs.exp.cttgta<-log(dn31.codon.pair$cttgta/cps.ctt.gta)
cps.obs.exp.cttgtc<-log(dn31.codon.pair$cttgtc/cps.ctt.gtc)
cps.obs.exp.cttgtg<-log(dn31.codon.pair$cttgtg/cps.ctt.gtg)
cps.obs.exp.cttgtt<-log(dn31.codon.pair$cttgtt/cps.ctt.gtt)

#cps.obs.exp.ctttaa<-log(dn31.codon.pair$ctttaa/cps.ctt.taa)
cps.obs.exp.ctttac<-log(dn31.codon.pair$ctttac/cps.ctt.tac)
#cps.obs.exp.ctttag<-log(dn31.codon.pair$ctttag/cps.ctt.tag)
cps.obs.exp.ctttat<-log(dn31.codon.pair$ctttat/cps.ctt.tat)
cps.obs.exp.ctttca<-log(dn31.codon.pair$ctttca/cps.ctt.tca)
cps.obs.exp.ctttcc<-log(dn31.codon.pair$ctttcc/cps.ctt.tcc)
cps.obs.exp.ctttcg<-log(dn31.codon.pair$ctttcg/cps.ctt.tcg)
cps.obs.exp.ctttct<-log(dn31.codon.pair$ctttct/cps.ctt.tct)
#cps.obs.exp.ctttga<-log(dn31.codon.pair$ctttga/cps.ctt.tga)
cps.obs.exp.ctttgc<-log(dn31.codon.pair$ctttgc/cps.ctt.tgc)
cps.obs.exp.ctttgg<-log(dn31.codon.pair$ctttgg/cps.ctt.tgg)
cps.obs.exp.ctttgt<-log(dn31.codon.pair$ctttgt/cps.ctt.tgt)
cps.obs.exp.ctttta<-log(dn31.codon.pair$ctttta/cps.ctt.tta)
cps.obs.exp.cttttc<-log(dn31.codon.pair$cttttc/cps.ctt.ttc)
cps.obs.exp.cttttg<-log(dn31.codon.pair$cttttg/cps.ctt.ttg)
cps.obs.exp.cttttt<-log(dn31.codon.pair$cttttt/cps.ctt.ttt)




















cps.obs.exp.gaaaaa<-log(dn31.codon.pair$gaaaaa/cps.gaa.aaa)
cps.obs.exp.gaaaac<-log(dn31.codon.pair$gaaaac/cps.gaa.aac)
cps.obs.exp.gaaaag<-log(dn31.codon.pair$gaaaag/cps.gaa.aag)
cps.obs.exp.gaaaat<-log(dn31.codon.pair$gaaaat/cps.gaa.aat)
cps.obs.exp.gaaaca<-log(dn31.codon.pair$gaaaca/cps.gaa.aca)
cps.obs.exp.gaaacc<-log(dn31.codon.pair$gaaacc/cps.gaa.acc)
cps.obs.exp.gaaacg<-log(dn31.codon.pair$gaaacg/cps.gaa.acg)
cps.obs.exp.gaaact<-log(dn31.codon.pair$gaaact/cps.gaa.act)
cps.obs.exp.gaaaga<-log(dn31.codon.pair$gaaaga/cps.gaa.aga)
cps.obs.exp.gaaagc<-log(dn31.codon.pair$gaaagc/cps.gaa.agc)
cps.obs.exp.gaaagg<-log(dn31.codon.pair$gaaagg/cps.gaa.agg)
cps.obs.exp.gaaagt<-log(dn31.codon.pair$gaaagt/cps.gaa.agt)
cps.obs.exp.gaaata<-log(dn31.codon.pair$gaaata/cps.gaa.ata)
cps.obs.exp.gaaatc<-log(dn31.codon.pair$gaaatc/cps.gaa.atc)
cps.obs.exp.gaaatg<-log(dn31.codon.pair$gaaatg/cps.gaa.atg)
cps.obs.exp.gaaatt<-log(dn31.codon.pair$gaaatt/cps.gaa.att)

cps.obs.exp.gaacaa<-log(dn31.codon.pair$gaacaa/cps.gaa.caa)
cps.obs.exp.gaacac<-log(dn31.codon.pair$gaacac/cps.gaa.cac)
cps.obs.exp.gaacag<-log(dn31.codon.pair$gaacag/cps.gaa.cag)
cps.obs.exp.gaacat<-log(dn31.codon.pair$gaacat/cps.gaa.cat)
cps.obs.exp.gaacca<-log(dn31.codon.pair$gaacca/cps.gaa.cca)
cps.obs.exp.gaaccc<-log(dn31.codon.pair$gaaccc/cps.gaa.ccc)
cps.obs.exp.gaaccg<-log(dn31.codon.pair$gaaccg/cps.gaa.ccg)
cps.obs.exp.gaacct<-log(dn31.codon.pair$gaacct/cps.gaa.cct)
cps.obs.exp.gaacga<-log(dn31.codon.pair$gaacga/cps.gaa.cga)
cps.obs.exp.gaacgc<-log(dn31.codon.pair$gaacgc/cps.gaa.cgc)
cps.obs.exp.gaacgg<-log(dn31.codon.pair$gaacgg/cps.gaa.cgg)
cps.obs.exp.gaacgt<-log(dn31.codon.pair$gaacgt/cps.gaa.cgt)
cps.obs.exp.gaacta<-log(dn31.codon.pair$gaacta/cps.gaa.cta)
cps.obs.exp.gaactc<-log(dn31.codon.pair$gaactc/cps.gaa.ctc)
cps.obs.exp.gaactg<-log(dn31.codon.pair$gaactg/cps.gaa.ctg)
cps.obs.exp.gaactt<-log(dn31.codon.pair$gaactt/cps.gaa.ctt)

cps.obs.exp.gaagaa<-log(dn31.codon.pair$gaagaa/cps.gaa.gaa)
cps.obs.exp.gaagac<-log(dn31.codon.pair$gaagac/cps.gaa.gac)
cps.obs.exp.gaagag<-log(dn31.codon.pair$gaagag/cps.gaa.gag)
cps.obs.exp.gaagat<-log(dn31.codon.pair$gaagat/cps.gaa.gat)
cps.obs.exp.gaagca<-log(dn31.codon.pair$gaagca/cps.gaa.gca)
cps.obs.exp.gaagcc<-log(dn31.codon.pair$gaagcc/cps.gaa.gcc)
cps.obs.exp.gaagcg<-log(dn31.codon.pair$gaagcg/cps.gaa.gcg)
cps.obs.exp.gaagct<-log(dn31.codon.pair$gaagct/cps.gaa.gct)
cps.obs.exp.gaagga<-log(dn31.codon.pair$gaagga/cps.gaa.gga)
cps.obs.exp.gaaggc<-log(dn31.codon.pair$gaaggc/cps.gaa.ggc)
cps.obs.exp.gaaggg<-log(dn31.codon.pair$gaaggg/cps.gaa.ggg)
cps.obs.exp.gaaggt<-log(dn31.codon.pair$gaaggt/cps.gaa.ggt)
cps.obs.exp.gaagta<-log(dn31.codon.pair$gaagta/cps.gaa.gta)
cps.obs.exp.gaagtc<-log(dn31.codon.pair$gaagtc/cps.gaa.gtc)
cps.obs.exp.gaagtg<-log(dn31.codon.pair$gaagtg/cps.gaa.gtg)
cps.obs.exp.gaagtt<-log(dn31.codon.pair$gaagtt/cps.gaa.gtt)

#cps.obs.exp.gaataa<-log(dn31.codon.pair$gaataa/cps.gaa.taa)
cps.obs.exp.gaatac<-log(dn31.codon.pair$gaatac/cps.gaa.tac)
#cps.obs.exp.gaatag<-log(dn31.codon.pair$gaatag/cps.gaa.tag)
cps.obs.exp.gaatat<-log(dn31.codon.pair$gaatat/cps.gaa.tat)
cps.obs.exp.gaatca<-log(dn31.codon.pair$gaatca/cps.gaa.tca)
cps.obs.exp.gaatcc<-log(dn31.codon.pair$gaatcc/cps.gaa.tcc)
cps.obs.exp.gaatcg<-log(dn31.codon.pair$gaatcg/cps.gaa.tcg)
cps.obs.exp.gaatct<-log(dn31.codon.pair$gaatct/cps.gaa.tct)
#cps.obs.exp.gaatga<-log(dn31.codon.pair$gaatga/cps.gaa.tga)
cps.obs.exp.gaatgc<-log(dn31.codon.pair$gaatgc/cps.gaa.tgc)
cps.obs.exp.gaatgg<-log(dn31.codon.pair$gaatgg/cps.gaa.tgg)
cps.obs.exp.gaatgt<-log(dn31.codon.pair$gaatgt/cps.gaa.tgt)
cps.obs.exp.gaatta<-log(dn31.codon.pair$gaatta/cps.gaa.tta)
cps.obs.exp.gaattc<-log(dn31.codon.pair$gaattc/cps.gaa.ttc)
cps.obs.exp.gaattg<-log(dn31.codon.pair$gaattg/cps.gaa.ttg)
cps.obs.exp.gaattt<-log(dn31.codon.pair$gaattt/cps.gaa.ttt)









cps.obs.exp.gacaaa<-log(dn31.codon.pair$gacaaa/cps.gac.aaa)
cps.obs.exp.gacaac<-log(dn31.codon.pair$gacaac/cps.gac.aac)
cps.obs.exp.gacaag<-log(dn31.codon.pair$gacaag/cps.gac.aag)
cps.obs.exp.gacaat<-log(dn31.codon.pair$gacaat/cps.gac.aat)
cps.obs.exp.gacaca<-log(dn31.codon.pair$gacaca/cps.gac.aca)
cps.obs.exp.gacacc<-log(dn31.codon.pair$gacacc/cps.gac.acc)
cps.obs.exp.gacacg<-log(dn31.codon.pair$gacacg/cps.gac.acg)
cps.obs.exp.gacact<-log(dn31.codon.pair$gacact/cps.gac.act)
cps.obs.exp.gacaga<-log(dn31.codon.pair$gacaga/cps.gac.aga)
cps.obs.exp.gacagc<-log(dn31.codon.pair$gacagc/cps.gac.agc)
cps.obs.exp.gacagg<-log(dn31.codon.pair$gacagg/cps.gac.agg)
cps.obs.exp.gacagt<-log(dn31.codon.pair$gacagt/cps.gac.agt)
cps.obs.exp.gacata<-log(dn31.codon.pair$gacata/cps.gac.ata)
cps.obs.exp.gacatc<-log(dn31.codon.pair$gacatc/cps.gac.atc)
cps.obs.exp.gacatg<-log(dn31.codon.pair$gacatg/cps.gac.atg)
cps.obs.exp.gacatt<-log(dn31.codon.pair$gacatt/cps.gac.att)

cps.obs.exp.gaccaa<-log(dn31.codon.pair$gaccaa/cps.gac.caa)
cps.obs.exp.gaccac<-log(dn31.codon.pair$gaccac/cps.gac.cac)
cps.obs.exp.gaccag<-log(dn31.codon.pair$gaccag/cps.gac.cag)
cps.obs.exp.gaccat<-log(dn31.codon.pair$gaccat/cps.gac.cat)
cps.obs.exp.gaccca<-log(dn31.codon.pair$gaccca/cps.gac.cca)
cps.obs.exp.gacccc<-log(dn31.codon.pair$gacccc/cps.gac.ccc)
cps.obs.exp.gacccg<-log(dn31.codon.pair$gacccg/cps.gac.ccg)
cps.obs.exp.gaccct<-log(dn31.codon.pair$gaccct/cps.gac.cct)
cps.obs.exp.gaccga<-log(dn31.codon.pair$gaccga/cps.gac.cga)
cps.obs.exp.gaccgc<-log(dn31.codon.pair$gaccgc/cps.gac.cgc)
cps.obs.exp.gaccgg<-log(dn31.codon.pair$gaccgg/cps.gac.cgg)
cps.obs.exp.gaccgt<-log(dn31.codon.pair$gaccgt/cps.gac.cgt)
cps.obs.exp.gaccta<-log(dn31.codon.pair$gaccta/cps.gac.cta)
cps.obs.exp.gacctc<-log(dn31.codon.pair$gacctc/cps.gac.ctc)
cps.obs.exp.gacctg<-log(dn31.codon.pair$gacctg/cps.gac.ctg)
cps.obs.exp.gacctt<-log(dn31.codon.pair$gacctt/cps.gac.ctt)

cps.obs.exp.gacgaa<-log(dn31.codon.pair$gacgaa/cps.gac.gaa)
cps.obs.exp.gacgac<-log(dn31.codon.pair$gacgac/cps.gac.gac)
cps.obs.exp.gacgag<-log(dn31.codon.pair$gacgag/cps.gac.gag)
cps.obs.exp.gacgat<-log(dn31.codon.pair$gacgat/cps.gac.gat)
cps.obs.exp.gacgca<-log(dn31.codon.pair$gacgca/cps.gac.gca)
cps.obs.exp.gacgcc<-log(dn31.codon.pair$gacgcc/cps.gac.gcc)
cps.obs.exp.gacgcg<-log(dn31.codon.pair$gacgcg/cps.gac.gcg)
cps.obs.exp.gacgct<-log(dn31.codon.pair$gacgct/cps.gac.gct)
cps.obs.exp.gacgga<-log(dn31.codon.pair$gacgga/cps.gac.gga)
cps.obs.exp.gacggc<-log(dn31.codon.pair$gacggc/cps.gac.ggc)
cps.obs.exp.gacggg<-log(dn31.codon.pair$gacggg/cps.gac.ggg)
cps.obs.exp.gacggt<-log(dn31.codon.pair$gacggt/cps.gac.ggt)
cps.obs.exp.gacgta<-log(dn31.codon.pair$gacgta/cps.gac.gta)
cps.obs.exp.gacgtc<-log(dn31.codon.pair$gacgtc/cps.gac.gtc)
cps.obs.exp.gacgtg<-log(dn31.codon.pair$gacgtg/cps.gac.gtg)
cps.obs.exp.gacgtt<-log(dn31.codon.pair$gacgtt/cps.gac.gtt)

#cps.obs.exp.gactaa<-log(dn31.codon.pair$gactaa/cps.gac.taa)
cps.obs.exp.gactac<-log(dn31.codon.pair$gactac/cps.gac.tac)
#cps.obs.exp.gactag<-log(dn31.codon.pair$gactag/cps.gac.tag)
cps.obs.exp.gactat<-log(dn31.codon.pair$gactat/cps.gac.tat)
cps.obs.exp.gactca<-log(dn31.codon.pair$gactca/cps.gac.tca)
cps.obs.exp.gactcc<-log(dn31.codon.pair$gactcc/cps.gac.tcc)
cps.obs.exp.gactcg<-log(dn31.codon.pair$gactcg/cps.gac.tcg)
cps.obs.exp.gactct<-log(dn31.codon.pair$gactct/cps.gac.tct)
#cps.obs.exp.gactga<-log(dn31.codon.pair$gactga/cps.gac.tga)
cps.obs.exp.gactgc<-log(dn31.codon.pair$gactgc/cps.gac.tgc)
cps.obs.exp.gactgg<-log(dn31.codon.pair$gactgg/cps.gac.tgg)
cps.obs.exp.gactgt<-log(dn31.codon.pair$gactgt/cps.gac.tgt)
cps.obs.exp.gactta<-log(dn31.codon.pair$gactta/cps.gac.tta)
cps.obs.exp.gacttc<-log(dn31.codon.pair$gacttc/cps.gac.ttc)
cps.obs.exp.gacttg<-log(dn31.codon.pair$gacttg/cps.gac.ttg)
cps.obs.exp.gacttt<-log(dn31.codon.pair$gacttt/cps.gac.ttt)










cps.obs.exp.gagaaa<-log(dn31.codon.pair$gagaaa/cps.gag.aaa)
cps.obs.exp.gagaac<-log(dn31.codon.pair$gagaac/cps.gag.aac)
cps.obs.exp.gagaag<-log(dn31.codon.pair$gagaag/cps.gag.aag)
cps.obs.exp.gagaat<-log(dn31.codon.pair$gagaat/cps.gag.aat)
cps.obs.exp.gagaca<-log(dn31.codon.pair$gagaca/cps.gag.aca)
cps.obs.exp.gagacc<-log(dn31.codon.pair$gagacc/cps.gag.acc)
cps.obs.exp.gagacg<-log(dn31.codon.pair$gagacg/cps.gag.acg)
cps.obs.exp.gagact<-log(dn31.codon.pair$gagact/cps.gag.act)
cps.obs.exp.gagaga<-log(dn31.codon.pair$gagaga/cps.gag.aga)
cps.obs.exp.gagagc<-log(dn31.codon.pair$gagagc/cps.gag.agc)
cps.obs.exp.gagagg<-log(dn31.codon.pair$gagagg/cps.gag.agg)
cps.obs.exp.gagagt<-log(dn31.codon.pair$gagagt/cps.gag.agt)
cps.obs.exp.gagata<-log(dn31.codon.pair$gagata/cps.gag.ata)
cps.obs.exp.gagatc<-log(dn31.codon.pair$gagatc/cps.gag.atc)
cps.obs.exp.gagatg<-log(dn31.codon.pair$gagatg/cps.gag.atg)
cps.obs.exp.gagatt<-log(dn31.codon.pair$gagatt/cps.gag.att)

cps.obs.exp.gagcaa<-log(dn31.codon.pair$gagcaa/cps.gag.caa)
cps.obs.exp.gagcac<-log(dn31.codon.pair$gagcac/cps.gag.cac)
cps.obs.exp.gagcag<-log(dn31.codon.pair$gagcag/cps.gag.cag)
cps.obs.exp.gagcat<-log(dn31.codon.pair$gagcat/cps.gag.cat)
cps.obs.exp.gagcca<-log(dn31.codon.pair$gagcca/cps.gag.cca)
cps.obs.exp.gagccc<-log(dn31.codon.pair$gagccc/cps.gag.ccc)
cps.obs.exp.gagccg<-log(dn31.codon.pair$gagccg/cps.gag.ccg)
cps.obs.exp.gagcct<-log(dn31.codon.pair$gagcct/cps.gag.cct)
cps.obs.exp.gagcga<-log(dn31.codon.pair$gagcga/cps.gag.cga)
cps.obs.exp.gagcgc<-log(dn31.codon.pair$gagcgc/cps.gag.cgc)
cps.obs.exp.gagcgg<-log(dn31.codon.pair$gagcgg/cps.gag.cgg)
cps.obs.exp.gagcgt<-log(dn31.codon.pair$gagcgt/cps.gag.cgt)
cps.obs.exp.gagcta<-log(dn31.codon.pair$gagcta/cps.gag.cta)
cps.obs.exp.gagctc<-log(dn31.codon.pair$gagctc/cps.gag.ctc)
cps.obs.exp.gagctg<-log(dn31.codon.pair$gagctg/cps.gag.ctg)
cps.obs.exp.gagctt<-log(dn31.codon.pair$gagctt/cps.gag.ctt)

cps.obs.exp.gaggaa<-log(dn31.codon.pair$gaggaa/cps.gag.gaa)
cps.obs.exp.gaggac<-log(dn31.codon.pair$gaggac/cps.gag.gac)
cps.obs.exp.gaggag<-log(dn31.codon.pair$gaggag/cps.gag.gag)
cps.obs.exp.gaggat<-log(dn31.codon.pair$gaggat/cps.gag.gat)
cps.obs.exp.gaggca<-log(dn31.codon.pair$gaggca/cps.gag.gca)
cps.obs.exp.gaggcc<-log(dn31.codon.pair$gaggcc/cps.gag.gcc)
cps.obs.exp.gaggcg<-log(dn31.codon.pair$gaggcg/cps.gag.gcg)
cps.obs.exp.gaggct<-log(dn31.codon.pair$gaggct/cps.gag.gct)
cps.obs.exp.gaggga<-log(dn31.codon.pair$gaggga/cps.gag.gga)
cps.obs.exp.gagggc<-log(dn31.codon.pair$gagggc/cps.gag.ggc)
cps.obs.exp.gagggg<-log(dn31.codon.pair$gagggg/cps.gag.ggg)
cps.obs.exp.gagggt<-log(dn31.codon.pair$gagggt/cps.gag.ggt)
cps.obs.exp.gaggta<-log(dn31.codon.pair$gaggta/cps.gag.gta)
cps.obs.exp.gaggtc<-log(dn31.codon.pair$gaggtc/cps.gag.gtc)
cps.obs.exp.gaggtg<-log(dn31.codon.pair$gaggtg/cps.gag.gtg)
cps.obs.exp.gaggtt<-log(dn31.codon.pair$gaggtt/cps.gag.gtt)

#cps.obs.exp.gagtaa<-log(dn31.codon.pair$gagtaa/cps.gag.taa)
cps.obs.exp.gagtac<-log(dn31.codon.pair$gagtac/cps.gag.tac)
#cps.obs.exp.gagtag<-log(dn31.codon.pair$gagtag/cps.gag.tag)
cps.obs.exp.gagtat<-log(dn31.codon.pair$gagtat/cps.gag.tat)
cps.obs.exp.gagtca<-log(dn31.codon.pair$gagtca/cps.gag.tca)
cps.obs.exp.gagtcc<-log(dn31.codon.pair$gagtcc/cps.gag.tcc)
cps.obs.exp.gagtcg<-log(dn31.codon.pair$gagtcg/cps.gag.tcg)
cps.obs.exp.gagtct<-log(dn31.codon.pair$gagtct/cps.gag.tct)
#cps.obs.exp.gagtga<-log(dn31.codon.pair$gagtga/cps.gag.tga)
cps.obs.exp.gagtgc<-log(dn31.codon.pair$gagtgc/cps.gag.tgc)
cps.obs.exp.gagtgg<-log(dn31.codon.pair$gagtgg/cps.gag.tgg)
cps.obs.exp.gagtgt<-log(dn31.codon.pair$gagtgt/cps.gag.tgt)
cps.obs.exp.gagtta<-log(dn31.codon.pair$gagtta/cps.gag.tta)
cps.obs.exp.gagttc<-log(dn31.codon.pair$gagttc/cps.gag.ttc)
cps.obs.exp.gagttg<-log(dn31.codon.pair$gagttg/cps.gag.ttg)
cps.obs.exp.gagttt<-log(dn31.codon.pair$gagttt/cps.gag.ttt)









cps.obs.exp.gataaa<-log(dn31.codon.pair$gataaa/cps.gat.aaa)
cps.obs.exp.gataac<-log(dn31.codon.pair$gataac/cps.gat.aac)
cps.obs.exp.gataag<-log(dn31.codon.pair$gataag/cps.gat.aag)
cps.obs.exp.gataat<-log(dn31.codon.pair$gataat/cps.gat.aat)
cps.obs.exp.gataca<-log(dn31.codon.pair$gataca/cps.gat.aca)
cps.obs.exp.gatacc<-log(dn31.codon.pair$gatacc/cps.gat.acc)
cps.obs.exp.gatacg<-log(dn31.codon.pair$gatacg/cps.gat.acg)
cps.obs.exp.gatact<-log(dn31.codon.pair$gatact/cps.gat.act)
cps.obs.exp.gataga<-log(dn31.codon.pair$gataga/cps.gat.aga)
cps.obs.exp.gatagc<-log(dn31.codon.pair$gatagc/cps.gat.agc)
cps.obs.exp.gatagg<-log(dn31.codon.pair$gatagg/cps.gat.agg)
cps.obs.exp.gatagt<-log(dn31.codon.pair$gatagt/cps.gat.agt)
cps.obs.exp.gatata<-log(dn31.codon.pair$gatata/cps.gat.ata)
cps.obs.exp.gatatc<-log(dn31.codon.pair$gatatc/cps.gat.atc)
cps.obs.exp.gatatg<-log(dn31.codon.pair$gatatg/cps.gat.atg)
cps.obs.exp.gatatt<-log(dn31.codon.pair$gatatt/cps.gat.att)

cps.obs.exp.gatcaa<-log(dn31.codon.pair$gatcaa/cps.gat.caa)
cps.obs.exp.gatcac<-log(dn31.codon.pair$gatcac/cps.gat.cac)
cps.obs.exp.gatcag<-log(dn31.codon.pair$gatcag/cps.gat.cag)
cps.obs.exp.gatcat<-log(dn31.codon.pair$gatcat/cps.gat.cat)
cps.obs.exp.gatcca<-log(dn31.codon.pair$gatcca/cps.gat.cca)
cps.obs.exp.gatccc<-log(dn31.codon.pair$gatccc/cps.gat.ccc)
cps.obs.exp.gatccg<-log(dn31.codon.pair$gatccg/cps.gat.ccg)
cps.obs.exp.gatcct<-log(dn31.codon.pair$gatcct/cps.gat.cct)
cps.obs.exp.gatcga<-log(dn31.codon.pair$gatcga/cps.gat.cga)
cps.obs.exp.gatcgc<-log(dn31.codon.pair$gatcgc/cps.gat.cgc)
cps.obs.exp.gatcgg<-log(dn31.codon.pair$gatcgg/cps.gat.cgg)
cps.obs.exp.gatcgt<-log(dn31.codon.pair$gatcgt/cps.gat.cgt)
cps.obs.exp.gatcta<-log(dn31.codon.pair$gatcta/cps.gat.cta)
cps.obs.exp.gatctc<-log(dn31.codon.pair$gatctc/cps.gat.ctc)
cps.obs.exp.gatctg<-log(dn31.codon.pair$gatctg/cps.gat.ctg)
cps.obs.exp.gatctt<-log(dn31.codon.pair$gatctt/cps.gat.ctt)

cps.obs.exp.gatgaa<-log(dn31.codon.pair$gatgaa/cps.gat.gaa)
cps.obs.exp.gatgac<-log(dn31.codon.pair$gatgac/cps.gat.gac)
cps.obs.exp.gatgag<-log(dn31.codon.pair$gatgag/cps.gat.gag)
cps.obs.exp.gatgat<-log(dn31.codon.pair$gatgat/cps.gat.gat)
cps.obs.exp.gatgca<-log(dn31.codon.pair$gatgca/cps.gat.gca)
cps.obs.exp.gatgcc<-log(dn31.codon.pair$gatgcc/cps.gat.gcc)
cps.obs.exp.gatgcg<-log(dn31.codon.pair$gatgcg/cps.gat.gcg)
cps.obs.exp.gatgct<-log(dn31.codon.pair$gatgct/cps.gat.gct)
cps.obs.exp.gatgga<-log(dn31.codon.pair$gatgga/cps.gat.gga)
cps.obs.exp.gatggc<-log(dn31.codon.pair$gatggc/cps.gat.ggc)
cps.obs.exp.gatggg<-log(dn31.codon.pair$gatggg/cps.gat.ggg)
cps.obs.exp.gatggt<-log(dn31.codon.pair$gatggt/cps.gat.ggt)
cps.obs.exp.gatgta<-log(dn31.codon.pair$gatgta/cps.gat.gta)
cps.obs.exp.gatgtc<-log(dn31.codon.pair$gatgtc/cps.gat.gtc)
cps.obs.exp.gatgtg<-log(dn31.codon.pair$gatgtg/cps.gat.gtg)
cps.obs.exp.gatgtt<-log(dn31.codon.pair$gatgtt/cps.gat.gtt)

#cps.obs.exp.gattaa<-log(dn31.codon.pair$gattaa/cps.gat.taa)
cps.obs.exp.gattac<-log(dn31.codon.pair$gattac/cps.gat.tac)
#cps.obs.exp.gattag<-log(dn31.codon.pair$gattag/cps.gat.tag)
cps.obs.exp.gattat<-log(dn31.codon.pair$gattat/cps.gat.tat)
cps.obs.exp.gattca<-log(dn31.codon.pair$gattca/cps.gat.tca)
cps.obs.exp.gattcc<-log(dn31.codon.pair$gattcc/cps.gat.tcc)
cps.obs.exp.gattcg<-log(dn31.codon.pair$gattcg/cps.gat.tcg)
cps.obs.exp.gattct<-log(dn31.codon.pair$gattct/cps.gat.tct)
#cps.obs.exp.gattga<-log(dn31.codon.pair$gattga/cps.gat.tga)
cps.obs.exp.gattgc<-log(dn31.codon.pair$gattgc/cps.gat.tgc)
cps.obs.exp.gattgg<-log(dn31.codon.pair$gattgg/cps.gat.tgg)
cps.obs.exp.gattgt<-log(dn31.codon.pair$gattgt/cps.gat.tgt)
cps.obs.exp.gattta<-log(dn31.codon.pair$gattta/cps.gat.tta)
cps.obs.exp.gatttc<-log(dn31.codon.pair$gatttc/cps.gat.ttc)
cps.obs.exp.gatttg<-log(dn31.codon.pair$gatttg/cps.gat.ttg)
cps.obs.exp.gatttt<-log(dn31.codon.pair$gatttt/cps.gat.ttt)

















cps.obs.exp.gcaaaa<-log(dn31.codon.pair$gcaaaa/cps.gca.aaa)
cps.obs.exp.gcaaac<-log(dn31.codon.pair$gcaaac/cps.gca.aac)
cps.obs.exp.gcaaag<-log(dn31.codon.pair$gcaaag/cps.gca.aag)
cps.obs.exp.gcaaat<-log(dn31.codon.pair$gcaaat/cps.gca.aat)
cps.obs.exp.gcaaca<-log(dn31.codon.pair$gcaaca/cps.gca.aca)
cps.obs.exp.gcaacc<-log(dn31.codon.pair$gcaacc/cps.gca.acc)
cps.obs.exp.gcaacg<-log(dn31.codon.pair$gcaacg/cps.gca.acg)
cps.obs.exp.gcaact<-log(dn31.codon.pair$gcaact/cps.gca.act)
cps.obs.exp.gcaaga<-log(dn31.codon.pair$gcaaga/cps.gca.aga)
cps.obs.exp.gcaagc<-log(dn31.codon.pair$gcaagc/cps.gca.agc)
cps.obs.exp.gcaagg<-log(dn31.codon.pair$gcaagg/cps.gca.agg)
cps.obs.exp.gcaagt<-log(dn31.codon.pair$gcaagt/cps.gca.agt)
cps.obs.exp.gcaata<-log(dn31.codon.pair$gcaata/cps.gca.ata)
cps.obs.exp.gcaatc<-log(dn31.codon.pair$gcaatc/cps.gca.atc)
cps.obs.exp.gcaatg<-log(dn31.codon.pair$gcaatg/cps.gca.atg)
cps.obs.exp.gcaatt<-log(dn31.codon.pair$gcaatt/cps.gca.att)

cps.obs.exp.gcacaa<-log(dn31.codon.pair$gcacaa/cps.gca.caa)
cps.obs.exp.gcacac<-log(dn31.codon.pair$gcacac/cps.gca.cac)
cps.obs.exp.gcacag<-log(dn31.codon.pair$gcacag/cps.gca.cag)
cps.obs.exp.gcacat<-log(dn31.codon.pair$gcacat/cps.gca.cat)
cps.obs.exp.gcacca<-log(dn31.codon.pair$gcacca/cps.gca.cca)
cps.obs.exp.gcaccc<-log(dn31.codon.pair$gcaccc/cps.gca.ccc)
cps.obs.exp.gcaccg<-log(dn31.codon.pair$gcaccg/cps.gca.ccg)
cps.obs.exp.gcacct<-log(dn31.codon.pair$gcacct/cps.gca.cct)
cps.obs.exp.gcacga<-log(dn31.codon.pair$gcacga/cps.gca.cga)
cps.obs.exp.gcacgc<-log(dn31.codon.pair$gcacgc/cps.gca.cgc)
cps.obs.exp.gcacgg<-log(dn31.codon.pair$gcacgg/cps.gca.cgg)
cps.obs.exp.gcacgt<-log(dn31.codon.pair$gcacgt/cps.gca.cgt)
cps.obs.exp.gcacta<-log(dn31.codon.pair$gcacta/cps.gca.cta)
cps.obs.exp.gcactc<-log(dn31.codon.pair$gcactc/cps.gca.ctc)
cps.obs.exp.gcactg<-log(dn31.codon.pair$gcactg/cps.gca.ctg)
cps.obs.exp.gcactt<-log(dn31.codon.pair$gcactt/cps.gca.ctt)

cps.obs.exp.gcagaa<-log(dn31.codon.pair$gcagaa/cps.gca.gaa)
cps.obs.exp.gcagac<-log(dn31.codon.pair$gcagac/cps.gca.gac)
cps.obs.exp.gcagag<-log(dn31.codon.pair$gcagag/cps.gca.gag)
cps.obs.exp.gcagat<-log(dn31.codon.pair$gcagat/cps.gca.gat)
cps.obs.exp.gcagca<-log(dn31.codon.pair$gcagca/cps.gca.gca)
cps.obs.exp.gcagcc<-log(dn31.codon.pair$gcagcc/cps.gca.gcc)
cps.obs.exp.gcagcg<-log(dn31.codon.pair$gcagcg/cps.gca.gcg)
cps.obs.exp.gcagct<-log(dn31.codon.pair$gcagct/cps.gca.gct)
cps.obs.exp.gcagga<-log(dn31.codon.pair$gcagga/cps.gca.gga)
cps.obs.exp.gcaggc<-log(dn31.codon.pair$gcaggc/cps.gca.ggc)
cps.obs.exp.gcaggg<-log(dn31.codon.pair$gcaggg/cps.gca.ggg)
cps.obs.exp.gcaggt<-log(dn31.codon.pair$gcaggt/cps.gca.ggt)
cps.obs.exp.gcagta<-log(dn31.codon.pair$gcagta/cps.gca.gta)
cps.obs.exp.gcagtc<-log(dn31.codon.pair$gcagtc/cps.gca.gtc)
cps.obs.exp.gcagtg<-log(dn31.codon.pair$gcagtg/cps.gca.gtg)
cps.obs.exp.gcagtt<-log(dn31.codon.pair$gcagtt/cps.gca.gtt)

#cps.obs.exp.gcataa<-log(dn31.codon.pair$gcataa/cps.gca.taa)
cps.obs.exp.gcatac<-log(dn31.codon.pair$gcatac/cps.gca.tac)
#cps.obs.exp.gcatag<-log(dn31.codon.pair$gcatag/cps.gca.tag)
cps.obs.exp.gcatat<-log(dn31.codon.pair$gcatat/cps.gca.tat)
cps.obs.exp.gcatca<-log(dn31.codon.pair$gcatca/cps.gca.tca)
cps.obs.exp.gcatcc<-log(dn31.codon.pair$gcatcc/cps.gca.tcc)
cps.obs.exp.gcatcg<-log(dn31.codon.pair$gcatcg/cps.gca.tcg)
cps.obs.exp.gcatct<-log(dn31.codon.pair$gcatct/cps.gca.tct)
#cps.obs.exp.gcatga<-log(dn31.codon.pair$gcatga/cps.gca.tga)
cps.obs.exp.gcatgc<-log(dn31.codon.pair$gcatgc/cps.gca.tgc)
cps.obs.exp.gcatgg<-log(dn31.codon.pair$gcatgg/cps.gca.tgg)
cps.obs.exp.gcatgt<-log(dn31.codon.pair$gcatgt/cps.gca.tgt)
cps.obs.exp.gcatta<-log(dn31.codon.pair$gcatta/cps.gca.tta)
cps.obs.exp.gcattc<-log(dn31.codon.pair$gcattc/cps.gca.ttc)
cps.obs.exp.gcattg<-log(dn31.codon.pair$gcattg/cps.gca.ttg)
cps.obs.exp.gcattt<-log(dn31.codon.pair$gcattt/cps.gca.ttt)









cps.obs.exp.gccaaa<-log(dn31.codon.pair$gccaaa/cps.gcc.aaa)
cps.obs.exp.gccaac<-log(dn31.codon.pair$gccaac/cps.gcc.aac)
cps.obs.exp.gccaag<-log(dn31.codon.pair$gccaag/cps.gcc.aag)
cps.obs.exp.gccaat<-log(dn31.codon.pair$gccaat/cps.gcc.aat)
cps.obs.exp.gccaca<-log(dn31.codon.pair$gccaca/cps.gcc.aca)
cps.obs.exp.gccacc<-log(dn31.codon.pair$gccacc/cps.gcc.acc)
cps.obs.exp.gccacg<-log(dn31.codon.pair$gccacg/cps.gcc.acg)
cps.obs.exp.gccact<-log(dn31.codon.pair$gccact/cps.gcc.act)
cps.obs.exp.gccaga<-log(dn31.codon.pair$gccaga/cps.gcc.aga)
cps.obs.exp.gccagc<-log(dn31.codon.pair$gccagc/cps.gcc.agc)
cps.obs.exp.gccagg<-log(dn31.codon.pair$gccagg/cps.gcc.agg)
cps.obs.exp.gccagt<-log(dn31.codon.pair$gccagt/cps.gcc.agt)
cps.obs.exp.gccata<-log(dn31.codon.pair$gccata/cps.gcc.ata)
cps.obs.exp.gccatc<-log(dn31.codon.pair$gccatc/cps.gcc.atc)
cps.obs.exp.gccatg<-log(dn31.codon.pair$gccatg/cps.gcc.atg)
cps.obs.exp.gccatt<-log(dn31.codon.pair$gccatt/cps.gcc.att)

cps.obs.exp.gcccaa<-log(dn31.codon.pair$gcccaa/cps.gcc.caa)
cps.obs.exp.gcccac<-log(dn31.codon.pair$gcccac/cps.gcc.cac)
cps.obs.exp.gcccag<-log(dn31.codon.pair$gcccag/cps.gcc.cag)
cps.obs.exp.gcccat<-log(dn31.codon.pair$gcccat/cps.gcc.cat)
cps.obs.exp.gcccca<-log(dn31.codon.pair$gcccca/cps.gcc.cca)
cps.obs.exp.gccccc<-log(dn31.codon.pair$gccccc/cps.gcc.ccc)
cps.obs.exp.gccccg<-log(dn31.codon.pair$gccccg/cps.gcc.ccg)
cps.obs.exp.gcccct<-log(dn31.codon.pair$gcccct/cps.gcc.cct)
cps.obs.exp.gcccga<-log(dn31.codon.pair$gcccga/cps.gcc.cga)
cps.obs.exp.gcccgc<-log(dn31.codon.pair$gcccgc/cps.gcc.cgc)
cps.obs.exp.gcccgg<-log(dn31.codon.pair$gcccgg/cps.gcc.cgg)
cps.obs.exp.gcccgt<-log(dn31.codon.pair$gcccgt/cps.gcc.cgt)
cps.obs.exp.gcccta<-log(dn31.codon.pair$gcccta/cps.gcc.cta)
cps.obs.exp.gccctc<-log(dn31.codon.pair$gccctc/cps.gcc.ctc)
cps.obs.exp.gccctg<-log(dn31.codon.pair$gccctg/cps.gcc.ctg)
cps.obs.exp.gccctt<-log(dn31.codon.pair$gccctt/cps.gcc.ctt)

cps.obs.exp.gccgaa<-log(dn31.codon.pair$gccgaa/cps.gcc.gaa)
cps.obs.exp.gccgac<-log(dn31.codon.pair$gccgac/cps.gcc.gac)
cps.obs.exp.gccgag<-log(dn31.codon.pair$gccgag/cps.gcc.gag)
cps.obs.exp.gccgat<-log(dn31.codon.pair$gccgat/cps.gcc.gat)
cps.obs.exp.gccgca<-log(dn31.codon.pair$gccgca/cps.gcc.gca)
cps.obs.exp.gccgcc<-log(dn31.codon.pair$gccgcc/cps.gcc.gcc)
cps.obs.exp.gccgcg<-log(dn31.codon.pair$gccgcg/cps.gcc.gcg)
cps.obs.exp.gccgct<-log(dn31.codon.pair$gccgct/cps.gcc.gct)
cps.obs.exp.gccgga<-log(dn31.codon.pair$gccgga/cps.gcc.gga)
cps.obs.exp.gccggc<-log(dn31.codon.pair$gccggc/cps.gcc.ggc)
cps.obs.exp.gccggg<-log(dn31.codon.pair$gccggg/cps.gcc.ggg)
cps.obs.exp.gccggt<-log(dn31.codon.pair$gccggt/cps.gcc.ggt)
cps.obs.exp.gccgta<-log(dn31.codon.pair$gccgta/cps.gcc.gta)
cps.obs.exp.gccgtc<-log(dn31.codon.pair$gccgtc/cps.gcc.gtc)
cps.obs.exp.gccgtg<-log(dn31.codon.pair$gccgtg/cps.gcc.gtg)
cps.obs.exp.gccgtt<-log(dn31.codon.pair$gccgtt/cps.gcc.gtt)

#cps.obs.exp.gcctaa<-log(dn31.codon.pair$gcctaa/cps.gcc.taa)
cps.obs.exp.gcctac<-log(dn31.codon.pair$gcctac/cps.gcc.tac)
#cps.obs.exp.gcctag<-log(dn31.codon.pair$gcctag/cps.gcc.tag)
cps.obs.exp.gcctat<-log(dn31.codon.pair$gcctat/cps.gcc.tat)
cps.obs.exp.gcctca<-log(dn31.codon.pair$gcctca/cps.gcc.tca)
cps.obs.exp.gcctcc<-log(dn31.codon.pair$gcctcc/cps.gcc.tcc)
cps.obs.exp.gcctcg<-log(dn31.codon.pair$gcctcg/cps.gcc.tcg)
cps.obs.exp.gcctct<-log(dn31.codon.pair$gcctct/cps.gcc.tct)
#cps.obs.exp.gcctga<-log(dn31.codon.pair$gcctga/cps.gcc.tga)
cps.obs.exp.gcctgc<-log(dn31.codon.pair$gcctgc/cps.gcc.tgc)
cps.obs.exp.gcctgg<-log(dn31.codon.pair$gcctgg/cps.gcc.tgg)
cps.obs.exp.gcctgt<-log(dn31.codon.pair$gcctgt/cps.gcc.tgt)
cps.obs.exp.gcctta<-log(dn31.codon.pair$gcctta/cps.gcc.tta)
cps.obs.exp.gccttc<-log(dn31.codon.pair$gccttc/cps.gcc.ttc)
cps.obs.exp.gccttg<-log(dn31.codon.pair$gccttg/cps.gcc.ttg)
cps.obs.exp.gccttt<-log(dn31.codon.pair$gccttt/cps.gcc.ttt)









cps.obs.exp.gcgaaa<-log(dn31.codon.pair$gcgaaa/cps.gcg.aaa)
cps.obs.exp.gcgaac<-log(dn31.codon.pair$gcgaac/cps.gcg.aac)
cps.obs.exp.gcgaag<-log(dn31.codon.pair$gcgaag/cps.gcg.aag)
cps.obs.exp.gcgaat<-log(dn31.codon.pair$gcgaat/cps.gcg.aat)
cps.obs.exp.gcgaca<-log(dn31.codon.pair$gcgaca/cps.gcg.aca)
cps.obs.exp.gcgacc<-log(dn31.codon.pair$gcgacc/cps.gcg.acc)
cps.obs.exp.gcgacg<-log(dn31.codon.pair$gcgacg/cps.gcg.acg)
cps.obs.exp.gcgact<-log(dn31.codon.pair$gcgact/cps.gcg.act)
cps.obs.exp.gcgaga<-log(dn31.codon.pair$gcgaga/cps.gcg.aga)
cps.obs.exp.gcgagc<-log(dn31.codon.pair$gcgagc/cps.gcg.agc)
cps.obs.exp.gcgagg<-log(dn31.codon.pair$gcgagg/cps.gcg.agg)
cps.obs.exp.gcgagt<-log(dn31.codon.pair$gcgagt/cps.gcg.agt)
cps.obs.exp.gcgata<-log(dn31.codon.pair$gcgata/cps.gcg.ata)
cps.obs.exp.gcgatc<-log(dn31.codon.pair$gcgatc/cps.gcg.atc)
cps.obs.exp.gcgatg<-log(dn31.codon.pair$gcgatg/cps.gcg.atg)
cps.obs.exp.gcgatt<-log(dn31.codon.pair$gcgatt/cps.gcg.att)

cps.obs.exp.gcgcaa<-log(dn31.codon.pair$gcgcaa/cps.gcg.caa)
cps.obs.exp.gcgcac<-log(dn31.codon.pair$gcgcac/cps.gcg.cac)
cps.obs.exp.gcgcag<-log(dn31.codon.pair$gcgcag/cps.gcg.cag)
cps.obs.exp.gcgcat<-log(dn31.codon.pair$gcgcat/cps.gcg.cat)
cps.obs.exp.gcgcca<-log(dn31.codon.pair$gcgcca/cps.gcg.cca)
cps.obs.exp.gcgccc<-log(dn31.codon.pair$gcgccc/cps.gcg.ccc)
cps.obs.exp.gcgccg<-log(dn31.codon.pair$gcgccg/cps.gcg.ccg)
cps.obs.exp.gcgcct<-log(dn31.codon.pair$gcgcct/cps.gcg.cct)
cps.obs.exp.gcgcga<-log(dn31.codon.pair$gcgcga/cps.gcg.cga)
cps.obs.exp.gcgcgc<-log(dn31.codon.pair$gcgcgc/cps.gcg.cgc)
cps.obs.exp.gcgcgg<-log(dn31.codon.pair$gcgcgg/cps.gcg.cgg)
cps.obs.exp.gcgcgt<-log(dn31.codon.pair$gcgcgt/cps.gcg.cgt)
cps.obs.exp.gcgcta<-log(dn31.codon.pair$gcgcta/cps.gcg.cta)
cps.obs.exp.gcgctc<-log(dn31.codon.pair$gcgctc/cps.gcg.ctc)
cps.obs.exp.gcgctg<-log(dn31.codon.pair$gcgctg/cps.gcg.ctg)
cps.obs.exp.gcgctt<-log(dn31.codon.pair$gcgctt/cps.gcg.ctt)

cps.obs.exp.gcggaa<-log(dn31.codon.pair$gcggaa/cps.gcg.gaa)
cps.obs.exp.gcggac<-log(dn31.codon.pair$gcggac/cps.gcg.gac)
cps.obs.exp.gcggag<-log(dn31.codon.pair$gcggag/cps.gcg.gag)
cps.obs.exp.gcggat<-log(dn31.codon.pair$gcggat/cps.gcg.gat)
cps.obs.exp.gcggca<-log(dn31.codon.pair$gcggca/cps.gcg.gca)
cps.obs.exp.gcggcc<-log(dn31.codon.pair$gcggcc/cps.gcg.gcc)
cps.obs.exp.gcggcg<-log(dn31.codon.pair$gcggcg/cps.gcg.gcg)
cps.obs.exp.gcggct<-log(dn31.codon.pair$gcggct/cps.gcg.gct)
cps.obs.exp.gcggga<-log(dn31.codon.pair$gcggga/cps.gcg.gga)
cps.obs.exp.gcgggc<-log(dn31.codon.pair$gcgggc/cps.gcg.ggc)
cps.obs.exp.gcgggg<-log(dn31.codon.pair$gcgggg/cps.gcg.ggg)
cps.obs.exp.gcgggt<-log(dn31.codon.pair$gcgggt/cps.gcg.ggt)
cps.obs.exp.gcggta<-log(dn31.codon.pair$gcggta/cps.gcg.gta)
cps.obs.exp.gcggtc<-log(dn31.codon.pair$gcggtc/cps.gcg.gtc)
cps.obs.exp.gcggtg<-log(dn31.codon.pair$gcggtg/cps.gcg.gtg)
cps.obs.exp.gcggtt<-log(dn31.codon.pair$gcggtt/cps.gcg.gtt)

#cps.obs.exp.gcgtaa<-log(dn31.codon.pair$gcgtaa/cps.gcg.taa)
cps.obs.exp.gcgtac<-log(dn31.codon.pair$gcgtac/cps.gcg.tac)
#cps.obs.exp.gcgtag<-log(dn31.codon.pair$gcgtag/cps.gcg.tag)
cps.obs.exp.gcgtat<-log(dn31.codon.pair$gcgtat/cps.gcg.tat)
cps.obs.exp.gcgtca<-log(dn31.codon.pair$gcgtca/cps.gcg.tca)
cps.obs.exp.gcgtcc<-log(dn31.codon.pair$gcgtcc/cps.gcg.tcc)
cps.obs.exp.gcgtcg<-log(dn31.codon.pair$gcgtcg/cps.gcg.tcg)
cps.obs.exp.gcgtct<-log(dn31.codon.pair$gcgtct/cps.gcg.tct)
#cps.obs.exp.gcgtga<-log(dn31.codon.pair$gcgtga/cps.gcg.tga)
cps.obs.exp.gcgtgc<-log(dn31.codon.pair$gcgtgc/cps.gcg.tgc)
cps.obs.exp.gcgtgg<-log(dn31.codon.pair$gcgtgg/cps.gcg.tgg)
cps.obs.exp.gcgtgt<-log(dn31.codon.pair$gcgtgt/cps.gcg.tgt)
cps.obs.exp.gcgtta<-log(dn31.codon.pair$gcgtta/cps.gcg.tta)
cps.obs.exp.gcgttc<-log(dn31.codon.pair$gcgttc/cps.gcg.ttc)
cps.obs.exp.gcgttg<-log(dn31.codon.pair$gcgttg/cps.gcg.ttg)
cps.obs.exp.gcgttt<-log(dn31.codon.pair$gcgttt/cps.gcg.ttt)









cps.obs.exp.gctaaa<-log(dn31.codon.pair$gctaaa/cps.gct.aaa)
cps.obs.exp.gctaac<-log(dn31.codon.pair$gctaac/cps.gct.aac)
cps.obs.exp.gctaag<-log(dn31.codon.pair$gctaag/cps.gct.aag)
cps.obs.exp.gctaat<-log(dn31.codon.pair$gctaat/cps.gct.aat)
cps.obs.exp.gctaca<-log(dn31.codon.pair$gctaca/cps.gct.aca)
cps.obs.exp.gctacc<-log(dn31.codon.pair$gctacc/cps.gct.acc)
cps.obs.exp.gctacg<-log(dn31.codon.pair$gctacg/cps.gct.acg)
cps.obs.exp.gctact<-log(dn31.codon.pair$gctact/cps.gct.act)
cps.obs.exp.gctaga<-log(dn31.codon.pair$gctaga/cps.gct.aga)
cps.obs.exp.gctagc<-log(dn31.codon.pair$gctagc/cps.gct.agc)
cps.obs.exp.gctagg<-log(dn31.codon.pair$gctagg/cps.gct.agg)
cps.obs.exp.gctagt<-log(dn31.codon.pair$gctagt/cps.gct.agt)
cps.obs.exp.gctata<-log(dn31.codon.pair$gctata/cps.gct.ata)
cps.obs.exp.gctatc<-log(dn31.codon.pair$gctatc/cps.gct.atc)
cps.obs.exp.gctatg<-log(dn31.codon.pair$gctatg/cps.gct.atg)
cps.obs.exp.gctatt<-log(dn31.codon.pair$gctatt/cps.gct.att)

cps.obs.exp.gctcaa<-log(dn31.codon.pair$gctcaa/cps.gct.caa)
cps.obs.exp.gctcac<-log(dn31.codon.pair$gctcac/cps.gct.cac)
cps.obs.exp.gctcag<-log(dn31.codon.pair$gctcag/cps.gct.cag)
cps.obs.exp.gctcat<-log(dn31.codon.pair$gctcat/cps.gct.cat)
cps.obs.exp.gctcca<-log(dn31.codon.pair$gctcca/cps.gct.cca)
cps.obs.exp.gctccc<-log(dn31.codon.pair$gctccc/cps.gct.ccc)
cps.obs.exp.gctccg<-log(dn31.codon.pair$gctccg/cps.gct.ccg)
cps.obs.exp.gctcct<-log(dn31.codon.pair$gctcct/cps.gct.cct)
cps.obs.exp.gctcga<-log(dn31.codon.pair$gctcga/cps.gct.cga)
cps.obs.exp.gctcgc<-log(dn31.codon.pair$gctcgc/cps.gct.cgc)
cps.obs.exp.gctcgg<-log(dn31.codon.pair$gctcgg/cps.gct.cgg)
cps.obs.exp.gctcgt<-log(dn31.codon.pair$gctcgt/cps.gct.cgt)
cps.obs.exp.gctcta<-log(dn31.codon.pair$gctcta/cps.gct.cta)
cps.obs.exp.gctctc<-log(dn31.codon.pair$gctctc/cps.gct.ctc)
cps.obs.exp.gctctg<-log(dn31.codon.pair$gctctg/cps.gct.ctg)
cps.obs.exp.gctctt<-log(dn31.codon.pair$gctctt/cps.gct.ctt)

cps.obs.exp.gctgaa<-log(dn31.codon.pair$gctgaa/cps.gct.gaa)
cps.obs.exp.gctgac<-log(dn31.codon.pair$gctgac/cps.gct.gac)
cps.obs.exp.gctgag<-log(dn31.codon.pair$gctgag/cps.gct.gag)
cps.obs.exp.gctgat<-log(dn31.codon.pair$gctgat/cps.gct.gat)
cps.obs.exp.gctgca<-log(dn31.codon.pair$gctgca/cps.gct.gca)
cps.obs.exp.gctgcc<-log(dn31.codon.pair$gctgcc/cps.gct.gcc)
cps.obs.exp.gctgcg<-log(dn31.codon.pair$gctgcg/cps.gct.gcg)
cps.obs.exp.gctgct<-log(dn31.codon.pair$gctgct/cps.gct.gct)
cps.obs.exp.gctgga<-log(dn31.codon.pair$gctgga/cps.gct.gga)
cps.obs.exp.gctggc<-log(dn31.codon.pair$gctggc/cps.gct.ggc)
cps.obs.exp.gctggg<-log(dn31.codon.pair$gctggg/cps.gct.ggg)
cps.obs.exp.gctggt<-log(dn31.codon.pair$gctggt/cps.gct.ggt)
cps.obs.exp.gctgta<-log(dn31.codon.pair$gctgta/cps.gct.gta)
cps.obs.exp.gctgtc<-log(dn31.codon.pair$gctgtc/cps.gct.gtc)
cps.obs.exp.gctgtg<-log(dn31.codon.pair$gctgtg/cps.gct.gtg)
cps.obs.exp.gctgtt<-log(dn31.codon.pair$gctgtt/cps.gct.gtt)

#cps.obs.exp.gcttaa<-log(dn31.codon.pair$gcttaa/cps.gct.taa)
cps.obs.exp.gcttac<-log(dn31.codon.pair$gcttac/cps.gct.tac)
#cps.obs.exp.gcttag<-log(dn31.codon.pair$gcttag/cps.gct.tag)
cps.obs.exp.gcttat<-log(dn31.codon.pair$gcttat/cps.gct.tat)
cps.obs.exp.gcttca<-log(dn31.codon.pair$gcttca/cps.gct.tca)
cps.obs.exp.gcttcc<-log(dn31.codon.pair$gcttcc/cps.gct.tcc)
cps.obs.exp.gcttcg<-log(dn31.codon.pair$gcttcg/cps.gct.tcg)
cps.obs.exp.gcttct<-log(dn31.codon.pair$gcttct/cps.gct.tct)
#cps.obs.exp.gcttga<-log(dn31.codon.pair$gcttga/cps.gct.tga)
cps.obs.exp.gcttgc<-log(dn31.codon.pair$gcttgc/cps.gct.tgc)
cps.obs.exp.gcttgg<-log(dn31.codon.pair$gcttgg/cps.gct.tgg)
cps.obs.exp.gcttgt<-log(dn31.codon.pair$gcttgt/cps.gct.tgt)
cps.obs.exp.gcttta<-log(dn31.codon.pair$gcttta/cps.gct.tta)
cps.obs.exp.gctttc<-log(dn31.codon.pair$gctttc/cps.gct.ttc)
cps.obs.exp.gctttg<-log(dn31.codon.pair$gctttg/cps.gct.ttg)
cps.obs.exp.gctttt<-log(dn31.codon.pair$gctttt/cps.gct.ttt)



















cps.obs.exp.ggaaaa<-log(dn31.codon.pair$ggaaaa/cps.gga.aaa)
cps.obs.exp.ggaaac<-log(dn31.codon.pair$ggaaac/cps.gga.aac)
cps.obs.exp.ggaaag<-log(dn31.codon.pair$ggaaag/cps.gga.aag)
cps.obs.exp.ggaaat<-log(dn31.codon.pair$ggaaat/cps.gga.aat)
cps.obs.exp.ggaaca<-log(dn31.codon.pair$ggaaca/cps.gga.aca)
cps.obs.exp.ggaacc<-log(dn31.codon.pair$ggaacc/cps.gga.acc)
cps.obs.exp.ggaacg<-log(dn31.codon.pair$ggaacg/cps.gga.acg)
cps.obs.exp.ggaact<-log(dn31.codon.pair$ggaact/cps.gga.act)
cps.obs.exp.ggaaga<-log(dn31.codon.pair$ggaaga/cps.gga.aga)
cps.obs.exp.ggaagc<-log(dn31.codon.pair$ggaagc/cps.gga.agc)
cps.obs.exp.ggaagg<-log(dn31.codon.pair$ggaagg/cps.gga.agg)
cps.obs.exp.ggaagt<-log(dn31.codon.pair$ggaagt/cps.gga.agt)
cps.obs.exp.ggaata<-log(dn31.codon.pair$ggaata/cps.gga.ata)
cps.obs.exp.ggaatc<-log(dn31.codon.pair$ggaatc/cps.gga.atc)
cps.obs.exp.ggaatg<-log(dn31.codon.pair$ggaatg/cps.gga.atg)
cps.obs.exp.ggaatt<-log(dn31.codon.pair$ggaatt/cps.gga.att)

cps.obs.exp.ggacaa<-log(dn31.codon.pair$ggacaa/cps.gga.caa)
cps.obs.exp.ggacac<-log(dn31.codon.pair$ggacac/cps.gga.cac)
cps.obs.exp.ggacag<-log(dn31.codon.pair$ggacag/cps.gga.cag)
cps.obs.exp.ggacat<-log(dn31.codon.pair$ggacat/cps.gga.cat)
cps.obs.exp.ggacca<-log(dn31.codon.pair$ggacca/cps.gga.cca)
cps.obs.exp.ggaccc<-log(dn31.codon.pair$ggaccc/cps.gga.ccc)
cps.obs.exp.ggaccg<-log(dn31.codon.pair$ggaccg/cps.gga.ccg)
cps.obs.exp.ggacct<-log(dn31.codon.pair$ggacct/cps.gga.cct)
cps.obs.exp.ggacga<-log(dn31.codon.pair$ggacga/cps.gga.cga)
cps.obs.exp.ggacgc<-log(dn31.codon.pair$ggacgc/cps.gga.cgc)
cps.obs.exp.ggacgg<-log(dn31.codon.pair$ggacgg/cps.gga.cgg)
cps.obs.exp.ggacgt<-log(dn31.codon.pair$ggacgt/cps.gga.cgt)
cps.obs.exp.ggacta<-log(dn31.codon.pair$ggacta/cps.gga.cta)
cps.obs.exp.ggactc<-log(dn31.codon.pair$ggactc/cps.gga.ctc)
cps.obs.exp.ggactg<-log(dn31.codon.pair$ggactg/cps.gga.ctg)
cps.obs.exp.ggactt<-log(dn31.codon.pair$ggactt/cps.gga.ctt)

cps.obs.exp.ggagaa<-log(dn31.codon.pair$ggagaa/cps.gga.gaa)
cps.obs.exp.ggagac<-log(dn31.codon.pair$ggagac/cps.gga.gac)
cps.obs.exp.ggagag<-log(dn31.codon.pair$ggagag/cps.gga.gag)
cps.obs.exp.ggagat<-log(dn31.codon.pair$ggagat/cps.gga.gat)
cps.obs.exp.ggagca<-log(dn31.codon.pair$ggagca/cps.gga.gca)
cps.obs.exp.ggagcc<-log(dn31.codon.pair$ggagcc/cps.gga.gcc)
cps.obs.exp.ggagcg<-log(dn31.codon.pair$ggagcg/cps.gga.gcg)
cps.obs.exp.ggagct<-log(dn31.codon.pair$ggagct/cps.gga.gct)
cps.obs.exp.ggagga<-log(dn31.codon.pair$ggagga/cps.gga.gga)
cps.obs.exp.ggaggc<-log(dn31.codon.pair$ggaggc/cps.gga.ggc)
cps.obs.exp.ggaggg<-log(dn31.codon.pair$ggaggg/cps.gga.ggg)
cps.obs.exp.ggaggt<-log(dn31.codon.pair$ggaggt/cps.gga.ggt)
cps.obs.exp.ggagta<-log(dn31.codon.pair$ggagta/cps.gga.gta)
cps.obs.exp.ggagtc<-log(dn31.codon.pair$ggagtc/cps.gga.gtc)
cps.obs.exp.ggagtg<-log(dn31.codon.pair$ggagtg/cps.gga.gtg)
cps.obs.exp.ggagtt<-log(dn31.codon.pair$ggagtt/cps.gga.gtt)

#cps.obs.exp.ggataa<-log(dn31.codon.pair$ggataa/cps.gga.taa)
cps.obs.exp.ggatac<-log(dn31.codon.pair$ggatac/cps.gga.tac)
#cps.obs.exp.ggatag<-log(dn31.codon.pair$ggatag/cps.gga.tag)
cps.obs.exp.ggatat<-log(dn31.codon.pair$ggatat/cps.gga.tat)
cps.obs.exp.ggatca<-log(dn31.codon.pair$ggatca/cps.gga.tca)
cps.obs.exp.ggatcc<-log(dn31.codon.pair$ggatcc/cps.gga.tcc)
cps.obs.exp.ggatcg<-log(dn31.codon.pair$ggatcg/cps.gga.tcg)
cps.obs.exp.ggatct<-log(dn31.codon.pair$ggatct/cps.gga.tct)
#cps.obs.exp.ggatga<-log(dn31.codon.pair$ggatga/cps.gga.tga)
cps.obs.exp.ggatgc<-log(dn31.codon.pair$ggatgc/cps.gga.tgc)
cps.obs.exp.ggatgg<-log(dn31.codon.pair$ggatgg/cps.gga.tgg)
cps.obs.exp.ggatgt<-log(dn31.codon.pair$ggatgt/cps.gga.tgt)
cps.obs.exp.ggatta<-log(dn31.codon.pair$ggatta/cps.gga.tta)
cps.obs.exp.ggattc<-log(dn31.codon.pair$ggattc/cps.gga.ttc)
cps.obs.exp.ggattg<-log(dn31.codon.pair$ggattg/cps.gga.ttg)
cps.obs.exp.ggattt<-log(dn31.codon.pair$ggattt/cps.gga.ttt)










cps.obs.exp.ggcaaa<-log(dn31.codon.pair$ggcaaa/cps.ggc.aaa)
cps.obs.exp.ggcaac<-log(dn31.codon.pair$ggcaac/cps.ggc.aac)
cps.obs.exp.ggcaag<-log(dn31.codon.pair$ggcaag/cps.ggc.aag)
cps.obs.exp.ggcaat<-log(dn31.codon.pair$ggcaat/cps.ggc.aat)
cps.obs.exp.ggcaca<-log(dn31.codon.pair$ggcaca/cps.ggc.aca)
cps.obs.exp.ggcacc<-log(dn31.codon.pair$ggcacc/cps.ggc.acc)
cps.obs.exp.ggcacg<-log(dn31.codon.pair$ggcacg/cps.ggc.acg)
cps.obs.exp.ggcact<-log(dn31.codon.pair$ggcact/cps.ggc.act)
cps.obs.exp.ggcaga<-log(dn31.codon.pair$ggcaga/cps.ggc.aga)
cps.obs.exp.ggcagc<-log(dn31.codon.pair$ggcagc/cps.ggc.agc)
cps.obs.exp.ggcagg<-log(dn31.codon.pair$ggcagg/cps.ggc.agg)
cps.obs.exp.ggcagt<-log(dn31.codon.pair$ggcagt/cps.ggc.agt)
cps.obs.exp.ggcata<-log(dn31.codon.pair$ggcata/cps.ggc.ata)
cps.obs.exp.ggcatc<-log(dn31.codon.pair$ggcatc/cps.ggc.atc)
cps.obs.exp.ggcatg<-log(dn31.codon.pair$ggcatg/cps.ggc.atg)
cps.obs.exp.ggcatt<-log(dn31.codon.pair$ggcatt/cps.ggc.att)

cps.obs.exp.ggccaa<-log(dn31.codon.pair$ggccaa/cps.ggc.caa)
cps.obs.exp.ggccac<-log(dn31.codon.pair$ggccac/cps.ggc.cac)
cps.obs.exp.ggccag<-log(dn31.codon.pair$ggccag/cps.ggc.cag)
cps.obs.exp.ggccat<-log(dn31.codon.pair$ggccat/cps.ggc.cat)
cps.obs.exp.ggccca<-log(dn31.codon.pair$ggccca/cps.ggc.cca)
cps.obs.exp.ggcccc<-log(dn31.codon.pair$ggcccc/cps.ggc.ccc)
cps.obs.exp.ggcccg<-log(dn31.codon.pair$ggcccg/cps.ggc.ccg)
cps.obs.exp.ggccct<-log(dn31.codon.pair$ggccct/cps.ggc.cct)
cps.obs.exp.ggccga<-log(dn31.codon.pair$ggccga/cps.ggc.cga)
cps.obs.exp.ggccgc<-log(dn31.codon.pair$ggccgc/cps.ggc.cgc)
cps.obs.exp.ggccgg<-log(dn31.codon.pair$ggccgg/cps.ggc.cgg)
cps.obs.exp.ggccgt<-log(dn31.codon.pair$ggccgt/cps.ggc.cgt)
cps.obs.exp.ggccta<-log(dn31.codon.pair$ggccta/cps.ggc.cta)
cps.obs.exp.ggcctc<-log(dn31.codon.pair$ggcctc/cps.ggc.ctc)
cps.obs.exp.ggcctg<-log(dn31.codon.pair$ggcctg/cps.ggc.ctg)
cps.obs.exp.ggcctt<-log(dn31.codon.pair$ggcctt/cps.ggc.ctt)

cps.obs.exp.ggcgaa<-log(dn31.codon.pair$ggcgaa/cps.ggc.gaa)
cps.obs.exp.ggcgac<-log(dn31.codon.pair$ggcgac/cps.ggc.gac)
cps.obs.exp.ggcgag<-log(dn31.codon.pair$ggcgag/cps.ggc.gag)
cps.obs.exp.ggcgat<-log(dn31.codon.pair$ggcgat/cps.ggc.gat)
cps.obs.exp.ggcgca<-log(dn31.codon.pair$ggcgca/cps.ggc.gca)
cps.obs.exp.ggcgcc<-log(dn31.codon.pair$ggcgcc/cps.ggc.gcc)
cps.obs.exp.ggcgcg<-log(dn31.codon.pair$ggcgcg/cps.ggc.gcg)
cps.obs.exp.ggcgct<-log(dn31.codon.pair$ggcgct/cps.ggc.gct)
cps.obs.exp.ggcgga<-log(dn31.codon.pair$ggcgga/cps.ggc.gga)
cps.obs.exp.ggcggc<-log(dn31.codon.pair$ggcggc/cps.ggc.ggc)
cps.obs.exp.ggcggg<-log(dn31.codon.pair$ggcggg/cps.ggc.ggg)
cps.obs.exp.ggcggt<-log(dn31.codon.pair$ggcggt/cps.ggc.ggt)
cps.obs.exp.ggcgta<-log(dn31.codon.pair$ggcgta/cps.ggc.gta)
cps.obs.exp.ggcgtc<-log(dn31.codon.pair$ggcgtc/cps.ggc.gtc)
cps.obs.exp.ggcgtg<-log(dn31.codon.pair$ggcgtg/cps.ggc.gtg)
cps.obs.exp.ggcgtt<-log(dn31.codon.pair$ggcgtt/cps.ggc.gtt)

#cps.obs.exp.ggctaa<-log(dn31.codon.pair$ggctaa/cps.ggc.taa)
cps.obs.exp.ggctac<-log(dn31.codon.pair$ggctac/cps.ggc.tac)
#cps.obs.exp.ggctag<-log(dn31.codon.pair$ggctag/cps.ggc.tag)
cps.obs.exp.ggctat<-log(dn31.codon.pair$ggctat/cps.ggc.tat)
cps.obs.exp.ggctca<-log(dn31.codon.pair$ggctca/cps.ggc.tca)
cps.obs.exp.ggctcc<-log(dn31.codon.pair$ggctcc/cps.ggc.tcc)
cps.obs.exp.ggctcg<-log(dn31.codon.pair$ggctcg/cps.ggc.tcg)
cps.obs.exp.ggctct<-log(dn31.codon.pair$ggctct/cps.ggc.tct)
#cps.obs.exp.ggctga<-log(dn31.codon.pair$ggctga/cps.ggc.tga)
cps.obs.exp.ggctgc<-log(dn31.codon.pair$ggctgc/cps.ggc.tgc)
cps.obs.exp.ggctgg<-log(dn31.codon.pair$ggctgg/cps.ggc.tgg)
cps.obs.exp.ggctgt<-log(dn31.codon.pair$ggctgt/cps.ggc.tgt)
cps.obs.exp.ggctta<-log(dn31.codon.pair$ggctta/cps.ggc.tta)
cps.obs.exp.ggcttc<-log(dn31.codon.pair$ggcttc/cps.ggc.ttc)
cps.obs.exp.ggcttg<-log(dn31.codon.pair$ggcttg/cps.ggc.ttg)
cps.obs.exp.ggcttt<-log(dn31.codon.pair$ggcttt/cps.ggc.ttt)









cps.obs.exp.gggaaa<-log(dn31.codon.pair$gggaaa/cps.ggg.aaa)
cps.obs.exp.gggaac<-log(dn31.codon.pair$gggaac/cps.ggg.aac)
cps.obs.exp.gggaag<-log(dn31.codon.pair$gggaag/cps.ggg.aag)
cps.obs.exp.gggaat<-log(dn31.codon.pair$gggaat/cps.ggg.aat)
cps.obs.exp.gggaca<-log(dn31.codon.pair$gggaca/cps.ggg.aca)
cps.obs.exp.gggacc<-log(dn31.codon.pair$gggacc/cps.ggg.acc)
cps.obs.exp.gggacg<-log(dn31.codon.pair$gggacg/cps.ggg.acg)
cps.obs.exp.gggact<-log(dn31.codon.pair$gggact/cps.ggg.act)
cps.obs.exp.gggaga<-log(dn31.codon.pair$gggaga/cps.ggg.aga)
cps.obs.exp.gggagc<-log(dn31.codon.pair$gggagc/cps.ggg.agc)
cps.obs.exp.gggagg<-log(dn31.codon.pair$gggagg/cps.ggg.agg)
cps.obs.exp.gggagt<-log(dn31.codon.pair$gggagt/cps.ggg.agt)
cps.obs.exp.gggata<-log(dn31.codon.pair$gggata/cps.ggg.ata)
cps.obs.exp.gggatc<-log(dn31.codon.pair$gggatc/cps.ggg.atc)
cps.obs.exp.gggatg<-log(dn31.codon.pair$gggatg/cps.ggg.atg)
cps.obs.exp.gggatt<-log(dn31.codon.pair$gggatt/cps.ggg.att)

cps.obs.exp.gggcaa<-log(dn31.codon.pair$gggcaa/cps.ggg.caa)
cps.obs.exp.gggcac<-log(dn31.codon.pair$gggcac/cps.ggg.cac)
cps.obs.exp.gggcag<-log(dn31.codon.pair$gggcag/cps.ggg.cag)
cps.obs.exp.gggcat<-log(dn31.codon.pair$gggcat/cps.ggg.cat)
cps.obs.exp.gggcca<-log(dn31.codon.pair$gggcca/cps.ggg.cca)
cps.obs.exp.gggccc<-log(dn31.codon.pair$gggccc/cps.ggg.ccc)
cps.obs.exp.gggccg<-log(dn31.codon.pair$gggccg/cps.ggg.ccg)
cps.obs.exp.gggcct<-log(dn31.codon.pair$gggcct/cps.ggg.cct)
cps.obs.exp.gggcga<-log(dn31.codon.pair$gggcga/cps.ggg.cga)
cps.obs.exp.gggcgc<-log(dn31.codon.pair$gggcgc/cps.ggg.cgc)
cps.obs.exp.gggcgg<-log(dn31.codon.pair$gggcgg/cps.ggg.cgg)
cps.obs.exp.gggcgt<-log(dn31.codon.pair$gggcgt/cps.ggg.cgt)
cps.obs.exp.gggcta<-log(dn31.codon.pair$gggcta/cps.ggg.cta)
cps.obs.exp.gggctc<-log(dn31.codon.pair$gggctc/cps.ggg.ctc)
cps.obs.exp.gggctg<-log(dn31.codon.pair$gggctg/cps.ggg.ctg)
cps.obs.exp.gggctt<-log(dn31.codon.pair$gggctt/cps.ggg.ctt)

cps.obs.exp.ggggaa<-log(dn31.codon.pair$ggggaa/cps.ggg.gaa)
cps.obs.exp.ggggac<-log(dn31.codon.pair$ggggac/cps.ggg.gac)
cps.obs.exp.ggggag<-log(dn31.codon.pair$ggggag/cps.ggg.gag)
cps.obs.exp.ggggat<-log(dn31.codon.pair$ggggat/cps.ggg.gat)
cps.obs.exp.ggggca<-log(dn31.codon.pair$ggggca/cps.ggg.gca)
cps.obs.exp.ggggcc<-log(dn31.codon.pair$ggggcc/cps.ggg.gcc)
cps.obs.exp.ggggcg<-log(dn31.codon.pair$ggggcg/cps.ggg.gcg)
cps.obs.exp.ggggct<-log(dn31.codon.pair$ggggct/cps.ggg.gct)
cps.obs.exp.ggggga<-log(dn31.codon.pair$ggggga/cps.ggg.gga)
cps.obs.exp.gggggc<-log(dn31.codon.pair$gggggc/cps.ggg.ggc)
cps.obs.exp.gggggg<-log(dn31.codon.pair$gggggg/cps.ggg.ggg)
cps.obs.exp.gggggt<-log(dn31.codon.pair$gggggt/cps.ggg.ggt)
cps.obs.exp.ggggta<-log(dn31.codon.pair$ggggta/cps.ggg.gta)
cps.obs.exp.ggggtc<-log(dn31.codon.pair$ggggtc/cps.ggg.gtc)
cps.obs.exp.ggggtg<-log(dn31.codon.pair$ggggtg/cps.ggg.gtg)
cps.obs.exp.ggggtt<-log(dn31.codon.pair$ggggtt/cps.ggg.gtt)

#cps.obs.exp.gggtaa<-log(dn31.codon.pair$gggtaa/cps.ggg.taa)
cps.obs.exp.gggtac<-log(dn31.codon.pair$gggtac/cps.ggg.tac)
#cps.obs.exp.gggtag<-log(dn31.codon.pair$gggtag/cps.ggg.tag)
cps.obs.exp.gggtat<-log(dn31.codon.pair$gggtat/cps.ggg.tat)
cps.obs.exp.gggtca<-log(dn31.codon.pair$gggtca/cps.ggg.tca)
cps.obs.exp.gggtcc<-log(dn31.codon.pair$gggtcc/cps.ggg.tcc)
cps.obs.exp.gggtcg<-log(dn31.codon.pair$gggtcg/cps.ggg.tcg)
cps.obs.exp.gggtct<-log(dn31.codon.pair$gggtct/cps.ggg.tct)
#cps.obs.exp.gggtga<-log(dn31.codon.pair$gggtga/cps.ggg.tga)
cps.obs.exp.gggtgc<-log(dn31.codon.pair$gggtgc/cps.ggg.tgc)
cps.obs.exp.gggtgg<-log(dn31.codon.pair$gggtgg/cps.ggg.tgg)
cps.obs.exp.gggtgt<-log(dn31.codon.pair$gggtgt/cps.ggg.tgt)
cps.obs.exp.gggtta<-log(dn31.codon.pair$gggtta/cps.ggg.tta)
cps.obs.exp.gggttc<-log(dn31.codon.pair$gggttc/cps.ggg.ttc)
cps.obs.exp.gggttg<-log(dn31.codon.pair$gggttg/cps.ggg.ttg)
cps.obs.exp.gggttt<-log(dn31.codon.pair$gggttt/cps.ggg.ttt)











cps.obs.exp.ggtaaa<-log(dn31.codon.pair$ggtaaa/cps.ggt.aaa)
cps.obs.exp.ggtaac<-log(dn31.codon.pair$ggtaac/cps.ggt.aac)
cps.obs.exp.ggtaag<-log(dn31.codon.pair$ggtaag/cps.ggt.aag)
cps.obs.exp.ggtaat<-log(dn31.codon.pair$ggtaat/cps.ggt.aat)
cps.obs.exp.ggtaca<-log(dn31.codon.pair$ggtaca/cps.ggt.aca)
cps.obs.exp.ggtacc<-log(dn31.codon.pair$ggtacc/cps.ggt.acc)
cps.obs.exp.ggtacg<-log(dn31.codon.pair$ggtacg/cps.ggt.acg)
cps.obs.exp.ggtact<-log(dn31.codon.pair$ggtact/cps.ggt.act)
cps.obs.exp.ggtaga<-log(dn31.codon.pair$ggtaga/cps.ggt.aga)
cps.obs.exp.ggtagc<-log(dn31.codon.pair$ggtagc/cps.ggt.agc)
cps.obs.exp.ggtagg<-log(dn31.codon.pair$ggtagg/cps.ggt.agg)
cps.obs.exp.ggtagt<-log(dn31.codon.pair$ggtagt/cps.ggt.agt)
cps.obs.exp.ggtata<-log(dn31.codon.pair$ggtata/cps.ggt.ata)
cps.obs.exp.ggtatc<-log(dn31.codon.pair$ggtatc/cps.ggt.atc)
cps.obs.exp.ggtatg<-log(dn31.codon.pair$ggtatg/cps.ggt.atg)
cps.obs.exp.ggtatt<-log(dn31.codon.pair$ggtatt/cps.ggt.att)

cps.obs.exp.ggtcaa<-log(dn31.codon.pair$ggtcaa/cps.ggt.caa)
cps.obs.exp.ggtcac<-log(dn31.codon.pair$ggtcac/cps.ggt.cac)
cps.obs.exp.ggtcag<-log(dn31.codon.pair$ggtcag/cps.ggt.cag)
cps.obs.exp.ggtcat<-log(dn31.codon.pair$ggtcat/cps.ggt.cat)
cps.obs.exp.ggtcca<-log(dn31.codon.pair$ggtcca/cps.ggt.cca)
cps.obs.exp.ggtccc<-log(dn31.codon.pair$ggtccc/cps.ggt.ccc)
cps.obs.exp.ggtccg<-log(dn31.codon.pair$ggtccg/cps.ggt.ccg)
cps.obs.exp.ggtcct<-log(dn31.codon.pair$ggtcct/cps.ggt.cct)
cps.obs.exp.ggtcga<-log(dn31.codon.pair$ggtcga/cps.ggt.cga)
cps.obs.exp.ggtcgc<-log(dn31.codon.pair$ggtcgc/cps.ggt.cgc)
cps.obs.exp.ggtcgg<-log(dn31.codon.pair$ggtcgg/cps.ggt.cgg)
cps.obs.exp.ggtcgt<-log(dn31.codon.pair$ggtcgt/cps.ggt.cgt)
cps.obs.exp.ggtcta<-log(dn31.codon.pair$ggtcta/cps.ggt.cta)
cps.obs.exp.ggtctc<-log(dn31.codon.pair$ggtctc/cps.ggt.ctc)
cps.obs.exp.ggtctg<-log(dn31.codon.pair$ggtctg/cps.ggt.ctg)
cps.obs.exp.ggtctt<-log(dn31.codon.pair$ggtctt/cps.ggt.ctt)

cps.obs.exp.ggtgaa<-log(dn31.codon.pair$ggtgaa/cps.ggt.gaa)
cps.obs.exp.ggtgac<-log(dn31.codon.pair$ggtgac/cps.ggt.gac)
cps.obs.exp.ggtgag<-log(dn31.codon.pair$ggtgag/cps.ggt.gag)
cps.obs.exp.ggtgat<-log(dn31.codon.pair$ggtgat/cps.ggt.gat)
cps.obs.exp.ggtgca<-log(dn31.codon.pair$ggtgca/cps.ggt.gca)
cps.obs.exp.ggtgcc<-log(dn31.codon.pair$ggtgcc/cps.ggt.gcc)
cps.obs.exp.ggtgcg<-log(dn31.codon.pair$ggtgcg/cps.ggt.gcg)
cps.obs.exp.ggtgct<-log(dn31.codon.pair$ggtgct/cps.ggt.gct)
cps.obs.exp.ggtgga<-log(dn31.codon.pair$ggtgga/cps.ggt.gga)
cps.obs.exp.ggtggc<-log(dn31.codon.pair$ggtggc/cps.ggt.ggc)
cps.obs.exp.ggtggg<-log(dn31.codon.pair$ggtggg/cps.ggt.ggg)
cps.obs.exp.ggtggt<-log(dn31.codon.pair$ggtggt/cps.ggt.ggt)
cps.obs.exp.ggtgta<-log(dn31.codon.pair$ggtgta/cps.ggt.gta)
cps.obs.exp.ggtgtc<-log(dn31.codon.pair$ggtgtc/cps.ggt.gtc)
cps.obs.exp.ggtgtg<-log(dn31.codon.pair$ggtgtg/cps.ggt.gtg)
cps.obs.exp.ggtgtt<-log(dn31.codon.pair$ggtgtt/cps.ggt.gtt)

#cps.obs.exp.ggttaa<-log(dn31.codon.pair$ggttaa/cps.ggt.taa)
cps.obs.exp.ggttac<-log(dn31.codon.pair$ggttac/cps.ggt.tac)
#cps.obs.exp.ggttag<-log(dn31.codon.pair$ggttag/cps.ggt.tag)
cps.obs.exp.ggttat<-log(dn31.codon.pair$ggttat/cps.ggt.tat)
cps.obs.exp.ggttca<-log(dn31.codon.pair$ggttca/cps.ggt.tca)
cps.obs.exp.ggttcc<-log(dn31.codon.pair$ggttcc/cps.ggt.tcc)
cps.obs.exp.ggttcg<-log(dn31.codon.pair$ggttcg/cps.ggt.tcg)
cps.obs.exp.ggttct<-log(dn31.codon.pair$ggttct/cps.ggt.tct)
#cps.obs.exp.ggttga<-log(dn31.codon.pair$ggttga/cps.ggt.tga)
cps.obs.exp.ggttgc<-log(dn31.codon.pair$ggttgc/cps.ggt.tgc)
cps.obs.exp.ggttgg<-log(dn31.codon.pair$ggttgg/cps.ggt.tgg)
cps.obs.exp.ggttgt<-log(dn31.codon.pair$ggttgt/cps.ggt.tgt)
cps.obs.exp.ggttta<-log(dn31.codon.pair$ggttta/cps.ggt.tta)
cps.obs.exp.ggtttc<-log(dn31.codon.pair$ggtttc/cps.ggt.ttc)
cps.obs.exp.ggtttg<-log(dn31.codon.pair$ggtttg/cps.ggt.ttg)
cps.obs.exp.ggtttt<-log(dn31.codon.pair$ggtttt/cps.ggt.ttt)




















cps.obs.exp.gtaaaa<-log(dn31.codon.pair$gtaaaa/cps.gta.aaa)
cps.obs.exp.gtaaac<-log(dn31.codon.pair$gtaaac/cps.gta.aac)
cps.obs.exp.gtaaag<-log(dn31.codon.pair$gtaaag/cps.gta.aag)
cps.obs.exp.gtaaat<-log(dn31.codon.pair$gtaaat/cps.gta.aat)
cps.obs.exp.gtaaca<-log(dn31.codon.pair$gtaaca/cps.gta.aca)
cps.obs.exp.gtaacc<-log(dn31.codon.pair$gtaacc/cps.gta.acc)
cps.obs.exp.gtaacg<-log(dn31.codon.pair$gtaacg/cps.gta.acg)
cps.obs.exp.gtaact<-log(dn31.codon.pair$gtaact/cps.gta.act)
cps.obs.exp.gtaaga<-log(dn31.codon.pair$gtaaga/cps.gta.aga)
cps.obs.exp.gtaagc<-log(dn31.codon.pair$gtaagc/cps.gta.agc)
cps.obs.exp.gtaagg<-log(dn31.codon.pair$gtaagg/cps.gta.agg)
cps.obs.exp.gtaagt<-log(dn31.codon.pair$gtaagt/cps.gta.agt)
cps.obs.exp.gtaata<-log(dn31.codon.pair$gtaata/cps.gta.ata)
cps.obs.exp.gtaatc<-log(dn31.codon.pair$gtaatc/cps.gta.atc)
cps.obs.exp.gtaatg<-log(dn31.codon.pair$gtaatg/cps.gta.atg)
cps.obs.exp.gtaatt<-log(dn31.codon.pair$gtaatt/cps.gta.att)

cps.obs.exp.gtacaa<-log(dn31.codon.pair$gtacaa/cps.gta.caa)
cps.obs.exp.gtacac<-log(dn31.codon.pair$gtacac/cps.gta.cac)
cps.obs.exp.gtacag<-log(dn31.codon.pair$gtacag/cps.gta.cag)
cps.obs.exp.gtacat<-log(dn31.codon.pair$gtacat/cps.gta.cat)
cps.obs.exp.gtacca<-log(dn31.codon.pair$gtacca/cps.gta.cca)
cps.obs.exp.gtaccc<-log(dn31.codon.pair$gtaccc/cps.gta.ccc)
cps.obs.exp.gtaccg<-log(dn31.codon.pair$gtaccg/cps.gta.ccg)
cps.obs.exp.gtacct<-log(dn31.codon.pair$gtacct/cps.gta.cct)
cps.obs.exp.gtacga<-log(dn31.codon.pair$gtacga/cps.gta.cga)
cps.obs.exp.gtacgc<-log(dn31.codon.pair$gtacgc/cps.gta.cgc)
cps.obs.exp.gtacgg<-log(dn31.codon.pair$gtacgg/cps.gta.cgg)
cps.obs.exp.gtacgt<-log(dn31.codon.pair$gtacgt/cps.gta.cgt)
cps.obs.exp.gtacta<-log(dn31.codon.pair$gtacta/cps.gta.cta)
cps.obs.exp.gtactc<-log(dn31.codon.pair$gtactc/cps.gta.ctc)
cps.obs.exp.gtactg<-log(dn31.codon.pair$gtactg/cps.gta.ctg)
cps.obs.exp.gtactt<-log(dn31.codon.pair$gtactt/cps.gta.ctt)

cps.obs.exp.gtagaa<-log(dn31.codon.pair$gtagaa/cps.gta.gaa)
cps.obs.exp.gtagac<-log(dn31.codon.pair$gtagac/cps.gta.gac)
cps.obs.exp.gtagag<-log(dn31.codon.pair$gtagag/cps.gta.gag)
cps.obs.exp.gtagat<-log(dn31.codon.pair$gtagat/cps.gta.gat)
cps.obs.exp.gtagca<-log(dn31.codon.pair$gtagca/cps.gta.gca)
cps.obs.exp.gtagcc<-log(dn31.codon.pair$gtagcc/cps.gta.gcc)
cps.obs.exp.gtagcg<-log(dn31.codon.pair$gtagcg/cps.gta.gcg)
cps.obs.exp.gtagct<-log(dn31.codon.pair$gtagct/cps.gta.gct)
cps.obs.exp.gtagga<-log(dn31.codon.pair$gtagga/cps.gta.gga)
cps.obs.exp.gtaggc<-log(dn31.codon.pair$gtaggc/cps.gta.ggc)
cps.obs.exp.gtaggg<-log(dn31.codon.pair$gtaggg/cps.gta.ggg)
cps.obs.exp.gtaggt<-log(dn31.codon.pair$gtaggt/cps.gta.ggt)
cps.obs.exp.gtagta<-log(dn31.codon.pair$gtagta/cps.gta.gta)
cps.obs.exp.gtagtc<-log(dn31.codon.pair$gtagtc/cps.gta.gtc)
cps.obs.exp.gtagtg<-log(dn31.codon.pair$gtagtg/cps.gta.gtg)
cps.obs.exp.gtagtt<-log(dn31.codon.pair$gtagtt/cps.gta.gtt)

#cps.obs.exp.gtataa<-log(dn31.codon.pair$gtataa/cps.gta.taa)
cps.obs.exp.gtatac<-log(dn31.codon.pair$gtatac/cps.gta.tac)
#cps.obs.exp.gtatag<-log(dn31.codon.pair$gtatag/cps.gta.tag)
cps.obs.exp.gtatat<-log(dn31.codon.pair$gtatat/cps.gta.tat)
cps.obs.exp.gtatca<-log(dn31.codon.pair$gtatca/cps.gta.tca)
cps.obs.exp.gtatcc<-log(dn31.codon.pair$gtatcc/cps.gta.tcc)
cps.obs.exp.gtatcg<-log(dn31.codon.pair$gtatcg/cps.gta.tcg)
cps.obs.exp.gtatct<-log(dn31.codon.pair$gtatct/cps.gta.tct)
#cps.obs.exp.gtatga<-log(dn31.codon.pair$gtatga/cps.gta.tga)
cps.obs.exp.gtatgc<-log(dn31.codon.pair$gtatgc/cps.gta.tgc)
cps.obs.exp.gtatgg<-log(dn31.codon.pair$gtatgg/cps.gta.tgg)
cps.obs.exp.gtatgt<-log(dn31.codon.pair$gtatgt/cps.gta.tgt)
cps.obs.exp.gtatta<-log(dn31.codon.pair$gtatta/cps.gta.tta)
cps.obs.exp.gtattc<-log(dn31.codon.pair$gtattc/cps.gta.ttc)
cps.obs.exp.gtattg<-log(dn31.codon.pair$gtattg/cps.gta.ttg)
cps.obs.exp.gtattt<-log(dn31.codon.pair$gtattt/cps.gta.ttt)









cps.obs.exp.gtcaaa<-log(dn31.codon.pair$gtcaaa/cps.gtc.aaa)
cps.obs.exp.gtcaac<-log(dn31.codon.pair$gtcaac/cps.gtc.aac)
cps.obs.exp.gtcaag<-log(dn31.codon.pair$gtcaag/cps.gtc.aag)
cps.obs.exp.gtcaat<-log(dn31.codon.pair$gtcaat/cps.gtc.aat)
cps.obs.exp.gtcaca<-log(dn31.codon.pair$gtcaca/cps.gtc.aca)
cps.obs.exp.gtcacc<-log(dn31.codon.pair$gtcacc/cps.gtc.acc)
cps.obs.exp.gtcacg<-log(dn31.codon.pair$gtcacg/cps.gtc.acg)
cps.obs.exp.gtcact<-log(dn31.codon.pair$gtcact/cps.gtc.act)
cps.obs.exp.gtcaga<-log(dn31.codon.pair$gtcaga/cps.gtc.aga)
cps.obs.exp.gtcagc<-log(dn31.codon.pair$gtcagc/cps.gtc.agc)
cps.obs.exp.gtcagg<-log(dn31.codon.pair$gtcagg/cps.gtc.agg)
cps.obs.exp.gtcagt<-log(dn31.codon.pair$gtcagt/cps.gtc.agt)
cps.obs.exp.gtcata<-log(dn31.codon.pair$gtcata/cps.gtc.ata)
cps.obs.exp.gtcatc<-log(dn31.codon.pair$gtcatc/cps.gtc.atc)
cps.obs.exp.gtcatg<-log(dn31.codon.pair$gtcatg/cps.gtc.atg)
cps.obs.exp.gtcatt<-log(dn31.codon.pair$gtcatt/cps.gtc.att)

cps.obs.exp.gtccaa<-log(dn31.codon.pair$gtccaa/cps.gtc.caa)
cps.obs.exp.gtccac<-log(dn31.codon.pair$gtccac/cps.gtc.cac)
cps.obs.exp.gtccag<-log(dn31.codon.pair$gtccag/cps.gtc.cag)
cps.obs.exp.gtccat<-log(dn31.codon.pair$gtccat/cps.gtc.cat)
cps.obs.exp.gtccca<-log(dn31.codon.pair$gtccca/cps.gtc.cca)
cps.obs.exp.gtcccc<-log(dn31.codon.pair$gtcccc/cps.gtc.ccc)
cps.obs.exp.gtcccg<-log(dn31.codon.pair$gtcccg/cps.gtc.ccg)
cps.obs.exp.gtccct<-log(dn31.codon.pair$gtccct/cps.gtc.cct)
cps.obs.exp.gtccga<-log(dn31.codon.pair$gtccga/cps.gtc.cga)
cps.obs.exp.gtccgc<-log(dn31.codon.pair$gtccgc/cps.gtc.cgc)
cps.obs.exp.gtccgg<-log(dn31.codon.pair$gtccgg/cps.gtc.cgg)
cps.obs.exp.gtccgt<-log(dn31.codon.pair$gtccgt/cps.gtc.cgt)
cps.obs.exp.gtccta<-log(dn31.codon.pair$gtccta/cps.gtc.cta)
cps.obs.exp.gtcctc<-log(dn31.codon.pair$gtcctc/cps.gtc.ctc)
cps.obs.exp.gtcctg<-log(dn31.codon.pair$gtcctg/cps.gtc.ctg)
cps.obs.exp.gtcctt<-log(dn31.codon.pair$gtcctt/cps.gtc.ctt)

cps.obs.exp.gtcgaa<-log(dn31.codon.pair$gtcgaa/cps.gtc.gaa)
cps.obs.exp.gtcgac<-log(dn31.codon.pair$gtcgac/cps.gtc.gac)
cps.obs.exp.gtcgag<-log(dn31.codon.pair$gtcgag/cps.gtc.gag)
cps.obs.exp.gtcgat<-log(dn31.codon.pair$gtcgat/cps.gtc.gat)
cps.obs.exp.gtcgca<-log(dn31.codon.pair$gtcgca/cps.gtc.gca)
cps.obs.exp.gtcgcc<-log(dn31.codon.pair$gtcgcc/cps.gtc.gcc)
cps.obs.exp.gtcgcg<-log(dn31.codon.pair$gtcgcg/cps.gtc.gcg)
cps.obs.exp.gtcgct<-log(dn31.codon.pair$gtcgct/cps.gtc.gct)
cps.obs.exp.gtcgga<-log(dn31.codon.pair$gtcgga/cps.gtc.gga)
cps.obs.exp.gtcggc<-log(dn31.codon.pair$gtcggc/cps.gtc.ggc)
cps.obs.exp.gtcggg<-log(dn31.codon.pair$gtcggg/cps.gtc.ggg)
cps.obs.exp.gtcggt<-log(dn31.codon.pair$gtcggt/cps.gtc.ggt)
cps.obs.exp.gtcgta<-log(dn31.codon.pair$gtcgta/cps.gtc.gta)
cps.obs.exp.gtcgtc<-log(dn31.codon.pair$gtcgtc/cps.gtc.gtc)
cps.obs.exp.gtcgtg<-log(dn31.codon.pair$gtcgtg/cps.gtc.gtg)
cps.obs.exp.gtcgtt<-log(dn31.codon.pair$gtcgtt/cps.gtc.gtt)

#cps.obs.exp.gtctaa<-log(dn31.codon.pair$gtctaa/cps.gtc.taa)
cps.obs.exp.gtctac<-log(dn31.codon.pair$gtctac/cps.gtc.tac)
#cps.obs.exp.gtctag<-log(dn31.codon.pair$gtctag/cps.gtc.tag)
cps.obs.exp.gtctat<-log(dn31.codon.pair$gtctat/cps.gtc.tat)
cps.obs.exp.gtctca<-log(dn31.codon.pair$gtctca/cps.gtc.tca)
cps.obs.exp.gtctcc<-log(dn31.codon.pair$gtctcc/cps.gtc.tcc)
cps.obs.exp.gtctcg<-log(dn31.codon.pair$gtctcg/cps.gtc.tcg)
cps.obs.exp.gtctct<-log(dn31.codon.pair$gtctct/cps.gtc.tct)
#cps.obs.exp.gtctga<-log(dn31.codon.pair$gtctga/cps.gtc.tga)
cps.obs.exp.gtctgc<-log(dn31.codon.pair$gtctgc/cps.gtc.tgc)
cps.obs.exp.gtctgg<-log(dn31.codon.pair$gtctgg/cps.gtc.tgg)
cps.obs.exp.gtctgt<-log(dn31.codon.pair$gtctgt/cps.gtc.tgt)
cps.obs.exp.gtctta<-log(dn31.codon.pair$gtctta/cps.gtc.tta)
cps.obs.exp.gtcttc<-log(dn31.codon.pair$gtcttc/cps.gtc.ttc)
cps.obs.exp.gtcttg<-log(dn31.codon.pair$gtcttg/cps.gtc.ttg)
cps.obs.exp.gtcttt<-log(dn31.codon.pair$gtcttt/cps.gtc.ttt)











cps.obs.exp.gtgaaa<-log(dn31.codon.pair$gtgaaa/cps.gtg.aaa)
cps.obs.exp.gtgaac<-log(dn31.codon.pair$gtgaac/cps.gtg.aac)
cps.obs.exp.gtgaag<-log(dn31.codon.pair$gtgaag/cps.gtg.aag)
cps.obs.exp.gtgaat<-log(dn31.codon.pair$gtgaat/cps.gtg.aat)
cps.obs.exp.gtgaca<-log(dn31.codon.pair$gtgaca/cps.gtg.aca)
cps.obs.exp.gtgacc<-log(dn31.codon.pair$gtgacc/cps.gtg.acc)
cps.obs.exp.gtgacg<-log(dn31.codon.pair$gtgacg/cps.gtg.acg)
cps.obs.exp.gtgact<-log(dn31.codon.pair$gtgact/cps.gtg.act)
cps.obs.exp.gtgaga<-log(dn31.codon.pair$gtgaga/cps.gtg.aga)
cps.obs.exp.gtgagc<-log(dn31.codon.pair$gtgagc/cps.gtg.agc)
cps.obs.exp.gtgagg<-log(dn31.codon.pair$gtgagg/cps.gtg.agg)
cps.obs.exp.gtgagt<-log(dn31.codon.pair$gtgagt/cps.gtg.agt)
cps.obs.exp.gtgata<-log(dn31.codon.pair$gtgata/cps.gtg.ata)
cps.obs.exp.gtgatc<-log(dn31.codon.pair$gtgatc/cps.gtg.atc)
cps.obs.exp.gtgatg<-log(dn31.codon.pair$gtgatg/cps.gtg.atg)
cps.obs.exp.gtgatt<-log(dn31.codon.pair$gtgatt/cps.gtg.att)

cps.obs.exp.gtgcaa<-log(dn31.codon.pair$gtgcaa/cps.gtg.caa)
cps.obs.exp.gtgcac<-log(dn31.codon.pair$gtgcac/cps.gtg.cac)
cps.obs.exp.gtgcag<-log(dn31.codon.pair$gtgcag/cps.gtg.cag)
cps.obs.exp.gtgcat<-log(dn31.codon.pair$gtgcat/cps.gtg.cat)
cps.obs.exp.gtgcca<-log(dn31.codon.pair$gtgcca/cps.gtg.cca)
cps.obs.exp.gtgccc<-log(dn31.codon.pair$gtgccc/cps.gtg.ccc)
cps.obs.exp.gtgccg<-log(dn31.codon.pair$gtgccg/cps.gtg.ccg)
cps.obs.exp.gtgcct<-log(dn31.codon.pair$gtgcct/cps.gtg.cct)
cps.obs.exp.gtgcga<-log(dn31.codon.pair$gtgcga/cps.gtg.cga)
cps.obs.exp.gtgcgc<-log(dn31.codon.pair$gtgcgc/cps.gtg.cgc)
cps.obs.exp.gtgcgg<-log(dn31.codon.pair$gtgcgg/cps.gtg.cgg)
cps.obs.exp.gtgcgt<-log(dn31.codon.pair$gtgcgt/cps.gtg.cgt)
cps.obs.exp.gtgcta<-log(dn31.codon.pair$gtgcta/cps.gtg.cta)
cps.obs.exp.gtgctc<-log(dn31.codon.pair$gtgctc/cps.gtg.ctc)
cps.obs.exp.gtgctg<-log(dn31.codon.pair$gtgctg/cps.gtg.ctg)
cps.obs.exp.gtgctt<-log(dn31.codon.pair$gtgctt/cps.gtg.ctt)

cps.obs.exp.gtggaa<-log(dn31.codon.pair$gtggaa/cps.gtg.gaa)
cps.obs.exp.gtggac<-log(dn31.codon.pair$gtggac/cps.gtg.gac)
cps.obs.exp.gtggag<-log(dn31.codon.pair$gtggag/cps.gtg.gag)
cps.obs.exp.gtggat<-log(dn31.codon.pair$gtggat/cps.gtg.gat)
cps.obs.exp.gtggca<-log(dn31.codon.pair$gtggca/cps.gtg.gca)
cps.obs.exp.gtggcc<-log(dn31.codon.pair$gtggcc/cps.gtg.gcc)
cps.obs.exp.gtggcg<-log(dn31.codon.pair$gtggcg/cps.gtg.gcg)
cps.obs.exp.gtggct<-log(dn31.codon.pair$gtggct/cps.gtg.gct)
cps.obs.exp.gtggga<-log(dn31.codon.pair$gtggga/cps.gtg.gga)
cps.obs.exp.gtgggc<-log(dn31.codon.pair$gtgggc/cps.gtg.ggc)
cps.obs.exp.gtgggg<-log(dn31.codon.pair$gtgggg/cps.gtg.ggg)
cps.obs.exp.gtgggt<-log(dn31.codon.pair$gtgggt/cps.gtg.ggt)
cps.obs.exp.gtggta<-log(dn31.codon.pair$gtggta/cps.gtg.gta)
cps.obs.exp.gtggtc<-log(dn31.codon.pair$gtggtc/cps.gtg.gtc)
cps.obs.exp.gtggtg<-log(dn31.codon.pair$gtggtg/cps.gtg.gtg)
cps.obs.exp.gtggtt<-log(dn31.codon.pair$gtggtt/cps.gtg.gtt)

#cps.obs.exp.gtgtaa<-log(dn31.codon.pair$gtgtaa/cps.gtg.taa)
cps.obs.exp.gtgtac<-log(dn31.codon.pair$gtgtac/cps.gtg.tac)
#cps.obs.exp.gtgtag<-log(dn31.codon.pair$gtgtag/cps.gtg.tag)
cps.obs.exp.gtgtat<-log(dn31.codon.pair$gtgtat/cps.gtg.tat)
cps.obs.exp.gtgtca<-log(dn31.codon.pair$gtgtca/cps.gtg.tca)
cps.obs.exp.gtgtcc<-log(dn31.codon.pair$gtgtcc/cps.gtg.tcc)
cps.obs.exp.gtgtcg<-log(dn31.codon.pair$gtgtcg/cps.gtg.tcg)
cps.obs.exp.gtgtct<-log(dn31.codon.pair$gtgtct/cps.gtg.tct)
#cps.obs.exp.gtgtga<-log(dn31.codon.pair$gtgtga/cps.gtg.tga)
cps.obs.exp.gtgtgc<-log(dn31.codon.pair$gtgtgc/cps.gtg.tgc)
cps.obs.exp.gtgtgg<-log(dn31.codon.pair$gtgtgg/cps.gtg.tgg)
cps.obs.exp.gtgtgt<-log(dn31.codon.pair$gtgtgt/cps.gtg.tgt)
cps.obs.exp.gtgtta<-log(dn31.codon.pair$gtgtta/cps.gtg.tta)
cps.obs.exp.gtgttc<-log(dn31.codon.pair$gtgttc/cps.gtg.ttc)
cps.obs.exp.gtgttg<-log(dn31.codon.pair$gtgttg/cps.gtg.ttg)
cps.obs.exp.gtgttt<-log(dn31.codon.pair$gtgttt/cps.gtg.ttt)








cps.obs.exp.gttaaa<-log(dn31.codon.pair$gttaaa/cps.gtt.aaa)
cps.obs.exp.gttaac<-log(dn31.codon.pair$gttaac/cps.gtt.aac)
cps.obs.exp.gttaag<-log(dn31.codon.pair$gttaag/cps.gtt.aag)
cps.obs.exp.gttaat<-log(dn31.codon.pair$gttaat/cps.gtt.aat)
cps.obs.exp.gttaca<-log(dn31.codon.pair$gttaca/cps.gtt.aca)
cps.obs.exp.gttacc<-log(dn31.codon.pair$gttacc/cps.gtt.acc)
cps.obs.exp.gttacg<-log(dn31.codon.pair$gttacg/cps.gtt.acg)
cps.obs.exp.gttact<-log(dn31.codon.pair$gttact/cps.gtt.act)
cps.obs.exp.gttaga<-log(dn31.codon.pair$gttaga/cps.gtt.aga)
cps.obs.exp.gttagc<-log(dn31.codon.pair$gttagc/cps.gtt.agc)
cps.obs.exp.gttagg<-log(dn31.codon.pair$gttagg/cps.gtt.agg)
cps.obs.exp.gttagt<-log(dn31.codon.pair$gttagt/cps.gtt.agt)
cps.obs.exp.gttata<-log(dn31.codon.pair$gttata/cps.gtt.ata)
cps.obs.exp.gttatc<-log(dn31.codon.pair$gttatc/cps.gtt.atc)
cps.obs.exp.gttatg<-log(dn31.codon.pair$gttatg/cps.gtt.atg)
cps.obs.exp.gttatt<-log(dn31.codon.pair$gttatt/cps.gtt.att)

cps.obs.exp.gttcaa<-log(dn31.codon.pair$gttcaa/cps.gtt.caa)
cps.obs.exp.gttcac<-log(dn31.codon.pair$gttcac/cps.gtt.cac)
cps.obs.exp.gttcag<-log(dn31.codon.pair$gttcag/cps.gtt.cag)
cps.obs.exp.gttcat<-log(dn31.codon.pair$gttcat/cps.gtt.cat)
cps.obs.exp.gttcca<-log(dn31.codon.pair$gttcca/cps.gtt.cca)
cps.obs.exp.gttccc<-log(dn31.codon.pair$gttccc/cps.gtt.ccc)
cps.obs.exp.gttccg<-log(dn31.codon.pair$gttccg/cps.gtt.ccg)
cps.obs.exp.gttcct<-log(dn31.codon.pair$gttcct/cps.gtt.cct)
cps.obs.exp.gttcga<-log(dn31.codon.pair$gttcga/cps.gtt.cga)
cps.obs.exp.gttcgc<-log(dn31.codon.pair$gttcgc/cps.gtt.cgc)
cps.obs.exp.gttcgg<-log(dn31.codon.pair$gttcgg/cps.gtt.cgg)
cps.obs.exp.gttcgt<-log(dn31.codon.pair$gttcgt/cps.gtt.cgt)
cps.obs.exp.gttcta<-log(dn31.codon.pair$gttcta/cps.gtt.cta)
cps.obs.exp.gttctc<-log(dn31.codon.pair$gttctc/cps.gtt.ctc)
cps.obs.exp.gttctg<-log(dn31.codon.pair$gttctg/cps.gtt.ctg)
cps.obs.exp.gttctt<-log(dn31.codon.pair$gttctt/cps.gtt.ctt)

cps.obs.exp.gttgaa<-log(dn31.codon.pair$gttgaa/cps.gtt.gaa)
cps.obs.exp.gttgac<-log(dn31.codon.pair$gttgac/cps.gtt.gac)
cps.obs.exp.gttgag<-log(dn31.codon.pair$gttgag/cps.gtt.gag)
cps.obs.exp.gttgat<-log(dn31.codon.pair$gttgat/cps.gtt.gat)
cps.obs.exp.gttgca<-log(dn31.codon.pair$gttgca/cps.gtt.gca)
cps.obs.exp.gttgcc<-log(dn31.codon.pair$gttgcc/cps.gtt.gcc)
cps.obs.exp.gttgcg<-log(dn31.codon.pair$gttgcg/cps.gtt.gcg)
cps.obs.exp.gttgct<-log(dn31.codon.pair$gttgct/cps.gtt.gct)
cps.obs.exp.gttgga<-log(dn31.codon.pair$gttgga/cps.gtt.gga)
cps.obs.exp.gttggc<-log(dn31.codon.pair$gttggc/cps.gtt.ggc)
cps.obs.exp.gttggg<-log(dn31.codon.pair$gttggg/cps.gtt.ggg)
cps.obs.exp.gttggt<-log(dn31.codon.pair$gttggt/cps.gtt.ggt)
cps.obs.exp.gttgta<-log(dn31.codon.pair$gttgta/cps.gtt.gta)
cps.obs.exp.gttgtc<-log(dn31.codon.pair$gttgtc/cps.gtt.gtc)
cps.obs.exp.gttgtg<-log(dn31.codon.pair$gttgtg/cps.gtt.gtg)
cps.obs.exp.gttgtt<-log(dn31.codon.pair$gttgtt/cps.gtt.gtt)

#cps.obs.exp.gtttaa<-log(dn31.codon.pair$gtttaa/cps.gtt.taa)
cps.obs.exp.gtttac<-log(dn31.codon.pair$gtttac/cps.gtt.tac)
#cps.obs.exp.gtttag<-log(dn31.codon.pair$gtttag/cps.gtt.tag)
cps.obs.exp.gtttat<-log(dn31.codon.pair$gtttat/cps.gtt.tat)
cps.obs.exp.gtttca<-log(dn31.codon.pair$gtttca/cps.gtt.tca)
cps.obs.exp.gtttcc<-log(dn31.codon.pair$gtttcc/cps.gtt.tcc)
cps.obs.exp.gtttcg<-log(dn31.codon.pair$gtttcg/cps.gtt.tcg)
cps.obs.exp.gtttct<-log(dn31.codon.pair$gtttct/cps.gtt.tct)
#cps.obs.exp.gtttga<-log(dn31.codon.pair$gtttga/cps.gtt.tga)
cps.obs.exp.gtttgc<-log(dn31.codon.pair$gtttgc/cps.gtt.tgc)
cps.obs.exp.gtttgg<-log(dn31.codon.pair$gtttgg/cps.gtt.tgg)
cps.obs.exp.gtttgt<-log(dn31.codon.pair$gtttgt/cps.gtt.tgt)
cps.obs.exp.gtttta<-log(dn31.codon.pair$gtttta/cps.gtt.tta)
cps.obs.exp.gttttc<-log(dn31.codon.pair$gttttc/cps.gtt.ttc)
cps.obs.exp.gttttg<-log(dn31.codon.pair$gttttg/cps.gtt.ttg)
cps.obs.exp.gttttt<-log(dn31.codon.pair$gttttt/cps.gtt.ttt)

















#Stop codon


#cps.obs.exp.taaaaa<-log(dn31.codon.pair$taaaaa/cps.taa.aaa)
#cps.obs.exp.taaaac<-log(dn31.codon.pair$taaaac/cps.taa.aac)
#cps.obs.exp.taaaag<-log(dn31.codon.pair$taaaag/cps.taa.aag)
#cps.obs.exp.taaaat<-log(dn31.codon.pair$taaaat/cps.taa.aat)
#cps.obs.exp.taaaca<-log(dn31.codon.pair$taaaca/cps.taa.aca)
#cps.obs.exp.taaacc<-log(dn31.codon.pair$taaacc/cps.taa.acc)
#cps.obs.exp.taaacg<-log(dn31.codon.pair$taaacg/cps.taa.acg)
#cps.obs.exp.taaact<-log(dn31.codon.pair$taaact/cps.taa.act)
#cps.obs.exp.taaaga<-log(dn31.codon.pair$taaaga/cps.taa.aga)
#cps.obs.exp.taaagc<-log(dn31.codon.pair$taaagc/cps.taa.agc)
#cps.obs.exp.taaagg<-log(dn31.codon.pair$taaagg/cps.taa.agg)
#cps.obs.exp.taaagt<-log(dn31.codon.pair$taaagt/cps.taa.agt)
#cps.obs.exp.taaata<-log(dn31.codon.pair$taaata/cps.taa.ata)
#cps.obs.exp.taaatc<-log(dn31.codon.pair$taaatc/cps.taa.atc)
#cps.obs.exp.taaatg<-log(dn31.codon.pair$taaatg/cps.taa.atg)
#cps.obs.exp.taaatt<-log(dn31.codon.pair$taaatt/cps.taa.att)

#cps.obs.exp.taacaa<-log(dn31.codon.pair$taacaa/cps.taa.caa)
#cps.obs.exp.taacac<-log(dn31.codon.pair$taacac/cps.taa.cac)
#cps.obs.exp.taacag<-log(dn31.codon.pair$taacag/cps.taa.cag)
#cps.obs.exp.taacat<-log(dn31.codon.pair$taacat/cps.taa.cat)
#cps.obs.exp.taacca<-log(dn31.codon.pair$taacca/cps.taa.cca)
#cps.obs.exp.taaccc<-log(dn31.codon.pair$taaccc/cps.taa.ccc)
#cps.obs.exp.taaccg<-log(dn31.codon.pair$taaccg/cps.taa.ccg)
#cps.obs.exp.taacct<-log(dn31.codon.pair$taacct/cps.taa.cct)
#cps.obs.exp.taacga<-log(dn31.codon.pair$taacga/cps.taa.cga)
#cps.obs.exp.taacgc<-log(dn31.codon.pair$taacgc/cps.taa.cgc)
#cps.obs.exp.taacgg<-log(dn31.codon.pair$taacgg/cps.taa.cgg)
#cps.obs.exp.taacgt<-log(dn31.codon.pair$taacgt/cps.taa.cgt)
#cps.obs.exp.taacta<-log(dn31.codon.pair$taacta/cps.taa.cta)
#cps.obs.exp.taactc<-log(dn31.codon.pair$taactc/cps.taa.ctc)
#cps.obs.exp.taactg<-log(dn31.codon.pair$taactg/cps.taa.ctg)
#cps.obs.exp.taactt<-log(dn31.codon.pair$taactt/cps.taa.ctt)

#cps.obs.exp.taagaa<-log(dn31.codon.pair$taagaa/cps.taa.gaa)
#cps.obs.exp.taagac<-log(dn31.codon.pair$taagac/cps.taa.gac)
#cps.obs.exp.taagag<-log(dn31.codon.pair$taagag/cps.taa.gag)
#cps.obs.exp.taagat<-log(dn31.codon.pair$taagat/cps.taa.gat)
#cps.obs.exp.taagca<-log(dn31.codon.pair$taagca/cps.taa.gca)
#cps.obs.exp.taagcc<-log(dn31.codon.pair$taagcc/cps.taa.gcc)
#cps.obs.exp.taagcg<-log(dn31.codon.pair$taagcg/cps.taa.gcg)
#cps.obs.exp.taagct<-log(dn31.codon.pair$taagct/cps.taa.gct)
#cps.obs.exp.taagga<-log(dn31.codon.pair$taagga/cps.taa.gga)
#cps.obs.exp.taaggc<-log(dn31.codon.pair$taaggc/cps.taa.ggc)
#cps.obs.exp.taaggg<-log(dn31.codon.pair$taaggg/cps.taa.ggg)
#cps.obs.exp.taaggt<-log(dn31.codon.pair$taaggt/cps.taa.ggt)
#cps.obs.exp.taagta<-log(dn31.codon.pair$taagta/cps.taa.gta)
#cps.obs.exp.taagtc<-log(dn31.codon.pair$taagtc/cps.taa.gtc)
#cps.obs.exp.taagtg<-log(dn31.codon.pair$taagtg/cps.taa.gtg)
#cps.obs.exp.taagtt<-log(dn31.codon.pair$taagtt/cps.taa.gtt)

#cps.obs.exp.taataa<-log(dn31.codon.pair$taataa/cps.taa.taa)
#cps.obs.exp.taatac<-log(dn31.codon.pair$taatac/cps.taa.tac)
#cps.obs.exp.taatag<-log(dn31.codon.pair$taatag/cps.taa.tag)
#cps.obs.exp.taatat<-log(dn31.codon.pair$taatat/cps.taa.tat)
#cps.obs.exp.taatca<-log(dn31.codon.pair$taatca/cps.taa.tca)
#cps.obs.exp.taatcc<-log(dn31.codon.pair$taatcc/cps.taa.tcc)
#cps.obs.exp.taatcg<-log(dn31.codon.pair$taatcg/cps.taa.tcg)
#cps.obs.exp.taatct<-log(dn31.codon.pair$taatct/cps.taa.tct)
#cps.obs.exp.taatga<-log(dn31.codon.pair$taatga/cps.taa.tga)
#cps.obs.exp.taatgc<-log(dn31.codon.pair$taatgc/cps.taa.tgc)
#cps.obs.exp.taatgg<-log(dn31.codon.pair$taatgg/cps.taa.tgg)
#cps.obs.exp.taatgt<-log(dn31.codon.pair$taatgt/cps.taa.tgt)
#cps.obs.exp.taatta<-log(dn31.codon.pair$taatta/cps.taa.tta)
#cps.obs.exp.taattc<-log(dn31.codon.pair$taattc/cps.taa.ttc)
#cps.obs.exp.taattg<-log(dn31.codon.pair$taattg/cps.taa.ttg)
#cps.obs.exp.taattt<-log(dn31.codon.pair$taattt/cps.taa.ttt)









cps.obs.exp.tacaaa<-log(dn31.codon.pair$tacaaa/cps.tac.aaa)
cps.obs.exp.tacaac<-log(dn31.codon.pair$tacaac/cps.tac.aac)
cps.obs.exp.tacaag<-log(dn31.codon.pair$tacaag/cps.tac.aag)
cps.obs.exp.tacaat<-log(dn31.codon.pair$tacaat/cps.tac.aat)
cps.obs.exp.tacaca<-log(dn31.codon.pair$tacaca/cps.tac.aca)
cps.obs.exp.tacacc<-log(dn31.codon.pair$tacacc/cps.tac.acc)
cps.obs.exp.tacacg<-log(dn31.codon.pair$tacacg/cps.tac.acg)
cps.obs.exp.tacact<-log(dn31.codon.pair$tacact/cps.tac.act)
cps.obs.exp.tacaga<-log(dn31.codon.pair$tacaga/cps.tac.aga)
cps.obs.exp.tacagc<-log(dn31.codon.pair$tacagc/cps.tac.agc)
cps.obs.exp.tacagg<-log(dn31.codon.pair$tacagg/cps.tac.agg)
cps.obs.exp.tacagt<-log(dn31.codon.pair$tacagt/cps.tac.agt)
cps.obs.exp.tacata<-log(dn31.codon.pair$tacata/cps.tac.ata)
cps.obs.exp.tacatc<-log(dn31.codon.pair$tacatc/cps.tac.atc)
cps.obs.exp.tacatg<-log(dn31.codon.pair$tacatg/cps.tac.atg)
cps.obs.exp.tacatt<-log(dn31.codon.pair$tacatt/cps.tac.att)

cps.obs.exp.taccaa<-log(dn31.codon.pair$taccaa/cps.tac.caa)
cps.obs.exp.taccac<-log(dn31.codon.pair$taccac/cps.tac.cac)
cps.obs.exp.taccag<-log(dn31.codon.pair$taccag/cps.tac.cag)
cps.obs.exp.taccat<-log(dn31.codon.pair$taccat/cps.tac.cat)
cps.obs.exp.taccca<-log(dn31.codon.pair$taccca/cps.tac.cca)
cps.obs.exp.tacccc<-log(dn31.codon.pair$tacccc/cps.tac.ccc)
cps.obs.exp.tacccg<-log(dn31.codon.pair$tacccg/cps.tac.ccg)
cps.obs.exp.taccct<-log(dn31.codon.pair$taccct/cps.tac.cct)
cps.obs.exp.taccga<-log(dn31.codon.pair$taccga/cps.tac.cga)
cps.obs.exp.taccgc<-log(dn31.codon.pair$taccgc/cps.tac.cgc)
cps.obs.exp.taccgg<-log(dn31.codon.pair$taccgg/cps.tac.cgg)
cps.obs.exp.taccgt<-log(dn31.codon.pair$taccgt/cps.tac.cgt)
cps.obs.exp.taccta<-log(dn31.codon.pair$taccta/cps.tac.cta)
cps.obs.exp.tacctc<-log(dn31.codon.pair$tacctc/cps.tac.ctc)
cps.obs.exp.tacctg<-log(dn31.codon.pair$tacctg/cps.tac.ctg)
cps.obs.exp.tacctt<-log(dn31.codon.pair$tacctt/cps.tac.ctt)

cps.obs.exp.tacgaa<-log(dn31.codon.pair$tacgaa/cps.tac.gaa)
cps.obs.exp.tacgac<-log(dn31.codon.pair$tacgac/cps.tac.gac)
cps.obs.exp.tacgag<-log(dn31.codon.pair$tacgag/cps.tac.gag)
cps.obs.exp.tacgat<-log(dn31.codon.pair$tacgat/cps.tac.gat)
cps.obs.exp.tacgca<-log(dn31.codon.pair$tacgca/cps.tac.gca)
cps.obs.exp.tacgcc<-log(dn31.codon.pair$tacgcc/cps.tac.gcc)
cps.obs.exp.tacgcg<-log(dn31.codon.pair$tacgcg/cps.tac.gcg)
cps.obs.exp.tacgct<-log(dn31.codon.pair$tacgct/cps.tac.gct)
cps.obs.exp.tacgga<-log(dn31.codon.pair$tacgga/cps.tac.gga)
cps.obs.exp.tacggc<-log(dn31.codon.pair$tacggc/cps.tac.ggc)
cps.obs.exp.tacggg<-log(dn31.codon.pair$tacggg/cps.tac.ggg)
cps.obs.exp.tacggt<-log(dn31.codon.pair$tacggt/cps.tac.ggt)
cps.obs.exp.tacgta<-log(dn31.codon.pair$tacgta/cps.tac.gta)
cps.obs.exp.tacgtc<-log(dn31.codon.pair$tacgtc/cps.tac.gtc)
cps.obs.exp.tacgtg<-log(dn31.codon.pair$tacgtg/cps.tac.gtg)
cps.obs.exp.tacgtt<-log(dn31.codon.pair$tacgtt/cps.tac.gtt)

#cps.obs.exp.tactaa<-log(dn31.codon.pair$tactaa/cps.tac.taa)
cps.obs.exp.tactac<-log(dn31.codon.pair$tactac/cps.tac.tac)
#cps.obs.exp.tactag<-log(dn31.codon.pair$tactag/cps.tac.tag)
cps.obs.exp.tactat<-log(dn31.codon.pair$tactat/cps.tac.tat)
cps.obs.exp.tactca<-log(dn31.codon.pair$tactca/cps.tac.tca)
cps.obs.exp.tactcc<-log(dn31.codon.pair$tactcc/cps.tac.tcc)
cps.obs.exp.tactcg<-log(dn31.codon.pair$tactcg/cps.tac.tcg)
cps.obs.exp.tactct<-log(dn31.codon.pair$tactct/cps.tac.tct)
#cps.obs.exp.tactga<-log(dn31.codon.pair$tactga/cps.tac.tga)
cps.obs.exp.tactgc<-log(dn31.codon.pair$tactgc/cps.tac.tgc)
cps.obs.exp.tactgg<-log(dn31.codon.pair$tactgg/cps.tac.tgg)
cps.obs.exp.tactgt<-log(dn31.codon.pair$tactgt/cps.tac.tgt)
cps.obs.exp.tactta<-log(dn31.codon.pair$tactta/cps.tac.tta)
cps.obs.exp.tacttc<-log(dn31.codon.pair$tacttc/cps.tac.ttc)
cps.obs.exp.tacttg<-log(dn31.codon.pair$tacttg/cps.tac.ttg)
cps.obs.exp.tacttt<-log(dn31.codon.pair$tacttt/cps.tac.ttt)







#Stop Codon 


#cps.obs.exp.tagaaa<-log(dn31.codon.pair$tagaaa/cps.tag.aaa)
#cps.obs.exp.tagaac<-log(dn31.codon.pair$tagaac/cps.tag.aac)
#cps.obs.exp.tagaag<-log(dn31.codon.pair$tagaag/cps.tag.aag)
#cps.obs.exp.tagaat<-log(dn31.codon.pair$tagaat/cps.tag.aat)
#cps.obs.exp.tagaca<-log(dn31.codon.pair$tagaca/cps.tag.aca)
#cps.obs.exp.tagacc<-log(dn31.codon.pair$tagacc/cps.tag.acc)
#cps.obs.exp.tagacg<-log(dn31.codon.pair$tagacg/cps.tag.acg)
#cps.obs.exp.tagact<-log(dn31.codon.pair$tagact/cps.tag.act)
#cps.obs.exp.tagaga<-log(dn31.codon.pair$tagaga/cps.tag.aga)
#cps.obs.exp.tagagc<-log(dn31.codon.pair$tagagc/cps.tag.agc)
#cps.obs.exp.tagagg<-log(dn31.codon.pair$tagagg/cps.tag.agg)
#cps.obs.exp.tagagt<-log(dn31.codon.pair$tagagt/cps.tag.agt)
#cps.obs.exp.tagata<-log(dn31.codon.pair$tagata/cps.tag.ata)
#cps.obs.exp.tagatc<-log(dn31.codon.pair$tagatc/cps.tag.atc)
#cps.obs.exp.tagatg<-log(dn31.codon.pair$tagatg/cps.tag.atg)
#cps.obs.exp.tagatt<-log(dn31.codon.pair$tagatt/cps.tag.att)

#cps.obs.exp.tagcaa<-log(dn31.codon.pair$tagcaa/cps.tag.caa)
#cps.obs.exp.tagcac<-log(dn31.codon.pair$tagcac/cps.tag.cac)
#cps.obs.exp.tagcag<-log(dn31.codon.pair$tagcag/cps.tag.cag)
#cps.obs.exp.tagcat<-log(dn31.codon.pair$tagcat/cps.tag.cat)
#cps.obs.exp.tagcca<-log(dn31.codon.pair$tagcca/cps.tag.cca)
#cps.obs.exp.tagccc<-log(dn31.codon.pair$tagccc/cps.tag.ccc)
#cps.obs.exp.tagccg<-log(dn31.codon.pair$tagccg/cps.tag.ccg)
#cps.obs.exp.tagcct<-log(dn31.codon.pair$tagcct/cps.tag.cct)
#cps.obs.exp.tagcga<-log(dn31.codon.pair$tagcga/cps.tag.cga)
#cps.obs.exp.tagcgc<-log(dn31.codon.pair$tagcgc/cps.tag.cgc)
#cps.obs.exp.tagcgg<-log(dn31.codon.pair$tagcgg/cps.tag.cgg)
#cps.obs.exp.tagcgt<-log(dn31.codon.pair$tagcgt/cps.tag.cgt)
#cps.obs.exp.tagcta<-log(dn31.codon.pair$tagcta/cps.tag.cta)
#cps.obs.exp.tagctc<-log(dn31.codon.pair$tagctc/cps.tag.ctc)
#cps.obs.exp.tagctg<-log(dn31.codon.pair$tagctg/cps.tag.ctg)
#cps.obs.exp.tagctt<-log(dn31.codon.pair$tagctt/cps.tag.ctt)

#cps.obs.exp.taggaa<-log(dn31.codon.pair$taggaa/cps.tag.gaa)
#cps.obs.exp.taggac<-log(dn31.codon.pair$taggac/cps.tag.gac)
#cps.obs.exp.taggag<-log(dn31.codon.pair$taggag/cps.tag.gag)
#cps.obs.exp.taggat<-log(dn31.codon.pair$taggat/cps.tag.gat)
#cps.obs.exp.taggca<-log(dn31.codon.pair$taggca/cps.tag.gca)
#cps.obs.exp.taggcc<-log(dn31.codon.pair$taggcc/cps.tag.gcc)
#cps.obs.exp.taggcg<-log(dn31.codon.pair$taggcg/cps.tag.gcg)
#cps.obs.exp.taggct<-log(dn31.codon.pair$taggct/cps.tag.gct)
#cps.obs.exp.taggga<-log(dn31.codon.pair$taggga/cps.tag.gga)
#cps.obs.exp.tagggc<-log(dn31.codon.pair$tagggc/cps.tag.ggc)
#cps.obs.exp.tagggg<-log(dn31.codon.pair$tagggg/cps.tag.ggg)
#cps.obs.exp.tagggt<-log(dn31.codon.pair$tagggt/cps.tag.ggt)
#cps.obs.exp.taggta<-log(dn31.codon.pair$taggta/cps.tag.gta)
#cps.obs.exp.taggtc<-log(dn31.codon.pair$taggtc/cps.tag.gtc)
#cps.obs.exp.taggtg<-log(dn31.codon.pair$taggtg/cps.tag.gtg)
#cps.obs.exp.taggtt<-log(dn31.codon.pair$taggtt/cps.tag.gtt)

#cps.obs.exp.tagtaa<-log(dn31.codon.pair$tagtaa/cps.tag.taa)
#cps.obs.exp.tagtac<-log(dn31.codon.pair$tagtac/cps.tag.tac)
#cps.obs.exp.tagtag<-log(dn31.codon.pair$tagtag/cps.tag.tag)
#cps.obs.exp.tagtat<-log(dn31.codon.pair$tagtat/cps.tag.tat)
#cps.obs.exp.tagtca<-log(dn31.codon.pair$tagtca/cps.tag.tca)
#cps.obs.exp.tagtcc<-log(dn31.codon.pair$tagtcc/cps.tag.tcc)
#cps.obs.exp.tagtcg<-log(dn31.codon.pair$tagtcg/cps.tag.tcg)
#cps.obs.exp.tagtct<-log(dn31.codon.pair$tagtct/cps.tag.tct)
#cps.obs.exp.tagtga<-log(dn31.codon.pair$tagtga/cps.tag.tga)
#cps.obs.exp.tagtgc<-log(dn31.codon.pair$tagtgc/cps.tag.tgc)
#cps.obs.exp.tagtgg<-log(dn31.codon.pair$tagtgg/cps.tag.tgg)
#cps.obs.exp.tagtgt<-log(dn31.codon.pair$tagtgt/cps.tag.tgt)
#cps.obs.exp.tagtta<-log(dn31.codon.pair$tagtta/cps.tag.tta)
#cps.obs.exp.tagttc<-log(dn31.codon.pair$tagttc/cps.tag.ttc)
#cps.obs.exp.tagttg<-log(dn31.codon.pair$tagttg/cps.tag.ttg)
#cps.obs.exp.tagttt<-log(dn31.codon.pair$tagttt/cps.tag.ttt)









cps.obs.exp.tataaa<-log(dn31.codon.pair$tataaa/cps.tat.aaa)
cps.obs.exp.tataac<-log(dn31.codon.pair$tataac/cps.tat.aac)
cps.obs.exp.tataag<-log(dn31.codon.pair$tataag/cps.tat.aag)
cps.obs.exp.tataat<-log(dn31.codon.pair$tataat/cps.tat.aat)
cps.obs.exp.tataca<-log(dn31.codon.pair$tataca/cps.tat.aca)
cps.obs.exp.tatacc<-log(dn31.codon.pair$tatacc/cps.tat.acc)
cps.obs.exp.tatacg<-log(dn31.codon.pair$tatacg/cps.tat.acg)
cps.obs.exp.tatact<-log(dn31.codon.pair$tatact/cps.tat.act)
cps.obs.exp.tataga<-log(dn31.codon.pair$tataga/cps.tat.aga)
cps.obs.exp.tatagc<-log(dn31.codon.pair$tatagc/cps.tat.agc)
cps.obs.exp.tatagg<-log(dn31.codon.pair$tatagg/cps.tat.agg)
cps.obs.exp.tatagt<-log(dn31.codon.pair$tatagt/cps.tat.agt)
cps.obs.exp.tatata<-log(dn31.codon.pair$tatata/cps.tat.ata)
cps.obs.exp.tatatc<-log(dn31.codon.pair$tatatc/cps.tat.atc)
cps.obs.exp.tatatg<-log(dn31.codon.pair$tatatg/cps.tat.atg)
cps.obs.exp.tatatt<-log(dn31.codon.pair$tatatt/cps.tat.att)

cps.obs.exp.tatcaa<-log(dn31.codon.pair$tatcaa/cps.tat.caa)
cps.obs.exp.tatcac<-log(dn31.codon.pair$tatcac/cps.tat.cac)
cps.obs.exp.tatcag<-log(dn31.codon.pair$tatcag/cps.tat.cag)
cps.obs.exp.tatcat<-log(dn31.codon.pair$tatcat/cps.tat.cat)
cps.obs.exp.tatcca<-log(dn31.codon.pair$tatcca/cps.tat.cca)
cps.obs.exp.tatccc<-log(dn31.codon.pair$tatccc/cps.tat.ccc)
cps.obs.exp.tatccg<-log(dn31.codon.pair$tatccg/cps.tat.ccg)
cps.obs.exp.tatcct<-log(dn31.codon.pair$tatcct/cps.tat.cct)
cps.obs.exp.tatcga<-log(dn31.codon.pair$tatcga/cps.tat.cga)
cps.obs.exp.tatcgc<-log(dn31.codon.pair$tatcgc/cps.tat.cgc)
cps.obs.exp.tatcgg<-log(dn31.codon.pair$tatcgg/cps.tat.cgg)
cps.obs.exp.tatcgt<-log(dn31.codon.pair$tatcgt/cps.tat.cgt)
cps.obs.exp.tatcta<-log(dn31.codon.pair$tatcta/cps.tat.cta)
cps.obs.exp.tatctc<-log(dn31.codon.pair$tatctc/cps.tat.ctc)
cps.obs.exp.tatctg<-log(dn31.codon.pair$tatctg/cps.tat.ctg)
cps.obs.exp.tatctt<-log(dn31.codon.pair$tatctt/cps.tat.ctt)

cps.obs.exp.tatgaa<-log(dn31.codon.pair$tatgaa/cps.tat.gaa)
cps.obs.exp.tatgac<-log(dn31.codon.pair$tatgac/cps.tat.gac)
cps.obs.exp.tatgag<-log(dn31.codon.pair$tatgag/cps.tat.gag)
cps.obs.exp.tatgat<-log(dn31.codon.pair$tatgat/cps.tat.gat)
cps.obs.exp.tatgca<-log(dn31.codon.pair$tatgca/cps.tat.gca)
cps.obs.exp.tatgcc<-log(dn31.codon.pair$tatgcc/cps.tat.gcc)
cps.obs.exp.tatgcg<-log(dn31.codon.pair$tatgcg/cps.tat.gcg)
cps.obs.exp.tatgct<-log(dn31.codon.pair$tatgct/cps.tat.gct)
cps.obs.exp.tatgga<-log(dn31.codon.pair$tatgga/cps.tat.gga)
cps.obs.exp.tatggc<-log(dn31.codon.pair$tatggc/cps.tat.ggc)
cps.obs.exp.tatggg<-log(dn31.codon.pair$tatggg/cps.tat.ggg)
cps.obs.exp.tatggt<-log(dn31.codon.pair$tatggt/cps.tat.ggt)
cps.obs.exp.tatgta<-log(dn31.codon.pair$tatgta/cps.tat.gta)
cps.obs.exp.tatgtc<-log(dn31.codon.pair$tatgtc/cps.tat.gtc)
cps.obs.exp.tatgtg<-log(dn31.codon.pair$tatgtg/cps.tat.gtg)
cps.obs.exp.tatgtt<-log(dn31.codon.pair$tatgtt/cps.tat.gtt)

#cps.obs.exp.tattaa<-log(dn31.codon.pair$tattaa/cps.tat.taa)
cps.obs.exp.tattac<-log(dn31.codon.pair$tattac/cps.tat.tac)
#cps.obs.exp.tattag<-log(dn31.codon.pair$tattag/cps.tat.tag)
cps.obs.exp.tattat<-log(dn31.codon.pair$tattat/cps.tat.tat)
cps.obs.exp.tattca<-log(dn31.codon.pair$tattca/cps.tat.tca)
cps.obs.exp.tattcc<-log(dn31.codon.pair$tattcc/cps.tat.tcc)
cps.obs.exp.tattcg<-log(dn31.codon.pair$tattcg/cps.tat.tcg)
cps.obs.exp.tattct<-log(dn31.codon.pair$tattct/cps.tat.tct)
#cps.obs.exp.tattga<-log(dn31.codon.pair$tattga/cps.tat.tga)
cps.obs.exp.tattgc<-log(dn31.codon.pair$tattgc/cps.tat.tgc)
cps.obs.exp.tattgg<-log(dn31.codon.pair$tattgg/cps.tat.tgg)
cps.obs.exp.tattgt<-log(dn31.codon.pair$tattgt/cps.tat.tgt)
cps.obs.exp.tattta<-log(dn31.codon.pair$tattta/cps.tat.tta)
cps.obs.exp.tatttc<-log(dn31.codon.pair$tatttc/cps.tat.ttc)
cps.obs.exp.tatttg<-log(dn31.codon.pair$tatttg/cps.tat.ttg)
cps.obs.exp.tatttt<-log(dn31.codon.pair$tatttt/cps.tat.ttt)

















cps.obs.exp.tcaaaa<-log(dn31.codon.pair$tcaaaa/cps.tca.aaa)
cps.obs.exp.tcaaac<-log(dn31.codon.pair$tcaaac/cps.tca.aac)
cps.obs.exp.tcaaag<-log(dn31.codon.pair$tcaaag/cps.tca.aag)
cps.obs.exp.tcaaat<-log(dn31.codon.pair$tcaaat/cps.tca.aat)
cps.obs.exp.tcaaca<-log(dn31.codon.pair$tcaaca/cps.tca.aca)
cps.obs.exp.tcaacc<-log(dn31.codon.pair$tcaacc/cps.tca.acc)
cps.obs.exp.tcaacg<-log(dn31.codon.pair$tcaacg/cps.tca.acg)
cps.obs.exp.tcaact<-log(dn31.codon.pair$tcaact/cps.tca.act)
cps.obs.exp.tcaaga<-log(dn31.codon.pair$tcaaga/cps.tca.aga)
cps.obs.exp.tcaagc<-log(dn31.codon.pair$tcaagc/cps.tca.agc)
cps.obs.exp.tcaagg<-log(dn31.codon.pair$tcaagg/cps.tca.agg)
cps.obs.exp.tcaagt<-log(dn31.codon.pair$tcaagt/cps.tca.agt)
cps.obs.exp.tcaata<-log(dn31.codon.pair$tcaata/cps.tca.ata)
cps.obs.exp.tcaatc<-log(dn31.codon.pair$tcaatc/cps.tca.atc)
cps.obs.exp.tcaatg<-log(dn31.codon.pair$tcaatg/cps.tca.atg)
cps.obs.exp.tcaatt<-log(dn31.codon.pair$tcaatt/cps.tca.att)

cps.obs.exp.tcacaa<-log(dn31.codon.pair$tcacaa/cps.tca.caa)
cps.obs.exp.tcacac<-log(dn31.codon.pair$tcacac/cps.tca.cac)
cps.obs.exp.tcacag<-log(dn31.codon.pair$tcacag/cps.tca.cag)
cps.obs.exp.tcacat<-log(dn31.codon.pair$tcacat/cps.tca.cat)
cps.obs.exp.tcacca<-log(dn31.codon.pair$tcacca/cps.tca.cca)
cps.obs.exp.tcaccc<-log(dn31.codon.pair$tcaccc/cps.tca.ccc)
cps.obs.exp.tcaccg<-log(dn31.codon.pair$tcaccg/cps.tca.ccg)
cps.obs.exp.tcacct<-log(dn31.codon.pair$tcacct/cps.tca.cct)
cps.obs.exp.tcacga<-log(dn31.codon.pair$tcacga/cps.tca.cga)
cps.obs.exp.tcacgc<-log(dn31.codon.pair$tcacgc/cps.tca.cgc)
cps.obs.exp.tcacgg<-log(dn31.codon.pair$tcacgg/cps.tca.cgg)
cps.obs.exp.tcacgt<-log(dn31.codon.pair$tcacgt/cps.tca.cgt)
cps.obs.exp.tcacta<-log(dn31.codon.pair$tcacta/cps.tca.cta)
cps.obs.exp.tcactc<-log(dn31.codon.pair$tcactc/cps.tca.ctc)
cps.obs.exp.tcactg<-log(dn31.codon.pair$tcactg/cps.tca.ctg)
cps.obs.exp.tcactt<-log(dn31.codon.pair$tcactt/cps.tca.ctt)

cps.obs.exp.tcagaa<-log(dn31.codon.pair$tcagaa/cps.tca.gaa)
cps.obs.exp.tcagac<-log(dn31.codon.pair$tcagac/cps.tca.gac)
cps.obs.exp.tcagag<-log(dn31.codon.pair$tcagag/cps.tca.gag)
cps.obs.exp.tcagat<-log(dn31.codon.pair$tcagat/cps.tca.gat)
cps.obs.exp.tcagca<-log(dn31.codon.pair$tcagca/cps.tca.gca)
cps.obs.exp.tcagcc<-log(dn31.codon.pair$tcagcc/cps.tca.gcc)
cps.obs.exp.tcagcg<-log(dn31.codon.pair$tcagcg/cps.tca.gcg)
cps.obs.exp.tcagct<-log(dn31.codon.pair$tcagct/cps.tca.gct)
cps.obs.exp.tcagga<-log(dn31.codon.pair$tcagga/cps.tca.gga)
cps.obs.exp.tcaggc<-log(dn31.codon.pair$tcaggc/cps.tca.ggc)
cps.obs.exp.tcaggg<-log(dn31.codon.pair$tcaggg/cps.tca.ggg)
cps.obs.exp.tcaggt<-log(dn31.codon.pair$tcaggt/cps.tca.ggt)
cps.obs.exp.tcagta<-log(dn31.codon.pair$tcagta/cps.tca.gta)
cps.obs.exp.tcagtc<-log(dn31.codon.pair$tcagtc/cps.tca.gtc)
cps.obs.exp.tcagtg<-log(dn31.codon.pair$tcagtg/cps.tca.gtg)
cps.obs.exp.tcagtt<-log(dn31.codon.pair$tcagtt/cps.tca.gtt)

#cps.obs.exp.tcataa<-log(dn31.codon.pair$tcataa/cps.tca.taa)
cps.obs.exp.tcatac<-log(dn31.codon.pair$tcatac/cps.tca.tac)
#cps.obs.exp.tcatag<-log(dn31.codon.pair$tcatag/cps.tca.tag)
cps.obs.exp.tcatat<-log(dn31.codon.pair$tcatat/cps.tca.tat)
cps.obs.exp.tcatca<-log(dn31.codon.pair$tcatca/cps.tca.tca)
cps.obs.exp.tcatcc<-log(dn31.codon.pair$tcatcc/cps.tca.tcc)
cps.obs.exp.tcatcg<-log(dn31.codon.pair$tcatcg/cps.tca.tcg)
cps.obs.exp.tcatct<-log(dn31.codon.pair$tcatct/cps.tca.tct)
#cps.obs.exp.tcatga<-log(dn31.codon.pair$tcatga/cps.tca.tga)
cps.obs.exp.tcatgc<-log(dn31.codon.pair$tcatgc/cps.tca.tgc)
cps.obs.exp.tcatgg<-log(dn31.codon.pair$tcatgg/cps.tca.tgg)
cps.obs.exp.tcatgt<-log(dn31.codon.pair$tcatgt/cps.tca.tgt)
cps.obs.exp.tcatta<-log(dn31.codon.pair$tcatta/cps.tca.tta)
cps.obs.exp.tcattc<-log(dn31.codon.pair$tcattc/cps.tca.ttc)
cps.obs.exp.tcattg<-log(dn31.codon.pair$tcattg/cps.tca.ttg)
cps.obs.exp.tcattt<-log(dn31.codon.pair$tcattt/cps.tca.ttt)








cps.obs.exp.tccaaa<-log(dn31.codon.pair$tccaaa/cps.tcc.aaa)
cps.obs.exp.tccaac<-log(dn31.codon.pair$tccaac/cps.tcc.aac)
cps.obs.exp.tccaag<-log(dn31.codon.pair$tccaag/cps.tcc.aag)
cps.obs.exp.tccaat<-log(dn31.codon.pair$tccaat/cps.tcc.aat)
cps.obs.exp.tccaca<-log(dn31.codon.pair$tccaca/cps.tcc.aca)
cps.obs.exp.tccacc<-log(dn31.codon.pair$tccacc/cps.tcc.acc)
cps.obs.exp.tccacg<-log(dn31.codon.pair$tccacg/cps.tcc.acg)
cps.obs.exp.tccact<-log(dn31.codon.pair$tccact/cps.tcc.act)
cps.obs.exp.tccaga<-log(dn31.codon.pair$tccaga/cps.tcc.aga)
cps.obs.exp.tccagc<-log(dn31.codon.pair$tccagc/cps.tcc.agc)
cps.obs.exp.tccagg<-log(dn31.codon.pair$tccagg/cps.tcc.agg)
cps.obs.exp.tccagt<-log(dn31.codon.pair$tccagt/cps.tcc.agt)
cps.obs.exp.tccata<-log(dn31.codon.pair$tccata/cps.tcc.ata)
cps.obs.exp.tccatc<-log(dn31.codon.pair$tccatc/cps.tcc.atc)
cps.obs.exp.tccatg<-log(dn31.codon.pair$tccatg/cps.tcc.atg)
cps.obs.exp.tccatt<-log(dn31.codon.pair$tccatt/cps.tcc.att)

cps.obs.exp.tcccaa<-log(dn31.codon.pair$tcccaa/cps.tcc.caa)
cps.obs.exp.tcccac<-log(dn31.codon.pair$tcccac/cps.tcc.cac)
cps.obs.exp.tcccag<-log(dn31.codon.pair$tcccag/cps.tcc.cag)
cps.obs.exp.tcccat<-log(dn31.codon.pair$tcccat/cps.tcc.cat)
cps.obs.exp.tcccca<-log(dn31.codon.pair$tcccca/cps.tcc.cca)
cps.obs.exp.tccccc<-log(dn31.codon.pair$tccccc/cps.tcc.ccc)
cps.obs.exp.tccccg<-log(dn31.codon.pair$tccccg/cps.tcc.ccg)
cps.obs.exp.tcccct<-log(dn31.codon.pair$tcccct/cps.tcc.cct)
cps.obs.exp.tcccga<-log(dn31.codon.pair$tcccga/cps.tcc.cga)
cps.obs.exp.tcccgc<-log(dn31.codon.pair$tcccgc/cps.tcc.cgc)
cps.obs.exp.tcccgg<-log(dn31.codon.pair$tcccgg/cps.tcc.cgg)
cps.obs.exp.tcccgt<-log(dn31.codon.pair$tcccgt/cps.tcc.cgt)
cps.obs.exp.tcccta<-log(dn31.codon.pair$tcccta/cps.tcc.cta)
cps.obs.exp.tccctc<-log(dn31.codon.pair$tccctc/cps.tcc.ctc)
cps.obs.exp.tccctg<-log(dn31.codon.pair$tccctg/cps.tcc.ctg)
cps.obs.exp.tccctt<-log(dn31.codon.pair$tccctt/cps.tcc.ctt)

cps.obs.exp.tccgaa<-log(dn31.codon.pair$tccgaa/cps.tcc.gaa)
cps.obs.exp.tccgac<-log(dn31.codon.pair$tccgac/cps.tcc.gac)
cps.obs.exp.tccgag<-log(dn31.codon.pair$tccgag/cps.tcc.gag)
cps.obs.exp.tccgat<-log(dn31.codon.pair$tccgat/cps.tcc.gat)
cps.obs.exp.tccgca<-log(dn31.codon.pair$tccgca/cps.tcc.gca)
cps.obs.exp.tccgcc<-log(dn31.codon.pair$tccgcc/cps.tcc.gcc)
cps.obs.exp.tccgcg<-log(dn31.codon.pair$tccgcg/cps.tcc.gcg)
cps.obs.exp.tccgct<-log(dn31.codon.pair$tccgct/cps.tcc.gct)
cps.obs.exp.tccgga<-log(dn31.codon.pair$tccgga/cps.tcc.gga)
cps.obs.exp.tccggc<-log(dn31.codon.pair$tccggc/cps.tcc.ggc)
cps.obs.exp.tccggg<-log(dn31.codon.pair$tccggg/cps.tcc.ggg)
cps.obs.exp.tccggt<-log(dn31.codon.pair$tccggt/cps.tcc.ggt)
cps.obs.exp.tccgta<-log(dn31.codon.pair$tccgta/cps.tcc.gta)
cps.obs.exp.tccgtc<-log(dn31.codon.pair$tccgtc/cps.tcc.gtc)
cps.obs.exp.tccgtg<-log(dn31.codon.pair$tccgtg/cps.tcc.gtg)
cps.obs.exp.tccgtt<-log(dn31.codon.pair$tccgtt/cps.tcc.gtt)

#cps.obs.exp.tcctaa<-log(dn31.codon.pair$tcctaa/cps.tcc.taa)
cps.obs.exp.tcctac<-log(dn31.codon.pair$tcctac/cps.tcc.tac)
#cps.obs.exp.tcctag<-log(dn31.codon.pair$tcctag/cps.tcc.tag)
cps.obs.exp.tcctat<-log(dn31.codon.pair$tcctat/cps.tcc.tat)
cps.obs.exp.tcctca<-log(dn31.codon.pair$tcctca/cps.tcc.tca)
cps.obs.exp.tcctcc<-log(dn31.codon.pair$tcctcc/cps.tcc.tcc)
cps.obs.exp.tcctcg<-log(dn31.codon.pair$tcctcg/cps.tcc.tcg)
cps.obs.exp.tcctct<-log(dn31.codon.pair$tcctct/cps.tcc.tct)
#cps.obs.exp.tcctga<-log(dn31.codon.pair$tcctga/cps.tcc.tga)
cps.obs.exp.tcctgc<-log(dn31.codon.pair$tcctgc/cps.tcc.tgc)
cps.obs.exp.tcctgg<-log(dn31.codon.pair$tcctgg/cps.tcc.tgg)
cps.obs.exp.tcctgt<-log(dn31.codon.pair$tcctgt/cps.tcc.tgt)
cps.obs.exp.tcctta<-log(dn31.codon.pair$tcctta/cps.tcc.tta)
cps.obs.exp.tccttc<-log(dn31.codon.pair$tccttc/cps.tcc.ttc)
cps.obs.exp.tccttg<-log(dn31.codon.pair$tccttg/cps.tcc.ttg)
cps.obs.exp.tccttt<-log(dn31.codon.pair$tccttt/cps.tcc.ttt)









cps.obs.exp.tcgaaa<-log(dn31.codon.pair$tcgaaa/cps.tcg.aaa)
cps.obs.exp.tcgaac<-log(dn31.codon.pair$tcgaac/cps.tcg.aac)
cps.obs.exp.tcgaag<-log(dn31.codon.pair$tcgaag/cps.tcg.aag)
cps.obs.exp.tcgaat<-log(dn31.codon.pair$tcgaat/cps.tcg.aat)
cps.obs.exp.tcgaca<-log(dn31.codon.pair$tcgaca/cps.tcg.aca)
cps.obs.exp.tcgacc<-log(dn31.codon.pair$tcgacc/cps.tcg.acc)
cps.obs.exp.tcgacg<-log(dn31.codon.pair$tcgacg/cps.tcg.acg)
cps.obs.exp.tcgact<-log(dn31.codon.pair$tcgact/cps.tcg.act)
cps.obs.exp.tcgaga<-log(dn31.codon.pair$tcgaga/cps.tcg.aga)
cps.obs.exp.tcgagc<-log(dn31.codon.pair$tcgagc/cps.tcg.agc)
cps.obs.exp.tcgagg<-log(dn31.codon.pair$tcgagg/cps.tcg.agg)
cps.obs.exp.tcgagt<-log(dn31.codon.pair$tcgagt/cps.tcg.agt)
cps.obs.exp.tcgata<-log(dn31.codon.pair$tcgata/cps.tcg.ata)
cps.obs.exp.tcgatc<-log(dn31.codon.pair$tcgatc/cps.tcg.atc)
cps.obs.exp.tcgatg<-log(dn31.codon.pair$tcgatg/cps.tcg.atg)
cps.obs.exp.tcgatt<-log(dn31.codon.pair$tcgatt/cps.tcg.att)

cps.obs.exp.tcgcaa<-log(dn31.codon.pair$tcgcaa/cps.tcg.caa)
cps.obs.exp.tcgcac<-log(dn31.codon.pair$tcgcac/cps.tcg.cac)
cps.obs.exp.tcgcag<-log(dn31.codon.pair$tcgcag/cps.tcg.cag)
cps.obs.exp.tcgcat<-log(dn31.codon.pair$tcgcat/cps.tcg.cat)
cps.obs.exp.tcgcca<-log(dn31.codon.pair$tcgcca/cps.tcg.cca)
cps.obs.exp.tcgccc<-log(dn31.codon.pair$tcgccc/cps.tcg.ccc)
cps.obs.exp.tcgccg<-log(dn31.codon.pair$tcgccg/cps.tcg.ccg)
cps.obs.exp.tcgcct<-log(dn31.codon.pair$tcgcct/cps.tcg.cct)
cps.obs.exp.tcgcga<-log(dn31.codon.pair$tcgcga/cps.tcg.cga)
cps.obs.exp.tcgcgc<-log(dn31.codon.pair$tcgcgc/cps.tcg.cgc)
cps.obs.exp.tcgcgg<-log(dn31.codon.pair$tcgcgg/cps.tcg.cgg)
cps.obs.exp.tcgcgt<-log(dn31.codon.pair$tcgcgt/cps.tcg.cgt)
cps.obs.exp.tcgcta<-log(dn31.codon.pair$tcgcta/cps.tcg.cta)
cps.obs.exp.tcgctc<-log(dn31.codon.pair$tcgctc/cps.tcg.ctc)
cps.obs.exp.tcgctg<-log(dn31.codon.pair$tcgctg/cps.tcg.ctg)
cps.obs.exp.tcgctt<-log(dn31.codon.pair$tcgctt/cps.tcg.ctt)

cps.obs.exp.tcggaa<-log(dn31.codon.pair$tcggaa/cps.tcg.gaa)
cps.obs.exp.tcggac<-log(dn31.codon.pair$tcggac/cps.tcg.gac)
cps.obs.exp.tcggag<-log(dn31.codon.pair$tcggag/cps.tcg.gag)
cps.obs.exp.tcggat<-log(dn31.codon.pair$tcggat/cps.tcg.gat)
cps.obs.exp.tcggca<-log(dn31.codon.pair$tcggca/cps.tcg.gca)
cps.obs.exp.tcggcc<-log(dn31.codon.pair$tcggcc/cps.tcg.gcc)
cps.obs.exp.tcggcg<-log(dn31.codon.pair$tcggcg/cps.tcg.gcg)
cps.obs.exp.tcggct<-log(dn31.codon.pair$tcggct/cps.tcg.gct)
cps.obs.exp.tcggga<-log(dn31.codon.pair$tcggga/cps.tcg.gga)
cps.obs.exp.tcgggc<-log(dn31.codon.pair$tcgggc/cps.tcg.ggc)
cps.obs.exp.tcgggg<-log(dn31.codon.pair$tcgggg/cps.tcg.ggg)
cps.obs.exp.tcgggt<-log(dn31.codon.pair$tcgggt/cps.tcg.ggt)
cps.obs.exp.tcggta<-log(dn31.codon.pair$tcggta/cps.tcg.gta)
cps.obs.exp.tcggtc<-log(dn31.codon.pair$tcggtc/cps.tcg.gtc)
cps.obs.exp.tcggtg<-log(dn31.codon.pair$tcggtg/cps.tcg.gtg)
cps.obs.exp.tcggtt<-log(dn31.codon.pair$tcggtt/cps.tcg.gtt)

#cps.obs.exp.tcgtaa<-log(dn31.codon.pair$tcgtaa/cps.tcg.taa)
cps.obs.exp.tcgtac<-log(dn31.codon.pair$tcgtac/cps.tcg.tac)
#cps.obs.exp.tcgtag<-log(dn31.codon.pair$tcgtag/cps.tcg.tag)
cps.obs.exp.tcgtat<-log(dn31.codon.pair$tcgtat/cps.tcg.tat)
cps.obs.exp.tcgtca<-log(dn31.codon.pair$tcgtca/cps.tcg.tca)
cps.obs.exp.tcgtcc<-log(dn31.codon.pair$tcgtcc/cps.tcg.tcc)
cps.obs.exp.tcgtcg<-log(dn31.codon.pair$tcgtcg/cps.tcg.tcg)
cps.obs.exp.tcgtct<-log(dn31.codon.pair$tcgtct/cps.tcg.tct)
#cps.obs.exp.tcgtga<-log(dn31.codon.pair$tcgtga/cps.tcg.tga)
cps.obs.exp.tcgtgc<-log(dn31.codon.pair$tcgtgc/cps.tcg.tgc)
cps.obs.exp.tcgtgg<-log(dn31.codon.pair$tcgtgg/cps.tcg.tgg)
cps.obs.exp.tcgtgt<-log(dn31.codon.pair$tcgtgt/cps.tcg.tgt)
cps.obs.exp.tcgtta<-log(dn31.codon.pair$tcgtta/cps.tcg.tta)
cps.obs.exp.tcgttc<-log(dn31.codon.pair$tcgttc/cps.tcg.ttc)
cps.obs.exp.tcgttg<-log(dn31.codon.pair$tcgttg/cps.tcg.ttg)
cps.obs.exp.tcgttt<-log(dn31.codon.pair$tcgttt/cps.tcg.ttt)









cps.obs.exp.tctaaa<-log(dn31.codon.pair$tctaaa/cps.tct.aaa)
cps.obs.exp.tctaac<-log(dn31.codon.pair$tctaac/cps.tct.aac)
cps.obs.exp.tctaag<-log(dn31.codon.pair$tctaag/cps.tct.aag)
cps.obs.exp.tctaat<-log(dn31.codon.pair$tctaat/cps.tct.aat)
cps.obs.exp.tctaca<-log(dn31.codon.pair$tctaca/cps.tct.aca)
cps.obs.exp.tctacc<-log(dn31.codon.pair$tctacc/cps.tct.acc)
cps.obs.exp.tctacg<-log(dn31.codon.pair$tctacg/cps.tct.acg)
cps.obs.exp.tctact<-log(dn31.codon.pair$tctact/cps.tct.act)
cps.obs.exp.tctaga<-log(dn31.codon.pair$tctaga/cps.tct.aga)
cps.obs.exp.tctagc<-log(dn31.codon.pair$tctagc/cps.tct.agc)
cps.obs.exp.tctagg<-log(dn31.codon.pair$tctagg/cps.tct.agg)
cps.obs.exp.tctagt<-log(dn31.codon.pair$tctagt/cps.tct.agt)
cps.obs.exp.tctata<-log(dn31.codon.pair$tctata/cps.tct.ata)
cps.obs.exp.tctatc<-log(dn31.codon.pair$tctatc/cps.tct.atc)
cps.obs.exp.tctatg<-log(dn31.codon.pair$tctatg/cps.tct.atg)
cps.obs.exp.tctatt<-log(dn31.codon.pair$tctatt/cps.tct.att)

cps.obs.exp.tctcaa<-log(dn31.codon.pair$tctcaa/cps.tct.caa)
cps.obs.exp.tctcac<-log(dn31.codon.pair$tctcac/cps.tct.cac)
cps.obs.exp.tctcag<-log(dn31.codon.pair$tctcag/cps.tct.cag)
cps.obs.exp.tctcat<-log(dn31.codon.pair$tctcat/cps.tct.cat)
cps.obs.exp.tctcca<-log(dn31.codon.pair$tctcca/cps.tct.cca)
cps.obs.exp.tctccc<-log(dn31.codon.pair$tctccc/cps.tct.ccc)
cps.obs.exp.tctccg<-log(dn31.codon.pair$tctccg/cps.tct.ccg)
cps.obs.exp.tctcct<-log(dn31.codon.pair$tctcct/cps.tct.cct)
cps.obs.exp.tctcga<-log(dn31.codon.pair$tctcga/cps.tct.cga)
cps.obs.exp.tctcgc<-log(dn31.codon.pair$tctcgc/cps.tct.cgc)
cps.obs.exp.tctcgg<-log(dn31.codon.pair$tctcgg/cps.tct.cgg)
cps.obs.exp.tctcgt<-log(dn31.codon.pair$tctcgt/cps.tct.cgt)
cps.obs.exp.tctcta<-log(dn31.codon.pair$tctcta/cps.tct.cta)
cps.obs.exp.tctctc<-log(dn31.codon.pair$tctctc/cps.tct.ctc)
cps.obs.exp.tctctg<-log(dn31.codon.pair$tctctg/cps.tct.ctg)
cps.obs.exp.tctctt<-log(dn31.codon.pair$tctctt/cps.tct.ctt)

cps.obs.exp.tctgaa<-log(dn31.codon.pair$tctgaa/cps.tct.gaa)
cps.obs.exp.tctgac<-log(dn31.codon.pair$tctgac/cps.tct.gac)
cps.obs.exp.tctgag<-log(dn31.codon.pair$tctgag/cps.tct.gag)
cps.obs.exp.tctgat<-log(dn31.codon.pair$tctgat/cps.tct.gat)
cps.obs.exp.tctgca<-log(dn31.codon.pair$tctgca/cps.tct.gca)
cps.obs.exp.tctgcc<-log(dn31.codon.pair$tctgcc/cps.tct.gcc)
cps.obs.exp.tctgcg<-log(dn31.codon.pair$tctgcg/cps.tct.gcg)
cps.obs.exp.tctgct<-log(dn31.codon.pair$tctgct/cps.tct.gct)
cps.obs.exp.tctgga<-log(dn31.codon.pair$tctgga/cps.tct.gga)
cps.obs.exp.tctggc<-log(dn31.codon.pair$tctggc/cps.tct.ggc)
cps.obs.exp.tctggg<-log(dn31.codon.pair$tctggg/cps.tct.ggg)
cps.obs.exp.tctggt<-log(dn31.codon.pair$tctggt/cps.tct.ggt)
cps.obs.exp.tctgta<-log(dn31.codon.pair$tctgta/cps.tct.gta)
cps.obs.exp.tctgtc<-log(dn31.codon.pair$tctgtc/cps.tct.gtc)
cps.obs.exp.tctgtg<-log(dn31.codon.pair$tctgtg/cps.tct.gtg)
cps.obs.exp.tctgtt<-log(dn31.codon.pair$tctgtt/cps.tct.gtt)

#cps.obs.exp.tcttaa<-log(dn31.codon.pair$tcttaa/cps.tct.taa)
cps.obs.exp.tcttac<-log(dn31.codon.pair$tcttac/cps.tct.tac)
#cps.obs.exp.tcttag<-log(dn31.codon.pair$tcttag/cps.tct.tag)
cps.obs.exp.tcttat<-log(dn31.codon.pair$tcttat/cps.tct.tat)
cps.obs.exp.tcttca<-log(dn31.codon.pair$tcttca/cps.tct.tca)
cps.obs.exp.tcttcc<-log(dn31.codon.pair$tcttcc/cps.tct.tcc)
cps.obs.exp.tcttcg<-log(dn31.codon.pair$tcttcg/cps.tct.tcg)
cps.obs.exp.tcttct<-log(dn31.codon.pair$tcttct/cps.tct.tct)
#cps.obs.exp.tcttga<-log(dn31.codon.pair$tcttga/cps.tct.tga)
cps.obs.exp.tcttgc<-log(dn31.codon.pair$tcttgc/cps.tct.tgc)
cps.obs.exp.tcttgg<-log(dn31.codon.pair$tcttgg/cps.tct.tgg)
cps.obs.exp.tcttgt<-log(dn31.codon.pair$tcttgt/cps.tct.tgt)
cps.obs.exp.tcttta<-log(dn31.codon.pair$tcttta/cps.tct.tta)
cps.obs.exp.tctttc<-log(dn31.codon.pair$tctttc/cps.tct.ttc)
cps.obs.exp.tctttg<-log(dn31.codon.pair$tctttg/cps.tct.ttg)
cps.obs.exp.tctttt<-log(dn31.codon.pair$tctttt/cps.tct.ttt)
















#Stop Codon 


#cps.obs.exp.tgaaaa<-log(dn31.codon.pair$tgaaaa/cps.tga.aaa)
#cps.obs.exp.tgaaac<-log(dn31.codon.pair$tgaaac/cps.tga.aac)
#cps.obs.exp.tgaaag<-log(dn31.codon.pair$tgaaag/cps.tga.aag)
#cps.obs.exp.tgaaat<-log(dn31.codon.pair$tgaaat/cps.tga.aat)
#cps.obs.exp.tgaaca<-log(dn31.codon.pair$tgaaca/cps.tga.aca)
#cps.obs.exp.tgaacc<-log(dn31.codon.pair$tgaacc/cps.tga.acc)
#cps.obs.exp.tgaacg<-log(dn31.codon.pair$tgaacg/cps.tga.acg)
#cps.obs.exp.tgaact<-log(dn31.codon.pair$tgaact/cps.tga.act)
#cps.obs.exp.tgaaga<-log(dn31.codon.pair$tgaaga/cps.tga.aga)
#cps.obs.exp.tgaagc<-log(dn31.codon.pair$tgaagc/cps.tga.agc)
#cps.obs.exp.tgaagg<-log(dn31.codon.pair$tgaagg/cps.tga.agg)
#cps.obs.exp.tgaagt<-log(dn31.codon.pair$tgaagt/cps.tga.agt)
#cps.obs.exp.tgaata<-log(dn31.codon.pair$tgaata/cps.tga.ata)
#cps.obs.exp.tgaatc<-log(dn31.codon.pair$tgaatc/cps.tga.atc)
#cps.obs.exp.tgaatg<-log(dn31.codon.pair$tgaatg/cps.tga.atg)
#cps.obs.exp.tgaatt<-log(dn31.codon.pair$tgaatt/cps.tga.att)

#cps.obs.exp.tgacaa<-log(dn31.codon.pair$tgacaa/cps.tga.caa)
#cps.obs.exp.tgacac<-log(dn31.codon.pair$tgacac/cps.tga.cac)
#cps.obs.exp.tgacag<-log(dn31.codon.pair$tgacag/cps.tga.cag)
#cps.obs.exp.tgacat<-log(dn31.codon.pair$tgacat/cps.tga.cat)
#cps.obs.exp.tgacca<-log(dn31.codon.pair$tgacca/cps.tga.cca)
#cps.obs.exp.tgaccc<-log(dn31.codon.pair$tgaccc/cps.tga.ccc)
#cps.obs.exp.tgaccg<-log(dn31.codon.pair$tgaccg/cps.tga.ccg)
#cps.obs.exp.tgacct<-log(dn31.codon.pair$tgacct/cps.tga.cct)
#cps.obs.exp.tgacga<-log(dn31.codon.pair$tgacga/cps.tga.cga)
#cps.obs.exp.tgacgc<-log(dn31.codon.pair$tgacgc/cps.tga.cgc)
#cps.obs.exp.tgacgg<-log(dn31.codon.pair$tgacgg/cps.tga.cgg)
#cps.obs.exp.tgacgt<-log(dn31.codon.pair$tgacgt/cps.tga.cgt)
#cps.obs.exp.tgacta<-log(dn31.codon.pair$tgacta/cps.tga.cta)
#cps.obs.exp.tgactc<-log(dn31.codon.pair$tgactc/cps.tga.ctc)
#cps.obs.exp.tgactg<-log(dn31.codon.pair$tgactg/cps.tga.ctg)
#cps.obs.exp.tgactt<-log(dn31.codon.pair$tgactt/cps.tga.ctt)

#cps.obs.exp.tgagaa<-log(dn31.codon.pair$tgagaa/cps.tga.gaa)
#cps.obs.exp.tgagac<-log(dn31.codon.pair$tgagac/cps.tga.gac)
#cps.obs.exp.tgagag<-log(dn31.codon.pair$tgagag/cps.tga.gag)
#cps.obs.exp.tgagat<-log(dn31.codon.pair$tgagat/cps.tga.gat)
#cps.obs.exp.tgagca<-log(dn31.codon.pair$tgagca/cps.tga.gca)
#cps.obs.exp.tgagcc<-log(dn31.codon.pair$tgagcc/cps.tga.gcc)
#cps.obs.exp.tgagcg<-log(dn31.codon.pair$tgagcg/cps.tga.gcg)
#cps.obs.exp.tgagct<-log(dn31.codon.pair$tgagct/cps.tga.gct)
#cps.obs.exp.tgagga<-log(dn31.codon.pair$tgagga/cps.tga.gga)
#cps.obs.exp.tgaggc<-log(dn31.codon.pair$tgaggc/cps.tga.ggc)
#cps.obs.exp.tgaggg<-log(dn31.codon.pair$tgaggg/cps.tga.ggg)
#cps.obs.exp.tgaggt<-log(dn31.codon.pair$tgaggt/cps.tga.ggt)
#cps.obs.exp.tgagta<-log(dn31.codon.pair$tgagta/cps.tga.gta)
#cps.obs.exp.tgagtc<-log(dn31.codon.pair$tgagtc/cps.tga.gtc)
#cps.obs.exp.tgagtg<-log(dn31.codon.pair$tgagtg/cps.tga.gtg)
#cps.obs.exp.tgagtt<-log(dn31.codon.pair$tgagtt/cps.tga.gtt)

#cps.obs.exp.tgataa<-log(dn31.codon.pair$tgataa/cps.tga.taa)
#cps.obs.exp.tgatac<-log(dn31.codon.pair$tgatac/cps.tga.tac)
#cps.obs.exp.tgatag<-log(dn31.codon.pair$tgatag/cps.tga.tag)
#cps.obs.exp.tgatat<-log(dn31.codon.pair$tgatat/cps.tga.tat)
#cps.obs.exp.tgatca<-log(dn31.codon.pair$tgatca/cps.tga.tca)
#cps.obs.exp.tgatcc<-log(dn31.codon.pair$tgatcc/cps.tga.tcc)
#cps.obs.exp.tgatcg<-log(dn31.codon.pair$tgatcg/cps.tga.tcg)
#cps.obs.exp.tgatct<-log(dn31.codon.pair$tgatct/cps.tga.tct)
#cps.obs.exp.tgatga<-log(dn31.codon.pair$tgatga/cps.tga.tga)
#cps.obs.exp.tgatgc<-log(dn31.codon.pair$tgatgc/cps.tga.tgc)
#cps.obs.exp.tgatgg<-log(dn31.codon.pair$tgatgg/cps.tga.tgg)
#cps.obs.exp.tgatgt<-log(dn31.codon.pair$tgatgt/cps.tga.tgt)
#cps.obs.exp.tgatta<-log(dn31.codon.pair$tgatta/cps.tga.tta)
#cps.obs.exp.tgattc<-log(dn31.codon.pair$tgattc/cps.tga.ttc)
#cps.obs.exp.tgattg<-log(dn31.codon.pair$tgattg/cps.tga.ttg)
#cps.obs.exp.tgattt<-log(dn31.codon.pair$tgattt/cps.tga.ttt)









cps.obs.exp.tgcaaa<-log(dn31.codon.pair$tgcaaa/cps.tgc.aaa)
cps.obs.exp.tgcaac<-log(dn31.codon.pair$tgcaac/cps.tgc.aac)
cps.obs.exp.tgcaag<-log(dn31.codon.pair$tgcaag/cps.tgc.aag)
cps.obs.exp.tgcaat<-log(dn31.codon.pair$tgcaat/cps.tgc.aat)
cps.obs.exp.tgcaca<-log(dn31.codon.pair$tgcaca/cps.tgc.aca)
cps.obs.exp.tgcacc<-log(dn31.codon.pair$tgcacc/cps.tgc.acc)
cps.obs.exp.tgcacg<-log(dn31.codon.pair$tgcacg/cps.tgc.acg)
cps.obs.exp.tgcact<-log(dn31.codon.pair$tgcact/cps.tgc.act)
cps.obs.exp.tgcaga<-log(dn31.codon.pair$tgcaga/cps.tgc.aga)
cps.obs.exp.tgcagc<-log(dn31.codon.pair$tgcagc/cps.tgc.agc)
cps.obs.exp.tgcagg<-log(dn31.codon.pair$tgcagg/cps.tgc.agg)
cps.obs.exp.tgcagt<-log(dn31.codon.pair$tgcagt/cps.tgc.agt)
cps.obs.exp.tgcata<-log(dn31.codon.pair$tgcata/cps.tgc.ata)
cps.obs.exp.tgcatc<-log(dn31.codon.pair$tgcatc/cps.tgc.atc)
cps.obs.exp.tgcatg<-log(dn31.codon.pair$tgcatg/cps.tgc.atg)
cps.obs.exp.tgcatt<-log(dn31.codon.pair$tgcatt/cps.tgc.att)

cps.obs.exp.tgccaa<-log(dn31.codon.pair$tgccaa/cps.tgc.caa)
cps.obs.exp.tgccac<-log(dn31.codon.pair$tgccac/cps.tgc.cac)
cps.obs.exp.tgccag<-log(dn31.codon.pair$tgccag/cps.tgc.cag)
cps.obs.exp.tgccat<-log(dn31.codon.pair$tgccat/cps.tgc.cat)
cps.obs.exp.tgccca<-log(dn31.codon.pair$tgccca/cps.tgc.cca)
cps.obs.exp.tgcccc<-log(dn31.codon.pair$tgcccc/cps.tgc.ccc)
cps.obs.exp.tgcccg<-log(dn31.codon.pair$tgcccg/cps.tgc.ccg)
cps.obs.exp.tgccct<-log(dn31.codon.pair$tgccct/cps.tgc.cct)
cps.obs.exp.tgccga<-log(dn31.codon.pair$tgccga/cps.tgc.cga)
cps.obs.exp.tgccgc<-log(dn31.codon.pair$tgccgc/cps.tgc.cgc)
cps.obs.exp.tgccgg<-log(dn31.codon.pair$tgccgg/cps.tgc.cgg)
cps.obs.exp.tgccgt<-log(dn31.codon.pair$tgccgt/cps.tgc.cgt)
cps.obs.exp.tgccta<-log(dn31.codon.pair$tgccta/cps.tgc.cta)
cps.obs.exp.tgcctc<-log(dn31.codon.pair$tgcctc/cps.tgc.ctc)
cps.obs.exp.tgcctg<-log(dn31.codon.pair$tgcctg/cps.tgc.ctg)
cps.obs.exp.tgcctt<-log(dn31.codon.pair$tgcctt/cps.tgc.ctt)

cps.obs.exp.tgcgaa<-log(dn31.codon.pair$tgcgaa/cps.tgc.gaa)
cps.obs.exp.tgcgac<-log(dn31.codon.pair$tgcgac/cps.tgc.gac)
cps.obs.exp.tgcgag<-log(dn31.codon.pair$tgcgag/cps.tgc.gag)
cps.obs.exp.tgcgat<-log(dn31.codon.pair$tgcgat/cps.tgc.gat)
cps.obs.exp.tgcgca<-log(dn31.codon.pair$tgcgca/cps.tgc.gca)
cps.obs.exp.tgcgcc<-log(dn31.codon.pair$tgcgcc/cps.tgc.gcc)
cps.obs.exp.tgcgcg<-log(dn31.codon.pair$tgcgcg/cps.tgc.gcg)
cps.obs.exp.tgcgct<-log(dn31.codon.pair$tgcgct/cps.tgc.gct)
cps.obs.exp.tgcgga<-log(dn31.codon.pair$tgcgga/cps.tgc.gga)
cps.obs.exp.tgcggc<-log(dn31.codon.pair$tgcggc/cps.tgc.ggc)
cps.obs.exp.tgcggg<-log(dn31.codon.pair$tgcggg/cps.tgc.ggg)
cps.obs.exp.tgcggt<-log(dn31.codon.pair$tgcggt/cps.tgc.ggt)
cps.obs.exp.tgcgta<-log(dn31.codon.pair$tgcgta/cps.tgc.gta)
cps.obs.exp.tgcgtc<-log(dn31.codon.pair$tgcgtc/cps.tgc.gtc)
cps.obs.exp.tgcgtg<-log(dn31.codon.pair$tgcgtg/cps.tgc.gtg)
cps.obs.exp.tgcgtt<-log(dn31.codon.pair$tgcgtt/cps.tgc.gtt)

#cps.obs.exp.tgctaa<-log(dn31.codon.pair$tgctaa/cps.tgc.taa)
cps.obs.exp.tgctac<-log(dn31.codon.pair$tgctac/cps.tgc.tac)
#cps.obs.exp.tgctag<-log(dn31.codon.pair$tgctag/cps.tgc.tag)
cps.obs.exp.tgctat<-log(dn31.codon.pair$tgctat/cps.tgc.tat)
cps.obs.exp.tgctca<-log(dn31.codon.pair$tgctca/cps.tgc.tca)
cps.obs.exp.tgctcc<-log(dn31.codon.pair$tgctcc/cps.tgc.tcc)
cps.obs.exp.tgctcg<-log(dn31.codon.pair$tgctcg/cps.tgc.tcg)
cps.obs.exp.tgctct<-log(dn31.codon.pair$tgctct/cps.tgc.tct)
#cps.obs.exp.tgctga<-log(dn31.codon.pair$tgctga/cps.tgc.tga)
cps.obs.exp.tgctgc<-log(dn31.codon.pair$tgctgc/cps.tgc.tgc)
cps.obs.exp.tgctgg<-log(dn31.codon.pair$tgctgg/cps.tgc.tgg)
cps.obs.exp.tgctgt<-log(dn31.codon.pair$tgctgt/cps.tgc.tgt)
cps.obs.exp.tgctta<-log(dn31.codon.pair$tgctta/cps.tgc.tta)
cps.obs.exp.tgcttc<-log(dn31.codon.pair$tgcttc/cps.tgc.ttc)
cps.obs.exp.tgcttg<-log(dn31.codon.pair$tgcttg/cps.tgc.ttg)
cps.obs.exp.tgcttt<-log(dn31.codon.pair$tgcttt/cps.tgc.ttt)










cps.obs.exp.tggaaa<-log(dn31.codon.pair$tggaaa/cps.tgg.aaa)
cps.obs.exp.tggaac<-log(dn31.codon.pair$tggaac/cps.tgg.aac)
cps.obs.exp.tggaag<-log(dn31.codon.pair$tggaag/cps.tgg.aag)
cps.obs.exp.tggaat<-log(dn31.codon.pair$tggaat/cps.tgg.aat)
cps.obs.exp.tggaca<-log(dn31.codon.pair$tggaca/cps.tgg.aca)
cps.obs.exp.tggacc<-log(dn31.codon.pair$tggacc/cps.tgg.acc)
cps.obs.exp.tggacg<-log(dn31.codon.pair$tggacg/cps.tgg.acg)
cps.obs.exp.tggact<-log(dn31.codon.pair$tggact/cps.tgg.act)
cps.obs.exp.tggaga<-log(dn31.codon.pair$tggaga/cps.tgg.aga)
cps.obs.exp.tggagc<-log(dn31.codon.pair$tggagc/cps.tgg.agc)
cps.obs.exp.tggagg<-log(dn31.codon.pair$tggagg/cps.tgg.agg)
cps.obs.exp.tggagt<-log(dn31.codon.pair$tggagt/cps.tgg.agt)
cps.obs.exp.tggata<-log(dn31.codon.pair$tggata/cps.tgg.ata)
cps.obs.exp.tggatc<-log(dn31.codon.pair$tggatc/cps.tgg.atc)
cps.obs.exp.tggatg<-log(dn31.codon.pair$tggatg/cps.tgg.atg)
cps.obs.exp.tggatt<-log(dn31.codon.pair$tggatt/cps.tgg.att)

cps.obs.exp.tggcaa<-log(dn31.codon.pair$tggcaa/cps.tgg.caa)
cps.obs.exp.tggcac<-log(dn31.codon.pair$tggcac/cps.tgg.cac)
cps.obs.exp.tggcag<-log(dn31.codon.pair$tggcag/cps.tgg.cag)
cps.obs.exp.tggcat<-log(dn31.codon.pair$tggcat/cps.tgg.cat)
cps.obs.exp.tggcca<-log(dn31.codon.pair$tggcca/cps.tgg.cca)
cps.obs.exp.tggccc<-log(dn31.codon.pair$tggccc/cps.tgg.ccc)
cps.obs.exp.tggccg<-log(dn31.codon.pair$tggccg/cps.tgg.ccg)
cps.obs.exp.tggcct<-log(dn31.codon.pair$tggcct/cps.tgg.cct)
cps.obs.exp.tggcga<-log(dn31.codon.pair$tggcga/cps.tgg.cga)
cps.obs.exp.tggcgc<-log(dn31.codon.pair$tggcgc/cps.tgg.cgc)
cps.obs.exp.tggcgg<-log(dn31.codon.pair$tggcgg/cps.tgg.cgg)
cps.obs.exp.tggcgt<-log(dn31.codon.pair$tggcgt/cps.tgg.cgt)
cps.obs.exp.tggcta<-log(dn31.codon.pair$tggcta/cps.tgg.cta)
cps.obs.exp.tggctc<-log(dn31.codon.pair$tggctc/cps.tgg.ctc)
cps.obs.exp.tggctg<-log(dn31.codon.pair$tggctg/cps.tgg.ctg)
cps.obs.exp.tggctt<-log(dn31.codon.pair$tggctt/cps.tgg.ctt)

cps.obs.exp.tgggaa<-log(dn31.codon.pair$tgggaa/cps.tgg.gaa)
cps.obs.exp.tgggac<-log(dn31.codon.pair$tgggac/cps.tgg.gac)
cps.obs.exp.tgggag<-log(dn31.codon.pair$tgggag/cps.tgg.gag)
cps.obs.exp.tgggat<-log(dn31.codon.pair$tgggat/cps.tgg.gat)
cps.obs.exp.tgggca<-log(dn31.codon.pair$tgggca/cps.tgg.gca)
cps.obs.exp.tgggcc<-log(dn31.codon.pair$tgggcc/cps.tgg.gcc)
cps.obs.exp.tgggcg<-log(dn31.codon.pair$tgggcg/cps.tgg.gcg)
cps.obs.exp.tgggct<-log(dn31.codon.pair$tgggct/cps.tgg.gct)
cps.obs.exp.tgggga<-log(dn31.codon.pair$tgggga/cps.tgg.gga)
cps.obs.exp.tggggc<-log(dn31.codon.pair$tggggc/cps.tgg.ggc)
cps.obs.exp.tggggg<-log(dn31.codon.pair$tggggg/cps.tgg.ggg)
cps.obs.exp.tggggt<-log(dn31.codon.pair$tggggt/cps.tgg.ggt)
cps.obs.exp.tgggta<-log(dn31.codon.pair$tgggta/cps.tgg.gta)
cps.obs.exp.tgggtc<-log(dn31.codon.pair$tgggtc/cps.tgg.gtc)
cps.obs.exp.tgggtg<-log(dn31.codon.pair$tgggtg/cps.tgg.gtg)
cps.obs.exp.tgggtt<-log(dn31.codon.pair$tgggtt/cps.tgg.gtt)

#cps.obs.exp.tggtaa<-log(dn31.codon.pair$tggtaa/cps.tgg.taa)
cps.obs.exp.tggtac<-log(dn31.codon.pair$tggtac/cps.tgg.tac)
#cps.obs.exp.tggtag<-log(dn31.codon.pair$tggtag/cps.tgg.tag)
cps.obs.exp.tggtat<-log(dn31.codon.pair$tggtat/cps.tgg.tat)
cps.obs.exp.tggtca<-log(dn31.codon.pair$tggtca/cps.tgg.tca)
cps.obs.exp.tggtcc<-log(dn31.codon.pair$tggtcc/cps.tgg.tcc)
cps.obs.exp.tggtcg<-log(dn31.codon.pair$tggtcg/cps.tgg.tcg)
cps.obs.exp.tggtct<-log(dn31.codon.pair$tggtct/cps.tgg.tct)
#cps.obs.exp.tggtga<-log(dn31.codon.pair$tggtga/cps.tgg.tga)
cps.obs.exp.tggtgc<-log(dn31.codon.pair$tggtgc/cps.tgg.tgc)
cps.obs.exp.tggtgg<-log(dn31.codon.pair$tggtgg/cps.tgg.tgg)
cps.obs.exp.tggtgt<-log(dn31.codon.pair$tggtgt/cps.tgg.tgt)
cps.obs.exp.tggtta<-log(dn31.codon.pair$tggtta/cps.tgg.tta)
cps.obs.exp.tggttc<-log(dn31.codon.pair$tggttc/cps.tgg.ttc)
cps.obs.exp.tggttg<-log(dn31.codon.pair$tggttg/cps.tgg.ttg)
cps.obs.exp.tggttt<-log(dn31.codon.pair$tggttt/cps.tgg.ttt)










cps.obs.exp.tgtaaa<-log(dn31.codon.pair$tgtaaa/cps.tgt.aaa)
cps.obs.exp.tgtaac<-log(dn31.codon.pair$tgtaac/cps.tgt.aac)
cps.obs.exp.tgtaag<-log(dn31.codon.pair$tgtaag/cps.tgt.aag)
cps.obs.exp.tgtaat<-log(dn31.codon.pair$tgtaat/cps.tgt.aat)
cps.obs.exp.tgtaca<-log(dn31.codon.pair$tgtaca/cps.tgt.aca)
cps.obs.exp.tgtacc<-log(dn31.codon.pair$tgtacc/cps.tgt.acc)
cps.obs.exp.tgtacg<-log(dn31.codon.pair$tgtacg/cps.tgt.acg)
cps.obs.exp.tgtact<-log(dn31.codon.pair$tgtact/cps.tgt.act)
cps.obs.exp.tgtaga<-log(dn31.codon.pair$tgtaga/cps.tgt.aga)
cps.obs.exp.tgtagc<-log(dn31.codon.pair$tgtagc/cps.tgt.agc)
cps.obs.exp.tgtagg<-log(dn31.codon.pair$tgtagg/cps.tgt.agg)
cps.obs.exp.tgtagt<-log(dn31.codon.pair$tgtagt/cps.tgt.agt)
cps.obs.exp.tgtata<-log(dn31.codon.pair$tgtata/cps.tgt.ata)
cps.obs.exp.tgtatc<-log(dn31.codon.pair$tgtatc/cps.tgt.atc)
cps.obs.exp.tgtatg<-log(dn31.codon.pair$tgtatg/cps.tgt.atg)
cps.obs.exp.tgtatt<-log(dn31.codon.pair$tgtatt/cps.tgt.att)

cps.obs.exp.tgtcaa<-log(dn31.codon.pair$tgtcaa/cps.tgt.caa)
cps.obs.exp.tgtcac<-log(dn31.codon.pair$tgtcac/cps.tgt.cac)
cps.obs.exp.tgtcag<-log(dn31.codon.pair$tgtcag/cps.tgt.cag)
cps.obs.exp.tgtcat<-log(dn31.codon.pair$tgtcat/cps.tgt.cat)
cps.obs.exp.tgtcca<-log(dn31.codon.pair$tgtcca/cps.tgt.cca)
cps.obs.exp.tgtccc<-log(dn31.codon.pair$tgtccc/cps.tgt.ccc)
cps.obs.exp.tgtccg<-log(dn31.codon.pair$tgtccg/cps.tgt.ccg)
cps.obs.exp.tgtcct<-log(dn31.codon.pair$tgtcct/cps.tgt.cct)
cps.obs.exp.tgtcga<-log(dn31.codon.pair$tgtcga/cps.tgt.cga)
cps.obs.exp.tgtcgc<-log(dn31.codon.pair$tgtcgc/cps.tgt.cgc)
cps.obs.exp.tgtcgg<-log(dn31.codon.pair$tgtcgg/cps.tgt.cgg)
cps.obs.exp.tgtcgt<-log(dn31.codon.pair$tgtcgt/cps.tgt.cgt)
cps.obs.exp.tgtcta<-log(dn31.codon.pair$tgtcta/cps.tgt.cta)
cps.obs.exp.tgtctc<-log(dn31.codon.pair$tgtctc/cps.tgt.ctc)
cps.obs.exp.tgtctg<-log(dn31.codon.pair$tgtctg/cps.tgt.ctg)
cps.obs.exp.tgtctt<-log(dn31.codon.pair$tgtctt/cps.tgt.ctt)

cps.obs.exp.tgtgaa<-log(dn31.codon.pair$tgtgaa/cps.tgt.gaa)
cps.obs.exp.tgtgac<-log(dn31.codon.pair$tgtgac/cps.tgt.gac)
cps.obs.exp.tgtgag<-log(dn31.codon.pair$tgtgag/cps.tgt.gag)
cps.obs.exp.tgtgat<-log(dn31.codon.pair$tgtgat/cps.tgt.gat)
cps.obs.exp.tgtgca<-log(dn31.codon.pair$tgtgca/cps.tgt.gca)
cps.obs.exp.tgtgcc<-log(dn31.codon.pair$tgtgcc/cps.tgt.gcc)
cps.obs.exp.tgtgcg<-log(dn31.codon.pair$tgtgcg/cps.tgt.gcg)
cps.obs.exp.tgtgct<-log(dn31.codon.pair$tgtgct/cps.tgt.gct)
cps.obs.exp.tgtgga<-log(dn31.codon.pair$tgtgga/cps.tgt.gga)
cps.obs.exp.tgtggc<-log(dn31.codon.pair$tgtggc/cps.tgt.ggc)
cps.obs.exp.tgtggg<-log(dn31.codon.pair$tgtggg/cps.tgt.ggg)
cps.obs.exp.tgtggt<-log(dn31.codon.pair$tgtggt/cps.tgt.ggt)
cps.obs.exp.tgtgta<-log(dn31.codon.pair$tgtgta/cps.tgt.gta)
cps.obs.exp.tgtgtc<-log(dn31.codon.pair$tgtgtc/cps.tgt.gtc)
cps.obs.exp.tgtgtg<-log(dn31.codon.pair$tgtgtg/cps.tgt.gtg)
cps.obs.exp.tgtgtt<-log(dn31.codon.pair$tgtgtt/cps.tgt.gtt)

#cps.obs.exp.tgttaa<-log(dn31.codon.pair$tgttaa/cps.tgt.taa)
cps.obs.exp.tgttac<-log(dn31.codon.pair$tgttac/cps.tgt.tac)
#cps.obs.exp.tgttag<-log(dn31.codon.pair$tgttag/cps.tgt.tag)
cps.obs.exp.tgttat<-log(dn31.codon.pair$tgttat/cps.tgt.tat)
cps.obs.exp.tgttca<-log(dn31.codon.pair$tgttca/cps.tgt.tca)
cps.obs.exp.tgttcc<-log(dn31.codon.pair$tgttcc/cps.tgt.tcc)
cps.obs.exp.tgttcg<-log(dn31.codon.pair$tgttcg/cps.tgt.tcg)
cps.obs.exp.tgttct<-log(dn31.codon.pair$tgttct/cps.tgt.tct)
#cps.obs.exp.tgttga<-log(dn31.codon.pair$tgttga/cps.tgt.tga)
cps.obs.exp.tgttgc<-log(dn31.codon.pair$tgttgc/cps.tgt.tgc)
cps.obs.exp.tgttgg<-log(dn31.codon.pair$tgttgg/cps.tgt.tgg)
cps.obs.exp.tgttgt<-log(dn31.codon.pair$tgttgt/cps.tgt.tgt)
cps.obs.exp.tgttta<-log(dn31.codon.pair$tgttta/cps.tgt.tta)
cps.obs.exp.tgtttc<-log(dn31.codon.pair$tgtttc/cps.tgt.ttc)
cps.obs.exp.tgtttg<-log(dn31.codon.pair$tgtttg/cps.tgt.ttg)
cps.obs.exp.tgtttt<-log(dn31.codon.pair$tgtttt/cps.tgt.ttt)




















cps.obs.exp.ttaaaa<-log(dn31.codon.pair$ttaaaa/cps.tta.aaa)
cps.obs.exp.ttaaac<-log(dn31.codon.pair$ttaaac/cps.tta.aac)
cps.obs.exp.ttaaag<-log(dn31.codon.pair$ttaaag/cps.tta.aag)
cps.obs.exp.ttaaat<-log(dn31.codon.pair$ttaaat/cps.tta.aat)
cps.obs.exp.ttaaca<-log(dn31.codon.pair$ttaaca/cps.tta.aca)
cps.obs.exp.ttaacc<-log(dn31.codon.pair$ttaacc/cps.tta.acc)
cps.obs.exp.ttaacg<-log(dn31.codon.pair$ttaacg/cps.tta.acg)
cps.obs.exp.ttaact<-log(dn31.codon.pair$ttaact/cps.tta.act)
cps.obs.exp.ttaaga<-log(dn31.codon.pair$ttaaga/cps.tta.aga)
cps.obs.exp.ttaagc<-log(dn31.codon.pair$ttaagc/cps.tta.agc)
cps.obs.exp.ttaagg<-log(dn31.codon.pair$ttaagg/cps.tta.agg)
cps.obs.exp.ttaagt<-log(dn31.codon.pair$ttaagt/cps.tta.agt)
cps.obs.exp.ttaata<-log(dn31.codon.pair$ttaata/cps.tta.ata)
cps.obs.exp.ttaatc<-log(dn31.codon.pair$ttaatc/cps.tta.atc)
cps.obs.exp.ttaatg<-log(dn31.codon.pair$ttaatg/cps.tta.atg)
cps.obs.exp.ttaatt<-log(dn31.codon.pair$ttaatt/cps.tta.att)

cps.obs.exp.ttacaa<-log(dn31.codon.pair$ttacaa/cps.tta.caa)
cps.obs.exp.ttacac<-log(dn31.codon.pair$ttacac/cps.tta.cac)
cps.obs.exp.ttacag<-log(dn31.codon.pair$ttacag/cps.tta.cag)
cps.obs.exp.ttacat<-log(dn31.codon.pair$ttacat/cps.tta.cat)
cps.obs.exp.ttacca<-log(dn31.codon.pair$ttacca/cps.tta.cca)
cps.obs.exp.ttaccc<-log(dn31.codon.pair$ttaccc/cps.tta.ccc)
cps.obs.exp.ttaccg<-log(dn31.codon.pair$ttaccg/cps.tta.ccg)
cps.obs.exp.ttacct<-log(dn31.codon.pair$ttacct/cps.tta.cct)
cps.obs.exp.ttacga<-log(dn31.codon.pair$ttacga/cps.tta.cga)
cps.obs.exp.ttacgc<-log(dn31.codon.pair$ttacgc/cps.tta.cgc)
cps.obs.exp.ttacgg<-log(dn31.codon.pair$ttacgg/cps.tta.cgg)
cps.obs.exp.ttacgt<-log(dn31.codon.pair$ttacgt/cps.tta.cgt)
cps.obs.exp.ttacta<-log(dn31.codon.pair$ttacta/cps.tta.cta)
cps.obs.exp.ttactc<-log(dn31.codon.pair$ttactc/cps.tta.ctc)
cps.obs.exp.ttactg<-log(dn31.codon.pair$ttactg/cps.tta.ctg)
cps.obs.exp.ttactt<-log(dn31.codon.pair$ttactt/cps.tta.ctt)

cps.obs.exp.ttagaa<-log(dn31.codon.pair$ttagaa/cps.tta.gaa)
cps.obs.exp.ttagac<-log(dn31.codon.pair$ttagac/cps.tta.gac)
cps.obs.exp.ttagag<-log(dn31.codon.pair$ttagag/cps.tta.gag)
cps.obs.exp.ttagat<-log(dn31.codon.pair$ttagat/cps.tta.gat)
cps.obs.exp.ttagca<-log(dn31.codon.pair$ttagca/cps.tta.gca)
cps.obs.exp.ttagcc<-log(dn31.codon.pair$ttagcc/cps.tta.gcc)
cps.obs.exp.ttagcg<-log(dn31.codon.pair$ttagcg/cps.tta.gcg)
cps.obs.exp.ttagct<-log(dn31.codon.pair$ttagct/cps.tta.gct)
cps.obs.exp.ttagga<-log(dn31.codon.pair$ttagga/cps.tta.gga)
cps.obs.exp.ttaggc<-log(dn31.codon.pair$ttaggc/cps.tta.ggc)
cps.obs.exp.ttaggg<-log(dn31.codon.pair$ttaggg/cps.tta.ggg)
cps.obs.exp.ttaggt<-log(dn31.codon.pair$ttaggt/cps.tta.ggt)
cps.obs.exp.ttagta<-log(dn31.codon.pair$ttagta/cps.tta.gta)
cps.obs.exp.ttagtc<-log(dn31.codon.pair$ttagtc/cps.tta.gtc)
cps.obs.exp.ttagtg<-log(dn31.codon.pair$ttagtg/cps.tta.gtg)
cps.obs.exp.ttagtt<-log(dn31.codon.pair$ttagtt/cps.tta.gtt)

#cps.obs.exp.ttataa<-log(dn31.codon.pair$ttataa/cps.tta.taa)
cps.obs.exp.ttatac<-log(dn31.codon.pair$ttatac/cps.tta.tac)
#cps.obs.exp.ttatag<-log(dn31.codon.pair$ttatag/cps.tta.tag)
cps.obs.exp.ttatat<-log(dn31.codon.pair$ttatat/cps.tta.tat)
cps.obs.exp.ttatca<-log(dn31.codon.pair$ttatca/cps.tta.tca)
cps.obs.exp.ttatcc<-log(dn31.codon.pair$ttatcc/cps.tta.tcc)
cps.obs.exp.ttatcg<-log(dn31.codon.pair$ttatcg/cps.tta.tcg)
cps.obs.exp.ttatct<-log(dn31.codon.pair$ttatct/cps.tta.tct)
#cps.obs.exp.ttatga<-log(dn31.codon.pair$ttatga/cps.tta.tga)
cps.obs.exp.ttatgc<-log(dn31.codon.pair$ttatgc/cps.tta.tgc)
cps.obs.exp.ttatgg<-log(dn31.codon.pair$ttatgg/cps.tta.tgg)
cps.obs.exp.ttatgt<-log(dn31.codon.pair$ttatgt/cps.tta.tgt)
cps.obs.exp.ttatta<-log(dn31.codon.pair$ttatta/cps.tta.tta)
cps.obs.exp.ttattc<-log(dn31.codon.pair$ttattc/cps.tta.ttc)
cps.obs.exp.ttattg<-log(dn31.codon.pair$ttattg/cps.tta.ttg)
cps.obs.exp.ttattt<-log(dn31.codon.pair$ttattt/cps.tta.ttt)









cps.obs.exp.ttcaaa<-log(dn31.codon.pair$ttcaaa/cps.ttc.aaa)
cps.obs.exp.ttcaac<-log(dn31.codon.pair$ttcaac/cps.ttc.aac)
cps.obs.exp.ttcaag<-log(dn31.codon.pair$ttcaag/cps.ttc.aag)
cps.obs.exp.ttcaat<-log(dn31.codon.pair$ttcaat/cps.ttc.aat)
cps.obs.exp.ttcaca<-log(dn31.codon.pair$ttcaca/cps.ttc.aca)
cps.obs.exp.ttcacc<-log(dn31.codon.pair$ttcacc/cps.ttc.acc)
cps.obs.exp.ttcacg<-log(dn31.codon.pair$ttcacg/cps.ttc.acg)
cps.obs.exp.ttcact<-log(dn31.codon.pair$ttcact/cps.ttc.act)
cps.obs.exp.ttcaga<-log(dn31.codon.pair$ttcaga/cps.ttc.aga)
cps.obs.exp.ttcagc<-log(dn31.codon.pair$ttcagc/cps.ttc.agc)
cps.obs.exp.ttcagg<-log(dn31.codon.pair$ttcagg/cps.ttc.agg)
cps.obs.exp.ttcagt<-log(dn31.codon.pair$ttcagt/cps.ttc.agt)
cps.obs.exp.ttcata<-log(dn31.codon.pair$ttcata/cps.ttc.ata)
cps.obs.exp.ttcatc<-log(dn31.codon.pair$ttcatc/cps.ttc.atc)
cps.obs.exp.ttcatg<-log(dn31.codon.pair$ttcatg/cps.ttc.atg)
cps.obs.exp.ttcatt<-log(dn31.codon.pair$ttcatt/cps.ttc.att)

cps.obs.exp.ttccaa<-log(dn31.codon.pair$ttccaa/cps.ttc.caa)
cps.obs.exp.ttccac<-log(dn31.codon.pair$ttccac/cps.ttc.cac)
cps.obs.exp.ttccag<-log(dn31.codon.pair$ttccag/cps.ttc.cag)
cps.obs.exp.ttccat<-log(dn31.codon.pair$ttccat/cps.ttc.cat)
cps.obs.exp.ttccca<-log(dn31.codon.pair$ttccca/cps.ttc.cca)
cps.obs.exp.ttcccc<-log(dn31.codon.pair$ttcccc/cps.ttc.ccc)
cps.obs.exp.ttcccg<-log(dn31.codon.pair$ttcccg/cps.ttc.ccg)
cps.obs.exp.ttccct<-log(dn31.codon.pair$ttccct/cps.ttc.cct)
cps.obs.exp.ttccga<-log(dn31.codon.pair$ttccga/cps.ttc.cga)
cps.obs.exp.ttccgc<-log(dn31.codon.pair$ttccgc/cps.ttc.cgc)
cps.obs.exp.ttccgg<-log(dn31.codon.pair$ttccgg/cps.ttc.cgg)
cps.obs.exp.ttccgt<-log(dn31.codon.pair$ttccgt/cps.ttc.cgt)
cps.obs.exp.ttccta<-log(dn31.codon.pair$ttccta/cps.ttc.cta)
cps.obs.exp.ttcctc<-log(dn31.codon.pair$ttcctc/cps.ttc.ctc)
cps.obs.exp.ttcctg<-log(dn31.codon.pair$ttcctg/cps.ttc.ctg)
cps.obs.exp.ttcctt<-log(dn31.codon.pair$ttcctt/cps.ttc.ctt)

cps.obs.exp.ttcgaa<-log(dn31.codon.pair$ttcgaa/cps.ttc.gaa)
cps.obs.exp.ttcgac<-log(dn31.codon.pair$ttcgac/cps.ttc.gac)
cps.obs.exp.ttcgag<-log(dn31.codon.pair$ttcgag/cps.ttc.gag)
cps.obs.exp.ttcgat<-log(dn31.codon.pair$ttcgat/cps.ttc.gat)
cps.obs.exp.ttcgca<-log(dn31.codon.pair$ttcgca/cps.ttc.gca)
cps.obs.exp.ttcgcc<-log(dn31.codon.pair$ttcgcc/cps.ttc.gcc)
cps.obs.exp.ttcgcg<-log(dn31.codon.pair$ttcgcg/cps.ttc.gcg)
cps.obs.exp.ttcgct<-log(dn31.codon.pair$ttcgct/cps.ttc.gct)
cps.obs.exp.ttcgga<-log(dn31.codon.pair$ttcgga/cps.ttc.gga)
cps.obs.exp.ttcggc<-log(dn31.codon.pair$ttcggc/cps.ttc.ggc)
cps.obs.exp.ttcggg<-log(dn31.codon.pair$ttcggg/cps.ttc.ggg)
cps.obs.exp.ttcggt<-log(dn31.codon.pair$ttcggt/cps.ttc.ggt)
cps.obs.exp.ttcgta<-log(dn31.codon.pair$ttcgta/cps.ttc.gta)
cps.obs.exp.ttcgtc<-log(dn31.codon.pair$ttcgtc/cps.ttc.gtc)
cps.obs.exp.ttcgtg<-log(dn31.codon.pair$ttcgtg/cps.ttc.gtg)
cps.obs.exp.ttcgtt<-log(dn31.codon.pair$ttcgtt/cps.ttc.gtt)

#cps.obs.exp.ttctaa<-log(dn31.codon.pair$ttctaa/cps.ttc.taa)
cps.obs.exp.ttctac<-log(dn31.codon.pair$ttctac/cps.ttc.tac)
#cps.obs.exp.ttctag<-log(dn31.codon.pair$ttctag/cps.ttc.tag)
cps.obs.exp.ttctat<-log(dn31.codon.pair$ttctat/cps.ttc.tat)
cps.obs.exp.ttctca<-log(dn31.codon.pair$ttctca/cps.ttc.tca)
cps.obs.exp.ttctcc<-log(dn31.codon.pair$ttctcc/cps.ttc.tcc)
cps.obs.exp.ttctcg<-log(dn31.codon.pair$ttctcg/cps.ttc.tcg)
cps.obs.exp.ttctct<-log(dn31.codon.pair$ttctct/cps.ttc.tct)
#cps.obs.exp.ttctga<-log(dn31.codon.pair$ttctga/cps.ttc.tga)
cps.obs.exp.ttctgc<-log(dn31.codon.pair$ttctgc/cps.ttc.tgc)
cps.obs.exp.ttctgg<-log(dn31.codon.pair$ttctgg/cps.ttc.tgg)
cps.obs.exp.ttctgt<-log(dn31.codon.pair$ttctgt/cps.ttc.tgt)
cps.obs.exp.ttctta<-log(dn31.codon.pair$ttctta/cps.ttc.tta)
cps.obs.exp.ttcttc<-log(dn31.codon.pair$ttcttc/cps.ttc.ttc)
cps.obs.exp.ttcttg<-log(dn31.codon.pair$ttcttg/cps.ttc.ttg)
cps.obs.exp.ttcttt<-log(dn31.codon.pair$ttcttt/cps.ttc.ttt)










cps.obs.exp.ttgaaa<-log(dn31.codon.pair$ttgaaa/cps.ttg.aaa)
cps.obs.exp.ttgaac<-log(dn31.codon.pair$ttgaac/cps.ttg.aac)
cps.obs.exp.ttgaag<-log(dn31.codon.pair$ttgaag/cps.ttg.aag)
cps.obs.exp.ttgaat<-log(dn31.codon.pair$ttgaat/cps.ttg.aat)
cps.obs.exp.ttgaca<-log(dn31.codon.pair$ttgaca/cps.ttg.aca)
cps.obs.exp.ttgacc<-log(dn31.codon.pair$ttgacc/cps.ttg.acc)
cps.obs.exp.ttgacg<-log(dn31.codon.pair$ttgacg/cps.ttg.acg)
cps.obs.exp.ttgact<-log(dn31.codon.pair$ttgact/cps.ttg.act)
cps.obs.exp.ttgaga<-log(dn31.codon.pair$ttgaga/cps.ttg.aga)
cps.obs.exp.ttgagc<-log(dn31.codon.pair$ttgagc/cps.ttg.agc)
cps.obs.exp.ttgagg<-log(dn31.codon.pair$ttgagg/cps.ttg.agg)
cps.obs.exp.ttgagt<-log(dn31.codon.pair$ttgagt/cps.ttg.agt)
cps.obs.exp.ttgata<-log(dn31.codon.pair$ttgata/cps.ttg.ata)
cps.obs.exp.ttgatc<-log(dn31.codon.pair$ttgatc/cps.ttg.atc)
cps.obs.exp.ttgatg<-log(dn31.codon.pair$ttgatg/cps.ttg.atg)
cps.obs.exp.ttgatt<-log(dn31.codon.pair$ttgatt/cps.ttg.att)

cps.obs.exp.ttgcaa<-log(dn31.codon.pair$ttgcaa/cps.ttg.caa)
cps.obs.exp.ttgcac<-log(dn31.codon.pair$ttgcac/cps.ttg.cac)
cps.obs.exp.ttgcag<-log(dn31.codon.pair$ttgcag/cps.ttg.cag)
cps.obs.exp.ttgcat<-log(dn31.codon.pair$ttgcat/cps.ttg.cat)
cps.obs.exp.ttgcca<-log(dn31.codon.pair$ttgcca/cps.ttg.cca)
cps.obs.exp.ttgccc<-log(dn31.codon.pair$ttgccc/cps.ttg.ccc)
cps.obs.exp.ttgccg<-log(dn31.codon.pair$ttgccg/cps.ttg.ccg)
cps.obs.exp.ttgcct<-log(dn31.codon.pair$ttgcct/cps.ttg.cct)
cps.obs.exp.ttgcga<-log(dn31.codon.pair$ttgcga/cps.ttg.cga)
cps.obs.exp.ttgcgc<-log(dn31.codon.pair$ttgcgc/cps.ttg.cgc)
cps.obs.exp.ttgcgg<-log(dn31.codon.pair$ttgcgg/cps.ttg.cgg)
cps.obs.exp.ttgcgt<-log(dn31.codon.pair$ttgcgt/cps.ttg.cgt)
cps.obs.exp.ttgcta<-log(dn31.codon.pair$ttgcta/cps.ttg.cta)
cps.obs.exp.ttgctc<-log(dn31.codon.pair$ttgctc/cps.ttg.ctc)
cps.obs.exp.ttgctg<-log(dn31.codon.pair$ttgctg/cps.ttg.ctg)
cps.obs.exp.ttgctt<-log(dn31.codon.pair$ttgctt/cps.ttg.ctt)

cps.obs.exp.ttggaa<-log(dn31.codon.pair$ttggaa/cps.ttg.gaa)
cps.obs.exp.ttggac<-log(dn31.codon.pair$ttggac/cps.ttg.gac)
cps.obs.exp.ttggag<-log(dn31.codon.pair$ttggag/cps.ttg.gag)
cps.obs.exp.ttggat<-log(dn31.codon.pair$ttggat/cps.ttg.gat)
cps.obs.exp.ttggca<-log(dn31.codon.pair$ttggca/cps.ttg.gca)
cps.obs.exp.ttggcc<-log(dn31.codon.pair$ttggcc/cps.ttg.gcc)
cps.obs.exp.ttggcg<-log(dn31.codon.pair$ttggcg/cps.ttg.gcg)
cps.obs.exp.ttggct<-log(dn31.codon.pair$ttggct/cps.ttg.gct)
cps.obs.exp.ttggga<-log(dn31.codon.pair$ttggga/cps.ttg.gga)
cps.obs.exp.ttgggc<-log(dn31.codon.pair$ttgggc/cps.ttg.ggc)
cps.obs.exp.ttgggg<-log(dn31.codon.pair$ttgggg/cps.ttg.ggg)
cps.obs.exp.ttgggt<-log(dn31.codon.pair$ttgggt/cps.ttg.ggt)
cps.obs.exp.ttggta<-log(dn31.codon.pair$ttggta/cps.ttg.gta)
cps.obs.exp.ttggtc<-log(dn31.codon.pair$ttggtc/cps.ttg.gtc)
cps.obs.exp.ttggtg<-log(dn31.codon.pair$ttggtg/cps.ttg.gtg)
cps.obs.exp.ttggtt<-log(dn31.codon.pair$ttggtt/cps.ttg.gtt)

#cps.obs.exp.ttgtaa<-log(dn31.codon.pair$ttgtaa/cps.ttg.taa)
cps.obs.exp.ttgtac<-log(dn31.codon.pair$ttgtac/cps.ttg.tac)
#cps.obs.exp.ttgtag<-log(dn31.codon.pair$ttgtag/cps.ttg.tag)
cps.obs.exp.ttgtat<-log(dn31.codon.pair$ttgtat/cps.ttg.tat)
cps.obs.exp.ttgtca<-log(dn31.codon.pair$ttgtca/cps.ttg.tca)
cps.obs.exp.ttgtcc<-log(dn31.codon.pair$ttgtcc/cps.ttg.tcc)
cps.obs.exp.ttgtcg<-log(dn31.codon.pair$ttgtcg/cps.ttg.tcg)
cps.obs.exp.ttgtct<-log(dn31.codon.pair$ttgtct/cps.ttg.tct)
#cps.obs.exp.ttgtga<-log(dn31.codon.pair$ttgtga/cps.ttg.tga)
cps.obs.exp.ttgtgc<-log(dn31.codon.pair$ttgtgc/cps.ttg.tgc)
cps.obs.exp.ttgtgg<-log(dn31.codon.pair$ttgtgg/cps.ttg.tgg)
cps.obs.exp.ttgtgt<-log(dn31.codon.pair$ttgtgt/cps.ttg.tgt)
cps.obs.exp.ttgtta<-log(dn31.codon.pair$ttgtta/cps.ttg.tta)
cps.obs.exp.ttgttc<-log(dn31.codon.pair$ttgttc/cps.ttg.ttc)
cps.obs.exp.ttgttg<-log(dn31.codon.pair$ttgttg/cps.ttg.ttg)
cps.obs.exp.ttgttt<-log(dn31.codon.pair$ttgttt/cps.ttg.ttt)









cps.obs.exp.tttaaa<-log(dn31.codon.pair$tttaaa/cps.ttt.aaa)
cps.obs.exp.tttaac<-log(dn31.codon.pair$tttaac/cps.ttt.aac)
cps.obs.exp.tttaag<-log(dn31.codon.pair$tttaag/cps.ttt.aag)
cps.obs.exp.tttaat<-log(dn31.codon.pair$tttaat/cps.ttt.aat)
cps.obs.exp.tttaca<-log(dn31.codon.pair$tttaca/cps.ttt.aca)
cps.obs.exp.tttacc<-log(dn31.codon.pair$tttacc/cps.ttt.acc)
cps.obs.exp.tttacg<-log(dn31.codon.pair$tttacg/cps.ttt.acg)
cps.obs.exp.tttact<-log(dn31.codon.pair$tttact/cps.ttt.act)
cps.obs.exp.tttaga<-log(dn31.codon.pair$tttaga/cps.ttt.aga)
cps.obs.exp.tttagc<-log(dn31.codon.pair$tttagc/cps.ttt.agc)
cps.obs.exp.tttagg<-log(dn31.codon.pair$tttagg/cps.ttt.agg)
cps.obs.exp.tttagt<-log(dn31.codon.pair$tttagt/cps.ttt.agt)
cps.obs.exp.tttata<-log(dn31.codon.pair$tttata/cps.ttt.ata)
cps.obs.exp.tttatc<-log(dn31.codon.pair$tttatc/cps.ttt.atc)
cps.obs.exp.tttatg<-log(dn31.codon.pair$tttatg/cps.ttt.atg)
cps.obs.exp.tttatt<-log(dn31.codon.pair$tttatt/cps.ttt.att)

cps.obs.exp.tttcaa<-log(dn31.codon.pair$tttcaa/cps.ttt.caa)
cps.obs.exp.tttcac<-log(dn31.codon.pair$tttcac/cps.ttt.cac)
cps.obs.exp.tttcag<-log(dn31.codon.pair$tttcag/cps.ttt.cag)
cps.obs.exp.tttcat<-log(dn31.codon.pair$tttcat/cps.ttt.cat)
cps.obs.exp.tttcca<-log(dn31.codon.pair$tttcca/cps.ttt.cca)
cps.obs.exp.tttccc<-log(dn31.codon.pair$tttccc/cps.ttt.ccc)
cps.obs.exp.tttccg<-log(dn31.codon.pair$tttccg/cps.ttt.ccg)
cps.obs.exp.tttcct<-log(dn31.codon.pair$tttcct/cps.ttt.cct)
cps.obs.exp.tttcga<-log(dn31.codon.pair$tttcga/cps.ttt.cga)
cps.obs.exp.tttcgc<-log(dn31.codon.pair$tttcgc/cps.ttt.cgc)
cps.obs.exp.tttcgg<-log(dn31.codon.pair$tttcgg/cps.ttt.cgg)
cps.obs.exp.tttcgt<-log(dn31.codon.pair$tttcgt/cps.ttt.cgt)
cps.obs.exp.tttcta<-log(dn31.codon.pair$tttcta/cps.ttt.cta)
cps.obs.exp.tttctc<-log(dn31.codon.pair$tttctc/cps.ttt.ctc)
cps.obs.exp.tttctg<-log(dn31.codon.pair$tttctg/cps.ttt.ctg)
cps.obs.exp.tttctt<-log(dn31.codon.pair$tttctt/cps.ttt.ctt)

cps.obs.exp.tttgaa<-log(dn31.codon.pair$tttgaa/cps.ttt.gaa)
cps.obs.exp.tttgac<-log(dn31.codon.pair$tttgac/cps.ttt.gac)
cps.obs.exp.tttgag<-log(dn31.codon.pair$tttgag/cps.ttt.gag)
cps.obs.exp.tttgat<-log(dn31.codon.pair$tttgat/cps.ttt.gat)
cps.obs.exp.tttgca<-log(dn31.codon.pair$tttgca/cps.ttt.gca)
cps.obs.exp.tttgcc<-log(dn31.codon.pair$tttgcc/cps.ttt.gcc)
cps.obs.exp.tttgcg<-log(dn31.codon.pair$tttgcg/cps.ttt.gcg)
cps.obs.exp.tttgct<-log(dn31.codon.pair$tttgct/cps.ttt.gct)
cps.obs.exp.tttgga<-log(dn31.codon.pair$tttgga/cps.ttt.gga)
cps.obs.exp.tttggc<-log(dn31.codon.pair$tttggc/cps.ttt.ggc)
cps.obs.exp.tttggg<-log(dn31.codon.pair$tttggg/cps.ttt.ggg)
cps.obs.exp.tttggt<-log(dn31.codon.pair$tttggt/cps.ttt.ggt)
cps.obs.exp.tttgta<-log(dn31.codon.pair$tttgta/cps.ttt.gta)
cps.obs.exp.tttgtc<-log(dn31.codon.pair$tttgtc/cps.ttt.gtc)
cps.obs.exp.tttgtg<-log(dn31.codon.pair$tttgtg/cps.ttt.gtg)
cps.obs.exp.tttgtt<-log(dn31.codon.pair$tttgtt/cps.ttt.gtt)

#cps.obs.exp.ttttaa<-log(dn31.codon.pair$ttttaa/cps.ttt.taa)
cps.obs.exp.ttttac<-log(dn31.codon.pair$ttttac/cps.ttt.tac)
#cps.obs.exp.ttttag<-log(dn31.codon.pair$ttttag/cps.ttt.tag)
cps.obs.exp.ttttat<-log(dn31.codon.pair$ttttat/cps.ttt.tat)
cps.obs.exp.ttttca<-log(dn31.codon.pair$ttttca/cps.ttt.tca)
cps.obs.exp.ttttcc<-log(dn31.codon.pair$ttttcc/cps.ttt.tcc)
cps.obs.exp.ttttcg<-log(dn31.codon.pair$ttttcg/cps.ttt.tcg)
cps.obs.exp.ttttct<-log(dn31.codon.pair$ttttct/cps.ttt.tct)
#cps.obs.exp.ttttga<-log(dn31.codon.pair$ttttga/cps.ttt.tga)
cps.obs.exp.ttttgc<-log(dn31.codon.pair$ttttgc/cps.ttt.tgc)
cps.obs.exp.ttttgg<-log(dn31.codon.pair$ttttgg/cps.ttt.tgg)
cps.obs.exp.ttttgt<-log(dn31.codon.pair$ttttgt/cps.ttt.tgt)
cps.obs.exp.ttttta<-log(dn31.codon.pair$ttttta/cps.ttt.tta)
cps.obs.exp.tttttc<-log(dn31.codon.pair$tttttc/cps.ttt.ttc)
cps.obs.exp.tttttg<-log(dn31.codon.pair$tttttg/cps.ttt.ttg)
cps.obs.exp.tttttt<-log(dn31.codon.pair$tttttt/cps.ttt.ttt)





#Remove cps.obs vector
#rm(list=ls(pattern="cps.obs"))
#Show cps vector
ls(pattern="cps.obs")
#Add cps names and value 
names.cps=ls(pattern="cps.obs")
#Size of cps
length(names.cps)
#Create a table with the right size
cps.df=data.frame(matrix(rep(NA,times=10001*length(names.cps)),ncol=length(names.cps)))
#Get all cps valeu and add in a data frame
for(i in 1:length(names.cps)){
cps.df[,i]=get(names.cps[i])}
#Get col names to data frame 
names(cps.df)=names.cps

#cps.df[is.nan(cps.df)] <- 0

cps.df.final<-do.call(data.frame,lapply(cps.df, function(x) replace(x, is.infinite(x),0)))

save(list=ls(), file="/home2/dmacedod/test/dn3_6deg.RData")


print("saving")



#Save cps value as a table
write.table(cps.df.final, file="/home2/dmacedod/test/dn3_6deg_original.txt")





#CPB score
#Original with stop codon  cpb.dn31<-cps.obs.exp.aaaaaa+cps.obs.exp.aaaaac+cps.obs.exp.aaaaag+cps.obs.exp.aaaaat+cps.obs.exp.aaaaca+cps.obs.exp.aaaacc+cps.obs.exp.aaaacg+cps.obs.exp.aaaact+cps.obs.exp.aaaaga+cps.obs.exp.aaaagc+cps.obs.exp.aaaagg+cps.obs.exp.aaaagt+cps.obs.exp.aaaata+cps.obs.exp.aaaatc+cps.obs.exp.aaaatg+cps.obs.exp.aaaatt+cps.obs.exp.aaacaa+cps.obs.exp.aaacac+cps.obs.exp.aaacag+cps.obs.exp.aaacat+cps.obs.exp.aaacca+cps.obs.exp.aaaccc+cps.obs.exp.aaaccg+cps.obs.exp.aaacct+cps.obs.exp.aaacga+cps.obs.exp.aaacgc+cps.obs.exp.aaacgg+cps.obs.exp.aaacgt+cps.obs.exp.aaacta+cps.obs.exp.aaactc+cps.obs.exp.aaactg+cps.obs.exp.aaactt+cps.obs.exp.aaagaa+cps.obs.exp.aaagac+cps.obs.exp.aaagag+cps.obs.exp.aaagat+cps.obs.exp.aaagca+cps.obs.exp.aaagcc+cps.obs.exp.aaagcg+cps.obs.exp.aaagct+cps.obs.exp.aaagga+cps.obs.exp.aaaggc+cps.obs.exp.aaaggg+cps.obs.exp.aaaggt+cps.obs.exp.aaagta+cps.obs.exp.aaagtc+cps.obs.exp.aaagtg+cps.obs.exp.aaagtt+cps.obs.exp.aaataa+cps.obs.exp.aaatac+cps.obs.exp.aaatag+cps.obs.exp.aaatat+cps.obs.exp.aaatca+cps.obs.exp.aaatcc+cps.obs.exp.aaatcg+cps.obs.exp.aaatct+cps.obs.exp.aaatga+cps.obs.exp.aaatgc+cps.obs.exp.aaatgg+cps.obs.exp.aaatgt+cps.obs.exp.aaatta+cps.obs.exp.aaattc+cps.obs.exp.aaattg+cps.obs.exp.aaattt+cps.obs.exp.aacaaa+cps.obs.exp.aacaac+cps.obs.exp.aacaag+cps.obs.exp.aacaat+cps.obs.exp.aacaca+cps.obs.exp.aacacc+cps.obs.exp.aacacg+cps.obs.exp.aacact+cps.obs.exp.aacaga+cps.obs.exp.aacagc+cps.obs.exp.aacagg+cps.obs.exp.aacagt+cps.obs.exp.aacata+cps.obs.exp.aacatc+cps.obs.exp.aacatg+cps.obs.exp.aacatt+cps.obs.exp.aaccaa+cps.obs.exp.aaccac+cps.obs.exp.aaccag+cps.obs.exp.aaccat+cps.obs.exp.aaccca+cps.obs.exp.aacccc+cps.obs.exp.aacccg+cps.obs.exp.aaccct+cps.obs.exp.aaccga+cps.obs.exp.aaccgc+cps.obs.exp.aaccgg+cps.obs.exp.aaccgt+cps.obs.exp.aaccta+cps.obs.exp.aacctc+cps.obs.exp.aacctg+cps.obs.exp.aacctt+cps.obs.exp.aacgaa+cps.obs.exp.aacgac+cps.obs.exp.aacgag+cps.obs.exp.aacgat+cps.obs.exp.aacgca+cps.obs.exp.aacgcc+cps.obs.exp.aacgcg+cps.obs.exp.aacgct+cps.obs.exp.aacgga+cps.obs.exp.aacggc+cps.obs.exp.aacggg+cps.obs.exp.aacggt+cps.obs.exp.aacgta+cps.obs.exp.aacgtc+cps.obs.exp.aacgtg+cps.obs.exp.aacgtt+cps.obs.exp.aactaa+cps.obs.exp.aactac+cps.obs.exp.aactag+cps.obs.exp.aactat+cps.obs.exp.aactca+cps.obs.exp.aactcc+cps.obs.exp.aactcg+cps.obs.exp.aactct+cps.obs.exp.aactga+cps.obs.exp.aactgc+cps.obs.exp.aactgg+cps.obs.exp.aactgt+cps.obs.exp.aactta+cps.obs.exp.aacttc+cps.obs.exp.aacttg+cps.obs.exp.aacttt+cps.obs.exp.aagaaa+cps.obs.exp.aagaac+cps.obs.exp.aagaag+cps.obs.exp.aagaat+cps.obs.exp.aagaca+cps.obs.exp.aagacc+cps.obs.exp.aagacg+cps.obs.exp.aagact+cps.obs.exp.aagaga+cps.obs.exp.aagagc+cps.obs.exp.aagagg+cps.obs.exp.aagagt+cps.obs.exp.aagata+cps.obs.exp.aagatc+cps.obs.exp.aagatg+cps.obs.exp.aagatt+cps.obs.exp.aagcaa+cps.obs.exp.aagcac+cps.obs.exp.aagcag+cps.obs.exp.aagcat+cps.obs.exp.aagcca+cps.obs.exp.aagccc+cps.obs.exp.aagccg+cps.obs.exp.aagcct+cps.obs.exp.aagcga+cps.obs.exp.aagcgc+cps.obs.exp.aagcgg+cps.obs.exp.aagcgt+cps.obs.exp.aagcta+cps.obs.exp.aagctc+cps.obs.exp.aagctg+cps.obs.exp.aagctt+cps.obs.exp.aaggaa+cps.obs.exp.aaggac+cps.obs.exp.aaggag+cps.obs.exp.aaggat+cps.obs.exp.aaggca+cps.obs.exp.aaggcc+cps.obs.exp.aaggcg+cps.obs.exp.aaggct+cps.obs.exp.aaggga+cps.obs.exp.aagggc+cps.obs.exp.aagggg+cps.obs.exp.aagggt+cps.obs.exp.aaggta+cps.obs.exp.aaggtc+cps.obs.exp.aaggtg+cps.obs.exp.aaggtt+cps.obs.exp.aagtaa+cps.obs.exp.aagtac+cps.obs.exp.aagtag+cps.obs.exp.aagtat+cps.obs.exp.aagtca+cps.obs.exp.aagtcc+cps.obs.exp.aagtcg+cps.obs.exp.aagtct+cps.obs.exp.aagtga+cps.obs.exp.aagtgc+cps.obs.exp.aagtgg+cps.obs.exp.aagtgt+cps.obs.exp.aagtta+cps.obs.exp.aagttc+cps.obs.exp.aagttg+cps.obs.exp.aagttt+cps.obs.exp.aataaa+cps.obs.exp.aataac+cps.obs.exp.aataag+cps.obs.exp.aataat+cps.obs.exp.aataca+cps.obs.exp.aatacc+cps.obs.exp.aatacg+cps.obs.exp.aatact+cps.obs.exp.aataga+cps.obs.exp.aatagc+cps.obs.exp.aatagg+cps.obs.exp.aatagt+cps.obs.exp.aatata+cps.obs.exp.aatatc+cps.obs.exp.aatatg+cps.obs.exp.aatatt+cps.obs.exp.aatcaa+cps.obs.exp.aatcac+cps.obs.exp.aatcag+cps.obs.exp.aatcat+cps.obs.exp.aatcca+cps.obs.exp.aatccc+cps.obs.exp.aatccg+cps.obs.exp.aatcct+cps.obs.exp.aatcga+cps.obs.exp.aatcgc+cps.obs.exp.aatcgg+cps.obs.exp.aatcgt+cps.obs.exp.aatcta+cps.obs.exp.aatctc+cps.obs.exp.aatctg+cps.obs.exp.aatctt+cps.obs.exp.aatgaa+cps.obs.exp.aatgac+cps.obs.exp.aatgag+cps.obs.exp.aatgat+cps.obs.exp.aatgca+cps.obs.exp.aatgcc+cps.obs.exp.aatgcg+cps.obs.exp.aatgct+cps.obs.exp.aatgga+cps.obs.exp.aatggc+cps.obs.exp.aatggg+cps.obs.exp.aatggt+cps.obs.exp.aatgta+cps.obs.exp.aatgtc+cps.obs.exp.aatgtg+cps.obs.exp.aatgtt+cps.obs.exp.aattaa+cps.obs.exp.aattac+cps.obs.exp.aattag+cps.obs.exp.aattat+cps.obs.exp.aattca+cps.obs.exp.aattcc+cps.obs.exp.aattcg+cps.obs.exp.aattct+cps.obs.exp.aattga+cps.obs.exp.aattgc+cps.obs.exp.aattgg+cps.obs.exp.aattgt+cps.obs.exp.aattta+cps.obs.exp.aatttc+cps.obs.exp.aatttg+cps.obs.exp.aatttt+cps.obs.exp.acaaaa+cps.obs.exp.acaaac+cps.obs.exp.acaaag+cps.obs.exp.acaaat+cps.obs.exp.acaaca+cps.obs.exp.acaacc+cps.obs.exp.acaacg+cps.obs.exp.acaact+cps.obs.exp.acaaga+cps.obs.exp.acaagc+cps.obs.exp.acaagg+cps.obs.exp.acaagt+cps.obs.exp.acaata+cps.obs.exp.acaatc+cps.obs.exp.acaatg+cps.obs.exp.acaatt+cps.obs.exp.acacaa+cps.obs.exp.acacac+cps.obs.exp.acacag+cps.obs.exp.acacat+cps.obs.exp.acacca+cps.obs.exp.acaccc+cps.obs.exp.acaccg+cps.obs.exp.acacct+cps.obs.exp.acacga+cps.obs.exp.acacgc+cps.obs.exp.acacgg+cps.obs.exp.acacgt+cps.obs.exp.acacta+cps.obs.exp.acactc+cps.obs.exp.acactg+cps.obs.exp.acactt+cps.obs.exp.acagaa+cps.obs.exp.acagac+cps.obs.exp.acagag+cps.obs.exp.acagat+cps.obs.exp.acagca+cps.obs.exp.acagcc+cps.obs.exp.acagcg+cps.obs.exp.acagct+cps.obs.exp.acagga+cps.obs.exp.acaggc+cps.obs.exp.acaggg+cps.obs.exp.acaggt+cps.obs.exp.acagta+cps.obs.exp.acagtc+cps.obs.exp.acagtg+cps.obs.exp.acagtt+cps.obs.exp.acataa+cps.obs.exp.acatac+cps.obs.exp.acatag+cps.obs.exp.acatat+cps.obs.exp.acatca+cps.obs.exp.acatcc+cps.obs.exp.acatcg+cps.obs.exp.acatct+cps.obs.exp.acatga+cps.obs.exp.acatgc+cps.obs.exp.acatgg+cps.obs.exp.acatgt+cps.obs.exp.acatta+cps.obs.exp.acattc+cps.obs.exp.acattg+cps.obs.exp.acattt+cps.obs.exp.accaaa+cps.obs.exp.accaac+cps.obs.exp.accaag+cps.obs.exp.accaat+cps.obs.exp.accaca+cps.obs.exp.accacc+cps.obs.exp.accacg+cps.obs.exp.accact+cps.obs.exp.accaga+cps.obs.exp.accagc+cps.obs.exp.accagg+cps.obs.exp.accagt+cps.obs.exp.accata+cps.obs.exp.accatc+cps.obs.exp.accatg+cps.obs.exp.accatt+cps.obs.exp.acccaa+cps.obs.exp.acccac+cps.obs.exp.acccag+cps.obs.exp.acccat+cps.obs.exp.acccca+cps.obs.exp.accccc+cps.obs.exp.accccg+cps.obs.exp.acccct+cps.obs.exp.acccga+cps.obs.exp.acccgc+cps.obs.exp.acccgg+cps.obs.exp.acccgt+cps.obs.exp.acccta+cps.obs.exp.accctc+cps.obs.exp.accctg+cps.obs.exp.accctt+cps.obs.exp.accgaa+cps.obs.exp.accgac+cps.obs.exp.accgag+cps.obs.exp.accgat+cps.obs.exp.accgca+cps.obs.exp.accgcc+cps.obs.exp.accgcg+cps.obs.exp.accgct+cps.obs.exp.accgga+cps.obs.exp.accggc+cps.obs.exp.accggg+cps.obs.exp.accggt+cps.obs.exp.accgta+cps.obs.exp.accgtc+cps.obs.exp.accgtg+cps.obs.exp.accgtt+cps.obs.exp.acctaa+cps.obs.exp.acctac+cps.obs.exp.acctag+cps.obs.exp.acctat+cps.obs.exp.acctca+cps.obs.exp.acctcc+cps.obs.exp.acctcg+cps.obs.exp.acctct+cps.obs.exp.acctga+cps.obs.exp.acctgc+cps.obs.exp.acctgg+cps.obs.exp.acctgt+cps.obs.exp.acctta+cps.obs.exp.accttc+cps.obs.exp.accttg+cps.obs.exp.accttt+cps.obs.exp.acgaaa+cps.obs.exp.acgaac+cps.obs.exp.acgaag+cps.obs.exp.acgaat+cps.obs.exp.acgaca+cps.obs.exp.acgacc+cps.obs.exp.acgacg+cps.obs.exp.acgact+cps.obs.exp.acgaga+cps.obs.exp.acgagc+cps.obs.exp.acgagg+cps.obs.exp.acgagt+cps.obs.exp.acgata+cps.obs.exp.acgatc+cps.obs.exp.acgatg+cps.obs.exp.acgatt+cps.obs.exp.acgcaa+cps.obs.exp.acgcac+cps.obs.exp.acgcag+cps.obs.exp.acgcat+cps.obs.exp.acgcca+cps.obs.exp.acgccc+cps.obs.exp.acgccg+cps.obs.exp.acgcct+cps.obs.exp.acgcga+cps.obs.exp.acgcgc+cps.obs.exp.acgcgg+cps.obs.exp.acgcgt+cps.obs.exp.acgcta+cps.obs.exp.acgctc+cps.obs.exp.acgctg+cps.obs.exp.acgctt+cps.obs.exp.acggaa+cps.obs.exp.acggac+cps.obs.exp.acggag+cps.obs.exp.acggat+cps.obs.exp.acggca+cps.obs.exp.acggcc+cps.obs.exp.acggcg+cps.obs.exp.acggct+cps.obs.exp.acggga+cps.obs.exp.acgggc+cps.obs.exp.acgggg+cps.obs.exp.acgggt+cps.obs.exp.acggta+cps.obs.exp.acggtc+cps.obs.exp.acggtg+cps.obs.exp.acggtt+cps.obs.exp.acgtaa+cps.obs.exp.acgtac+cps.obs.exp.acgtag+cps.obs.exp.acgtat+cps.obs.exp.acgtca+cps.obs.exp.acgtcc+cps.obs.exp.acgtcg+cps.obs.exp.acgtct+cps.obs.exp.acgtga+cps.obs.exp.acgtgc+cps.obs.exp.acgtgg+cps.obs.exp.acgtgt+cps.obs.exp.acgtta+cps.obs.exp.acgttc+cps.obs.exp.acgttg+cps.obs.exp.acgttt+cps.obs.exp.actaaa+cps.obs.exp.actaac+cps.obs.exp.actaag+cps.obs.exp.actaat+cps.obs.exp.actaca+cps.obs.exp.actacc+cps.obs.exp.actacg+cps.obs.exp.actact+cps.obs.exp.actaga+cps.obs.exp.actagc+cps.obs.exp.actagg+cps.obs.exp.actagt+cps.obs.exp.actata+cps.obs.exp.actatc+cps.obs.exp.actatg+cps.obs.exp.actatt+cps.obs.exp.actcaa+cps.obs.exp.actcac+cps.obs.exp.actcag+cps.obs.exp.actcat+cps.obs.exp.actcca+cps.obs.exp.actccc+cps.obs.exp.actccg+cps.obs.exp.actcct+cps.obs.exp.actcga+cps.obs.exp.actcgc+cps.obs.exp.actcgg+cps.obs.exp.actcgt+cps.obs.exp.actcta+cps.obs.exp.actctc+cps.obs.exp.actctg+cps.obs.exp.actctt+cps.obs.exp.actgaa+cps.obs.exp.actgac+cps.obs.exp.actgag+cps.obs.exp.actgat+cps.obs.exp.actgca+cps.obs.exp.actgcc+cps.obs.exp.actgcg+cps.obs.exp.actgct+cps.obs.exp.actgga+cps.obs.exp.actggc+cps.obs.exp.actggg+cps.obs.exp.actggt+cps.obs.exp.actgta+cps.obs.exp.actgtc+cps.obs.exp.actgtg+cps.obs.exp.actgtt+cps.obs.exp.acttaa+cps.obs.exp.acttac+cps.obs.exp.acttag+cps.obs.exp.acttat+cps.obs.exp.acttca+cps.obs.exp.acttcc+cps.obs.exp.acttcg+cps.obs.exp.acttct+cps.obs.exp.acttga+cps.obs.exp.acttgc+cps.obs.exp.acttgg+cps.obs.exp.acttgt+cps.obs.exp.acttta+cps.obs.exp.actttc+cps.obs.exp.actttg+cps.obs.exp.actttt+cps.obs.exp.agaaaa+cps.obs.exp.agaaac+cps.obs.exp.agaaag+cps.obs.exp.agaaat+cps.obs.exp.agaaca+cps.obs.exp.agaacc+cps.obs.exp.agaacg+cps.obs.exp.agaact+cps.obs.exp.agaaga+cps.obs.exp.agaagc+cps.obs.exp.agaagg+cps.obs.exp.agaagt+cps.obs.exp.agaata+cps.obs.exp.agaatc+cps.obs.exp.agaatg+cps.obs.exp.agaatt+cps.obs.exp.agacaa+cps.obs.exp.agacac+cps.obs.exp.agacag+cps.obs.exp.agacat+cps.obs.exp.agacca+cps.obs.exp.agaccc+cps.obs.exp.agaccg+cps.obs.exp.agacct+cps.obs.exp.agacga+cps.obs.exp.agacgc+cps.obs.exp.agacgg+cps.obs.exp.agacgt+cps.obs.exp.agacta+cps.obs.exp.agactc+cps.obs.exp.agactg+cps.obs.exp.agactt+cps.obs.exp.agagaa+cps.obs.exp.agagac+cps.obs.exp.agagag+cps.obs.exp.agagat+cps.obs.exp.agagca+cps.obs.exp.agagcc+cps.obs.exp.agagcg+cps.obs.exp.agagct+cps.obs.exp.agagga+cps.obs.exp.agaggc+cps.obs.exp.agaggg+cps.obs.exp.agaggt+cps.obs.exp.agagta+cps.obs.exp.agagtc+cps.obs.exp.agagtg+cps.obs.exp.agagtt+cps.obs.exp.agataa+cps.obs.exp.agatac+cps.obs.exp.agatag+cps.obs.exp.agatat+cps.obs.exp.agatca+cps.obs.exp.agatcc+cps.obs.exp.agatcg+cps.obs.exp.agatct+cps.obs.exp.agatga+cps.obs.exp.agatgc+cps.obs.exp.agatgg+cps.obs.exp.agatgt+cps.obs.exp.agatta+cps.obs.exp.agattc+cps.obs.exp.agattg+cps.obs.exp.agattt+cps.obs.exp.agcaaa+cps.obs.exp.agcaac+cps.obs.exp.agcaag+cps.obs.exp.agcaat+cps.obs.exp.agcaca+cps.obs.exp.agcacc+cps.obs.exp.agcacg+cps.obs.exp.agcact+cps.obs.exp.agcaga+cps.obs.exp.agcagc+cps.obs.exp.agcagg+cps.obs.exp.agcagt+cps.obs.exp.agcata+cps.obs.exp.agcatc+cps.obs.exp.agcatg+cps.obs.exp.agcatt+cps.obs.exp.agccaa+cps.obs.exp.agccac+cps.obs.exp.agccag+cps.obs.exp.agccat+cps.obs.exp.agccca+cps.obs.exp.agcccc+cps.obs.exp.agcccg+cps.obs.exp.agccct+cps.obs.exp.agccga+cps.obs.exp.agccgc+cps.obs.exp.agccgg+cps.obs.exp.agccgt+cps.obs.exp.agccta+cps.obs.exp.agcctc+cps.obs.exp.agcctg+cps.obs.exp.agcctt+cps.obs.exp.agcgaa+cps.obs.exp.agcgac+cps.obs.exp.agcgag+cps.obs.exp.agcgat+cps.obs.exp.agcgca+cps.obs.exp.agcgcc+cps.obs.exp.agcgcg+cps.obs.exp.agcgct+cps.obs.exp.agcgga+cps.obs.exp.agcggc+cps.obs.exp.agcggg+cps.obs.exp.agcggt+cps.obs.exp.agcgta+cps.obs.exp.agcgtc+cps.obs.exp.agcgtg+cps.obs.exp.agcgtt+cps.obs.exp.agctaa+cps.obs.exp.agctac+cps.obs.exp.agctag+cps.obs.exp.agctat+cps.obs.exp.agctca+cps.obs.exp.agctcc+cps.obs.exp.agctcg+cps.obs.exp.agctct+cps.obs.exp.agctga+cps.obs.exp.agctgc+cps.obs.exp.agctgg+cps.obs.exp.agctgt+cps.obs.exp.agctta+cps.obs.exp.agcttc+cps.obs.exp.agcttg+cps.obs.exp.agcttt+cps.obs.exp.aggaaa+cps.obs.exp.aggaac+cps.obs.exp.aggaag+cps.obs.exp.aggaat+cps.obs.exp.aggaca+cps.obs.exp.aggacc+cps.obs.exp.aggacg+cps.obs.exp.aggact+cps.obs.exp.aggaga+cps.obs.exp.aggagc+cps.obs.exp.aggagg+cps.obs.exp.aggagt+cps.obs.exp.aggata+cps.obs.exp.aggatc+cps.obs.exp.aggatg+cps.obs.exp.aggatt+cps.obs.exp.aggcaa+cps.obs.exp.aggcac+cps.obs.exp.aggcag+cps.obs.exp.aggcat+cps.obs.exp.aggcca+cps.obs.exp.aggccc+cps.obs.exp.aggccg+cps.obs.exp.aggcct+cps.obs.exp.aggcga+cps.obs.exp.aggcgc+cps.obs.exp.aggcgg+cps.obs.exp.aggcgt+cps.obs.exp.aggcta+cps.obs.exp.aggctc+cps.obs.exp.aggctg+cps.obs.exp.aggctt+cps.obs.exp.agggaa+cps.obs.exp.agggac+cps.obs.exp.agggag+cps.obs.exp.agggat+cps.obs.exp.agggca+cps.obs.exp.agggcc+cps.obs.exp.agggcg+cps.obs.exp.agggct+cps.obs.exp.agggga+cps.obs.exp.aggggc+cps.obs.exp.aggggg+cps.obs.exp.aggggt+cps.obs.exp.agggta+cps.obs.exp.agggtc+cps.obs.exp.agggtg+cps.obs.exp.agggtt+cps.obs.exp.aggtaa+cps.obs.exp.aggtac+cps.obs.exp.aggtag+cps.obs.exp.aggtat+cps.obs.exp.aggtca+cps.obs.exp.aggtcc+cps.obs.exp.aggtcg+cps.obs.exp.aggtct+cps.obs.exp.aggtga+cps.obs.exp.aggtgc+cps.obs.exp.aggtgg+cps.obs.exp.aggtgt+cps.obs.exp.aggtta+cps.obs.exp.aggttc+cps.obs.exp.aggttg+cps.obs.exp.aggttt+cps.obs.exp.agtaaa+cps.obs.exp.agtaac+cps.obs.exp.agtaag+cps.obs.exp.agtaat+cps.obs.exp.agtaca+cps.obs.exp.agtacc+cps.obs.exp.agtacg+cps.obs.exp.agtact+cps.obs.exp.agtaga+cps.obs.exp.agtagc+cps.obs.exp.agtagg+cps.obs.exp.agtagt+cps.obs.exp.agtata+cps.obs.exp.agtatc+cps.obs.exp.agtatg+cps.obs.exp.agtatt+cps.obs.exp.agtcaa+cps.obs.exp.agtcac+cps.obs.exp.agtcag+cps.obs.exp.agtcat+cps.obs.exp.agtcca+cps.obs.exp.agtccc+cps.obs.exp.agtccg+cps.obs.exp.agtcct+cps.obs.exp.agtcga+cps.obs.exp.agtcgc+cps.obs.exp.agtcgg+cps.obs.exp.agtcgt+cps.obs.exp.agtcta+cps.obs.exp.agtctc+cps.obs.exp.agtctg+cps.obs.exp.agtctt+cps.obs.exp.agtgaa+cps.obs.exp.agtgac+cps.obs.exp.agtgag+cps.obs.exp.agtgat+cps.obs.exp.agtgca+cps.obs.exp.agtgcc+cps.obs.exp.agtgcg+cps.obs.exp.agtgct+cps.obs.exp.agtgga+cps.obs.exp.agtggc+cps.obs.exp.agtggg+cps.obs.exp.agtggt+cps.obs.exp.agtgta+cps.obs.exp.agtgtc+cps.obs.exp.agtgtg+cps.obs.exp.agtgtt+cps.obs.exp.agttaa+cps.obs.exp.agttac+cps.obs.exp.agttag+cps.obs.exp.agttat+cps.obs.exp.agttca+cps.obs.exp.agttcc+cps.obs.exp.agttcg+cps.obs.exp.agttct+cps.obs.exp.agttga+cps.obs.exp.agttgc+cps.obs.exp.agttgg+cps.obs.exp.agttgt+cps.obs.exp.agttta+cps.obs.exp.agtttc+cps.obs.exp.agtttg+cps.obs.exp.agtttt+cps.obs.exp.ataaaa+cps.obs.exp.ataaac+cps.obs.exp.ataaag+cps.obs.exp.ataaat+cps.obs.exp.ataaca+cps.obs.exp.ataacc+cps.obs.exp.ataacg+cps.obs.exp.ataact+cps.obs.exp.ataaga+cps.obs.exp.ataagc+cps.obs.exp.ataagg+cps.obs.exp.ataagt+cps.obs.exp.ataata+cps.obs.exp.ataatc+cps.obs.exp.ataatg+cps.obs.exp.ataatt+cps.obs.exp.atacaa+cps.obs.exp.atacac+cps.obs.exp.atacag+cps.obs.exp.atacat+cps.obs.exp.atacca+cps.obs.exp.ataccc+cps.obs.exp.ataccg+cps.obs.exp.atacct+cps.obs.exp.atacga+cps.obs.exp.atacgc+cps.obs.exp.atacgg+cps.obs.exp.atacgt+cps.obs.exp.atacta+cps.obs.exp.atactc+cps.obs.exp.atactg+cps.obs.exp.atactt+cps.obs.exp.atagaa+cps.obs.exp.atagac+cps.obs.exp.atagag+cps.obs.exp.atagat+cps.obs.exp.atagca+cps.obs.exp.atagcc+cps.obs.exp.atagcg+cps.obs.exp.atagct+cps.obs.exp.atagga+cps.obs.exp.ataggc+cps.obs.exp.ataggg+cps.obs.exp.ataggt+cps.obs.exp.atagta+cps.obs.exp.atagtc+cps.obs.exp.atagtg+cps.obs.exp.atagtt+cps.obs.exp.atataa+cps.obs.exp.atatac+cps.obs.exp.atatag+cps.obs.exp.atatat+cps.obs.exp.atatca+cps.obs.exp.atatcc+cps.obs.exp.atatcg+cps.obs.exp.atatct+cps.obs.exp.atatga+cps.obs.exp.atatgc+cps.obs.exp.atatgg+cps.obs.exp.atatgt+cps.obs.exp.atatta+cps.obs.exp.atattc+cps.obs.exp.atattg+cps.obs.exp.atattt+cps.obs.exp.atcaaa+cps.obs.exp.atcaac+cps.obs.exp.atcaag+cps.obs.exp.atcaat+cps.obs.exp.atcaca+cps.obs.exp.atcacc+cps.obs.exp.atcacg+cps.obs.exp.atcact+cps.obs.exp.atcaga+cps.obs.exp.atcagc+cps.obs.exp.atcagg+cps.obs.exp.atcagt+cps.obs.exp.atcata+cps.obs.exp.atcatc+cps.obs.exp.atcatg+cps.obs.exp.atcatt+cps.obs.exp.atccaa+cps.obs.exp.atccac+cps.obs.exp.atccag+cps.obs.exp.atccat+cps.obs.exp.atccca+cps.obs.exp.atcccc+cps.obs.exp.atcccg+cps.obs.exp.atccct+cps.obs.exp.atccga+cps.obs.exp.atccgc+cps.obs.exp.atccgg+cps.obs.exp.atccgt+cps.obs.exp.atccta+cps.obs.exp.atcctc+cps.obs.exp.atcctg+cps.obs.exp.atcctt+cps.obs.exp.atcgaa+cps.obs.exp.atcgac+cps.obs.exp.atcgag+cps.obs.exp.atcgat+cps.obs.exp.atcgca+cps.obs.exp.atcgcc+cps.obs.exp.atcgcg+cps.obs.exp.atcgct+cps.obs.exp.atcgga+cps.obs.exp.atcggc+cps.obs.exp.atcggg+cps.obs.exp.atcggt+cps.obs.exp.atcgta+cps.obs.exp.atcgtc+cps.obs.exp.atcgtg+cps.obs.exp.atcgtt+cps.obs.exp.atctaa+cps.obs.exp.atctac+cps.obs.exp.atctag+cps.obs.exp.atctat+cps.obs.exp.atctca+cps.obs.exp.atctcc+cps.obs.exp.atctcg+cps.obs.exp.atctct+cps.obs.exp.atctga+cps.obs.exp.atctgc+cps.obs.exp.atctgg+cps.obs.exp.atctgt+cps.obs.exp.atctta+cps.obs.exp.atcttc+cps.obs.exp.atcttg+cps.obs.exp.atcttt+cps.obs.exp.atgaaa+cps.obs.exp.atgaac+cps.obs.exp.atgaag+cps.obs.exp.atgaat+cps.obs.exp.atgaca+cps.obs.exp.atgacc+cps.obs.exp.atgacg+cps.obs.exp.atgact+cps.obs.exp.atgaga+cps.obs.exp.atgagc+cps.obs.exp.atgagg+cps.obs.exp.atgagt+cps.obs.exp.atgata+cps.obs.exp.atgatc+cps.obs.exp.atgatg+cps.obs.exp.atgatt+cps.obs.exp.atgcaa+cps.obs.exp.atgcac+cps.obs.exp.atgcag+cps.obs.exp.atgcat+cps.obs.exp.atgcca+cps.obs.exp.atgccc+cps.obs.exp.atgccg+cps.obs.exp.atgcct+cps.obs.exp.atgcga+cps.obs.exp.atgcgc+cps.obs.exp.atgcgg+cps.obs.exp.atgcgt+cps.obs.exp.atgcta+cps.obs.exp.atgctc+cps.obs.exp.atgctg+cps.obs.exp.atgctt+cps.obs.exp.atggaa+cps.obs.exp.atggac+cps.obs.exp.atggag+cps.obs.exp.atggat+cps.obs.exp.atggca+cps.obs.exp.atggcc+cps.obs.exp.atggcg+cps.obs.exp.atggct+cps.obs.exp.atggga+cps.obs.exp.atgggc+cps.obs.exp.atgggg+cps.obs.exp.atgggt+cps.obs.exp.atggta+cps.obs.exp.atggtc+cps.obs.exp.atggtg+cps.obs.exp.atggtt+cps.obs.exp.atgtaa+cps.obs.exp.atgtac+cps.obs.exp.atgtag+cps.obs.exp.atgtat+cps.obs.exp.atgtca+cps.obs.exp.atgtcc+cps.obs.exp.atgtcg+cps.obs.exp.atgtct+cps.obs.exp.atgtga+cps.obs.exp.atgtgc+cps.obs.exp.atgtgg+cps.obs.exp.atgtgt+cps.obs.exp.atgtta+cps.obs.exp.atgttc+cps.obs.exp.atgttg+cps.obs.exp.atgttt+cps.obs.exp.attaaa+cps.obs.exp.attaac+cps.obs.exp.attaag+cps.obs.exp.attaat+cps.obs.exp.attaca+cps.obs.exp.attacc+cps.obs.exp.attacg+cps.obs.exp.attact+cps.obs.exp.attaga+cps.obs.exp.attagc+cps.obs.exp.attagg+cps.obs.exp.attagt+cps.obs.exp.attata+cps.obs.exp.attatc+cps.obs.exp.attatg+cps.obs.exp.attatt+cps.obs.exp.attcaa+cps.obs.exp.attcac+cps.obs.exp.attcag+cps.obs.exp.attcat+cps.obs.exp.attcca+cps.obs.exp.attccc+cps.obs.exp.attccg+cps.obs.exp.attcct+cps.obs.exp.attcga+cps.obs.exp.attcgc+cps.obs.exp.attcgg+cps.obs.exp.attcgt+cps.obs.exp.attcta+cps.obs.exp.attctc+cps.obs.exp.attctg+cps.obs.exp.attctt+cps.obs.exp.attgaa+cps.obs.exp.attgac+cps.obs.exp.attgag+cps.obs.exp.attgat+cps.obs.exp.attgca+cps.obs.exp.attgcc+cps.obs.exp.attgcg+cps.obs.exp.attgct+cps.obs.exp.attgga+cps.obs.exp.attggc+cps.obs.exp.attggg+cps.obs.exp.attggt+cps.obs.exp.attgta+cps.obs.exp.attgtc+cps.obs.exp.attgtg+cps.obs.exp.attgtt+cps.obs.exp.atttaa+cps.obs.exp.atttac+cps.obs.exp.atttag+cps.obs.exp.atttat+cps.obs.exp.atttca+cps.obs.exp.atttcc+cps.obs.exp.atttcg+cps.obs.exp.atttct+cps.obs.exp.atttga+cps.obs.exp.atttgc+cps.obs.exp.atttgg+cps.obs.exp.atttgt+cps.obs.exp.atttta+cps.obs.exp.attttc+cps.obs.exp.attttg+cps.obs.exp.attttt+cps.obs.exp.caaaaa+cps.obs.exp.caaaac+cps.obs.exp.caaaag+cps.obs.exp.caaaat+cps.obs.exp.caaaca+cps.obs.exp.caaacc+cps.obs.exp.caaacg+cps.obs.exp.caaact+cps.obs.exp.caaaga+cps.obs.exp.caaagc+cps.obs.exp.caaagg+cps.obs.exp.caaagt+cps.obs.exp.caaata+cps.obs.exp.caaatc+cps.obs.exp.caaatg+cps.obs.exp.caaatt+cps.obs.exp.caacaa+cps.obs.exp.caacac+cps.obs.exp.caacag+cps.obs.exp.caacat+cps.obs.exp.caacca+cps.obs.exp.caaccc+cps.obs.exp.caaccg+cps.obs.exp.caacct+cps.obs.exp.caacga+cps.obs.exp.caacgc+cps.obs.exp.caacgg+cps.obs.exp.caacgt+cps.obs.exp.caacta+cps.obs.exp.caactc+cps.obs.exp.caactg+cps.obs.exp.caactt+cps.obs.exp.caagaa+cps.obs.exp.caagac+cps.obs.exp.caagag+cps.obs.exp.caagat+cps.obs.exp.caagca+cps.obs.exp.caagcc+cps.obs.exp.caagcg+cps.obs.exp.caagct+cps.obs.exp.caagga+cps.obs.exp.caaggc+cps.obs.exp.caaggg+cps.obs.exp.caaggt+cps.obs.exp.caagta+cps.obs.exp.caagtc+cps.obs.exp.caagtg+cps.obs.exp.caagtt+cps.obs.exp.caataa+cps.obs.exp.caatac+cps.obs.exp.caatag+cps.obs.exp.caatat+cps.obs.exp.caatca+cps.obs.exp.caatcc+cps.obs.exp.caatcg+cps.obs.exp.caatct+cps.obs.exp.caatga+cps.obs.exp.caatgc+cps.obs.exp.caatgg+cps.obs.exp.caatgt+cps.obs.exp.caatta+cps.obs.exp.caattc+cps.obs.exp.caattg+cps.obs.exp.caattt+cps.obs.exp.cacaaa+cps.obs.exp.cacaac+cps.obs.exp.cacaag+cps.obs.exp.cacaat+cps.obs.exp.cacaca+cps.obs.exp.cacacc+cps.obs.exp.cacacg+cps.obs.exp.cacact+cps.obs.exp.cacaga+cps.obs.exp.cacagc+cps.obs.exp.cacagg+cps.obs.exp.cacagt+cps.obs.exp.cacata+cps.obs.exp.cacatc+cps.obs.exp.cacatg+cps.obs.exp.cacatt+cps.obs.exp.caccaa+cps.obs.exp.caccac+cps.obs.exp.caccag+cps.obs.exp.caccat+cps.obs.exp.caccca+cps.obs.exp.cacccc+cps.obs.exp.cacccg+cps.obs.exp.caccct+cps.obs.exp.caccga+cps.obs.exp.caccgc+cps.obs.exp.caccgg+cps.obs.exp.caccgt+cps.obs.exp.caccta+cps.obs.exp.cacctc+cps.obs.exp.cacctg+cps.obs.exp.cacctt+cps.obs.exp.cacgaa+cps.obs.exp.cacgac+cps.obs.exp.cacgag+cps.obs.exp.cacgat+cps.obs.exp.cacgca+cps.obs.exp.cacgcc+cps.obs.exp.cacgcg+cps.obs.exp.cacgct+cps.obs.exp.cacgga+cps.obs.exp.cacggc+cps.obs.exp.cacggg+cps.obs.exp.cacggt+cps.obs.exp.cacgta+cps.obs.exp.cacgtc+cps.obs.exp.cacgtg+cps.obs.exp.cacgtt+cps.obs.exp.cactaa+cps.obs.exp.cactac+cps.obs.exp.cactag+cps.obs.exp.cactat+cps.obs.exp.cactca+cps.obs.exp.cactcc+cps.obs.exp.cactcg+cps.obs.exp.cactct+cps.obs.exp.cactga+cps.obs.exp.cactgc+cps.obs.exp.cactgg+cps.obs.exp.cactgt+cps.obs.exp.cactta+cps.obs.exp.cacttc+cps.obs.exp.cacttg+cps.obs.exp.cacttt+cps.obs.exp.cagaaa+cps.obs.exp.cagaac+cps.obs.exp.cagaag+cps.obs.exp.cagaat+cps.obs.exp.cagaca+cps.obs.exp.cagacc+cps.obs.exp.cagacg+cps.obs.exp.cagact+cps.obs.exp.cagaga+cps.obs.exp.cagagc+cps.obs.exp.cagagg+cps.obs.exp.cagagt+cps.obs.exp.cagata+cps.obs.exp.cagatc+cps.obs.exp.cagatg+cps.obs.exp.cagatt+cps.obs.exp.cagcaa+cps.obs.exp.cagcac+cps.obs.exp.cagcag+cps.obs.exp.cagcat+cps.obs.exp.cagcca+cps.obs.exp.cagccc+cps.obs.exp.cagccg+cps.obs.exp.cagcct+cps.obs.exp.cagcga+cps.obs.exp.cagcgc+cps.obs.exp.cagcgg+cps.obs.exp.cagcgt+cps.obs.exp.cagcta+cps.obs.exp.cagctc+cps.obs.exp.cagctg+cps.obs.exp.cagctt+cps.obs.exp.caggaa+cps.obs.exp.caggac+cps.obs.exp.caggag+cps.obs.exp.caggat+cps.obs.exp.caggca+cps.obs.exp.caggcc+cps.obs.exp.caggcg+cps.obs.exp.caggct+cps.obs.exp.caggga+cps.obs.exp.cagggc+cps.obs.exp.cagggg+cps.obs.exp.cagggt+cps.obs.exp.caggta+cps.obs.exp.caggtc+cps.obs.exp.caggtg+cps.obs.exp.caggtt+cps.obs.exp.cagtaa+cps.obs.exp.cagtac+cps.obs.exp.cagtag+cps.obs.exp.cagtat+cps.obs.exp.cagtca+cps.obs.exp.cagtcc+cps.obs.exp.cagtcg+cps.obs.exp.cagtct+cps.obs.exp.cagtga+cps.obs.exp.cagtgc+cps.obs.exp.cagtgg+cps.obs.exp.cagtgt+cps.obs.exp.cagtta+cps.obs.exp.cagttc+cps.obs.exp.cagttg+cps.obs.exp.cagttt+cps.obs.exp.cataaa+cps.obs.exp.cataac+cps.obs.exp.cataag+cps.obs.exp.cataat+cps.obs.exp.cataca+cps.obs.exp.catacc+cps.obs.exp.catacg+cps.obs.exp.catact+cps.obs.exp.cataga+cps.obs.exp.catagc+cps.obs.exp.catagg+cps.obs.exp.catagt+cps.obs.exp.catata+cps.obs.exp.catatc+cps.obs.exp.catatg+cps.obs.exp.catatt+cps.obs.exp.catcaa+cps.obs.exp.catcac+cps.obs.exp.catcag+cps.obs.exp.catcat+cps.obs.exp.catcca+cps.obs.exp.catccc+cps.obs.exp.catccg+cps.obs.exp.catcct+cps.obs.exp.catcga+cps.obs.exp.catcgc+cps.obs.exp.catcgg+cps.obs.exp.catcgt+cps.obs.exp.catcta+cps.obs.exp.catctc+cps.obs.exp.catctg+cps.obs.exp.catctt+cps.obs.exp.catgaa+cps.obs.exp.catgac+cps.obs.exp.catgag+cps.obs.exp.catgat+cps.obs.exp.catgca+cps.obs.exp.catgcc+cps.obs.exp.catgcg+cps.obs.exp.catgct+cps.obs.exp.catgga+cps.obs.exp.catggc+cps.obs.exp.catggg+cps.obs.exp.catggt+cps.obs.exp.catgta+cps.obs.exp.catgtc+cps.obs.exp.catgtg+cps.obs.exp.catgtt+cps.obs.exp.cattaa+cps.obs.exp.cattac+cps.obs.exp.cattag+cps.obs.exp.cattat+cps.obs.exp.cattca+cps.obs.exp.cattcc+cps.obs.exp.cattcg+cps.obs.exp.cattct+cps.obs.exp.cattga+cps.obs.exp.cattgc+cps.obs.exp.cattgg+cps.obs.exp.cattgt+cps.obs.exp.cattta+cps.obs.exp.catttc+cps.obs.exp.catttg+cps.obs.exp.catttt+cps.obs.exp.ccaaaa+cps.obs.exp.ccaaac+cps.obs.exp.ccaaag+cps.obs.exp.ccaaat+cps.obs.exp.ccaaca+cps.obs.exp.ccaacc+cps.obs.exp.ccaacg+cps.obs.exp.ccaact+cps.obs.exp.ccaaga+cps.obs.exp.ccaagc+cps.obs.exp.ccaagg+cps.obs.exp.ccaagt+cps.obs.exp.ccaata+cps.obs.exp.ccaatc+cps.obs.exp.ccaatg+cps.obs.exp.ccaatt+cps.obs.exp.ccacaa+cps.obs.exp.ccacac+cps.obs.exp.ccacag+cps.obs.exp.ccacat+cps.obs.exp.ccacca+cps.obs.exp.ccaccc+cps.obs.exp.ccaccg+cps.obs.exp.ccacct+cps.obs.exp.ccacga+cps.obs.exp.ccacgc+cps.obs.exp.ccacgg+cps.obs.exp.ccacgt+cps.obs.exp.ccacta+cps.obs.exp.ccactc+cps.obs.exp.ccactg+cps.obs.exp.ccactt+cps.obs.exp.ccagaa+cps.obs.exp.ccagac+cps.obs.exp.ccagag+cps.obs.exp.ccagat+cps.obs.exp.ccagca+cps.obs.exp.ccagcc+cps.obs.exp.ccagcg+cps.obs.exp.ccagct+cps.obs.exp.ccagga+cps.obs.exp.ccaggc+cps.obs.exp.ccaggg+cps.obs.exp.ccaggt+cps.obs.exp.ccagta+cps.obs.exp.ccagtc+cps.obs.exp.ccagtg+cps.obs.exp.ccagtt+cps.obs.exp.ccataa+cps.obs.exp.ccatac+cps.obs.exp.ccatag+cps.obs.exp.ccatat+cps.obs.exp.ccatca+cps.obs.exp.ccatcc+cps.obs.exp.ccatcg+cps.obs.exp.ccatct+cps.obs.exp.ccatga+cps.obs.exp.ccatgc+cps.obs.exp.ccatgg+cps.obs.exp.ccatgt+cps.obs.exp.ccatta+cps.obs.exp.ccattc+cps.obs.exp.ccattg+cps.obs.exp.ccattt+cps.obs.exp.cccaaa+cps.obs.exp.cccaac+cps.obs.exp.cccaag+cps.obs.exp.cccaat+cps.obs.exp.cccaca+cps.obs.exp.cccacc+cps.obs.exp.cccacg+cps.obs.exp.cccact+cps.obs.exp.cccaga+cps.obs.exp.cccagc+cps.obs.exp.cccagg+cps.obs.exp.cccagt+cps.obs.exp.cccata+cps.obs.exp.cccatc+cps.obs.exp.cccatg+cps.obs.exp.cccatt+cps.obs.exp.ccccaa+cps.obs.exp.ccccac+cps.obs.exp.ccccag+cps.obs.exp.ccccat+cps.obs.exp.ccccca+cps.obs.exp.cccccc+cps.obs.exp.cccccg+cps.obs.exp.ccccct+cps.obs.exp.ccccga+cps.obs.exp.ccccgc+cps.obs.exp.ccccgg+cps.obs.exp.ccccgt+cps.obs.exp.ccccta+cps.obs.exp.cccctc+cps.obs.exp.cccctg+cps.obs.exp.cccctt+cps.obs.exp.cccgaa+cps.obs.exp.cccgac+cps.obs.exp.cccgag+cps.obs.exp.cccgat+cps.obs.exp.cccgca+cps.obs.exp.cccgcc+cps.obs.exp.cccgcg+cps.obs.exp.cccgct+cps.obs.exp.cccgga+cps.obs.exp.cccggc+cps.obs.exp.cccggg+cps.obs.exp.cccggt+cps.obs.exp.cccgta+cps.obs.exp.cccgtc+cps.obs.exp.cccgtg+cps.obs.exp.cccgtt+cps.obs.exp.ccctaa+cps.obs.exp.ccctac+cps.obs.exp.ccctag+cps.obs.exp.ccctat+cps.obs.exp.ccctca+cps.obs.exp.ccctcc+cps.obs.exp.ccctcg+cps.obs.exp.ccctct+cps.obs.exp.ccctga+cps.obs.exp.ccctgc+cps.obs.exp.ccctgg+cps.obs.exp.ccctgt+cps.obs.exp.ccctta+cps.obs.exp.cccttc+cps.obs.exp.cccttg+cps.obs.exp.cccttt+cps.obs.exp.ccgaaa+cps.obs.exp.ccgaac+cps.obs.exp.ccgaag+cps.obs.exp.ccgaat+cps.obs.exp.ccgaca+cps.obs.exp.ccgacc+cps.obs.exp.ccgacg+cps.obs.exp.ccgact+cps.obs.exp.ccgaga+cps.obs.exp.ccgagc+cps.obs.exp.ccgagg+cps.obs.exp.ccgagt+cps.obs.exp.ccgata+cps.obs.exp.ccgatc+cps.obs.exp.ccgatg+cps.obs.exp.ccgatt+cps.obs.exp.ccgcaa+cps.obs.exp.ccgcac+cps.obs.exp.ccgcag+cps.obs.exp.ccgcat+cps.obs.exp.ccgcca+cps.obs.exp.ccgccc+cps.obs.exp.ccgccg+cps.obs.exp.ccgcct+cps.obs.exp.ccgcga+cps.obs.exp.ccgcgc+cps.obs.exp.ccgcgg+cps.obs.exp.ccgcgt+cps.obs.exp.ccgcta+cps.obs.exp.ccgctc+cps.obs.exp.ccgctg+cps.obs.exp.ccgctt+cps.obs.exp.ccggaa+cps.obs.exp.ccggac+cps.obs.exp.ccggag+cps.obs.exp.ccggat+cps.obs.exp.ccggca+cps.obs.exp.ccggcc+cps.obs.exp.ccggcg+cps.obs.exp.ccggct+cps.obs.exp.ccggga+cps.obs.exp.ccgggc+cps.obs.exp.ccgggg+cps.obs.exp.ccgggt+cps.obs.exp.ccggta+cps.obs.exp.ccggtc+cps.obs.exp.ccggtg+cps.obs.exp.ccggtt+cps.obs.exp.ccgtaa+cps.obs.exp.ccgtac+cps.obs.exp.ccgtag+cps.obs.exp.ccgtat+cps.obs.exp.ccgtca+cps.obs.exp.ccgtcc+cps.obs.exp.ccgtcg+cps.obs.exp.ccgtct+cps.obs.exp.ccgtga+cps.obs.exp.ccgtgc+cps.obs.exp.ccgtgg+cps.obs.exp.ccgtgt+cps.obs.exp.ccgtta+cps.obs.exp.ccgttc+cps.obs.exp.ccgttg+cps.obs.exp.ccgttt+cps.obs.exp.cctaaa+cps.obs.exp.cctaac+cps.obs.exp.cctaag+cps.obs.exp.cctaat+cps.obs.exp.cctaca+cps.obs.exp.cctacc+cps.obs.exp.cctacg+cps.obs.exp.cctact+cps.obs.exp.cctaga+cps.obs.exp.cctagc+cps.obs.exp.cctagg+cps.obs.exp.cctagt+cps.obs.exp.cctata+cps.obs.exp.cctatc+cps.obs.exp.cctatg+cps.obs.exp.cctatt+cps.obs.exp.cctcaa+cps.obs.exp.cctcac+cps.obs.exp.cctcag+cps.obs.exp.cctcat+cps.obs.exp.cctcca+cps.obs.exp.cctccc+cps.obs.exp.cctccg+cps.obs.exp.cctcct+cps.obs.exp.cctcga+cps.obs.exp.cctcgc+cps.obs.exp.cctcgg+cps.obs.exp.cctcgt+cps.obs.exp.cctcta+cps.obs.exp.cctctc+cps.obs.exp.cctctg+cps.obs.exp.cctctt+cps.obs.exp.cctgaa+cps.obs.exp.cctgac+cps.obs.exp.cctgag+cps.obs.exp.cctgat+cps.obs.exp.cctgca+cps.obs.exp.cctgcc+cps.obs.exp.cctgcg+cps.obs.exp.cctgct+cps.obs.exp.cctgga+cps.obs.exp.cctggc+cps.obs.exp.cctggg+cps.obs.exp.cctggt+cps.obs.exp.cctgta+cps.obs.exp.cctgtc+cps.obs.exp.cctgtg+cps.obs.exp.cctgtt+cps.obs.exp.ccttaa+cps.obs.exp.ccttac+cps.obs.exp.ccttag+cps.obs.exp.ccttat+cps.obs.exp.ccttca+cps.obs.exp.ccttcc+cps.obs.exp.ccttcg+cps.obs.exp.ccttct+cps.obs.exp.ccttga+cps.obs.exp.ccttgc+cps.obs.exp.ccttgg+cps.obs.exp.ccttgt+cps.obs.exp.ccttta+cps.obs.exp.cctttc+cps.obs.exp.cctttg+cps.obs.exp.cctttt+cps.obs.exp.cgaaaa+cps.obs.exp.cgaaac+cps.obs.exp.cgaaag+cps.obs.exp.cgaaat+cps.obs.exp.cgaaca+cps.obs.exp.cgaacc+cps.obs.exp.cgaacg+cps.obs.exp.cgaact+cps.obs.exp.cgaaga+cps.obs.exp.cgaagc+cps.obs.exp.cgaagg+cps.obs.exp.cgaagt+cps.obs.exp.cgaata+cps.obs.exp.cgaatc+cps.obs.exp.cgaatg+cps.obs.exp.cgaatt+cps.obs.exp.cgacaa+cps.obs.exp.cgacac+cps.obs.exp.cgacag+cps.obs.exp.cgacat+cps.obs.exp.cgacca+cps.obs.exp.cgaccc+cps.obs.exp.cgaccg+cps.obs.exp.cgacct+cps.obs.exp.cgacga+cps.obs.exp.cgacgc+cps.obs.exp.cgacgg+cps.obs.exp.cgacgt+cps.obs.exp.cgacta+cps.obs.exp.cgactc+cps.obs.exp.cgactg+cps.obs.exp.cgactt+cps.obs.exp.cgagaa+cps.obs.exp.cgagac+cps.obs.exp.cgagag+cps.obs.exp.cgagat+cps.obs.exp.cgagca+cps.obs.exp.cgagcc+cps.obs.exp.cgagcg+cps.obs.exp.cgagct+cps.obs.exp.cgagga+cps.obs.exp.cgaggc+cps.obs.exp.cgaggg+cps.obs.exp.cgaggt+cps.obs.exp.cgagta+cps.obs.exp.cgagtc+cps.obs.exp.cgagtg+cps.obs.exp.cgagtt+cps.obs.exp.cgataa+cps.obs.exp.cgatac+cps.obs.exp.cgatag+cps.obs.exp.cgatat+cps.obs.exp.cgatca+cps.obs.exp.cgatcc+cps.obs.exp.cgatcg+cps.obs.exp.cgatct+cps.obs.exp.cgatga+cps.obs.exp.cgatgc+cps.obs.exp.cgatgg+cps.obs.exp.cgatgt+cps.obs.exp.cgatta+cps.obs.exp.cgattc+cps.obs.exp.cgattg+cps.obs.exp.cgattt+cps.obs.exp.cgcaaa+cps.obs.exp.cgcaac+cps.obs.exp.cgcaag+cps.obs.exp.cgcaat+cps.obs.exp.cgcaca+cps.obs.exp.cgcacc+cps.obs.exp.cgcacg+cps.obs.exp.cgcact+cps.obs.exp.cgcaga+cps.obs.exp.cgcagc+cps.obs.exp.cgcagg+cps.obs.exp.cgcagt+cps.obs.exp.cgcata+cps.obs.exp.cgcatc+cps.obs.exp.cgcatg+cps.obs.exp.cgcatt+cps.obs.exp.cgccaa+cps.obs.exp.cgccac+cps.obs.exp.cgccag+cps.obs.exp.cgccat+cps.obs.exp.cgccca+cps.obs.exp.cgcccc+cps.obs.exp.cgcccg+cps.obs.exp.cgccct+cps.obs.exp.cgccga+cps.obs.exp.cgccgc+cps.obs.exp.cgccgg+cps.obs.exp.cgccgt+cps.obs.exp.cgccta+cps.obs.exp.cgcctc+cps.obs.exp.cgcctg+cps.obs.exp.cgcctt+cps.obs.exp.cgcgaa+cps.obs.exp.cgcgac+cps.obs.exp.cgcgag+cps.obs.exp.cgcgat+cps.obs.exp.cgcgca+cps.obs.exp.cgcgcc+cps.obs.exp.cgcgcg+cps.obs.exp.cgcgct+cps.obs.exp.cgcgga+cps.obs.exp.cgcggc+cps.obs.exp.cgcggg+cps.obs.exp.cgcggt+cps.obs.exp.cgcgta+cps.obs.exp.cgcgtc+cps.obs.exp.cgcgtg+cps.obs.exp.cgcgtt+cps.obs.exp.cgctaa+cps.obs.exp.cgctac+cps.obs.exp.cgctag+cps.obs.exp.cgctat+cps.obs.exp.cgctca+cps.obs.exp.cgctcc+cps.obs.exp.cgctcg+cps.obs.exp.cgctct+cps.obs.exp.cgctga+cps.obs.exp.cgctgc+cps.obs.exp.cgctgg+cps.obs.exp.cgctgt+cps.obs.exp.cgctta+cps.obs.exp.cgcttc+cps.obs.exp.cgcttg+cps.obs.exp.cgcttt+cps.obs.exp.cggaaa+cps.obs.exp.cggaac+cps.obs.exp.cggaag+cps.obs.exp.cggaat+cps.obs.exp.cggaca+cps.obs.exp.cggacc+cps.obs.exp.cggacg+cps.obs.exp.cggact+cps.obs.exp.cggaga+cps.obs.exp.cggagc+cps.obs.exp.cggagg+cps.obs.exp.cggagt+cps.obs.exp.cggata+cps.obs.exp.cggatc+cps.obs.exp.cggatg+cps.obs.exp.cggatt+cps.obs.exp.cggcaa+cps.obs.exp.cggcac+cps.obs.exp.cggcag+cps.obs.exp.cggcat+cps.obs.exp.cggcca+cps.obs.exp.cggccc+cps.obs.exp.cggccg+cps.obs.exp.cggcct+cps.obs.exp.cggcga+cps.obs.exp.cggcgc+cps.obs.exp.cggcgg+cps.obs.exp.cggcgt+cps.obs.exp.cggcta+cps.obs.exp.cggctc+cps.obs.exp.cggctg+cps.obs.exp.cggctt+cps.obs.exp.cgggaa+cps.obs.exp.cgggac+cps.obs.exp.cgggag+cps.obs.exp.cgggat+cps.obs.exp.cgggca+cps.obs.exp.cgggcc+cps.obs.exp.cgggcg+cps.obs.exp.cgggct+cps.obs.exp.cgggga+cps.obs.exp.cggggc+cps.obs.exp.cggggg+cps.obs.exp.cggggt+cps.obs.exp.cgggta+cps.obs.exp.cgggtc+cps.obs.exp.cgggtg+cps.obs.exp.cgggtt+cps.obs.exp.cggtaa+cps.obs.exp.cggtac+cps.obs.exp.cggtag+cps.obs.exp.cggtat+cps.obs.exp.cggtca+cps.obs.exp.cggtcc+cps.obs.exp.cggtcg+cps.obs.exp.cggtct+cps.obs.exp.cggtga+cps.obs.exp.cggtgc+cps.obs.exp.cggtgg+cps.obs.exp.cggtgt+cps.obs.exp.cggtta+cps.obs.exp.cggttc+cps.obs.exp.cggttg+cps.obs.exp.cggttt+cps.obs.exp.cgtaaa+cps.obs.exp.cgtaac+cps.obs.exp.cgtaag+cps.obs.exp.cgtaat+cps.obs.exp.cgtaca+cps.obs.exp.cgtacc+cps.obs.exp.cgtacg+cps.obs.exp.cgtact+cps.obs.exp.cgtaga+cps.obs.exp.cgtagc+cps.obs.exp.cgtagg+cps.obs.exp.cgtagt+cps.obs.exp.cgtata+cps.obs.exp.cgtatc+cps.obs.exp.cgtatg+cps.obs.exp.cgtatt+cps.obs.exp.cgtcaa+cps.obs.exp.cgtcac+cps.obs.exp.cgtcag+cps.obs.exp.cgtcat+cps.obs.exp.cgtcca+cps.obs.exp.cgtccc+cps.obs.exp.cgtccg+cps.obs.exp.cgtcct+cps.obs.exp.cgtcga+cps.obs.exp.cgtcgc+cps.obs.exp.cgtcgg+cps.obs.exp.cgtcgt+cps.obs.exp.cgtcta+cps.obs.exp.cgtctc+cps.obs.exp.cgtctg+cps.obs.exp.cgtctt+cps.obs.exp.cgtgaa+cps.obs.exp.cgtgac+cps.obs.exp.cgtgag+cps.obs.exp.cgtgat+cps.obs.exp.cgtgca+cps.obs.exp.cgtgcc+cps.obs.exp.cgtgcg+cps.obs.exp.cgtgct+cps.obs.exp.cgtgga+cps.obs.exp.cgtggc+cps.obs.exp.cgtggg+cps.obs.exp.cgtggt+cps.obs.exp.cgtgta+cps.obs.exp.cgtgtc+cps.obs.exp.cgtgtg+cps.obs.exp.cgtgtt+cps.obs.exp.cgttaa+cps.obs.exp.cgttac+cps.obs.exp.cgttag+cps.obs.exp.cgttat+cps.obs.exp.cgttca+cps.obs.exp.cgttcc+cps.obs.exp.cgttcg+cps.obs.exp.cgttct+cps.obs.exp.cgttga+cps.obs.exp.cgttgc+cps.obs.exp.cgttgg+cps.obs.exp.cgttgt+cps.obs.exp.cgttta+cps.obs.exp.cgtttc+cps.obs.exp.cgtttg+cps.obs.exp.cgtttt+cps.obs.exp.ctaaaa+cps.obs.exp.ctaaac+cps.obs.exp.ctaaag+cps.obs.exp.ctaaat+cps.obs.exp.ctaaca+cps.obs.exp.ctaacc+cps.obs.exp.ctaacg+cps.obs.exp.ctaact+cps.obs.exp.ctaaga+cps.obs.exp.ctaagc+cps.obs.exp.ctaagg+cps.obs.exp.ctaagt+cps.obs.exp.ctaata+cps.obs.exp.ctaatc+cps.obs.exp.ctaatg+cps.obs.exp.ctaatt+cps.obs.exp.ctacaa+cps.obs.exp.ctacac+cps.obs.exp.ctacag+cps.obs.exp.ctacat+cps.obs.exp.ctacca+cps.obs.exp.ctaccc+cps.obs.exp.ctaccg+cps.obs.exp.ctacct+cps.obs.exp.ctacga+cps.obs.exp.ctacgc+cps.obs.exp.ctacgg+cps.obs.exp.ctacgt+cps.obs.exp.ctacta+cps.obs.exp.ctactc+cps.obs.exp.ctactg+cps.obs.exp.ctactt+cps.obs.exp.ctagaa+cps.obs.exp.ctagac+cps.obs.exp.ctagag+cps.obs.exp.ctagat+cps.obs.exp.ctagca+cps.obs.exp.ctagcc+cps.obs.exp.ctagcg+cps.obs.exp.ctagct+cps.obs.exp.ctagga+cps.obs.exp.ctaggc+cps.obs.exp.ctaggg+cps.obs.exp.ctaggt+cps.obs.exp.ctagta+cps.obs.exp.ctagtc+cps.obs.exp.ctagtg+cps.obs.exp.ctagtt+cps.obs.exp.ctataa+cps.obs.exp.ctatac+cps.obs.exp.ctatag+cps.obs.exp.ctatat+cps.obs.exp.ctatca+cps.obs.exp.ctatcc+cps.obs.exp.ctatcg+cps.obs.exp.ctatct+cps.obs.exp.ctatga+cps.obs.exp.ctatgc+cps.obs.exp.ctatgg+cps.obs.exp.ctatgt+cps.obs.exp.ctatta+cps.obs.exp.ctattc+cps.obs.exp.ctattg+cps.obs.exp.ctattt+cps.obs.exp.ctcaaa+cps.obs.exp.ctcaac+cps.obs.exp.ctcaag+cps.obs.exp.ctcaat+cps.obs.exp.ctcaca+cps.obs.exp.ctcacc+cps.obs.exp.ctcacg+cps.obs.exp.ctcact+cps.obs.exp.ctcaga+cps.obs.exp.ctcagc+cps.obs.exp.ctcagg+cps.obs.exp.ctcagt+cps.obs.exp.ctcata+cps.obs.exp.ctcatc+cps.obs.exp.ctcatg+cps.obs.exp.ctcatt+cps.obs.exp.ctccaa+cps.obs.exp.ctccac+cps.obs.exp.ctccag+cps.obs.exp.ctccat+cps.obs.exp.ctccca+cps.obs.exp.ctcccc+cps.obs.exp.ctcccg+cps.obs.exp.ctccct+cps.obs.exp.ctccga+cps.obs.exp.ctccgc+cps.obs.exp.ctccgg+cps.obs.exp.ctccgt+cps.obs.exp.ctccta+cps.obs.exp.ctcctc+cps.obs.exp.ctcctg+cps.obs.exp.ctcctt+cps.obs.exp.ctcgaa+cps.obs.exp.ctcgac+cps.obs.exp.ctcgag+cps.obs.exp.ctcgat+cps.obs.exp.ctcgca+cps.obs.exp.ctcgcc+cps.obs.exp.ctcgcg+cps.obs.exp.ctcgct+cps.obs.exp.ctcgga+cps.obs.exp.ctcggc+cps.obs.exp.ctcggg+cps.obs.exp.ctcggt+cps.obs.exp.ctcgta+cps.obs.exp.ctcgtc+cps.obs.exp.ctcgtg+cps.obs.exp.ctcgtt+cps.obs.exp.ctctaa+cps.obs.exp.ctctac+cps.obs.exp.ctctag+cps.obs.exp.ctctat+cps.obs.exp.ctctca+cps.obs.exp.ctctcc+cps.obs.exp.ctctcg+cps.obs.exp.ctctct+cps.obs.exp.ctctga+cps.obs.exp.ctctgc+cps.obs.exp.ctctgg+cps.obs.exp.ctctgt+cps.obs.exp.ctctta+cps.obs.exp.ctcttc+cps.obs.exp.ctcttg+cps.obs.exp.ctcttt+cps.obs.exp.ctgaaa+cps.obs.exp.ctgaac+cps.obs.exp.ctgaag+cps.obs.exp.ctgaat+cps.obs.exp.ctgaca+cps.obs.exp.ctgacc+cps.obs.exp.ctgacg+cps.obs.exp.ctgact+cps.obs.exp.ctgaga+cps.obs.exp.ctgagc+cps.obs.exp.ctgagg+cps.obs.exp.ctgagt+cps.obs.exp.ctgata+cps.obs.exp.ctgatc+cps.obs.exp.ctgatg+cps.obs.exp.ctgatt+cps.obs.exp.ctgcaa+cps.obs.exp.ctgcac+cps.obs.exp.ctgcag+cps.obs.exp.ctgcat+cps.obs.exp.ctgcca+cps.obs.exp.ctgccc+cps.obs.exp.ctgccg+cps.obs.exp.ctgcct+cps.obs.exp.ctgcga+cps.obs.exp.ctgcgc+cps.obs.exp.ctgcgg+cps.obs.exp.ctgcgt+cps.obs.exp.ctgcta+cps.obs.exp.ctgctc+cps.obs.exp.ctgctg+cps.obs.exp.ctgctt+cps.obs.exp.ctggaa+cps.obs.exp.ctggac+cps.obs.exp.ctggag+cps.obs.exp.ctggat+cps.obs.exp.ctggca+cps.obs.exp.ctggcc+cps.obs.exp.ctggcg+cps.obs.exp.ctggct+cps.obs.exp.ctggga+cps.obs.exp.ctgggc+cps.obs.exp.ctgggg+cps.obs.exp.ctgggt+cps.obs.exp.ctggta+cps.obs.exp.ctggtc+cps.obs.exp.ctggtg+cps.obs.exp.ctggtt+cps.obs.exp.ctgtaa+cps.obs.exp.ctgtac+cps.obs.exp.ctgtag+cps.obs.exp.ctgtat+cps.obs.exp.ctgtca+cps.obs.exp.ctgtcc+cps.obs.exp.ctgtcg+cps.obs.exp.ctgtct+cps.obs.exp.ctgtga+cps.obs.exp.ctgtgc+cps.obs.exp.ctgtgg+cps.obs.exp.ctgtgt+cps.obs.exp.ctgtta+cps.obs.exp.ctgttc+cps.obs.exp.ctgttg+cps.obs.exp.ctgttt+cps.obs.exp.cttaaa+cps.obs.exp.cttaac+cps.obs.exp.cttaag+cps.obs.exp.cttaat+cps.obs.exp.cttaca+cps.obs.exp.cttacc+cps.obs.exp.cttacg+cps.obs.exp.cttact+cps.obs.exp.cttaga+cps.obs.exp.cttagc+cps.obs.exp.cttagg+cps.obs.exp.cttagt+cps.obs.exp.cttata+cps.obs.exp.cttatc+cps.obs.exp.cttatg+cps.obs.exp.cttatt+cps.obs.exp.cttcaa+cps.obs.exp.cttcac+cps.obs.exp.cttcag+cps.obs.exp.cttcat+cps.obs.exp.cttcca+cps.obs.exp.cttccc+cps.obs.exp.cttccg+cps.obs.exp.cttcct+cps.obs.exp.cttcga+cps.obs.exp.cttcgc+cps.obs.exp.cttcgg+cps.obs.exp.cttcgt+cps.obs.exp.cttcta+cps.obs.exp.cttctc+cps.obs.exp.cttctg+cps.obs.exp.cttctt+cps.obs.exp.cttgaa+cps.obs.exp.cttgac+cps.obs.exp.cttgag+cps.obs.exp.cttgat+cps.obs.exp.cttgca+cps.obs.exp.cttgcc+cps.obs.exp.cttgcg+cps.obs.exp.cttgct+cps.obs.exp.cttgga+cps.obs.exp.cttggc+cps.obs.exp.cttggg+cps.obs.exp.cttggt+cps.obs.exp.cttgta+cps.obs.exp.cttgtc+cps.obs.exp.cttgtg+cps.obs.exp.cttgtt+cps.obs.exp.ctttaa+cps.obs.exp.ctttac+cps.obs.exp.ctttag+cps.obs.exp.ctttat+cps.obs.exp.ctttca+cps.obs.exp.ctttcc+cps.obs.exp.ctttcg+cps.obs.exp.ctttct+cps.obs.exp.ctttga+cps.obs.exp.ctttgc+cps.obs.exp.ctttgg+cps.obs.exp.ctttgt+cps.obs.exp.ctttta+cps.obs.exp.cttttc+cps.obs.exp.cttttg+cps.obs.exp.cttttt+cps.obs.exp.gaaaaa+cps.obs.exp.gaaaac+cps.obs.exp.gaaaag+cps.obs.exp.gaaaat+cps.obs.exp.gaaaca+cps.obs.exp.gaaacc+cps.obs.exp.gaaacg+cps.obs.exp.gaaact+cps.obs.exp.gaaaga+cps.obs.exp.gaaagc+cps.obs.exp.gaaagg+cps.obs.exp.gaaagt+cps.obs.exp.gaaata+cps.obs.exp.gaaatc+cps.obs.exp.gaaatg+cps.obs.exp.gaaatt+cps.obs.exp.gaacaa+cps.obs.exp.gaacac+cps.obs.exp.gaacag+cps.obs.exp.gaacat+cps.obs.exp.gaacca+cps.obs.exp.gaaccc+cps.obs.exp.gaaccg+cps.obs.exp.gaacct+cps.obs.exp.gaacga+cps.obs.exp.gaacgc+cps.obs.exp.gaacgg+cps.obs.exp.gaacgt+cps.obs.exp.gaacta+cps.obs.exp.gaactc+cps.obs.exp.gaactg+cps.obs.exp.gaactt+cps.obs.exp.gaagaa+cps.obs.exp.gaagac+cps.obs.exp.gaagag+cps.obs.exp.gaagat+cps.obs.exp.gaagca+cps.obs.exp.gaagcc+cps.obs.exp.gaagcg+cps.obs.exp.gaagct+cps.obs.exp.gaagga+cps.obs.exp.gaaggc+cps.obs.exp.gaaggg+cps.obs.exp.gaaggt+cps.obs.exp.gaagta+cps.obs.exp.gaagtc+cps.obs.exp.gaagtg+cps.obs.exp.gaagtt+cps.obs.exp.gaataa+cps.obs.exp.gaatac+cps.obs.exp.gaatag+cps.obs.exp.gaatat+cps.obs.exp.gaatca+cps.obs.exp.gaatcc+cps.obs.exp.gaatcg+cps.obs.exp.gaatct+cps.obs.exp.gaatga+cps.obs.exp.gaatgc+cps.obs.exp.gaatgg+cps.obs.exp.gaatgt+cps.obs.exp.gaatta+cps.obs.exp.gaattc+cps.obs.exp.gaattg+cps.obs.exp.gaattt+cps.obs.exp.gacaaa+cps.obs.exp.gacaac+cps.obs.exp.gacaag+cps.obs.exp.gacaat+cps.obs.exp.gacaca+cps.obs.exp.gacacc+cps.obs.exp.gacacg+cps.obs.exp.gacact+cps.obs.exp.gacaga+cps.obs.exp.gacagc+cps.obs.exp.gacagg+cps.obs.exp.gacagt+cps.obs.exp.gacata+cps.obs.exp.gacatc+cps.obs.exp.gacatg+cps.obs.exp.gacatt+cps.obs.exp.gaccaa+cps.obs.exp.gaccac+cps.obs.exp.gaccag+cps.obs.exp.gaccat+cps.obs.exp.gaccca+cps.obs.exp.gacccc+cps.obs.exp.gacccg+cps.obs.exp.gaccct+cps.obs.exp.gaccga+cps.obs.exp.gaccgc+cps.obs.exp.gaccgg+cps.obs.exp.gaccgt+cps.obs.exp.gaccta+cps.obs.exp.gacctc+cps.obs.exp.gacctg+cps.obs.exp.gacctt+cps.obs.exp.gacgaa+cps.obs.exp.gacgac+cps.obs.exp.gacgag+cps.obs.exp.gacgat+cps.obs.exp.gacgca+cps.obs.exp.gacgcc+cps.obs.exp.gacgcg+cps.obs.exp.gacgct+cps.obs.exp.gacgga+cps.obs.exp.gacggc+cps.obs.exp.gacggg+cps.obs.exp.gacggt+cps.obs.exp.gacgta+cps.obs.exp.gacgtc+cps.obs.exp.gacgtg+cps.obs.exp.gacgtt+cps.obs.exp.gactaa+cps.obs.exp.gactac+cps.obs.exp.gactag+cps.obs.exp.gactat+cps.obs.exp.gactca+cps.obs.exp.gactcc+cps.obs.exp.gactcg+cps.obs.exp.gactct+cps.obs.exp.gactga+cps.obs.exp.gactgc+cps.obs.exp.gactgg+cps.obs.exp.gactgt+cps.obs.exp.gactta+cps.obs.exp.gacttc+cps.obs.exp.gacttg+cps.obs.exp.gacttt+cps.obs.exp.gagaaa+cps.obs.exp.gagaac+cps.obs.exp.gagaag+cps.obs.exp.gagaat+cps.obs.exp.gagaca+cps.obs.exp.gagacc+cps.obs.exp.gagacg+cps.obs.exp.gagact+cps.obs.exp.gagaga+cps.obs.exp.gagagc+cps.obs.exp.gagagg+cps.obs.exp.gagagt+cps.obs.exp.gagata+cps.obs.exp.gagatc+cps.obs.exp.gagatg+cps.obs.exp.gagatt+cps.obs.exp.gagcaa+cps.obs.exp.gagcac+cps.obs.exp.gagcag+cps.obs.exp.gagcat+cps.obs.exp.gagcca+cps.obs.exp.gagccc+cps.obs.exp.gagccg+cps.obs.exp.gagcct+cps.obs.exp.gagcga+cps.obs.exp.gagcgc+cps.obs.exp.gagcgg+cps.obs.exp.gagcgt+cps.obs.exp.gagcta+cps.obs.exp.gagctc+cps.obs.exp.gagctg+cps.obs.exp.gagctt+cps.obs.exp.gaggaa+cps.obs.exp.gaggac+cps.obs.exp.gaggag+cps.obs.exp.gaggat+cps.obs.exp.gaggca+cps.obs.exp.gaggcc+cps.obs.exp.gaggcg+cps.obs.exp.gaggct+cps.obs.exp.gaggga+cps.obs.exp.gagggc+cps.obs.exp.gagggg+cps.obs.exp.gagggt+cps.obs.exp.gaggta+cps.obs.exp.gaggtc+cps.obs.exp.gaggtg+cps.obs.exp.gaggtt+cps.obs.exp.gagtaa+cps.obs.exp.gagtac+cps.obs.exp.gagtag+cps.obs.exp.gagtat+cps.obs.exp.gagtca+cps.obs.exp.gagtcc+cps.obs.exp.gagtcg+cps.obs.exp.gagtct+cps.obs.exp.gagtga+cps.obs.exp.gagtgc+cps.obs.exp.gagtgg+cps.obs.exp.gagtgt+cps.obs.exp.gagtta+cps.obs.exp.gagttc+cps.obs.exp.gagttg+cps.obs.exp.gagttt+cps.obs.exp.gataaa+cps.obs.exp.gataac+cps.obs.exp.gataag+cps.obs.exp.gataat+cps.obs.exp.gataca+cps.obs.exp.gatacc+cps.obs.exp.gatacg+cps.obs.exp.gatact+cps.obs.exp.gataga+cps.obs.exp.gatagc+cps.obs.exp.gatagg+cps.obs.exp.gatagt+cps.obs.exp.gatata+cps.obs.exp.gatatc+cps.obs.exp.gatatg+cps.obs.exp.gatatt+cps.obs.exp.gatcaa+cps.obs.exp.gatcac+cps.obs.exp.gatcag+cps.obs.exp.gatcat+cps.obs.exp.gatcca+cps.obs.exp.gatccc+cps.obs.exp.gatccg+cps.obs.exp.gatcct+cps.obs.exp.gatcga+cps.obs.exp.gatcgc+cps.obs.exp.gatcgg+cps.obs.exp.gatcgt+cps.obs.exp.gatcta+cps.obs.exp.gatctc+cps.obs.exp.gatctg+cps.obs.exp.gatctt+cps.obs.exp.gatgaa+cps.obs.exp.gatgac+cps.obs.exp.gatgag+cps.obs.exp.gatgat+cps.obs.exp.gatgca+cps.obs.exp.gatgcc+cps.obs.exp.gatgcg+cps.obs.exp.gatgct+cps.obs.exp.gatgga+cps.obs.exp.gatggc+cps.obs.exp.gatggg+cps.obs.exp.gatggt+cps.obs.exp.gatgta+cps.obs.exp.gatgtc+cps.obs.exp.gatgtg+cps.obs.exp.gatgtt+cps.obs.exp.gattaa+cps.obs.exp.gattac+cps.obs.exp.gattag+cps.obs.exp.gattat+cps.obs.exp.gattca+cps.obs.exp.gattcc+cps.obs.exp.gattcg+cps.obs.exp.gattct+cps.obs.exp.gattga+cps.obs.exp.gattgc+cps.obs.exp.gattgg+cps.obs.exp.gattgt+cps.obs.exp.gattta+cps.obs.exp.gatttc+cps.obs.exp.gatttg+cps.obs.exp.gatttt+cps.obs.exp.gcaaaa+cps.obs.exp.gcaaac+cps.obs.exp.gcaaag+cps.obs.exp.gcaaat+cps.obs.exp.gcaaca+cps.obs.exp.gcaacc+cps.obs.exp.gcaacg+cps.obs.exp.gcaact+cps.obs.exp.gcaaga+cps.obs.exp.gcaagc+cps.obs.exp.gcaagg+cps.obs.exp.gcaagt+cps.obs.exp.gcaata+cps.obs.exp.gcaatc+cps.obs.exp.gcaatg+cps.obs.exp.gcaatt+cps.obs.exp.gcacaa+cps.obs.exp.gcacac+cps.obs.exp.gcacag+cps.obs.exp.gcacat+cps.obs.exp.gcacca+cps.obs.exp.gcaccc+cps.obs.exp.gcaccg+cps.obs.exp.gcacct+cps.obs.exp.gcacga+cps.obs.exp.gcacgc+cps.obs.exp.gcacgg+cps.obs.exp.gcacgt+cps.obs.exp.gcacta+cps.obs.exp.gcactc+cps.obs.exp.gcactg+cps.obs.exp.gcactt+cps.obs.exp.gcagaa+cps.obs.exp.gcagac+cps.obs.exp.gcagag+cps.obs.exp.gcagat+cps.obs.exp.gcagca+cps.obs.exp.gcagcc+cps.obs.exp.gcagcg+cps.obs.exp.gcagct+cps.obs.exp.gcagga+cps.obs.exp.gcaggc+cps.obs.exp.gcaggg+cps.obs.exp.gcaggt+cps.obs.exp.gcagta+cps.obs.exp.gcagtc+cps.obs.exp.gcagtg+cps.obs.exp.gcagtt+cps.obs.exp.gcataa+cps.obs.exp.gcatac+cps.obs.exp.gcatag+cps.obs.exp.gcatat+cps.obs.exp.gcatca+cps.obs.exp.gcatcc+cps.obs.exp.gcatcg+cps.obs.exp.gcatct+cps.obs.exp.gcatga+cps.obs.exp.gcatgc+cps.obs.exp.gcatgg+cps.obs.exp.gcatgt+cps.obs.exp.gcatta+cps.obs.exp.gcattc+cps.obs.exp.gcattg+cps.obs.exp.gcattt+cps.obs.exp.gccaaa+cps.obs.exp.gccaac+cps.obs.exp.gccaag+cps.obs.exp.gccaat+cps.obs.exp.gccaca+cps.obs.exp.gccacc+cps.obs.exp.gccacg+cps.obs.exp.gccact+cps.obs.exp.gccaga+cps.obs.exp.gccagc+cps.obs.exp.gccagg+cps.obs.exp.gccagt+cps.obs.exp.gccata+cps.obs.exp.gccatc+cps.obs.exp.gccatg+cps.obs.exp.gccatt+cps.obs.exp.gcccaa+cps.obs.exp.gcccac+cps.obs.exp.gcccag+cps.obs.exp.gcccat+cps.obs.exp.gcccca+cps.obs.exp.gccccc+cps.obs.exp.gccccg+cps.obs.exp.gcccct+cps.obs.exp.gcccga+cps.obs.exp.gcccgc+cps.obs.exp.gcccgg+cps.obs.exp.gcccgt+cps.obs.exp.gcccta+cps.obs.exp.gccctc+cps.obs.exp.gccctg+cps.obs.exp.gccctt+cps.obs.exp.gccgaa+cps.obs.exp.gccgac+cps.obs.exp.gccgag+cps.obs.exp.gccgat+cps.obs.exp.gccgca+cps.obs.exp.gccgcc+cps.obs.exp.gccgcg+cps.obs.exp.gccgct+cps.obs.exp.gccgga+cps.obs.exp.gccggc+cps.obs.exp.gccggg+cps.obs.exp.gccggt+cps.obs.exp.gccgta+cps.obs.exp.gccgtc+cps.obs.exp.gccgtg+cps.obs.exp.gccgtt+cps.obs.exp.gcctaa+cps.obs.exp.gcctac+cps.obs.exp.gcctag+cps.obs.exp.gcctat+cps.obs.exp.gcctca+cps.obs.exp.gcctcc+cps.obs.exp.gcctcg+cps.obs.exp.gcctct+cps.obs.exp.gcctga+cps.obs.exp.gcctgc+cps.obs.exp.gcctgg+cps.obs.exp.gcctgt+cps.obs.exp.gcctta+cps.obs.exp.gccttc+cps.obs.exp.gccttg+cps.obs.exp.gccttt+cps.obs.exp.gcgaaa+cps.obs.exp.gcgaac+cps.obs.exp.gcgaag+cps.obs.exp.gcgaat+cps.obs.exp.gcgaca+cps.obs.exp.gcgacc+cps.obs.exp.gcgacg+cps.obs.exp.gcgact+cps.obs.exp.gcgaga+cps.obs.exp.gcgagc+cps.obs.exp.gcgagg+cps.obs.exp.gcgagt+cps.obs.exp.gcgata+cps.obs.exp.gcgatc+cps.obs.exp.gcgatg+cps.obs.exp.gcgatt+cps.obs.exp.gcgcaa+cps.obs.exp.gcgcac+cps.obs.exp.gcgcag+cps.obs.exp.gcgcat+cps.obs.exp.gcgcca+cps.obs.exp.gcgccc+cps.obs.exp.gcgccg+cps.obs.exp.gcgcct+cps.obs.exp.gcgcga+cps.obs.exp.gcgcgc+cps.obs.exp.gcgcgg+cps.obs.exp.gcgcgt+cps.obs.exp.gcgcta+cps.obs.exp.gcgctc+cps.obs.exp.gcgctg+cps.obs.exp.gcgctt+cps.obs.exp.gcggaa+cps.obs.exp.gcggac+cps.obs.exp.gcggag+cps.obs.exp.gcggat+cps.obs.exp.gcggca+cps.obs.exp.gcggcc+cps.obs.exp.gcggcg+cps.obs.exp.gcggct+cps.obs.exp.gcggga+cps.obs.exp.gcgggc+cps.obs.exp.gcgggg+cps.obs.exp.gcgggt+cps.obs.exp.gcggta+cps.obs.exp.gcggtc+cps.obs.exp.gcggtg+cps.obs.exp.gcggtt+cps.obs.exp.gcgtaa+cps.obs.exp.gcgtac+cps.obs.exp.gcgtag+cps.obs.exp.gcgtat+cps.obs.exp.gcgtca+cps.obs.exp.gcgtcc+cps.obs.exp.gcgtcg+cps.obs.exp.gcgtct+cps.obs.exp.gcgtga+cps.obs.exp.gcgtgc+cps.obs.exp.gcgtgg+cps.obs.exp.gcgtgt+cps.obs.exp.gcgtta+cps.obs.exp.gcgttc+cps.obs.exp.gcgttg+cps.obs.exp.gcgttt+cps.obs.exp.gctaaa+cps.obs.exp.gctaac+cps.obs.exp.gctaag+cps.obs.exp.gctaat+cps.obs.exp.gctaca+cps.obs.exp.gctacc+cps.obs.exp.gctacg+cps.obs.exp.gctact+cps.obs.exp.gctaga+cps.obs.exp.gctagc+cps.obs.exp.gctagg+cps.obs.exp.gctagt+cps.obs.exp.gctata+cps.obs.exp.gctatc+cps.obs.exp.gctatg+cps.obs.exp.gctatt+cps.obs.exp.gctcaa+cps.obs.exp.gctcac+cps.obs.exp.gctcag+cps.obs.exp.gctcat+cps.obs.exp.gctcca+cps.obs.exp.gctccc+cps.obs.exp.gctccg+cps.obs.exp.gctcct+cps.obs.exp.gctcga+cps.obs.exp.gctcgc+cps.obs.exp.gctcgg+cps.obs.exp.gctcgt+cps.obs.exp.gctcta+cps.obs.exp.gctctc+cps.obs.exp.gctctg+cps.obs.exp.gctctt+cps.obs.exp.gctgaa+cps.obs.exp.gctgac+cps.obs.exp.gctgag+cps.obs.exp.gctgat+cps.obs.exp.gctgca+cps.obs.exp.gctgcc+cps.obs.exp.gctgcg+cps.obs.exp.gctgct+cps.obs.exp.gctgga+cps.obs.exp.gctggc+cps.obs.exp.gctggg+cps.obs.exp.gctggt+cps.obs.exp.gctgta+cps.obs.exp.gctgtc+cps.obs.exp.gctgtg+cps.obs.exp.gctgtt+cps.obs.exp.gcttaa+cps.obs.exp.gcttac+cps.obs.exp.gcttag+cps.obs.exp.gcttat+cps.obs.exp.gcttca+cps.obs.exp.gcttcc+cps.obs.exp.gcttcg+cps.obs.exp.gcttct+cps.obs.exp.gcttga+cps.obs.exp.gcttgc+cps.obs.exp.gcttgg+cps.obs.exp.gcttgt+cps.obs.exp.gcttta+cps.obs.exp.gctttc+cps.obs.exp.gctttg+cps.obs.exp.gctttt+cps.obs.exp.ggaaaa+cps.obs.exp.ggaaac+cps.obs.exp.ggaaag+cps.obs.exp.ggaaat+cps.obs.exp.ggaaca+cps.obs.exp.ggaacc+cps.obs.exp.ggaacg+cps.obs.exp.ggaact+cps.obs.exp.ggaaga+cps.obs.exp.ggaagc+cps.obs.exp.ggaagg+cps.obs.exp.ggaagt+cps.obs.exp.ggaata+cps.obs.exp.ggaatc+cps.obs.exp.ggaatg+cps.obs.exp.ggaatt+cps.obs.exp.ggacaa+cps.obs.exp.ggacac+cps.obs.exp.ggacag+cps.obs.exp.ggacat+cps.obs.exp.ggacca+cps.obs.exp.ggaccc+cps.obs.exp.ggaccg+cps.obs.exp.ggacct+cps.obs.exp.ggacga+cps.obs.exp.ggacgc+cps.obs.exp.ggacgg+cps.obs.exp.ggacgt+cps.obs.exp.ggacta+cps.obs.exp.ggactc+cps.obs.exp.ggactg+cps.obs.exp.ggactt+cps.obs.exp.ggagaa+cps.obs.exp.ggagac+cps.obs.exp.ggagag+cps.obs.exp.ggagat+cps.obs.exp.ggagca+cps.obs.exp.ggagcc+cps.obs.exp.ggagcg+cps.obs.exp.ggagct+cps.obs.exp.ggagga+cps.obs.exp.ggaggc+cps.obs.exp.ggaggg+cps.obs.exp.ggaggt+cps.obs.exp.ggagta+cps.obs.exp.ggagtc+cps.obs.exp.ggagtg+cps.obs.exp.ggagtt+cps.obs.exp.ggataa+cps.obs.exp.ggatac+cps.obs.exp.ggatag+cps.obs.exp.ggatat+cps.obs.exp.ggatca+cps.obs.exp.ggatcc+cps.obs.exp.ggatcg+cps.obs.exp.ggatct+cps.obs.exp.ggatga+cps.obs.exp.ggatgc+cps.obs.exp.ggatgg+cps.obs.exp.ggatgt+cps.obs.exp.ggatta+cps.obs.exp.ggattc+cps.obs.exp.ggattg+cps.obs.exp.ggattt+cps.obs.exp.ggcaaa+cps.obs.exp.ggcaac+cps.obs.exp.ggcaag+cps.obs.exp.ggcaat+cps.obs.exp.ggcaca+cps.obs.exp.ggcacc+cps.obs.exp.ggcacg+cps.obs.exp.ggcact+cps.obs.exp.ggcaga+cps.obs.exp.ggcagc+cps.obs.exp.ggcagg+cps.obs.exp.ggcagt+cps.obs.exp.ggcata+cps.obs.exp.ggcatc+cps.obs.exp.ggcatg+cps.obs.exp.ggcatt+cps.obs.exp.ggccaa+cps.obs.exp.ggccac+cps.obs.exp.ggccag+cps.obs.exp.ggccat+cps.obs.exp.ggccca+cps.obs.exp.ggcccc+cps.obs.exp.ggcccg+cps.obs.exp.ggccct+cps.obs.exp.ggccga+cps.obs.exp.ggccgc+cps.obs.exp.ggccgg+cps.obs.exp.ggccgt+cps.obs.exp.ggccta+cps.obs.exp.ggcctc+cps.obs.exp.ggcctg+cps.obs.exp.ggcctt+cps.obs.exp.ggcgaa+cps.obs.exp.ggcgac+cps.obs.exp.ggcgag+cps.obs.exp.ggcgat+cps.obs.exp.ggcgca+cps.obs.exp.ggcgcc+cps.obs.exp.ggcgcg+cps.obs.exp.ggcgct+cps.obs.exp.ggcgga+cps.obs.exp.ggcggc+cps.obs.exp.ggcggg+cps.obs.exp.ggcggt+cps.obs.exp.ggcgta+cps.obs.exp.ggcgtc+cps.obs.exp.ggcgtg+cps.obs.exp.ggcgtt+cps.obs.exp.ggctaa+cps.obs.exp.ggctac+cps.obs.exp.ggctag+cps.obs.exp.ggctat+cps.obs.exp.ggctca+cps.obs.exp.ggctcc+cps.obs.exp.ggctcg+cps.obs.exp.ggctct+cps.obs.exp.ggctga+cps.obs.exp.ggctgc+cps.obs.exp.ggctgg+cps.obs.exp.ggctgt+cps.obs.exp.ggctta+cps.obs.exp.ggcttc+cps.obs.exp.ggcttg+cps.obs.exp.ggcttt+cps.obs.exp.gggaaa+cps.obs.exp.gggaac+cps.obs.exp.gggaag+cps.obs.exp.gggaat+cps.obs.exp.gggaca+cps.obs.exp.gggacc+cps.obs.exp.gggacg+cps.obs.exp.gggact+cps.obs.exp.gggaga+cps.obs.exp.gggagc+cps.obs.exp.gggagg+cps.obs.exp.gggagt+cps.obs.exp.gggata+cps.obs.exp.gggatc+cps.obs.exp.gggatg+cps.obs.exp.gggatt+cps.obs.exp.gggcaa+cps.obs.exp.gggcac+cps.obs.exp.gggcag+cps.obs.exp.gggcat+cps.obs.exp.gggcca+cps.obs.exp.gggccc+cps.obs.exp.gggccg+cps.obs.exp.gggcct+cps.obs.exp.gggcga+cps.obs.exp.gggcgc+cps.obs.exp.gggcgg+cps.obs.exp.gggcgt+cps.obs.exp.gggcta+cps.obs.exp.gggctc+cps.obs.exp.gggctg+cps.obs.exp.gggctt+cps.obs.exp.ggggaa+cps.obs.exp.ggggac+cps.obs.exp.ggggag+cps.obs.exp.ggggat+cps.obs.exp.ggggca+cps.obs.exp.ggggcc+cps.obs.exp.ggggcg+cps.obs.exp.ggggct+cps.obs.exp.ggggga+cps.obs.exp.gggggc+cps.obs.exp.gggggg+cps.obs.exp.gggggt+cps.obs.exp.ggggta+cps.obs.exp.ggggtc+cps.obs.exp.ggggtg+cps.obs.exp.ggggtt+cps.obs.exp.gggtaa+cps.obs.exp.gggtac+cps.obs.exp.gggtag+cps.obs.exp.gggtat+cps.obs.exp.gggtca+cps.obs.exp.gggtcc+cps.obs.exp.gggtcg+cps.obs.exp.gggtct+cps.obs.exp.gggtga+cps.obs.exp.gggtgc+cps.obs.exp.gggtgg+cps.obs.exp.gggtgt+cps.obs.exp.gggtta+cps.obs.exp.gggttc+cps.obs.exp.gggttg+cps.obs.exp.gggttt+cps.obs.exp.ggtaaa+cps.obs.exp.ggtaac+cps.obs.exp.ggtaag+cps.obs.exp.ggtaat+cps.obs.exp.ggtaca+cps.obs.exp.ggtacc+cps.obs.exp.ggtacg+cps.obs.exp.ggtact+cps.obs.exp.ggtaga+cps.obs.exp.ggtagc+cps.obs.exp.ggtagg+cps.obs.exp.ggtagt+cps.obs.exp.ggtata+cps.obs.exp.ggtatc+cps.obs.exp.ggtatg+cps.obs.exp.ggtatt+cps.obs.exp.ggtcaa+cps.obs.exp.ggtcac+cps.obs.exp.ggtcag+cps.obs.exp.ggtcat+cps.obs.exp.ggtcca+cps.obs.exp.ggtccc+cps.obs.exp.ggtccg+cps.obs.exp.ggtcct+cps.obs.exp.ggtcga+cps.obs.exp.ggtcgc+cps.obs.exp.ggtcgg+cps.obs.exp.ggtcgt+cps.obs.exp.ggtcta+cps.obs.exp.ggtctc+cps.obs.exp.ggtctg+cps.obs.exp.ggtctt+cps.obs.exp.ggtgaa+cps.obs.exp.ggtgac+cps.obs.exp.ggtgag+cps.obs.exp.ggtgat+cps.obs.exp.ggtgca+cps.obs.exp.ggtgcc+cps.obs.exp.ggtgcg+cps.obs.exp.ggtgct+cps.obs.exp.ggtgga+cps.obs.exp.ggtggc+cps.obs.exp.ggtggg+cps.obs.exp.ggtggt+cps.obs.exp.ggtgta+cps.obs.exp.ggtgtc+cps.obs.exp.ggtgtg+cps.obs.exp.ggtgtt+cps.obs.exp.ggttaa+cps.obs.exp.ggttac+cps.obs.exp.ggttag+cps.obs.exp.ggttat+cps.obs.exp.ggttca+cps.obs.exp.ggttcc+cps.obs.exp.ggttcg+cps.obs.exp.ggttct+cps.obs.exp.ggttga+cps.obs.exp.ggttgc+cps.obs.exp.ggttgg+cps.obs.exp.ggttgt+cps.obs.exp.ggttta+cps.obs.exp.ggtttc+cps.obs.exp.ggtttg+cps.obs.exp.ggtttt+cps.obs.exp.gtaaaa+cps.obs.exp.gtaaac+cps.obs.exp.gtaaag+cps.obs.exp.gtaaat+cps.obs.exp.gtaaca+cps.obs.exp.gtaacc+cps.obs.exp.gtaacg+cps.obs.exp.gtaact+cps.obs.exp.gtaaga+cps.obs.exp.gtaagc+cps.obs.exp.gtaagg+cps.obs.exp.gtaagt+cps.obs.exp.gtaata+cps.obs.exp.gtaatc+cps.obs.exp.gtaatg+cps.obs.exp.gtaatt+cps.obs.exp.gtacaa+cps.obs.exp.gtacac+cps.obs.exp.gtacag+cps.obs.exp.gtacat+cps.obs.exp.gtacca+cps.obs.exp.gtaccc+cps.obs.exp.gtaccg+cps.obs.exp.gtacct+cps.obs.exp.gtacga+cps.obs.exp.gtacgc+cps.obs.exp.gtacgg+cps.obs.exp.gtacgt+cps.obs.exp.gtacta+cps.obs.exp.gtactc+cps.obs.exp.gtactg+cps.obs.exp.gtactt+cps.obs.exp.gtagaa+cps.obs.exp.gtagac+cps.obs.exp.gtagag+cps.obs.exp.gtagat+cps.obs.exp.gtagca+cps.obs.exp.gtagcc+cps.obs.exp.gtagcg+cps.obs.exp.gtagct+cps.obs.exp.gtagga+cps.obs.exp.gtaggc+cps.obs.exp.gtaggg+cps.obs.exp.gtaggt+cps.obs.exp.gtagta+cps.obs.exp.gtagtc+cps.obs.exp.gtagtg+cps.obs.exp.gtagtt+cps.obs.exp.gtataa+cps.obs.exp.gtatac+cps.obs.exp.gtatag+cps.obs.exp.gtatat+cps.obs.exp.gtatca+cps.obs.exp.gtatcc+cps.obs.exp.gtatcg+cps.obs.exp.gtatct+cps.obs.exp.gtatga+cps.obs.exp.gtatgc+cps.obs.exp.gtatgg+cps.obs.exp.gtatgt+cps.obs.exp.gtatta+cps.obs.exp.gtattc+cps.obs.exp.gtattg+cps.obs.exp.gtattt+cps.obs.exp.gtcaaa+cps.obs.exp.gtcaac+cps.obs.exp.gtcaag+cps.obs.exp.gtcaat+cps.obs.exp.gtcaca+cps.obs.exp.gtcacc+cps.obs.exp.gtcacg+cps.obs.exp.gtcact+cps.obs.exp.gtcaga+cps.obs.exp.gtcagc+cps.obs.exp.gtcagg+cps.obs.exp.gtcagt+cps.obs.exp.gtcata+cps.obs.exp.gtcatc+cps.obs.exp.gtcatg+cps.obs.exp.gtcatt+cps.obs.exp.gtccaa+cps.obs.exp.gtccac+cps.obs.exp.gtccag+cps.obs.exp.gtccat+cps.obs.exp.gtccca+cps.obs.exp.gtcccc+cps.obs.exp.gtcccg+cps.obs.exp.gtccct+cps.obs.exp.gtccga+cps.obs.exp.gtccgc+cps.obs.exp.gtccgg+cps.obs.exp.gtccgt+cps.obs.exp.gtccta+cps.obs.exp.gtcctc+cps.obs.exp.gtcctg+cps.obs.exp.gtcctt+cps.obs.exp.gtcgaa+cps.obs.exp.gtcgac+cps.obs.exp.gtcgag+cps.obs.exp.gtcgat+cps.obs.exp.gtcgca+cps.obs.exp.gtcgcc+cps.obs.exp.gtcgcg+cps.obs.exp.gtcgct+cps.obs.exp.gtcgga+cps.obs.exp.gtcggc+cps.obs.exp.gtcggg+cps.obs.exp.gtcggt+cps.obs.exp.gtcgta+cps.obs.exp.gtcgtc+cps.obs.exp.gtcgtg+cps.obs.exp.gtcgtt+cps.obs.exp.gtctaa+cps.obs.exp.gtctac+cps.obs.exp.gtctag+cps.obs.exp.gtctat+cps.obs.exp.gtctca+cps.obs.exp.gtctcc+cps.obs.exp.gtctcg+cps.obs.exp.gtctct+cps.obs.exp.gtctga+cps.obs.exp.gtctgc+cps.obs.exp.gtctgg+cps.obs.exp.gtctgt+cps.obs.exp.gtctta+cps.obs.exp.gtcttc+cps.obs.exp.gtcttg+cps.obs.exp.gtcttt+cps.obs.exp.gtgaaa+cps.obs.exp.gtgaac+cps.obs.exp.gtgaag+cps.obs.exp.gtgaat+cps.obs.exp.gtgaca+cps.obs.exp.gtgacc+cps.obs.exp.gtgacg+cps.obs.exp.gtgact+cps.obs.exp.gtgaga+cps.obs.exp.gtgagc+cps.obs.exp.gtgagg+cps.obs.exp.gtgagt+cps.obs.exp.gtgata+cps.obs.exp.gtgatc+cps.obs.exp.gtgatg+cps.obs.exp.gtgatt+cps.obs.exp.gtgcaa+cps.obs.exp.gtgcac+cps.obs.exp.gtgcag+cps.obs.exp.gtgcat+cps.obs.exp.gtgcca+cps.obs.exp.gtgccc+cps.obs.exp.gtgccg+cps.obs.exp.gtgcct+cps.obs.exp.gtgcga+cps.obs.exp.gtgcgc+cps.obs.exp.gtgcgg+cps.obs.exp.gtgcgt+cps.obs.exp.gtgcta+cps.obs.exp.gtgctc+cps.obs.exp.gtgctg+cps.obs.exp.gtgctt+cps.obs.exp.gtggaa+cps.obs.exp.gtggac+cps.obs.exp.gtggag+cps.obs.exp.gtggat+cps.obs.exp.gtggca+cps.obs.exp.gtggcc+cps.obs.exp.gtggcg+cps.obs.exp.gtggct+cps.obs.exp.gtggga+cps.obs.exp.gtgggc+cps.obs.exp.gtgggg+cps.obs.exp.gtgggt+cps.obs.exp.gtggta+cps.obs.exp.gtggtc+cps.obs.exp.gtggtg+cps.obs.exp.gtggtt+cps.obs.exp.gtgtaa+cps.obs.exp.gtgtac+cps.obs.exp.gtgtag+cps.obs.exp.gtgtat+cps.obs.exp.gtgtca+cps.obs.exp.gtgtcc+cps.obs.exp.gtgtcg+cps.obs.exp.gtgtct+cps.obs.exp.gtgtga+cps.obs.exp.gtgtgc+cps.obs.exp.gtgtgg+cps.obs.exp.gtgtgt+cps.obs.exp.gtgtta+cps.obs.exp.gtgttc+cps.obs.exp.gtgttg+cps.obs.exp.gtgttt+cps.obs.exp.gttaaa+cps.obs.exp.gttaac+cps.obs.exp.gttaag+cps.obs.exp.gttaat+cps.obs.exp.gttaca+cps.obs.exp.gttacc+cps.obs.exp.gttacg+cps.obs.exp.gttact+cps.obs.exp.gttaga+cps.obs.exp.gttagc+cps.obs.exp.gttagg+cps.obs.exp.gttagt+cps.obs.exp.gttata+cps.obs.exp.gttatc+cps.obs.exp.gttatg+cps.obs.exp.gttatt+cps.obs.exp.gttcaa+cps.obs.exp.gttcac+cps.obs.exp.gttcag+cps.obs.exp.gttcat+cps.obs.exp.gttcca+cps.obs.exp.gttccc+cps.obs.exp.gttccg+cps.obs.exp.gttcct+cps.obs.exp.gttcga+cps.obs.exp.gttcgc+cps.obs.exp.gttcgg+cps.obs.exp.gttcgt+cps.obs.exp.gttcta+cps.obs.exp.gttctc+cps.obs.exp.gttctg+cps.obs.exp.gttctt+cps.obs.exp.gttgaa+cps.obs.exp.gttgac+cps.obs.exp.gttgag+cps.obs.exp.gttgat+cps.obs.exp.gttgca+cps.obs.exp.gttgcc+cps.obs.exp.gttgcg+cps.obs.exp.gttgct+cps.obs.exp.gttgga+cps.obs.exp.gttggc+cps.obs.exp.gttggg+cps.obs.exp.gttggt+cps.obs.exp.gttgta+cps.obs.exp.gttgtc+cps.obs.exp.gttgtg+cps.obs.exp.gttgtt+cps.obs.exp.gtttaa+cps.obs.exp.gtttac+cps.obs.exp.gtttag+cps.obs.exp.gtttat+cps.obs.exp.gtttca+cps.obs.exp.gtttcc+cps.obs.exp.gtttcg+cps.obs.exp.gtttct+cps.obs.exp.gtttga+cps.obs.exp.gtttgc+cps.obs.exp.gtttgg+cps.obs.exp.gtttgt+cps.obs.exp.gtttta+cps.obs.exp.gttttc+cps.obs.exp.gttttg+cps.obs.exp.gttttt+cps.obs.exp.taaaaa+cps.obs.exp.taaaac+cps.obs.exp.taaaag+cps.obs.exp.taaaat+cps.obs.exp.taaaca+cps.obs.exp.taaacc+cps.obs.exp.taaacg+cps.obs.exp.taaact+cps.obs.exp.taaaga+cps.obs.exp.taaagc+cps.obs.exp.taaagg+cps.obs.exp.taaagt+cps.obs.exp.taaata+cps.obs.exp.taaatc+cps.obs.exp.taaatg+cps.obs.exp.taaatt+cps.obs.exp.taacaa+cps.obs.exp.taacac+cps.obs.exp.taacag+cps.obs.exp.taacat+cps.obs.exp.taacca+cps.obs.exp.taaccc+cps.obs.exp.taaccg+cps.obs.exp.taacct+cps.obs.exp.taacga+cps.obs.exp.taacgc+cps.obs.exp.taacgg+cps.obs.exp.taacgt+cps.obs.exp.taacta+cps.obs.exp.taactc+cps.obs.exp.taactg+cps.obs.exp.taactt+cps.obs.exp.taagaa+cps.obs.exp.taagac+cps.obs.exp.taagag+cps.obs.exp.taagat+cps.obs.exp.taagca+cps.obs.exp.taagcc+cps.obs.exp.taagcg+cps.obs.exp.taagct+cps.obs.exp.taagga+cps.obs.exp.taaggc+cps.obs.exp.taaggg+cps.obs.exp.taaggt+cps.obs.exp.taagta+cps.obs.exp.taagtc+cps.obs.exp.taagtg+cps.obs.exp.taagtt+cps.obs.exp.taataa+cps.obs.exp.taatac+cps.obs.exp.taatag+cps.obs.exp.taatat+cps.obs.exp.taatca+cps.obs.exp.taatcc+cps.obs.exp.taatcg+cps.obs.exp.taatct+cps.obs.exp.taatga+cps.obs.exp.taatgc+cps.obs.exp.taatgg+cps.obs.exp.taatgt+cps.obs.exp.taatta+cps.obs.exp.taattc+cps.obs.exp.taattg+cps.obs.exp.taattt+cps.obs.exp.tacaaa+cps.obs.exp.tacaac+cps.obs.exp.tacaag+cps.obs.exp.tacaat+cps.obs.exp.tacaca+cps.obs.exp.tacacc+cps.obs.exp.tacacg+cps.obs.exp.tacact+cps.obs.exp.tacaga+cps.obs.exp.tacagc+cps.obs.exp.tacagg+cps.obs.exp.tacagt+cps.obs.exp.tacata+cps.obs.exp.tacatc+cps.obs.exp.tacatg+cps.obs.exp.tacatt+cps.obs.exp.taccaa+cps.obs.exp.taccac+cps.obs.exp.taccag+cps.obs.exp.taccat+cps.obs.exp.taccca+cps.obs.exp.tacccc+cps.obs.exp.tacccg+cps.obs.exp.taccct+cps.obs.exp.taccga+cps.obs.exp.taccgc+cps.obs.exp.taccgg+cps.obs.exp.taccgt+cps.obs.exp.taccta+cps.obs.exp.tacctc+cps.obs.exp.tacctg+cps.obs.exp.tacctt+cps.obs.exp.tacgaa+cps.obs.exp.tacgac+cps.obs.exp.tacgag+cps.obs.exp.tacgat+cps.obs.exp.tacgca+cps.obs.exp.tacgcc+cps.obs.exp.tacgcg+cps.obs.exp.tacgct+cps.obs.exp.tacgga+cps.obs.exp.tacggc+cps.obs.exp.tacggg+cps.obs.exp.tacggt+cps.obs.exp.tacgta+cps.obs.exp.tacgtc+cps.obs.exp.tacgtg+cps.obs.exp.tacgtt+cps.obs.exp.tactaa+cps.obs.exp.tactac+cps.obs.exp.tactag+cps.obs.exp.tactat+cps.obs.exp.tactca+cps.obs.exp.tactcc+cps.obs.exp.tactcg+cps.obs.exp.tactct+cps.obs.exp.tactga+cps.obs.exp.tactgc+cps.obs.exp.tactgg+cps.obs.exp.tactgt+cps.obs.exp.tactta+cps.obs.exp.tacttc+cps.obs.exp.tacttg+cps.obs.exp.tacttt+cps.obs.exp.tagaaa+cps.obs.exp.tagaac+cps.obs.exp.tagaag+cps.obs.exp.tagaat+cps.obs.exp.tagaca+cps.obs.exp.tagacc+cps.obs.exp.tagacg+cps.obs.exp.tagact+cps.obs.exp.tagaga+cps.obs.exp.tagagc+cps.obs.exp.tagagg+cps.obs.exp.tagagt+cps.obs.exp.tagata+cps.obs.exp.tagatc+cps.obs.exp.tagatg+cps.obs.exp.tagatt+cps.obs.exp.tagcaa+cps.obs.exp.tagcac+cps.obs.exp.tagcag+cps.obs.exp.tagcat+cps.obs.exp.tagcca+cps.obs.exp.tagccc+cps.obs.exp.tagccg+cps.obs.exp.tagcct+cps.obs.exp.tagcga+cps.obs.exp.tagcgc+cps.obs.exp.tagcgg+cps.obs.exp.tagcgt+cps.obs.exp.tagcta+cps.obs.exp.tagctc+cps.obs.exp.tagctg+cps.obs.exp.tagctt+cps.obs.exp.taggaa+cps.obs.exp.taggac+cps.obs.exp.taggag+cps.obs.exp.taggat+cps.obs.exp.taggca+cps.obs.exp.taggcc+cps.obs.exp.taggcg+cps.obs.exp.taggct+cps.obs.exp.taggga+cps.obs.exp.tagggc+cps.obs.exp.tagggg+cps.obs.exp.tagggt+cps.obs.exp.taggta+cps.obs.exp.taggtc+cps.obs.exp.taggtg+cps.obs.exp.taggtt+cps.obs.exp.tagtaa+cps.obs.exp.tagtac+cps.obs.exp.tagtag+cps.obs.exp.tagtat+cps.obs.exp.tagtca+cps.obs.exp.tagtcc+cps.obs.exp.tagtcg+cps.obs.exp.tagtct+cps.obs.exp.tagtga+cps.obs.exp.tagtgc+cps.obs.exp.tagtgg+cps.obs.exp.tagtgt+cps.obs.exp.tagtta+cps.obs.exp.tagttc+cps.obs.exp.tagttg+cps.obs.exp.tagttt+cps.obs.exp.tataaa+cps.obs.exp.tataac+cps.obs.exp.tataag+cps.obs.exp.tataat+cps.obs.exp.tataca+cps.obs.exp.tatacc+cps.obs.exp.tatacg+cps.obs.exp.tatact+cps.obs.exp.tataga+cps.obs.exp.tatagc+cps.obs.exp.tatagg+cps.obs.exp.tatagt+cps.obs.exp.tatata+cps.obs.exp.tatatc+cps.obs.exp.tatatg+cps.obs.exp.tatatt+cps.obs.exp.tatcaa+cps.obs.exp.tatcac+cps.obs.exp.tatcag+cps.obs.exp.tatcat+cps.obs.exp.tatcca+cps.obs.exp.tatccc+cps.obs.exp.tatccg+cps.obs.exp.tatcct+cps.obs.exp.tatcga+cps.obs.exp.tatcgc+cps.obs.exp.tatcgg+cps.obs.exp.tatcgt+cps.obs.exp.tatcta+cps.obs.exp.tatctc+cps.obs.exp.tatctg+cps.obs.exp.tatctt+cps.obs.exp.tatgaa+cps.obs.exp.tatgac+cps.obs.exp.tatgag+cps.obs.exp.tatgat+cps.obs.exp.tatgca+cps.obs.exp.tatgcc+cps.obs.exp.tatgcg+cps.obs.exp.tatgct+cps.obs.exp.tatgga+cps.obs.exp.tatggc+cps.obs.exp.tatggg+cps.obs.exp.tatggt+cps.obs.exp.tatgta+cps.obs.exp.tatgtc+cps.obs.exp.tatgtg+cps.obs.exp.tatgtt+cps.obs.exp.tattaa+cps.obs.exp.tattac+cps.obs.exp.tattag+cps.obs.exp.tattat+cps.obs.exp.tattca+cps.obs.exp.tattcc+cps.obs.exp.tattcg+cps.obs.exp.tattct+cps.obs.exp.tattga+cps.obs.exp.tattgc+cps.obs.exp.tattgg+cps.obs.exp.tattgt+cps.obs.exp.tattta+cps.obs.exp.tatttc+cps.obs.exp.tatttg+cps.obs.exp.tatttt+cps.obs.exp.tcaaaa+cps.obs.exp.tcaaac+cps.obs.exp.tcaaag+cps.obs.exp.tcaaat+cps.obs.exp.tcaaca+cps.obs.exp.tcaacc+cps.obs.exp.tcaacg+cps.obs.exp.tcaact+cps.obs.exp.tcaaga+cps.obs.exp.tcaagc+cps.obs.exp.tcaagg+cps.obs.exp.tcaagt+cps.obs.exp.tcaata+cps.obs.exp.tcaatc+cps.obs.exp.tcaatg+cps.obs.exp.tcaatt+cps.obs.exp.tcacaa+cps.obs.exp.tcacac+cps.obs.exp.tcacag+cps.obs.exp.tcacat+cps.obs.exp.tcacca+cps.obs.exp.tcaccc+cps.obs.exp.tcaccg+cps.obs.exp.tcacct+cps.obs.exp.tcacga+cps.obs.exp.tcacgc+cps.obs.exp.tcacgg+cps.obs.exp.tcacgt+cps.obs.exp.tcacta+cps.obs.exp.tcactc+cps.obs.exp.tcactg+cps.obs.exp.tcactt+cps.obs.exp.tcagaa+cps.obs.exp.tcagac+cps.obs.exp.tcagag+cps.obs.exp.tcagat+cps.obs.exp.tcagca+cps.obs.exp.tcagcc+cps.obs.exp.tcagcg+cps.obs.exp.tcagct+cps.obs.exp.tcagga+cps.obs.exp.tcaggc+cps.obs.exp.tcaggg+cps.obs.exp.tcaggt+cps.obs.exp.tcagta+cps.obs.exp.tcagtc+cps.obs.exp.tcagtg+cps.obs.exp.tcagtt+cps.obs.exp.tcataa+cps.obs.exp.tcatac+cps.obs.exp.tcatag+cps.obs.exp.tcatat+cps.obs.exp.tcatca+cps.obs.exp.tcatcc+cps.obs.exp.tcatcg+cps.obs.exp.tcatct+cps.obs.exp.tcatga+cps.obs.exp.tcatgc+cps.obs.exp.tcatgg+cps.obs.exp.tcatgt+cps.obs.exp.tcatta+cps.obs.exp.tcattc+cps.obs.exp.tcattg+cps.obs.exp.tcattt+cps.obs.exp.tccaaa+cps.obs.exp.tccaac+cps.obs.exp.tccaag+cps.obs.exp.tccaat+cps.obs.exp.tccaca+cps.obs.exp.tccacc+cps.obs.exp.tccacg+cps.obs.exp.tccact+cps.obs.exp.tccaga+cps.obs.exp.tccagc+cps.obs.exp.tccagg+cps.obs.exp.tccagt+cps.obs.exp.tccata+cps.obs.exp.tccatc+cps.obs.exp.tccatg+cps.obs.exp.tccatt+cps.obs.exp.tcccaa+cps.obs.exp.tcccac+cps.obs.exp.tcccag+cps.obs.exp.tcccat+cps.obs.exp.tcccca+cps.obs.exp.tccccc+cps.obs.exp.tccccg+cps.obs.exp.tcccct+cps.obs.exp.tcccga+cps.obs.exp.tcccgc+cps.obs.exp.tcccgg+cps.obs.exp.tcccgt+cps.obs.exp.tcccta+cps.obs.exp.tccctc+cps.obs.exp.tccctg+cps.obs.exp.tccctt+cps.obs.exp.tccgaa+cps.obs.exp.tccgac+cps.obs.exp.tccgag+cps.obs.exp.tccgat+cps.obs.exp.tccgca+cps.obs.exp.tccgcc+cps.obs.exp.tccgcg+cps.obs.exp.tccgct+cps.obs.exp.tccgga+cps.obs.exp.tccggc+cps.obs.exp.tccggg+cps.obs.exp.tccggt+cps.obs.exp.tccgta+cps.obs.exp.tccgtc+cps.obs.exp.tccgtg+cps.obs.exp.tccgtt+cps.obs.exp.tcctaa+cps.obs.exp.tcctac+cps.obs.exp.tcctag+cps.obs.exp.tcctat+cps.obs.exp.tcctca+cps.obs.exp.tcctcc+cps.obs.exp.tcctcg+cps.obs.exp.tcctct+cps.obs.exp.tcctga+cps.obs.exp.tcctgc+cps.obs.exp.tcctgg+cps.obs.exp.tcctgt+cps.obs.exp.tcctta+cps.obs.exp.tccttc+cps.obs.exp.tccttg+cps.obs.exp.tccttt+cps.obs.exp.tcgaaa+cps.obs.exp.tcgaac+cps.obs.exp.tcgaag+cps.obs.exp.tcgaat+cps.obs.exp.tcgaca+cps.obs.exp.tcgacc+cps.obs.exp.tcgacg+cps.obs.exp.tcgact+cps.obs.exp.tcgaga+cps.obs.exp.tcgagc+cps.obs.exp.tcgagg+cps.obs.exp.tcgagt+cps.obs.exp.tcgata+cps.obs.exp.tcgatc+cps.obs.exp.tcgatg+cps.obs.exp.tcgatt+cps.obs.exp.tcgcaa+cps.obs.exp.tcgcac+cps.obs.exp.tcgcag+cps.obs.exp.tcgcat+cps.obs.exp.tcgcca+cps.obs.exp.tcgccc+cps.obs.exp.tcgccg+cps.obs.exp.tcgcct+cps.obs.exp.tcgcga+cps.obs.exp.tcgcgc+cps.obs.exp.tcgcgg+cps.obs.exp.tcgcgt+cps.obs.exp.tcgcta+cps.obs.exp.tcgctc+cps.obs.exp.tcgctg+cps.obs.exp.tcgctt+cps.obs.exp.tcggaa+cps.obs.exp.tcggac+cps.obs.exp.tcggag+cps.obs.exp.tcggat+cps.obs.exp.tcggca+cps.obs.exp.tcggcc+cps.obs.exp.tcggcg+cps.obs.exp.tcggct+cps.obs.exp.tcggga+cps.obs.exp.tcgggc+cps.obs.exp.tcgggg+cps.obs.exp.tcgggt+cps.obs.exp.tcggta+cps.obs.exp.tcggtc+cps.obs.exp.tcggtg+cps.obs.exp.tcggtt+cps.obs.exp.tcgtaa+cps.obs.exp.tcgtac+cps.obs.exp.tcgtag+cps.obs.exp.tcgtat+cps.obs.exp.tcgtca+cps.obs.exp.tcgtcc+cps.obs.exp.tcgtcg+cps.obs.exp.tcgtct+cps.obs.exp.tcgtga+cps.obs.exp.tcgtgc+cps.obs.exp.tcgtgg+cps.obs.exp.tcgtgt+cps.obs.exp.tcgtta+cps.obs.exp.tcgttc+cps.obs.exp.tcgttg+cps.obs.exp.tcgttt+cps.obs.exp.tctaaa+cps.obs.exp.tctaac+cps.obs.exp.tctaag+cps.obs.exp.tctaat+cps.obs.exp.tctaca+cps.obs.exp.tctacc+cps.obs.exp.tctacg+cps.obs.exp.tctact+cps.obs.exp.tctaga+cps.obs.exp.tctagc+cps.obs.exp.tctagg+cps.obs.exp.tctagt+cps.obs.exp.tctata+cps.obs.exp.tctatc+cps.obs.exp.tctatg+cps.obs.exp.tctatt+cps.obs.exp.tctcaa+cps.obs.exp.tctcac+cps.obs.exp.tctcag+cps.obs.exp.tctcat+cps.obs.exp.tctcca+cps.obs.exp.tctccc+cps.obs.exp.tctccg+cps.obs.exp.tctcct+cps.obs.exp.tctcga+cps.obs.exp.tctcgc+cps.obs.exp.tctcgg+cps.obs.exp.tctcgt+cps.obs.exp.tctcta+cps.obs.exp.tctctc+cps.obs.exp.tctctg+cps.obs.exp.tctctt+cps.obs.exp.tctgaa+cps.obs.exp.tctgac+cps.obs.exp.tctgag+cps.obs.exp.tctgat+cps.obs.exp.tctgca+cps.obs.exp.tctgcc+cps.obs.exp.tctgcg+cps.obs.exp.tctgct+cps.obs.exp.tctgga+cps.obs.exp.tctggc+cps.obs.exp.tctggg+cps.obs.exp.tctggt+cps.obs.exp.tctgta+cps.obs.exp.tctgtc+cps.obs.exp.tctgtg+cps.obs.exp.tctgtt+cps.obs.exp.tcttaa+cps.obs.exp.tcttac+cps.obs.exp.tcttag+cps.obs.exp.tcttat+cps.obs.exp.tcttca+cps.obs.exp.tcttcc+cps.obs.exp.tcttcg+cps.obs.exp.tcttct+cps.obs.exp.tcttga+cps.obs.exp.tcttgc+cps.obs.exp.tcttgg+cps.obs.exp.tcttgt+cps.obs.exp.tcttta+cps.obs.exp.tctttc+cps.obs.exp.tctttg+cps.obs.exp.tctttt+cps.obs.exp.tgaaaa+cps.obs.exp.tgaaac+cps.obs.exp.tgaaag+cps.obs.exp.tgaaat+cps.obs.exp.tgaaca+cps.obs.exp.tgaacc+cps.obs.exp.tgaacg+cps.obs.exp.tgaact+cps.obs.exp.tgaaga+cps.obs.exp.tgaagc+cps.obs.exp.tgaagg+cps.obs.exp.tgaagt+cps.obs.exp.tgaata+cps.obs.exp.tgaatc+cps.obs.exp.tgaatg+cps.obs.exp.tgaatt+cps.obs.exp.tgacaa+cps.obs.exp.tgacac+cps.obs.exp.tgacag+cps.obs.exp.tgacat+cps.obs.exp.tgacca+cps.obs.exp.tgaccc+cps.obs.exp.tgaccg+cps.obs.exp.tgacct+cps.obs.exp.tgacga+cps.obs.exp.tgacgc+cps.obs.exp.tgacgg+cps.obs.exp.tgacgt+cps.obs.exp.tgacta+cps.obs.exp.tgactc+cps.obs.exp.tgactg+cps.obs.exp.tgactt+cps.obs.exp.tgagaa+cps.obs.exp.tgagac+cps.obs.exp.tgagag+cps.obs.exp.tgagat+cps.obs.exp.tgagca+cps.obs.exp.tgagcc+cps.obs.exp.tgagcg+cps.obs.exp.tgagct+cps.obs.exp.tgagga+cps.obs.exp.tgaggc+cps.obs.exp.tgaggg+cps.obs.exp.tgaggt+cps.obs.exp.tgagta+cps.obs.exp.tgagtc+cps.obs.exp.tgagtg+cps.obs.exp.tgagtt+cps.obs.exp.tgataa+cps.obs.exp.tgatac+cps.obs.exp.tgatag+cps.obs.exp.tgatat+cps.obs.exp.tgatca+cps.obs.exp.tgatcc+cps.obs.exp.tgatcg+cps.obs.exp.tgatct+cps.obs.exp.tgatga+cps.obs.exp.tgatgc+cps.obs.exp.tgatgg+cps.obs.exp.tgatgt+cps.obs.exp.tgatta+cps.obs.exp.tgattc+cps.obs.exp.tgattg+cps.obs.exp.tgattt+cps.obs.exp.tgcaaa+cps.obs.exp.tgcaac+cps.obs.exp.tgcaag+cps.obs.exp.tgcaat+cps.obs.exp.tgcaca+cps.obs.exp.tgcacc+cps.obs.exp.tgcacg+cps.obs.exp.tgcact+cps.obs.exp.tgcaga+cps.obs.exp.tgcagc+cps.obs.exp.tgcagg+cps.obs.exp.tgcagt+cps.obs.exp.tgcata+cps.obs.exp.tgcatc+cps.obs.exp.tgcatg+cps.obs.exp.tgcatt+cps.obs.exp.tgccaa+cps.obs.exp.tgccac+cps.obs.exp.tgccag+cps.obs.exp.tgccat+cps.obs.exp.tgccca+cps.obs.exp.tgcccc+cps.obs.exp.tgcccg+cps.obs.exp.tgccct+cps.obs.exp.tgccga+cps.obs.exp.tgccgc+cps.obs.exp.tgccgg+cps.obs.exp.tgccgt+cps.obs.exp.tgccta+cps.obs.exp.tgcctc+cps.obs.exp.tgcctg+cps.obs.exp.tgcctt+cps.obs.exp.tgcgaa+cps.obs.exp.tgcgac+cps.obs.exp.tgcgag+cps.obs.exp.tgcgat+cps.obs.exp.tgcgca+cps.obs.exp.tgcgcc+cps.obs.exp.tgcgcg+cps.obs.exp.tgcgct+cps.obs.exp.tgcgga+cps.obs.exp.tgcggc+cps.obs.exp.tgcggg+cps.obs.exp.tgcggt+cps.obs.exp.tgcgta+cps.obs.exp.tgcgtc+cps.obs.exp.tgcgtg+cps.obs.exp.tgcgtt+cps.obs.exp.tgctaa+cps.obs.exp.tgctac+cps.obs.exp.tgctag+cps.obs.exp.tgctat+cps.obs.exp.tgctca+cps.obs.exp.tgctcc+cps.obs.exp.tgctcg+cps.obs.exp.tgctct+cps.obs.exp.tgctga+cps.obs.exp.tgctgc+cps.obs.exp.tgctgg+cps.obs.exp.tgctgt+cps.obs.exp.tgctta+cps.obs.exp.tgcttc+cps.obs.exp.tgcttg+cps.obs.exp.tgcttt+cps.obs.exp.tggaaa+cps.obs.exp.tggaac+cps.obs.exp.tggaag+cps.obs.exp.tggaat+cps.obs.exp.tggaca+cps.obs.exp.tggacc+cps.obs.exp.tggacg+cps.obs.exp.tggact+cps.obs.exp.tggaga+cps.obs.exp.tggagc+cps.obs.exp.tggagg+cps.obs.exp.tggagt+cps.obs.exp.tggata+cps.obs.exp.tggatc+cps.obs.exp.tggatg+cps.obs.exp.tggatt+cps.obs.exp.tggcaa+cps.obs.exp.tggcac+cps.obs.exp.tggcag+cps.obs.exp.tggcat+cps.obs.exp.tggcca+cps.obs.exp.tggccc+cps.obs.exp.tggccg+cps.obs.exp.tggcct+cps.obs.exp.tggcga+cps.obs.exp.tggcgc+cps.obs.exp.tggcgg+cps.obs.exp.tggcgt+cps.obs.exp.tggcta+cps.obs.exp.tggctc+cps.obs.exp.tggctg+cps.obs.exp.tggctt+cps.obs.exp.tgggaa+cps.obs.exp.tgggac+cps.obs.exp.tgggag+cps.obs.exp.tgggat+cps.obs.exp.tgggca+cps.obs.exp.tgggcc+cps.obs.exp.tgggcg+cps.obs.exp.tgggct+cps.obs.exp.tgggga+cps.obs.exp.tggggc+cps.obs.exp.tggggg+cps.obs.exp.tggggt+cps.obs.exp.tgggta+cps.obs.exp.tgggtc+cps.obs.exp.tgggtg+cps.obs.exp.tgggtt+cps.obs.exp.tggtaa+cps.obs.exp.tggtac+cps.obs.exp.tggtag+cps.obs.exp.tggtat+cps.obs.exp.tggtca+cps.obs.exp.tggtcc+cps.obs.exp.tggtcg+cps.obs.exp.tggtct+cps.obs.exp.tggtga+cps.obs.exp.tggtgc+cps.obs.exp.tggtgg+cps.obs.exp.tggtgt+cps.obs.exp.tggtta+cps.obs.exp.tggttc+cps.obs.exp.tggttg+cps.obs.exp.tggttt+cps.obs.exp.tgtaaa+cps.obs.exp.tgtaac+cps.obs.exp.tgtaag+cps.obs.exp.tgtaat+cps.obs.exp.tgtaca+cps.obs.exp.tgtacc+cps.obs.exp.tgtacg+cps.obs.exp.tgtact+cps.obs.exp.tgtaga+cps.obs.exp.tgtagc+cps.obs.exp.tgtagg+cps.obs.exp.tgtagt+cps.obs.exp.tgtata+cps.obs.exp.tgtatc+cps.obs.exp.tgtatg+cps.obs.exp.tgtatt+cps.obs.exp.tgtcaa+cps.obs.exp.tgtcac+cps.obs.exp.tgtcag+cps.obs.exp.tgtcat+cps.obs.exp.tgtcca+cps.obs.exp.tgtccc+cps.obs.exp.tgtccg+cps.obs.exp.tgtcct+cps.obs.exp.tgtcga+cps.obs.exp.tgtcgc+cps.obs.exp.tgtcgg+cps.obs.exp.tgtcgt+cps.obs.exp.tgtcta+cps.obs.exp.tgtctc+cps.obs.exp.tgtctg+cps.obs.exp.tgtctt+cps.obs.exp.tgtgaa+cps.obs.exp.tgtgac+cps.obs.exp.tgtgag+cps.obs.exp.tgtgat+cps.obs.exp.tgtgca+cps.obs.exp.tgtgcc+cps.obs.exp.tgtgcg+cps.obs.exp.tgtgct+cps.obs.exp.tgtgga+cps.obs.exp.tgtggc+cps.obs.exp.tgtggg+cps.obs.exp.tgtggt+cps.obs.exp.tgtgta+cps.obs.exp.tgtgtc+cps.obs.exp.tgtgtg+cps.obs.exp.tgtgtt+cps.obs.exp.tgttaa+cps.obs.exp.tgttac+cps.obs.exp.tgttag+cps.obs.exp.tgttat+cps.obs.exp.tgttca+cps.obs.exp.tgttcc+cps.obs.exp.tgttcg+cps.obs.exp.tgttct+cps.obs.exp.tgttga+cps.obs.exp.tgttgc+cps.obs.exp.tgttgg+cps.obs.exp.tgttgt+cps.obs.exp.tgttta+cps.obs.exp.tgtttc+cps.obs.exp.tgtttg+cps.obs.exp.tgtttt+cps.obs.exp.ttaaaa+cps.obs.exp.ttaaac+cps.obs.exp.ttaaag+cps.obs.exp.ttaaat+cps.obs.exp.ttaaca+cps.obs.exp.ttaacc+cps.obs.exp.ttaacg+cps.obs.exp.ttaact+cps.obs.exp.ttaaga+cps.obs.exp.ttaagc+cps.obs.exp.ttaagg+cps.obs.exp.ttaagt+cps.obs.exp.ttaata+cps.obs.exp.ttaatc+cps.obs.exp.ttaatg+cps.obs.exp.ttaatt+cps.obs.exp.ttacaa+cps.obs.exp.ttacac+cps.obs.exp.ttacag+cps.obs.exp.ttacat+cps.obs.exp.ttacca+cps.obs.exp.ttaccc+cps.obs.exp.ttaccg+cps.obs.exp.ttacct+cps.obs.exp.ttacga+cps.obs.exp.ttacgc+cps.obs.exp.ttacgg+cps.obs.exp.ttacgt+cps.obs.exp.ttacta+cps.obs.exp.ttactc+cps.obs.exp.ttactg+cps.obs.exp.ttactt+cps.obs.exp.ttagaa+cps.obs.exp.ttagac+cps.obs.exp.ttagag+cps.obs.exp.ttagat+cps.obs.exp.ttagca+cps.obs.exp.ttagcc+cps.obs.exp.ttagcg+cps.obs.exp.ttagct+cps.obs.exp.ttagga+cps.obs.exp.ttaggc+cps.obs.exp.ttaggg+cps.obs.exp.ttaggt+cps.obs.exp.ttagta+cps.obs.exp.ttagtc+cps.obs.exp.ttagtg+cps.obs.exp.ttagtt+cps.obs.exp.ttataa+cps.obs.exp.ttatac+cps.obs.exp.ttatag+cps.obs.exp.ttatat+cps.obs.exp.ttatca+cps.obs.exp.ttatcc+cps.obs.exp.ttatcg+cps.obs.exp.ttatct+cps.obs.exp.ttatga+cps.obs.exp.ttatgc+cps.obs.exp.ttatgg+cps.obs.exp.ttatgt+cps.obs.exp.ttatta+cps.obs.exp.ttattc+cps.obs.exp.ttattg+cps.obs.exp.ttattt+cps.obs.exp.ttcaaa+cps.obs.exp.ttcaac+cps.obs.exp.ttcaag+cps.obs.exp.ttcaat+cps.obs.exp.ttcaca+cps.obs.exp.ttcacc+cps.obs.exp.ttcacg+cps.obs.exp.ttcact+cps.obs.exp.ttcaga+cps.obs.exp.ttcagc+cps.obs.exp.ttcagg+cps.obs.exp.ttcagt+cps.obs.exp.ttcata+cps.obs.exp.ttcatc+cps.obs.exp.ttcatg+cps.obs.exp.ttcatt+cps.obs.exp.ttccaa+cps.obs.exp.ttccac+cps.obs.exp.ttccag+cps.obs.exp.ttccat+cps.obs.exp.ttccca+cps.obs.exp.ttcccc+cps.obs.exp.ttcccg+cps.obs.exp.ttccct+cps.obs.exp.ttccga+cps.obs.exp.ttccgc+cps.obs.exp.ttccgg+cps.obs.exp.ttccgt+cps.obs.exp.ttccta+cps.obs.exp.ttcctc+cps.obs.exp.ttcctg+cps.obs.exp.ttcctt+cps.obs.exp.ttcgaa+cps.obs.exp.ttcgac+cps.obs.exp.ttcgag+cps.obs.exp.ttcgat+cps.obs.exp.ttcgca+cps.obs.exp.ttcgcc+cps.obs.exp.ttcgcg+cps.obs.exp.ttcgct+cps.obs.exp.ttcgga+cps.obs.exp.ttcggc+cps.obs.exp.ttcggg+cps.obs.exp.ttcggt+cps.obs.exp.ttcgta+cps.obs.exp.ttcgtc+cps.obs.exp.ttcgtg+cps.obs.exp.ttcgtt+cps.obs.exp.ttctaa+cps.obs.exp.ttctac+cps.obs.exp.ttctag+cps.obs.exp.ttctat+cps.obs.exp.ttctca+cps.obs.exp.ttctcc+cps.obs.exp.ttctcg+cps.obs.exp.ttctct+cps.obs.exp.ttctga+cps.obs.exp.ttctgc+cps.obs.exp.ttctgg+cps.obs.exp.ttctgt+cps.obs.exp.ttctta+cps.obs.exp.ttcttc+cps.obs.exp.ttcttg+cps.obs.exp.ttcttt+cps.obs.exp.ttgaaa+cps.obs.exp.ttgaac+cps.obs.exp.ttgaag+cps.obs.exp.ttgaat+cps.obs.exp.ttgaca+cps.obs.exp.ttgacc+cps.obs.exp.ttgacg+cps.obs.exp.ttgact+cps.obs.exp.ttgaga+cps.obs.exp.ttgagc+cps.obs.exp.ttgagg+cps.obs.exp.ttgagt+cps.obs.exp.ttgata+cps.obs.exp.ttgatc+cps.obs.exp.ttgatg+cps.obs.exp.ttgatt+cps.obs.exp.ttgcaa+cps.obs.exp.ttgcac+cps.obs.exp.ttgcag+cps.obs.exp.ttgcat+cps.obs.exp.ttgcca+cps.obs.exp.ttgccc+cps.obs.exp.ttgccg+cps.obs.exp.ttgcct+cps.obs.exp.ttgcga+cps.obs.exp.ttgcgc+cps.obs.exp.ttgcgg+cps.obs.exp.ttgcgt+cps.obs.exp.ttgcta+cps.obs.exp.ttgctc+cps.obs.exp.ttgctg+cps.obs.exp.ttgctt+cps.obs.exp.ttggaa+cps.obs.exp.ttggac+cps.obs.exp.ttggag+cps.obs.exp.ttggat+cps.obs.exp.ttggca+cps.obs.exp.ttggcc+cps.obs.exp.ttggcg+cps.obs.exp.ttggct+cps.obs.exp.ttggga+cps.obs.exp.ttgggc+cps.obs.exp.ttgggg+cps.obs.exp.ttgggt+cps.obs.exp.ttggta+cps.obs.exp.ttggtc+cps.obs.exp.ttggtg+cps.obs.exp.ttggtt+cps.obs.exp.ttgtaa+cps.obs.exp.ttgtac+cps.obs.exp.ttgtag+cps.obs.exp.ttgtat+cps.obs.exp.ttgtca+cps.obs.exp.ttgtcc+cps.obs.exp.ttgtcg+cps.obs.exp.ttgtct+cps.obs.exp.ttgtga+cps.obs.exp.ttgtgc+cps.obs.exp.ttgtgg+cps.obs.exp.ttgtgt+cps.obs.exp.ttgtta+cps.obs.exp.ttgttc+cps.obs.exp.ttgttg+cps.obs.exp.ttgttt+cps.obs.exp.tttaaa+cps.obs.exp.tttaac+cps.obs.exp.tttaag+cps.obs.exp.tttaat+cps.obs.exp.tttaca+cps.obs.exp.tttacc+cps.obs.exp.tttacg+cps.obs.exp.tttact+cps.obs.exp.tttaga+cps.obs.exp.tttagc+cps.obs.exp.tttagg+cps.obs.exp.tttagt+cps.obs.exp.tttata+cps.obs.exp.tttatc+cps.obs.exp.tttatg+cps.obs.exp.tttatt+cps.obs.exp.tttcaa+cps.obs.exp.tttcac+cps.obs.exp.tttcag+cps.obs.exp.tttcat+cps.obs.exp.tttcca+cps.obs.exp.tttccc+cps.obs.exp.tttccg+cps.obs.exp.tttcct+cps.obs.exp.tttcga+cps.obs.exp.tttcgc+cps.obs.exp.tttcgg+cps.obs.exp.tttcgt+cps.obs.exp.tttcta+cps.obs.exp.tttctc+cps.obs.exp.tttctg+cps.obs.exp.tttctt+cps.obs.exp.tttgaa+cps.obs.exp.tttgac+cps.obs.exp.tttgag+cps.obs.exp.tttgat+cps.obs.exp.tttgca+cps.obs.exp.tttgcc+cps.obs.exp.tttgcg+cps.obs.exp.tttgct+cps.obs.exp.tttgga+cps.obs.exp.tttggc+cps.obs.exp.tttggg+cps.obs.exp.tttggt+cps.obs.exp.tttgta+cps.obs.exp.tttgtc+cps.obs.exp.tttgtg+cps.obs.exp.tttgtt+cps.obs.exp.ttttaa+cps.obs.exp.ttttac+cps.obs.exp.ttttag+cps.obs.exp.ttttat+cps.obs.exp.ttttca+cps.obs.exp.ttttcc+cps.obs.exp.ttttcg+cps.obs.exp.ttttct+cps.obs.exp.ttttga+cps.obs.exp.ttttgc+cps.obs.exp.ttttgg+cps.obs.exp.ttttgt+cps.obs.exp.ttttta+cps.obs.exp.tttttc+cps.obs.exp.tttttg+cps.obs.exp.tttttt
 
#cpb.dn31<-cps.obs.exp.aaaaaa+cps.obs.exp.aaaaac+cps.obs.exp.aaaaag+cps.obs.exp.aaaaat+cps.obs.exp.aaaaca+cps.obs.exp.aaaacc+cps.obs.exp.aaaacg+cps.obs.exp.aaaact+cps.obs.exp.aaaaga+cps.obs.exp.aaaagc+cps.obs.exp.aaaagg+cps.obs.exp.aaaagt+cps.obs.exp.aaaata+cps.obs.exp.aaaatc+cps.obs.exp.aaaatg+cps.obs.exp.aaaatt+cps.obs.exp.aaacaa+cps.obs.exp.aaacac+cps.obs.exp.aaacag+cps.obs.exp.aaacat+cps.obs.exp.aaacca+cps.obs.exp.aaaccc+cps.obs.exp.aaaccg+cps.obs.exp.aaacct+cps.obs.exp.aaacga+cps.obs.exp.aaacgc+cps.obs.exp.aaacgg+cps.obs.exp.aaacgt+cps.obs.exp.aaacta+cps.obs.exp.aaactc+cps.obs.exp.aaactg+cps.obs.exp.aaactt+cps.obs.exp.aaagaa+cps.obs.exp.aaagac+cps.obs.exp.aaagag+cps.obs.exp.aaagat+cps.obs.exp.aaagca+cps.obs.exp.aaagcc+cps.obs.exp.aaagcg+cps.obs.exp.aaagct+cps.obs.exp.aaagga+cps.obs.exp.aaaggc+cps.obs.exp.aaaggg+cps.obs.exp.aaaggt+cps.obs.exp.aaagta+cps.obs.exp.aaagtc+cps.obs.exp.aaagtg+cps.obs.exp.aaagtt+cps.obs.exp.aaatac+cps.obs.exp.aaatat+cps.obs.exp.aaatca+cps.obs.exp.aaatcc+cps.obs.exp.aaatcg+cps.obs.exp.aaatct+cps.obs.exp.aaatgc+cps.obs.exp.aaatgg+cps.obs.exp.aaatgt+cps.obs.exp.aaatta+cps.obs.exp.aaattc+cps.obs.exp.aaattg+cps.obs.exp.aaattt+cps.obs.exp.aacaaa+cps.obs.exp.aacaac+cps.obs.exp.aacaag+cps.obs.exp.aacaat+cps.obs.exp.aacaca+cps.obs.exp.aacacc+cps.obs.exp.aacacg+cps.obs.exp.aacact+cps.obs.exp.aacaga+cps.obs.exp.aacagc+cps.obs.exp.aacagg+cps.obs.exp.aacagt+cps.obs.exp.aacata+cps.obs.exp.aacatc+cps.obs.exp.aacatg+cps.obs.exp.aacatt+cps.obs.exp.aaccaa+cps.obs.exp.aaccac+cps.obs.exp.aaccag+cps.obs.exp.aaccat+cps.obs.exp.aaccca+cps.obs.exp.aacccc+cps.obs.exp.aacccg+cps.obs.exp.aaccct+cps.obs.exp.aaccga+cps.obs.exp.aaccgc+cps.obs.exp.aaccgg+cps.obs.exp.aaccgt+cps.obs.exp.aaccta+cps.obs.exp.aacctc+cps.obs.exp.aacctg+cps.obs.exp.aacctt+cps.obs.exp.aacgaa+cps.obs.exp.aacgac+cps.obs.exp.aacgag+cps.obs.exp.aacgat+cps.obs.exp.aacgca+cps.obs.exp.aacgcc+cps.obs.exp.aacgcg+cps.obs.exp.aacgct+cps.obs.exp.aacgga+cps.obs.exp.aacggc+cps.obs.exp.aacggg+cps.obs.exp.aacggt+cps.obs.exp.aacgta+cps.obs.exp.aacgtc+cps.obs.exp.aacgtg+cps.obs.exp.aacgtt+cps.obs.exp.aactac+cps.obs.exp.aactat+cps.obs.exp.aactca+cps.obs.exp.aactcc+cps.obs.exp.aactcg+cps.obs.exp.aactct+cps.obs.exp.aactgc+cps.obs.exp.aactgg+cps.obs.exp.aactgt+cps.obs.exp.aactta+cps.obs.exp.aacttc+cps.obs.exp.aacttg+cps.obs.exp.aacttt+cps.obs.exp.aagaaa+cps.obs.exp.aagaac+cps.obs.exp.aagaag+cps.obs.exp.aagaat+cps.obs.exp.aagaca+cps.obs.exp.aagacc+cps.obs.exp.aagacg+cps.obs.exp.aagact+cps.obs.exp.aagaga+cps.obs.exp.aagagc+cps.obs.exp.aagagg+cps.obs.exp.aagagt+cps.obs.exp.aagata+cps.obs.exp.aagatc+cps.obs.exp.aagatg+cps.obs.exp.aagatt+cps.obs.exp.aagcaa+cps.obs.exp.aagcac+cps.obs.exp.aagcag+cps.obs.exp.aagcat+cps.obs.exp.aagcca+cps.obs.exp.aagccc+cps.obs.exp.aagccg+cps.obs.exp.aagcct+cps.obs.exp.aagcga+cps.obs.exp.aagcgc+cps.obs.exp.aagcgg+cps.obs.exp.aagcgt+cps.obs.exp.aagcta+cps.obs.exp.aagctc+cps.obs.exp.aagctg+cps.obs.exp.aagctt+cps.obs.exp.aaggaa+cps.obs.exp.aaggac+cps.obs.exp.aaggag+cps.obs.exp.aaggat+cps.obs.exp.aaggca+cps.obs.exp.aaggcc+cps.obs.exp.aaggcg+cps.obs.exp.aaggct+cps.obs.exp.aaggga+cps.obs.exp.aagggc+cps.obs.exp.aagggg+cps.obs.exp.aagggt+cps.obs.exp.aaggta+cps.obs.exp.aaggtc+cps.obs.exp.aaggtg+cps.obs.exp.aaggtt+cps.obs.exp.aagtac+cps.obs.exp.aagtat+cps.obs.exp.aagtca+cps.obs.exp.aagtcc+cps.obs.exp.aagtcg+cps.obs.exp.aagtct+cps.obs.exp.aagtgc+cps.obs.exp.aagtgg+cps.obs.exp.aagtgt+cps.obs.exp.aagtta+cps.obs.exp.aagttc+cps.obs.exp.aagttg+cps.obs.exp.aagttt+cps.obs.exp.aataaa+cps.obs.exp.aataac+cps.obs.exp.aataag+cps.obs.exp.aataat+cps.obs.exp.aataca+cps.obs.exp.aatacc+cps.obs.exp.aatacg+cps.obs.exp.aatact+cps.obs.exp.aataga+cps.obs.exp.aatagc+cps.obs.exp.aatagg+cps.obs.exp.aatagt+cps.obs.exp.aatata+cps.obs.exp.aatatc+cps.obs.exp.aatatg+cps.obs.exp.aatatt+cps.obs.exp.aatcaa+cps.obs.exp.aatcac+cps.obs.exp.aatcag+cps.obs.exp.aatcat+cps.obs.exp.aatcca+cps.obs.exp.aatccc+cps.obs.exp.aatccg+cps.obs.exp.aatcct+cps.obs.exp.aatcga+cps.obs.exp.aatcgc+cps.obs.exp.aatcgg+cps.obs.exp.aatcgt+cps.obs.exp.aatcta+cps.obs.exp.aatctc+cps.obs.exp.aatctg+cps.obs.exp.aatctt+cps.obs.exp.aatgaa+cps.obs.exp.aatgac+cps.obs.exp.aatgag+cps.obs.exp.aatgat+cps.obs.exp.aatgca+cps.obs.exp.aatgcc+cps.obs.exp.aatgcg+cps.obs.exp.aatgct+cps.obs.exp.aatgga+cps.obs.exp.aatggc+cps.obs.exp.aatggg+cps.obs.exp.aatggt+cps.obs.exp.aatgta+cps.obs.exp.aatgtc+cps.obs.exp.aatgtg+cps.obs.exp.aatgtt+cps.obs.exp.aattac+cps.obs.exp.aattat+cps.obs.exp.aattca+cps.obs.exp.aattcc+cps.obs.exp.aattcg+cps.obs.exp.aattct+cps.obs.exp.aattgc+cps.obs.exp.aattgg+cps.obs.exp.aattgt+cps.obs.exp.aattta+cps.obs.exp.aatttc+cps.obs.exp.aatttg+cps.obs.exp.aatttt+cps.obs.exp.acaaaa+cps.obs.exp.acaaac+cps.obs.exp.acaaag+cps.obs.exp.acaaat+cps.obs.exp.acaaca+cps.obs.exp.acaacc+cps.obs.exp.acaacg+cps.obs.exp.acaact+cps.obs.exp.acaaga+cps.obs.exp.acaagc+cps.obs.exp.acaagg+cps.obs.exp.acaagt+cps.obs.exp.acaata+cps.obs.exp.acaatc+cps.obs.exp.acaatg+cps.obs.exp.acaatt+cps.obs.exp.acacaa+cps.obs.exp.acacac+cps.obs.exp.acacag+cps.obs.exp.acacat+cps.obs.exp.acacca+cps.obs.exp.acaccc+cps.obs.exp.acaccg+cps.obs.exp.acacct+cps.obs.exp.acacga+cps.obs.exp.acacgc+cps.obs.exp.acacgg+cps.obs.exp.acacgt+cps.obs.exp.acacta+cps.obs.exp.acactc+cps.obs.exp.acactg+cps.obs.exp.acactt+cps.obs.exp.acagaa+cps.obs.exp.acagac+cps.obs.exp.acagag+cps.obs.exp.acagat+cps.obs.exp.acagca+cps.obs.exp.acagcc+cps.obs.exp.acagcg+cps.obs.exp.acagct+cps.obs.exp.acagga+cps.obs.exp.acaggc+cps.obs.exp.acaggg+cps.obs.exp.acaggt+cps.obs.exp.acagta+cps.obs.exp.acagtc+cps.obs.exp.acagtg+cps.obs.exp.acagtt+cps.obs.exp.acatac+cps.obs.exp.acatat+cps.obs.exp.acatca+cps.obs.exp.acatcc+cps.obs.exp.acatcg+cps.obs.exp.acatct+cps.obs.exp.acatgc+cps.obs.exp.acatgg+cps.obs.exp.acatgt+cps.obs.exp.acatta+cps.obs.exp.acattc+cps.obs.exp.acattg+cps.obs.exp.acattt+cps.obs.exp.accaaa+cps.obs.exp.accaac+cps.obs.exp.accaag+cps.obs.exp.accaat+cps.obs.exp.accaca+cps.obs.exp.accacc+cps.obs.exp.accacg+cps.obs.exp.accact+cps.obs.exp.accaga+cps.obs.exp.accagc+cps.obs.exp.accagg+cps.obs.exp.accagt+cps.obs.exp.accata+cps.obs.exp.accatc+cps.obs.exp.accatg+cps.obs.exp.accatt+cps.obs.exp.acccaa+cps.obs.exp.acccac+cps.obs.exp.acccag+cps.obs.exp.acccat+cps.obs.exp.acccca+cps.obs.exp.accccc+cps.obs.exp.accccg+cps.obs.exp.acccct+cps.obs.exp.acccga+cps.obs.exp.acccgc+cps.obs.exp.acccgg+cps.obs.exp.acccgt+cps.obs.exp.acccta+cps.obs.exp.accctc+cps.obs.exp.accctg+cps.obs.exp.accctt+cps.obs.exp.accgaa+cps.obs.exp.accgac+cps.obs.exp.accgag+cps.obs.exp.accgat+cps.obs.exp.accgca+cps.obs.exp.accgcc+cps.obs.exp.accgcg+cps.obs.exp.accgct+cps.obs.exp.accgga+cps.obs.exp.accggc+cps.obs.exp.accggg+cps.obs.exp.accggt+cps.obs.exp.accgta+cps.obs.exp.accgtc+cps.obs.exp.accgtg+cps.obs.exp.accgtt+cps.obs.exp.acctac+cps.obs.exp.acctat+cps.obs.exp.acctca+cps.obs.exp.acctcc+cps.obs.exp.acctcg+cps.obs.exp.acctct+cps.obs.exp.acctgc+cps.obs.exp.acctgg+cps.obs.exp.acctgt+cps.obs.exp.acctta+cps.obs.exp.accttc+cps.obs.exp.accttg+cps.obs.exp.accttt+cps.obs.exp.acgaaa+cps.obs.exp.acgaac+cps.obs.exp.acgaag+cps.obs.exp.acgaat+cps.obs.exp.acgaca+cps.obs.exp.acgacc+cps.obs.exp.acgacg+cps.obs.exp.acgact+cps.obs.exp.acgaga+cps.obs.exp.acgagc+cps.obs.exp.acgagg+cps.obs.exp.acgagt+cps.obs.exp.acgata+cps.obs.exp.acgatc+cps.obs.exp.acgatg+cps.obs.exp.acgatt+cps.obs.exp.acgcaa+cps.obs.exp.acgcac+cps.obs.exp.acgcag+cps.obs.exp.acgcat+cps.obs.exp.acgcca+cps.obs.exp.acgccc+cps.obs.exp.acgccg+cps.obs.exp.acgcct+cps.obs.exp.acgcga+cps.obs.exp.acgcgc+cps.obs.exp.acgcgg+cps.obs.exp.acgcgt+cps.obs.exp.acgcta+cps.obs.exp.acgctc+cps.obs.exp.acgctg+cps.obs.exp.acgctt+cps.obs.exp.acggaa+cps.obs.exp.acggac+cps.obs.exp.acggag+cps.obs.exp.acggat+cps.obs.exp.acggca+cps.obs.exp.acggcc+cps.obs.exp.acggcg+cps.obs.exp.acggct+cps.obs.exp.acggga+cps.obs.exp.acgggc+cps.obs.exp.acgggg+cps.obs.exp.acgggt+cps.obs.exp.acggta+cps.obs.exp.acggtc+cps.obs.exp.acggtg+cps.obs.exp.acggtt+cps.obs.exp.acgtac+cps.obs.exp.acgtat+cps.obs.exp.acgtca+cps.obs.exp.acgtcc+cps.obs.exp.acgtcg+cps.obs.exp.acgtct+cps.obs.exp.acgtgc+cps.obs.exp.acgtgg+cps.obs.exp.acgtgt+cps.obs.exp.acgtta+cps.obs.exp.acgttc+cps.obs.exp.acgttg+cps.obs.exp.acgttt+cps.obs.exp.actaaa+cps.obs.exp.actaac+cps.obs.exp.actaag+cps.obs.exp.actaat+cps.obs.exp.actaca+cps.obs.exp.actacc+cps.obs.exp.actacg+cps.obs.exp.actact+cps.obs.exp.actaga+cps.obs.exp.actagc+cps.obs.exp.actagg+cps.obs.exp.actagt+cps.obs.exp.actata+cps.obs.exp.actatc+cps.obs.exp.actatg+cps.obs.exp.actatt+cps.obs.exp.actcaa+cps.obs.exp.actcac+cps.obs.exp.actcag+cps.obs.exp.actcat+cps.obs.exp.actcca+cps.obs.exp.actccc+cps.obs.exp.actccg+cps.obs.exp.actcct+cps.obs.exp.actcga+cps.obs.exp.actcgc+cps.obs.exp.actcgg+cps.obs.exp.actcgt+cps.obs.exp.actcta+cps.obs.exp.actctc+cps.obs.exp.actctg+cps.obs.exp.actctt+cps.obs.exp.actgaa+cps.obs.exp.actgac+cps.obs.exp.actgag+cps.obs.exp.actgat+cps.obs.exp.actgca+cps.obs.exp.actgcc+cps.obs.exp.actgcg+cps.obs.exp.actgct+cps.obs.exp.actgga+cps.obs.exp.actggc+cps.obs.exp.actggg+cps.obs.exp.actggt+cps.obs.exp.actgta+cps.obs.exp.actgtc+cps.obs.exp.actgtg+cps.obs.exp.actgtt+cps.obs.exp.acttac+cps.obs.exp.acttat+cps.obs.exp.acttca+cps.obs.exp.acttcc+cps.obs.exp.acttcg+cps.obs.exp.acttct+cps.obs.exp.acttgc+cps.obs.exp.acttgg+cps.obs.exp.acttgt+cps.obs.exp.acttta+cps.obs.exp.actttc+cps.obs.exp.actttg+cps.obs.exp.actttt+cps.obs.exp.agaaaa+cps.obs.exp.agaaac+cps.obs.exp.agaaag+cps.obs.exp.agaaat+cps.obs.exp.agaaca+cps.obs.exp.agaacc+cps.obs.exp.agaacg+cps.obs.exp.agaact+cps.obs.exp.agaaga+cps.obs.exp.agaagc+cps.obs.exp.agaagg+cps.obs.exp.agaagt+cps.obs.exp.agaata+cps.obs.exp.agaatc+cps.obs.exp.agaatg+cps.obs.exp.agaatt+cps.obs.exp.agacaa+cps.obs.exp.agacac+cps.obs.exp.agacag+cps.obs.exp.agacat+cps.obs.exp.agacca+cps.obs.exp.agaccc+cps.obs.exp.agaccg+cps.obs.exp.agacct+cps.obs.exp.agacga+cps.obs.exp.agacgc+cps.obs.exp.agacgg+cps.obs.exp.agacgt+cps.obs.exp.agacta+cps.obs.exp.agactc+cps.obs.exp.agactg+cps.obs.exp.agactt+cps.obs.exp.agagaa+cps.obs.exp.agagac+cps.obs.exp.agagag+cps.obs.exp.agagat+cps.obs.exp.agagca+cps.obs.exp.agagcc+cps.obs.exp.agagcg+cps.obs.exp.agagct+cps.obs.exp.agagga+cps.obs.exp.agaggc+cps.obs.exp.agaggg+cps.obs.exp.agaggt+cps.obs.exp.agagta+cps.obs.exp.agagtc+cps.obs.exp.agagtg+cps.obs.exp.agagtt+cps.obs.exp.agatac+cps.obs.exp.agatat+cps.obs.exp.agatca+cps.obs.exp.agatcc+cps.obs.exp.agatcg+cps.obs.exp.agatct+cps.obs.exp.agatgc+cps.obs.exp.agatgg+cps.obs.exp.agatgt+cps.obs.exp.agatta+cps.obs.exp.agattc+cps.obs.exp.agattg+cps.obs.exp.agattt+cps.obs.exp.agcaaa+cps.obs.exp.agcaac+cps.obs.exp.agcaag+cps.obs.exp.agcaat+cps.obs.exp.agcaca+cps.obs.exp.agcacc+cps.obs.exp.agcacg+cps.obs.exp.agcact+cps.obs.exp.agcaga+cps.obs.exp.agcagc+cps.obs.exp.agcagg+cps.obs.exp.agcagt+cps.obs.exp.agcata+cps.obs.exp.agcatc+cps.obs.exp.agcatg+cps.obs.exp.agcatt+cps.obs.exp.agccaa+cps.obs.exp.agccac+cps.obs.exp.agccag+cps.obs.exp.agccat+cps.obs.exp.agccca+cps.obs.exp.agcccc+cps.obs.exp.agcccg+cps.obs.exp.agccct+cps.obs.exp.agccga+cps.obs.exp.agccgc+cps.obs.exp.agccgg+cps.obs.exp.agccgt+cps.obs.exp.agccta+cps.obs.exp.agcctc+cps.obs.exp.agcctg+cps.obs.exp.agcctt+cps.obs.exp.agcgaa+cps.obs.exp.agcgac+cps.obs.exp.agcgag+cps.obs.exp.agcgat+cps.obs.exp.agcgca+cps.obs.exp.agcgcc+cps.obs.exp.agcgcg+cps.obs.exp.agcgct+cps.obs.exp.agcgga+cps.obs.exp.agcggc+cps.obs.exp.agcggg+cps.obs.exp.agcggt+cps.obs.exp.agcgta+cps.obs.exp.agcgtc+cps.obs.exp.agcgtg+cps.obs.exp.agcgtt+cps.obs.exp.agctac+cps.obs.exp.agctat+cps.obs.exp.agctca+cps.obs.exp.agctcc+cps.obs.exp.agctcg+cps.obs.exp.agctct+cps.obs.exp.agctgc+cps.obs.exp.agctgg+cps.obs.exp.agctgt+cps.obs.exp.agctta+cps.obs.exp.agcttc+cps.obs.exp.agcttg+cps.obs.exp.agcttt+cps.obs.exp.aggaaa+cps.obs.exp.aggaac+cps.obs.exp.aggaag+cps.obs.exp.aggaat+cps.obs.exp.aggaca+cps.obs.exp.aggacc+cps.obs.exp.aggacg+cps.obs.exp.aggact+cps.obs.exp.aggaga+cps.obs.exp.aggagc+cps.obs.exp.aggagg+cps.obs.exp.aggagt+cps.obs.exp.aggata+cps.obs.exp.aggatc+cps.obs.exp.aggatg+cps.obs.exp.aggatt+cps.obs.exp.aggcaa+cps.obs.exp.aggcac+cps.obs.exp.aggcag+cps.obs.exp.aggcat+cps.obs.exp.aggcca+cps.obs.exp.aggccc+cps.obs.exp.aggccg+cps.obs.exp.aggcct+cps.obs.exp.aggcga+cps.obs.exp.aggcgc+cps.obs.exp.aggcgg+cps.obs.exp.aggcgt+cps.obs.exp.aggcta+cps.obs.exp.aggctc+cps.obs.exp.aggctg+cps.obs.exp.aggctt+cps.obs.exp.agggaa+cps.obs.exp.agggac+cps.obs.exp.agggag+cps.obs.exp.agggat+cps.obs.exp.agggca+cps.obs.exp.agggcc+cps.obs.exp.agggcg+cps.obs.exp.agggct+cps.obs.exp.agggga+cps.obs.exp.aggggc+cps.obs.exp.aggggg+cps.obs.exp.aggggt+cps.obs.exp.agggta+cps.obs.exp.agggtc+cps.obs.exp.agggtg+cps.obs.exp.agggtt+cps.obs.exp.aggtac+cps.obs.exp.aggtat+cps.obs.exp.aggtca+cps.obs.exp.aggtcc+cps.obs.exp.aggtcg+cps.obs.exp.aggtct+cps.obs.exp.aggtgc+cps.obs.exp.aggtgg+cps.obs.exp.aggtgt+cps.obs.exp.aggtta+cps.obs.exp.aggttc+cps.obs.exp.aggttg+cps.obs.exp.aggttt+cps.obs.exp.agtaaa+cps.obs.exp.agtaac+cps.obs.exp.agtaag+cps.obs.exp.agtaat+cps.obs.exp.agtaca+cps.obs.exp.agtacc+cps.obs.exp.agtacg+cps.obs.exp.agtact+cps.obs.exp.agtaga+cps.obs.exp.agtagc+cps.obs.exp.agtagg+cps.obs.exp.agtagt+cps.obs.exp.agtata+cps.obs.exp.agtatc+cps.obs.exp.agtatg+cps.obs.exp.agtatt+cps.obs.exp.agtcaa+cps.obs.exp.agtcac+cps.obs.exp.agtcag+cps.obs.exp.agtcat+cps.obs.exp.agtcca+cps.obs.exp.agtccc+cps.obs.exp.agtccg+cps.obs.exp.agtcct+cps.obs.exp.agtcga+cps.obs.exp.agtcgc+cps.obs.exp.agtcgg+cps.obs.exp.agtcgt+cps.obs.exp.agtcta+cps.obs.exp.agtctc+cps.obs.exp.agtctg+cps.obs.exp.agtctt+cps.obs.exp.agtgaa+cps.obs.exp.agtgac+cps.obs.exp.agtgag+cps.obs.exp.agtgat+cps.obs.exp.agtgca+cps.obs.exp.agtgcc+cps.obs.exp.agtgcg+cps.obs.exp.agtgct+cps.obs.exp.agtgga+cps.obs.exp.agtggc+cps.obs.exp.agtggg+cps.obs.exp.agtggt+cps.obs.exp.agtgta+cps.obs.exp.agtgtc+cps.obs.exp.agtgtg+cps.obs.exp.agtgtt+cps.obs.exp.agttac+cps.obs.exp.agttat+cps.obs.exp.agttca+cps.obs.exp.agttcc+cps.obs.exp.agttcg+cps.obs.exp.agttct+cps.obs.exp.agttgc+cps.obs.exp.agttgg+cps.obs.exp.agttgt+cps.obs.exp.agttta+cps.obs.exp.agtttc+cps.obs.exp.agtttg+cps.obs.exp.agtttt+cps.obs.exp.ataaaa+cps.obs.exp.ataaac+cps.obs.exp.ataaag+cps.obs.exp.ataaat+cps.obs.exp.ataaca+cps.obs.exp.ataacc+cps.obs.exp.ataacg+cps.obs.exp.ataact+cps.obs.exp.ataaga+cps.obs.exp.ataagc+cps.obs.exp.ataagg+cps.obs.exp.ataagt+cps.obs.exp.ataata+cps.obs.exp.ataatc+cps.obs.exp.ataatg+cps.obs.exp.ataatt+cps.obs.exp.atacaa+cps.obs.exp.atacac+cps.obs.exp.atacag+cps.obs.exp.atacat+cps.obs.exp.atacca+cps.obs.exp.ataccc+cps.obs.exp.ataccg+cps.obs.exp.atacct+cps.obs.exp.atacga+cps.obs.exp.atacgc+cps.obs.exp.atacgg+cps.obs.exp.atacgt+cps.obs.exp.atacta+cps.obs.exp.atactc+cps.obs.exp.atactg+cps.obs.exp.atactt+cps.obs.exp.atagaa+cps.obs.exp.atagac+cps.obs.exp.atagag+cps.obs.exp.atagat+cps.obs.exp.atagca+cps.obs.exp.atagcc+cps.obs.exp.atagcg+cps.obs.exp.atagct+cps.obs.exp.atagga+cps.obs.exp.ataggc+cps.obs.exp.ataggg+cps.obs.exp.ataggt+cps.obs.exp.atagta+cps.obs.exp.atagtc+cps.obs.exp.atagtg+cps.obs.exp.atagtt+cps.obs.exp.atatac+cps.obs.exp.atatat+cps.obs.exp.atatca+cps.obs.exp.atatcc+cps.obs.exp.atatcg+cps.obs.exp.atatct+cps.obs.exp.atatgc+cps.obs.exp.atatgg+cps.obs.exp.atatgt+cps.obs.exp.atatta+cps.obs.exp.atattc+cps.obs.exp.atattg+cps.obs.exp.atattt+cps.obs.exp.atcaaa+cps.obs.exp.atcaac+cps.obs.exp.atcaag+cps.obs.exp.atcaat+cps.obs.exp.atcaca+cps.obs.exp.atcacc+cps.obs.exp.atcacg+cps.obs.exp.atcact+cps.obs.exp.atcaga+cps.obs.exp.atcagc+cps.obs.exp.atcagg+cps.obs.exp.atcagt+cps.obs.exp.atcata+cps.obs.exp.atcatc+cps.obs.exp.atcatg+cps.obs.exp.atcatt+cps.obs.exp.atccaa+cps.obs.exp.atccac+cps.obs.exp.atccag+cps.obs.exp.atccat+cps.obs.exp.atccca+cps.obs.exp.atcccc+cps.obs.exp.atcccg+cps.obs.exp.atccct+cps.obs.exp.atccga+cps.obs.exp.atccgc+cps.obs.exp.atccgg+cps.obs.exp.atccgt+cps.obs.exp.atccta+cps.obs.exp.atcctc+cps.obs.exp.atcctg+cps.obs.exp.atcctt+cps.obs.exp.atcgaa+cps.obs.exp.atcgac+cps.obs.exp.atcgag+cps.obs.exp.atcgat+cps.obs.exp.atcgca+cps.obs.exp.atcgcc+cps.obs.exp.atcgcg+cps.obs.exp.atcgct+cps.obs.exp.atcgga+cps.obs.exp.atcggc+cps.obs.exp.atcggg+cps.obs.exp.atcggt+cps.obs.exp.atcgta+cps.obs.exp.atcgtc+cps.obs.exp.atcgtg+cps.obs.exp.atcgtt+cps.obs.exp.atctac+cps.obs.exp.atctat+cps.obs.exp.atctca+cps.obs.exp.atctcc+cps.obs.exp.atctcg+cps.obs.exp.atctct+cps.obs.exp.atctgc+cps.obs.exp.atctgg+cps.obs.exp.atctgt+cps.obs.exp.atctta+cps.obs.exp.atcttc+cps.obs.exp.atcttg+cps.obs.exp.atcttt+cps.obs.exp.atgaaa+cps.obs.exp.atgaac+cps.obs.exp.atgaag+cps.obs.exp.atgaat+cps.obs.exp.atgaca+cps.obs.exp.atgacc+cps.obs.exp.atgacg+cps.obs.exp.atgact+cps.obs.exp.atgaga+cps.obs.exp.atgagc+cps.obs.exp.atgagg+cps.obs.exp.atgagt+cps.obs.exp.atgata+cps.obs.exp.atgatc+cps.obs.exp.atgatg+cps.obs.exp.atgatt+cps.obs.exp.atgcaa+cps.obs.exp.atgcac+cps.obs.exp.atgcag+cps.obs.exp.atgcat+cps.obs.exp.atgcca+cps.obs.exp.atgccc+cps.obs.exp.atgccg+cps.obs.exp.atgcct+cps.obs.exp.atgcga+cps.obs.exp.atgcgc+cps.obs.exp.atgcgg+cps.obs.exp.atgcgt+cps.obs.exp.atgcta+cps.obs.exp.atgctc+cps.obs.exp.atgctg+cps.obs.exp.atgctt+cps.obs.exp.atggaa+cps.obs.exp.atggac+cps.obs.exp.atggag+cps.obs.exp.atggat+cps.obs.exp.atggca+cps.obs.exp.atggcc+cps.obs.exp.atggcg+cps.obs.exp.atggct+cps.obs.exp.atggga+cps.obs.exp.atgggc+cps.obs.exp.atgggg+cps.obs.exp.atgggt+cps.obs.exp.atggta+cps.obs.exp.atggtc+cps.obs.exp.atggtg+cps.obs.exp.atggtt+cps.obs.exp.atgtac+cps.obs.exp.atgtat+cps.obs.exp.atgtca+cps.obs.exp.atgtcc+cps.obs.exp.atgtcg+cps.obs.exp.atgtct+cps.obs.exp.atgtgc+cps.obs.exp.atgtgg+cps.obs.exp.atgtgt+cps.obs.exp.atgtta+cps.obs.exp.atgttc+cps.obs.exp.atgttg+cps.obs.exp.atgttt+cps.obs.exp.attaaa+cps.obs.exp.attaac+cps.obs.exp.attaag+cps.obs.exp.attaat+cps.obs.exp.attaca+cps.obs.exp.attacc+cps.obs.exp.attacg+cps.obs.exp.attact+cps.obs.exp.attaga+cps.obs.exp.attagc+cps.obs.exp.attagg+cps.obs.exp.attagt+cps.obs.exp.attata+cps.obs.exp.attatc+cps.obs.exp.attatg+cps.obs.exp.attatt+cps.obs.exp.attcaa+cps.obs.exp.attcac+cps.obs.exp.attcag+cps.obs.exp.attcat+cps.obs.exp.attcca+cps.obs.exp.attccc+cps.obs.exp.attccg+cps.obs.exp.attcct+cps.obs.exp.attcga+cps.obs.exp.attcgc+cps.obs.exp.attcgg+cps.obs.exp.attcgt+cps.obs.exp.attcta+cps.obs.exp.attctc+cps.obs.exp.attctg+cps.obs.exp.attctt+cps.obs.exp.attgaa+cps.obs.exp.attgac+cps.obs.exp.attgag+cps.obs.exp.attgat+cps.obs.exp.attgca+cps.obs.exp.attgcc+cps.obs.exp.attgcg+cps.obs.exp.attgct+cps.obs.exp.attgga+cps.obs.exp.attggc+cps.obs.exp.attggg+cps.obs.exp.attggt+cps.obs.exp.attgta+cps.obs.exp.attgtc+cps.obs.exp.attgtg+cps.obs.exp.attgtt+cps.obs.exp.atttac+cps.obs.exp.atttat+cps.obs.exp.atttca+cps.obs.exp.atttcc+cps.obs.exp.atttcg+cps.obs.exp.atttct+cps.obs.exp.atttgc+cps.obs.exp.atttgg+cps.obs.exp.atttgt+cps.obs.exp.atttta+cps.obs.exp.attttc+cps.obs.exp.attttg+cps.obs.exp.attttt+cps.obs.exp.caaaaa+cps.obs.exp.caaaac+cps.obs.exp.caaaag+cps.obs.exp.caaaat+cps.obs.exp.caaaca+cps.obs.exp.caaacc+cps.obs.exp.caaacg+cps.obs.exp.caaact+cps.obs.exp.caaaga+cps.obs.exp.caaagc+cps.obs.exp.caaagg+cps.obs.exp.caaagt+cps.obs.exp.caaata+cps.obs.exp.caaatc+cps.obs.exp.caaatg+cps.obs.exp.caaatt+cps.obs.exp.caacaa+cps.obs.exp.caacac+cps.obs.exp.caacag+cps.obs.exp.caacat+cps.obs.exp.caacca+cps.obs.exp.caaccc+cps.obs.exp.caaccg+cps.obs.exp.caacct+cps.obs.exp.caacga+cps.obs.exp.caacgc+cps.obs.exp.caacgg+cps.obs.exp.caacgt+cps.obs.exp.caacta+cps.obs.exp.caactc+cps.obs.exp.caactg+cps.obs.exp.caactt+cps.obs.exp.caagaa+cps.obs.exp.caagac+cps.obs.exp.caagag+cps.obs.exp.caagat+cps.obs.exp.caagca+cps.obs.exp.caagcc+cps.obs.exp.caagcg+cps.obs.exp.caagct+cps.obs.exp.caagga+cps.obs.exp.caaggc+cps.obs.exp.caaggg+cps.obs.exp.caaggt+cps.obs.exp.caagta+cps.obs.exp.caagtc+cps.obs.exp.caagtg+cps.obs.exp.caagtt+cps.obs.exp.caatac+cps.obs.exp.caatat+cps.obs.exp.caatca+cps.obs.exp.caatcc+cps.obs.exp.caatcg+cps.obs.exp.caatct+cps.obs.exp.caatgc+cps.obs.exp.caatgg+cps.obs.exp.caatgt+cps.obs.exp.caatta+cps.obs.exp.caattc+cps.obs.exp.caattg+cps.obs.exp.caattt+cps.obs.exp.cacaaa+cps.obs.exp.cacaac+cps.obs.exp.cacaag+cps.obs.exp.cacaat+cps.obs.exp.cacaca+cps.obs.exp.cacacc+cps.obs.exp.cacacg+cps.obs.exp.cacact+cps.obs.exp.cacaga+cps.obs.exp.cacagc+cps.obs.exp.cacagg+cps.obs.exp.cacagt+cps.obs.exp.cacata+cps.obs.exp.cacatc+cps.obs.exp.cacatg+cps.obs.exp.cacatt+cps.obs.exp.caccaa+cps.obs.exp.caccac+cps.obs.exp.caccag+cps.obs.exp.caccat+cps.obs.exp.caccca+cps.obs.exp.cacccc+cps.obs.exp.cacccg+cps.obs.exp.caccct+cps.obs.exp.caccga+cps.obs.exp.caccgc+cps.obs.exp.caccgg+cps.obs.exp.caccgt+cps.obs.exp.caccta+cps.obs.exp.cacctc+cps.obs.exp.cacctg+cps.obs.exp.cacctt+cps.obs.exp.cacgaa+cps.obs.exp.cacgac+cps.obs.exp.cacgag+cps.obs.exp.cacgat+cps.obs.exp.cacgca+cps.obs.exp.cacgcc+cps.obs.exp.cacgcg+cps.obs.exp.cacgct+cps.obs.exp.cacgga+cps.obs.exp.cacggc+cps.obs.exp.cacggg+cps.obs.exp.cacggt+cps.obs.exp.cacgta+cps.obs.exp.cacgtc+cps.obs.exp.cacgtg+cps.obs.exp.cacgtt+cps.obs.exp.cactac+cps.obs.exp.cactat+cps.obs.exp.cactca+cps.obs.exp.cactcc+cps.obs.exp.cactcg+cps.obs.exp.cactct+cps.obs.exp.cactgc+cps.obs.exp.cactgg+cps.obs.exp.cactgt+cps.obs.exp.cactta+cps.obs.exp.cacttc+cps.obs.exp.cacttg+cps.obs.exp.cacttt+cps.obs.exp.cagaaa+cps.obs.exp.cagaac+cps.obs.exp.cagaag+cps.obs.exp.cagaat+cps.obs.exp.cagaca+cps.obs.exp.cagacc+cps.obs.exp.cagacg+cps.obs.exp.cagact+cps.obs.exp.cagaga+cps.obs.exp.cagagc+cps.obs.exp.cagagg+cps.obs.exp.cagagt+cps.obs.exp.cagata+cps.obs.exp.cagatc+cps.obs.exp.cagatg+cps.obs.exp.cagatt+cps.obs.exp.cagcaa+cps.obs.exp.cagcac+cps.obs.exp.cagcag+cps.obs.exp.cagcat+cps.obs.exp.cagcca+cps.obs.exp.cagccc+cps.obs.exp.cagccg+cps.obs.exp.cagcct+cps.obs.exp.cagcga+cps.obs.exp.cagcgc+cps.obs.exp.cagcgg+cps.obs.exp.cagcgt+cps.obs.exp.cagcta+cps.obs.exp.cagctc+cps.obs.exp.cagctg+cps.obs.exp.cagctt+cps.obs.exp.caggaa+cps.obs.exp.caggac+cps.obs.exp.caggag+cps.obs.exp.caggat+cps.obs.exp.caggca+cps.obs.exp.caggcc+cps.obs.exp.caggcg+cps.obs.exp.caggct+cps.obs.exp.caggga+cps.obs.exp.cagggc+cps.obs.exp.cagggg+cps.obs.exp.cagggt+cps.obs.exp.caggta+cps.obs.exp.caggtc+cps.obs.exp.caggtg+cps.obs.exp.caggtt+cps.obs.exp.cagtac+cps.obs.exp.cagtat+cps.obs.exp.cagtca+cps.obs.exp.cagtcc+cps.obs.exp.cagtcg+cps.obs.exp.cagtct+cps.obs.exp.cagtgc+cps.obs.exp.cagtgg+cps.obs.exp.cagtgt+cps.obs.exp.cagtta+cps.obs.exp.cagttc+cps.obs.exp.cagttg+cps.obs.exp.cagttt+cps.obs.exp.cataaa+cps.obs.exp.cataac+cps.obs.exp.cataag+cps.obs.exp.cataat+cps.obs.exp.cataca+cps.obs.exp.catacc+cps.obs.exp.catacg+cps.obs.exp.catact+cps.obs.exp.cataga+cps.obs.exp.catagc+cps.obs.exp.catagg+cps.obs.exp.catagt+cps.obs.exp.catata+cps.obs.exp.catatc+cps.obs.exp.catatg+cps.obs.exp.catatt+cps.obs.exp.catcaa+cps.obs.exp.catcac+cps.obs.exp.catcag+cps.obs.exp.catcat+cps.obs.exp.catcca+cps.obs.exp.catccc+cps.obs.exp.catccg+cps.obs.exp.catcct+cps.obs.exp.catcga+cps.obs.exp.catcgc+cps.obs.exp.catcgg+cps.obs.exp.catcgt+cps.obs.exp.catcta+cps.obs.exp.catctc+cps.obs.exp.catctg+cps.obs.exp.catctt+cps.obs.exp.catgaa+cps.obs.exp.catgac+cps.obs.exp.catgag+cps.obs.exp.catgat+cps.obs.exp.catgca+cps.obs.exp.catgcc+cps.obs.exp.catgcg+cps.obs.exp.catgct+cps.obs.exp.catgga+cps.obs.exp.catggc+cps.obs.exp.catggg+cps.obs.exp.catggt+cps.obs.exp.catgta+cps.obs.exp.catgtc+cps.obs.exp.catgtg+cps.obs.exp.catgtt+cps.obs.exp.cattac+cps.obs.exp.cattat+cps.obs.exp.cattca+cps.obs.exp.cattcc+cps.obs.exp.cattcg+cps.obs.exp.cattct+cps.obs.exp.cattgc+cps.obs.exp.cattgg+cps.obs.exp.cattgt+cps.obs.exp.cattta+cps.obs.exp.catttc+cps.obs.exp.catttg+cps.obs.exp.catttt+cps.obs.exp.ccaaaa+cps.obs.exp.ccaaac+cps.obs.exp.ccaaag+cps.obs.exp.ccaaat+cps.obs.exp.ccaaca+cps.obs.exp.ccaacc+cps.obs.exp.ccaacg+cps.obs.exp.ccaact+cps.obs.exp.ccaaga+cps.obs.exp.ccaagc+cps.obs.exp.ccaagg+cps.obs.exp.ccaagt+cps.obs.exp.ccaata+cps.obs.exp.ccaatc+cps.obs.exp.ccaatg+cps.obs.exp.ccaatt+cps.obs.exp.ccacaa+cps.obs.exp.ccacac+cps.obs.exp.ccacag+cps.obs.exp.ccacat+cps.obs.exp.ccacca+cps.obs.exp.ccaccc+cps.obs.exp.ccaccg+cps.obs.exp.ccacct+cps.obs.exp.ccacga+cps.obs.exp.ccacgc+cps.obs.exp.ccacgg+cps.obs.exp.ccacgt+cps.obs.exp.ccacta+cps.obs.exp.ccactc+cps.obs.exp.ccactg+cps.obs.exp.ccactt+cps.obs.exp.ccagaa+cps.obs.exp.ccagac+cps.obs.exp.ccagag+cps.obs.exp.ccagat+cps.obs.exp.ccagca+cps.obs.exp.ccagcc+cps.obs.exp.ccagcg+cps.obs.exp.ccagct+cps.obs.exp.ccagga+cps.obs.exp.ccaggc+cps.obs.exp.ccaggg+cps.obs.exp.ccaggt+cps.obs.exp.ccagta+cps.obs.exp.ccagtc+cps.obs.exp.ccagtg+cps.obs.exp.ccagtt+cps.obs.exp.ccatac+cps.obs.exp.ccatat+cps.obs.exp.ccatca+cps.obs.exp.ccatcc+cps.obs.exp.ccatcg+cps.obs.exp.ccatct+cps.obs.exp.ccatgc+cps.obs.exp.ccatgg+cps.obs.exp.ccatgt+cps.obs.exp.ccatta+cps.obs.exp.ccattc+cps.obs.exp.ccattg+cps.obs.exp.ccattt+cps.obs.exp.cccaaa+cps.obs.exp.cccaac+cps.obs.exp.cccaag+cps.obs.exp.cccaat+cps.obs.exp.cccaca+cps.obs.exp.cccacc+cps.obs.exp.cccacg+cps.obs.exp.cccact+cps.obs.exp.cccaga+cps.obs.exp.cccagc+cps.obs.exp.cccagg+cps.obs.exp.cccagt+cps.obs.exp.cccata+cps.obs.exp.cccatc+cps.obs.exp.cccatg+cps.obs.exp.cccatt+cps.obs.exp.ccccaa+cps.obs.exp.ccccac+cps.obs.exp.ccccag+cps.obs.exp.ccccat+cps.obs.exp.ccccca+cps.obs.exp.cccccc+cps.obs.exp.cccccg+cps.obs.exp.ccccct+cps.obs.exp.ccccga+cps.obs.exp.ccccgc+cps.obs.exp.ccccgg+cps.obs.exp.ccccgt+cps.obs.exp.ccccta+cps.obs.exp.cccctc+cps.obs.exp.cccctg+cps.obs.exp.cccctt+cps.obs.exp.cccgaa+cps.obs.exp.cccgac+cps.obs.exp.cccgag+cps.obs.exp.cccgat+cps.obs.exp.cccgca+cps.obs.exp.cccgcc+cps.obs.exp.cccgcg+cps.obs.exp.cccgct+cps.obs.exp.cccgga+cps.obs.exp.cccggc+cps.obs.exp.cccggg+cps.obs.exp.cccggt+cps.obs.exp.cccgta+cps.obs.exp.cccgtc+cps.obs.exp.cccgtg+cps.obs.exp.cccgtt+cps.obs.exp.ccctac+cps.obs.exp.ccctat+cps.obs.exp.ccctca+cps.obs.exp.ccctcc+cps.obs.exp.ccctcg+cps.obs.exp.ccctct+cps.obs.exp.ccctgc+cps.obs.exp.ccctgg+cps.obs.exp.ccctgt+cps.obs.exp.ccctta+cps.obs.exp.cccttc+cps.obs.exp.cccttg+cps.obs.exp.cccttt+cps.obs.exp.ccgaaa+cps.obs.exp.ccgaac+cps.obs.exp.ccgaag+cps.obs.exp.ccgaat+cps.obs.exp.ccgaca+cps.obs.exp.ccgacc+cps.obs.exp.ccgacg+cps.obs.exp.ccgact+cps.obs.exp.ccgaga+cps.obs.exp.ccgagc+cps.obs.exp.ccgagg+cps.obs.exp.ccgagt+cps.obs.exp.ccgata+cps.obs.exp.ccgatc+cps.obs.exp.ccgatg+cps.obs.exp.ccgatt+cps.obs.exp.ccgcaa+cps.obs.exp.ccgcac+cps.obs.exp.ccgcag+cps.obs.exp.ccgcat+cps.obs.exp.ccgcca+cps.obs.exp.ccgccc+cps.obs.exp.ccgccg+cps.obs.exp.ccgcct+cps.obs.exp.ccgcga+cps.obs.exp.ccgcgc+cps.obs.exp.ccgcgg+cps.obs.exp.ccgcgt+cps.obs.exp.ccgcta+cps.obs.exp.ccgctc+cps.obs.exp.ccgctg+cps.obs.exp.ccgctt+cps.obs.exp.ccggaa+cps.obs.exp.ccggac+cps.obs.exp.ccggag+cps.obs.exp.ccggat+cps.obs.exp.ccggca+cps.obs.exp.ccggcc+cps.obs.exp.ccggcg+cps.obs.exp.ccggct+cps.obs.exp.ccggga+cps.obs.exp.ccgggc+cps.obs.exp.ccgggg+cps.obs.exp.ccgggt+cps.obs.exp.ccggta+cps.obs.exp.ccggtc+cps.obs.exp.ccggtg+cps.obs.exp.ccggtt+cps.obs.exp.ccgtac+cps.obs.exp.ccgtat+cps.obs.exp.ccgtca+cps.obs.exp.ccgtcc+cps.obs.exp.ccgtcg+cps.obs.exp.ccgtct+cps.obs.exp.ccgtgc+cps.obs.exp.ccgtgg+cps.obs.exp.ccgtgt+cps.obs.exp.ccgtta+cps.obs.exp.ccgttc+cps.obs.exp.ccgttg+cps.obs.exp.ccgttt+cps.obs.exp.cctaaa+cps.obs.exp.cctaac+cps.obs.exp.cctaag+cps.obs.exp.cctaat+cps.obs.exp.cctaca+cps.obs.exp.cctacc+cps.obs.exp.cctacg+cps.obs.exp.cctact+cps.obs.exp.cctaga+cps.obs.exp.cctagc+cps.obs.exp.cctagg+cps.obs.exp.cctagt+cps.obs.exp.cctata+cps.obs.exp.cctatc+cps.obs.exp.cctatg+cps.obs.exp.cctatt+cps.obs.exp.cctcaa+cps.obs.exp.cctcac+cps.obs.exp.cctcag+cps.obs.exp.cctcat+cps.obs.exp.cctcca+cps.obs.exp.cctccc+cps.obs.exp.cctccg+cps.obs.exp.cctcct+cps.obs.exp.cctcga+cps.obs.exp.cctcgc+cps.obs.exp.cctcgg+cps.obs.exp.cctcgt+cps.obs.exp.cctcta+cps.obs.exp.cctctc+cps.obs.exp.cctctg+cps.obs.exp.cctctt+cps.obs.exp.cctgaa+cps.obs.exp.cctgac+cps.obs.exp.cctgag+cps.obs.exp.cctgat+cps.obs.exp.cctgca+cps.obs.exp.cctgcc+cps.obs.exp.cctgcg+cps.obs.exp.cctgct+cps.obs.exp.cctgga+cps.obs.exp.cctggc+cps.obs.exp.cctggg+cps.obs.exp.cctggt+cps.obs.exp.cctgta+cps.obs.exp.cctgtc+cps.obs.exp.cctgtg+cps.obs.exp.cctgtt+cps.obs.exp.ccttac+cps.obs.exp.ccttat+cps.obs.exp.ccttca+cps.obs.exp.ccttcc+cps.obs.exp.ccttcg+cps.obs.exp.ccttct+cps.obs.exp.ccttgc+cps.obs.exp.ccttgg+cps.obs.exp.ccttgt+cps.obs.exp.ccttta+cps.obs.exp.cctttc+cps.obs.exp.cctttg+cps.obs.exp.cctttt+cps.obs.exp.cgaaaa+cps.obs.exp.cgaaac+cps.obs.exp.cgaaag+cps.obs.exp.cgaaat+cps.obs.exp.cgaaca+cps.obs.exp.cgaacc+cps.obs.exp.cgaacg+cps.obs.exp.cgaact+cps.obs.exp.cgaaga+cps.obs.exp.cgaagc+cps.obs.exp.cgaagg+cps.obs.exp.cgaagt+cps.obs.exp.cgaata+cps.obs.exp.cgaatc+cps.obs.exp.cgaatg+cps.obs.exp.cgaatt+cps.obs.exp.cgacaa+cps.obs.exp.cgacac+cps.obs.exp.cgacag+cps.obs.exp.cgacat+cps.obs.exp.cgacca+cps.obs.exp.cgaccc+cps.obs.exp.cgaccg+cps.obs.exp.cgacct+cps.obs.exp.cgacga+cps.obs.exp.cgacgc+cps.obs.exp.cgacgg+cps.obs.exp.cgacgt+cps.obs.exp.cgacta+cps.obs.exp.cgactc+cps.obs.exp.cgactg+cps.obs.exp.cgactt+cps.obs.exp.cgagaa+cps.obs.exp.cgagac+cps.obs.exp.cgagag+cps.obs.exp.cgagat+cps.obs.exp.cgagca+cps.obs.exp.cgagcc+cps.obs.exp.cgagcg+cps.obs.exp.cgagct+cps.obs.exp.cgagga+cps.obs.exp.cgaggc+cps.obs.exp.cgaggg+cps.obs.exp.cgaggt+cps.obs.exp.cgagta+cps.obs.exp.cgagtc+cps.obs.exp.cgagtg+cps.obs.exp.cgagtt+cps.obs.exp.cgatac+cps.obs.exp.cgatat+cps.obs.exp.cgatca+cps.obs.exp.cgatcc+cps.obs.exp.cgatcg+cps.obs.exp.cgatct+cps.obs.exp.cgatgc+cps.obs.exp.cgatgg+cps.obs.exp.cgatgt+cps.obs.exp.cgatta+cps.obs.exp.cgattc+cps.obs.exp.cgattg+cps.obs.exp.cgattt+cps.obs.exp.cgcaaa+cps.obs.exp.cgcaac+cps.obs.exp.cgcaag+cps.obs.exp.cgcaat+cps.obs.exp.cgcaca+cps.obs.exp.cgcacc+cps.obs.exp.cgcacg+cps.obs.exp.cgcact+cps.obs.exp.cgcaga+cps.obs.exp.cgcagc+cps.obs.exp.cgcagg+cps.obs.exp.cgcagt+cps.obs.exp.cgcata+cps.obs.exp.cgcatc+cps.obs.exp.cgcatg+cps.obs.exp.cgcatt+cps.obs.exp.cgccaa+cps.obs.exp.cgccac+cps.obs.exp.cgccag+cps.obs.exp.cgccat+cps.obs.exp.cgccca+cps.obs.exp.cgcccc+cps.obs.exp.cgcccg+cps.obs.exp.cgccct+cps.obs.exp.cgccga+cps.obs.exp.cgccgc+cps.obs.exp.cgccgg+cps.obs.exp.cgccgt+cps.obs.exp.cgccta+cps.obs.exp.cgcctc+cps.obs.exp.cgcctg+cps.obs.exp.cgcctt+cps.obs.exp.cgcgaa+cps.obs.exp.cgcgac+cps.obs.exp.cgcgag+cps.obs.exp.cgcgat+cps.obs.exp.cgcgca+cps.obs.exp.cgcgcc+cps.obs.exp.cgcgcg+cps.obs.exp.cgcgct+cps.obs.exp.cgcgga+cps.obs.exp.cgcggc+cps.obs.exp.cgcggg+cps.obs.exp.cgcggt+cps.obs.exp.cgcgta+cps.obs.exp.cgcgtc+cps.obs.exp.cgcgtg+cps.obs.exp.cgcgtt+cps.obs.exp.cgctac+cps.obs.exp.cgctat+cps.obs.exp.cgctca+cps.obs.exp.cgctcc+cps.obs.exp.cgctcg+cps.obs.exp.cgctct+cps.obs.exp.cgctgc+cps.obs.exp.cgctgg+cps.obs.exp.cgctgt+cps.obs.exp.cgctta+cps.obs.exp.cgcttc+cps.obs.exp.cgcttg+cps.obs.exp.cgcttt+cps.obs.exp.cggaaa+cps.obs.exp.cggaac+cps.obs.exp.cggaag+cps.obs.exp.cggaat+cps.obs.exp.cggaca+cps.obs.exp.cggacc+cps.obs.exp.cggacg+cps.obs.exp.cggact+cps.obs.exp.cggaga+cps.obs.exp.cggagc+cps.obs.exp.cggagg+cps.obs.exp.cggagt+cps.obs.exp.cggata+cps.obs.exp.cggatc+cps.obs.exp.cggatg+cps.obs.exp.cggatt+cps.obs.exp.cggcaa+cps.obs.exp.cggcac+cps.obs.exp.cggcag+cps.obs.exp.cggcat+cps.obs.exp.cggcca+cps.obs.exp.cggccc+cps.obs.exp.cggccg+cps.obs.exp.cggcct+cps.obs.exp.cggcga+cps.obs.exp.cggcgc+cps.obs.exp.cggcgg+cps.obs.exp.cggcgt+cps.obs.exp.cggcta+cps.obs.exp.cggctc+cps.obs.exp.cggctg+cps.obs.exp.cggctt+cps.obs.exp.cgggaa+cps.obs.exp.cgggac+cps.obs.exp.cgggag+cps.obs.exp.cgggat+cps.obs.exp.cgggca+cps.obs.exp.cgggcc+cps.obs.exp.cgggcg+cps.obs.exp.cgggct+cps.obs.exp.cgggga+cps.obs.exp.cggggc+cps.obs.exp.cggggg+cps.obs.exp.cggggt+cps.obs.exp.cgggta+cps.obs.exp.cgggtc+cps.obs.exp.cgggtg+cps.obs.exp.cgggtt+cps.obs.exp.cggtac+cps.obs.exp.cggtat+cps.obs.exp.cggtca+cps.obs.exp.cggtcc+cps.obs.exp.cggtcg+cps.obs.exp.cggtct+cps.obs.exp.cggtgc+cps.obs.exp.cggtgg+cps.obs.exp.cggtgt+cps.obs.exp.cggtta+cps.obs.exp.cggttc+cps.obs.exp.cggttg+cps.obs.exp.cggttt+cps.obs.exp.cgtaaa+cps.obs.exp.cgtaac+cps.obs.exp.cgtaag+cps.obs.exp.cgtaat+cps.obs.exp.cgtaca+cps.obs.exp.cgtacc+cps.obs.exp.cgtacg+cps.obs.exp.cgtact+cps.obs.exp.cgtaga+cps.obs.exp.cgtagc+cps.obs.exp.cgtagg+cps.obs.exp.cgtagt+cps.obs.exp.cgtata+cps.obs.exp.cgtatc+cps.obs.exp.cgtatg+cps.obs.exp.cgtatt+cps.obs.exp.cgtcaa+cps.obs.exp.cgtcac+cps.obs.exp.cgtcag+cps.obs.exp.cgtcat+cps.obs.exp.cgtcca+cps.obs.exp.cgtccc+cps.obs.exp.cgtccg+cps.obs.exp.cgtcct+cps.obs.exp.cgtcga+cps.obs.exp.cgtcgc+cps.obs.exp.cgtcgg+cps.obs.exp.cgtcgt+cps.obs.exp.cgtcta+cps.obs.exp.cgtctc+cps.obs.exp.cgtctg+cps.obs.exp.cgtctt+cps.obs.exp.cgtgaa+cps.obs.exp.cgtgac+cps.obs.exp.cgtgag+cps.obs.exp.cgtgat+cps.obs.exp.cgtgca+cps.obs.exp.cgtgcc+cps.obs.exp.cgtgcg+cps.obs.exp.cgtgct+cps.obs.exp.cgtgga+cps.obs.exp.cgtggc+cps.obs.exp.cgtggg+cps.obs.exp.cgtggt+cps.obs.exp.cgtgta+cps.obs.exp.cgtgtc+cps.obs.exp.cgtgtg+cps.obs.exp.cgtgtt+cps.obs.exp.cgttac+cps.obs.exp.cgttat+cps.obs.exp.cgttca+cps.obs.exp.cgttcc+cps.obs.exp.cgttcg+cps.obs.exp.cgttct+cps.obs.exp.cgttgc+cps.obs.exp.cgttgg+cps.obs.exp.cgttgt+cps.obs.exp.cgttta+cps.obs.exp.cgtttc+cps.obs.exp.cgtttg+cps.obs.exp.cgtttt+cps.obs.exp.ctaaaa+cps.obs.exp.ctaaac+cps.obs.exp.ctaaag+cps.obs.exp.ctaaat+cps.obs.exp.ctaaca+cps.obs.exp.ctaacc+cps.obs.exp.ctaacg+cps.obs.exp.ctaact+cps.obs.exp.ctaaga+cps.obs.exp.ctaagc+cps.obs.exp.ctaagg+cps.obs.exp.ctaagt+cps.obs.exp.ctaata+cps.obs.exp.ctaatc+cps.obs.exp.ctaatg+cps.obs.exp.ctaatt+cps.obs.exp.ctacaa+cps.obs.exp.ctacac+cps.obs.exp.ctacag+cps.obs.exp.ctacat+cps.obs.exp.ctacca+cps.obs.exp.ctaccc+cps.obs.exp.ctaccg+cps.obs.exp.ctacct+cps.obs.exp.ctacga+cps.obs.exp.ctacgc+cps.obs.exp.ctacgg+cps.obs.exp.ctacgt+cps.obs.exp.ctacta+cps.obs.exp.ctactc+cps.obs.exp.ctactg+cps.obs.exp.ctactt+cps.obs.exp.ctagaa+cps.obs.exp.ctagac+cps.obs.exp.ctagag+cps.obs.exp.ctagat+cps.obs.exp.ctagca+cps.obs.exp.ctagcc+cps.obs.exp.ctagcg+cps.obs.exp.ctagct+cps.obs.exp.ctagga+cps.obs.exp.ctaggc+cps.obs.exp.ctaggg+cps.obs.exp.ctaggt+cps.obs.exp.ctagta+cps.obs.exp.ctagtc+cps.obs.exp.ctagtg+cps.obs.exp.ctagtt+cps.obs.exp.ctatac+cps.obs.exp.ctatat+cps.obs.exp.ctatca+cps.obs.exp.ctatcc+cps.obs.exp.ctatcg+cps.obs.exp.ctatct+cps.obs.exp.ctatgc+cps.obs.exp.ctatgg+cps.obs.exp.ctatgt+cps.obs.exp.ctatta+cps.obs.exp.ctattc+cps.obs.exp.ctattg+cps.obs.exp.ctattt+cps.obs.exp.ctcaaa+cps.obs.exp.ctcaac+cps.obs.exp.ctcaag+cps.obs.exp.ctcaat+cps.obs.exp.ctcaca+cps.obs.exp.ctcacc+cps.obs.exp.ctcacg+cps.obs.exp.ctcact+cps.obs.exp.ctcaga+cps.obs.exp.ctcagc+cps.obs.exp.ctcagg+cps.obs.exp.ctcagt+cps.obs.exp.ctcata+cps.obs.exp.ctcatc+cps.obs.exp.ctcatg+cps.obs.exp.ctcatt+cps.obs.exp.ctccaa+cps.obs.exp.ctccac+cps.obs.exp.ctccag+cps.obs.exp.ctccat+cps.obs.exp.ctccca+cps.obs.exp.ctcccc+cps.obs.exp.ctcccg+cps.obs.exp.ctccct+cps.obs.exp.ctccga+cps.obs.exp.ctccgc+cps.obs.exp.ctccgg+cps.obs.exp.ctccgt+cps.obs.exp.ctccta+cps.obs.exp.ctcctc+cps.obs.exp.ctcctg+cps.obs.exp.ctcctt+cps.obs.exp.ctcgaa+cps.obs.exp.ctcgac+cps.obs.exp.ctcgag+cps.obs.exp.ctcgat+cps.obs.exp.ctcgca+cps.obs.exp.ctcgcc+cps.obs.exp.ctcgcg+cps.obs.exp.ctcgct+cps.obs.exp.ctcgga+cps.obs.exp.ctcggc+cps.obs.exp.ctcggg+cps.obs.exp.ctcggt+cps.obs.exp.ctcgta+cps.obs.exp.ctcgtc+cps.obs.exp.ctcgtg+cps.obs.exp.ctcgtt+cps.obs.exp.ctctac+cps.obs.exp.ctctat+cps.obs.exp.ctctca+cps.obs.exp.ctctcc+cps.obs.exp.ctctcg+cps.obs.exp.ctctct+cps.obs.exp.ctctgc+cps.obs.exp.ctctgg+cps.obs.exp.ctctgt+cps.obs.exp.ctctta+cps.obs.exp.ctcttc+cps.obs.exp.ctcttg+cps.obs.exp.ctcttt+cps.obs.exp.ctgaaa+cps.obs.exp.ctgaac+cps.obs.exp.ctgaag+cps.obs.exp.ctgaat+cps.obs.exp.ctgaca+cps.obs.exp.ctgacc+cps.obs.exp.ctgacg+cps.obs.exp.ctgact+cps.obs.exp.ctgaga+cps.obs.exp.ctgagc+cps.obs.exp.ctgagg+cps.obs.exp.ctgagt+cps.obs.exp.ctgata+cps.obs.exp.ctgatc+cps.obs.exp.ctgatg+cps.obs.exp.ctgatt+cps.obs.exp.ctgcaa+cps.obs.exp.ctgcac+cps.obs.exp.ctgcag+cps.obs.exp.ctgcat+cps.obs.exp.ctgcca+cps.obs.exp.ctgccc+cps.obs.exp.ctgccg+cps.obs.exp.ctgcct+cps.obs.exp.ctgcga+cps.obs.exp.ctgcgc+cps.obs.exp.ctgcgg+cps.obs.exp.ctgcgt+cps.obs.exp.ctgcta+cps.obs.exp.ctgctc+cps.obs.exp.ctgctg+cps.obs.exp.ctgctt+cps.obs.exp.ctggaa+cps.obs.exp.ctggac+cps.obs.exp.ctggag+cps.obs.exp.ctggat+cps.obs.exp.ctggca+cps.obs.exp.ctggcc+cps.obs.exp.ctggcg+cps.obs.exp.ctggct+cps.obs.exp.ctggga+cps.obs.exp.ctgggc+cps.obs.exp.ctgggg+cps.obs.exp.ctgggt+cps.obs.exp.ctggta+cps.obs.exp.ctggtc+cps.obs.exp.ctggtg+cps.obs.exp.ctggtt+cps.obs.exp.ctgtac+cps.obs.exp.ctgtat+cps.obs.exp.ctgtca+cps.obs.exp.ctgtcc+cps.obs.exp.ctgtcg+cps.obs.exp.ctgtct+cps.obs.exp.ctgtgc+cps.obs.exp.ctgtgg+cps.obs.exp.ctgtgt+cps.obs.exp.ctgtta+cps.obs.exp.ctgttc+cps.obs.exp.ctgttg+cps.obs.exp.ctgttt+cps.obs.exp.cttaaa+cps.obs.exp.cttaac+cps.obs.exp.cttaag+cps.obs.exp.cttaat+cps.obs.exp.cttaca+cps.obs.exp.cttacc+cps.obs.exp.cttacg+cps.obs.exp.cttact+cps.obs.exp.cttaga+cps.obs.exp.cttagc+cps.obs.exp.cttagg+cps.obs.exp.cttagt+cps.obs.exp.cttata+cps.obs.exp.cttatc+cps.obs.exp.cttatg+cps.obs.exp.cttatt+cps.obs.exp.cttcaa+cps.obs.exp.cttcac+cps.obs.exp.cttcag+cps.obs.exp.cttcat+cps.obs.exp.cttcca+cps.obs.exp.cttccc+cps.obs.exp.cttccg+cps.obs.exp.cttcct+cps.obs.exp.cttcga+cps.obs.exp.cttcgc+cps.obs.exp.cttcgg+cps.obs.exp.cttcgt+cps.obs.exp.cttcta+cps.obs.exp.cttctc+cps.obs.exp.cttctg+cps.obs.exp.cttctt+cps.obs.exp.cttgaa+cps.obs.exp.cttgac+cps.obs.exp.cttgag+cps.obs.exp.cttgat+cps.obs.exp.cttgca+cps.obs.exp.cttgcc+cps.obs.exp.cttgcg+cps.obs.exp.cttgct+cps.obs.exp.cttgga+cps.obs.exp.cttggc+cps.obs.exp.cttggg+cps.obs.exp.cttggt+cps.obs.exp.cttgta+cps.obs.exp.cttgtc+cps.obs.exp.cttgtg+cps.obs.exp.cttgtt+cps.obs.exp.ctttac+cps.obs.exp.ctttat+cps.obs.exp.ctttca+cps.obs.exp.ctttcc+cps.obs.exp.ctttcg+cps.obs.exp.ctttct+cps.obs.exp.ctttgc+cps.obs.exp.ctttgg+cps.obs.exp.ctttgt+cps.obs.exp.ctttta+cps.obs.exp.cttttc+cps.obs.exp.cttttg+cps.obs.exp.cttttt+cps.obs.exp.gaaaaa+cps.obs.exp.gaaaac+cps.obs.exp.gaaaag+cps.obs.exp.gaaaat+cps.obs.exp.gaaaca+cps.obs.exp.gaaacc+cps.obs.exp.gaaacg+cps.obs.exp.gaaact+cps.obs.exp.gaaaga+cps.obs.exp.gaaagc+cps.obs.exp.gaaagg+cps.obs.exp.gaaagt+cps.obs.exp.gaaata+cps.obs.exp.gaaatc+cps.obs.exp.gaaatg+cps.obs.exp.gaaatt+cps.obs.exp.gaacaa+cps.obs.exp.gaacac+cps.obs.exp.gaacag+cps.obs.exp.gaacat+cps.obs.exp.gaacca+cps.obs.exp.gaaccc+cps.obs.exp.gaaccg+cps.obs.exp.gaacct+cps.obs.exp.gaacga+cps.obs.exp.gaacgc+cps.obs.exp.gaacgg+cps.obs.exp.gaacgt+cps.obs.exp.gaacta+cps.obs.exp.gaactc+cps.obs.exp.gaactg+cps.obs.exp.gaactt+cps.obs.exp.gaagaa+cps.obs.exp.gaagac+cps.obs.exp.gaagag+cps.obs.exp.gaagat+cps.obs.exp.gaagca+cps.obs.exp.gaagcc+cps.obs.exp.gaagcg+cps.obs.exp.gaagct+cps.obs.exp.gaagga+cps.obs.exp.gaaggc+cps.obs.exp.gaaggg+cps.obs.exp.gaaggt+cps.obs.exp.gaagta+cps.obs.exp.gaagtc+cps.obs.exp.gaagtg+cps.obs.exp.gaagtt+cps.obs.exp.gaatac+cps.obs.exp.gaatat+cps.obs.exp.gaatca+cps.obs.exp.gaatcc+cps.obs.exp.gaatcg+cps.obs.exp.gaatct+cps.obs.exp.gaatgc+cps.obs.exp.gaatgg+cps.obs.exp.gaatgt+cps.obs.exp.gaatta+cps.obs.exp.gaattc+cps.obs.exp.gaattg+cps.obs.exp.gaattt+cps.obs.exp.gacaaa+cps.obs.exp.gacaac+cps.obs.exp.gacaag+cps.obs.exp.gacaat+cps.obs.exp.gacaca+cps.obs.exp.gacacc+cps.obs.exp.gacacg+cps.obs.exp.gacact+cps.obs.exp.gacaga+cps.obs.exp.gacagc+cps.obs.exp.gacagg+cps.obs.exp.gacagt+cps.obs.exp.gacata+cps.obs.exp.gacatc+cps.obs.exp.gacatg+cps.obs.exp.gacatt+cps.obs.exp.gaccaa+cps.obs.exp.gaccac+cps.obs.exp.gaccag+cps.obs.exp.gaccat+cps.obs.exp.gaccca+cps.obs.exp.gacccc+cps.obs.exp.gacccg+cps.obs.exp.gaccct+cps.obs.exp.gaccga+cps.obs.exp.gaccgc+cps.obs.exp.gaccgg+cps.obs.exp.gaccgt+cps.obs.exp.gaccta+cps.obs.exp.gacctc+cps.obs.exp.gacctg+cps.obs.exp.gacctt+cps.obs.exp.gacgaa+cps.obs.exp.gacgac+cps.obs.exp.gacgag+cps.obs.exp.gacgat+cps.obs.exp.gacgca+cps.obs.exp.gacgcc+cps.obs.exp.gacgcg+cps.obs.exp.gacgct+cps.obs.exp.gacgga+cps.obs.exp.gacggc+cps.obs.exp.gacggg+cps.obs.exp.gacggt+cps.obs.exp.gacgta+cps.obs.exp.gacgtc+cps.obs.exp.gacgtg+cps.obs.exp.gacgtt+cps.obs.exp.gactac+cps.obs.exp.gactat+cps.obs.exp.gactca+cps.obs.exp.gactcc+cps.obs.exp.gactcg+cps.obs.exp.gactct+cps.obs.exp.gactgc+cps.obs.exp.gactgg+cps.obs.exp.gactgt+cps.obs.exp.gactta+cps.obs.exp.gacttc+cps.obs.exp.gacttg+cps.obs.exp.gacttt+cps.obs.exp.gagaaa+cps.obs.exp.gagaac+cps.obs.exp.gagaag+cps.obs.exp.gagaat+cps.obs.exp.gagaca+cps.obs.exp.gagacc+cps.obs.exp.gagacg+cps.obs.exp.gagact+cps.obs.exp.gagaga+cps.obs.exp.gagagc+cps.obs.exp.gagagg+cps.obs.exp.gagagt+cps.obs.exp.gagata+cps.obs.exp.gagatc+cps.obs.exp.gagatg+cps.obs.exp.gagatt+cps.obs.exp.gagcaa+cps.obs.exp.gagcac+cps.obs.exp.gagcag+cps.obs.exp.gagcat+cps.obs.exp.gagcca+cps.obs.exp.gagccc+cps.obs.exp.gagccg+cps.obs.exp.gagcct+cps.obs.exp.gagcga+cps.obs.exp.gagcgc+cps.obs.exp.gagcgg+cps.obs.exp.gagcgt+cps.obs.exp.gagcta+cps.obs.exp.gagctc+cps.obs.exp.gagctg+cps.obs.exp.gagctt+cps.obs.exp.gaggaa+cps.obs.exp.gaggac+cps.obs.exp.gaggag+cps.obs.exp.gaggat+cps.obs.exp.gaggca+cps.obs.exp.gaggcc+cps.obs.exp.gaggcg+cps.obs.exp.gaggct+cps.obs.exp.gaggga+cps.obs.exp.gagggc+cps.obs.exp.gagggg+cps.obs.exp.gagggt+cps.obs.exp.gaggta+cps.obs.exp.gaggtc+cps.obs.exp.gaggtg+cps.obs.exp.gaggtt+cps.obs.exp.gagtac+cps.obs.exp.gagtat+cps.obs.exp.gagtca+cps.obs.exp.gagtcc+cps.obs.exp.gagtcg+cps.obs.exp.gagtct+cps.obs.exp.gagtgc+cps.obs.exp.gagtgg+cps.obs.exp.gagtgt+cps.obs.exp.gagtta+cps.obs.exp.gagttc+cps.obs.exp.gagttg+cps.obs.exp.gagttt+cps.obs.exp.gataaa+cps.obs.exp.gataac+cps.obs.exp.gataag+cps.obs.exp.gataat+cps.obs.exp.gataca+cps.obs.exp.gatacc+cps.obs.exp.gatacg+cps.obs.exp.gatact+cps.obs.exp.gataga+cps.obs.exp.gatagc+cps.obs.exp.gatagg+cps.obs.exp.gatagt+cps.obs.exp.gatata+cps.obs.exp.gatatc+cps.obs.exp.gatatg+cps.obs.exp.gatatt+cps.obs.exp.gatcaa+cps.obs.exp.gatcac+cps.obs.exp.gatcag+cps.obs.exp.gatcat+cps.obs.exp.gatcca+cps.obs.exp.gatccc+cps.obs.exp.gatccg+cps.obs.exp.gatcct+cps.obs.exp.gatcga+cps.obs.exp.gatcgc+cps.obs.exp.gatcgg+cps.obs.exp.gatcgt+cps.obs.exp.gatcta+cps.obs.exp.gatctc+cps.obs.exp.gatctg+cps.obs.exp.gatctt+cps.obs.exp.gatgaa+cps.obs.exp.gatgac+cps.obs.exp.gatgag+cps.obs.exp.gatgat+cps.obs.exp.gatgca+cps.obs.exp.gatgcc+cps.obs.exp.gatgcg+cps.obs.exp.gatgct+cps.obs.exp.gatgga+cps.obs.exp.gatggc+cps.obs.exp.gatggg+cps.obs.exp.gatggt+cps.obs.exp.gatgta+cps.obs.exp.gatgtc+cps.obs.exp.gatgtg+cps.obs.exp.gatgtt+cps.obs.exp.gattac+cps.obs.exp.gattat+cps.obs.exp.gattca+cps.obs.exp.gattcc+cps.obs.exp.gattcg+cps.obs.exp.gattct+cps.obs.exp.gattgc+cps.obs.exp.gattgg+cps.obs.exp.gattgt+cps.obs.exp.gattta+cps.obs.exp.gatttc+cps.obs.exp.gatttg+cps.obs.exp.gatttt+cps.obs.exp.gcaaaa+cps.obs.exp.gcaaac+cps.obs.exp.gcaaag+cps.obs.exp.gcaaat+cps.obs.exp.gcaaca+cps.obs.exp.gcaacc+cps.obs.exp.gcaacg+cps.obs.exp.gcaact+cps.obs.exp.gcaaga+cps.obs.exp.gcaagc+cps.obs.exp.gcaagg+cps.obs.exp.gcaagt+cps.obs.exp.gcaata+cps.obs.exp.gcaatc+cps.obs.exp.gcaatg+cps.obs.exp.gcaatt+cps.obs.exp.gcacaa+cps.obs.exp.gcacac+cps.obs.exp.gcacag+cps.obs.exp.gcacat+cps.obs.exp.gcacca+cps.obs.exp.gcaccc+cps.obs.exp.gcaccg+cps.obs.exp.gcacct+cps.obs.exp.gcacga+cps.obs.exp.gcacgc+cps.obs.exp.gcacgg+cps.obs.exp.gcacgt+cps.obs.exp.gcacta+cps.obs.exp.gcactc+cps.obs.exp.gcactg+cps.obs.exp.gcactt+cps.obs.exp.gcagaa+cps.obs.exp.gcagac+cps.obs.exp.gcagag+cps.obs.exp.gcagat+cps.obs.exp.gcagca+cps.obs.exp.gcagcc+cps.obs.exp.gcagcg+cps.obs.exp.gcagct+cps.obs.exp.gcagga+cps.obs.exp.gcaggc+cps.obs.exp.gcaggg+cps.obs.exp.gcaggt+cps.obs.exp.gcagta+cps.obs.exp.gcagtc+cps.obs.exp.gcagtg+cps.obs.exp.gcagtt+cps.obs.exp.gcatac+cps.obs.exp.gcatat+cps.obs.exp.gcatca+cps.obs.exp.gcatcc+cps.obs.exp.gcatcg+cps.obs.exp.gcatct+cps.obs.exp.gcatgc+cps.obs.exp.gcatgg+cps.obs.exp.gcatgt+cps.obs.exp.gcatta+cps.obs.exp.gcattc+cps.obs.exp.gcattg+cps.obs.exp.gcattt+cps.obs.exp.gccaaa+cps.obs.exp.gccaac+cps.obs.exp.gccaag+cps.obs.exp.gccaat+cps.obs.exp.gccaca+cps.obs.exp.gccacc+cps.obs.exp.gccacg+cps.obs.exp.gccact+cps.obs.exp.gccaga+cps.obs.exp.gccagc+cps.obs.exp.gccagg+cps.obs.exp.gccagt+cps.obs.exp.gccata+cps.obs.exp.gccatc+cps.obs.exp.gccatg+cps.obs.exp.gccatt+cps.obs.exp.gcccaa+cps.obs.exp.gcccac+cps.obs.exp.gcccag+cps.obs.exp.gcccat+cps.obs.exp.gcccca+cps.obs.exp.gccccc+cps.obs.exp.gccccg+cps.obs.exp.gcccct+cps.obs.exp.gcccga+cps.obs.exp.gcccgc+cps.obs.exp.gcccgg+cps.obs.exp.gcccgt+cps.obs.exp.gcccta+cps.obs.exp.gccctc+cps.obs.exp.gccctg+cps.obs.exp.gccctt+cps.obs.exp.gccgaa+cps.obs.exp.gccgac+cps.obs.exp.gccgag+cps.obs.exp.gccgat+cps.obs.exp.gccgca+cps.obs.exp.gccgcc+cps.obs.exp.gccgcg+cps.obs.exp.gccgct+cps.obs.exp.gccgga+cps.obs.exp.gccggc+cps.obs.exp.gccggg+cps.obs.exp.gccggt+cps.obs.exp.gccgta+cps.obs.exp.gccgtc+cps.obs.exp.gccgtg+cps.obs.exp.gccgtt+cps.obs.exp.gcctac+cps.obs.exp.gcctat+cps.obs.exp.gcctca+cps.obs.exp.gcctcc+cps.obs.exp.gcctcg+cps.obs.exp.gcctct+cps.obs.exp.gcctgc+cps.obs.exp.gcctgg+cps.obs.exp.gcctgt+cps.obs.exp.gcctta+cps.obs.exp.gccttc+cps.obs.exp.gccttg+cps.obs.exp.gccttt+cps.obs.exp.gcgaaa+cps.obs.exp.gcgaac+cps.obs.exp.gcgaag+cps.obs.exp.gcgaat+cps.obs.exp.gcgaca+cps.obs.exp.gcgacc+cps.obs.exp.gcgacg+cps.obs.exp.gcgact+cps.obs.exp.gcgaga+cps.obs.exp.gcgagc+cps.obs.exp.gcgagg+cps.obs.exp.gcgagt+cps.obs.exp.gcgata+cps.obs.exp.gcgatc+cps.obs.exp.gcgatg+cps.obs.exp.gcgatt+cps.obs.exp.gcgcaa+cps.obs.exp.gcgcac+cps.obs.exp.gcgcag+cps.obs.exp.gcgcat+cps.obs.exp.gcgcca+cps.obs.exp.gcgccc+cps.obs.exp.gcgccg+cps.obs.exp.gcgcct+cps.obs.exp.gcgcga+cps.obs.exp.gcgcgc+cps.obs.exp.gcgcgg+cps.obs.exp.gcgcgt+cps.obs.exp.gcgcta+cps.obs.exp.gcgctc+cps.obs.exp.gcgctg+cps.obs.exp.gcgctt+cps.obs.exp.gcggaa+cps.obs.exp.gcggac+cps.obs.exp.gcggag+cps.obs.exp.gcggat+cps.obs.exp.gcggca+cps.obs.exp.gcggcc+cps.obs.exp.gcggcg+cps.obs.exp.gcggct+cps.obs.exp.gcggga+cps.obs.exp.gcgggc+cps.obs.exp.gcgggg+cps.obs.exp.gcgggt+cps.obs.exp.gcggta+cps.obs.exp.gcggtc+cps.obs.exp.gcggtg+cps.obs.exp.gcggtt+cps.obs.exp.gcgtac+cps.obs.exp.gcgtat+cps.obs.exp.gcgtca+cps.obs.exp.gcgtcc+cps.obs.exp.gcgtcg+cps.obs.exp.gcgtct+cps.obs.exp.gcgtgc+cps.obs.exp.gcgtgg+cps.obs.exp.gcgtgt+cps.obs.exp.gcgtta+cps.obs.exp.gcgttc+cps.obs.exp.gcgttg+cps.obs.exp.gcgttt+cps.obs.exp.gctaaa+cps.obs.exp.gctaac+cps.obs.exp.gctaag+cps.obs.exp.gctaat+cps.obs.exp.gctaca+cps.obs.exp.gctacc+cps.obs.exp.gctacg+cps.obs.exp.gctact+cps.obs.exp.gctaga+cps.obs.exp.gctagc+cps.obs.exp.gctagg+cps.obs.exp.gctagt+cps.obs.exp.gctata+cps.obs.exp.gctatc+cps.obs.exp.gctatg+cps.obs.exp.gctatt+cps.obs.exp.gctcaa+cps.obs.exp.gctcac+cps.obs.exp.gctcag+cps.obs.exp.gctcat+cps.obs.exp.gctcca+cps.obs.exp.gctccc+cps.obs.exp.gctccg+cps.obs.exp.gctcct+cps.obs.exp.gctcga+cps.obs.exp.gctcgc+cps.obs.exp.gctcgg+cps.obs.exp.gctcgt+cps.obs.exp.gctcta+cps.obs.exp.gctctc+cps.obs.exp.gctctg+cps.obs.exp.gctctt+cps.obs.exp.gctgaa+cps.obs.exp.gctgac+cps.obs.exp.gctgag+cps.obs.exp.gctgat+cps.obs.exp.gctgca+cps.obs.exp.gctgcc+cps.obs.exp.gctgcg+cps.obs.exp.gctgct+cps.obs.exp.gctgga+cps.obs.exp.gctggc+cps.obs.exp.gctggg+cps.obs.exp.gctggt+cps.obs.exp.gctgta+cps.obs.exp.gctgtc+cps.obs.exp.gctgtg+cps.obs.exp.gctgtt+cps.obs.exp.gcttac+cps.obs.exp.gcttat+cps.obs.exp.gcttca+cps.obs.exp.gcttcc+cps.obs.exp.gcttcg+cps.obs.exp.gcttct+cps.obs.exp.gcttgc+cps.obs.exp.gcttgg+cps.obs.exp.gcttgt+cps.obs.exp.gcttta+cps.obs.exp.gctttc+cps.obs.exp.gctttg+cps.obs.exp.gctttt+cps.obs.exp.ggaaaa+cps.obs.exp.ggaaac+cps.obs.exp.ggaaag+cps.obs.exp.ggaaat+cps.obs.exp.ggaaca+cps.obs.exp.ggaacc+cps.obs.exp.ggaacg+cps.obs.exp.ggaact+cps.obs.exp.ggaaga+cps.obs.exp.ggaagc+cps.obs.exp.ggaagg+cps.obs.exp.ggaagt+cps.obs.exp.ggaata+cps.obs.exp.ggaatc+cps.obs.exp.ggaatg+cps.obs.exp.ggaatt+cps.obs.exp.ggacaa+cps.obs.exp.ggacac+cps.obs.exp.ggacag+cps.obs.exp.ggacat+cps.obs.exp.ggacca+cps.obs.exp.ggaccc+cps.obs.exp.ggaccg+cps.obs.exp.ggacct+cps.obs.exp.ggacga+cps.obs.exp.ggacgc+cps.obs.exp.ggacgg+cps.obs.exp.ggacgt+cps.obs.exp.ggacta+cps.obs.exp.ggactc+cps.obs.exp.ggactg+cps.obs.exp.ggactt+cps.obs.exp.ggagaa+cps.obs.exp.ggagac+cps.obs.exp.ggagag+cps.obs.exp.ggagat+cps.obs.exp.ggagca+cps.obs.exp.ggagcc+cps.obs.exp.ggagcg+cps.obs.exp.ggagct+cps.obs.exp.ggagga+cps.obs.exp.ggaggc+cps.obs.exp.ggaggg+cps.obs.exp.ggaggt+cps.obs.exp.ggagta+cps.obs.exp.ggagtc+cps.obs.exp.ggagtg+cps.obs.exp.ggagtt+cps.obs.exp.ggatac+cps.obs.exp.ggatat+cps.obs.exp.ggatca+cps.obs.exp.ggatcc+cps.obs.exp.ggatcg+cps.obs.exp.ggatct+cps.obs.exp.ggatgc+cps.obs.exp.ggatgg+cps.obs.exp.ggatgt+cps.obs.exp.ggatta+cps.obs.exp.ggattc+cps.obs.exp.ggattg+cps.obs.exp.ggattt+cps.obs.exp.ggcaaa+cps.obs.exp.ggcaac+cps.obs.exp.ggcaag+cps.obs.exp.ggcaat+cps.obs.exp.ggcaca+cps.obs.exp.ggcacc+cps.obs.exp.ggcacg+cps.obs.exp.ggcact+cps.obs.exp.ggcaga+cps.obs.exp.ggcagc+cps.obs.exp.ggcagg+cps.obs.exp.ggcagt+cps.obs.exp.ggcata+cps.obs.exp.ggcatc+cps.obs.exp.ggcatg+cps.obs.exp.ggcatt+cps.obs.exp.ggccaa+cps.obs.exp.ggccac+cps.obs.exp.ggccag+cps.obs.exp.ggccat+cps.obs.exp.ggccca+cps.obs.exp.ggcccc+cps.obs.exp.ggcccg+cps.obs.exp.ggccct+cps.obs.exp.ggccga+cps.obs.exp.ggccgc+cps.obs.exp.ggccgg+cps.obs.exp.ggccgt+cps.obs.exp.ggccta+cps.obs.exp.ggcctc+cps.obs.exp.ggcctg+cps.obs.exp.ggcctt+cps.obs.exp.ggcgaa+cps.obs.exp.ggcgac+cps.obs.exp.ggcgag+cps.obs.exp.ggcgat+cps.obs.exp.ggcgca+cps.obs.exp.ggcgcc+cps.obs.exp.ggcgcg+cps.obs.exp.ggcgct+cps.obs.exp.ggcgga+cps.obs.exp.ggcggc+cps.obs.exp.ggcggg+cps.obs.exp.ggcggt+cps.obs.exp.ggcgta+cps.obs.exp.ggcgtc+cps.obs.exp.ggcgtg+cps.obs.exp.ggcgtt+cps.obs.exp.ggctac+cps.obs.exp.ggctat+cps.obs.exp.ggctca+cps.obs.exp.ggctcc+cps.obs.exp.ggctcg+cps.obs.exp.ggctct+cps.obs.exp.ggctgc+cps.obs.exp.ggctgg+cps.obs.exp.ggctgt+cps.obs.exp.ggctta+cps.obs.exp.ggcttc+cps.obs.exp.ggcttg+cps.obs.exp.ggcttt+cps.obs.exp.gggaaa+cps.obs.exp.gggaac+cps.obs.exp.gggaag+cps.obs.exp.gggaat+cps.obs.exp.gggaca+cps.obs.exp.gggacc+cps.obs.exp.gggacg+cps.obs.exp.gggact+cps.obs.exp.gggaga+cps.obs.exp.gggagc+cps.obs.exp.gggagg+cps.obs.exp.gggagt+cps.obs.exp.gggata+cps.obs.exp.gggatc+cps.obs.exp.gggatg+cps.obs.exp.gggatt+cps.obs.exp.gggcaa+cps.obs.exp.gggcac+cps.obs.exp.gggcag+cps.obs.exp.gggcat+cps.obs.exp.gggcca+cps.obs.exp.gggccc+cps.obs.exp.gggccg+cps.obs.exp.gggcct+cps.obs.exp.gggcga+cps.obs.exp.gggcgc+cps.obs.exp.gggcgg+cps.obs.exp.gggcgt+cps.obs.exp.gggcta+cps.obs.exp.gggctc+cps.obs.exp.gggctg+cps.obs.exp.gggctt+cps.obs.exp.ggggaa+cps.obs.exp.ggggac+cps.obs.exp.ggggag+cps.obs.exp.ggggat+cps.obs.exp.ggggca+cps.obs.exp.ggggcc+cps.obs.exp.ggggcg+cps.obs.exp.ggggct+cps.obs.exp.ggggga+cps.obs.exp.gggggc+cps.obs.exp.gggggg+cps.obs.exp.gggggt+cps.obs.exp.ggggta+cps.obs.exp.ggggtc+cps.obs.exp.ggggtg+cps.obs.exp.ggggtt+cps.obs.exp.gggtac+cps.obs.exp.gggtat+cps.obs.exp.gggtca+cps.obs.exp.gggtcc+cps.obs.exp.gggtcg+cps.obs.exp.gggtct+cps.obs.exp.gggtgc+cps.obs.exp.gggtgg+cps.obs.exp.gggtgt+cps.obs.exp.gggtta+cps.obs.exp.gggttc+cps.obs.exp.gggttg+cps.obs.exp.gggttt+cps.obs.exp.ggtaaa+cps.obs.exp.ggtaac+cps.obs.exp.ggtaag+cps.obs.exp.ggtaat+cps.obs.exp.ggtaca+cps.obs.exp.ggtacc+cps.obs.exp.ggtacg+cps.obs.exp.ggtact+cps.obs.exp.ggtaga+cps.obs.exp.ggtagc+cps.obs.exp.ggtagg+cps.obs.exp.ggtagt+cps.obs.exp.ggtata+cps.obs.exp.ggtatc+cps.obs.exp.ggtatg+cps.obs.exp.ggtatt+cps.obs.exp.ggtcaa+cps.obs.exp.ggtcac+cps.obs.exp.ggtcag+cps.obs.exp.ggtcat+cps.obs.exp.ggtcca+cps.obs.exp.ggtccc+cps.obs.exp.ggtccg+cps.obs.exp.ggtcct+cps.obs.exp.ggtcga+cps.obs.exp.ggtcgc+cps.obs.exp.ggtcgg+cps.obs.exp.ggtcgt+cps.obs.exp.ggtcta+cps.obs.exp.ggtctc+cps.obs.exp.ggtctg+cps.obs.exp.ggtctt+cps.obs.exp.ggtgaa+cps.obs.exp.ggtgac+cps.obs.exp.ggtgag+cps.obs.exp.ggtgat+cps.obs.exp.ggtgca+cps.obs.exp.ggtgcc+cps.obs.exp.ggtgcg+cps.obs.exp.ggtgct+cps.obs.exp.ggtgga+cps.obs.exp.ggtggc+cps.obs.exp.ggtggg+cps.obs.exp.ggtggt+cps.obs.exp.ggtgta+cps.obs.exp.ggtgtc+cps.obs.exp.ggtgtg+cps.obs.exp.ggtgtt+cps.obs.exp.ggttac+cps.obs.exp.ggttat+cps.obs.exp.ggttca+cps.obs.exp.ggttcc+cps.obs.exp.ggttcg+cps.obs.exp.ggttct+cps.obs.exp.ggttgc+cps.obs.exp.ggttgg+cps.obs.exp.ggttgt+cps.obs.exp.ggttta+cps.obs.exp.ggtttc+cps.obs.exp.ggtttg+cps.obs.exp.ggtttt+cps.obs.exp.gtaaaa+cps.obs.exp.gtaaac+cps.obs.exp.gtaaag+cps.obs.exp.gtaaat+cps.obs.exp.gtaaca+cps.obs.exp.gtaacc+cps.obs.exp.gtaacg+cps.obs.exp.gtaact+cps.obs.exp.gtaaga+cps.obs.exp.gtaagc+cps.obs.exp.gtaagg+cps.obs.exp.gtaagt+cps.obs.exp.gtaata+cps.obs.exp.gtaatc+cps.obs.exp.gtaatg+cps.obs.exp.gtaatt+cps.obs.exp.gtacaa+cps.obs.exp.gtacac+cps.obs.exp.gtacag+cps.obs.exp.gtacat+cps.obs.exp.gtacca+cps.obs.exp.gtaccc+cps.obs.exp.gtaccg+cps.obs.exp.gtacct+cps.obs.exp.gtacga+cps.obs.exp.gtacgc+cps.obs.exp.gtacgg+cps.obs.exp.gtacgt+cps.obs.exp.gtacta+cps.obs.exp.gtactc+cps.obs.exp.gtactg+cps.obs.exp.gtactt+cps.obs.exp.gtagaa+cps.obs.exp.gtagac+cps.obs.exp.gtagag+cps.obs.exp.gtagat+cps.obs.exp.gtagca+cps.obs.exp.gtagcc+cps.obs.exp.gtagcg+cps.obs.exp.gtagct+cps.obs.exp.gtagga+cps.obs.exp.gtaggc+cps.obs.exp.gtaggg+cps.obs.exp.gtaggt+cps.obs.exp.gtagta+cps.obs.exp.gtagtc+cps.obs.exp.gtagtg+cps.obs.exp.gtagtt+cps.obs.exp.gtatac+cps.obs.exp.gtatat+cps.obs.exp.gtatca+cps.obs.exp.gtatcc+cps.obs.exp.gtatcg+cps.obs.exp.gtatct+cps.obs.exp.gtatgc+cps.obs.exp.gtatgg+cps.obs.exp.gtatgt+cps.obs.exp.gtatta+cps.obs.exp.gtattc+cps.obs.exp.gtattg+cps.obs.exp.gtattt+cps.obs.exp.gtcaaa+cps.obs.exp.gtcaac+cps.obs.exp.gtcaag+cps.obs.exp.gtcaat+cps.obs.exp.gtcaca+cps.obs.exp.gtcacc+cps.obs.exp.gtcacg+cps.obs.exp.gtcact+cps.obs.exp.gtcaga+cps.obs.exp.gtcagc+cps.obs.exp.gtcagg+cps.obs.exp.gtcagt+cps.obs.exp.gtcata+cps.obs.exp.gtcatc+cps.obs.exp.gtcatg+cps.obs.exp.gtcatt+cps.obs.exp.gtccaa+cps.obs.exp.gtccac+cps.obs.exp.gtccag+cps.obs.exp.gtccat+cps.obs.exp.gtccca+cps.obs.exp.gtcccc+cps.obs.exp.gtcccg+cps.obs.exp.gtccct+cps.obs.exp.gtccga+cps.obs.exp.gtccgc+cps.obs.exp.gtccgg+cps.obs.exp.gtccgt+cps.obs.exp.gtccta+cps.obs.exp.gtcctc+cps.obs.exp.gtcctg+cps.obs.exp.gtcctt+cps.obs.exp.gtcgaa+cps.obs.exp.gtcgac+cps.obs.exp.gtcgag+cps.obs.exp.gtcgat+cps.obs.exp.gtcgca+cps.obs.exp.gtcgcc+cps.obs.exp.gtcgcg+cps.obs.exp.gtcgct+cps.obs.exp.gtcgga+cps.obs.exp.gtcggc+cps.obs.exp.gtcggg+cps.obs.exp.gtcggt+cps.obs.exp.gtcgta+cps.obs.exp.gtcgtc+cps.obs.exp.gtcgtg+cps.obs.exp.gtcgtt+cps.obs.exp.gtctac+cps.obs.exp.gtctat+cps.obs.exp.gtctca+cps.obs.exp.gtctcc+cps.obs.exp.gtctcg+cps.obs.exp.gtctct+cps.obs.exp.gtctgc+cps.obs.exp.gtctgg+cps.obs.exp.gtctgt+cps.obs.exp.gtctta+cps.obs.exp.gtcttc+cps.obs.exp.gtcttg+cps.obs.exp.gtcttt+cps.obs.exp.gtgaaa+cps.obs.exp.gtgaac+cps.obs.exp.gtgaag+cps.obs.exp.gtgaat+cps.obs.exp.gtgaca+cps.obs.exp.gtgacc+cps.obs.exp.gtgacg+cps.obs.exp.gtgact+cps.obs.exp.gtgaga+cps.obs.exp.gtgagc+cps.obs.exp.gtgagg+cps.obs.exp.gtgagt+cps.obs.exp.gtgata+cps.obs.exp.gtgatc+cps.obs.exp.gtgatg+cps.obs.exp.gtgatt+cps.obs.exp.gtgcaa+cps.obs.exp.gtgcac+cps.obs.exp.gtgcag+cps.obs.exp.gtgcat+cps.obs.exp.gtgcca+cps.obs.exp.gtgccc+cps.obs.exp.gtgccg+cps.obs.exp.gtgcct+cps.obs.exp.gtgcga+cps.obs.exp.gtgcgc+cps.obs.exp.gtgcgg+cps.obs.exp.gtgcgt+cps.obs.exp.gtgcta+cps.obs.exp.gtgctc+cps.obs.exp.gtgctg+cps.obs.exp.gtgctt+cps.obs.exp.gtggaa+cps.obs.exp.gtggac+cps.obs.exp.gtggag+cps.obs.exp.gtggat+cps.obs.exp.gtggca+cps.obs.exp.gtggcc+cps.obs.exp.gtggcg+cps.obs.exp.gtggct+cps.obs.exp.gtggga+cps.obs.exp.gtgggc+cps.obs.exp.gtgggg+cps.obs.exp.gtgggt+cps.obs.exp.gtggta+cps.obs.exp.gtggtc+cps.obs.exp.gtggtg+cps.obs.exp.gtggtt+cps.obs.exp.gtgtac+cps.obs.exp.gtgtat+cps.obs.exp.gtgtca+cps.obs.exp.gtgtcc+cps.obs.exp.gtgtcg+cps.obs.exp.gtgtct+cps.obs.exp.gtgtgc+cps.obs.exp.gtgtgg+cps.obs.exp.gtgtgt+cps.obs.exp.gtgtta+cps.obs.exp.gtgttc+cps.obs.exp.gtgttg+cps.obs.exp.gtgttt+cps.obs.exp.gttaaa+cps.obs.exp.gttaac+cps.obs.exp.gttaag+cps.obs.exp.gttaat+cps.obs.exp.gttaca+cps.obs.exp.gttacc+cps.obs.exp.gttacg+cps.obs.exp.gttact+cps.obs.exp.gttaga+cps.obs.exp.gttagc+cps.obs.exp.gttagg+cps.obs.exp.gttagt+cps.obs.exp.gttata+cps.obs.exp.gttatc+cps.obs.exp.gttatg+cps.obs.exp.gttatt+cps.obs.exp.gttcaa+cps.obs.exp.gttcac+cps.obs.exp.gttcag+cps.obs.exp.gttcat+cps.obs.exp.gttcca+cps.obs.exp.gttccc+cps.obs.exp.gttccg+cps.obs.exp.gttcct+cps.obs.exp.gttcga+cps.obs.exp.gttcgc+cps.obs.exp.gttcgg+cps.obs.exp.gttcgt+cps.obs.exp.gttcta+cps.obs.exp.gttctc+cps.obs.exp.gttctg+cps.obs.exp.gttctt+cps.obs.exp.gttgaa+cps.obs.exp.gttgac+cps.obs.exp.gttgag+cps.obs.exp.gttgat+cps.obs.exp.gttgca+cps.obs.exp.gttgcc+cps.obs.exp.gttgcg+cps.obs.exp.gttgct+cps.obs.exp.gttgga+cps.obs.exp.gttggc+cps.obs.exp.gttggg+cps.obs.exp.gttggt+cps.obs.exp.gttgta+cps.obs.exp.gttgtc+cps.obs.exp.gttgtg+cps.obs.exp.gttgtt+cps.obs.exp.gtttac+cps.obs.exp.gtttat+cps.obs.exp.gtttca+cps.obs.exp.gtttcc+cps.obs.exp.gtttcg+cps.obs.exp.gtttct+cps.obs.exp.gtttgc+cps.obs.exp.gtttgg+cps.obs.exp.gtttgt+cps.obs.exp.gtttta+cps.obs.exp.gttttc+cps.obs.exp.gttttg+cps.obs.exp.gttttt+cps.obs.exp.tacaaa+cps.obs.exp.tacaac+cps.obs.exp.tacaag+cps.obs.exp.tacaat+cps.obs.exp.tacaca+cps.obs.exp.tacacc+cps.obs.exp.tacacg+cps.obs.exp.tacact+cps.obs.exp.tacaga+cps.obs.exp.tacagc+cps.obs.exp.tacagg+cps.obs.exp.tacagt+cps.obs.exp.tacata+cps.obs.exp.tacatc+cps.obs.exp.tacatg+cps.obs.exp.tacatt+cps.obs.exp.taccaa+cps.obs.exp.taccac+cps.obs.exp.taccag+cps.obs.exp.taccat+cps.obs.exp.taccca+cps.obs.exp.tacccc+cps.obs.exp.tacccg+cps.obs.exp.taccct+cps.obs.exp.taccga+cps.obs.exp.taccgc+cps.obs.exp.taccgg+cps.obs.exp.taccgt+cps.obs.exp.taccta+cps.obs.exp.tacctc+cps.obs.exp.tacctg+cps.obs.exp.tacctt+cps.obs.exp.tacgaa+cps.obs.exp.tacgac+cps.obs.exp.tacgag+cps.obs.exp.tacgat+cps.obs.exp.tacgca+cps.obs.exp.tacgcc+cps.obs.exp.tacgcg+cps.obs.exp.tacgct+cps.obs.exp.tacgga+cps.obs.exp.tacggc+cps.obs.exp.tacggg+cps.obs.exp.tacggt+cps.obs.exp.tacgta+cps.obs.exp.tacgtc+cps.obs.exp.tacgtg+cps.obs.exp.tacgtt+cps.obs.exp.tactac+cps.obs.exp.tactat+cps.obs.exp.tactca+cps.obs.exp.tactcc+cps.obs.exp.tactcg+cps.obs.exp.tactct+cps.obs.exp.tactgc+cps.obs.exp.tactgg+cps.obs.exp.tactgt+cps.obs.exp.tactta+cps.obs.exp.tacttc+cps.obs.exp.tacttg+cps.obs.exp.tacttt+cps.obs.exp.tataaa+cps.obs.exp.tataac+cps.obs.exp.tataag+cps.obs.exp.tataat+cps.obs.exp.tataca+cps.obs.exp.tatacc+cps.obs.exp.tatacg+cps.obs.exp.tatact+cps.obs.exp.tataga+cps.obs.exp.tatagc+cps.obs.exp.tatagg+cps.obs.exp.tatagt+cps.obs.exp.tatata+cps.obs.exp.tatatc+cps.obs.exp.tatatg+cps.obs.exp.tatatt+cps.obs.exp.tatcaa+cps.obs.exp.tatcac+cps.obs.exp.tatcag+cps.obs.exp.tatcat+cps.obs.exp.tatcca+cps.obs.exp.tatccc+cps.obs.exp.tatccg+cps.obs.exp.tatcct+cps.obs.exp.tatcga+cps.obs.exp.tatcgc+cps.obs.exp.tatcgg+cps.obs.exp.tatcgt+cps.obs.exp.tatcta+cps.obs.exp.tatctc+cps.obs.exp.tatctg+cps.obs.exp.tatctt+cps.obs.exp.tatgaa+cps.obs.exp.tatgac+cps.obs.exp.tatgag+cps.obs.exp.tatgat+cps.obs.exp.tatgca+cps.obs.exp.tatgcc+cps.obs.exp.tatgcg+cps.obs.exp.tatgct+cps.obs.exp.tatgga+cps.obs.exp.tatggc+cps.obs.exp.tatggg+cps.obs.exp.tatggt+cps.obs.exp.tatgta+cps.obs.exp.tatgtc+cps.obs.exp.tatgtg+cps.obs.exp.tatgtt+cps.obs.exp.tattac+cps.obs.exp.tattat+cps.obs.exp.tattca+cps.obs.exp.tattcc+cps.obs.exp.tattcg+cps.obs.exp.tattct+cps.obs.exp.tattgc+cps.obs.exp.tattgg+cps.obs.exp.tattgt+cps.obs.exp.tattta+cps.obs.exp.tatttc+cps.obs.exp.tatttg+cps.obs.exp.tatttt+cps.obs.exp.tcaaaa+cps.obs.exp.tcaaac+cps.obs.exp.tcaaag+cps.obs.exp.tcaaat+cps.obs.exp.tcaaca+cps.obs.exp.tcaacc+cps.obs.exp.tcaacg+cps.obs.exp.tcaact+cps.obs.exp.tcaaga+cps.obs.exp.tcaagc+cps.obs.exp.tcaagg+cps.obs.exp.tcaagt+cps.obs.exp.tcaata+cps.obs.exp.tcaatc+cps.obs.exp.tcaatg+cps.obs.exp.tcaatt+cps.obs.exp.tcacaa+cps.obs.exp.tcacac+cps.obs.exp.tcacag+cps.obs.exp.tcacat+cps.obs.exp.tcacca+cps.obs.exp.tcaccc+cps.obs.exp.tcaccg+cps.obs.exp.tcacct+cps.obs.exp.tcacga+cps.obs.exp.tcacgc+cps.obs.exp.tcacgg+cps.obs.exp.tcacgt+cps.obs.exp.tcacta+cps.obs.exp.tcactc+cps.obs.exp.tcactg+cps.obs.exp.tcactt+cps.obs.exp.tcagaa+cps.obs.exp.tcagac+cps.obs.exp.tcagag+cps.obs.exp.tcagat+cps.obs.exp.tcagca+cps.obs.exp.tcagcc+cps.obs.exp.tcagcg+cps.obs.exp.tcagct+cps.obs.exp.tcagga+cps.obs.exp.tcaggc+cps.obs.exp.tcaggg+cps.obs.exp.tcaggt+cps.obs.exp.tcagta+cps.obs.exp.tcagtc+cps.obs.exp.tcagtg+cps.obs.exp.tcagtt+cps.obs.exp.tcatac+cps.obs.exp.tcatat+cps.obs.exp.tcatca+cps.obs.exp.tcatcc+cps.obs.exp.tcatcg+cps.obs.exp.tcatct+cps.obs.exp.tcatgc+cps.obs.exp.tcatgg+cps.obs.exp.tcatgt+cps.obs.exp.tcatta+cps.obs.exp.tcattc+cps.obs.exp.tcattg+cps.obs.exp.tcattt+cps.obs.exp.tccaaa+cps.obs.exp.tccaac+cps.obs.exp.tccaag+cps.obs.exp.tccaat+cps.obs.exp.tccaca+cps.obs.exp.tccacc+cps.obs.exp.tccacg+cps.obs.exp.tccact+cps.obs.exp.tccaga+cps.obs.exp.tccagc+cps.obs.exp.tccagg+cps.obs.exp.tccagt+cps.obs.exp.tccata+cps.obs.exp.tccatc+cps.obs.exp.tccatg+cps.obs.exp.tccatt+cps.obs.exp.tcccaa+cps.obs.exp.tcccac+cps.obs.exp.tcccag+cps.obs.exp.tcccat+cps.obs.exp.tcccca+cps.obs.exp.tccccc+cps.obs.exp.tccccg+cps.obs.exp.tcccct+cps.obs.exp.tcccga+cps.obs.exp.tcccgc+cps.obs.exp.tcccgg+cps.obs.exp.tcccgt+cps.obs.exp.tcccta+cps.obs.exp.tccctc+cps.obs.exp.tccctg+cps.obs.exp.tccctt+cps.obs.exp.tccgaa+cps.obs.exp.tccgac+cps.obs.exp.tccgag+cps.obs.exp.tccgat+cps.obs.exp.tccgca+cps.obs.exp.tccgcc+cps.obs.exp.tccgcg+cps.obs.exp.tccgct+cps.obs.exp.tccgga+cps.obs.exp.tccggc+cps.obs.exp.tccggg+cps.obs.exp.tccggt+cps.obs.exp.tccgta+cps.obs.exp.tccgtc+cps.obs.exp.tccgtg+cps.obs.exp.tccgtt+cps.obs.exp.tcctac+cps.obs.exp.tcctat+cps.obs.exp.tcctca+cps.obs.exp.tcctcc+cps.obs.exp.tcctcg+cps.obs.exp.tcctct+cps.obs.exp.tcctgc+cps.obs.exp.tcctgg+cps.obs.exp.tcctgt+cps.obs.exp.tcctta+cps.obs.exp.tccttc+cps.obs.exp.tccttg+cps.obs.exp.tccttt+cps.obs.exp.tcgaaa+cps.obs.exp.tcgaac+cps.obs.exp.tcgaag+cps.obs.exp.tcgaat+cps.obs.exp.tcgaca+cps.obs.exp.tcgacc+cps.obs.exp.tcgacg+cps.obs.exp.tcgact+cps.obs.exp.tcgaga+cps.obs.exp.tcgagc+cps.obs.exp.tcgagg+cps.obs.exp.tcgagt+cps.obs.exp.tcgata+cps.obs.exp.tcgatc+cps.obs.exp.tcgatg+cps.obs.exp.tcgatt+cps.obs.exp.tcgcaa+cps.obs.exp.tcgcac+cps.obs.exp.tcgcag+cps.obs.exp.tcgcat+cps.obs.exp.tcgcca+cps.obs.exp.tcgccc+cps.obs.exp.tcgccg+cps.obs.exp.tcgcct+cps.obs.exp.tcgcga+cps.obs.exp.tcgcgc+cps.obs.exp.tcgcgg+cps.obs.exp.tcgcgt+cps.obs.exp.tcgcta+cps.obs.exp.tcgctc+cps.obs.exp.tcgctg+cps.obs.exp.tcgctt+cps.obs.exp.tcggaa+cps.obs.exp.tcggac+cps.obs.exp.tcggag+cps.obs.exp.tcggat+cps.obs.exp.tcggca+cps.obs.exp.tcggcc+cps.obs.exp.tcggcg+cps.obs.exp.tcggct+cps.obs.exp.tcggga+cps.obs.exp.tcgggc+cps.obs.exp.tcgggg+cps.obs.exp.tcgggt+cps.obs.exp.tcggta+cps.obs.exp.tcggtc+cps.obs.exp.tcggtg+cps.obs.exp.tcggtt+cps.obs.exp.tcgtac+cps.obs.exp.tcgtat+cps.obs.exp.tcgtca+cps.obs.exp.tcgtcc+cps.obs.exp.tcgtcg+cps.obs.exp.tcgtct+cps.obs.exp.tcgtgc+cps.obs.exp.tcgtgg+cps.obs.exp.tcgtgt+cps.obs.exp.tcgtta+cps.obs.exp.tcgttc+cps.obs.exp.tcgttg+cps.obs.exp.tcgttt+cps.obs.exp.tctaaa+cps.obs.exp.tctaac+cps.obs.exp.tctaag+cps.obs.exp.tctaat+cps.obs.exp.tctaca+cps.obs.exp.tctacc+cps.obs.exp.tctacg+cps.obs.exp.tctact+cps.obs.exp.tctaga+cps.obs.exp.tctagc+cps.obs.exp.tctagg+cps.obs.exp.tctagt+cps.obs.exp.tctata+cps.obs.exp.tctatc+cps.obs.exp.tctatg+cps.obs.exp.tctatt+cps.obs.exp.tctcaa+cps.obs.exp.tctcac+cps.obs.exp.tctcag+cps.obs.exp.tctcat+cps.obs.exp.tctcca+cps.obs.exp.tctccc+cps.obs.exp.tctccg+cps.obs.exp.tctcct+cps.obs.exp.tctcga+cps.obs.exp.tctcgc+cps.obs.exp.tctcgg+cps.obs.exp.tctcgt+cps.obs.exp.tctcta+cps.obs.exp.tctctc+cps.obs.exp.tctctg+cps.obs.exp.tctctt+cps.obs.exp.tctgaa+cps.obs.exp.tctgac+cps.obs.exp.tctgag+cps.obs.exp.tctgat+cps.obs.exp.tctgca+cps.obs.exp.tctgcc+cps.obs.exp.tctgcg+cps.obs.exp.tctgct+cps.obs.exp.tctgga+cps.obs.exp.tctggc+cps.obs.exp.tctggg+cps.obs.exp.tctggt+cps.obs.exp.tctgta+cps.obs.exp.tctgtc+cps.obs.exp.tctgtg+cps.obs.exp.tctgtt+cps.obs.exp.tcttac+cps.obs.exp.tcttat+cps.obs.exp.tcttca+cps.obs.exp.tcttcc+cps.obs.exp.tcttcg+cps.obs.exp.tcttct+cps.obs.exp.tcttgc+cps.obs.exp.tcttgg+cps.obs.exp.tcttgt+cps.obs.exp.tcttta+cps.obs.exp.tctttc+cps.obs.exp.tctttg+cps.obs.exp.tctttt+cps.obs.exp.tgcaaa+cps.obs.exp.tgcaac+cps.obs.exp.tgcaag+cps.obs.exp.tgcaat+cps.obs.exp.tgcaca+cps.obs.exp.tgcacc+cps.obs.exp.tgcacg+cps.obs.exp.tgcact+cps.obs.exp.tgcaga+cps.obs.exp.tgcagc+cps.obs.exp.tgcagg+cps.obs.exp.tgcagt+cps.obs.exp.tgcata+cps.obs.exp.tgcatc+cps.obs.exp.tgcatg+cps.obs.exp.tgcatt+cps.obs.exp.tgccaa+cps.obs.exp.tgccac+cps.obs.exp.tgccag+cps.obs.exp.tgccat+cps.obs.exp.tgccca+cps.obs.exp.tgcccc+cps.obs.exp.tgcccg+cps.obs.exp.tgccct+cps.obs.exp.tgccga+cps.obs.exp.tgccgc+cps.obs.exp.tgccgg+cps.obs.exp.tgccgt+cps.obs.exp.tgccta+cps.obs.exp.tgcctc+cps.obs.exp.tgcctg+cps.obs.exp.tgcctt+cps.obs.exp.tgcgaa+cps.obs.exp.tgcgac+cps.obs.exp.tgcgag+cps.obs.exp.tgcgat+cps.obs.exp.tgcgca+cps.obs.exp.tgcgcc+cps.obs.exp.tgcgcg+cps.obs.exp.tgcgct+cps.obs.exp.tgcgga+cps.obs.exp.tgcggc+cps.obs.exp.tgcggg+cps.obs.exp.tgcggt+cps.obs.exp.tgcgta+cps.obs.exp.tgcgtc+cps.obs.exp.tgcgtg+cps.obs.exp.tgcgtt+cps.obs.exp.tgctac+cps.obs.exp.tgctat+cps.obs.exp.tgctca+cps.obs.exp.tgctcc+cps.obs.exp.tgctcg+cps.obs.exp.tgctct+cps.obs.exp.tgctgc+cps.obs.exp.tgctgg+cps.obs.exp.tgctgt+cps.obs.exp.tgctta+cps.obs.exp.tgcttc+cps.obs.exp.tgcttg+cps.obs.exp.tgcttt+cps.obs.exp.tggaaa+cps.obs.exp.tggaac+cps.obs.exp.tggaag+cps.obs.exp.tggaat+cps.obs.exp.tggaca+cps.obs.exp.tggacc+cps.obs.exp.tggacg+cps.obs.exp.tggact+cps.obs.exp.tggaga+cps.obs.exp.tggagc+cps.obs.exp.tggagg+cps.obs.exp.tggagt+cps.obs.exp.tggata+cps.obs.exp.tggatc+cps.obs.exp.tggatg+cps.obs.exp.tggatt+cps.obs.exp.tggcaa+cps.obs.exp.tggcac+cps.obs.exp.tggcag+cps.obs.exp.tggcat+cps.obs.exp.tggcca+cps.obs.exp.tggccc+cps.obs.exp.tggccg+cps.obs.exp.tggcct+cps.obs.exp.tggcga+cps.obs.exp.tggcgc+cps.obs.exp.tggcgg+cps.obs.exp.tggcgt+cps.obs.exp.tggcta+cps.obs.exp.tggctc+cps.obs.exp.tggctg+cps.obs.exp.tggctt+cps.obs.exp.tgggaa+cps.obs.exp.tgggac+cps.obs.exp.tgggag+cps.obs.exp.tgggat+cps.obs.exp.tgggca+cps.obs.exp.tgggcc+cps.obs.exp.tgggcg+cps.obs.exp.tgggct+cps.obs.exp.tgggga+cps.obs.exp.tggggc+cps.obs.exp.tggggg+cps.obs.exp.tggggt+cps.obs.exp.tgggta+cps.obs.exp.tgggtc+cps.obs.exp.tgggtg+cps.obs.exp.tgggtt+cps.obs.exp.tggtac+cps.obs.exp.tggtat+cps.obs.exp.tggtca+cps.obs.exp.tggtcc+cps.obs.exp.tggtcg+cps.obs.exp.tggtct+cps.obs.exp.tggtgc+cps.obs.exp.tggtgg+cps.obs.exp.tggtgt+cps.obs.exp.tggtta+cps.obs.exp.tggttc+cps.obs.exp.tggttg+cps.obs.exp.tggttt+cps.obs.exp.tgtaaa+cps.obs.exp.tgtaac+cps.obs.exp.tgtaag+cps.obs.exp.tgtaat+cps.obs.exp.tgtaca+cps.obs.exp.tgtacc+cps.obs.exp.tgtacg+cps.obs.exp.tgtact+cps.obs.exp.tgtaga+cps.obs.exp.tgtagc+cps.obs.exp.tgtagg+cps.obs.exp.tgtagt+cps.obs.exp.tgtata+cps.obs.exp.tgtatc+cps.obs.exp.tgtatg+cps.obs.exp.tgtatt+cps.obs.exp.tgtcaa+cps.obs.exp.tgtcac+cps.obs.exp.tgtcag+cps.obs.exp.tgtcat+cps.obs.exp.tgtcca+cps.obs.exp.tgtccc+cps.obs.exp.tgtccg+cps.obs.exp.tgtcct+cps.obs.exp.tgtcga+cps.obs.exp.tgtcgc+cps.obs.exp.tgtcgg+cps.obs.exp.tgtcgt+cps.obs.exp.tgtcta+cps.obs.exp.tgtctc+cps.obs.exp.tgtctg+cps.obs.exp.tgtctt+cps.obs.exp.tgtgaa+cps.obs.exp.tgtgac+cps.obs.exp.tgtgag+cps.obs.exp.tgtgat+cps.obs.exp.tgtgca+cps.obs.exp.tgtgcc+cps.obs.exp.tgtgcg+cps.obs.exp.tgtgct+cps.obs.exp.tgtgga+cps.obs.exp.tgtggc+cps.obs.exp.tgtggg+cps.obs.exp.tgtggt+cps.obs.exp.tgtgta+cps.obs.exp.tgtgtc+cps.obs.exp.tgtgtg+cps.obs.exp.tgtgtt+cps.obs.exp.tgttac+cps.obs.exp.tgttat+cps.obs.exp.tgttca+cps.obs.exp.tgttcc+cps.obs.exp.tgttcg+cps.obs.exp.tgttct+cps.obs.exp.tgttgc+cps.obs.exp.tgttgg+cps.obs.exp.tgttgt+cps.obs.exp.tgttta+cps.obs.exp.tgtttc+cps.obs.exp.tgtttg+cps.obs.exp.tgtttt+cps.obs.exp.ttaaaa+cps.obs.exp.ttaaac+cps.obs.exp.ttaaag+cps.obs.exp.ttaaat+cps.obs.exp.ttaaca+cps.obs.exp.ttaacc+cps.obs.exp.ttaacg+cps.obs.exp.ttaact+cps.obs.exp.ttaaga+cps.obs.exp.ttaagc+cps.obs.exp.ttaagg+cps.obs.exp.ttaagt+cps.obs.exp.ttaata+cps.obs.exp.ttaatc+cps.obs.exp.ttaatg+cps.obs.exp.ttaatt+cps.obs.exp.ttacaa+cps.obs.exp.ttacac+cps.obs.exp.ttacag+cps.obs.exp.ttacat+cps.obs.exp.ttacca+cps.obs.exp.ttaccc+cps.obs.exp.ttaccg+cps.obs.exp.ttacct+cps.obs.exp.ttacga+cps.obs.exp.ttacgc+cps.obs.exp.ttacgg+cps.obs.exp.ttacgt+cps.obs.exp.ttacta+cps.obs.exp.ttactc+cps.obs.exp.ttactg+cps.obs.exp.ttactt+cps.obs.exp.ttagaa+cps.obs.exp.ttagac+cps.obs.exp.ttagag+cps.obs.exp.ttagat+cps.obs.exp.ttagca+cps.obs.exp.ttagcc+cps.obs.exp.ttagcg+cps.obs.exp.ttagct+cps.obs.exp.ttagga+cps.obs.exp.ttaggc+cps.obs.exp.ttaggg+cps.obs.exp.ttaggt+cps.obs.exp.ttagta+cps.obs.exp.ttagtc+cps.obs.exp.ttagtg+cps.obs.exp.ttagtt+cps.obs.exp.ttatac+cps.obs.exp.ttatat+cps.obs.exp.ttatca+cps.obs.exp.ttatcc+cps.obs.exp.ttatcg+cps.obs.exp.ttatct+cps.obs.exp.ttatgc+cps.obs.exp.ttatgg+cps.obs.exp.ttatgt+cps.obs.exp.ttatta+cps.obs.exp.ttattc+cps.obs.exp.ttattg+cps.obs.exp.ttattt+cps.obs.exp.ttcaaa+cps.obs.exp.ttcaac+cps.obs.exp.ttcaag+cps.obs.exp.ttcaat+cps.obs.exp.ttcaca+cps.obs.exp.ttcacc+cps.obs.exp.ttcacg+cps.obs.exp.ttcact+cps.obs.exp.ttcaga+cps.obs.exp.ttcagc+cps.obs.exp.ttcagg+cps.obs.exp.ttcagt+cps.obs.exp.ttcata+cps.obs.exp.ttcatc+cps.obs.exp.ttcatg+cps.obs.exp.ttcatt+cps.obs.exp.ttccaa+cps.obs.exp.ttccac+cps.obs.exp.ttccag+cps.obs.exp.ttccat+cps.obs.exp.ttccca+cps.obs.exp.ttcccc+cps.obs.exp.ttcccg+cps.obs.exp.ttccct+cps.obs.exp.ttccga+cps.obs.exp.ttccgc+cps.obs.exp.ttccgg+cps.obs.exp.ttccgt+cps.obs.exp.ttccta+cps.obs.exp.ttcctc+cps.obs.exp.ttcctg+cps.obs.exp.ttcctt+cps.obs.exp.ttcgaa+cps.obs.exp.ttcgac+cps.obs.exp.ttcgag+cps.obs.exp.ttcgat+cps.obs.exp.ttcgca+cps.obs.exp.ttcgcc+cps.obs.exp.ttcgcg+cps.obs.exp.ttcgct+cps.obs.exp.ttcgga+cps.obs.exp.ttcggc+cps.obs.exp.ttcggg+cps.obs.exp.ttcggt+cps.obs.exp.ttcgta+cps.obs.exp.ttcgtc+cps.obs.exp.ttcgtg+cps.obs.exp.ttcgtt+cps.obs.exp.ttctac+cps.obs.exp.ttctat+cps.obs.exp.ttctca+cps.obs.exp.ttctcc+cps.obs.exp.ttctcg+cps.obs.exp.ttctct+cps.obs.exp.ttctgc+cps.obs.exp.ttctgg+cps.obs.exp.ttctgt+cps.obs.exp.ttctta+cps.obs.exp.ttcttc+cps.obs.exp.ttcttg+cps.obs.exp.ttcttt+cps.obs.exp.ttgaaa+cps.obs.exp.ttgaac+cps.obs.exp.ttgaag+cps.obs.exp.ttgaat+cps.obs.exp.ttgaca+cps.obs.exp.ttgacc+cps.obs.exp.ttgacg+cps.obs.exp.ttgact+cps.obs.exp.ttgaga+cps.obs.exp.ttgagc+cps.obs.exp.ttgagg+cps.obs.exp.ttgagt+cps.obs.exp.ttgata+cps.obs.exp.ttgatc+cps.obs.exp.ttgatg+cps.obs.exp.ttgatt+cps.obs.exp.ttgcaa+cps.obs.exp.ttgcac+cps.obs.exp.ttgcag+cps.obs.exp.ttgcat+cps.obs.exp.ttgcca+cps.obs.exp.ttgccc+cps.obs.exp.ttgccg+cps.obs.exp.ttgcct+cps.obs.exp.ttgcga+cps.obs.exp.ttgcgc+cps.obs.exp.ttgcgg+cps.obs.exp.ttgcgt+cps.obs.exp.ttgcta+cps.obs.exp.ttgctc+cps.obs.exp.ttgctg+cps.obs.exp.ttgctt+cps.obs.exp.ttggaa+cps.obs.exp.ttggac+cps.obs.exp.ttggag+cps.obs.exp.ttggat+cps.obs.exp.ttggca+cps.obs.exp.ttggcc+cps.obs.exp.ttggcg+cps.obs.exp.ttggct+cps.obs.exp.ttggga+cps.obs.exp.ttgggc+cps.obs.exp.ttgggg+cps.obs.exp.ttgggt+cps.obs.exp.ttggta+cps.obs.exp.ttggtc+cps.obs.exp.ttggtg+cps.obs.exp.ttggtt+cps.obs.exp.ttgtac+cps.obs.exp.ttgtat+cps.obs.exp.ttgtca+cps.obs.exp.ttgtcc+cps.obs.exp.ttgtcg+cps.obs.exp.ttgtct+cps.obs.exp.ttgtgc+cps.obs.exp.ttgtgg+cps.obs.exp.ttgtgt+cps.obs.exp.ttgtta+cps.obs.exp.ttgttc+cps.obs.exp.ttgttg+cps.obs.exp.ttgttt+cps.obs.exp.tttaaa+cps.obs.exp.tttaac+cps.obs.exp.tttaag+cps.obs.exp.tttaat+cps.obs.exp.tttaca+cps.obs.exp.tttacc+cps.obs.exp.tttacg+cps.obs.exp.tttact+cps.obs.exp.tttaga+cps.obs.exp.tttagc+cps.obs.exp.tttagg+cps.obs.exp.tttagt+cps.obs.exp.tttata+cps.obs.exp.tttatc+cps.obs.exp.tttatg+cps.obs.exp.tttatt+cps.obs.exp.tttcaa+cps.obs.exp.tttcac+cps.obs.exp.tttcag+cps.obs.exp.tttcat+cps.obs.exp.tttcca+cps.obs.exp.tttccc+cps.obs.exp.tttccg+cps.obs.exp.tttcct+cps.obs.exp.tttcga+cps.obs.exp.tttcgc+cps.obs.exp.tttcgg+cps.obs.exp.tttcgt+cps.obs.exp.tttcta+cps.obs.exp.tttctc+cps.obs.exp.tttctg+cps.obs.exp.tttctt+cps.obs.exp.tttgaa+cps.obs.exp.tttgac+cps.obs.exp.tttgag+cps.obs.exp.tttgat+cps.obs.exp.tttgca+cps.obs.exp.tttgcc+cps.obs.exp.tttgcg+cps.obs.exp.tttgct+cps.obs.exp.tttgga+cps.obs.exp.tttggc+cps.obs.exp.tttggg+cps.obs.exp.tttggt+cps.obs.exp.tttgta+cps.obs.exp.tttgtc+cps.obs.exp.tttgtg+cps.obs.exp.tttgtt+cps.obs.exp.ttttac+cps.obs.exp.ttttat+cps.obs.exp.ttttca+cps.obs.exp.ttttcc+cps.obs.exp.ttttcg+cps.obs.exp.ttttct+cps.obs.exp.ttttgc+cps.obs.exp.ttttgg+cps.obs.exp.ttttgt+cps.obs.exp.ttttta+cps.obs.exp.tttttc+cps.obs.exp.tttttg+cps.obs.exp.tttttt

#Read CPS table without Inf
#cps.df.mod<-read.table(file="/home2/dmacedod/test/dn3_6deg_original_mod.txt", header=TRUE)

#Generate sum of all codon pair per sequence
cps.df.mod.sum<-apply(cps.df.final, 1, function(x) sum(x))

#Calculate CPB
cpb.dn31<-cps.df.mod.sum/3693

#Write as a table (CPB)
write.table(cpb.dn31, file="/home2/dmacedod/test/dn3_6deg_cpb_result.txt")

#Graph 
#png("CPB/dn3_6deg_cpb.png")
#hist(cpb.dn31, main="CPB (dn3 6deg)", xlab=NULL)
#abline(v=cpb.dn31[1], col=2,lty=2)
#text (0.65, 250, "WT", cex=0.5, col=2)
#dev.off()


               

         
               
               
