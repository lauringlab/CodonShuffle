library("plyr")
library("seqinr")

#Not counting stop codons


#Nucleotide
#dn231 <- read.fasta(file = "3dn_6deg/P1_3dn_6deg_with_original_seq.fas")                     
n3 <- read.fasta(file = "3N/P1_3n_with_original_seq.fas")                     
#dn23 <- read.fasta(file = "dN23/P1_dn23_with_original_seq.fas")                     
#dn31 <- read.fasta(file = "dN31/P1_dn31_with_original_seq.fas")                     

#AA
#dn231.aa <- read.fasta(file = "3dn_6deg/P1_3dn_6deg_with_original_seq_aa")                     
n3.aa <- read.fasta(file = "3N/P1_3n_with_original_seq_aa.fas")                     
#dn23.aa <- read.fasta(file = "dN23/P1_dn23_with_original_seq_aa.fas")                     
#dn31.aa <- read.fasta(file = "dN31/P1_dn31_with_original_seq_aa.fas")

#Count trinucleotide frequency in the sequence
#lapply(n3, function(x) count(x,3))
#Generate a Table with trinucleotide frequency 
#dn231.trinuc<-ldply(dn231, function(x) count(x,3)) 
n3.trinuc<-ldply(n3, function(x) count(x,3)) 
#dn23.trinuc<-ldply(dn23, function(x) count(x,3)) 
#dn31.trinuc<-ldply(dn31, function(x) count(x,3)) 

#Count frequency of codon pair(two trinucleotide) 
#dn231.codon.pair<-ldply(dn231, function(x) count(x,6)) 
n3.codon.pair<-ldply(n3, function(x) count(x,6)) 
#dn23.codon.pair<-ldply(dn23, function(x) count(x,6)) 
#dn31.codon.pair<-ldply(dn31, function(x) count(x,6))

#Generate a Table with AA frequency 
#dn231.count.aa<-ldply(dn231.aa, function(x) count(x,1, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 
n3.count.aa<-ldply(n3.aa, function(x) count(x,1, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 
#dn23.count.aa<-ldply(dn23.aa, function(x) count(x,1, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 
#dn31.count.aa<-ldply(dn31.aa, function(x) count(x,1, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 

#Generate a Table with diaa frequency 
#dn231.count.diaa<-ldply(dn231.aa, function(x) count(x,2, alphabet = c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v"))) 
n3.count.diaa<-ldply(n3.aa, function(x) count(x,2)) 
#dn23.count.diaa<-ldply(dn23.aa, function(x) count(x,2)) 
#dn31.count.diaa<-ldply(dn31.aa, function(x) count(x,2)) 

#CPS expected
#((F(A)xF(B) expected frequency)/F(X)xF(Y))xF(XY)
cps.aaa.aaa<-((n3.trinuc$aaa*n3.trinuc$aaa)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aaa.aac<-((n3.trinuc$aaa*n3.trinuc$aac)/(n3.count.aa$k*n3.count.aa$n))/n3.count.diaa$kn
cps.aaa.aag<-((n3.trinuc$aaa*n3.trinuc$aag)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aaa.aat<-((n3.trinuc$aaa*n3.trinuc$aat)/(n3.count.aa$k*n3.count.aa$n))/n3.count.diaa$kn

cps.aaa.aca<-((n3.trinuc$aaa*n3.trinuc$aca)/(n3.count.aa$k*n3.count.aa$t))/n3.count.diaa$kt
cps.aaa.acc<-((n3.trinuc$aaa*n3.trinuc$acc)/(n3.count.aa$k*n3.count.aa$t))/n3.count.diaa$kt
cps.aaa.acg<-((n3.trinuc$aaa*n3.trinuc$acg)/(n3.count.aa$k*n3.count.aa$t))/n3.count.diaa$kt
cps.aaa.act<-((n3.trinuc$aaa*n3.trinuc$act)/(n3.count.aa$k*n3.count.aa$t))/n3.count.diaa$kt

cps.aaa.aga<-((n3.trinuc$aaa*n3.trinuc$aga)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aaa.agc<-((n3.trinuc$aaa*n3.trinuc$agc)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks
cps.aaa.agg<-((n3.trinuc$aaa*n3.trinuc$agg)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aaa.agt<-((n3.trinuc$aaa*n3.trinuc$agt)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks

cps.aaa.ata<-((n3.trinuc$aaa*n3.trinuc$ata)/(n3.count.aa$k*n3.count.aa$i))/n3.count.diaa$ki
cps.aaa.atc<-((n3.trinuc$aaa*n3.trinuc$atc)/(n3.count.aa$k*n3.count.aa$i))/n3.count.diaa$ki
cps.aaa.atg<-((n3.trinuc$aaa*n3.trinuc$atg)/(n3.count.aa$k*n3.count.aa$m))/n3.count.diaa$km
cps.aaa.att<-((n3.trinuc$aaa*n3.trinuc$att)/(n3.count.aa$k*n3.count.aa$i))/n3.count.diaa$ki

cps.aaa.caa<-((n3.trinuc$aaa*n3.trinuc$caa)/(n3.count.aa$k*n3.count.aa$q))/n3.count.diaa$kq
cps.aaa.cac<-((n3.trinuc$aaa*n3.trinuc$cac)/(n3.count.aa$k*n3.count.aa$h))/n3.count.diaa$kh
cps.aaa.cag<-((n3.trinuc$aaa*n3.trinuc$cag)/(n3.count.aa$k*n3.count.aa$q))/n3.count.diaa$kq
cps.aaa.cat<-((n3.trinuc$aaa*n3.trinuc$cat)/(n3.count.aa$k*n3.count.aa$h))/n3.count.diaa$kh

cps.aaa.cca<-((n3.trinuc$aaa*n3.trinuc$cca)/(n3.count.aa$k*n3.count.aa$p))/n3.count.diaa$kp
cps.aaa.ccc<-((n3.trinuc$aaa*n3.trinuc$ccc)/(n3.count.aa$k*n3.count.aa$p))/n3.count.diaa$kp
cps.aaa.ccg<-((n3.trinuc$aaa*n3.trinuc$ccg)/(n3.count.aa$k*n3.count.aa$p))/n3.count.diaa$kp
cps.aaa.cct<-((n3.trinuc$aaa*n3.trinuc$cct)/(n3.count.aa$k*n3.count.aa$p))/n3.count.diaa$kp

cps.aaa.cga<-((n3.trinuc$aaa*n3.trinuc$cga)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aaa.cgc<-((n3.trinuc$aaa*n3.trinuc$cgc)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aaa.cgg<-((n3.trinuc$aaa*n3.trinuc$cgg)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aaa.cgt<-((n3.trinuc$aaa*n3.trinuc$cgt)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr

cps.aaa.cta<-((n3.trinuc$aaa*n3.trinuc$cta)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aaa.ctc<-((n3.trinuc$aaa*n3.trinuc$ctc)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aaa.ctg<-((n3.trinuc$aaa*n3.trinuc$ctg)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aaa.ctt<-((n3.trinuc$aaa*n3.trinuc$ctt)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl

cps.aaa.gaa<-((n3.trinuc$aaa*n3.trinuc$gaa)/(n3.count.aa$k*n3.count.aa$e))/n3.count.diaa$ke
cps.aaa.gac<-((n3.trinuc$aaa*n3.trinuc$gac)/(n3.count.aa$k*n3.count.aa$d))/n3.count.diaa$kd
cps.aaa.gag<-((n3.trinuc$aaa*n3.trinuc$gag)/(n3.count.aa$k*n3.count.aa$e))/n3.count.diaa$ke
cps.aaa.gat<-((n3.trinuc$aaa*n3.trinuc$gat)/(n3.count.aa$k*n3.count.aa$d))/n3.count.diaa$kd

cps.aaa.gca<-((n3.trinuc$aaa*n3.trinuc$gca)/(n3.count.aa$k*n3.count.aa$a))/n3.count.diaa$ka
cps.aaa.gcc<-((n3.trinuc$aaa*n3.trinuc$gcc)/(n3.count.aa$k*n3.count.aa$a))/n3.count.diaa$ka
cps.aaa.gcg<-((n3.trinuc$aaa*n3.trinuc$gcg)/(n3.count.aa$k*n3.count.aa$a))/n3.count.diaa$ka
cps.aaa.gct<-((n3.trinuc$aaa*n3.trinuc$gct)/(n3.count.aa$k*n3.count.aa$a))/n3.count.diaa$ka

cps.aaa.gga<-((n3.trinuc$aaa*n3.trinuc$gga)/(n3.count.aa$k*n3.count.aa$g))/n3.count.diaa$kg
cps.aaa.ggc<-((n3.trinuc$aaa*n3.trinuc$ggc)/(n3.count.aa$k*n3.count.aa$g))/n3.count.diaa$kg
cps.aaa.ggg<-((n3.trinuc$aaa*n3.trinuc$ggg)/(n3.count.aa$k*n3.count.aa$g))/n3.count.diaa$kg
cps.aaa.ggt<-((n3.trinuc$aaa*n3.trinuc$ggt)/(n3.count.aa$k*n3.count.aa$g))/n3.count.diaa$kg

cps.aaa.gta<-((n3.trinuc$aaa*n3.trinuc$gta)/(n3.count.aa$k*n3.count.aa$v))/n3.count.diaa$kv
cps.aaa.gtc<-((n3.trinuc$aaa*n3.trinuc$gtc)/(n3.count.aa$k*n3.count.aa$v))/n3.count.diaa$kv
cps.aaa.gtg<-((n3.trinuc$aaa*n3.trinuc$gtg)/(n3.count.aa$k*n3.count.aa$v))/n3.count.diaa$kv
cps.aaa.gtt<-((n3.trinuc$aaa*n3.trinuc$gtt)/(n3.count.aa$k*n3.count.aa$v))/n3.count.diaa$kv

#stop codon
#cps.aaa.taa<-((n3.trinuc$aaa*n3.trinuc$taa)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aaa.tac<-((n3.trinuc$aaa*n3.trinuc$tac)/(n3.count.aa$k*n3.count.aa$y))/n3.count.diaa$ky
#stop codon
#cps.aaa.tag<-((n3.trinuc$aaa*n3.trinuc$tag)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aaa.tat<-((n3.trinuc$aaa*n3.trinuc$tat)/(n3.count.aa$k*n3.count.aa$y))/n3.count.diaa$ky

cps.aaa.tca<-((n3.trinuc$aaa*n3.trinuc$tca)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks
cps.aaa.tcc<-((n3.trinuc$aaa*n3.trinuc$tcc)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks
cps.aaa.tcg<-((n3.trinuc$aaa*n3.trinuc$tcg)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks
cps.aaa.tct<-((n3.trinuc$aaa*n3.trinuc$tct)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks

#Stop codon
#cps.aaa.tga<-((n3.trinuc$aaa*n3.trinuc$tga)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aaa.tgc<-((n3.trinuc$aaa*n3.trinuc$tgc)/(n3.count.aa$k*n3.count.aa$c))/n3.count.diaa$kc
cps.aaa.tgg<-((n3.trinuc$aaa*n3.trinuc$tgg)/(n3.count.aa$k*n3.count.aa$w))/n3.count.diaa$kw
cps.aaa.tgt<-((n3.trinuc$aaa*n3.trinuc$tgt)/(n3.count.aa$k*n3.count.aa$c))/n3.count.diaa$kc

cps.aaa.tta<-((n3.trinuc$aaa*n3.trinuc$tta)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aaa.ttc<-((n3.trinuc$aaa*n3.trinuc$ttc)/(n3.count.aa$k*n3.count.aa$f))/n3.count.diaa$kf
cps.aaa.ttg<-((n3.trinuc$aaa*n3.trinuc$ttg)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aaa.ttt<-((n3.trinuc$aaa*n3.trinuc$ttt)/(n3.count.aa$k*n3.count.aa$f))/n3.count.diaa$kf




cps.aac.aaa<-((n3.trinuc$aac*n3.trinuc$aaa)/(n3.count.aa$n*n3.count.aa$k))/n3.count.diaa$nk
cps.aac.aac<-((n3.trinuc$aac*n3.trinuc$aac)/(n3.count.aa$n*n3.count.aa$n))/n3.count.diaa$nn
cps.aac.aag<-((n3.trinuc$aac*n3.trinuc$aag)/(n3.count.aa$n*n3.count.aa$k))/n3.count.diaa$nk
cps.aac.aat<-((n3.trinuc$aac*n3.trinuc$aat)/(n3.count.aa$n*n3.count.aa$n))/n3.count.diaa$nn

cps.aac.aca<-((n3.trinuc$aac*n3.trinuc$aca)/(n3.count.aa$n*n3.count.aa$t))/n3.count.diaa$nt
cps.aac.acc<-((n3.trinuc$aac*n3.trinuc$acc)/(n3.count.aa$n*n3.count.aa$t))/n3.count.diaa$nt
cps.aac.acg<-((n3.trinuc$aac*n3.trinuc$acg)/(n3.count.aa$n*n3.count.aa$t))/n3.count.diaa$nt
cps.aac.act<-((n3.trinuc$aac*n3.trinuc$act)/(n3.count.aa$n*n3.count.aa$t))/n3.count.diaa$nt

cps.aac.aga<-((n3.trinuc$aac*n3.trinuc$aga)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aac.agc<-((n3.trinuc$aac*n3.trinuc$agc)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns
cps.aac.agg<-((n3.trinuc$aac*n3.trinuc$agg)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aac.agt<-((n3.trinuc$aac*n3.trinuc$agt)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns

cps.aac.ata<-((n3.trinuc$aac*n3.trinuc$ata)/(n3.count.aa$n*n3.count.aa$i))/n3.count.diaa$ni
cps.aac.atc<-((n3.trinuc$aac*n3.trinuc$atc)/(n3.count.aa$n*n3.count.aa$i))/n3.count.diaa$ni
cps.aac.atg<-((n3.trinuc$aac*n3.trinuc$atg)/(n3.count.aa$n*n3.count.aa$m))/n3.count.diaa$nm
cps.aac.att<-((n3.trinuc$aac*n3.trinuc$att)/(n3.count.aa$n*n3.count.aa$i))/n3.count.diaa$ni

cps.aac.caa<-((n3.trinuc$aac*n3.trinuc$caa)/(n3.count.aa$n*n3.count.aa$q))/n3.count.diaa$nq
cps.aac.cac<-((n3.trinuc$aac*n3.trinuc$cac)/(n3.count.aa$n*n3.count.aa$h))/n3.count.diaa$nh
cps.aac.cag<-((n3.trinuc$aac*n3.trinuc$cag)/(n3.count.aa$n*n3.count.aa$q))/n3.count.diaa$nq
cps.aac.cat<-((n3.trinuc$aac*n3.trinuc$cat)/(n3.count.aa$n*n3.count.aa$h))/n3.count.diaa$nh

cps.aac.cca<-((n3.trinuc$aac*n3.trinuc$cca)/(n3.count.aa$n*n3.count.aa$p))/n3.count.diaa$np
cps.aac.ccc<-((n3.trinuc$aac*n3.trinuc$ccc)/(n3.count.aa$n*n3.count.aa$p))/n3.count.diaa$np
cps.aac.ccg<-((n3.trinuc$aac*n3.trinuc$ccg)/(n3.count.aa$n*n3.count.aa$p))/n3.count.diaa$np
cps.aac.cct<-((n3.trinuc$aac*n3.trinuc$cct)/(n3.count.aa$n*n3.count.aa$p))/n3.count.diaa$np

cps.aac.cga<-((n3.trinuc$aac*n3.trinuc$cga)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aac.cgc<-((n3.trinuc$aac*n3.trinuc$cgc)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aac.cgg<-((n3.trinuc$aac*n3.trinuc$cgg)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aac.cgt<-((n3.trinuc$aac*n3.trinuc$cgt)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr

cps.aac.cta<-((n3.trinuc$aac*n3.trinuc$cta)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aac.ctc<-((n3.trinuc$aac*n3.trinuc$ctc)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aac.ctg<-((n3.trinuc$aac*n3.trinuc$ctg)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aac.ctt<-((n3.trinuc$aac*n3.trinuc$ctt)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl

cps.aac.gaa<-((n3.trinuc$aac*n3.trinuc$gaa)/(n3.count.aa$n*n3.count.aa$e))/n3.count.diaa$ne
cps.aac.gac<-((n3.trinuc$aac*n3.trinuc$gac)/(n3.count.aa$n*n3.count.aa$d))/n3.count.diaa$nd
cps.aac.gag<-((n3.trinuc$aac*n3.trinuc$gag)/(n3.count.aa$n*n3.count.aa$e))/n3.count.diaa$ne
cps.aac.gat<-((n3.trinuc$aac*n3.trinuc$gat)/(n3.count.aa$n*n3.count.aa$d))/n3.count.diaa$nd

cps.aac.gca<-((n3.trinuc$aac*n3.trinuc$gca)/(n3.count.aa$n*n3.count.aa$a))/n3.count.diaa$na
cps.aac.gcc<-((n3.trinuc$aac*n3.trinuc$gcc)/(n3.count.aa$n*n3.count.aa$a))/n3.count.diaa$na
cps.aac.gcg<-((n3.trinuc$aac*n3.trinuc$gcg)/(n3.count.aa$n*n3.count.aa$a))/n3.count.diaa$na
cps.aac.gct<-((n3.trinuc$aac*n3.trinuc$gct)/(n3.count.aa$n*n3.count.aa$a))/n3.count.diaa$na

cps.aac.gga<-((n3.trinuc$aac*n3.trinuc$gga)/(n3.count.aa$n*n3.count.aa$g))/n3.count.diaa$ng
cps.aac.ggc<-((n3.trinuc$aac*n3.trinuc$ggc)/(n3.count.aa$n*n3.count.aa$g))/n3.count.diaa$ng
cps.aac.ggg<-((n3.trinuc$aac*n3.trinuc$ggg)/(n3.count.aa$n*n3.count.aa$g))/n3.count.diaa$ng
cps.aac.ggt<-((n3.trinuc$aac*n3.trinuc$ggt)/(n3.count.aa$n*n3.count.aa$g))/n3.count.diaa$ng

cps.aac.gta<-((n3.trinuc$aac*n3.trinuc$gta)/(n3.count.aa$n*n3.count.aa$v))/n3.count.diaa$nv
cps.aac.gtc<-((n3.trinuc$aac*n3.trinuc$gtc)/(n3.count.aa$n*n3.count.aa$v))/n3.count.diaa$nv
cps.aac.gtg<-((n3.trinuc$aac*n3.trinuc$gtg)/(n3.count.aa$n*n3.count.aa$v))/n3.count.diaa$nv
cps.aac.gtt<-((n3.trinuc$aac*n3.trinuc$gtt)/(n3.count.aa$n*n3.count.aa$v))/n3.count.diaa$nv

#Stop codon
#cps.aac.taa<-((n3.trinuc$aac*n3.trinuc$taa)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aac.tac<-((n3.trinuc$aac*n3.trinuc$tac)/(n3.count.aa$n*n3.count.aa$y))/n3.count.diaa$ny
#Stop codon
#cps.aac.tag<-((n3.trinuc$aac*n3.trinuc$tag)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aac.tat<-((n3.trinuc$aac*n3.trinuc$tat)/(n3.count.aa$n*n3.count.aa$y))/n3.count.diaa$ny

cps.aac.tca<-((n3.trinuc$aac*n3.trinuc$tca)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns
cps.aac.tcc<-((n3.trinuc$aac*n3.trinuc$tcc)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns
cps.aac.tcg<-((n3.trinuc$aac*n3.trinuc$tcg)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns
cps.aac.tct<-((n3.trinuc$aac*n3.trinuc$tct)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns

#Stop codon
#cps.aac.tga<-((n3.trinuc$aac*n3.trinuc$tga)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aac.tgc<-((n3.trinuc$aac*n3.trinuc$tgc)/(n3.count.aa$n*n3.count.aa$c))/n3.count.diaa$nc
cps.aac.tgg<-((n3.trinuc$aac*n3.trinuc$tgg)/(n3.count.aa$n*n3.count.aa$w))/n3.count.diaa$nw
cps.aac.tgt<-((n3.trinuc$aac*n3.trinuc$tgt)/(n3.count.aa$n*n3.count.aa$c))/n3.count.diaa$nc

cps.aac.tta<-((n3.trinuc$aac*n3.trinuc$tta)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aac.ttc<-((n3.trinuc$aac*n3.trinuc$ttc)/(n3.count.aa$n*n3.count.aa$f))/n3.count.diaa$nf
cps.aac.ttg<-((n3.trinuc$aac*n3.trinuc$ttg)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aac.ttt<-((n3.trinuc$aac*n3.trinuc$ttt)/(n3.count.aa$n*n3.count.aa$f))/n3.count.diaa$nf





cps.aag.aaa<-((n3.trinuc$aag*n3.trinuc$aaa)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aag.aac<-((n3.trinuc$aag*n3.trinuc$aac)/(n3.count.aa$k*n3.count.aa$n))/n3.count.diaa$kn
cps.aag.aag<-((n3.trinuc$aag*n3.trinuc$aag)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aag.aat<-((n3.trinuc$aag*n3.trinuc$aat)/(n3.count.aa$k*n3.count.aa$n))/n3.count.diaa$kn

cps.aag.aca<-((n3.trinuc$aag*n3.trinuc$aca)/(n3.count.aa$k*n3.count.aa$t))/n3.count.diaa$kt
cps.aag.acc<-((n3.trinuc$aag*n3.trinuc$acc)/(n3.count.aa$k*n3.count.aa$t))/n3.count.diaa$kt
cps.aag.acg<-((n3.trinuc$aag*n3.trinuc$acg)/(n3.count.aa$k*n3.count.aa$t))/n3.count.diaa$kt
cps.aag.act<-((n3.trinuc$aag*n3.trinuc$act)/(n3.count.aa$k*n3.count.aa$t))/n3.count.diaa$kt

cps.aag.aga<-((n3.trinuc$aag*n3.trinuc$aga)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aag.agc<-((n3.trinuc$aag*n3.trinuc$agc)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks
cps.aag.agg<-((n3.trinuc$aag*n3.trinuc$agg)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aag.agt<-((n3.trinuc$aag*n3.trinuc$agt)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks

cps.aag.ata<-((n3.trinuc$aag*n3.trinuc$ata)/(n3.count.aa$k*n3.count.aa$i))/n3.count.diaa$ki
cps.aag.atc<-((n3.trinuc$aag*n3.trinuc$atc)/(n3.count.aa$k*n3.count.aa$i))/n3.count.diaa$ki
cps.aag.atg<-((n3.trinuc$aag*n3.trinuc$atg)/(n3.count.aa$k*n3.count.aa$m))/n3.count.diaa$km
cps.aag.att<-((n3.trinuc$aag*n3.trinuc$att)/(n3.count.aa$k*n3.count.aa$i))/n3.count.diaa$ki

cps.aag.caa<-((n3.trinuc$aag*n3.trinuc$caa)/(n3.count.aa$k*n3.count.aa$q))/n3.count.diaa$kq
cps.aag.cac<-((n3.trinuc$aag*n3.trinuc$cac)/(n3.count.aa$k*n3.count.aa$h))/n3.count.diaa$kh
cps.aag.cag<-((n3.trinuc$aag*n3.trinuc$cag)/(n3.count.aa$k*n3.count.aa$q))/n3.count.diaa$kq
cps.aag.cat<-((n3.trinuc$aag*n3.trinuc$cat)/(n3.count.aa$k*n3.count.aa$h))/n3.count.diaa$kh

cps.aag.cca<-((n3.trinuc$aag*n3.trinuc$cca)/(n3.count.aa$k*n3.count.aa$p))/n3.count.diaa$kp
cps.aag.ccc<-((n3.trinuc$aag*n3.trinuc$ccc)/(n3.count.aa$k*n3.count.aa$p))/n3.count.diaa$kp
cps.aag.ccg<-((n3.trinuc$aag*n3.trinuc$ccg)/(n3.count.aa$k*n3.count.aa$p))/n3.count.diaa$kp
cps.aag.cct<-((n3.trinuc$aag*n3.trinuc$cct)/(n3.count.aa$k*n3.count.aa$p))/n3.count.diaa$kp

cps.aag.cga<-((n3.trinuc$aag*n3.trinuc$cga)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aag.cgc<-((n3.trinuc$aag*n3.trinuc$cgc)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aag.cgg<-((n3.trinuc$aag*n3.trinuc$cgg)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr
cps.aag.cgt<-((n3.trinuc$aag*n3.trinuc$cgt)/(n3.count.aa$k*n3.count.aa$r))/n3.count.diaa$kr

cps.aag.cta<-((n3.trinuc$aag*n3.trinuc$cta)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aag.ctc<-((n3.trinuc$aag*n3.trinuc$ctc)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aag.ctg<-((n3.trinuc$aag*n3.trinuc$ctg)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aag.ctt<-((n3.trinuc$aag*n3.trinuc$ctt)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl

cps.aag.gaa<-((n3.trinuc$aag*n3.trinuc$gaa)/(n3.count.aa$k*n3.count.aa$e))/n3.count.diaa$ke
cps.aag.gac<-((n3.trinuc$aag*n3.trinuc$gac)/(n3.count.aa$k*n3.count.aa$d))/n3.count.diaa$kd
cps.aag.gag<-((n3.trinuc$aag*n3.trinuc$gag)/(n3.count.aa$k*n3.count.aa$e))/n3.count.diaa$ke
cps.aag.gat<-((n3.trinuc$aag*n3.trinuc$gat)/(n3.count.aa$k*n3.count.aa$d))/n3.count.diaa$kd

cps.aag.gca<-((n3.trinuc$aag*n3.trinuc$gca)/(n3.count.aa$k*n3.count.aa$a))/n3.count.diaa$ka
cps.aag.gcc<-((n3.trinuc$aag*n3.trinuc$gcc)/(n3.count.aa$k*n3.count.aa$a))/n3.count.diaa$ka
cps.aag.gcg<-((n3.trinuc$aag*n3.trinuc$gcg)/(n3.count.aa$k*n3.count.aa$a))/n3.count.diaa$ka
cps.aag.gct<-((n3.trinuc$aag*n3.trinuc$gct)/(n3.count.aa$k*n3.count.aa$a))/n3.count.diaa$ka

cps.aag.gga<-((n3.trinuc$aag*n3.trinuc$gga)/(n3.count.aa$k*n3.count.aa$g))/n3.count.diaa$kg
cps.aag.ggc<-((n3.trinuc$aag*n3.trinuc$ggc)/(n3.count.aa$k*n3.count.aa$g))/n3.count.diaa$kg
cps.aag.ggg<-((n3.trinuc$aag*n3.trinuc$ggg)/(n3.count.aa$k*n3.count.aa$g))/n3.count.diaa$kg
cps.aag.ggt<-((n3.trinuc$aag*n3.trinuc$ggt)/(n3.count.aa$k*n3.count.aa$g))/n3.count.diaa$kg

cps.aag.gta<-((n3.trinuc$aag*n3.trinuc$gta)/(n3.count.aa$k*n3.count.aa$v))/n3.count.diaa$kv
cps.aag.gtc<-((n3.trinuc$aag*n3.trinuc$gtc)/(n3.count.aa$k*n3.count.aa$v))/n3.count.diaa$kv
cps.aag.gtg<-((n3.trinuc$aag*n3.trinuc$gtg)/(n3.count.aa$k*n3.count.aa$v))/n3.count.diaa$kv
cps.aag.gtt<-((n3.trinuc$aag*n3.trinuc$gtt)/(n3.count.aa$k*n3.count.aa$v))/n3.count.diaa$kv

#Stop codon
#cps.aag.taa<-((n3.trinuc$aag*n3.trinuc$taa)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aag.tac<-((n3.trinuc$aag*n3.trinuc$tac)/(n3.count.aa$k*n3.count.aa$y))/n3.count.diaa$ky
#Stop codon
#cps.aag.tag<-((n3.trinuc$aag*n3.trinuc$tag)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aag.tat<-((n3.trinuc$aag*n3.trinuc$tat)/(n3.count.aa$k*n3.count.aa$y))/n3.count.diaa$ky

cps.aag.tca<-((n3.trinuc$aag*n3.trinuc$tca)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks
cps.aag.tcc<-((n3.trinuc$aag*n3.trinuc$tcc)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks
cps.aag.tcg<-((n3.trinuc$aag*n3.trinuc$tcg)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks
cps.aag.tct<-((n3.trinuc$aag*n3.trinuc$tct)/(n3.count.aa$k*n3.count.aa$s))/n3.count.diaa$ks

#Stop codon
#cps.aag.tga<-((n3.trinuc$aag*n3.trinuc$tga)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aag.tgc<-((n3.trinuc$aag*n3.trinuc$tgc)/(n3.count.aa$k*n3.count.aa$c))/n3.count.diaa$kc
cps.aag.tgg<-((n3.trinuc$aag*n3.trinuc$tgg)/(n3.count.aa$k*n3.count.aa$w))/n3.count.diaa$kw
cps.aag.tgt<-((n3.trinuc$aag*n3.trinuc$tgt)/(n3.count.aa$k*n3.count.aa$c))/n3.count.diaa$kc

cps.aag.tta<-((n3.trinuc$aag*n3.trinuc$tta)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aag.ttc<-((n3.trinuc$aag*n3.trinuc$ttc)/(n3.count.aa$k*n3.count.aa$f))/n3.count.diaa$kf
cps.aag.ttg<-((n3.trinuc$aag*n3.trinuc$ttg)/(n3.count.aa$k*n3.count.aa$l))/n3.count.diaa$kl
cps.aag.ttt<-((n3.trinuc$aag*n3.trinuc$ttt)/(n3.count.aa$k*n3.count.aa$f))/n3.count.diaa$kf








cps.aat.aaa<-((n3.trinuc$aat*n3.trinuc$aaa)/(n3.count.aa$n*n3.count.aa$k))/n3.count.diaa$nk
cps.aat.aac<-((n3.trinuc$aat*n3.trinuc$aac)/(n3.count.aa$n*n3.count.aa$n))/n3.count.diaa$nn
cps.aat.aag<-((n3.trinuc$aat*n3.trinuc$aag)/(n3.count.aa$n*n3.count.aa$k))/n3.count.diaa$nk
cps.aat.aat<-((n3.trinuc$aat*n3.trinuc$aat)/(n3.count.aa$n*n3.count.aa$n))/n3.count.diaa$nn

cps.aat.aca<-((n3.trinuc$aat*n3.trinuc$aca)/(n3.count.aa$n*n3.count.aa$t))/n3.count.diaa$nt
cps.aat.acc<-((n3.trinuc$aat*n3.trinuc$acc)/(n3.count.aa$n*n3.count.aa$t))/n3.count.diaa$nt
cps.aat.acg<-((n3.trinuc$aat*n3.trinuc$acg)/(n3.count.aa$n*n3.count.aa$t))/n3.count.diaa$nt
cps.aat.act<-((n3.trinuc$aat*n3.trinuc$act)/(n3.count.aa$n*n3.count.aa$t))/n3.count.diaa$nt

cps.aat.aga<-((n3.trinuc$aat*n3.trinuc$aga)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aat.agc<-((n3.trinuc$aat*n3.trinuc$agc)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns
cps.aat.agg<-((n3.trinuc$aat*n3.trinuc$agg)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aat.agt<-((n3.trinuc$aat*n3.trinuc$agt)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns

cps.aat.ata<-((n3.trinuc$aat*n3.trinuc$ata)/(n3.count.aa$n*n3.count.aa$i))/n3.count.diaa$ni
cps.aat.atc<-((n3.trinuc$aat*n3.trinuc$atc)/(n3.count.aa$n*n3.count.aa$i))/n3.count.diaa$ni
cps.aat.atg<-((n3.trinuc$aat*n3.trinuc$atg)/(n3.count.aa$n*n3.count.aa$m))/n3.count.diaa$nm
cps.aat.att<-((n3.trinuc$aat*n3.trinuc$att)/(n3.count.aa$n*n3.count.aa$i))/n3.count.diaa$ni

cps.aat.caa<-((n3.trinuc$aat*n3.trinuc$caa)/(n3.count.aa$n*n3.count.aa$q))/n3.count.diaa$nq
cps.aat.cac<-((n3.trinuc$aat*n3.trinuc$cac)/(n3.count.aa$n*n3.count.aa$h))/n3.count.diaa$nh
cps.aat.cag<-((n3.trinuc$aat*n3.trinuc$cag)/(n3.count.aa$n*n3.count.aa$q))/n3.count.diaa$nq
cps.aat.cat<-((n3.trinuc$aat*n3.trinuc$cat)/(n3.count.aa$n*n3.count.aa$h))/n3.count.diaa$nh

cps.aat.cca<-((n3.trinuc$aat*n3.trinuc$cca)/(n3.count.aa$n*n3.count.aa$p))/n3.count.diaa$np
cps.aat.ccc<-((n3.trinuc$aat*n3.trinuc$ccc)/(n3.count.aa$n*n3.count.aa$p))/n3.count.diaa$np
cps.aat.ccg<-((n3.trinuc$aat*n3.trinuc$ccg)/(n3.count.aa$n*n3.count.aa$p))/n3.count.diaa$np
cps.aat.cct<-((n3.trinuc$aat*n3.trinuc$cct)/(n3.count.aa$n*n3.count.aa$p))/n3.count.diaa$np

cps.aat.cga<-((n3.trinuc$aat*n3.trinuc$cga)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aat.cgc<-((n3.trinuc$aat*n3.trinuc$cgc)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aat.cgg<-((n3.trinuc$aat*n3.trinuc$cgg)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr
cps.aat.cgt<-((n3.trinuc$aat*n3.trinuc$cgt)/(n3.count.aa$n*n3.count.aa$r))/n3.count.diaa$nr

cps.aat.cta<-((n3.trinuc$aat*n3.trinuc$cta)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aat.ctc<-((n3.trinuc$aat*n3.trinuc$ctc)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aat.ctg<-((n3.trinuc$aat*n3.trinuc$ctg)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aat.ctt<-((n3.trinuc$aat*n3.trinuc$ctt)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl

cps.aat.gaa<-((n3.trinuc$aat*n3.trinuc$gaa)/(n3.count.aa$n*n3.count.aa$e))/n3.count.diaa$ne
cps.aat.gac<-((n3.trinuc$aat*n3.trinuc$gac)/(n3.count.aa$n*n3.count.aa$d))/n3.count.diaa$nd
cps.aat.gag<-((n3.trinuc$aat*n3.trinuc$gag)/(n3.count.aa$n*n3.count.aa$e))/n3.count.diaa$ne
cps.aat.gat<-((n3.trinuc$aat*n3.trinuc$gat)/(n3.count.aa$n*n3.count.aa$d))/n3.count.diaa$nd

cps.aat.gca<-((n3.trinuc$aat*n3.trinuc$gca)/(n3.count.aa$n*n3.count.aa$a))/n3.count.diaa$na
cps.aat.gcc<-((n3.trinuc$aat*n3.trinuc$gcc)/(n3.count.aa$n*n3.count.aa$a))/n3.count.diaa$na
cps.aat.gcg<-((n3.trinuc$aat*n3.trinuc$gcg)/(n3.count.aa$n*n3.count.aa$a))/n3.count.diaa$na
cps.aat.gct<-((n3.trinuc$aat*n3.trinuc$gct)/(n3.count.aa$n*n3.count.aa$a))/n3.count.diaa$na

cps.aat.gga<-((n3.trinuc$aat*n3.trinuc$gga)/(n3.count.aa$n*n3.count.aa$g))/n3.count.diaa$ng
cps.aat.ggc<-((n3.trinuc$aat*n3.trinuc$ggc)/(n3.count.aa$n*n3.count.aa$g))/n3.count.diaa$ng
cps.aat.ggg<-((n3.trinuc$aat*n3.trinuc$ggg)/(n3.count.aa$n*n3.count.aa$g))/n3.count.diaa$ng
cps.aat.ggt<-((n3.trinuc$aat*n3.trinuc$ggt)/(n3.count.aa$n*n3.count.aa$g))/n3.count.diaa$ng

cps.aat.gta<-((n3.trinuc$aat*n3.trinuc$gta)/(n3.count.aa$n*n3.count.aa$v))/n3.count.diaa$nv
cps.aat.gtc<-((n3.trinuc$aat*n3.trinuc$gtc)/(n3.count.aa$n*n3.count.aa$v))/n3.count.diaa$nv
cps.aat.gtg<-((n3.trinuc$aat*n3.trinuc$gtg)/(n3.count.aa$n*n3.count.aa$v))/n3.count.diaa$nv
cps.aat.gtt<-((n3.trinuc$aat*n3.trinuc$gtt)/(n3.count.aa$n*n3.count.aa$v))/n3.count.diaa$nv

#Stop codon
#cps.aat.taa<-((n3.trinuc$aat*n3.trinuc$taa)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aat.tac<-((n3.trinuc$aat*n3.trinuc$tac)/(n3.count.aa$n*n3.count.aa$y))/n3.count.diaa$ny
#Stop codon
#cps.aat.tag<-((n3.trinuc$aat*n3.trinuc$tag)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aat.tat<-((n3.trinuc$aat*n3.trinuc$tat)/(n3.count.aa$n*n3.count.aa$y))/n3.count.diaa$ny

cps.aat.tca<-((n3.trinuc$aat*n3.trinuc$tca)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns
cps.aat.tcc<-((n3.trinuc$aat*n3.trinuc$tcc)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns
cps.aat.tcg<-((n3.trinuc$aat*n3.trinuc$tcg)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns
cps.aat.tct<-((n3.trinuc$aat*n3.trinuc$tct)/(n3.count.aa$n*n3.count.aa$s))/n3.count.diaa$ns

#Stop codon
#cps.aat.tga<-((n3.trinuc$aat*n3.trinuc$tga)/(n3.count.aa$k*n3.count.aa$k))/n3.count.diaa$kk
cps.aat.tgc<-((n3.trinuc$aat*n3.trinuc$tgc)/(n3.count.aa$n*n3.count.aa$c))/n3.count.diaa$nc
cps.aat.tgg<-((n3.trinuc$aat*n3.trinuc$tgg)/(n3.count.aa$n*n3.count.aa$w))/n3.count.diaa$nw
cps.aat.tgt<-((n3.trinuc$aat*n3.trinuc$tgt)/(n3.count.aa$n*n3.count.aa$c))/n3.count.diaa$nc

cps.aat.tta<-((n3.trinuc$aat*n3.trinuc$tta)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aat.ttc<-((n3.trinuc$aat*n3.trinuc$ttc)/(n3.count.aa$n*n3.count.aa$f))/n3.count.diaa$nf
cps.aat.ttg<-((n3.trinuc$aat*n3.trinuc$ttg)/(n3.count.aa$n*n3.count.aa$l))/n3.count.diaa$nl
cps.aat.ttt<-((n3.trinuc$aat*n3.trinuc$ttt)/(n3.count.aa$n*n3.count.aa$f))/n3.count.diaa$nf






cps.aca.aaa<-((n3.trinuc$aca*n3.trinuc$aaa)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$tk
cps.aca.aac<-((n3.trinuc$aca*n3.trinuc$aac)/(n3.count.aa$t*n3.count.aa$n))/n3.count.diaa$tn
cps.aca.aag<-((n3.trinuc$aca*n3.trinuc$aag)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$tk
cps.aca.aat<-((n3.trinuc$aca*n3.trinuc$aat)/(n3.count.aa$t*n3.count.aa$n))/n3.count.diaa$tn

cps.aca.aca<-((n3.trinuc$aca*n3.trinuc$aca)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.aca.acc<-((n3.trinuc$aca*n3.trinuc$acc)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.aca.acg<-((n3.trinuc$aca*n3.trinuc$acg)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.aca.act<-((n3.trinuc$aca*n3.trinuc$act)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt

cps.aca.aga<-((n3.trinuc$aca*n3.trinuc$aga)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.aca.agc<-((n3.trinuc$aca*n3.trinuc$agc)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.aca.agg<-((n3.trinuc$aca*n3.trinuc$agg)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.aca.agt<-((n3.trinuc$aca*n3.trinuc$agt)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts

cps.aca.ata<-((n3.trinuc$aca*n3.trinuc$ata)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti
cps.aca.atc<-((n3.trinuc$aca*n3.trinuc$atc)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti
cps.aca.atg<-((n3.trinuc$aca*n3.trinuc$atg)/(n3.count.aa$t*n3.count.aa$m))/n3.count.diaa$tm
cps.aca.att<-((n3.trinuc$aca*n3.trinuc$att)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti

cps.aca.caa<-((n3.trinuc$aca*n3.trinuc$caa)/(n3.count.aa$t*n3.count.aa$q))/n3.count.diaa$tq
cps.aca.cac<-((n3.trinuc$aca*n3.trinuc$cac)/(n3.count.aa$t*n3.count.aa$h))/n3.count.diaa$th
cps.aca.cag<-((n3.trinuc$aca*n3.trinuc$cag)/(n3.count.aa$t*n3.count.aa$q))/n3.count.diaa$tq
cps.aca.cat<-((n3.trinuc$aca*n3.trinuc$cat)/(n3.count.aa$t*n3.count.aa$h))/n3.count.diaa$th

cps.aca.cca<-((n3.trinuc$aca*n3.trinuc$cca)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.aca.ccc<-((n3.trinuc$aca*n3.trinuc$ccc)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.aca.ccg<-((n3.trinuc$aca*n3.trinuc$ccg)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.aca.cct<-((n3.trinuc$aca*n3.trinuc$cct)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp

cps.aca.cga<-((n3.trinuc$aca*n3.trinuc$cga)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.aca.cgc<-((n3.trinuc$aca*n3.trinuc$cgc)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.aca.cgg<-((n3.trinuc$aca*n3.trinuc$cgg)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.aca.cgt<-((n3.trinuc$aca*n3.trinuc$cgt)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr

cps.aca.cta<-((n3.trinuc$aca*n3.trinuc$cta)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.aca.ctc<-((n3.trinuc$aca*n3.trinuc$ctc)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.aca.ctg<-((n3.trinuc$aca*n3.trinuc$ctg)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.aca.ctt<-((n3.trinuc$aca*n3.trinuc$ctt)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl

cps.aca.gaa<-((n3.trinuc$aca*n3.trinuc$gaa)/(n3.count.aa$t*n3.count.aa$e))/n3.count.diaa$te
cps.aca.gac<-((n3.trinuc$aca*n3.trinuc$gac)/(n3.count.aa$t*n3.count.aa$d))/n3.count.diaa$td
cps.aca.gag<-((n3.trinuc$aca*n3.trinuc$gag)/(n3.count.aa$t*n3.count.aa$e))/n3.count.diaa$te
cps.aca.gat<-((n3.trinuc$aca*n3.trinuc$gat)/(n3.count.aa$t*n3.count.aa$d))/n3.count.diaa$td

cps.aca.gca<-((n3.trinuc$aca*n3.trinuc$gca)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.aca.gcc<-((n3.trinuc$aca*n3.trinuc$gcc)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.aca.gcg<-((n3.trinuc$aca*n3.trinuc$gcg)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.aca.gct<-((n3.trinuc$aca*n3.trinuc$gct)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta

cps.aca.gga<-((n3.trinuc$aca*n3.trinuc$gga)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.aca.ggc<-((n3.trinuc$aca*n3.trinuc$ggc)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.aca.ggg<-((n3.trinuc$aca*n3.trinuc$ggg)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.aca.ggt<-((n3.trinuc$aca*n3.trinuc$ggt)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg

cps.aca.gta<-((n3.trinuc$aca*n3.trinuc$gta)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.aca.gtc<-((n3.trinuc$aca*n3.trinuc$gtc)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.aca.gtg<-((n3.trinuc$aca*n3.trinuc$gtg)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.aca.gtt<-((n3.trinuc$aca*n3.trinuc$gtt)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv

#Stop codon
#cps.aca.taa<-((n3.trinuc$aca*n3.trinuc$taa)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.aca.tac<-((n3.trinuc$aca*n3.trinuc$tac)/(n3.count.aa$t*n3.count.aa$y))/n3.count.diaa$ty
#Stop codon
#cps.aca.tag<-((n3.trinuc$aca*n3.trinuc$tag)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.aca.tat<-((n3.trinuc$aca*n3.trinuc$tat)/(n3.count.aa$t*n3.count.aa$y))/n3.count.diaa$ty

cps.aca.tca<-((n3.trinuc$aca*n3.trinuc$tca)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.aca.tcc<-((n3.trinuc$aca*n3.trinuc$tcc)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.aca.tcg<-((n3.trinuc$aca*n3.trinuc$tcg)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.aca.tct<-((n3.trinuc$aca*n3.trinuc$tct)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts

#Stop codon
#cps.aca.tga<-((n3.trinuc$aca*n3.trinuc$tga)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.aca.tgc<-((n3.trinuc$aca*n3.trinuc$tgc)/(n3.count.aa$t*n3.count.aa$c))/n3.count.diaa$tc
cps.aca.tgg<-((n3.trinuc$aca*n3.trinuc$tgg)/(n3.count.aa$t*n3.count.aa$w))/n3.count.diaa$tw
cps.aca.tgt<-((n3.trinuc$aca*n3.trinuc$tgt)/(n3.count.aa$t*n3.count.aa$c))/n3.count.diaa$tc

cps.aca.tta<-((n3.trinuc$aca*n3.trinuc$tta)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.aca.ttc<-((n3.trinuc$aca*n3.trinuc$ttc)/(n3.count.aa$t*n3.count.aa$f))/n3.count.diaa$tf
cps.aca.ttg<-((n3.trinuc$aca*n3.trinuc$ttg)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.aca.ttt<-((n3.trinuc$aca*n3.trinuc$ttt)/(n3.count.aa$t*n3.count.aa$f))/n3.count.diaa$tf







cps.acc.aaa<-((n3.trinuc$acc*n3.trinuc$aaa)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$tk
cps.acc.aac<-((n3.trinuc$acc*n3.trinuc$aac)/(n3.count.aa$t*n3.count.aa$n))/n3.count.diaa$tn
cps.acc.aag<-((n3.trinuc$acc*n3.trinuc$aag)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$tk
cps.acc.aat<-((n3.trinuc$acc*n3.trinuc$aat)/(n3.count.aa$t*n3.count.aa$n))/n3.count.diaa$tn

cps.acc.aca<-((n3.trinuc$acc*n3.trinuc$aca)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.acc.acc<-((n3.trinuc$acc*n3.trinuc$acc)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.acc.acg<-((n3.trinuc$acc*n3.trinuc$acg)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.acc.act<-((n3.trinuc$acc*n3.trinuc$act)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt

cps.acc.aga<-((n3.trinuc$acc*n3.trinuc$aga)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acc.agc<-((n3.trinuc$acc*n3.trinuc$agc)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.acc.agg<-((n3.trinuc$acc*n3.trinuc$agg)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acc.agt<-((n3.trinuc$acc*n3.trinuc$agt)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts

cps.acc.ata<-((n3.trinuc$acc*n3.trinuc$ata)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti
cps.acc.atc<-((n3.trinuc$acc*n3.trinuc$atc)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti
cps.acc.atg<-((n3.trinuc$acc*n3.trinuc$atg)/(n3.count.aa$t*n3.count.aa$m))/n3.count.diaa$tm
cps.acc.att<-((n3.trinuc$acc*n3.trinuc$att)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti

cps.acc.caa<-((n3.trinuc$acc*n3.trinuc$caa)/(n3.count.aa$t*n3.count.aa$q))/n3.count.diaa$tq
cps.acc.cac<-((n3.trinuc$acc*n3.trinuc$cac)/(n3.count.aa$t*n3.count.aa$h))/n3.count.diaa$th
cps.acc.cag<-((n3.trinuc$acc*n3.trinuc$cag)/(n3.count.aa$t*n3.count.aa$q))/n3.count.diaa$tq
cps.acc.cat<-((n3.trinuc$acc*n3.trinuc$cat)/(n3.count.aa$t*n3.count.aa$h))/n3.count.diaa$th

cps.acc.cca<-((n3.trinuc$acc*n3.trinuc$cca)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.acc.ccc<-((n3.trinuc$acc*n3.trinuc$ccc)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.acc.ccg<-((n3.trinuc$acc*n3.trinuc$ccg)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.acc.cct<-((n3.trinuc$acc*n3.trinuc$cct)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp

cps.acc.cga<-((n3.trinuc$acc*n3.trinuc$cga)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acc.cgc<-((n3.trinuc$acc*n3.trinuc$cgc)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acc.cgg<-((n3.trinuc$acc*n3.trinuc$cgg)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acc.cgt<-((n3.trinuc$acc*n3.trinuc$cgt)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr

cps.acc.cta<-((n3.trinuc$acc*n3.trinuc$cta)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acc.ctc<-((n3.trinuc$acc*n3.trinuc$ctc)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acc.ctg<-((n3.trinuc$acc*n3.trinuc$ctg)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acc.ctt<-((n3.trinuc$acc*n3.trinuc$ctt)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl

cps.acc.gaa<-((n3.trinuc$acc*n3.trinuc$gaa)/(n3.count.aa$t*n3.count.aa$e))/n3.count.diaa$te
cps.acc.gac<-((n3.trinuc$acc*n3.trinuc$gac)/(n3.count.aa$t*n3.count.aa$d))/n3.count.diaa$td
cps.acc.gag<-((n3.trinuc$acc*n3.trinuc$gag)/(n3.count.aa$t*n3.count.aa$e))/n3.count.diaa$te
cps.acc.gat<-((n3.trinuc$acc*n3.trinuc$gat)/(n3.count.aa$t*n3.count.aa$d))/n3.count.diaa$td

cps.acc.gca<-((n3.trinuc$acc*n3.trinuc$gca)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.acc.gcc<-((n3.trinuc$acc*n3.trinuc$gcc)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.acc.gcg<-((n3.trinuc$acc*n3.trinuc$gcg)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.acc.gct<-((n3.trinuc$acc*n3.trinuc$gct)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta

cps.acc.gga<-((n3.trinuc$acc*n3.trinuc$gga)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.acc.ggc<-((n3.trinuc$acc*n3.trinuc$ggc)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.acc.ggg<-((n3.trinuc$acc*n3.trinuc$ggg)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.acc.ggt<-((n3.trinuc$acc*n3.trinuc$ggt)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg

cps.acc.gta<-((n3.trinuc$acc*n3.trinuc$gta)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.acc.gtc<-((n3.trinuc$acc*n3.trinuc$gtc)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.acc.gtg<-((n3.trinuc$acc*n3.trinuc$gtg)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.acc.gtt<-((n3.trinuc$acc*n3.trinuc$gtt)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv

#Stop codon
#cps.acc.taa<-((n3.trinuc$acc*n3.trinuc$taa)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.acc.tac<-((n3.trinuc$acc*n3.trinuc$tac)/(n3.count.aa$t*n3.count.aa$y))/n3.count.diaa$ty
#stop codon
#cps.acc.tag<-((n3.trinuc$acc*n3.trinuc$tag)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.acc.tat<-((n3.trinuc$acc*n3.trinuc$tat)/(n3.count.aa$t*n3.count.aa$y))/n3.count.diaa$ty

cps.acc.tca<-((n3.trinuc$acc*n3.trinuc$tca)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.acc.tcc<-((n3.trinuc$acc*n3.trinuc$tcc)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.acc.tcg<-((n3.trinuc$acc*n3.trinuc$tcg)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.acc.tct<-((n3.trinuc$acc*n3.trinuc$tct)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts

#Stop codon
#cps.acc.tga<-((n3.trinuc$acc*n3.trinuc$tga)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.acc.tgc<-((n3.trinuc$acc*n3.trinuc$tgc)/(n3.count.aa$t*n3.count.aa$c))/n3.count.diaa$tc
cps.acc.tgg<-((n3.trinuc$acc*n3.trinuc$tgg)/(n3.count.aa$t*n3.count.aa$w))/n3.count.diaa$tw
cps.acc.tgt<-((n3.trinuc$acc*n3.trinuc$tgt)/(n3.count.aa$t*n3.count.aa$c))/n3.count.diaa$tc

cps.acc.tta<-((n3.trinuc$acc*n3.trinuc$tta)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acc.ttc<-((n3.trinuc$acc*n3.trinuc$ttc)/(n3.count.aa$t*n3.count.aa$f))/n3.count.diaa$tf
cps.acc.ttg<-((n3.trinuc$acc*n3.trinuc$ttg)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acc.ttt<-((n3.trinuc$acc*n3.trinuc$ttt)/(n3.count.aa$t*n3.count.aa$f))/n3.count.diaa$tf







cps.acg.aaa<-((n3.trinuc$acg*n3.trinuc$aaa)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$tk
cps.acg.aac<-((n3.trinuc$acg*n3.trinuc$aac)/(n3.count.aa$t*n3.count.aa$n))/n3.count.diaa$tn
cps.acg.aag<-((n3.trinuc$acg*n3.trinuc$aag)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$tk
cps.acg.aat<-((n3.trinuc$acg*n3.trinuc$aat)/(n3.count.aa$t*n3.count.aa$n))/n3.count.diaa$tn

cps.acg.aca<-((n3.trinuc$acg*n3.trinuc$aca)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.acg.acc<-((n3.trinuc$acg*n3.trinuc$acc)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.acg.acg<-((n3.trinuc$acg*n3.trinuc$acg)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.acg.act<-((n3.trinuc$acg*n3.trinuc$act)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt

cps.acg.aga<-((n3.trinuc$acg*n3.trinuc$aga)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acg.agc<-((n3.trinuc$acg*n3.trinuc$agc)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.acg.agg<-((n3.trinuc$acg*n3.trinuc$agg)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acg.agt<-((n3.trinuc$acg*n3.trinuc$agt)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts

cps.acg.ata<-((n3.trinuc$acg*n3.trinuc$ata)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti
cps.acg.atc<-((n3.trinuc$acg*n3.trinuc$atc)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti
cps.acg.atg<-((n3.trinuc$acg*n3.trinuc$atg)/(n3.count.aa$t*n3.count.aa$m))/n3.count.diaa$tm
cps.acg.att<-((n3.trinuc$acg*n3.trinuc$att)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti

cps.acg.caa<-((n3.trinuc$acg*n3.trinuc$caa)/(n3.count.aa$t*n3.count.aa$q))/n3.count.diaa$tq
cps.acg.cac<-((n3.trinuc$acg*n3.trinuc$cac)/(n3.count.aa$t*n3.count.aa$h))/n3.count.diaa$th
cps.acg.cag<-((n3.trinuc$acg*n3.trinuc$cag)/(n3.count.aa$t*n3.count.aa$q))/n3.count.diaa$tq
cps.acg.cat<-((n3.trinuc$acg*n3.trinuc$cat)/(n3.count.aa$t*n3.count.aa$h))/n3.count.diaa$th

cps.acg.cca<-((n3.trinuc$acg*n3.trinuc$cca)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.acg.ccc<-((n3.trinuc$acg*n3.trinuc$ccc)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.acg.ccg<-((n3.trinuc$acg*n3.trinuc$ccg)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.acg.cct<-((n3.trinuc$acg*n3.trinuc$cct)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp

cps.acg.cga<-((n3.trinuc$acg*n3.trinuc$cga)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acg.cgc<-((n3.trinuc$acg*n3.trinuc$cgc)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acg.cgg<-((n3.trinuc$acg*n3.trinuc$cgg)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.acg.cgt<-((n3.trinuc$acg*n3.trinuc$cgt)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr

cps.acg.cta<-((n3.trinuc$acg*n3.trinuc$cta)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acg.ctc<-((n3.trinuc$acg*n3.trinuc$ctc)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acg.ctg<-((n3.trinuc$acg*n3.trinuc$ctg)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acg.ctt<-((n3.trinuc$acg*n3.trinuc$ctt)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl

cps.acg.gaa<-((n3.trinuc$acg*n3.trinuc$gaa)/(n3.count.aa$t*n3.count.aa$e))/n3.count.diaa$te
cps.acg.gac<-((n3.trinuc$acg*n3.trinuc$gac)/(n3.count.aa$t*n3.count.aa$d))/n3.count.diaa$td
cps.acg.gag<-((n3.trinuc$acg*n3.trinuc$gag)/(n3.count.aa$t*n3.count.aa$e))/n3.count.diaa$te
cps.acg.gat<-((n3.trinuc$acg*n3.trinuc$gat)/(n3.count.aa$t*n3.count.aa$d))/n3.count.diaa$td

cps.acg.gca<-((n3.trinuc$acg*n3.trinuc$gca)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.acg.gcc<-((n3.trinuc$acg*n3.trinuc$gcc)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.acg.gcg<-((n3.trinuc$acg*n3.trinuc$gcg)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.acg.gct<-((n3.trinuc$acg*n3.trinuc$gct)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta

cps.acg.gga<-((n3.trinuc$acg*n3.trinuc$gga)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.acg.ggc<-((n3.trinuc$acg*n3.trinuc$ggc)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.acg.ggg<-((n3.trinuc$acg*n3.trinuc$ggg)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.acg.ggt<-((n3.trinuc$acg*n3.trinuc$ggt)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg

cps.acg.gta<-((n3.trinuc$acg*n3.trinuc$gta)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.acg.gtc<-((n3.trinuc$acg*n3.trinuc$gtc)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.acg.gtg<-((n3.trinuc$acg*n3.trinuc$gtg)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.acg.gtt<-((n3.trinuc$acg*n3.trinuc$gtt)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv

#Stop codon
#cps.acg.taa<-((n3.trinuc$acg*n3.trinuc$taa)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.acg.tac<-((n3.trinuc$acg*n3.trinuc$tac)/(n3.count.aa$t*n3.count.aa$y))/n3.count.diaa$ty
#Stop codon
#cps.acg.tag<-((n3.trinuc$acg*n3.trinuc$tag)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.acg.tat<-((n3.trinuc$acg*n3.trinuc$tat)/(n3.count.aa$t*n3.count.aa$y))/n3.count.diaa$ty

cps.acg.tca<-((n3.trinuc$acg*n3.trinuc$tca)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.acg.tcc<-((n3.trinuc$acg*n3.trinuc$tcc)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.acg.tcg<-((n3.trinuc$acg*n3.trinuc$tcg)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.acg.tct<-((n3.trinuc$acg*n3.trinuc$tct)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts

#Stop codon
#cps.acg.tga<-((n3.trinuc$acg*n3.trinuc$tga)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.acg.tgc<-((n3.trinuc$acg*n3.trinuc$tgc)/(n3.count.aa$t*n3.count.aa$c))/n3.count.diaa$tc
cps.acg.tgg<-((n3.trinuc$acg*n3.trinuc$tgg)/(n3.count.aa$t*n3.count.aa$w))/n3.count.diaa$tw
cps.acg.tgt<-((n3.trinuc$acg*n3.trinuc$tgt)/(n3.count.aa$t*n3.count.aa$c))/n3.count.diaa$tc

cps.acg.tta<-((n3.trinuc$acg*n3.trinuc$tta)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acg.ttc<-((n3.trinuc$acg*n3.trinuc$ttc)/(n3.count.aa$t*n3.count.aa$f))/n3.count.diaa$tf
cps.acg.ttg<-((n3.trinuc$acg*n3.trinuc$ttg)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.acg.ttt<-((n3.trinuc$acg*n3.trinuc$ttt)/(n3.count.aa$t*n3.count.aa$f))/n3.count.diaa$tf







cps.act.aaa<-((n3.trinuc$act*n3.trinuc$aaa)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$tk
cps.act.aac<-((n3.trinuc$act*n3.trinuc$aac)/(n3.count.aa$t*n3.count.aa$n))/n3.count.diaa$tn
cps.act.aag<-((n3.trinuc$act*n3.trinuc$aag)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$tk
cps.act.aat<-((n3.trinuc$act*n3.trinuc$aat)/(n3.count.aa$t*n3.count.aa$n))/n3.count.diaa$tn

cps.act.aca<-((n3.trinuc$act*n3.trinuc$aca)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.act.acc<-((n3.trinuc$act*n3.trinuc$acc)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.act.acg<-((n3.trinuc$act*n3.trinuc$acg)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt
cps.act.act<-((n3.trinuc$act*n3.trinuc$act)/(n3.count.aa$t*n3.count.aa$t))/n3.count.diaa$tt

cps.act.aga<-((n3.trinuc$act*n3.trinuc$aga)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.act.agc<-((n3.trinuc$act*n3.trinuc$agc)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.act.agg<-((n3.trinuc$act*n3.trinuc$agg)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.act.agt<-((n3.trinuc$act*n3.trinuc$agt)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts

cps.act.ata<-((n3.trinuc$act*n3.trinuc$ata)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti
cps.act.atc<-((n3.trinuc$act*n3.trinuc$atc)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti
cps.act.atg<-((n3.trinuc$act*n3.trinuc$atg)/(n3.count.aa$t*n3.count.aa$m))/n3.count.diaa$tw
cps.act.att<-((n3.trinuc$act*n3.trinuc$att)/(n3.count.aa$t*n3.count.aa$i))/n3.count.diaa$ti

cps.act.caa<-((n3.trinuc$act*n3.trinuc$caa)/(n3.count.aa$t*n3.count.aa$q))/n3.count.diaa$tq
cps.act.cac<-((n3.trinuc$act*n3.trinuc$cac)/(n3.count.aa$t*n3.count.aa$h))/n3.count.diaa$th
cps.act.cag<-((n3.trinuc$act*n3.trinuc$cag)/(n3.count.aa$t*n3.count.aa$q))/n3.count.diaa$tq
cps.act.cat<-((n3.trinuc$act*n3.trinuc$cat)/(n3.count.aa$t*n3.count.aa$h))/n3.count.diaa$th

cps.act.cca<-((n3.trinuc$act*n3.trinuc$cca)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.act.ccc<-((n3.trinuc$act*n3.trinuc$ccc)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.act.ccg<-((n3.trinuc$act*n3.trinuc$ccg)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp
cps.act.cct<-((n3.trinuc$act*n3.trinuc$cct)/(n3.count.aa$t*n3.count.aa$p))/n3.count.diaa$tp

cps.act.cga<-((n3.trinuc$act*n3.trinuc$cga)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.act.cgc<-((n3.trinuc$act*n3.trinuc$cgc)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.act.cgg<-((n3.trinuc$act*n3.trinuc$cgg)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr
cps.act.cgt<-((n3.trinuc$act*n3.trinuc$cgt)/(n3.count.aa$t*n3.count.aa$r))/n3.count.diaa$tr

cps.act.cta<-((n3.trinuc$act*n3.trinuc$cta)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.act.ctc<-((n3.trinuc$act*n3.trinuc$ctc)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.act.ctg<-((n3.trinuc$act*n3.trinuc$ctg)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl
cps.act.ctt<-((n3.trinuc$act*n3.trinuc$ctt)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$tl

cps.act.gaa<-((n3.trinuc$act*n3.trinuc$gaa)/(n3.count.aa$t*n3.count.aa$e))/n3.count.diaa$te
cps.act.gac<-((n3.trinuc$act*n3.trinuc$gac)/(n3.count.aa$t*n3.count.aa$d))/n3.count.diaa$td
cps.act.gag<-((n3.trinuc$act*n3.trinuc$gag)/(n3.count.aa$t*n3.count.aa$e))/n3.count.diaa$te
cps.act.gat<-((n3.trinuc$act*n3.trinuc$gat)/(n3.count.aa$t*n3.count.aa$d))/n3.count.diaa$td

cps.act.gca<-((n3.trinuc$act*n3.trinuc$gca)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.act.gcc<-((n3.trinuc$act*n3.trinuc$gcc)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.act.gcg<-((n3.trinuc$act*n3.trinuc$gcg)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta
cps.act.gct<-((n3.trinuc$act*n3.trinuc$gct)/(n3.count.aa$t*n3.count.aa$a))/n3.count.diaa$ta

cps.act.gga<-((n3.trinuc$act*n3.trinuc$gga)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.act.ggc<-((n3.trinuc$act*n3.trinuc$ggc)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.act.ggg<-((n3.trinuc$act*n3.trinuc$ggg)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg
cps.act.ggt<-((n3.trinuc$act*n3.trinuc$ggt)/(n3.count.aa$t*n3.count.aa$g))/n3.count.diaa$tg

cps.act.gta<-((n3.trinuc$act*n3.trinuc$gta)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.act.gtc<-((n3.trinuc$act*n3.trinuc$gtc)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.act.gtg<-((n3.trinuc$act*n3.trinuc$gtg)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv
cps.act.gtt<-((n3.trinuc$act*n3.trinuc$gtt)/(n3.count.aa$t*n3.count.aa$v))/n3.count.diaa$tv

#Stop codon
#cps.act.taa<-((n3.trinuc$act*n3.trinuc$taa)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.act.tac<-((n3.trinuc$act*n3.trinuc$tac)/(n3.count.aa$t*n3.count.aa$y))/n3.count.diaa$ty
#Stop codon
#cps.act.tag<-((n3.trinuc$act*n3.trinuc$tag)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.act.tat<-((n3.trinuc$act*n3.trinuc$tat)/(n3.count.aa$t*n3.count.aa$y))/n3.count.diaa$ty

cps.act.tca<-((n3.trinuc$act*n3.trinuc$tca)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.act.tcc<-((n3.trinuc$act*n3.trinuc$tcc)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.act.tcg<-((n3.trinuc$act*n3.trinuc$tcg)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts
cps.act.tct<-((n3.trinuc$act*n3.trinuc$tct)/(n3.count.aa$t*n3.count.aa$s))/n3.count.diaa$ts

#Stop codon
#cps.act.tga<-((n3.trinuc$act*n3.trinuc$tga)/(n3.count.aa$t*n3.count.aa$k))/n3.count.diaa$kk
cps.act.tgc<-((n3.trinuc$act*n3.trinuc$tgc)/(n3.count.aa$t*n3.count.aa$c))/n3.count.diaa$ts
cps.act.tgg<-((n3.trinuc$act*n3.trinuc$tgg)/(n3.count.aa$t*n3.count.aa$w))/n3.count.diaa$ts
cps.act.tgt<-((n3.trinuc$act*n3.trinuc$tgt)/(n3.count.aa$t*n3.count.aa$c))/n3.count.diaa$ts

cps.act.tta<-((n3.trinuc$act*n3.trinuc$tta)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$ts
cps.act.ttc<-((n3.trinuc$act*n3.trinuc$ttc)/(n3.count.aa$t*n3.count.aa$f))/n3.count.diaa$ts
cps.act.ttg<-((n3.trinuc$act*n3.trinuc$ttg)/(n3.count.aa$t*n3.count.aa$l))/n3.count.diaa$ts
cps.act.ttt<-((n3.trinuc$act*n3.trinuc$ttt)/(n3.count.aa$t*n3.count.aa$f))/n3.count.diaa$ts













cps.aga.aaa<-((n3.trinuc$aga*n3.trinuc$aaa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.aga.aac<-((n3.trinuc$aga*n3.trinuc$aac)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn
cps.aga.aag<-((n3.trinuc$aga*n3.trinuc$aag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.aga.aat<-((n3.trinuc$aga*n3.trinuc$aat)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn

cps.aga.aca<-((n3.trinuc$aga*n3.trinuc$aca)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.aga.acc<-((n3.trinuc$aga*n3.trinuc$acc)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.aga.acg<-((n3.trinuc$aga*n3.trinuc$acg)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.aga.act<-((n3.trinuc$aga*n3.trinuc$act)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt

cps.aga.aga<-((n3.trinuc$aga*n3.trinuc$aga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.aga.agc<-((n3.trinuc$aga*n3.trinuc$agc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.aga.agg<-((n3.trinuc$aga*n3.trinuc$agg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.aga.agt<-((n3.trinuc$aga*n3.trinuc$agt)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

cps.aga.ata<-((n3.trinuc$aga*n3.trinuc$ata)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.aga.atc<-((n3.trinuc$aga*n3.trinuc$atc)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.aga.atg<-((n3.trinuc$aga*n3.trinuc$atg)/(n3.count.aa$r*n3.count.aa$m))/n3.count.diaa$rm
cps.aga.att<-((n3.trinuc$aga*n3.trinuc$att)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri

cps.aga.caa<-((n3.trinuc$aga*n3.trinuc$caa)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.aga.cac<-((n3.trinuc$aga*n3.trinuc$cac)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh
cps.aga.cag<-((n3.trinuc$aga*n3.trinuc$cag)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.aga.cat<-((n3.trinuc$aga*n3.trinuc$cat)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh

cps.aga.cca<-((n3.trinuc$aga*n3.trinuc$cca)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.aga.ccc<-((n3.trinuc$aga*n3.trinuc$ccc)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.aga.ccg<-((n3.trinuc$aga*n3.trinuc$ccg)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.aga.cct<-((n3.trinuc$aga*n3.trinuc$cct)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp

cps.aga.cga<-((n3.trinuc$aga*n3.trinuc$cga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.aga.cgc<-((n3.trinuc$aga*n3.trinuc$cgc)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.aga.cgg<-((n3.trinuc$aga*n3.trinuc$cgg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.aga.cgt<-((n3.trinuc$aga*n3.trinuc$cgt)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr

cps.aga.cta<-((n3.trinuc$aga*n3.trinuc$cta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.aga.ctc<-((n3.trinuc$aga*n3.trinuc$ctc)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.aga.ctg<-((n3.trinuc$aga*n3.trinuc$ctg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.aga.ctt<-((n3.trinuc$aga*n3.trinuc$ctt)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl

cps.aga.gaa<-((n3.trinuc$aga*n3.trinuc$gaa)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.aga.gac<-((n3.trinuc$aga*n3.trinuc$gac)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd
cps.aga.gag<-((n3.trinuc$aga*n3.trinuc$gag)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.aga.gat<-((n3.trinuc$aga*n3.trinuc$gat)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd

cps.aga.gca<-((n3.trinuc$aga*n3.trinuc$gca)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.aga.gcc<-((n3.trinuc$aga*n3.trinuc$gcc)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.aga.gcg<-((n3.trinuc$aga*n3.trinuc$gcg)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.aga.gct<-((n3.trinuc$aga*n3.trinuc$gct)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra

cps.aga.gga<-((n3.trinuc$aga*n3.trinuc$gga)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.aga.ggc<-((n3.trinuc$aga*n3.trinuc$ggc)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.aga.ggg<-((n3.trinuc$aga*n3.trinuc$ggg)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.aga.ggt<-((n3.trinuc$aga*n3.trinuc$ggt)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg

cps.aga.gta<-((n3.trinuc$aga*n3.trinuc$gta)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.aga.gtc<-((n3.trinuc$aga*n3.trinuc$gtc)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.aga.gtg<-((n3.trinuc$aga*n3.trinuc$gtg)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.aga.gtt<-((n3.trinuc$aga*n3.trinuc$gtt)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv

#Stop codon
#cps.aga.taa<-((n3.trinuc$aga*n3.trinuc$taa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.aga.tac<-((n3.trinuc$aga*n3.trinuc$tac)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry
#Stop codon
#cps.aga.tag<-((n3.trinuc$aga*n3.trinuc$tag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.aga.tat<-((n3.trinuc$aga*n3.trinuc$tat)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry

cps.aga.tca<-((n3.trinuc$aga*n3.trinuc$tca)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.aga.tcc<-((n3.trinuc$aga*n3.trinuc$tcc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.aga.tcg<-((n3.trinuc$aga*n3.trinuc$tcg)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.aga.tct<-((n3.trinuc$aga*n3.trinuc$tct)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

#Stop codon
#cps.aga.tga<-((n3.trinuc$aga*n3.trinuc$tga)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.aga.tgc<-((n3.trinuc$aga*n3.trinuc$tgc)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc
cps.aga.tgg<-((n3.trinuc$aga*n3.trinuc$tgg)/(n3.count.aa$r*n3.count.aa$w))/n3.count.diaa$rw
cps.aga.tgt<-((n3.trinuc$aga*n3.trinuc$tgt)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc

cps.aga.tta<-((n3.trinuc$aga*n3.trinuc$tta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.aga.ttc<-((n3.trinuc$aga*n3.trinuc$ttc)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf
cps.aga.ttg<-((n3.trinuc$aga*n3.trinuc$ttg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.aga.ttt<-((n3.trinuc$aga*n3.trinuc$ttt)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf







cps.agc.aaa<-((n3.trinuc$agc*n3.trinuc$aaa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.agc.aac<-((n3.trinuc$acg*n3.trinuc$aac)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn
cps.agc.aag<-((n3.trinuc$agc*n3.trinuc$aag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.agc.aat<-((n3.trinuc$agc*n3.trinuc$aat)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn

cps.agc.aca<-((n3.trinuc$agc*n3.trinuc$aca)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.agc.acc<-((n3.trinuc$agc*n3.trinuc$acc)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.agc.acg<-((n3.trinuc$agc*n3.trinuc$acg)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.agc.act<-((n3.trinuc$agc*n3.trinuc$act)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st

cps.agc.aga<-((n3.trinuc$agc*n3.trinuc$aga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agc.agc<-((n3.trinuc$agc*n3.trinuc$agc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.agc.agg<-((n3.trinuc$agc*n3.trinuc$agg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agc.agt<-((n3.trinuc$agc*n3.trinuc$agt)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

cps.agc.ata<-((n3.trinuc$agc*n3.trinuc$ata)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.agc.atc<-((n3.trinuc$agc*n3.trinuc$atc)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.agc.atg<-((n3.trinuc$agc*n3.trinuc$atg)/(n3.count.aa$s*n3.count.aa$m))/n3.count.diaa$sm
cps.agc.att<-((n3.trinuc$agc*n3.trinuc$att)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si

cps.agc.caa<-((n3.trinuc$agc*n3.trinuc$caa)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.agc.cac<-((n3.trinuc$agc*n3.trinuc$cac)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh
cps.agc.cag<-((n3.trinuc$agc*n3.trinuc$cag)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.agc.cat<-((n3.trinuc$agc*n3.trinuc$cat)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh

cps.agc.cca<-((n3.trinuc$agc*n3.trinuc$cca)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.agc.ccc<-((n3.trinuc$agc*n3.trinuc$ccc)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.agc.ccg<-((n3.trinuc$agc*n3.trinuc$ccg)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.agc.cct<-((n3.trinuc$agc*n3.trinuc$cct)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp

cps.agc.cga<-((n3.trinuc$agc*n3.trinuc$cga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agc.cgc<-((n3.trinuc$agc*n3.trinuc$cgc)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agc.cgg<-((n3.trinuc$agc*n3.trinuc$cgg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agc.cgt<-((n3.trinuc$agc*n3.trinuc$cgt)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr

cps.agc.cta<-((n3.trinuc$agc*n3.trinuc$cta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agc.ctc<-((n3.trinuc$agc*n3.trinuc$ctc)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agc.ctg<-((n3.trinuc$agc*n3.trinuc$ctg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agc.ctt<-((n3.trinuc$agc*n3.trinuc$ctt)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl

cps.agc.gaa<-((n3.trinuc$agc*n3.trinuc$gaa)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.agc.gac<-((n3.trinuc$agc*n3.trinuc$gac)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd
cps.agc.gag<-((n3.trinuc$agc*n3.trinuc$gag)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.agc.gat<-((n3.trinuc$agc*n3.trinuc$gat)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd

cps.agc.gca<-((n3.trinuc$agc*n3.trinuc$gca)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.agc.gcc<-((n3.trinuc$agc*n3.trinuc$gcc)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.agc.gcg<-((n3.trinuc$agc*n3.trinuc$gcg)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.agc.gct<-((n3.trinuc$agc*n3.trinuc$gct)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa

cps.agc.gga<-((n3.trinuc$agc*n3.trinuc$gga)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.agc.ggc<-((n3.trinuc$agc*n3.trinuc$ggc)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.agc.ggg<-((n3.trinuc$agc*n3.trinuc$ggg)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.agc.ggt<-((n3.trinuc$agc*n3.trinuc$ggt)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg

cps.agc.gta<-((n3.trinuc$agc*n3.trinuc$gta)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.agc.gtc<-((n3.trinuc$agc*n3.trinuc$gtc)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.agc.gtg<-((n3.trinuc$agc*n3.trinuc$gtg)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.agc.gtt<-((n3.trinuc$agc*n3.trinuc$gtt)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv

#Stop codon
#cps.agc.taa<-((n3.trinuc$agc*n3.trinuc$taa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$kk
cps.agc.tac<-((n3.trinuc$agc*n3.trinuc$tac)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy
#Stop codon
#cps.agc.tag<-((n3.trinuc$agc*n3.trinuc$tag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$kk
cps.agc.tat<-((n3.trinuc$agc*n3.trinuc$tat)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy

cps.agc.tca<-((n3.trinuc$agc*n3.trinuc$tca)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.agc.tcc<-((n3.trinuc$agc*n3.trinuc$tcc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.agc.tcg<-((n3.trinuc$agc*n3.trinuc$tcg)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.agc.tct<-((n3.trinuc$agc*n3.trinuc$tct)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

#Stop codon
#cps.agc.tga<-((n3.trinuc$agc*n3.trinuc$tga)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$kk
cps.agc.tgc<-((n3.trinuc$agc*n3.trinuc$tgc)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc
cps.agc.tgg<-((n3.trinuc$agc*n3.trinuc$tgg)/(n3.count.aa$s*n3.count.aa$w))/n3.count.diaa$sw
cps.agc.tgt<-((n3.trinuc$acg*n3.trinuc$tgt)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc

cps.agc.tta<-((n3.trinuc$agc*n3.trinuc$tta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agc.ttc<-((n3.trinuc$agc*n3.trinuc$ttc)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf
cps.agc.ttg<-((n3.trinuc$agc*n3.trinuc$ttg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agc.ttt<-((n3.trinuc$agc*n3.trinuc$ttt)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf







cps.agg.aaa<-((n3.trinuc$agg*n3.trinuc$aaa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.agg.aac<-((n3.trinuc$agg*n3.trinuc$aac)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn
cps.agg.aag<-((n3.trinuc$agg*n3.trinuc$aag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.agg.aat<-((n3.trinuc$agg*n3.trinuc$aat)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn

cps.agg.aca<-((n3.trinuc$agg*n3.trinuc$aca)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.agg.acc<-((n3.trinuc$agg*n3.trinuc$acc)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.agg.acg<-((n3.trinuc$agg*n3.trinuc$acg)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.agg.act<-((n3.trinuc$agg*n3.trinuc$act)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt

cps.agg.aga<-((n3.trinuc$agg*n3.trinuc$aga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.agg.agc<-((n3.trinuc$agg*n3.trinuc$agc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.agg.agg<-((n3.trinuc$agg*n3.trinuc$agg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.agg.agt<-((n3.trinuc$agg*n3.trinuc$agt)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

cps.agg.ata<-((n3.trinuc$agg*n3.trinuc$ata)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.agg.atc<-((n3.trinuc$agg*n3.trinuc$atc)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.agg.atg<-((n3.trinuc$agg*n3.trinuc$atg)/(n3.count.aa$r*n3.count.aa$m))/n3.count.diaa$rm
cps.agg.att<-((n3.trinuc$agg*n3.trinuc$att)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri

cps.agg.caa<-((n3.trinuc$agg*n3.trinuc$caa)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.agg.cac<-((n3.trinuc$agg*n3.trinuc$cac)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh
cps.agg.cag<-((n3.trinuc$agg*n3.trinuc$cag)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.agg.cat<-((n3.trinuc$agg*n3.trinuc$cat)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh

cps.agg.cca<-((n3.trinuc$agg*n3.trinuc$cca)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.agg.ccc<-((n3.trinuc$agg*n3.trinuc$ccc)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.agg.ccg<-((n3.trinuc$agg*n3.trinuc$ccg)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.agg.cct<-((n3.trinuc$agg*n3.trinuc$cct)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp

cps.agg.cga<-((n3.trinuc$agg*n3.trinuc$cga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.agg.cgc<-((n3.trinuc$agg*n3.trinuc$cgc)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.agg.cgg<-((n3.trinuc$agg*n3.trinuc$cgg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.agg.cgt<-((n3.trinuc$agg*n3.trinuc$cgt)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr

cps.agg.cta<-((n3.trinuc$agg*n3.trinuc$cta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.agg.ctc<-((n3.trinuc$agg*n3.trinuc$ctc)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.agg.ctg<-((n3.trinuc$agg*n3.trinuc$ctg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.agg.ctt<-((n3.trinuc$agg*n3.trinuc$ctt)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl

cps.agg.gaa<-((n3.trinuc$agg*n3.trinuc$gaa)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.agg.gac<-((n3.trinuc$agg*n3.trinuc$gac)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd
cps.agg.gag<-((n3.trinuc$agg*n3.trinuc$gag)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.agg.gat<-((n3.trinuc$agg*n3.trinuc$gat)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd

cps.agg.gca<-((n3.trinuc$agg*n3.trinuc$gca)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.agg.gcc<-((n3.trinuc$agg*n3.trinuc$gcc)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.agg.gcg<-((n3.trinuc$agg*n3.trinuc$gcg)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.agg.gct<-((n3.trinuc$agg*n3.trinuc$gct)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra

cps.agg.gga<-((n3.trinuc$agg*n3.trinuc$gga)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.agg.ggc<-((n3.trinuc$agg*n3.trinuc$ggc)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.agg.ggg<-((n3.trinuc$agg*n3.trinuc$ggg)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.agg.ggt<-((n3.trinuc$agg*n3.trinuc$ggt)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg

cps.agg.gta<-((n3.trinuc$agg*n3.trinuc$gta)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.agg.gtc<-((n3.trinuc$agg*n3.trinuc$gtc)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.agg.gtg<-((n3.trinuc$agg*n3.trinuc$gtg)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.agg.gtt<-((n3.trinuc$agg*n3.trinuc$gtt)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv

#Stop codon
#cps.agg.taa<-((n3.trinuc$agg*n3.trinuc$taa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.agg.tac<-((n3.trinuc$agg*n3.trinuc$tac)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry
#Stop codon
#cps.agg.tag<-((n3.trinuc$agg*n3.trinuc$tag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.agg.tat<-((n3.trinuc$agg*n3.trinuc$tat)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry

cps.agg.tca<-((n3.trinuc$agg*n3.trinuc$tca)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.agg.tcc<-((n3.trinuc$agg*n3.trinuc$tcc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.agg.tcg<-((n3.trinuc$agg*n3.trinuc$tcg)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.agg.tct<-((n3.trinuc$agg*n3.trinuc$tct)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

#Stop codon
#cps.agg.tga<-((n3.trinuc$agg*n3.trinuc$tga)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.agg.tgc<-((n3.trinuc$agg*n3.trinuc$tgc)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc
cps.agg.tgg<-((n3.trinuc$agg*n3.trinuc$tgg)/(n3.count.aa$r*n3.count.aa$w))/n3.count.diaa$rw
cps.agg.tgt<-((n3.trinuc$agg*n3.trinuc$tgt)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc

cps.agg.tta<-((n3.trinuc$agg*n3.trinuc$tta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.agg.ttc<-((n3.trinuc$agg*n3.trinuc$ttc)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf
cps.agg.ttg<-((n3.trinuc$agg*n3.trinuc$ttg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.agg.ttt<-((n3.trinuc$agg*n3.trinuc$ttt)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf







cps.agt.aaa<-((n3.trinuc$agt*n3.trinuc$aaa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.agt.aac<-((n3.trinuc$agt*n3.trinuc$aac)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn
cps.agt.aag<-((n3.trinuc$agt*n3.trinuc$aag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.agt.aat<-((n3.trinuc$agt*n3.trinuc$aat)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn

cps.agt.aca<-((n3.trinuc$agt*n3.trinuc$aca)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.agt.acc<-((n3.trinuc$agt*n3.trinuc$acc)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.agt.acg<-((n3.trinuc$agt*n3.trinuc$acg)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.agt.act<-((n3.trinuc$agt*n3.trinuc$act)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st

cps.agt.aga<-((n3.trinuc$agt*n3.trinuc$aga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agt.agc<-((n3.trinuc$agt*n3.trinuc$agc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.agt.agg<-((n3.trinuc$agt*n3.trinuc$agg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agt.agt<-((n3.trinuc$agt*n3.trinuc$agt)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

cps.agt.ata<-((n3.trinuc$agt*n3.trinuc$ata)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.agt.atc<-((n3.trinuc$agt*n3.trinuc$atc)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.agt.atg<-((n3.trinuc$agt*n3.trinuc$atg)/(n3.count.aa$s*n3.count.aa$m))/n3.count.diaa$sm
cps.agt.att<-((n3.trinuc$agt*n3.trinuc$att)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si

cps.agt.caa<-((n3.trinuc$agt*n3.trinuc$caa)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.agt.cac<-((n3.trinuc$agt*n3.trinuc$cac)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh
cps.agt.cag<-((n3.trinuc$agt*n3.trinuc$cag)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.agt.cat<-((n3.trinuc$agt*n3.trinuc$cat)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh

cps.agt.cca<-((n3.trinuc$agt*n3.trinuc$cca)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.agt.ccc<-((n3.trinuc$agt*n3.trinuc$ccc)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.agt.ccg<-((n3.trinuc$agt*n3.trinuc$ccg)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.agt.cct<-((n3.trinuc$agt*n3.trinuc$cct)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp

cps.agt.cga<-((n3.trinuc$agt*n3.trinuc$cga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agt.cgc<-((n3.trinuc$agt*n3.trinuc$cgc)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agt.cgg<-((n3.trinuc$agt*n3.trinuc$cgg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.agt.cgt<-((n3.trinuc$agt*n3.trinuc$cgt)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr

cps.agt.cta<-((n3.trinuc$agt*n3.trinuc$cta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agt.ctc<-((n3.trinuc$agt*n3.trinuc$ctc)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agt.ctg<-((n3.trinuc$agt*n3.trinuc$ctg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agt.ctt<-((n3.trinuc$agt*n3.trinuc$ctt)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl

cps.agt.gaa<-((n3.trinuc$agt*n3.trinuc$gaa)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.agt.gac<-((n3.trinuc$agt*n3.trinuc$gac)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd
cps.agt.gag<-((n3.trinuc$agt*n3.trinuc$gag)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.agt.gat<-((n3.trinuc$agt*n3.trinuc$gat)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd

cps.agt.gca<-((n3.trinuc$agt*n3.trinuc$gca)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.agt.gcc<-((n3.trinuc$agt*n3.trinuc$gcc)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.agt.gcg<-((n3.trinuc$agt*n3.trinuc$gcg)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.agt.gct<-((n3.trinuc$agt*n3.trinuc$gct)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa

cps.agt.gga<-((n3.trinuc$agt*n3.trinuc$gga)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.agt.ggc<-((n3.trinuc$agt*n3.trinuc$ggc)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.agt.ggg<-((n3.trinuc$agt*n3.trinuc$ggg)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.agt.ggt<-((n3.trinuc$agt*n3.trinuc$ggt)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg

cps.agt.gta<-((n3.trinuc$agt*n3.trinuc$gta)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.agt.gtc<-((n3.trinuc$agt*n3.trinuc$gtc)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.agt.gtg<-((n3.trinuc$agt*n3.trinuc$gtg)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.agt.gtt<-((n3.trinuc$agt*n3.trinuc$gtt)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv

#Stop codon
#cps.agt.taa<-((n3.trinuc$agt*n3.trinuc$taa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$kk
cps.agt.tac<-((n3.trinuc$agt*n3.trinuc$tac)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy
#Stop codon
#cps.agt.tag<-((n3.trinuc$agt*n3.trinuc$tag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$kk
cps.agt.tat<-((n3.trinuc$agt*n3.trinuc$tat)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy

cps.agt.tca<-((n3.trinuc$agt*n3.trinuc$tca)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.agt.tcc<-((n3.trinuc$agt*n3.trinuc$tcc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.agt.tcg<-((n3.trinuc$agt*n3.trinuc$tcg)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.agt.tct<-((n3.trinuc$agt*n3.trinuc$tct)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

#Stop codon
#cps.agt.tga<-((n3.trinuc$agt*n3.trinuc$tga)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$kk
cps.agt.tgc<-((n3.trinuc$agt*n3.trinuc$tgc)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc
cps.agt.tgg<-((n3.trinuc$agt*n3.trinuc$tgg)/(n3.count.aa$s*n3.count.aa$w))/n3.count.diaa$sw
cps.agt.tgt<-((n3.trinuc$agt*n3.trinuc$tgt)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc

cps.agt.tta<-((n3.trinuc$agt*n3.trinuc$tta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agt.ttc<-((n3.trinuc$agt*n3.trinuc$ttc)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf
cps.agt.ttg<-((n3.trinuc$agt*n3.trinuc$ttg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.agt.ttt<-((n3.trinuc$agt*n3.trinuc$ttt)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf













cps.ata.aaa<-((n3.trinuc$ata*n3.trinuc$aaa)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.ata.aac<-((n3.trinuc$ata*n3.trinuc$aac)/(n3.count.aa$i*n3.count.aa$n))/n3.count.diaa$"in"
cps.ata.aag<-((n3.trinuc$ata*n3.trinuc$aag)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.ata.aat<-((n3.trinuc$ata*n3.trinuc$aat)/(n3.count.aa$i*n3.count.aa$n))/n3.count.diaa$"in"

cps.ata.aca<-((n3.trinuc$ata*n3.trinuc$aca)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.ata.acc<-((n3.trinuc$ata*n3.trinuc$acc)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.ata.acg<-((n3.trinuc$ata*n3.trinuc$acg)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.ata.act<-((n3.trinuc$ata*n3.trinuc$act)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it

cps.ata.aga<-((n3.trinuc$ata*n3.trinuc$aga)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.ata.agc<-((n3.trinuc$ata*n3.trinuc$agc)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.ata.agg<-((n3.trinuc$ata*n3.trinuc$agg)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.ata.agt<-((n3.trinuc$ata*n3.trinuc$agt)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is

cps.ata.ata<-((n3.trinuc$ata*n3.trinuc$ata)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii
cps.ata.atc<-((n3.trinuc$ata*n3.trinuc$atc)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii
cps.ata.atg<-((n3.trinuc$ata*n3.trinuc$atg)/(n3.count.aa$i*n3.count.aa$m))/n3.count.diaa$im
cps.ata.att<-((n3.trinuc$ata*n3.trinuc$att)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii

cps.ata.caa<-((n3.trinuc$ata*n3.trinuc$caa)/(n3.count.aa$i*n3.count.aa$q))/n3.count.diaa$iq
cps.ata.cac<-((n3.trinuc$ata*n3.trinuc$cac)/(n3.count.aa$i*n3.count.aa$h))/n3.count.diaa$ih
cps.ata.cag<-((n3.trinuc$ata*n3.trinuc$cag)/(n3.count.aa$i*n3.count.aa$q))/n3.count.diaa$iq
cps.ata.cat<-((n3.trinuc$ata*n3.trinuc$cat)/(n3.count.aa$i*n3.count.aa$h))/n3.count.diaa$ih

cps.ata.cca<-((n3.trinuc$ata*n3.trinuc$cca)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.ata.ccc<-((n3.trinuc$ata*n3.trinuc$ccc)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.ata.ccg<-((n3.trinuc$ata*n3.trinuc$ccg)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.ata.cct<-((n3.trinuc$ata*n3.trinuc$cct)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip

cps.ata.cga<-((n3.trinuc$ata*n3.trinuc$cga)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.ata.cgc<-((n3.trinuc$ata*n3.trinuc$cgc)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.ata.cgg<-((n3.trinuc$ata*n3.trinuc$cgg)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.ata.cgt<-((n3.trinuc$ata*n3.trinuc$cgt)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir

cps.ata.cta<-((n3.trinuc$ata*n3.trinuc$cta)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.ata.ctc<-((n3.trinuc$ata*n3.trinuc$ctc)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.ata.ctg<-((n3.trinuc$ata*n3.trinuc$ctg)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.ata.ctt<-((n3.trinuc$ata*n3.trinuc$ctt)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il

cps.ata.gaa<-((n3.trinuc$ata*n3.trinuc$gaa)/(n3.count.aa$i*n3.count.aa$e))/n3.count.diaa$ie
cps.ata.gac<-((n3.trinuc$ata*n3.trinuc$gac)/(n3.count.aa$i*n3.count.aa$d))/n3.count.diaa$id
cps.ata.gag<-((n3.trinuc$ata*n3.trinuc$gag)/(n3.count.aa$i*n3.count.aa$e))/n3.count.diaa$ie
cps.ata.gat<-((n3.trinuc$ata*n3.trinuc$gat)/(n3.count.aa$i*n3.count.aa$d))/n3.count.diaa$id

cps.ata.gca<-((n3.trinuc$ata*n3.trinuc$gca)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.ata.gcc<-((n3.trinuc$ata*n3.trinuc$gcc)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.ata.gcg<-((n3.trinuc$ata*n3.trinuc$gcg)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.ata.gct<-((n3.trinuc$ata*n3.trinuc$gct)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia

cps.ata.gga<-((n3.trinuc$ata*n3.trinuc$gga)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.ata.ggc<-((n3.trinuc$ata*n3.trinuc$ggc)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.ata.ggg<-((n3.trinuc$ata*n3.trinuc$ggg)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.ata.ggt<-((n3.trinuc$ata*n3.trinuc$ggt)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig

cps.ata.gta<-((n3.trinuc$ata*n3.trinuc$gta)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.ata.gtc<-((n3.trinuc$ata*n3.trinuc$gtc)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.ata.gtg<-((n3.trinuc$ata*n3.trinuc$gtg)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.ata.gtt<-((n3.trinuc$ata*n3.trinuc$gtt)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv

#Stop codon
#cps.ata.taa<-((n3.trinuc$ata*n3.trinuc$taa)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.ata.tac<-((n3.trinuc$ata*n3.trinuc$tac)/(n3.count.aa$i*n3.count.aa$y))/n3.count.diaa$iy
#Stop codon
#cps.ata.tag<-((n3.trinuc$ata*n3.trinuc$tag)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.ata.tat<-((n3.trinuc$ata*n3.trinuc$tat)/(n3.count.aa$i*n3.count.aa$y))/n3.count.diaa$iy

cps.ata.tca<-((n3.trinuc$ata*n3.trinuc$tca)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.ata.tcc<-((n3.trinuc$ata*n3.trinuc$tcc)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.ata.tcg<-((n3.trinuc$ata*n3.trinuc$tcg)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.ata.tct<-((n3.trinuc$ata*n3.trinuc$tct)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is

#Stop codon
#cps.ata.tga<-((n3.trinuc$ata*n3.trinuc$tga)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.ata.tgc<-((n3.trinuc$ata*n3.trinuc$tgc)/(n3.count.aa$i*n3.count.aa$c))/n3.count.diaa$ic
cps.ata.tgg<-((n3.trinuc$ata*n3.trinuc$tgg)/(n3.count.aa$i*n3.count.aa$w))/n3.count.diaa$iw
cps.ata.tgt<-((n3.trinuc$ata*n3.trinuc$tgt)/(n3.count.aa$i*n3.count.aa$c))/n3.count.diaa$ic

cps.ata.tta<-((n3.trinuc$ata*n3.trinuc$tta)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.ata.ttc<-((n3.trinuc$ata*n3.trinuc$ttc)/(n3.count.aa$i*n3.count.aa$f))/n3.count.diaa$"if"
cps.ata.ttg<-((n3.trinuc$ata*n3.trinuc$ttg)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.ata.ttt<-((n3.trinuc$ata*n3.trinuc$ttt)/(n3.count.aa$i*n3.count.aa$f))/n3.count.diaa$"if"










cps.atc.aaa<-((n3.trinuc$atc*n3.trinuc$aaa)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.atc.aac<-((n3.trinuc$atc*n3.trinuc$aac)/(n3.count.aa$i*n3.count.aa$n))/n3.count.diaa$"in"
cps.atc.aag<-((n3.trinuc$atc*n3.trinuc$aag)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.atc.aat<-((n3.trinuc$atc*n3.trinuc$aat)/(n3.count.aa$i*n3.count.aa$n))/n3.count.diaa$"in"

cps.atc.aca<-((n3.trinuc$atc*n3.trinuc$aca)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.atc.acc<-((n3.trinuc$atc*n3.trinuc$acc)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.atc.acg<-((n3.trinuc$atc*n3.trinuc$acg)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.atc.act<-((n3.trinuc$atc*n3.trinuc$act)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it

cps.atc.aga<-((n3.trinuc$atc*n3.trinuc$aga)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.atc.agc<-((n3.trinuc$atc*n3.trinuc$agc)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.atc.agg<-((n3.trinuc$atc*n3.trinuc$agg)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.atc.agt<-((n3.trinuc$atc*n3.trinuc$agt)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is

cps.atc.ata<-((n3.trinuc$atc*n3.trinuc$ata)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii
cps.atc.atc<-((n3.trinuc$atc*n3.trinuc$atc)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii
cps.atc.atg<-((n3.trinuc$atc*n3.trinuc$atg)/(n3.count.aa$i*n3.count.aa$m))/n3.count.diaa$im
cps.atc.att<-((n3.trinuc$atc*n3.trinuc$att)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii

cps.atc.caa<-((n3.trinuc$atc*n3.trinuc$caa)/(n3.count.aa$i*n3.count.aa$q))/n3.count.diaa$iq
cps.atc.cac<-((n3.trinuc$atc*n3.trinuc$cac)/(n3.count.aa$i*n3.count.aa$h))/n3.count.diaa$ih
cps.atc.cag<-((n3.trinuc$atc*n3.trinuc$cag)/(n3.count.aa$i*n3.count.aa$q))/n3.count.diaa$iq
cps.atc.cat<-((n3.trinuc$atc*n3.trinuc$cat)/(n3.count.aa$i*n3.count.aa$h))/n3.count.diaa$ih

cps.atc.cca<-((n3.trinuc$atc*n3.trinuc$cca)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.atc.ccc<-((n3.trinuc$atc*n3.trinuc$ccc)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.atc.ccg<-((n3.trinuc$atc*n3.trinuc$ccg)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.atc.cct<-((n3.trinuc$atc*n3.trinuc$cct)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip

cps.atc.cga<-((n3.trinuc$atc*n3.trinuc$cga)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.atc.cgc<-((n3.trinuc$atc*n3.trinuc$cgc)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.atc.cgg<-((n3.trinuc$atc*n3.trinuc$cgg)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.atc.cgt<-((n3.trinuc$atc*n3.trinuc$cgt)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir

cps.atc.cta<-((n3.trinuc$atc*n3.trinuc$cta)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.atc.ctc<-((n3.trinuc$atc*n3.trinuc$ctc)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.atc.ctg<-((n3.trinuc$atc*n3.trinuc$ctg)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.atc.ctt<-((n3.trinuc$atc*n3.trinuc$ctt)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il

cps.atc.gaa<-((n3.trinuc$atc*n3.trinuc$gaa)/(n3.count.aa$i*n3.count.aa$e))/n3.count.diaa$ie
cps.atc.gac<-((n3.trinuc$atc*n3.trinuc$gac)/(n3.count.aa$i*n3.count.aa$d))/n3.count.diaa$id
cps.atc.gag<-((n3.trinuc$atc*n3.trinuc$gag)/(n3.count.aa$i*n3.count.aa$e))/n3.count.diaa$ie
cps.atc.gat<-((n3.trinuc$atc*n3.trinuc$gat)/(n3.count.aa$i*n3.count.aa$d))/n3.count.diaa$id

cps.atc.gca<-((n3.trinuc$atc*n3.trinuc$gca)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.atc.gcc<-((n3.trinuc$atc*n3.trinuc$gcc)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.atc.gcg<-((n3.trinuc$atc*n3.trinuc$gcg)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.atc.gct<-((n3.trinuc$atc*n3.trinuc$gct)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia

cps.atc.gga<-((n3.trinuc$atc*n3.trinuc$gga)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.atc.ggc<-((n3.trinuc$atc*n3.trinuc$ggc)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.atc.ggg<-((n3.trinuc$atc*n3.trinuc$ggg)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.atc.ggt<-((n3.trinuc$atc*n3.trinuc$ggt)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig

cps.atc.gta<-((n3.trinuc$atc*n3.trinuc$gta)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.atc.gtc<-((n3.trinuc$atc*n3.trinuc$gtc)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.atc.gtg<-((n3.trinuc$atc*n3.trinuc$gtg)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.atc.gtt<-((n3.trinuc$atc*n3.trinuc$gtt)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv

#Stop codon
#cps.atc.taa<-((n3.trinuc$atc*n3.trinuc$taa)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.atc.tac<-((n3.trinuc$atc*n3.trinuc$tac)/(n3.count.aa$i*n3.count.aa$y))/n3.count.diaa$iy
#Stop codon
#cps.atc.tag<-((n3.trinuc$atc*n3.trinuc$tag)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.atc.tat<-((n3.trinuc$atc*n3.trinuc$tat)/(n3.count.aa$i*n3.count.aa$y))/n3.count.diaa$iy

cps.atc.tca<-((n3.trinuc$atc*n3.trinuc$tca)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.atc.tcc<-((n3.trinuc$atc*n3.trinuc$tcc)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.atc.tcg<-((n3.trinuc$atc*n3.trinuc$tcg)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.atc.tct<-((n3.trinuc$atc*n3.trinuc$tct)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is

#Stop codon
#cps.atc.tga<-((n3.trinuc$atc*n3.trinuc$tga)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.atc.tgc<-((n3.trinuc$atc*n3.trinuc$tgc)/(n3.count.aa$i*n3.count.aa$c))/n3.count.diaa$ic
cps.atc.tgg<-((n3.trinuc$atc*n3.trinuc$tgg)/(n3.count.aa$i*n3.count.aa$w))/n3.count.diaa$iw
cps.atc.tgt<-((n3.trinuc$atc*n3.trinuc$tgt)/(n3.count.aa$i*n3.count.aa$c))/n3.count.diaa$ic

cps.atc.tta<-((n3.trinuc$atc*n3.trinuc$tta)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.atc.ttc<-((n3.trinuc$atc*n3.trinuc$ttc)/(n3.count.aa$i*n3.count.aa$f))/n3.count.diaa$"if"
cps.atc.ttg<-((n3.trinuc$atc*n3.trinuc$ttg)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.atc.ttt<-((n3.trinuc$atc*n3.trinuc$ttt)/(n3.count.aa$i*n3.count.aa$f))/n3.count.diaa$"if"










cps.atg.aaa<-((n3.trinuc$atg*n3.trinuc$aaa)/(n3.count.aa$m*n3.count.aa$k))/n3.count.diaa$mk
cps.atg.aac<-((n3.trinuc$atg*n3.trinuc$aac)/(n3.count.aa$m*n3.count.aa$n))/n3.count.diaa$mn
cps.atg.aag<-((n3.trinuc$atg*n3.trinuc$aag)/(n3.count.aa$m*n3.count.aa$k))/n3.count.diaa$mk
cps.atg.aat<-((n3.trinuc$atg*n3.trinuc$aat)/(n3.count.aa$m*n3.count.aa$n))/n3.count.diaa$mn

cps.atg.aca<-((n3.trinuc$atg*n3.trinuc$aca)/(n3.count.aa$m*n3.count.aa$t))/n3.count.diaa$mt
cps.atg.acc<-((n3.trinuc$atg*n3.trinuc$acc)/(n3.count.aa$m*n3.count.aa$t))/n3.count.diaa$mt
cps.atg.acg<-((n3.trinuc$atg*n3.trinuc$acg)/(n3.count.aa$m*n3.count.aa$t))/n3.count.diaa$mt
cps.atg.act<-((n3.trinuc$atg*n3.trinuc$act)/(n3.count.aa$m*n3.count.aa$t))/n3.count.diaa$mt

cps.atg.aga<-((n3.trinuc$atg*n3.trinuc$aga)/(n3.count.aa$m*n3.count.aa$r))/n3.count.diaa$mr
cps.atg.agc<-((n3.trinuc$atg*n3.trinuc$agc)/(n3.count.aa$m*n3.count.aa$s))/n3.count.diaa$ms
cps.atg.agg<-((n3.trinuc$atg*n3.trinuc$agg)/(n3.count.aa$m*n3.count.aa$r))/n3.count.diaa$mr
cps.atg.agt<-((n3.trinuc$atg*n3.trinuc$agt)/(n3.count.aa$m*n3.count.aa$s))/n3.count.diaa$ms

cps.atg.ata<-((n3.trinuc$atg*n3.trinuc$ata)/(n3.count.aa$m*n3.count.aa$i))/n3.count.diaa$mi
cps.atg.atc<-((n3.trinuc$atg*n3.trinuc$atc)/(n3.count.aa$m*n3.count.aa$i))/n3.count.diaa$mi
cps.atg.atg<-((n3.trinuc$atg*n3.trinuc$atg)/(n3.count.aa$m*n3.count.aa$m))/n3.count.diaa$mm
cps.atg.att<-((n3.trinuc$atg*n3.trinuc$att)/(n3.count.aa$m*n3.count.aa$i))/n3.count.diaa$mi

cps.atg.caa<-((n3.trinuc$atg*n3.trinuc$caa)/(n3.count.aa$m*n3.count.aa$q))/n3.count.diaa$mq
cps.atg.cac<-((n3.trinuc$atg*n3.trinuc$cac)/(n3.count.aa$m*n3.count.aa$h))/n3.count.diaa$mh
cps.atg.cag<-((n3.trinuc$atg*n3.trinuc$cag)/(n3.count.aa$m*n3.count.aa$q))/n3.count.diaa$mq
cps.atg.cat<-((n3.trinuc$atg*n3.trinuc$cat)/(n3.count.aa$m*n3.count.aa$h))/n3.count.diaa$mh

cps.atg.cca<-((n3.trinuc$atg*n3.trinuc$cca)/(n3.count.aa$m*n3.count.aa$p))/n3.count.diaa$mp
cps.atg.ccc<-((n3.trinuc$atg*n3.trinuc$ccc)/(n3.count.aa$m*n3.count.aa$p))/n3.count.diaa$mp
cps.atg.ccg<-((n3.trinuc$atg*n3.trinuc$ccg)/(n3.count.aa$m*n3.count.aa$p))/n3.count.diaa$mp
cps.atg.cct<-((n3.trinuc$atg*n3.trinuc$cct)/(n3.count.aa$m*n3.count.aa$p))/n3.count.diaa$mp

cps.atg.cga<-((n3.trinuc$atg*n3.trinuc$cga)/(n3.count.aa$m*n3.count.aa$r))/n3.count.diaa$mr
cps.atg.cgc<-((n3.trinuc$atg*n3.trinuc$cgc)/(n3.count.aa$m*n3.count.aa$r))/n3.count.diaa$mr
cps.atg.cgg<-((n3.trinuc$atg*n3.trinuc$cgg)/(n3.count.aa$m*n3.count.aa$r))/n3.count.diaa$mr
cps.atg.cgt<-((n3.trinuc$atg*n3.trinuc$cgt)/(n3.count.aa$m*n3.count.aa$r))/n3.count.diaa$mr

cps.atg.cta<-((n3.trinuc$atg*n3.trinuc$cta)/(n3.count.aa$m*n3.count.aa$l))/n3.count.diaa$ml
cps.atg.ctc<-((n3.trinuc$atg*n3.trinuc$ctc)/(n3.count.aa$m*n3.count.aa$l))/n3.count.diaa$ml
cps.atg.ctg<-((n3.trinuc$atg*n3.trinuc$ctg)/(n3.count.aa$m*n3.count.aa$l))/n3.count.diaa$ml
cps.atg.ctt<-((n3.trinuc$atg*n3.trinuc$ctt)/(n3.count.aa$m*n3.count.aa$l))/n3.count.diaa$ml

cps.atg.gaa<-((n3.trinuc$atg*n3.trinuc$gaa)/(n3.count.aa$m*n3.count.aa$e))/n3.count.diaa$me
cps.atg.gac<-((n3.trinuc$atg*n3.trinuc$gac)/(n3.count.aa$m*n3.count.aa$d))/n3.count.diaa$md
cps.atg.gag<-((n3.trinuc$atg*n3.trinuc$gag)/(n3.count.aa$m*n3.count.aa$e))/n3.count.diaa$me
cps.atg.gat<-((n3.trinuc$atg*n3.trinuc$gat)/(n3.count.aa$m*n3.count.aa$d))/n3.count.diaa$md

cps.atg.gca<-((n3.trinuc$atg*n3.trinuc$gca)/(n3.count.aa$m*n3.count.aa$a))/n3.count.diaa$ma
cps.atg.gcc<-((n3.trinuc$atg*n3.trinuc$gcc)/(n3.count.aa$m*n3.count.aa$a))/n3.count.diaa$ma
cps.atg.gcg<-((n3.trinuc$atg*n3.trinuc$gcg)/(n3.count.aa$m*n3.count.aa$a))/n3.count.diaa$ma
cps.atg.gct<-((n3.trinuc$atg*n3.trinuc$gct)/(n3.count.aa$m*n3.count.aa$a))/n3.count.diaa$ma

cps.atg.gga<-((n3.trinuc$atg*n3.trinuc$gga)/(n3.count.aa$m*n3.count.aa$g))/n3.count.diaa$mg
cps.atg.ggc<-((n3.trinuc$atg*n3.trinuc$ggc)/(n3.count.aa$m*n3.count.aa$g))/n3.count.diaa$mg
cps.atg.ggg<-((n3.trinuc$atg*n3.trinuc$ggg)/(n3.count.aa$m*n3.count.aa$g))/n3.count.diaa$mg
cps.atg.ggt<-((n3.trinuc$atg*n3.trinuc$ggt)/(n3.count.aa$m*n3.count.aa$g))/n3.count.diaa$mg

cps.atg.gta<-((n3.trinuc$atg*n3.trinuc$gta)/(n3.count.aa$m*n3.count.aa$v))/n3.count.diaa$mv
cps.atg.gtc<-((n3.trinuc$atg*n3.trinuc$gtc)/(n3.count.aa$m*n3.count.aa$v))/n3.count.diaa$mv
cps.atg.gtg<-((n3.trinuc$atg*n3.trinuc$gtg)/(n3.count.aa$m*n3.count.aa$v))/n3.count.diaa$mv
cps.atg.gtt<-((n3.trinuc$atg*n3.trinuc$gtt)/(n3.count.aa$m*n3.count.aa$v))/n3.count.diaa$mv

#Stop codon
#cps.atg.taa<-((n3.trinuc$atg*n3.trinuc$taa)/(n3.count.aa$m*n3.count.aa$k))/n3.count.diaa$mk
cps.atg.tac<-((n3.trinuc$atg*n3.trinuc$tac)/(n3.count.aa$m*n3.count.aa$y))/n3.count.diaa$my
#Stop codon
#cps.atg.tag<-((n3.trinuc$atg*n3.trinuc$tag)/(n3.count.aa$m*n3.count.aa$k))/n3.count.diaa$mk
cps.atg.tat<-((n3.trinuc$atg*n3.trinuc$tat)/(n3.count.aa$m*n3.count.aa$y))/n3.count.diaa$my

cps.atg.tca<-((n3.trinuc$atg*n3.trinuc$tca)/(n3.count.aa$m*n3.count.aa$s))/n3.count.diaa$ms
cps.atg.tcc<-((n3.trinuc$atg*n3.trinuc$tcc)/(n3.count.aa$m*n3.count.aa$s))/n3.count.diaa$ms
cps.atg.tcg<-((n3.trinuc$atg*n3.trinuc$tcg)/(n3.count.aa$m*n3.count.aa$s))/n3.count.diaa$ms
cps.atg.tct<-((n3.trinuc$atg*n3.trinuc$tct)/(n3.count.aa$m*n3.count.aa$s))/n3.count.diaa$ms

#Stop codon
#cps.atg.tga<-((n3.trinuc$atg*n3.trinuc$tga)/(n3.count.aa$m*n3.count.aa$k))/n3.count.diaa$mk
cps.atg.tgc<-((n3.trinuc$atg*n3.trinuc$tgc)/(n3.count.aa$m*n3.count.aa$c))/n3.count.diaa$mc
cps.atg.tgg<-((n3.trinuc$atg*n3.trinuc$tgg)/(n3.count.aa$m*n3.count.aa$w))/n3.count.diaa$mw
cps.atg.tgt<-((n3.trinuc$atg*n3.trinuc$tgt)/(n3.count.aa$m*n3.count.aa$c))/n3.count.diaa$mc

cps.atg.tta<-((n3.trinuc$atg*n3.trinuc$tta)/(n3.count.aa$m*n3.count.aa$l))/n3.count.diaa$ml
cps.atg.ttc<-((n3.trinuc$atg*n3.trinuc$ttc)/(n3.count.aa$m*n3.count.aa$f))/n3.count.diaa$mf
cps.atg.ttg<-((n3.trinuc$atg*n3.trinuc$ttg)/(n3.count.aa$m*n3.count.aa$l))/n3.count.diaa$ml
cps.atg.ttt<-((n3.trinuc$atg*n3.trinuc$ttt)/(n3.count.aa$m*n3.count.aa$f))/n3.count.diaa$mf








cps.att.aaa<-((n3.trinuc$att*n3.trinuc$aaa)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.att.aac<-((n3.trinuc$att*n3.trinuc$aac)/(n3.count.aa$i*n3.count.aa$n))/n3.count.diaa$"in"
cps.att.aag<-((n3.trinuc$att*n3.trinuc$aag)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.att.aat<-((n3.trinuc$att*n3.trinuc$aat)/(n3.count.aa$i*n3.count.aa$n))/n3.count.diaa$"in"

cps.att.aca<-((n3.trinuc$att*n3.trinuc$aca)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.att.acc<-((n3.trinuc$att*n3.trinuc$acc)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.att.acg<-((n3.trinuc$att*n3.trinuc$acg)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it
cps.att.act<-((n3.trinuc$att*n3.trinuc$act)/(n3.count.aa$i*n3.count.aa$t))/n3.count.diaa$it

cps.att.aga<-((n3.trinuc$att*n3.trinuc$aga)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.att.agc<-((n3.trinuc$att*n3.trinuc$agc)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.att.agg<-((n3.trinuc$att*n3.trinuc$agg)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.att.agt<-((n3.trinuc$att*n3.trinuc$agt)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is

cps.att.ata<-((n3.trinuc$att*n3.trinuc$ata)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii
cps.att.atc<-((n3.trinuc$att*n3.trinuc$atc)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii
cps.att.atg<-((n3.trinuc$att*n3.trinuc$atg)/(n3.count.aa$i*n3.count.aa$m))/n3.count.diaa$im
cps.att.att<-((n3.trinuc$att*n3.trinuc$att)/(n3.count.aa$i*n3.count.aa$i))/n3.count.diaa$ii

cps.att.caa<-((n3.trinuc$att*n3.trinuc$caa)/(n3.count.aa$i*n3.count.aa$q))/n3.count.diaa$iq
cps.att.cac<-((n3.trinuc$att*n3.trinuc$cac)/(n3.count.aa$i*n3.count.aa$h))/n3.count.diaa$ih
cps.att.cag<-((n3.trinuc$att*n3.trinuc$cag)/(n3.count.aa$i*n3.count.aa$q))/n3.count.diaa$iq
cps.att.cat<-((n3.trinuc$att*n3.trinuc$cat)/(n3.count.aa$i*n3.count.aa$h))/n3.count.diaa$ih

cps.att.cca<-((n3.trinuc$att*n3.trinuc$cca)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.att.ccc<-((n3.trinuc$att*n3.trinuc$ccc)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.att.ccg<-((n3.trinuc$att*n3.trinuc$ccg)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip
cps.att.cct<-((n3.trinuc$att*n3.trinuc$cct)/(n3.count.aa$i*n3.count.aa$p))/n3.count.diaa$ip

cps.att.cga<-((n3.trinuc$att*n3.trinuc$cga)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.att.cgc<-((n3.trinuc$att*n3.trinuc$cgc)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.att.cgg<-((n3.trinuc$att*n3.trinuc$cgg)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir
cps.att.cgt<-((n3.trinuc$att*n3.trinuc$cgt)/(n3.count.aa$i*n3.count.aa$r))/n3.count.diaa$ir

cps.att.cta<-((n3.trinuc$att*n3.trinuc$cta)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.att.ctc<-((n3.trinuc$att*n3.trinuc$ctc)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.att.ctg<-((n3.trinuc$att*n3.trinuc$ctg)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.att.ctt<-((n3.trinuc$att*n3.trinuc$ctt)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il

cps.att.gaa<-((n3.trinuc$att*n3.trinuc$gaa)/(n3.count.aa$i*n3.count.aa$e))/n3.count.diaa$ie
cps.att.gac<-((n3.trinuc$att*n3.trinuc$gac)/(n3.count.aa$i*n3.count.aa$d))/n3.count.diaa$id
cps.att.gag<-((n3.trinuc$att*n3.trinuc$gag)/(n3.count.aa$i*n3.count.aa$e))/n3.count.diaa$ie
cps.att.gat<-((n3.trinuc$att*n3.trinuc$gat)/(n3.count.aa$i*n3.count.aa$d))/n3.count.diaa$id

cps.att.gca<-((n3.trinuc$att*n3.trinuc$gca)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.att.gcc<-((n3.trinuc$att*n3.trinuc$gcc)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.att.gcg<-((n3.trinuc$att*n3.trinuc$gcg)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia
cps.att.gct<-((n3.trinuc$att*n3.trinuc$gct)/(n3.count.aa$i*n3.count.aa$a))/n3.count.diaa$ia

cps.att.gga<-((n3.trinuc$att*n3.trinuc$gga)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.att.ggc<-((n3.trinuc$att*n3.trinuc$ggc)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.att.ggg<-((n3.trinuc$att*n3.trinuc$ggg)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig
cps.att.ggt<-((n3.trinuc$att*n3.trinuc$ggt)/(n3.count.aa$i*n3.count.aa$g))/n3.count.diaa$ig

cps.att.gta<-((n3.trinuc$att*n3.trinuc$gta)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.att.gtc<-((n3.trinuc$att*n3.trinuc$gtc)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.att.gtg<-((n3.trinuc$att*n3.trinuc$gtg)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv
cps.att.gtt<-((n3.trinuc$att*n3.trinuc$gtt)/(n3.count.aa$i*n3.count.aa$v))/n3.count.diaa$iv

#Stop codon
#cps.att.taa<-((n3.trinuc$att*n3.trinuc$taa)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.att.tac<-((n3.trinuc$att*n3.trinuc$tac)/(n3.count.aa$i*n3.count.aa$y))/n3.count.diaa$iy
#Stop codon
#cps.att.tag<-((n3.trinuc$att*n3.trinuc$tag)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.att.tat<-((n3.trinuc$att*n3.trinuc$tat)/(n3.count.aa$i*n3.count.aa$y))/n3.count.diaa$iy

cps.att.tca<-((n3.trinuc$att*n3.trinuc$tca)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.att.tcc<-((n3.trinuc$att*n3.trinuc$tcc)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.att.tcg<-((n3.trinuc$att*n3.trinuc$tcg)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is
cps.att.tct<-((n3.trinuc$att*n3.trinuc$tct)/(n3.count.aa$i*n3.count.aa$s))/n3.count.diaa$is

#Stop codon
#cps.att.tga<-((n3.trinuc$att*n3.trinuc$tga)/(n3.count.aa$i*n3.count.aa$k))/n3.count.diaa$ik
cps.att.tgc<-((n3.trinuc$att*n3.trinuc$tgc)/(n3.count.aa$i*n3.count.aa$c))/n3.count.diaa$ic
cps.att.tgg<-((n3.trinuc$att*n3.trinuc$tgg)/(n3.count.aa$i*n3.count.aa$w))/n3.count.diaa$iw
cps.att.tgt<-((n3.trinuc$att*n3.trinuc$tgt)/(n3.count.aa$i*n3.count.aa$c))/n3.count.diaa$ic

cps.att.tta<-((n3.trinuc$att*n3.trinuc$tta)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.att.ttc<-((n3.trinuc$att*n3.trinuc$ttc)/(n3.count.aa$i*n3.count.aa$f))/n3.count.diaa$"if"
cps.att.ttg<-((n3.trinuc$att*n3.trinuc$ttg)/(n3.count.aa$i*n3.count.aa$l))/n3.count.diaa$il
cps.att.ttt<-((n3.trinuc$att*n3.trinuc$ttt)/(n3.count.aa$i*n3.count.aa$f))/n3.count.diaa$"if"








cps.caa.aaa<-((n3.trinuc$caa*n3.trinuc$aaa)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.caa.aac<-((n3.trinuc$caa*n3.trinuc$aac)/(n3.count.aa$q*n3.count.aa$n))/n3.count.diaa$qn
cps.caa.aag<-((n3.trinuc$caa*n3.trinuc$aag)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.caa.aat<-((n3.trinuc$caa*n3.trinuc$aat)/(n3.count.aa$q*n3.count.aa$n))/n3.count.diaa$qn

cps.caa.aca<-((n3.trinuc$caa*n3.trinuc$aca)/(n3.count.aa$q*n3.count.aa$t))/n3.count.diaa$qt
cps.caa.acc<-((n3.trinuc$caa*n3.trinuc$acc)/(n3.count.aa$q*n3.count.aa$t))/n3.count.diaa$qt
cps.caa.acg<-((n3.trinuc$caa*n3.trinuc$acg)/(n3.count.aa$q*n3.count.aa$t))/n3.count.diaa$qt
cps.caa.act<-((n3.trinuc$caa*n3.trinuc$act)/(n3.count.aa$q*n3.count.aa$t))/n3.count.diaa$qt

cps.caa.aga<-((n3.trinuc$caa*n3.trinuc$aga)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.caa.agc<-((n3.trinuc$caa*n3.trinuc$agc)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs
cps.caa.agg<-((n3.trinuc$caa*n3.trinuc$agg)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.caa.agt<-((n3.trinuc$caa*n3.trinuc$agt)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs

cps.caa.ata<-((n3.trinuc$caa*n3.trinuc$ata)/(n3.count.aa$q*n3.count.aa$i))/n3.count.diaa$qi
cps.caa.atc<-((n3.trinuc$caa*n3.trinuc$atc)/(n3.count.aa$q*n3.count.aa$i))/n3.count.diaa$qi
cps.caa.atg<-((n3.trinuc$caa*n3.trinuc$atg)/(n3.count.aa$q*n3.count.aa$m))/n3.count.diaa$qm
cps.caa.att<-((n3.trinuc$caa*n3.trinuc$att)/(n3.count.aa$q*n3.count.aa$i))/n3.count.diaa$qi

cps.caa.caa<-((n3.trinuc$caa*n3.trinuc$caa)/(n3.count.aa$q*n3.count.aa$q))/n3.count.diaa$qq
cps.caa.cac<-((n3.trinuc$caa*n3.trinuc$cac)/(n3.count.aa$q*n3.count.aa$h))/n3.count.diaa$qh
cps.caa.cag<-((n3.trinuc$caa*n3.trinuc$cag)/(n3.count.aa$q*n3.count.aa$q))/n3.count.diaa$qq
cps.caa.cat<-((n3.trinuc$caa*n3.trinuc$cat)/(n3.count.aa$q*n3.count.aa$h))/n3.count.diaa$qh

cps.caa.cca<-((n3.trinuc$caa*n3.trinuc$cca)/(n3.count.aa$q*n3.count.aa$p))/n3.count.diaa$qp
cps.caa.ccc<-((n3.trinuc$caa*n3.trinuc$ccc)/(n3.count.aa$q*n3.count.aa$p))/n3.count.diaa$qp
cps.caa.ccg<-((n3.trinuc$caa*n3.trinuc$ccg)/(n3.count.aa$q*n3.count.aa$p))/n3.count.diaa$qp
cps.caa.cct<-((n3.trinuc$caa*n3.trinuc$cct)/(n3.count.aa$q*n3.count.aa$p))/n3.count.diaa$qp

cps.caa.cga<-((n3.trinuc$caa*n3.trinuc$cga)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.caa.cgc<-((n3.trinuc$caa*n3.trinuc$cgc)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.caa.cgg<-((n3.trinuc$caa*n3.trinuc$cgg)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.caa.cgt<-((n3.trinuc$caa*n3.trinuc$cgt)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr

cps.caa.cta<-((n3.trinuc$caa*n3.trinuc$cta)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.caa.ctc<-((n3.trinuc$caa*n3.trinuc$ctc)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.caa.ctg<-((n3.trinuc$caa*n3.trinuc$ctg)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.caa.ctt<-((n3.trinuc$caa*n3.trinuc$ctt)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql

cps.caa.gaa<-((n3.trinuc$caa*n3.trinuc$gaa)/(n3.count.aa$q*n3.count.aa$e))/n3.count.diaa$qe
cps.caa.gac<-((n3.trinuc$caa*n3.trinuc$gac)/(n3.count.aa$q*n3.count.aa$d))/n3.count.diaa$qd
cps.caa.gag<-((n3.trinuc$caa*n3.trinuc$gag)/(n3.count.aa$q*n3.count.aa$e))/n3.count.diaa$qe
cps.caa.gat<-((n3.trinuc$caa*n3.trinuc$gat)/(n3.count.aa$q*n3.count.aa$d))/n3.count.diaa$qd

cps.caa.gca<-((n3.trinuc$caa*n3.trinuc$gca)/(n3.count.aa$q*n3.count.aa$a))/n3.count.diaa$qa
cps.caa.gcc<-((n3.trinuc$caa*n3.trinuc$gcc)/(n3.count.aa$q*n3.count.aa$a))/n3.count.diaa$qa
cps.caa.gcg<-((n3.trinuc$caa*n3.trinuc$gcg)/(n3.count.aa$q*n3.count.aa$a))/n3.count.diaa$qa
cps.caa.gct<-((n3.trinuc$caa*n3.trinuc$gct)/(n3.count.aa$q*n3.count.aa$a))/n3.count.diaa$qa

cps.caa.gga<-((n3.trinuc$caa*n3.trinuc$gga)/(n3.count.aa$q*n3.count.aa$g))/n3.count.diaa$qg
cps.caa.ggc<-((n3.trinuc$caa*n3.trinuc$ggc)/(n3.count.aa$q*n3.count.aa$g))/n3.count.diaa$qg
cps.caa.ggg<-((n3.trinuc$caa*n3.trinuc$ggg)/(n3.count.aa$q*n3.count.aa$g))/n3.count.diaa$qg
cps.caa.ggt<-((n3.trinuc$caa*n3.trinuc$ggt)/(n3.count.aa$q*n3.count.aa$g))/n3.count.diaa$qg

cps.caa.gta<-((n3.trinuc$caa*n3.trinuc$gta)/(n3.count.aa$q*n3.count.aa$v))/n3.count.diaa$qv
cps.caa.gtc<-((n3.trinuc$caa*n3.trinuc$gtc)/(n3.count.aa$q*n3.count.aa$v))/n3.count.diaa$qv
cps.caa.gtg<-((n3.trinuc$caa*n3.trinuc$gtg)/(n3.count.aa$q*n3.count.aa$v))/n3.count.diaa$qv
cps.caa.gtt<-((n3.trinuc$caa*n3.trinuc$gtt)/(n3.count.aa$q*n3.count.aa$v))/n3.count.diaa$qv

#Stop codon
#cps.caa.taa<-((n3.trinuc$caa*n3.trinuc$taa)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.caa.tac<-((n3.trinuc$caa*n3.trinuc$tac)/(n3.count.aa$q*n3.count.aa$y))/n3.count.diaa$qy
#Stop codon
#cps.caa.tag<-((n3.trinuc$caa*n3.trinuc$tag)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.caa.tat<-((n3.trinuc$caa*n3.trinuc$tat)/(n3.count.aa$q*n3.count.aa$y))/n3.count.diaa$qy

cps.caa.tca<-((n3.trinuc$caa*n3.trinuc$tca)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs
cps.caa.tcc<-((n3.trinuc$caa*n3.trinuc$tcc)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs
cps.caa.tcg<-((n3.trinuc$caa*n3.trinuc$tcg)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs
cps.caa.tct<-((n3.trinuc$caa*n3.trinuc$tct)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs

#Stop codon
#cps.caa.tga<-((n3.trinuc$caa*n3.trinuc$tga)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.caa.tgc<-((n3.trinuc$caa*n3.trinuc$tgc)/(n3.count.aa$q*n3.count.aa$c))/n3.count.diaa$qc
cps.caa.tgg<-((n3.trinuc$caa*n3.trinuc$tgg)/(n3.count.aa$q*n3.count.aa$w))/n3.count.diaa$qw
cps.caa.tgt<-((n3.trinuc$caa*n3.trinuc$tgt)/(n3.count.aa$q*n3.count.aa$c))/n3.count.diaa$qc

cps.caa.tta<-((n3.trinuc$caa*n3.trinuc$tta)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.caa.ttc<-((n3.trinuc$caa*n3.trinuc$ttc)/(n3.count.aa$q*n3.count.aa$f))/n3.count.diaa$qf
cps.caa.ttg<-((n3.trinuc$caa*n3.trinuc$ttg)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.caa.ttt<-((n3.trinuc$caa*n3.trinuc$ttt)/(n3.count.aa$q*n3.count.aa$f))/n3.count.diaa$qf








cps.cac.aaa<-((n3.trinuc$cac*n3.trinuc$aaa)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cac.aac<-((n3.trinuc$cac*n3.trinuc$aac)/(n3.count.aa$h*n3.count.aa$n))/n3.count.diaa$hn
cps.cac.aag<-((n3.trinuc$cac*n3.trinuc$aag)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cac.aat<-((n3.trinuc$cac*n3.trinuc$aat)/(n3.count.aa$h*n3.count.aa$n))/n3.count.diaa$hn

cps.cac.aca<-((n3.trinuc$cac*n3.trinuc$aca)/(n3.count.aa$h*n3.count.aa$t))/n3.count.diaa$ht
cps.cac.acc<-((n3.trinuc$cac*n3.trinuc$acc)/(n3.count.aa$h*n3.count.aa$t))/n3.count.diaa$ht
cps.cac.acg<-((n3.trinuc$cac*n3.trinuc$acg)/(n3.count.aa$h*n3.count.aa$t))/n3.count.diaa$ht
cps.cac.act<-((n3.trinuc$cac*n3.trinuc$act)/(n3.count.aa$h*n3.count.aa$t))/n3.count.diaa$ht

cps.cac.aga<-((n3.trinuc$cac*n3.trinuc$aga)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cac.agc<-((n3.trinuc$cac*n3.trinuc$agc)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs
cps.cac.agg<-((n3.trinuc$cac*n3.trinuc$agg)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cac.agt<-((n3.trinuc$cac*n3.trinuc$agt)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs

cps.cac.ata<-((n3.trinuc$cac*n3.trinuc$ata)/(n3.count.aa$h*n3.count.aa$i))/n3.count.diaa$hi
cps.cac.atc<-((n3.trinuc$cac*n3.trinuc$atc)/(n3.count.aa$h*n3.count.aa$i))/n3.count.diaa$hi
cps.cac.atg<-((n3.trinuc$cac*n3.trinuc$atg)/(n3.count.aa$h*n3.count.aa$m))/n3.count.diaa$hm
cps.cac.att<-((n3.trinuc$cac*n3.trinuc$att)/(n3.count.aa$h*n3.count.aa$i))/n3.count.diaa$hi

cps.cac.caa<-((n3.trinuc$cac*n3.trinuc$caa)/(n3.count.aa$h*n3.count.aa$q))/n3.count.diaa$hq
cps.cac.cac<-((n3.trinuc$cac*n3.trinuc$cac)/(n3.count.aa$h*n3.count.aa$h))/n3.count.diaa$hh
cps.cac.cag<-((n3.trinuc$cac*n3.trinuc$cag)/(n3.count.aa$h*n3.count.aa$q))/n3.count.diaa$hq
cps.cac.cat<-((n3.trinuc$cac*n3.trinuc$cat)/(n3.count.aa$h*n3.count.aa$h))/n3.count.diaa$hh

cps.cac.cca<-((n3.trinuc$cac*n3.trinuc$cca)/(n3.count.aa$h*n3.count.aa$p))/n3.count.diaa$hp
cps.cac.ccc<-((n3.trinuc$cac*n3.trinuc$ccc)/(n3.count.aa$h*n3.count.aa$p))/n3.count.diaa$hp
cps.cac.ccg<-((n3.trinuc$cac*n3.trinuc$ccg)/(n3.count.aa$h*n3.count.aa$p))/n3.count.diaa$hp
cps.cac.cct<-((n3.trinuc$cac*n3.trinuc$cct)/(n3.count.aa$h*n3.count.aa$p))/n3.count.diaa$hp

cps.cac.cga<-((n3.trinuc$cac*n3.trinuc$cga)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cac.cgc<-((n3.trinuc$cac*n3.trinuc$cgc)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cac.cgg<-((n3.trinuc$cac*n3.trinuc$cgg)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cac.cgt<-((n3.trinuc$cac*n3.trinuc$cgt)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr

cps.cac.cta<-((n3.trinuc$cac*n3.trinuc$cta)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cac.ctc<-((n3.trinuc$cac*n3.trinuc$ctc)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cac.ctg<-((n3.trinuc$cac*n3.trinuc$ctg)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cac.ctt<-((n3.trinuc$cac*n3.trinuc$ctt)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl

cps.cac.gaa<-((n3.trinuc$cac*n3.trinuc$gaa)/(n3.count.aa$h*n3.count.aa$e))/n3.count.diaa$he
cps.cac.gac<-((n3.trinuc$cac*n3.trinuc$gac)/(n3.count.aa$h*n3.count.aa$d))/n3.count.diaa$hd
cps.cac.gag<-((n3.trinuc$cac*n3.trinuc$gag)/(n3.count.aa$h*n3.count.aa$e))/n3.count.diaa$he
cps.cac.gat<-((n3.trinuc$cac*n3.trinuc$gat)/(n3.count.aa$h*n3.count.aa$d))/n3.count.diaa$hd

cps.cac.gca<-((n3.trinuc$cac*n3.trinuc$gca)/(n3.count.aa$h*n3.count.aa$a))/n3.count.diaa$ha
cps.cac.gcc<-((n3.trinuc$cac*n3.trinuc$gcc)/(n3.count.aa$h*n3.count.aa$a))/n3.count.diaa$ha
cps.cac.gcg<-((n3.trinuc$cac*n3.trinuc$gcg)/(n3.count.aa$h*n3.count.aa$a))/n3.count.diaa$ha
cps.cac.gct<-((n3.trinuc$cac*n3.trinuc$gct)/(n3.count.aa$h*n3.count.aa$a))/n3.count.diaa$ha

cps.cac.gga<-((n3.trinuc$cac*n3.trinuc$gga)/(n3.count.aa$h*n3.count.aa$g))/n3.count.diaa$hg
cps.cac.ggc<-((n3.trinuc$cac*n3.trinuc$ggc)/(n3.count.aa$h*n3.count.aa$g))/n3.count.diaa$hg
cps.cac.ggg<-((n3.trinuc$cac*n3.trinuc$ggg)/(n3.count.aa$h*n3.count.aa$g))/n3.count.diaa$hg
cps.cac.ggt<-((n3.trinuc$cac*n3.trinuc$ggt)/(n3.count.aa$h*n3.count.aa$g))/n3.count.diaa$hg

cps.cac.gta<-((n3.trinuc$cac*n3.trinuc$gta)/(n3.count.aa$h*n3.count.aa$v))/n3.count.diaa$hv
cps.cac.gtc<-((n3.trinuc$cac*n3.trinuc$gtc)/(n3.count.aa$h*n3.count.aa$v))/n3.count.diaa$hv
cps.cac.gtg<-((n3.trinuc$cac*n3.trinuc$gtg)/(n3.count.aa$h*n3.count.aa$v))/n3.count.diaa$hv
cps.cac.gtt<-((n3.trinuc$cac*n3.trinuc$gtt)/(n3.count.aa$h*n3.count.aa$v))/n3.count.diaa$hv

#Stop codon
#cps.cac.taa<-((n3.trinuc$cac*n3.trinuc$taa)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cac.tac<-((n3.trinuc$cac*n3.trinuc$tac)/(n3.count.aa$h*n3.count.aa$y))/n3.count.diaa$hy
#Stop codon
#cps.cac.tag<-((n3.trinuc$cac*n3.trinuc$tag)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cac.tat<-((n3.trinuc$cac*n3.trinuc$tat)/(n3.count.aa$h*n3.count.aa$y))/n3.count.diaa$hy

cps.cac.tca<-((n3.trinuc$cac*n3.trinuc$tca)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs
cps.cac.tcc<-((n3.trinuc$cac*n3.trinuc$tcc)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs
cps.cac.tcg<-((n3.trinuc$cac*n3.trinuc$tcg)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs
cps.cac.tct<-((n3.trinuc$cac*n3.trinuc$tct)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs

#Stop codon
#cps.cac.tga<-((n3.trinuc$cac*n3.trinuc$tga)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cac.tgc<-((n3.trinuc$cac*n3.trinuc$tgc)/(n3.count.aa$h*n3.count.aa$c))/n3.count.diaa$hc
cps.cac.tgg<-((n3.trinuc$cac*n3.trinuc$tgg)/(n3.count.aa$h*n3.count.aa$w))/n3.count.diaa$hw
cps.cac.tgt<-((n3.trinuc$cac*n3.trinuc$tgt)/(n3.count.aa$h*n3.count.aa$c))/n3.count.diaa$hc

cps.cac.tta<-((n3.trinuc$cac*n3.trinuc$tta)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cac.ttc<-((n3.trinuc$cac*n3.trinuc$ttc)/(n3.count.aa$h*n3.count.aa$f))/n3.count.diaa$hf
cps.cac.ttg<-((n3.trinuc$cac*n3.trinuc$ttg)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cac.ttt<-((n3.trinuc$cac*n3.trinuc$ttt)/(n3.count.aa$h*n3.count.aa$f))/n3.count.diaa$hf








cps.cag.aaa<-((n3.trinuc$cag*n3.trinuc$aaa)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.cag.aac<-((n3.trinuc$cag*n3.trinuc$aac)/(n3.count.aa$q*n3.count.aa$n))/n3.count.diaa$qn
cps.cag.aag<-((n3.trinuc$cag*n3.trinuc$aag)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.cag.aat<-((n3.trinuc$cag*n3.trinuc$aat)/(n3.count.aa$q*n3.count.aa$n))/n3.count.diaa$qn

cps.cag.aca<-((n3.trinuc$cag*n3.trinuc$aca)/(n3.count.aa$q*n3.count.aa$t))/n3.count.diaa$qt
cps.cag.acc<-((n3.trinuc$cag*n3.trinuc$acc)/(n3.count.aa$q*n3.count.aa$t))/n3.count.diaa$qt
cps.cag.acg<-((n3.trinuc$cag*n3.trinuc$acg)/(n3.count.aa$q*n3.count.aa$t))/n3.count.diaa$qt
cps.cag.act<-((n3.trinuc$cag*n3.trinuc$act)/(n3.count.aa$q*n3.count.aa$t))/n3.count.diaa$qt

cps.cag.aga<-((n3.trinuc$cag*n3.trinuc$aga)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.cag.agc<-((n3.trinuc$cag*n3.trinuc$agc)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs
cps.cag.agg<-((n3.trinuc$cag*n3.trinuc$agg)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.cag.agt<-((n3.trinuc$cag*n3.trinuc$agt)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs

cps.cag.ata<-((n3.trinuc$cag*n3.trinuc$ata)/(n3.count.aa$q*n3.count.aa$i))/n3.count.diaa$qi
cps.cag.atc<-((n3.trinuc$cag*n3.trinuc$atc)/(n3.count.aa$q*n3.count.aa$i))/n3.count.diaa$qi
cps.cag.atg<-((n3.trinuc$cag*n3.trinuc$atg)/(n3.count.aa$q*n3.count.aa$m))/n3.count.diaa$qm
cps.cag.att<-((n3.trinuc$cag*n3.trinuc$att)/(n3.count.aa$q*n3.count.aa$i))/n3.count.diaa$qi

cps.cag.caa<-((n3.trinuc$cag*n3.trinuc$caa)/(n3.count.aa$q*n3.count.aa$q))/n3.count.diaa$qq
cps.cag.cac<-((n3.trinuc$cag*n3.trinuc$cac)/(n3.count.aa$q*n3.count.aa$h))/n3.count.diaa$qh
cps.cag.cag<-((n3.trinuc$cag*n3.trinuc$cag)/(n3.count.aa$q*n3.count.aa$q))/n3.count.diaa$qq
cps.cag.cat<-((n3.trinuc$cag*n3.trinuc$cat)/(n3.count.aa$q*n3.count.aa$h))/n3.count.diaa$qh

cps.cag.cca<-((n3.trinuc$cag*n3.trinuc$cca)/(n3.count.aa$q*n3.count.aa$p))/n3.count.diaa$qp
cps.cag.ccc<-((n3.trinuc$cag*n3.trinuc$ccc)/(n3.count.aa$q*n3.count.aa$p))/n3.count.diaa$qp
cps.cag.ccg<-((n3.trinuc$cag*n3.trinuc$ccg)/(n3.count.aa$q*n3.count.aa$p))/n3.count.diaa$qp
cps.cag.cct<-((n3.trinuc$cag*n3.trinuc$cct)/(n3.count.aa$q*n3.count.aa$p))/n3.count.diaa$qp

cps.cag.cga<-((n3.trinuc$cag*n3.trinuc$cga)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.cag.cgc<-((n3.trinuc$cag*n3.trinuc$cgc)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.cag.cgg<-((n3.trinuc$cag*n3.trinuc$cgg)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr
cps.cag.cgt<-((n3.trinuc$cag*n3.trinuc$cgt)/(n3.count.aa$q*n3.count.aa$r))/n3.count.diaa$qr

cps.cag.cta<-((n3.trinuc$cag*n3.trinuc$cta)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.cag.ctc<-((n3.trinuc$cag*n3.trinuc$ctc)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.cag.ctg<-((n3.trinuc$cag*n3.trinuc$ctg)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.cag.ctt<-((n3.trinuc$cag*n3.trinuc$ctt)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql

cps.cag.gaa<-((n3.trinuc$cag*n3.trinuc$gaa)/(n3.count.aa$q*n3.count.aa$e))/n3.count.diaa$qe
cps.cag.gac<-((n3.trinuc$cag*n3.trinuc$gac)/(n3.count.aa$q*n3.count.aa$d))/n3.count.diaa$qd
cps.cag.gag<-((n3.trinuc$cag*n3.trinuc$gag)/(n3.count.aa$q*n3.count.aa$e))/n3.count.diaa$qe
cps.cag.gat<-((n3.trinuc$cag*n3.trinuc$gat)/(n3.count.aa$q*n3.count.aa$d))/n3.count.diaa$qd

cps.cag.gca<-((n3.trinuc$cag*n3.trinuc$gca)/(n3.count.aa$q*n3.count.aa$a))/n3.count.diaa$qa
cps.cag.gcc<-((n3.trinuc$cag*n3.trinuc$gcc)/(n3.count.aa$q*n3.count.aa$a))/n3.count.diaa$qa
cps.cag.gcg<-((n3.trinuc$cag*n3.trinuc$gcg)/(n3.count.aa$q*n3.count.aa$a))/n3.count.diaa$qa
cps.cag.gct<-((n3.trinuc$cag*n3.trinuc$gct)/(n3.count.aa$q*n3.count.aa$a))/n3.count.diaa$qa

cps.cag.gga<-((n3.trinuc$cag*n3.trinuc$gga)/(n3.count.aa$q*n3.count.aa$g))/n3.count.diaa$qg
cps.cag.ggc<-((n3.trinuc$cag*n3.trinuc$ggc)/(n3.count.aa$q*n3.count.aa$g))/n3.count.diaa$qg
cps.cag.ggg<-((n3.trinuc$cag*n3.trinuc$ggg)/(n3.count.aa$q*n3.count.aa$g))/n3.count.diaa$qg
cps.cag.ggt<-((n3.trinuc$cag*n3.trinuc$ggt)/(n3.count.aa$q*n3.count.aa$g))/n3.count.diaa$qg

cps.cag.gta<-((n3.trinuc$cag*n3.trinuc$gta)/(n3.count.aa$q*n3.count.aa$v))/n3.count.diaa$qv
cps.cag.gtc<-((n3.trinuc$cag*n3.trinuc$gtc)/(n3.count.aa$q*n3.count.aa$v))/n3.count.diaa$qv
cps.cag.gtg<-((n3.trinuc$cag*n3.trinuc$gtg)/(n3.count.aa$q*n3.count.aa$v))/n3.count.diaa$qv
cps.cag.gtt<-((n3.trinuc$cag*n3.trinuc$gtt)/(n3.count.aa$q*n3.count.aa$v))/n3.count.diaa$qv

#Stop codon
#cps.cag.taa<-((n3.trinuc$cag*n3.trinuc$taa)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.cag.tac<-((n3.trinuc$cag*n3.trinuc$tac)/(n3.count.aa$q*n3.count.aa$y))/n3.count.diaa$qy
#Stop codon
#cps.cag.tag<-((n3.trinuc$cag*n3.trinuc$tag)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.cag.tat<-((n3.trinuc$cag*n3.trinuc$tat)/(n3.count.aa$q*n3.count.aa$y))/n3.count.diaa$qy

cps.cag.tca<-((n3.trinuc$cag*n3.trinuc$tca)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs
cps.cag.tcc<-((n3.trinuc$cag*n3.trinuc$tcc)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs
cps.cag.tcg<-((n3.trinuc$cag*n3.trinuc$tcg)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs
cps.cag.tct<-((n3.trinuc$cag*n3.trinuc$tct)/(n3.count.aa$q*n3.count.aa$s))/n3.count.diaa$qs

#Stop codon
#cps.cag.tga<-((n3.trinuc$cag*n3.trinuc$tga)/(n3.count.aa$q*n3.count.aa$k))/n3.count.diaa$qk
cps.cag.tgc<-((n3.trinuc$cag*n3.trinuc$tgc)/(n3.count.aa$q*n3.count.aa$c))/n3.count.diaa$qc
cps.cag.tgg<-((n3.trinuc$cag*n3.trinuc$tgg)/(n3.count.aa$q*n3.count.aa$w))/n3.count.diaa$qw
cps.cag.tgt<-((n3.trinuc$cag*n3.trinuc$tgt)/(n3.count.aa$q*n3.count.aa$c))/n3.count.diaa$qc

cps.cag.tta<-((n3.trinuc$cag*n3.trinuc$tta)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.cag.ttc<-((n3.trinuc$cag*n3.trinuc$ttc)/(n3.count.aa$q*n3.count.aa$f))/n3.count.diaa$qf
cps.cag.ttg<-((n3.trinuc$cag*n3.trinuc$ttg)/(n3.count.aa$q*n3.count.aa$l))/n3.count.diaa$ql
cps.cag.ttt<-((n3.trinuc$cag*n3.trinuc$ttt)/(n3.count.aa$q*n3.count.aa$f))/n3.count.diaa$qf








cps.cat.aaa<-((n3.trinuc$cat*n3.trinuc$aaa)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cat.aac<-((n3.trinuc$cat*n3.trinuc$aac)/(n3.count.aa$h*n3.count.aa$n))/n3.count.diaa$hn
cps.cat.aag<-((n3.trinuc$cat*n3.trinuc$aag)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cat.aat<-((n3.trinuc$cat*n3.trinuc$aat)/(n3.count.aa$h*n3.count.aa$n))/n3.count.diaa$hn

cps.cat.aca<-((n3.trinuc$cat*n3.trinuc$aca)/(n3.count.aa$h*n3.count.aa$t))/n3.count.diaa$ht
cps.cat.acc<-((n3.trinuc$cat*n3.trinuc$acc)/(n3.count.aa$h*n3.count.aa$t))/n3.count.diaa$ht
cps.cat.acg<-((n3.trinuc$cat*n3.trinuc$acg)/(n3.count.aa$h*n3.count.aa$t))/n3.count.diaa$ht
cps.cat.act<-((n3.trinuc$cat*n3.trinuc$act)/(n3.count.aa$h*n3.count.aa$t))/n3.count.diaa$ht

cps.cat.aga<-((n3.trinuc$cat*n3.trinuc$aga)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cat.agc<-((n3.trinuc$cat*n3.trinuc$agc)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs
cps.cat.agg<-((n3.trinuc$cat*n3.trinuc$agg)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cat.agt<-((n3.trinuc$cat*n3.trinuc$agt)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs

cps.cat.ata<-((n3.trinuc$cat*n3.trinuc$ata)/(n3.count.aa$h*n3.count.aa$i))/n3.count.diaa$hi
cps.cat.atc<-((n3.trinuc$cat*n3.trinuc$atc)/(n3.count.aa$h*n3.count.aa$i))/n3.count.diaa$hi
cps.cat.atg<-((n3.trinuc$cat*n3.trinuc$atg)/(n3.count.aa$h*n3.count.aa$m))/n3.count.diaa$hm
cps.cat.att<-((n3.trinuc$cat*n3.trinuc$att)/(n3.count.aa$h*n3.count.aa$i))/n3.count.diaa$hi

cps.cat.caa<-((n3.trinuc$cat*n3.trinuc$caa)/(n3.count.aa$h*n3.count.aa$q))/n3.count.diaa$hq
cps.cat.cac<-((n3.trinuc$cat*n3.trinuc$cac)/(n3.count.aa$h*n3.count.aa$h))/n3.count.diaa$hh
cps.cat.cag<-((n3.trinuc$cat*n3.trinuc$cag)/(n3.count.aa$h*n3.count.aa$q))/n3.count.diaa$hq
cps.cat.cat<-((n3.trinuc$cat*n3.trinuc$cat)/(n3.count.aa$h*n3.count.aa$h))/n3.count.diaa$hh

cps.cat.cca<-((n3.trinuc$cat*n3.trinuc$cca)/(n3.count.aa$h*n3.count.aa$p))/n3.count.diaa$hp
cps.cat.ccc<-((n3.trinuc$cat*n3.trinuc$ccc)/(n3.count.aa$h*n3.count.aa$p))/n3.count.diaa$hp
cps.cat.ccg<-((n3.trinuc$cat*n3.trinuc$ccg)/(n3.count.aa$h*n3.count.aa$p))/n3.count.diaa$hp
cps.cat.cct<-((n3.trinuc$cat*n3.trinuc$cct)/(n3.count.aa$h*n3.count.aa$p))/n3.count.diaa$hp

cps.cat.cga<-((n3.trinuc$cat*n3.trinuc$cga)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cat.cgc<-((n3.trinuc$cat*n3.trinuc$cgc)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cat.cgg<-((n3.trinuc$cat*n3.trinuc$cgg)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr
cps.cat.cgt<-((n3.trinuc$cat*n3.trinuc$cgt)/(n3.count.aa$h*n3.count.aa$r))/n3.count.diaa$hr

cps.cat.cta<-((n3.trinuc$cat*n3.trinuc$cta)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cat.ctc<-((n3.trinuc$cat*n3.trinuc$ctc)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cat.ctg<-((n3.trinuc$cat*n3.trinuc$ctg)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cat.ctt<-((n3.trinuc$cat*n3.trinuc$ctt)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl

cps.cat.gaa<-((n3.trinuc$cat*n3.trinuc$gaa)/(n3.count.aa$h*n3.count.aa$e))/n3.count.diaa$he
cps.cat.gac<-((n3.trinuc$cat*n3.trinuc$gac)/(n3.count.aa$h*n3.count.aa$d))/n3.count.diaa$hd
cps.cat.gag<-((n3.trinuc$cat*n3.trinuc$gag)/(n3.count.aa$h*n3.count.aa$e))/n3.count.diaa$he
cps.cat.gat<-((n3.trinuc$cat*n3.trinuc$gat)/(n3.count.aa$h*n3.count.aa$d))/n3.count.diaa$hd

cps.cat.gca<-((n3.trinuc$cat*n3.trinuc$gca)/(n3.count.aa$h*n3.count.aa$a))/n3.count.diaa$ha
cps.cat.gcc<-((n3.trinuc$cat*n3.trinuc$gcc)/(n3.count.aa$h*n3.count.aa$a))/n3.count.diaa$ha
cps.cat.gcg<-((n3.trinuc$cat*n3.trinuc$gcg)/(n3.count.aa$h*n3.count.aa$a))/n3.count.diaa$ha
cps.cat.gct<-((n3.trinuc$cat*n3.trinuc$gct)/(n3.count.aa$h*n3.count.aa$a))/n3.count.diaa$ha

cps.cat.gga<-((n3.trinuc$cat*n3.trinuc$gga)/(n3.count.aa$h*n3.count.aa$g))/n3.count.diaa$hg
cps.cat.ggc<-((n3.trinuc$cat*n3.trinuc$ggc)/(n3.count.aa$h*n3.count.aa$g))/n3.count.diaa$hg
cps.cat.ggg<-((n3.trinuc$cat*n3.trinuc$ggg)/(n3.count.aa$h*n3.count.aa$g))/n3.count.diaa$hg
cps.cat.ggt<-((n3.trinuc$cat*n3.trinuc$ggt)/(n3.count.aa$h*n3.count.aa$g))/n3.count.diaa$hg

cps.cat.gta<-((n3.trinuc$cat*n3.trinuc$gta)/(n3.count.aa$h*n3.count.aa$v))/n3.count.diaa$hv
cps.cat.gtc<-((n3.trinuc$cat*n3.trinuc$gtc)/(n3.count.aa$h*n3.count.aa$v))/n3.count.diaa$hv
cps.cat.gtg<-((n3.trinuc$cat*n3.trinuc$gtg)/(n3.count.aa$h*n3.count.aa$v))/n3.count.diaa$hv
cps.cat.gtt<-((n3.trinuc$cat*n3.trinuc$gtt)/(n3.count.aa$h*n3.count.aa$v))/n3.count.diaa$hv

#Stop codon
#cps.cat.taa<-((n3.trinuc$cat*n3.trinuc$taa)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cat.tac<-((n3.trinuc$cat*n3.trinuc$tac)/(n3.count.aa$h*n3.count.aa$y))/n3.count.diaa$hy
#Stop codon
#cps.cat.tag<-((n3.trinuc$cat*n3.trinuc$tag)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cat.tat<-((n3.trinuc$cat*n3.trinuc$tat)/(n3.count.aa$h*n3.count.aa$y))/n3.count.diaa$hy

cps.cat.tca<-((n3.trinuc$cat*n3.trinuc$tca)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs
cps.cat.tcc<-((n3.trinuc$cat*n3.trinuc$tcc)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs
cps.cat.tcg<-((n3.trinuc$cat*n3.trinuc$tcg)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs
cps.cat.tct<-((n3.trinuc$cat*n3.trinuc$tct)/(n3.count.aa$h*n3.count.aa$s))/n3.count.diaa$hs

#Stop codon
#cps.cat.tga<-((n3.trinuc$cat*n3.trinuc$tga)/(n3.count.aa$h*n3.count.aa$k))/n3.count.diaa$hk
cps.cat.tgc<-((n3.trinuc$cat*n3.trinuc$tgc)/(n3.count.aa$h*n3.count.aa$c))/n3.count.diaa$hc
cps.cat.tgg<-((n3.trinuc$cat*n3.trinuc$tgg)/(n3.count.aa$h*n3.count.aa$w))/n3.count.diaa$hw
cps.cat.tgt<-((n3.trinuc$cat*n3.trinuc$tgt)/(n3.count.aa$h*n3.count.aa$c))/n3.count.diaa$hc

cps.cat.tta<-((n3.trinuc$cat*n3.trinuc$tta)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cat.ttc<-((n3.trinuc$cat*n3.trinuc$ttc)/(n3.count.aa$h*n3.count.aa$f))/n3.count.diaa$hf
cps.cat.ttg<-((n3.trinuc$cat*n3.trinuc$ttg)/(n3.count.aa$h*n3.count.aa$l))/n3.count.diaa$hl
cps.cat.ttt<-((n3.trinuc$cat*n3.trinuc$ttt)/(n3.count.aa$h*n3.count.aa$f))/n3.count.diaa$hf








cps.cca.aaa<-((n3.trinuc$cca*n3.trinuc$aaa)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cca.aac<-((n3.trinuc$cca*n3.trinuc$aac)/(n3.count.aa$p*n3.count.aa$n))/n3.count.diaa$pn
cps.cca.aag<-((n3.trinuc$cca*n3.trinuc$aag)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cca.aat<-((n3.trinuc$cca*n3.trinuc$aat)/(n3.count.aa$p*n3.count.aa$n))/n3.count.diaa$pn

cps.cca.aca<-((n3.trinuc$cca*n3.trinuc$aca)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.cca.acc<-((n3.trinuc$cca*n3.trinuc$acc)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.cca.acg<-((n3.trinuc$cca*n3.trinuc$acg)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.cca.act<-((n3.trinuc$cca*n3.trinuc$act)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt

cps.cca.aga<-((n3.trinuc$cca*n3.trinuc$aga)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cca.agc<-((n3.trinuc$cca*n3.trinuc$agc)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.cca.agg<-((n3.trinuc$cca*n3.trinuc$agg)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cca.agt<-((n3.trinuc$cca*n3.trinuc$agt)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps

cps.cca.ata<-((n3.trinuc$cca*n3.trinuc$ata)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi
cps.cca.atc<-((n3.trinuc$cca*n3.trinuc$atc)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi
cps.cca.atg<-((n3.trinuc$cca*n3.trinuc$atg)/(n3.count.aa$p*n3.count.aa$m))/n3.count.diaa$pm
cps.cca.att<-((n3.trinuc$cca*n3.trinuc$att)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi

cps.cca.caa<-((n3.trinuc$cca*n3.trinuc$caa)/(n3.count.aa$p*n3.count.aa$q))/n3.count.diaa$pq
cps.cca.cac<-((n3.trinuc$cca*n3.trinuc$cac)/(n3.count.aa$p*n3.count.aa$h))/n3.count.diaa$ph
cps.cca.cag<-((n3.trinuc$cca*n3.trinuc$cag)/(n3.count.aa$p*n3.count.aa$q))/n3.count.diaa$pq
cps.cca.cat<-((n3.trinuc$cca*n3.trinuc$cat)/(n3.count.aa$p*n3.count.aa$h))/n3.count.diaa$ph

cps.cca.cca<-((n3.trinuc$cca*n3.trinuc$cca)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.cca.ccc<-((n3.trinuc$cca*n3.trinuc$ccc)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.cca.ccg<-((n3.trinuc$cca*n3.trinuc$ccg)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.cca.cct<-((n3.trinuc$cca*n3.trinuc$cct)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp

cps.cca.cga<-((n3.trinuc$cca*n3.trinuc$cga)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cca.cgc<-((n3.trinuc$cca*n3.trinuc$cgc)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cca.cgg<-((n3.trinuc$cca*n3.trinuc$cgg)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cca.cgt<-((n3.trinuc$cca*n3.trinuc$cgt)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr

cps.cca.cta<-((n3.trinuc$cca*n3.trinuc$cta)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cca.ctc<-((n3.trinuc$cca*n3.trinuc$ctc)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cca.ctg<-((n3.trinuc$cca*n3.trinuc$ctg)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cca.ctt<-((n3.trinuc$cca*n3.trinuc$ctt)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl

cps.cca.gaa<-((n3.trinuc$cca*n3.trinuc$gaa)/(n3.count.aa$p*n3.count.aa$e))/n3.count.diaa$pe
cps.cca.gac<-((n3.trinuc$cca*n3.trinuc$gac)/(n3.count.aa$p*n3.count.aa$d))/n3.count.diaa$pd
cps.cca.gag<-((n3.trinuc$cca*n3.trinuc$gag)/(n3.count.aa$p*n3.count.aa$e))/n3.count.diaa$pe
cps.cca.gat<-((n3.trinuc$cca*n3.trinuc$gat)/(n3.count.aa$p*n3.count.aa$d))/n3.count.diaa$pd

cps.cca.gca<-((n3.trinuc$cca*n3.trinuc$gca)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.cca.gcc<-((n3.trinuc$cca*n3.trinuc$gcc)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.cca.gcg<-((n3.trinuc$cca*n3.trinuc$gcg)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.cca.gct<-((n3.trinuc$cca*n3.trinuc$gct)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa

cps.cca.gga<-((n3.trinuc$cca*n3.trinuc$gga)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.cca.ggc<-((n3.trinuc$cca*n3.trinuc$ggc)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.cca.ggg<-((n3.trinuc$cca*n3.trinuc$ggg)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.cca.ggt<-((n3.trinuc$cca*n3.trinuc$ggt)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg

cps.cca.gta<-((n3.trinuc$cca*n3.trinuc$gta)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.cca.gtc<-((n3.trinuc$cca*n3.trinuc$gtc)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.cca.gtg<-((n3.trinuc$cca*n3.trinuc$gtg)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.cca.gtt<-((n3.trinuc$cca*n3.trinuc$gtt)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv

#Stop codon
#cps.cca.taa<-((n3.trinuc$cca*n3.trinuc$taa)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cca.tac<-((n3.trinuc$cca*n3.trinuc$tac)/(n3.count.aa$p*n3.count.aa$y))/n3.count.diaa$py
#Stop codon
#cps.cca.tag<-((n3.trinuc$cca*n3.trinuc$tag)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cca.tat<-((n3.trinuc$cca*n3.trinuc$tat)/(n3.count.aa$p*n3.count.aa$y))/n3.count.diaa$py

cps.cca.tca<-((n3.trinuc$cca*n3.trinuc$tca)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.cca.tcc<-((n3.trinuc$cca*n3.trinuc$tcc)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.cca.tcg<-((n3.trinuc$cca*n3.trinuc$tcg)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.cca.tct<-((n3.trinuc$cca*n3.trinuc$tct)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps

#Stop codon
#cps.cca.tga<-((n3.trinuc$cca*n3.trinuc$tga)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cca.tgc<-((n3.trinuc$cca*n3.trinuc$tgc)/(n3.count.aa$p*n3.count.aa$c))/n3.count.diaa$pc
cps.cca.tgg<-((n3.trinuc$cca*n3.trinuc$tgg)/(n3.count.aa$p*n3.count.aa$w))/n3.count.diaa$pw
cps.cca.tgt<-((n3.trinuc$cca*n3.trinuc$tgt)/(n3.count.aa$p*n3.count.aa$c))/n3.count.diaa$pc

cps.cca.tta<-((n3.trinuc$cca*n3.trinuc$tta)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cca.ttc<-((n3.trinuc$cca*n3.trinuc$ttc)/(n3.count.aa$p*n3.count.aa$f))/n3.count.diaa$pf
cps.cca.ttg<-((n3.trinuc$cca*n3.trinuc$ttg)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cca.ttt<-((n3.trinuc$cca*n3.trinuc$ttt)/(n3.count.aa$p*n3.count.aa$f))/n3.count.diaa$pf








cps.ccc.aaa<-((n3.trinuc$ccc*n3.trinuc$aaa)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccc.aac<-((n3.trinuc$ccc*n3.trinuc$aac)/(n3.count.aa$p*n3.count.aa$n))/n3.count.diaa$pn
cps.ccc.aag<-((n3.trinuc$ccc*n3.trinuc$aag)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccc.aat<-((n3.trinuc$ccc*n3.trinuc$aat)/(n3.count.aa$p*n3.count.aa$n))/n3.count.diaa$pn

cps.ccc.aca<-((n3.trinuc$ccc*n3.trinuc$aca)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.ccc.acc<-((n3.trinuc$ccc*n3.trinuc$acc)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.ccc.acg<-((n3.trinuc$ccc*n3.trinuc$acg)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.ccc.act<-((n3.trinuc$ccc*n3.trinuc$act)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt

cps.ccc.aga<-((n3.trinuc$ccc*n3.trinuc$aga)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccc.agc<-((n3.trinuc$ccc*n3.trinuc$agc)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.ccc.agg<-((n3.trinuc$ccc*n3.trinuc$agg)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccc.agt<-((n3.trinuc$ccc*n3.trinuc$agt)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps

cps.ccc.ata<-((n3.trinuc$ccc*n3.trinuc$ata)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi
cps.ccc.atc<-((n3.trinuc$ccc*n3.trinuc$atc)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi
cps.ccc.atg<-((n3.trinuc$ccc*n3.trinuc$atg)/(n3.count.aa$p*n3.count.aa$m))/n3.count.diaa$pm
cps.ccc.att<-((n3.trinuc$ccc*n3.trinuc$att)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi

cps.ccc.caa<-((n3.trinuc$ccc*n3.trinuc$caa)/(n3.count.aa$p*n3.count.aa$q))/n3.count.diaa$pq
cps.ccc.cac<-((n3.trinuc$ccc*n3.trinuc$cac)/(n3.count.aa$p*n3.count.aa$h))/n3.count.diaa$ph
cps.ccc.cag<-((n3.trinuc$ccc*n3.trinuc$cag)/(n3.count.aa$p*n3.count.aa$q))/n3.count.diaa$pq
cps.ccc.cat<-((n3.trinuc$ccc*n3.trinuc$cat)/(n3.count.aa$p*n3.count.aa$h))/n3.count.diaa$ph

cps.ccc.cca<-((n3.trinuc$ccc*n3.trinuc$cca)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.ccc.ccc<-((n3.trinuc$ccc*n3.trinuc$ccc)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.ccc.ccg<-((n3.trinuc$ccc*n3.trinuc$ccg)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.ccc.cct<-((n3.trinuc$ccc*n3.trinuc$cct)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp

cps.ccc.cga<-((n3.trinuc$ccc*n3.trinuc$cga)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccc.cgc<-((n3.trinuc$ccc*n3.trinuc$cgc)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccc.cgg<-((n3.trinuc$ccc*n3.trinuc$cgg)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccc.cgt<-((n3.trinuc$ccc*n3.trinuc$cgt)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr

cps.ccc.cta<-((n3.trinuc$ccc*n3.trinuc$cta)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccc.ctc<-((n3.trinuc$ccc*n3.trinuc$ctc)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccc.ctg<-((n3.trinuc$ccc*n3.trinuc$ctg)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccc.ctt<-((n3.trinuc$ccc*n3.trinuc$ctt)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl

cps.ccc.gaa<-((n3.trinuc$ccc*n3.trinuc$gaa)/(n3.count.aa$p*n3.count.aa$e))/n3.count.diaa$pe
cps.ccc.gac<-((n3.trinuc$ccc*n3.trinuc$gac)/(n3.count.aa$p*n3.count.aa$d))/n3.count.diaa$pd
cps.ccc.gag<-((n3.trinuc$ccc*n3.trinuc$gag)/(n3.count.aa$p*n3.count.aa$e))/n3.count.diaa$pe
cps.ccc.gat<-((n3.trinuc$ccc*n3.trinuc$gat)/(n3.count.aa$p*n3.count.aa$d))/n3.count.diaa$pd

cps.ccc.gca<-((n3.trinuc$ccc*n3.trinuc$gca)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.ccc.gcc<-((n3.trinuc$ccc*n3.trinuc$gcc)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.ccc.gcg<-((n3.trinuc$ccc*n3.trinuc$gcg)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.ccc.gct<-((n3.trinuc$ccc*n3.trinuc$gct)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa

cps.ccc.gga<-((n3.trinuc$ccc*n3.trinuc$gga)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.ccc.ggc<-((n3.trinuc$ccc*n3.trinuc$ggc)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.ccc.ggg<-((n3.trinuc$ccc*n3.trinuc$ggg)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.ccc.ggt<-((n3.trinuc$ccc*n3.trinuc$ggt)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg

cps.ccc.gta<-((n3.trinuc$ccc*n3.trinuc$gta)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.ccc.gtc<-((n3.trinuc$ccc*n3.trinuc$gtc)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.ccc.gtg<-((n3.trinuc$ccc*n3.trinuc$gtg)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.ccc.gtt<-((n3.trinuc$ccc*n3.trinuc$gtt)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv

#Stop codon
#cps.ccc.taa<-((n3.trinuc$ccc*n3.trinuc$taa)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccc.tac<-((n3.trinuc$ccc*n3.trinuc$tac)/(n3.count.aa$p*n3.count.aa$y))/n3.count.diaa$py
#Stop codon
#cps.ccc.tag<-((n3.trinuc$ccc*n3.trinuc$tag)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccc.tat<-((n3.trinuc$ccc*n3.trinuc$tat)/(n3.count.aa$p*n3.count.aa$y))/n3.count.diaa$py

cps.ccc.tca<-((n3.trinuc$ccc*n3.trinuc$tca)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.ccc.tcc<-((n3.trinuc$ccc*n3.trinuc$tcc)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.ccc.tcg<-((n3.trinuc$ccc*n3.trinuc$tcg)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.ccc.tct<-((n3.trinuc$ccc*n3.trinuc$tct)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps

#Stop codon
#cps.ccc.tga<-((n3.trinuc$ccc*n3.trinuc$tga)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccc.tgc<-((n3.trinuc$ccc*n3.trinuc$tgc)/(n3.count.aa$p*n3.count.aa$c))/n3.count.diaa$pc
cps.ccc.tgg<-((n3.trinuc$ccc*n3.trinuc$tgg)/(n3.count.aa$p*n3.count.aa$w))/n3.count.diaa$pw
cps.ccc.tgt<-((n3.trinuc$ccc*n3.trinuc$tgt)/(n3.count.aa$p*n3.count.aa$c))/n3.count.diaa$pc

cps.ccc.tta<-((n3.trinuc$ccc*n3.trinuc$tta)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccc.ttc<-((n3.trinuc$ccc*n3.trinuc$ttc)/(n3.count.aa$p*n3.count.aa$f))/n3.count.diaa$pf
cps.ccc.ttg<-((n3.trinuc$ccc*n3.trinuc$ttg)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccc.ttt<-((n3.trinuc$ccc*n3.trinuc$ttt)/(n3.count.aa$p*n3.count.aa$f))/n3.count.diaa$pf








cps.ccg.aaa<-((n3.trinuc$ccg*n3.trinuc$aaa)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccg.aac<-((n3.trinuc$ccg*n3.trinuc$aac)/(n3.count.aa$p*n3.count.aa$n))/n3.count.diaa$pn
cps.ccg.aag<-((n3.trinuc$ccg*n3.trinuc$aag)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccg.aat<-((n3.trinuc$ccg*n3.trinuc$aat)/(n3.count.aa$p*n3.count.aa$n))/n3.count.diaa$pn

cps.ccg.aca<-((n3.trinuc$ccg*n3.trinuc$aca)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.ccg.acc<-((n3.trinuc$ccg*n3.trinuc$acc)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.ccg.acg<-((n3.trinuc$ccg*n3.trinuc$acg)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.ccg.act<-((n3.trinuc$ccg*n3.trinuc$act)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt

cps.ccg.aga<-((n3.trinuc$ccg*n3.trinuc$aga)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccg.agc<-((n3.trinuc$ccg*n3.trinuc$agc)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.ccg.agg<-((n3.trinuc$ccg*n3.trinuc$agg)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccg.agt<-((n3.trinuc$ccg*n3.trinuc$agt)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps

cps.ccg.ata<-((n3.trinuc$ccg*n3.trinuc$ata)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi
cps.ccg.atc<-((n3.trinuc$ccg*n3.trinuc$atc)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi
cps.ccg.atg<-((n3.trinuc$ccg*n3.trinuc$atg)/(n3.count.aa$p*n3.count.aa$m))/n3.count.diaa$pm
cps.ccg.att<-((n3.trinuc$ccg*n3.trinuc$att)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi

cps.ccg.caa<-((n3.trinuc$ccg*n3.trinuc$caa)/(n3.count.aa$p*n3.count.aa$q))/n3.count.diaa$pq
cps.ccg.cac<-((n3.trinuc$ccg*n3.trinuc$cac)/(n3.count.aa$p*n3.count.aa$h))/n3.count.diaa$ph
cps.ccg.cag<-((n3.trinuc$ccg*n3.trinuc$cag)/(n3.count.aa$p*n3.count.aa$q))/n3.count.diaa$pq
cps.ccg.cat<-((n3.trinuc$ccg*n3.trinuc$cat)/(n3.count.aa$p*n3.count.aa$h))/n3.count.diaa$ph

cps.ccg.cca<-((n3.trinuc$ccg*n3.trinuc$cca)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.ccg.ccc<-((n3.trinuc$ccg*n3.trinuc$ccc)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.ccg.ccg<-((n3.trinuc$ccg*n3.trinuc$ccg)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.ccg.cct<-((n3.trinuc$ccg*n3.trinuc$cct)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp

cps.ccg.cga<-((n3.trinuc$ccg*n3.trinuc$cga)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccg.cgc<-((n3.trinuc$ccg*n3.trinuc$cgc)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccg.cgg<-((n3.trinuc$ccg*n3.trinuc$cgg)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.ccg.cgt<-((n3.trinuc$ccg*n3.trinuc$cgt)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr

cps.ccg.cta<-((n3.trinuc$ccg*n3.trinuc$cta)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccg.ctc<-((n3.trinuc$ccg*n3.trinuc$ctc)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccg.ctg<-((n3.trinuc$ccg*n3.trinuc$ctg)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccg.ctt<-((n3.trinuc$ccg*n3.trinuc$ctt)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl

cps.ccg.gaa<-((n3.trinuc$ccg*n3.trinuc$gaa)/(n3.count.aa$p*n3.count.aa$e))/n3.count.diaa$pe
cps.ccg.gac<-((n3.trinuc$ccg*n3.trinuc$gac)/(n3.count.aa$p*n3.count.aa$d))/n3.count.diaa$pd
cps.ccg.gag<-((n3.trinuc$ccg*n3.trinuc$gag)/(n3.count.aa$p*n3.count.aa$e))/n3.count.diaa$pe
cps.ccg.gat<-((n3.trinuc$ccg*n3.trinuc$gat)/(n3.count.aa$p*n3.count.aa$d))/n3.count.diaa$pd

cps.ccg.gca<-((n3.trinuc$ccg*n3.trinuc$gca)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.ccg.gcc<-((n3.trinuc$ccg*n3.trinuc$gcc)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.ccg.gcg<-((n3.trinuc$ccg*n3.trinuc$gcg)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.ccg.gct<-((n3.trinuc$ccg*n3.trinuc$gct)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa

cps.ccg.gga<-((n3.trinuc$ccg*n3.trinuc$gga)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.ccg.ggc<-((n3.trinuc$ccg*n3.trinuc$ggc)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.ccg.ggg<-((n3.trinuc$ccg*n3.trinuc$ggg)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.ccg.ggt<-((n3.trinuc$ccg*n3.trinuc$ggt)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg

cps.ccg.gta<-((n3.trinuc$ccg*n3.trinuc$gta)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.ccg.gtc<-((n3.trinuc$ccg*n3.trinuc$gtc)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.ccg.gtg<-((n3.trinuc$ccg*n3.trinuc$gtg)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.ccg.gtt<-((n3.trinuc$ccg*n3.trinuc$gtt)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv

#Stop codon
#cps.ccg.taa<-((n3.trinuc$ccg*n3.trinuc$taa)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccg.tac<-((n3.trinuc$ccg*n3.trinuc$tac)/(n3.count.aa$p*n3.count.aa$y))/n3.count.diaa$py
#Stop codon
#cps.ccg.tag<-((n3.trinuc$ccg*n3.trinuc$tag)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccg.tat<-((n3.trinuc$ccg*n3.trinuc$tat)/(n3.count.aa$p*n3.count.aa$y))/n3.count.diaa$py

cps.ccg.tca<-((n3.trinuc$ccg*n3.trinuc$tca)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.ccg.tcc<-((n3.trinuc$ccg*n3.trinuc$tcc)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.ccg.tcg<-((n3.trinuc$ccg*n3.trinuc$tcg)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.ccg.tct<-((n3.trinuc$ccg*n3.trinuc$tct)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps

#Stop codon
#cps.ccg.tga<-((n3.trinuc$ccg*n3.trinuc$tga)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.ccg.tgc<-((n3.trinuc$ccg*n3.trinuc$tgc)/(n3.count.aa$p*n3.count.aa$c))/n3.count.diaa$pc
cps.ccg.tgg<-((n3.trinuc$ccg*n3.trinuc$tgg)/(n3.count.aa$p*n3.count.aa$w))/n3.count.diaa$pw
cps.ccg.tgt<-((n3.trinuc$ccg*n3.trinuc$tgt)/(n3.count.aa$p*n3.count.aa$c))/n3.count.diaa$pc

cps.ccg.tta<-((n3.trinuc$ccg*n3.trinuc$tta)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccg.ttc<-((n3.trinuc$ccg*n3.trinuc$ttc)/(n3.count.aa$p*n3.count.aa$f))/n3.count.diaa$pf
cps.ccg.ttg<-((n3.trinuc$ccg*n3.trinuc$ttg)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.ccg.ttt<-((n3.trinuc$ccg*n3.trinuc$ttt)/(n3.count.aa$p*n3.count.aa$f))/n3.count.diaa$pf








cps.cct.aaa<-((n3.trinuc$cct*n3.trinuc$aaa)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cct.aac<-((n3.trinuc$cct*n3.trinuc$aac)/(n3.count.aa$p*n3.count.aa$n))/n3.count.diaa$pn
cps.cct.aag<-((n3.trinuc$cct*n3.trinuc$aag)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cct.aat<-((n3.trinuc$cct*n3.trinuc$aat)/(n3.count.aa$p*n3.count.aa$n))/n3.count.diaa$pn

cps.cct.aca<-((n3.trinuc$cct*n3.trinuc$aca)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.cct.acc<-((n3.trinuc$cct*n3.trinuc$acc)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.cct.acg<-((n3.trinuc$cct*n3.trinuc$acg)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt
cps.cct.act<-((n3.trinuc$cct*n3.trinuc$act)/(n3.count.aa$p*n3.count.aa$t))/n3.count.diaa$pt

cps.cct.aga<-((n3.trinuc$cct*n3.trinuc$aga)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cct.agc<-((n3.trinuc$cct*n3.trinuc$agc)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.cct.agg<-((n3.trinuc$cct*n3.trinuc$agg)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cct.agt<-((n3.trinuc$cct*n3.trinuc$agt)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps

cps.cct.ata<-((n3.trinuc$cct*n3.trinuc$ata)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi
cps.cct.atc<-((n3.trinuc$cct*n3.trinuc$atc)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi
cps.cct.atg<-((n3.trinuc$cct*n3.trinuc$atg)/(n3.count.aa$p*n3.count.aa$m))/n3.count.diaa$pm
cps.cct.att<-((n3.trinuc$cct*n3.trinuc$att)/(n3.count.aa$p*n3.count.aa$i))/n3.count.diaa$pi

cps.cct.caa<-((n3.trinuc$cct*n3.trinuc$caa)/(n3.count.aa$p*n3.count.aa$q))/n3.count.diaa$pq
cps.cct.cac<-((n3.trinuc$cct*n3.trinuc$cac)/(n3.count.aa$p*n3.count.aa$h))/n3.count.diaa$ph
cps.cct.cag<-((n3.trinuc$cct*n3.trinuc$cag)/(n3.count.aa$p*n3.count.aa$q))/n3.count.diaa$pq
cps.cct.cat<-((n3.trinuc$cct*n3.trinuc$cat)/(n3.count.aa$p*n3.count.aa$h))/n3.count.diaa$ph

cps.cct.cca<-((n3.trinuc$cct*n3.trinuc$cca)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.cct.ccc<-((n3.trinuc$cct*n3.trinuc$ccc)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.cct.ccg<-((n3.trinuc$cct*n3.trinuc$ccg)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp
cps.cct.cct<-((n3.trinuc$cct*n3.trinuc$cct)/(n3.count.aa$p*n3.count.aa$p))/n3.count.diaa$pp

cps.cct.cga<-((n3.trinuc$cct*n3.trinuc$cga)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cct.cgc<-((n3.trinuc$cct*n3.trinuc$cgc)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cct.cgg<-((n3.trinuc$cct*n3.trinuc$cgg)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr
cps.cct.cgt<-((n3.trinuc$cct*n3.trinuc$cgt)/(n3.count.aa$p*n3.count.aa$r))/n3.count.diaa$pr

cps.cct.cta<-((n3.trinuc$cct*n3.trinuc$cta)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cct.ctc<-((n3.trinuc$cct*n3.trinuc$ctc)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cct.ctg<-((n3.trinuc$cct*n3.trinuc$ctg)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cct.ctt<-((n3.trinuc$cct*n3.trinuc$ctt)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl

cps.cct.gaa<-((n3.trinuc$cct*n3.trinuc$gaa)/(n3.count.aa$p*n3.count.aa$e))/n3.count.diaa$pe
cps.cct.gac<-((n3.trinuc$cct*n3.trinuc$gac)/(n3.count.aa$p*n3.count.aa$d))/n3.count.diaa$pd
cps.cct.gag<-((n3.trinuc$cct*n3.trinuc$gag)/(n3.count.aa$p*n3.count.aa$e))/n3.count.diaa$pe
cps.cct.gat<-((n3.trinuc$cct*n3.trinuc$gat)/(n3.count.aa$p*n3.count.aa$d))/n3.count.diaa$pd

cps.cct.gca<-((n3.trinuc$cct*n3.trinuc$gca)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.cct.gcc<-((n3.trinuc$cct*n3.trinuc$gcc)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.cct.gcg<-((n3.trinuc$cct*n3.trinuc$gcg)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa
cps.cct.gct<-((n3.trinuc$cct*n3.trinuc$gct)/(n3.count.aa$p*n3.count.aa$a))/n3.count.diaa$pa

cps.cct.gga<-((n3.trinuc$cct*n3.trinuc$gga)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.cct.ggc<-((n3.trinuc$cct*n3.trinuc$ggc)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.cct.ggg<-((n3.trinuc$cct*n3.trinuc$ggg)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg
cps.cct.ggt<-((n3.trinuc$cct*n3.trinuc$ggt)/(n3.count.aa$p*n3.count.aa$g))/n3.count.diaa$pg

cps.cct.gta<-((n3.trinuc$cct*n3.trinuc$gta)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.cct.gtc<-((n3.trinuc$cct*n3.trinuc$gtc)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.cct.gtg<-((n3.trinuc$cct*n3.trinuc$gtg)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv
cps.cct.gtt<-((n3.trinuc$cct*n3.trinuc$gtt)/(n3.count.aa$p*n3.count.aa$v))/n3.count.diaa$pv

#Stop codon
#cps.cct.taa<-((n3.trinuc$cct*n3.trinuc$taa)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cct.tac<-((n3.trinuc$cct*n3.trinuc$tac)/(n3.count.aa$p*n3.count.aa$y))/n3.count.diaa$py
#Stop codon
#cps.cct.tag<-((n3.trinuc$cct*n3.trinuc$tag)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cct.tat<-((n3.trinuc$cct*n3.trinuc$tat)/(n3.count.aa$p*n3.count.aa$y))/n3.count.diaa$py

cps.cct.tca<-((n3.trinuc$cct*n3.trinuc$tca)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.cct.tcc<-((n3.trinuc$cct*n3.trinuc$tcc)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.cct.tcg<-((n3.trinuc$cct*n3.trinuc$tcg)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps
cps.cct.tct<-((n3.trinuc$cct*n3.trinuc$tct)/(n3.count.aa$p*n3.count.aa$s))/n3.count.diaa$ps

#Stop codon
#cps.cct.tga<-((n3.trinuc$cct*n3.trinuc$tga)/(n3.count.aa$p*n3.count.aa$k))/n3.count.diaa$pk
cps.cct.tgc<-((n3.trinuc$cct*n3.trinuc$tgc)/(n3.count.aa$p*n3.count.aa$c))/n3.count.diaa$pc
cps.cct.tgg<-((n3.trinuc$cct*n3.trinuc$tgg)/(n3.count.aa$p*n3.count.aa$w))/n3.count.diaa$pw
cps.cct.tgt<-((n3.trinuc$cct*n3.trinuc$tgt)/(n3.count.aa$p*n3.count.aa$c))/n3.count.diaa$pc

cps.cct.tta<-((n3.trinuc$cct*n3.trinuc$tta)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cct.ttc<-((n3.trinuc$cct*n3.trinuc$ttc)/(n3.count.aa$p*n3.count.aa$f))/n3.count.diaa$pf
cps.cct.ttg<-((n3.trinuc$cct*n3.trinuc$ttg)/(n3.count.aa$p*n3.count.aa$l))/n3.count.diaa$pl
cps.cct.ttt<-((n3.trinuc$cct*n3.trinuc$ttt)/(n3.count.aa$p*n3.count.aa$f))/n3.count.diaa$pf








cps.cga.aaa<-((n3.trinuc$cga*n3.trinuc$aaa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.cga.aac<-((n3.trinuc$cga*n3.trinuc$aac)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn
cps.cga.aag<-((n3.trinuc$cga*n3.trinuc$aag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.cga.aat<-((n3.trinuc$cga*n3.trinuc$aat)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn

cps.cga.aca<-((n3.trinuc$cga*n3.trinuc$aca)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cga.acc<-((n3.trinuc$cga*n3.trinuc$acc)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cga.acg<-((n3.trinuc$cga*n3.trinuc$acg)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cga.act<-((n3.trinuc$cga*n3.trinuc$act)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt

cps.cga.aga<-((n3.trinuc$cga*n3.trinuc$aga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cga.agc<-((n3.trinuc$cga*n3.trinuc$agc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cga.agg<-((n3.trinuc$cga*n3.trinuc$agg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cga.agt<-((n3.trinuc$cga*n3.trinuc$agt)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

cps.cga.ata<-((n3.trinuc$cga*n3.trinuc$ata)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.cga.atc<-((n3.trinuc$cga*n3.trinuc$atc)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.cga.atg<-((n3.trinuc$cga*n3.trinuc$atg)/(n3.count.aa$r*n3.count.aa$m))/n3.count.diaa$rm
cps.cga.att<-((n3.trinuc$cga*n3.trinuc$att)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri

cps.cga.caa<-((n3.trinuc$cga*n3.trinuc$caa)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.cga.cac<-((n3.trinuc$cga*n3.trinuc$cac)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh
cps.cga.cag<-((n3.trinuc$cga*n3.trinuc$cag)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.cga.cat<-((n3.trinuc$cga*n3.trinuc$cat)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh

cps.cga.cca<-((n3.trinuc$cga*n3.trinuc$cca)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cga.ccc<-((n3.trinuc$cga*n3.trinuc$ccc)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cga.ccg<-((n3.trinuc$cga*n3.trinuc$ccg)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cga.cct<-((n3.trinuc$cga*n3.trinuc$cct)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp

cps.cga.cga<-((n3.trinuc$cga*n3.trinuc$cga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cga.cgc<-((n3.trinuc$cga*n3.trinuc$cgc)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cga.cgg<-((n3.trinuc$cga*n3.trinuc$cgg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cga.cgt<-((n3.trinuc$cga*n3.trinuc$cgt)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr

cps.cga.cta<-((n3.trinuc$cga*n3.trinuc$cta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cga.ctc<-((n3.trinuc$cga*n3.trinuc$ctc)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cga.ctg<-((n3.trinuc$cga*n3.trinuc$ctg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cga.ctt<-((n3.trinuc$cga*n3.trinuc$ctt)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl

cps.cga.gaa<-((n3.trinuc$cga*n3.trinuc$gaa)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.cga.gac<-((n3.trinuc$cga*n3.trinuc$gac)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd
cps.cga.gag<-((n3.trinuc$cga*n3.trinuc$gag)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.cga.gat<-((n3.trinuc$cga*n3.trinuc$gat)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd

cps.cga.gca<-((n3.trinuc$cga*n3.trinuc$gca)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cga.gcc<-((n3.trinuc$cga*n3.trinuc$gcc)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cga.gcg<-((n3.trinuc$cga*n3.trinuc$gcg)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cga.gct<-((n3.trinuc$cga*n3.trinuc$gct)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra

cps.cga.gga<-((n3.trinuc$cga*n3.trinuc$gga)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cga.ggc<-((n3.trinuc$cga*n3.trinuc$ggc)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cga.ggg<-((n3.trinuc$cga*n3.trinuc$ggg)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cga.ggt<-((n3.trinuc$cga*n3.trinuc$ggt)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg

cps.cga.gta<-((n3.trinuc$cga*n3.trinuc$gta)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cga.gtc<-((n3.trinuc$cga*n3.trinuc$gtc)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cga.gtg<-((n3.trinuc$cga*n3.trinuc$gtg)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cga.gtt<-((n3.trinuc$cga*n3.trinuc$gtt)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv

#Stop codon
#cps.cga.taa<-((n3.trinuc$cga*n3.trinuc$taa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cga.tac<-((n3.trinuc$cga*n3.trinuc$tac)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry
#Stop codon
#cps.cga.tag<-((n3.trinuc$cga*n3.trinuc$tag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cga.tat<-((n3.trinuc$cga*n3.trinuc$tat)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry

cps.cga.tca<-((n3.trinuc$cga*n3.trinuc$tca)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cga.tcc<-((n3.trinuc$cga*n3.trinuc$tcc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cga.tcg<-((n3.trinuc$cga*n3.trinuc$tcg)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cga.tct<-((n3.trinuc$cga*n3.trinuc$tct)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

#Stop codon
#cps.cga.tga<-((n3.trinuc$cga*n3.trinuc$tga)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cga.tgc<-((n3.trinuc$cga*n3.trinuc$tgc)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc
cps.cga.tgg<-((n3.trinuc$cga*n3.trinuc$tgg)/(n3.count.aa$r*n3.count.aa$w))/n3.count.diaa$rw
cps.cga.tgt<-((n3.trinuc$cga*n3.trinuc$tgt)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc

cps.cga.tta<-((n3.trinuc$cga*n3.trinuc$tta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cga.ttc<-((n3.trinuc$cga*n3.trinuc$ttc)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf
cps.cga.ttg<-((n3.trinuc$cga*n3.trinuc$ttg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cga.ttt<-((n3.trinuc$cga*n3.trinuc$ttt)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf








cps.cgc.aaa<-((n3.trinuc$cgc*n3.trinuc$aaa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.cgc.aac<-((n3.trinuc$cgc*n3.trinuc$aac)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn
cps.cgc.aag<-((n3.trinuc$cgc*n3.trinuc$aag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.cgc.aat<-((n3.trinuc$cgc*n3.trinuc$aat)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn

cps.cgc.aca<-((n3.trinuc$cgc*n3.trinuc$aca)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgc.acc<-((n3.trinuc$cgc*n3.trinuc$acc)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgc.acg<-((n3.trinuc$cgc*n3.trinuc$acg)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgc.act<-((n3.trinuc$cgc*n3.trinuc$act)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt

cps.cgc.aga<-((n3.trinuc$cgc*n3.trinuc$aga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgc.agc<-((n3.trinuc$cgc*n3.trinuc$agc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgc.agg<-((n3.trinuc$cgc*n3.trinuc$agg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgc.agt<-((n3.trinuc$cgc*n3.trinuc$agt)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

cps.cgc.ata<-((n3.trinuc$cgc*n3.trinuc$ata)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.cgc.atc<-((n3.trinuc$cgc*n3.trinuc$atc)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.cgc.atg<-((n3.trinuc$cgc*n3.trinuc$atg)/(n3.count.aa$r*n3.count.aa$m))/n3.count.diaa$rm
cps.cgc.att<-((n3.trinuc$cgc*n3.trinuc$att)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri

cps.cgc.caa<-((n3.trinuc$cgc*n3.trinuc$caa)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.cgc.cac<-((n3.trinuc$cgc*n3.trinuc$cac)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh
cps.cgc.cag<-((n3.trinuc$cgc*n3.trinuc$cag)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.cgc.cat<-((n3.trinuc$cgc*n3.trinuc$cat)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh

cps.cgc.cca<-((n3.trinuc$cgc*n3.trinuc$cca)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgc.ccc<-((n3.trinuc$cgc*n3.trinuc$ccc)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgc.ccg<-((n3.trinuc$cgc*n3.trinuc$ccg)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgc.cct<-((n3.trinuc$cgc*n3.trinuc$cct)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp

cps.cgc.cga<-((n3.trinuc$cgc*n3.trinuc$cga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgc.cgc<-((n3.trinuc$cgc*n3.trinuc$cgc)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgc.cgg<-((n3.trinuc$cgc*n3.trinuc$cgg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgc.cgt<-((n3.trinuc$cgc*n3.trinuc$cgt)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr

cps.cgc.cta<-((n3.trinuc$cgc*n3.trinuc$cta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgc.ctc<-((n3.trinuc$cgc*n3.trinuc$ctc)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgc.ctg<-((n3.trinuc$cgc*n3.trinuc$ctg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgc.ctt<-((n3.trinuc$cgc*n3.trinuc$ctt)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl

cps.cgc.gaa<-((n3.trinuc$cgc*n3.trinuc$gaa)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.cgc.gac<-((n3.trinuc$cgc*n3.trinuc$gac)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd
cps.cgc.gag<-((n3.trinuc$cgc*n3.trinuc$gag)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.cgc.gat<-((n3.trinuc$cgc*n3.trinuc$gat)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd

cps.cgc.gca<-((n3.trinuc$cgc*n3.trinuc$gca)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgc.gcc<-((n3.trinuc$cgc*n3.trinuc$gcc)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgc.gcg<-((n3.trinuc$cgc*n3.trinuc$gcg)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgc.gct<-((n3.trinuc$cgc*n3.trinuc$gct)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra

cps.cgc.gga<-((n3.trinuc$cgc*n3.trinuc$gga)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgc.ggc<-((n3.trinuc$cgc*n3.trinuc$ggc)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgc.ggg<-((n3.trinuc$cgc*n3.trinuc$ggg)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgc.ggt<-((n3.trinuc$cgc*n3.trinuc$ggt)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg

cps.cgc.gta<-((n3.trinuc$cgc*n3.trinuc$gta)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgc.gtc<-((n3.trinuc$cgc*n3.trinuc$gtc)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgc.gtg<-((n3.trinuc$cgc*n3.trinuc$gtg)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgc.gtt<-((n3.trinuc$cgc*n3.trinuc$gtt)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv

#Stop codon
#cps.cgc.taa<-((n3.trinuc$cgc*n3.trinuc$taa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgc.tac<-((n3.trinuc$cgc*n3.trinuc$tac)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry
#Stop codon
#cps.cgc.tag<-((n3.trinuc$cgc*n3.trinuc$tag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgc.tat<-((n3.trinuc$cgc*n3.trinuc$tat)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry

cps.cgc.tca<-((n3.trinuc$cgc*n3.trinuc$tca)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgc.tcc<-((n3.trinuc$cgc*n3.trinuc$tcc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgc.tcg<-((n3.trinuc$cgc*n3.trinuc$tcg)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgc.tct<-((n3.trinuc$cgc*n3.trinuc$tct)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

#Stop codon
#cps.cgc.tga<-((n3.trinuc$cgc*n3.trinuc$tga)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgc.tgc<-((n3.trinuc$cgc*n3.trinuc$tgc)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc
cps.cgc.tgg<-((n3.trinuc$cgc*n3.trinuc$tgg)/(n3.count.aa$r*n3.count.aa$w))/n3.count.diaa$rw
cps.cgc.tgt<-((n3.trinuc$cgc*n3.trinuc$tgt)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc

cps.cgc.tta<-((n3.trinuc$cgc*n3.trinuc$tta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgc.ttc<-((n3.trinuc$cgc*n3.trinuc$ttc)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf
cps.cgc.ttg<-((n3.trinuc$cgc*n3.trinuc$ttg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgc.ttt<-((n3.trinuc$cgc*n3.trinuc$ttt)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf








cps.cgg.aaa<-((n3.trinuc$cgg*n3.trinuc$aaa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.cgg.aac<-((n3.trinuc$cgg*n3.trinuc$aac)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn
cps.cgg.aag<-((n3.trinuc$cgg*n3.trinuc$aag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.cgg.aat<-((n3.trinuc$cgg*n3.trinuc$aat)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn

cps.cgg.aca<-((n3.trinuc$cgg*n3.trinuc$aca)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgg.acc<-((n3.trinuc$cgg*n3.trinuc$acc)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgg.acg<-((n3.trinuc$cgg*n3.trinuc$acg)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgg.act<-((n3.trinuc$cgg*n3.trinuc$act)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt

cps.cgg.aga<-((n3.trinuc$cgg*n3.trinuc$aga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgg.agc<-((n3.trinuc$cgg*n3.trinuc$agc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgg.agg<-((n3.trinuc$cgg*n3.trinuc$agg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgg.agt<-((n3.trinuc$cgg*n3.trinuc$agt)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

cps.cgg.ata<-((n3.trinuc$cgg*n3.trinuc$ata)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.cgg.atc<-((n3.trinuc$cgg*n3.trinuc$atc)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.cgg.atg<-((n3.trinuc$cgg*n3.trinuc$atg)/(n3.count.aa$r*n3.count.aa$m))/n3.count.diaa$rm
cps.cgg.att<-((n3.trinuc$cgg*n3.trinuc$att)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri

cps.cgg.caa<-((n3.trinuc$cgg*n3.trinuc$caa)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.cgg.cac<-((n3.trinuc$cgg*n3.trinuc$cac)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh
cps.cgg.cag<-((n3.trinuc$cgg*n3.trinuc$cag)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.cgg.cat<-((n3.trinuc$cgg*n3.trinuc$cat)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh

cps.cgg.cca<-((n3.trinuc$cgg*n3.trinuc$cca)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgg.ccc<-((n3.trinuc$cgg*n3.trinuc$ccc)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgg.ccg<-((n3.trinuc$cgg*n3.trinuc$ccg)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgg.cct<-((n3.trinuc$cgg*n3.trinuc$cct)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp

cps.cgg.cga<-((n3.trinuc$cgg*n3.trinuc$cga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgg.cgc<-((n3.trinuc$cgg*n3.trinuc$cgc)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgg.cgg<-((n3.trinuc$cgg*n3.trinuc$cgg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgg.cgt<-((n3.trinuc$cgg*n3.trinuc$cgt)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr

cps.cgg.cta<-((n3.trinuc$cgg*n3.trinuc$cta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgg.ctc<-((n3.trinuc$cgg*n3.trinuc$ctc)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgg.ctg<-((n3.trinuc$cgg*n3.trinuc$ctg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgg.ctt<-((n3.trinuc$cgg*n3.trinuc$ctt)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl

cps.cgg.gaa<-((n3.trinuc$cgg*n3.trinuc$gaa)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.cgg.gac<-((n3.trinuc$cgg*n3.trinuc$gac)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd
cps.cgg.gag<-((n3.trinuc$cgg*n3.trinuc$gag)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.cgg.gat<-((n3.trinuc$cgg*n3.trinuc$gat)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd

cps.cgg.gca<-((n3.trinuc$cgg*n3.trinuc$gca)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgg.gcc<-((n3.trinuc$cgg*n3.trinuc$gcc)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgg.gcg<-((n3.trinuc$cgg*n3.trinuc$gcg)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgg.gct<-((n3.trinuc$cgg*n3.trinuc$gct)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra

cps.cgg.gga<-((n3.trinuc$cgg*n3.trinuc$gga)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgg.ggc<-((n3.trinuc$cgg*n3.trinuc$ggc)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgg.ggg<-((n3.trinuc$cgg*n3.trinuc$ggg)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgg.ggt<-((n3.trinuc$cgg*n3.trinuc$ggt)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg

cps.cgg.gta<-((n3.trinuc$cgg*n3.trinuc$gta)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgg.gtc<-((n3.trinuc$cgg*n3.trinuc$gtc)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgg.gtg<-((n3.trinuc$cgg*n3.trinuc$gtg)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgg.gtt<-((n3.trinuc$cgg*n3.trinuc$gtt)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv

#Stop codon
#cps.cgg.taa<-((n3.trinuc$cgg*n3.trinuc$taa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgg.tac<-((n3.trinuc$cgg*n3.trinuc$tac)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry
#Stop codon
#cps.cgg.tag<-((n3.trinuc$cgg*n3.trinuc$tag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgg.tat<-((n3.trinuc$cgg*n3.trinuc$tat)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry

cps.cgg.tca<-((n3.trinuc$cgg*n3.trinuc$tca)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgg.tcc<-((n3.trinuc$cgg*n3.trinuc$tcc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgg.tcg<-((n3.trinuc$cgg*n3.trinuc$tcg)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgg.tct<-((n3.trinuc$cgg*n3.trinuc$tct)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

#Stop codon
#cps.cgg.tga<-((n3.trinuc$cgg*n3.trinuc$tga)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgg.tgc<-((n3.trinuc$cgg*n3.trinuc$tgc)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc
cps.cgg.tgg<-((n3.trinuc$cgg*n3.trinuc$tgg)/(n3.count.aa$r*n3.count.aa$w))/n3.count.diaa$rw
cps.cgg.tgt<-((n3.trinuc$cgg*n3.trinuc$tgt)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc

cps.cgg.tta<-((n3.trinuc$cgg*n3.trinuc$tta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgg.ttc<-((n3.trinuc$cgg*n3.trinuc$ttc)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf
cps.cgg.ttg<-((n3.trinuc$cgg*n3.trinuc$ttg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgg.ttt<-((n3.trinuc$cgg*n3.trinuc$ttt)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf








cps.cgt.aaa<-((n3.trinuc$cgt*n3.trinuc$aaa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.cgt.aac<-((n3.trinuc$cgt*n3.trinuc$aac)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn
cps.cgt.aag<-((n3.trinuc$cgt*n3.trinuc$aag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$rk
cps.cgt.aat<-((n3.trinuc$cgt*n3.trinuc$aat)/(n3.count.aa$r*n3.count.aa$n))/n3.count.diaa$rn

cps.cgt.aca<-((n3.trinuc$cgt*n3.trinuc$aca)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgt.acc<-((n3.trinuc$cgt*n3.trinuc$acc)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgt.acg<-((n3.trinuc$cgt*n3.trinuc$acg)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt
cps.cgt.act<-((n3.trinuc$cgt*n3.trinuc$act)/(n3.count.aa$r*n3.count.aa$t))/n3.count.diaa$rt

cps.cgt.aga<-((n3.trinuc$cgt*n3.trinuc$aga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgt.agc<-((n3.trinuc$cgt*n3.trinuc$agc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgt.agg<-((n3.trinuc$cgt*n3.trinuc$agg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgt.agt<-((n3.trinuc$cgt*n3.trinuc$agt)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

cps.cgt.ata<-((n3.trinuc$cgt*n3.trinuc$ata)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.cgt.atc<-((n3.trinuc$cgt*n3.trinuc$atc)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri
cps.cgt.atg<-((n3.trinuc$cgt*n3.trinuc$atg)/(n3.count.aa$r*n3.count.aa$m))/n3.count.diaa$rm
cps.cgt.att<-((n3.trinuc$cgt*n3.trinuc$att)/(n3.count.aa$r*n3.count.aa$i))/n3.count.diaa$ri

cps.cgt.caa<-((n3.trinuc$cgt*n3.trinuc$caa)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.cgt.cac<-((n3.trinuc$cgt*n3.trinuc$cac)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh
cps.cgt.cag<-((n3.trinuc$cgt*n3.trinuc$cag)/(n3.count.aa$r*n3.count.aa$q))/n3.count.diaa$rq
cps.cgt.cat<-((n3.trinuc$cgt*n3.trinuc$cat)/(n3.count.aa$r*n3.count.aa$h))/n3.count.diaa$rh

cps.cgt.cca<-((n3.trinuc$cgt*n3.trinuc$cca)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgt.ccc<-((n3.trinuc$cgt*n3.trinuc$ccc)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgt.ccg<-((n3.trinuc$cgt*n3.trinuc$ccg)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp
cps.cgt.cct<-((n3.trinuc$cgt*n3.trinuc$cct)/(n3.count.aa$r*n3.count.aa$p))/n3.count.diaa$rp

cps.cgt.cga<-((n3.trinuc$cgt*n3.trinuc$cga)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgt.cgc<-((n3.trinuc$cgt*n3.trinuc$cgc)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgt.cgg<-((n3.trinuc$cgt*n3.trinuc$cgg)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr
cps.cgt.cgt<-((n3.trinuc$cgt*n3.trinuc$cgt)/(n3.count.aa$r*n3.count.aa$r))/n3.count.diaa$rr

cps.cgt.cta<-((n3.trinuc$cgt*n3.trinuc$cta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgt.ctc<-((n3.trinuc$cgt*n3.trinuc$ctc)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgt.ctg<-((n3.trinuc$cgt*n3.trinuc$ctg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgt.ctt<-((n3.trinuc$cgt*n3.trinuc$ctt)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl

cps.cgt.gaa<-((n3.trinuc$cgt*n3.trinuc$gaa)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.cgt.gac<-((n3.trinuc$cgt*n3.trinuc$gac)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd
cps.cgt.gag<-((n3.trinuc$cgt*n3.trinuc$gag)/(n3.count.aa$r*n3.count.aa$e))/n3.count.diaa$re
cps.cgt.gat<-((n3.trinuc$cgt*n3.trinuc$gat)/(n3.count.aa$r*n3.count.aa$d))/n3.count.diaa$rd

cps.cgt.gca<-((n3.trinuc$cgt*n3.trinuc$gca)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgt.gcc<-((n3.trinuc$cgt*n3.trinuc$gcc)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgt.gcg<-((n3.trinuc$cgt*n3.trinuc$gcg)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra
cps.cgt.gct<-((n3.trinuc$cgt*n3.trinuc$gct)/(n3.count.aa$r*n3.count.aa$a))/n3.count.diaa$ra

cps.cgt.gga<-((n3.trinuc$cgt*n3.trinuc$gga)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgt.ggc<-((n3.trinuc$cgt*n3.trinuc$ggc)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgt.ggg<-((n3.trinuc$cgt*n3.trinuc$ggg)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg
cps.cgt.ggt<-((n3.trinuc$cgt*n3.trinuc$ggt)/(n3.count.aa$r*n3.count.aa$g))/n3.count.diaa$rg

cps.cgt.gta<-((n3.trinuc$cgt*n3.trinuc$gta)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgt.gtc<-((n3.trinuc$cgt*n3.trinuc$gtc)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgt.gtg<-((n3.trinuc$cgt*n3.trinuc$gtg)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv
cps.cgt.gtt<-((n3.trinuc$cgt*n3.trinuc$gtt)/(n3.count.aa$r*n3.count.aa$v))/n3.count.diaa$rv

#Stop codon
#cps.cgt.taa<-((n3.trinuc$cgt*n3.trinuc$taa)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgt.tac<-((n3.trinuc$cgt*n3.trinuc$tac)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry
#Stop codon
#cps.cgt.tag<-((n3.trinuc$cgt*n3.trinuc$tag)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgt.tat<-((n3.trinuc$cgt*n3.trinuc$tat)/(n3.count.aa$r*n3.count.aa$y))/n3.count.diaa$ry

cps.cgt.tca<-((n3.trinuc$cgt*n3.trinuc$tca)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgt.tcc<-((n3.trinuc$cgt*n3.trinuc$tcc)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgt.tcg<-((n3.trinuc$cgt*n3.trinuc$tcg)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs
cps.cgt.tct<-((n3.trinuc$cgt*n3.trinuc$tct)/(n3.count.aa$r*n3.count.aa$s))/n3.count.diaa$rs

#Stop codon
#cps.cgt.tga<-((n3.trinuc$cgt*n3.trinuc$tga)/(n3.count.aa$r*n3.count.aa$k))/n3.count.diaa$kk
cps.cgt.tgc<-((n3.trinuc$cgt*n3.trinuc$tgc)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc
cps.cgt.tgg<-((n3.trinuc$cgt*n3.trinuc$tgg)/(n3.count.aa$r*n3.count.aa$w))/n3.count.diaa$rw
cps.cgt.tgt<-((n3.trinuc$cgt*n3.trinuc$tgt)/(n3.count.aa$r*n3.count.aa$c))/n3.count.diaa$rc

cps.cgt.tta<-((n3.trinuc$cgt*n3.trinuc$tta)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgt.ttc<-((n3.trinuc$cgt*n3.trinuc$ttc)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf
cps.cgt.ttg<-((n3.trinuc$cgt*n3.trinuc$ttg)/(n3.count.aa$r*n3.count.aa$l))/n3.count.diaa$rl
cps.cgt.ttt<-((n3.trinuc$cgt*n3.trinuc$ttt)/(n3.count.aa$r*n3.count.aa$f))/n3.count.diaa$rf








cps.cta.aaa<-((n3.trinuc$cta*n3.trinuc$aaa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.cta.aac<-((n3.trinuc$cta*n3.trinuc$aac)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln
cps.cta.aag<-((n3.trinuc$cta*n3.trinuc$aag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.cta.aat<-((n3.trinuc$cta*n3.trinuc$aat)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln

cps.cta.aca<-((n3.trinuc$cta*n3.trinuc$aca)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.cta.acc<-((n3.trinuc$cta*n3.trinuc$acc)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.cta.acg<-((n3.trinuc$cta*n3.trinuc$acg)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.cta.act<-((n3.trinuc$cta*n3.trinuc$act)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt

cps.cta.aga<-((n3.trinuc$cta*n3.trinuc$aga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.cta.agc<-((n3.trinuc$cta*n3.trinuc$agc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.cta.agg<-((n3.trinuc$cta*n3.trinuc$agg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.cta.agt<-((n3.trinuc$cta*n3.trinuc$agt)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

cps.cta.ata<-((n3.trinuc$cta*n3.trinuc$ata)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.cta.atc<-((n3.trinuc$cta*n3.trinuc$atc)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.cta.atg<-((n3.trinuc$cta*n3.trinuc$atg)/(n3.count.aa$l*n3.count.aa$m))/n3.count.diaa$lm
cps.cta.att<-((n3.trinuc$cta*n3.trinuc$att)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li

cps.cta.caa<-((n3.trinuc$cta*n3.trinuc$caa)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.cta.cac<-((n3.trinuc$cta*n3.trinuc$cac)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh
cps.cta.cag<-((n3.trinuc$cta*n3.trinuc$cag)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.cta.cat<-((n3.trinuc$cta*n3.trinuc$cat)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh

cps.cta.cca<-((n3.trinuc$cta*n3.trinuc$cca)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.cta.ccc<-((n3.trinuc$cta*n3.trinuc$ccc)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.cta.ccg<-((n3.trinuc$cta*n3.trinuc$ccg)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.cta.cct<-((n3.trinuc$cta*n3.trinuc$cct)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp

cps.cta.cga<-((n3.trinuc$cta*n3.trinuc$cga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.cta.cgc<-((n3.trinuc$cta*n3.trinuc$cgc)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.cta.cgg<-((n3.trinuc$cta*n3.trinuc$cgg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.cta.cgt<-((n3.trinuc$cta*n3.trinuc$cgt)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr

cps.cta.cta<-((n3.trinuc$cta*n3.trinuc$cta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.cta.ctc<-((n3.trinuc$cta*n3.trinuc$ctc)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.cta.ctg<-((n3.trinuc$cta*n3.trinuc$ctg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.cta.ctt<-((n3.trinuc$cta*n3.trinuc$ctt)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll

cps.cta.gaa<-((n3.trinuc$cta*n3.trinuc$gaa)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.cta.gac<-((n3.trinuc$cta*n3.trinuc$gac)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld
cps.cta.gag<-((n3.trinuc$cta*n3.trinuc$gag)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.cta.gat<-((n3.trinuc$cta*n3.trinuc$gat)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld

cps.cta.gca<-((n3.trinuc$cta*n3.trinuc$gca)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.cta.gcc<-((n3.trinuc$cta*n3.trinuc$gcc)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.cta.gcg<-((n3.trinuc$cta*n3.trinuc$gcg)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.cta.gct<-((n3.trinuc$cta*n3.trinuc$gct)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la

cps.cta.gga<-((n3.trinuc$cta*n3.trinuc$gga)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.cta.ggc<-((n3.trinuc$cta*n3.trinuc$ggc)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.cta.ggg<-((n3.trinuc$cta*n3.trinuc$ggg)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.cta.ggt<-((n3.trinuc$cta*n3.trinuc$ggt)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg

cps.cta.gta<-((n3.trinuc$cta*n3.trinuc$gta)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.cta.gtc<-((n3.trinuc$cta*n3.trinuc$gtc)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.cta.gtg<-((n3.trinuc$cta*n3.trinuc$gtg)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.cta.gtt<-((n3.trinuc$cta*n3.trinuc$gtt)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv

#Stop codon
#cps.cta.taa<-((n3.trinuc$cta*n3.trinuc$taa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.cta.tac<-((n3.trinuc$cta*n3.trinuc$tac)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly
#Stop codon
#cps.cta.tag<-((n3.trinuc$cta*n3.trinuc$tag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.cta.tat<-((n3.trinuc$cta*n3.trinuc$tat)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly

cps.cta.tca<-((n3.trinuc$cta*n3.trinuc$tca)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.cta.tcc<-((n3.trinuc$cta*n3.trinuc$tcc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.cta.tcg<-((n3.trinuc$cta*n3.trinuc$tcg)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.cta.tct<-((n3.trinuc$cta*n3.trinuc$tct)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

#Stop codon
#cps.cta.tga<-((n3.trinuc$cta*n3.trinuc$tga)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.cta.tgc<-((n3.trinuc$cta*n3.trinuc$tgc)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc
cps.cta.tgg<-((n3.trinuc$cta*n3.trinuc$tgg)/(n3.count.aa$l*n3.count.aa$w))/n3.count.diaa$lw
cps.cta.tgt<-((n3.trinuc$cta*n3.trinuc$tgt)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc

cps.cta.tta<-((n3.trinuc$cta*n3.trinuc$tta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.cta.ttc<-((n3.trinuc$cta*n3.trinuc$ttc)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf
cps.cta.ttg<-((n3.trinuc$cta*n3.trinuc$ttg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.cta.ttt<-((n3.trinuc$cta*n3.trinuc$ttt)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf








cps.ctc.aaa<-((n3.trinuc$ctc*n3.trinuc$aaa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctc.aac<-((n3.trinuc$ctc*n3.trinuc$aac)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln
cps.ctc.aag<-((n3.trinuc$ctc*n3.trinuc$aag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctc.aat<-((n3.trinuc$ctc*n3.trinuc$aat)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln

cps.ctc.aca<-((n3.trinuc$ctc*n3.trinuc$aca)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctc.acc<-((n3.trinuc$ctc*n3.trinuc$acc)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctc.acg<-((n3.trinuc$ctc*n3.trinuc$acg)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctc.act<-((n3.trinuc$ctc*n3.trinuc$act)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt

cps.ctc.aga<-((n3.trinuc$ctc*n3.trinuc$aga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctc.agc<-((n3.trinuc$ctc*n3.trinuc$agc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctc.agg<-((n3.trinuc$ctc*n3.trinuc$agg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctc.agt<-((n3.trinuc$ctc*n3.trinuc$agt)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

cps.ctc.ata<-((n3.trinuc$ctc*n3.trinuc$ata)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.ctc.atc<-((n3.trinuc$ctc*n3.trinuc$atc)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.ctc.atg<-((n3.trinuc$ctc*n3.trinuc$atg)/(n3.count.aa$l*n3.count.aa$m))/n3.count.diaa$lm
cps.ctc.att<-((n3.trinuc$ctc*n3.trinuc$att)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li

cps.ctc.caa<-((n3.trinuc$ctc*n3.trinuc$caa)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.ctc.cac<-((n3.trinuc$ctc*n3.trinuc$cac)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh
cps.ctc.cag<-((n3.trinuc$ctc*n3.trinuc$cag)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.ctc.cat<-((n3.trinuc$ctc*n3.trinuc$cat)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh

cps.ctc.cca<-((n3.trinuc$ctc*n3.trinuc$cca)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctc.ccc<-((n3.trinuc$ctc*n3.trinuc$ccc)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctc.ccg<-((n3.trinuc$ctc*n3.trinuc$ccg)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctc.cct<-((n3.trinuc$ctc*n3.trinuc$cct)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp

cps.ctc.cga<-((n3.trinuc$ctc*n3.trinuc$cga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctc.cgc<-((n3.trinuc$ctc*n3.trinuc$cgc)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctc.cgg<-((n3.trinuc$ctc*n3.trinuc$cgg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctc.cgt<-((n3.trinuc$ctc*n3.trinuc$cgt)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr

cps.ctc.cta<-((n3.trinuc$ctc*n3.trinuc$cta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctc.ctc<-((n3.trinuc$ctc*n3.trinuc$ctc)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctc.ctg<-((n3.trinuc$ctc*n3.trinuc$ctg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctc.ctt<-((n3.trinuc$ctc*n3.trinuc$ctt)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll

cps.ctc.gaa<-((n3.trinuc$ctc*n3.trinuc$gaa)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.ctc.gac<-((n3.trinuc$ctc*n3.trinuc$gac)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld
cps.ctc.gag<-((n3.trinuc$ctc*n3.trinuc$gag)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.ctc.gat<-((n3.trinuc$ctc*n3.trinuc$gat)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld

cps.ctc.gca<-((n3.trinuc$ctc*n3.trinuc$gca)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctc.gcc<-((n3.trinuc$ctc*n3.trinuc$gcc)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctc.gcg<-((n3.trinuc$ctc*n3.trinuc$gcg)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctc.gct<-((n3.trinuc$ctc*n3.trinuc$gct)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la

cps.ctc.gga<-((n3.trinuc$ctc*n3.trinuc$gga)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctc.ggc<-((n3.trinuc$ctc*n3.trinuc$ggc)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctc.ggg<-((n3.trinuc$ctc*n3.trinuc$ggg)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctc.ggt<-((n3.trinuc$ctc*n3.trinuc$ggt)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg

cps.ctc.gta<-((n3.trinuc$ctc*n3.trinuc$gta)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctc.gtc<-((n3.trinuc$ctc*n3.trinuc$gtc)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctc.gtg<-((n3.trinuc$ctc*n3.trinuc$gtg)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctc.gtt<-((n3.trinuc$ctc*n3.trinuc$gtt)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv

#Stop codon
#cps.ctc.taa<-((n3.trinuc$ctc*n3.trinuc$taa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctc.tac<-((n3.trinuc$ctc*n3.trinuc$tac)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly
#Stop codon
#cps.ctc.tag<-((n3.trinuc$ctc*n3.trinuc$tag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctc.tat<-((n3.trinuc$ctc*n3.trinuc$tat)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly

cps.ctc.tca<-((n3.trinuc$ctc*n3.trinuc$tca)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctc.tcc<-((n3.trinuc$ctc*n3.trinuc$tcc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctc.tcg<-((n3.trinuc$ctc*n3.trinuc$tcg)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctc.tct<-((n3.trinuc$ctc*n3.trinuc$tct)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

#Stop codon
#cps.ctc.tga<-((n3.trinuc$ctc*n3.trinuc$tga)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctc.tgc<-((n3.trinuc$ctc*n3.trinuc$tgc)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc
cps.ctc.tgg<-((n3.trinuc$ctc*n3.trinuc$tgg)/(n3.count.aa$l*n3.count.aa$w))/n3.count.diaa$lw
cps.ctc.tgt<-((n3.trinuc$ctc*n3.trinuc$tgt)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc

cps.ctc.tta<-((n3.trinuc$ctc*n3.trinuc$tta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctc.ttc<-((n3.trinuc$ctc*n3.trinuc$ttc)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf
cps.ctc.ttg<-((n3.trinuc$ctc*n3.trinuc$ttg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctc.ttt<-((n3.trinuc$ctc*n3.trinuc$ttt)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf








cps.ctg.aaa<-((n3.trinuc$ctg*n3.trinuc$aaa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctg.aac<-((n3.trinuc$ctg*n3.trinuc$aac)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln
cps.ctg.aag<-((n3.trinuc$ctg*n3.trinuc$aag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctg.aat<-((n3.trinuc$ctg*n3.trinuc$aat)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln

cps.ctg.aca<-((n3.trinuc$ctg*n3.trinuc$aca)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctg.acc<-((n3.trinuc$ctg*n3.trinuc$acc)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctg.acg<-((n3.trinuc$ctg*n3.trinuc$acg)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctg.act<-((n3.trinuc$ctg*n3.trinuc$act)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt

cps.ctg.aga<-((n3.trinuc$ctg*n3.trinuc$aga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctg.agc<-((n3.trinuc$ctg*n3.trinuc$agc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctg.agg<-((n3.trinuc$ctg*n3.trinuc$agg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctg.agt<-((n3.trinuc$ctg*n3.trinuc$agt)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

cps.ctg.ata<-((n3.trinuc$ctg*n3.trinuc$ata)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.ctg.atc<-((n3.trinuc$ctg*n3.trinuc$atc)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.ctg.atg<-((n3.trinuc$ctg*n3.trinuc$atg)/(n3.count.aa$l*n3.count.aa$m))/n3.count.diaa$lm
cps.ctg.att<-((n3.trinuc$ctg*n3.trinuc$att)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li

cps.ctg.caa<-((n3.trinuc$ctg*n3.trinuc$caa)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.ctg.cac<-((n3.trinuc$ctg*n3.trinuc$cac)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh
cps.ctg.cag<-((n3.trinuc$ctg*n3.trinuc$cag)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.ctg.cat<-((n3.trinuc$ctg*n3.trinuc$cat)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh

cps.ctg.cca<-((n3.trinuc$ctg*n3.trinuc$cca)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctg.ccc<-((n3.trinuc$ctg*n3.trinuc$ccc)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctg.ccg<-((n3.trinuc$ctg*n3.trinuc$ccg)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctg.cct<-((n3.trinuc$ctg*n3.trinuc$cct)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp

cps.ctg.cga<-((n3.trinuc$ctg*n3.trinuc$cga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctg.cgc<-((n3.trinuc$ctg*n3.trinuc$cgc)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctg.cgg<-((n3.trinuc$ctg*n3.trinuc$cgg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctg.cgt<-((n3.trinuc$ctg*n3.trinuc$cgt)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr

cps.ctg.cta<-((n3.trinuc$ctg*n3.trinuc$cta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctg.ctc<-((n3.trinuc$ctg*n3.trinuc$ctc)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctg.ctg<-((n3.trinuc$ctg*n3.trinuc$ctg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctg.ctt<-((n3.trinuc$ctg*n3.trinuc$ctt)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll

cps.ctg.gaa<-((n3.trinuc$ctg*n3.trinuc$gaa)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.ctg.gac<-((n3.trinuc$ctg*n3.trinuc$gac)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld
cps.ctg.gag<-((n3.trinuc$ctg*n3.trinuc$gag)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.ctg.gat<-((n3.trinuc$ctg*n3.trinuc$gat)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld

cps.ctg.gca<-((n3.trinuc$ctg*n3.trinuc$gca)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctg.gcc<-((n3.trinuc$ctg*n3.trinuc$gcc)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctg.gcg<-((n3.trinuc$ctg*n3.trinuc$gcg)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctg.gct<-((n3.trinuc$ctg*n3.trinuc$gct)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la

cps.ctg.gga<-((n3.trinuc$ctg*n3.trinuc$gga)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctg.ggc<-((n3.trinuc$ctg*n3.trinuc$ggc)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctg.ggg<-((n3.trinuc$ctg*n3.trinuc$ggg)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctg.ggt<-((n3.trinuc$ctg*n3.trinuc$ggt)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg

cps.ctg.gta<-((n3.trinuc$ctg*n3.trinuc$gta)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctg.gtc<-((n3.trinuc$ctg*n3.trinuc$gtc)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctg.gtg<-((n3.trinuc$ctg*n3.trinuc$gtg)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctg.gtt<-((n3.trinuc$ctg*n3.trinuc$gtt)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv

#Stop codon
#cps.ctg.taa<-((n3.trinuc$ctg*n3.trinuc$taa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctg.tac<-((n3.trinuc$ctg*n3.trinuc$tac)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly
#Stop codon
#cps.ctg.tag<-((n3.trinuc$ctg*n3.trinuc$tag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctg.tat<-((n3.trinuc$ctg*n3.trinuc$tat)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly

cps.ctg.tca<-((n3.trinuc$ctg*n3.trinuc$tca)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctg.tcc<-((n3.trinuc$ctg*n3.trinuc$tcc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctg.tcg<-((n3.trinuc$ctg*n3.trinuc$tcg)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctg.tct<-((n3.trinuc$ctg*n3.trinuc$tct)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

#Stop codon
#cps.ctg.tga<-((n3.trinuc$ctg*n3.trinuc$tga)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctg.tgc<-((n3.trinuc$ctg*n3.trinuc$tgc)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc
cps.ctg.tgg<-((n3.trinuc$ctg*n3.trinuc$tgg)/(n3.count.aa$l*n3.count.aa$w))/n3.count.diaa$lw
cps.ctg.tgt<-((n3.trinuc$ctg*n3.trinuc$tgt)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc

cps.ctg.tta<-((n3.trinuc$ctg*n3.trinuc$tta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctg.ttc<-((n3.trinuc$ctg*n3.trinuc$ttc)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf
cps.ctg.ttg<-((n3.trinuc$ctg*n3.trinuc$ttg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctg.ttt<-((n3.trinuc$ctg*n3.trinuc$ttt)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf








cps.ctt.aaa<-((n3.trinuc$ctt*n3.trinuc$aaa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctt.aac<-((n3.trinuc$ctt*n3.trinuc$aac)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln
cps.ctt.aag<-((n3.trinuc$ctt*n3.trinuc$aag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctt.aat<-((n3.trinuc$ctt*n3.trinuc$aat)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln

cps.ctt.aca<-((n3.trinuc$ctt*n3.trinuc$aca)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctt.acc<-((n3.trinuc$ctt*n3.trinuc$acc)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctt.acg<-((n3.trinuc$ctt*n3.trinuc$acg)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ctt.act<-((n3.trinuc$ctt*n3.trinuc$act)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt

cps.ctt.aga<-((n3.trinuc$ctt*n3.trinuc$aga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctt.agc<-((n3.trinuc$ctt*n3.trinuc$agc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctt.agg<-((n3.trinuc$ctt*n3.trinuc$agg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctt.agt<-((n3.trinuc$ctt*n3.trinuc$agt)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

cps.ctt.ata<-((n3.trinuc$ctt*n3.trinuc$ata)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.ctt.atc<-((n3.trinuc$ctt*n3.trinuc$atc)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.ctt.atg<-((n3.trinuc$ctt*n3.trinuc$atg)/(n3.count.aa$l*n3.count.aa$m))/n3.count.diaa$lm
cps.ctt.att<-((n3.trinuc$ctt*n3.trinuc$att)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li

cps.ctt.caa<-((n3.trinuc$ctt*n3.trinuc$caa)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.ctt.cac<-((n3.trinuc$ctt*n3.trinuc$cac)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh
cps.ctt.cag<-((n3.trinuc$ctt*n3.trinuc$cag)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.ctt.cat<-((n3.trinuc$ctt*n3.trinuc$cat)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh

cps.ctt.cca<-((n3.trinuc$ctt*n3.trinuc$cca)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctt.ccc<-((n3.trinuc$ctt*n3.trinuc$ccc)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctt.ccg<-((n3.trinuc$ctt*n3.trinuc$ccg)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ctt.cct<-((n3.trinuc$ctt*n3.trinuc$cct)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp

cps.ctt.cga<-((n3.trinuc$ctt*n3.trinuc$cga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctt.cgc<-((n3.trinuc$ctt*n3.trinuc$cgc)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctt.cgg<-((n3.trinuc$ctt*n3.trinuc$cgg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ctt.cgt<-((n3.trinuc$ctt*n3.trinuc$cgt)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr

cps.ctt.cta<-((n3.trinuc$ctt*n3.trinuc$cta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctt.ctc<-((n3.trinuc$ctt*n3.trinuc$ctc)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctt.ctg<-((n3.trinuc$ctt*n3.trinuc$ctg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctt.ctt<-((n3.trinuc$ctt*n3.trinuc$ctt)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll

cps.ctt.gaa<-((n3.trinuc$ctt*n3.trinuc$gaa)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.ctt.gac<-((n3.trinuc$ctt*n3.trinuc$gac)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld
cps.ctt.gag<-((n3.trinuc$ctt*n3.trinuc$gag)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.ctt.gat<-((n3.trinuc$ctt*n3.trinuc$gat)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld

cps.ctt.gca<-((n3.trinuc$ctt*n3.trinuc$gca)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctt.gcc<-((n3.trinuc$ctt*n3.trinuc$gcc)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctt.gcg<-((n3.trinuc$ctt*n3.trinuc$gcg)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ctt.gct<-((n3.trinuc$ctt*n3.trinuc$gct)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la

cps.ctt.gga<-((n3.trinuc$ctt*n3.trinuc$gga)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctt.ggc<-((n3.trinuc$ctt*n3.trinuc$ggc)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctt.ggg<-((n3.trinuc$ctt*n3.trinuc$ggg)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ctt.ggt<-((n3.trinuc$ctt*n3.trinuc$ggt)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg

cps.ctt.gta<-((n3.trinuc$ctt*n3.trinuc$gta)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctt.gtc<-((n3.trinuc$ctt*n3.trinuc$gtc)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctt.gtg<-((n3.trinuc$ctt*n3.trinuc$gtg)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ctt.gtt<-((n3.trinuc$ctt*n3.trinuc$gtt)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv

#Stop codon
#cps.ctt.taa<-((n3.trinuc$ctt*n3.trinuc$taa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctt.tac<-((n3.trinuc$ctt*n3.trinuc$tac)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly
#Stop codon
#cps.ctt.tag<-((n3.trinuc$ctt*n3.trinuc$tag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctt.tat<-((n3.trinuc$ctt*n3.trinuc$tat)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly

cps.ctt.tca<-((n3.trinuc$ctt*n3.trinuc$tca)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctt.tcc<-((n3.trinuc$ctt*n3.trinuc$tcc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctt.tcg<-((n3.trinuc$ctt*n3.trinuc$tcg)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ctt.tct<-((n3.trinuc$ctt*n3.trinuc$tct)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

#Stop codon
#cps.ctt.tga<-((n3.trinuc$ctt*n3.trinuc$tga)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ctt.tgc<-((n3.trinuc$ctt*n3.trinuc$tgc)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc
cps.ctt.tgg<-((n3.trinuc$ctt*n3.trinuc$tgg)/(n3.count.aa$l*n3.count.aa$w))/n3.count.diaa$lw
cps.ctt.tgt<-((n3.trinuc$ctt*n3.trinuc$tgt)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc

cps.ctt.tta<-((n3.trinuc$ctt*n3.trinuc$tta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctt.ttc<-((n3.trinuc$ctt*n3.trinuc$ttc)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf
cps.ctt.ttg<-((n3.trinuc$ctt*n3.trinuc$ttg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ctt.ttt<-((n3.trinuc$ctt*n3.trinuc$ttt)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf








cps.gaa.aaa<-((n3.trinuc$gaa*n3.trinuc$aaa)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gaa.aac<-((n3.trinuc$gaa*n3.trinuc$aac)/(n3.count.aa$e*n3.count.aa$n))/n3.count.diaa$en
cps.gaa.aag<-((n3.trinuc$gaa*n3.trinuc$aag)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gaa.aat<-((n3.trinuc$gaa*n3.trinuc$aat)/(n3.count.aa$e*n3.count.aa$n))/n3.count.diaa$en

cps.gaa.aca<-((n3.trinuc$gaa*n3.trinuc$aca)/(n3.count.aa$e*n3.count.aa$t))/n3.count.diaa$et
cps.gaa.acc<-((n3.trinuc$gaa*n3.trinuc$acc)/(n3.count.aa$e*n3.count.aa$t))/n3.count.diaa$et
cps.gaa.acg<-((n3.trinuc$gaa*n3.trinuc$acg)/(n3.count.aa$e*n3.count.aa$t))/n3.count.diaa$et
cps.gaa.act<-((n3.trinuc$gaa*n3.trinuc$act)/(n3.count.aa$e*n3.count.aa$t))/n3.count.diaa$et

cps.gaa.aga<-((n3.trinuc$gaa*n3.trinuc$aga)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gaa.agc<-((n3.trinuc$gaa*n3.trinuc$agc)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es
cps.gaa.agg<-((n3.trinuc$gaa*n3.trinuc$agg)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gaa.agt<-((n3.trinuc$gaa*n3.trinuc$agt)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es

cps.gaa.ata<-((n3.trinuc$gaa*n3.trinuc$ata)/(n3.count.aa$e*n3.count.aa$i))/n3.count.diaa$ei
cps.gaa.atc<-((n3.trinuc$gaa*n3.trinuc$atc)/(n3.count.aa$e*n3.count.aa$i))/n3.count.diaa$ei
cps.gaa.atg<-((n3.trinuc$gaa*n3.trinuc$atg)/(n3.count.aa$e*n3.count.aa$m))/n3.count.diaa$em
cps.gaa.att<-((n3.trinuc$gaa*n3.trinuc$att)/(n3.count.aa$e*n3.count.aa$i))/n3.count.diaa$ei

cps.gaa.caa<-((n3.trinuc$gaa*n3.trinuc$caa)/(n3.count.aa$e*n3.count.aa$q))/n3.count.diaa$eq
cps.gaa.cac<-((n3.trinuc$gaa*n3.trinuc$cac)/(n3.count.aa$e*n3.count.aa$h))/n3.count.diaa$eh
cps.gaa.cag<-((n3.trinuc$gaa*n3.trinuc$cag)/(n3.count.aa$e*n3.count.aa$q))/n3.count.diaa$eq
cps.gaa.cat<-((n3.trinuc$gaa*n3.trinuc$cat)/(n3.count.aa$e*n3.count.aa$h))/n3.count.diaa$eh

cps.gaa.cca<-((n3.trinuc$gaa*n3.trinuc$cca)/(n3.count.aa$e*n3.count.aa$p))/n3.count.diaa$ep
cps.gaa.ccc<-((n3.trinuc$gaa*n3.trinuc$ccc)/(n3.count.aa$e*n3.count.aa$p))/n3.count.diaa$ep
cps.gaa.ccg<-((n3.trinuc$gaa*n3.trinuc$ccg)/(n3.count.aa$e*n3.count.aa$p))/n3.count.diaa$ep
cps.gaa.cct<-((n3.trinuc$gaa*n3.trinuc$cct)/(n3.count.aa$e*n3.count.aa$p))/n3.count.diaa$ep

cps.gaa.cga<-((n3.trinuc$gaa*n3.trinuc$cga)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gaa.cgc<-((n3.trinuc$gaa*n3.trinuc$cgc)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gaa.cgg<-((n3.trinuc$gaa*n3.trinuc$cgg)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gaa.cgt<-((n3.trinuc$gaa*n3.trinuc$cgt)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er

cps.gaa.cta<-((n3.trinuc$gaa*n3.trinuc$cta)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gaa.ctc<-((n3.trinuc$gaa*n3.trinuc$ctc)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gaa.ctg<-((n3.trinuc$gaa*n3.trinuc$ctg)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gaa.ctt<-((n3.trinuc$gaa*n3.trinuc$ctt)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el

cps.gaa.gaa<-((n3.trinuc$gaa*n3.trinuc$gaa)/(n3.count.aa$e*n3.count.aa$e))/n3.count.diaa$ee
cps.gaa.gac<-((n3.trinuc$gaa*n3.trinuc$gac)/(n3.count.aa$e*n3.count.aa$d))/n3.count.diaa$ed
cps.gaa.gag<-((n3.trinuc$gaa*n3.trinuc$gag)/(n3.count.aa$e*n3.count.aa$e))/n3.count.diaa$ee
cps.gaa.gat<-((n3.trinuc$gaa*n3.trinuc$gat)/(n3.count.aa$e*n3.count.aa$d))/n3.count.diaa$ed

cps.gaa.gca<-((n3.trinuc$gaa*n3.trinuc$gca)/(n3.count.aa$e*n3.count.aa$a))/n3.count.diaa$ea
cps.gaa.gcc<-((n3.trinuc$gaa*n3.trinuc$gcc)/(n3.count.aa$e*n3.count.aa$a))/n3.count.diaa$ea
cps.gaa.gcg<-((n3.trinuc$gaa*n3.trinuc$gcg)/(n3.count.aa$e*n3.count.aa$a))/n3.count.diaa$ea
cps.gaa.gct<-((n3.trinuc$gaa*n3.trinuc$gct)/(n3.count.aa$e*n3.count.aa$a))/n3.count.diaa$ea

cps.gaa.gga<-((n3.trinuc$gaa*n3.trinuc$gga)/(n3.count.aa$e*n3.count.aa$g))/n3.count.diaa$eg
cps.gaa.ggc<-((n3.trinuc$gaa*n3.trinuc$ggc)/(n3.count.aa$e*n3.count.aa$g))/n3.count.diaa$eg
cps.gaa.ggg<-((n3.trinuc$gaa*n3.trinuc$ggg)/(n3.count.aa$e*n3.count.aa$g))/n3.count.diaa$eg
cps.gaa.ggt<-((n3.trinuc$gaa*n3.trinuc$ggt)/(n3.count.aa$e*n3.count.aa$g))/n3.count.diaa$eg

cps.gaa.gta<-((n3.trinuc$gaa*n3.trinuc$gta)/(n3.count.aa$e*n3.count.aa$v))/n3.count.diaa$ev
cps.gaa.gtc<-((n3.trinuc$gaa*n3.trinuc$gtc)/(n3.count.aa$e*n3.count.aa$v))/n3.count.diaa$ev
cps.gaa.gtg<-((n3.trinuc$gaa*n3.trinuc$gtg)/(n3.count.aa$e*n3.count.aa$v))/n3.count.diaa$ev
cps.gaa.gtt<-((n3.trinuc$gaa*n3.trinuc$gtt)/(n3.count.aa$e*n3.count.aa$v))/n3.count.diaa$ev

#Stop codon
#cps.gaa.taa<-((n3.trinuc$gaa*n3.trinuc$taa)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gaa.tac<-((n3.trinuc$gaa*n3.trinuc$tac)/(n3.count.aa$e*n3.count.aa$y))/n3.count.diaa$ey
#Stop codon
#cps.gaa.tag<-((n3.trinuc$gaa*n3.trinuc$tag)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gaa.tat<-((n3.trinuc$gaa*n3.trinuc$tat)/(n3.count.aa$e*n3.count.aa$y))/n3.count.diaa$ey

cps.gaa.tca<-((n3.trinuc$gaa*n3.trinuc$tca)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es
cps.gaa.tcc<-((n3.trinuc$gaa*n3.trinuc$tcc)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es
cps.gaa.tcg<-((n3.trinuc$gaa*n3.trinuc$tcg)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es
cps.gaa.tct<-((n3.trinuc$gaa*n3.trinuc$tct)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es

#Stop codon
#cps.gaa.tga<-((n3.trinuc$gaa*n3.trinuc$tga)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gaa.tgc<-((n3.trinuc$gaa*n3.trinuc$tgc)/(n3.count.aa$e*n3.count.aa$c))/n3.count.diaa$ec
cps.gaa.tgg<-((n3.trinuc$gaa*n3.trinuc$tgg)/(n3.count.aa$e*n3.count.aa$w))/n3.count.diaa$ew
cps.gaa.tgt<-((n3.trinuc$gaa*n3.trinuc$tgt)/(n3.count.aa$e*n3.count.aa$c))/n3.count.diaa$ec

cps.gaa.tta<-((n3.trinuc$gaa*n3.trinuc$tta)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gaa.ttc<-((n3.trinuc$gaa*n3.trinuc$ttc)/(n3.count.aa$e*n3.count.aa$f))/n3.count.diaa$ef
cps.gaa.ttg<-((n3.trinuc$gaa*n3.trinuc$ttg)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gaa.ttt<-((n3.trinuc$gaa*n3.trinuc$ttt)/(n3.count.aa$e*n3.count.aa$f))/n3.count.diaa$ef








cps.gac.aaa<-((n3.trinuc$gac*n3.trinuc$aaa)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gac.aac<-((n3.trinuc$gac*n3.trinuc$aac)/(n3.count.aa$d*n3.count.aa$n))/n3.count.diaa$dn
cps.gac.aag<-((n3.trinuc$gac*n3.trinuc$aag)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gac.aat<-((n3.trinuc$gac*n3.trinuc$aat)/(n3.count.aa$d*n3.count.aa$n))/n3.count.diaa$dn

cps.gac.aca<-((n3.trinuc$gac*n3.trinuc$aca)/(n3.count.aa$d*n3.count.aa$t))/n3.count.diaa$dt
cps.gac.acc<-((n3.trinuc$gac*n3.trinuc$acc)/(n3.count.aa$d*n3.count.aa$t))/n3.count.diaa$dt
cps.gac.acg<-((n3.trinuc$gac*n3.trinuc$acg)/(n3.count.aa$d*n3.count.aa$t))/n3.count.diaa$dt
cps.gac.act<-((n3.trinuc$gac*n3.trinuc$act)/(n3.count.aa$d*n3.count.aa$t))/n3.count.diaa$dt

cps.gac.aga<-((n3.trinuc$gac*n3.trinuc$aga)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gac.agc<-((n3.trinuc$gac*n3.trinuc$agc)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds
cps.gac.agg<-((n3.trinuc$gac*n3.trinuc$agg)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gac.agt<-((n3.trinuc$gac*n3.trinuc$agt)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds

cps.gac.ata<-((n3.trinuc$gac*n3.trinuc$ata)/(n3.count.aa$d*n3.count.aa$i))/n3.count.diaa$di
cps.gac.atc<-((n3.trinuc$gac*n3.trinuc$atc)/(n3.count.aa$d*n3.count.aa$i))/n3.count.diaa$di
cps.gac.atg<-((n3.trinuc$gac*n3.trinuc$atg)/(n3.count.aa$d*n3.count.aa$m))/n3.count.diaa$dm
cps.gac.att<-((n3.trinuc$gac*n3.trinuc$att)/(n3.count.aa$d*n3.count.aa$i))/n3.count.diaa$di

cps.gac.caa<-((n3.trinuc$gac*n3.trinuc$caa)/(n3.count.aa$d*n3.count.aa$q))/n3.count.diaa$dq
cps.gac.cac<-((n3.trinuc$gac*n3.trinuc$cac)/(n3.count.aa$d*n3.count.aa$h))/n3.count.diaa$dh
cps.gac.cag<-((n3.trinuc$gac*n3.trinuc$cag)/(n3.count.aa$d*n3.count.aa$q))/n3.count.diaa$dq
cps.gac.cat<-((n3.trinuc$gac*n3.trinuc$cat)/(n3.count.aa$d*n3.count.aa$h))/n3.count.diaa$dh

cps.gac.cca<-((n3.trinuc$gac*n3.trinuc$cca)/(n3.count.aa$d*n3.count.aa$p))/n3.count.diaa$dp
cps.gac.ccc<-((n3.trinuc$gac*n3.trinuc$ccc)/(n3.count.aa$d*n3.count.aa$p))/n3.count.diaa$dp
cps.gac.ccg<-((n3.trinuc$gac*n3.trinuc$ccg)/(n3.count.aa$d*n3.count.aa$p))/n3.count.diaa$dp
cps.gac.cct<-((n3.trinuc$gac*n3.trinuc$cct)/(n3.count.aa$d*n3.count.aa$p))/n3.count.diaa$dp

cps.gac.cga<-((n3.trinuc$gac*n3.trinuc$cga)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gac.cgc<-((n3.trinuc$gac*n3.trinuc$cgc)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gac.cgg<-((n3.trinuc$gac*n3.trinuc$cgg)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gac.cgt<-((n3.trinuc$gac*n3.trinuc$cgt)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr

cps.gac.cta<-((n3.trinuc$gac*n3.trinuc$cta)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gac.ctc<-((n3.trinuc$gac*n3.trinuc$ctc)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gac.ctg<-((n3.trinuc$gac*n3.trinuc$ctg)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gac.ctt<-((n3.trinuc$gac*n3.trinuc$ctt)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl

cps.gac.gaa<-((n3.trinuc$gac*n3.trinuc$gaa)/(n3.count.aa$d*n3.count.aa$e))/n3.count.diaa$de
cps.gac.gac<-((n3.trinuc$gac*n3.trinuc$gac)/(n3.count.aa$d*n3.count.aa$d))/n3.count.diaa$dd
cps.gac.gag<-((n3.trinuc$gac*n3.trinuc$gag)/(n3.count.aa$d*n3.count.aa$e))/n3.count.diaa$de
cps.gac.gat<-((n3.trinuc$gac*n3.trinuc$gat)/(n3.count.aa$d*n3.count.aa$d))/n3.count.diaa$dd

cps.gac.gca<-((n3.trinuc$gac*n3.trinuc$gca)/(n3.count.aa$d*n3.count.aa$a))/n3.count.diaa$da
cps.gac.gcc<-((n3.trinuc$gac*n3.trinuc$gcc)/(n3.count.aa$d*n3.count.aa$a))/n3.count.diaa$da
cps.gac.gcg<-((n3.trinuc$gac*n3.trinuc$gcg)/(n3.count.aa$d*n3.count.aa$a))/n3.count.diaa$da
cps.gac.gct<-((n3.trinuc$gac*n3.trinuc$gct)/(n3.count.aa$d*n3.count.aa$a))/n3.count.diaa$da

cps.gac.gga<-((n3.trinuc$gac*n3.trinuc$gga)/(n3.count.aa$d*n3.count.aa$g))/n3.count.diaa$dg
cps.gac.ggc<-((n3.trinuc$gac*n3.trinuc$ggc)/(n3.count.aa$d*n3.count.aa$g))/n3.count.diaa$dg
cps.gac.ggg<-((n3.trinuc$gac*n3.trinuc$ggg)/(n3.count.aa$d*n3.count.aa$g))/n3.count.diaa$dg
cps.gac.ggt<-((n3.trinuc$gac*n3.trinuc$ggt)/(n3.count.aa$d*n3.count.aa$g))/n3.count.diaa$dg

cps.gac.gta<-((n3.trinuc$gac*n3.trinuc$gta)/(n3.count.aa$d*n3.count.aa$v))/n3.count.diaa$dv
cps.gac.gtc<-((n3.trinuc$gac*n3.trinuc$gtc)/(n3.count.aa$d*n3.count.aa$v))/n3.count.diaa$dv
cps.gac.gtg<-((n3.trinuc$gac*n3.trinuc$gtg)/(n3.count.aa$d*n3.count.aa$v))/n3.count.diaa$dv
cps.gac.gtt<-((n3.trinuc$gac*n3.trinuc$gtt)/(n3.count.aa$d*n3.count.aa$v))/n3.count.diaa$dv

#Stop codon
#cps.gac.taa<-((n3.trinuc$gac*n3.trinuc$taa)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gac.tac<-((n3.trinuc$gac*n3.trinuc$tac)/(n3.count.aa$d*n3.count.aa$y))/n3.count.diaa$dy
#Stop codon
#cps.gac.tag<-((n3.trinuc$gac*n3.trinuc$tag)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gac.tat<-((n3.trinuc$gac*n3.trinuc$tat)/(n3.count.aa$d*n3.count.aa$y))/n3.count.diaa$dy

cps.gac.tca<-((n3.trinuc$gac*n3.trinuc$tca)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds
cps.gac.tcc<-((n3.trinuc$gac*n3.trinuc$tcc)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds
cps.gac.tcg<-((n3.trinuc$gac*n3.trinuc$tcg)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds
cps.gac.tct<-((n3.trinuc$gac*n3.trinuc$tct)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds

#Stop codon
#cps.gac.tga<-((n3.trinuc$gac*n3.trinuc$tga)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gac.tgc<-((n3.trinuc$gac*n3.trinuc$tgc)/(n3.count.aa$d*n3.count.aa$c))/n3.count.diaa$dc
cps.gac.tgg<-((n3.trinuc$gac*n3.trinuc$tgg)/(n3.count.aa$d*n3.count.aa$w))/n3.count.diaa$dw
cps.gac.tgt<-((n3.trinuc$gac*n3.trinuc$tgt)/(n3.count.aa$d*n3.count.aa$c))/n3.count.diaa$dc

cps.gac.tta<-((n3.trinuc$gac*n3.trinuc$tta)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gac.ttc<-((n3.trinuc$gac*n3.trinuc$ttc)/(n3.count.aa$d*n3.count.aa$f))/n3.count.diaa$df
cps.gac.ttg<-((n3.trinuc$gac*n3.trinuc$ttg)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gac.ttt<-((n3.trinuc$gac*n3.trinuc$ttt)/(n3.count.aa$d*n3.count.aa$f))/n3.count.diaa$df








cps.gag.aaa<-((n3.trinuc$gag*n3.trinuc$aaa)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gag.aac<-((n3.trinuc$gag*n3.trinuc$aac)/(n3.count.aa$e*n3.count.aa$n))/n3.count.diaa$en
cps.gag.aag<-((n3.trinuc$gag*n3.trinuc$aag)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gag.aat<-((n3.trinuc$gag*n3.trinuc$aat)/(n3.count.aa$e*n3.count.aa$n))/n3.count.diaa$en

cps.gag.aca<-((n3.trinuc$gag*n3.trinuc$aca)/(n3.count.aa$e*n3.count.aa$t))/n3.count.diaa$et
cps.gag.acc<-((n3.trinuc$gag*n3.trinuc$acc)/(n3.count.aa$e*n3.count.aa$t))/n3.count.diaa$et
cps.gag.acg<-((n3.trinuc$gag*n3.trinuc$acg)/(n3.count.aa$e*n3.count.aa$t))/n3.count.diaa$et
cps.gag.act<-((n3.trinuc$gag*n3.trinuc$act)/(n3.count.aa$e*n3.count.aa$t))/n3.count.diaa$et

cps.gag.aga<-((n3.trinuc$gag*n3.trinuc$aga)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gag.agc<-((n3.trinuc$gag*n3.trinuc$agc)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es
cps.gag.agg<-((n3.trinuc$gag*n3.trinuc$agg)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gag.agt<-((n3.trinuc$gag*n3.trinuc$agt)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es

cps.gag.ata<-((n3.trinuc$gag*n3.trinuc$ata)/(n3.count.aa$e*n3.count.aa$i))/n3.count.diaa$ei
cps.gag.atc<-((n3.trinuc$gag*n3.trinuc$atc)/(n3.count.aa$e*n3.count.aa$i))/n3.count.diaa$ei
cps.gag.atg<-((n3.trinuc$gag*n3.trinuc$atg)/(n3.count.aa$e*n3.count.aa$m))/n3.count.diaa$em
cps.gag.att<-((n3.trinuc$gag*n3.trinuc$att)/(n3.count.aa$e*n3.count.aa$i))/n3.count.diaa$ei

cps.gag.caa<-((n3.trinuc$gag*n3.trinuc$caa)/(n3.count.aa$e*n3.count.aa$q))/n3.count.diaa$eq
cps.gag.cac<-((n3.trinuc$gag*n3.trinuc$cac)/(n3.count.aa$e*n3.count.aa$h))/n3.count.diaa$eh
cps.gag.cag<-((n3.trinuc$gag*n3.trinuc$cag)/(n3.count.aa$e*n3.count.aa$q))/n3.count.diaa$eq
cps.gag.cat<-((n3.trinuc$gag*n3.trinuc$cat)/(n3.count.aa$e*n3.count.aa$h))/n3.count.diaa$eh

cps.gag.cca<-((n3.trinuc$gag*n3.trinuc$cca)/(n3.count.aa$e*n3.count.aa$p))/n3.count.diaa$ep
cps.gag.ccc<-((n3.trinuc$gag*n3.trinuc$ccc)/(n3.count.aa$e*n3.count.aa$p))/n3.count.diaa$ep
cps.gag.ccg<-((n3.trinuc$gag*n3.trinuc$ccg)/(n3.count.aa$e*n3.count.aa$p))/n3.count.diaa$ep
cps.gag.cct<-((n3.trinuc$gag*n3.trinuc$cct)/(n3.count.aa$e*n3.count.aa$p))/n3.count.diaa$ep

cps.gag.cga<-((n3.trinuc$gag*n3.trinuc$cga)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gag.cgc<-((n3.trinuc$gag*n3.trinuc$cgc)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gag.cgg<-((n3.trinuc$gag*n3.trinuc$cgg)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er
cps.gag.cgt<-((n3.trinuc$gag*n3.trinuc$cgt)/(n3.count.aa$e*n3.count.aa$r))/n3.count.diaa$er

cps.gag.cta<-((n3.trinuc$gag*n3.trinuc$cta)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gag.ctc<-((n3.trinuc$gag*n3.trinuc$ctc)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gag.ctg<-((n3.trinuc$gag*n3.trinuc$ctg)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gag.ctt<-((n3.trinuc$gag*n3.trinuc$ctt)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el

cps.gag.gaa<-((n3.trinuc$gag*n3.trinuc$gaa)/(n3.count.aa$e*n3.count.aa$e))/n3.count.diaa$ee
cps.gag.gac<-((n3.trinuc$gag*n3.trinuc$gac)/(n3.count.aa$e*n3.count.aa$d))/n3.count.diaa$ed
cps.gag.gag<-((n3.trinuc$gag*n3.trinuc$gag)/(n3.count.aa$e*n3.count.aa$e))/n3.count.diaa$ee
cps.gag.gat<-((n3.trinuc$gag*n3.trinuc$gat)/(n3.count.aa$e*n3.count.aa$d))/n3.count.diaa$ed

cps.gag.gca<-((n3.trinuc$gag*n3.trinuc$gca)/(n3.count.aa$e*n3.count.aa$a))/n3.count.diaa$ea
cps.gag.gcc<-((n3.trinuc$gag*n3.trinuc$gcc)/(n3.count.aa$e*n3.count.aa$a))/n3.count.diaa$ea
cps.gag.gcg<-((n3.trinuc$gag*n3.trinuc$gcg)/(n3.count.aa$e*n3.count.aa$a))/n3.count.diaa$ea
cps.gag.gct<-((n3.trinuc$gag*n3.trinuc$gct)/(n3.count.aa$e*n3.count.aa$a))/n3.count.diaa$ea

cps.gag.gga<-((n3.trinuc$gag*n3.trinuc$gga)/(n3.count.aa$e*n3.count.aa$g))/n3.count.diaa$eg
cps.gag.ggc<-((n3.trinuc$gag*n3.trinuc$ggc)/(n3.count.aa$e*n3.count.aa$g))/n3.count.diaa$eg
cps.gag.ggg<-((n3.trinuc$gag*n3.trinuc$ggg)/(n3.count.aa$e*n3.count.aa$g))/n3.count.diaa$eg
cps.gag.ggt<-((n3.trinuc$gag*n3.trinuc$ggt)/(n3.count.aa$e*n3.count.aa$g))/n3.count.diaa$eg

cps.gag.gta<-((n3.trinuc$gag*n3.trinuc$gta)/(n3.count.aa$e*n3.count.aa$v))/n3.count.diaa$ev
cps.gag.gtc<-((n3.trinuc$gag*n3.trinuc$gtc)/(n3.count.aa$e*n3.count.aa$v))/n3.count.diaa$ev
cps.gag.gtg<-((n3.trinuc$gag*n3.trinuc$gtg)/(n3.count.aa$e*n3.count.aa$v))/n3.count.diaa$ev
cps.gag.gtt<-((n3.trinuc$gag*n3.trinuc$gtt)/(n3.count.aa$e*n3.count.aa$v))/n3.count.diaa$ev

#Stop codon
#cps.gag.taa<-((n3.trinuc$gag*n3.trinuc$taa)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gag.tac<-((n3.trinuc$gag*n3.trinuc$tac)/(n3.count.aa$e*n3.count.aa$y))/n3.count.diaa$ey
#Stop codon
#cps.gag.tag<-((n3.trinuc$gag*n3.trinuc$tag)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gag.tat<-((n3.trinuc$gag*n3.trinuc$tat)/(n3.count.aa$e*n3.count.aa$y))/n3.count.diaa$ey

cps.gag.tca<-((n3.trinuc$gag*n3.trinuc$tca)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es
cps.gag.tcc<-((n3.trinuc$gag*n3.trinuc$tcc)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es
cps.gag.tcg<-((n3.trinuc$gag*n3.trinuc$tcg)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es
cps.gag.tct<-((n3.trinuc$gag*n3.trinuc$tct)/(n3.count.aa$e*n3.count.aa$s))/n3.count.diaa$es

#Stop codon
#cps.gag.tga<-((n3.trinuc$gag*n3.trinuc$tga)/(n3.count.aa$e*n3.count.aa$k))/n3.count.diaa$ek
cps.gag.tgc<-((n3.trinuc$gag*n3.trinuc$tgc)/(n3.count.aa$e*n3.count.aa$c))/n3.count.diaa$ec
cps.gag.tgg<-((n3.trinuc$gag*n3.trinuc$tgg)/(n3.count.aa$e*n3.count.aa$w))/n3.count.diaa$ew
cps.gag.tgt<-((n3.trinuc$gag*n3.trinuc$tgt)/(n3.count.aa$e*n3.count.aa$c))/n3.count.diaa$ec

cps.gag.tta<-((n3.trinuc$gag*n3.trinuc$tta)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gag.ttc<-((n3.trinuc$gag*n3.trinuc$ttc)/(n3.count.aa$e*n3.count.aa$f))/n3.count.diaa$ef
cps.gag.ttg<-((n3.trinuc$gag*n3.trinuc$ttg)/(n3.count.aa$e*n3.count.aa$l))/n3.count.diaa$el
cps.gag.ttt<-((n3.trinuc$gag*n3.trinuc$ttt)/(n3.count.aa$e*n3.count.aa$f))/n3.count.diaa$ef








cps.gat.aaa<-((n3.trinuc$gat*n3.trinuc$aaa)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gat.aac<-((n3.trinuc$gat*n3.trinuc$aac)/(n3.count.aa$d*n3.count.aa$n))/n3.count.diaa$dn
cps.gat.aag<-((n3.trinuc$gat*n3.trinuc$aag)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gat.aat<-((n3.trinuc$gat*n3.trinuc$aat)/(n3.count.aa$d*n3.count.aa$n))/n3.count.diaa$dn

cps.gat.aca<-((n3.trinuc$gat*n3.trinuc$aca)/(n3.count.aa$d*n3.count.aa$t))/n3.count.diaa$dt
cps.gat.acc<-((n3.trinuc$gat*n3.trinuc$acc)/(n3.count.aa$d*n3.count.aa$t))/n3.count.diaa$dt
cps.gat.acg<-((n3.trinuc$gat*n3.trinuc$acg)/(n3.count.aa$d*n3.count.aa$t))/n3.count.diaa$dt
cps.gat.act<-((n3.trinuc$gat*n3.trinuc$act)/(n3.count.aa$d*n3.count.aa$t))/n3.count.diaa$dt

cps.gat.aga<-((n3.trinuc$gat*n3.trinuc$aga)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gat.agc<-((n3.trinuc$gat*n3.trinuc$agc)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds
cps.gat.agg<-((n3.trinuc$gat*n3.trinuc$agg)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gat.agt<-((n3.trinuc$gat*n3.trinuc$agt)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds

cps.gat.ata<-((n3.trinuc$gat*n3.trinuc$ata)/(n3.count.aa$d*n3.count.aa$i))/n3.count.diaa$di
cps.gat.atc<-((n3.trinuc$gat*n3.trinuc$atc)/(n3.count.aa$d*n3.count.aa$i))/n3.count.diaa$di
cps.gat.atg<-((n3.trinuc$gat*n3.trinuc$atg)/(n3.count.aa$d*n3.count.aa$m))/n3.count.diaa$dm
cps.gat.att<-((n3.trinuc$gat*n3.trinuc$att)/(n3.count.aa$d*n3.count.aa$i))/n3.count.diaa$di

cps.gat.caa<-((n3.trinuc$gat*n3.trinuc$caa)/(n3.count.aa$d*n3.count.aa$q))/n3.count.diaa$dq
cps.gat.cac<-((n3.trinuc$gat*n3.trinuc$cac)/(n3.count.aa$d*n3.count.aa$h))/n3.count.diaa$dh
cps.gat.cag<-((n3.trinuc$gat*n3.trinuc$cag)/(n3.count.aa$d*n3.count.aa$q))/n3.count.diaa$dq
cps.gat.cat<-((n3.trinuc$gat*n3.trinuc$cat)/(n3.count.aa$d*n3.count.aa$h))/n3.count.diaa$dh

cps.gat.cca<-((n3.trinuc$gat*n3.trinuc$cca)/(n3.count.aa$d*n3.count.aa$p))/n3.count.diaa$dp
cps.gat.ccc<-((n3.trinuc$gat*n3.trinuc$ccc)/(n3.count.aa$d*n3.count.aa$p))/n3.count.diaa$dp
cps.gat.ccg<-((n3.trinuc$gat*n3.trinuc$ccg)/(n3.count.aa$d*n3.count.aa$p))/n3.count.diaa$dp
cps.gat.cct<-((n3.trinuc$gat*n3.trinuc$cct)/(n3.count.aa$d*n3.count.aa$p))/n3.count.diaa$dp

cps.gat.cga<-((n3.trinuc$gat*n3.trinuc$cga)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gat.cgc<-((n3.trinuc$gat*n3.trinuc$cgc)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gat.cgg<-((n3.trinuc$gat*n3.trinuc$cgg)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr
cps.gat.cgt<-((n3.trinuc$gat*n3.trinuc$cgt)/(n3.count.aa$d*n3.count.aa$r))/n3.count.diaa$dr

cps.gat.cta<-((n3.trinuc$gat*n3.trinuc$cta)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gat.ctc<-((n3.trinuc$gat*n3.trinuc$ctc)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gat.ctg<-((n3.trinuc$gat*n3.trinuc$ctg)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gat.ctt<-((n3.trinuc$gat*n3.trinuc$ctt)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl

cps.gat.gaa<-((n3.trinuc$gat*n3.trinuc$gaa)/(n3.count.aa$d*n3.count.aa$e))/n3.count.diaa$de
cps.gat.gac<-((n3.trinuc$gat*n3.trinuc$gac)/(n3.count.aa$d*n3.count.aa$d))/n3.count.diaa$dd
cps.gat.gag<-((n3.trinuc$gat*n3.trinuc$gag)/(n3.count.aa$d*n3.count.aa$e))/n3.count.diaa$de
cps.gat.gat<-((n3.trinuc$gat*n3.trinuc$gat)/(n3.count.aa$d*n3.count.aa$d))/n3.count.diaa$dd

cps.gat.gca<-((n3.trinuc$gat*n3.trinuc$gca)/(n3.count.aa$d*n3.count.aa$a))/n3.count.diaa$da
cps.gat.gcc<-((n3.trinuc$gat*n3.trinuc$gcc)/(n3.count.aa$d*n3.count.aa$a))/n3.count.diaa$da
cps.gat.gcg<-((n3.trinuc$gat*n3.trinuc$gcg)/(n3.count.aa$d*n3.count.aa$a))/n3.count.diaa$da
cps.gat.gct<-((n3.trinuc$gat*n3.trinuc$gct)/(n3.count.aa$d*n3.count.aa$a))/n3.count.diaa$da

cps.gat.gga<-((n3.trinuc$gat*n3.trinuc$gga)/(n3.count.aa$d*n3.count.aa$g))/n3.count.diaa$dg
cps.gat.ggc<-((n3.trinuc$gat*n3.trinuc$ggc)/(n3.count.aa$d*n3.count.aa$g))/n3.count.diaa$dg
cps.gat.ggg<-((n3.trinuc$gat*n3.trinuc$ggg)/(n3.count.aa$d*n3.count.aa$g))/n3.count.diaa$dg
cps.gat.ggt<-((n3.trinuc$gat*n3.trinuc$ggt)/(n3.count.aa$d*n3.count.aa$g))/n3.count.diaa$dg

cps.gat.gta<-((n3.trinuc$gat*n3.trinuc$gta)/(n3.count.aa$d*n3.count.aa$v))/n3.count.diaa$dv
cps.gat.gtc<-((n3.trinuc$gat*n3.trinuc$gtc)/(n3.count.aa$d*n3.count.aa$v))/n3.count.diaa$dv
cps.gat.gtg<-((n3.trinuc$gat*n3.trinuc$gtg)/(n3.count.aa$d*n3.count.aa$v))/n3.count.diaa$dv
cps.gat.gtt<-((n3.trinuc$gat*n3.trinuc$gtt)/(n3.count.aa$d*n3.count.aa$v))/n3.count.diaa$dv

#Stop codon
#cps.gat.taa<-((n3.trinuc$gat*n3.trinuc$taa)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gat.tac<-((n3.trinuc$gat*n3.trinuc$tac)/(n3.count.aa$d*n3.count.aa$y))/n3.count.diaa$dy
#Stop codon
#cps.gat.tag<-((n3.trinuc$gat*n3.trinuc$tag)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gat.tat<-((n3.trinuc$gat*n3.trinuc$tat)/(n3.count.aa$d*n3.count.aa$y))/n3.count.diaa$dy

cps.gat.tca<-((n3.trinuc$gat*n3.trinuc$tca)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds
cps.gat.tcc<-((n3.trinuc$gat*n3.trinuc$tcc)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds
cps.gat.tcg<-((n3.trinuc$gat*n3.trinuc$tcg)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds
cps.gat.tct<-((n3.trinuc$gat*n3.trinuc$tct)/(n3.count.aa$d*n3.count.aa$s))/n3.count.diaa$ds

#Stop codon
#cps.gat.tga<-((n3.trinuc$gat*n3.trinuc$tga)/(n3.count.aa$d*n3.count.aa$k))/n3.count.diaa$dk
cps.gat.tgc<-((n3.trinuc$gat*n3.trinuc$tgc)/(n3.count.aa$d*n3.count.aa$c))/n3.count.diaa$dc
cps.gat.tgg<-((n3.trinuc$gat*n3.trinuc$tgg)/(n3.count.aa$d*n3.count.aa$w))/n3.count.diaa$dw
cps.gat.tgt<-((n3.trinuc$gat*n3.trinuc$tgt)/(n3.count.aa$d*n3.count.aa$c))/n3.count.diaa$dc

cps.gat.tta<-((n3.trinuc$gat*n3.trinuc$tta)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gat.ttc<-((n3.trinuc$gat*n3.trinuc$ttc)/(n3.count.aa$d*n3.count.aa$f))/n3.count.diaa$df
cps.gat.ttg<-((n3.trinuc$gat*n3.trinuc$ttg)/(n3.count.aa$d*n3.count.aa$l))/n3.count.diaa$dl
cps.gat.ttt<-((n3.trinuc$gat*n3.trinuc$ttt)/(n3.count.aa$d*n3.count.aa$f))/n3.count.diaa$df








cps.gca.aaa<-((n3.trinuc$gca*n3.trinuc$aaa)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gca.aac<-((n3.trinuc$gca*n3.trinuc$aac)/(n3.count.aa$a*n3.count.aa$n))/n3.count.diaa$an
cps.gca.aag<-((n3.trinuc$gca*n3.trinuc$aag)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gca.aat<-((n3.trinuc$gca*n3.trinuc$aat)/(n3.count.aa$a*n3.count.aa$n))/n3.count.diaa$an

cps.gca.aca<-((n3.trinuc$gca*n3.trinuc$aca)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gca.acc<-((n3.trinuc$gca*n3.trinuc$acc)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gca.acg<-((n3.trinuc$gca*n3.trinuc$acg)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gca.act<-((n3.trinuc$gca*n3.trinuc$act)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at

cps.gca.aga<-((n3.trinuc$gca*n3.trinuc$aga)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gca.agc<-((n3.trinuc$gca*n3.trinuc$agc)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gca.agg<-((n3.trinuc$gca*n3.trinuc$agg)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gca.agt<-((n3.trinuc$gca*n3.trinuc$agt)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as

cps.gca.ata<-((n3.trinuc$gca*n3.trinuc$ata)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai
cps.gca.atc<-((n3.trinuc$gca*n3.trinuc$atc)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai
cps.gca.atg<-((n3.trinuc$gca*n3.trinuc$atg)/(n3.count.aa$a*n3.count.aa$m))/n3.count.diaa$am
cps.gca.att<-((n3.trinuc$gca*n3.trinuc$att)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai

cps.gca.caa<-((n3.trinuc$gca*n3.trinuc$caa)/(n3.count.aa$a*n3.count.aa$q))/n3.count.diaa$aq
cps.gca.cac<-((n3.trinuc$gca*n3.trinuc$cac)/(n3.count.aa$a*n3.count.aa$h))/n3.count.diaa$ah
cps.gca.cag<-((n3.trinuc$gca*n3.trinuc$cag)/(n3.count.aa$a*n3.count.aa$q))/n3.count.diaa$aq
cps.gca.cat<-((n3.trinuc$gca*n3.trinuc$cat)/(n3.count.aa$a*n3.count.aa$h))/n3.count.diaa$ah

cps.gca.cca<-((n3.trinuc$gca*n3.trinuc$cca)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gca.ccc<-((n3.trinuc$gca*n3.trinuc$ccc)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gca.ccg<-((n3.trinuc$gca*n3.trinuc$ccg)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gca.cct<-((n3.trinuc$gca*n3.trinuc$cct)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap

cps.gca.cga<-((n3.trinuc$gca*n3.trinuc$cga)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gca.cgc<-((n3.trinuc$gca*n3.trinuc$cgc)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gca.cgg<-((n3.trinuc$gca*n3.trinuc$cgg)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gca.cgt<-((n3.trinuc$gca*n3.trinuc$cgt)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar

cps.gca.cta<-((n3.trinuc$gca*n3.trinuc$cta)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gca.ctc<-((n3.trinuc$gca*n3.trinuc$ctc)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gca.ctg<-((n3.trinuc$gca*n3.trinuc$ctg)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gca.ctt<-((n3.trinuc$gca*n3.trinuc$ctt)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al

cps.gca.gaa<-((n3.trinuc$gca*n3.trinuc$gaa)/(n3.count.aa$a*n3.count.aa$e))/n3.count.diaa$ae
cps.gca.gac<-((n3.trinuc$gca*n3.trinuc$gac)/(n3.count.aa$a*n3.count.aa$d))/n3.count.diaa$ad
cps.gca.gag<-((n3.trinuc$gca*n3.trinuc$gag)/(n3.count.aa$a*n3.count.aa$e))/n3.count.diaa$ae
cps.gca.gat<-((n3.trinuc$gca*n3.trinuc$gat)/(n3.count.aa$a*n3.count.aa$d))/n3.count.diaa$ad

cps.gca.gca<-((n3.trinuc$gca*n3.trinuc$gca)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gca.gcc<-((n3.trinuc$gca*n3.trinuc$gcc)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gca.gcg<-((n3.trinuc$gca*n3.trinuc$gcg)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gca.gct<-((n3.trinuc$gca*n3.trinuc$gct)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa

cps.gca.gga<-((n3.trinuc$gca*n3.trinuc$gga)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gca.ggc<-((n3.trinuc$gca*n3.trinuc$ggc)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gca.ggg<-((n3.trinuc$gca*n3.trinuc$ggg)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gca.ggt<-((n3.trinuc$gca*n3.trinuc$ggt)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag

cps.gca.gta<-((n3.trinuc$gca*n3.trinuc$gta)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gca.gtc<-((n3.trinuc$gca*n3.trinuc$gtc)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gca.gtg<-((n3.trinuc$gca*n3.trinuc$gtg)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gca.gtt<-((n3.trinuc$gca*n3.trinuc$gtt)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av

#Stop codon
#cps.gca.taa<-((n3.trinuc$gca*n3.trinuc$taa)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gca.tac<-((n3.trinuc$gca*n3.trinuc$tac)/(n3.count.aa$a*n3.count.aa$y))/n3.count.diaa$ay
#Stop codon
#cps.gca.tag<-((n3.trinuc$gca*n3.trinuc$tag)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gca.tat<-((n3.trinuc$gca*n3.trinuc$tat)/(n3.count.aa$a*n3.count.aa$y))/n3.count.diaa$ay

cps.gca.tca<-((n3.trinuc$gca*n3.trinuc$tca)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gca.tcc<-((n3.trinuc$gca*n3.trinuc$tcc)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gca.tcg<-((n3.trinuc$gca*n3.trinuc$tcg)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gca.tct<-((n3.trinuc$gca*n3.trinuc$tct)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as

#Stop codon
#cps.gca.tga<-((n3.trinuc$gca*n3.trinuc$tga)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gca.tgc<-((n3.trinuc$gca*n3.trinuc$tgc)/(n3.count.aa$a*n3.count.aa$c))/n3.count.diaa$ac
cps.gca.tgg<-((n3.trinuc$gca*n3.trinuc$tgg)/(n3.count.aa$a*n3.count.aa$w))/n3.count.diaa$aw
cps.gca.tgt<-((n3.trinuc$gca*n3.trinuc$tgt)/(n3.count.aa$a*n3.count.aa$c))/n3.count.diaa$ac

cps.gca.tta<-((n3.trinuc$gca*n3.trinuc$tta)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gca.ttc<-((n3.trinuc$gca*n3.trinuc$ttc)/(n3.count.aa$a*n3.count.aa$f))/n3.count.diaa$af
cps.gca.ttg<-((n3.trinuc$gca*n3.trinuc$ttg)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gca.ttt<-((n3.trinuc$gca*n3.trinuc$ttt)/(n3.count.aa$a*n3.count.aa$f))/n3.count.diaa$af








cps.gcc.aaa<-((n3.trinuc$gcc*n3.trinuc$aaa)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcc.aac<-((n3.trinuc$gcc*n3.trinuc$aac)/(n3.count.aa$a*n3.count.aa$n))/n3.count.diaa$an
cps.gcc.aag<-((n3.trinuc$gcc*n3.trinuc$aag)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcc.aat<-((n3.trinuc$gcc*n3.trinuc$aat)/(n3.count.aa$a*n3.count.aa$n))/n3.count.diaa$an

cps.gcc.aca<-((n3.trinuc$gcc*n3.trinuc$aca)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gcc.acc<-((n3.trinuc$gcc*n3.trinuc$acc)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gcc.acg<-((n3.trinuc$gcc*n3.trinuc$acg)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gcc.act<-((n3.trinuc$gcc*n3.trinuc$act)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at

cps.gcc.aga<-((n3.trinuc$gcc*n3.trinuc$aga)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcc.agc<-((n3.trinuc$gcc*n3.trinuc$agc)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gcc.agg<-((n3.trinuc$gcc*n3.trinuc$agg)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcc.agt<-((n3.trinuc$gcc*n3.trinuc$agt)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as

cps.gcc.ata<-((n3.trinuc$gcc*n3.trinuc$ata)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai
cps.gcc.atc<-((n3.trinuc$gcc*n3.trinuc$atc)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai
cps.gcc.atg<-((n3.trinuc$gcc*n3.trinuc$atg)/(n3.count.aa$a*n3.count.aa$m))/n3.count.diaa$am
cps.gcc.att<-((n3.trinuc$gcc*n3.trinuc$att)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai

cps.gcc.caa<-((n3.trinuc$gcc*n3.trinuc$caa)/(n3.count.aa$a*n3.count.aa$q))/n3.count.diaa$aq
cps.gcc.cac<-((n3.trinuc$gcc*n3.trinuc$cac)/(n3.count.aa$a*n3.count.aa$h))/n3.count.diaa$ah
cps.gcc.cag<-((n3.trinuc$gcc*n3.trinuc$cag)/(n3.count.aa$a*n3.count.aa$q))/n3.count.diaa$aq
cps.gcc.cat<-((n3.trinuc$gcc*n3.trinuc$cat)/(n3.count.aa$a*n3.count.aa$h))/n3.count.diaa$ah

cps.gcc.cca<-((n3.trinuc$gcc*n3.trinuc$cca)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gcc.ccc<-((n3.trinuc$gcc*n3.trinuc$ccc)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gcc.ccg<-((n3.trinuc$gcc*n3.trinuc$ccg)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gcc.cct<-((n3.trinuc$gcc*n3.trinuc$cct)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap

cps.gcc.cga<-((n3.trinuc$gcc*n3.trinuc$cga)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcc.cgc<-((n3.trinuc$gcc*n3.trinuc$cgc)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcc.cgg<-((n3.trinuc$gcc*n3.trinuc$cgg)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcc.cgt<-((n3.trinuc$gcc*n3.trinuc$cgt)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar

cps.gcc.cta<-((n3.trinuc$gcc*n3.trinuc$cta)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcc.ctc<-((n3.trinuc$gcc*n3.trinuc$ctc)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcc.ctg<-((n3.trinuc$gcc*n3.trinuc$ctg)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcc.ctt<-((n3.trinuc$gcc*n3.trinuc$ctt)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al

cps.gcc.gaa<-((n3.trinuc$gcc*n3.trinuc$gaa)/(n3.count.aa$a*n3.count.aa$e))/n3.count.diaa$ae
cps.gcc.gac<-((n3.trinuc$gcc*n3.trinuc$gac)/(n3.count.aa$a*n3.count.aa$d))/n3.count.diaa$ad
cps.gcc.gag<-((n3.trinuc$gcc*n3.trinuc$gag)/(n3.count.aa$a*n3.count.aa$e))/n3.count.diaa$ae
cps.gcc.gat<-((n3.trinuc$gcc*n3.trinuc$gat)/(n3.count.aa$a*n3.count.aa$d))/n3.count.diaa$ad

cps.gcc.gca<-((n3.trinuc$gcc*n3.trinuc$gca)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gcc.gcc<-((n3.trinuc$gcc*n3.trinuc$gcc)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gcc.gcg<-((n3.trinuc$gcc*n3.trinuc$gcg)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gcc.gct<-((n3.trinuc$gcc*n3.trinuc$gct)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa

cps.gcc.gga<-((n3.trinuc$gcc*n3.trinuc$gga)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gcc.ggc<-((n3.trinuc$gcc*n3.trinuc$ggc)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gcc.ggg<-((n3.trinuc$gcc*n3.trinuc$ggg)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gcc.ggt<-((n3.trinuc$gcc*n3.trinuc$ggt)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag

cps.gcc.gta<-((n3.trinuc$gcc*n3.trinuc$gta)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gcc.gtc<-((n3.trinuc$gcc*n3.trinuc$gtc)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gcc.gtg<-((n3.trinuc$gcc*n3.trinuc$gtg)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gcc.gtt<-((n3.trinuc$gcc*n3.trinuc$gtt)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av

#Stop codon
#cps.gcc.taa<-((n3.trinuc$gcc*n3.trinuc$taa)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcc.tac<-((n3.trinuc$gcc*n3.trinuc$tac)/(n3.count.aa$a*n3.count.aa$y))/n3.count.diaa$ay
#Stop codon
#cps.gcc.tag<-((n3.trinuc$gcc*n3.trinuc$tag)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcc.tat<-((n3.trinuc$gcc*n3.trinuc$tat)/(n3.count.aa$a*n3.count.aa$y))/n3.count.diaa$ay

cps.gcc.tca<-((n3.trinuc$gcc*n3.trinuc$tca)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gcc.tcc<-((n3.trinuc$gcc*n3.trinuc$tcc)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gcc.tcg<-((n3.trinuc$gcc*n3.trinuc$tcg)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gcc.tct<-((n3.trinuc$gcc*n3.trinuc$tct)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as

#Stop codon
#cps.gcc.tga<-((n3.trinuc$gcc*n3.trinuc$tga)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcc.tgc<-((n3.trinuc$gcc*n3.trinuc$tgc)/(n3.count.aa$a*n3.count.aa$c))/n3.count.diaa$ac
cps.gcc.tgg<-((n3.trinuc$gcc*n3.trinuc$tgg)/(n3.count.aa$a*n3.count.aa$w))/n3.count.diaa$aw
cps.gcc.tgt<-((n3.trinuc$gcc*n3.trinuc$tgt)/(n3.count.aa$a*n3.count.aa$c))/n3.count.diaa$ac

cps.gcc.tta<-((n3.trinuc$gcc*n3.trinuc$tta)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcc.ttc<-((n3.trinuc$gcc*n3.trinuc$ttc)/(n3.count.aa$a*n3.count.aa$f))/n3.count.diaa$af
cps.gcc.ttg<-((n3.trinuc$gcc*n3.trinuc$ttg)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcc.ttt<-((n3.trinuc$gcc*n3.trinuc$ttt)/(n3.count.aa$a*n3.count.aa$f))/n3.count.diaa$af








cps.gcg.aaa<-((n3.trinuc$gcg*n3.trinuc$aaa)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcg.aac<-((n3.trinuc$gcg*n3.trinuc$aac)/(n3.count.aa$a*n3.count.aa$n))/n3.count.diaa$an
cps.gcg.aag<-((n3.trinuc$gcg*n3.trinuc$aag)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcg.aat<-((n3.trinuc$gcg*n3.trinuc$aat)/(n3.count.aa$a*n3.count.aa$n))/n3.count.diaa$an

cps.gcg.aca<-((n3.trinuc$gcg*n3.trinuc$aca)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gcg.acc<-((n3.trinuc$gcg*n3.trinuc$acc)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gcg.acg<-((n3.trinuc$gcg*n3.trinuc$acg)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gcg.act<-((n3.trinuc$gcg*n3.trinuc$act)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at

cps.gcg.aga<-((n3.trinuc$gcg*n3.trinuc$aga)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcg.agc<-((n3.trinuc$gcg*n3.trinuc$agc)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gcg.agg<-((n3.trinuc$gcg*n3.trinuc$agg)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcg.agt<-((n3.trinuc$gcg*n3.trinuc$agt)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as

cps.gcg.ata<-((n3.trinuc$gcg*n3.trinuc$ata)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai
cps.gcg.atc<-((n3.trinuc$gcg*n3.trinuc$atc)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai
cps.gcg.atg<-((n3.trinuc$gcg*n3.trinuc$atg)/(n3.count.aa$a*n3.count.aa$m))/n3.count.diaa$am
cps.gcg.att<-((n3.trinuc$gcg*n3.trinuc$att)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai

cps.gcg.caa<-((n3.trinuc$gcg*n3.trinuc$caa)/(n3.count.aa$a*n3.count.aa$q))/n3.count.diaa$aq
cps.gcg.cac<-((n3.trinuc$gcg*n3.trinuc$cac)/(n3.count.aa$a*n3.count.aa$h))/n3.count.diaa$ah
cps.gcg.cag<-((n3.trinuc$gcg*n3.trinuc$cag)/(n3.count.aa$a*n3.count.aa$q))/n3.count.diaa$aq
cps.gcg.cat<-((n3.trinuc$gcg*n3.trinuc$cat)/(n3.count.aa$a*n3.count.aa$h))/n3.count.diaa$ah

cps.gcg.cca<-((n3.trinuc$gcg*n3.trinuc$cca)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gcg.ccc<-((n3.trinuc$gcg*n3.trinuc$ccc)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gcg.ccg<-((n3.trinuc$gcg*n3.trinuc$ccg)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gcg.cct<-((n3.trinuc$gcg*n3.trinuc$cct)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap

cps.gcg.cga<-((n3.trinuc$gcg*n3.trinuc$cga)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcg.cgc<-((n3.trinuc$gcg*n3.trinuc$cgc)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcg.cgg<-((n3.trinuc$gcg*n3.trinuc$cgg)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gcg.cgt<-((n3.trinuc$gcg*n3.trinuc$cgt)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar

cps.gcg.cta<-((n3.trinuc$gcg*n3.trinuc$cta)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcg.ctc<-((n3.trinuc$gcg*n3.trinuc$ctc)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcg.ctg<-((n3.trinuc$gcg*n3.trinuc$ctg)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcg.ctt<-((n3.trinuc$gcg*n3.trinuc$ctt)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al

cps.gcg.gaa<-((n3.trinuc$gcg*n3.trinuc$gaa)/(n3.count.aa$a*n3.count.aa$e))/n3.count.diaa$ae
cps.gcg.gac<-((n3.trinuc$gcg*n3.trinuc$gac)/(n3.count.aa$a*n3.count.aa$d))/n3.count.diaa$ad
cps.gcg.gag<-((n3.trinuc$gcg*n3.trinuc$gag)/(n3.count.aa$a*n3.count.aa$e))/n3.count.diaa$ae
cps.gcg.gat<-((n3.trinuc$gcg*n3.trinuc$gat)/(n3.count.aa$a*n3.count.aa$d))/n3.count.diaa$ad

cps.gcg.gca<-((n3.trinuc$gcg*n3.trinuc$gca)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gcg.gcc<-((n3.trinuc$gcg*n3.trinuc$gcc)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gcg.gcg<-((n3.trinuc$gcg*n3.trinuc$gcg)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gcg.gct<-((n3.trinuc$gcg*n3.trinuc$gct)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa

cps.gcg.gga<-((n3.trinuc$gcg*n3.trinuc$gga)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gcg.ggc<-((n3.trinuc$gcg*n3.trinuc$ggc)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gcg.ggg<-((n3.trinuc$gcg*n3.trinuc$ggg)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gcg.ggt<-((n3.trinuc$gcg*n3.trinuc$ggt)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag

cps.gcg.gta<-((n3.trinuc$gcg*n3.trinuc$gta)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gcg.gtc<-((n3.trinuc$gcg*n3.trinuc$gtc)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gcg.gtg<-((n3.trinuc$gcg*n3.trinuc$gtg)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gcg.gtt<-((n3.trinuc$gcg*n3.trinuc$gtt)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av

#Stop codon
#cps.gcg.taa<-((n3.trinuc$gcg*n3.trinuc$taa)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcg.tac<-((n3.trinuc$gcg*n3.trinuc$tac)/(n3.count.aa$a*n3.count.aa$y))/n3.count.diaa$ay
#Stop codon
#cps.gcg.tag<-((n3.trinuc$gcg*n3.trinuc$tag)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcg.tat<-((n3.trinuc$gcg*n3.trinuc$tat)/(n3.count.aa$a*n3.count.aa$y))/n3.count.diaa$ay

cps.gcg.tca<-((n3.trinuc$gcg*n3.trinuc$tca)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gcg.tcc<-((n3.trinuc$gcg*n3.trinuc$tcc)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gcg.tcg<-((n3.trinuc$gcg*n3.trinuc$tcg)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gcg.tct<-((n3.trinuc$gcg*n3.trinuc$tct)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as

#Stop codon
#cps.gcg.tga<-((n3.trinuc$gcg*n3.trinuc$tga)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gcg.tgc<-((n3.trinuc$gcg*n3.trinuc$tgc)/(n3.count.aa$a*n3.count.aa$c))/n3.count.diaa$ac
cps.gcg.tgg<-((n3.trinuc$gcg*n3.trinuc$tgg)/(n3.count.aa$a*n3.count.aa$w))/n3.count.diaa$aw
cps.gcg.tgt<-((n3.trinuc$gcg*n3.trinuc$tgt)/(n3.count.aa$a*n3.count.aa$c))/n3.count.diaa$ac

cps.gcg.tta<-((n3.trinuc$gcg*n3.trinuc$tta)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcg.ttc<-((n3.trinuc$gcg*n3.trinuc$ttc)/(n3.count.aa$a*n3.count.aa$f))/n3.count.diaa$af
cps.gcg.ttg<-((n3.trinuc$gcg*n3.trinuc$ttg)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gcg.ttt<-((n3.trinuc$gcg*n3.trinuc$ttt)/(n3.count.aa$a*n3.count.aa$f))/n3.count.diaa$af








cps.gct.aaa<-((n3.trinuc$gct*n3.trinuc$aaa)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gct.aac<-((n3.trinuc$gct*n3.trinuc$aac)/(n3.count.aa$a*n3.count.aa$n))/n3.count.diaa$an
cps.gct.aag<-((n3.trinuc$gct*n3.trinuc$aag)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gct.aat<-((n3.trinuc$gct*n3.trinuc$aat)/(n3.count.aa$a*n3.count.aa$n))/n3.count.diaa$an

cps.gct.aca<-((n3.trinuc$gct*n3.trinuc$aca)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gct.acc<-((n3.trinuc$gct*n3.trinuc$acc)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gct.acg<-((n3.trinuc$gct*n3.trinuc$acg)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at
cps.gct.act<-((n3.trinuc$gct*n3.trinuc$act)/(n3.count.aa$a*n3.count.aa$t))/n3.count.diaa$at

cps.gct.aga<-((n3.trinuc$gct*n3.trinuc$aga)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gct.agc<-((n3.trinuc$gct*n3.trinuc$agc)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gct.agg<-((n3.trinuc$gct*n3.trinuc$agg)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gct.agt<-((n3.trinuc$gct*n3.trinuc$agt)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as

cps.gct.ata<-((n3.trinuc$gct*n3.trinuc$ata)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai
cps.gct.atc<-((n3.trinuc$gct*n3.trinuc$atc)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai
cps.gct.atg<-((n3.trinuc$gct*n3.trinuc$atg)/(n3.count.aa$a*n3.count.aa$m))/n3.count.diaa$am
cps.gct.att<-((n3.trinuc$gct*n3.trinuc$att)/(n3.count.aa$a*n3.count.aa$i))/n3.count.diaa$ai

cps.gct.caa<-((n3.trinuc$gct*n3.trinuc$caa)/(n3.count.aa$a*n3.count.aa$q))/n3.count.diaa$aq
cps.gct.cac<-((n3.trinuc$gct*n3.trinuc$cac)/(n3.count.aa$a*n3.count.aa$h))/n3.count.diaa$ah
cps.gct.cag<-((n3.trinuc$gct*n3.trinuc$cag)/(n3.count.aa$a*n3.count.aa$q))/n3.count.diaa$aq
cps.gct.cat<-((n3.trinuc$gct*n3.trinuc$cat)/(n3.count.aa$a*n3.count.aa$h))/n3.count.diaa$ah

cps.gct.cca<-((n3.trinuc$gct*n3.trinuc$cca)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gct.ccc<-((n3.trinuc$gct*n3.trinuc$ccc)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gct.ccg<-((n3.trinuc$gct*n3.trinuc$ccg)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap
cps.gct.cct<-((n3.trinuc$gct*n3.trinuc$cct)/(n3.count.aa$a*n3.count.aa$p))/n3.count.diaa$ap

cps.gct.cga<-((n3.trinuc$gct*n3.trinuc$cga)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gct.cgc<-((n3.trinuc$gct*n3.trinuc$cgc)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gct.cgg<-((n3.trinuc$gct*n3.trinuc$cgg)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar
cps.gct.cgt<-((n3.trinuc$gct*n3.trinuc$cgt)/(n3.count.aa$a*n3.count.aa$r))/n3.count.diaa$ar

cps.gct.cta<-((n3.trinuc$gct*n3.trinuc$cta)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gct.ctc<-((n3.trinuc$gct*n3.trinuc$ctc)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gct.ctg<-((n3.trinuc$gct*n3.trinuc$ctg)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gct.ctt<-((n3.trinuc$gct*n3.trinuc$ctt)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al

cps.gct.gaa<-((n3.trinuc$gct*n3.trinuc$gaa)/(n3.count.aa$a*n3.count.aa$e))/n3.count.diaa$ae
cps.gct.gac<-((n3.trinuc$gct*n3.trinuc$gac)/(n3.count.aa$a*n3.count.aa$d))/n3.count.diaa$ad
cps.gct.gag<-((n3.trinuc$gct*n3.trinuc$gag)/(n3.count.aa$a*n3.count.aa$e))/n3.count.diaa$ae
cps.gct.gat<-((n3.trinuc$gct*n3.trinuc$gat)/(n3.count.aa$a*n3.count.aa$d))/n3.count.diaa$ad

cps.gct.gca<-((n3.trinuc$gct*n3.trinuc$gca)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gct.gcc<-((n3.trinuc$gct*n3.trinuc$gcc)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gct.gcg<-((n3.trinuc$gct*n3.trinuc$gcg)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa
cps.gct.gct<-((n3.trinuc$gct*n3.trinuc$gct)/(n3.count.aa$a*n3.count.aa$a))/n3.count.diaa$aa

cps.gct.gga<-((n3.trinuc$gct*n3.trinuc$gga)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gct.ggc<-((n3.trinuc$gct*n3.trinuc$ggc)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gct.ggg<-((n3.trinuc$gct*n3.trinuc$ggg)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag
cps.gct.ggt<-((n3.trinuc$gct*n3.trinuc$ggt)/(n3.count.aa$a*n3.count.aa$g))/n3.count.diaa$ag

cps.gct.gta<-((n3.trinuc$gct*n3.trinuc$gta)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gct.gtc<-((n3.trinuc$gct*n3.trinuc$gtc)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gct.gtg<-((n3.trinuc$gct*n3.trinuc$gtg)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av
cps.gct.gtt<-((n3.trinuc$gct*n3.trinuc$gtt)/(n3.count.aa$a*n3.count.aa$v))/n3.count.diaa$av

#Stop codon
#cps.gct.taa<-((n3.trinuc$gct*n3.trinuc$taa)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gct.tac<-((n3.trinuc$gct*n3.trinuc$tac)/(n3.count.aa$a*n3.count.aa$y))/n3.count.diaa$ay
#Stop codon
#cps.gct.tag<-((n3.trinuc$gct*n3.trinuc$tag)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gct.tat<-((n3.trinuc$gct*n3.trinuc$tat)/(n3.count.aa$a*n3.count.aa$y))/n3.count.diaa$ay

cps.gct.tca<-((n3.trinuc$gct*n3.trinuc$tca)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gct.tcc<-((n3.trinuc$gct*n3.trinuc$tcc)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gct.tcg<-((n3.trinuc$gct*n3.trinuc$tcg)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as
cps.gct.tct<-((n3.trinuc$gct*n3.trinuc$tct)/(n3.count.aa$a*n3.count.aa$s))/n3.count.diaa$as

#Stop codon
#cps.gct.tga<-((n3.trinuc$gct*n3.trinuc$tga)/(n3.count.aa$a*n3.count.aa$k))/n3.count.diaa$ak
cps.gct.tgc<-((n3.trinuc$gct*n3.trinuc$tgc)/(n3.count.aa$a*n3.count.aa$c))/n3.count.diaa$ac
cps.gct.tgg<-((n3.trinuc$gct*n3.trinuc$tgg)/(n3.count.aa$a*n3.count.aa$w))/n3.count.diaa$aw
cps.gct.tgt<-((n3.trinuc$gct*n3.trinuc$tgt)/(n3.count.aa$a*n3.count.aa$c))/n3.count.diaa$ac

cps.gct.tta<-((n3.trinuc$gct*n3.trinuc$tta)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gct.ttc<-((n3.trinuc$gct*n3.trinuc$ttc)/(n3.count.aa$a*n3.count.aa$f))/n3.count.diaa$af
cps.gct.ttg<-((n3.trinuc$gct*n3.trinuc$ttg)/(n3.count.aa$a*n3.count.aa$l))/n3.count.diaa$al
cps.gct.ttt<-((n3.trinuc$gct*n3.trinuc$ttt)/(n3.count.aa$a*n3.count.aa$f))/n3.count.diaa$af








cps.gga.aaa<-((n3.trinuc$gga*n3.trinuc$aaa)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.gga.aac<-((n3.trinuc$gga*n3.trinuc$aac)/(n3.count.aa$g*n3.count.aa$n))/n3.count.diaa$gn
cps.gga.aag<-((n3.trinuc$gga*n3.trinuc$aag)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.gga.aat<-((n3.trinuc$gga*n3.trinuc$aat)/(n3.count.aa$g*n3.count.aa$n))/n3.count.diaa$gn

cps.gga.aca<-((n3.trinuc$gga*n3.trinuc$aca)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.gga.acc<-((n3.trinuc$gga*n3.trinuc$acc)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.gga.acg<-((n3.trinuc$gga*n3.trinuc$acg)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.gga.act<-((n3.trinuc$gga*n3.trinuc$act)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt

cps.gga.aga<-((n3.trinuc$gga*n3.trinuc$aga)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.gga.agc<-((n3.trinuc$gga*n3.trinuc$agc)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.gga.agg<-((n3.trinuc$gga*n3.trinuc$agg)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.gga.agt<-((n3.trinuc$gga*n3.trinuc$agt)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs

cps.gga.ata<-((n3.trinuc$gga*n3.trinuc$ata)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi
cps.gga.atc<-((n3.trinuc$gga*n3.trinuc$atc)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi
cps.gga.atg<-((n3.trinuc$gga*n3.trinuc$atg)/(n3.count.aa$g*n3.count.aa$m))/n3.count.diaa$gm
cps.gga.att<-((n3.trinuc$gga*n3.trinuc$att)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi

cps.gga.caa<-((n3.trinuc$gga*n3.trinuc$caa)/(n3.count.aa$g*n3.count.aa$q))/n3.count.diaa$gq
cps.gga.cac<-((n3.trinuc$gga*n3.trinuc$cac)/(n3.count.aa$g*n3.count.aa$h))/n3.count.diaa$gh
cps.gga.cag<-((n3.trinuc$gga*n3.trinuc$cag)/(n3.count.aa$g*n3.count.aa$q))/n3.count.diaa$gq
cps.gga.cat<-((n3.trinuc$gga*n3.trinuc$cat)/(n3.count.aa$g*n3.count.aa$h))/n3.count.diaa$gh

cps.gga.cca<-((n3.trinuc$gga*n3.trinuc$cca)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.gga.ccc<-((n3.trinuc$gga*n3.trinuc$ccc)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.gga.ccg<-((n3.trinuc$gga*n3.trinuc$ccg)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.gga.cct<-((n3.trinuc$gga*n3.trinuc$cct)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp

cps.gga.cga<-((n3.trinuc$gga*n3.trinuc$cga)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.gga.cgc<-((n3.trinuc$gga*n3.trinuc$cgc)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.gga.cgg<-((n3.trinuc$gga*n3.trinuc$cgg)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.gga.cgt<-((n3.trinuc$gga*n3.trinuc$cgt)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr

cps.gga.cta<-((n3.trinuc$gga*n3.trinuc$cta)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.gga.ctc<-((n3.trinuc$gga*n3.trinuc$ctc)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.gga.ctg<-((n3.trinuc$gga*n3.trinuc$ctg)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.gga.ctt<-((n3.trinuc$gga*n3.trinuc$ctt)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl

cps.gga.gaa<-((n3.trinuc$gga*n3.trinuc$gaa)/(n3.count.aa$g*n3.count.aa$e))/n3.count.diaa$ge
cps.gga.gac<-((n3.trinuc$gga*n3.trinuc$gac)/(n3.count.aa$g*n3.count.aa$d))/n3.count.diaa$gd
cps.gga.gag<-((n3.trinuc$gga*n3.trinuc$gag)/(n3.count.aa$g*n3.count.aa$e))/n3.count.diaa$ge
cps.gga.gat<-((n3.trinuc$gga*n3.trinuc$gat)/(n3.count.aa$g*n3.count.aa$d))/n3.count.diaa$gd

cps.gga.gca<-((n3.trinuc$gga*n3.trinuc$gca)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.gga.gcc<-((n3.trinuc$gga*n3.trinuc$gcc)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.gga.gcg<-((n3.trinuc$gga*n3.trinuc$gcg)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.gga.gct<-((n3.trinuc$gga*n3.trinuc$gct)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga

cps.gga.gga<-((n3.trinuc$gga*n3.trinuc$gga)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.gga.ggc<-((n3.trinuc$gga*n3.trinuc$ggc)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.gga.ggg<-((n3.trinuc$gga*n3.trinuc$ggg)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.gga.ggt<-((n3.trinuc$gga*n3.trinuc$ggt)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg

cps.gga.gta<-((n3.trinuc$gga*n3.trinuc$gta)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.gga.gtc<-((n3.trinuc$gga*n3.trinuc$gtc)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.gga.gtg<-((n3.trinuc$gga*n3.trinuc$gtg)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.gga.gtt<-((n3.trinuc$gga*n3.trinuc$gtt)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv

#Stop codon
#cps.gga.taa<-((n3.trinuc$gga*n3.trinuc$taa)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.gga.tac<-((n3.trinuc$gga*n3.trinuc$tac)/(n3.count.aa$g*n3.count.aa$y))/n3.count.diaa$gy
#Stop codon
#cps.gga.tag<-((n3.trinuc$gga*n3.trinuc$tag)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.gga.tat<-((n3.trinuc$gga*n3.trinuc$tat)/(n3.count.aa$g*n3.count.aa$y))/n3.count.diaa$gy

cps.gga.tca<-((n3.trinuc$gga*n3.trinuc$tca)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.gga.tcc<-((n3.trinuc$gga*n3.trinuc$tcc)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.gga.tcg<-((n3.trinuc$gga*n3.trinuc$tcg)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.gga.tct<-((n3.trinuc$gga*n3.trinuc$tct)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs

#Stop codon
#cps.gga.tga<-((n3.trinuc$gga*n3.trinuc$tga)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.gga.tgc<-((n3.trinuc$gga*n3.trinuc$tgc)/(n3.count.aa$g*n3.count.aa$c))/n3.count.diaa$gc
cps.gga.tgg<-((n3.trinuc$gga*n3.trinuc$tgg)/(n3.count.aa$g*n3.count.aa$w))/n3.count.diaa$gw
cps.gga.tgt<-((n3.trinuc$gga*n3.trinuc$tgt)/(n3.count.aa$g*n3.count.aa$c))/n3.count.diaa$gc

cps.gga.tta<-((n3.trinuc$gga*n3.trinuc$tta)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.gga.ttc<-((n3.trinuc$gga*n3.trinuc$ttc)/(n3.count.aa$g*n3.count.aa$f))/n3.count.diaa$gf
cps.gga.ttg<-((n3.trinuc$gga*n3.trinuc$ttg)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.gga.ttt<-((n3.trinuc$gga*n3.trinuc$ttt)/(n3.count.aa$g*n3.count.aa$f))/n3.count.diaa$gf








cps.ggc.aaa<-((n3.trinuc$ggc*n3.trinuc$aaa)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggc.aac<-((n3.trinuc$ggc*n3.trinuc$aac)/(n3.count.aa$g*n3.count.aa$n))/n3.count.diaa$gn
cps.ggc.aag<-((n3.trinuc$ggc*n3.trinuc$aag)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggc.aat<-((n3.trinuc$ggc*n3.trinuc$aat)/(n3.count.aa$g*n3.count.aa$n))/n3.count.diaa$gn

cps.ggc.aca<-((n3.trinuc$ggc*n3.trinuc$aca)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggc.acc<-((n3.trinuc$ggc*n3.trinuc$acc)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggc.acg<-((n3.trinuc$ggc*n3.trinuc$acg)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggc.act<-((n3.trinuc$ggc*n3.trinuc$act)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt

cps.ggc.aga<-((n3.trinuc$ggc*n3.trinuc$aga)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggc.agc<-((n3.trinuc$ggc*n3.trinuc$agc)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggc.agg<-((n3.trinuc$ggc*n3.trinuc$agg)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggc.agt<-((n3.trinuc$ggc*n3.trinuc$agt)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs

cps.ggc.ata<-((n3.trinuc$ggc*n3.trinuc$ata)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi
cps.ggc.atc<-((n3.trinuc$ggc*n3.trinuc$atc)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi
cps.ggc.atg<-((n3.trinuc$ggc*n3.trinuc$atg)/(n3.count.aa$g*n3.count.aa$m))/n3.count.diaa$gm
cps.ggc.att<-((n3.trinuc$ggc*n3.trinuc$att)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi

cps.ggc.caa<-((n3.trinuc$ggc*n3.trinuc$caa)/(n3.count.aa$g*n3.count.aa$q))/n3.count.diaa$gq
cps.ggc.cac<-((n3.trinuc$ggc*n3.trinuc$cac)/(n3.count.aa$g*n3.count.aa$h))/n3.count.diaa$gh
cps.ggc.cag<-((n3.trinuc$ggc*n3.trinuc$cag)/(n3.count.aa$g*n3.count.aa$q))/n3.count.diaa$gq
cps.ggc.cat<-((n3.trinuc$ggc*n3.trinuc$cat)/(n3.count.aa$g*n3.count.aa$h))/n3.count.diaa$gh

cps.ggc.cca<-((n3.trinuc$ggc*n3.trinuc$cca)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggc.ccc<-((n3.trinuc$ggc*n3.trinuc$ccc)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggc.ccg<-((n3.trinuc$ggc*n3.trinuc$ccg)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggc.cct<-((n3.trinuc$ggc*n3.trinuc$cct)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp

cps.ggc.cga<-((n3.trinuc$ggc*n3.trinuc$cga)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggc.cgc<-((n3.trinuc$ggc*n3.trinuc$cgc)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggc.cgg<-((n3.trinuc$ggc*n3.trinuc$cgg)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggc.cgt<-((n3.trinuc$ggc*n3.trinuc$cgt)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr

cps.ggc.cta<-((n3.trinuc$ggc*n3.trinuc$cta)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggc.ctc<-((n3.trinuc$ggc*n3.trinuc$ctc)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggc.ctg<-((n3.trinuc$ggc*n3.trinuc$ctg)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggc.ctt<-((n3.trinuc$ggc*n3.trinuc$ctt)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl

cps.ggc.gaa<-((n3.trinuc$ggc*n3.trinuc$gaa)/(n3.count.aa$g*n3.count.aa$e))/n3.count.diaa$ge
cps.ggc.gac<-((n3.trinuc$ggc*n3.trinuc$gac)/(n3.count.aa$g*n3.count.aa$d))/n3.count.diaa$gd
cps.ggc.gag<-((n3.trinuc$ggc*n3.trinuc$gag)/(n3.count.aa$g*n3.count.aa$e))/n3.count.diaa$ge
cps.ggc.gat<-((n3.trinuc$ggc*n3.trinuc$gat)/(n3.count.aa$g*n3.count.aa$d))/n3.count.diaa$gd

cps.ggc.gca<-((n3.trinuc$ggc*n3.trinuc$gca)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggc.gcc<-((n3.trinuc$ggc*n3.trinuc$gcc)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggc.gcg<-((n3.trinuc$ggc*n3.trinuc$gcg)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggc.gct<-((n3.trinuc$ggc*n3.trinuc$gct)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga

cps.ggc.gga<-((n3.trinuc$ggc*n3.trinuc$gga)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggc.ggc<-((n3.trinuc$ggc*n3.trinuc$ggc)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggc.ggg<-((n3.trinuc$ggc*n3.trinuc$ggg)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggc.ggt<-((n3.trinuc$ggc*n3.trinuc$ggt)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg

cps.ggc.gta<-((n3.trinuc$ggc*n3.trinuc$gta)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggc.gtc<-((n3.trinuc$ggc*n3.trinuc$gtc)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggc.gtg<-((n3.trinuc$ggc*n3.trinuc$gtg)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggc.gtt<-((n3.trinuc$ggc*n3.trinuc$gtt)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv

#Stop codon
#cps.ggc.taa<-((n3.trinuc$ggc*n3.trinuc$taa)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggc.tac<-((n3.trinuc$ggc*n3.trinuc$tac)/(n3.count.aa$g*n3.count.aa$y))/n3.count.diaa$gy
#Stop codon
#cps.ggc.tag<-((n3.trinuc$ggc*n3.trinuc$tag)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggc.tat<-((n3.trinuc$ggc*n3.trinuc$tat)/(n3.count.aa$g*n3.count.aa$y))/n3.count.diaa$gy

cps.ggc.tca<-((n3.trinuc$ggc*n3.trinuc$tca)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggc.tcc<-((n3.trinuc$ggc*n3.trinuc$tcc)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggc.tcg<-((n3.trinuc$ggc*n3.trinuc$tcg)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggc.tct<-((n3.trinuc$ggc*n3.trinuc$tct)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs

#Stop codon
#cps.ggc.tga<-((n3.trinuc$ggc*n3.trinuc$tga)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggc.tgc<-((n3.trinuc$ggc*n3.trinuc$tgc)/(n3.count.aa$g*n3.count.aa$c))/n3.count.diaa$gc
cps.ggc.tgg<-((n3.trinuc$ggc*n3.trinuc$tgg)/(n3.count.aa$g*n3.count.aa$w))/n3.count.diaa$gw
cps.ggc.tgt<-((n3.trinuc$ggc*n3.trinuc$tgt)/(n3.count.aa$g*n3.count.aa$c))/n3.count.diaa$gc

cps.ggc.tta<-((n3.trinuc$ggc*n3.trinuc$tta)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggc.ttc<-((n3.trinuc$ggc*n3.trinuc$ttc)/(n3.count.aa$g*n3.count.aa$f))/n3.count.diaa$gf
cps.ggc.ttg<-((n3.trinuc$ggc*n3.trinuc$ttg)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggc.ttt<-((n3.trinuc$ggc*n3.trinuc$ttt)/(n3.count.aa$g*n3.count.aa$f))/n3.count.diaa$gf








cps.ggg.aaa<-((n3.trinuc$ggg*n3.trinuc$aaa)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggg.aac<-((n3.trinuc$ggg*n3.trinuc$aac)/(n3.count.aa$g*n3.count.aa$n))/n3.count.diaa$gn
cps.ggg.aag<-((n3.trinuc$ggg*n3.trinuc$aag)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggg.aat<-((n3.trinuc$ggg*n3.trinuc$aat)/(n3.count.aa$g*n3.count.aa$n))/n3.count.diaa$gn

cps.ggg.aca<-((n3.trinuc$ggg*n3.trinuc$aca)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggg.acc<-((n3.trinuc$ggg*n3.trinuc$acc)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggg.acg<-((n3.trinuc$ggg*n3.trinuc$acg)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggg.act<-((n3.trinuc$ggg*n3.trinuc$act)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt

cps.ggg.aga<-((n3.trinuc$ggg*n3.trinuc$aga)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggg.agc<-((n3.trinuc$ggg*n3.trinuc$agc)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggg.agg<-((n3.trinuc$ggg*n3.trinuc$agg)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggg.agt<-((n3.trinuc$ggg*n3.trinuc$agt)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs

cps.ggg.ata<-((n3.trinuc$ggg*n3.trinuc$ata)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi
cps.ggg.atc<-((n3.trinuc$ggg*n3.trinuc$atc)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi
cps.ggg.atg<-((n3.trinuc$ggg*n3.trinuc$atg)/(n3.count.aa$g*n3.count.aa$m))/n3.count.diaa$gm
cps.ggg.att<-((n3.trinuc$ggg*n3.trinuc$att)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi

cps.ggg.caa<-((n3.trinuc$ggg*n3.trinuc$caa)/(n3.count.aa$g*n3.count.aa$q))/n3.count.diaa$gq
cps.ggg.cac<-((n3.trinuc$ggg*n3.trinuc$cac)/(n3.count.aa$g*n3.count.aa$h))/n3.count.diaa$gh
cps.ggg.cag<-((n3.trinuc$ggg*n3.trinuc$cag)/(n3.count.aa$g*n3.count.aa$q))/n3.count.diaa$gq
cps.ggg.cat<-((n3.trinuc$ggg*n3.trinuc$cat)/(n3.count.aa$g*n3.count.aa$h))/n3.count.diaa$gh

cps.ggg.cca<-((n3.trinuc$ggg*n3.trinuc$cca)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggg.ccc<-((n3.trinuc$ggg*n3.trinuc$ccc)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggg.ccg<-((n3.trinuc$ggg*n3.trinuc$ccg)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggg.cct<-((n3.trinuc$ggg*n3.trinuc$cct)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp

cps.ggg.cga<-((n3.trinuc$ggg*n3.trinuc$cga)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggg.cgc<-((n3.trinuc$ggg*n3.trinuc$cgc)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggg.cgg<-((n3.trinuc$ggg*n3.trinuc$cgg)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggg.cgt<-((n3.trinuc$ggg*n3.trinuc$cgt)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr

cps.ggg.cta<-((n3.trinuc$ggg*n3.trinuc$cta)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggg.ctc<-((n3.trinuc$ggg*n3.trinuc$ctc)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggg.ctg<-((n3.trinuc$ggg*n3.trinuc$ctg)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggg.ctt<-((n3.trinuc$ggg*n3.trinuc$ctt)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl

cps.ggg.gaa<-((n3.trinuc$ggg*n3.trinuc$gaa)/(n3.count.aa$g*n3.count.aa$e))/n3.count.diaa$ge
cps.ggg.gac<-((n3.trinuc$ggg*n3.trinuc$gac)/(n3.count.aa$g*n3.count.aa$d))/n3.count.diaa$gd
cps.ggg.gag<-((n3.trinuc$ggg*n3.trinuc$gag)/(n3.count.aa$g*n3.count.aa$e))/n3.count.diaa$ge
cps.ggg.gat<-((n3.trinuc$ggg*n3.trinuc$gat)/(n3.count.aa$g*n3.count.aa$d))/n3.count.diaa$gd

cps.ggg.gca<-((n3.trinuc$ggg*n3.trinuc$gca)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggg.gcc<-((n3.trinuc$ggg*n3.trinuc$gcc)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggg.gcg<-((n3.trinuc$ggg*n3.trinuc$gcg)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggg.gct<-((n3.trinuc$ggg*n3.trinuc$gct)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga

cps.ggg.gga<-((n3.trinuc$ggg*n3.trinuc$gga)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggg.ggc<-((n3.trinuc$ggg*n3.trinuc$ggc)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggg.ggg<-((n3.trinuc$ggg*n3.trinuc$ggg)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggg.ggt<-((n3.trinuc$ggg*n3.trinuc$ggt)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg

cps.ggg.gta<-((n3.trinuc$ggg*n3.trinuc$gta)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggg.gtc<-((n3.trinuc$ggg*n3.trinuc$gtc)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggg.gtg<-((n3.trinuc$ggg*n3.trinuc$gtg)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggg.gtt<-((n3.trinuc$ggg*n3.trinuc$gtt)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv

#Stop codon
#cps.ggg.taa<-((n3.trinuc$ggg*n3.trinuc$taa)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggg.tac<-((n3.trinuc$ggg*n3.trinuc$tac)/(n3.count.aa$g*n3.count.aa$y))/n3.count.diaa$gy
#Stop codon
#cps.ggg.tag<-((n3.trinuc$ggg*n3.trinuc$tag)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggg.tat<-((n3.trinuc$ggg*n3.trinuc$tat)/(n3.count.aa$g*n3.count.aa$y))/n3.count.diaa$gy

cps.ggg.tca<-((n3.trinuc$ggg*n3.trinuc$tca)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggg.tcc<-((n3.trinuc$ggg*n3.trinuc$tcc)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggg.tcg<-((n3.trinuc$ggg*n3.trinuc$tcg)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggg.tct<-((n3.trinuc$ggg*n3.trinuc$tct)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs

#Stop codon
#cps.ggg.tga<-((n3.trinuc$ggg*n3.trinuc$tga)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggg.tgc<-((n3.trinuc$ggg*n3.trinuc$tgc)/(n3.count.aa$g*n3.count.aa$c))/n3.count.diaa$gc
cps.ggg.tgg<-((n3.trinuc$ggg*n3.trinuc$tgg)/(n3.count.aa$g*n3.count.aa$w))/n3.count.diaa$gw
cps.ggg.tgt<-((n3.trinuc$ggg*n3.trinuc$tgt)/(n3.count.aa$g*n3.count.aa$c))/n3.count.diaa$gc

cps.ggg.tta<-((n3.trinuc$ggg*n3.trinuc$tta)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggg.ttc<-((n3.trinuc$ggg*n3.trinuc$ttc)/(n3.count.aa$g*n3.count.aa$f))/n3.count.diaa$gf
cps.ggg.ttg<-((n3.trinuc$ggg*n3.trinuc$ttg)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggg.ttt<-((n3.trinuc$ggg*n3.trinuc$ttt)/(n3.count.aa$g*n3.count.aa$f))/n3.count.diaa$gf








cps.ggt.aaa<-((n3.trinuc$ggt*n3.trinuc$aaa)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggt.aac<-((n3.trinuc$ggt*n3.trinuc$aac)/(n3.count.aa$g*n3.count.aa$n))/n3.count.diaa$gn
cps.ggt.aag<-((n3.trinuc$ggt*n3.trinuc$aag)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggt.aat<-((n3.trinuc$ggt*n3.trinuc$aat)/(n3.count.aa$g*n3.count.aa$n))/n3.count.diaa$gn

cps.ggt.aca<-((n3.trinuc$ggt*n3.trinuc$aca)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggt.acc<-((n3.trinuc$ggt*n3.trinuc$acc)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggt.acg<-((n3.trinuc$ggt*n3.trinuc$acg)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt
cps.ggt.act<-((n3.trinuc$ggt*n3.trinuc$act)/(n3.count.aa$g*n3.count.aa$t))/n3.count.diaa$gt

cps.ggt.aga<-((n3.trinuc$ggt*n3.trinuc$aga)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggt.agc<-((n3.trinuc$ggt*n3.trinuc$agc)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggt.agg<-((n3.trinuc$ggt*n3.trinuc$agg)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggt.agt<-((n3.trinuc$ggt*n3.trinuc$agt)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs

cps.ggt.ata<-((n3.trinuc$ggt*n3.trinuc$ata)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi
cps.ggt.atc<-((n3.trinuc$ggt*n3.trinuc$atc)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi
cps.ggt.atg<-((n3.trinuc$ggt*n3.trinuc$atg)/(n3.count.aa$g*n3.count.aa$m))/n3.count.diaa$gm
cps.ggt.att<-((n3.trinuc$ggt*n3.trinuc$att)/(n3.count.aa$g*n3.count.aa$i))/n3.count.diaa$gi

cps.ggt.caa<-((n3.trinuc$ggt*n3.trinuc$caa)/(n3.count.aa$g*n3.count.aa$q))/n3.count.diaa$gq
cps.ggt.cac<-((n3.trinuc$ggt*n3.trinuc$cac)/(n3.count.aa$g*n3.count.aa$h))/n3.count.diaa$gh
cps.ggt.cag<-((n3.trinuc$ggt*n3.trinuc$cag)/(n3.count.aa$g*n3.count.aa$q))/n3.count.diaa$gq
cps.ggt.cat<-((n3.trinuc$ggt*n3.trinuc$cat)/(n3.count.aa$g*n3.count.aa$h))/n3.count.diaa$gh

cps.ggt.cca<-((n3.trinuc$ggt*n3.trinuc$cca)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggt.ccc<-((n3.trinuc$ggt*n3.trinuc$ccc)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggt.ccg<-((n3.trinuc$ggt*n3.trinuc$ccg)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp
cps.ggt.cct<-((n3.trinuc$ggt*n3.trinuc$cct)/(n3.count.aa$g*n3.count.aa$p))/n3.count.diaa$gp

cps.ggt.cga<-((n3.trinuc$ggt*n3.trinuc$cga)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggt.cgc<-((n3.trinuc$ggt*n3.trinuc$cgc)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggt.cgg<-((n3.trinuc$ggt*n3.trinuc$cgg)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr
cps.ggt.cgt<-((n3.trinuc$ggt*n3.trinuc$cgt)/(n3.count.aa$g*n3.count.aa$r))/n3.count.diaa$gr

cps.ggt.cta<-((n3.trinuc$ggt*n3.trinuc$cta)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggt.ctc<-((n3.trinuc$ggt*n3.trinuc$ctc)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggt.ctg<-((n3.trinuc$ggt*n3.trinuc$ctg)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggt.ctt<-((n3.trinuc$ggt*n3.trinuc$ctt)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl

cps.ggt.gaa<-((n3.trinuc$ggt*n3.trinuc$gaa)/(n3.count.aa$g*n3.count.aa$e))/n3.count.diaa$ge
cps.ggt.gac<-((n3.trinuc$ggt*n3.trinuc$gac)/(n3.count.aa$g*n3.count.aa$d))/n3.count.diaa$gd
cps.ggt.gag<-((n3.trinuc$ggt*n3.trinuc$gag)/(n3.count.aa$g*n3.count.aa$e))/n3.count.diaa$ge
cps.ggt.gat<-((n3.trinuc$ggt*n3.trinuc$gat)/(n3.count.aa$g*n3.count.aa$d))/n3.count.diaa$gd

cps.ggt.gca<-((n3.trinuc$ggt*n3.trinuc$gca)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggt.gcc<-((n3.trinuc$ggt*n3.trinuc$gcc)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggt.gcg<-((n3.trinuc$ggt*n3.trinuc$gcg)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga
cps.ggt.gct<-((n3.trinuc$ggt*n3.trinuc$gct)/(n3.count.aa$g*n3.count.aa$a))/n3.count.diaa$ga

cps.ggt.gga<-((n3.trinuc$ggt*n3.trinuc$gga)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggt.ggc<-((n3.trinuc$ggt*n3.trinuc$ggc)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggt.ggg<-((n3.trinuc$ggt*n3.trinuc$ggg)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg
cps.ggt.ggt<-((n3.trinuc$ggt*n3.trinuc$ggt)/(n3.count.aa$g*n3.count.aa$g))/n3.count.diaa$gg

cps.ggt.gta<-((n3.trinuc$ggt*n3.trinuc$gta)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggt.gtc<-((n3.trinuc$ggt*n3.trinuc$gtc)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggt.gtg<-((n3.trinuc$ggt*n3.trinuc$gtg)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv
cps.ggt.gtt<-((n3.trinuc$ggt*n3.trinuc$gtt)/(n3.count.aa$g*n3.count.aa$v))/n3.count.diaa$gv

#Stop codon
#cps.ggt.taa<-((n3.trinuc$ggt*n3.trinuc$taa)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggt.tac<-((n3.trinuc$ggt*n3.trinuc$tac)/(n3.count.aa$g*n3.count.aa$y))/n3.count.diaa$gy
#Stop codon
#cps.ggt.tag<-((n3.trinuc$ggt*n3.trinuc$tag)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggt.tat<-((n3.trinuc$ggt*n3.trinuc$tat)/(n3.count.aa$g*n3.count.aa$y))/n3.count.diaa$gy

cps.ggt.tca<-((n3.trinuc$ggt*n3.trinuc$tca)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggt.tcc<-((n3.trinuc$ggt*n3.trinuc$tcc)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggt.tcg<-((n3.trinuc$ggt*n3.trinuc$tcg)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs
cps.ggt.tct<-((n3.trinuc$ggt*n3.trinuc$tct)/(n3.count.aa$g*n3.count.aa$s))/n3.count.diaa$gs

#Stop codon
#cps.ggt.tga<-((n3.trinuc$ggt*n3.trinuc$tga)/(n3.count.aa$g*n3.count.aa$k))/n3.count.diaa$gk
cps.ggt.tgc<-((n3.trinuc$ggt*n3.trinuc$tgc)/(n3.count.aa$g*n3.count.aa$c))/n3.count.diaa$gc
cps.ggt.tgg<-((n3.trinuc$ggt*n3.trinuc$tgg)/(n3.count.aa$g*n3.count.aa$w))/n3.count.diaa$gw
cps.ggt.tgt<-((n3.trinuc$ggt*n3.trinuc$tgt)/(n3.count.aa$g*n3.count.aa$c))/n3.count.diaa$gc

cps.ggt.tta<-((n3.trinuc$ggt*n3.trinuc$tta)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggt.ttc<-((n3.trinuc$ggt*n3.trinuc$ttc)/(n3.count.aa$g*n3.count.aa$f))/n3.count.diaa$gf
cps.ggt.ttg<-((n3.trinuc$ggt*n3.trinuc$ttg)/(n3.count.aa$g*n3.count.aa$l))/n3.count.diaa$gl
cps.ggt.ttt<-((n3.trinuc$ggt*n3.trinuc$ttt)/(n3.count.aa$g*n3.count.aa$f))/n3.count.diaa$gf








cps.gta.aaa<-((n3.trinuc$gta*n3.trinuc$aaa)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gta.aac<-((n3.trinuc$gta*n3.trinuc$aac)/(n3.count.aa$v*n3.count.aa$n))/n3.count.diaa$vn
cps.gta.aag<-((n3.trinuc$gta*n3.trinuc$aag)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gta.aat<-((n3.trinuc$gta*n3.trinuc$aat)/(n3.count.aa$v*n3.count.aa$n))/n3.count.diaa$vn

cps.gta.aca<-((n3.trinuc$gta*n3.trinuc$aca)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gta.acc<-((n3.trinuc$gta*n3.trinuc$acc)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gta.acg<-((n3.trinuc$gta*n3.trinuc$acg)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gta.act<-((n3.trinuc$gta*n3.trinuc$act)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt

cps.gta.aga<-((n3.trinuc$gta*n3.trinuc$aga)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gta.agc<-((n3.trinuc$gta*n3.trinuc$agc)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gta.agg<-((n3.trinuc$gta*n3.trinuc$agg)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gta.agt<-((n3.trinuc$gta*n3.trinuc$agt)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs

cps.gta.ata<-((n3.trinuc$gta*n3.trinuc$ata)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi
cps.gta.atc<-((n3.trinuc$gta*n3.trinuc$atc)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi
cps.gta.atg<-((n3.trinuc$gta*n3.trinuc$atg)/(n3.count.aa$v*n3.count.aa$m))/n3.count.diaa$vm
cps.gta.att<-((n3.trinuc$gta*n3.trinuc$att)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi

cps.gta.caa<-((n3.trinuc$gta*n3.trinuc$caa)/(n3.count.aa$v*n3.count.aa$q))/n3.count.diaa$vq
cps.gta.cac<-((n3.trinuc$gta*n3.trinuc$cac)/(n3.count.aa$v*n3.count.aa$h))/n3.count.diaa$vh
cps.gta.cag<-((n3.trinuc$gta*n3.trinuc$cag)/(n3.count.aa$v*n3.count.aa$q))/n3.count.diaa$vq
cps.gta.cat<-((n3.trinuc$gta*n3.trinuc$cat)/(n3.count.aa$v*n3.count.aa$h))/n3.count.diaa$vh

cps.gta.cca<-((n3.trinuc$gta*n3.trinuc$cca)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gta.ccc<-((n3.trinuc$gta*n3.trinuc$ccc)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gta.ccg<-((n3.trinuc$gta*n3.trinuc$ccg)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gta.cct<-((n3.trinuc$gta*n3.trinuc$cct)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp

cps.gta.cga<-((n3.trinuc$gta*n3.trinuc$cga)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gta.cgc<-((n3.trinuc$gta*n3.trinuc$cgc)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gta.cgg<-((n3.trinuc$gta*n3.trinuc$cgg)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gta.cgt<-((n3.trinuc$gta*n3.trinuc$cgt)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr

cps.gta.cta<-((n3.trinuc$gta*n3.trinuc$cta)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gta.ctc<-((n3.trinuc$gta*n3.trinuc$ctc)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gta.ctg<-((n3.trinuc$gta*n3.trinuc$ctg)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gta.ctt<-((n3.trinuc$gta*n3.trinuc$ctt)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl

cps.gta.gaa<-((n3.trinuc$gta*n3.trinuc$gaa)/(n3.count.aa$v*n3.count.aa$e))/n3.count.diaa$ve
cps.gta.gac<-((n3.trinuc$gta*n3.trinuc$gac)/(n3.count.aa$v*n3.count.aa$d))/n3.count.diaa$vd
cps.gta.gag<-((n3.trinuc$gta*n3.trinuc$gag)/(n3.count.aa$v*n3.count.aa$e))/n3.count.diaa$ve
cps.gta.gat<-((n3.trinuc$gta*n3.trinuc$gat)/(n3.count.aa$v*n3.count.aa$d))/n3.count.diaa$vd

cps.gta.gca<-((n3.trinuc$gta*n3.trinuc$gca)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gta.gcc<-((n3.trinuc$gta*n3.trinuc$gcc)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gta.gcg<-((n3.trinuc$gta*n3.trinuc$gcg)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gta.gct<-((n3.trinuc$gta*n3.trinuc$gct)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va

cps.gta.gga<-((n3.trinuc$gta*n3.trinuc$gga)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gta.ggc<-((n3.trinuc$gta*n3.trinuc$ggc)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gta.ggg<-((n3.trinuc$gta*n3.trinuc$ggg)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gta.ggt<-((n3.trinuc$gta*n3.trinuc$ggt)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg

cps.gta.gta<-((n3.trinuc$gta*n3.trinuc$gta)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gta.gtc<-((n3.trinuc$gta*n3.trinuc$gtc)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gta.gtg<-((n3.trinuc$gta*n3.trinuc$gtg)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gta.gtt<-((n3.trinuc$gta*n3.trinuc$gtt)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv

#Stop codon
#cps.gta.taa<-((n3.trinuc$gta*n3.trinuc$taa)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gta.tac<-((n3.trinuc$gta*n3.trinuc$tac)/(n3.count.aa$v*n3.count.aa$y))/n3.count.diaa$vy
#Stop codon
#cps.gta.tag<-((n3.trinuc$gta*n3.trinuc$tag)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gta.tat<-((n3.trinuc$gta*n3.trinuc$tat)/(n3.count.aa$v*n3.count.aa$y))/n3.count.diaa$vy

cps.gta.tca<-((n3.trinuc$gta*n3.trinuc$tca)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gta.tcc<-((n3.trinuc$gta*n3.trinuc$tcc)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gta.tcg<-((n3.trinuc$gta*n3.trinuc$tcg)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gta.tct<-((n3.trinuc$gta*n3.trinuc$tct)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs

#Stop codon
#cps.gta.tga<-((n3.trinuc$gta*n3.trinuc$tga)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gta.tgc<-((n3.trinuc$gta*n3.trinuc$tgc)/(n3.count.aa$v*n3.count.aa$c))/n3.count.diaa$vc
cps.gta.tgg<-((n3.trinuc$gta*n3.trinuc$tgg)/(n3.count.aa$v*n3.count.aa$w))/n3.count.diaa$vw
cps.gta.tgt<-((n3.trinuc$gta*n3.trinuc$tgt)/(n3.count.aa$v*n3.count.aa$c))/n3.count.diaa$vc

cps.gta.tta<-((n3.trinuc$gta*n3.trinuc$tta)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gta.ttc<-((n3.trinuc$gta*n3.trinuc$ttc)/(n3.count.aa$v*n3.count.aa$f))/n3.count.diaa$vf
cps.gta.ttg<-((n3.trinuc$gta*n3.trinuc$ttg)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gta.ttt<-((n3.trinuc$gta*n3.trinuc$ttt)/(n3.count.aa$v*n3.count.aa$f))/n3.count.diaa$vf








cps.gtc.aaa<-((n3.trinuc$gtc*n3.trinuc$aaa)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtc.aac<-((n3.trinuc$gtc*n3.trinuc$aac)/(n3.count.aa$v*n3.count.aa$n))/n3.count.diaa$vn
cps.gtc.aag<-((n3.trinuc$gtc*n3.trinuc$aag)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtc.aat<-((n3.trinuc$gtc*n3.trinuc$aat)/(n3.count.aa$v*n3.count.aa$n))/n3.count.diaa$vn

cps.gtc.aca<-((n3.trinuc$gtc*n3.trinuc$aca)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtc.acc<-((n3.trinuc$gtc*n3.trinuc$acc)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtc.acg<-((n3.trinuc$gtc*n3.trinuc$acg)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtc.act<-((n3.trinuc$gtc*n3.trinuc$act)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt

cps.gtc.aga<-((n3.trinuc$gtc*n3.trinuc$aga)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtc.agc<-((n3.trinuc$gtc*n3.trinuc$agc)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtc.agg<-((n3.trinuc$gtc*n3.trinuc$agg)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtc.agt<-((n3.trinuc$gtc*n3.trinuc$agt)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs

cps.gtc.ata<-((n3.trinuc$gtc*n3.trinuc$ata)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi
cps.gtc.atc<-((n3.trinuc$gtc*n3.trinuc$atc)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi
cps.gtc.atg<-((n3.trinuc$gtc*n3.trinuc$atg)/(n3.count.aa$v*n3.count.aa$m))/n3.count.diaa$vm
cps.gtc.att<-((n3.trinuc$gtc*n3.trinuc$att)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi

cps.gtc.caa<-((n3.trinuc$gtc*n3.trinuc$caa)/(n3.count.aa$v*n3.count.aa$q))/n3.count.diaa$vq
cps.gtc.cac<-((n3.trinuc$gtc*n3.trinuc$cac)/(n3.count.aa$v*n3.count.aa$h))/n3.count.diaa$vh
cps.gtc.cag<-((n3.trinuc$gtc*n3.trinuc$cag)/(n3.count.aa$v*n3.count.aa$q))/n3.count.diaa$vq
cps.gtc.cat<-((n3.trinuc$gtc*n3.trinuc$cat)/(n3.count.aa$v*n3.count.aa$h))/n3.count.diaa$vh

cps.gtc.cca<-((n3.trinuc$gtc*n3.trinuc$cca)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtc.ccc<-((n3.trinuc$gtc*n3.trinuc$ccc)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtc.ccg<-((n3.trinuc$gtc*n3.trinuc$ccg)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtc.cct<-((n3.trinuc$gtc*n3.trinuc$cct)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp

cps.gtc.cga<-((n3.trinuc$gtc*n3.trinuc$cga)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtc.cgc<-((n3.trinuc$gtc*n3.trinuc$cgc)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtc.cgg<-((n3.trinuc$gtc*n3.trinuc$cgg)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtc.cgt<-((n3.trinuc$gtc*n3.trinuc$cgt)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr

cps.gtc.cta<-((n3.trinuc$gtc*n3.trinuc$cta)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtc.ctc<-((n3.trinuc$gtc*n3.trinuc$ctc)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtc.ctg<-((n3.trinuc$gtc*n3.trinuc$ctg)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtc.ctt<-((n3.trinuc$gtc*n3.trinuc$ctt)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl

cps.gtc.gaa<-((n3.trinuc$gtc*n3.trinuc$gaa)/(n3.count.aa$v*n3.count.aa$e))/n3.count.diaa$ve
cps.gtc.gac<-((n3.trinuc$gtc*n3.trinuc$gac)/(n3.count.aa$v*n3.count.aa$d))/n3.count.diaa$vd
cps.gtc.gag<-((n3.trinuc$gtc*n3.trinuc$gag)/(n3.count.aa$v*n3.count.aa$e))/n3.count.diaa$ve
cps.gtc.gat<-((n3.trinuc$gtc*n3.trinuc$gat)/(n3.count.aa$v*n3.count.aa$d))/n3.count.diaa$vd

cps.gtc.gca<-((n3.trinuc$gtc*n3.trinuc$gca)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtc.gcc<-((n3.trinuc$gtc*n3.trinuc$gcc)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtc.gcg<-((n3.trinuc$gtc*n3.trinuc$gcg)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtc.gct<-((n3.trinuc$gtc*n3.trinuc$gct)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va

cps.gtc.gga<-((n3.trinuc$gtc*n3.trinuc$gga)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtc.ggc<-((n3.trinuc$gtc*n3.trinuc$ggc)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtc.ggg<-((n3.trinuc$gtc*n3.trinuc$ggg)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtc.ggt<-((n3.trinuc$gtc*n3.trinuc$ggt)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg

cps.gtc.gta<-((n3.trinuc$gtc*n3.trinuc$gta)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtc.gtc<-((n3.trinuc$gtc*n3.trinuc$gtc)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtc.gtg<-((n3.trinuc$gtc*n3.trinuc$gtg)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtc.gtt<-((n3.trinuc$gtc*n3.trinuc$gtt)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv

#Stop codon
#cps.gtc.taa<-((n3.trinuc$gtc*n3.trinuc$taa)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtc.tac<-((n3.trinuc$gtc*n3.trinuc$tac)/(n3.count.aa$v*n3.count.aa$y))/n3.count.diaa$vy
#Stop codon
#cps.gtc.tag<-((n3.trinuc$gtc*n3.trinuc$tag)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtc.tat<-((n3.trinuc$gtc*n3.trinuc$tat)/(n3.count.aa$v*n3.count.aa$y))/n3.count.diaa$vy

cps.gtc.tca<-((n3.trinuc$gtc*n3.trinuc$tca)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtc.tcc<-((n3.trinuc$gtc*n3.trinuc$tcc)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtc.tcg<-((n3.trinuc$gtc*n3.trinuc$tcg)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtc.tct<-((n3.trinuc$gtc*n3.trinuc$tct)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs

#Stop codon
#cps.gtc.tga<-((n3.trinuc$gtc*n3.trinuc$tga)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtc.tgc<-((n3.trinuc$gtc*n3.trinuc$tgc)/(n3.count.aa$v*n3.count.aa$c))/n3.count.diaa$vc
cps.gtc.tgg<-((n3.trinuc$gtc*n3.trinuc$tgg)/(n3.count.aa$v*n3.count.aa$w))/n3.count.diaa$vw
cps.gtc.tgt<-((n3.trinuc$gtc*n3.trinuc$tgt)/(n3.count.aa$v*n3.count.aa$c))/n3.count.diaa$vc

cps.gtc.tta<-((n3.trinuc$gtc*n3.trinuc$tta)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtc.ttc<-((n3.trinuc$gtc*n3.trinuc$ttc)/(n3.count.aa$v*n3.count.aa$f))/n3.count.diaa$vf
cps.gtc.ttg<-((n3.trinuc$gtc*n3.trinuc$ttg)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtc.ttt<-((n3.trinuc$gtc*n3.trinuc$ttt)/(n3.count.aa$v*n3.count.aa$f))/n3.count.diaa$vf








cps.gtg.aaa<-((n3.trinuc$gtg*n3.trinuc$aaa)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtg.aac<-((n3.trinuc$gtg*n3.trinuc$aac)/(n3.count.aa$v*n3.count.aa$n))/n3.count.diaa$vn
cps.gtg.aag<-((n3.trinuc$gtg*n3.trinuc$aag)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtg.aat<-((n3.trinuc$gtg*n3.trinuc$aat)/(n3.count.aa$v*n3.count.aa$n))/n3.count.diaa$vn

cps.gtg.aca<-((n3.trinuc$gtg*n3.trinuc$aca)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtg.acc<-((n3.trinuc$gtg*n3.trinuc$acc)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtg.acg<-((n3.trinuc$gtg*n3.trinuc$acg)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtg.act<-((n3.trinuc$gtg*n3.trinuc$act)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt

cps.gtg.aga<-((n3.trinuc$gtg*n3.trinuc$aga)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtg.agc<-((n3.trinuc$gtg*n3.trinuc$agc)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtg.agg<-((n3.trinuc$gtg*n3.trinuc$agg)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtg.agt<-((n3.trinuc$gtg*n3.trinuc$agt)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs

cps.gtg.ata<-((n3.trinuc$gtg*n3.trinuc$ata)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi
cps.gtg.atc<-((n3.trinuc$gtg*n3.trinuc$atc)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi
cps.gtg.atg<-((n3.trinuc$gtg*n3.trinuc$atg)/(n3.count.aa$v*n3.count.aa$m))/n3.count.diaa$vm
cps.gtg.att<-((n3.trinuc$gtg*n3.trinuc$att)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi

cps.gtg.caa<-((n3.trinuc$gtg*n3.trinuc$caa)/(n3.count.aa$v*n3.count.aa$q))/n3.count.diaa$vq
cps.gtg.cac<-((n3.trinuc$gtg*n3.trinuc$cac)/(n3.count.aa$v*n3.count.aa$h))/n3.count.diaa$vh
cps.gtg.cag<-((n3.trinuc$gtg*n3.trinuc$cag)/(n3.count.aa$v*n3.count.aa$q))/n3.count.diaa$vq
cps.gtg.cat<-((n3.trinuc$gtg*n3.trinuc$cat)/(n3.count.aa$v*n3.count.aa$h))/n3.count.diaa$vh

cps.gtg.cca<-((n3.trinuc$gtg*n3.trinuc$cca)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtg.ccc<-((n3.trinuc$gtg*n3.trinuc$ccc)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtg.ccg<-((n3.trinuc$gtg*n3.trinuc$ccg)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtg.cct<-((n3.trinuc$gtg*n3.trinuc$cct)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp

cps.gtg.cga<-((n3.trinuc$gtg*n3.trinuc$cga)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtg.cgc<-((n3.trinuc$gtg*n3.trinuc$cgc)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtg.cgg<-((n3.trinuc$gtg*n3.trinuc$cgg)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtg.cgt<-((n3.trinuc$gtg*n3.trinuc$cgt)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr

cps.gtg.cta<-((n3.trinuc$gtg*n3.trinuc$cta)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtg.ctc<-((n3.trinuc$gtg*n3.trinuc$ctc)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtg.ctg<-((n3.trinuc$gtg*n3.trinuc$ctg)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtg.ctt<-((n3.trinuc$gtg*n3.trinuc$ctt)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl

cps.gtg.gaa<-((n3.trinuc$gtg*n3.trinuc$gaa)/(n3.count.aa$v*n3.count.aa$e))/n3.count.diaa$ve
cps.gtg.gac<-((n3.trinuc$gtg*n3.trinuc$gac)/(n3.count.aa$v*n3.count.aa$d))/n3.count.diaa$vd
cps.gtg.gag<-((n3.trinuc$gtg*n3.trinuc$gag)/(n3.count.aa$v*n3.count.aa$e))/n3.count.diaa$ve
cps.gtg.gat<-((n3.trinuc$gtg*n3.trinuc$gat)/(n3.count.aa$v*n3.count.aa$d))/n3.count.diaa$vd

cps.gtg.gca<-((n3.trinuc$gtg*n3.trinuc$gca)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtg.gcc<-((n3.trinuc$gtg*n3.trinuc$gcc)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtg.gcg<-((n3.trinuc$gtg*n3.trinuc$gcg)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtg.gct<-((n3.trinuc$gtg*n3.trinuc$gct)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va

cps.gtg.gga<-((n3.trinuc$gtg*n3.trinuc$gga)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtg.ggc<-((n3.trinuc$gtg*n3.trinuc$ggc)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtg.ggg<-((n3.trinuc$gtg*n3.trinuc$ggg)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtg.ggt<-((n3.trinuc$gtg*n3.trinuc$ggt)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg

cps.gtg.gta<-((n3.trinuc$gtg*n3.trinuc$gta)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtg.gtc<-((n3.trinuc$gtg*n3.trinuc$gtc)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtg.gtg<-((n3.trinuc$gtg*n3.trinuc$gtg)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtg.gtt<-((n3.trinuc$gtg*n3.trinuc$gtt)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv

#Stop codon
#cps.gtg.taa<-((n3.trinuc$gtg*n3.trinuc$taa)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtg.tac<-((n3.trinuc$gtg*n3.trinuc$tac)/(n3.count.aa$v*n3.count.aa$y))/n3.count.diaa$vy
#Stop codon
#cps.gtg.tag<-((n3.trinuc$gtg*n3.trinuc$tag)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtg.tat<-((n3.trinuc$gtg*n3.trinuc$tat)/(n3.count.aa$v*n3.count.aa$y))/n3.count.diaa$vy

cps.gtg.tca<-((n3.trinuc$gtg*n3.trinuc$tca)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtg.tcc<-((n3.trinuc$gtg*n3.trinuc$tcc)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtg.tcg<-((n3.trinuc$gtg*n3.trinuc$tcg)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtg.tct<-((n3.trinuc$gtg*n3.trinuc$tct)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs

#Stop codon
#cps.gtg.tga<-((n3.trinuc$gtg*n3.trinuc$tga)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtg.tgc<-((n3.trinuc$gtg*n3.trinuc$tgc)/(n3.count.aa$v*n3.count.aa$c))/n3.count.diaa$vc
cps.gtg.tgg<-((n3.trinuc$gtg*n3.trinuc$tgg)/(n3.count.aa$v*n3.count.aa$w))/n3.count.diaa$vw
cps.gtg.tgt<-((n3.trinuc$gtg*n3.trinuc$tgt)/(n3.count.aa$v*n3.count.aa$c))/n3.count.diaa$vc

cps.gtg.tta<-((n3.trinuc$gtg*n3.trinuc$tta)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtg.ttc<-((n3.trinuc$gtg*n3.trinuc$ttc)/(n3.count.aa$v*n3.count.aa$f))/n3.count.diaa$vf
cps.gtg.ttg<-((n3.trinuc$gtg*n3.trinuc$ttg)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtg.ttt<-((n3.trinuc$gtg*n3.trinuc$ttt)/(n3.count.aa$v*n3.count.aa$f))/n3.count.diaa$vf








cps.gtt.aaa<-((n3.trinuc$gtt*n3.trinuc$aaa)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtt.aac<-((n3.trinuc$gtt*n3.trinuc$aac)/(n3.count.aa$v*n3.count.aa$n))/n3.count.diaa$vn
cps.gtt.aag<-((n3.trinuc$gtt*n3.trinuc$aag)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtt.aat<-((n3.trinuc$gtt*n3.trinuc$aat)/(n3.count.aa$v*n3.count.aa$n))/n3.count.diaa$vn

cps.gtt.aca<-((n3.trinuc$gtt*n3.trinuc$aca)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtt.acc<-((n3.trinuc$gtt*n3.trinuc$acc)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtt.acg<-((n3.trinuc$gtt*n3.trinuc$acg)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt
cps.gtt.act<-((n3.trinuc$gtt*n3.trinuc$act)/(n3.count.aa$v*n3.count.aa$t))/n3.count.diaa$vt

cps.gtt.aga<-((n3.trinuc$gtt*n3.trinuc$aga)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtt.agc<-((n3.trinuc$gtt*n3.trinuc$agc)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtt.agg<-((n3.trinuc$gtt*n3.trinuc$agg)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtt.agt<-((n3.trinuc$gtt*n3.trinuc$agt)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs

cps.gtt.ata<-((n3.trinuc$gtt*n3.trinuc$ata)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi
cps.gtt.atc<-((n3.trinuc$gtt*n3.trinuc$atc)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi
cps.gtt.atg<-((n3.trinuc$gtt*n3.trinuc$atg)/(n3.count.aa$v*n3.count.aa$m))/n3.count.diaa$vm
cps.gtt.att<-((n3.trinuc$gtt*n3.trinuc$att)/(n3.count.aa$v*n3.count.aa$i))/n3.count.diaa$vi

cps.gtt.caa<-((n3.trinuc$gtt*n3.trinuc$caa)/(n3.count.aa$v*n3.count.aa$q))/n3.count.diaa$vq
cps.gtt.cac<-((n3.trinuc$gtt*n3.trinuc$cac)/(n3.count.aa$v*n3.count.aa$h))/n3.count.diaa$vh
cps.gtt.cag<-((n3.trinuc$gtt*n3.trinuc$cag)/(n3.count.aa$v*n3.count.aa$q))/n3.count.diaa$vq
cps.gtt.cat<-((n3.trinuc$gtt*n3.trinuc$cat)/(n3.count.aa$v*n3.count.aa$h))/n3.count.diaa$vh

cps.gtt.cca<-((n3.trinuc$gtt*n3.trinuc$cca)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtt.ccc<-((n3.trinuc$gtt*n3.trinuc$ccc)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtt.ccg<-((n3.trinuc$gtt*n3.trinuc$ccg)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp
cps.gtt.cct<-((n3.trinuc$gtt*n3.trinuc$cct)/(n3.count.aa$v*n3.count.aa$p))/n3.count.diaa$vp

cps.gtt.cga<-((n3.trinuc$gtt*n3.trinuc$cga)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtt.cgc<-((n3.trinuc$gtt*n3.trinuc$cgc)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtt.cgg<-((n3.trinuc$gtt*n3.trinuc$cgg)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr
cps.gtt.cgt<-((n3.trinuc$gtt*n3.trinuc$cgt)/(n3.count.aa$v*n3.count.aa$r))/n3.count.diaa$vr

cps.gtt.cta<-((n3.trinuc$gtt*n3.trinuc$cta)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtt.ctc<-((n3.trinuc$gtt*n3.trinuc$ctc)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtt.ctg<-((n3.trinuc$gtt*n3.trinuc$ctg)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtt.ctt<-((n3.trinuc$gtt*n3.trinuc$ctt)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl

cps.gtt.gaa<-((n3.trinuc$gtt*n3.trinuc$gaa)/(n3.count.aa$v*n3.count.aa$e))/n3.count.diaa$ve
cps.gtt.gac<-((n3.trinuc$gtt*n3.trinuc$gac)/(n3.count.aa$v*n3.count.aa$d))/n3.count.diaa$vd
cps.gtt.gag<-((n3.trinuc$gtt*n3.trinuc$gag)/(n3.count.aa$v*n3.count.aa$e))/n3.count.diaa$ve
cps.gtt.gat<-((n3.trinuc$gtt*n3.trinuc$gat)/(n3.count.aa$v*n3.count.aa$d))/n3.count.diaa$vd

cps.gtt.gca<-((n3.trinuc$gtt*n3.trinuc$gca)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtt.gcc<-((n3.trinuc$gtt*n3.trinuc$gcc)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtt.gcg<-((n3.trinuc$gtt*n3.trinuc$gcg)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va
cps.gtt.gct<-((n3.trinuc$gtt*n3.trinuc$gct)/(n3.count.aa$v*n3.count.aa$a))/n3.count.diaa$va

cps.gtt.gga<-((n3.trinuc$gtt*n3.trinuc$gga)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtt.ggc<-((n3.trinuc$gtt*n3.trinuc$ggc)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtt.ggg<-((n3.trinuc$gtt*n3.trinuc$ggg)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg
cps.gtt.ggt<-((n3.trinuc$gtt*n3.trinuc$ggt)/(n3.count.aa$v*n3.count.aa$g))/n3.count.diaa$vg

cps.gtt.gta<-((n3.trinuc$gtt*n3.trinuc$gta)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtt.gtc<-((n3.trinuc$gtt*n3.trinuc$gtc)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtt.gtg<-((n3.trinuc$gtt*n3.trinuc$gtg)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv
cps.gtt.gtt<-((n3.trinuc$gtt*n3.trinuc$gtt)/(n3.count.aa$v*n3.count.aa$v))/n3.count.diaa$vv

#Stop codon
#cps.gtt.taa<-((n3.trinuc$gtt*n3.trinuc$taa)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtt.tac<-((n3.trinuc$gtt*n3.trinuc$tac)/(n3.count.aa$v*n3.count.aa$y))/n3.count.diaa$vy
#Stop codon
#cps.gtt.tag<-((n3.trinuc$gtt*n3.trinuc$tag)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtt.tat<-((n3.trinuc$gtt*n3.trinuc$tat)/(n3.count.aa$v*n3.count.aa$y))/n3.count.diaa$vy

cps.gtt.tca<-((n3.trinuc$gtt*n3.trinuc$tca)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtt.tcc<-((n3.trinuc$gtt*n3.trinuc$tcc)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtt.tcg<-((n3.trinuc$gtt*n3.trinuc$tcg)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs
cps.gtt.tct<-((n3.trinuc$gtt*n3.trinuc$tct)/(n3.count.aa$v*n3.count.aa$s))/n3.count.diaa$vs

#Stop codon
#cps.gtt.tga<-((n3.trinuc$gtt*n3.trinuc$tga)/(n3.count.aa$v*n3.count.aa$k))/n3.count.diaa$vk
cps.gtt.tgc<-((n3.trinuc$gtt*n3.trinuc$tgc)/(n3.count.aa$v*n3.count.aa$c))/n3.count.diaa$vc
cps.gtt.tgg<-((n3.trinuc$gtt*n3.trinuc$tgg)/(n3.count.aa$v*n3.count.aa$w))/n3.count.diaa$vw
cps.gtt.tgt<-((n3.trinuc$gtt*n3.trinuc$tgt)/(n3.count.aa$v*n3.count.aa$c))/n3.count.diaa$vc

cps.gtt.tta<-((n3.trinuc$gtt*n3.trinuc$tta)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtt.ttc<-((n3.trinuc$gtt*n3.trinuc$ttc)/(n3.count.aa$v*n3.count.aa$f))/n3.count.diaa$vf
cps.gtt.ttg<-((n3.trinuc$gtt*n3.trinuc$ttg)/(n3.count.aa$v*n3.count.aa$l))/n3.count.diaa$vl
cps.gtt.ttt<-((n3.trinuc$gtt*n3.trinuc$ttt)/(n3.count.aa$v*n3.count.aa$f))/n3.count.diaa$vf










#Stop codon
#cps.taa.aaa<-((n3.trinuc$taa*n3.trinuc$aaa)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.taa.aac<-((n3.trinuc$taa*n3.trinuc$aac)/(n3.count.aa$y*n3.count.aa$n))/n3.count.diaa$yn
#cps.taa.aag<-((n3.trinuc$taa*n3.trinuc$aag)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.taa.aat<-((n3.trinuc$taa*n3.trinuc$aat)/(n3.count.aa$y*n3.count.aa$n))/n3.count.diaa$yn

#cps.taa.aca<-((n3.trinuc$taa*n3.trinuc$aca)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
#cps.taa.acc<-((n3.trinuc$taa*n3.trinuc$acc)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
#cps.taa.acg<-((n3.trinuc$taa*n3.trinuc$acg)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
#cps.taa.act<-((n3.trinuc$taa*n3.trinuc$act)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt

#cps.taa.aga<-((n3.trinuc$taa*n3.trinuc$aga)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.taa.agc<-((n3.trinuc$taa*n3.trinuc$agc)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
#cps.taa.agg<-((n3.trinuc$taa*n3.trinuc$agg)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.taa.agt<-((n3.trinuc$taa*n3.trinuc$agt)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys

#cps.taa.ata<-((n3.trinuc$taa*n3.trinuc$ata)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi
#cps.taa.atc<-((n3.trinuc$taa*n3.trinuc$atc)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi
#cps.taa.atg<-((n3.trinuc$taa*n3.trinuc$atg)/(n3.count.aa$y*n3.count.aa$m))/n3.count.diaa$ym
#cps.taa.att<-((n3.trinuc$taa*n3.trinuc$att)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi

#cps.taa.caa<-((n3.trinuc$taa*n3.trinuc$caa)/(n3.count.aa$y*n3.count.aa$q))/n3.count.diaa$yq
#cps.taa.cac<-((n3.trinuc$taa*n3.trinuc$cac)/(n3.count.aa$y*n3.count.aa$h))/n3.count.diaa$yh
#cps.taa.cag<-((n3.trinuc$taa*n3.trinuc$cag)/(n3.count.aa$y*n3.count.aa$q))/n3.count.diaa$yq
#cps.taa.cat<-((n3.trinuc$taa*n3.trinuc$cat)/(n3.count.aa$y*n3.count.aa$h))/n3.count.diaa$yh

#cps.taa.cca<-((n3.trinuc$taa*n3.trinuc$cca)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
#cps.taa.ccc<-((n3.trinuc$taa*n3.trinuc$ccc)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
#cps.taa.ccg<-((n3.trinuc$taa*n3.trinuc$ccg)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
#cps.taa.cct<-((n3.trinuc$taa*n3.trinuc$cct)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp

#cps.taa.cga<-((n3.trinuc$taa*n3.trinuc$cga)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.taa.cgc<-((n3.trinuc$taa*n3.trinuc$cgc)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.taa.cgg<-((n3.trinuc$taa*n3.trinuc$cgg)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.taa.cgt<-((n3.trinuc$taa*n3.trinuc$cgt)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr

#cps.taa.cta<-((n3.trinuc$taa*n3.trinuc$cta)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.taa.ctc<-((n3.trinuc$taa*n3.trinuc$ctc)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.taa.ctg<-((n3.trinuc$taa*n3.trinuc$ctg)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.taa.ctt<-((n3.trinuc$taa*n3.trinuc$ctt)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl

#cps.taa.gaa<-((n3.trinuc$taa*n3.trinuc$gaa)/(n3.count.aa$y*n3.count.aa$e))/n3.count.diaa$ye
#cps.taa.gac<-((n3.trinuc$taa*n3.trinuc$gac)/(n3.count.aa$y*n3.count.aa$d))/n3.count.diaa$yd
#cps.taa.gag<-((n3.trinuc$taa*n3.trinuc$gag)/(n3.count.aa$y*n3.count.aa$e))/n3.count.diaa$ye
#cps.taa.gat<-((n3.trinuc$taa*n3.trinuc$gat)/(n3.count.aa$y*n3.count.aa$d))/n3.count.diaa$yd

#cps.taa.gca<-((n3.trinuc$taa*n3.trinuc$gca)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
#cps.taa.gcc<-((n3.trinuc$taa*n3.trinuc$gcc)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
#cps.taa.gcg<-((n3.trinuc$taa*n3.trinuc$gcg)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
#cps.taa.gct<-((n3.trinuc$taa*n3.trinuc$gct)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya

#cps.taa.gga<-((n3.trinuc$taa*n3.trinuc$gga)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
#cps.taa.ggc<-((n3.trinuc$taa*n3.trinuc$ggc)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
#cps.taa.ggg<-((n3.trinuc$taa*n3.trinuc$ggg)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
#cps.taa.ggt<-((n3.trinuc$taa*n3.trinuc$ggt)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg

#cps.taa.gta<-((n3.trinuc$taa*n3.trinuc$gta)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
#cps.taa.gtc<-((n3.trinuc$taa*n3.trinuc$gtc)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
#cps.taa.gtg<-((n3.trinuc$taa*n3.trinuc$gtg)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
#cps.taa.gtt<-((n3.trinuc$taa*n3.trinuc$gtt)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv

#Stop codon
#cps.taa.taa<-((n3.trinuc$taa*n3.trinuc$taa)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.taa.tac<-((n3.trinuc$taa*n3.trinuc$tac)/(n3.count.aa$y*n3.count.aa$y))/n3.count.diaa$yy
#Stop codon
#cps.taa.tag<-((n3.trinuc$taa*n3.trinuc$tag)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.taa.tat<-((n3.trinuc$taa*n3.trinuc$tat)/(n3.count.aa$y*n3.count.aa$y))/n3.count.diaa$yy

#cps.taa.tca<-((n3.trinuc$taa*n3.trinuc$tca)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
#cps.taa.tcc<-((n3.trinuc$taa*n3.trinuc$tcc)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
#cps.taa.tcg<-((n3.trinuc$taa*n3.trinuc$tcg)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
#cps.taa.tct<-((n3.trinuc$taa*n3.trinuc$tct)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys

#Stop codon
#cps.taa.tga<-((n3.trinuc$taa*n3.trinuc$tga)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.taa.tgc<-((n3.trinuc$taa*n3.trinuc$tgc)/(n3.count.aa$y*n3.count.aa$c))/n3.count.diaa$yc
#cps.taa.tgg<-((n3.trinuc$taa*n3.trinuc$tgg)/(n3.count.aa$y*n3.count.aa$w))/n3.count.diaa$yw
#cps.taa.tgt<-((n3.trinuc$taa*n3.trinuc$tgt)/(n3.count.aa$y*n3.count.aa$c))/n3.count.diaa$yc

#cps.taa.tta<-((n3.trinuc$taa*n3.trinuc$tta)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.taa.ttc<-((n3.trinuc$taa*n3.trinuc$ttc)/(n3.count.aa$y*n3.count.aa$f))/n3.count.diaa$yf
#cps.taa.ttg<-((n3.trinuc$taa*n3.trinuc$ttg)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.taa.ttt<-((n3.trinuc$taa*n3.trinuc$ttt)/(n3.count.aa$y*n3.count.aa$f))/n3.count.diaa$yf



















cps.tac.aaa<-((n3.trinuc$tac*n3.trinuc$aaa)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tac.aac<-((n3.trinuc$tac*n3.trinuc$aac)/(n3.count.aa$y*n3.count.aa$n))/n3.count.diaa$yn
cps.tac.aag<-((n3.trinuc$tac*n3.trinuc$aag)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tac.aat<-((n3.trinuc$tac*n3.trinuc$aat)/(n3.count.aa$y*n3.count.aa$n))/n3.count.diaa$yn

cps.tac.aca<-((n3.trinuc$tac*n3.trinuc$aca)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
cps.tac.acc<-((n3.trinuc$tac*n3.trinuc$acc)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
cps.tac.acg<-((n3.trinuc$tac*n3.trinuc$acg)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
cps.tac.act<-((n3.trinuc$tac*n3.trinuc$act)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt

cps.tac.aga<-((n3.trinuc$tac*n3.trinuc$aga)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tac.agc<-((n3.trinuc$tac*n3.trinuc$agc)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
cps.tac.agg<-((n3.trinuc$tac*n3.trinuc$agg)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tac.agt<-((n3.trinuc$tac*n3.trinuc$agt)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys

cps.tac.ata<-((n3.trinuc$tac*n3.trinuc$ata)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi
cps.tac.atc<-((n3.trinuc$tac*n3.trinuc$atc)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi
cps.tac.atg<-((n3.trinuc$tac*n3.trinuc$atg)/(n3.count.aa$y*n3.count.aa$m))/n3.count.diaa$ym
cps.tac.att<-((n3.trinuc$tac*n3.trinuc$att)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi

cps.tac.caa<-((n3.trinuc$tac*n3.trinuc$caa)/(n3.count.aa$y*n3.count.aa$q))/n3.count.diaa$yq
cps.tac.cac<-((n3.trinuc$tac*n3.trinuc$cac)/(n3.count.aa$y*n3.count.aa$h))/n3.count.diaa$yh
cps.tac.cag<-((n3.trinuc$tac*n3.trinuc$cag)/(n3.count.aa$y*n3.count.aa$q))/n3.count.diaa$yq
cps.tac.cat<-((n3.trinuc$tac*n3.trinuc$cat)/(n3.count.aa$y*n3.count.aa$h))/n3.count.diaa$yh

cps.tac.cca<-((n3.trinuc$tac*n3.trinuc$cca)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
cps.tac.ccc<-((n3.trinuc$tac*n3.trinuc$ccc)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
cps.tac.ccg<-((n3.trinuc$tac*n3.trinuc$ccg)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
cps.tac.cct<-((n3.trinuc$tac*n3.trinuc$cct)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp

cps.tac.cga<-((n3.trinuc$tac*n3.trinuc$cga)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tac.cgc<-((n3.trinuc$tac*n3.trinuc$cgc)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tac.cgg<-((n3.trinuc$tac*n3.trinuc$cgg)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tac.cgt<-((n3.trinuc$tac*n3.trinuc$cgt)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr

cps.tac.cta<-((n3.trinuc$tac*n3.trinuc$cta)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tac.ctc<-((n3.trinuc$tac*n3.trinuc$ctc)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tac.ctg<-((n3.trinuc$tac*n3.trinuc$ctg)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tac.ctt<-((n3.trinuc$tac*n3.trinuc$ctt)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl

cps.tac.gaa<-((n3.trinuc$tac*n3.trinuc$gaa)/(n3.count.aa$y*n3.count.aa$e))/n3.count.diaa$ye
cps.tac.gac<-((n3.trinuc$tac*n3.trinuc$gac)/(n3.count.aa$y*n3.count.aa$d))/n3.count.diaa$yd
cps.tac.gag<-((n3.trinuc$tac*n3.trinuc$gag)/(n3.count.aa$y*n3.count.aa$e))/n3.count.diaa$ye
cps.tac.gat<-((n3.trinuc$tac*n3.trinuc$gat)/(n3.count.aa$y*n3.count.aa$d))/n3.count.diaa$yd

cps.tac.gca<-((n3.trinuc$tac*n3.trinuc$gca)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
cps.tac.gcc<-((n3.trinuc$tac*n3.trinuc$gcc)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
cps.tac.gcg<-((n3.trinuc$tac*n3.trinuc$gcg)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
cps.tac.gct<-((n3.trinuc$tac*n3.trinuc$gct)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya

cps.tac.gga<-((n3.trinuc$tac*n3.trinuc$gga)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
cps.tac.ggc<-((n3.trinuc$tac*n3.trinuc$ggc)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
cps.tac.ggg<-((n3.trinuc$tac*n3.trinuc$ggg)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
cps.tac.ggt<-((n3.trinuc$tac*n3.trinuc$ggt)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg

cps.tac.gta<-((n3.trinuc$tac*n3.trinuc$gta)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
cps.tac.gtc<-((n3.trinuc$tac*n3.trinuc$gtc)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
cps.tac.gtg<-((n3.trinuc$tac*n3.trinuc$gtg)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
cps.tac.gtt<-((n3.trinuc$tac*n3.trinuc$gtt)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv

#Stop codon
#cps.tac.taa<-((n3.trinuc$tac*n3.trinuc$taa)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tac.tac<-((n3.trinuc$tac*n3.trinuc$tac)/(n3.count.aa$y*n3.count.aa$y))/n3.count.diaa$yy
#Stop codon
#cps.tac.tag<-((n3.trinuc$tac*n3.trinuc$tag)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tac.tat<-((n3.trinuc$tac*n3.trinuc$tat)/(n3.count.aa$y*n3.count.aa$y))/n3.count.diaa$yy

cps.tac.tca<-((n3.trinuc$tac*n3.trinuc$tca)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
cps.tac.tcc<-((n3.trinuc$tac*n3.trinuc$tcc)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
cps.tac.tcg<-((n3.trinuc$tac*n3.trinuc$tcg)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
cps.tac.tct<-((n3.trinuc$tac*n3.trinuc$tct)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys

#Stop codon
#cps.tac.tga<-((n3.trinuc$tac*n3.trinuc$tga)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tac.tgc<-((n3.trinuc$tac*n3.trinuc$tgc)/(n3.count.aa$y*n3.count.aa$c))/n3.count.diaa$yc
cps.tac.tgg<-((n3.trinuc$tac*n3.trinuc$tgg)/(n3.count.aa$y*n3.count.aa$w))/n3.count.diaa$yw
cps.tac.tgt<-((n3.trinuc$tac*n3.trinuc$tgt)/(n3.count.aa$y*n3.count.aa$c))/n3.count.diaa$yc

cps.tac.tta<-((n3.trinuc$tac*n3.trinuc$tta)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tac.ttc<-((n3.trinuc$tac*n3.trinuc$ttc)/(n3.count.aa$y*n3.count.aa$f))/n3.count.diaa$yf
cps.tac.ttg<-((n3.trinuc$tac*n3.trinuc$ttg)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tac.ttt<-((n3.trinuc$tac*n3.trinuc$ttt)/(n3.count.aa$y*n3.count.aa$f))/n3.count.diaa$yf








#Stop codon

#cps.tag.aaa<-((n3.trinuc$tag*n3.trinuc$aaa)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.tag.aac<-((n3.trinuc$tag*n3.trinuc$aac)/(n3.count.aa$y*n3.count.aa$n))/n3.count.diaa$yn
#cps.tag.aag<-((n3.trinuc$tag*n3.trinuc$aag)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.tag.aat<-((n3.trinuc$tag*n3.trinuc$aat)/(n3.count.aa$y*n3.count.aa$n))/n3.count.diaa$yn

#cps.tag.aca<-((n3.trinuc$tag*n3.trinuc$aca)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
#cps.tag.acc<-((n3.trinuc$tag*n3.trinuc$acc)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
#cps.tag.acg<-((n3.trinuc$tag*n3.trinuc$acg)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
#cps.tag.act<-((n3.trinuc$tag*n3.trinuc$act)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt

#cps.tag.aga<-((n3.trinuc$tag*n3.trinuc$aga)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.tag.agc<-((n3.trinuc$tag*n3.trinuc$agc)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
#cps.tag.agg<-((n3.trinuc$tag*n3.trinuc$agg)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.tag.agt<-((n3.trinuc$tag*n3.trinuc$agt)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys

#cps.tag.ata<-((n3.trinuc$tag*n3.trinuc$ata)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi
#cps.tag.atc<-((n3.trinuc$tag*n3.trinuc$atc)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi
#cps.tag.atg<-((n3.trinuc$tag*n3.trinuc$atg)/(n3.count.aa$y*n3.count.aa$m))/n3.count.diaa$ym
#cps.tag.att<-((n3.trinuc$tag*n3.trinuc$att)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi

#cps.tag.caa<-((n3.trinuc$tag*n3.trinuc$caa)/(n3.count.aa$y*n3.count.aa$q))/n3.count.diaa$yq
#cps.tag.cac<-((n3.trinuc$tag*n3.trinuc$cac)/(n3.count.aa$y*n3.count.aa$h))/n3.count.diaa$yh
#cps.tag.cag<-((n3.trinuc$tag*n3.trinuc$cag)/(n3.count.aa$y*n3.count.aa$q))/n3.count.diaa$yq
#cps.tag.cat<-((n3.trinuc$tag*n3.trinuc$cat)/(n3.count.aa$y*n3.count.aa$h))/n3.count.diaa$yh

#cps.tag.cca<-((n3.trinuc$tag*n3.trinuc$cca)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
#cps.tag.ccc<-((n3.trinuc$tag*n3.trinuc$ccc)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
#cps.tag.ccg<-((n3.trinuc$tag*n3.trinuc$ccg)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
#cps.tag.cct<-((n3.trinuc$tag*n3.trinuc$cct)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp

#cps.tag.cga<-((n3.trinuc$tag*n3.trinuc$cga)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.tag.cgc<-((n3.trinuc$tag*n3.trinuc$cgc)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.tag.cgg<-((n3.trinuc$tag*n3.trinuc$cgg)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
#cps.tag.cgt<-((n3.trinuc$tag*n3.trinuc$cgt)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr

#cps.tag.cta<-((n3.trinuc$tag*n3.trinuc$cta)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.tag.ctc<-((n3.trinuc$tag*n3.trinuc$ctc)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.tag.ctg<-((n3.trinuc$tag*n3.trinuc$ctg)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.tag.ctt<-((n3.trinuc$tag*n3.trinuc$ctt)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl

#cps.tag.gaa<-((n3.trinuc$tag*n3.trinuc$gaa)/(n3.count.aa$y*n3.count.aa$e))/n3.count.diaa$ye
#cps.tag.gac<-((n3.trinuc$tag*n3.trinuc$gac)/(n3.count.aa$y*n3.count.aa$d))/n3.count.diaa$yd
#cps.tag.gag<-((n3.trinuc$tag*n3.trinuc$gag)/(n3.count.aa$y*n3.count.aa$e))/n3.count.diaa$ye
#cps.tag.gat<-((n3.trinuc$tag*n3.trinuc$gat)/(n3.count.aa$y*n3.count.aa$d))/n3.count.diaa$yd

#cps.tag.gca<-((n3.trinuc$tag*n3.trinuc$gca)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
#cps.tag.gcc<-((n3.trinuc$tag*n3.trinuc$gcc)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
#cps.tag.gcg<-((n3.trinuc$tag*n3.trinuc$gcg)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
#cps.tag.gct<-((n3.trinuc$tag*n3.trinuc$gct)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya

#cps.tag.gga<-((n3.trinuc$tag*n3.trinuc$gga)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
#cps.tag.ggc<-((n3.trinuc$tag*n3.trinuc$ggc)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
#cps.tag.ggg<-((n3.trinuc$tag*n3.trinuc$ggg)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
#cps.tag.ggt<-((n3.trinuc$tag*n3.trinuc$ggt)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg

#cps.tag.gta<-((n3.trinuc$tag*n3.trinuc$gta)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
#cps.tag.gtc<-((n3.trinuc$tag*n3.trinuc$gtc)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
#cps.tag.gtg<-((n3.trinuc$tag*n3.trinuc$gtg)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
#cps.tag.gtt<-((n3.trinuc$tag*n3.trinuc$gtt)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv

#Stop codon
#cps.tag.taa<-((n3.trinuc$tag*n3.trinuc$taa)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.tag.tac<-((n3.trinuc$tag*n3.trinuc$tac)/(n3.count.aa$y*n3.count.aa$y))/n3.count.diaa$yy
#Stop codon
#cps.tag.tag<-((n3.trinuc$tag*n3.trinuc$tag)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.tag.tat<-((n3.trinuc$tag*n3.trinuc$tat)/(n3.count.aa$y*n3.count.aa$y))/n3.count.diaa$yy

#cps.tag.tca<-((n3.trinuc$tag*n3.trinuc$tca)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
#cps.tag.tcc<-((n3.trinuc$tag*n3.trinuc$tcc)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
#cps.tag.tcg<-((n3.trinuc$tag*n3.trinuc$tcg)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
#cps.tag.tct<-((n3.trinuc$tag*n3.trinuc$tct)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys

#Stop codon
#cps.tag.tga<-((n3.trinuc$tag*n3.trinuc$tga)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
#cps.tag.tgc<-((n3.trinuc$tag*n3.trinuc$tgc)/(n3.count.aa$y*n3.count.aa$c))/n3.count.diaa$yc
#cps.tag.tgg<-((n3.trinuc$tag*n3.trinuc$tgg)/(n3.count.aa$y*n3.count.aa$w))/n3.count.diaa$yw
#cps.tag.tgt<-((n3.trinuc$tag*n3.trinuc$tgt)/(n3.count.aa$y*n3.count.aa$c))/n3.count.diaa$yc

#cps.tag.tta<-((n3.trinuc$tag*n3.trinuc$tta)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.tag.ttc<-((n3.trinuc$tag*n3.trinuc$ttc)/(n3.count.aa$y*n3.count.aa$f))/n3.count.diaa$yf
#cps.tag.ttg<-((n3.trinuc$tag*n3.trinuc$ttg)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
#cps.tag.ttt<-((n3.trinuc$tag*n3.trinuc$ttt)/(n3.count.aa$y*n3.count.aa$f))/n3.count.diaa$yf

















cps.tat.aaa<-((n3.trinuc$tat*n3.trinuc$aaa)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tat.aac<-((n3.trinuc$tat*n3.trinuc$aac)/(n3.count.aa$y*n3.count.aa$n))/n3.count.diaa$yn
cps.tat.aag<-((n3.trinuc$tat*n3.trinuc$aag)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tat.aat<-((n3.trinuc$tat*n3.trinuc$aat)/(n3.count.aa$y*n3.count.aa$n))/n3.count.diaa$yn

cps.tat.aca<-((n3.trinuc$tat*n3.trinuc$aca)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
cps.tat.acc<-((n3.trinuc$tat*n3.trinuc$acc)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
cps.tat.acg<-((n3.trinuc$tat*n3.trinuc$acg)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt
cps.tat.act<-((n3.trinuc$tat*n3.trinuc$act)/(n3.count.aa$y*n3.count.aa$t))/n3.count.diaa$yt

cps.tat.aga<-((n3.trinuc$tat*n3.trinuc$aga)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tat.agc<-((n3.trinuc$tat*n3.trinuc$agc)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
cps.tat.agg<-((n3.trinuc$tat*n3.trinuc$agg)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tat.agt<-((n3.trinuc$tat*n3.trinuc$agt)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys

cps.tat.ata<-((n3.trinuc$tat*n3.trinuc$ata)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi
cps.tat.atc<-((n3.trinuc$tat*n3.trinuc$atc)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi
cps.tat.atg<-((n3.trinuc$tat*n3.trinuc$atg)/(n3.count.aa$y*n3.count.aa$m))/n3.count.diaa$ym
cps.tat.att<-((n3.trinuc$tat*n3.trinuc$att)/(n3.count.aa$y*n3.count.aa$i))/n3.count.diaa$yi

cps.tat.caa<-((n3.trinuc$tat*n3.trinuc$caa)/(n3.count.aa$y*n3.count.aa$q))/n3.count.diaa$yq
cps.tat.cac<-((n3.trinuc$tat*n3.trinuc$cac)/(n3.count.aa$y*n3.count.aa$h))/n3.count.diaa$yh
cps.tat.cag<-((n3.trinuc$tat*n3.trinuc$cag)/(n3.count.aa$y*n3.count.aa$q))/n3.count.diaa$yq
cps.tat.cat<-((n3.trinuc$tat*n3.trinuc$cat)/(n3.count.aa$y*n3.count.aa$h))/n3.count.diaa$yh

cps.tat.cca<-((n3.trinuc$tat*n3.trinuc$cca)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
cps.tat.ccc<-((n3.trinuc$tat*n3.trinuc$ccc)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
cps.tat.ccg<-((n3.trinuc$tat*n3.trinuc$ccg)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp
cps.tat.cct<-((n3.trinuc$tat*n3.trinuc$cct)/(n3.count.aa$y*n3.count.aa$p))/n3.count.diaa$yp

cps.tat.cga<-((n3.trinuc$tat*n3.trinuc$cga)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tat.cgc<-((n3.trinuc$tat*n3.trinuc$cgc)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tat.cgg<-((n3.trinuc$tat*n3.trinuc$cgg)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr
cps.tat.cgt<-((n3.trinuc$tat*n3.trinuc$cgt)/(n3.count.aa$y*n3.count.aa$r))/n3.count.diaa$yr

cps.tat.cta<-((n3.trinuc$tat*n3.trinuc$cta)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tat.ctc<-((n3.trinuc$tat*n3.trinuc$ctc)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tat.ctg<-((n3.trinuc$tat*n3.trinuc$ctg)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tat.ctt<-((n3.trinuc$tat*n3.trinuc$ctt)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl

cps.tat.gaa<-((n3.trinuc$tat*n3.trinuc$gaa)/(n3.count.aa$y*n3.count.aa$e))/n3.count.diaa$ye
cps.tat.gac<-((n3.trinuc$tat*n3.trinuc$gac)/(n3.count.aa$y*n3.count.aa$d))/n3.count.diaa$yd
cps.tat.gag<-((n3.trinuc$tat*n3.trinuc$gag)/(n3.count.aa$y*n3.count.aa$e))/n3.count.diaa$ye
cps.tat.gat<-((n3.trinuc$tat*n3.trinuc$gat)/(n3.count.aa$y*n3.count.aa$d))/n3.count.diaa$yd

cps.tat.gca<-((n3.trinuc$tat*n3.trinuc$gca)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
cps.tat.gcc<-((n3.trinuc$tat*n3.trinuc$gcc)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
cps.tat.gcg<-((n3.trinuc$tat*n3.trinuc$gcg)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya
cps.tat.gct<-((n3.trinuc$tat*n3.trinuc$gct)/(n3.count.aa$y*n3.count.aa$a))/n3.count.diaa$ya

cps.tat.gga<-((n3.trinuc$tat*n3.trinuc$gga)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
cps.tat.ggc<-((n3.trinuc$tat*n3.trinuc$ggc)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
cps.tat.ggg<-((n3.trinuc$tat*n3.trinuc$ggg)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg
cps.tat.ggt<-((n3.trinuc$tat*n3.trinuc$ggt)/(n3.count.aa$y*n3.count.aa$g))/n3.count.diaa$yg

cps.tat.gta<-((n3.trinuc$tat*n3.trinuc$gta)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
cps.tat.gtc<-((n3.trinuc$tat*n3.trinuc$gtc)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
cps.tat.gtg<-((n3.trinuc$tat*n3.trinuc$gtg)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv
cps.tat.gtt<-((n3.trinuc$tat*n3.trinuc$gtt)/(n3.count.aa$y*n3.count.aa$v))/n3.count.diaa$yv

#Stop codon
#cps.tat.taa<-((n3.trinuc$tat*n3.trinuc$taa)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tat.tac<-((n3.trinuc$tat*n3.trinuc$tac)/(n3.count.aa$y*n3.count.aa$y))/n3.count.diaa$yy
#Stop codon
#cps.tat.tag<-((n3.trinuc$tat*n3.trinuc$tag)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tat.tat<-((n3.trinuc$tat*n3.trinuc$tat)/(n3.count.aa$y*n3.count.aa$y))/n3.count.diaa$yy

cps.tat.tca<-((n3.trinuc$tat*n3.trinuc$tca)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
cps.tat.tcc<-((n3.trinuc$tat*n3.trinuc$tcc)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
cps.tat.tcg<-((n3.trinuc$tat*n3.trinuc$tcg)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys
cps.tat.tct<-((n3.trinuc$tat*n3.trinuc$tct)/(n3.count.aa$y*n3.count.aa$s))/n3.count.diaa$ys

#Stop codon
#cps.tat.tga<-((n3.trinuc$tat*n3.trinuc$tga)/(n3.count.aa$y*n3.count.aa$k))/n3.count.diaa$yk
cps.tat.tgc<-((n3.trinuc$tat*n3.trinuc$tgc)/(n3.count.aa$y*n3.count.aa$c))/n3.count.diaa$yc
cps.tat.tgg<-((n3.trinuc$tat*n3.trinuc$tgg)/(n3.count.aa$y*n3.count.aa$w))/n3.count.diaa$yw
cps.tat.tgt<-((n3.trinuc$tat*n3.trinuc$tgt)/(n3.count.aa$y*n3.count.aa$c))/n3.count.diaa$yc

cps.tat.tta<-((n3.trinuc$tat*n3.trinuc$tta)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tat.ttc<-((n3.trinuc$tat*n3.trinuc$ttc)/(n3.count.aa$y*n3.count.aa$f))/n3.count.diaa$yf
cps.tat.ttg<-((n3.trinuc$tat*n3.trinuc$ttg)/(n3.count.aa$y*n3.count.aa$l))/n3.count.diaa$yl
cps.tat.ttt<-((n3.trinuc$tat*n3.trinuc$ttt)/(n3.count.aa$y*n3.count.aa$f))/n3.count.diaa$yf







cps.tca.aaa<-((n3.trinuc$tca*n3.trinuc$aaa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tca.aac<-((n3.trinuc$tca*n3.trinuc$aac)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn
cps.tca.aag<-((n3.trinuc$tca*n3.trinuc$aag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tca.aat<-((n3.trinuc$tca*n3.trinuc$aat)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn

cps.tca.aca<-((n3.trinuc$tca*n3.trinuc$aca)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tca.acc<-((n3.trinuc$tca*n3.trinuc$acc)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tca.acg<-((n3.trinuc$tca*n3.trinuc$acg)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tca.act<-((n3.trinuc$tca*n3.trinuc$act)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st

cps.tca.aga<-((n3.trinuc$tca*n3.trinuc$aga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tca.agc<-((n3.trinuc$tca*n3.trinuc$agc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tca.agg<-((n3.trinuc$tca*n3.trinuc$agg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tca.agt<-((n3.trinuc$tca*n3.trinuc$agt)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

cps.tca.ata<-((n3.trinuc$tca*n3.trinuc$ata)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.tca.atc<-((n3.trinuc$tca*n3.trinuc$atc)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.tca.atg<-((n3.trinuc$tca*n3.trinuc$atg)/(n3.count.aa$s*n3.count.aa$m))/n3.count.diaa$sm
cps.tca.att<-((n3.trinuc$tca*n3.trinuc$att)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si

cps.tca.caa<-((n3.trinuc$tca*n3.trinuc$caa)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.tca.cac<-((n3.trinuc$tca*n3.trinuc$cac)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh
cps.tca.cag<-((n3.trinuc$tca*n3.trinuc$cag)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.tca.cat<-((n3.trinuc$tca*n3.trinuc$cat)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh

cps.tca.cca<-((n3.trinuc$tca*n3.trinuc$cca)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tca.ccc<-((n3.trinuc$tca*n3.trinuc$ccc)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tca.ccg<-((n3.trinuc$tca*n3.trinuc$ccg)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tca.cct<-((n3.trinuc$tca*n3.trinuc$cct)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp

cps.tca.cga<-((n3.trinuc$tca*n3.trinuc$cga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tca.cgc<-((n3.trinuc$tca*n3.trinuc$cgc)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tca.cgg<-((n3.trinuc$tca*n3.trinuc$cgg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tca.cgt<-((n3.trinuc$tca*n3.trinuc$cgt)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr

cps.tca.cta<-((n3.trinuc$tca*n3.trinuc$cta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tca.ctc<-((n3.trinuc$tca*n3.trinuc$ctc)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tca.ctg<-((n3.trinuc$tca*n3.trinuc$ctg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tca.ctt<-((n3.trinuc$tca*n3.trinuc$ctt)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl

cps.tca.gaa<-((n3.trinuc$tca*n3.trinuc$gaa)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.tca.gac<-((n3.trinuc$tca*n3.trinuc$gac)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd
cps.tca.gag<-((n3.trinuc$tca*n3.trinuc$gag)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.tca.gat<-((n3.trinuc$tca*n3.trinuc$gat)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd

cps.tca.gca<-((n3.trinuc$tca*n3.trinuc$gca)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tca.gcc<-((n3.trinuc$tca*n3.trinuc$gcc)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tca.gcg<-((n3.trinuc$tca*n3.trinuc$gcg)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tca.gct<-((n3.trinuc$tca*n3.trinuc$gct)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa

cps.tca.gga<-((n3.trinuc$tca*n3.trinuc$gga)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tca.ggc<-((n3.trinuc$tca*n3.trinuc$ggc)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tca.ggg<-((n3.trinuc$tca*n3.trinuc$ggg)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tca.ggt<-((n3.trinuc$tca*n3.trinuc$ggt)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg

cps.tca.gta<-((n3.trinuc$tca*n3.trinuc$gta)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tca.gtc<-((n3.trinuc$tca*n3.trinuc$gtc)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tca.gtg<-((n3.trinuc$tca*n3.trinuc$gtg)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tca.gtt<-((n3.trinuc$tca*n3.trinuc$gtt)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv

#Stop codon
#cps.tca.taa<-((n3.trinuc$tca*n3.trinuc$taa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tca.tac<-((n3.trinuc$tca*n3.trinuc$tac)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy
#Stop codon
#cps.tca.tag<-((n3.trinuc$tca*n3.trinuc$tag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tca.tat<-((n3.trinuc$tca*n3.trinuc$tat)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy

cps.tca.tca<-((n3.trinuc$tca*n3.trinuc$tca)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tca.tcc<-((n3.trinuc$tca*n3.trinuc$tcc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tca.tcg<-((n3.trinuc$tca*n3.trinuc$tcg)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tca.tct<-((n3.trinuc$tca*n3.trinuc$tct)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

#Stop codon
#cps.tca.tga<-((n3.trinuc$tca*n3.trinuc$tga)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tca.tgc<-((n3.trinuc$tca*n3.trinuc$tgc)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc
cps.tca.tgg<-((n3.trinuc$tca*n3.trinuc$tgg)/(n3.count.aa$s*n3.count.aa$w))/n3.count.diaa$sw
cps.tca.tgt<-((n3.trinuc$tca*n3.trinuc$tgt)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc

cps.tca.tta<-((n3.trinuc$tca*n3.trinuc$tta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tca.ttc<-((n3.trinuc$tca*n3.trinuc$ttc)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf
cps.tca.ttg<-((n3.trinuc$tca*n3.trinuc$ttg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tca.ttt<-((n3.trinuc$tca*n3.trinuc$ttt)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf








cps.tcc.aaa<-((n3.trinuc$tcc*n3.trinuc$aaa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcc.aac<-((n3.trinuc$tcc*n3.trinuc$aac)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn
cps.tcc.aag<-((n3.trinuc$tcc*n3.trinuc$aag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcc.aat<-((n3.trinuc$tcc*n3.trinuc$aat)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn

cps.tcc.aca<-((n3.trinuc$tcc*n3.trinuc$aca)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tcc.acc<-((n3.trinuc$tcc*n3.trinuc$acc)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tcc.acg<-((n3.trinuc$tcc*n3.trinuc$acg)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tcc.act<-((n3.trinuc$tcc*n3.trinuc$act)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st

cps.tcc.aga<-((n3.trinuc$tcc*n3.trinuc$aga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcc.agc<-((n3.trinuc$tcc*n3.trinuc$agc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tcc.agg<-((n3.trinuc$tcc*n3.trinuc$agg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcc.agt<-((n3.trinuc$tcc*n3.trinuc$agt)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

cps.tcc.ata<-((n3.trinuc$tcc*n3.trinuc$ata)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.tcc.atc<-((n3.trinuc$tcc*n3.trinuc$atc)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.tcc.atg<-((n3.trinuc$tcc*n3.trinuc$atg)/(n3.count.aa$s*n3.count.aa$m))/n3.count.diaa$sm
cps.tcc.att<-((n3.trinuc$tcc*n3.trinuc$att)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si

cps.tcc.caa<-((n3.trinuc$tcc*n3.trinuc$caa)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.tcc.cac<-((n3.trinuc$tcc*n3.trinuc$cac)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh
cps.tcc.cag<-((n3.trinuc$tcc*n3.trinuc$cag)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.tcc.cat<-((n3.trinuc$tcc*n3.trinuc$cat)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh

cps.tcc.cca<-((n3.trinuc$tcc*n3.trinuc$cca)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tcc.ccc<-((n3.trinuc$tcc*n3.trinuc$ccc)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tcc.ccg<-((n3.trinuc$tcc*n3.trinuc$ccg)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tcc.cct<-((n3.trinuc$tcc*n3.trinuc$cct)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp

cps.tcc.cga<-((n3.trinuc$tcc*n3.trinuc$cga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcc.cgc<-((n3.trinuc$tcc*n3.trinuc$cgc)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcc.cgg<-((n3.trinuc$tcc*n3.trinuc$cgg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcc.cgt<-((n3.trinuc$tcc*n3.trinuc$cgt)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr

cps.tcc.cta<-((n3.trinuc$tcc*n3.trinuc$cta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcc.ctc<-((n3.trinuc$tcc*n3.trinuc$ctc)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcc.ctg<-((n3.trinuc$tcc*n3.trinuc$ctg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcc.ctt<-((n3.trinuc$tcc*n3.trinuc$ctt)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl

cps.tcc.gaa<-((n3.trinuc$tcc*n3.trinuc$gaa)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.tcc.gac<-((n3.trinuc$tcc*n3.trinuc$gac)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd
cps.tcc.gag<-((n3.trinuc$tcc*n3.trinuc$gag)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.tcc.gat<-((n3.trinuc$tcc*n3.trinuc$gat)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd

cps.tcc.gca<-((n3.trinuc$tcc*n3.trinuc$gca)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tcc.gcc<-((n3.trinuc$tcc*n3.trinuc$gcc)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tcc.gcg<-((n3.trinuc$tcc*n3.trinuc$gcg)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tcc.gct<-((n3.trinuc$tcc*n3.trinuc$gct)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa

cps.tcc.gga<-((n3.trinuc$tcc*n3.trinuc$gga)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tcc.ggc<-((n3.trinuc$tcc*n3.trinuc$ggc)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tcc.ggg<-((n3.trinuc$tcc*n3.trinuc$ggg)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tcc.ggt<-((n3.trinuc$tcc*n3.trinuc$ggt)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg

cps.tcc.gta<-((n3.trinuc$tcc*n3.trinuc$gta)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tcc.gtc<-((n3.trinuc$tcc*n3.trinuc$gtc)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tcc.gtg<-((n3.trinuc$tcc*n3.trinuc$gtg)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tcc.gtt<-((n3.trinuc$tcc*n3.trinuc$gtt)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv

#Stop codon
#cps.tcc.taa<-((n3.trinuc$tcc*n3.trinuc$taa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcc.tac<-((n3.trinuc$tcc*n3.trinuc$tac)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy
#Stop codon
#cps.tcc.tag<-((n3.trinuc$tcc*n3.trinuc$tag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcc.tat<-((n3.trinuc$tcc*n3.trinuc$tat)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy

cps.tcc.tca<-((n3.trinuc$tcc*n3.trinuc$tca)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tcc.tcc<-((n3.trinuc$tcc*n3.trinuc$tcc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tcc.tcg<-((n3.trinuc$tcc*n3.trinuc$tcg)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tcc.tct<-((n3.trinuc$tcc*n3.trinuc$tct)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

#Stop codon
#cps.tcc.tga<-((n3.trinuc$tcc*n3.trinuc$tga)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcc.tgc<-((n3.trinuc$tcc*n3.trinuc$tgc)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc
cps.tcc.tgg<-((n3.trinuc$tcc*n3.trinuc$tgg)/(n3.count.aa$s*n3.count.aa$w))/n3.count.diaa$sw
cps.tcc.tgt<-((n3.trinuc$tcc*n3.trinuc$tgt)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc

cps.tcc.tta<-((n3.trinuc$tcc*n3.trinuc$tta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcc.ttc<-((n3.trinuc$tcc*n3.trinuc$ttc)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf
cps.tcc.ttg<-((n3.trinuc$tcc*n3.trinuc$ttg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcc.ttt<-((n3.trinuc$tcc*n3.trinuc$ttt)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf








cps.tcg.aaa<-((n3.trinuc$tcg*n3.trinuc$aaa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcg.aac<-((n3.trinuc$tcg*n3.trinuc$aac)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn
cps.tcg.aag<-((n3.trinuc$tcg*n3.trinuc$aag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcg.aat<-((n3.trinuc$tcg*n3.trinuc$aat)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn

cps.tcg.aca<-((n3.trinuc$tcg*n3.trinuc$aca)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tcg.acc<-((n3.trinuc$tcg*n3.trinuc$acc)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tcg.acg<-((n3.trinuc$tcg*n3.trinuc$acg)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tcg.act<-((n3.trinuc$tcg*n3.trinuc$act)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st

cps.tcg.aga<-((n3.trinuc$tcg*n3.trinuc$aga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcg.agc<-((n3.trinuc$tcg*n3.trinuc$agc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tcg.agg<-((n3.trinuc$tcg*n3.trinuc$agg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcg.agt<-((n3.trinuc$tcg*n3.trinuc$agt)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

cps.tcg.ata<-((n3.trinuc$tcg*n3.trinuc$ata)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.tcg.atc<-((n3.trinuc$tcg*n3.trinuc$atc)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.tcg.atg<-((n3.trinuc$tcg*n3.trinuc$atg)/(n3.count.aa$s*n3.count.aa$m))/n3.count.diaa$sm
cps.tcg.att<-((n3.trinuc$tcg*n3.trinuc$att)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si

cps.tcg.caa<-((n3.trinuc$tcg*n3.trinuc$caa)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.tcg.cac<-((n3.trinuc$tcg*n3.trinuc$cac)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh
cps.tcg.cag<-((n3.trinuc$tcg*n3.trinuc$cag)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.tcg.cat<-((n3.trinuc$tcg*n3.trinuc$cat)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh

cps.tcg.cca<-((n3.trinuc$tcg*n3.trinuc$cca)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tcg.ccc<-((n3.trinuc$tcg*n3.trinuc$ccc)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tcg.ccg<-((n3.trinuc$tcg*n3.trinuc$ccg)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tcg.cct<-((n3.trinuc$tcg*n3.trinuc$cct)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp

cps.tcg.cga<-((n3.trinuc$tcg*n3.trinuc$cga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcg.cgc<-((n3.trinuc$tcg*n3.trinuc$cgc)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcg.cgg<-((n3.trinuc$tcg*n3.trinuc$cgg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tcg.cgt<-((n3.trinuc$tcg*n3.trinuc$cgt)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr

cps.tcg.cta<-((n3.trinuc$tcg*n3.trinuc$cta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcg.ctc<-((n3.trinuc$tcg*n3.trinuc$ctc)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcg.ctg<-((n3.trinuc$tcg*n3.trinuc$ctg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcg.ctt<-((n3.trinuc$tcg*n3.trinuc$ctt)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl

cps.tcg.gaa<-((n3.trinuc$tcg*n3.trinuc$gaa)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.tcg.gac<-((n3.trinuc$tcg*n3.trinuc$gac)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd
cps.tcg.gag<-((n3.trinuc$tcg*n3.trinuc$gag)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.tcg.gat<-((n3.trinuc$tcg*n3.trinuc$gat)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd

cps.tcg.gca<-((n3.trinuc$tcg*n3.trinuc$gca)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tcg.gcc<-((n3.trinuc$tcg*n3.trinuc$gcc)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tcg.gcg<-((n3.trinuc$tcg*n3.trinuc$gcg)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tcg.gct<-((n3.trinuc$tcg*n3.trinuc$gct)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa

cps.tcg.gga<-((n3.trinuc$tcg*n3.trinuc$gga)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tcg.ggc<-((n3.trinuc$tcg*n3.trinuc$ggc)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tcg.ggg<-((n3.trinuc$tcg*n3.trinuc$ggg)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tcg.ggt<-((n3.trinuc$tcg*n3.trinuc$ggt)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg

cps.tcg.gta<-((n3.trinuc$tcg*n3.trinuc$gta)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tcg.gtc<-((n3.trinuc$tcg*n3.trinuc$gtc)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tcg.gtg<-((n3.trinuc$tcg*n3.trinuc$gtg)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tcg.gtt<-((n3.trinuc$tcg*n3.trinuc$gtt)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv

#Stop codon
#cps.tcg.taa<-((n3.trinuc$tcg*n3.trinuc$taa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcg.tac<-((n3.trinuc$tcg*n3.trinuc$tac)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy
#Stop codon
#cps.tcg.tag<-((n3.trinuc$tcg*n3.trinuc$tag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcg.tat<-((n3.trinuc$tcg*n3.trinuc$tat)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy

cps.tcg.tca<-((n3.trinuc$tcg*n3.trinuc$tca)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tcg.tcc<-((n3.trinuc$tcg*n3.trinuc$tcc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tcg.tcg<-((n3.trinuc$tcg*n3.trinuc$tcg)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tcg.tct<-((n3.trinuc$tcg*n3.trinuc$tct)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

#Stop codon
#cps.tcg.tga<-((n3.trinuc$tcg*n3.trinuc$tga)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tcg.tgc<-((n3.trinuc$tcg*n3.trinuc$tgc)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc
cps.tcg.tgg<-((n3.trinuc$tcg*n3.trinuc$tgg)/(n3.count.aa$s*n3.count.aa$w))/n3.count.diaa$sw
cps.tcg.tgt<-((n3.trinuc$tcg*n3.trinuc$tgt)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc

cps.tcg.tta<-((n3.trinuc$tcg*n3.trinuc$tta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcg.ttc<-((n3.trinuc$tcg*n3.trinuc$ttc)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf
cps.tcg.ttg<-((n3.trinuc$tcg*n3.trinuc$ttg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tcg.ttt<-((n3.trinuc$tcg*n3.trinuc$ttt)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf








cps.tct.aaa<-((n3.trinuc$tct*n3.trinuc$aaa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tct.aac<-((n3.trinuc$tct*n3.trinuc$aac)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn
cps.tct.aag<-((n3.trinuc$tct*n3.trinuc$aag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tct.aat<-((n3.trinuc$tct*n3.trinuc$aat)/(n3.count.aa$s*n3.count.aa$n))/n3.count.diaa$sn

cps.tct.aca<-((n3.trinuc$tct*n3.trinuc$aca)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tct.acc<-((n3.trinuc$tct*n3.trinuc$acc)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tct.acg<-((n3.trinuc$tct*n3.trinuc$acg)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st
cps.tct.act<-((n3.trinuc$tct*n3.trinuc$act)/(n3.count.aa$s*n3.count.aa$t))/n3.count.diaa$st

cps.tct.aga<-((n3.trinuc$tct*n3.trinuc$aga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tct.agc<-((n3.trinuc$tct*n3.trinuc$agc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tct.agg<-((n3.trinuc$tct*n3.trinuc$agg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tct.agt<-((n3.trinuc$tct*n3.trinuc$agt)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

cps.tct.ata<-((n3.trinuc$tct*n3.trinuc$ata)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.tct.atc<-((n3.trinuc$tct*n3.trinuc$atc)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si
cps.tct.atg<-((n3.trinuc$tct*n3.trinuc$atg)/(n3.count.aa$s*n3.count.aa$m))/n3.count.diaa$sm
cps.tct.att<-((n3.trinuc$tct*n3.trinuc$att)/(n3.count.aa$s*n3.count.aa$i))/n3.count.diaa$si

cps.tct.caa<-((n3.trinuc$tct*n3.trinuc$caa)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.tct.cac<-((n3.trinuc$tct*n3.trinuc$cac)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh
cps.tct.cag<-((n3.trinuc$tct*n3.trinuc$cag)/(n3.count.aa$s*n3.count.aa$q))/n3.count.diaa$sq
cps.tct.cat<-((n3.trinuc$tct*n3.trinuc$cat)/(n3.count.aa$s*n3.count.aa$h))/n3.count.diaa$sh

cps.tct.cca<-((n3.trinuc$tct*n3.trinuc$cca)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tct.ccc<-((n3.trinuc$tct*n3.trinuc$ccc)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tct.ccg<-((n3.trinuc$tct*n3.trinuc$ccg)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp
cps.tct.cct<-((n3.trinuc$tct*n3.trinuc$cct)/(n3.count.aa$s*n3.count.aa$p))/n3.count.diaa$sp

cps.tct.cga<-((n3.trinuc$tct*n3.trinuc$cga)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tct.cgc<-((n3.trinuc$tct*n3.trinuc$cgc)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tct.cgg<-((n3.trinuc$tct*n3.trinuc$cgg)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr
cps.tct.cgt<-((n3.trinuc$tct*n3.trinuc$cgt)/(n3.count.aa$s*n3.count.aa$r))/n3.count.diaa$sr

cps.tct.cta<-((n3.trinuc$tct*n3.trinuc$cta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tct.ctc<-((n3.trinuc$tct*n3.trinuc$ctc)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tct.ctg<-((n3.trinuc$tct*n3.trinuc$ctg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tct.ctt<-((n3.trinuc$tct*n3.trinuc$ctt)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl

cps.tct.gaa<-((n3.trinuc$tct*n3.trinuc$gaa)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.tct.gac<-((n3.trinuc$tct*n3.trinuc$gac)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd
cps.tct.gag<-((n3.trinuc$tct*n3.trinuc$gag)/(n3.count.aa$s*n3.count.aa$e))/n3.count.diaa$se
cps.tct.gat<-((n3.trinuc$tct*n3.trinuc$gat)/(n3.count.aa$s*n3.count.aa$d))/n3.count.diaa$sd

cps.tct.gca<-((n3.trinuc$tct*n3.trinuc$gca)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tct.gcc<-((n3.trinuc$tct*n3.trinuc$gcc)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tct.gcg<-((n3.trinuc$tct*n3.trinuc$gcg)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa
cps.tct.gct<-((n3.trinuc$tct*n3.trinuc$gct)/(n3.count.aa$s*n3.count.aa$a))/n3.count.diaa$sa

cps.tct.gga<-((n3.trinuc$tct*n3.trinuc$gga)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tct.ggc<-((n3.trinuc$tct*n3.trinuc$ggc)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tct.ggg<-((n3.trinuc$tct*n3.trinuc$ggg)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg
cps.tct.ggt<-((n3.trinuc$tct*n3.trinuc$ggt)/(n3.count.aa$s*n3.count.aa$g))/n3.count.diaa$sg

cps.tct.gta<-((n3.trinuc$tct*n3.trinuc$gta)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tct.gtc<-((n3.trinuc$tct*n3.trinuc$gtc)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tct.gtg<-((n3.trinuc$tct*n3.trinuc$gtg)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv
cps.tct.gtt<-((n3.trinuc$tct*n3.trinuc$gtt)/(n3.count.aa$s*n3.count.aa$v))/n3.count.diaa$sv

#Stop codon
#cps.tct.taa<-((n3.trinuc$tct*n3.trinuc$taa)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tct.tac<-((n3.trinuc$tct*n3.trinuc$tac)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy
#Stop codon
#cps.tct.tag<-((n3.trinuc$tct*n3.trinuc$tag)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tct.tat<-((n3.trinuc$tct*n3.trinuc$tat)/(n3.count.aa$s*n3.count.aa$y))/n3.count.diaa$sy

cps.tct.tca<-((n3.trinuc$tct*n3.trinuc$tca)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tct.tcc<-((n3.trinuc$tct*n3.trinuc$tcc)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tct.tcg<-((n3.trinuc$tct*n3.trinuc$tcg)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss
cps.tct.tct<-((n3.trinuc$tct*n3.trinuc$tct)/(n3.count.aa$s*n3.count.aa$s))/n3.count.diaa$ss

#Stop codon
#cps.tct.tga<-((n3.trinuc$tct*n3.trinuc$tga)/(n3.count.aa$s*n3.count.aa$k))/n3.count.diaa$sk
cps.tct.tgc<-((n3.trinuc$tct*n3.trinuc$tgc)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc
cps.tct.tgg<-((n3.trinuc$tct*n3.trinuc$tgg)/(n3.count.aa$s*n3.count.aa$w))/n3.count.diaa$sw
cps.tct.tgt<-((n3.trinuc$tct*n3.trinuc$tgt)/(n3.count.aa$s*n3.count.aa$c))/n3.count.diaa$sc

cps.tct.tta<-((n3.trinuc$tct*n3.trinuc$tta)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tct.ttc<-((n3.trinuc$tct*n3.trinuc$ttc)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf
cps.tct.ttg<-((n3.trinuc$tct*n3.trinuc$ttg)/(n3.count.aa$s*n3.count.aa$l))/n3.count.diaa$sl
cps.tct.ttt<-((n3.trinuc$tct*n3.trinuc$ttt)/(n3.count.aa$s*n3.count.aa$f))/n3.count.diaa$sf











#Stop codon
#cps.tga.aaa<-((n3.trinuc$tga*n3.trinuc$aaa)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
#cps.tga.aac<-((n3.trinuc$tga*n3.trinuc$aac)/(n3.count.aa$c*n3.count.aa$n))/n3.count.diaa$cn
#cps.tga.aag<-((n3.trinuc$tga*n3.trinuc$aag)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
#cps.tga.aat<-((n3.trinuc$tga*n3.trinuc$aat)/(n3.count.aa$c*n3.count.aa$n))/n3.count.diaa$cn

#cps.tga.aca<-((n3.trinuc$tga*n3.trinuc$aca)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
#cps.tga.acc<-((n3.trinuc$tga*n3.trinuc$acc)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
#cps.tga.acg<-((n3.trinuc$tga*n3.trinuc$acg)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
#cps.tga.act<-((n3.trinuc$tga*n3.trinuc$act)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct

#cps.tga.aga<-((n3.trinuc$tga*n3.trinuc$aga)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
#cps.tga.agc<-((n3.trinuc$tga*n3.trinuc$agc)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
#cps.tga.agg<-((n3.trinuc$tga*n3.trinuc$agg)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
#cps.tga.agt<-((n3.trinuc$tga*n3.trinuc$agt)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs

#cps.tga.ata<-((n3.trinuc$tga*n3.trinuc$ata)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci
#cps.tga.atc<-((n3.trinuc$tga*n3.trinuc$atc)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci
#cps.tga.atg<-((n3.trinuc$tga*n3.trinuc$atg)/(n3.count.aa$c*n3.count.aa$m))/n3.count.diaa$cm
#cps.tga.att<-((n3.trinuc$tga*n3.trinuc$att)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci

#cps.tga.caa<-((n3.trinuc$tga*n3.trinuc$caa)/(n3.count.aa$c*n3.count.aa$q))/n3.count.diaa$cq
#cps.tga.cac<-((n3.trinuc$tga*n3.trinuc$cac)/(n3.count.aa$c*n3.count.aa$h))/n3.count.diaa$ch
#cps.tga.cag<-((n3.trinuc$tga*n3.trinuc$cag)/(n3.count.aa$c*n3.count.aa$q))/n3.count.diaa$cq
#cps.tga.cat<-((n3.trinuc$tga*n3.trinuc$cat)/(n3.count.aa$c*n3.count.aa$h))/n3.count.diaa$ch

#cps.tga.cca<-((n3.trinuc$tga*n3.trinuc$cca)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
#cps.tga.ccc<-((n3.trinuc$tga*n3.trinuc$ccc)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
#cps.tga.ccg<-((n3.trinuc$tga*n3.trinuc$ccg)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
#cps.tga.cct<-((n3.trinuc$tga*n3.trinuc$cct)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp

#cps.tga.cga<-((n3.trinuc$tga*n3.trinuc$cga)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
#cps.tga.cgc<-((n3.trinuc$tga*n3.trinuc$cgc)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
#cps.tga.cgg<-((n3.trinuc$tga*n3.trinuc$cgg)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
#cps.tga.cgt<-((n3.trinuc$tga*n3.trinuc$cgt)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr

#cps.tga.cta<-((n3.trinuc$tga*n3.trinuc$cta)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
#cps.tga.ctc<-((n3.trinuc$tga*n3.trinuc$ctc)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
#cps.tga.ctg<-((n3.trinuc$tga*n3.trinuc$ctg)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
#cps.tga.ctt<-((n3.trinuc$tga*n3.trinuc$ctt)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl

#cps.tga.gaa<-((n3.trinuc$tga*n3.trinuc$gaa)/(n3.count.aa$c*n3.count.aa$e))/n3.count.diaa$ce
#cps.tga.gac<-((n3.trinuc$tga*n3.trinuc$gac)/(n3.count.aa$c*n3.count.aa$d))/n3.count.diaa$cd
#cps.tga.gag<-((n3.trinuc$tga*n3.trinuc$gag)/(n3.count.aa$c*n3.count.aa$e))/n3.count.diaa$ce
#cps.tga.gat<-((n3.trinuc$tga*n3.trinuc$gat)/(n3.count.aa$c*n3.count.aa$d))/n3.count.diaa$cd

#cps.tga.gca<-((n3.trinuc$tga*n3.trinuc$gca)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
#cps.tga.gcc<-((n3.trinuc$tga*n3.trinuc$gcc)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
#cps.tga.gcg<-((n3.trinuc$tga*n3.trinuc$gcg)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
#cps.tga.gct<-((n3.trinuc$tga*n3.trinuc$gct)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca

#cps.tga.gga<-((n3.trinuc$tga*n3.trinuc$gga)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
#cps.tga.ggc<-((n3.trinuc$tga*n3.trinuc$ggc)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
#cps.tga.ggg<-((n3.trinuc$tga*n3.trinuc$ggg)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
#cps.tga.ggt<-((n3.trinuc$tga*n3.trinuc$ggt)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg

#cps.tga.gta<-((n3.trinuc$tga*n3.trinuc$gta)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
#cps.tga.gtc<-((n3.trinuc$tga*n3.trinuc$gtc)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
#cps.tga.gtg<-((n3.trinuc$tga*n3.trinuc$gtg)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
#cps.tga.gtt<-((n3.trinuc$tga*n3.trinuc$gtt)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv

#Stop codon
#cps.tga.taa<-((n3.trinuc$tga*n3.trinuc$taa)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
#cps.tga.tac<-((n3.trinuc$tga*n3.trinuc$tac)/(n3.count.aa$c*n3.count.aa$y))/n3.count.diaa$cy
#Stop codon
#cps.tga.tag<-((n3.trinuc$tga*n3.trinuc$tag)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
#cps.tga.tat<-((n3.trinuc$tga*n3.trinuc$tat)/(n3.count.aa$c*n3.count.aa$y))/n3.count.diaa$cy

#cps.tga.tca<-((n3.trinuc$tga*n3.trinuc$tca)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
#cps.tga.tcc<-((n3.trinuc$tga*n3.trinuc$tcc)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
#cps.tga.tcg<-((n3.trinuc$tga*n3.trinuc$tcg)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
#cps.tga.tct<-((n3.trinuc$tga*n3.trinuc$tct)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs

#Stop codon
#cps.tga.tga<-((n3.trinuc$tga*n3.trinuc$tga)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
#cps.tga.tgc<-((n3.trinuc$tga*n3.trinuc$tgc)/(n3.count.aa$c*n3.count.aa$c))/n3.count.diaa$cc
#cps.tga.tgg<-((n3.trinuc$tga*n3.trinuc$tgg)/(n3.count.aa$c*n3.count.aa$w))/n3.count.diaa$cw
#cps.tga.tgt<-((n3.trinuc$tga*n3.trinuc$tgt)/(n3.count.aa$c*n3.count.aa$c))/n3.count.diaa$cc

#cps.tga.tta<-((n3.trinuc$tga*n3.trinuc$tta)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
#cps.tga.ttc<-((n3.trinuc$tga*n3.trinuc$ttc)/(n3.count.aa$c*n3.count.aa$f))/n3.count.diaa$cf
#cps.tga.ttg<-((n3.trinuc$tga*n3.trinuc$ttg)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
#cps.tga.ttt<-((n3.trinuc$tga*n3.trinuc$ttt)/(n3.count.aa$c*n3.count.aa$f))/n3.count.diaa$cf












cps.tgc.aaa<-((n3.trinuc$tgc*n3.trinuc$aaa)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgc.aac<-((n3.trinuc$tgc*n3.trinuc$aac)/(n3.count.aa$c*n3.count.aa$n))/n3.count.diaa$cn
cps.tgc.aag<-((n3.trinuc$tgc*n3.trinuc$aag)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgc.aat<-((n3.trinuc$tgc*n3.trinuc$aat)/(n3.count.aa$c*n3.count.aa$n))/n3.count.diaa$cn

cps.tgc.aca<-((n3.trinuc$tgc*n3.trinuc$aca)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
cps.tgc.acc<-((n3.trinuc$tgc*n3.trinuc$acc)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
cps.tgc.acg<-((n3.trinuc$tgc*n3.trinuc$acg)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
cps.tgc.act<-((n3.trinuc$tgc*n3.trinuc$act)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct

cps.tgc.aga<-((n3.trinuc$tgc*n3.trinuc$aga)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgc.agc<-((n3.trinuc$tgc*n3.trinuc$agc)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
cps.tgc.agg<-((n3.trinuc$tgc*n3.trinuc$agg)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgc.agt<-((n3.trinuc$tgc*n3.trinuc$agt)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs

cps.tgc.ata<-((n3.trinuc$tgc*n3.trinuc$ata)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci
cps.tgc.atc<-((n3.trinuc$tgc*n3.trinuc$atc)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci
cps.tgc.atg<-((n3.trinuc$tgc*n3.trinuc$atg)/(n3.count.aa$c*n3.count.aa$m))/n3.count.diaa$cm
cps.tgc.att<-((n3.trinuc$tgc*n3.trinuc$att)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci

cps.tgc.caa<-((n3.trinuc$tgc*n3.trinuc$caa)/(n3.count.aa$c*n3.count.aa$q))/n3.count.diaa$cq
cps.tgc.cac<-((n3.trinuc$tgc*n3.trinuc$cac)/(n3.count.aa$c*n3.count.aa$h))/n3.count.diaa$ch
cps.tgc.cag<-((n3.trinuc$tgc*n3.trinuc$cag)/(n3.count.aa$c*n3.count.aa$q))/n3.count.diaa$cq
cps.tgc.cat<-((n3.trinuc$tgc*n3.trinuc$cat)/(n3.count.aa$c*n3.count.aa$h))/n3.count.diaa$ch

cps.tgc.cca<-((n3.trinuc$tgc*n3.trinuc$cca)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
cps.tgc.ccc<-((n3.trinuc$tgc*n3.trinuc$ccc)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
cps.tgc.ccg<-((n3.trinuc$tgc*n3.trinuc$ccg)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
cps.tgc.cct<-((n3.trinuc$tgc*n3.trinuc$cct)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp

cps.tgc.cga<-((n3.trinuc$tgc*n3.trinuc$cga)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgc.cgc<-((n3.trinuc$tgc*n3.trinuc$cgc)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgc.cgg<-((n3.trinuc$tgc*n3.trinuc$cgg)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgc.cgt<-((n3.trinuc$tgc*n3.trinuc$cgt)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr

cps.tgc.cta<-((n3.trinuc$tgc*n3.trinuc$cta)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgc.ctc<-((n3.trinuc$tgc*n3.trinuc$ctc)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgc.ctg<-((n3.trinuc$tgc*n3.trinuc$ctg)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgc.ctt<-((n3.trinuc$tgc*n3.trinuc$ctt)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl

cps.tgc.gaa<-((n3.trinuc$tgc*n3.trinuc$gaa)/(n3.count.aa$c*n3.count.aa$e))/n3.count.diaa$ce
cps.tgc.gac<-((n3.trinuc$tgc*n3.trinuc$gac)/(n3.count.aa$c*n3.count.aa$d))/n3.count.diaa$cd
cps.tgc.gag<-((n3.trinuc$tgc*n3.trinuc$gag)/(n3.count.aa$c*n3.count.aa$e))/n3.count.diaa$ce
cps.tgc.gat<-((n3.trinuc$tgc*n3.trinuc$gat)/(n3.count.aa$c*n3.count.aa$d))/n3.count.diaa$cd

cps.tgc.gca<-((n3.trinuc$tgc*n3.trinuc$gca)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
cps.tgc.gcc<-((n3.trinuc$tgc*n3.trinuc$gcc)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
cps.tgc.gcg<-((n3.trinuc$tgc*n3.trinuc$gcg)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
cps.tgc.gct<-((n3.trinuc$tgc*n3.trinuc$gct)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca

cps.tgc.gga<-((n3.trinuc$tgc*n3.trinuc$gga)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
cps.tgc.ggc<-((n3.trinuc$tgc*n3.trinuc$ggc)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
cps.tgc.ggg<-((n3.trinuc$tgc*n3.trinuc$ggg)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
cps.tgc.ggt<-((n3.trinuc$tgc*n3.trinuc$ggt)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg

cps.tgc.gta<-((n3.trinuc$tgc*n3.trinuc$gta)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
cps.tgc.gtc<-((n3.trinuc$tgc*n3.trinuc$gtc)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
cps.tgc.gtg<-((n3.trinuc$tgc*n3.trinuc$gtg)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
cps.tgc.gtt<-((n3.trinuc$tgc*n3.trinuc$gtt)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv

#Stop codon
#cps.tgc.taa<-((n3.trinuc$tgc*n3.trinuc$taa)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgc.tac<-((n3.trinuc$tgc*n3.trinuc$tac)/(n3.count.aa$c*n3.count.aa$y))/n3.count.diaa$cy
#Stop codon
#cps.tgc.tag<-((n3.trinuc$tgc*n3.trinuc$tag)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgc.tat<-((n3.trinuc$tgc*n3.trinuc$tat)/(n3.count.aa$c*n3.count.aa$y))/n3.count.diaa$cy

cps.tgc.tca<-((n3.trinuc$tgc*n3.trinuc$tca)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
cps.tgc.tcc<-((n3.trinuc$tgc*n3.trinuc$tcc)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
cps.tgc.tcg<-((n3.trinuc$tgc*n3.trinuc$tcg)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
cps.tgc.tct<-((n3.trinuc$tgc*n3.trinuc$tct)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs

#Stop codon
#cps.tgc.tga<-((n3.trinuc$tgc*n3.trinuc$tga)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgc.tgc<-((n3.trinuc$tgc*n3.trinuc$tgc)/(n3.count.aa$c*n3.count.aa$c))/n3.count.diaa$cc
cps.tgc.tgg<-((n3.trinuc$tgc*n3.trinuc$tgg)/(n3.count.aa$c*n3.count.aa$w))/n3.count.diaa$cw
cps.tgc.tgt<-((n3.trinuc$tgc*n3.trinuc$tgt)/(n3.count.aa$c*n3.count.aa$c))/n3.count.diaa$cc

cps.tgc.tta<-((n3.trinuc$tgc*n3.trinuc$tta)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgc.ttc<-((n3.trinuc$tgc*n3.trinuc$ttc)/(n3.count.aa$c*n3.count.aa$f))/n3.count.diaa$cf
cps.tgc.ttg<-((n3.trinuc$tgc*n3.trinuc$ttg)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgc.ttt<-((n3.trinuc$tgc*n3.trinuc$ttt)/(n3.count.aa$c*n3.count.aa$f))/n3.count.diaa$cf












cps.tgg.aaa<-((n3.trinuc$tgg*n3.trinuc$aaa)/(n3.count.aa$w*n3.count.aa$k))/n3.count.diaa$wk
cps.tgg.aac<-((n3.trinuc$tgg*n3.trinuc$aac)/(n3.count.aa$w*n3.count.aa$n))/n3.count.diaa$wn
cps.tgg.aag<-((n3.trinuc$tgg*n3.trinuc$aag)/(n3.count.aa$w*n3.count.aa$k))/n3.count.diaa$wk
cps.tgg.aat<-((n3.trinuc$tgg*n3.trinuc$aat)/(n3.count.aa$w*n3.count.aa$n))/n3.count.diaa$wn

cps.tgg.aca<-((n3.trinuc$tgg*n3.trinuc$aca)/(n3.count.aa$w*n3.count.aa$t))/n3.count.diaa$wt
cps.tgg.acc<-((n3.trinuc$tgg*n3.trinuc$acc)/(n3.count.aa$w*n3.count.aa$t))/n3.count.diaa$wt
cps.tgg.acg<-((n3.trinuc$tgg*n3.trinuc$acg)/(n3.count.aa$w*n3.count.aa$t))/n3.count.diaa$wt
cps.tgg.act<-((n3.trinuc$tgg*n3.trinuc$act)/(n3.count.aa$w*n3.count.aa$t))/n3.count.diaa$wt

cps.tgg.aga<-((n3.trinuc$tgg*n3.trinuc$aga)/(n3.count.aa$w*n3.count.aa$r))/n3.count.diaa$wr
cps.tgg.agc<-((n3.trinuc$tgg*n3.trinuc$agc)/(n3.count.aa$w*n3.count.aa$s))/n3.count.diaa$ws
cps.tgg.agg<-((n3.trinuc$tgg*n3.trinuc$agg)/(n3.count.aa$w*n3.count.aa$r))/n3.count.diaa$wr
cps.tgg.agt<-((n3.trinuc$tgg*n3.trinuc$agt)/(n3.count.aa$w*n3.count.aa$s))/n3.count.diaa$ws

cps.tgg.ata<-((n3.trinuc$tgg*n3.trinuc$ata)/(n3.count.aa$w*n3.count.aa$i))/n3.count.diaa$wi
cps.tgg.atc<-((n3.trinuc$tgg*n3.trinuc$atc)/(n3.count.aa$w*n3.count.aa$i))/n3.count.diaa$wi
cps.tgg.atg<-((n3.trinuc$tgg*n3.trinuc$atg)/(n3.count.aa$w*n3.count.aa$m))/n3.count.diaa$wm
cps.tgg.att<-((n3.trinuc$tgg*n3.trinuc$att)/(n3.count.aa$w*n3.count.aa$i))/n3.count.diaa$wi

cps.tgg.caa<-((n3.trinuc$tgg*n3.trinuc$caa)/(n3.count.aa$w*n3.count.aa$q))/n3.count.diaa$wq
cps.tgg.cac<-((n3.trinuc$tgg*n3.trinuc$cac)/(n3.count.aa$w*n3.count.aa$h))/n3.count.diaa$wh
cps.tgg.cag<-((n3.trinuc$tgg*n3.trinuc$cag)/(n3.count.aa$w*n3.count.aa$q))/n3.count.diaa$wq
cps.tgg.cat<-((n3.trinuc$tgg*n3.trinuc$cat)/(n3.count.aa$w*n3.count.aa$h))/n3.count.diaa$wh

cps.tgg.cca<-((n3.trinuc$tgg*n3.trinuc$cca)/(n3.count.aa$w*n3.count.aa$p))/n3.count.diaa$wp
cps.tgg.ccc<-((n3.trinuc$tgg*n3.trinuc$ccc)/(n3.count.aa$w*n3.count.aa$p))/n3.count.diaa$wp
cps.tgg.ccg<-((n3.trinuc$tgg*n3.trinuc$ccg)/(n3.count.aa$w*n3.count.aa$p))/n3.count.diaa$wp
cps.tgg.cct<-((n3.trinuc$tgg*n3.trinuc$cct)/(n3.count.aa$w*n3.count.aa$p))/n3.count.diaa$wp

cps.tgg.cga<-((n3.trinuc$tgg*n3.trinuc$cga)/(n3.count.aa$w*n3.count.aa$r))/n3.count.diaa$wr
cps.tgg.cgc<-((n3.trinuc$tgg*n3.trinuc$cgc)/(n3.count.aa$w*n3.count.aa$r))/n3.count.diaa$wr
cps.tgg.cgg<-((n3.trinuc$tgg*n3.trinuc$cgg)/(n3.count.aa$w*n3.count.aa$r))/n3.count.diaa$wr
cps.tgg.cgt<-((n3.trinuc$tgg*n3.trinuc$cgt)/(n3.count.aa$w*n3.count.aa$r))/n3.count.diaa$wr

cps.tgg.cta<-((n3.trinuc$tgg*n3.trinuc$cta)/(n3.count.aa$w*n3.count.aa$l))/n3.count.diaa$wl
cps.tgg.ctc<-((n3.trinuc$tgg*n3.trinuc$ctc)/(n3.count.aa$w*n3.count.aa$l))/n3.count.diaa$wl
cps.tgg.ctg<-((n3.trinuc$tgg*n3.trinuc$ctg)/(n3.count.aa$w*n3.count.aa$l))/n3.count.diaa$wl
cps.tgg.ctt<-((n3.trinuc$tgg*n3.trinuc$ctt)/(n3.count.aa$w*n3.count.aa$l))/n3.count.diaa$wl

cps.tgg.gaa<-((n3.trinuc$tgg*n3.trinuc$gaa)/(n3.count.aa$w*n3.count.aa$e))/n3.count.diaa$we
cps.tgg.gac<-((n3.trinuc$tgg*n3.trinuc$gac)/(n3.count.aa$w*n3.count.aa$d))/n3.count.diaa$wd
cps.tgg.gag<-((n3.trinuc$tgg*n3.trinuc$gag)/(n3.count.aa$w*n3.count.aa$e))/n3.count.diaa$we
cps.tgg.gat<-((n3.trinuc$tgg*n3.trinuc$gat)/(n3.count.aa$w*n3.count.aa$d))/n3.count.diaa$wd

cps.tgg.gca<-((n3.trinuc$tgg*n3.trinuc$gca)/(n3.count.aa$w*n3.count.aa$a))/n3.count.diaa$wa
cps.tgg.gcc<-((n3.trinuc$tgg*n3.trinuc$gcc)/(n3.count.aa$w*n3.count.aa$a))/n3.count.diaa$wa
cps.tgg.gcg<-((n3.trinuc$tgg*n3.trinuc$gcg)/(n3.count.aa$w*n3.count.aa$a))/n3.count.diaa$wa
cps.tgg.gct<-((n3.trinuc$tgg*n3.trinuc$gct)/(n3.count.aa$w*n3.count.aa$a))/n3.count.diaa$wa

cps.tgg.gga<-((n3.trinuc$tgg*n3.trinuc$gga)/(n3.count.aa$w*n3.count.aa$g))/n3.count.diaa$wg
cps.tgg.ggc<-((n3.trinuc$tgg*n3.trinuc$ggc)/(n3.count.aa$w*n3.count.aa$g))/n3.count.diaa$wg
cps.tgg.ggg<-((n3.trinuc$tgg*n3.trinuc$ggg)/(n3.count.aa$w*n3.count.aa$g))/n3.count.diaa$wg
cps.tgg.ggt<-((n3.trinuc$tgg*n3.trinuc$ggt)/(n3.count.aa$w*n3.count.aa$g))/n3.count.diaa$wg

cps.tgg.gta<-((n3.trinuc$tgg*n3.trinuc$gta)/(n3.count.aa$w*n3.count.aa$v))/n3.count.diaa$wv
cps.tgg.gtc<-((n3.trinuc$tgg*n3.trinuc$gtc)/(n3.count.aa$w*n3.count.aa$v))/n3.count.diaa$wv
cps.tgg.gtg<-((n3.trinuc$tgg*n3.trinuc$gtg)/(n3.count.aa$w*n3.count.aa$v))/n3.count.diaa$wv
cps.tgg.gtt<-((n3.trinuc$tgg*n3.trinuc$gtt)/(n3.count.aa$w*n3.count.aa$v))/n3.count.diaa$wv

#Stop codon
#cps.tgg.taa<-((n3.trinuc$tgg*n3.trinuc$taa)/(n3.count.aa$w*n3.count.aa$k))/n3.count.diaa$wk
cps.tgg.tac<-((n3.trinuc$tgg*n3.trinuc$tac)/(n3.count.aa$w*n3.count.aa$y))/n3.count.diaa$wy
#Stop codon
#cps.tgg.tag<-((n3.trinuc$tgg*n3.trinuc$tag)/(n3.count.aa$w*n3.count.aa$k))/n3.count.diaa$wk
cps.tgg.tat<-((n3.trinuc$tgg*n3.trinuc$tat)/(n3.count.aa$w*n3.count.aa$y))/n3.count.diaa$wy

cps.tgg.tca<-((n3.trinuc$tgg*n3.trinuc$tca)/(n3.count.aa$w*n3.count.aa$s))/n3.count.diaa$ws
cps.tgg.tcc<-((n3.trinuc$tgg*n3.trinuc$tcc)/(n3.count.aa$w*n3.count.aa$s))/n3.count.diaa$ws
cps.tgg.tcg<-((n3.trinuc$tgg*n3.trinuc$tcg)/(n3.count.aa$w*n3.count.aa$s))/n3.count.diaa$ws
cps.tgg.tct<-((n3.trinuc$tgg*n3.trinuc$tct)/(n3.count.aa$w*n3.count.aa$s))/n3.count.diaa$ws

#Stop codon
#cps.tgg.tga<-((n3.trinuc$tgg*n3.trinuc$tga)/(n3.count.aa$w*n3.count.aa$k))/n3.count.diaa$wk
cps.tgg.tgc<-((n3.trinuc$tgg*n3.trinuc$tgc)/(n3.count.aa$w*n3.count.aa$c))/n3.count.diaa$wc
cps.tgg.tgg<-((n3.trinuc$tgg*n3.trinuc$tgg)/(n3.count.aa$w*n3.count.aa$w))/n3.count.diaa$ww
cps.tgg.tgt<-((n3.trinuc$tgg*n3.trinuc$tgt)/(n3.count.aa$w*n3.count.aa$c))/n3.count.diaa$wc

cps.tgg.tta<-((n3.trinuc$tgg*n3.trinuc$tta)/(n3.count.aa$w*n3.count.aa$l))/n3.count.diaa$wl
cps.tgg.ttc<-((n3.trinuc$tgg*n3.trinuc$ttc)/(n3.count.aa$w*n3.count.aa$f))/n3.count.diaa$wf
cps.tgg.ttg<-((n3.trinuc$tgg*n3.trinuc$ttg)/(n3.count.aa$w*n3.count.aa$l))/n3.count.diaa$wl
cps.tgg.ttt<-((n3.trinuc$tgg*n3.trinuc$ttt)/(n3.count.aa$w*n3.count.aa$f))/n3.count.diaa$wf








cps.tgt.aaa<-((n3.trinuc$tgt*n3.trinuc$aaa)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgt.aac<-((n3.trinuc$tgt*n3.trinuc$aac)/(n3.count.aa$c*n3.count.aa$n))/n3.count.diaa$cn
cps.tgt.aag<-((n3.trinuc$tgt*n3.trinuc$aag)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgt.aat<-((n3.trinuc$tgt*n3.trinuc$aat)/(n3.count.aa$c*n3.count.aa$n))/n3.count.diaa$cn

cps.tgt.aca<-((n3.trinuc$tgt*n3.trinuc$aca)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
cps.tgt.acc<-((n3.trinuc$tgt*n3.trinuc$acc)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
cps.tgt.acg<-((n3.trinuc$tgt*n3.trinuc$acg)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct
cps.tgt.act<-((n3.trinuc$tgt*n3.trinuc$act)/(n3.count.aa$c*n3.count.aa$t))/n3.count.diaa$ct

cps.tgt.aga<-((n3.trinuc$tgt*n3.trinuc$aga)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgt.agc<-((n3.trinuc$tgt*n3.trinuc$agc)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
cps.tgt.agg<-((n3.trinuc$tgt*n3.trinuc$agg)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgt.agt<-((n3.trinuc$tgt*n3.trinuc$agt)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs

cps.tgt.ata<-((n3.trinuc$tgt*n3.trinuc$ata)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci
cps.tgt.atc<-((n3.trinuc$tgt*n3.trinuc$atc)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci
cps.tgt.atg<-((n3.trinuc$tgt*n3.trinuc$atg)/(n3.count.aa$c*n3.count.aa$m))/n3.count.diaa$cm
cps.tgt.att<-((n3.trinuc$tgt*n3.trinuc$att)/(n3.count.aa$c*n3.count.aa$i))/n3.count.diaa$ci

cps.tgt.caa<-((n3.trinuc$tgt*n3.trinuc$caa)/(n3.count.aa$c*n3.count.aa$q))/n3.count.diaa$cq
cps.tgt.cac<-((n3.trinuc$tgt*n3.trinuc$cac)/(n3.count.aa$c*n3.count.aa$h))/n3.count.diaa$ch
cps.tgt.cag<-((n3.trinuc$tgt*n3.trinuc$cag)/(n3.count.aa$c*n3.count.aa$q))/n3.count.diaa$cq
cps.tgt.cat<-((n3.trinuc$tgt*n3.trinuc$cat)/(n3.count.aa$c*n3.count.aa$h))/n3.count.diaa$ch

cps.tgt.cca<-((n3.trinuc$tgt*n3.trinuc$cca)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
cps.tgt.ccc<-((n3.trinuc$tgt*n3.trinuc$ccc)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
cps.tgt.ccg<-((n3.trinuc$tgt*n3.trinuc$ccg)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp
cps.tgt.cct<-((n3.trinuc$tgt*n3.trinuc$cct)/(n3.count.aa$c*n3.count.aa$p))/n3.count.diaa$cp

cps.tgt.cga<-((n3.trinuc$tgt*n3.trinuc$cga)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgt.cgc<-((n3.trinuc$tgt*n3.trinuc$cgc)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgt.cgg<-((n3.trinuc$tgt*n3.trinuc$cgg)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr
cps.tgt.cgt<-((n3.trinuc$tgt*n3.trinuc$cgt)/(n3.count.aa$c*n3.count.aa$r))/n3.count.diaa$cr

cps.tgt.cta<-((n3.trinuc$tgt*n3.trinuc$cta)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgt.ctc<-((n3.trinuc$tgt*n3.trinuc$ctc)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgt.ctg<-((n3.trinuc$tgt*n3.trinuc$ctg)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgt.ctt<-((n3.trinuc$tgt*n3.trinuc$ctt)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl

cps.tgt.gaa<-((n3.trinuc$tgt*n3.trinuc$gaa)/(n3.count.aa$c*n3.count.aa$e))/n3.count.diaa$ce
cps.tgt.gac<-((n3.trinuc$tgt*n3.trinuc$gac)/(n3.count.aa$c*n3.count.aa$d))/n3.count.diaa$cd
cps.tgt.gag<-((n3.trinuc$tgt*n3.trinuc$gag)/(n3.count.aa$c*n3.count.aa$e))/n3.count.diaa$ce
cps.tgt.gat<-((n3.trinuc$tgt*n3.trinuc$gat)/(n3.count.aa$c*n3.count.aa$d))/n3.count.diaa$cd

cps.tgt.gca<-((n3.trinuc$tgt*n3.trinuc$gca)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
cps.tgt.gcc<-((n3.trinuc$tgt*n3.trinuc$gcc)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
cps.tgt.gcg<-((n3.trinuc$tgt*n3.trinuc$gcg)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca
cps.tgt.gct<-((n3.trinuc$tgt*n3.trinuc$gct)/(n3.count.aa$c*n3.count.aa$a))/n3.count.diaa$ca

cps.tgt.gga<-((n3.trinuc$tgt*n3.trinuc$gga)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
cps.tgt.ggc<-((n3.trinuc$tgt*n3.trinuc$ggc)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
cps.tgt.ggg<-((n3.trinuc$tgt*n3.trinuc$ggg)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg
cps.tgt.ggt<-((n3.trinuc$tgt*n3.trinuc$ggt)/(n3.count.aa$c*n3.count.aa$g))/n3.count.diaa$cg

cps.tgt.gta<-((n3.trinuc$tgt*n3.trinuc$gta)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
cps.tgt.gtc<-((n3.trinuc$tgt*n3.trinuc$gtc)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
cps.tgt.gtg<-((n3.trinuc$tgt*n3.trinuc$gtg)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv
cps.tgt.gtt<-((n3.trinuc$tgt*n3.trinuc$gtt)/(n3.count.aa$c*n3.count.aa$v))/n3.count.diaa$cv

#Stop codon
#cps.tgt.taa<-((n3.trinuc$tgt*n3.trinuc$taa)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgt.tac<-((n3.trinuc$tgt*n3.trinuc$tac)/(n3.count.aa$c*n3.count.aa$y))/n3.count.diaa$cy
#Stop codon
#cps.tgt.tag<-((n3.trinuc$tgt*n3.trinuc$tag)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgt.tat<-((n3.trinuc$tgt*n3.trinuc$tat)/(n3.count.aa$c*n3.count.aa$y))/n3.count.diaa$cy

cps.tgt.tca<-((n3.trinuc$tgt*n3.trinuc$tca)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
cps.tgt.tcc<-((n3.trinuc$tgt*n3.trinuc$tcc)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
cps.tgt.tcg<-((n3.trinuc$tgt*n3.trinuc$tcg)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs
cps.tgt.tct<-((n3.trinuc$tgt*n3.trinuc$tct)/(n3.count.aa$c*n3.count.aa$s))/n3.count.diaa$cs

#Stop codon
#cps.tgt.tga<-((n3.trinuc$tgt*n3.trinuc$tga)/(n3.count.aa$c*n3.count.aa$k))/n3.count.diaa$ck
cps.tgt.tgc<-((n3.trinuc$tgt*n3.trinuc$tgc)/(n3.count.aa$c*n3.count.aa$c))/n3.count.diaa$cc
cps.tgt.tgg<-((n3.trinuc$tgt*n3.trinuc$tgg)/(n3.count.aa$c*n3.count.aa$w))/n3.count.diaa$cw
cps.tgt.tgt<-((n3.trinuc$tgt*n3.trinuc$tgt)/(n3.count.aa$c*n3.count.aa$c))/n3.count.diaa$cc

cps.tgt.tta<-((n3.trinuc$tgt*n3.trinuc$tta)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgt.ttc<-((n3.trinuc$tgt*n3.trinuc$ttc)/(n3.count.aa$c*n3.count.aa$f))/n3.count.diaa$cf
cps.tgt.ttg<-((n3.trinuc$tgt*n3.trinuc$ttg)/(n3.count.aa$c*n3.count.aa$l))/n3.count.diaa$cl
cps.tgt.ttt<-((n3.trinuc$tgt*n3.trinuc$ttt)/(n3.count.aa$c*n3.count.aa$f))/n3.count.diaa$cf








cps.tta.aaa<-((n3.trinuc$tta*n3.trinuc$aaa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.tta.aac<-((n3.trinuc$tta*n3.trinuc$aac)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln
cps.tta.aag<-((n3.trinuc$tta*n3.trinuc$aag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.tta.aat<-((n3.trinuc$tta*n3.trinuc$aat)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln

cps.tta.aca<-((n3.trinuc$tta*n3.trinuc$aca)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.tta.acc<-((n3.trinuc$tta*n3.trinuc$acc)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.tta.acg<-((n3.trinuc$tta*n3.trinuc$acg)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.tta.act<-((n3.trinuc$tta*n3.trinuc$act)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt

cps.tta.aga<-((n3.trinuc$tta*n3.trinuc$aga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.tta.agc<-((n3.trinuc$tta*n3.trinuc$agc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.tta.agg<-((n3.trinuc$tta*n3.trinuc$agg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.tta.agt<-((n3.trinuc$tta*n3.trinuc$agt)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

cps.tta.ata<-((n3.trinuc$tta*n3.trinuc$ata)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.tta.atc<-((n3.trinuc$tta*n3.trinuc$atc)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.tta.atg<-((n3.trinuc$tta*n3.trinuc$atg)/(n3.count.aa$l*n3.count.aa$m))/n3.count.diaa$lm
cps.tta.att<-((n3.trinuc$tta*n3.trinuc$att)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li

cps.tta.caa<-((n3.trinuc$tta*n3.trinuc$caa)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.tta.cac<-((n3.trinuc$tta*n3.trinuc$cac)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh
cps.tta.cag<-((n3.trinuc$tta*n3.trinuc$cag)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.tta.cat<-((n3.trinuc$tta*n3.trinuc$cat)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh

cps.tta.cca<-((n3.trinuc$tta*n3.trinuc$cca)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.tta.ccc<-((n3.trinuc$tta*n3.trinuc$ccc)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.tta.ccg<-((n3.trinuc$tta*n3.trinuc$ccg)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.tta.cct<-((n3.trinuc$tta*n3.trinuc$cct)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp

cps.tta.cga<-((n3.trinuc$tta*n3.trinuc$cga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.tta.cgc<-((n3.trinuc$tta*n3.trinuc$cgc)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.tta.cgg<-((n3.trinuc$tta*n3.trinuc$cgg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.tta.cgt<-((n3.trinuc$tta*n3.trinuc$cgt)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr

cps.tta.cta<-((n3.trinuc$tta*n3.trinuc$cta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.tta.ctc<-((n3.trinuc$tta*n3.trinuc$ctc)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.tta.ctg<-((n3.trinuc$tta*n3.trinuc$ctg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.tta.ctt<-((n3.trinuc$tta*n3.trinuc$ctt)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll

cps.tta.gaa<-((n3.trinuc$tta*n3.trinuc$gaa)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.tta.gac<-((n3.trinuc$tta*n3.trinuc$gac)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld
cps.tta.gag<-((n3.trinuc$tta*n3.trinuc$gag)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.tta.gat<-((n3.trinuc$tta*n3.trinuc$gat)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld

cps.tta.gca<-((n3.trinuc$tta*n3.trinuc$gca)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.tta.gcc<-((n3.trinuc$tta*n3.trinuc$gcc)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.tta.gcg<-((n3.trinuc$tta*n3.trinuc$gcg)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.tta.gct<-((n3.trinuc$tta*n3.trinuc$gct)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la

cps.tta.gga<-((n3.trinuc$tta*n3.trinuc$gga)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.tta.ggc<-((n3.trinuc$tta*n3.trinuc$ggc)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.tta.ggg<-((n3.trinuc$tta*n3.trinuc$ggg)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.tta.ggt<-((n3.trinuc$tta*n3.trinuc$ggt)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg

cps.tta.gta<-((n3.trinuc$tta*n3.trinuc$gta)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.tta.gtc<-((n3.trinuc$tta*n3.trinuc$gtc)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.tta.gtg<-((n3.trinuc$tta*n3.trinuc$gtg)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.tta.gtt<-((n3.trinuc$tta*n3.trinuc$gtt)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv

#Stop codon
#cps.tta.taa<-((n3.trinuc$tta*n3.trinuc$taa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.tta.tac<-((n3.trinuc$tta*n3.trinuc$tac)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly
#Stop codon
#cps.tta.tag<-((n3.trinuc$tta*n3.trinuc$tag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.tta.tat<-((n3.trinuc$tta*n3.trinuc$tat)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly

cps.tta.tca<-((n3.trinuc$tta*n3.trinuc$tca)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.tta.tcc<-((n3.trinuc$tta*n3.trinuc$tcc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.tta.tcg<-((n3.trinuc$tta*n3.trinuc$tcg)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.tta.tct<-((n3.trinuc$tta*n3.trinuc$tct)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

#Stop codon
#cps.tta.tga<-((n3.trinuc$tta*n3.trinuc$tga)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.tta.tgc<-((n3.trinuc$tta*n3.trinuc$tgc)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc
cps.tta.tgg<-((n3.trinuc$tta*n3.trinuc$tgg)/(n3.count.aa$l*n3.count.aa$w))/n3.count.diaa$lw
cps.tta.tgt<-((n3.trinuc$tta*n3.trinuc$tgt)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc

cps.tta.tta<-((n3.trinuc$tta*n3.trinuc$tta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.tta.ttc<-((n3.trinuc$tta*n3.trinuc$ttc)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf
cps.tta.ttg<-((n3.trinuc$tta*n3.trinuc$ttg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.tta.ttt<-((n3.trinuc$tta*n3.trinuc$ttt)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf








cps.ttc.aaa<-((n3.trinuc$ttc*n3.trinuc$aaa)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttc.aac<-((n3.trinuc$ttc*n3.trinuc$aac)/(n3.count.aa$f*n3.count.aa$n))/n3.count.diaa$fn
cps.ttc.aag<-((n3.trinuc$ttc*n3.trinuc$aag)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttc.aat<-((n3.trinuc$ttc*n3.trinuc$aat)/(n3.count.aa$f*n3.count.aa$n))/n3.count.diaa$fn

cps.ttc.aca<-((n3.trinuc$ttc*n3.trinuc$aca)/(n3.count.aa$f*n3.count.aa$t))/n3.count.diaa$ft
cps.ttc.acc<-((n3.trinuc$ttc*n3.trinuc$acc)/(n3.count.aa$f*n3.count.aa$t))/n3.count.diaa$ft
cps.ttc.acg<-((n3.trinuc$ttc*n3.trinuc$acg)/(n3.count.aa$f*n3.count.aa$t))/n3.count.diaa$ft
cps.ttc.act<-((n3.trinuc$ttc*n3.trinuc$act)/(n3.count.aa$f*n3.count.aa$t))/n3.count.diaa$ft

cps.ttc.aga<-((n3.trinuc$ttc*n3.trinuc$aga)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttc.agc<-((n3.trinuc$ttc*n3.trinuc$agc)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs
cps.ttc.agg<-((n3.trinuc$ttc*n3.trinuc$agg)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttc.agt<-((n3.trinuc$ttc*n3.trinuc$agt)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs

cps.ttc.ata<-((n3.trinuc$ttc*n3.trinuc$ata)/(n3.count.aa$f*n3.count.aa$i))/n3.count.diaa$fi
cps.ttc.atc<-((n3.trinuc$ttc*n3.trinuc$atc)/(n3.count.aa$f*n3.count.aa$i))/n3.count.diaa$fi
cps.ttc.atg<-((n3.trinuc$ttc*n3.trinuc$atg)/(n3.count.aa$f*n3.count.aa$m))/n3.count.diaa$fm
cps.ttc.att<-((n3.trinuc$ttc*n3.trinuc$att)/(n3.count.aa$f*n3.count.aa$i))/n3.count.diaa$fi

cps.ttc.caa<-((n3.trinuc$ttc*n3.trinuc$caa)/(n3.count.aa$f*n3.count.aa$q))/n3.count.diaa$fq
cps.ttc.cac<-((n3.trinuc$ttc*n3.trinuc$cac)/(n3.count.aa$f*n3.count.aa$h))/n3.count.diaa$fh
cps.ttc.cag<-((n3.trinuc$ttc*n3.trinuc$cag)/(n3.count.aa$f*n3.count.aa$q))/n3.count.diaa$fq
cps.ttc.cat<-((n3.trinuc$ttc*n3.trinuc$cat)/(n3.count.aa$f*n3.count.aa$h))/n3.count.diaa$fh

cps.ttc.cca<-((n3.trinuc$ttc*n3.trinuc$cca)/(n3.count.aa$f*n3.count.aa$p))/n3.count.diaa$fp
cps.ttc.ccc<-((n3.trinuc$ttc*n3.trinuc$ccc)/(n3.count.aa$f*n3.count.aa$p))/n3.count.diaa$fp
cps.ttc.ccg<-((n3.trinuc$ttc*n3.trinuc$ccg)/(n3.count.aa$f*n3.count.aa$p))/n3.count.diaa$fp
cps.ttc.cct<-((n3.trinuc$ttc*n3.trinuc$cct)/(n3.count.aa$f*n3.count.aa$p))/n3.count.diaa$fp

cps.ttc.cga<-((n3.trinuc$ttc*n3.trinuc$cga)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttc.cgc<-((n3.trinuc$ttc*n3.trinuc$cgc)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttc.cgg<-((n3.trinuc$ttc*n3.trinuc$cgg)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttc.cgt<-((n3.trinuc$ttc*n3.trinuc$cgt)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr

cps.ttc.cta<-((n3.trinuc$ttc*n3.trinuc$cta)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttc.ctc<-((n3.trinuc$ttc*n3.trinuc$ctc)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttc.ctg<-((n3.trinuc$ttc*n3.trinuc$ctg)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttc.ctt<-((n3.trinuc$ttc*n3.trinuc$ctt)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl

cps.ttc.gaa<-((n3.trinuc$ttc*n3.trinuc$gaa)/(n3.count.aa$f*n3.count.aa$e))/n3.count.diaa$fe
cps.ttc.gac<-((n3.trinuc$ttc*n3.trinuc$gac)/(n3.count.aa$f*n3.count.aa$d))/n3.count.diaa$fd
cps.ttc.gag<-((n3.trinuc$ttc*n3.trinuc$gag)/(n3.count.aa$f*n3.count.aa$e))/n3.count.diaa$fe
cps.ttc.gat<-((n3.trinuc$ttc*n3.trinuc$gat)/(n3.count.aa$f*n3.count.aa$d))/n3.count.diaa$fd

cps.ttc.gca<-((n3.trinuc$ttc*n3.trinuc$gca)/(n3.count.aa$f*n3.count.aa$a))/n3.count.diaa$fa
cps.ttc.gcc<-((n3.trinuc$ttc*n3.trinuc$gcc)/(n3.count.aa$f*n3.count.aa$a))/n3.count.diaa$fa
cps.ttc.gcg<-((n3.trinuc$ttc*n3.trinuc$gcg)/(n3.count.aa$f*n3.count.aa$a))/n3.count.diaa$fa
cps.ttc.gct<-((n3.trinuc$ttc*n3.trinuc$gct)/(n3.count.aa$f*n3.count.aa$a))/n3.count.diaa$fa

cps.ttc.gga<-((n3.trinuc$ttc*n3.trinuc$gga)/(n3.count.aa$f*n3.count.aa$g))/n3.count.diaa$fg
cps.ttc.ggc<-((n3.trinuc$ttc*n3.trinuc$ggc)/(n3.count.aa$f*n3.count.aa$g))/n3.count.diaa$fg
cps.ttc.ggg<-((n3.trinuc$ttc*n3.trinuc$ggg)/(n3.count.aa$f*n3.count.aa$g))/n3.count.diaa$fg
cps.ttc.ggt<-((n3.trinuc$ttc*n3.trinuc$ggt)/(n3.count.aa$f*n3.count.aa$g))/n3.count.diaa$fg

cps.ttc.gta<-((n3.trinuc$ttc*n3.trinuc$gta)/(n3.count.aa$f*n3.count.aa$v))/n3.count.diaa$fv
cps.ttc.gtc<-((n3.trinuc$ttc*n3.trinuc$gtc)/(n3.count.aa$f*n3.count.aa$v))/n3.count.diaa$fv
cps.ttc.gtg<-((n3.trinuc$ttc*n3.trinuc$gtg)/(n3.count.aa$f*n3.count.aa$v))/n3.count.diaa$fv
cps.ttc.gtt<-((n3.trinuc$ttc*n3.trinuc$gtt)/(n3.count.aa$f*n3.count.aa$v))/n3.count.diaa$fv

#Stop codon
#cps.ttc.taa<-((n3.trinuc$ttc*n3.trinuc$taa)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttc.tac<-((n3.trinuc$ttc*n3.trinuc$tac)/(n3.count.aa$f*n3.count.aa$y))/n3.count.diaa$fy
#Stop codon
#cps.ttc.tag<-((n3.trinuc$ttc*n3.trinuc$tag)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttc.tat<-((n3.trinuc$ttc*n3.trinuc$tat)/(n3.count.aa$f*n3.count.aa$y))/n3.count.diaa$fy

cps.ttc.tca<-((n3.trinuc$ttc*n3.trinuc$tca)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs
cps.ttc.tcc<-((n3.trinuc$ttc*n3.trinuc$tcc)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs
cps.ttc.tcg<-((n3.trinuc$ttc*n3.trinuc$tcg)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs
cps.ttc.tct<-((n3.trinuc$ttc*n3.trinuc$tct)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs

#Stop codon
#cps.ttc.tga<-((n3.trinuc$ttc*n3.trinuc$tga)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttc.tgc<-((n3.trinuc$ttc*n3.trinuc$tgc)/(n3.count.aa$f*n3.count.aa$c))/n3.count.diaa$fc
cps.ttc.tgg<-((n3.trinuc$ttc*n3.trinuc$tgg)/(n3.count.aa$f*n3.count.aa$w))/n3.count.diaa$fw
cps.ttc.tgt<-((n3.trinuc$ttc*n3.trinuc$tgt)/(n3.count.aa$f*n3.count.aa$c))/n3.count.diaa$fc

cps.ttc.tta<-((n3.trinuc$ttc*n3.trinuc$tta)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttc.ttc<-((n3.trinuc$ttc*n3.trinuc$ttc)/(n3.count.aa$f*n3.count.aa$f))/n3.count.diaa$ff
cps.ttc.ttg<-((n3.trinuc$ttc*n3.trinuc$ttg)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttc.ttt<-((n3.trinuc$ttc*n3.trinuc$ttt)/(n3.count.aa$f*n3.count.aa$f))/n3.count.diaa$ff








cps.ttg.aaa<-((n3.trinuc$ttg*n3.trinuc$aaa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ttg.aac<-((n3.trinuc$ttg*n3.trinuc$aac)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln
cps.ttg.aag<-((n3.trinuc$ttg*n3.trinuc$aag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ttg.aat<-((n3.trinuc$ttg*n3.trinuc$aat)/(n3.count.aa$l*n3.count.aa$n))/n3.count.diaa$ln

cps.ttg.aca<-((n3.trinuc$ttg*n3.trinuc$aca)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ttg.acc<-((n3.trinuc$ttg*n3.trinuc$acc)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ttg.acg<-((n3.trinuc$ttg*n3.trinuc$acg)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt
cps.ttg.act<-((n3.trinuc$ttg*n3.trinuc$act)/(n3.count.aa$l*n3.count.aa$t))/n3.count.diaa$lt

cps.ttg.aga<-((n3.trinuc$ttg*n3.trinuc$aga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ttg.agc<-((n3.trinuc$ttg*n3.trinuc$agc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ttg.agg<-((n3.trinuc$ttg*n3.trinuc$agg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ttg.agt<-((n3.trinuc$ttg*n3.trinuc$agt)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

cps.ttg.ata<-((n3.trinuc$ttg*n3.trinuc$ata)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.ttg.atc<-((n3.trinuc$ttg*n3.trinuc$atc)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li
cps.ttg.atg<-((n3.trinuc$ttg*n3.trinuc$atg)/(n3.count.aa$l*n3.count.aa$m))/n3.count.diaa$lm
cps.ttg.att<-((n3.trinuc$ttg*n3.trinuc$att)/(n3.count.aa$l*n3.count.aa$i))/n3.count.diaa$li

cps.ttg.caa<-((n3.trinuc$ttg*n3.trinuc$caa)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.ttg.cac<-((n3.trinuc$ttg*n3.trinuc$cac)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh
cps.ttg.cag<-((n3.trinuc$ttg*n3.trinuc$cag)/(n3.count.aa$l*n3.count.aa$q))/n3.count.diaa$lq
cps.ttg.cat<-((n3.trinuc$ttg*n3.trinuc$cat)/(n3.count.aa$l*n3.count.aa$h))/n3.count.diaa$lh

cps.ttg.cca<-((n3.trinuc$ttg*n3.trinuc$cca)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ttg.ccc<-((n3.trinuc$ttg*n3.trinuc$ccc)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ttg.ccg<-((n3.trinuc$ttg*n3.trinuc$ccg)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp
cps.ttg.cct<-((n3.trinuc$ttg*n3.trinuc$cct)/(n3.count.aa$l*n3.count.aa$p))/n3.count.diaa$lp

cps.ttg.cga<-((n3.trinuc$ttg*n3.trinuc$cga)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ttg.cgc<-((n3.trinuc$ttg*n3.trinuc$cgc)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ttg.cgg<-((n3.trinuc$ttg*n3.trinuc$cgg)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr
cps.ttg.cgt<-((n3.trinuc$ttg*n3.trinuc$cgt)/(n3.count.aa$l*n3.count.aa$r))/n3.count.diaa$lr

cps.ttg.cta<-((n3.trinuc$ttg*n3.trinuc$cta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ttg.ctc<-((n3.trinuc$ttg*n3.trinuc$ctc)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ttg.ctg<-((n3.trinuc$ttg*n3.trinuc$ctg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ttg.ctt<-((n3.trinuc$ttg*n3.trinuc$ctt)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll

cps.ttg.gaa<-((n3.trinuc$ttg*n3.trinuc$gaa)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.ttg.gac<-((n3.trinuc$ttg*n3.trinuc$gac)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld
cps.ttg.gag<-((n3.trinuc$ttg*n3.trinuc$gag)/(n3.count.aa$l*n3.count.aa$e))/n3.count.diaa$le
cps.ttg.gat<-((n3.trinuc$ttg*n3.trinuc$gat)/(n3.count.aa$l*n3.count.aa$d))/n3.count.diaa$ld

cps.ttg.gca<-((n3.trinuc$ttg*n3.trinuc$gca)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ttg.gcc<-((n3.trinuc$ttg*n3.trinuc$gcc)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ttg.gcg<-((n3.trinuc$ttg*n3.trinuc$gcg)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la
cps.ttg.gct<-((n3.trinuc$ttg*n3.trinuc$gct)/(n3.count.aa$l*n3.count.aa$a))/n3.count.diaa$la

cps.ttg.gga<-((n3.trinuc$ttg*n3.trinuc$gga)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ttg.ggc<-((n3.trinuc$ttg*n3.trinuc$ggc)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ttg.ggg<-((n3.trinuc$ttg*n3.trinuc$ggg)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg
cps.ttg.ggt<-((n3.trinuc$ttg*n3.trinuc$ggt)/(n3.count.aa$l*n3.count.aa$g))/n3.count.diaa$lg

cps.ttg.gta<-((n3.trinuc$ttg*n3.trinuc$gta)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ttg.gtc<-((n3.trinuc$ttg*n3.trinuc$gtc)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ttg.gtg<-((n3.trinuc$ttg*n3.trinuc$gtg)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv
cps.ttg.gtt<-((n3.trinuc$ttg*n3.trinuc$gtt)/(n3.count.aa$l*n3.count.aa$v))/n3.count.diaa$lv

#Stop codon
#cps.ttg.taa<-((n3.trinuc$ttg*n3.trinuc$taa)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ttg.tac<-((n3.trinuc$ttg*n3.trinuc$tac)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly
#Stop codon
#cps.ttg.tag<-((n3.trinuc$ttg*n3.trinuc$tag)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ttg.tat<-((n3.trinuc$ttg*n3.trinuc$tat)/(n3.count.aa$l*n3.count.aa$y))/n3.count.diaa$ly

cps.ttg.tca<-((n3.trinuc$ttg*n3.trinuc$tca)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ttg.tcc<-((n3.trinuc$ttg*n3.trinuc$tcc)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ttg.tcg<-((n3.trinuc$ttg*n3.trinuc$tcg)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls
cps.ttg.tct<-((n3.trinuc$ttg*n3.trinuc$tct)/(n3.count.aa$l*n3.count.aa$s))/n3.count.diaa$ls

#Stop codon
#cps.ttg.tga<-((n3.trinuc$ttg*n3.trinuc$tga)/(n3.count.aa$l*n3.count.aa$k))/n3.count.diaa$lk
cps.ttg.tgc<-((n3.trinuc$ttg*n3.trinuc$tgc)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc
cps.ttg.tgg<-((n3.trinuc$ttg*n3.trinuc$tgg)/(n3.count.aa$l*n3.count.aa$w))/n3.count.diaa$lw
cps.ttg.tgt<-((n3.trinuc$ttg*n3.trinuc$tgt)/(n3.count.aa$l*n3.count.aa$c))/n3.count.diaa$lc

cps.ttg.tta<-((n3.trinuc$ttg*n3.trinuc$tta)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ttg.ttc<-((n3.trinuc$ttg*n3.trinuc$ttc)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf
cps.ttg.ttg<-((n3.trinuc$ttg*n3.trinuc$ttg)/(n3.count.aa$l*n3.count.aa$l))/n3.count.diaa$ll
cps.ttg.ttt<-((n3.trinuc$ttg*n3.trinuc$ttt)/(n3.count.aa$l*n3.count.aa$f))/n3.count.diaa$lf








cps.ttt.aaa<-((n3.trinuc$ttt*n3.trinuc$aaa)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttt.aac<-((n3.trinuc$ttt*n3.trinuc$aac)/(n3.count.aa$f*n3.count.aa$n))/n3.count.diaa$fn
cps.ttt.aag<-((n3.trinuc$ttt*n3.trinuc$aag)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttt.aat<-((n3.trinuc$ttt*n3.trinuc$aat)/(n3.count.aa$f*n3.count.aa$n))/n3.count.diaa$fn

cps.ttt.aca<-((n3.trinuc$ttt*n3.trinuc$aca)/(n3.count.aa$f*n3.count.aa$t))/n3.count.diaa$ft
cps.ttt.acc<-((n3.trinuc$ttt*n3.trinuc$acc)/(n3.count.aa$f*n3.count.aa$t))/n3.count.diaa$ft
cps.ttt.acg<-((n3.trinuc$ttt*n3.trinuc$acg)/(n3.count.aa$f*n3.count.aa$t))/n3.count.diaa$ft
cps.ttt.act<-((n3.trinuc$ttt*n3.trinuc$act)/(n3.count.aa$f*n3.count.aa$t))/n3.count.diaa$ft

cps.ttt.aga<-((n3.trinuc$ttt*n3.trinuc$aga)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttt.agc<-((n3.trinuc$ttt*n3.trinuc$agc)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs
cps.ttt.agg<-((n3.trinuc$ttt*n3.trinuc$agg)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttt.agt<-((n3.trinuc$ttt*n3.trinuc$agt)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs

cps.ttt.ata<-((n3.trinuc$ttt*n3.trinuc$ata)/(n3.count.aa$f*n3.count.aa$i))/n3.count.diaa$fi
cps.ttt.atc<-((n3.trinuc$ttt*n3.trinuc$atc)/(n3.count.aa$f*n3.count.aa$i))/n3.count.diaa$fi
cps.ttt.atg<-((n3.trinuc$ttt*n3.trinuc$atg)/(n3.count.aa$f*n3.count.aa$m))/n3.count.diaa$fm
cps.ttt.att<-((n3.trinuc$ttt*n3.trinuc$att)/(n3.count.aa$f*n3.count.aa$i))/n3.count.diaa$fi

cps.ttt.caa<-((n3.trinuc$ttt*n3.trinuc$caa)/(n3.count.aa$f*n3.count.aa$q))/n3.count.diaa$fq
cps.ttt.cac<-((n3.trinuc$ttt*n3.trinuc$cac)/(n3.count.aa$f*n3.count.aa$h))/n3.count.diaa$fh
cps.ttt.cag<-((n3.trinuc$ttt*n3.trinuc$cag)/(n3.count.aa$f*n3.count.aa$q))/n3.count.diaa$fq
cps.ttt.cat<-((n3.trinuc$ttt*n3.trinuc$cat)/(n3.count.aa$f*n3.count.aa$h))/n3.count.diaa$fh

cps.ttt.cca<-((n3.trinuc$ttt*n3.trinuc$cca)/(n3.count.aa$f*n3.count.aa$p))/n3.count.diaa$fp
cps.ttt.ccc<-((n3.trinuc$ttt*n3.trinuc$ccc)/(n3.count.aa$f*n3.count.aa$p))/n3.count.diaa$fp
cps.ttt.ccg<-((n3.trinuc$ttt*n3.trinuc$ccg)/(n3.count.aa$f*n3.count.aa$p))/n3.count.diaa$fp
cps.ttt.cct<-((n3.trinuc$ttt*n3.trinuc$cct)/(n3.count.aa$f*n3.count.aa$p))/n3.count.diaa$fp

cps.ttt.cga<-((n3.trinuc$ttt*n3.trinuc$cga)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttt.cgc<-((n3.trinuc$ttt*n3.trinuc$cgc)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttt.cgg<-((n3.trinuc$ttt*n3.trinuc$cgg)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr
cps.ttt.cgt<-((n3.trinuc$ttt*n3.trinuc$cgt)/(n3.count.aa$f*n3.count.aa$r))/n3.count.diaa$fr

cps.ttt.cta<-((n3.trinuc$ttt*n3.trinuc$cta)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttt.ctc<-((n3.trinuc$ttt*n3.trinuc$ctc)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttt.ctg<-((n3.trinuc$ttt*n3.trinuc$ctg)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttt.ctt<-((n3.trinuc$ttt*n3.trinuc$ctt)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl

cps.ttt.gaa<-((n3.trinuc$ttt*n3.trinuc$gaa)/(n3.count.aa$f*n3.count.aa$e))/n3.count.diaa$fe
cps.ttt.gac<-((n3.trinuc$ttt*n3.trinuc$gac)/(n3.count.aa$f*n3.count.aa$d))/n3.count.diaa$fd
cps.ttt.gag<-((n3.trinuc$ttt*n3.trinuc$gag)/(n3.count.aa$f*n3.count.aa$e))/n3.count.diaa$fe
cps.ttt.gat<-((n3.trinuc$ttt*n3.trinuc$gat)/(n3.count.aa$f*n3.count.aa$d))/n3.count.diaa$fd

cps.ttt.gca<-((n3.trinuc$ttt*n3.trinuc$gca)/(n3.count.aa$f*n3.count.aa$a))/n3.count.diaa$fa
cps.ttt.gcc<-((n3.trinuc$ttt*n3.trinuc$gcc)/(n3.count.aa$f*n3.count.aa$a))/n3.count.diaa$fa
cps.ttt.gcg<-((n3.trinuc$ttt*n3.trinuc$gcg)/(n3.count.aa$f*n3.count.aa$a))/n3.count.diaa$fa
cps.ttt.gct<-((n3.trinuc$ttt*n3.trinuc$gct)/(n3.count.aa$f*n3.count.aa$a))/n3.count.diaa$fa

cps.ttt.gga<-((n3.trinuc$ttt*n3.trinuc$gga)/(n3.count.aa$f*n3.count.aa$g))/n3.count.diaa$fg
cps.ttt.ggc<-((n3.trinuc$ttt*n3.trinuc$ggc)/(n3.count.aa$f*n3.count.aa$g))/n3.count.diaa$fg
cps.ttt.ggg<-((n3.trinuc$ttt*n3.trinuc$ggg)/(n3.count.aa$f*n3.count.aa$g))/n3.count.diaa$fg
cps.ttt.ggt<-((n3.trinuc$ttt*n3.trinuc$ggt)/(n3.count.aa$f*n3.count.aa$g))/n3.count.diaa$fg

cps.ttt.gta<-((n3.trinuc$ttt*n3.trinuc$gta)/(n3.count.aa$f*n3.count.aa$v))/n3.count.diaa$fv
cps.ttt.gtc<-((n3.trinuc$ttt*n3.trinuc$gtc)/(n3.count.aa$f*n3.count.aa$v))/n3.count.diaa$fv
cps.ttt.gtg<-((n3.trinuc$ttt*n3.trinuc$gtg)/(n3.count.aa$f*n3.count.aa$v))/n3.count.diaa$fv
cps.ttt.gtt<-((n3.trinuc$ttt*n3.trinuc$gtt)/(n3.count.aa$f*n3.count.aa$v))/n3.count.diaa$fv

#Stop codon
#cps.ttt.taa<-((n3.trinuc$ttt*n3.trinuc$taa)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttt.tac<-((n3.trinuc$ttt*n3.trinuc$tac)/(n3.count.aa$f*n3.count.aa$y))/n3.count.diaa$fy
#Stop codon
#cps.ttt.tag<-((n3.trinuc$ttt*n3.trinuc$tag)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttt.tat<-((n3.trinuc$ttt*n3.trinuc$tat)/(n3.count.aa$f*n3.count.aa$y))/n3.count.diaa$fy

cps.ttt.tca<-((n3.trinuc$ttt*n3.trinuc$tca)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs
cps.ttt.tcc<-((n3.trinuc$ttt*n3.trinuc$tcc)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs
cps.ttt.tcg<-((n3.trinuc$ttt*n3.trinuc$tcg)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs
cps.ttt.tct<-((n3.trinuc$ttt*n3.trinuc$tct)/(n3.count.aa$f*n3.count.aa$s))/n3.count.diaa$fs

#Stop codon
#cps.ttt.tga<-((n3.trinuc$ttt*n3.trinuc$tga)/(n3.count.aa$f*n3.count.aa$k))/n3.count.diaa$fk
cps.ttt.tgc<-((n3.trinuc$ttt*n3.trinuc$tgc)/(n3.count.aa$f*n3.count.aa$c))/n3.count.diaa$fc
cps.ttt.tgg<-((n3.trinuc$ttt*n3.trinuc$tgg)/(n3.count.aa$f*n3.count.aa$w))/n3.count.diaa$fw
cps.ttt.tgt<-((n3.trinuc$ttt*n3.trinuc$tgt)/(n3.count.aa$f*n3.count.aa$c))/n3.count.diaa$fc

cps.ttt.tta<-((n3.trinuc$ttt*n3.trinuc$tta)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttt.ttc<-((n3.trinuc$ttt*n3.trinuc$ttc)/(n3.count.aa$f*n3.count.aa$f))/n3.count.diaa$ff
cps.ttt.ttg<-((n3.trinuc$ttt*n3.trinuc$ttg)/(n3.count.aa$f*n3.count.aa$l))/n3.count.diaa$fl
cps.ttt.ttt<-((n3.trinuc$ttt*n3.trinuc$ttt)/(n3.count.aa$f*n3.count.aa$f))/n3.count.diaa$ff


























































#CPS and CPS observed
#cps.observed/cps.expected
cps.obs.exp.aaaaaa<-log(n3.codon.pair$aaaaaa/cps.aaa.aaa)
cps.obs.exp.aaaaac<-log(n3.codon.pair$aaaaac/cps.aaa.aac)
cps.obs.exp.aaaaag<-log(n3.codon.pair$aaaaag/cps.aaa.aag)
cps.obs.exp.aaaaat<-log(n3.codon.pair$aaaaat/cps.aaa.aat)
cps.obs.exp.aaaaca<-log(n3.codon.pair$aaaaca/cps.aaa.aca)
cps.obs.exp.aaaacc<-log(n3.codon.pair$aaaacc/cps.aaa.acc)
cps.obs.exp.aaaacg<-log(n3.codon.pair$aaaacg/cps.aaa.acg)
cps.obs.exp.aaaact<-log(n3.codon.pair$aaaact/cps.aaa.act)
cps.obs.exp.aaaaga<-log(n3.codon.pair$aaaaga/cps.aaa.aga)
cps.obs.exp.aaaagc<-log(n3.codon.pair$aaaagc/cps.aaa.agc)
cps.obs.exp.aaaagg<-log(n3.codon.pair$aaaagg/cps.aaa.agg)
cps.obs.exp.aaaagt<-log(n3.codon.pair$aaaagt/cps.aaa.agt)
cps.obs.exp.aaaata<-log(n3.codon.pair$aaaata/cps.aaa.ata)
cps.obs.exp.aaaatc<-log(n3.codon.pair$aaaatc/cps.aaa.atc)
cps.obs.exp.aaaatg<-log(n3.codon.pair$aaaatg/cps.aaa.atg)
cps.obs.exp.aaaatt<-log(n3.codon.pair$aaaatt/cps.aaa.att)

cps.obs.exp.aaacaa<-log(n3.codon.pair$aaacaa/cps.aaa.caa)
cps.obs.exp.aaacac<-log(n3.codon.pair$aaacac/cps.aaa.cac)
cps.obs.exp.aaacag<-log(n3.codon.pair$aaacag/cps.aaa.cag)
cps.obs.exp.aaacat<-log(n3.codon.pair$aaacat/cps.aaa.cat)
cps.obs.exp.aaacca<-log(n3.codon.pair$aaacca/cps.aaa.cca)
cps.obs.exp.aaaccc<-log(n3.codon.pair$aaaccc/cps.aaa.ccc)
cps.obs.exp.aaaccg<-log(n3.codon.pair$aaaccg/cps.aaa.ccg)
cps.obs.exp.aaacct<-log(n3.codon.pair$aaacct/cps.aaa.cct)
cps.obs.exp.aaacga<-log(n3.codon.pair$aaacga/cps.aaa.cga)
cps.obs.exp.aaacgc<-log(n3.codon.pair$aaacgc/cps.aaa.cgc)
cps.obs.exp.aaacgg<-log(n3.codon.pair$aaacgg/cps.aaa.cgg)
cps.obs.exp.aaacgt<-log(n3.codon.pair$aaacgt/cps.aaa.cgt)
cps.obs.exp.aaacta<-log(n3.codon.pair$aaacta/cps.aaa.cta)
cps.obs.exp.aaactc<-log(n3.codon.pair$aaactc/cps.aaa.ctc)
cps.obs.exp.aaactg<-log(n3.codon.pair$aaactg/cps.aaa.ctg)
cps.obs.exp.aaactt<-log(n3.codon.pair$aaactt/cps.aaa.ctt)

cps.obs.exp.aaagaa<-log(n3.codon.pair$aaagaa/cps.aaa.gaa)
cps.obs.exp.aaagac<-log(n3.codon.pair$aaagac/cps.aaa.gac)
cps.obs.exp.aaagag<-log(n3.codon.pair$aaagag/cps.aaa.gag)
cps.obs.exp.aaagat<-log(n3.codon.pair$aaagat/cps.aaa.gat)
cps.obs.exp.aaagca<-log(n3.codon.pair$aaagca/cps.aaa.gca)
cps.obs.exp.aaagcc<-log(n3.codon.pair$aaagcc/cps.aaa.gcc)
cps.obs.exp.aaagcg<-log(n3.codon.pair$aaagcg/cps.aaa.gcg)
cps.obs.exp.aaagct<-log(n3.codon.pair$aaagct/cps.aaa.gct)
cps.obs.exp.aaagga<-log(n3.codon.pair$aaagga/cps.aaa.gga)
cps.obs.exp.aaaggc<-log(n3.codon.pair$aaaggc/cps.aaa.ggc)
cps.obs.exp.aaaggg<-log(n3.codon.pair$aaaggg/cps.aaa.ggg)
cps.obs.exp.aaaggt<-log(n3.codon.pair$aaaggt/cps.aaa.ggt)
cps.obs.exp.aaagta<-log(n3.codon.pair$aaagta/cps.aaa.gta)
cps.obs.exp.aaagtc<-log(n3.codon.pair$aaagtc/cps.aaa.gtc)
cps.obs.exp.aaagtg<-log(n3.codon.pair$aaagtg/cps.aaa.gtg)
cps.obs.exp.aaagtt<-log(n3.codon.pair$aaagtt/cps.aaa.gtt)

#cps.obs.exp.aaataa<-log(n3.codon.pair$aaataa/cps.aaa.taa)
cps.obs.exp.aaatac<-log(n3.codon.pair$aaatac/cps.aaa.tac)
#cps.obs.exp.aaatag<-log(n3.codon.pair$aaatag/cps.aaa.tag)
cps.obs.exp.aaatat<-log(n3.codon.pair$aaatat/cps.aaa.tat)
cps.obs.exp.aaatca<-log(n3.codon.pair$aaatca/cps.aaa.tca)
cps.obs.exp.aaatcc<-log(n3.codon.pair$aaatcc/cps.aaa.tcc)
cps.obs.exp.aaatcg<-log(n3.codon.pair$aaatcg/cps.aaa.tcg)
cps.obs.exp.aaatct<-log(n3.codon.pair$aaatct/cps.aaa.tct)
#cps.obs.exp.aaatga<-log(n3.codon.pair$aaatga/cps.aaa.tga)
cps.obs.exp.aaatgc<-log(n3.codon.pair$aaatgc/cps.aaa.tgc)
cps.obs.exp.aaatgg<-log(n3.codon.pair$aaatgg/cps.aaa.tgg)
cps.obs.exp.aaatgt<-log(n3.codon.pair$aaatgt/cps.aaa.tgt)
cps.obs.exp.aaatta<-log(n3.codon.pair$aaatta/cps.aaa.tta)
cps.obs.exp.aaattc<-log(n3.codon.pair$aaattc/cps.aaa.ttc)
cps.obs.exp.aaattg<-log(n3.codon.pair$aaattg/cps.aaa.ttg)
cps.obs.exp.aaattt<-log(n3.codon.pair$aaattt/cps.aaa.ttt)




cps.obs.exp.aacaaa<-log(n3.codon.pair$aacaaa/cps.aac.aaa)
cps.obs.exp.aacaac<-log(n3.codon.pair$aacaac/cps.aac.aac)
cps.obs.exp.aacaag<-log(n3.codon.pair$aacaag/cps.aac.aag)
cps.obs.exp.aacaat<-log(n3.codon.pair$aacaat/cps.aac.aat)
cps.obs.exp.aacaca<-log(n3.codon.pair$aacaca/cps.aac.aca)
cps.obs.exp.aacacc<-log(n3.codon.pair$aacacc/cps.aac.acc)
cps.obs.exp.aacacg<-log(n3.codon.pair$aacacg/cps.aac.acg)
cps.obs.exp.aacact<-log(n3.codon.pair$aacact/cps.aac.act)
cps.obs.exp.aacaga<-log(n3.codon.pair$aacaga/cps.aac.aga)
cps.obs.exp.aacagc<-log(n3.codon.pair$aacagc/cps.aac.agc)
cps.obs.exp.aacagg<-log(n3.codon.pair$aacagg/cps.aac.agg)
cps.obs.exp.aacagt<-log(n3.codon.pair$aacagt/cps.aac.agt)
cps.obs.exp.aacata<-log(n3.codon.pair$aacata/cps.aac.ata)
cps.obs.exp.aacatc<-log(n3.codon.pair$aacatc/cps.aac.atc)
cps.obs.exp.aacatg<-log(n3.codon.pair$aacatg/cps.aac.atg)
cps.obs.exp.aacatt<-log(n3.codon.pair$aacatt/cps.aac.att)

cps.obs.exp.aaccaa<-log(n3.codon.pair$aaccaa/cps.aac.caa)
cps.obs.exp.aaccac<-log(n3.codon.pair$aaccac/cps.aac.cac)
cps.obs.exp.aaccag<-log(n3.codon.pair$aaccag/cps.aac.cag)
cps.obs.exp.aaccat<-log(n3.codon.pair$aaccat/cps.aac.cat)
cps.obs.exp.aaccca<-log(n3.codon.pair$aaccca/cps.aac.cca)
cps.obs.exp.aacccc<-log(n3.codon.pair$aacccc/cps.aac.ccc)
cps.obs.exp.aacccg<-log(n3.codon.pair$aacccg/cps.aac.ccg)
cps.obs.exp.aaccct<-log(n3.codon.pair$aaccct/cps.aac.cct)
cps.obs.exp.aaccga<-log(n3.codon.pair$aaccga/cps.aac.cga)
cps.obs.exp.aaccgc<-log(n3.codon.pair$aaccgc/cps.aac.cgc)
cps.obs.exp.aaccgg<-log(n3.codon.pair$aaccgg/cps.aac.cgg)
cps.obs.exp.aaccgt<-log(n3.codon.pair$aaccgt/cps.aac.cgt)
cps.obs.exp.aaccta<-log(n3.codon.pair$aaccta/cps.aac.cta)
cps.obs.exp.aacctc<-log(n3.codon.pair$aacctc/cps.aac.ctc)
cps.obs.exp.aacctg<-log(n3.codon.pair$aacctg/cps.aac.ctg)
cps.obs.exp.aacctt<-log(n3.codon.pair$aacctt/cps.aac.ctt)

cps.obs.exp.aacgaa<-log(n3.codon.pair$aacgaa/cps.aac.gaa)
cps.obs.exp.aacgac<-log(n3.codon.pair$aacgac/cps.aac.gac)
cps.obs.exp.aacgag<-log(n3.codon.pair$aacgag/cps.aac.gag)
cps.obs.exp.aacgat<-log(n3.codon.pair$aacgat/cps.aac.gat)
cps.obs.exp.aacgca<-log(n3.codon.pair$aacgca/cps.aac.gca)
cps.obs.exp.aacgcc<-log(n3.codon.pair$aacgcc/cps.aac.gcc)
cps.obs.exp.aacgcg<-log(n3.codon.pair$aacgcg/cps.aac.gcg)
cps.obs.exp.aacgct<-log(n3.codon.pair$aacgct/cps.aac.gct)
cps.obs.exp.aacgga<-log(n3.codon.pair$aacgga/cps.aac.gga)
cps.obs.exp.aacggc<-log(n3.codon.pair$aacggc/cps.aac.ggc)
cps.obs.exp.aacggg<-log(n3.codon.pair$aacggg/cps.aac.ggg)
cps.obs.exp.aacggt<-log(n3.codon.pair$aacggt/cps.aac.ggt)
cps.obs.exp.aacgta<-log(n3.codon.pair$aacgta/cps.aac.gta)
cps.obs.exp.aacgtc<-log(n3.codon.pair$aacgtc/cps.aac.gtc)
cps.obs.exp.aacgtg<-log(n3.codon.pair$aacgtg/cps.aac.gtg)
cps.obs.exp.aacgtt<-log(n3.codon.pair$aacgtt/cps.aac.gtt)

#cps.obs.exp.aactaa<-log(n3.codon.pair$aactaa/cps.aac.taa)
cps.obs.exp.aactac<-log(n3.codon.pair$aactac/cps.aac.tac)
#cps.obs.exp.aactag<-log(n3.codon.pair$aactag/cps.aac.tag)
cps.obs.exp.aactat<-log(n3.codon.pair$aactat/cps.aac.tat)
cps.obs.exp.aactca<-log(n3.codon.pair$aactca/cps.aac.tca)
cps.obs.exp.aactcc<-log(n3.codon.pair$aactcc/cps.aac.tcc)
cps.obs.exp.aactcg<-log(n3.codon.pair$aactcg/cps.aac.tcg)
cps.obs.exp.aactct<-log(n3.codon.pair$aactct/cps.aac.tct)
#cps.obs.exp.aactga<-log(n3.codon.pair$aactga/cps.aac.tga)
cps.obs.exp.aactgc<-log(n3.codon.pair$aactgc/cps.aac.tgc)
cps.obs.exp.aactgg<-log(n3.codon.pair$aactgg/cps.aac.tgg)
cps.obs.exp.aactgt<-log(n3.codon.pair$aactgt/cps.aac.tgt)
cps.obs.exp.aactta<-log(n3.codon.pair$aactta/cps.aac.tta)
cps.obs.exp.aacttc<-log(n3.codon.pair$aacttc/cps.aac.ttc)
cps.obs.exp.aacttg<-log(n3.codon.pair$aacttg/cps.aac.ttg)
cps.obs.exp.aacttt<-log(n3.codon.pair$aacttt/cps.aac.ttt)








cps.obs.exp.aagaaa<-log(n3.codon.pair$aagaaa/cps.aag.aaa)
cps.obs.exp.aagaac<-log(n3.codon.pair$aagaac/cps.aag.aac)
cps.obs.exp.aagaag<-log(n3.codon.pair$aagaag/cps.aag.aag)
cps.obs.exp.aagaat<-log(n3.codon.pair$aagaat/cps.aag.aat)
cps.obs.exp.aagaca<-log(n3.codon.pair$aagaca/cps.aag.aca)
cps.obs.exp.aagacc<-log(n3.codon.pair$aagacc/cps.aag.acc)
cps.obs.exp.aagacg<-log(n3.codon.pair$aagacg/cps.aag.acg)
cps.obs.exp.aagact<-log(n3.codon.pair$aagact/cps.aag.act)
cps.obs.exp.aagaga<-log(n3.codon.pair$aagaga/cps.aag.aga)
cps.obs.exp.aagagc<-log(n3.codon.pair$aagagc/cps.aag.agc)
cps.obs.exp.aagagg<-log(n3.codon.pair$aagagg/cps.aag.agg)
cps.obs.exp.aagagt<-log(n3.codon.pair$aagagt/cps.aag.agt)
cps.obs.exp.aagata<-log(n3.codon.pair$aagata/cps.aag.ata)
cps.obs.exp.aagatc<-log(n3.codon.pair$aagatc/cps.aag.atc)
cps.obs.exp.aagatg<-log(n3.codon.pair$aagatg/cps.aag.atg)
cps.obs.exp.aagatt<-log(n3.codon.pair$aagatt/cps.aag.att)

cps.obs.exp.aagcaa<-log(n3.codon.pair$aagcaa/cps.aag.caa)
cps.obs.exp.aagcac<-log(n3.codon.pair$aagcac/cps.aag.cac)
cps.obs.exp.aagcag<-log(n3.codon.pair$aagcag/cps.aag.cag)
cps.obs.exp.aagcat<-log(n3.codon.pair$aagcat/cps.aag.cat)
cps.obs.exp.aagcca<-log(n3.codon.pair$aagcca/cps.aag.cca)
cps.obs.exp.aagccc<-log(n3.codon.pair$aagccc/cps.aag.ccc)
cps.obs.exp.aagccg<-log(n3.codon.pair$aagccg/cps.aag.ccg)
cps.obs.exp.aagcct<-log(n3.codon.pair$aagcct/cps.aag.cct)
cps.obs.exp.aagcga<-log(n3.codon.pair$aagcga/cps.aag.cga)
cps.obs.exp.aagcgc<-log(n3.codon.pair$aagcgc/cps.aag.cgc)
cps.obs.exp.aagcgg<-log(n3.codon.pair$aagcgg/cps.aag.cgg)
cps.obs.exp.aagcgt<-log(n3.codon.pair$aagcgt/cps.aag.cgt)
cps.obs.exp.aagcta<-log(n3.codon.pair$aagcta/cps.aag.cta)
cps.obs.exp.aagctc<-log(n3.codon.pair$aagctc/cps.aag.ctc)
cps.obs.exp.aagctg<-log(n3.codon.pair$aagctg/cps.aag.ctg)
cps.obs.exp.aagctt<-log(n3.codon.pair$aagctt/cps.aag.ctt)

cps.obs.exp.aaggaa<-log(n3.codon.pair$aaggaa/cps.aag.gaa)
cps.obs.exp.aaggac<-log(n3.codon.pair$aaggac/cps.aag.gac)
cps.obs.exp.aaggag<-log(n3.codon.pair$aaggag/cps.aag.gag)
cps.obs.exp.aaggat<-log(n3.codon.pair$aaggat/cps.aag.gat)
cps.obs.exp.aaggca<-log(n3.codon.pair$aaggca/cps.aag.gca)
cps.obs.exp.aaggcc<-log(n3.codon.pair$aaggcc/cps.aag.gcc)
cps.obs.exp.aaggcg<-log(n3.codon.pair$aaggcg/cps.aag.gcg)
cps.obs.exp.aaggct<-log(n3.codon.pair$aaggct/cps.aag.gct)
cps.obs.exp.aaggga<-log(n3.codon.pair$aaggga/cps.aag.gga)
cps.obs.exp.aagggc<-log(n3.codon.pair$aagggc/cps.aag.ggc)
cps.obs.exp.aagggg<-log(n3.codon.pair$aagggg/cps.aag.ggg)
cps.obs.exp.aagggt<-log(n3.codon.pair$aagggt/cps.aag.ggt)
cps.obs.exp.aaggta<-log(n3.codon.pair$aaggta/cps.aag.gta)
cps.obs.exp.aaggtc<-log(n3.codon.pair$aaggtc/cps.aag.gtc)
cps.obs.exp.aaggtg<-log(n3.codon.pair$aaggtg/cps.aag.gtg)
cps.obs.exp.aaggtt<-log(n3.codon.pair$aaggtt/cps.aag.gtt)

#cps.obs.exp.aagtaa<-log(n3.codon.pair$aagtaa/cps.aag.taa)
cps.obs.exp.aagtac<-log(n3.codon.pair$aagtac/cps.aag.tac)
#cps.obs.exp.aagtag<-log(n3.codon.pair$aagtag/cps.aag.tag)
cps.obs.exp.aagtat<-log(n3.codon.pair$aagtat/cps.aag.tat)
cps.obs.exp.aagtca<-log(n3.codon.pair$aagtca/cps.aag.tca)
cps.obs.exp.aagtcc<-log(n3.codon.pair$aagtcc/cps.aag.tcc)
cps.obs.exp.aagtcg<-log(n3.codon.pair$aagtcg/cps.aag.tcg)
cps.obs.exp.aagtct<-log(n3.codon.pair$aagtct/cps.aag.tct)
#cps.obs.exp.aagtga<-log(n3.codon.pair$aagtga/cps.aag.tga)
cps.obs.exp.aagtgc<-log(n3.codon.pair$aagtgc/cps.aag.tgc)
cps.obs.exp.aagtgg<-log(n3.codon.pair$aagtgg/cps.aag.tgg)
cps.obs.exp.aagtgt<-log(n3.codon.pair$aagtgt/cps.aag.tgt)
cps.obs.exp.aagtta<-log(n3.codon.pair$aagtta/cps.aag.tta)
cps.obs.exp.aagttc<-log(n3.codon.pair$aagttc/cps.aag.ttc)
cps.obs.exp.aagttg<-log(n3.codon.pair$aagttg/cps.aag.ttg)
cps.obs.exp.aagttt<-log(n3.codon.pair$aagttt/cps.aag.ttt)







cps.obs.exp.aataaa<-log(n3.codon.pair$aataaa/cps.aat.aaa)
cps.obs.exp.aataac<-log(n3.codon.pair$aataac/cps.aat.aac)
cps.obs.exp.aataag<-log(n3.codon.pair$aataag/cps.aat.aag)
cps.obs.exp.aataat<-log(n3.codon.pair$aataat/cps.aat.aat)
cps.obs.exp.aataca<-log(n3.codon.pair$aataca/cps.aat.aca)
cps.obs.exp.aatacc<-log(n3.codon.pair$aatacc/cps.aat.acc)
cps.obs.exp.aatacg<-log(n3.codon.pair$aatacg/cps.aat.acg)
cps.obs.exp.aatact<-log(n3.codon.pair$aatact/cps.aat.act)
cps.obs.exp.aataga<-log(n3.codon.pair$aataga/cps.aat.aga)
cps.obs.exp.aatagc<-log(n3.codon.pair$aatagc/cps.aat.agc)
cps.obs.exp.aatagg<-log(n3.codon.pair$aatagg/cps.aat.agg)
cps.obs.exp.aatagt<-log(n3.codon.pair$aatagt/cps.aat.agt)
cps.obs.exp.aatata<-log(n3.codon.pair$aatata/cps.aat.ata)
cps.obs.exp.aatatc<-log(n3.codon.pair$aatatc/cps.aat.atc)
cps.obs.exp.aatatg<-log(n3.codon.pair$aatatg/cps.aat.atg)
cps.obs.exp.aatatt<-log(n3.codon.pair$aatatt/cps.aat.att)

cps.obs.exp.aatcaa<-log(n3.codon.pair$aatcaa/cps.aat.caa)
cps.obs.exp.aatcac<-log(n3.codon.pair$aatcac/cps.aat.cac)
cps.obs.exp.aatcag<-log(n3.codon.pair$aatcag/cps.aat.cag)
cps.obs.exp.aatcat<-log(n3.codon.pair$aatcat/cps.aat.cat)
cps.obs.exp.aatcca<-log(n3.codon.pair$aatcca/cps.aat.cca)
cps.obs.exp.aatccc<-log(n3.codon.pair$aatccc/cps.aat.ccc)
cps.obs.exp.aatccg<-log(n3.codon.pair$aatccg/cps.aat.ccg)
cps.obs.exp.aatcct<-log(n3.codon.pair$aatcct/cps.aat.cct)
cps.obs.exp.aatcga<-log(n3.codon.pair$aatcga/cps.aat.cga)
cps.obs.exp.aatcgc<-log(n3.codon.pair$aatcgc/cps.aat.cgc)
cps.obs.exp.aatcgg<-log(n3.codon.pair$aatcgg/cps.aat.cgg)
cps.obs.exp.aatcgt<-log(n3.codon.pair$aatcgt/cps.aat.cgt)
cps.obs.exp.aatcta<-log(n3.codon.pair$aatcta/cps.aat.cta)
cps.obs.exp.aatctc<-log(n3.codon.pair$aatctc/cps.aat.ctc)
cps.obs.exp.aatctg<-log(n3.codon.pair$aatctg/cps.aat.ctg)
cps.obs.exp.aatctt<-log(n3.codon.pair$aatctt/cps.aat.ctt)

cps.obs.exp.aatgaa<-log(n3.codon.pair$aatgaa/cps.aat.gaa)
cps.obs.exp.aatgac<-log(n3.codon.pair$aatgac/cps.aat.gac)
cps.obs.exp.aatgag<-log(n3.codon.pair$aatgag/cps.aat.gag)
cps.obs.exp.aatgat<-log(n3.codon.pair$aatgat/cps.aat.gat)
cps.obs.exp.aatgca<-log(n3.codon.pair$aatgca/cps.aat.gca)
cps.obs.exp.aatgcc<-log(n3.codon.pair$aatgcc/cps.aat.gcc)
cps.obs.exp.aatgcg<-log(n3.codon.pair$aatgcg/cps.aat.gcg)
cps.obs.exp.aatgct<-log(n3.codon.pair$aatgct/cps.aat.gct)
cps.obs.exp.aatgga<-log(n3.codon.pair$aatgga/cps.aat.gga)
cps.obs.exp.aatggc<-log(n3.codon.pair$aatggc/cps.aat.ggc)
cps.obs.exp.aatggg<-log(n3.codon.pair$aatggg/cps.aat.ggg)
cps.obs.exp.aatggt<-log(n3.codon.pair$aatggt/cps.aat.ggt)
cps.obs.exp.aatgta<-log(n3.codon.pair$aatgta/cps.aat.gta)
cps.obs.exp.aatgtc<-log(n3.codon.pair$aatgtc/cps.aat.gtc)
cps.obs.exp.aatgtg<-log(n3.codon.pair$aatgtg/cps.aat.gtg)
cps.obs.exp.aatgtt<-log(n3.codon.pair$aatgtt/cps.aat.gtt)

#cps.obs.exp.aattaa<-log(n3.codon.pair$aattaa/cps.aat.taa)
cps.obs.exp.aattac<-log(n3.codon.pair$aattac/cps.aat.tac)
#cps.obs.exp.aattag<-log(n3.codon.pair$aattag/cps.aat.tag)
cps.obs.exp.aattat<-log(n3.codon.pair$aattat/cps.aat.tat)
cps.obs.exp.aattca<-log(n3.codon.pair$aattca/cps.aat.tca)
cps.obs.exp.aattcc<-log(n3.codon.pair$aattcc/cps.aat.tcc)
cps.obs.exp.aattcg<-log(n3.codon.pair$aattcg/cps.aat.tcg)
cps.obs.exp.aattct<-log(n3.codon.pair$aattct/cps.aat.tct)
#cps.obs.exp.aattga<-log(n3.codon.pair$aattga/cps.aat.tga)
cps.obs.exp.aattgc<-log(n3.codon.pair$aattgc/cps.aat.tgc)
cps.obs.exp.aattgg<-log(n3.codon.pair$aattgg/cps.aat.tgg)
cps.obs.exp.aattgt<-log(n3.codon.pair$aattgt/cps.aat.tgt)
cps.obs.exp.aattta<-log(n3.codon.pair$aattta/cps.aat.tta)
cps.obs.exp.aatttc<-log(n3.codon.pair$aatttc/cps.aat.ttc)
cps.obs.exp.aatttg<-log(n3.codon.pair$aatttg/cps.aat.ttg)
cps.obs.exp.aatttt<-log(n3.codon.pair$aatttt/cps.aat.ttt)

















cps.obs.exp.acaaaa<-log(n3.codon.pair$acaaaa/cps.aca.aaa)
cps.obs.exp.acaaac<-log(n3.codon.pair$acaaac/cps.aca.aac)
cps.obs.exp.acaaag<-log(n3.codon.pair$acaaag/cps.aca.aag)
cps.obs.exp.acaaat<-log(n3.codon.pair$acaaat/cps.aca.aat)
cps.obs.exp.acaaca<-log(n3.codon.pair$acaaca/cps.aca.aca)
cps.obs.exp.acaacc<-log(n3.codon.pair$acaacc/cps.aca.acc)
cps.obs.exp.acaacg<-log(n3.codon.pair$acaacg/cps.aca.acg)
cps.obs.exp.acaact<-log(n3.codon.pair$acaact/cps.aca.act)
cps.obs.exp.acaaga<-log(n3.codon.pair$acaaga/cps.aca.aga)
cps.obs.exp.acaagc<-log(n3.codon.pair$acaagc/cps.aca.agc)
cps.obs.exp.acaagg<-log(n3.codon.pair$acaagg/cps.aca.agg)
cps.obs.exp.acaagt<-log(n3.codon.pair$acaagt/cps.aca.agt)
cps.obs.exp.acaata<-log(n3.codon.pair$acaata/cps.aca.ata)
cps.obs.exp.acaatc<-log(n3.codon.pair$acaatc/cps.aca.atc)
cps.obs.exp.acaatg<-log(n3.codon.pair$acaatg/cps.aca.atg)
cps.obs.exp.acaatt<-log(n3.codon.pair$acaatt/cps.aca.att)

cps.obs.exp.acacaa<-log(n3.codon.pair$acacaa/cps.aca.caa)
cps.obs.exp.acacac<-log(n3.codon.pair$acacac/cps.aca.cac)
cps.obs.exp.acacag<-log(n3.codon.pair$acacag/cps.aca.cag)
cps.obs.exp.acacat<-log(n3.codon.pair$acacat/cps.aca.cat)
cps.obs.exp.acacca<-log(n3.codon.pair$acacca/cps.aca.cca)
cps.obs.exp.acaccc<-log(n3.codon.pair$acaccc/cps.aca.ccc)
cps.obs.exp.acaccg<-log(n3.codon.pair$acaccg/cps.aca.ccg)
cps.obs.exp.acacct<-log(n3.codon.pair$acacct/cps.aca.cct)
cps.obs.exp.acacga<-log(n3.codon.pair$acacga/cps.aca.cga)
cps.obs.exp.acacgc<-log(n3.codon.pair$acacgc/cps.aca.cgc)
cps.obs.exp.acacgg<-log(n3.codon.pair$acacgg/cps.aca.cgg)
cps.obs.exp.acacgt<-log(n3.codon.pair$acacgt/cps.aca.cgt)
cps.obs.exp.acacta<-log(n3.codon.pair$acacta/cps.aca.cta)
cps.obs.exp.acactc<-log(n3.codon.pair$acactc/cps.aca.ctc)
cps.obs.exp.acactg<-log(n3.codon.pair$acactg/cps.aca.ctg)
cps.obs.exp.acactt<-log(n3.codon.pair$acactt/cps.aca.ctt)

cps.obs.exp.acagaa<-log(n3.codon.pair$acagaa/cps.aca.gaa)
cps.obs.exp.acagac<-log(n3.codon.pair$acagac/cps.aca.gac)
cps.obs.exp.acagag<-log(n3.codon.pair$acagag/cps.aca.gag)
cps.obs.exp.acagat<-log(n3.codon.pair$acagat/cps.aca.gat)
cps.obs.exp.acagca<-log(n3.codon.pair$acagca/cps.aca.gca)
cps.obs.exp.acagcc<-log(n3.codon.pair$acagcc/cps.aca.gcc)
cps.obs.exp.acagcg<-log(n3.codon.pair$acagcg/cps.aca.gcg)
cps.obs.exp.acagct<-log(n3.codon.pair$acagct/cps.aca.gct)
cps.obs.exp.acagga<-log(n3.codon.pair$acagga/cps.aca.gga)
cps.obs.exp.acaggc<-log(n3.codon.pair$acaggc/cps.aca.ggc)
cps.obs.exp.acaggg<-log(n3.codon.pair$acaggg/cps.aca.ggg)
cps.obs.exp.acaggt<-log(n3.codon.pair$acaggt/cps.aca.ggt)
cps.obs.exp.acagta<-log(n3.codon.pair$acagta/cps.aca.gta)
cps.obs.exp.acagtc<-log(n3.codon.pair$acagtc/cps.aca.gtc)
cps.obs.exp.acagtg<-log(n3.codon.pair$acagtg/cps.aca.gtg)
cps.obs.exp.acagtt<-log(n3.codon.pair$acagtt/cps.aca.gtt)

#cps.obs.exp.acataa<-log(n3.codon.pair$acataa/cps.aca.taa)
cps.obs.exp.acatac<-log(n3.codon.pair$acatac/cps.aca.tac)
#cps.obs.exp.acatag<-log(n3.codon.pair$acatag/cps.aca.tag)
cps.obs.exp.acatat<-log(n3.codon.pair$acatat/cps.aca.tat)
cps.obs.exp.acatca<-log(n3.codon.pair$acatca/cps.aca.tca)
cps.obs.exp.acatcc<-log(n3.codon.pair$acatcc/cps.aca.tcc)
cps.obs.exp.acatcg<-log(n3.codon.pair$acatcg/cps.aca.tcg)
cps.obs.exp.acatct<-log(n3.codon.pair$acatct/cps.aca.tct)
#cps.obs.exp.acatga<-log(n3.codon.pair$acatga/cps.aca.tga)
cps.obs.exp.acatgc<-log(n3.codon.pair$acatgc/cps.aca.tgc)
cps.obs.exp.acatgg<-log(n3.codon.pair$acatgg/cps.aca.tgg)
cps.obs.exp.acatgt<-log(n3.codon.pair$acatgt/cps.aca.tgt)
cps.obs.exp.acatta<-log(n3.codon.pair$acatta/cps.aca.tta)
cps.obs.exp.acattc<-log(n3.codon.pair$acattc/cps.aca.ttc)
cps.obs.exp.acattg<-log(n3.codon.pair$acattg/cps.aca.ttg)
cps.obs.exp.acattt<-log(n3.codon.pair$acattt/cps.aca.ttt)









cps.obs.exp.accaaa<-log(n3.codon.pair$accaaa/cps.acc.aaa)
cps.obs.exp.accaac<-log(n3.codon.pair$accaac/cps.acc.aac)
cps.obs.exp.accaag<-log(n3.codon.pair$accaag/cps.acc.aag)
cps.obs.exp.accaat<-log(n3.codon.pair$accaat/cps.acc.aat)
cps.obs.exp.accaca<-log(n3.codon.pair$accaca/cps.acc.aca)
cps.obs.exp.accacc<-log(n3.codon.pair$accacc/cps.acc.acc)
cps.obs.exp.accacg<-log(n3.codon.pair$accacg/cps.acc.acg)
cps.obs.exp.accact<-log(n3.codon.pair$accact/cps.acc.act)
cps.obs.exp.accaga<-log(n3.codon.pair$accaga/cps.acc.aga)
cps.obs.exp.accagc<-log(n3.codon.pair$accagc/cps.acc.agc)
cps.obs.exp.accagg<-log(n3.codon.pair$accagg/cps.acc.agg)
cps.obs.exp.accagt<-log(n3.codon.pair$accagt/cps.acc.agt)
cps.obs.exp.accata<-log(n3.codon.pair$accata/cps.acc.ata)
cps.obs.exp.accatc<-log(n3.codon.pair$accatc/cps.acc.atc)
cps.obs.exp.accatg<-log(n3.codon.pair$accatg/cps.acc.atg)
cps.obs.exp.accatt<-log(n3.codon.pair$accatt/cps.acc.att)

cps.obs.exp.acccaa<-log(n3.codon.pair$acccaa/cps.acc.caa)
cps.obs.exp.acccac<-log(n3.codon.pair$acccac/cps.acc.cac)
cps.obs.exp.acccag<-log(n3.codon.pair$acccag/cps.acc.cag)
cps.obs.exp.acccat<-log(n3.codon.pair$acccat/cps.acc.cat)
cps.obs.exp.acccca<-log(n3.codon.pair$acccca/cps.acc.cca)
cps.obs.exp.accccc<-log(n3.codon.pair$accccc/cps.acc.ccc)
cps.obs.exp.accccg<-log(n3.codon.pair$accccg/cps.acc.ccg)
cps.obs.exp.acccct<-log(n3.codon.pair$acccct/cps.acc.cct)
cps.obs.exp.acccga<-log(n3.codon.pair$acccga/cps.acc.cga)
cps.obs.exp.acccgc<-log(n3.codon.pair$acccgc/cps.acc.cgc)
cps.obs.exp.acccgg<-log(n3.codon.pair$acccgg/cps.acc.cgg)
cps.obs.exp.acccgt<-log(n3.codon.pair$acccgt/cps.acc.cgt)
cps.obs.exp.acccta<-log(n3.codon.pair$acccta/cps.acc.cta)
cps.obs.exp.accctc<-log(n3.codon.pair$accctc/cps.acc.ctc)
cps.obs.exp.accctg<-log(n3.codon.pair$accctg/cps.acc.ctg)
cps.obs.exp.accctt<-log(n3.codon.pair$accctt/cps.acc.ctt)

cps.obs.exp.accgaa<-log(n3.codon.pair$accgaa/cps.acc.gaa)
cps.obs.exp.accgac<-log(n3.codon.pair$accgac/cps.acc.gac)
cps.obs.exp.accgag<-log(n3.codon.pair$accgag/cps.acc.gag)
cps.obs.exp.accgat<-log(n3.codon.pair$accgat/cps.acc.gat)
cps.obs.exp.accgca<-log(n3.codon.pair$accgca/cps.acc.gca)
cps.obs.exp.accgcc<-log(n3.codon.pair$accgcc/cps.acc.gcc)
cps.obs.exp.accgcg<-log(n3.codon.pair$accgcg/cps.acc.gcg)
cps.obs.exp.accgct<-log(n3.codon.pair$accgct/cps.acc.gct)
cps.obs.exp.accgga<-log(n3.codon.pair$accgga/cps.acc.gga)
cps.obs.exp.accggc<-log(n3.codon.pair$accggc/cps.acc.ggc)
cps.obs.exp.accggg<-log(n3.codon.pair$accggg/cps.acc.ggg)
cps.obs.exp.accggt<-log(n3.codon.pair$accggt/cps.acc.ggt)
cps.obs.exp.accgta<-log(n3.codon.pair$accgta/cps.acc.gta)
cps.obs.exp.accgtc<-log(n3.codon.pair$accgtc/cps.acc.gtc)
cps.obs.exp.accgtg<-log(n3.codon.pair$accgtg/cps.acc.gtg)
cps.obs.exp.accgtt<-log(n3.codon.pair$accgtt/cps.acc.gtt)

#cps.obs.exp.acctaa<-log(n3.codon.pair$acctaa/cps.acc.taa)
cps.obs.exp.acctac<-log(n3.codon.pair$acctac/cps.acc.tac)
#cps.obs.exp.acctag<-log(n3.codon.pair$acctag/cps.acc.tag)
cps.obs.exp.acctat<-log(n3.codon.pair$acctat/cps.acc.tat)
cps.obs.exp.acctca<-log(n3.codon.pair$acctca/cps.acc.tca)
cps.obs.exp.acctcc<-log(n3.codon.pair$acctcc/cps.acc.tcc)
cps.obs.exp.acctcg<-log(n3.codon.pair$acctcg/cps.acc.tcg)
cps.obs.exp.acctct<-log(n3.codon.pair$acctct/cps.acc.tct)
#cps.obs.exp.acctga<-log(n3.codon.pair$acctga/cps.acc.tga)
cps.obs.exp.acctgc<-log(n3.codon.pair$acctgc/cps.acc.tgc)
cps.obs.exp.acctgg<-log(n3.codon.pair$acctgg/cps.acc.tgg)
cps.obs.exp.acctgt<-log(n3.codon.pair$acctgt/cps.acc.tgt)
cps.obs.exp.acctta<-log(n3.codon.pair$acctta/cps.acc.tta)
cps.obs.exp.accttc<-log(n3.codon.pair$accttc/cps.acc.ttc)
cps.obs.exp.accttg<-log(n3.codon.pair$accttg/cps.acc.ttg)
cps.obs.exp.accttt<-log(n3.codon.pair$accttt/cps.acc.ttt)








cps.obs.exp.acgaaa<-log(n3.codon.pair$acgaaa/cps.acg.aaa)
cps.obs.exp.acgaac<-log(n3.codon.pair$acgaac/cps.acg.aac)
cps.obs.exp.acgaag<-log(n3.codon.pair$acgaag/cps.acg.aag)
cps.obs.exp.acgaat<-log(n3.codon.pair$acgaat/cps.acg.aat)
cps.obs.exp.acgaca<-log(n3.codon.pair$acgaca/cps.acg.aca)
cps.obs.exp.acgacc<-log(n3.codon.pair$acgacc/cps.acg.acc)
cps.obs.exp.acgacg<-log(n3.codon.pair$acgacg/cps.acg.acg)
cps.obs.exp.acgact<-log(n3.codon.pair$acgact/cps.acg.act)
cps.obs.exp.acgaga<-log(n3.codon.pair$acgaga/cps.acg.aga)
cps.obs.exp.acgagc<-log(n3.codon.pair$acgagc/cps.acg.agc)
cps.obs.exp.acgagg<-log(n3.codon.pair$acgagg/cps.acg.agg)
cps.obs.exp.acgagt<-log(n3.codon.pair$acgagt/cps.acg.agt)
cps.obs.exp.acgata<-log(n3.codon.pair$acgata/cps.acg.ata)
cps.obs.exp.acgatc<-log(n3.codon.pair$acgatc/cps.acg.atc)
cps.obs.exp.acgatg<-log(n3.codon.pair$acgatg/cps.acg.atg)
cps.obs.exp.acgatt<-log(n3.codon.pair$acgatt/cps.acg.att)

cps.obs.exp.acgcaa<-log(n3.codon.pair$acgcaa/cps.acg.caa)
cps.obs.exp.acgcac<-log(n3.codon.pair$acgcac/cps.acg.cac)
cps.obs.exp.acgcag<-log(n3.codon.pair$acgcag/cps.acg.cag)
cps.obs.exp.acgcat<-log(n3.codon.pair$acgcat/cps.acg.cat)
cps.obs.exp.acgcca<-log(n3.codon.pair$acgcca/cps.acg.cca)
cps.obs.exp.acgccc<-log(n3.codon.pair$acgccc/cps.acg.ccc)
cps.obs.exp.acgccg<-log(n3.codon.pair$acgccg/cps.acg.ccg)
cps.obs.exp.acgcct<-log(n3.codon.pair$acgcct/cps.acg.cct)
cps.obs.exp.acgcga<-log(n3.codon.pair$acgcga/cps.acg.cga)
cps.obs.exp.acgcgc<-log(n3.codon.pair$acgcgc/cps.acg.cgc)
cps.obs.exp.acgcgg<-log(n3.codon.pair$acgcgg/cps.acg.cgg)
cps.obs.exp.acgcgt<-log(n3.codon.pair$acgcgt/cps.acg.cgt)
cps.obs.exp.acgcta<-log(n3.codon.pair$acgcta/cps.acg.cta)
cps.obs.exp.acgctc<-log(n3.codon.pair$acgctc/cps.acg.ctc)
cps.obs.exp.acgctg<-log(n3.codon.pair$acgctg/cps.acg.ctg)
cps.obs.exp.acgctt<-log(n3.codon.pair$acgctt/cps.acg.ctt)

cps.obs.exp.acggaa<-log(n3.codon.pair$acggaa/cps.acg.gaa)
cps.obs.exp.acggac<-log(n3.codon.pair$acggac/cps.acg.gac)
cps.obs.exp.acggag<-log(n3.codon.pair$acggag/cps.acg.gag)
cps.obs.exp.acggat<-log(n3.codon.pair$acggat/cps.acg.gat)
cps.obs.exp.acggca<-log(n3.codon.pair$acggca/cps.acg.gca)
cps.obs.exp.acggcc<-log(n3.codon.pair$acggcc/cps.acg.gcc)
cps.obs.exp.acggcg<-log(n3.codon.pair$acggcg/cps.acg.gcg)
cps.obs.exp.acggct<-log(n3.codon.pair$acggct/cps.acg.gct)
cps.obs.exp.acggga<-log(n3.codon.pair$acggga/cps.acg.gga)
cps.obs.exp.acgggc<-log(n3.codon.pair$acgggc/cps.acg.ggc)
cps.obs.exp.acgggg<-log(n3.codon.pair$acgggg/cps.acg.ggg)
cps.obs.exp.acgggt<-log(n3.codon.pair$acgggt/cps.acg.ggt)
cps.obs.exp.acggta<-log(n3.codon.pair$acggta/cps.acg.gta)
cps.obs.exp.acggtc<-log(n3.codon.pair$acggtc/cps.acg.gtc)
cps.obs.exp.acggtg<-log(n3.codon.pair$acggtg/cps.acg.gtg)
cps.obs.exp.acggtt<-log(n3.codon.pair$acggtt/cps.acg.gtt)

#cps.obs.exp.acgtaa<-log(n3.codon.pair$acgtaa/cps.acg.taa)
cps.obs.exp.acgtac<-log(n3.codon.pair$acgtac/cps.acg.tac)
#cps.obs.exp.acgtag<-log(n3.codon.pair$acgtag/cps.acg.tag)
cps.obs.exp.acgtat<-log(n3.codon.pair$acgtat/cps.acg.tat)
cps.obs.exp.acgtca<-log(n3.codon.pair$acgtca/cps.acg.tca)
cps.obs.exp.acgtcc<-log(n3.codon.pair$acgtcc/cps.acg.tcc)
cps.obs.exp.acgtcg<-log(n3.codon.pair$acgtcg/cps.acg.tcg)
cps.obs.exp.acgtct<-log(n3.codon.pair$acgtct/cps.acg.tct)
#cps.obs.exp.acgtga<-log(n3.codon.pair$acgtga/cps.acg.tga)
cps.obs.exp.acgtgc<-log(n3.codon.pair$acgtgc/cps.acg.tgc)
cps.obs.exp.acgtgg<-log(n3.codon.pair$acgtgg/cps.acg.tgg)
cps.obs.exp.acgtgt<-log(n3.codon.pair$acgtgt/cps.acg.tgt)
cps.obs.exp.acgtta<-log(n3.codon.pair$acgtta/cps.acg.tta)
cps.obs.exp.acgttc<-log(n3.codon.pair$acgttc/cps.acg.ttc)
cps.obs.exp.acgttg<-log(n3.codon.pair$acgttg/cps.acg.ttg)
cps.obs.exp.acgttt<-log(n3.codon.pair$acgttt/cps.acg.ttt)








cps.obs.exp.actaaa<-log(n3.codon.pair$actaaa/cps.act.aaa)
cps.obs.exp.actaac<-log(n3.codon.pair$actaac/cps.act.aac)
cps.obs.exp.actaag<-log(n3.codon.pair$actaag/cps.act.aag)
cps.obs.exp.actaat<-log(n3.codon.pair$actaat/cps.act.aat)
cps.obs.exp.actaca<-log(n3.codon.pair$actaca/cps.act.aca)
cps.obs.exp.actacc<-log(n3.codon.pair$actacc/cps.act.acc)
cps.obs.exp.actacg<-log(n3.codon.pair$actacg/cps.act.acg)
cps.obs.exp.actact<-log(n3.codon.pair$actact/cps.act.act)
cps.obs.exp.actaga<-log(n3.codon.pair$actaga/cps.act.aga)
cps.obs.exp.actagc<-log(n3.codon.pair$actagc/cps.act.agc)
cps.obs.exp.actagg<-log(n3.codon.pair$actagg/cps.act.agg)
cps.obs.exp.actagt<-log(n3.codon.pair$actagt/cps.act.agt)
cps.obs.exp.actata<-log(n3.codon.pair$actata/cps.act.ata)
cps.obs.exp.actatc<-log(n3.codon.pair$actatc/cps.act.atc)
cps.obs.exp.actatg<-log(n3.codon.pair$actatg/cps.act.atg)
cps.obs.exp.actatt<-log(n3.codon.pair$actatt/cps.act.att)

cps.obs.exp.actcaa<-log(n3.codon.pair$actcaa/cps.act.caa)
cps.obs.exp.actcac<-log(n3.codon.pair$actcac/cps.act.cac)
cps.obs.exp.actcag<-log(n3.codon.pair$actcag/cps.act.cag)
cps.obs.exp.actcat<-log(n3.codon.pair$actcat/cps.act.cat)
cps.obs.exp.actcca<-log(n3.codon.pair$actcca/cps.act.cca)
cps.obs.exp.actccc<-log(n3.codon.pair$actccc/cps.act.ccc)
cps.obs.exp.actccg<-log(n3.codon.pair$actccg/cps.act.ccg)
cps.obs.exp.actcct<-log(n3.codon.pair$actcct/cps.act.cct)
cps.obs.exp.actcga<-log(n3.codon.pair$actcga/cps.act.cga)
cps.obs.exp.actcgc<-log(n3.codon.pair$actcgc/cps.act.cgc)
cps.obs.exp.actcgg<-log(n3.codon.pair$actcgg/cps.act.cgg)
cps.obs.exp.actcgt<-log(n3.codon.pair$actcgt/cps.act.cgt)
cps.obs.exp.actcta<-log(n3.codon.pair$actcta/cps.act.cta)
cps.obs.exp.actctc<-log(n3.codon.pair$actctc/cps.act.ctc)
cps.obs.exp.actctg<-log(n3.codon.pair$actctg/cps.act.ctg)
cps.obs.exp.actctt<-log(n3.codon.pair$actctt/cps.act.ctt)

cps.obs.exp.actgaa<-log(n3.codon.pair$actgaa/cps.act.gaa)
cps.obs.exp.actgac<-log(n3.codon.pair$actgac/cps.act.gac)
cps.obs.exp.actgag<-log(n3.codon.pair$actgag/cps.act.gag)
cps.obs.exp.actgat<-log(n3.codon.pair$actgat/cps.act.gat)
cps.obs.exp.actgca<-log(n3.codon.pair$actgca/cps.act.gca)
cps.obs.exp.actgcc<-log(n3.codon.pair$actgcc/cps.act.gcc)
cps.obs.exp.actgcg<-log(n3.codon.pair$actgcg/cps.act.gcg)
cps.obs.exp.actgct<-log(n3.codon.pair$actgct/cps.act.gct)
cps.obs.exp.actgga<-log(n3.codon.pair$actgga/cps.act.gga)
cps.obs.exp.actggc<-log(n3.codon.pair$actggc/cps.act.ggc)
cps.obs.exp.actggg<-log(n3.codon.pair$actggg/cps.act.ggg)
cps.obs.exp.actggt<-log(n3.codon.pair$actggt/cps.act.ggt)
cps.obs.exp.actgta<-log(n3.codon.pair$actgta/cps.act.gta)
cps.obs.exp.actgtc<-log(n3.codon.pair$actgtc/cps.act.gtc)
cps.obs.exp.actgtg<-log(n3.codon.pair$actgtg/cps.act.gtg)
cps.obs.exp.actgtt<-log(n3.codon.pair$actgtt/cps.act.gtt)

#cps.obs.exp.acttaa<-log(n3.codon.pair$acttaa/cps.act.taa)
cps.obs.exp.acttac<-log(n3.codon.pair$acttac/cps.act.tac)
#cps.obs.exp.acttag<-log(n3.codon.pair$acttag/cps.act.tag)
cps.obs.exp.acttat<-log(n3.codon.pair$acttat/cps.act.tat)
cps.obs.exp.acttca<-log(n3.codon.pair$acttca/cps.act.tca)
cps.obs.exp.acttcc<-log(n3.codon.pair$acttcc/cps.act.tcc)
cps.obs.exp.acttcg<-log(n3.codon.pair$acttcg/cps.act.tcg)
cps.obs.exp.acttct<-log(n3.codon.pair$acttct/cps.act.tct)
#cps.obs.exp.acttga<-log(n3.codon.pair$acttga/cps.act.tga)
cps.obs.exp.acttgc<-log(n3.codon.pair$acttgc/cps.act.tgc)
cps.obs.exp.acttgg<-log(n3.codon.pair$acttgg/cps.act.tgg)
cps.obs.exp.acttgt<-log(n3.codon.pair$acttgt/cps.act.tgt)
cps.obs.exp.acttta<-log(n3.codon.pair$acttta/cps.act.tta)
cps.obs.exp.actttc<-log(n3.codon.pair$actttc/cps.act.ttc)
cps.obs.exp.actttg<-log(n3.codon.pair$actttg/cps.act.ttg)
cps.obs.exp.actttt<-log(n3.codon.pair$actttt/cps.act.ttt)



















cps.obs.exp.agaaaa<-log(n3.codon.pair$agaaaa/cps.aga.aaa)
cps.obs.exp.agaaac<-log(n3.codon.pair$agaaac/cps.aga.aac)
cps.obs.exp.agaaag<-log(n3.codon.pair$agaaag/cps.aga.aag)
cps.obs.exp.agaaat<-log(n3.codon.pair$agaaat/cps.aga.aat)
cps.obs.exp.agaaca<-log(n3.codon.pair$agaaca/cps.aga.aca)
cps.obs.exp.agaacc<-log(n3.codon.pair$agaacc/cps.aga.acc)
cps.obs.exp.agaacg<-log(n3.codon.pair$agaacg/cps.aga.acg)
cps.obs.exp.agaact<-log(n3.codon.pair$agaact/cps.aga.act)
cps.obs.exp.agaaga<-log(n3.codon.pair$agaaga/cps.aga.aga)
cps.obs.exp.agaagc<-log(n3.codon.pair$agaagc/cps.aga.agc)
cps.obs.exp.agaagg<-log(n3.codon.pair$agaagg/cps.aga.agg)
cps.obs.exp.agaagt<-log(n3.codon.pair$agaagt/cps.aga.agt)
cps.obs.exp.agaata<-log(n3.codon.pair$agaata/cps.aga.ata)
cps.obs.exp.agaatc<-log(n3.codon.pair$agaatc/cps.aga.atc)
cps.obs.exp.agaatg<-log(n3.codon.pair$agaatg/cps.aga.atg)
cps.obs.exp.agaatt<-log(n3.codon.pair$agaatt/cps.aga.att)

cps.obs.exp.agacaa<-log(n3.codon.pair$agacaa/cps.aga.caa)
cps.obs.exp.agacac<-log(n3.codon.pair$agacac/cps.aga.cac)
cps.obs.exp.agacag<-log(n3.codon.pair$agacag/cps.aga.cag)
cps.obs.exp.agacat<-log(n3.codon.pair$agacat/cps.aga.cat)
cps.obs.exp.agacca<-log(n3.codon.pair$agacca/cps.aga.cca)
cps.obs.exp.agaccc<-log(n3.codon.pair$agaccc/cps.aga.ccc)
cps.obs.exp.agaccg<-log(n3.codon.pair$agaccg/cps.aga.ccg)
cps.obs.exp.agacct<-log(n3.codon.pair$agacct/cps.aga.cct)
cps.obs.exp.agacga<-log(n3.codon.pair$agacga/cps.aga.cga)
cps.obs.exp.agacgc<-log(n3.codon.pair$agacgc/cps.aga.cgc)
cps.obs.exp.agacgg<-log(n3.codon.pair$agacgg/cps.aga.cgg)
cps.obs.exp.agacgt<-log(n3.codon.pair$agacgt/cps.aga.cgt)
cps.obs.exp.agacta<-log(n3.codon.pair$agacta/cps.aga.cta)
cps.obs.exp.agactc<-log(n3.codon.pair$agactc/cps.aga.ctc)
cps.obs.exp.agactg<-log(n3.codon.pair$agactg/cps.aga.ctg)
cps.obs.exp.agactt<-log(n3.codon.pair$agactt/cps.aga.ctt)

cps.obs.exp.agagaa<-log(n3.codon.pair$agagaa/cps.aga.gaa)
cps.obs.exp.agagac<-log(n3.codon.pair$agagac/cps.aga.gac)
cps.obs.exp.agagag<-log(n3.codon.pair$agagag/cps.aga.gag)
cps.obs.exp.agagat<-log(n3.codon.pair$agagat/cps.aga.gat)
cps.obs.exp.agagca<-log(n3.codon.pair$agagca/cps.aga.gca)
cps.obs.exp.agagcc<-log(n3.codon.pair$agagcc/cps.aga.gcc)
cps.obs.exp.agagcg<-log(n3.codon.pair$agagcg/cps.aga.gcg)
cps.obs.exp.agagct<-log(n3.codon.pair$agagct/cps.aga.gct)
cps.obs.exp.agagga<-log(n3.codon.pair$agagga/cps.aga.gga)
cps.obs.exp.agaggc<-log(n3.codon.pair$agaggc/cps.aga.ggc)
cps.obs.exp.agaggg<-log(n3.codon.pair$agaggg/cps.aga.ggg)
cps.obs.exp.agaggt<-log(n3.codon.pair$agaggt/cps.aga.ggt)
cps.obs.exp.agagta<-log(n3.codon.pair$agagta/cps.aga.gta)
cps.obs.exp.agagtc<-log(n3.codon.pair$agagtc/cps.aga.gtc)
cps.obs.exp.agagtg<-log(n3.codon.pair$agagtg/cps.aga.gtg)
cps.obs.exp.agagtt<-log(n3.codon.pair$agagtt/cps.aga.gtt)

#cps.obs.exp.agataa<-log(n3.codon.pair$agataa/cps.aga.taa)
cps.obs.exp.agatac<-log(n3.codon.pair$agatac/cps.aga.tac)
#cps.obs.exp.agatag<-log(n3.codon.pair$agatag/cps.aga.tag)
cps.obs.exp.agatat<-log(n3.codon.pair$agatat/cps.aga.tat)
cps.obs.exp.agatca<-log(n3.codon.pair$agatca/cps.aga.tca)
cps.obs.exp.agatcc<-log(n3.codon.pair$agatcc/cps.aga.tcc)
cps.obs.exp.agatcg<-log(n3.codon.pair$agatcg/cps.aga.tcg)
cps.obs.exp.agatct<-log(n3.codon.pair$agatct/cps.aga.tct)
#cps.obs.exp.agatga<-log(n3.codon.pair$agatga/cps.aga.tga)
cps.obs.exp.agatgc<-log(n3.codon.pair$agatgc/cps.aga.tgc)
cps.obs.exp.agatgg<-log(n3.codon.pair$agatgg/cps.aga.tgg)
cps.obs.exp.agatgt<-log(n3.codon.pair$agatgt/cps.aga.tgt)
cps.obs.exp.agatta<-log(n3.codon.pair$agatta/cps.aga.tta)
cps.obs.exp.agattc<-log(n3.codon.pair$agattc/cps.aga.ttc)
cps.obs.exp.agattg<-log(n3.codon.pair$agattg/cps.aga.ttg)
cps.obs.exp.agattt<-log(n3.codon.pair$agattt/cps.aga.ttt)







cps.obs.exp.agcaaa<-log(n3.codon.pair$agcaaa/cps.agc.aaa)
cps.obs.exp.agcaac<-log(n3.codon.pair$agcaac/cps.agc.aac)
cps.obs.exp.agcaag<-log(n3.codon.pair$agcaag/cps.agc.aag)
cps.obs.exp.agcaat<-log(n3.codon.pair$agcaat/cps.agc.aat)
cps.obs.exp.agcaca<-log(n3.codon.pair$agcaca/cps.agc.aca)
cps.obs.exp.agcacc<-log(n3.codon.pair$agcacc/cps.agc.acc)
cps.obs.exp.agcacg<-log(n3.codon.pair$agcacg/cps.agc.acg)
cps.obs.exp.agcact<-log(n3.codon.pair$agcact/cps.agc.act)
cps.obs.exp.agcaga<-log(n3.codon.pair$agcaga/cps.agc.aga)
cps.obs.exp.agcagc<-log(n3.codon.pair$agcagc/cps.agc.agc)
cps.obs.exp.agcagg<-log(n3.codon.pair$agcagg/cps.agc.agg)
cps.obs.exp.agcagt<-log(n3.codon.pair$agcagt/cps.agc.agt)
cps.obs.exp.agcata<-log(n3.codon.pair$agcata/cps.agc.ata)
cps.obs.exp.agcatc<-log(n3.codon.pair$agcatc/cps.agc.atc)
cps.obs.exp.agcatg<-log(n3.codon.pair$agcatg/cps.agc.atg)
cps.obs.exp.agcatt<-log(n3.codon.pair$agcatt/cps.agc.att)

cps.obs.exp.agccaa<-log(n3.codon.pair$agccaa/cps.agc.caa)
cps.obs.exp.agccac<-log(n3.codon.pair$agccac/cps.agc.cac)
cps.obs.exp.agccag<-log(n3.codon.pair$agccag/cps.agc.cag)
cps.obs.exp.agccat<-log(n3.codon.pair$agccat/cps.agc.cat)
cps.obs.exp.agccca<-log(n3.codon.pair$agccca/cps.agc.cca)
cps.obs.exp.agcccc<-log(n3.codon.pair$agcccc/cps.agc.ccc)
cps.obs.exp.agcccg<-log(n3.codon.pair$agcccg/cps.agc.ccg)
cps.obs.exp.agccct<-log(n3.codon.pair$agccct/cps.agc.cct)
cps.obs.exp.agccga<-log(n3.codon.pair$agccga/cps.agc.cga)
cps.obs.exp.agccgc<-log(n3.codon.pair$agccgc/cps.agc.cgc)
cps.obs.exp.agccgg<-log(n3.codon.pair$agccgg/cps.agc.cgg)
cps.obs.exp.agccgt<-log(n3.codon.pair$agccgt/cps.agc.cgt)
cps.obs.exp.agccta<-log(n3.codon.pair$agccta/cps.agc.cta)
cps.obs.exp.agcctc<-log(n3.codon.pair$agcctc/cps.agc.ctc)
cps.obs.exp.agcctg<-log(n3.codon.pair$agcctg/cps.agc.ctg)
cps.obs.exp.agcctt<-log(n3.codon.pair$agcctt/cps.agc.ctt)

cps.obs.exp.agcgaa<-log(n3.codon.pair$agcgaa/cps.agc.gaa)
cps.obs.exp.agcgac<-log(n3.codon.pair$agcgac/cps.agc.gac)
cps.obs.exp.agcgag<-log(n3.codon.pair$agcgag/cps.agc.gag)
cps.obs.exp.agcgat<-log(n3.codon.pair$agcgat/cps.agc.gat)
cps.obs.exp.agcgca<-log(n3.codon.pair$agcgca/cps.agc.gca)
cps.obs.exp.agcgcc<-log(n3.codon.pair$agcgcc/cps.agc.gcc)
cps.obs.exp.agcgcg<-log(n3.codon.pair$agcgcg/cps.agc.gcg)
cps.obs.exp.agcgct<-log(n3.codon.pair$agcgct/cps.agc.gct)
cps.obs.exp.agcgga<-log(n3.codon.pair$agcgga/cps.agc.gga)
cps.obs.exp.agcggc<-log(n3.codon.pair$agcggc/cps.agc.ggc)
cps.obs.exp.agcggg<-log(n3.codon.pair$agcggg/cps.agc.ggg)
cps.obs.exp.agcggt<-log(n3.codon.pair$agcggt/cps.agc.ggt)
cps.obs.exp.agcgta<-log(n3.codon.pair$agcgta/cps.agc.gta)
cps.obs.exp.agcgtc<-log(n3.codon.pair$agcgtc/cps.agc.gtc)
cps.obs.exp.agcgtg<-log(n3.codon.pair$agcgtg/cps.agc.gtg)
cps.obs.exp.agcgtt<-log(n3.codon.pair$agcgtt/cps.agc.gtt)

#cps.obs.exp.agctaa<-log(n3.codon.pair$agctaa/cps.agc.taa)
cps.obs.exp.agctac<-log(n3.codon.pair$agctac/cps.agc.tac)
#cps.obs.exp.agctag<-log(n3.codon.pair$agctag/cps.agc.tag)
cps.obs.exp.agctat<-log(n3.codon.pair$agctat/cps.agc.tat)
cps.obs.exp.agctca<-log(n3.codon.pair$agctca/cps.agc.tca)
cps.obs.exp.agctcc<-log(n3.codon.pair$agctcc/cps.agc.tcc)
cps.obs.exp.agctcg<-log(n3.codon.pair$agctcg/cps.agc.tcg)
cps.obs.exp.agctct<-log(n3.codon.pair$agctct/cps.agc.tct)
#cps.obs.exp.agctga<-log(n3.codon.pair$agctga/cps.agc.tga)
cps.obs.exp.agctgc<-log(n3.codon.pair$agctgc/cps.agc.tgc)
cps.obs.exp.agctgg<-log(n3.codon.pair$agctgg/cps.agc.tgg)
cps.obs.exp.agctgt<-log(n3.codon.pair$agctgt/cps.agc.tgt)
cps.obs.exp.agctta<-log(n3.codon.pair$agctta/cps.agc.tta)
cps.obs.exp.agcttc<-log(n3.codon.pair$agcttc/cps.agc.ttc)
cps.obs.exp.agcttg<-log(n3.codon.pair$agcttg/cps.agc.ttg)
cps.obs.exp.agcttt<-log(n3.codon.pair$agcttt/cps.agc.ttt)









cps.obs.exp.aggaaa<-log(n3.codon.pair$aggaaa/cps.agg.aaa)
cps.obs.exp.aggaac<-log(n3.codon.pair$aggaac/cps.agg.aac)
cps.obs.exp.aggaag<-log(n3.codon.pair$aggaag/cps.agg.aag)
cps.obs.exp.aggaat<-log(n3.codon.pair$aggaat/cps.agg.aat)
cps.obs.exp.aggaca<-log(n3.codon.pair$aggaca/cps.agg.aca)
cps.obs.exp.aggacc<-log(n3.codon.pair$aggacc/cps.agg.acc)
cps.obs.exp.aggacg<-log(n3.codon.pair$aggacg/cps.agg.acg)
cps.obs.exp.aggact<-log(n3.codon.pair$aggact/cps.agg.act)
cps.obs.exp.aggaga<-log(n3.codon.pair$aggaga/cps.agg.aga)
cps.obs.exp.aggagc<-log(n3.codon.pair$aggagc/cps.agg.agc)
cps.obs.exp.aggagg<-log(n3.codon.pair$aggagg/cps.agg.agg)
cps.obs.exp.aggagt<-log(n3.codon.pair$aggagt/cps.agg.agt)
cps.obs.exp.aggata<-log(n3.codon.pair$aggata/cps.agg.ata)
cps.obs.exp.aggatc<-log(n3.codon.pair$aggatc/cps.agg.atc)
cps.obs.exp.aggatg<-log(n3.codon.pair$aggatg/cps.agg.atg)
cps.obs.exp.aggatt<-log(n3.codon.pair$aggatt/cps.agg.att)

cps.obs.exp.aggcaa<-log(n3.codon.pair$aggcaa/cps.agg.caa)
cps.obs.exp.aggcac<-log(n3.codon.pair$aggcac/cps.agg.cac)
cps.obs.exp.aggcag<-log(n3.codon.pair$aggcag/cps.agg.cag)
cps.obs.exp.aggcat<-log(n3.codon.pair$aggcat/cps.agg.cat)
cps.obs.exp.aggcca<-log(n3.codon.pair$aggcca/cps.agg.cca)
cps.obs.exp.aggccc<-log(n3.codon.pair$aggccc/cps.agg.ccc)
cps.obs.exp.aggccg<-log(n3.codon.pair$aggccg/cps.agg.ccg)
cps.obs.exp.aggcct<-log(n3.codon.pair$aggcct/cps.agg.cct)
cps.obs.exp.aggcga<-log(n3.codon.pair$aggcga/cps.agg.cga)
cps.obs.exp.aggcgc<-log(n3.codon.pair$aggcgc/cps.agg.cgc)
cps.obs.exp.aggcgg<-log(n3.codon.pair$aggcgg/cps.agg.cgg)
cps.obs.exp.aggcgt<-log(n3.codon.pair$aggcgt/cps.agg.cgt)
cps.obs.exp.aggcta<-log(n3.codon.pair$aggcta/cps.agg.cta)
cps.obs.exp.aggctc<-log(n3.codon.pair$aggctc/cps.agg.ctc)
cps.obs.exp.aggctg<-log(n3.codon.pair$aggctg/cps.agg.ctg)
cps.obs.exp.aggctt<-log(n3.codon.pair$aggctt/cps.agg.ctt)

cps.obs.exp.agggaa<-log(n3.codon.pair$agggaa/cps.agg.gaa)
cps.obs.exp.agggac<-log(n3.codon.pair$agggac/cps.agg.gac)
cps.obs.exp.agggag<-log(n3.codon.pair$agggag/cps.agg.gag)
cps.obs.exp.agggat<-log(n3.codon.pair$agggat/cps.agg.gat)
cps.obs.exp.agggca<-log(n3.codon.pair$agggca/cps.agg.gca)
cps.obs.exp.agggcc<-log(n3.codon.pair$agggcc/cps.agg.gcc)
cps.obs.exp.agggcg<-log(n3.codon.pair$agggcg/cps.agg.gcg)
cps.obs.exp.agggct<-log(n3.codon.pair$agggct/cps.agg.gct)
cps.obs.exp.agggga<-log(n3.codon.pair$agggga/cps.agg.gga)
cps.obs.exp.aggggc<-log(n3.codon.pair$aggggc/cps.agg.ggc)
cps.obs.exp.aggggg<-log(n3.codon.pair$aggggg/cps.agg.ggg)
cps.obs.exp.aggggt<-log(n3.codon.pair$aggggt/cps.agg.ggt)
cps.obs.exp.agggta<-log(n3.codon.pair$agggta/cps.agg.gta)
cps.obs.exp.agggtc<-log(n3.codon.pair$agggtc/cps.agg.gtc)
cps.obs.exp.agggtg<-log(n3.codon.pair$agggtg/cps.agg.gtg)
cps.obs.exp.agggtt<-log(n3.codon.pair$agggtt/cps.agg.gtt)

#cps.obs.exp.aggtaa<-log(n3.codon.pair$aggtaa/cps.agg.taa)
cps.obs.exp.aggtac<-log(n3.codon.pair$aggtac/cps.agg.tac)
#cps.obs.exp.aggtag<-log(n3.codon.pair$aggtag/cps.agg.tag)
cps.obs.exp.aggtat<-log(n3.codon.pair$aggtat/cps.agg.tat)
cps.obs.exp.aggtca<-log(n3.codon.pair$aggtca/cps.agg.tca)
cps.obs.exp.aggtcc<-log(n3.codon.pair$aggtcc/cps.agg.tcc)
cps.obs.exp.aggtcg<-log(n3.codon.pair$aggtcg/cps.agg.tcg)
cps.obs.exp.aggtct<-log(n3.codon.pair$aggtct/cps.agg.tct)
#cps.obs.exp.aggtga<-log(n3.codon.pair$aggtga/cps.agg.tga)
cps.obs.exp.aggtgc<-log(n3.codon.pair$aggtgc/cps.agg.tgc)
cps.obs.exp.aggtgg<-log(n3.codon.pair$aggtgg/cps.agg.tgg)
cps.obs.exp.aggtgt<-log(n3.codon.pair$aggtgt/cps.agg.tgt)
cps.obs.exp.aggtta<-log(n3.codon.pair$aggtta/cps.agg.tta)
cps.obs.exp.aggttc<-log(n3.codon.pair$aggttc/cps.agg.ttc)
cps.obs.exp.aggttg<-log(n3.codon.pair$aggttg/cps.agg.ttg)
cps.obs.exp.aggttt<-log(n3.codon.pair$aggttt/cps.agg.ttt)







cps.obs.exp.agtaaa<-log(n3.codon.pair$agtaaa/cps.agt.aaa)
cps.obs.exp.agtaac<-log(n3.codon.pair$agtaac/cps.agt.aac)
cps.obs.exp.agtaag<-log(n3.codon.pair$agtaag/cps.agt.aag)
cps.obs.exp.agtaat<-log(n3.codon.pair$agtaat/cps.agt.aat)
cps.obs.exp.agtaca<-log(n3.codon.pair$agtaca/cps.agt.aca)
cps.obs.exp.agtacc<-log(n3.codon.pair$agtacc/cps.agt.acc)
cps.obs.exp.agtacg<-log(n3.codon.pair$agtacg/cps.agt.acg)
cps.obs.exp.agtact<-log(n3.codon.pair$agtact/cps.agt.act)
cps.obs.exp.agtaga<-log(n3.codon.pair$agtaga/cps.agt.aga)
cps.obs.exp.agtagc<-log(n3.codon.pair$agtagc/cps.agt.agc)
cps.obs.exp.agtagg<-log(n3.codon.pair$agtagg/cps.agt.agg)
cps.obs.exp.agtagt<-log(n3.codon.pair$agtagt/cps.agt.agt)
cps.obs.exp.agtata<-log(n3.codon.pair$agtata/cps.agt.ata)
cps.obs.exp.agtatc<-log(n3.codon.pair$agtatc/cps.agt.atc)
cps.obs.exp.agtatg<-log(n3.codon.pair$agtatg/cps.agt.atg)
cps.obs.exp.agtatt<-log(n3.codon.pair$agtatt/cps.agt.att)

cps.obs.exp.agtcaa<-log(n3.codon.pair$agtcaa/cps.agt.caa)
cps.obs.exp.agtcac<-log(n3.codon.pair$agtcac/cps.agt.cac)
cps.obs.exp.agtcag<-log(n3.codon.pair$agtcag/cps.agt.cag)
cps.obs.exp.agtcat<-log(n3.codon.pair$agtcat/cps.agt.cat)
cps.obs.exp.agtcca<-log(n3.codon.pair$agtcca/cps.agt.cca)
cps.obs.exp.agtccc<-log(n3.codon.pair$agtccc/cps.agt.ccc)
cps.obs.exp.agtccg<-log(n3.codon.pair$agtccg/cps.agt.ccg)
cps.obs.exp.agtcct<-log(n3.codon.pair$agtcct/cps.agt.cct)
cps.obs.exp.agtcga<-log(n3.codon.pair$agtcga/cps.agt.cga)
cps.obs.exp.agtcgc<-log(n3.codon.pair$agtcgc/cps.agt.cgc)
cps.obs.exp.agtcgg<-log(n3.codon.pair$agtcgg/cps.agt.cgg)
cps.obs.exp.agtcgt<-log(n3.codon.pair$agtcgt/cps.agt.cgt)
cps.obs.exp.agtcta<-log(n3.codon.pair$agtcta/cps.agt.cta)
cps.obs.exp.agtctc<-log(n3.codon.pair$agtctc/cps.agt.ctc)
cps.obs.exp.agtctg<-log(n3.codon.pair$agtctg/cps.agt.ctg)
cps.obs.exp.agtctt<-log(n3.codon.pair$agtctt/cps.agt.ctt)

cps.obs.exp.agtgaa<-log(n3.codon.pair$agtgaa/cps.agt.gaa)
cps.obs.exp.agtgac<-log(n3.codon.pair$agtgac/cps.agt.gac)
cps.obs.exp.agtgag<-log(n3.codon.pair$agtgag/cps.agt.gag)
cps.obs.exp.agtgat<-log(n3.codon.pair$agtgat/cps.agt.gat)
cps.obs.exp.agtgca<-log(n3.codon.pair$agtgca/cps.agt.gca)
cps.obs.exp.agtgcc<-log(n3.codon.pair$agtgcc/cps.agt.gcc)
cps.obs.exp.agtgcg<-log(n3.codon.pair$agtgcg/cps.agt.gcg)
cps.obs.exp.agtgct<-log(n3.codon.pair$agtgct/cps.agt.gct)
cps.obs.exp.agtgga<-log(n3.codon.pair$agtgga/cps.agt.gga)
cps.obs.exp.agtggc<-log(n3.codon.pair$agtggc/cps.agt.ggc)
cps.obs.exp.agtggg<-log(n3.codon.pair$agtggg/cps.agt.ggg)
cps.obs.exp.agtggt<-log(n3.codon.pair$agtggt/cps.agt.ggt)
cps.obs.exp.agtgta<-log(n3.codon.pair$agtgta/cps.agt.gta)
cps.obs.exp.agtgtc<-log(n3.codon.pair$agtgtc/cps.agt.gtc)
cps.obs.exp.agtgtg<-log(n3.codon.pair$agtgtg/cps.agt.gtg)
cps.obs.exp.agtgtt<-log(n3.codon.pair$agtgtt/cps.agt.gtt)

#cps.obs.exp.agttaa<-log(n3.codon.pair$agttaa/cps.agt.taa)
cps.obs.exp.agttac<-log(n3.codon.pair$agttac/cps.agt.tac)
#cps.obs.exp.agttag<-log(n3.codon.pair$agttag/cps.agt.tag)
cps.obs.exp.agttat<-log(n3.codon.pair$agttat/cps.agt.tat)
cps.obs.exp.agttca<-log(n3.codon.pair$agttca/cps.agt.tca)
cps.obs.exp.agttcc<-log(n3.codon.pair$agttcc/cps.agt.tcc)
cps.obs.exp.agttcg<-log(n3.codon.pair$agttcg/cps.agt.tcg)
cps.obs.exp.agttct<-log(n3.codon.pair$agttct/cps.agt.tct)
#cps.obs.exp.agttga<-log(n3.codon.pair$agttga/cps.agt.tga)
cps.obs.exp.agttgc<-log(n3.codon.pair$agttgc/cps.agt.tgc)
cps.obs.exp.agttgg<-log(n3.codon.pair$agttgg/cps.agt.tgg)
cps.obs.exp.agttgt<-log(n3.codon.pair$agttgt/cps.agt.tgt)
cps.obs.exp.agttta<-log(n3.codon.pair$agttta/cps.agt.tta)
cps.obs.exp.agtttc<-log(n3.codon.pair$agtttc/cps.agt.ttc)
cps.obs.exp.agtttg<-log(n3.codon.pair$agtttg/cps.agt.ttg)
cps.obs.exp.agtttt<-log(n3.codon.pair$agtttt/cps.agt.ttt)




















cps.obs.exp.ataaaa<-log(n3.codon.pair$ataaaa/cps.ata.aaa)
cps.obs.exp.ataaac<-log(n3.codon.pair$ataaac/cps.ata.aac)
cps.obs.exp.ataaag<-log(n3.codon.pair$ataaag/cps.ata.aag)
cps.obs.exp.ataaat<-log(n3.codon.pair$ataaat/cps.ata.aat)
cps.obs.exp.ataaca<-log(n3.codon.pair$ataaca/cps.ata.aca)
cps.obs.exp.ataacc<-log(n3.codon.pair$ataacc/cps.ata.acc)
cps.obs.exp.ataacg<-log(n3.codon.pair$ataacg/cps.ata.acg)
cps.obs.exp.ataact<-log(n3.codon.pair$ataact/cps.ata.act)
cps.obs.exp.ataaga<-log(n3.codon.pair$ataaga/cps.ata.aga)
cps.obs.exp.ataagc<-log(n3.codon.pair$ataagc/cps.ata.agc)
cps.obs.exp.ataagg<-log(n3.codon.pair$ataagg/cps.ata.agg)
cps.obs.exp.ataagt<-log(n3.codon.pair$ataagt/cps.ata.agt)
cps.obs.exp.ataata<-log(n3.codon.pair$ataata/cps.ata.ata)
cps.obs.exp.ataatc<-log(n3.codon.pair$ataatc/cps.ata.atc)
cps.obs.exp.ataatg<-log(n3.codon.pair$ataatg/cps.ata.atg)
cps.obs.exp.ataatt<-log(n3.codon.pair$ataatt/cps.ata.att)

cps.obs.exp.atacaa<-log(n3.codon.pair$atacaa/cps.ata.caa)
cps.obs.exp.atacac<-log(n3.codon.pair$atacac/cps.ata.cac)
cps.obs.exp.atacag<-log(n3.codon.pair$atacag/cps.ata.cag)
cps.obs.exp.atacat<-log(n3.codon.pair$atacat/cps.ata.cat)
cps.obs.exp.atacca<-log(n3.codon.pair$atacca/cps.ata.cca)
cps.obs.exp.ataccc<-log(n3.codon.pair$ataccc/cps.ata.ccc)
cps.obs.exp.ataccg<-log(n3.codon.pair$ataccg/cps.ata.ccg)
cps.obs.exp.atacct<-log(n3.codon.pair$atacct/cps.ata.cct)
cps.obs.exp.atacga<-log(n3.codon.pair$atacga/cps.ata.cga)
cps.obs.exp.atacgc<-log(n3.codon.pair$atacgc/cps.ata.cgc)
cps.obs.exp.atacgg<-log(n3.codon.pair$atacgg/cps.ata.cgg)
cps.obs.exp.atacgt<-log(n3.codon.pair$atacgt/cps.ata.cgt)
cps.obs.exp.atacta<-log(n3.codon.pair$atacta/cps.ata.cta)
cps.obs.exp.atactc<-log(n3.codon.pair$atactc/cps.ata.ctc)
cps.obs.exp.atactg<-log(n3.codon.pair$atactg/cps.ata.ctg)
cps.obs.exp.atactt<-log(n3.codon.pair$atactt/cps.ata.ctt)

cps.obs.exp.atagaa<-log(n3.codon.pair$atagaa/cps.ata.gaa)
cps.obs.exp.atagac<-log(n3.codon.pair$atagac/cps.ata.gac)
cps.obs.exp.atagag<-log(n3.codon.pair$atagag/cps.ata.gag)
cps.obs.exp.atagat<-log(n3.codon.pair$atagat/cps.ata.gat)
cps.obs.exp.atagca<-log(n3.codon.pair$atagca/cps.ata.gca)
cps.obs.exp.atagcc<-log(n3.codon.pair$atagcc/cps.ata.gcc)
cps.obs.exp.atagcg<-log(n3.codon.pair$atagcg/cps.ata.gcg)
cps.obs.exp.atagct<-log(n3.codon.pair$atagct/cps.ata.gct)
cps.obs.exp.atagga<-log(n3.codon.pair$atagga/cps.ata.gga)
cps.obs.exp.ataggc<-log(n3.codon.pair$ataggc/cps.ata.ggc)
cps.obs.exp.ataggg<-log(n3.codon.pair$ataggg/cps.ata.ggg)
cps.obs.exp.ataggt<-log(n3.codon.pair$ataggt/cps.ata.ggt)
cps.obs.exp.atagta<-log(n3.codon.pair$atagta/cps.ata.gta)
cps.obs.exp.atagtc<-log(n3.codon.pair$atagtc/cps.ata.gtc)
cps.obs.exp.atagtg<-log(n3.codon.pair$atagtg/cps.ata.gtg)
cps.obs.exp.atagtt<-log(n3.codon.pair$atagtt/cps.ata.gtt)

#cps.obs.exp.atataa<-log(n3.codon.pair$atataa/cps.ata.taa)
cps.obs.exp.atatac<-log(n3.codon.pair$atatac/cps.ata.tac)
#cps.obs.exp.atatag<-log(n3.codon.pair$atatag/cps.ata.tag)
cps.obs.exp.atatat<-log(n3.codon.pair$atatat/cps.ata.tat)
cps.obs.exp.atatca<-log(n3.codon.pair$atatca/cps.ata.tca)
cps.obs.exp.atatcc<-log(n3.codon.pair$atatcc/cps.ata.tcc)
cps.obs.exp.atatcg<-log(n3.codon.pair$atatcg/cps.ata.tcg)
cps.obs.exp.atatct<-log(n3.codon.pair$atatct/cps.ata.tct)
#cps.obs.exp.atatga<-log(n3.codon.pair$atatga/cps.ata.tga)
cps.obs.exp.atatgc<-log(n3.codon.pair$atatgc/cps.ata.tgc)
cps.obs.exp.atatgg<-log(n3.codon.pair$atatgg/cps.ata.tgg)
cps.obs.exp.atatgt<-log(n3.codon.pair$atatgt/cps.ata.tgt)
cps.obs.exp.atatta<-log(n3.codon.pair$atatta/cps.ata.tta)
cps.obs.exp.atattc<-log(n3.codon.pair$atattc/cps.ata.ttc)
cps.obs.exp.atattg<-log(n3.codon.pair$atattg/cps.ata.ttg)
cps.obs.exp.atattt<-log(n3.codon.pair$atattt/cps.ata.ttt)







cps.obs.exp.atcaaa<-log(n3.codon.pair$atcaaa/cps.atc.aaa)
cps.obs.exp.atcaac<-log(n3.codon.pair$atcaac/cps.atc.aac)
cps.obs.exp.atcaag<-log(n3.codon.pair$atcaag/cps.atc.aag)
cps.obs.exp.atcaat<-log(n3.codon.pair$atcaat/cps.atc.aat)
cps.obs.exp.atcaca<-log(n3.codon.pair$atcaca/cps.atc.aca)
cps.obs.exp.atcacc<-log(n3.codon.pair$atcacc/cps.atc.acc)
cps.obs.exp.atcacg<-log(n3.codon.pair$atcacg/cps.atc.acg)
cps.obs.exp.atcact<-log(n3.codon.pair$atcact/cps.atc.act)
cps.obs.exp.atcaga<-log(n3.codon.pair$atcaga/cps.atc.aga)
cps.obs.exp.atcagc<-log(n3.codon.pair$atcagc/cps.atc.agc)
cps.obs.exp.atcagg<-log(n3.codon.pair$atcagg/cps.atc.agg)
cps.obs.exp.atcagt<-log(n3.codon.pair$atcagt/cps.atc.agt)
cps.obs.exp.atcata<-log(n3.codon.pair$atcata/cps.atc.ata)
cps.obs.exp.atcatc<-log(n3.codon.pair$atcatc/cps.atc.atc)
cps.obs.exp.atcatg<-log(n3.codon.pair$atcatg/cps.atc.atg)
cps.obs.exp.atcatt<-log(n3.codon.pair$atcatt/cps.atc.att)

cps.obs.exp.atccaa<-log(n3.codon.pair$atccaa/cps.atc.caa)
cps.obs.exp.atccac<-log(n3.codon.pair$atccac/cps.atc.cac)
cps.obs.exp.atccag<-log(n3.codon.pair$atccag/cps.atc.cag)
cps.obs.exp.atccat<-log(n3.codon.pair$atccat/cps.atc.cat)
cps.obs.exp.atccca<-log(n3.codon.pair$atccca/cps.atc.cca)
cps.obs.exp.atcccc<-log(n3.codon.pair$atcccc/cps.atc.ccc)
cps.obs.exp.atcccg<-log(n3.codon.pair$atcccg/cps.atc.ccg)
cps.obs.exp.atccct<-log(n3.codon.pair$atccct/cps.atc.cct)
cps.obs.exp.atccga<-log(n3.codon.pair$atccga/cps.atc.cga)
cps.obs.exp.atccgc<-log(n3.codon.pair$atccgc/cps.atc.cgc)
cps.obs.exp.atccgg<-log(n3.codon.pair$atccgg/cps.atc.cgg)
cps.obs.exp.atccgt<-log(n3.codon.pair$atccgt/cps.atc.cgt)
cps.obs.exp.atccta<-log(n3.codon.pair$atccta/cps.atc.cta)
cps.obs.exp.atcctc<-log(n3.codon.pair$atcctc/cps.atc.ctc)
cps.obs.exp.atcctg<-log(n3.codon.pair$atcctg/cps.atc.ctg)
cps.obs.exp.atcctt<-log(n3.codon.pair$atcctt/cps.atc.ctt)

cps.obs.exp.atcgaa<-log(n3.codon.pair$atcgaa/cps.atc.gaa)
cps.obs.exp.atcgac<-log(n3.codon.pair$atcgac/cps.atc.gac)
cps.obs.exp.atcgag<-log(n3.codon.pair$atcgag/cps.atc.gag)
cps.obs.exp.atcgat<-log(n3.codon.pair$atcgat/cps.atc.gat)
cps.obs.exp.atcgca<-log(n3.codon.pair$atcgca/cps.atc.gca)
cps.obs.exp.atcgcc<-log(n3.codon.pair$atcgcc/cps.atc.gcc)
cps.obs.exp.atcgcg<-log(n3.codon.pair$atcgcg/cps.atc.gcg)
cps.obs.exp.atcgct<-log(n3.codon.pair$atcgct/cps.atc.gct)
cps.obs.exp.atcgga<-log(n3.codon.pair$atcgga/cps.atc.gga)
cps.obs.exp.atcggc<-log(n3.codon.pair$atcggc/cps.atc.ggc)
cps.obs.exp.atcggg<-log(n3.codon.pair$atcggg/cps.atc.ggg)
cps.obs.exp.atcggt<-log(n3.codon.pair$atcggt/cps.atc.ggt)
cps.obs.exp.atcgta<-log(n3.codon.pair$atcgta/cps.atc.gta)
cps.obs.exp.atcgtc<-log(n3.codon.pair$atcgtc/cps.atc.gtc)
cps.obs.exp.atcgtg<-log(n3.codon.pair$atcgtg/cps.atc.gtg)
cps.obs.exp.atcgtt<-log(n3.codon.pair$atcgtt/cps.atc.gtt)

#cps.obs.exp.atctaa<-log(n3.codon.pair$atctaa/cps.atc.taa)
cps.obs.exp.atctac<-log(n3.codon.pair$atctac/cps.atc.tac)
#cps.obs.exp.atctag<-log(n3.codon.pair$atctag/cps.atc.tag)
cps.obs.exp.atctat<-log(n3.codon.pair$atctat/cps.atc.tat)
cps.obs.exp.atctca<-log(n3.codon.pair$atctca/cps.atc.tca)
cps.obs.exp.atctcc<-log(n3.codon.pair$atctcc/cps.atc.tcc)
cps.obs.exp.atctcg<-log(n3.codon.pair$atctcg/cps.atc.tcg)
cps.obs.exp.atctct<-log(n3.codon.pair$atctct/cps.atc.tct)
#cps.obs.exp.atctga<-log(n3.codon.pair$atctga/cps.atc.tga)
cps.obs.exp.atctgc<-log(n3.codon.pair$atctgc/cps.atc.tgc)
cps.obs.exp.atctgg<-log(n3.codon.pair$atctgg/cps.atc.tgg)
cps.obs.exp.atctgt<-log(n3.codon.pair$atctgt/cps.atc.tgt)
cps.obs.exp.atctta<-log(n3.codon.pair$atctta/cps.atc.tta)
cps.obs.exp.atcttc<-log(n3.codon.pair$atcttc/cps.atc.ttc)
cps.obs.exp.atcttg<-log(n3.codon.pair$atcttg/cps.atc.ttg)
cps.obs.exp.atcttt<-log(n3.codon.pair$atcttt/cps.atc.ttt)








cps.obs.exp.atgaaa<-log(n3.codon.pair$atgaaa/cps.atg.aaa)
cps.obs.exp.atgaac<-log(n3.codon.pair$atgaac/cps.atg.aac)
cps.obs.exp.atgaag<-log(n3.codon.pair$atgaag/cps.atg.aag)
cps.obs.exp.atgaat<-log(n3.codon.pair$atgaat/cps.atg.aat)
cps.obs.exp.atgaca<-log(n3.codon.pair$atgaca/cps.atg.aca)
cps.obs.exp.atgacc<-log(n3.codon.pair$atgacc/cps.atg.acc)
cps.obs.exp.atgacg<-log(n3.codon.pair$atgacg/cps.atg.acg)
cps.obs.exp.atgact<-log(n3.codon.pair$atgact/cps.atg.act)
cps.obs.exp.atgaga<-log(n3.codon.pair$atgaga/cps.atg.aga)
cps.obs.exp.atgagc<-log(n3.codon.pair$atgagc/cps.atg.agc)
cps.obs.exp.atgagg<-log(n3.codon.pair$atgagg/cps.atg.agg)
cps.obs.exp.atgagt<-log(n3.codon.pair$atgagt/cps.atg.agt)
cps.obs.exp.atgata<-log(n3.codon.pair$atgata/cps.atg.ata)
cps.obs.exp.atgatc<-log(n3.codon.pair$atgatc/cps.atg.atc)
cps.obs.exp.atgatg<-log(n3.codon.pair$atgatg/cps.atg.atg)
cps.obs.exp.atgatt<-log(n3.codon.pair$atgatt/cps.atg.att)

cps.obs.exp.atgcaa<-log(n3.codon.pair$atgcaa/cps.atg.caa)
cps.obs.exp.atgcac<-log(n3.codon.pair$atgcac/cps.atg.cac)
cps.obs.exp.atgcag<-log(n3.codon.pair$atgcag/cps.atg.cag)
cps.obs.exp.atgcat<-log(n3.codon.pair$atgcat/cps.atg.cat)
cps.obs.exp.atgcca<-log(n3.codon.pair$atgcca/cps.atg.cca)
cps.obs.exp.atgccc<-log(n3.codon.pair$atgccc/cps.atg.ccc)
cps.obs.exp.atgccg<-log(n3.codon.pair$atgccg/cps.atg.ccg)
cps.obs.exp.atgcct<-log(n3.codon.pair$atgcct/cps.atg.cct)
cps.obs.exp.atgcga<-log(n3.codon.pair$atgcga/cps.atg.cga)
cps.obs.exp.atgcgc<-log(n3.codon.pair$atgcgc/cps.atg.cgc)
cps.obs.exp.atgcgg<-log(n3.codon.pair$atgcgg/cps.atg.cgg)
cps.obs.exp.atgcgt<-log(n3.codon.pair$atgcgt/cps.atg.cgt)
cps.obs.exp.atgcta<-log(n3.codon.pair$atgcta/cps.atg.cta)
cps.obs.exp.atgctc<-log(n3.codon.pair$atgctc/cps.atg.ctc)
cps.obs.exp.atgctg<-log(n3.codon.pair$atgctg/cps.atg.ctg)
cps.obs.exp.atgctt<-log(n3.codon.pair$atgctt/cps.atg.ctt)

cps.obs.exp.atggaa<-log(n3.codon.pair$atggaa/cps.atg.gaa)
cps.obs.exp.atggac<-log(n3.codon.pair$atggac/cps.atg.gac)
cps.obs.exp.atggag<-log(n3.codon.pair$atggag/cps.atg.gag)
cps.obs.exp.atggat<-log(n3.codon.pair$atggat/cps.atg.gat)
cps.obs.exp.atggca<-log(n3.codon.pair$atggca/cps.atg.gca)
cps.obs.exp.atggcc<-log(n3.codon.pair$atggcc/cps.atg.gcc)
cps.obs.exp.atggcg<-log(n3.codon.pair$atggcg/cps.atg.gcg)
cps.obs.exp.atggct<-log(n3.codon.pair$atggct/cps.atg.gct)
cps.obs.exp.atggga<-log(n3.codon.pair$atggga/cps.atg.gga)
cps.obs.exp.atgggc<-log(n3.codon.pair$atgggc/cps.atg.ggc)
cps.obs.exp.atgggg<-log(n3.codon.pair$atgggg/cps.atg.ggg)
cps.obs.exp.atgggt<-log(n3.codon.pair$atgggt/cps.atg.ggt)
cps.obs.exp.atggta<-log(n3.codon.pair$atggta/cps.atg.gta)
cps.obs.exp.atggtc<-log(n3.codon.pair$atggtc/cps.atg.gtc)
cps.obs.exp.atggtg<-log(n3.codon.pair$atggtg/cps.atg.gtg)
cps.obs.exp.atggtt<-log(n3.codon.pair$atggtt/cps.atg.gtt)

#cps.obs.exp.atgtaa<-log(n3.codon.pair$atgtaa/cps.atg.taa)
cps.obs.exp.atgtac<-log(n3.codon.pair$atgtac/cps.atg.tac)
#cps.obs.exp.atgtag<-log(n3.codon.pair$atgtag/cps.atg.tag)
cps.obs.exp.atgtat<-log(n3.codon.pair$atgtat/cps.atg.tat)
cps.obs.exp.atgtca<-log(n3.codon.pair$atgtca/cps.atg.tca)
cps.obs.exp.atgtcc<-log(n3.codon.pair$atgtcc/cps.atg.tcc)
cps.obs.exp.atgtcg<-log(n3.codon.pair$atgtcg/cps.atg.tcg)
cps.obs.exp.atgtct<-log(n3.codon.pair$atgtct/cps.atg.tct)
#cps.obs.exp.atgtga<-log(n3.codon.pair$atgtga/cps.atg.tga)
cps.obs.exp.atgtgc<-log(n3.codon.pair$atgtgc/cps.atg.tgc)
cps.obs.exp.atgtgg<-log(n3.codon.pair$atgtgg/cps.atg.tgg)
cps.obs.exp.atgtgt<-log(n3.codon.pair$atgtgt/cps.atg.tgt)
cps.obs.exp.atgtta<-log(n3.codon.pair$atgtta/cps.atg.tta)
cps.obs.exp.atgttc<-log(n3.codon.pair$atgttc/cps.atg.ttc)
cps.obs.exp.atgttg<-log(n3.codon.pair$atgttg/cps.atg.ttg)
cps.obs.exp.atgttt<-log(n3.codon.pair$atgttt/cps.atg.ttt)







cps.obs.exp.attaaa<-log(n3.codon.pair$attaaa/cps.att.aaa)
cps.obs.exp.attaac<-log(n3.codon.pair$attaac/cps.att.aac)
cps.obs.exp.attaag<-log(n3.codon.pair$attaag/cps.att.aag)
cps.obs.exp.attaat<-log(n3.codon.pair$attaat/cps.att.aat)
cps.obs.exp.attaca<-log(n3.codon.pair$attaca/cps.att.aca)
cps.obs.exp.attacc<-log(n3.codon.pair$attacc/cps.att.acc)
cps.obs.exp.attacg<-log(n3.codon.pair$attacg/cps.att.acg)
cps.obs.exp.attact<-log(n3.codon.pair$attact/cps.att.act)
cps.obs.exp.attaga<-log(n3.codon.pair$attaga/cps.att.aga)
cps.obs.exp.attagc<-log(n3.codon.pair$attagc/cps.att.agc)
cps.obs.exp.attagg<-log(n3.codon.pair$attagg/cps.att.agg)
cps.obs.exp.attagt<-log(n3.codon.pair$attagt/cps.att.agt)
cps.obs.exp.attata<-log(n3.codon.pair$attata/cps.att.ata)
cps.obs.exp.attatc<-log(n3.codon.pair$attatc/cps.att.atc)
cps.obs.exp.attatg<-log(n3.codon.pair$attatg/cps.att.atg)
cps.obs.exp.attatt<-log(n3.codon.pair$attatt/cps.att.att)

cps.obs.exp.attcaa<-log(n3.codon.pair$attcaa/cps.att.caa)
cps.obs.exp.attcac<-log(n3.codon.pair$attcac/cps.att.cac)
cps.obs.exp.attcag<-log(n3.codon.pair$attcag/cps.att.cag)
cps.obs.exp.attcat<-log(n3.codon.pair$attcat/cps.att.cat)
cps.obs.exp.attcca<-log(n3.codon.pair$attcca/cps.att.cca)
cps.obs.exp.attccc<-log(n3.codon.pair$attccc/cps.att.ccc)
cps.obs.exp.attccg<-log(n3.codon.pair$attccg/cps.att.ccg)
cps.obs.exp.attcct<-log(n3.codon.pair$attcct/cps.att.cct)
cps.obs.exp.attcga<-log(n3.codon.pair$attcga/cps.att.cga)
cps.obs.exp.attcgc<-log(n3.codon.pair$attcgc/cps.att.cgc)
cps.obs.exp.attcgg<-log(n3.codon.pair$attcgg/cps.att.cgg)
cps.obs.exp.attcgt<-log(n3.codon.pair$attcgt/cps.att.cgt)
cps.obs.exp.attcta<-log(n3.codon.pair$attcta/cps.att.cta)
cps.obs.exp.attctc<-log(n3.codon.pair$attctc/cps.att.ctc)
cps.obs.exp.attctg<-log(n3.codon.pair$attctg/cps.att.ctg)
cps.obs.exp.attctt<-log(n3.codon.pair$attctt/cps.att.ctt)

cps.obs.exp.attgaa<-log(n3.codon.pair$attgaa/cps.att.gaa)
cps.obs.exp.attgac<-log(n3.codon.pair$attgac/cps.att.gac)
cps.obs.exp.attgag<-log(n3.codon.pair$attgag/cps.att.gag)
cps.obs.exp.attgat<-log(n3.codon.pair$attgat/cps.att.gat)
cps.obs.exp.attgca<-log(n3.codon.pair$attgca/cps.att.gca)
cps.obs.exp.attgcc<-log(n3.codon.pair$attgcc/cps.att.gcc)
cps.obs.exp.attgcg<-log(n3.codon.pair$attgcg/cps.att.gcg)
cps.obs.exp.attgct<-log(n3.codon.pair$attgct/cps.att.gct)
cps.obs.exp.attgga<-log(n3.codon.pair$attgga/cps.att.gga)
cps.obs.exp.attggc<-log(n3.codon.pair$attggc/cps.att.ggc)
cps.obs.exp.attggg<-log(n3.codon.pair$attggg/cps.att.ggg)
cps.obs.exp.attggt<-log(n3.codon.pair$attggt/cps.att.ggt)
cps.obs.exp.attgta<-log(n3.codon.pair$attgta/cps.att.gta)
cps.obs.exp.attgtc<-log(n3.codon.pair$attgtc/cps.att.gtc)
cps.obs.exp.attgtg<-log(n3.codon.pair$attgtg/cps.att.gtg)
cps.obs.exp.attgtt<-log(n3.codon.pair$attgtt/cps.att.gtt)

#cps.obs.exp.atttaa<-log(n3.codon.pair$atttaa/cps.att.taa)
cps.obs.exp.atttac<-log(n3.codon.pair$atttac/cps.att.tac)
#cps.obs.exp.atttag<-log(n3.codon.pair$atttag/cps.att.tag)
cps.obs.exp.atttat<-log(n3.codon.pair$atttat/cps.att.tat)
cps.obs.exp.atttca<-log(n3.codon.pair$atttca/cps.att.tca)
cps.obs.exp.atttcc<-log(n3.codon.pair$atttcc/cps.att.tcc)
cps.obs.exp.atttcg<-log(n3.codon.pair$atttcg/cps.att.tcg)
cps.obs.exp.atttct<-log(n3.codon.pair$atttct/cps.att.tct)
#cps.obs.exp.atttga<-log(n3.codon.pair$atttga/cps.att.tga)
cps.obs.exp.atttgc<-log(n3.codon.pair$atttgc/cps.att.tgc)
cps.obs.exp.atttgg<-log(n3.codon.pair$atttgg/cps.att.tgg)
cps.obs.exp.atttgt<-log(n3.codon.pair$atttgt/cps.att.tgt)
cps.obs.exp.atttta<-log(n3.codon.pair$atttta/cps.att.tta)
cps.obs.exp.attttc<-log(n3.codon.pair$attttc/cps.att.ttc)
cps.obs.exp.attttg<-log(n3.codon.pair$attttg/cps.att.ttg)
cps.obs.exp.attttt<-log(n3.codon.pair$attttt/cps.att.ttt)



















cps.obs.exp.caaaaa<-log(n3.codon.pair$caaaaa/cps.caa.aaa)
cps.obs.exp.caaaac<-log(n3.codon.pair$caaaac/cps.caa.aac)
cps.obs.exp.caaaag<-log(n3.codon.pair$caaaag/cps.caa.aag)
cps.obs.exp.caaaat<-log(n3.codon.pair$caaaat/cps.caa.aat)
cps.obs.exp.caaaca<-log(n3.codon.pair$caaaca/cps.caa.aca)
cps.obs.exp.caaacc<-log(n3.codon.pair$caaacc/cps.caa.acc)
cps.obs.exp.caaacg<-log(n3.codon.pair$caaacg/cps.caa.acg)
cps.obs.exp.caaact<-log(n3.codon.pair$caaact/cps.caa.act)
cps.obs.exp.caaaga<-log(n3.codon.pair$caaaga/cps.caa.aga)
cps.obs.exp.caaagc<-log(n3.codon.pair$caaagc/cps.caa.agc)
cps.obs.exp.caaagg<-log(n3.codon.pair$caaagg/cps.caa.agg)
cps.obs.exp.caaagt<-log(n3.codon.pair$caaagt/cps.caa.agt)
cps.obs.exp.caaata<-log(n3.codon.pair$caaata/cps.caa.ata)
cps.obs.exp.caaatc<-log(n3.codon.pair$caaatc/cps.caa.atc)
cps.obs.exp.caaatg<-log(n3.codon.pair$caaatg/cps.caa.atg)
cps.obs.exp.caaatt<-log(n3.codon.pair$caaatt/cps.caa.att)

cps.obs.exp.caacaa<-log(n3.codon.pair$caacaa/cps.caa.caa)
cps.obs.exp.caacac<-log(n3.codon.pair$caacac/cps.caa.cac)
cps.obs.exp.caacag<-log(n3.codon.pair$caacag/cps.caa.cag)
cps.obs.exp.caacat<-log(n3.codon.pair$caacat/cps.caa.cat)
cps.obs.exp.caacca<-log(n3.codon.pair$caacca/cps.caa.cca)
cps.obs.exp.caaccc<-log(n3.codon.pair$caaccc/cps.caa.ccc)
cps.obs.exp.caaccg<-log(n3.codon.pair$caaccg/cps.caa.ccg)
cps.obs.exp.caacct<-log(n3.codon.pair$caacct/cps.caa.cct)
cps.obs.exp.caacga<-log(n3.codon.pair$caacga/cps.caa.cga)
cps.obs.exp.caacgc<-log(n3.codon.pair$caacgc/cps.caa.cgc)
cps.obs.exp.caacgg<-log(n3.codon.pair$caacgg/cps.caa.cgg)
cps.obs.exp.caacgt<-log(n3.codon.pair$caacgt/cps.caa.cgt)
cps.obs.exp.caacta<-log(n3.codon.pair$caacta/cps.caa.cta)
cps.obs.exp.caactc<-log(n3.codon.pair$caactc/cps.caa.ctc)
cps.obs.exp.caactg<-log(n3.codon.pair$caactg/cps.caa.ctg)
cps.obs.exp.caactt<-log(n3.codon.pair$caactt/cps.caa.ctt)

cps.obs.exp.caagaa<-log(n3.codon.pair$caagaa/cps.caa.gaa)
cps.obs.exp.caagac<-log(n3.codon.pair$caagac/cps.caa.gac)
cps.obs.exp.caagag<-log(n3.codon.pair$caagag/cps.caa.gag)
cps.obs.exp.caagat<-log(n3.codon.pair$caagat/cps.caa.gat)
cps.obs.exp.caagca<-log(n3.codon.pair$caagca/cps.caa.gca)
cps.obs.exp.caagcc<-log(n3.codon.pair$caagcc/cps.caa.gcc)
cps.obs.exp.caagcg<-log(n3.codon.pair$caagcg/cps.caa.gcg)
cps.obs.exp.caagct<-log(n3.codon.pair$caagct/cps.caa.gct)
cps.obs.exp.caagga<-log(n3.codon.pair$caagga/cps.caa.gga)
cps.obs.exp.caaggc<-log(n3.codon.pair$caaggc/cps.caa.ggc)
cps.obs.exp.caaggg<-log(n3.codon.pair$caaggg/cps.caa.ggg)
cps.obs.exp.caaggt<-log(n3.codon.pair$caaggt/cps.caa.ggt)
cps.obs.exp.caagta<-log(n3.codon.pair$caagta/cps.caa.gta)
cps.obs.exp.caagtc<-log(n3.codon.pair$caagtc/cps.caa.gtc)
cps.obs.exp.caagtg<-log(n3.codon.pair$caagtg/cps.caa.gtg)
cps.obs.exp.caagtt<-log(n3.codon.pair$caagtt/cps.caa.gtt)

#cps.obs.exp.caataa<-log(n3.codon.pair$caataa/cps.caa.taa)
cps.obs.exp.caatac<-log(n3.codon.pair$caatac/cps.caa.tac)
#cps.obs.exp.caatag<-log(n3.codon.pair$caatag/cps.caa.tag)
cps.obs.exp.caatat<-log(n3.codon.pair$caatat/cps.caa.tat)
cps.obs.exp.caatca<-log(n3.codon.pair$caatca/cps.caa.tca)
cps.obs.exp.caatcc<-log(n3.codon.pair$caatcc/cps.caa.tcc)
cps.obs.exp.caatcg<-log(n3.codon.pair$caatcg/cps.caa.tcg)
cps.obs.exp.caatct<-log(n3.codon.pair$caatct/cps.caa.tct)
#cps.obs.exp.caatga<-log(n3.codon.pair$caatga/cps.caa.tga)
cps.obs.exp.caatgc<-log(n3.codon.pair$caatgc/cps.caa.tgc)
cps.obs.exp.caatgg<-log(n3.codon.pair$caatgg/cps.caa.tgg)
cps.obs.exp.caatgt<-log(n3.codon.pair$caatgt/cps.caa.tgt)
cps.obs.exp.caatta<-log(n3.codon.pair$caatta/cps.caa.tta)
cps.obs.exp.caattc<-log(n3.codon.pair$caattc/cps.caa.ttc)
cps.obs.exp.caattg<-log(n3.codon.pair$caattg/cps.caa.ttg)
cps.obs.exp.caattt<-log(n3.codon.pair$caattt/cps.caa.ttt)







cps.obs.exp.cacaaa<-log(n3.codon.pair$cacaaa/cps.cac.aaa)
cps.obs.exp.cacaac<-log(n3.codon.pair$cacaac/cps.cac.aac)
cps.obs.exp.cacaag<-log(n3.codon.pair$cacaag/cps.cac.aag)
cps.obs.exp.cacaat<-log(n3.codon.pair$cacaat/cps.cac.aat)
cps.obs.exp.cacaca<-log(n3.codon.pair$cacaca/cps.cac.aca)
cps.obs.exp.cacacc<-log(n3.codon.pair$cacacc/cps.cac.acc)
cps.obs.exp.cacacg<-log(n3.codon.pair$cacacg/cps.cac.acg)
cps.obs.exp.cacact<-log(n3.codon.pair$cacact/cps.cac.act)
cps.obs.exp.cacaga<-log(n3.codon.pair$cacaga/cps.cac.aga)
cps.obs.exp.cacagc<-log(n3.codon.pair$cacagc/cps.cac.agc)
cps.obs.exp.cacagg<-log(n3.codon.pair$cacagg/cps.cac.agg)
cps.obs.exp.cacagt<-log(n3.codon.pair$cacagt/cps.cac.agt)
cps.obs.exp.cacata<-log(n3.codon.pair$cacata/cps.cac.ata)
cps.obs.exp.cacatc<-log(n3.codon.pair$cacatc/cps.cac.atc)
cps.obs.exp.cacatg<-log(n3.codon.pair$cacatg/cps.cac.atg)
cps.obs.exp.cacatt<-log(n3.codon.pair$cacatt/cps.cac.att)

cps.obs.exp.caccaa<-log(n3.codon.pair$caccaa/cps.cac.caa)
cps.obs.exp.caccac<-log(n3.codon.pair$caccac/cps.cac.cac)
cps.obs.exp.caccag<-log(n3.codon.pair$caccag/cps.cac.cag)
cps.obs.exp.caccat<-log(n3.codon.pair$caccat/cps.cac.cat)
cps.obs.exp.caccca<-log(n3.codon.pair$caccca/cps.cac.cca)
cps.obs.exp.cacccc<-log(n3.codon.pair$cacccc/cps.cac.ccc)
cps.obs.exp.cacccg<-log(n3.codon.pair$cacccg/cps.cac.ccg)
cps.obs.exp.caccct<-log(n3.codon.pair$caccct/cps.cac.cct)
cps.obs.exp.caccga<-log(n3.codon.pair$caccga/cps.cac.cga)
cps.obs.exp.caccgc<-log(n3.codon.pair$caccgc/cps.cac.cgc)
cps.obs.exp.caccgg<-log(n3.codon.pair$caccgg/cps.cac.cgg)
cps.obs.exp.caccgt<-log(n3.codon.pair$caccgt/cps.cac.cgt)
cps.obs.exp.caccta<-log(n3.codon.pair$caccta/cps.cac.cta)
cps.obs.exp.cacctc<-log(n3.codon.pair$cacctc/cps.cac.ctc)
cps.obs.exp.cacctg<-log(n3.codon.pair$cacctg/cps.cac.ctg)
cps.obs.exp.cacctt<-log(n3.codon.pair$cacctt/cps.cac.ctt)

cps.obs.exp.cacgaa<-log(n3.codon.pair$cacgaa/cps.cac.gaa)
cps.obs.exp.cacgac<-log(n3.codon.pair$cacgac/cps.cac.gac)
cps.obs.exp.cacgag<-log(n3.codon.pair$cacgag/cps.cac.gag)
cps.obs.exp.cacgat<-log(n3.codon.pair$cacgat/cps.cac.gat)
cps.obs.exp.cacgca<-log(n3.codon.pair$cacgca/cps.cac.gca)
cps.obs.exp.cacgcc<-log(n3.codon.pair$cacgcc/cps.cac.gcc)
cps.obs.exp.cacgcg<-log(n3.codon.pair$cacgcg/cps.cac.gcg)
cps.obs.exp.cacgct<-log(n3.codon.pair$cacgct/cps.cac.gct)
cps.obs.exp.cacgga<-log(n3.codon.pair$cacgga/cps.cac.gga)
cps.obs.exp.cacggc<-log(n3.codon.pair$cacggc/cps.cac.ggc)
cps.obs.exp.cacggg<-log(n3.codon.pair$cacggg/cps.cac.ggg)
cps.obs.exp.cacggt<-log(n3.codon.pair$cacggt/cps.cac.ggt)
cps.obs.exp.cacgta<-log(n3.codon.pair$cacgta/cps.cac.gta)
cps.obs.exp.cacgtc<-log(n3.codon.pair$cacgtc/cps.cac.gtc)
cps.obs.exp.cacgtg<-log(n3.codon.pair$cacgtg/cps.cac.gtg)
cps.obs.exp.cacgtt<-log(n3.codon.pair$cacgtt/cps.cac.gtt)

#cps.obs.exp.cactaa<-log(n3.codon.pair$cactaa/cps.cac.taa)
cps.obs.exp.cactac<-log(n3.codon.pair$cactac/cps.cac.tac)
#cps.obs.exp.cactag<-log(n3.codon.pair$cactag/cps.cac.tag)
cps.obs.exp.cactat<-log(n3.codon.pair$cactat/cps.cac.tat)
cps.obs.exp.cactca<-log(n3.codon.pair$cactca/cps.cac.tca)
cps.obs.exp.cactcc<-log(n3.codon.pair$cactcc/cps.cac.tcc)
cps.obs.exp.cactcg<-log(n3.codon.pair$cactcg/cps.cac.tcg)
cps.obs.exp.cactct<-log(n3.codon.pair$cactct/cps.cac.tct)
#cps.obs.exp.cactga<-log(n3.codon.pair$cactga/cps.cac.tga)
cps.obs.exp.cactgc<-log(n3.codon.pair$cactgc/cps.cac.tgc)
cps.obs.exp.cactgg<-log(n3.codon.pair$cactgg/cps.cac.tgg)
cps.obs.exp.cactgt<-log(n3.codon.pair$cactgt/cps.cac.tgt)
cps.obs.exp.cactta<-log(n3.codon.pair$cactta/cps.cac.tta)
cps.obs.exp.cacttc<-log(n3.codon.pair$cacttc/cps.cac.ttc)
cps.obs.exp.cacttg<-log(n3.codon.pair$cacttg/cps.cac.ttg)
cps.obs.exp.cacttt<-log(n3.codon.pair$cacttt/cps.cac.ttt)









cps.obs.exp.cagaaa<-log(n3.codon.pair$cagaaa/cps.cag.aaa)
cps.obs.exp.cagaac<-log(n3.codon.pair$cagaac/cps.cag.aac)
cps.obs.exp.cagaag<-log(n3.codon.pair$cagaag/cps.cag.aag)
cps.obs.exp.cagaat<-log(n3.codon.pair$cagaat/cps.cag.aat)
cps.obs.exp.cagaca<-log(n3.codon.pair$cagaca/cps.cag.aca)
cps.obs.exp.cagacc<-log(n3.codon.pair$cagacc/cps.cag.acc)
cps.obs.exp.cagacg<-log(n3.codon.pair$cagacg/cps.cag.acg)
cps.obs.exp.cagact<-log(n3.codon.pair$cagact/cps.cag.act)
cps.obs.exp.cagaga<-log(n3.codon.pair$cagaga/cps.cag.aga)
cps.obs.exp.cagagc<-log(n3.codon.pair$cagagc/cps.cag.agc)
cps.obs.exp.cagagg<-log(n3.codon.pair$cagagg/cps.cag.agg)
cps.obs.exp.cagagt<-log(n3.codon.pair$cagagt/cps.cag.agt)
cps.obs.exp.cagata<-log(n3.codon.pair$cagata/cps.cag.ata)
cps.obs.exp.cagatc<-log(n3.codon.pair$cagatc/cps.cag.atc)
cps.obs.exp.cagatg<-log(n3.codon.pair$cagatg/cps.cag.atg)
cps.obs.exp.cagatt<-log(n3.codon.pair$cagatt/cps.cag.att)

cps.obs.exp.cagcaa<-log(n3.codon.pair$cagcaa/cps.cag.caa)
cps.obs.exp.cagcac<-log(n3.codon.pair$cagcac/cps.cag.cac)
cps.obs.exp.cagcag<-log(n3.codon.pair$cagcag/cps.cag.cag)
cps.obs.exp.cagcat<-log(n3.codon.pair$cagcat/cps.cag.cat)
cps.obs.exp.cagcca<-log(n3.codon.pair$cagcca/cps.cag.cca)
cps.obs.exp.cagccc<-log(n3.codon.pair$cagccc/cps.cag.ccc)
cps.obs.exp.cagccg<-log(n3.codon.pair$cagccg/cps.cag.ccg)
cps.obs.exp.cagcct<-log(n3.codon.pair$cagcct/cps.cag.cct)
cps.obs.exp.cagcga<-log(n3.codon.pair$cagcga/cps.cag.cga)
cps.obs.exp.cagcgc<-log(n3.codon.pair$cagcgc/cps.cag.cgc)
cps.obs.exp.cagcgg<-log(n3.codon.pair$cagcgg/cps.cag.cgg)
cps.obs.exp.cagcgt<-log(n3.codon.pair$cagcgt/cps.cag.cgt)
cps.obs.exp.cagcta<-log(n3.codon.pair$cagcta/cps.cag.cta)
cps.obs.exp.cagctc<-log(n3.codon.pair$cagctc/cps.cag.ctc)
cps.obs.exp.cagctg<-log(n3.codon.pair$cagctg/cps.cag.ctg)
cps.obs.exp.cagctt<-log(n3.codon.pair$cagctt/cps.cag.ctt)

cps.obs.exp.caggaa<-log(n3.codon.pair$caggaa/cps.cag.gaa)
cps.obs.exp.caggac<-log(n3.codon.pair$caggac/cps.cag.gac)
cps.obs.exp.caggag<-log(n3.codon.pair$caggag/cps.cag.gag)
cps.obs.exp.caggat<-log(n3.codon.pair$caggat/cps.cag.gat)
cps.obs.exp.caggca<-log(n3.codon.pair$caggca/cps.cag.gca)
cps.obs.exp.caggcc<-log(n3.codon.pair$caggcc/cps.cag.gcc)
cps.obs.exp.caggcg<-log(n3.codon.pair$caggcg/cps.cag.gcg)
cps.obs.exp.caggct<-log(n3.codon.pair$caggct/cps.cag.gct)
cps.obs.exp.caggga<-log(n3.codon.pair$caggga/cps.cag.gga)
cps.obs.exp.cagggc<-log(n3.codon.pair$cagggc/cps.cag.ggc)
cps.obs.exp.cagggg<-log(n3.codon.pair$cagggg/cps.cag.ggg)
cps.obs.exp.cagggt<-log(n3.codon.pair$cagggt/cps.cag.ggt)
cps.obs.exp.caggta<-log(n3.codon.pair$caggta/cps.cag.gta)
cps.obs.exp.caggtc<-log(n3.codon.pair$caggtc/cps.cag.gtc)
cps.obs.exp.caggtg<-log(n3.codon.pair$caggtg/cps.cag.gtg)
cps.obs.exp.caggtt<-log(n3.codon.pair$caggtt/cps.cag.gtt)

#cps.obs.exp.cagtaa<-log(n3.codon.pair$cagtaa/cps.cag.taa)
cps.obs.exp.cagtac<-log(n3.codon.pair$cagtac/cps.cag.tac)
#cps.obs.exp.cagtag<-log(n3.codon.pair$cagtag/cps.cag.tag)
cps.obs.exp.cagtat<-log(n3.codon.pair$cagtat/cps.cag.tat)
cps.obs.exp.cagtca<-log(n3.codon.pair$cagtca/cps.cag.tca)
cps.obs.exp.cagtcc<-log(n3.codon.pair$cagtcc/cps.cag.tcc)
cps.obs.exp.cagtcg<-log(n3.codon.pair$cagtcg/cps.cag.tcg)
cps.obs.exp.cagtct<-log(n3.codon.pair$cagtct/cps.cag.tct)
#cps.obs.exp.cagtga<-log(n3.codon.pair$cagtga/cps.cag.tga)
cps.obs.exp.cagtgc<-log(n3.codon.pair$cagtgc/cps.cag.tgc)
cps.obs.exp.cagtgg<-log(n3.codon.pair$cagtgg/cps.cag.tgg)
cps.obs.exp.cagtgt<-log(n3.codon.pair$cagtgt/cps.cag.tgt)
cps.obs.exp.cagtta<-log(n3.codon.pair$cagtta/cps.cag.tta)
cps.obs.exp.cagttc<-log(n3.codon.pair$cagttc/cps.cag.ttc)
cps.obs.exp.cagttg<-log(n3.codon.pair$cagttg/cps.cag.ttg)
cps.obs.exp.cagttt<-log(n3.codon.pair$cagttt/cps.cag.ttt)








cps.obs.exp.cataaa<-log(n3.codon.pair$cataaa/cps.cat.aaa)
cps.obs.exp.cataac<-log(n3.codon.pair$cataac/cps.cat.aac)
cps.obs.exp.cataag<-log(n3.codon.pair$cataag/cps.cat.aag)
cps.obs.exp.cataat<-log(n3.codon.pair$cataat/cps.cat.aat)
cps.obs.exp.cataca<-log(n3.codon.pair$cataca/cps.cat.aca)
cps.obs.exp.catacc<-log(n3.codon.pair$catacc/cps.cat.acc)
cps.obs.exp.catacg<-log(n3.codon.pair$catacg/cps.cat.acg)
cps.obs.exp.catact<-log(n3.codon.pair$catact/cps.cat.act)
cps.obs.exp.cataga<-log(n3.codon.pair$cataga/cps.cat.aga)
cps.obs.exp.catagc<-log(n3.codon.pair$catagc/cps.cat.agc)
cps.obs.exp.catagg<-log(n3.codon.pair$catagg/cps.cat.agg)
cps.obs.exp.catagt<-log(n3.codon.pair$catagt/cps.cat.agt)
cps.obs.exp.catata<-log(n3.codon.pair$catata/cps.cat.ata)
cps.obs.exp.catatc<-log(n3.codon.pair$catatc/cps.cat.atc)
cps.obs.exp.catatg<-log(n3.codon.pair$catatg/cps.cat.atg)
cps.obs.exp.catatt<-log(n3.codon.pair$catatt/cps.cat.att)

cps.obs.exp.catcaa<-log(n3.codon.pair$catcaa/cps.cat.caa)
cps.obs.exp.catcac<-log(n3.codon.pair$catcac/cps.cat.cac)
cps.obs.exp.catcag<-log(n3.codon.pair$catcag/cps.cat.cag)
cps.obs.exp.catcat<-log(n3.codon.pair$catcat/cps.cat.cat)
cps.obs.exp.catcca<-log(n3.codon.pair$catcca/cps.cat.cca)
cps.obs.exp.catccc<-log(n3.codon.pair$catccc/cps.cat.ccc)
cps.obs.exp.catccg<-log(n3.codon.pair$catccg/cps.cat.ccg)
cps.obs.exp.catcct<-log(n3.codon.pair$catcct/cps.cat.cct)
cps.obs.exp.catcga<-log(n3.codon.pair$catcga/cps.cat.cga)
cps.obs.exp.catcgc<-log(n3.codon.pair$catcgc/cps.cat.cgc)
cps.obs.exp.catcgg<-log(n3.codon.pair$catcgg/cps.cat.cgg)
cps.obs.exp.catcgt<-log(n3.codon.pair$catcgt/cps.cat.cgt)
cps.obs.exp.catcta<-log(n3.codon.pair$catcta/cps.cat.cta)
cps.obs.exp.catctc<-log(n3.codon.pair$catctc/cps.cat.ctc)
cps.obs.exp.catctg<-log(n3.codon.pair$catctg/cps.cat.ctg)
cps.obs.exp.catctt<-log(n3.codon.pair$catctt/cps.cat.ctt)

cps.obs.exp.catgaa<-log(n3.codon.pair$catgaa/cps.cat.gaa)
cps.obs.exp.catgac<-log(n3.codon.pair$catgac/cps.cat.gac)
cps.obs.exp.catgag<-log(n3.codon.pair$catgag/cps.cat.gag)
cps.obs.exp.catgat<-log(n3.codon.pair$catgat/cps.cat.gat)
cps.obs.exp.catgca<-log(n3.codon.pair$catgca/cps.cat.gca)
cps.obs.exp.catgcc<-log(n3.codon.pair$catgcc/cps.cat.gcc)
cps.obs.exp.catgcg<-log(n3.codon.pair$catgcg/cps.cat.gcg)
cps.obs.exp.catgct<-log(n3.codon.pair$catgct/cps.cat.gct)
cps.obs.exp.catgga<-log(n3.codon.pair$catgga/cps.cat.gga)
cps.obs.exp.catggc<-log(n3.codon.pair$catggc/cps.cat.ggc)
cps.obs.exp.catggg<-log(n3.codon.pair$catggg/cps.cat.ggg)
cps.obs.exp.catggt<-log(n3.codon.pair$catggt/cps.cat.ggt)
cps.obs.exp.catgta<-log(n3.codon.pair$catgta/cps.cat.gta)
cps.obs.exp.catgtc<-log(n3.codon.pair$catgtc/cps.cat.gtc)
cps.obs.exp.catgtg<-log(n3.codon.pair$catgtg/cps.cat.gtg)
cps.obs.exp.catgtt<-log(n3.codon.pair$catgtt/cps.cat.gtt)

#cps.obs.exp.cattaa<-log(n3.codon.pair$cattaa/cps.cat.taa)
cps.obs.exp.cattac<-log(n3.codon.pair$cattac/cps.cat.tac)
#cps.obs.exp.cattag<-log(n3.codon.pair$cattag/cps.cat.tag)
cps.obs.exp.cattat<-log(n3.codon.pair$cattat/cps.cat.tat)
cps.obs.exp.cattca<-log(n3.codon.pair$cattca/cps.cat.tca)
cps.obs.exp.cattcc<-log(n3.codon.pair$cattcc/cps.cat.tcc)
cps.obs.exp.cattcg<-log(n3.codon.pair$cattcg/cps.cat.tcg)
cps.obs.exp.cattct<-log(n3.codon.pair$cattct/cps.cat.tct)
#cps.obs.exp.cattga<-log(n3.codon.pair$cattga/cps.cat.tga)
cps.obs.exp.cattgc<-log(n3.codon.pair$cattgc/cps.cat.tgc)
cps.obs.exp.cattgg<-log(n3.codon.pair$cattgg/cps.cat.tgg)
cps.obs.exp.cattgt<-log(n3.codon.pair$cattgt/cps.cat.tgt)
cps.obs.exp.cattta<-log(n3.codon.pair$cattta/cps.cat.tta)
cps.obs.exp.catttc<-log(n3.codon.pair$catttc/cps.cat.ttc)
cps.obs.exp.catttg<-log(n3.codon.pair$catttg/cps.cat.ttg)
cps.obs.exp.catttt<-log(n3.codon.pair$catttt/cps.cat.ttt)

















cps.obs.exp.ccaaaa<-log(n3.codon.pair$ccaaaa/cps.cca.aaa)
cps.obs.exp.ccaaac<-log(n3.codon.pair$ccaaac/cps.cca.aac)
cps.obs.exp.ccaaag<-log(n3.codon.pair$ccaaag/cps.cca.aag)
cps.obs.exp.ccaaat<-log(n3.codon.pair$ccaaat/cps.cca.aat)
cps.obs.exp.ccaaca<-log(n3.codon.pair$ccaaca/cps.cca.aca)
cps.obs.exp.ccaacc<-log(n3.codon.pair$ccaacc/cps.cca.acc)
cps.obs.exp.ccaacg<-log(n3.codon.pair$ccaacg/cps.cca.acg)
cps.obs.exp.ccaact<-log(n3.codon.pair$ccaact/cps.cca.act)
cps.obs.exp.ccaaga<-log(n3.codon.pair$ccaaga/cps.cca.aga)
cps.obs.exp.ccaagc<-log(n3.codon.pair$ccaagc/cps.cca.agc)
cps.obs.exp.ccaagg<-log(n3.codon.pair$ccaagg/cps.cca.agg)
cps.obs.exp.ccaagt<-log(n3.codon.pair$ccaagt/cps.cca.agt)
cps.obs.exp.ccaata<-log(n3.codon.pair$ccaata/cps.cca.ata)
cps.obs.exp.ccaatc<-log(n3.codon.pair$ccaatc/cps.cca.atc)
cps.obs.exp.ccaatg<-log(n3.codon.pair$ccaatg/cps.cca.atg)
cps.obs.exp.ccaatt<-log(n3.codon.pair$ccaatt/cps.cca.att)

cps.obs.exp.ccacaa<-log(n3.codon.pair$ccacaa/cps.cca.caa)
cps.obs.exp.ccacac<-log(n3.codon.pair$ccacac/cps.cca.cac)
cps.obs.exp.ccacag<-log(n3.codon.pair$ccacag/cps.cca.cag)
cps.obs.exp.ccacat<-log(n3.codon.pair$ccacat/cps.cca.cat)
cps.obs.exp.ccacca<-log(n3.codon.pair$ccacca/cps.cca.cca)
cps.obs.exp.ccaccc<-log(n3.codon.pair$ccaccc/cps.cca.ccc)
cps.obs.exp.ccaccg<-log(n3.codon.pair$ccaccg/cps.cca.ccg)
cps.obs.exp.ccacct<-log(n3.codon.pair$ccacct/cps.cca.cct)
cps.obs.exp.ccacga<-log(n3.codon.pair$ccacga/cps.cca.cga)
cps.obs.exp.ccacgc<-log(n3.codon.pair$ccacgc/cps.cca.cgc)
cps.obs.exp.ccacgg<-log(n3.codon.pair$ccacgg/cps.cca.cgg)
cps.obs.exp.ccacgt<-log(n3.codon.pair$ccacgt/cps.cca.cgt)
cps.obs.exp.ccacta<-log(n3.codon.pair$ccacta/cps.cca.cta)
cps.obs.exp.ccactc<-log(n3.codon.pair$ccactc/cps.cca.ctc)
cps.obs.exp.ccactg<-log(n3.codon.pair$ccactg/cps.cca.ctg)
cps.obs.exp.ccactt<-log(n3.codon.pair$ccactt/cps.cca.ctt)

cps.obs.exp.ccagaa<-log(n3.codon.pair$ccagaa/cps.cca.gaa)
cps.obs.exp.ccagac<-log(n3.codon.pair$ccagac/cps.cca.gac)
cps.obs.exp.ccagag<-log(n3.codon.pair$ccagag/cps.cca.gag)
cps.obs.exp.ccagat<-log(n3.codon.pair$ccagat/cps.cca.gat)
cps.obs.exp.ccagca<-log(n3.codon.pair$ccagca/cps.cca.gca)
cps.obs.exp.ccagcc<-log(n3.codon.pair$ccagcc/cps.cca.gcc)
cps.obs.exp.ccagcg<-log(n3.codon.pair$ccagcg/cps.cca.gcg)
cps.obs.exp.ccagct<-log(n3.codon.pair$ccagct/cps.cca.gct)
cps.obs.exp.ccagga<-log(n3.codon.pair$ccagga/cps.cca.gga)
cps.obs.exp.ccaggc<-log(n3.codon.pair$ccaggc/cps.cca.ggc)
cps.obs.exp.ccaggg<-log(n3.codon.pair$ccaggg/cps.cca.ggg)
cps.obs.exp.ccaggt<-log(n3.codon.pair$ccaggt/cps.cca.ggt)
cps.obs.exp.ccagta<-log(n3.codon.pair$ccagta/cps.cca.gta)
cps.obs.exp.ccagtc<-log(n3.codon.pair$ccagtc/cps.cca.gtc)
cps.obs.exp.ccagtg<-log(n3.codon.pair$ccagtg/cps.cca.gtg)
cps.obs.exp.ccagtt<-log(n3.codon.pair$ccagtt/cps.cca.gtt)

#cps.obs.exp.ccataa<-log(n3.codon.pair$ccataa/cps.cca.taa)
cps.obs.exp.ccatac<-log(n3.codon.pair$ccatac/cps.cca.tac)
#cps.obs.exp.ccatag<-log(n3.codon.pair$ccatag/cps.cca.tag)
cps.obs.exp.ccatat<-log(n3.codon.pair$ccatat/cps.cca.tat)
cps.obs.exp.ccatca<-log(n3.codon.pair$ccatca/cps.cca.tca)
cps.obs.exp.ccatcc<-log(n3.codon.pair$ccatcc/cps.cca.tcc)
cps.obs.exp.ccatcg<-log(n3.codon.pair$ccatcg/cps.cca.tcg)
cps.obs.exp.ccatct<-log(n3.codon.pair$ccatct/cps.cca.tct)
#cps.obs.exp.ccatga<-log(n3.codon.pair$ccatga/cps.cca.tga)
cps.obs.exp.ccatgc<-log(n3.codon.pair$ccatgc/cps.cca.tgc)
cps.obs.exp.ccatgg<-log(n3.codon.pair$ccatgg/cps.cca.tgg)
cps.obs.exp.ccatgt<-log(n3.codon.pair$ccatgt/cps.cca.tgt)
cps.obs.exp.ccatta<-log(n3.codon.pair$ccatta/cps.cca.tta)
cps.obs.exp.ccattc<-log(n3.codon.pair$ccattc/cps.cca.ttc)
cps.obs.exp.ccattg<-log(n3.codon.pair$ccattg/cps.cca.ttg)
cps.obs.exp.ccattt<-log(n3.codon.pair$ccattt/cps.cca.ttt)







cps.obs.exp.cccaaa<-log(n3.codon.pair$cccaaa/cps.ccc.aaa)
cps.obs.exp.cccaac<-log(n3.codon.pair$cccaac/cps.ccc.aac)
cps.obs.exp.cccaag<-log(n3.codon.pair$cccaag/cps.ccc.aag)
cps.obs.exp.cccaat<-log(n3.codon.pair$cccaat/cps.ccc.aat)
cps.obs.exp.cccaca<-log(n3.codon.pair$cccaca/cps.ccc.aca)
cps.obs.exp.cccacc<-log(n3.codon.pair$cccacc/cps.ccc.acc)
cps.obs.exp.cccacg<-log(n3.codon.pair$cccacg/cps.ccc.acg)
cps.obs.exp.cccact<-log(n3.codon.pair$cccact/cps.ccc.act)
cps.obs.exp.cccaga<-log(n3.codon.pair$cccaga/cps.ccc.aga)
cps.obs.exp.cccagc<-log(n3.codon.pair$cccagc/cps.ccc.agc)
cps.obs.exp.cccagg<-log(n3.codon.pair$cccagg/cps.ccc.agg)
cps.obs.exp.cccagt<-log(n3.codon.pair$cccagt/cps.ccc.agt)
cps.obs.exp.cccata<-log(n3.codon.pair$cccata/cps.ccc.ata)
cps.obs.exp.cccatc<-log(n3.codon.pair$cccatc/cps.ccc.atc)
cps.obs.exp.cccatg<-log(n3.codon.pair$cccatg/cps.ccc.atg)
cps.obs.exp.cccatt<-log(n3.codon.pair$cccatt/cps.ccc.att)

cps.obs.exp.ccccaa<-log(n3.codon.pair$ccccaa/cps.ccc.caa)
cps.obs.exp.ccccac<-log(n3.codon.pair$ccccac/cps.ccc.cac)
cps.obs.exp.ccccag<-log(n3.codon.pair$ccccag/cps.ccc.cag)
cps.obs.exp.ccccat<-log(n3.codon.pair$ccccat/cps.ccc.cat)
cps.obs.exp.ccccca<-log(n3.codon.pair$ccccca/cps.ccc.cca)
cps.obs.exp.cccccc<-log(n3.codon.pair$cccccc/cps.ccc.ccc)
cps.obs.exp.cccccg<-log(n3.codon.pair$cccccg/cps.ccc.ccg)
cps.obs.exp.ccccct<-log(n3.codon.pair$ccccct/cps.ccc.cct)
cps.obs.exp.ccccga<-log(n3.codon.pair$ccccga/cps.ccc.cga)
cps.obs.exp.ccccgc<-log(n3.codon.pair$ccccgc/cps.ccc.cgc)
cps.obs.exp.ccccgg<-log(n3.codon.pair$ccccgg/cps.ccc.cgg)
cps.obs.exp.ccccgt<-log(n3.codon.pair$ccccgt/cps.ccc.cgt)
cps.obs.exp.ccccta<-log(n3.codon.pair$ccccta/cps.ccc.cta)
cps.obs.exp.cccctc<-log(n3.codon.pair$cccctc/cps.ccc.ctc)
cps.obs.exp.cccctg<-log(n3.codon.pair$cccctg/cps.ccc.ctg)
cps.obs.exp.cccctt<-log(n3.codon.pair$cccctt/cps.ccc.ctt)

cps.obs.exp.cccgaa<-log(n3.codon.pair$cccgaa/cps.ccc.gaa)
cps.obs.exp.cccgac<-log(n3.codon.pair$cccgac/cps.ccc.gac)
cps.obs.exp.cccgag<-log(n3.codon.pair$cccgag/cps.ccc.gag)
cps.obs.exp.cccgat<-log(n3.codon.pair$cccgat/cps.ccc.gat)
cps.obs.exp.cccgca<-log(n3.codon.pair$cccgca/cps.ccc.gca)
cps.obs.exp.cccgcc<-log(n3.codon.pair$cccgcc/cps.ccc.gcc)
cps.obs.exp.cccgcg<-log(n3.codon.pair$cccgcg/cps.ccc.gcg)
cps.obs.exp.cccgct<-log(n3.codon.pair$cccgct/cps.ccc.gct)
cps.obs.exp.cccgga<-log(n3.codon.pair$cccgga/cps.ccc.gga)
cps.obs.exp.cccggc<-log(n3.codon.pair$cccggc/cps.ccc.ggc)
cps.obs.exp.cccggg<-log(n3.codon.pair$cccggg/cps.ccc.ggg)
cps.obs.exp.cccggt<-log(n3.codon.pair$cccggt/cps.ccc.ggt)
cps.obs.exp.cccgta<-log(n3.codon.pair$cccgta/cps.ccc.gta)
cps.obs.exp.cccgtc<-log(n3.codon.pair$cccgtc/cps.ccc.gtc)
cps.obs.exp.cccgtg<-log(n3.codon.pair$cccgtg/cps.ccc.gtg)
cps.obs.exp.cccgtt<-log(n3.codon.pair$cccgtt/cps.ccc.gtt)

#cps.obs.exp.ccctaa<-log(n3.codon.pair$ccctaa/cps.ccc.taa)
cps.obs.exp.ccctac<-log(n3.codon.pair$ccctac/cps.ccc.tac)
#cps.obs.exp.ccctag<-log(n3.codon.pair$ccctag/cps.ccc.tag)
cps.obs.exp.ccctat<-log(n3.codon.pair$ccctat/cps.ccc.tat)
cps.obs.exp.ccctca<-log(n3.codon.pair$ccctca/cps.ccc.tca)
cps.obs.exp.ccctcc<-log(n3.codon.pair$ccctcc/cps.ccc.tcc)
cps.obs.exp.ccctcg<-log(n3.codon.pair$ccctcg/cps.ccc.tcg)
cps.obs.exp.ccctct<-log(n3.codon.pair$ccctct/cps.ccc.tct)
#cps.obs.exp.ccctga<-log(n3.codon.pair$ccctga/cps.ccc.tga)
cps.obs.exp.ccctgc<-log(n3.codon.pair$ccctgc/cps.ccc.tgc)
cps.obs.exp.ccctgg<-log(n3.codon.pair$ccctgg/cps.ccc.tgg)
cps.obs.exp.ccctgt<-log(n3.codon.pair$ccctgt/cps.ccc.tgt)
cps.obs.exp.ccctta<-log(n3.codon.pair$ccctta/cps.ccc.tta)
cps.obs.exp.cccttc<-log(n3.codon.pair$cccttc/cps.ccc.ttc)
cps.obs.exp.cccttg<-log(n3.codon.pair$cccttg/cps.ccc.ttg)
cps.obs.exp.cccttt<-log(n3.codon.pair$cccttt/cps.ccc.ttt)









cps.obs.exp.ccgaaa<-log(n3.codon.pair$ccgaaa/cps.ccg.aaa)
cps.obs.exp.ccgaac<-log(n3.codon.pair$ccgaac/cps.ccg.aac)
cps.obs.exp.ccgaag<-log(n3.codon.pair$ccgaag/cps.ccg.aag)
cps.obs.exp.ccgaat<-log(n3.codon.pair$ccgaat/cps.ccg.aat)
cps.obs.exp.ccgaca<-log(n3.codon.pair$ccgaca/cps.ccg.aca)
cps.obs.exp.ccgacc<-log(n3.codon.pair$ccgacc/cps.ccg.acc)
cps.obs.exp.ccgacg<-log(n3.codon.pair$ccgacg/cps.ccg.acg)
cps.obs.exp.ccgact<-log(n3.codon.pair$ccgact/cps.ccg.act)
cps.obs.exp.ccgaga<-log(n3.codon.pair$ccgaga/cps.ccg.aga)
cps.obs.exp.ccgagc<-log(n3.codon.pair$ccgagc/cps.ccg.agc)
cps.obs.exp.ccgagg<-log(n3.codon.pair$ccgagg/cps.ccg.agg)
cps.obs.exp.ccgagt<-log(n3.codon.pair$ccgagt/cps.ccg.agt)
cps.obs.exp.ccgata<-log(n3.codon.pair$ccgata/cps.ccg.ata)
cps.obs.exp.ccgatc<-log(n3.codon.pair$ccgatc/cps.ccg.atc)
cps.obs.exp.ccgatg<-log(n3.codon.pair$ccgatg/cps.ccg.atg)
cps.obs.exp.ccgatt<-log(n3.codon.pair$ccgatt/cps.ccg.att)

cps.obs.exp.ccgcaa<-log(n3.codon.pair$ccgcaa/cps.ccg.caa)
cps.obs.exp.ccgcac<-log(n3.codon.pair$ccgcac/cps.ccg.cac)
cps.obs.exp.ccgcag<-log(n3.codon.pair$ccgcag/cps.ccg.cag)
cps.obs.exp.ccgcat<-log(n3.codon.pair$ccgcat/cps.ccg.cat)
cps.obs.exp.ccgcca<-log(n3.codon.pair$ccgcca/cps.ccg.cca)
cps.obs.exp.ccgccc<-log(n3.codon.pair$ccgccc/cps.ccg.ccc)
cps.obs.exp.ccgccg<-log(n3.codon.pair$ccgccg/cps.ccg.ccg)
cps.obs.exp.ccgcct<-log(n3.codon.pair$ccgcct/cps.ccg.cct)
cps.obs.exp.ccgcga<-log(n3.codon.pair$ccgcga/cps.ccg.cga)
cps.obs.exp.ccgcgc<-log(n3.codon.pair$ccgcgc/cps.ccg.cgc)
cps.obs.exp.ccgcgg<-log(n3.codon.pair$ccgcgg/cps.ccg.cgg)
cps.obs.exp.ccgcgt<-log(n3.codon.pair$ccgcgt/cps.ccg.cgt)
cps.obs.exp.ccgcta<-log(n3.codon.pair$ccgcta/cps.ccg.cta)
cps.obs.exp.ccgctc<-log(n3.codon.pair$ccgctc/cps.ccg.ctc)
cps.obs.exp.ccgctg<-log(n3.codon.pair$ccgctg/cps.ccg.ctg)
cps.obs.exp.ccgctt<-log(n3.codon.pair$ccgctt/cps.ccg.ctt)

cps.obs.exp.ccggaa<-log(n3.codon.pair$ccggaa/cps.ccg.gaa)
cps.obs.exp.ccggac<-log(n3.codon.pair$ccggac/cps.ccg.gac)
cps.obs.exp.ccggag<-log(n3.codon.pair$ccggag/cps.ccg.gag)
cps.obs.exp.ccggat<-log(n3.codon.pair$ccggat/cps.ccg.gat)
cps.obs.exp.ccggca<-log(n3.codon.pair$ccggca/cps.ccg.gca)
cps.obs.exp.ccggcc<-log(n3.codon.pair$ccggcc/cps.ccg.gcc)
cps.obs.exp.ccggcg<-log(n3.codon.pair$ccggcg/cps.ccg.gcg)
cps.obs.exp.ccggct<-log(n3.codon.pair$ccggct/cps.ccg.gct)
cps.obs.exp.ccggga<-log(n3.codon.pair$ccggga/cps.ccg.gga)
cps.obs.exp.ccgggc<-log(n3.codon.pair$ccgggc/cps.ccg.ggc)
cps.obs.exp.ccgggg<-log(n3.codon.pair$ccgggg/cps.ccg.ggg)
cps.obs.exp.ccgggt<-log(n3.codon.pair$ccgggt/cps.ccg.ggt)
cps.obs.exp.ccggta<-log(n3.codon.pair$ccggta/cps.ccg.gta)
cps.obs.exp.ccggtc<-log(n3.codon.pair$ccggtc/cps.ccg.gtc)
cps.obs.exp.ccggtg<-log(n3.codon.pair$ccggtg/cps.ccg.gtg)
cps.obs.exp.ccggtt<-log(n3.codon.pair$ccggtt/cps.ccg.gtt)

#cps.obs.exp.ccgtaa<-log(n3.codon.pair$ccgtaa/cps.ccg.taa)
cps.obs.exp.ccgtac<-log(n3.codon.pair$ccgtac/cps.ccg.tac)
#cps.obs.exp.ccgtag<-log(n3.codon.pair$ccgtag/cps.ccg.tag)
cps.obs.exp.ccgtat<-log(n3.codon.pair$ccgtat/cps.ccg.tat)
cps.obs.exp.ccgtca<-log(n3.codon.pair$ccgtca/cps.ccg.tca)
cps.obs.exp.ccgtcc<-log(n3.codon.pair$ccgtcc/cps.ccg.tcc)
cps.obs.exp.ccgtcg<-log(n3.codon.pair$ccgtcg/cps.ccg.tcg)
cps.obs.exp.ccgtct<-log(n3.codon.pair$ccgtct/cps.ccg.tct)
#cps.obs.exp.ccgtga<-log(n3.codon.pair$ccgtga/cps.ccg.tga)
cps.obs.exp.ccgtgc<-log(n3.codon.pair$ccgtgc/cps.ccg.tgc)
cps.obs.exp.ccgtgg<-log(n3.codon.pair$ccgtgg/cps.ccg.tgg)
cps.obs.exp.ccgtgt<-log(n3.codon.pair$ccgtgt/cps.ccg.tgt)
cps.obs.exp.ccgtta<-log(n3.codon.pair$ccgtta/cps.ccg.tta)
cps.obs.exp.ccgttc<-log(n3.codon.pair$ccgttc/cps.ccg.ttc)
cps.obs.exp.ccgttg<-log(n3.codon.pair$ccgttg/cps.ccg.ttg)
cps.obs.exp.ccgttt<-log(n3.codon.pair$ccgttt/cps.ccg.ttt)







cps.obs.exp.cctaaa<-log(n3.codon.pair$cctaaa/cps.cct.aaa)
cps.obs.exp.cctaac<-log(n3.codon.pair$cctaac/cps.cct.aac)
cps.obs.exp.cctaag<-log(n3.codon.pair$cctaag/cps.cct.aag)
cps.obs.exp.cctaat<-log(n3.codon.pair$cctaat/cps.cct.aat)
cps.obs.exp.cctaca<-log(n3.codon.pair$cctaca/cps.cct.aca)
cps.obs.exp.cctacc<-log(n3.codon.pair$cctacc/cps.cct.acc)
cps.obs.exp.cctacg<-log(n3.codon.pair$cctacg/cps.cct.acg)
cps.obs.exp.cctact<-log(n3.codon.pair$cctact/cps.cct.act)
cps.obs.exp.cctaga<-log(n3.codon.pair$cctaga/cps.cct.aga)
cps.obs.exp.cctagc<-log(n3.codon.pair$cctagc/cps.cct.agc)
cps.obs.exp.cctagg<-log(n3.codon.pair$cctagg/cps.cct.agg)
cps.obs.exp.cctagt<-log(n3.codon.pair$cctagt/cps.cct.agt)
cps.obs.exp.cctata<-log(n3.codon.pair$cctata/cps.cct.ata)
cps.obs.exp.cctatc<-log(n3.codon.pair$cctatc/cps.cct.atc)
cps.obs.exp.cctatg<-log(n3.codon.pair$cctatg/cps.cct.atg)
cps.obs.exp.cctatt<-log(n3.codon.pair$cctatt/cps.cct.att)

cps.obs.exp.cctcaa<-log(n3.codon.pair$cctcaa/cps.cct.caa)
cps.obs.exp.cctcac<-log(n3.codon.pair$cctcac/cps.cct.cac)
cps.obs.exp.cctcag<-log(n3.codon.pair$cctcag/cps.cct.cag)
cps.obs.exp.cctcat<-log(n3.codon.pair$cctcat/cps.cct.cat)
cps.obs.exp.cctcca<-log(n3.codon.pair$cctcca/cps.cct.cca)
cps.obs.exp.cctccc<-log(n3.codon.pair$cctccc/cps.cct.ccc)
cps.obs.exp.cctccg<-log(n3.codon.pair$cctccg/cps.cct.ccg)
cps.obs.exp.cctcct<-log(n3.codon.pair$cctcct/cps.cct.cct)
cps.obs.exp.cctcga<-log(n3.codon.pair$cctcga/cps.cct.cga)
cps.obs.exp.cctcgc<-log(n3.codon.pair$cctcgc/cps.cct.cgc)
cps.obs.exp.cctcgg<-log(n3.codon.pair$cctcgg/cps.cct.cgg)
cps.obs.exp.cctcgt<-log(n3.codon.pair$cctcgt/cps.cct.cgt)
cps.obs.exp.cctcta<-log(n3.codon.pair$cctcta/cps.cct.cta)
cps.obs.exp.cctctc<-log(n3.codon.pair$cctctc/cps.cct.ctc)
cps.obs.exp.cctctg<-log(n3.codon.pair$cctctg/cps.cct.ctg)
cps.obs.exp.cctctt<-log(n3.codon.pair$cctctt/cps.cct.ctt)

cps.obs.exp.cctgaa<-log(n3.codon.pair$cctgaa/cps.cct.gaa)
cps.obs.exp.cctgac<-log(n3.codon.pair$cctgac/cps.cct.gac)
cps.obs.exp.cctgag<-log(n3.codon.pair$cctgag/cps.cct.gag)
cps.obs.exp.cctgat<-log(n3.codon.pair$cctgat/cps.cct.gat)
cps.obs.exp.cctgca<-log(n3.codon.pair$cctgca/cps.cct.gca)
cps.obs.exp.cctgcc<-log(n3.codon.pair$cctgcc/cps.cct.gcc)
cps.obs.exp.cctgcg<-log(n3.codon.pair$cctgcg/cps.cct.gcg)
cps.obs.exp.cctgct<-log(n3.codon.pair$cctgct/cps.cct.gct)
cps.obs.exp.cctgga<-log(n3.codon.pair$cctgga/cps.cct.gga)
cps.obs.exp.cctggc<-log(n3.codon.pair$cctggc/cps.cct.ggc)
cps.obs.exp.cctggg<-log(n3.codon.pair$cctggg/cps.cct.ggg)
cps.obs.exp.cctggt<-log(n3.codon.pair$cctggt/cps.cct.ggt)
cps.obs.exp.cctgta<-log(n3.codon.pair$cctgta/cps.cct.gta)
cps.obs.exp.cctgtc<-log(n3.codon.pair$cctgtc/cps.cct.gtc)
cps.obs.exp.cctgtg<-log(n3.codon.pair$cctgtg/cps.cct.gtg)
cps.obs.exp.cctgtt<-log(n3.codon.pair$cctgtt/cps.cct.gtt)

#cps.obs.exp.ccttaa<-log(n3.codon.pair$ccttaa/cps.cct.taa)
cps.obs.exp.ccttac<-log(n3.codon.pair$ccttac/cps.cct.tac)
#cps.obs.exp.ccttag<-log(n3.codon.pair$ccttag/cps.cct.tag)
cps.obs.exp.ccttat<-log(n3.codon.pair$ccttat/cps.cct.tat)
cps.obs.exp.ccttca<-log(n3.codon.pair$ccttca/cps.cct.tca)
cps.obs.exp.ccttcc<-log(n3.codon.pair$ccttcc/cps.cct.tcc)
cps.obs.exp.ccttcg<-log(n3.codon.pair$ccttcg/cps.cct.tcg)
cps.obs.exp.ccttct<-log(n3.codon.pair$ccttct/cps.cct.tct)
#cps.obs.exp.ccttga<-log(n3.codon.pair$ccttga/cps.cct.tga)
cps.obs.exp.ccttgc<-log(n3.codon.pair$ccttgc/cps.cct.tgc)
cps.obs.exp.ccttgg<-log(n3.codon.pair$ccttgg/cps.cct.tgg)
cps.obs.exp.ccttgt<-log(n3.codon.pair$ccttgt/cps.cct.tgt)
cps.obs.exp.ccttta<-log(n3.codon.pair$ccttta/cps.cct.tta)
cps.obs.exp.cctttc<-log(n3.codon.pair$cctttc/cps.cct.ttc)
cps.obs.exp.cctttg<-log(n3.codon.pair$cctttg/cps.cct.ttg)
cps.obs.exp.cctttt<-log(n3.codon.pair$cctttt/cps.cct.ttt)



















cps.obs.exp.cgaaaa<-log(n3.codon.pair$cgaaaa/cps.cga.aaa)
cps.obs.exp.cgaaac<-log(n3.codon.pair$cgaaac/cps.cga.aac)
cps.obs.exp.cgaaag<-log(n3.codon.pair$cgaaag/cps.cga.aag)
cps.obs.exp.cgaaat<-log(n3.codon.pair$cgaaat/cps.cga.aat)
cps.obs.exp.cgaaca<-log(n3.codon.pair$cgaaca/cps.cga.aca)
cps.obs.exp.cgaacc<-log(n3.codon.pair$cgaacc/cps.cga.acc)
cps.obs.exp.cgaacg<-log(n3.codon.pair$cgaacg/cps.cga.acg)
cps.obs.exp.cgaact<-log(n3.codon.pair$cgaact/cps.cga.act)
cps.obs.exp.cgaaga<-log(n3.codon.pair$cgaaga/cps.cga.aga)
cps.obs.exp.cgaagc<-log(n3.codon.pair$cgaagc/cps.cga.agc)
cps.obs.exp.cgaagg<-log(n3.codon.pair$cgaagg/cps.cga.agg)
cps.obs.exp.cgaagt<-log(n3.codon.pair$cgaagt/cps.cga.agt)
cps.obs.exp.cgaata<-log(n3.codon.pair$cgaata/cps.cga.ata)
cps.obs.exp.cgaatc<-log(n3.codon.pair$cgaatc/cps.cga.atc)
cps.obs.exp.cgaatg<-log(n3.codon.pair$cgaatg/cps.cga.atg)
cps.obs.exp.cgaatt<-log(n3.codon.pair$cgaatt/cps.cga.att)

cps.obs.exp.cgacaa<-log(n3.codon.pair$cgacaa/cps.cga.caa)
cps.obs.exp.cgacac<-log(n3.codon.pair$cgacac/cps.cga.cac)
cps.obs.exp.cgacag<-log(n3.codon.pair$cgacag/cps.cga.cag)
cps.obs.exp.cgacat<-log(n3.codon.pair$cgacat/cps.cga.cat)
cps.obs.exp.cgacca<-log(n3.codon.pair$cgacca/cps.cga.cca)
cps.obs.exp.cgaccc<-log(n3.codon.pair$cgaccc/cps.cga.ccc)
cps.obs.exp.cgaccg<-log(n3.codon.pair$cgaccg/cps.cga.ccg)
cps.obs.exp.cgacct<-log(n3.codon.pair$cgacct/cps.cga.cct)
cps.obs.exp.cgacga<-log(n3.codon.pair$cgacga/cps.cga.cga)
cps.obs.exp.cgacgc<-log(n3.codon.pair$cgacgc/cps.cga.cgc)
cps.obs.exp.cgacgg<-log(n3.codon.pair$cgacgg/cps.cga.cgg)
cps.obs.exp.cgacgt<-log(n3.codon.pair$cgacgt/cps.cga.cgt)
cps.obs.exp.cgacta<-log(n3.codon.pair$cgacta/cps.cga.cta)
cps.obs.exp.cgactc<-log(n3.codon.pair$cgactc/cps.cga.ctc)
cps.obs.exp.cgactg<-log(n3.codon.pair$cgactg/cps.cga.ctg)
cps.obs.exp.cgactt<-log(n3.codon.pair$cgactt/cps.cga.ctt)

cps.obs.exp.cgagaa<-log(n3.codon.pair$cgagaa/cps.cga.gaa)
cps.obs.exp.cgagac<-log(n3.codon.pair$cgagac/cps.cga.gac)
cps.obs.exp.cgagag<-log(n3.codon.pair$cgagag/cps.cga.gag)
cps.obs.exp.cgagat<-log(n3.codon.pair$cgagat/cps.cga.gat)
cps.obs.exp.cgagca<-log(n3.codon.pair$cgagca/cps.cga.gca)
cps.obs.exp.cgagcc<-log(n3.codon.pair$cgagcc/cps.cga.gcc)
cps.obs.exp.cgagcg<-log(n3.codon.pair$cgagcg/cps.cga.gcg)
cps.obs.exp.cgagct<-log(n3.codon.pair$cgagct/cps.cga.gct)
cps.obs.exp.cgagga<-log(n3.codon.pair$cgagga/cps.cga.gga)
cps.obs.exp.cgaggc<-log(n3.codon.pair$cgaggc/cps.cga.ggc)
cps.obs.exp.cgaggg<-log(n3.codon.pair$cgaggg/cps.cga.ggg)
cps.obs.exp.cgaggt<-log(n3.codon.pair$cgaggt/cps.cga.ggt)
cps.obs.exp.cgagta<-log(n3.codon.pair$cgagta/cps.cga.gta)
cps.obs.exp.cgagtc<-log(n3.codon.pair$cgagtc/cps.cga.gtc)
cps.obs.exp.cgagtg<-log(n3.codon.pair$cgagtg/cps.cga.gtg)
cps.obs.exp.cgagtt<-log(n3.codon.pair$cgagtt/cps.cga.gtt)

#cps.obs.exp.cgataa<-log(n3.codon.pair$cgataa/cps.cga.taa)
cps.obs.exp.cgatac<-log(n3.codon.pair$cgatac/cps.cga.tac)
#cps.obs.exp.cgatag<-log(n3.codon.pair$cgatag/cps.cga.tag)
cps.obs.exp.cgatat<-log(n3.codon.pair$cgatat/cps.cga.tat)
cps.obs.exp.cgatca<-log(n3.codon.pair$cgatca/cps.cga.tca)
cps.obs.exp.cgatcc<-log(n3.codon.pair$cgatcc/cps.cga.tcc)
cps.obs.exp.cgatcg<-log(n3.codon.pair$cgatcg/cps.cga.tcg)
cps.obs.exp.cgatct<-log(n3.codon.pair$cgatct/cps.cga.tct)
#cps.obs.exp.cgatga<-log(n3.codon.pair$cgatga/cps.cga.tga)
cps.obs.exp.cgatgc<-log(n3.codon.pair$cgatgc/cps.cga.tgc)
cps.obs.exp.cgatgg<-log(n3.codon.pair$cgatgg/cps.cga.tgg)
cps.obs.exp.cgatgt<-log(n3.codon.pair$cgatgt/cps.cga.tgt)
cps.obs.exp.cgatta<-log(n3.codon.pair$cgatta/cps.cga.tta)
cps.obs.exp.cgattc<-log(n3.codon.pair$cgattc/cps.cga.ttc)
cps.obs.exp.cgattg<-log(n3.codon.pair$cgattg/cps.cga.ttg)
cps.obs.exp.cgattt<-log(n3.codon.pair$cgattt/cps.cga.ttt)








cps.obs.exp.cgcaaa<-log(n3.codon.pair$cgcaaa/cps.cgc.aaa)
cps.obs.exp.cgcaac<-log(n3.codon.pair$cgcaac/cps.cgc.aac)
cps.obs.exp.cgcaag<-log(n3.codon.pair$cgcaag/cps.cgc.aag)
cps.obs.exp.cgcaat<-log(n3.codon.pair$cgcaat/cps.cgc.aat)
cps.obs.exp.cgcaca<-log(n3.codon.pair$cgcaca/cps.cgc.aca)
cps.obs.exp.cgcacc<-log(n3.codon.pair$cgcacc/cps.cgc.acc)
cps.obs.exp.cgcacg<-log(n3.codon.pair$cgcacg/cps.cgc.acg)
cps.obs.exp.cgcact<-log(n3.codon.pair$cgcact/cps.cgc.act)
cps.obs.exp.cgcaga<-log(n3.codon.pair$cgcaga/cps.cgc.aga)
cps.obs.exp.cgcagc<-log(n3.codon.pair$cgcagc/cps.cgc.agc)
cps.obs.exp.cgcagg<-log(n3.codon.pair$cgcagg/cps.cgc.agg)
cps.obs.exp.cgcagt<-log(n3.codon.pair$cgcagt/cps.cgc.agt)
cps.obs.exp.cgcata<-log(n3.codon.pair$cgcata/cps.cgc.ata)
cps.obs.exp.cgcatc<-log(n3.codon.pair$cgcatc/cps.cgc.atc)
cps.obs.exp.cgcatg<-log(n3.codon.pair$cgcatg/cps.cgc.atg)
cps.obs.exp.cgcatt<-log(n3.codon.pair$cgcatt/cps.cgc.att)

cps.obs.exp.cgccaa<-log(n3.codon.pair$cgccaa/cps.cgc.caa)
cps.obs.exp.cgccac<-log(n3.codon.pair$cgccac/cps.cgc.cac)
cps.obs.exp.cgccag<-log(n3.codon.pair$cgccag/cps.cgc.cag)
cps.obs.exp.cgccat<-log(n3.codon.pair$cgccat/cps.cgc.cat)
cps.obs.exp.cgccca<-log(n3.codon.pair$cgccca/cps.cgc.cca)
cps.obs.exp.cgcccc<-log(n3.codon.pair$cgcccc/cps.cgc.ccc)
cps.obs.exp.cgcccg<-log(n3.codon.pair$cgcccg/cps.cgc.ccg)
cps.obs.exp.cgccct<-log(n3.codon.pair$cgccct/cps.cgc.cct)
cps.obs.exp.cgccga<-log(n3.codon.pair$cgccga/cps.cgc.cga)
cps.obs.exp.cgccgc<-log(n3.codon.pair$cgccgc/cps.cgc.cgc)
cps.obs.exp.cgccgg<-log(n3.codon.pair$cgccgg/cps.cgc.cgg)
cps.obs.exp.cgccgt<-log(n3.codon.pair$cgccgt/cps.cgc.cgt)
cps.obs.exp.cgccta<-log(n3.codon.pair$cgccta/cps.cgc.cta)
cps.obs.exp.cgcctc<-log(n3.codon.pair$cgcctc/cps.cgc.ctc)
cps.obs.exp.cgcctg<-log(n3.codon.pair$cgcctg/cps.cgc.ctg)
cps.obs.exp.cgcctt<-log(n3.codon.pair$cgcctt/cps.cgc.ctt)

cps.obs.exp.cgcgaa<-log(n3.codon.pair$cgcgaa/cps.cgc.gaa)
cps.obs.exp.cgcgac<-log(n3.codon.pair$cgcgac/cps.cgc.gac)
cps.obs.exp.cgcgag<-log(n3.codon.pair$cgcgag/cps.cgc.gag)
cps.obs.exp.cgcgat<-log(n3.codon.pair$cgcgat/cps.cgc.gat)
cps.obs.exp.cgcgca<-log(n3.codon.pair$cgcgca/cps.cgc.gca)
cps.obs.exp.cgcgcc<-log(n3.codon.pair$cgcgcc/cps.cgc.gcc)
cps.obs.exp.cgcgcg<-log(n3.codon.pair$cgcgcg/cps.cgc.gcg)
cps.obs.exp.cgcgct<-log(n3.codon.pair$cgcgct/cps.cgc.gct)
cps.obs.exp.cgcgga<-log(n3.codon.pair$cgcgga/cps.cgc.gga)
cps.obs.exp.cgcggc<-log(n3.codon.pair$cgcggc/cps.cgc.ggc)
cps.obs.exp.cgcggg<-log(n3.codon.pair$cgcggg/cps.cgc.ggg)
cps.obs.exp.cgcggt<-log(n3.codon.pair$cgcggt/cps.cgc.ggt)
cps.obs.exp.cgcgta<-log(n3.codon.pair$cgcgta/cps.cgc.gta)
cps.obs.exp.cgcgtc<-log(n3.codon.pair$cgcgtc/cps.cgc.gtc)
cps.obs.exp.cgcgtg<-log(n3.codon.pair$cgcgtg/cps.cgc.gtg)
cps.obs.exp.cgcgtt<-log(n3.codon.pair$cgcgtt/cps.cgc.gtt)

#cps.obs.exp.cgctaa<-log(n3.codon.pair$cgctaa/cps.cgc.taa)
cps.obs.exp.cgctac<-log(n3.codon.pair$cgctac/cps.cgc.tac)
#cps.obs.exp.cgctag<-log(n3.codon.pair$cgctag/cps.cgc.tag)
cps.obs.exp.cgctat<-log(n3.codon.pair$cgctat/cps.cgc.tat)
cps.obs.exp.cgctca<-log(n3.codon.pair$cgctca/cps.cgc.tca)
cps.obs.exp.cgctcc<-log(n3.codon.pair$cgctcc/cps.cgc.tcc)
cps.obs.exp.cgctcg<-log(n3.codon.pair$cgctcg/cps.cgc.tcg)
cps.obs.exp.cgctct<-log(n3.codon.pair$cgctct/cps.cgc.tct)
#cps.obs.exp.cgctga<-log(n3.codon.pair$cgctga/cps.cgc.tga)
cps.obs.exp.cgctgc<-log(n3.codon.pair$cgctgc/cps.cgc.tgc)
cps.obs.exp.cgctgg<-log(n3.codon.pair$cgctgg/cps.cgc.tgg)
cps.obs.exp.cgctgt<-log(n3.codon.pair$cgctgt/cps.cgc.tgt)
cps.obs.exp.cgctta<-log(n3.codon.pair$cgctta/cps.cgc.tta)
cps.obs.exp.cgcttc<-log(n3.codon.pair$cgcttc/cps.cgc.ttc)
cps.obs.exp.cgcttg<-log(n3.codon.pair$cgcttg/cps.cgc.ttg)
cps.obs.exp.cgcttt<-log(n3.codon.pair$cgcttt/cps.cgc.ttt)










cps.obs.exp.cggaaa<-log(n3.codon.pair$cggaaa/cps.cgg.aaa)
cps.obs.exp.cggaac<-log(n3.codon.pair$cggaac/cps.cgg.aac)
cps.obs.exp.cggaag<-log(n3.codon.pair$cggaag/cps.cgg.aag)
cps.obs.exp.cggaat<-log(n3.codon.pair$cggaat/cps.cgg.aat)
cps.obs.exp.cggaca<-log(n3.codon.pair$cggaca/cps.cgg.aca)
cps.obs.exp.cggacc<-log(n3.codon.pair$cggacc/cps.cgg.acc)
cps.obs.exp.cggacg<-log(n3.codon.pair$cggacg/cps.cgg.acg)
cps.obs.exp.cggact<-log(n3.codon.pair$cggact/cps.cgg.act)
cps.obs.exp.cggaga<-log(n3.codon.pair$cggaga/cps.cgg.aga)
cps.obs.exp.cggagc<-log(n3.codon.pair$cggagc/cps.cgg.agc)
cps.obs.exp.cggagg<-log(n3.codon.pair$cggagg/cps.cgg.agg)
cps.obs.exp.cggagt<-log(n3.codon.pair$cggagt/cps.cgg.agt)
cps.obs.exp.cggata<-log(n3.codon.pair$cggata/cps.cgg.ata)
cps.obs.exp.cggatc<-log(n3.codon.pair$cggatc/cps.cgg.atc)
cps.obs.exp.cggatg<-log(n3.codon.pair$cggatg/cps.cgg.atg)
cps.obs.exp.cggatt<-log(n3.codon.pair$cggatt/cps.cgg.att)

cps.obs.exp.cggcaa<-log(n3.codon.pair$cggcaa/cps.cgg.caa)
cps.obs.exp.cggcac<-log(n3.codon.pair$cggcac/cps.cgg.cac)
cps.obs.exp.cggcag<-log(n3.codon.pair$cggcag/cps.cgg.cag)
cps.obs.exp.cggcat<-log(n3.codon.pair$cggcat/cps.cgg.cat)
cps.obs.exp.cggcca<-log(n3.codon.pair$cggcca/cps.cgg.cca)
cps.obs.exp.cggccc<-log(n3.codon.pair$cggccc/cps.cgg.ccc)
cps.obs.exp.cggccg<-log(n3.codon.pair$cggccg/cps.cgg.ccg)
cps.obs.exp.cggcct<-log(n3.codon.pair$cggcct/cps.cgg.cct)
cps.obs.exp.cggcga<-log(n3.codon.pair$cggcga/cps.cgg.cga)
cps.obs.exp.cggcgc<-log(n3.codon.pair$cggcgc/cps.cgg.cgc)
cps.obs.exp.cggcgg<-log(n3.codon.pair$cggcgg/cps.cgg.cgg)
cps.obs.exp.cggcgt<-log(n3.codon.pair$cggcgt/cps.cgg.cgt)
cps.obs.exp.cggcta<-log(n3.codon.pair$cggcta/cps.cgg.cta)
cps.obs.exp.cggctc<-log(n3.codon.pair$cggctc/cps.cgg.ctc)
cps.obs.exp.cggctg<-log(n3.codon.pair$cggctg/cps.cgg.ctg)
cps.obs.exp.cggctt<-log(n3.codon.pair$cggctt/cps.cgg.ctt)

cps.obs.exp.cgggaa<-log(n3.codon.pair$cgggaa/cps.cgg.gaa)
cps.obs.exp.cgggac<-log(n3.codon.pair$cgggac/cps.cgg.gac)
cps.obs.exp.cgggag<-log(n3.codon.pair$cgggag/cps.cgg.gag)
cps.obs.exp.cgggat<-log(n3.codon.pair$cgggat/cps.cgg.gat)
cps.obs.exp.cgggca<-log(n3.codon.pair$cgggca/cps.cgg.gca)
cps.obs.exp.cgggcc<-log(n3.codon.pair$cgggcc/cps.cgg.gcc)
cps.obs.exp.cgggcg<-log(n3.codon.pair$cgggcg/cps.cgg.gcg)
cps.obs.exp.cgggct<-log(n3.codon.pair$cgggct/cps.cgg.gct)
cps.obs.exp.cgggga<-log(n3.codon.pair$cgggga/cps.cgg.gga)
cps.obs.exp.cggggc<-log(n3.codon.pair$cggggc/cps.cgg.ggc)
cps.obs.exp.cggggg<-log(n3.codon.pair$cggggg/cps.cgg.ggg)
cps.obs.exp.cggggt<-log(n3.codon.pair$cggggt/cps.cgg.ggt)
cps.obs.exp.cgggta<-log(n3.codon.pair$cgggta/cps.cgg.gta)
cps.obs.exp.cgggtc<-log(n3.codon.pair$cgggtc/cps.cgg.gtc)
cps.obs.exp.cgggtg<-log(n3.codon.pair$cgggtg/cps.cgg.gtg)
cps.obs.exp.cgggtt<-log(n3.codon.pair$cgggtt/cps.cgg.gtt)

#cps.obs.exp.cggtaa<-log(n3.codon.pair$cggtaa/cps.cgg.taa)
cps.obs.exp.cggtac<-log(n3.codon.pair$cggtac/cps.cgg.tac)
#cps.obs.exp.cggtag<-log(n3.codon.pair$cggtag/cps.cgg.tag)
cps.obs.exp.cggtat<-log(n3.codon.pair$cggtat/cps.cgg.tat)
cps.obs.exp.cggtca<-log(n3.codon.pair$cggtca/cps.cgg.tca)
cps.obs.exp.cggtcc<-log(n3.codon.pair$cggtcc/cps.cgg.tcc)
cps.obs.exp.cggtcg<-log(n3.codon.pair$cggtcg/cps.cgg.tcg)
cps.obs.exp.cggtct<-log(n3.codon.pair$cggtct/cps.cgg.tct)
#cps.obs.exp.cggtga<-log(n3.codon.pair$cggtga/cps.cgg.tga)
cps.obs.exp.cggtgc<-log(n3.codon.pair$cggtgc/cps.cgg.tgc)
cps.obs.exp.cggtgg<-log(n3.codon.pair$cggtgg/cps.cgg.tgg)
cps.obs.exp.cggtgt<-log(n3.codon.pair$cggtgt/cps.cgg.tgt)
cps.obs.exp.cggtta<-log(n3.codon.pair$cggtta/cps.cgg.tta)
cps.obs.exp.cggttc<-log(n3.codon.pair$cggttc/cps.cgg.ttc)
cps.obs.exp.cggttg<-log(n3.codon.pair$cggttg/cps.cgg.ttg)
cps.obs.exp.cggttt<-log(n3.codon.pair$cggttt/cps.cgg.ttt)








cps.obs.exp.cgtaaa<-log(n3.codon.pair$cgtaaa/cps.cgt.aaa)
cps.obs.exp.cgtaac<-log(n3.codon.pair$cgtaac/cps.cgt.aac)
cps.obs.exp.cgtaag<-log(n3.codon.pair$cgtaag/cps.cgt.aag)
cps.obs.exp.cgtaat<-log(n3.codon.pair$cgtaat/cps.cgt.aat)
cps.obs.exp.cgtaca<-log(n3.codon.pair$cgtaca/cps.cgt.aca)
cps.obs.exp.cgtacc<-log(n3.codon.pair$cgtacc/cps.cgt.acc)
cps.obs.exp.cgtacg<-log(n3.codon.pair$cgtacg/cps.cgt.acg)
cps.obs.exp.cgtact<-log(n3.codon.pair$cgtact/cps.cgt.act)
cps.obs.exp.cgtaga<-log(n3.codon.pair$cgtaga/cps.cgt.aga)
cps.obs.exp.cgtagc<-log(n3.codon.pair$cgtagc/cps.cgt.agc)
cps.obs.exp.cgtagg<-log(n3.codon.pair$cgtagg/cps.cgt.agg)
cps.obs.exp.cgtagt<-log(n3.codon.pair$cgtagt/cps.cgt.agt)
cps.obs.exp.cgtata<-log(n3.codon.pair$cgtata/cps.cgt.ata)
cps.obs.exp.cgtatc<-log(n3.codon.pair$cgtatc/cps.cgt.atc)
cps.obs.exp.cgtatg<-log(n3.codon.pair$cgtatg/cps.cgt.atg)
cps.obs.exp.cgtatt<-log(n3.codon.pair$cgtatt/cps.cgt.att)

cps.obs.exp.cgtcaa<-log(n3.codon.pair$cgtcaa/cps.cgt.caa)
cps.obs.exp.cgtcac<-log(n3.codon.pair$cgtcac/cps.cgt.cac)
cps.obs.exp.cgtcag<-log(n3.codon.pair$cgtcag/cps.cgt.cag)
cps.obs.exp.cgtcat<-log(n3.codon.pair$cgtcat/cps.cgt.cat)
cps.obs.exp.cgtcca<-log(n3.codon.pair$cgtcca/cps.cgt.cca)
cps.obs.exp.cgtccc<-log(n3.codon.pair$cgtccc/cps.cgt.ccc)
cps.obs.exp.cgtccg<-log(n3.codon.pair$cgtccg/cps.cgt.ccg)
cps.obs.exp.cgtcct<-log(n3.codon.pair$cgtcct/cps.cgt.cct)
cps.obs.exp.cgtcga<-log(n3.codon.pair$cgtcga/cps.cgt.cga)
cps.obs.exp.cgtcgc<-log(n3.codon.pair$cgtcgc/cps.cgt.cgc)
cps.obs.exp.cgtcgg<-log(n3.codon.pair$cgtcgg/cps.cgt.cgg)
cps.obs.exp.cgtcgt<-log(n3.codon.pair$cgtcgt/cps.cgt.cgt)
cps.obs.exp.cgtcta<-log(n3.codon.pair$cgtcta/cps.cgt.cta)
cps.obs.exp.cgtctc<-log(n3.codon.pair$cgtctc/cps.cgt.ctc)
cps.obs.exp.cgtctg<-log(n3.codon.pair$cgtctg/cps.cgt.ctg)
cps.obs.exp.cgtctt<-log(n3.codon.pair$cgtctt/cps.cgt.ctt)

cps.obs.exp.cgtgaa<-log(n3.codon.pair$cgtgaa/cps.cgt.gaa)
cps.obs.exp.cgtgac<-log(n3.codon.pair$cgtgac/cps.cgt.gac)
cps.obs.exp.cgtgag<-log(n3.codon.pair$cgtgag/cps.cgt.gag)
cps.obs.exp.cgtgat<-log(n3.codon.pair$cgtgat/cps.cgt.gat)
cps.obs.exp.cgtgca<-log(n3.codon.pair$cgtgca/cps.cgt.gca)
cps.obs.exp.cgtgcc<-log(n3.codon.pair$cgtgcc/cps.cgt.gcc)
cps.obs.exp.cgtgcg<-log(n3.codon.pair$cgtgcg/cps.cgt.gcg)
cps.obs.exp.cgtgct<-log(n3.codon.pair$cgtgct/cps.cgt.gct)
cps.obs.exp.cgtgga<-log(n3.codon.pair$cgtgga/cps.cgt.gga)
cps.obs.exp.cgtggc<-log(n3.codon.pair$cgtggc/cps.cgt.ggc)
cps.obs.exp.cgtggg<-log(n3.codon.pair$cgtggg/cps.cgt.ggg)
cps.obs.exp.cgtggt<-log(n3.codon.pair$cgtggt/cps.cgt.ggt)
cps.obs.exp.cgtgta<-log(n3.codon.pair$cgtgta/cps.cgt.gta)
cps.obs.exp.cgtgtc<-log(n3.codon.pair$cgtgtc/cps.cgt.gtc)
cps.obs.exp.cgtgtg<-log(n3.codon.pair$cgtgtg/cps.cgt.gtg)
cps.obs.exp.cgtgtt<-log(n3.codon.pair$cgtgtt/cps.cgt.gtt)

#cps.obs.exp.cgttaa<-log(n3.codon.pair$cgttaa/cps.cgt.taa)
cps.obs.exp.cgttac<-log(n3.codon.pair$cgttac/cps.cgt.tac)
#cps.obs.exp.cgttag<-log(n3.codon.pair$cgttag/cps.cgt.tag)
cps.obs.exp.cgttat<-log(n3.codon.pair$cgttat/cps.cgt.tat)
cps.obs.exp.cgttca<-log(n3.codon.pair$cgttca/cps.cgt.tca)
cps.obs.exp.cgttcc<-log(n3.codon.pair$cgttcc/cps.cgt.tcc)
cps.obs.exp.cgttcg<-log(n3.codon.pair$cgttcg/cps.cgt.tcg)
cps.obs.exp.cgttct<-log(n3.codon.pair$cgttct/cps.cgt.tct)
#cps.obs.exp.cgttga<-log(n3.codon.pair$cgttga/cps.cgt.tga)
cps.obs.exp.cgttgc<-log(n3.codon.pair$cgttgc/cps.cgt.tgc)
cps.obs.exp.cgttgg<-log(n3.codon.pair$cgttgg/cps.cgt.tgg)
cps.obs.exp.cgttgt<-log(n3.codon.pair$cgttgt/cps.cgt.tgt)
cps.obs.exp.cgttta<-log(n3.codon.pair$cgttta/cps.cgt.tta)
cps.obs.exp.cgtttc<-log(n3.codon.pair$cgtttc/cps.cgt.ttc)
cps.obs.exp.cgtttg<-log(n3.codon.pair$cgtttg/cps.cgt.ttg)
cps.obs.exp.cgtttt<-log(n3.codon.pair$cgtttt/cps.cgt.ttt)




















cps.obs.exp.ctaaaa<-log(n3.codon.pair$ctaaaa/cps.cta.aaa)
cps.obs.exp.ctaaac<-log(n3.codon.pair$ctaaac/cps.cta.aac)
cps.obs.exp.ctaaag<-log(n3.codon.pair$ctaaag/cps.cta.aag)
cps.obs.exp.ctaaat<-log(n3.codon.pair$ctaaat/cps.cta.aat)
cps.obs.exp.ctaaca<-log(n3.codon.pair$ctaaca/cps.cta.aca)
cps.obs.exp.ctaacc<-log(n3.codon.pair$ctaacc/cps.cta.acc)
cps.obs.exp.ctaacg<-log(n3.codon.pair$ctaacg/cps.cta.acg)
cps.obs.exp.ctaact<-log(n3.codon.pair$ctaact/cps.cta.act)
cps.obs.exp.ctaaga<-log(n3.codon.pair$ctaaga/cps.cta.aga)
cps.obs.exp.ctaagc<-log(n3.codon.pair$ctaagc/cps.cta.agc)
cps.obs.exp.ctaagg<-log(n3.codon.pair$ctaagg/cps.cta.agg)
cps.obs.exp.ctaagt<-log(n3.codon.pair$ctaagt/cps.cta.agt)
cps.obs.exp.ctaata<-log(n3.codon.pair$ctaata/cps.cta.ata)
cps.obs.exp.ctaatc<-log(n3.codon.pair$ctaatc/cps.cta.atc)
cps.obs.exp.ctaatg<-log(n3.codon.pair$ctaatg/cps.cta.atg)
cps.obs.exp.ctaatt<-log(n3.codon.pair$ctaatt/cps.cta.att)

cps.obs.exp.ctacaa<-log(n3.codon.pair$ctacaa/cps.cta.caa)
cps.obs.exp.ctacac<-log(n3.codon.pair$ctacac/cps.cta.cac)
cps.obs.exp.ctacag<-log(n3.codon.pair$ctacag/cps.cta.cag)
cps.obs.exp.ctacat<-log(n3.codon.pair$ctacat/cps.cta.cat)
cps.obs.exp.ctacca<-log(n3.codon.pair$ctacca/cps.cta.cca)
cps.obs.exp.ctaccc<-log(n3.codon.pair$ctaccc/cps.cta.ccc)
cps.obs.exp.ctaccg<-log(n3.codon.pair$ctaccg/cps.cta.ccg)
cps.obs.exp.ctacct<-log(n3.codon.pair$ctacct/cps.cta.cct)
cps.obs.exp.ctacga<-log(n3.codon.pair$ctacga/cps.cta.cga)
cps.obs.exp.ctacgc<-log(n3.codon.pair$ctacgc/cps.cta.cgc)
cps.obs.exp.ctacgg<-log(n3.codon.pair$ctacgg/cps.cta.cgg)
cps.obs.exp.ctacgt<-log(n3.codon.pair$ctacgt/cps.cta.cgt)
cps.obs.exp.ctacta<-log(n3.codon.pair$ctacta/cps.cta.cta)
cps.obs.exp.ctactc<-log(n3.codon.pair$ctactc/cps.cta.ctc)
cps.obs.exp.ctactg<-log(n3.codon.pair$ctactg/cps.cta.ctg)
cps.obs.exp.ctactt<-log(n3.codon.pair$ctactt/cps.cta.ctt)

cps.obs.exp.ctagaa<-log(n3.codon.pair$ctagaa/cps.cta.gaa)
cps.obs.exp.ctagac<-log(n3.codon.pair$ctagac/cps.cta.gac)
cps.obs.exp.ctagag<-log(n3.codon.pair$ctagag/cps.cta.gag)
cps.obs.exp.ctagat<-log(n3.codon.pair$ctagat/cps.cta.gat)
cps.obs.exp.ctagca<-log(n3.codon.pair$ctagca/cps.cta.gca)
cps.obs.exp.ctagcc<-log(n3.codon.pair$ctagcc/cps.cta.gcc)
cps.obs.exp.ctagcg<-log(n3.codon.pair$ctagcg/cps.cta.gcg)
cps.obs.exp.ctagct<-log(n3.codon.pair$ctagct/cps.cta.gct)
cps.obs.exp.ctagga<-log(n3.codon.pair$ctagga/cps.cta.gga)
cps.obs.exp.ctaggc<-log(n3.codon.pair$ctaggc/cps.cta.ggc)
cps.obs.exp.ctaggg<-log(n3.codon.pair$ctaggg/cps.cta.ggg)
cps.obs.exp.ctaggt<-log(n3.codon.pair$ctaggt/cps.cta.ggt)
cps.obs.exp.ctagta<-log(n3.codon.pair$ctagta/cps.cta.gta)
cps.obs.exp.ctagtc<-log(n3.codon.pair$ctagtc/cps.cta.gtc)
cps.obs.exp.ctagtg<-log(n3.codon.pair$ctagtg/cps.cta.gtg)
cps.obs.exp.ctagtt<-log(n3.codon.pair$ctagtt/cps.cta.gtt)

#cps.obs.exp.ctataa<-log(n3.codon.pair$ctataa/cps.cta.taa)
cps.obs.exp.ctatac<-log(n3.codon.pair$ctatac/cps.cta.tac)
#cps.obs.exp.ctatag<-log(n3.codon.pair$ctatag/cps.cta.tag)
cps.obs.exp.ctatat<-log(n3.codon.pair$ctatat/cps.cta.tat)
cps.obs.exp.ctatca<-log(n3.codon.pair$ctatca/cps.cta.tca)
cps.obs.exp.ctatcc<-log(n3.codon.pair$ctatcc/cps.cta.tcc)
cps.obs.exp.ctatcg<-log(n3.codon.pair$ctatcg/cps.cta.tcg)
cps.obs.exp.ctatct<-log(n3.codon.pair$ctatct/cps.cta.tct)
#cps.obs.exp.ctatga<-log(n3.codon.pair$ctatga/cps.cta.tga)
cps.obs.exp.ctatgc<-log(n3.codon.pair$ctatgc/cps.cta.tgc)
cps.obs.exp.ctatgg<-log(n3.codon.pair$ctatgg/cps.cta.tgg)
cps.obs.exp.ctatgt<-log(n3.codon.pair$ctatgt/cps.cta.tgt)
cps.obs.exp.ctatta<-log(n3.codon.pair$ctatta/cps.cta.tta)
cps.obs.exp.ctattc<-log(n3.codon.pair$ctattc/cps.cta.ttc)
cps.obs.exp.ctattg<-log(n3.codon.pair$ctattg/cps.cta.ttg)
cps.obs.exp.ctattt<-log(n3.codon.pair$ctattt/cps.cta.ttt)









cps.obs.exp.ctcaaa<-log(n3.codon.pair$ctcaaa/cps.ctc.aaa)
cps.obs.exp.ctcaac<-log(n3.codon.pair$ctcaac/cps.ctc.aac)
cps.obs.exp.ctcaag<-log(n3.codon.pair$ctcaag/cps.ctc.aag)
cps.obs.exp.ctcaat<-log(n3.codon.pair$ctcaat/cps.ctc.aat)
cps.obs.exp.ctcaca<-log(n3.codon.pair$ctcaca/cps.ctc.aca)
cps.obs.exp.ctcacc<-log(n3.codon.pair$ctcacc/cps.ctc.acc)
cps.obs.exp.ctcacg<-log(n3.codon.pair$ctcacg/cps.ctc.acg)
cps.obs.exp.ctcact<-log(n3.codon.pair$ctcact/cps.ctc.act)
cps.obs.exp.ctcaga<-log(n3.codon.pair$ctcaga/cps.ctc.aga)
cps.obs.exp.ctcagc<-log(n3.codon.pair$ctcagc/cps.ctc.agc)
cps.obs.exp.ctcagg<-log(n3.codon.pair$ctcagg/cps.ctc.agg)
cps.obs.exp.ctcagt<-log(n3.codon.pair$ctcagt/cps.ctc.agt)
cps.obs.exp.ctcata<-log(n3.codon.pair$ctcata/cps.ctc.ata)
cps.obs.exp.ctcatc<-log(n3.codon.pair$ctcatc/cps.ctc.atc)
cps.obs.exp.ctcatg<-log(n3.codon.pair$ctcatg/cps.ctc.atg)
cps.obs.exp.ctcatt<-log(n3.codon.pair$ctcatt/cps.ctc.att)

cps.obs.exp.ctccaa<-log(n3.codon.pair$ctccaa/cps.ctc.caa)
cps.obs.exp.ctccac<-log(n3.codon.pair$ctccac/cps.ctc.cac)
cps.obs.exp.ctccag<-log(n3.codon.pair$ctccag/cps.ctc.cag)
cps.obs.exp.ctccat<-log(n3.codon.pair$ctccat/cps.ctc.cat)
cps.obs.exp.ctccca<-log(n3.codon.pair$ctccca/cps.ctc.cca)
cps.obs.exp.ctcccc<-log(n3.codon.pair$ctcccc/cps.ctc.ccc)
cps.obs.exp.ctcccg<-log(n3.codon.pair$ctcccg/cps.ctc.ccg)
cps.obs.exp.ctccct<-log(n3.codon.pair$ctccct/cps.ctc.cct)
cps.obs.exp.ctccga<-log(n3.codon.pair$ctccga/cps.ctc.cga)
cps.obs.exp.ctccgc<-log(n3.codon.pair$ctccgc/cps.ctc.cgc)
cps.obs.exp.ctccgg<-log(n3.codon.pair$ctccgg/cps.ctc.cgg)
cps.obs.exp.ctccgt<-log(n3.codon.pair$ctccgt/cps.ctc.cgt)
cps.obs.exp.ctccta<-log(n3.codon.pair$ctccta/cps.ctc.cta)
cps.obs.exp.ctcctc<-log(n3.codon.pair$ctcctc/cps.ctc.ctc)
cps.obs.exp.ctcctg<-log(n3.codon.pair$ctcctg/cps.ctc.ctg)
cps.obs.exp.ctcctt<-log(n3.codon.pair$ctcctt/cps.ctc.ctt)

cps.obs.exp.ctcgaa<-log(n3.codon.pair$ctcgaa/cps.ctc.gaa)
cps.obs.exp.ctcgac<-log(n3.codon.pair$ctcgac/cps.ctc.gac)
cps.obs.exp.ctcgag<-log(n3.codon.pair$ctcgag/cps.ctc.gag)
cps.obs.exp.ctcgat<-log(n3.codon.pair$ctcgat/cps.ctc.gat)
cps.obs.exp.ctcgca<-log(n3.codon.pair$ctcgca/cps.ctc.gca)
cps.obs.exp.ctcgcc<-log(n3.codon.pair$ctcgcc/cps.ctc.gcc)
cps.obs.exp.ctcgcg<-log(n3.codon.pair$ctcgcg/cps.ctc.gcg)
cps.obs.exp.ctcgct<-log(n3.codon.pair$ctcgct/cps.ctc.gct)
cps.obs.exp.ctcgga<-log(n3.codon.pair$ctcgga/cps.ctc.gga)
cps.obs.exp.ctcggc<-log(n3.codon.pair$ctcggc/cps.ctc.ggc)
cps.obs.exp.ctcggg<-log(n3.codon.pair$ctcggg/cps.ctc.ggg)
cps.obs.exp.ctcggt<-log(n3.codon.pair$ctcggt/cps.ctc.ggt)
cps.obs.exp.ctcgta<-log(n3.codon.pair$ctcgta/cps.ctc.gta)
cps.obs.exp.ctcgtc<-log(n3.codon.pair$ctcgtc/cps.ctc.gtc)
cps.obs.exp.ctcgtg<-log(n3.codon.pair$ctcgtg/cps.ctc.gtg)
cps.obs.exp.ctcgtt<-log(n3.codon.pair$ctcgtt/cps.ctc.gtt)

#cps.obs.exp.ctctaa<-log(n3.codon.pair$ctctaa/cps.ctc.taa)
cps.obs.exp.ctctac<-log(n3.codon.pair$ctctac/cps.ctc.tac)
#cps.obs.exp.ctctag<-log(n3.codon.pair$ctctag/cps.ctc.tag)
cps.obs.exp.ctctat<-log(n3.codon.pair$ctctat/cps.ctc.tat)
cps.obs.exp.ctctca<-log(n3.codon.pair$ctctca/cps.ctc.tca)
cps.obs.exp.ctctcc<-log(n3.codon.pair$ctctcc/cps.ctc.tcc)
cps.obs.exp.ctctcg<-log(n3.codon.pair$ctctcg/cps.ctc.tcg)
cps.obs.exp.ctctct<-log(n3.codon.pair$ctctct/cps.ctc.tct)
#cps.obs.exp.ctctga<-log(n3.codon.pair$ctctga/cps.ctc.tga)
cps.obs.exp.ctctgc<-log(n3.codon.pair$ctctgc/cps.ctc.tgc)
cps.obs.exp.ctctgg<-log(n3.codon.pair$ctctgg/cps.ctc.tgg)
cps.obs.exp.ctctgt<-log(n3.codon.pair$ctctgt/cps.ctc.tgt)
cps.obs.exp.ctctta<-log(n3.codon.pair$ctctta/cps.ctc.tta)
cps.obs.exp.ctcttc<-log(n3.codon.pair$ctcttc/cps.ctc.ttc)
cps.obs.exp.ctcttg<-log(n3.codon.pair$ctcttg/cps.ctc.ttg)
cps.obs.exp.ctcttt<-log(n3.codon.pair$ctcttt/cps.ctc.ttt)









cps.obs.exp.ctgaaa<-log(n3.codon.pair$ctgaaa/cps.ctg.aaa)
cps.obs.exp.ctgaac<-log(n3.codon.pair$ctgaac/cps.ctg.aac)
cps.obs.exp.ctgaag<-log(n3.codon.pair$ctgaag/cps.ctg.aag)
cps.obs.exp.ctgaat<-log(n3.codon.pair$ctgaat/cps.ctg.aat)
cps.obs.exp.ctgaca<-log(n3.codon.pair$ctgaca/cps.ctg.aca)
cps.obs.exp.ctgacc<-log(n3.codon.pair$ctgacc/cps.ctg.acc)
cps.obs.exp.ctgacg<-log(n3.codon.pair$ctgacg/cps.ctg.acg)
cps.obs.exp.ctgact<-log(n3.codon.pair$ctgact/cps.ctg.act)
cps.obs.exp.ctgaga<-log(n3.codon.pair$ctgaga/cps.ctg.aga)
cps.obs.exp.ctgagc<-log(n3.codon.pair$ctgagc/cps.ctg.agc)
cps.obs.exp.ctgagg<-log(n3.codon.pair$ctgagg/cps.ctg.agg)
cps.obs.exp.ctgagt<-log(n3.codon.pair$ctgagt/cps.ctg.agt)
cps.obs.exp.ctgata<-log(n3.codon.pair$ctgata/cps.ctg.ata)
cps.obs.exp.ctgatc<-log(n3.codon.pair$ctgatc/cps.ctg.atc)
cps.obs.exp.ctgatg<-log(n3.codon.pair$ctgatg/cps.ctg.atg)
cps.obs.exp.ctgatt<-log(n3.codon.pair$ctgatt/cps.ctg.att)

cps.obs.exp.ctgcaa<-log(n3.codon.pair$ctgcaa/cps.ctg.caa)
cps.obs.exp.ctgcac<-log(n3.codon.pair$ctgcac/cps.ctg.cac)
cps.obs.exp.ctgcag<-log(n3.codon.pair$ctgcag/cps.ctg.cag)
cps.obs.exp.ctgcat<-log(n3.codon.pair$ctgcat/cps.ctg.cat)
cps.obs.exp.ctgcca<-log(n3.codon.pair$ctgcca/cps.ctg.cca)
cps.obs.exp.ctgccc<-log(n3.codon.pair$ctgccc/cps.ctg.ccc)
cps.obs.exp.ctgccg<-log(n3.codon.pair$ctgccg/cps.ctg.ccg)
cps.obs.exp.ctgcct<-log(n3.codon.pair$ctgcct/cps.ctg.cct)
cps.obs.exp.ctgcga<-log(n3.codon.pair$ctgcga/cps.ctg.cga)
cps.obs.exp.ctgcgc<-log(n3.codon.pair$ctgcgc/cps.ctg.cgc)
cps.obs.exp.ctgcgg<-log(n3.codon.pair$ctgcgg/cps.ctg.cgg)
cps.obs.exp.ctgcgt<-log(n3.codon.pair$ctgcgt/cps.ctg.cgt)
cps.obs.exp.ctgcta<-log(n3.codon.pair$ctgcta/cps.ctg.cta)
cps.obs.exp.ctgctc<-log(n3.codon.pair$ctgctc/cps.ctg.ctc)
cps.obs.exp.ctgctg<-log(n3.codon.pair$ctgctg/cps.ctg.ctg)
cps.obs.exp.ctgctt<-log(n3.codon.pair$ctgctt/cps.ctg.ctt)

cps.obs.exp.ctggaa<-log(n3.codon.pair$ctggaa/cps.ctg.gaa)
cps.obs.exp.ctggac<-log(n3.codon.pair$ctggac/cps.ctg.gac)
cps.obs.exp.ctggag<-log(n3.codon.pair$ctggag/cps.ctg.gag)
cps.obs.exp.ctggat<-log(n3.codon.pair$ctggat/cps.ctg.gat)
cps.obs.exp.ctggca<-log(n3.codon.pair$ctggca/cps.ctg.gca)
cps.obs.exp.ctggcc<-log(n3.codon.pair$ctggcc/cps.ctg.gcc)
cps.obs.exp.ctggcg<-log(n3.codon.pair$ctggcg/cps.ctg.gcg)
cps.obs.exp.ctggct<-log(n3.codon.pair$ctggct/cps.ctg.gct)
cps.obs.exp.ctggga<-log(n3.codon.pair$ctggga/cps.ctg.gga)
cps.obs.exp.ctgggc<-log(n3.codon.pair$ctgggc/cps.ctg.ggc)
cps.obs.exp.ctgggg<-log(n3.codon.pair$ctgggg/cps.ctg.ggg)
cps.obs.exp.ctgggt<-log(n3.codon.pair$ctgggt/cps.ctg.ggt)
cps.obs.exp.ctggta<-log(n3.codon.pair$ctggta/cps.ctg.gta)
cps.obs.exp.ctggtc<-log(n3.codon.pair$ctggtc/cps.ctg.gtc)
cps.obs.exp.ctggtg<-log(n3.codon.pair$ctggtg/cps.ctg.gtg)
cps.obs.exp.ctggtt<-log(n3.codon.pair$ctggtt/cps.ctg.gtt)

#cps.obs.exp.ctgtaa<-log(n3.codon.pair$ctgtaa/cps.ctg.taa)
cps.obs.exp.ctgtac<-log(n3.codon.pair$ctgtac/cps.ctg.tac)
#cps.obs.exp.ctgtag<-log(n3.codon.pair$ctgtag/cps.ctg.tag)
cps.obs.exp.ctgtat<-log(n3.codon.pair$ctgtat/cps.ctg.tat)
cps.obs.exp.ctgtca<-log(n3.codon.pair$ctgtca/cps.ctg.tca)
cps.obs.exp.ctgtcc<-log(n3.codon.pair$ctgtcc/cps.ctg.tcc)
cps.obs.exp.ctgtcg<-log(n3.codon.pair$ctgtcg/cps.ctg.tcg)
cps.obs.exp.ctgtct<-log(n3.codon.pair$ctgtct/cps.ctg.tct)
#cps.obs.exp.ctgtga<-log(n3.codon.pair$ctgtga/cps.ctg.tga)
cps.obs.exp.ctgtgc<-log(n3.codon.pair$ctgtgc/cps.ctg.tgc)
cps.obs.exp.ctgtgg<-log(n3.codon.pair$ctgtgg/cps.ctg.tgg)
cps.obs.exp.ctgtgt<-log(n3.codon.pair$ctgtgt/cps.ctg.tgt)
cps.obs.exp.ctgtta<-log(n3.codon.pair$ctgtta/cps.ctg.tta)
cps.obs.exp.ctgttc<-log(n3.codon.pair$ctgttc/cps.ctg.ttc)
cps.obs.exp.ctgttg<-log(n3.codon.pair$ctgttg/cps.ctg.ttg)
cps.obs.exp.ctgttt<-log(n3.codon.pair$ctgttt/cps.ctg.ttt)








cps.obs.exp.cttaaa<-log(n3.codon.pair$cttaaa/cps.ctt.aaa)
cps.obs.exp.cttaac<-log(n3.codon.pair$cttaac/cps.ctt.aac)
cps.obs.exp.cttaag<-log(n3.codon.pair$cttaag/cps.ctt.aag)
cps.obs.exp.cttaat<-log(n3.codon.pair$cttaat/cps.ctt.aat)
cps.obs.exp.cttaca<-log(n3.codon.pair$cttaca/cps.ctt.aca)
cps.obs.exp.cttacc<-log(n3.codon.pair$cttacc/cps.ctt.acc)
cps.obs.exp.cttacg<-log(n3.codon.pair$cttacg/cps.ctt.acg)
cps.obs.exp.cttact<-log(n3.codon.pair$cttact/cps.ctt.act)
cps.obs.exp.cttaga<-log(n3.codon.pair$cttaga/cps.ctt.aga)
cps.obs.exp.cttagc<-log(n3.codon.pair$cttagc/cps.ctt.agc)
cps.obs.exp.cttagg<-log(n3.codon.pair$cttagg/cps.ctt.agg)
cps.obs.exp.cttagt<-log(n3.codon.pair$cttagt/cps.ctt.agt)
cps.obs.exp.cttata<-log(n3.codon.pair$cttata/cps.ctt.ata)
cps.obs.exp.cttatc<-log(n3.codon.pair$cttatc/cps.ctt.atc)
cps.obs.exp.cttatg<-log(n3.codon.pair$cttatg/cps.ctt.atg)
cps.obs.exp.cttatt<-log(n3.codon.pair$cttatt/cps.ctt.att)

cps.obs.exp.cttcaa<-log(n3.codon.pair$cttcaa/cps.ctt.caa)
cps.obs.exp.cttcac<-log(n3.codon.pair$cttcac/cps.ctt.cac)
cps.obs.exp.cttcag<-log(n3.codon.pair$cttcag/cps.ctt.cag)
cps.obs.exp.cttcat<-log(n3.codon.pair$cttcat/cps.ctt.cat)
cps.obs.exp.cttcca<-log(n3.codon.pair$cttcca/cps.ctt.cca)
cps.obs.exp.cttccc<-log(n3.codon.pair$cttccc/cps.ctt.ccc)
cps.obs.exp.cttccg<-log(n3.codon.pair$cttccg/cps.ctt.ccg)
cps.obs.exp.cttcct<-log(n3.codon.pair$cttcct/cps.ctt.cct)
cps.obs.exp.cttcga<-log(n3.codon.pair$cttcga/cps.ctt.cga)
cps.obs.exp.cttcgc<-log(n3.codon.pair$cttcgc/cps.ctt.cgc)
cps.obs.exp.cttcgg<-log(n3.codon.pair$cttcgg/cps.ctt.cgg)
cps.obs.exp.cttcgt<-log(n3.codon.pair$cttcgt/cps.ctt.cgt)
cps.obs.exp.cttcta<-log(n3.codon.pair$cttcta/cps.ctt.cta)
cps.obs.exp.cttctc<-log(n3.codon.pair$cttctc/cps.ctt.ctc)
cps.obs.exp.cttctg<-log(n3.codon.pair$cttctg/cps.ctt.ctg)
cps.obs.exp.cttctt<-log(n3.codon.pair$cttctt/cps.ctt.ctt)

cps.obs.exp.cttgaa<-log(n3.codon.pair$cttgaa/cps.ctt.gaa)
cps.obs.exp.cttgac<-log(n3.codon.pair$cttgac/cps.ctt.gac)
cps.obs.exp.cttgag<-log(n3.codon.pair$cttgag/cps.ctt.gag)
cps.obs.exp.cttgat<-log(n3.codon.pair$cttgat/cps.ctt.gat)
cps.obs.exp.cttgca<-log(n3.codon.pair$cttgca/cps.ctt.gca)
cps.obs.exp.cttgcc<-log(n3.codon.pair$cttgcc/cps.ctt.gcc)
cps.obs.exp.cttgcg<-log(n3.codon.pair$cttgcg/cps.ctt.gcg)
cps.obs.exp.cttgct<-log(n3.codon.pair$cttgct/cps.ctt.gct)
cps.obs.exp.cttgga<-log(n3.codon.pair$cttgga/cps.ctt.gga)
cps.obs.exp.cttggc<-log(n3.codon.pair$cttggc/cps.ctt.ggc)
cps.obs.exp.cttggg<-log(n3.codon.pair$cttggg/cps.ctt.ggg)
cps.obs.exp.cttggt<-log(n3.codon.pair$cttggt/cps.ctt.ggt)
cps.obs.exp.cttgta<-log(n3.codon.pair$cttgta/cps.ctt.gta)
cps.obs.exp.cttgtc<-log(n3.codon.pair$cttgtc/cps.ctt.gtc)
cps.obs.exp.cttgtg<-log(n3.codon.pair$cttgtg/cps.ctt.gtg)
cps.obs.exp.cttgtt<-log(n3.codon.pair$cttgtt/cps.ctt.gtt)

#cps.obs.exp.ctttaa<-log(n3.codon.pair$ctttaa/cps.ctt.taa)
cps.obs.exp.ctttac<-log(n3.codon.pair$ctttac/cps.ctt.tac)
#cps.obs.exp.ctttag<-log(n3.codon.pair$ctttag/cps.ctt.tag)
cps.obs.exp.ctttat<-log(n3.codon.pair$ctttat/cps.ctt.tat)
cps.obs.exp.ctttca<-log(n3.codon.pair$ctttca/cps.ctt.tca)
cps.obs.exp.ctttcc<-log(n3.codon.pair$ctttcc/cps.ctt.tcc)
cps.obs.exp.ctttcg<-log(n3.codon.pair$ctttcg/cps.ctt.tcg)
cps.obs.exp.ctttct<-log(n3.codon.pair$ctttct/cps.ctt.tct)
#cps.obs.exp.ctttga<-log(n3.codon.pair$ctttga/cps.ctt.tga)
cps.obs.exp.ctttgc<-log(n3.codon.pair$ctttgc/cps.ctt.tgc)
cps.obs.exp.ctttgg<-log(n3.codon.pair$ctttgg/cps.ctt.tgg)
cps.obs.exp.ctttgt<-log(n3.codon.pair$ctttgt/cps.ctt.tgt)
cps.obs.exp.ctttta<-log(n3.codon.pair$ctttta/cps.ctt.tta)
cps.obs.exp.cttttc<-log(n3.codon.pair$cttttc/cps.ctt.ttc)
cps.obs.exp.cttttg<-log(n3.codon.pair$cttttg/cps.ctt.ttg)
cps.obs.exp.cttttt<-log(n3.codon.pair$cttttt/cps.ctt.ttt)




















cps.obs.exp.gaaaaa<-log(n3.codon.pair$gaaaaa/cps.gaa.aaa)
cps.obs.exp.gaaaac<-log(n3.codon.pair$gaaaac/cps.gaa.aac)
cps.obs.exp.gaaaag<-log(n3.codon.pair$gaaaag/cps.gaa.aag)
cps.obs.exp.gaaaat<-log(n3.codon.pair$gaaaat/cps.gaa.aat)
cps.obs.exp.gaaaca<-log(n3.codon.pair$gaaaca/cps.gaa.aca)
cps.obs.exp.gaaacc<-log(n3.codon.pair$gaaacc/cps.gaa.acc)
cps.obs.exp.gaaacg<-log(n3.codon.pair$gaaacg/cps.gaa.acg)
cps.obs.exp.gaaact<-log(n3.codon.pair$gaaact/cps.gaa.act)
cps.obs.exp.gaaaga<-log(n3.codon.pair$gaaaga/cps.gaa.aga)
cps.obs.exp.gaaagc<-log(n3.codon.pair$gaaagc/cps.gaa.agc)
cps.obs.exp.gaaagg<-log(n3.codon.pair$gaaagg/cps.gaa.agg)
cps.obs.exp.gaaagt<-log(n3.codon.pair$gaaagt/cps.gaa.agt)
cps.obs.exp.gaaata<-log(n3.codon.pair$gaaata/cps.gaa.ata)
cps.obs.exp.gaaatc<-log(n3.codon.pair$gaaatc/cps.gaa.atc)
cps.obs.exp.gaaatg<-log(n3.codon.pair$gaaatg/cps.gaa.atg)
cps.obs.exp.gaaatt<-log(n3.codon.pair$gaaatt/cps.gaa.att)

cps.obs.exp.gaacaa<-log(n3.codon.pair$gaacaa/cps.gaa.caa)
cps.obs.exp.gaacac<-log(n3.codon.pair$gaacac/cps.gaa.cac)
cps.obs.exp.gaacag<-log(n3.codon.pair$gaacag/cps.gaa.cag)
cps.obs.exp.gaacat<-log(n3.codon.pair$gaacat/cps.gaa.cat)
cps.obs.exp.gaacca<-log(n3.codon.pair$gaacca/cps.gaa.cca)
cps.obs.exp.gaaccc<-log(n3.codon.pair$gaaccc/cps.gaa.ccc)
cps.obs.exp.gaaccg<-log(n3.codon.pair$gaaccg/cps.gaa.ccg)
cps.obs.exp.gaacct<-log(n3.codon.pair$gaacct/cps.gaa.cct)
cps.obs.exp.gaacga<-log(n3.codon.pair$gaacga/cps.gaa.cga)
cps.obs.exp.gaacgc<-log(n3.codon.pair$gaacgc/cps.gaa.cgc)
cps.obs.exp.gaacgg<-log(n3.codon.pair$gaacgg/cps.gaa.cgg)
cps.obs.exp.gaacgt<-log(n3.codon.pair$gaacgt/cps.gaa.cgt)
cps.obs.exp.gaacta<-log(n3.codon.pair$gaacta/cps.gaa.cta)
cps.obs.exp.gaactc<-log(n3.codon.pair$gaactc/cps.gaa.ctc)
cps.obs.exp.gaactg<-log(n3.codon.pair$gaactg/cps.gaa.ctg)
cps.obs.exp.gaactt<-log(n3.codon.pair$gaactt/cps.gaa.ctt)

cps.obs.exp.gaagaa<-log(n3.codon.pair$gaagaa/cps.gaa.gaa)
cps.obs.exp.gaagac<-log(n3.codon.pair$gaagac/cps.gaa.gac)
cps.obs.exp.gaagag<-log(n3.codon.pair$gaagag/cps.gaa.gag)
cps.obs.exp.gaagat<-log(n3.codon.pair$gaagat/cps.gaa.gat)
cps.obs.exp.gaagca<-log(n3.codon.pair$gaagca/cps.gaa.gca)
cps.obs.exp.gaagcc<-log(n3.codon.pair$gaagcc/cps.gaa.gcc)
cps.obs.exp.gaagcg<-log(n3.codon.pair$gaagcg/cps.gaa.gcg)
cps.obs.exp.gaagct<-log(n3.codon.pair$gaagct/cps.gaa.gct)
cps.obs.exp.gaagga<-log(n3.codon.pair$gaagga/cps.gaa.gga)
cps.obs.exp.gaaggc<-log(n3.codon.pair$gaaggc/cps.gaa.ggc)
cps.obs.exp.gaaggg<-log(n3.codon.pair$gaaggg/cps.gaa.ggg)
cps.obs.exp.gaaggt<-log(n3.codon.pair$gaaggt/cps.gaa.ggt)
cps.obs.exp.gaagta<-log(n3.codon.pair$gaagta/cps.gaa.gta)
cps.obs.exp.gaagtc<-log(n3.codon.pair$gaagtc/cps.gaa.gtc)
cps.obs.exp.gaagtg<-log(n3.codon.pair$gaagtg/cps.gaa.gtg)
cps.obs.exp.gaagtt<-log(n3.codon.pair$gaagtt/cps.gaa.gtt)

#cps.obs.exp.gaataa<-log(n3.codon.pair$gaataa/cps.gaa.taa)
cps.obs.exp.gaatac<-log(n3.codon.pair$gaatac/cps.gaa.tac)
#cps.obs.exp.gaatag<-log(n3.codon.pair$gaatag/cps.gaa.tag)
cps.obs.exp.gaatat<-log(n3.codon.pair$gaatat/cps.gaa.tat)
cps.obs.exp.gaatca<-log(n3.codon.pair$gaatca/cps.gaa.tca)
cps.obs.exp.gaatcc<-log(n3.codon.pair$gaatcc/cps.gaa.tcc)
cps.obs.exp.gaatcg<-log(n3.codon.pair$gaatcg/cps.gaa.tcg)
cps.obs.exp.gaatct<-log(n3.codon.pair$gaatct/cps.gaa.tct)
#cps.obs.exp.gaatga<-log(n3.codon.pair$gaatga/cps.gaa.tga)
cps.obs.exp.gaatgc<-log(n3.codon.pair$gaatgc/cps.gaa.tgc)
cps.obs.exp.gaatgg<-log(n3.codon.pair$gaatgg/cps.gaa.tgg)
cps.obs.exp.gaatgt<-log(n3.codon.pair$gaatgt/cps.gaa.tgt)
cps.obs.exp.gaatta<-log(n3.codon.pair$gaatta/cps.gaa.tta)
cps.obs.exp.gaattc<-log(n3.codon.pair$gaattc/cps.gaa.ttc)
cps.obs.exp.gaattg<-log(n3.codon.pair$gaattg/cps.gaa.ttg)
cps.obs.exp.gaattt<-log(n3.codon.pair$gaattt/cps.gaa.ttt)









cps.obs.exp.gacaaa<-log(n3.codon.pair$gacaaa/cps.gac.aaa)
cps.obs.exp.gacaac<-log(n3.codon.pair$gacaac/cps.gac.aac)
cps.obs.exp.gacaag<-log(n3.codon.pair$gacaag/cps.gac.aag)
cps.obs.exp.gacaat<-log(n3.codon.pair$gacaat/cps.gac.aat)
cps.obs.exp.gacaca<-log(n3.codon.pair$gacaca/cps.gac.aca)
cps.obs.exp.gacacc<-log(n3.codon.pair$gacacc/cps.gac.acc)
cps.obs.exp.gacacg<-log(n3.codon.pair$gacacg/cps.gac.acg)
cps.obs.exp.gacact<-log(n3.codon.pair$gacact/cps.gac.act)
cps.obs.exp.gacaga<-log(n3.codon.pair$gacaga/cps.gac.aga)
cps.obs.exp.gacagc<-log(n3.codon.pair$gacagc/cps.gac.agc)
cps.obs.exp.gacagg<-log(n3.codon.pair$gacagg/cps.gac.agg)
cps.obs.exp.gacagt<-log(n3.codon.pair$gacagt/cps.gac.agt)
cps.obs.exp.gacata<-log(n3.codon.pair$gacata/cps.gac.ata)
cps.obs.exp.gacatc<-log(n3.codon.pair$gacatc/cps.gac.atc)
cps.obs.exp.gacatg<-log(n3.codon.pair$gacatg/cps.gac.atg)
cps.obs.exp.gacatt<-log(n3.codon.pair$gacatt/cps.gac.att)

cps.obs.exp.gaccaa<-log(n3.codon.pair$gaccaa/cps.gac.caa)
cps.obs.exp.gaccac<-log(n3.codon.pair$gaccac/cps.gac.cac)
cps.obs.exp.gaccag<-log(n3.codon.pair$gaccag/cps.gac.cag)
cps.obs.exp.gaccat<-log(n3.codon.pair$gaccat/cps.gac.cat)
cps.obs.exp.gaccca<-log(n3.codon.pair$gaccca/cps.gac.cca)
cps.obs.exp.gacccc<-log(n3.codon.pair$gacccc/cps.gac.ccc)
cps.obs.exp.gacccg<-log(n3.codon.pair$gacccg/cps.gac.ccg)
cps.obs.exp.gaccct<-log(n3.codon.pair$gaccct/cps.gac.cct)
cps.obs.exp.gaccga<-log(n3.codon.pair$gaccga/cps.gac.cga)
cps.obs.exp.gaccgc<-log(n3.codon.pair$gaccgc/cps.gac.cgc)
cps.obs.exp.gaccgg<-log(n3.codon.pair$gaccgg/cps.gac.cgg)
cps.obs.exp.gaccgt<-log(n3.codon.pair$gaccgt/cps.gac.cgt)
cps.obs.exp.gaccta<-log(n3.codon.pair$gaccta/cps.gac.cta)
cps.obs.exp.gacctc<-log(n3.codon.pair$gacctc/cps.gac.ctc)
cps.obs.exp.gacctg<-log(n3.codon.pair$gacctg/cps.gac.ctg)
cps.obs.exp.gacctt<-log(n3.codon.pair$gacctt/cps.gac.ctt)

cps.obs.exp.gacgaa<-log(n3.codon.pair$gacgaa/cps.gac.gaa)
cps.obs.exp.gacgac<-log(n3.codon.pair$gacgac/cps.gac.gac)
cps.obs.exp.gacgag<-log(n3.codon.pair$gacgag/cps.gac.gag)
cps.obs.exp.gacgat<-log(n3.codon.pair$gacgat/cps.gac.gat)
cps.obs.exp.gacgca<-log(n3.codon.pair$gacgca/cps.gac.gca)
cps.obs.exp.gacgcc<-log(n3.codon.pair$gacgcc/cps.gac.gcc)
cps.obs.exp.gacgcg<-log(n3.codon.pair$gacgcg/cps.gac.gcg)
cps.obs.exp.gacgct<-log(n3.codon.pair$gacgct/cps.gac.gct)
cps.obs.exp.gacgga<-log(n3.codon.pair$gacgga/cps.gac.gga)
cps.obs.exp.gacggc<-log(n3.codon.pair$gacggc/cps.gac.ggc)
cps.obs.exp.gacggg<-log(n3.codon.pair$gacggg/cps.gac.ggg)
cps.obs.exp.gacggt<-log(n3.codon.pair$gacggt/cps.gac.ggt)
cps.obs.exp.gacgta<-log(n3.codon.pair$gacgta/cps.gac.gta)
cps.obs.exp.gacgtc<-log(n3.codon.pair$gacgtc/cps.gac.gtc)
cps.obs.exp.gacgtg<-log(n3.codon.pair$gacgtg/cps.gac.gtg)
cps.obs.exp.gacgtt<-log(n3.codon.pair$gacgtt/cps.gac.gtt)

#cps.obs.exp.gactaa<-log(n3.codon.pair$gactaa/cps.gac.taa)
cps.obs.exp.gactac<-log(n3.codon.pair$gactac/cps.gac.tac)
#cps.obs.exp.gactag<-log(n3.codon.pair$gactag/cps.gac.tag)
cps.obs.exp.gactat<-log(n3.codon.pair$gactat/cps.gac.tat)
cps.obs.exp.gactca<-log(n3.codon.pair$gactca/cps.gac.tca)
cps.obs.exp.gactcc<-log(n3.codon.pair$gactcc/cps.gac.tcc)
cps.obs.exp.gactcg<-log(n3.codon.pair$gactcg/cps.gac.tcg)
cps.obs.exp.gactct<-log(n3.codon.pair$gactct/cps.gac.tct)
#cps.obs.exp.gactga<-log(n3.codon.pair$gactga/cps.gac.tga)
cps.obs.exp.gactgc<-log(n3.codon.pair$gactgc/cps.gac.tgc)
cps.obs.exp.gactgg<-log(n3.codon.pair$gactgg/cps.gac.tgg)
cps.obs.exp.gactgt<-log(n3.codon.pair$gactgt/cps.gac.tgt)
cps.obs.exp.gactta<-log(n3.codon.pair$gactta/cps.gac.tta)
cps.obs.exp.gacttc<-log(n3.codon.pair$gacttc/cps.gac.ttc)
cps.obs.exp.gacttg<-log(n3.codon.pair$gacttg/cps.gac.ttg)
cps.obs.exp.gacttt<-log(n3.codon.pair$gacttt/cps.gac.ttt)










cps.obs.exp.gagaaa<-log(n3.codon.pair$gagaaa/cps.gag.aaa)
cps.obs.exp.gagaac<-log(n3.codon.pair$gagaac/cps.gag.aac)
cps.obs.exp.gagaag<-log(n3.codon.pair$gagaag/cps.gag.aag)
cps.obs.exp.gagaat<-log(n3.codon.pair$gagaat/cps.gag.aat)
cps.obs.exp.gagaca<-log(n3.codon.pair$gagaca/cps.gag.aca)
cps.obs.exp.gagacc<-log(n3.codon.pair$gagacc/cps.gag.acc)
cps.obs.exp.gagacg<-log(n3.codon.pair$gagacg/cps.gag.acg)
cps.obs.exp.gagact<-log(n3.codon.pair$gagact/cps.gag.act)
cps.obs.exp.gagaga<-log(n3.codon.pair$gagaga/cps.gag.aga)
cps.obs.exp.gagagc<-log(n3.codon.pair$gagagc/cps.gag.agc)
cps.obs.exp.gagagg<-log(n3.codon.pair$gagagg/cps.gag.agg)
cps.obs.exp.gagagt<-log(n3.codon.pair$gagagt/cps.gag.agt)
cps.obs.exp.gagata<-log(n3.codon.pair$gagata/cps.gag.ata)
cps.obs.exp.gagatc<-log(n3.codon.pair$gagatc/cps.gag.atc)
cps.obs.exp.gagatg<-log(n3.codon.pair$gagatg/cps.gag.atg)
cps.obs.exp.gagatt<-log(n3.codon.pair$gagatt/cps.gag.att)

cps.obs.exp.gagcaa<-log(n3.codon.pair$gagcaa/cps.gag.caa)
cps.obs.exp.gagcac<-log(n3.codon.pair$gagcac/cps.gag.cac)
cps.obs.exp.gagcag<-log(n3.codon.pair$gagcag/cps.gag.cag)
cps.obs.exp.gagcat<-log(n3.codon.pair$gagcat/cps.gag.cat)
cps.obs.exp.gagcca<-log(n3.codon.pair$gagcca/cps.gag.cca)
cps.obs.exp.gagccc<-log(n3.codon.pair$gagccc/cps.gag.ccc)
cps.obs.exp.gagccg<-log(n3.codon.pair$gagccg/cps.gag.ccg)
cps.obs.exp.gagcct<-log(n3.codon.pair$gagcct/cps.gag.cct)
cps.obs.exp.gagcga<-log(n3.codon.pair$gagcga/cps.gag.cga)
cps.obs.exp.gagcgc<-log(n3.codon.pair$gagcgc/cps.gag.cgc)
cps.obs.exp.gagcgg<-log(n3.codon.pair$gagcgg/cps.gag.cgg)
cps.obs.exp.gagcgt<-log(n3.codon.pair$gagcgt/cps.gag.cgt)
cps.obs.exp.gagcta<-log(n3.codon.pair$gagcta/cps.gag.cta)
cps.obs.exp.gagctc<-log(n3.codon.pair$gagctc/cps.gag.ctc)
cps.obs.exp.gagctg<-log(n3.codon.pair$gagctg/cps.gag.ctg)
cps.obs.exp.gagctt<-log(n3.codon.pair$gagctt/cps.gag.ctt)

cps.obs.exp.gaggaa<-log(n3.codon.pair$gaggaa/cps.gag.gaa)
cps.obs.exp.gaggac<-log(n3.codon.pair$gaggac/cps.gag.gac)
cps.obs.exp.gaggag<-log(n3.codon.pair$gaggag/cps.gag.gag)
cps.obs.exp.gaggat<-log(n3.codon.pair$gaggat/cps.gag.gat)
cps.obs.exp.gaggca<-log(n3.codon.pair$gaggca/cps.gag.gca)
cps.obs.exp.gaggcc<-log(n3.codon.pair$gaggcc/cps.gag.gcc)
cps.obs.exp.gaggcg<-log(n3.codon.pair$gaggcg/cps.gag.gcg)
cps.obs.exp.gaggct<-log(n3.codon.pair$gaggct/cps.gag.gct)
cps.obs.exp.gaggga<-log(n3.codon.pair$gaggga/cps.gag.gga)
cps.obs.exp.gagggc<-log(n3.codon.pair$gagggc/cps.gag.ggc)
cps.obs.exp.gagggg<-log(n3.codon.pair$gagggg/cps.gag.ggg)
cps.obs.exp.gagggt<-log(n3.codon.pair$gagggt/cps.gag.ggt)
cps.obs.exp.gaggta<-log(n3.codon.pair$gaggta/cps.gag.gta)
cps.obs.exp.gaggtc<-log(n3.codon.pair$gaggtc/cps.gag.gtc)
cps.obs.exp.gaggtg<-log(n3.codon.pair$gaggtg/cps.gag.gtg)
cps.obs.exp.gaggtt<-log(n3.codon.pair$gaggtt/cps.gag.gtt)

#cps.obs.exp.gagtaa<-log(n3.codon.pair$gagtaa/cps.gag.taa)
cps.obs.exp.gagtac<-log(n3.codon.pair$gagtac/cps.gag.tac)
#cps.obs.exp.gagtag<-log(n3.codon.pair$gagtag/cps.gag.tag)
cps.obs.exp.gagtat<-log(n3.codon.pair$gagtat/cps.gag.tat)
cps.obs.exp.gagtca<-log(n3.codon.pair$gagtca/cps.gag.tca)
cps.obs.exp.gagtcc<-log(n3.codon.pair$gagtcc/cps.gag.tcc)
cps.obs.exp.gagtcg<-log(n3.codon.pair$gagtcg/cps.gag.tcg)
cps.obs.exp.gagtct<-log(n3.codon.pair$gagtct/cps.gag.tct)
#cps.obs.exp.gagtga<-log(n3.codon.pair$gagtga/cps.gag.tga)
cps.obs.exp.gagtgc<-log(n3.codon.pair$gagtgc/cps.gag.tgc)
cps.obs.exp.gagtgg<-log(n3.codon.pair$gagtgg/cps.gag.tgg)
cps.obs.exp.gagtgt<-log(n3.codon.pair$gagtgt/cps.gag.tgt)
cps.obs.exp.gagtta<-log(n3.codon.pair$gagtta/cps.gag.tta)
cps.obs.exp.gagttc<-log(n3.codon.pair$gagttc/cps.gag.ttc)
cps.obs.exp.gagttg<-log(n3.codon.pair$gagttg/cps.gag.ttg)
cps.obs.exp.gagttt<-log(n3.codon.pair$gagttt/cps.gag.ttt)









cps.obs.exp.gataaa<-log(n3.codon.pair$gataaa/cps.gat.aaa)
cps.obs.exp.gataac<-log(n3.codon.pair$gataac/cps.gat.aac)
cps.obs.exp.gataag<-log(n3.codon.pair$gataag/cps.gat.aag)
cps.obs.exp.gataat<-log(n3.codon.pair$gataat/cps.gat.aat)
cps.obs.exp.gataca<-log(n3.codon.pair$gataca/cps.gat.aca)
cps.obs.exp.gatacc<-log(n3.codon.pair$gatacc/cps.gat.acc)
cps.obs.exp.gatacg<-log(n3.codon.pair$gatacg/cps.gat.acg)
cps.obs.exp.gatact<-log(n3.codon.pair$gatact/cps.gat.act)
cps.obs.exp.gataga<-log(n3.codon.pair$gataga/cps.gat.aga)
cps.obs.exp.gatagc<-log(n3.codon.pair$gatagc/cps.gat.agc)
cps.obs.exp.gatagg<-log(n3.codon.pair$gatagg/cps.gat.agg)
cps.obs.exp.gatagt<-log(n3.codon.pair$gatagt/cps.gat.agt)
cps.obs.exp.gatata<-log(n3.codon.pair$gatata/cps.gat.ata)
cps.obs.exp.gatatc<-log(n3.codon.pair$gatatc/cps.gat.atc)
cps.obs.exp.gatatg<-log(n3.codon.pair$gatatg/cps.gat.atg)
cps.obs.exp.gatatt<-log(n3.codon.pair$gatatt/cps.gat.att)

cps.obs.exp.gatcaa<-log(n3.codon.pair$gatcaa/cps.gat.caa)
cps.obs.exp.gatcac<-log(n3.codon.pair$gatcac/cps.gat.cac)
cps.obs.exp.gatcag<-log(n3.codon.pair$gatcag/cps.gat.cag)
cps.obs.exp.gatcat<-log(n3.codon.pair$gatcat/cps.gat.cat)
cps.obs.exp.gatcca<-log(n3.codon.pair$gatcca/cps.gat.cca)
cps.obs.exp.gatccc<-log(n3.codon.pair$gatccc/cps.gat.ccc)
cps.obs.exp.gatccg<-log(n3.codon.pair$gatccg/cps.gat.ccg)
cps.obs.exp.gatcct<-log(n3.codon.pair$gatcct/cps.gat.cct)
cps.obs.exp.gatcga<-log(n3.codon.pair$gatcga/cps.gat.cga)
cps.obs.exp.gatcgc<-log(n3.codon.pair$gatcgc/cps.gat.cgc)
cps.obs.exp.gatcgg<-log(n3.codon.pair$gatcgg/cps.gat.cgg)
cps.obs.exp.gatcgt<-log(n3.codon.pair$gatcgt/cps.gat.cgt)
cps.obs.exp.gatcta<-log(n3.codon.pair$gatcta/cps.gat.cta)
cps.obs.exp.gatctc<-log(n3.codon.pair$gatctc/cps.gat.ctc)
cps.obs.exp.gatctg<-log(n3.codon.pair$gatctg/cps.gat.ctg)
cps.obs.exp.gatctt<-log(n3.codon.pair$gatctt/cps.gat.ctt)

cps.obs.exp.gatgaa<-log(n3.codon.pair$gatgaa/cps.gat.gaa)
cps.obs.exp.gatgac<-log(n3.codon.pair$gatgac/cps.gat.gac)
cps.obs.exp.gatgag<-log(n3.codon.pair$gatgag/cps.gat.gag)
cps.obs.exp.gatgat<-log(n3.codon.pair$gatgat/cps.gat.gat)
cps.obs.exp.gatgca<-log(n3.codon.pair$gatgca/cps.gat.gca)
cps.obs.exp.gatgcc<-log(n3.codon.pair$gatgcc/cps.gat.gcc)
cps.obs.exp.gatgcg<-log(n3.codon.pair$gatgcg/cps.gat.gcg)
cps.obs.exp.gatgct<-log(n3.codon.pair$gatgct/cps.gat.gct)
cps.obs.exp.gatgga<-log(n3.codon.pair$gatgga/cps.gat.gga)
cps.obs.exp.gatggc<-log(n3.codon.pair$gatggc/cps.gat.ggc)
cps.obs.exp.gatggg<-log(n3.codon.pair$gatggg/cps.gat.ggg)
cps.obs.exp.gatggt<-log(n3.codon.pair$gatggt/cps.gat.ggt)
cps.obs.exp.gatgta<-log(n3.codon.pair$gatgta/cps.gat.gta)
cps.obs.exp.gatgtc<-log(n3.codon.pair$gatgtc/cps.gat.gtc)
cps.obs.exp.gatgtg<-log(n3.codon.pair$gatgtg/cps.gat.gtg)
cps.obs.exp.gatgtt<-log(n3.codon.pair$gatgtt/cps.gat.gtt)

#cps.obs.exp.gattaa<-log(n3.codon.pair$gattaa/cps.gat.taa)
cps.obs.exp.gattac<-log(n3.codon.pair$gattac/cps.gat.tac)
#cps.obs.exp.gattag<-log(n3.codon.pair$gattag/cps.gat.tag)
cps.obs.exp.gattat<-log(n3.codon.pair$gattat/cps.gat.tat)
cps.obs.exp.gattca<-log(n3.codon.pair$gattca/cps.gat.tca)
cps.obs.exp.gattcc<-log(n3.codon.pair$gattcc/cps.gat.tcc)
cps.obs.exp.gattcg<-log(n3.codon.pair$gattcg/cps.gat.tcg)
cps.obs.exp.gattct<-log(n3.codon.pair$gattct/cps.gat.tct)
#cps.obs.exp.gattga<-log(n3.codon.pair$gattga/cps.gat.tga)
cps.obs.exp.gattgc<-log(n3.codon.pair$gattgc/cps.gat.tgc)
cps.obs.exp.gattgg<-log(n3.codon.pair$gattgg/cps.gat.tgg)
cps.obs.exp.gattgt<-log(n3.codon.pair$gattgt/cps.gat.tgt)
cps.obs.exp.gattta<-log(n3.codon.pair$gattta/cps.gat.tta)
cps.obs.exp.gatttc<-log(n3.codon.pair$gatttc/cps.gat.ttc)
cps.obs.exp.gatttg<-log(n3.codon.pair$gatttg/cps.gat.ttg)
cps.obs.exp.gatttt<-log(n3.codon.pair$gatttt/cps.gat.ttt)

















cps.obs.exp.gcaaaa<-log(n3.codon.pair$gcaaaa/cps.gca.aaa)
cps.obs.exp.gcaaac<-log(n3.codon.pair$gcaaac/cps.gca.aac)
cps.obs.exp.gcaaag<-log(n3.codon.pair$gcaaag/cps.gca.aag)
cps.obs.exp.gcaaat<-log(n3.codon.pair$gcaaat/cps.gca.aat)
cps.obs.exp.gcaaca<-log(n3.codon.pair$gcaaca/cps.gca.aca)
cps.obs.exp.gcaacc<-log(n3.codon.pair$gcaacc/cps.gca.acc)
cps.obs.exp.gcaacg<-log(n3.codon.pair$gcaacg/cps.gca.acg)
cps.obs.exp.gcaact<-log(n3.codon.pair$gcaact/cps.gca.act)
cps.obs.exp.gcaaga<-log(n3.codon.pair$gcaaga/cps.gca.aga)
cps.obs.exp.gcaagc<-log(n3.codon.pair$gcaagc/cps.gca.agc)
cps.obs.exp.gcaagg<-log(n3.codon.pair$gcaagg/cps.gca.agg)
cps.obs.exp.gcaagt<-log(n3.codon.pair$gcaagt/cps.gca.agt)
cps.obs.exp.gcaata<-log(n3.codon.pair$gcaata/cps.gca.ata)
cps.obs.exp.gcaatc<-log(n3.codon.pair$gcaatc/cps.gca.atc)
cps.obs.exp.gcaatg<-log(n3.codon.pair$gcaatg/cps.gca.atg)
cps.obs.exp.gcaatt<-log(n3.codon.pair$gcaatt/cps.gca.att)

cps.obs.exp.gcacaa<-log(n3.codon.pair$gcacaa/cps.gca.caa)
cps.obs.exp.gcacac<-log(n3.codon.pair$gcacac/cps.gca.cac)
cps.obs.exp.gcacag<-log(n3.codon.pair$gcacag/cps.gca.cag)
cps.obs.exp.gcacat<-log(n3.codon.pair$gcacat/cps.gca.cat)
cps.obs.exp.gcacca<-log(n3.codon.pair$gcacca/cps.gca.cca)
cps.obs.exp.gcaccc<-log(n3.codon.pair$gcaccc/cps.gca.ccc)
cps.obs.exp.gcaccg<-log(n3.codon.pair$gcaccg/cps.gca.ccg)
cps.obs.exp.gcacct<-log(n3.codon.pair$gcacct/cps.gca.cct)
cps.obs.exp.gcacga<-log(n3.codon.pair$gcacga/cps.gca.cga)
cps.obs.exp.gcacgc<-log(n3.codon.pair$gcacgc/cps.gca.cgc)
cps.obs.exp.gcacgg<-log(n3.codon.pair$gcacgg/cps.gca.cgg)
cps.obs.exp.gcacgt<-log(n3.codon.pair$gcacgt/cps.gca.cgt)
cps.obs.exp.gcacta<-log(n3.codon.pair$gcacta/cps.gca.cta)
cps.obs.exp.gcactc<-log(n3.codon.pair$gcactc/cps.gca.ctc)
cps.obs.exp.gcactg<-log(n3.codon.pair$gcactg/cps.gca.ctg)
cps.obs.exp.gcactt<-log(n3.codon.pair$gcactt/cps.gca.ctt)

cps.obs.exp.gcagaa<-log(n3.codon.pair$gcagaa/cps.gca.gaa)
cps.obs.exp.gcagac<-log(n3.codon.pair$gcagac/cps.gca.gac)
cps.obs.exp.gcagag<-log(n3.codon.pair$gcagag/cps.gca.gag)
cps.obs.exp.gcagat<-log(n3.codon.pair$gcagat/cps.gca.gat)
cps.obs.exp.gcagca<-log(n3.codon.pair$gcagca/cps.gca.gca)
cps.obs.exp.gcagcc<-log(n3.codon.pair$gcagcc/cps.gca.gcc)
cps.obs.exp.gcagcg<-log(n3.codon.pair$gcagcg/cps.gca.gcg)
cps.obs.exp.gcagct<-log(n3.codon.pair$gcagct/cps.gca.gct)
cps.obs.exp.gcagga<-log(n3.codon.pair$gcagga/cps.gca.gga)
cps.obs.exp.gcaggc<-log(n3.codon.pair$gcaggc/cps.gca.ggc)
cps.obs.exp.gcaggg<-log(n3.codon.pair$gcaggg/cps.gca.ggg)
cps.obs.exp.gcaggt<-log(n3.codon.pair$gcaggt/cps.gca.ggt)
cps.obs.exp.gcagta<-log(n3.codon.pair$gcagta/cps.gca.gta)
cps.obs.exp.gcagtc<-log(n3.codon.pair$gcagtc/cps.gca.gtc)
cps.obs.exp.gcagtg<-log(n3.codon.pair$gcagtg/cps.gca.gtg)
cps.obs.exp.gcagtt<-log(n3.codon.pair$gcagtt/cps.gca.gtt)

#cps.obs.exp.gcataa<-log(n3.codon.pair$gcataa/cps.gca.taa)
cps.obs.exp.gcatac<-log(n3.codon.pair$gcatac/cps.gca.tac)
#cps.obs.exp.gcatag<-log(n3.codon.pair$gcatag/cps.gca.tag)
cps.obs.exp.gcatat<-log(n3.codon.pair$gcatat/cps.gca.tat)
cps.obs.exp.gcatca<-log(n3.codon.pair$gcatca/cps.gca.tca)
cps.obs.exp.gcatcc<-log(n3.codon.pair$gcatcc/cps.gca.tcc)
cps.obs.exp.gcatcg<-log(n3.codon.pair$gcatcg/cps.gca.tcg)
cps.obs.exp.gcatct<-log(n3.codon.pair$gcatct/cps.gca.tct)
#cps.obs.exp.gcatga<-log(n3.codon.pair$gcatga/cps.gca.tga)
cps.obs.exp.gcatgc<-log(n3.codon.pair$gcatgc/cps.gca.tgc)
cps.obs.exp.gcatgg<-log(n3.codon.pair$gcatgg/cps.gca.tgg)
cps.obs.exp.gcatgt<-log(n3.codon.pair$gcatgt/cps.gca.tgt)
cps.obs.exp.gcatta<-log(n3.codon.pair$gcatta/cps.gca.tta)
cps.obs.exp.gcattc<-log(n3.codon.pair$gcattc/cps.gca.ttc)
cps.obs.exp.gcattg<-log(n3.codon.pair$gcattg/cps.gca.ttg)
cps.obs.exp.gcattt<-log(n3.codon.pair$gcattt/cps.gca.ttt)









cps.obs.exp.gccaaa<-log(n3.codon.pair$gccaaa/cps.gcc.aaa)
cps.obs.exp.gccaac<-log(n3.codon.pair$gccaac/cps.gcc.aac)
cps.obs.exp.gccaag<-log(n3.codon.pair$gccaag/cps.gcc.aag)
cps.obs.exp.gccaat<-log(n3.codon.pair$gccaat/cps.gcc.aat)
cps.obs.exp.gccaca<-log(n3.codon.pair$gccaca/cps.gcc.aca)
cps.obs.exp.gccacc<-log(n3.codon.pair$gccacc/cps.gcc.acc)
cps.obs.exp.gccacg<-log(n3.codon.pair$gccacg/cps.gcc.acg)
cps.obs.exp.gccact<-log(n3.codon.pair$gccact/cps.gcc.act)
cps.obs.exp.gccaga<-log(n3.codon.pair$gccaga/cps.gcc.aga)
cps.obs.exp.gccagc<-log(n3.codon.pair$gccagc/cps.gcc.agc)
cps.obs.exp.gccagg<-log(n3.codon.pair$gccagg/cps.gcc.agg)
cps.obs.exp.gccagt<-log(n3.codon.pair$gccagt/cps.gcc.agt)
cps.obs.exp.gccata<-log(n3.codon.pair$gccata/cps.gcc.ata)
cps.obs.exp.gccatc<-log(n3.codon.pair$gccatc/cps.gcc.atc)
cps.obs.exp.gccatg<-log(n3.codon.pair$gccatg/cps.gcc.atg)
cps.obs.exp.gccatt<-log(n3.codon.pair$gccatt/cps.gcc.att)

cps.obs.exp.gcccaa<-log(n3.codon.pair$gcccaa/cps.gcc.caa)
cps.obs.exp.gcccac<-log(n3.codon.pair$gcccac/cps.gcc.cac)
cps.obs.exp.gcccag<-log(n3.codon.pair$gcccag/cps.gcc.cag)
cps.obs.exp.gcccat<-log(n3.codon.pair$gcccat/cps.gcc.cat)
cps.obs.exp.gcccca<-log(n3.codon.pair$gcccca/cps.gcc.cca)
cps.obs.exp.gccccc<-log(n3.codon.pair$gccccc/cps.gcc.ccc)
cps.obs.exp.gccccg<-log(n3.codon.pair$gccccg/cps.gcc.ccg)
cps.obs.exp.gcccct<-log(n3.codon.pair$gcccct/cps.gcc.cct)
cps.obs.exp.gcccga<-log(n3.codon.pair$gcccga/cps.gcc.cga)
cps.obs.exp.gcccgc<-log(n3.codon.pair$gcccgc/cps.gcc.cgc)
cps.obs.exp.gcccgg<-log(n3.codon.pair$gcccgg/cps.gcc.cgg)
cps.obs.exp.gcccgt<-log(n3.codon.pair$gcccgt/cps.gcc.cgt)
cps.obs.exp.gcccta<-log(n3.codon.pair$gcccta/cps.gcc.cta)
cps.obs.exp.gccctc<-log(n3.codon.pair$gccctc/cps.gcc.ctc)
cps.obs.exp.gccctg<-log(n3.codon.pair$gccctg/cps.gcc.ctg)
cps.obs.exp.gccctt<-log(n3.codon.pair$gccctt/cps.gcc.ctt)

cps.obs.exp.gccgaa<-log(n3.codon.pair$gccgaa/cps.gcc.gaa)
cps.obs.exp.gccgac<-log(n3.codon.pair$gccgac/cps.gcc.gac)
cps.obs.exp.gccgag<-log(n3.codon.pair$gccgag/cps.gcc.gag)
cps.obs.exp.gccgat<-log(n3.codon.pair$gccgat/cps.gcc.gat)
cps.obs.exp.gccgca<-log(n3.codon.pair$gccgca/cps.gcc.gca)
cps.obs.exp.gccgcc<-log(n3.codon.pair$gccgcc/cps.gcc.gcc)
cps.obs.exp.gccgcg<-log(n3.codon.pair$gccgcg/cps.gcc.gcg)
cps.obs.exp.gccgct<-log(n3.codon.pair$gccgct/cps.gcc.gct)
cps.obs.exp.gccgga<-log(n3.codon.pair$gccgga/cps.gcc.gga)
cps.obs.exp.gccggc<-log(n3.codon.pair$gccggc/cps.gcc.ggc)
cps.obs.exp.gccggg<-log(n3.codon.pair$gccggg/cps.gcc.ggg)
cps.obs.exp.gccggt<-log(n3.codon.pair$gccggt/cps.gcc.ggt)
cps.obs.exp.gccgta<-log(n3.codon.pair$gccgta/cps.gcc.gta)
cps.obs.exp.gccgtc<-log(n3.codon.pair$gccgtc/cps.gcc.gtc)
cps.obs.exp.gccgtg<-log(n3.codon.pair$gccgtg/cps.gcc.gtg)
cps.obs.exp.gccgtt<-log(n3.codon.pair$gccgtt/cps.gcc.gtt)

#cps.obs.exp.gcctaa<-log(n3.codon.pair$gcctaa/cps.gcc.taa)
cps.obs.exp.gcctac<-log(n3.codon.pair$gcctac/cps.gcc.tac)
#cps.obs.exp.gcctag<-log(n3.codon.pair$gcctag/cps.gcc.tag)
cps.obs.exp.gcctat<-log(n3.codon.pair$gcctat/cps.gcc.tat)
cps.obs.exp.gcctca<-log(n3.codon.pair$gcctca/cps.gcc.tca)
cps.obs.exp.gcctcc<-log(n3.codon.pair$gcctcc/cps.gcc.tcc)
cps.obs.exp.gcctcg<-log(n3.codon.pair$gcctcg/cps.gcc.tcg)
cps.obs.exp.gcctct<-log(n3.codon.pair$gcctct/cps.gcc.tct)
#cps.obs.exp.gcctga<-log(n3.codon.pair$gcctga/cps.gcc.tga)
cps.obs.exp.gcctgc<-log(n3.codon.pair$gcctgc/cps.gcc.tgc)
cps.obs.exp.gcctgg<-log(n3.codon.pair$gcctgg/cps.gcc.tgg)
cps.obs.exp.gcctgt<-log(n3.codon.pair$gcctgt/cps.gcc.tgt)
cps.obs.exp.gcctta<-log(n3.codon.pair$gcctta/cps.gcc.tta)
cps.obs.exp.gccttc<-log(n3.codon.pair$gccttc/cps.gcc.ttc)
cps.obs.exp.gccttg<-log(n3.codon.pair$gccttg/cps.gcc.ttg)
cps.obs.exp.gccttt<-log(n3.codon.pair$gccttt/cps.gcc.ttt)









cps.obs.exp.gcgaaa<-log(n3.codon.pair$gcgaaa/cps.gcg.aaa)
cps.obs.exp.gcgaac<-log(n3.codon.pair$gcgaac/cps.gcg.aac)
cps.obs.exp.gcgaag<-log(n3.codon.pair$gcgaag/cps.gcg.aag)
cps.obs.exp.gcgaat<-log(n3.codon.pair$gcgaat/cps.gcg.aat)
cps.obs.exp.gcgaca<-log(n3.codon.pair$gcgaca/cps.gcg.aca)
cps.obs.exp.gcgacc<-log(n3.codon.pair$gcgacc/cps.gcg.acc)
cps.obs.exp.gcgacg<-log(n3.codon.pair$gcgacg/cps.gcg.acg)
cps.obs.exp.gcgact<-log(n3.codon.pair$gcgact/cps.gcg.act)
cps.obs.exp.gcgaga<-log(n3.codon.pair$gcgaga/cps.gcg.aga)
cps.obs.exp.gcgagc<-log(n3.codon.pair$gcgagc/cps.gcg.agc)
cps.obs.exp.gcgagg<-log(n3.codon.pair$gcgagg/cps.gcg.agg)
cps.obs.exp.gcgagt<-log(n3.codon.pair$gcgagt/cps.gcg.agt)
cps.obs.exp.gcgata<-log(n3.codon.pair$gcgata/cps.gcg.ata)
cps.obs.exp.gcgatc<-log(n3.codon.pair$gcgatc/cps.gcg.atc)
cps.obs.exp.gcgatg<-log(n3.codon.pair$gcgatg/cps.gcg.atg)
cps.obs.exp.gcgatt<-log(n3.codon.pair$gcgatt/cps.gcg.att)

cps.obs.exp.gcgcaa<-log(n3.codon.pair$gcgcaa/cps.gcg.caa)
cps.obs.exp.gcgcac<-log(n3.codon.pair$gcgcac/cps.gcg.cac)
cps.obs.exp.gcgcag<-log(n3.codon.pair$gcgcag/cps.gcg.cag)
cps.obs.exp.gcgcat<-log(n3.codon.pair$gcgcat/cps.gcg.cat)
cps.obs.exp.gcgcca<-log(n3.codon.pair$gcgcca/cps.gcg.cca)
cps.obs.exp.gcgccc<-log(n3.codon.pair$gcgccc/cps.gcg.ccc)
cps.obs.exp.gcgccg<-log(n3.codon.pair$gcgccg/cps.gcg.ccg)
cps.obs.exp.gcgcct<-log(n3.codon.pair$gcgcct/cps.gcg.cct)
cps.obs.exp.gcgcga<-log(n3.codon.pair$gcgcga/cps.gcg.cga)
cps.obs.exp.gcgcgc<-log(n3.codon.pair$gcgcgc/cps.gcg.cgc)
cps.obs.exp.gcgcgg<-log(n3.codon.pair$gcgcgg/cps.gcg.cgg)
cps.obs.exp.gcgcgt<-log(n3.codon.pair$gcgcgt/cps.gcg.cgt)
cps.obs.exp.gcgcta<-log(n3.codon.pair$gcgcta/cps.gcg.cta)
cps.obs.exp.gcgctc<-log(n3.codon.pair$gcgctc/cps.gcg.ctc)
cps.obs.exp.gcgctg<-log(n3.codon.pair$gcgctg/cps.gcg.ctg)
cps.obs.exp.gcgctt<-log(n3.codon.pair$gcgctt/cps.gcg.ctt)

cps.obs.exp.gcggaa<-log(n3.codon.pair$gcggaa/cps.gcg.gaa)
cps.obs.exp.gcggac<-log(n3.codon.pair$gcggac/cps.gcg.gac)
cps.obs.exp.gcggag<-log(n3.codon.pair$gcggag/cps.gcg.gag)
cps.obs.exp.gcggat<-log(n3.codon.pair$gcggat/cps.gcg.gat)
cps.obs.exp.gcggca<-log(n3.codon.pair$gcggca/cps.gcg.gca)
cps.obs.exp.gcggcc<-log(n3.codon.pair$gcggcc/cps.gcg.gcc)
cps.obs.exp.gcggcg<-log(n3.codon.pair$gcggcg/cps.gcg.gcg)
cps.obs.exp.gcggct<-log(n3.codon.pair$gcggct/cps.gcg.gct)
cps.obs.exp.gcggga<-log(n3.codon.pair$gcggga/cps.gcg.gga)
cps.obs.exp.gcgggc<-log(n3.codon.pair$gcgggc/cps.gcg.ggc)
cps.obs.exp.gcgggg<-log(n3.codon.pair$gcgggg/cps.gcg.ggg)
cps.obs.exp.gcgggt<-log(n3.codon.pair$gcgggt/cps.gcg.ggt)
cps.obs.exp.gcggta<-log(n3.codon.pair$gcggta/cps.gcg.gta)
cps.obs.exp.gcggtc<-log(n3.codon.pair$gcggtc/cps.gcg.gtc)
cps.obs.exp.gcggtg<-log(n3.codon.pair$gcggtg/cps.gcg.gtg)
cps.obs.exp.gcggtt<-log(n3.codon.pair$gcggtt/cps.gcg.gtt)

#cps.obs.exp.gcgtaa<-log(n3.codon.pair$gcgtaa/cps.gcg.taa)
cps.obs.exp.gcgtac<-log(n3.codon.pair$gcgtac/cps.gcg.tac)
#cps.obs.exp.gcgtag<-log(n3.codon.pair$gcgtag/cps.gcg.tag)
cps.obs.exp.gcgtat<-log(n3.codon.pair$gcgtat/cps.gcg.tat)
cps.obs.exp.gcgtca<-log(n3.codon.pair$gcgtca/cps.gcg.tca)
cps.obs.exp.gcgtcc<-log(n3.codon.pair$gcgtcc/cps.gcg.tcc)
cps.obs.exp.gcgtcg<-log(n3.codon.pair$gcgtcg/cps.gcg.tcg)
cps.obs.exp.gcgtct<-log(n3.codon.pair$gcgtct/cps.gcg.tct)
#cps.obs.exp.gcgtga<-log(n3.codon.pair$gcgtga/cps.gcg.tga)
cps.obs.exp.gcgtgc<-log(n3.codon.pair$gcgtgc/cps.gcg.tgc)
cps.obs.exp.gcgtgg<-log(n3.codon.pair$gcgtgg/cps.gcg.tgg)
cps.obs.exp.gcgtgt<-log(n3.codon.pair$gcgtgt/cps.gcg.tgt)
cps.obs.exp.gcgtta<-log(n3.codon.pair$gcgtta/cps.gcg.tta)
cps.obs.exp.gcgttc<-log(n3.codon.pair$gcgttc/cps.gcg.ttc)
cps.obs.exp.gcgttg<-log(n3.codon.pair$gcgttg/cps.gcg.ttg)
cps.obs.exp.gcgttt<-log(n3.codon.pair$gcgttt/cps.gcg.ttt)









cps.obs.exp.gctaaa<-log(n3.codon.pair$gctaaa/cps.gct.aaa)
cps.obs.exp.gctaac<-log(n3.codon.pair$gctaac/cps.gct.aac)
cps.obs.exp.gctaag<-log(n3.codon.pair$gctaag/cps.gct.aag)
cps.obs.exp.gctaat<-log(n3.codon.pair$gctaat/cps.gct.aat)
cps.obs.exp.gctaca<-log(n3.codon.pair$gctaca/cps.gct.aca)
cps.obs.exp.gctacc<-log(n3.codon.pair$gctacc/cps.gct.acc)
cps.obs.exp.gctacg<-log(n3.codon.pair$gctacg/cps.gct.acg)
cps.obs.exp.gctact<-log(n3.codon.pair$gctact/cps.gct.act)
cps.obs.exp.gctaga<-log(n3.codon.pair$gctaga/cps.gct.aga)
cps.obs.exp.gctagc<-log(n3.codon.pair$gctagc/cps.gct.agc)
cps.obs.exp.gctagg<-log(n3.codon.pair$gctagg/cps.gct.agg)
cps.obs.exp.gctagt<-log(n3.codon.pair$gctagt/cps.gct.agt)
cps.obs.exp.gctata<-log(n3.codon.pair$gctata/cps.gct.ata)
cps.obs.exp.gctatc<-log(n3.codon.pair$gctatc/cps.gct.atc)
cps.obs.exp.gctatg<-log(n3.codon.pair$gctatg/cps.gct.atg)
cps.obs.exp.gctatt<-log(n3.codon.pair$gctatt/cps.gct.att)

cps.obs.exp.gctcaa<-log(n3.codon.pair$gctcaa/cps.gct.caa)
cps.obs.exp.gctcac<-log(n3.codon.pair$gctcac/cps.gct.cac)
cps.obs.exp.gctcag<-log(n3.codon.pair$gctcag/cps.gct.cag)
cps.obs.exp.gctcat<-log(n3.codon.pair$gctcat/cps.gct.cat)
cps.obs.exp.gctcca<-log(n3.codon.pair$gctcca/cps.gct.cca)
cps.obs.exp.gctccc<-log(n3.codon.pair$gctccc/cps.gct.ccc)
cps.obs.exp.gctccg<-log(n3.codon.pair$gctccg/cps.gct.ccg)
cps.obs.exp.gctcct<-log(n3.codon.pair$gctcct/cps.gct.cct)
cps.obs.exp.gctcga<-log(n3.codon.pair$gctcga/cps.gct.cga)
cps.obs.exp.gctcgc<-log(n3.codon.pair$gctcgc/cps.gct.cgc)
cps.obs.exp.gctcgg<-log(n3.codon.pair$gctcgg/cps.gct.cgg)
cps.obs.exp.gctcgt<-log(n3.codon.pair$gctcgt/cps.gct.cgt)
cps.obs.exp.gctcta<-log(n3.codon.pair$gctcta/cps.gct.cta)
cps.obs.exp.gctctc<-log(n3.codon.pair$gctctc/cps.gct.ctc)
cps.obs.exp.gctctg<-log(n3.codon.pair$gctctg/cps.gct.ctg)
cps.obs.exp.gctctt<-log(n3.codon.pair$gctctt/cps.gct.ctt)

cps.obs.exp.gctgaa<-log(n3.codon.pair$gctgaa/cps.gct.gaa)
cps.obs.exp.gctgac<-log(n3.codon.pair$gctgac/cps.gct.gac)
cps.obs.exp.gctgag<-log(n3.codon.pair$gctgag/cps.gct.gag)
cps.obs.exp.gctgat<-log(n3.codon.pair$gctgat/cps.gct.gat)
cps.obs.exp.gctgca<-log(n3.codon.pair$gctgca/cps.gct.gca)
cps.obs.exp.gctgcc<-log(n3.codon.pair$gctgcc/cps.gct.gcc)
cps.obs.exp.gctgcg<-log(n3.codon.pair$gctgcg/cps.gct.gcg)
cps.obs.exp.gctgct<-log(n3.codon.pair$gctgct/cps.gct.gct)
cps.obs.exp.gctgga<-log(n3.codon.pair$gctgga/cps.gct.gga)
cps.obs.exp.gctggc<-log(n3.codon.pair$gctggc/cps.gct.ggc)
cps.obs.exp.gctggg<-log(n3.codon.pair$gctggg/cps.gct.ggg)
cps.obs.exp.gctggt<-log(n3.codon.pair$gctggt/cps.gct.ggt)
cps.obs.exp.gctgta<-log(n3.codon.pair$gctgta/cps.gct.gta)
cps.obs.exp.gctgtc<-log(n3.codon.pair$gctgtc/cps.gct.gtc)
cps.obs.exp.gctgtg<-log(n3.codon.pair$gctgtg/cps.gct.gtg)
cps.obs.exp.gctgtt<-log(n3.codon.pair$gctgtt/cps.gct.gtt)

#cps.obs.exp.gcttaa<-log(n3.codon.pair$gcttaa/cps.gct.taa)
cps.obs.exp.gcttac<-log(n3.codon.pair$gcttac/cps.gct.tac)
#cps.obs.exp.gcttag<-log(n3.codon.pair$gcttag/cps.gct.tag)
cps.obs.exp.gcttat<-log(n3.codon.pair$gcttat/cps.gct.tat)
cps.obs.exp.gcttca<-log(n3.codon.pair$gcttca/cps.gct.tca)
cps.obs.exp.gcttcc<-log(n3.codon.pair$gcttcc/cps.gct.tcc)
cps.obs.exp.gcttcg<-log(n3.codon.pair$gcttcg/cps.gct.tcg)
cps.obs.exp.gcttct<-log(n3.codon.pair$gcttct/cps.gct.tct)
#cps.obs.exp.gcttga<-log(n3.codon.pair$gcttga/cps.gct.tga)
cps.obs.exp.gcttgc<-log(n3.codon.pair$gcttgc/cps.gct.tgc)
cps.obs.exp.gcttgg<-log(n3.codon.pair$gcttgg/cps.gct.tgg)
cps.obs.exp.gcttgt<-log(n3.codon.pair$gcttgt/cps.gct.tgt)
cps.obs.exp.gcttta<-log(n3.codon.pair$gcttta/cps.gct.tta)
cps.obs.exp.gctttc<-log(n3.codon.pair$gctttc/cps.gct.ttc)
cps.obs.exp.gctttg<-log(n3.codon.pair$gctttg/cps.gct.ttg)
cps.obs.exp.gctttt<-log(n3.codon.pair$gctttt/cps.gct.ttt)



















cps.obs.exp.ggaaaa<-log(n3.codon.pair$ggaaaa/cps.gga.aaa)
cps.obs.exp.ggaaac<-log(n3.codon.pair$ggaaac/cps.gga.aac)
cps.obs.exp.ggaaag<-log(n3.codon.pair$ggaaag/cps.gga.aag)
cps.obs.exp.ggaaat<-log(n3.codon.pair$ggaaat/cps.gga.aat)
cps.obs.exp.ggaaca<-log(n3.codon.pair$ggaaca/cps.gga.aca)
cps.obs.exp.ggaacc<-log(n3.codon.pair$ggaacc/cps.gga.acc)
cps.obs.exp.ggaacg<-log(n3.codon.pair$ggaacg/cps.gga.acg)
cps.obs.exp.ggaact<-log(n3.codon.pair$ggaact/cps.gga.act)
cps.obs.exp.ggaaga<-log(n3.codon.pair$ggaaga/cps.gga.aga)
cps.obs.exp.ggaagc<-log(n3.codon.pair$ggaagc/cps.gga.agc)
cps.obs.exp.ggaagg<-log(n3.codon.pair$ggaagg/cps.gga.agg)
cps.obs.exp.ggaagt<-log(n3.codon.pair$ggaagt/cps.gga.agt)
cps.obs.exp.ggaata<-log(n3.codon.pair$ggaata/cps.gga.ata)
cps.obs.exp.ggaatc<-log(n3.codon.pair$ggaatc/cps.gga.atc)
cps.obs.exp.ggaatg<-log(n3.codon.pair$ggaatg/cps.gga.atg)
cps.obs.exp.ggaatt<-log(n3.codon.pair$ggaatt/cps.gga.att)

cps.obs.exp.ggacaa<-log(n3.codon.pair$ggacaa/cps.gga.caa)
cps.obs.exp.ggacac<-log(n3.codon.pair$ggacac/cps.gga.cac)
cps.obs.exp.ggacag<-log(n3.codon.pair$ggacag/cps.gga.cag)
cps.obs.exp.ggacat<-log(n3.codon.pair$ggacat/cps.gga.cat)
cps.obs.exp.ggacca<-log(n3.codon.pair$ggacca/cps.gga.cca)
cps.obs.exp.ggaccc<-log(n3.codon.pair$ggaccc/cps.gga.ccc)
cps.obs.exp.ggaccg<-log(n3.codon.pair$ggaccg/cps.gga.ccg)
cps.obs.exp.ggacct<-log(n3.codon.pair$ggacct/cps.gga.cct)
cps.obs.exp.ggacga<-log(n3.codon.pair$ggacga/cps.gga.cga)
cps.obs.exp.ggacgc<-log(n3.codon.pair$ggacgc/cps.gga.cgc)
cps.obs.exp.ggacgg<-log(n3.codon.pair$ggacgg/cps.gga.cgg)
cps.obs.exp.ggacgt<-log(n3.codon.pair$ggacgt/cps.gga.cgt)
cps.obs.exp.ggacta<-log(n3.codon.pair$ggacta/cps.gga.cta)
cps.obs.exp.ggactc<-log(n3.codon.pair$ggactc/cps.gga.ctc)
cps.obs.exp.ggactg<-log(n3.codon.pair$ggactg/cps.gga.ctg)
cps.obs.exp.ggactt<-log(n3.codon.pair$ggactt/cps.gga.ctt)

cps.obs.exp.ggagaa<-log(n3.codon.pair$ggagaa/cps.gga.gaa)
cps.obs.exp.ggagac<-log(n3.codon.pair$ggagac/cps.gga.gac)
cps.obs.exp.ggagag<-log(n3.codon.pair$ggagag/cps.gga.gag)
cps.obs.exp.ggagat<-log(n3.codon.pair$ggagat/cps.gga.gat)
cps.obs.exp.ggagca<-log(n3.codon.pair$ggagca/cps.gga.gca)
cps.obs.exp.ggagcc<-log(n3.codon.pair$ggagcc/cps.gga.gcc)
cps.obs.exp.ggagcg<-log(n3.codon.pair$ggagcg/cps.gga.gcg)
cps.obs.exp.ggagct<-log(n3.codon.pair$ggagct/cps.gga.gct)
cps.obs.exp.ggagga<-log(n3.codon.pair$ggagga/cps.gga.gga)
cps.obs.exp.ggaggc<-log(n3.codon.pair$ggaggc/cps.gga.ggc)
cps.obs.exp.ggaggg<-log(n3.codon.pair$ggaggg/cps.gga.ggg)
cps.obs.exp.ggaggt<-log(n3.codon.pair$ggaggt/cps.gga.ggt)
cps.obs.exp.ggagta<-log(n3.codon.pair$ggagta/cps.gga.gta)
cps.obs.exp.ggagtc<-log(n3.codon.pair$ggagtc/cps.gga.gtc)
cps.obs.exp.ggagtg<-log(n3.codon.pair$ggagtg/cps.gga.gtg)
cps.obs.exp.ggagtt<-log(n3.codon.pair$ggagtt/cps.gga.gtt)

#cps.obs.exp.ggataa<-log(n3.codon.pair$ggataa/cps.gga.taa)
cps.obs.exp.ggatac<-log(n3.codon.pair$ggatac/cps.gga.tac)
#cps.obs.exp.ggatag<-log(n3.codon.pair$ggatag/cps.gga.tag)
cps.obs.exp.ggatat<-log(n3.codon.pair$ggatat/cps.gga.tat)
cps.obs.exp.ggatca<-log(n3.codon.pair$ggatca/cps.gga.tca)
cps.obs.exp.ggatcc<-log(n3.codon.pair$ggatcc/cps.gga.tcc)
cps.obs.exp.ggatcg<-log(n3.codon.pair$ggatcg/cps.gga.tcg)
cps.obs.exp.ggatct<-log(n3.codon.pair$ggatct/cps.gga.tct)
#cps.obs.exp.ggatga<-log(n3.codon.pair$ggatga/cps.gga.tga)
cps.obs.exp.ggatgc<-log(n3.codon.pair$ggatgc/cps.gga.tgc)
cps.obs.exp.ggatgg<-log(n3.codon.pair$ggatgg/cps.gga.tgg)
cps.obs.exp.ggatgt<-log(n3.codon.pair$ggatgt/cps.gga.tgt)
cps.obs.exp.ggatta<-log(n3.codon.pair$ggatta/cps.gga.tta)
cps.obs.exp.ggattc<-log(n3.codon.pair$ggattc/cps.gga.ttc)
cps.obs.exp.ggattg<-log(n3.codon.pair$ggattg/cps.gga.ttg)
cps.obs.exp.ggattt<-log(n3.codon.pair$ggattt/cps.gga.ttt)










cps.obs.exp.ggcaaa<-log(n3.codon.pair$ggcaaa/cps.ggc.aaa)
cps.obs.exp.ggcaac<-log(n3.codon.pair$ggcaac/cps.ggc.aac)
cps.obs.exp.ggcaag<-log(n3.codon.pair$ggcaag/cps.ggc.aag)
cps.obs.exp.ggcaat<-log(n3.codon.pair$ggcaat/cps.ggc.aat)
cps.obs.exp.ggcaca<-log(n3.codon.pair$ggcaca/cps.ggc.aca)
cps.obs.exp.ggcacc<-log(n3.codon.pair$ggcacc/cps.ggc.acc)
cps.obs.exp.ggcacg<-log(n3.codon.pair$ggcacg/cps.ggc.acg)
cps.obs.exp.ggcact<-log(n3.codon.pair$ggcact/cps.ggc.act)
cps.obs.exp.ggcaga<-log(n3.codon.pair$ggcaga/cps.ggc.aga)
cps.obs.exp.ggcagc<-log(n3.codon.pair$ggcagc/cps.ggc.agc)
cps.obs.exp.ggcagg<-log(n3.codon.pair$ggcagg/cps.ggc.agg)
cps.obs.exp.ggcagt<-log(n3.codon.pair$ggcagt/cps.ggc.agt)
cps.obs.exp.ggcata<-log(n3.codon.pair$ggcata/cps.ggc.ata)
cps.obs.exp.ggcatc<-log(n3.codon.pair$ggcatc/cps.ggc.atc)
cps.obs.exp.ggcatg<-log(n3.codon.pair$ggcatg/cps.ggc.atg)
cps.obs.exp.ggcatt<-log(n3.codon.pair$ggcatt/cps.ggc.att)

cps.obs.exp.ggccaa<-log(n3.codon.pair$ggccaa/cps.ggc.caa)
cps.obs.exp.ggccac<-log(n3.codon.pair$ggccac/cps.ggc.cac)
cps.obs.exp.ggccag<-log(n3.codon.pair$ggccag/cps.ggc.cag)
cps.obs.exp.ggccat<-log(n3.codon.pair$ggccat/cps.ggc.cat)
cps.obs.exp.ggccca<-log(n3.codon.pair$ggccca/cps.ggc.cca)
cps.obs.exp.ggcccc<-log(n3.codon.pair$ggcccc/cps.ggc.ccc)
cps.obs.exp.ggcccg<-log(n3.codon.pair$ggcccg/cps.ggc.ccg)
cps.obs.exp.ggccct<-log(n3.codon.pair$ggccct/cps.ggc.cct)
cps.obs.exp.ggccga<-log(n3.codon.pair$ggccga/cps.ggc.cga)
cps.obs.exp.ggccgc<-log(n3.codon.pair$ggccgc/cps.ggc.cgc)
cps.obs.exp.ggccgg<-log(n3.codon.pair$ggccgg/cps.ggc.cgg)
cps.obs.exp.ggccgt<-log(n3.codon.pair$ggccgt/cps.ggc.cgt)
cps.obs.exp.ggccta<-log(n3.codon.pair$ggccta/cps.ggc.cta)
cps.obs.exp.ggcctc<-log(n3.codon.pair$ggcctc/cps.ggc.ctc)
cps.obs.exp.ggcctg<-log(n3.codon.pair$ggcctg/cps.ggc.ctg)
cps.obs.exp.ggcctt<-log(n3.codon.pair$ggcctt/cps.ggc.ctt)

cps.obs.exp.ggcgaa<-log(n3.codon.pair$ggcgaa/cps.ggc.gaa)
cps.obs.exp.ggcgac<-log(n3.codon.pair$ggcgac/cps.ggc.gac)
cps.obs.exp.ggcgag<-log(n3.codon.pair$ggcgag/cps.ggc.gag)
cps.obs.exp.ggcgat<-log(n3.codon.pair$ggcgat/cps.ggc.gat)
cps.obs.exp.ggcgca<-log(n3.codon.pair$ggcgca/cps.ggc.gca)
cps.obs.exp.ggcgcc<-log(n3.codon.pair$ggcgcc/cps.ggc.gcc)
cps.obs.exp.ggcgcg<-log(n3.codon.pair$ggcgcg/cps.ggc.gcg)
cps.obs.exp.ggcgct<-log(n3.codon.pair$ggcgct/cps.ggc.gct)
cps.obs.exp.ggcgga<-log(n3.codon.pair$ggcgga/cps.ggc.gga)
cps.obs.exp.ggcggc<-log(n3.codon.pair$ggcggc/cps.ggc.ggc)
cps.obs.exp.ggcggg<-log(n3.codon.pair$ggcggg/cps.ggc.ggg)
cps.obs.exp.ggcggt<-log(n3.codon.pair$ggcggt/cps.ggc.ggt)
cps.obs.exp.ggcgta<-log(n3.codon.pair$ggcgta/cps.ggc.gta)
cps.obs.exp.ggcgtc<-log(n3.codon.pair$ggcgtc/cps.ggc.gtc)
cps.obs.exp.ggcgtg<-log(n3.codon.pair$ggcgtg/cps.ggc.gtg)
cps.obs.exp.ggcgtt<-log(n3.codon.pair$ggcgtt/cps.ggc.gtt)

#cps.obs.exp.ggctaa<-log(n3.codon.pair$ggctaa/cps.ggc.taa)
cps.obs.exp.ggctac<-log(n3.codon.pair$ggctac/cps.ggc.tac)
#cps.obs.exp.ggctag<-log(n3.codon.pair$ggctag/cps.ggc.tag)
cps.obs.exp.ggctat<-log(n3.codon.pair$ggctat/cps.ggc.tat)
cps.obs.exp.ggctca<-log(n3.codon.pair$ggctca/cps.ggc.tca)
cps.obs.exp.ggctcc<-log(n3.codon.pair$ggctcc/cps.ggc.tcc)
cps.obs.exp.ggctcg<-log(n3.codon.pair$ggctcg/cps.ggc.tcg)
cps.obs.exp.ggctct<-log(n3.codon.pair$ggctct/cps.ggc.tct)
#cps.obs.exp.ggctga<-log(n3.codon.pair$ggctga/cps.ggc.tga)
cps.obs.exp.ggctgc<-log(n3.codon.pair$ggctgc/cps.ggc.tgc)
cps.obs.exp.ggctgg<-log(n3.codon.pair$ggctgg/cps.ggc.tgg)
cps.obs.exp.ggctgt<-log(n3.codon.pair$ggctgt/cps.ggc.tgt)
cps.obs.exp.ggctta<-log(n3.codon.pair$ggctta/cps.ggc.tta)
cps.obs.exp.ggcttc<-log(n3.codon.pair$ggcttc/cps.ggc.ttc)
cps.obs.exp.ggcttg<-log(n3.codon.pair$ggcttg/cps.ggc.ttg)
cps.obs.exp.ggcttt<-log(n3.codon.pair$ggcttt/cps.ggc.ttt)









cps.obs.exp.gggaaa<-log(n3.codon.pair$gggaaa/cps.ggg.aaa)
cps.obs.exp.gggaac<-log(n3.codon.pair$gggaac/cps.ggg.aac)
cps.obs.exp.gggaag<-log(n3.codon.pair$gggaag/cps.ggg.aag)
cps.obs.exp.gggaat<-log(n3.codon.pair$gggaat/cps.ggg.aat)
cps.obs.exp.gggaca<-log(n3.codon.pair$gggaca/cps.ggg.aca)
cps.obs.exp.gggacc<-log(n3.codon.pair$gggacc/cps.ggg.acc)
cps.obs.exp.gggacg<-log(n3.codon.pair$gggacg/cps.ggg.acg)
cps.obs.exp.gggact<-log(n3.codon.pair$gggact/cps.ggg.act)
cps.obs.exp.gggaga<-log(n3.codon.pair$gggaga/cps.ggg.aga)
cps.obs.exp.gggagc<-log(n3.codon.pair$gggagc/cps.ggg.agc)
cps.obs.exp.gggagg<-log(n3.codon.pair$gggagg/cps.ggg.agg)
cps.obs.exp.gggagt<-log(n3.codon.pair$gggagt/cps.ggg.agt)
cps.obs.exp.gggata<-log(n3.codon.pair$gggata/cps.ggg.ata)
cps.obs.exp.gggatc<-log(n3.codon.pair$gggatc/cps.ggg.atc)
cps.obs.exp.gggatg<-log(n3.codon.pair$gggatg/cps.ggg.atg)
cps.obs.exp.gggatt<-log(n3.codon.pair$gggatt/cps.ggg.att)

cps.obs.exp.gggcaa<-log(n3.codon.pair$gggcaa/cps.ggg.caa)
cps.obs.exp.gggcac<-log(n3.codon.pair$gggcac/cps.ggg.cac)
cps.obs.exp.gggcag<-log(n3.codon.pair$gggcag/cps.ggg.cag)
cps.obs.exp.gggcat<-log(n3.codon.pair$gggcat/cps.ggg.cat)
cps.obs.exp.gggcca<-log(n3.codon.pair$gggcca/cps.ggg.cca)
cps.obs.exp.gggccc<-log(n3.codon.pair$gggccc/cps.ggg.ccc)
cps.obs.exp.gggccg<-log(n3.codon.pair$gggccg/cps.ggg.ccg)
cps.obs.exp.gggcct<-log(n3.codon.pair$gggcct/cps.ggg.cct)
cps.obs.exp.gggcga<-log(n3.codon.pair$gggcga/cps.ggg.cga)
cps.obs.exp.gggcgc<-log(n3.codon.pair$gggcgc/cps.ggg.cgc)
cps.obs.exp.gggcgg<-log(n3.codon.pair$gggcgg/cps.ggg.cgg)
cps.obs.exp.gggcgt<-log(n3.codon.pair$gggcgt/cps.ggg.cgt)
cps.obs.exp.gggcta<-log(n3.codon.pair$gggcta/cps.ggg.cta)
cps.obs.exp.gggctc<-log(n3.codon.pair$gggctc/cps.ggg.ctc)
cps.obs.exp.gggctg<-log(n3.codon.pair$gggctg/cps.ggg.ctg)
cps.obs.exp.gggctt<-log(n3.codon.pair$gggctt/cps.ggg.ctt)

cps.obs.exp.ggggaa<-log(n3.codon.pair$ggggaa/cps.ggg.gaa)
cps.obs.exp.ggggac<-log(n3.codon.pair$ggggac/cps.ggg.gac)
cps.obs.exp.ggggag<-log(n3.codon.pair$ggggag/cps.ggg.gag)
cps.obs.exp.ggggat<-log(n3.codon.pair$ggggat/cps.ggg.gat)
cps.obs.exp.ggggca<-log(n3.codon.pair$ggggca/cps.ggg.gca)
cps.obs.exp.ggggcc<-log(n3.codon.pair$ggggcc/cps.ggg.gcc)
cps.obs.exp.ggggcg<-log(n3.codon.pair$ggggcg/cps.ggg.gcg)
cps.obs.exp.ggggct<-log(n3.codon.pair$ggggct/cps.ggg.gct)
cps.obs.exp.ggggga<-log(n3.codon.pair$ggggga/cps.ggg.gga)
cps.obs.exp.gggggc<-log(n3.codon.pair$gggggc/cps.ggg.ggc)
cps.obs.exp.gggggg<-log(n3.codon.pair$gggggg/cps.ggg.ggg)
cps.obs.exp.gggggt<-log(n3.codon.pair$gggggt/cps.ggg.ggt)
cps.obs.exp.ggggta<-log(n3.codon.pair$ggggta/cps.ggg.gta)
cps.obs.exp.ggggtc<-log(n3.codon.pair$ggggtc/cps.ggg.gtc)
cps.obs.exp.ggggtg<-log(n3.codon.pair$ggggtg/cps.ggg.gtg)
cps.obs.exp.ggggtt<-log(n3.codon.pair$ggggtt/cps.ggg.gtt)

#cps.obs.exp.gggtaa<-log(n3.codon.pair$gggtaa/cps.ggg.taa)
cps.obs.exp.gggtac<-log(n3.codon.pair$gggtac/cps.ggg.tac)
#cps.obs.exp.gggtag<-log(n3.codon.pair$gggtag/cps.ggg.tag)
cps.obs.exp.gggtat<-log(n3.codon.pair$gggtat/cps.ggg.tat)
cps.obs.exp.gggtca<-log(n3.codon.pair$gggtca/cps.ggg.tca)
cps.obs.exp.gggtcc<-log(n3.codon.pair$gggtcc/cps.ggg.tcc)
cps.obs.exp.gggtcg<-log(n3.codon.pair$gggtcg/cps.ggg.tcg)
cps.obs.exp.gggtct<-log(n3.codon.pair$gggtct/cps.ggg.tct)
#cps.obs.exp.gggtga<-log(n3.codon.pair$gggtga/cps.ggg.tga)
cps.obs.exp.gggtgc<-log(n3.codon.pair$gggtgc/cps.ggg.tgc)
cps.obs.exp.gggtgg<-log(n3.codon.pair$gggtgg/cps.ggg.tgg)
cps.obs.exp.gggtgt<-log(n3.codon.pair$gggtgt/cps.ggg.tgt)
cps.obs.exp.gggtta<-log(n3.codon.pair$gggtta/cps.ggg.tta)
cps.obs.exp.gggttc<-log(n3.codon.pair$gggttc/cps.ggg.ttc)
cps.obs.exp.gggttg<-log(n3.codon.pair$gggttg/cps.ggg.ttg)
cps.obs.exp.gggttt<-log(n3.codon.pair$gggttt/cps.ggg.ttt)











cps.obs.exp.ggtaaa<-log(n3.codon.pair$ggtaaa/cps.ggt.aaa)
cps.obs.exp.ggtaac<-log(n3.codon.pair$ggtaac/cps.ggt.aac)
cps.obs.exp.ggtaag<-log(n3.codon.pair$ggtaag/cps.ggt.aag)
cps.obs.exp.ggtaat<-log(n3.codon.pair$ggtaat/cps.ggt.aat)
cps.obs.exp.ggtaca<-log(n3.codon.pair$ggtaca/cps.ggt.aca)
cps.obs.exp.ggtacc<-log(n3.codon.pair$ggtacc/cps.ggt.acc)
cps.obs.exp.ggtacg<-log(n3.codon.pair$ggtacg/cps.ggt.acg)
cps.obs.exp.ggtact<-log(n3.codon.pair$ggtact/cps.ggt.act)
cps.obs.exp.ggtaga<-log(n3.codon.pair$ggtaga/cps.ggt.aga)
cps.obs.exp.ggtagc<-log(n3.codon.pair$ggtagc/cps.ggt.agc)
cps.obs.exp.ggtagg<-log(n3.codon.pair$ggtagg/cps.ggt.agg)
cps.obs.exp.ggtagt<-log(n3.codon.pair$ggtagt/cps.ggt.agt)
cps.obs.exp.ggtata<-log(n3.codon.pair$ggtata/cps.ggt.ata)
cps.obs.exp.ggtatc<-log(n3.codon.pair$ggtatc/cps.ggt.atc)
cps.obs.exp.ggtatg<-log(n3.codon.pair$ggtatg/cps.ggt.atg)
cps.obs.exp.ggtatt<-log(n3.codon.pair$ggtatt/cps.ggt.att)

cps.obs.exp.ggtcaa<-log(n3.codon.pair$ggtcaa/cps.ggt.caa)
cps.obs.exp.ggtcac<-log(n3.codon.pair$ggtcac/cps.ggt.cac)
cps.obs.exp.ggtcag<-log(n3.codon.pair$ggtcag/cps.ggt.cag)
cps.obs.exp.ggtcat<-log(n3.codon.pair$ggtcat/cps.ggt.cat)
cps.obs.exp.ggtcca<-log(n3.codon.pair$ggtcca/cps.ggt.cca)
cps.obs.exp.ggtccc<-log(n3.codon.pair$ggtccc/cps.ggt.ccc)
cps.obs.exp.ggtccg<-log(n3.codon.pair$ggtccg/cps.ggt.ccg)
cps.obs.exp.ggtcct<-log(n3.codon.pair$ggtcct/cps.ggt.cct)
cps.obs.exp.ggtcga<-log(n3.codon.pair$ggtcga/cps.ggt.cga)
cps.obs.exp.ggtcgc<-log(n3.codon.pair$ggtcgc/cps.ggt.cgc)
cps.obs.exp.ggtcgg<-log(n3.codon.pair$ggtcgg/cps.ggt.cgg)
cps.obs.exp.ggtcgt<-log(n3.codon.pair$ggtcgt/cps.ggt.cgt)
cps.obs.exp.ggtcta<-log(n3.codon.pair$ggtcta/cps.ggt.cta)
cps.obs.exp.ggtctc<-log(n3.codon.pair$ggtctc/cps.ggt.ctc)
cps.obs.exp.ggtctg<-log(n3.codon.pair$ggtctg/cps.ggt.ctg)
cps.obs.exp.ggtctt<-log(n3.codon.pair$ggtctt/cps.ggt.ctt)

cps.obs.exp.ggtgaa<-log(n3.codon.pair$ggtgaa/cps.ggt.gaa)
cps.obs.exp.ggtgac<-log(n3.codon.pair$ggtgac/cps.ggt.gac)
cps.obs.exp.ggtgag<-log(n3.codon.pair$ggtgag/cps.ggt.gag)
cps.obs.exp.ggtgat<-log(n3.codon.pair$ggtgat/cps.ggt.gat)
cps.obs.exp.ggtgca<-log(n3.codon.pair$ggtgca/cps.ggt.gca)
cps.obs.exp.ggtgcc<-log(n3.codon.pair$ggtgcc/cps.ggt.gcc)
cps.obs.exp.ggtgcg<-log(n3.codon.pair$ggtgcg/cps.ggt.gcg)
cps.obs.exp.ggtgct<-log(n3.codon.pair$ggtgct/cps.ggt.gct)
cps.obs.exp.ggtgga<-log(n3.codon.pair$ggtgga/cps.ggt.gga)
cps.obs.exp.ggtggc<-log(n3.codon.pair$ggtggc/cps.ggt.ggc)
cps.obs.exp.ggtggg<-log(n3.codon.pair$ggtggg/cps.ggt.ggg)
cps.obs.exp.ggtggt<-log(n3.codon.pair$ggtggt/cps.ggt.ggt)
cps.obs.exp.ggtgta<-log(n3.codon.pair$ggtgta/cps.ggt.gta)
cps.obs.exp.ggtgtc<-log(n3.codon.pair$ggtgtc/cps.ggt.gtc)
cps.obs.exp.ggtgtg<-log(n3.codon.pair$ggtgtg/cps.ggt.gtg)
cps.obs.exp.ggtgtt<-log(n3.codon.pair$ggtgtt/cps.ggt.gtt)

#cps.obs.exp.ggttaa<-log(n3.codon.pair$ggttaa/cps.ggt.taa)
cps.obs.exp.ggttac<-log(n3.codon.pair$ggttac/cps.ggt.tac)
#cps.obs.exp.ggttag<-log(n3.codon.pair$ggttag/cps.ggt.tag)
cps.obs.exp.ggttat<-log(n3.codon.pair$ggttat/cps.ggt.tat)
cps.obs.exp.ggttca<-log(n3.codon.pair$ggttca/cps.ggt.tca)
cps.obs.exp.ggttcc<-log(n3.codon.pair$ggttcc/cps.ggt.tcc)
cps.obs.exp.ggttcg<-log(n3.codon.pair$ggttcg/cps.ggt.tcg)
cps.obs.exp.ggttct<-log(n3.codon.pair$ggttct/cps.ggt.tct)
#cps.obs.exp.ggttga<-log(n3.codon.pair$ggttga/cps.ggt.tga)
cps.obs.exp.ggttgc<-log(n3.codon.pair$ggttgc/cps.ggt.tgc)
cps.obs.exp.ggttgg<-log(n3.codon.pair$ggttgg/cps.ggt.tgg)
cps.obs.exp.ggttgt<-log(n3.codon.pair$ggttgt/cps.ggt.tgt)
cps.obs.exp.ggttta<-log(n3.codon.pair$ggttta/cps.ggt.tta)
cps.obs.exp.ggtttc<-log(n3.codon.pair$ggtttc/cps.ggt.ttc)
cps.obs.exp.ggtttg<-log(n3.codon.pair$ggtttg/cps.ggt.ttg)
cps.obs.exp.ggtttt<-log(n3.codon.pair$ggtttt/cps.ggt.ttt)




















cps.obs.exp.gtaaaa<-log(n3.codon.pair$gtaaaa/cps.gta.aaa)
cps.obs.exp.gtaaac<-log(n3.codon.pair$gtaaac/cps.gta.aac)
cps.obs.exp.gtaaag<-log(n3.codon.pair$gtaaag/cps.gta.aag)
cps.obs.exp.gtaaat<-log(n3.codon.pair$gtaaat/cps.gta.aat)
cps.obs.exp.gtaaca<-log(n3.codon.pair$gtaaca/cps.gta.aca)
cps.obs.exp.gtaacc<-log(n3.codon.pair$gtaacc/cps.gta.acc)
cps.obs.exp.gtaacg<-log(n3.codon.pair$gtaacg/cps.gta.acg)
cps.obs.exp.gtaact<-log(n3.codon.pair$gtaact/cps.gta.act)
cps.obs.exp.gtaaga<-log(n3.codon.pair$gtaaga/cps.gta.aga)
cps.obs.exp.gtaagc<-log(n3.codon.pair$gtaagc/cps.gta.agc)
cps.obs.exp.gtaagg<-log(n3.codon.pair$gtaagg/cps.gta.agg)
cps.obs.exp.gtaagt<-log(n3.codon.pair$gtaagt/cps.gta.agt)
cps.obs.exp.gtaata<-log(n3.codon.pair$gtaata/cps.gta.ata)
cps.obs.exp.gtaatc<-log(n3.codon.pair$gtaatc/cps.gta.atc)
cps.obs.exp.gtaatg<-log(n3.codon.pair$gtaatg/cps.gta.atg)
cps.obs.exp.gtaatt<-log(n3.codon.pair$gtaatt/cps.gta.att)

cps.obs.exp.gtacaa<-log(n3.codon.pair$gtacaa/cps.gta.caa)
cps.obs.exp.gtacac<-log(n3.codon.pair$gtacac/cps.gta.cac)
cps.obs.exp.gtacag<-log(n3.codon.pair$gtacag/cps.gta.cag)
cps.obs.exp.gtacat<-log(n3.codon.pair$gtacat/cps.gta.cat)
cps.obs.exp.gtacca<-log(n3.codon.pair$gtacca/cps.gta.cca)
cps.obs.exp.gtaccc<-log(n3.codon.pair$gtaccc/cps.gta.ccc)
cps.obs.exp.gtaccg<-log(n3.codon.pair$gtaccg/cps.gta.ccg)
cps.obs.exp.gtacct<-log(n3.codon.pair$gtacct/cps.gta.cct)
cps.obs.exp.gtacga<-log(n3.codon.pair$gtacga/cps.gta.cga)
cps.obs.exp.gtacgc<-log(n3.codon.pair$gtacgc/cps.gta.cgc)
cps.obs.exp.gtacgg<-log(n3.codon.pair$gtacgg/cps.gta.cgg)
cps.obs.exp.gtacgt<-log(n3.codon.pair$gtacgt/cps.gta.cgt)
cps.obs.exp.gtacta<-log(n3.codon.pair$gtacta/cps.gta.cta)
cps.obs.exp.gtactc<-log(n3.codon.pair$gtactc/cps.gta.ctc)
cps.obs.exp.gtactg<-log(n3.codon.pair$gtactg/cps.gta.ctg)
cps.obs.exp.gtactt<-log(n3.codon.pair$gtactt/cps.gta.ctt)

cps.obs.exp.gtagaa<-log(n3.codon.pair$gtagaa/cps.gta.gaa)
cps.obs.exp.gtagac<-log(n3.codon.pair$gtagac/cps.gta.gac)
cps.obs.exp.gtagag<-log(n3.codon.pair$gtagag/cps.gta.gag)
cps.obs.exp.gtagat<-log(n3.codon.pair$gtagat/cps.gta.gat)
cps.obs.exp.gtagca<-log(n3.codon.pair$gtagca/cps.gta.gca)
cps.obs.exp.gtagcc<-log(n3.codon.pair$gtagcc/cps.gta.gcc)
cps.obs.exp.gtagcg<-log(n3.codon.pair$gtagcg/cps.gta.gcg)
cps.obs.exp.gtagct<-log(n3.codon.pair$gtagct/cps.gta.gct)
cps.obs.exp.gtagga<-log(n3.codon.pair$gtagga/cps.gta.gga)
cps.obs.exp.gtaggc<-log(n3.codon.pair$gtaggc/cps.gta.ggc)
cps.obs.exp.gtaggg<-log(n3.codon.pair$gtaggg/cps.gta.ggg)
cps.obs.exp.gtaggt<-log(n3.codon.pair$gtaggt/cps.gta.ggt)
cps.obs.exp.gtagta<-log(n3.codon.pair$gtagta/cps.gta.gta)
cps.obs.exp.gtagtc<-log(n3.codon.pair$gtagtc/cps.gta.gtc)
cps.obs.exp.gtagtg<-log(n3.codon.pair$gtagtg/cps.gta.gtg)
cps.obs.exp.gtagtt<-log(n3.codon.pair$gtagtt/cps.gta.gtt)

#cps.obs.exp.gtataa<-log(n3.codon.pair$gtataa/cps.gta.taa)
cps.obs.exp.gtatac<-log(n3.codon.pair$gtatac/cps.gta.tac)
#cps.obs.exp.gtatag<-log(n3.codon.pair$gtatag/cps.gta.tag)
cps.obs.exp.gtatat<-log(n3.codon.pair$gtatat/cps.gta.tat)
cps.obs.exp.gtatca<-log(n3.codon.pair$gtatca/cps.gta.tca)
cps.obs.exp.gtatcc<-log(n3.codon.pair$gtatcc/cps.gta.tcc)
cps.obs.exp.gtatcg<-log(n3.codon.pair$gtatcg/cps.gta.tcg)
cps.obs.exp.gtatct<-log(n3.codon.pair$gtatct/cps.gta.tct)
#cps.obs.exp.gtatga<-log(n3.codon.pair$gtatga/cps.gta.tga)
cps.obs.exp.gtatgc<-log(n3.codon.pair$gtatgc/cps.gta.tgc)
cps.obs.exp.gtatgg<-log(n3.codon.pair$gtatgg/cps.gta.tgg)
cps.obs.exp.gtatgt<-log(n3.codon.pair$gtatgt/cps.gta.tgt)
cps.obs.exp.gtatta<-log(n3.codon.pair$gtatta/cps.gta.tta)
cps.obs.exp.gtattc<-log(n3.codon.pair$gtattc/cps.gta.ttc)
cps.obs.exp.gtattg<-log(n3.codon.pair$gtattg/cps.gta.ttg)
cps.obs.exp.gtattt<-log(n3.codon.pair$gtattt/cps.gta.ttt)









cps.obs.exp.gtcaaa<-log(n3.codon.pair$gtcaaa/cps.gtc.aaa)
cps.obs.exp.gtcaac<-log(n3.codon.pair$gtcaac/cps.gtc.aac)
cps.obs.exp.gtcaag<-log(n3.codon.pair$gtcaag/cps.gtc.aag)
cps.obs.exp.gtcaat<-log(n3.codon.pair$gtcaat/cps.gtc.aat)
cps.obs.exp.gtcaca<-log(n3.codon.pair$gtcaca/cps.gtc.aca)
cps.obs.exp.gtcacc<-log(n3.codon.pair$gtcacc/cps.gtc.acc)
cps.obs.exp.gtcacg<-log(n3.codon.pair$gtcacg/cps.gtc.acg)
cps.obs.exp.gtcact<-log(n3.codon.pair$gtcact/cps.gtc.act)
cps.obs.exp.gtcaga<-log(n3.codon.pair$gtcaga/cps.gtc.aga)
cps.obs.exp.gtcagc<-log(n3.codon.pair$gtcagc/cps.gtc.agc)
cps.obs.exp.gtcagg<-log(n3.codon.pair$gtcagg/cps.gtc.agg)
cps.obs.exp.gtcagt<-log(n3.codon.pair$gtcagt/cps.gtc.agt)
cps.obs.exp.gtcata<-log(n3.codon.pair$gtcata/cps.gtc.ata)
cps.obs.exp.gtcatc<-log(n3.codon.pair$gtcatc/cps.gtc.atc)
cps.obs.exp.gtcatg<-log(n3.codon.pair$gtcatg/cps.gtc.atg)
cps.obs.exp.gtcatt<-log(n3.codon.pair$gtcatt/cps.gtc.att)

cps.obs.exp.gtccaa<-log(n3.codon.pair$gtccaa/cps.gtc.caa)
cps.obs.exp.gtccac<-log(n3.codon.pair$gtccac/cps.gtc.cac)
cps.obs.exp.gtccag<-log(n3.codon.pair$gtccag/cps.gtc.cag)
cps.obs.exp.gtccat<-log(n3.codon.pair$gtccat/cps.gtc.cat)
cps.obs.exp.gtccca<-log(n3.codon.pair$gtccca/cps.gtc.cca)
cps.obs.exp.gtcccc<-log(n3.codon.pair$gtcccc/cps.gtc.ccc)
cps.obs.exp.gtcccg<-log(n3.codon.pair$gtcccg/cps.gtc.ccg)
cps.obs.exp.gtccct<-log(n3.codon.pair$gtccct/cps.gtc.cct)
cps.obs.exp.gtccga<-log(n3.codon.pair$gtccga/cps.gtc.cga)
cps.obs.exp.gtccgc<-log(n3.codon.pair$gtccgc/cps.gtc.cgc)
cps.obs.exp.gtccgg<-log(n3.codon.pair$gtccgg/cps.gtc.cgg)
cps.obs.exp.gtccgt<-log(n3.codon.pair$gtccgt/cps.gtc.cgt)
cps.obs.exp.gtccta<-log(n3.codon.pair$gtccta/cps.gtc.cta)
cps.obs.exp.gtcctc<-log(n3.codon.pair$gtcctc/cps.gtc.ctc)
cps.obs.exp.gtcctg<-log(n3.codon.pair$gtcctg/cps.gtc.ctg)
cps.obs.exp.gtcctt<-log(n3.codon.pair$gtcctt/cps.gtc.ctt)

cps.obs.exp.gtcgaa<-log(n3.codon.pair$gtcgaa/cps.gtc.gaa)
cps.obs.exp.gtcgac<-log(n3.codon.pair$gtcgac/cps.gtc.gac)
cps.obs.exp.gtcgag<-log(n3.codon.pair$gtcgag/cps.gtc.gag)
cps.obs.exp.gtcgat<-log(n3.codon.pair$gtcgat/cps.gtc.gat)
cps.obs.exp.gtcgca<-log(n3.codon.pair$gtcgca/cps.gtc.gca)
cps.obs.exp.gtcgcc<-log(n3.codon.pair$gtcgcc/cps.gtc.gcc)
cps.obs.exp.gtcgcg<-log(n3.codon.pair$gtcgcg/cps.gtc.gcg)
cps.obs.exp.gtcgct<-log(n3.codon.pair$gtcgct/cps.gtc.gct)
cps.obs.exp.gtcgga<-log(n3.codon.pair$gtcgga/cps.gtc.gga)
cps.obs.exp.gtcggc<-log(n3.codon.pair$gtcggc/cps.gtc.ggc)
cps.obs.exp.gtcggg<-log(n3.codon.pair$gtcggg/cps.gtc.ggg)
cps.obs.exp.gtcggt<-log(n3.codon.pair$gtcggt/cps.gtc.ggt)
cps.obs.exp.gtcgta<-log(n3.codon.pair$gtcgta/cps.gtc.gta)
cps.obs.exp.gtcgtc<-log(n3.codon.pair$gtcgtc/cps.gtc.gtc)
cps.obs.exp.gtcgtg<-log(n3.codon.pair$gtcgtg/cps.gtc.gtg)
cps.obs.exp.gtcgtt<-log(n3.codon.pair$gtcgtt/cps.gtc.gtt)

#cps.obs.exp.gtctaa<-log(n3.codon.pair$gtctaa/cps.gtc.taa)
cps.obs.exp.gtctac<-log(n3.codon.pair$gtctac/cps.gtc.tac)
#cps.obs.exp.gtctag<-log(n3.codon.pair$gtctag/cps.gtc.tag)
cps.obs.exp.gtctat<-log(n3.codon.pair$gtctat/cps.gtc.tat)
cps.obs.exp.gtctca<-log(n3.codon.pair$gtctca/cps.gtc.tca)
cps.obs.exp.gtctcc<-log(n3.codon.pair$gtctcc/cps.gtc.tcc)
cps.obs.exp.gtctcg<-log(n3.codon.pair$gtctcg/cps.gtc.tcg)
cps.obs.exp.gtctct<-log(n3.codon.pair$gtctct/cps.gtc.tct)
#cps.obs.exp.gtctga<-log(n3.codon.pair$gtctga/cps.gtc.tga)
cps.obs.exp.gtctgc<-log(n3.codon.pair$gtctgc/cps.gtc.tgc)
cps.obs.exp.gtctgg<-log(n3.codon.pair$gtctgg/cps.gtc.tgg)
cps.obs.exp.gtctgt<-log(n3.codon.pair$gtctgt/cps.gtc.tgt)
cps.obs.exp.gtctta<-log(n3.codon.pair$gtctta/cps.gtc.tta)
cps.obs.exp.gtcttc<-log(n3.codon.pair$gtcttc/cps.gtc.ttc)
cps.obs.exp.gtcttg<-log(n3.codon.pair$gtcttg/cps.gtc.ttg)
cps.obs.exp.gtcttt<-log(n3.codon.pair$gtcttt/cps.gtc.ttt)











cps.obs.exp.gtgaaa<-log(n3.codon.pair$gtgaaa/cps.gtg.aaa)
cps.obs.exp.gtgaac<-log(n3.codon.pair$gtgaac/cps.gtg.aac)
cps.obs.exp.gtgaag<-log(n3.codon.pair$gtgaag/cps.gtg.aag)
cps.obs.exp.gtgaat<-log(n3.codon.pair$gtgaat/cps.gtg.aat)
cps.obs.exp.gtgaca<-log(n3.codon.pair$gtgaca/cps.gtg.aca)
cps.obs.exp.gtgacc<-log(n3.codon.pair$gtgacc/cps.gtg.acc)
cps.obs.exp.gtgacg<-log(n3.codon.pair$gtgacg/cps.gtg.acg)
cps.obs.exp.gtgact<-log(n3.codon.pair$gtgact/cps.gtg.act)
cps.obs.exp.gtgaga<-log(n3.codon.pair$gtgaga/cps.gtg.aga)
cps.obs.exp.gtgagc<-log(n3.codon.pair$gtgagc/cps.gtg.agc)
cps.obs.exp.gtgagg<-log(n3.codon.pair$gtgagg/cps.gtg.agg)
cps.obs.exp.gtgagt<-log(n3.codon.pair$gtgagt/cps.gtg.agt)
cps.obs.exp.gtgata<-log(n3.codon.pair$gtgata/cps.gtg.ata)
cps.obs.exp.gtgatc<-log(n3.codon.pair$gtgatc/cps.gtg.atc)
cps.obs.exp.gtgatg<-log(n3.codon.pair$gtgatg/cps.gtg.atg)
cps.obs.exp.gtgatt<-log(n3.codon.pair$gtgatt/cps.gtg.att)

cps.obs.exp.gtgcaa<-log(n3.codon.pair$gtgcaa/cps.gtg.caa)
cps.obs.exp.gtgcac<-log(n3.codon.pair$gtgcac/cps.gtg.cac)
cps.obs.exp.gtgcag<-log(n3.codon.pair$gtgcag/cps.gtg.cag)
cps.obs.exp.gtgcat<-log(n3.codon.pair$gtgcat/cps.gtg.cat)
cps.obs.exp.gtgcca<-log(n3.codon.pair$gtgcca/cps.gtg.cca)
cps.obs.exp.gtgccc<-log(n3.codon.pair$gtgccc/cps.gtg.ccc)
cps.obs.exp.gtgccg<-log(n3.codon.pair$gtgccg/cps.gtg.ccg)
cps.obs.exp.gtgcct<-log(n3.codon.pair$gtgcct/cps.gtg.cct)
cps.obs.exp.gtgcga<-log(n3.codon.pair$gtgcga/cps.gtg.cga)
cps.obs.exp.gtgcgc<-log(n3.codon.pair$gtgcgc/cps.gtg.cgc)
cps.obs.exp.gtgcgg<-log(n3.codon.pair$gtgcgg/cps.gtg.cgg)
cps.obs.exp.gtgcgt<-log(n3.codon.pair$gtgcgt/cps.gtg.cgt)
cps.obs.exp.gtgcta<-log(n3.codon.pair$gtgcta/cps.gtg.cta)
cps.obs.exp.gtgctc<-log(n3.codon.pair$gtgctc/cps.gtg.ctc)
cps.obs.exp.gtgctg<-log(n3.codon.pair$gtgctg/cps.gtg.ctg)
cps.obs.exp.gtgctt<-log(n3.codon.pair$gtgctt/cps.gtg.ctt)

cps.obs.exp.gtggaa<-log(n3.codon.pair$gtggaa/cps.gtg.gaa)
cps.obs.exp.gtggac<-log(n3.codon.pair$gtggac/cps.gtg.gac)
cps.obs.exp.gtggag<-log(n3.codon.pair$gtggag/cps.gtg.gag)
cps.obs.exp.gtggat<-log(n3.codon.pair$gtggat/cps.gtg.gat)
cps.obs.exp.gtggca<-log(n3.codon.pair$gtggca/cps.gtg.gca)
cps.obs.exp.gtggcc<-log(n3.codon.pair$gtggcc/cps.gtg.gcc)
cps.obs.exp.gtggcg<-log(n3.codon.pair$gtggcg/cps.gtg.gcg)
cps.obs.exp.gtggct<-log(n3.codon.pair$gtggct/cps.gtg.gct)
cps.obs.exp.gtggga<-log(n3.codon.pair$gtggga/cps.gtg.gga)
cps.obs.exp.gtgggc<-log(n3.codon.pair$gtgggc/cps.gtg.ggc)
cps.obs.exp.gtgggg<-log(n3.codon.pair$gtgggg/cps.gtg.ggg)
cps.obs.exp.gtgggt<-log(n3.codon.pair$gtgggt/cps.gtg.ggt)
cps.obs.exp.gtggta<-log(n3.codon.pair$gtggta/cps.gtg.gta)
cps.obs.exp.gtggtc<-log(n3.codon.pair$gtggtc/cps.gtg.gtc)
cps.obs.exp.gtggtg<-log(n3.codon.pair$gtggtg/cps.gtg.gtg)
cps.obs.exp.gtggtt<-log(n3.codon.pair$gtggtt/cps.gtg.gtt)

#cps.obs.exp.gtgtaa<-log(n3.codon.pair$gtgtaa/cps.gtg.taa)
cps.obs.exp.gtgtac<-log(n3.codon.pair$gtgtac/cps.gtg.tac)
#cps.obs.exp.gtgtag<-log(n3.codon.pair$gtgtag/cps.gtg.tag)
cps.obs.exp.gtgtat<-log(n3.codon.pair$gtgtat/cps.gtg.tat)
cps.obs.exp.gtgtca<-log(n3.codon.pair$gtgtca/cps.gtg.tca)
cps.obs.exp.gtgtcc<-log(n3.codon.pair$gtgtcc/cps.gtg.tcc)
cps.obs.exp.gtgtcg<-log(n3.codon.pair$gtgtcg/cps.gtg.tcg)
cps.obs.exp.gtgtct<-log(n3.codon.pair$gtgtct/cps.gtg.tct)
#cps.obs.exp.gtgtga<-log(n3.codon.pair$gtgtga/cps.gtg.tga)
cps.obs.exp.gtgtgc<-log(n3.codon.pair$gtgtgc/cps.gtg.tgc)
cps.obs.exp.gtgtgg<-log(n3.codon.pair$gtgtgg/cps.gtg.tgg)
cps.obs.exp.gtgtgt<-log(n3.codon.pair$gtgtgt/cps.gtg.tgt)
cps.obs.exp.gtgtta<-log(n3.codon.pair$gtgtta/cps.gtg.tta)
cps.obs.exp.gtgttc<-log(n3.codon.pair$gtgttc/cps.gtg.ttc)
cps.obs.exp.gtgttg<-log(n3.codon.pair$gtgttg/cps.gtg.ttg)
cps.obs.exp.gtgttt<-log(n3.codon.pair$gtgttt/cps.gtg.ttt)








cps.obs.exp.gttaaa<-log(n3.codon.pair$gttaaa/cps.gtt.aaa)
cps.obs.exp.gttaac<-log(n3.codon.pair$gttaac/cps.gtt.aac)
cps.obs.exp.gttaag<-log(n3.codon.pair$gttaag/cps.gtt.aag)
cps.obs.exp.gttaat<-log(n3.codon.pair$gttaat/cps.gtt.aat)
cps.obs.exp.gttaca<-log(n3.codon.pair$gttaca/cps.gtt.aca)
cps.obs.exp.gttacc<-log(n3.codon.pair$gttacc/cps.gtt.acc)
cps.obs.exp.gttacg<-log(n3.codon.pair$gttacg/cps.gtt.acg)
cps.obs.exp.gttact<-log(n3.codon.pair$gttact/cps.gtt.act)
cps.obs.exp.gttaga<-log(n3.codon.pair$gttaga/cps.gtt.aga)
cps.obs.exp.gttagc<-log(n3.codon.pair$gttagc/cps.gtt.agc)
cps.obs.exp.gttagg<-log(n3.codon.pair$gttagg/cps.gtt.agg)
cps.obs.exp.gttagt<-log(n3.codon.pair$gttagt/cps.gtt.agt)
cps.obs.exp.gttata<-log(n3.codon.pair$gttata/cps.gtt.ata)
cps.obs.exp.gttatc<-log(n3.codon.pair$gttatc/cps.gtt.atc)
cps.obs.exp.gttatg<-log(n3.codon.pair$gttatg/cps.gtt.atg)
cps.obs.exp.gttatt<-log(n3.codon.pair$gttatt/cps.gtt.att)

cps.obs.exp.gttcaa<-log(n3.codon.pair$gttcaa/cps.gtt.caa)
cps.obs.exp.gttcac<-log(n3.codon.pair$gttcac/cps.gtt.cac)
cps.obs.exp.gttcag<-log(n3.codon.pair$gttcag/cps.gtt.cag)
cps.obs.exp.gttcat<-log(n3.codon.pair$gttcat/cps.gtt.cat)
cps.obs.exp.gttcca<-log(n3.codon.pair$gttcca/cps.gtt.cca)
cps.obs.exp.gttccc<-log(n3.codon.pair$gttccc/cps.gtt.ccc)
cps.obs.exp.gttccg<-log(n3.codon.pair$gttccg/cps.gtt.ccg)
cps.obs.exp.gttcct<-log(n3.codon.pair$gttcct/cps.gtt.cct)
cps.obs.exp.gttcga<-log(n3.codon.pair$gttcga/cps.gtt.cga)
cps.obs.exp.gttcgc<-log(n3.codon.pair$gttcgc/cps.gtt.cgc)
cps.obs.exp.gttcgg<-log(n3.codon.pair$gttcgg/cps.gtt.cgg)
cps.obs.exp.gttcgt<-log(n3.codon.pair$gttcgt/cps.gtt.cgt)
cps.obs.exp.gttcta<-log(n3.codon.pair$gttcta/cps.gtt.cta)
cps.obs.exp.gttctc<-log(n3.codon.pair$gttctc/cps.gtt.ctc)
cps.obs.exp.gttctg<-log(n3.codon.pair$gttctg/cps.gtt.ctg)
cps.obs.exp.gttctt<-log(n3.codon.pair$gttctt/cps.gtt.ctt)

cps.obs.exp.gttgaa<-log(n3.codon.pair$gttgaa/cps.gtt.gaa)
cps.obs.exp.gttgac<-log(n3.codon.pair$gttgac/cps.gtt.gac)
cps.obs.exp.gttgag<-log(n3.codon.pair$gttgag/cps.gtt.gag)
cps.obs.exp.gttgat<-log(n3.codon.pair$gttgat/cps.gtt.gat)
cps.obs.exp.gttgca<-log(n3.codon.pair$gttgca/cps.gtt.gca)
cps.obs.exp.gttgcc<-log(n3.codon.pair$gttgcc/cps.gtt.gcc)
cps.obs.exp.gttgcg<-log(n3.codon.pair$gttgcg/cps.gtt.gcg)
cps.obs.exp.gttgct<-log(n3.codon.pair$gttgct/cps.gtt.gct)
cps.obs.exp.gttgga<-log(n3.codon.pair$gttgga/cps.gtt.gga)
cps.obs.exp.gttggc<-log(n3.codon.pair$gttggc/cps.gtt.ggc)
cps.obs.exp.gttggg<-log(n3.codon.pair$gttggg/cps.gtt.ggg)
cps.obs.exp.gttggt<-log(n3.codon.pair$gttggt/cps.gtt.ggt)
cps.obs.exp.gttgta<-log(n3.codon.pair$gttgta/cps.gtt.gta)
cps.obs.exp.gttgtc<-log(n3.codon.pair$gttgtc/cps.gtt.gtc)
cps.obs.exp.gttgtg<-log(n3.codon.pair$gttgtg/cps.gtt.gtg)
cps.obs.exp.gttgtt<-log(n3.codon.pair$gttgtt/cps.gtt.gtt)

#cps.obs.exp.gtttaa<-log(n3.codon.pair$gtttaa/cps.gtt.taa)
cps.obs.exp.gtttac<-log(n3.codon.pair$gtttac/cps.gtt.tac)
#cps.obs.exp.gtttag<-log(n3.codon.pair$gtttag/cps.gtt.tag)
cps.obs.exp.gtttat<-log(n3.codon.pair$gtttat/cps.gtt.tat)
cps.obs.exp.gtttca<-log(n3.codon.pair$gtttca/cps.gtt.tca)
cps.obs.exp.gtttcc<-log(n3.codon.pair$gtttcc/cps.gtt.tcc)
cps.obs.exp.gtttcg<-log(n3.codon.pair$gtttcg/cps.gtt.tcg)
cps.obs.exp.gtttct<-log(n3.codon.pair$gtttct/cps.gtt.tct)
#cps.obs.exp.gtttga<-log(n3.codon.pair$gtttga/cps.gtt.tga)
cps.obs.exp.gtttgc<-log(n3.codon.pair$gtttgc/cps.gtt.tgc)
cps.obs.exp.gtttgg<-log(n3.codon.pair$gtttgg/cps.gtt.tgg)
cps.obs.exp.gtttgt<-log(n3.codon.pair$gtttgt/cps.gtt.tgt)
cps.obs.exp.gtttta<-log(n3.codon.pair$gtttta/cps.gtt.tta)
cps.obs.exp.gttttc<-log(n3.codon.pair$gttttc/cps.gtt.ttc)
cps.obs.exp.gttttg<-log(n3.codon.pair$gttttg/cps.gtt.ttg)
cps.obs.exp.gttttt<-log(n3.codon.pair$gttttt/cps.gtt.ttt)

















#Stop codon


#cps.obs.exp.taaaaa<-log(n3.codon.pair$taaaaa/cps.taa.aaa)
#cps.obs.exp.taaaac<-log(n3.codon.pair$taaaac/cps.taa.aac)
#cps.obs.exp.taaaag<-log(n3.codon.pair$taaaag/cps.taa.aag)
#cps.obs.exp.taaaat<-log(n3.codon.pair$taaaat/cps.taa.aat)
#cps.obs.exp.taaaca<-log(n3.codon.pair$taaaca/cps.taa.aca)
#cps.obs.exp.taaacc<-log(n3.codon.pair$taaacc/cps.taa.acc)
#cps.obs.exp.taaacg<-log(n3.codon.pair$taaacg/cps.taa.acg)
#cps.obs.exp.taaact<-log(n3.codon.pair$taaact/cps.taa.act)
#cps.obs.exp.taaaga<-log(n3.codon.pair$taaaga/cps.taa.aga)
#cps.obs.exp.taaagc<-log(n3.codon.pair$taaagc/cps.taa.agc)
#cps.obs.exp.taaagg<-log(n3.codon.pair$taaagg/cps.taa.agg)
#cps.obs.exp.taaagt<-log(n3.codon.pair$taaagt/cps.taa.agt)
#cps.obs.exp.taaata<-log(n3.codon.pair$taaata/cps.taa.ata)
#cps.obs.exp.taaatc<-log(n3.codon.pair$taaatc/cps.taa.atc)
#cps.obs.exp.taaatg<-log(n3.codon.pair$taaatg/cps.taa.atg)
#cps.obs.exp.taaatt<-log(n3.codon.pair$taaatt/cps.taa.att)

#cps.obs.exp.taacaa<-log(n3.codon.pair$taacaa/cps.taa.caa)
#cps.obs.exp.taacac<-log(n3.codon.pair$taacac/cps.taa.cac)
#cps.obs.exp.taacag<-log(n3.codon.pair$taacag/cps.taa.cag)
#cps.obs.exp.taacat<-log(n3.codon.pair$taacat/cps.taa.cat)
#cps.obs.exp.taacca<-log(n3.codon.pair$taacca/cps.taa.cca)
#cps.obs.exp.taaccc<-log(n3.codon.pair$taaccc/cps.taa.ccc)
#cps.obs.exp.taaccg<-log(n3.codon.pair$taaccg/cps.taa.ccg)
#cps.obs.exp.taacct<-log(n3.codon.pair$taacct/cps.taa.cct)
#cps.obs.exp.taacga<-log(n3.codon.pair$taacga/cps.taa.cga)
#cps.obs.exp.taacgc<-log(n3.codon.pair$taacgc/cps.taa.cgc)
#cps.obs.exp.taacgg<-log(n3.codon.pair$taacgg/cps.taa.cgg)
#cps.obs.exp.taacgt<-log(n3.codon.pair$taacgt/cps.taa.cgt)
#cps.obs.exp.taacta<-log(n3.codon.pair$taacta/cps.taa.cta)
#cps.obs.exp.taactc<-log(n3.codon.pair$taactc/cps.taa.ctc)
#cps.obs.exp.taactg<-log(n3.codon.pair$taactg/cps.taa.ctg)
#cps.obs.exp.taactt<-log(n3.codon.pair$taactt/cps.taa.ctt)

#cps.obs.exp.taagaa<-log(n3.codon.pair$taagaa/cps.taa.gaa)
#cps.obs.exp.taagac<-log(n3.codon.pair$taagac/cps.taa.gac)
#cps.obs.exp.taagag<-log(n3.codon.pair$taagag/cps.taa.gag)
#cps.obs.exp.taagat<-log(n3.codon.pair$taagat/cps.taa.gat)
#cps.obs.exp.taagca<-log(n3.codon.pair$taagca/cps.taa.gca)
#cps.obs.exp.taagcc<-log(n3.codon.pair$taagcc/cps.taa.gcc)
#cps.obs.exp.taagcg<-log(n3.codon.pair$taagcg/cps.taa.gcg)
#cps.obs.exp.taagct<-log(n3.codon.pair$taagct/cps.taa.gct)
#cps.obs.exp.taagga<-log(n3.codon.pair$taagga/cps.taa.gga)
#cps.obs.exp.taaggc<-log(n3.codon.pair$taaggc/cps.taa.ggc)
#cps.obs.exp.taaggg<-log(n3.codon.pair$taaggg/cps.taa.ggg)
#cps.obs.exp.taaggt<-log(n3.codon.pair$taaggt/cps.taa.ggt)
#cps.obs.exp.taagta<-log(n3.codon.pair$taagta/cps.taa.gta)
#cps.obs.exp.taagtc<-log(n3.codon.pair$taagtc/cps.taa.gtc)
#cps.obs.exp.taagtg<-log(n3.codon.pair$taagtg/cps.taa.gtg)
#cps.obs.exp.taagtt<-log(n3.codon.pair$taagtt/cps.taa.gtt)

#cps.obs.exp.taataa<-log(n3.codon.pair$taataa/cps.taa.taa)
#cps.obs.exp.taatac<-log(n3.codon.pair$taatac/cps.taa.tac)
#cps.obs.exp.taatag<-log(n3.codon.pair$taatag/cps.taa.tag)
#cps.obs.exp.taatat<-log(n3.codon.pair$taatat/cps.taa.tat)
#cps.obs.exp.taatca<-log(n3.codon.pair$taatca/cps.taa.tca)
#cps.obs.exp.taatcc<-log(n3.codon.pair$taatcc/cps.taa.tcc)
#cps.obs.exp.taatcg<-log(n3.codon.pair$taatcg/cps.taa.tcg)
#cps.obs.exp.taatct<-log(n3.codon.pair$taatct/cps.taa.tct)
#cps.obs.exp.taatga<-log(n3.codon.pair$taatga/cps.taa.tga)
#cps.obs.exp.taatgc<-log(n3.codon.pair$taatgc/cps.taa.tgc)
#cps.obs.exp.taatgg<-log(n3.codon.pair$taatgg/cps.taa.tgg)
#cps.obs.exp.taatgt<-log(n3.codon.pair$taatgt/cps.taa.tgt)
#cps.obs.exp.taatta<-log(n3.codon.pair$taatta/cps.taa.tta)
#cps.obs.exp.taattc<-log(n3.codon.pair$taattc/cps.taa.ttc)
#cps.obs.exp.taattg<-log(n3.codon.pair$taattg/cps.taa.ttg)
#cps.obs.exp.taattt<-log(n3.codon.pair$taattt/cps.taa.ttt)









cps.obs.exp.tacaaa<-log(n3.codon.pair$tacaaa/cps.tac.aaa)
cps.obs.exp.tacaac<-log(n3.codon.pair$tacaac/cps.tac.aac)
cps.obs.exp.tacaag<-log(n3.codon.pair$tacaag/cps.tac.aag)
cps.obs.exp.tacaat<-log(n3.codon.pair$tacaat/cps.tac.aat)
cps.obs.exp.tacaca<-log(n3.codon.pair$tacaca/cps.tac.aca)
cps.obs.exp.tacacc<-log(n3.codon.pair$tacacc/cps.tac.acc)
cps.obs.exp.tacacg<-log(n3.codon.pair$tacacg/cps.tac.acg)
cps.obs.exp.tacact<-log(n3.codon.pair$tacact/cps.tac.act)
cps.obs.exp.tacaga<-log(n3.codon.pair$tacaga/cps.tac.aga)
cps.obs.exp.tacagc<-log(n3.codon.pair$tacagc/cps.tac.agc)
cps.obs.exp.tacagg<-log(n3.codon.pair$tacagg/cps.tac.agg)
cps.obs.exp.tacagt<-log(n3.codon.pair$tacagt/cps.tac.agt)
cps.obs.exp.tacata<-log(n3.codon.pair$tacata/cps.tac.ata)
cps.obs.exp.tacatc<-log(n3.codon.pair$tacatc/cps.tac.atc)
cps.obs.exp.tacatg<-log(n3.codon.pair$tacatg/cps.tac.atg)
cps.obs.exp.tacatt<-log(n3.codon.pair$tacatt/cps.tac.att)

cps.obs.exp.taccaa<-log(n3.codon.pair$taccaa/cps.tac.caa)
cps.obs.exp.taccac<-log(n3.codon.pair$taccac/cps.tac.cac)
cps.obs.exp.taccag<-log(n3.codon.pair$taccag/cps.tac.cag)
cps.obs.exp.taccat<-log(n3.codon.pair$taccat/cps.tac.cat)
cps.obs.exp.taccca<-log(n3.codon.pair$taccca/cps.tac.cca)
cps.obs.exp.tacccc<-log(n3.codon.pair$tacccc/cps.tac.ccc)
cps.obs.exp.tacccg<-log(n3.codon.pair$tacccg/cps.tac.ccg)
cps.obs.exp.taccct<-log(n3.codon.pair$taccct/cps.tac.cct)
cps.obs.exp.taccga<-log(n3.codon.pair$taccga/cps.tac.cga)
cps.obs.exp.taccgc<-log(n3.codon.pair$taccgc/cps.tac.cgc)
cps.obs.exp.taccgg<-log(n3.codon.pair$taccgg/cps.tac.cgg)
cps.obs.exp.taccgt<-log(n3.codon.pair$taccgt/cps.tac.cgt)
cps.obs.exp.taccta<-log(n3.codon.pair$taccta/cps.tac.cta)
cps.obs.exp.tacctc<-log(n3.codon.pair$tacctc/cps.tac.ctc)
cps.obs.exp.tacctg<-log(n3.codon.pair$tacctg/cps.tac.ctg)
cps.obs.exp.tacctt<-log(n3.codon.pair$tacctt/cps.tac.ctt)

cps.obs.exp.tacgaa<-log(n3.codon.pair$tacgaa/cps.tac.gaa)
cps.obs.exp.tacgac<-log(n3.codon.pair$tacgac/cps.tac.gac)
cps.obs.exp.tacgag<-log(n3.codon.pair$tacgag/cps.tac.gag)
cps.obs.exp.tacgat<-log(n3.codon.pair$tacgat/cps.tac.gat)
cps.obs.exp.tacgca<-log(n3.codon.pair$tacgca/cps.tac.gca)
cps.obs.exp.tacgcc<-log(n3.codon.pair$tacgcc/cps.tac.gcc)
cps.obs.exp.tacgcg<-log(n3.codon.pair$tacgcg/cps.tac.gcg)
cps.obs.exp.tacgct<-log(n3.codon.pair$tacgct/cps.tac.gct)
cps.obs.exp.tacgga<-log(n3.codon.pair$tacgga/cps.tac.gga)
cps.obs.exp.tacggc<-log(n3.codon.pair$tacggc/cps.tac.ggc)
cps.obs.exp.tacggg<-log(n3.codon.pair$tacggg/cps.tac.ggg)
cps.obs.exp.tacggt<-log(n3.codon.pair$tacggt/cps.tac.ggt)
cps.obs.exp.tacgta<-log(n3.codon.pair$tacgta/cps.tac.gta)
cps.obs.exp.tacgtc<-log(n3.codon.pair$tacgtc/cps.tac.gtc)
cps.obs.exp.tacgtg<-log(n3.codon.pair$tacgtg/cps.tac.gtg)
cps.obs.exp.tacgtt<-log(n3.codon.pair$tacgtt/cps.tac.gtt)

#cps.obs.exp.tactaa<-log(n3.codon.pair$tactaa/cps.tac.taa)
cps.obs.exp.tactac<-log(n3.codon.pair$tactac/cps.tac.tac)
#cps.obs.exp.tactag<-log(n3.codon.pair$tactag/cps.tac.tag)
cps.obs.exp.tactat<-log(n3.codon.pair$tactat/cps.tac.tat)
cps.obs.exp.tactca<-log(n3.codon.pair$tactca/cps.tac.tca)
cps.obs.exp.tactcc<-log(n3.codon.pair$tactcc/cps.tac.tcc)
cps.obs.exp.tactcg<-log(n3.codon.pair$tactcg/cps.tac.tcg)
cps.obs.exp.tactct<-log(n3.codon.pair$tactct/cps.tac.tct)
#cps.obs.exp.tactga<-log(n3.codon.pair$tactga/cps.tac.tga)
cps.obs.exp.tactgc<-log(n3.codon.pair$tactgc/cps.tac.tgc)
cps.obs.exp.tactgg<-log(n3.codon.pair$tactgg/cps.tac.tgg)
cps.obs.exp.tactgt<-log(n3.codon.pair$tactgt/cps.tac.tgt)
cps.obs.exp.tactta<-log(n3.codon.pair$tactta/cps.tac.tta)
cps.obs.exp.tacttc<-log(n3.codon.pair$tacttc/cps.tac.ttc)
cps.obs.exp.tacttg<-log(n3.codon.pair$tacttg/cps.tac.ttg)
cps.obs.exp.tacttt<-log(n3.codon.pair$tacttt/cps.tac.ttt)







#Stop Codon 


#cps.obs.exp.tagaaa<-log(n3.codon.pair$tagaaa/cps.tag.aaa)
#cps.obs.exp.tagaac<-log(n3.codon.pair$tagaac/cps.tag.aac)
#cps.obs.exp.tagaag<-log(n3.codon.pair$tagaag/cps.tag.aag)
#cps.obs.exp.tagaat<-log(n3.codon.pair$tagaat/cps.tag.aat)
#cps.obs.exp.tagaca<-log(n3.codon.pair$tagaca/cps.tag.aca)
#cps.obs.exp.tagacc<-log(n3.codon.pair$tagacc/cps.tag.acc)
#cps.obs.exp.tagacg<-log(n3.codon.pair$tagacg/cps.tag.acg)
#cps.obs.exp.tagact<-log(n3.codon.pair$tagact/cps.tag.act)
#cps.obs.exp.tagaga<-log(n3.codon.pair$tagaga/cps.tag.aga)
#cps.obs.exp.tagagc<-log(n3.codon.pair$tagagc/cps.tag.agc)
#cps.obs.exp.tagagg<-log(n3.codon.pair$tagagg/cps.tag.agg)
#cps.obs.exp.tagagt<-log(n3.codon.pair$tagagt/cps.tag.agt)
#cps.obs.exp.tagata<-log(n3.codon.pair$tagata/cps.tag.ata)
#cps.obs.exp.tagatc<-log(n3.codon.pair$tagatc/cps.tag.atc)
#cps.obs.exp.tagatg<-log(n3.codon.pair$tagatg/cps.tag.atg)
#cps.obs.exp.tagatt<-log(n3.codon.pair$tagatt/cps.tag.att)

#cps.obs.exp.tagcaa<-log(n3.codon.pair$tagcaa/cps.tag.caa)
#cps.obs.exp.tagcac<-log(n3.codon.pair$tagcac/cps.tag.cac)
#cps.obs.exp.tagcag<-log(n3.codon.pair$tagcag/cps.tag.cag)
#cps.obs.exp.tagcat<-log(n3.codon.pair$tagcat/cps.tag.cat)
#cps.obs.exp.tagcca<-log(n3.codon.pair$tagcca/cps.tag.cca)
#cps.obs.exp.tagccc<-log(n3.codon.pair$tagccc/cps.tag.ccc)
#cps.obs.exp.tagccg<-log(n3.codon.pair$tagccg/cps.tag.ccg)
#cps.obs.exp.tagcct<-log(n3.codon.pair$tagcct/cps.tag.cct)
#cps.obs.exp.tagcga<-log(n3.codon.pair$tagcga/cps.tag.cga)
#cps.obs.exp.tagcgc<-log(n3.codon.pair$tagcgc/cps.tag.cgc)
#cps.obs.exp.tagcgg<-log(n3.codon.pair$tagcgg/cps.tag.cgg)
#cps.obs.exp.tagcgt<-log(n3.codon.pair$tagcgt/cps.tag.cgt)
#cps.obs.exp.tagcta<-log(n3.codon.pair$tagcta/cps.tag.cta)
#cps.obs.exp.tagctc<-log(n3.codon.pair$tagctc/cps.tag.ctc)
#cps.obs.exp.tagctg<-log(n3.codon.pair$tagctg/cps.tag.ctg)
#cps.obs.exp.tagctt<-log(n3.codon.pair$tagctt/cps.tag.ctt)

#cps.obs.exp.taggaa<-log(n3.codon.pair$taggaa/cps.tag.gaa)
#cps.obs.exp.taggac<-log(n3.codon.pair$taggac/cps.tag.gac)
#cps.obs.exp.taggag<-log(n3.codon.pair$taggag/cps.tag.gag)
#cps.obs.exp.taggat<-log(n3.codon.pair$taggat/cps.tag.gat)
#cps.obs.exp.taggca<-log(n3.codon.pair$taggca/cps.tag.gca)
#cps.obs.exp.taggcc<-log(n3.codon.pair$taggcc/cps.tag.gcc)
#cps.obs.exp.taggcg<-log(n3.codon.pair$taggcg/cps.tag.gcg)
#cps.obs.exp.taggct<-log(n3.codon.pair$taggct/cps.tag.gct)
#cps.obs.exp.taggga<-log(n3.codon.pair$taggga/cps.tag.gga)
#cps.obs.exp.tagggc<-log(n3.codon.pair$tagggc/cps.tag.ggc)
#cps.obs.exp.tagggg<-log(n3.codon.pair$tagggg/cps.tag.ggg)
#cps.obs.exp.tagggt<-log(n3.codon.pair$tagggt/cps.tag.ggt)
#cps.obs.exp.taggta<-log(n3.codon.pair$taggta/cps.tag.gta)
#cps.obs.exp.taggtc<-log(n3.codon.pair$taggtc/cps.tag.gtc)
#cps.obs.exp.taggtg<-log(n3.codon.pair$taggtg/cps.tag.gtg)
#cps.obs.exp.taggtt<-log(n3.codon.pair$taggtt/cps.tag.gtt)

#cps.obs.exp.tagtaa<-log(n3.codon.pair$tagtaa/cps.tag.taa)
#cps.obs.exp.tagtac<-log(n3.codon.pair$tagtac/cps.tag.tac)
#cps.obs.exp.tagtag<-log(n3.codon.pair$tagtag/cps.tag.tag)
#cps.obs.exp.tagtat<-log(n3.codon.pair$tagtat/cps.tag.tat)
#cps.obs.exp.tagtca<-log(n3.codon.pair$tagtca/cps.tag.tca)
#cps.obs.exp.tagtcc<-log(n3.codon.pair$tagtcc/cps.tag.tcc)
#cps.obs.exp.tagtcg<-log(n3.codon.pair$tagtcg/cps.tag.tcg)
#cps.obs.exp.tagtct<-log(n3.codon.pair$tagtct/cps.tag.tct)
#cps.obs.exp.tagtga<-log(n3.codon.pair$tagtga/cps.tag.tga)
#cps.obs.exp.tagtgc<-log(n3.codon.pair$tagtgc/cps.tag.tgc)
#cps.obs.exp.tagtgg<-log(n3.codon.pair$tagtgg/cps.tag.tgg)
#cps.obs.exp.tagtgt<-log(n3.codon.pair$tagtgt/cps.tag.tgt)
#cps.obs.exp.tagtta<-log(n3.codon.pair$tagtta/cps.tag.tta)
#cps.obs.exp.tagttc<-log(n3.codon.pair$tagttc/cps.tag.ttc)
#cps.obs.exp.tagttg<-log(n3.codon.pair$tagttg/cps.tag.ttg)
#cps.obs.exp.tagttt<-log(n3.codon.pair$tagttt/cps.tag.ttt)









cps.obs.exp.tataaa<-log(n3.codon.pair$tataaa/cps.tat.aaa)
cps.obs.exp.tataac<-log(n3.codon.pair$tataac/cps.tat.aac)
cps.obs.exp.tataag<-log(n3.codon.pair$tataag/cps.tat.aag)
cps.obs.exp.tataat<-log(n3.codon.pair$tataat/cps.tat.aat)
cps.obs.exp.tataca<-log(n3.codon.pair$tataca/cps.tat.aca)
cps.obs.exp.tatacc<-log(n3.codon.pair$tatacc/cps.tat.acc)
cps.obs.exp.tatacg<-log(n3.codon.pair$tatacg/cps.tat.acg)
cps.obs.exp.tatact<-log(n3.codon.pair$tatact/cps.tat.act)
cps.obs.exp.tataga<-log(n3.codon.pair$tataga/cps.tat.aga)
cps.obs.exp.tatagc<-log(n3.codon.pair$tatagc/cps.tat.agc)
cps.obs.exp.tatagg<-log(n3.codon.pair$tatagg/cps.tat.agg)
cps.obs.exp.tatagt<-log(n3.codon.pair$tatagt/cps.tat.agt)
cps.obs.exp.tatata<-log(n3.codon.pair$tatata/cps.tat.ata)
cps.obs.exp.tatatc<-log(n3.codon.pair$tatatc/cps.tat.atc)
cps.obs.exp.tatatg<-log(n3.codon.pair$tatatg/cps.tat.atg)
cps.obs.exp.tatatt<-log(n3.codon.pair$tatatt/cps.tat.att)

cps.obs.exp.tatcaa<-log(n3.codon.pair$tatcaa/cps.tat.caa)
cps.obs.exp.tatcac<-log(n3.codon.pair$tatcac/cps.tat.cac)
cps.obs.exp.tatcag<-log(n3.codon.pair$tatcag/cps.tat.cag)
cps.obs.exp.tatcat<-log(n3.codon.pair$tatcat/cps.tat.cat)
cps.obs.exp.tatcca<-log(n3.codon.pair$tatcca/cps.tat.cca)
cps.obs.exp.tatccc<-log(n3.codon.pair$tatccc/cps.tat.ccc)
cps.obs.exp.tatccg<-log(n3.codon.pair$tatccg/cps.tat.ccg)
cps.obs.exp.tatcct<-log(n3.codon.pair$tatcct/cps.tat.cct)
cps.obs.exp.tatcga<-log(n3.codon.pair$tatcga/cps.tat.cga)
cps.obs.exp.tatcgc<-log(n3.codon.pair$tatcgc/cps.tat.cgc)
cps.obs.exp.tatcgg<-log(n3.codon.pair$tatcgg/cps.tat.cgg)
cps.obs.exp.tatcgt<-log(n3.codon.pair$tatcgt/cps.tat.cgt)
cps.obs.exp.tatcta<-log(n3.codon.pair$tatcta/cps.tat.cta)
cps.obs.exp.tatctc<-log(n3.codon.pair$tatctc/cps.tat.ctc)
cps.obs.exp.tatctg<-log(n3.codon.pair$tatctg/cps.tat.ctg)
cps.obs.exp.tatctt<-log(n3.codon.pair$tatctt/cps.tat.ctt)

cps.obs.exp.tatgaa<-log(n3.codon.pair$tatgaa/cps.tat.gaa)
cps.obs.exp.tatgac<-log(n3.codon.pair$tatgac/cps.tat.gac)
cps.obs.exp.tatgag<-log(n3.codon.pair$tatgag/cps.tat.gag)
cps.obs.exp.tatgat<-log(n3.codon.pair$tatgat/cps.tat.gat)
cps.obs.exp.tatgca<-log(n3.codon.pair$tatgca/cps.tat.gca)
cps.obs.exp.tatgcc<-log(n3.codon.pair$tatgcc/cps.tat.gcc)
cps.obs.exp.tatgcg<-log(n3.codon.pair$tatgcg/cps.tat.gcg)
cps.obs.exp.tatgct<-log(n3.codon.pair$tatgct/cps.tat.gct)
cps.obs.exp.tatgga<-log(n3.codon.pair$tatgga/cps.tat.gga)
cps.obs.exp.tatggc<-log(n3.codon.pair$tatggc/cps.tat.ggc)
cps.obs.exp.tatggg<-log(n3.codon.pair$tatggg/cps.tat.ggg)
cps.obs.exp.tatggt<-log(n3.codon.pair$tatggt/cps.tat.ggt)
cps.obs.exp.tatgta<-log(n3.codon.pair$tatgta/cps.tat.gta)
cps.obs.exp.tatgtc<-log(n3.codon.pair$tatgtc/cps.tat.gtc)
cps.obs.exp.tatgtg<-log(n3.codon.pair$tatgtg/cps.tat.gtg)
cps.obs.exp.tatgtt<-log(n3.codon.pair$tatgtt/cps.tat.gtt)

#cps.obs.exp.tattaa<-log(n3.codon.pair$tattaa/cps.tat.taa)
cps.obs.exp.tattac<-log(n3.codon.pair$tattac/cps.tat.tac)
#cps.obs.exp.tattag<-log(n3.codon.pair$tattag/cps.tat.tag)
cps.obs.exp.tattat<-log(n3.codon.pair$tattat/cps.tat.tat)
cps.obs.exp.tattca<-log(n3.codon.pair$tattca/cps.tat.tca)
cps.obs.exp.tattcc<-log(n3.codon.pair$tattcc/cps.tat.tcc)
cps.obs.exp.tattcg<-log(n3.codon.pair$tattcg/cps.tat.tcg)
cps.obs.exp.tattct<-log(n3.codon.pair$tattct/cps.tat.tct)
#cps.obs.exp.tattga<-log(n3.codon.pair$tattga/cps.tat.tga)
cps.obs.exp.tattgc<-log(n3.codon.pair$tattgc/cps.tat.tgc)
cps.obs.exp.tattgg<-log(n3.codon.pair$tattgg/cps.tat.tgg)
cps.obs.exp.tattgt<-log(n3.codon.pair$tattgt/cps.tat.tgt)
cps.obs.exp.tattta<-log(n3.codon.pair$tattta/cps.tat.tta)
cps.obs.exp.tatttc<-log(n3.codon.pair$tatttc/cps.tat.ttc)
cps.obs.exp.tatttg<-log(n3.codon.pair$tatttg/cps.tat.ttg)
cps.obs.exp.tatttt<-log(n3.codon.pair$tatttt/cps.tat.ttt)

















cps.obs.exp.tcaaaa<-log(n3.codon.pair$tcaaaa/cps.tca.aaa)
cps.obs.exp.tcaaac<-log(n3.codon.pair$tcaaac/cps.tca.aac)
cps.obs.exp.tcaaag<-log(n3.codon.pair$tcaaag/cps.tca.aag)
cps.obs.exp.tcaaat<-log(n3.codon.pair$tcaaat/cps.tca.aat)
cps.obs.exp.tcaaca<-log(n3.codon.pair$tcaaca/cps.tca.aca)
cps.obs.exp.tcaacc<-log(n3.codon.pair$tcaacc/cps.tca.acc)
cps.obs.exp.tcaacg<-log(n3.codon.pair$tcaacg/cps.tca.acg)
cps.obs.exp.tcaact<-log(n3.codon.pair$tcaact/cps.tca.act)
cps.obs.exp.tcaaga<-log(n3.codon.pair$tcaaga/cps.tca.aga)
cps.obs.exp.tcaagc<-log(n3.codon.pair$tcaagc/cps.tca.agc)
cps.obs.exp.tcaagg<-log(n3.codon.pair$tcaagg/cps.tca.agg)
cps.obs.exp.tcaagt<-log(n3.codon.pair$tcaagt/cps.tca.agt)
cps.obs.exp.tcaata<-log(n3.codon.pair$tcaata/cps.tca.ata)
cps.obs.exp.tcaatc<-log(n3.codon.pair$tcaatc/cps.tca.atc)
cps.obs.exp.tcaatg<-log(n3.codon.pair$tcaatg/cps.tca.atg)
cps.obs.exp.tcaatt<-log(n3.codon.pair$tcaatt/cps.tca.att)

cps.obs.exp.tcacaa<-log(n3.codon.pair$tcacaa/cps.tca.caa)
cps.obs.exp.tcacac<-log(n3.codon.pair$tcacac/cps.tca.cac)
cps.obs.exp.tcacag<-log(n3.codon.pair$tcacag/cps.tca.cag)
cps.obs.exp.tcacat<-log(n3.codon.pair$tcacat/cps.tca.cat)
cps.obs.exp.tcacca<-log(n3.codon.pair$tcacca/cps.tca.cca)
cps.obs.exp.tcaccc<-log(n3.codon.pair$tcaccc/cps.tca.ccc)
cps.obs.exp.tcaccg<-log(n3.codon.pair$tcaccg/cps.tca.ccg)
cps.obs.exp.tcacct<-log(n3.codon.pair$tcacct/cps.tca.cct)
cps.obs.exp.tcacga<-log(n3.codon.pair$tcacga/cps.tca.cga)
cps.obs.exp.tcacgc<-log(n3.codon.pair$tcacgc/cps.tca.cgc)
cps.obs.exp.tcacgg<-log(n3.codon.pair$tcacgg/cps.tca.cgg)
cps.obs.exp.tcacgt<-log(n3.codon.pair$tcacgt/cps.tca.cgt)
cps.obs.exp.tcacta<-log(n3.codon.pair$tcacta/cps.tca.cta)
cps.obs.exp.tcactc<-log(n3.codon.pair$tcactc/cps.tca.ctc)
cps.obs.exp.tcactg<-log(n3.codon.pair$tcactg/cps.tca.ctg)
cps.obs.exp.tcactt<-log(n3.codon.pair$tcactt/cps.tca.ctt)

cps.obs.exp.tcagaa<-log(n3.codon.pair$tcagaa/cps.tca.gaa)
cps.obs.exp.tcagac<-log(n3.codon.pair$tcagac/cps.tca.gac)
cps.obs.exp.tcagag<-log(n3.codon.pair$tcagag/cps.tca.gag)
cps.obs.exp.tcagat<-log(n3.codon.pair$tcagat/cps.tca.gat)
cps.obs.exp.tcagca<-log(n3.codon.pair$tcagca/cps.tca.gca)
cps.obs.exp.tcagcc<-log(n3.codon.pair$tcagcc/cps.tca.gcc)
cps.obs.exp.tcagcg<-log(n3.codon.pair$tcagcg/cps.tca.gcg)
cps.obs.exp.tcagct<-log(n3.codon.pair$tcagct/cps.tca.gct)
cps.obs.exp.tcagga<-log(n3.codon.pair$tcagga/cps.tca.gga)
cps.obs.exp.tcaggc<-log(n3.codon.pair$tcaggc/cps.tca.ggc)
cps.obs.exp.tcaggg<-log(n3.codon.pair$tcaggg/cps.tca.ggg)
cps.obs.exp.tcaggt<-log(n3.codon.pair$tcaggt/cps.tca.ggt)
cps.obs.exp.tcagta<-log(n3.codon.pair$tcagta/cps.tca.gta)
cps.obs.exp.tcagtc<-log(n3.codon.pair$tcagtc/cps.tca.gtc)
cps.obs.exp.tcagtg<-log(n3.codon.pair$tcagtg/cps.tca.gtg)
cps.obs.exp.tcagtt<-log(n3.codon.pair$tcagtt/cps.tca.gtt)

#cps.obs.exp.tcataa<-log(n3.codon.pair$tcataa/cps.tca.taa)
cps.obs.exp.tcatac<-log(n3.codon.pair$tcatac/cps.tca.tac)
#cps.obs.exp.tcatag<-log(n3.codon.pair$tcatag/cps.tca.tag)
cps.obs.exp.tcatat<-log(n3.codon.pair$tcatat/cps.tca.tat)
cps.obs.exp.tcatca<-log(n3.codon.pair$tcatca/cps.tca.tca)
cps.obs.exp.tcatcc<-log(n3.codon.pair$tcatcc/cps.tca.tcc)
cps.obs.exp.tcatcg<-log(n3.codon.pair$tcatcg/cps.tca.tcg)
cps.obs.exp.tcatct<-log(n3.codon.pair$tcatct/cps.tca.tct)
#cps.obs.exp.tcatga<-log(n3.codon.pair$tcatga/cps.tca.tga)
cps.obs.exp.tcatgc<-log(n3.codon.pair$tcatgc/cps.tca.tgc)
cps.obs.exp.tcatgg<-log(n3.codon.pair$tcatgg/cps.tca.tgg)
cps.obs.exp.tcatgt<-log(n3.codon.pair$tcatgt/cps.tca.tgt)
cps.obs.exp.tcatta<-log(n3.codon.pair$tcatta/cps.tca.tta)
cps.obs.exp.tcattc<-log(n3.codon.pair$tcattc/cps.tca.ttc)
cps.obs.exp.tcattg<-log(n3.codon.pair$tcattg/cps.tca.ttg)
cps.obs.exp.tcattt<-log(n3.codon.pair$tcattt/cps.tca.ttt)








cps.obs.exp.tccaaa<-log(n3.codon.pair$tccaaa/cps.tcc.aaa)
cps.obs.exp.tccaac<-log(n3.codon.pair$tccaac/cps.tcc.aac)
cps.obs.exp.tccaag<-log(n3.codon.pair$tccaag/cps.tcc.aag)
cps.obs.exp.tccaat<-log(n3.codon.pair$tccaat/cps.tcc.aat)
cps.obs.exp.tccaca<-log(n3.codon.pair$tccaca/cps.tcc.aca)
cps.obs.exp.tccacc<-log(n3.codon.pair$tccacc/cps.tcc.acc)
cps.obs.exp.tccacg<-log(n3.codon.pair$tccacg/cps.tcc.acg)
cps.obs.exp.tccact<-log(n3.codon.pair$tccact/cps.tcc.act)
cps.obs.exp.tccaga<-log(n3.codon.pair$tccaga/cps.tcc.aga)
cps.obs.exp.tccagc<-log(n3.codon.pair$tccagc/cps.tcc.agc)
cps.obs.exp.tccagg<-log(n3.codon.pair$tccagg/cps.tcc.agg)
cps.obs.exp.tccagt<-log(n3.codon.pair$tccagt/cps.tcc.agt)
cps.obs.exp.tccata<-log(n3.codon.pair$tccata/cps.tcc.ata)
cps.obs.exp.tccatc<-log(n3.codon.pair$tccatc/cps.tcc.atc)
cps.obs.exp.tccatg<-log(n3.codon.pair$tccatg/cps.tcc.atg)
cps.obs.exp.tccatt<-log(n3.codon.pair$tccatt/cps.tcc.att)

cps.obs.exp.tcccaa<-log(n3.codon.pair$tcccaa/cps.tcc.caa)
cps.obs.exp.tcccac<-log(n3.codon.pair$tcccac/cps.tcc.cac)
cps.obs.exp.tcccag<-log(n3.codon.pair$tcccag/cps.tcc.cag)
cps.obs.exp.tcccat<-log(n3.codon.pair$tcccat/cps.tcc.cat)
cps.obs.exp.tcccca<-log(n3.codon.pair$tcccca/cps.tcc.cca)
cps.obs.exp.tccccc<-log(n3.codon.pair$tccccc/cps.tcc.ccc)
cps.obs.exp.tccccg<-log(n3.codon.pair$tccccg/cps.tcc.ccg)
cps.obs.exp.tcccct<-log(n3.codon.pair$tcccct/cps.tcc.cct)
cps.obs.exp.tcccga<-log(n3.codon.pair$tcccga/cps.tcc.cga)
cps.obs.exp.tcccgc<-log(n3.codon.pair$tcccgc/cps.tcc.cgc)
cps.obs.exp.tcccgg<-log(n3.codon.pair$tcccgg/cps.tcc.cgg)
cps.obs.exp.tcccgt<-log(n3.codon.pair$tcccgt/cps.tcc.cgt)
cps.obs.exp.tcccta<-log(n3.codon.pair$tcccta/cps.tcc.cta)
cps.obs.exp.tccctc<-log(n3.codon.pair$tccctc/cps.tcc.ctc)
cps.obs.exp.tccctg<-log(n3.codon.pair$tccctg/cps.tcc.ctg)
cps.obs.exp.tccctt<-log(n3.codon.pair$tccctt/cps.tcc.ctt)

cps.obs.exp.tccgaa<-log(n3.codon.pair$tccgaa/cps.tcc.gaa)
cps.obs.exp.tccgac<-log(n3.codon.pair$tccgac/cps.tcc.gac)
cps.obs.exp.tccgag<-log(n3.codon.pair$tccgag/cps.tcc.gag)
cps.obs.exp.tccgat<-log(n3.codon.pair$tccgat/cps.tcc.gat)
cps.obs.exp.tccgca<-log(n3.codon.pair$tccgca/cps.tcc.gca)
cps.obs.exp.tccgcc<-log(n3.codon.pair$tccgcc/cps.tcc.gcc)
cps.obs.exp.tccgcg<-log(n3.codon.pair$tccgcg/cps.tcc.gcg)
cps.obs.exp.tccgct<-log(n3.codon.pair$tccgct/cps.tcc.gct)
cps.obs.exp.tccgga<-log(n3.codon.pair$tccgga/cps.tcc.gga)
cps.obs.exp.tccggc<-log(n3.codon.pair$tccggc/cps.tcc.ggc)
cps.obs.exp.tccggg<-log(n3.codon.pair$tccggg/cps.tcc.ggg)
cps.obs.exp.tccggt<-log(n3.codon.pair$tccggt/cps.tcc.ggt)
cps.obs.exp.tccgta<-log(n3.codon.pair$tccgta/cps.tcc.gta)
cps.obs.exp.tccgtc<-log(n3.codon.pair$tccgtc/cps.tcc.gtc)
cps.obs.exp.tccgtg<-log(n3.codon.pair$tccgtg/cps.tcc.gtg)
cps.obs.exp.tccgtt<-log(n3.codon.pair$tccgtt/cps.tcc.gtt)

#cps.obs.exp.tcctaa<-log(n3.codon.pair$tcctaa/cps.tcc.taa)
cps.obs.exp.tcctac<-log(n3.codon.pair$tcctac/cps.tcc.tac)
#cps.obs.exp.tcctag<-log(n3.codon.pair$tcctag/cps.tcc.tag)
cps.obs.exp.tcctat<-log(n3.codon.pair$tcctat/cps.tcc.tat)
cps.obs.exp.tcctca<-log(n3.codon.pair$tcctca/cps.tcc.tca)
cps.obs.exp.tcctcc<-log(n3.codon.pair$tcctcc/cps.tcc.tcc)
cps.obs.exp.tcctcg<-log(n3.codon.pair$tcctcg/cps.tcc.tcg)
cps.obs.exp.tcctct<-log(n3.codon.pair$tcctct/cps.tcc.tct)
#cps.obs.exp.tcctga<-log(n3.codon.pair$tcctga/cps.tcc.tga)
cps.obs.exp.tcctgc<-log(n3.codon.pair$tcctgc/cps.tcc.tgc)
cps.obs.exp.tcctgg<-log(n3.codon.pair$tcctgg/cps.tcc.tgg)
cps.obs.exp.tcctgt<-log(n3.codon.pair$tcctgt/cps.tcc.tgt)
cps.obs.exp.tcctta<-log(n3.codon.pair$tcctta/cps.tcc.tta)
cps.obs.exp.tccttc<-log(n3.codon.pair$tccttc/cps.tcc.ttc)
cps.obs.exp.tccttg<-log(n3.codon.pair$tccttg/cps.tcc.ttg)
cps.obs.exp.tccttt<-log(n3.codon.pair$tccttt/cps.tcc.ttt)









cps.obs.exp.tcgaaa<-log(n3.codon.pair$tcgaaa/cps.tcg.aaa)
cps.obs.exp.tcgaac<-log(n3.codon.pair$tcgaac/cps.tcg.aac)
cps.obs.exp.tcgaag<-log(n3.codon.pair$tcgaag/cps.tcg.aag)
cps.obs.exp.tcgaat<-log(n3.codon.pair$tcgaat/cps.tcg.aat)
cps.obs.exp.tcgaca<-log(n3.codon.pair$tcgaca/cps.tcg.aca)
cps.obs.exp.tcgacc<-log(n3.codon.pair$tcgacc/cps.tcg.acc)
cps.obs.exp.tcgacg<-log(n3.codon.pair$tcgacg/cps.tcg.acg)
cps.obs.exp.tcgact<-log(n3.codon.pair$tcgact/cps.tcg.act)
cps.obs.exp.tcgaga<-log(n3.codon.pair$tcgaga/cps.tcg.aga)
cps.obs.exp.tcgagc<-log(n3.codon.pair$tcgagc/cps.tcg.agc)
cps.obs.exp.tcgagg<-log(n3.codon.pair$tcgagg/cps.tcg.agg)
cps.obs.exp.tcgagt<-log(n3.codon.pair$tcgagt/cps.tcg.agt)
cps.obs.exp.tcgata<-log(n3.codon.pair$tcgata/cps.tcg.ata)
cps.obs.exp.tcgatc<-log(n3.codon.pair$tcgatc/cps.tcg.atc)
cps.obs.exp.tcgatg<-log(n3.codon.pair$tcgatg/cps.tcg.atg)
cps.obs.exp.tcgatt<-log(n3.codon.pair$tcgatt/cps.tcg.att)

cps.obs.exp.tcgcaa<-log(n3.codon.pair$tcgcaa/cps.tcg.caa)
cps.obs.exp.tcgcac<-log(n3.codon.pair$tcgcac/cps.tcg.cac)
cps.obs.exp.tcgcag<-log(n3.codon.pair$tcgcag/cps.tcg.cag)
cps.obs.exp.tcgcat<-log(n3.codon.pair$tcgcat/cps.tcg.cat)
cps.obs.exp.tcgcca<-log(n3.codon.pair$tcgcca/cps.tcg.cca)
cps.obs.exp.tcgccc<-log(n3.codon.pair$tcgccc/cps.tcg.ccc)
cps.obs.exp.tcgccg<-log(n3.codon.pair$tcgccg/cps.tcg.ccg)
cps.obs.exp.tcgcct<-log(n3.codon.pair$tcgcct/cps.tcg.cct)
cps.obs.exp.tcgcga<-log(n3.codon.pair$tcgcga/cps.tcg.cga)
cps.obs.exp.tcgcgc<-log(n3.codon.pair$tcgcgc/cps.tcg.cgc)
cps.obs.exp.tcgcgg<-log(n3.codon.pair$tcgcgg/cps.tcg.cgg)
cps.obs.exp.tcgcgt<-log(n3.codon.pair$tcgcgt/cps.tcg.cgt)
cps.obs.exp.tcgcta<-log(n3.codon.pair$tcgcta/cps.tcg.cta)
cps.obs.exp.tcgctc<-log(n3.codon.pair$tcgctc/cps.tcg.ctc)
cps.obs.exp.tcgctg<-log(n3.codon.pair$tcgctg/cps.tcg.ctg)
cps.obs.exp.tcgctt<-log(n3.codon.pair$tcgctt/cps.tcg.ctt)

cps.obs.exp.tcggaa<-log(n3.codon.pair$tcggaa/cps.tcg.gaa)
cps.obs.exp.tcggac<-log(n3.codon.pair$tcggac/cps.tcg.gac)
cps.obs.exp.tcggag<-log(n3.codon.pair$tcggag/cps.tcg.gag)
cps.obs.exp.tcggat<-log(n3.codon.pair$tcggat/cps.tcg.gat)
cps.obs.exp.tcggca<-log(n3.codon.pair$tcggca/cps.tcg.gca)
cps.obs.exp.tcggcc<-log(n3.codon.pair$tcggcc/cps.tcg.gcc)
cps.obs.exp.tcggcg<-log(n3.codon.pair$tcggcg/cps.tcg.gcg)
cps.obs.exp.tcggct<-log(n3.codon.pair$tcggct/cps.tcg.gct)
cps.obs.exp.tcggga<-log(n3.codon.pair$tcggga/cps.tcg.gga)
cps.obs.exp.tcgggc<-log(n3.codon.pair$tcgggc/cps.tcg.ggc)
cps.obs.exp.tcgggg<-log(n3.codon.pair$tcgggg/cps.tcg.ggg)
cps.obs.exp.tcgggt<-log(n3.codon.pair$tcgggt/cps.tcg.ggt)
cps.obs.exp.tcggta<-log(n3.codon.pair$tcggta/cps.tcg.gta)
cps.obs.exp.tcggtc<-log(n3.codon.pair$tcggtc/cps.tcg.gtc)
cps.obs.exp.tcggtg<-log(n3.codon.pair$tcggtg/cps.tcg.gtg)
cps.obs.exp.tcggtt<-log(n3.codon.pair$tcggtt/cps.tcg.gtt)

#cps.obs.exp.tcgtaa<-log(n3.codon.pair$tcgtaa/cps.tcg.taa)
cps.obs.exp.tcgtac<-log(n3.codon.pair$tcgtac/cps.tcg.tac)
#cps.obs.exp.tcgtag<-log(n3.codon.pair$tcgtag/cps.tcg.tag)
cps.obs.exp.tcgtat<-log(n3.codon.pair$tcgtat/cps.tcg.tat)
cps.obs.exp.tcgtca<-log(n3.codon.pair$tcgtca/cps.tcg.tca)
cps.obs.exp.tcgtcc<-log(n3.codon.pair$tcgtcc/cps.tcg.tcc)
cps.obs.exp.tcgtcg<-log(n3.codon.pair$tcgtcg/cps.tcg.tcg)
cps.obs.exp.tcgtct<-log(n3.codon.pair$tcgtct/cps.tcg.tct)
#cps.obs.exp.tcgtga<-log(n3.codon.pair$tcgtga/cps.tcg.tga)
cps.obs.exp.tcgtgc<-log(n3.codon.pair$tcgtgc/cps.tcg.tgc)
cps.obs.exp.tcgtgg<-log(n3.codon.pair$tcgtgg/cps.tcg.tgg)
cps.obs.exp.tcgtgt<-log(n3.codon.pair$tcgtgt/cps.tcg.tgt)
cps.obs.exp.tcgtta<-log(n3.codon.pair$tcgtta/cps.tcg.tta)
cps.obs.exp.tcgttc<-log(n3.codon.pair$tcgttc/cps.tcg.ttc)
cps.obs.exp.tcgttg<-log(n3.codon.pair$tcgttg/cps.tcg.ttg)
cps.obs.exp.tcgttt<-log(n3.codon.pair$tcgttt/cps.tcg.ttt)









cps.obs.exp.tctaaa<-log(n3.codon.pair$tctaaa/cps.tct.aaa)
cps.obs.exp.tctaac<-log(n3.codon.pair$tctaac/cps.tct.aac)
cps.obs.exp.tctaag<-log(n3.codon.pair$tctaag/cps.tct.aag)
cps.obs.exp.tctaat<-log(n3.codon.pair$tctaat/cps.tct.aat)
cps.obs.exp.tctaca<-log(n3.codon.pair$tctaca/cps.tct.aca)
cps.obs.exp.tctacc<-log(n3.codon.pair$tctacc/cps.tct.acc)
cps.obs.exp.tctacg<-log(n3.codon.pair$tctacg/cps.tct.acg)
cps.obs.exp.tctact<-log(n3.codon.pair$tctact/cps.tct.act)
cps.obs.exp.tctaga<-log(n3.codon.pair$tctaga/cps.tct.aga)
cps.obs.exp.tctagc<-log(n3.codon.pair$tctagc/cps.tct.agc)
cps.obs.exp.tctagg<-log(n3.codon.pair$tctagg/cps.tct.agg)
cps.obs.exp.tctagt<-log(n3.codon.pair$tctagt/cps.tct.agt)
cps.obs.exp.tctata<-log(n3.codon.pair$tctata/cps.tct.ata)
cps.obs.exp.tctatc<-log(n3.codon.pair$tctatc/cps.tct.atc)
cps.obs.exp.tctatg<-log(n3.codon.pair$tctatg/cps.tct.atg)
cps.obs.exp.tctatt<-log(n3.codon.pair$tctatt/cps.tct.att)

cps.obs.exp.tctcaa<-log(n3.codon.pair$tctcaa/cps.tct.caa)
cps.obs.exp.tctcac<-log(n3.codon.pair$tctcac/cps.tct.cac)
cps.obs.exp.tctcag<-log(n3.codon.pair$tctcag/cps.tct.cag)
cps.obs.exp.tctcat<-log(n3.codon.pair$tctcat/cps.tct.cat)
cps.obs.exp.tctcca<-log(n3.codon.pair$tctcca/cps.tct.cca)
cps.obs.exp.tctccc<-log(n3.codon.pair$tctccc/cps.tct.ccc)
cps.obs.exp.tctccg<-log(n3.codon.pair$tctccg/cps.tct.ccg)
cps.obs.exp.tctcct<-log(n3.codon.pair$tctcct/cps.tct.cct)
cps.obs.exp.tctcga<-log(n3.codon.pair$tctcga/cps.tct.cga)
cps.obs.exp.tctcgc<-log(n3.codon.pair$tctcgc/cps.tct.cgc)
cps.obs.exp.tctcgg<-log(n3.codon.pair$tctcgg/cps.tct.cgg)
cps.obs.exp.tctcgt<-log(n3.codon.pair$tctcgt/cps.tct.cgt)
cps.obs.exp.tctcta<-log(n3.codon.pair$tctcta/cps.tct.cta)
cps.obs.exp.tctctc<-log(n3.codon.pair$tctctc/cps.tct.ctc)
cps.obs.exp.tctctg<-log(n3.codon.pair$tctctg/cps.tct.ctg)
cps.obs.exp.tctctt<-log(n3.codon.pair$tctctt/cps.tct.ctt)

cps.obs.exp.tctgaa<-log(n3.codon.pair$tctgaa/cps.tct.gaa)
cps.obs.exp.tctgac<-log(n3.codon.pair$tctgac/cps.tct.gac)
cps.obs.exp.tctgag<-log(n3.codon.pair$tctgag/cps.tct.gag)
cps.obs.exp.tctgat<-log(n3.codon.pair$tctgat/cps.tct.gat)
cps.obs.exp.tctgca<-log(n3.codon.pair$tctgca/cps.tct.gca)
cps.obs.exp.tctgcc<-log(n3.codon.pair$tctgcc/cps.tct.gcc)
cps.obs.exp.tctgcg<-log(n3.codon.pair$tctgcg/cps.tct.gcg)
cps.obs.exp.tctgct<-log(n3.codon.pair$tctgct/cps.tct.gct)
cps.obs.exp.tctgga<-log(n3.codon.pair$tctgga/cps.tct.gga)
cps.obs.exp.tctggc<-log(n3.codon.pair$tctggc/cps.tct.ggc)
cps.obs.exp.tctggg<-log(n3.codon.pair$tctggg/cps.tct.ggg)
cps.obs.exp.tctggt<-log(n3.codon.pair$tctggt/cps.tct.ggt)
cps.obs.exp.tctgta<-log(n3.codon.pair$tctgta/cps.tct.gta)
cps.obs.exp.tctgtc<-log(n3.codon.pair$tctgtc/cps.tct.gtc)
cps.obs.exp.tctgtg<-log(n3.codon.pair$tctgtg/cps.tct.gtg)
cps.obs.exp.tctgtt<-log(n3.codon.pair$tctgtt/cps.tct.gtt)

#cps.obs.exp.tcttaa<-log(n3.codon.pair$tcttaa/cps.tct.taa)
cps.obs.exp.tcttac<-log(n3.codon.pair$tcttac/cps.tct.tac)
#cps.obs.exp.tcttag<-log(n3.codon.pair$tcttag/cps.tct.tag)
cps.obs.exp.tcttat<-log(n3.codon.pair$tcttat/cps.tct.tat)
cps.obs.exp.tcttca<-log(n3.codon.pair$tcttca/cps.tct.tca)
cps.obs.exp.tcttcc<-log(n3.codon.pair$tcttcc/cps.tct.tcc)
cps.obs.exp.tcttcg<-log(n3.codon.pair$tcttcg/cps.tct.tcg)
cps.obs.exp.tcttct<-log(n3.codon.pair$tcttct/cps.tct.tct)
#cps.obs.exp.tcttga<-log(n3.codon.pair$tcttga/cps.tct.tga)
cps.obs.exp.tcttgc<-log(n3.codon.pair$tcttgc/cps.tct.tgc)
cps.obs.exp.tcttgg<-log(n3.codon.pair$tcttgg/cps.tct.tgg)
cps.obs.exp.tcttgt<-log(n3.codon.pair$tcttgt/cps.tct.tgt)
cps.obs.exp.tcttta<-log(n3.codon.pair$tcttta/cps.tct.tta)
cps.obs.exp.tctttc<-log(n3.codon.pair$tctttc/cps.tct.ttc)
cps.obs.exp.tctttg<-log(n3.codon.pair$tctttg/cps.tct.ttg)
cps.obs.exp.tctttt<-log(n3.codon.pair$tctttt/cps.tct.ttt)
















#Stop Codon 


#cps.obs.exp.tgaaaa<-log(n3.codon.pair$tgaaaa/cps.tga.aaa)
#cps.obs.exp.tgaaac<-log(n3.codon.pair$tgaaac/cps.tga.aac)
#cps.obs.exp.tgaaag<-log(n3.codon.pair$tgaaag/cps.tga.aag)
#cps.obs.exp.tgaaat<-log(n3.codon.pair$tgaaat/cps.tga.aat)
#cps.obs.exp.tgaaca<-log(n3.codon.pair$tgaaca/cps.tga.aca)
#cps.obs.exp.tgaacc<-log(n3.codon.pair$tgaacc/cps.tga.acc)
#cps.obs.exp.tgaacg<-log(n3.codon.pair$tgaacg/cps.tga.acg)
#cps.obs.exp.tgaact<-log(n3.codon.pair$tgaact/cps.tga.act)
#cps.obs.exp.tgaaga<-log(n3.codon.pair$tgaaga/cps.tga.aga)
#cps.obs.exp.tgaagc<-log(n3.codon.pair$tgaagc/cps.tga.agc)
#cps.obs.exp.tgaagg<-log(n3.codon.pair$tgaagg/cps.tga.agg)
#cps.obs.exp.tgaagt<-log(n3.codon.pair$tgaagt/cps.tga.agt)
#cps.obs.exp.tgaata<-log(n3.codon.pair$tgaata/cps.tga.ata)
#cps.obs.exp.tgaatc<-log(n3.codon.pair$tgaatc/cps.tga.atc)
#cps.obs.exp.tgaatg<-log(n3.codon.pair$tgaatg/cps.tga.atg)
#cps.obs.exp.tgaatt<-log(n3.codon.pair$tgaatt/cps.tga.att)

#cps.obs.exp.tgacaa<-log(n3.codon.pair$tgacaa/cps.tga.caa)
#cps.obs.exp.tgacac<-log(n3.codon.pair$tgacac/cps.tga.cac)
#cps.obs.exp.tgacag<-log(n3.codon.pair$tgacag/cps.tga.cag)
#cps.obs.exp.tgacat<-log(n3.codon.pair$tgacat/cps.tga.cat)
#cps.obs.exp.tgacca<-log(n3.codon.pair$tgacca/cps.tga.cca)
#cps.obs.exp.tgaccc<-log(n3.codon.pair$tgaccc/cps.tga.ccc)
#cps.obs.exp.tgaccg<-log(n3.codon.pair$tgaccg/cps.tga.ccg)
#cps.obs.exp.tgacct<-log(n3.codon.pair$tgacct/cps.tga.cct)
#cps.obs.exp.tgacga<-log(n3.codon.pair$tgacga/cps.tga.cga)
#cps.obs.exp.tgacgc<-log(n3.codon.pair$tgacgc/cps.tga.cgc)
#cps.obs.exp.tgacgg<-log(n3.codon.pair$tgacgg/cps.tga.cgg)
#cps.obs.exp.tgacgt<-log(n3.codon.pair$tgacgt/cps.tga.cgt)
#cps.obs.exp.tgacta<-log(n3.codon.pair$tgacta/cps.tga.cta)
#cps.obs.exp.tgactc<-log(n3.codon.pair$tgactc/cps.tga.ctc)
#cps.obs.exp.tgactg<-log(n3.codon.pair$tgactg/cps.tga.ctg)
#cps.obs.exp.tgactt<-log(n3.codon.pair$tgactt/cps.tga.ctt)

#cps.obs.exp.tgagaa<-log(n3.codon.pair$tgagaa/cps.tga.gaa)
#cps.obs.exp.tgagac<-log(n3.codon.pair$tgagac/cps.tga.gac)
#cps.obs.exp.tgagag<-log(n3.codon.pair$tgagag/cps.tga.gag)
#cps.obs.exp.tgagat<-log(n3.codon.pair$tgagat/cps.tga.gat)
#cps.obs.exp.tgagca<-log(n3.codon.pair$tgagca/cps.tga.gca)
#cps.obs.exp.tgagcc<-log(n3.codon.pair$tgagcc/cps.tga.gcc)
#cps.obs.exp.tgagcg<-log(n3.codon.pair$tgagcg/cps.tga.gcg)
#cps.obs.exp.tgagct<-log(n3.codon.pair$tgagct/cps.tga.gct)
#cps.obs.exp.tgagga<-log(n3.codon.pair$tgagga/cps.tga.gga)
#cps.obs.exp.tgaggc<-log(n3.codon.pair$tgaggc/cps.tga.ggc)
#cps.obs.exp.tgaggg<-log(n3.codon.pair$tgaggg/cps.tga.ggg)
#cps.obs.exp.tgaggt<-log(n3.codon.pair$tgaggt/cps.tga.ggt)
#cps.obs.exp.tgagta<-log(n3.codon.pair$tgagta/cps.tga.gta)
#cps.obs.exp.tgagtc<-log(n3.codon.pair$tgagtc/cps.tga.gtc)
#cps.obs.exp.tgagtg<-log(n3.codon.pair$tgagtg/cps.tga.gtg)
#cps.obs.exp.tgagtt<-log(n3.codon.pair$tgagtt/cps.tga.gtt)

#cps.obs.exp.tgataa<-log(n3.codon.pair$tgataa/cps.tga.taa)
#cps.obs.exp.tgatac<-log(n3.codon.pair$tgatac/cps.tga.tac)
#cps.obs.exp.tgatag<-log(n3.codon.pair$tgatag/cps.tga.tag)
#cps.obs.exp.tgatat<-log(n3.codon.pair$tgatat/cps.tga.tat)
#cps.obs.exp.tgatca<-log(n3.codon.pair$tgatca/cps.tga.tca)
#cps.obs.exp.tgatcc<-log(n3.codon.pair$tgatcc/cps.tga.tcc)
#cps.obs.exp.tgatcg<-log(n3.codon.pair$tgatcg/cps.tga.tcg)
#cps.obs.exp.tgatct<-log(n3.codon.pair$tgatct/cps.tga.tct)
#cps.obs.exp.tgatga<-log(n3.codon.pair$tgatga/cps.tga.tga)
#cps.obs.exp.tgatgc<-log(n3.codon.pair$tgatgc/cps.tga.tgc)
#cps.obs.exp.tgatgg<-log(n3.codon.pair$tgatgg/cps.tga.tgg)
#cps.obs.exp.tgatgt<-log(n3.codon.pair$tgatgt/cps.tga.tgt)
#cps.obs.exp.tgatta<-log(n3.codon.pair$tgatta/cps.tga.tta)
#cps.obs.exp.tgattc<-log(n3.codon.pair$tgattc/cps.tga.ttc)
#cps.obs.exp.tgattg<-log(n3.codon.pair$tgattg/cps.tga.ttg)
#cps.obs.exp.tgattt<-log(n3.codon.pair$tgattt/cps.tga.ttt)









cps.obs.exp.tgcaaa<-log(n3.codon.pair$tgcaaa/cps.tgc.aaa)
cps.obs.exp.tgcaac<-log(n3.codon.pair$tgcaac/cps.tgc.aac)
cps.obs.exp.tgcaag<-log(n3.codon.pair$tgcaag/cps.tgc.aag)
cps.obs.exp.tgcaat<-log(n3.codon.pair$tgcaat/cps.tgc.aat)
cps.obs.exp.tgcaca<-log(n3.codon.pair$tgcaca/cps.tgc.aca)
cps.obs.exp.tgcacc<-log(n3.codon.pair$tgcacc/cps.tgc.acc)
cps.obs.exp.tgcacg<-log(n3.codon.pair$tgcacg/cps.tgc.acg)
cps.obs.exp.tgcact<-log(n3.codon.pair$tgcact/cps.tgc.act)
cps.obs.exp.tgcaga<-log(n3.codon.pair$tgcaga/cps.tgc.aga)
cps.obs.exp.tgcagc<-log(n3.codon.pair$tgcagc/cps.tgc.agc)
cps.obs.exp.tgcagg<-log(n3.codon.pair$tgcagg/cps.tgc.agg)
cps.obs.exp.tgcagt<-log(n3.codon.pair$tgcagt/cps.tgc.agt)
cps.obs.exp.tgcata<-log(n3.codon.pair$tgcata/cps.tgc.ata)
cps.obs.exp.tgcatc<-log(n3.codon.pair$tgcatc/cps.tgc.atc)
cps.obs.exp.tgcatg<-log(n3.codon.pair$tgcatg/cps.tgc.atg)
cps.obs.exp.tgcatt<-log(n3.codon.pair$tgcatt/cps.tgc.att)

cps.obs.exp.tgccaa<-log(n3.codon.pair$tgccaa/cps.tgc.caa)
cps.obs.exp.tgccac<-log(n3.codon.pair$tgccac/cps.tgc.cac)
cps.obs.exp.tgccag<-log(n3.codon.pair$tgccag/cps.tgc.cag)
cps.obs.exp.tgccat<-log(n3.codon.pair$tgccat/cps.tgc.cat)
cps.obs.exp.tgccca<-log(n3.codon.pair$tgccca/cps.tgc.cca)
cps.obs.exp.tgcccc<-log(n3.codon.pair$tgcccc/cps.tgc.ccc)
cps.obs.exp.tgcccg<-log(n3.codon.pair$tgcccg/cps.tgc.ccg)
cps.obs.exp.tgccct<-log(n3.codon.pair$tgccct/cps.tgc.cct)
cps.obs.exp.tgccga<-log(n3.codon.pair$tgccga/cps.tgc.cga)
cps.obs.exp.tgccgc<-log(n3.codon.pair$tgccgc/cps.tgc.cgc)
cps.obs.exp.tgccgg<-log(n3.codon.pair$tgccgg/cps.tgc.cgg)
cps.obs.exp.tgccgt<-log(n3.codon.pair$tgccgt/cps.tgc.cgt)
cps.obs.exp.tgccta<-log(n3.codon.pair$tgccta/cps.tgc.cta)
cps.obs.exp.tgcctc<-log(n3.codon.pair$tgcctc/cps.tgc.ctc)
cps.obs.exp.tgcctg<-log(n3.codon.pair$tgcctg/cps.tgc.ctg)
cps.obs.exp.tgcctt<-log(n3.codon.pair$tgcctt/cps.tgc.ctt)

cps.obs.exp.tgcgaa<-log(n3.codon.pair$tgcgaa/cps.tgc.gaa)
cps.obs.exp.tgcgac<-log(n3.codon.pair$tgcgac/cps.tgc.gac)
cps.obs.exp.tgcgag<-log(n3.codon.pair$tgcgag/cps.tgc.gag)
cps.obs.exp.tgcgat<-log(n3.codon.pair$tgcgat/cps.tgc.gat)
cps.obs.exp.tgcgca<-log(n3.codon.pair$tgcgca/cps.tgc.gca)
cps.obs.exp.tgcgcc<-log(n3.codon.pair$tgcgcc/cps.tgc.gcc)
cps.obs.exp.tgcgcg<-log(n3.codon.pair$tgcgcg/cps.tgc.gcg)
cps.obs.exp.tgcgct<-log(n3.codon.pair$tgcgct/cps.tgc.gct)
cps.obs.exp.tgcgga<-log(n3.codon.pair$tgcgga/cps.tgc.gga)
cps.obs.exp.tgcggc<-log(n3.codon.pair$tgcggc/cps.tgc.ggc)
cps.obs.exp.tgcggg<-log(n3.codon.pair$tgcggg/cps.tgc.ggg)
cps.obs.exp.tgcggt<-log(n3.codon.pair$tgcggt/cps.tgc.ggt)
cps.obs.exp.tgcgta<-log(n3.codon.pair$tgcgta/cps.tgc.gta)
cps.obs.exp.tgcgtc<-log(n3.codon.pair$tgcgtc/cps.tgc.gtc)
cps.obs.exp.tgcgtg<-log(n3.codon.pair$tgcgtg/cps.tgc.gtg)
cps.obs.exp.tgcgtt<-log(n3.codon.pair$tgcgtt/cps.tgc.gtt)

#cps.obs.exp.tgctaa<-log(n3.codon.pair$tgctaa/cps.tgc.taa)
cps.obs.exp.tgctac<-log(n3.codon.pair$tgctac/cps.tgc.tac)
#cps.obs.exp.tgctag<-log(n3.codon.pair$tgctag/cps.tgc.tag)
cps.obs.exp.tgctat<-log(n3.codon.pair$tgctat/cps.tgc.tat)
cps.obs.exp.tgctca<-log(n3.codon.pair$tgctca/cps.tgc.tca)
cps.obs.exp.tgctcc<-log(n3.codon.pair$tgctcc/cps.tgc.tcc)
cps.obs.exp.tgctcg<-log(n3.codon.pair$tgctcg/cps.tgc.tcg)
cps.obs.exp.tgctct<-log(n3.codon.pair$tgctct/cps.tgc.tct)
#cps.obs.exp.tgctga<-log(n3.codon.pair$tgctga/cps.tgc.tga)
cps.obs.exp.tgctgc<-log(n3.codon.pair$tgctgc/cps.tgc.tgc)
cps.obs.exp.tgctgg<-log(n3.codon.pair$tgctgg/cps.tgc.tgg)
cps.obs.exp.tgctgt<-log(n3.codon.pair$tgctgt/cps.tgc.tgt)
cps.obs.exp.tgctta<-log(n3.codon.pair$tgctta/cps.tgc.tta)
cps.obs.exp.tgcttc<-log(n3.codon.pair$tgcttc/cps.tgc.ttc)
cps.obs.exp.tgcttg<-log(n3.codon.pair$tgcttg/cps.tgc.ttg)
cps.obs.exp.tgcttt<-log(n3.codon.pair$tgcttt/cps.tgc.ttt)










cps.obs.exp.tggaaa<-log(n3.codon.pair$tggaaa/cps.tgg.aaa)
cps.obs.exp.tggaac<-log(n3.codon.pair$tggaac/cps.tgg.aac)
cps.obs.exp.tggaag<-log(n3.codon.pair$tggaag/cps.tgg.aag)
cps.obs.exp.tggaat<-log(n3.codon.pair$tggaat/cps.tgg.aat)
cps.obs.exp.tggaca<-log(n3.codon.pair$tggaca/cps.tgg.aca)
cps.obs.exp.tggacc<-log(n3.codon.pair$tggacc/cps.tgg.acc)
cps.obs.exp.tggacg<-log(n3.codon.pair$tggacg/cps.tgg.acg)
cps.obs.exp.tggact<-log(n3.codon.pair$tggact/cps.tgg.act)
cps.obs.exp.tggaga<-log(n3.codon.pair$tggaga/cps.tgg.aga)
cps.obs.exp.tggagc<-log(n3.codon.pair$tggagc/cps.tgg.agc)
cps.obs.exp.tggagg<-log(n3.codon.pair$tggagg/cps.tgg.agg)
cps.obs.exp.tggagt<-log(n3.codon.pair$tggagt/cps.tgg.agt)
cps.obs.exp.tggata<-log(n3.codon.pair$tggata/cps.tgg.ata)
cps.obs.exp.tggatc<-log(n3.codon.pair$tggatc/cps.tgg.atc)
cps.obs.exp.tggatg<-log(n3.codon.pair$tggatg/cps.tgg.atg)
cps.obs.exp.tggatt<-log(n3.codon.pair$tggatt/cps.tgg.att)

cps.obs.exp.tggcaa<-log(n3.codon.pair$tggcaa/cps.tgg.caa)
cps.obs.exp.tggcac<-log(n3.codon.pair$tggcac/cps.tgg.cac)
cps.obs.exp.tggcag<-log(n3.codon.pair$tggcag/cps.tgg.cag)
cps.obs.exp.tggcat<-log(n3.codon.pair$tggcat/cps.tgg.cat)
cps.obs.exp.tggcca<-log(n3.codon.pair$tggcca/cps.tgg.cca)
cps.obs.exp.tggccc<-log(n3.codon.pair$tggccc/cps.tgg.ccc)
cps.obs.exp.tggccg<-log(n3.codon.pair$tggccg/cps.tgg.ccg)
cps.obs.exp.tggcct<-log(n3.codon.pair$tggcct/cps.tgg.cct)
cps.obs.exp.tggcga<-log(n3.codon.pair$tggcga/cps.tgg.cga)
cps.obs.exp.tggcgc<-log(n3.codon.pair$tggcgc/cps.tgg.cgc)
cps.obs.exp.tggcgg<-log(n3.codon.pair$tggcgg/cps.tgg.cgg)
cps.obs.exp.tggcgt<-log(n3.codon.pair$tggcgt/cps.tgg.cgt)
cps.obs.exp.tggcta<-log(n3.codon.pair$tggcta/cps.tgg.cta)
cps.obs.exp.tggctc<-log(n3.codon.pair$tggctc/cps.tgg.ctc)
cps.obs.exp.tggctg<-log(n3.codon.pair$tggctg/cps.tgg.ctg)
cps.obs.exp.tggctt<-log(n3.codon.pair$tggctt/cps.tgg.ctt)

cps.obs.exp.tgggaa<-log(n3.codon.pair$tgggaa/cps.tgg.gaa)
cps.obs.exp.tgggac<-log(n3.codon.pair$tgggac/cps.tgg.gac)
cps.obs.exp.tgggag<-log(n3.codon.pair$tgggag/cps.tgg.gag)
cps.obs.exp.tgggat<-log(n3.codon.pair$tgggat/cps.tgg.gat)
cps.obs.exp.tgggca<-log(n3.codon.pair$tgggca/cps.tgg.gca)
cps.obs.exp.tgggcc<-log(n3.codon.pair$tgggcc/cps.tgg.gcc)
cps.obs.exp.tgggcg<-log(n3.codon.pair$tgggcg/cps.tgg.gcg)
cps.obs.exp.tgggct<-log(n3.codon.pair$tgggct/cps.tgg.gct)
cps.obs.exp.tgggga<-log(n3.codon.pair$tgggga/cps.tgg.gga)
cps.obs.exp.tggggc<-log(n3.codon.pair$tggggc/cps.tgg.ggc)
cps.obs.exp.tggggg<-log(n3.codon.pair$tggggg/cps.tgg.ggg)
cps.obs.exp.tggggt<-log(n3.codon.pair$tggggt/cps.tgg.ggt)
cps.obs.exp.tgggta<-log(n3.codon.pair$tgggta/cps.tgg.gta)
cps.obs.exp.tgggtc<-log(n3.codon.pair$tgggtc/cps.tgg.gtc)
cps.obs.exp.tgggtg<-log(n3.codon.pair$tgggtg/cps.tgg.gtg)
cps.obs.exp.tgggtt<-log(n3.codon.pair$tgggtt/cps.tgg.gtt)

#cps.obs.exp.tggtaa<-log(n3.codon.pair$tggtaa/cps.tgg.taa)
cps.obs.exp.tggtac<-log(n3.codon.pair$tggtac/cps.tgg.tac)
#cps.obs.exp.tggtag<-log(n3.codon.pair$tggtag/cps.tgg.tag)
cps.obs.exp.tggtat<-log(n3.codon.pair$tggtat/cps.tgg.tat)
cps.obs.exp.tggtca<-log(n3.codon.pair$tggtca/cps.tgg.tca)
cps.obs.exp.tggtcc<-log(n3.codon.pair$tggtcc/cps.tgg.tcc)
cps.obs.exp.tggtcg<-log(n3.codon.pair$tggtcg/cps.tgg.tcg)
cps.obs.exp.tggtct<-log(n3.codon.pair$tggtct/cps.tgg.tct)
#cps.obs.exp.tggtga<-log(n3.codon.pair$tggtga/cps.tgg.tga)
cps.obs.exp.tggtgc<-log(n3.codon.pair$tggtgc/cps.tgg.tgc)
cps.obs.exp.tggtgg<-log(n3.codon.pair$tggtgg/cps.tgg.tgg)
cps.obs.exp.tggtgt<-log(n3.codon.pair$tggtgt/cps.tgg.tgt)
cps.obs.exp.tggtta<-log(n3.codon.pair$tggtta/cps.tgg.tta)
cps.obs.exp.tggttc<-log(n3.codon.pair$tggttc/cps.tgg.ttc)
cps.obs.exp.tggttg<-log(n3.codon.pair$tggttg/cps.tgg.ttg)
cps.obs.exp.tggttt<-log(n3.codon.pair$tggttt/cps.tgg.ttt)










cps.obs.exp.tgtaaa<-log(n3.codon.pair$tgtaaa/cps.tgt.aaa)
cps.obs.exp.tgtaac<-log(n3.codon.pair$tgtaac/cps.tgt.aac)
cps.obs.exp.tgtaag<-log(n3.codon.pair$tgtaag/cps.tgt.aag)
cps.obs.exp.tgtaat<-log(n3.codon.pair$tgtaat/cps.tgt.aat)
cps.obs.exp.tgtaca<-log(n3.codon.pair$tgtaca/cps.tgt.aca)
cps.obs.exp.tgtacc<-log(n3.codon.pair$tgtacc/cps.tgt.acc)
cps.obs.exp.tgtacg<-log(n3.codon.pair$tgtacg/cps.tgt.acg)
cps.obs.exp.tgtact<-log(n3.codon.pair$tgtact/cps.tgt.act)
cps.obs.exp.tgtaga<-log(n3.codon.pair$tgtaga/cps.tgt.aga)
cps.obs.exp.tgtagc<-log(n3.codon.pair$tgtagc/cps.tgt.agc)
cps.obs.exp.tgtagg<-log(n3.codon.pair$tgtagg/cps.tgt.agg)
cps.obs.exp.tgtagt<-log(n3.codon.pair$tgtagt/cps.tgt.agt)
cps.obs.exp.tgtata<-log(n3.codon.pair$tgtata/cps.tgt.ata)
cps.obs.exp.tgtatc<-log(n3.codon.pair$tgtatc/cps.tgt.atc)
cps.obs.exp.tgtatg<-log(n3.codon.pair$tgtatg/cps.tgt.atg)
cps.obs.exp.tgtatt<-log(n3.codon.pair$tgtatt/cps.tgt.att)

cps.obs.exp.tgtcaa<-log(n3.codon.pair$tgtcaa/cps.tgt.caa)
cps.obs.exp.tgtcac<-log(n3.codon.pair$tgtcac/cps.tgt.cac)
cps.obs.exp.tgtcag<-log(n3.codon.pair$tgtcag/cps.tgt.cag)
cps.obs.exp.tgtcat<-log(n3.codon.pair$tgtcat/cps.tgt.cat)
cps.obs.exp.tgtcca<-log(n3.codon.pair$tgtcca/cps.tgt.cca)
cps.obs.exp.tgtccc<-log(n3.codon.pair$tgtccc/cps.tgt.ccc)
cps.obs.exp.tgtccg<-log(n3.codon.pair$tgtccg/cps.tgt.ccg)
cps.obs.exp.tgtcct<-log(n3.codon.pair$tgtcct/cps.tgt.cct)
cps.obs.exp.tgtcga<-log(n3.codon.pair$tgtcga/cps.tgt.cga)
cps.obs.exp.tgtcgc<-log(n3.codon.pair$tgtcgc/cps.tgt.cgc)
cps.obs.exp.tgtcgg<-log(n3.codon.pair$tgtcgg/cps.tgt.cgg)
cps.obs.exp.tgtcgt<-log(n3.codon.pair$tgtcgt/cps.tgt.cgt)
cps.obs.exp.tgtcta<-log(n3.codon.pair$tgtcta/cps.tgt.cta)
cps.obs.exp.tgtctc<-log(n3.codon.pair$tgtctc/cps.tgt.ctc)
cps.obs.exp.tgtctg<-log(n3.codon.pair$tgtctg/cps.tgt.ctg)
cps.obs.exp.tgtctt<-log(n3.codon.pair$tgtctt/cps.tgt.ctt)

cps.obs.exp.tgtgaa<-log(n3.codon.pair$tgtgaa/cps.tgt.gaa)
cps.obs.exp.tgtgac<-log(n3.codon.pair$tgtgac/cps.tgt.gac)
cps.obs.exp.tgtgag<-log(n3.codon.pair$tgtgag/cps.tgt.gag)
cps.obs.exp.tgtgat<-log(n3.codon.pair$tgtgat/cps.tgt.gat)
cps.obs.exp.tgtgca<-log(n3.codon.pair$tgtgca/cps.tgt.gca)
cps.obs.exp.tgtgcc<-log(n3.codon.pair$tgtgcc/cps.tgt.gcc)
cps.obs.exp.tgtgcg<-log(n3.codon.pair$tgtgcg/cps.tgt.gcg)
cps.obs.exp.tgtgct<-log(n3.codon.pair$tgtgct/cps.tgt.gct)
cps.obs.exp.tgtgga<-log(n3.codon.pair$tgtgga/cps.tgt.gga)
cps.obs.exp.tgtggc<-log(n3.codon.pair$tgtggc/cps.tgt.ggc)
cps.obs.exp.tgtggg<-log(n3.codon.pair$tgtggg/cps.tgt.ggg)
cps.obs.exp.tgtggt<-log(n3.codon.pair$tgtggt/cps.tgt.ggt)
cps.obs.exp.tgtgta<-log(n3.codon.pair$tgtgta/cps.tgt.gta)
cps.obs.exp.tgtgtc<-log(n3.codon.pair$tgtgtc/cps.tgt.gtc)
cps.obs.exp.tgtgtg<-log(n3.codon.pair$tgtgtg/cps.tgt.gtg)
cps.obs.exp.tgtgtt<-log(n3.codon.pair$tgtgtt/cps.tgt.gtt)

#cps.obs.exp.tgttaa<-log(n3.codon.pair$tgttaa/cps.tgt.taa)
cps.obs.exp.tgttac<-log(n3.codon.pair$tgttac/cps.tgt.tac)
#cps.obs.exp.tgttag<-log(n3.codon.pair$tgttag/cps.tgt.tag)
cps.obs.exp.tgttat<-log(n3.codon.pair$tgttat/cps.tgt.tat)
cps.obs.exp.tgttca<-log(n3.codon.pair$tgttca/cps.tgt.tca)
cps.obs.exp.tgttcc<-log(n3.codon.pair$tgttcc/cps.tgt.tcc)
cps.obs.exp.tgttcg<-log(n3.codon.pair$tgttcg/cps.tgt.tcg)
cps.obs.exp.tgttct<-log(n3.codon.pair$tgttct/cps.tgt.tct)
#cps.obs.exp.tgttga<-log(n3.codon.pair$tgttga/cps.tgt.tga)
cps.obs.exp.tgttgc<-log(n3.codon.pair$tgttgc/cps.tgt.tgc)
cps.obs.exp.tgttgg<-log(n3.codon.pair$tgttgg/cps.tgt.tgg)
cps.obs.exp.tgttgt<-log(n3.codon.pair$tgttgt/cps.tgt.tgt)
cps.obs.exp.tgttta<-log(n3.codon.pair$tgttta/cps.tgt.tta)
cps.obs.exp.tgtttc<-log(n3.codon.pair$tgtttc/cps.tgt.ttc)
cps.obs.exp.tgtttg<-log(n3.codon.pair$tgtttg/cps.tgt.ttg)
cps.obs.exp.tgtttt<-log(n3.codon.pair$tgtttt/cps.tgt.ttt)




















cps.obs.exp.ttaaaa<-log(n3.codon.pair$ttaaaa/cps.tta.aaa)
cps.obs.exp.ttaaac<-log(n3.codon.pair$ttaaac/cps.tta.aac)
cps.obs.exp.ttaaag<-log(n3.codon.pair$ttaaag/cps.tta.aag)
cps.obs.exp.ttaaat<-log(n3.codon.pair$ttaaat/cps.tta.aat)
cps.obs.exp.ttaaca<-log(n3.codon.pair$ttaaca/cps.tta.aca)
cps.obs.exp.ttaacc<-log(n3.codon.pair$ttaacc/cps.tta.acc)
cps.obs.exp.ttaacg<-log(n3.codon.pair$ttaacg/cps.tta.acg)
cps.obs.exp.ttaact<-log(n3.codon.pair$ttaact/cps.tta.act)
cps.obs.exp.ttaaga<-log(n3.codon.pair$ttaaga/cps.tta.aga)
cps.obs.exp.ttaagc<-log(n3.codon.pair$ttaagc/cps.tta.agc)
cps.obs.exp.ttaagg<-log(n3.codon.pair$ttaagg/cps.tta.agg)
cps.obs.exp.ttaagt<-log(n3.codon.pair$ttaagt/cps.tta.agt)
cps.obs.exp.ttaata<-log(n3.codon.pair$ttaata/cps.tta.ata)
cps.obs.exp.ttaatc<-log(n3.codon.pair$ttaatc/cps.tta.atc)
cps.obs.exp.ttaatg<-log(n3.codon.pair$ttaatg/cps.tta.atg)
cps.obs.exp.ttaatt<-log(n3.codon.pair$ttaatt/cps.tta.att)

cps.obs.exp.ttacaa<-log(n3.codon.pair$ttacaa/cps.tta.caa)
cps.obs.exp.ttacac<-log(n3.codon.pair$ttacac/cps.tta.cac)
cps.obs.exp.ttacag<-log(n3.codon.pair$ttacag/cps.tta.cag)
cps.obs.exp.ttacat<-log(n3.codon.pair$ttacat/cps.tta.cat)
cps.obs.exp.ttacca<-log(n3.codon.pair$ttacca/cps.tta.cca)
cps.obs.exp.ttaccc<-log(n3.codon.pair$ttaccc/cps.tta.ccc)
cps.obs.exp.ttaccg<-log(n3.codon.pair$ttaccg/cps.tta.ccg)
cps.obs.exp.ttacct<-log(n3.codon.pair$ttacct/cps.tta.cct)
cps.obs.exp.ttacga<-log(n3.codon.pair$ttacga/cps.tta.cga)
cps.obs.exp.ttacgc<-log(n3.codon.pair$ttacgc/cps.tta.cgc)
cps.obs.exp.ttacgg<-log(n3.codon.pair$ttacgg/cps.tta.cgg)
cps.obs.exp.ttacgt<-log(n3.codon.pair$ttacgt/cps.tta.cgt)
cps.obs.exp.ttacta<-log(n3.codon.pair$ttacta/cps.tta.cta)
cps.obs.exp.ttactc<-log(n3.codon.pair$ttactc/cps.tta.ctc)
cps.obs.exp.ttactg<-log(n3.codon.pair$ttactg/cps.tta.ctg)
cps.obs.exp.ttactt<-log(n3.codon.pair$ttactt/cps.tta.ctt)

cps.obs.exp.ttagaa<-log(n3.codon.pair$ttagaa/cps.tta.gaa)
cps.obs.exp.ttagac<-log(n3.codon.pair$ttagac/cps.tta.gac)
cps.obs.exp.ttagag<-log(n3.codon.pair$ttagag/cps.tta.gag)
cps.obs.exp.ttagat<-log(n3.codon.pair$ttagat/cps.tta.gat)
cps.obs.exp.ttagca<-log(n3.codon.pair$ttagca/cps.tta.gca)
cps.obs.exp.ttagcc<-log(n3.codon.pair$ttagcc/cps.tta.gcc)
cps.obs.exp.ttagcg<-log(n3.codon.pair$ttagcg/cps.tta.gcg)
cps.obs.exp.ttagct<-log(n3.codon.pair$ttagct/cps.tta.gct)
cps.obs.exp.ttagga<-log(n3.codon.pair$ttagga/cps.tta.gga)
cps.obs.exp.ttaggc<-log(n3.codon.pair$ttaggc/cps.tta.ggc)
cps.obs.exp.ttaggg<-log(n3.codon.pair$ttaggg/cps.tta.ggg)
cps.obs.exp.ttaggt<-log(n3.codon.pair$ttaggt/cps.tta.ggt)
cps.obs.exp.ttagta<-log(n3.codon.pair$ttagta/cps.tta.gta)
cps.obs.exp.ttagtc<-log(n3.codon.pair$ttagtc/cps.tta.gtc)
cps.obs.exp.ttagtg<-log(n3.codon.pair$ttagtg/cps.tta.gtg)
cps.obs.exp.ttagtt<-log(n3.codon.pair$ttagtt/cps.tta.gtt)

#cps.obs.exp.ttataa<-log(n3.codon.pair$ttataa/cps.tta.taa)
cps.obs.exp.ttatac<-log(n3.codon.pair$ttatac/cps.tta.tac)
#cps.obs.exp.ttatag<-log(n3.codon.pair$ttatag/cps.tta.tag)
cps.obs.exp.ttatat<-log(n3.codon.pair$ttatat/cps.tta.tat)
cps.obs.exp.ttatca<-log(n3.codon.pair$ttatca/cps.tta.tca)
cps.obs.exp.ttatcc<-log(n3.codon.pair$ttatcc/cps.tta.tcc)
cps.obs.exp.ttatcg<-log(n3.codon.pair$ttatcg/cps.tta.tcg)
cps.obs.exp.ttatct<-log(n3.codon.pair$ttatct/cps.tta.tct)
#cps.obs.exp.ttatga<-log(n3.codon.pair$ttatga/cps.tta.tga)
cps.obs.exp.ttatgc<-log(n3.codon.pair$ttatgc/cps.tta.tgc)
cps.obs.exp.ttatgg<-log(n3.codon.pair$ttatgg/cps.tta.tgg)
cps.obs.exp.ttatgt<-log(n3.codon.pair$ttatgt/cps.tta.tgt)
cps.obs.exp.ttatta<-log(n3.codon.pair$ttatta/cps.tta.tta)
cps.obs.exp.ttattc<-log(n3.codon.pair$ttattc/cps.tta.ttc)
cps.obs.exp.ttattg<-log(n3.codon.pair$ttattg/cps.tta.ttg)
cps.obs.exp.ttattt<-log(n3.codon.pair$ttattt/cps.tta.ttt)









cps.obs.exp.ttcaaa<-log(n3.codon.pair$ttcaaa/cps.ttc.aaa)
cps.obs.exp.ttcaac<-log(n3.codon.pair$ttcaac/cps.ttc.aac)
cps.obs.exp.ttcaag<-log(n3.codon.pair$ttcaag/cps.ttc.aag)
cps.obs.exp.ttcaat<-log(n3.codon.pair$ttcaat/cps.ttc.aat)
cps.obs.exp.ttcaca<-log(n3.codon.pair$ttcaca/cps.ttc.aca)
cps.obs.exp.ttcacc<-log(n3.codon.pair$ttcacc/cps.ttc.acc)
cps.obs.exp.ttcacg<-log(n3.codon.pair$ttcacg/cps.ttc.acg)
cps.obs.exp.ttcact<-log(n3.codon.pair$ttcact/cps.ttc.act)
cps.obs.exp.ttcaga<-log(n3.codon.pair$ttcaga/cps.ttc.aga)
cps.obs.exp.ttcagc<-log(n3.codon.pair$ttcagc/cps.ttc.agc)
cps.obs.exp.ttcagg<-log(n3.codon.pair$ttcagg/cps.ttc.agg)
cps.obs.exp.ttcagt<-log(n3.codon.pair$ttcagt/cps.ttc.agt)
cps.obs.exp.ttcata<-log(n3.codon.pair$ttcata/cps.ttc.ata)
cps.obs.exp.ttcatc<-log(n3.codon.pair$ttcatc/cps.ttc.atc)
cps.obs.exp.ttcatg<-log(n3.codon.pair$ttcatg/cps.ttc.atg)
cps.obs.exp.ttcatt<-log(n3.codon.pair$ttcatt/cps.ttc.att)

cps.obs.exp.ttccaa<-log(n3.codon.pair$ttccaa/cps.ttc.caa)
cps.obs.exp.ttccac<-log(n3.codon.pair$ttccac/cps.ttc.cac)
cps.obs.exp.ttccag<-log(n3.codon.pair$ttccag/cps.ttc.cag)
cps.obs.exp.ttccat<-log(n3.codon.pair$ttccat/cps.ttc.cat)
cps.obs.exp.ttccca<-log(n3.codon.pair$ttccca/cps.ttc.cca)
cps.obs.exp.ttcccc<-log(n3.codon.pair$ttcccc/cps.ttc.ccc)
cps.obs.exp.ttcccg<-log(n3.codon.pair$ttcccg/cps.ttc.ccg)
cps.obs.exp.ttccct<-log(n3.codon.pair$ttccct/cps.ttc.cct)
cps.obs.exp.ttccga<-log(n3.codon.pair$ttccga/cps.ttc.cga)
cps.obs.exp.ttccgc<-log(n3.codon.pair$ttccgc/cps.ttc.cgc)
cps.obs.exp.ttccgg<-log(n3.codon.pair$ttccgg/cps.ttc.cgg)
cps.obs.exp.ttccgt<-log(n3.codon.pair$ttccgt/cps.ttc.cgt)
cps.obs.exp.ttccta<-log(n3.codon.pair$ttccta/cps.ttc.cta)
cps.obs.exp.ttcctc<-log(n3.codon.pair$ttcctc/cps.ttc.ctc)
cps.obs.exp.ttcctg<-log(n3.codon.pair$ttcctg/cps.ttc.ctg)
cps.obs.exp.ttcctt<-log(n3.codon.pair$ttcctt/cps.ttc.ctt)

cps.obs.exp.ttcgaa<-log(n3.codon.pair$ttcgaa/cps.ttc.gaa)
cps.obs.exp.ttcgac<-log(n3.codon.pair$ttcgac/cps.ttc.gac)
cps.obs.exp.ttcgag<-log(n3.codon.pair$ttcgag/cps.ttc.gag)
cps.obs.exp.ttcgat<-log(n3.codon.pair$ttcgat/cps.ttc.gat)
cps.obs.exp.ttcgca<-log(n3.codon.pair$ttcgca/cps.ttc.gca)
cps.obs.exp.ttcgcc<-log(n3.codon.pair$ttcgcc/cps.ttc.gcc)
cps.obs.exp.ttcgcg<-log(n3.codon.pair$ttcgcg/cps.ttc.gcg)
cps.obs.exp.ttcgct<-log(n3.codon.pair$ttcgct/cps.ttc.gct)
cps.obs.exp.ttcgga<-log(n3.codon.pair$ttcgga/cps.ttc.gga)
cps.obs.exp.ttcggc<-log(n3.codon.pair$ttcggc/cps.ttc.ggc)
cps.obs.exp.ttcggg<-log(n3.codon.pair$ttcggg/cps.ttc.ggg)
cps.obs.exp.ttcggt<-log(n3.codon.pair$ttcggt/cps.ttc.ggt)
cps.obs.exp.ttcgta<-log(n3.codon.pair$ttcgta/cps.ttc.gta)
cps.obs.exp.ttcgtc<-log(n3.codon.pair$ttcgtc/cps.ttc.gtc)
cps.obs.exp.ttcgtg<-log(n3.codon.pair$ttcgtg/cps.ttc.gtg)
cps.obs.exp.ttcgtt<-log(n3.codon.pair$ttcgtt/cps.ttc.gtt)

#cps.obs.exp.ttctaa<-log(n3.codon.pair$ttctaa/cps.ttc.taa)
cps.obs.exp.ttctac<-log(n3.codon.pair$ttctac/cps.ttc.tac)
#cps.obs.exp.ttctag<-log(n3.codon.pair$ttctag/cps.ttc.tag)
cps.obs.exp.ttctat<-log(n3.codon.pair$ttctat/cps.ttc.tat)
cps.obs.exp.ttctca<-log(n3.codon.pair$ttctca/cps.ttc.tca)
cps.obs.exp.ttctcc<-log(n3.codon.pair$ttctcc/cps.ttc.tcc)
cps.obs.exp.ttctcg<-log(n3.codon.pair$ttctcg/cps.ttc.tcg)
cps.obs.exp.ttctct<-log(n3.codon.pair$ttctct/cps.ttc.tct)
#cps.obs.exp.ttctga<-log(n3.codon.pair$ttctga/cps.ttc.tga)
cps.obs.exp.ttctgc<-log(n3.codon.pair$ttctgc/cps.ttc.tgc)
cps.obs.exp.ttctgg<-log(n3.codon.pair$ttctgg/cps.ttc.tgg)
cps.obs.exp.ttctgt<-log(n3.codon.pair$ttctgt/cps.ttc.tgt)
cps.obs.exp.ttctta<-log(n3.codon.pair$ttctta/cps.ttc.tta)
cps.obs.exp.ttcttc<-log(n3.codon.pair$ttcttc/cps.ttc.ttc)
cps.obs.exp.ttcttg<-log(n3.codon.pair$ttcttg/cps.ttc.ttg)
cps.obs.exp.ttcttt<-log(n3.codon.pair$ttcttt/cps.ttc.ttt)










cps.obs.exp.ttgaaa<-log(n3.codon.pair$ttgaaa/cps.ttg.aaa)
cps.obs.exp.ttgaac<-log(n3.codon.pair$ttgaac/cps.ttg.aac)
cps.obs.exp.ttgaag<-log(n3.codon.pair$ttgaag/cps.ttg.aag)
cps.obs.exp.ttgaat<-log(n3.codon.pair$ttgaat/cps.ttg.aat)
cps.obs.exp.ttgaca<-log(n3.codon.pair$ttgaca/cps.ttg.aca)
cps.obs.exp.ttgacc<-log(n3.codon.pair$ttgacc/cps.ttg.acc)
cps.obs.exp.ttgacg<-log(n3.codon.pair$ttgacg/cps.ttg.acg)
cps.obs.exp.ttgact<-log(n3.codon.pair$ttgact/cps.ttg.act)
cps.obs.exp.ttgaga<-log(n3.codon.pair$ttgaga/cps.ttg.aga)
cps.obs.exp.ttgagc<-log(n3.codon.pair$ttgagc/cps.ttg.agc)
cps.obs.exp.ttgagg<-log(n3.codon.pair$ttgagg/cps.ttg.agg)
cps.obs.exp.ttgagt<-log(n3.codon.pair$ttgagt/cps.ttg.agt)
cps.obs.exp.ttgata<-log(n3.codon.pair$ttgata/cps.ttg.ata)
cps.obs.exp.ttgatc<-log(n3.codon.pair$ttgatc/cps.ttg.atc)
cps.obs.exp.ttgatg<-log(n3.codon.pair$ttgatg/cps.ttg.atg)
cps.obs.exp.ttgatt<-log(n3.codon.pair$ttgatt/cps.ttg.att)

cps.obs.exp.ttgcaa<-log(n3.codon.pair$ttgcaa/cps.ttg.caa)
cps.obs.exp.ttgcac<-log(n3.codon.pair$ttgcac/cps.ttg.cac)
cps.obs.exp.ttgcag<-log(n3.codon.pair$ttgcag/cps.ttg.cag)
cps.obs.exp.ttgcat<-log(n3.codon.pair$ttgcat/cps.ttg.cat)
cps.obs.exp.ttgcca<-log(n3.codon.pair$ttgcca/cps.ttg.cca)
cps.obs.exp.ttgccc<-log(n3.codon.pair$ttgccc/cps.ttg.ccc)
cps.obs.exp.ttgccg<-log(n3.codon.pair$ttgccg/cps.ttg.ccg)
cps.obs.exp.ttgcct<-log(n3.codon.pair$ttgcct/cps.ttg.cct)
cps.obs.exp.ttgcga<-log(n3.codon.pair$ttgcga/cps.ttg.cga)
cps.obs.exp.ttgcgc<-log(n3.codon.pair$ttgcgc/cps.ttg.cgc)
cps.obs.exp.ttgcgg<-log(n3.codon.pair$ttgcgg/cps.ttg.cgg)
cps.obs.exp.ttgcgt<-log(n3.codon.pair$ttgcgt/cps.ttg.cgt)
cps.obs.exp.ttgcta<-log(n3.codon.pair$ttgcta/cps.ttg.cta)
cps.obs.exp.ttgctc<-log(n3.codon.pair$ttgctc/cps.ttg.ctc)
cps.obs.exp.ttgctg<-log(n3.codon.pair$ttgctg/cps.ttg.ctg)
cps.obs.exp.ttgctt<-log(n3.codon.pair$ttgctt/cps.ttg.ctt)

cps.obs.exp.ttggaa<-log(n3.codon.pair$ttggaa/cps.ttg.gaa)
cps.obs.exp.ttggac<-log(n3.codon.pair$ttggac/cps.ttg.gac)
cps.obs.exp.ttggag<-log(n3.codon.pair$ttggag/cps.ttg.gag)
cps.obs.exp.ttggat<-log(n3.codon.pair$ttggat/cps.ttg.gat)
cps.obs.exp.ttggca<-log(n3.codon.pair$ttggca/cps.ttg.gca)
cps.obs.exp.ttggcc<-log(n3.codon.pair$ttggcc/cps.ttg.gcc)
cps.obs.exp.ttggcg<-log(n3.codon.pair$ttggcg/cps.ttg.gcg)
cps.obs.exp.ttggct<-log(n3.codon.pair$ttggct/cps.ttg.gct)
cps.obs.exp.ttggga<-log(n3.codon.pair$ttggga/cps.ttg.gga)
cps.obs.exp.ttgggc<-log(n3.codon.pair$ttgggc/cps.ttg.ggc)
cps.obs.exp.ttgggg<-log(n3.codon.pair$ttgggg/cps.ttg.ggg)
cps.obs.exp.ttgggt<-log(n3.codon.pair$ttgggt/cps.ttg.ggt)
cps.obs.exp.ttggta<-log(n3.codon.pair$ttggta/cps.ttg.gta)
cps.obs.exp.ttggtc<-log(n3.codon.pair$ttggtc/cps.ttg.gtc)
cps.obs.exp.ttggtg<-log(n3.codon.pair$ttggtg/cps.ttg.gtg)
cps.obs.exp.ttggtt<-log(n3.codon.pair$ttggtt/cps.ttg.gtt)

#cps.obs.exp.ttgtaa<-log(n3.codon.pair$ttgtaa/cps.ttg.taa)
cps.obs.exp.ttgtac<-log(n3.codon.pair$ttgtac/cps.ttg.tac)
#cps.obs.exp.ttgtag<-log(n3.codon.pair$ttgtag/cps.ttg.tag)
cps.obs.exp.ttgtat<-log(n3.codon.pair$ttgtat/cps.ttg.tat)
cps.obs.exp.ttgtca<-log(n3.codon.pair$ttgtca/cps.ttg.tca)
cps.obs.exp.ttgtcc<-log(n3.codon.pair$ttgtcc/cps.ttg.tcc)
cps.obs.exp.ttgtcg<-log(n3.codon.pair$ttgtcg/cps.ttg.tcg)
cps.obs.exp.ttgtct<-log(n3.codon.pair$ttgtct/cps.ttg.tct)
#cps.obs.exp.ttgtga<-log(n3.codon.pair$ttgtga/cps.ttg.tga)
cps.obs.exp.ttgtgc<-log(n3.codon.pair$ttgtgc/cps.ttg.tgc)
cps.obs.exp.ttgtgg<-log(n3.codon.pair$ttgtgg/cps.ttg.tgg)
cps.obs.exp.ttgtgt<-log(n3.codon.pair$ttgtgt/cps.ttg.tgt)
cps.obs.exp.ttgtta<-log(n3.codon.pair$ttgtta/cps.ttg.tta)
cps.obs.exp.ttgttc<-log(n3.codon.pair$ttgttc/cps.ttg.ttc)
cps.obs.exp.ttgttg<-log(n3.codon.pair$ttgttg/cps.ttg.ttg)
cps.obs.exp.ttgttt<-log(n3.codon.pair$ttgttt/cps.ttg.ttt)









cps.obs.exp.tttaaa<-log(n3.codon.pair$tttaaa/cps.ttt.aaa)
cps.obs.exp.tttaac<-log(n3.codon.pair$tttaac/cps.ttt.aac)
cps.obs.exp.tttaag<-log(n3.codon.pair$tttaag/cps.ttt.aag)
cps.obs.exp.tttaat<-log(n3.codon.pair$tttaat/cps.ttt.aat)
cps.obs.exp.tttaca<-log(n3.codon.pair$tttaca/cps.ttt.aca)
cps.obs.exp.tttacc<-log(n3.codon.pair$tttacc/cps.ttt.acc)
cps.obs.exp.tttacg<-log(n3.codon.pair$tttacg/cps.ttt.acg)
cps.obs.exp.tttact<-log(n3.codon.pair$tttact/cps.ttt.act)
cps.obs.exp.tttaga<-log(n3.codon.pair$tttaga/cps.ttt.aga)
cps.obs.exp.tttagc<-log(n3.codon.pair$tttagc/cps.ttt.agc)
cps.obs.exp.tttagg<-log(n3.codon.pair$tttagg/cps.ttt.agg)
cps.obs.exp.tttagt<-log(n3.codon.pair$tttagt/cps.ttt.agt)
cps.obs.exp.tttata<-log(n3.codon.pair$tttata/cps.ttt.ata)
cps.obs.exp.tttatc<-log(n3.codon.pair$tttatc/cps.ttt.atc)
cps.obs.exp.tttatg<-log(n3.codon.pair$tttatg/cps.ttt.atg)
cps.obs.exp.tttatt<-log(n3.codon.pair$tttatt/cps.ttt.att)

cps.obs.exp.tttcaa<-log(n3.codon.pair$tttcaa/cps.ttt.caa)
cps.obs.exp.tttcac<-log(n3.codon.pair$tttcac/cps.ttt.cac)
cps.obs.exp.tttcag<-log(n3.codon.pair$tttcag/cps.ttt.cag)
cps.obs.exp.tttcat<-log(n3.codon.pair$tttcat/cps.ttt.cat)
cps.obs.exp.tttcca<-log(n3.codon.pair$tttcca/cps.ttt.cca)
cps.obs.exp.tttccc<-log(n3.codon.pair$tttccc/cps.ttt.ccc)
cps.obs.exp.tttccg<-log(n3.codon.pair$tttccg/cps.ttt.ccg)
cps.obs.exp.tttcct<-log(n3.codon.pair$tttcct/cps.ttt.cct)
cps.obs.exp.tttcga<-log(n3.codon.pair$tttcga/cps.ttt.cga)
cps.obs.exp.tttcgc<-log(n3.codon.pair$tttcgc/cps.ttt.cgc)
cps.obs.exp.tttcgg<-log(n3.codon.pair$tttcgg/cps.ttt.cgg)
cps.obs.exp.tttcgt<-log(n3.codon.pair$tttcgt/cps.ttt.cgt)
cps.obs.exp.tttcta<-log(n3.codon.pair$tttcta/cps.ttt.cta)
cps.obs.exp.tttctc<-log(n3.codon.pair$tttctc/cps.ttt.ctc)
cps.obs.exp.tttctg<-log(n3.codon.pair$tttctg/cps.ttt.ctg)
cps.obs.exp.tttctt<-log(n3.codon.pair$tttctt/cps.ttt.ctt)

cps.obs.exp.tttgaa<-log(n3.codon.pair$tttgaa/cps.ttt.gaa)
cps.obs.exp.tttgac<-log(n3.codon.pair$tttgac/cps.ttt.gac)
cps.obs.exp.tttgag<-log(n3.codon.pair$tttgag/cps.ttt.gag)
cps.obs.exp.tttgat<-log(n3.codon.pair$tttgat/cps.ttt.gat)
cps.obs.exp.tttgca<-log(n3.codon.pair$tttgca/cps.ttt.gca)
cps.obs.exp.tttgcc<-log(n3.codon.pair$tttgcc/cps.ttt.gcc)
cps.obs.exp.tttgcg<-log(n3.codon.pair$tttgcg/cps.ttt.gcg)
cps.obs.exp.tttgct<-log(n3.codon.pair$tttgct/cps.ttt.gct)
cps.obs.exp.tttgga<-log(n3.codon.pair$tttgga/cps.ttt.gga)
cps.obs.exp.tttggc<-log(n3.codon.pair$tttggc/cps.ttt.ggc)
cps.obs.exp.tttggg<-log(n3.codon.pair$tttggg/cps.ttt.ggg)
cps.obs.exp.tttggt<-log(n3.codon.pair$tttggt/cps.ttt.ggt)
cps.obs.exp.tttgta<-log(n3.codon.pair$tttgta/cps.ttt.gta)
cps.obs.exp.tttgtc<-log(n3.codon.pair$tttgtc/cps.ttt.gtc)
cps.obs.exp.tttgtg<-log(n3.codon.pair$tttgtg/cps.ttt.gtg)
cps.obs.exp.tttgtt<-log(n3.codon.pair$tttgtt/cps.ttt.gtt)

#cps.obs.exp.ttttaa<-log(n3.codon.pair$ttttaa/cps.ttt.taa)
cps.obs.exp.ttttac<-log(n3.codon.pair$ttttac/cps.ttt.tac)
#cps.obs.exp.ttttag<-log(n3.codon.pair$ttttag/cps.ttt.tag)
cps.obs.exp.ttttat<-log(n3.codon.pair$ttttat/cps.ttt.tat)
cps.obs.exp.ttttca<-log(n3.codon.pair$ttttca/cps.ttt.tca)
cps.obs.exp.ttttcc<-log(n3.codon.pair$ttttcc/cps.ttt.tcc)
cps.obs.exp.ttttcg<-log(n3.codon.pair$ttttcg/cps.ttt.tcg)
cps.obs.exp.ttttct<-log(n3.codon.pair$ttttct/cps.ttt.tct)
#cps.obs.exp.ttttga<-log(n3.codon.pair$ttttga/cps.ttt.tga)
cps.obs.exp.ttttgc<-log(n3.codon.pair$ttttgc/cps.ttt.tgc)
cps.obs.exp.ttttgg<-log(n3.codon.pair$ttttgg/cps.ttt.tgg)
cps.obs.exp.ttttgt<-log(n3.codon.pair$ttttgt/cps.ttt.tgt)
cps.obs.exp.ttttta<-log(n3.codon.pair$ttttta/cps.ttt.tta)
cps.obs.exp.tttttc<-log(n3.codon.pair$tttttc/cps.ttt.ttc)
cps.obs.exp.tttttg<-log(n3.codon.pair$tttttg/cps.ttt.ttg)
cps.obs.exp.tttttt<-log(n3.codon.pair$tttttt/cps.ttt.ttt)





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
#Original with stop codon  cpb.n3<-cps.obs.exp.aaaaaa+cps.obs.exp.aaaaac+cps.obs.exp.aaaaag+cps.obs.exp.aaaaat+cps.obs.exp.aaaaca+cps.obs.exp.aaaacc+cps.obs.exp.aaaacg+cps.obs.exp.aaaact+cps.obs.exp.aaaaga+cps.obs.exp.aaaagc+cps.obs.exp.aaaagg+cps.obs.exp.aaaagt+cps.obs.exp.aaaata+cps.obs.exp.aaaatc+cps.obs.exp.aaaatg+cps.obs.exp.aaaatt+cps.obs.exp.aaacaa+cps.obs.exp.aaacac+cps.obs.exp.aaacag+cps.obs.exp.aaacat+cps.obs.exp.aaacca+cps.obs.exp.aaaccc+cps.obs.exp.aaaccg+cps.obs.exp.aaacct+cps.obs.exp.aaacga+cps.obs.exp.aaacgc+cps.obs.exp.aaacgg+cps.obs.exp.aaacgt+cps.obs.exp.aaacta+cps.obs.exp.aaactc+cps.obs.exp.aaactg+cps.obs.exp.aaactt+cps.obs.exp.aaagaa+cps.obs.exp.aaagac+cps.obs.exp.aaagag+cps.obs.exp.aaagat+cps.obs.exp.aaagca+cps.obs.exp.aaagcc+cps.obs.exp.aaagcg+cps.obs.exp.aaagct+cps.obs.exp.aaagga+cps.obs.exp.aaaggc+cps.obs.exp.aaaggg+cps.obs.exp.aaaggt+cps.obs.exp.aaagta+cps.obs.exp.aaagtc+cps.obs.exp.aaagtg+cps.obs.exp.aaagtt+cps.obs.exp.aaataa+cps.obs.exp.aaatac+cps.obs.exp.aaatag+cps.obs.exp.aaatat+cps.obs.exp.aaatca+cps.obs.exp.aaatcc+cps.obs.exp.aaatcg+cps.obs.exp.aaatct+cps.obs.exp.aaatga+cps.obs.exp.aaatgc+cps.obs.exp.aaatgg+cps.obs.exp.aaatgt+cps.obs.exp.aaatta+cps.obs.exp.aaattc+cps.obs.exp.aaattg+cps.obs.exp.aaattt+cps.obs.exp.aacaaa+cps.obs.exp.aacaac+cps.obs.exp.aacaag+cps.obs.exp.aacaat+cps.obs.exp.aacaca+cps.obs.exp.aacacc+cps.obs.exp.aacacg+cps.obs.exp.aacact+cps.obs.exp.aacaga+cps.obs.exp.aacagc+cps.obs.exp.aacagg+cps.obs.exp.aacagt+cps.obs.exp.aacata+cps.obs.exp.aacatc+cps.obs.exp.aacatg+cps.obs.exp.aacatt+cps.obs.exp.aaccaa+cps.obs.exp.aaccac+cps.obs.exp.aaccag+cps.obs.exp.aaccat+cps.obs.exp.aaccca+cps.obs.exp.aacccc+cps.obs.exp.aacccg+cps.obs.exp.aaccct+cps.obs.exp.aaccga+cps.obs.exp.aaccgc+cps.obs.exp.aaccgg+cps.obs.exp.aaccgt+cps.obs.exp.aaccta+cps.obs.exp.aacctc+cps.obs.exp.aacctg+cps.obs.exp.aacctt+cps.obs.exp.aacgaa+cps.obs.exp.aacgac+cps.obs.exp.aacgag+cps.obs.exp.aacgat+cps.obs.exp.aacgca+cps.obs.exp.aacgcc+cps.obs.exp.aacgcg+cps.obs.exp.aacgct+cps.obs.exp.aacgga+cps.obs.exp.aacggc+cps.obs.exp.aacggg+cps.obs.exp.aacggt+cps.obs.exp.aacgta+cps.obs.exp.aacgtc+cps.obs.exp.aacgtg+cps.obs.exp.aacgtt+cps.obs.exp.aactaa+cps.obs.exp.aactac+cps.obs.exp.aactag+cps.obs.exp.aactat+cps.obs.exp.aactca+cps.obs.exp.aactcc+cps.obs.exp.aactcg+cps.obs.exp.aactct+cps.obs.exp.aactga+cps.obs.exp.aactgc+cps.obs.exp.aactgg+cps.obs.exp.aactgt+cps.obs.exp.aactta+cps.obs.exp.aacttc+cps.obs.exp.aacttg+cps.obs.exp.aacttt+cps.obs.exp.aagaaa+cps.obs.exp.aagaac+cps.obs.exp.aagaag+cps.obs.exp.aagaat+cps.obs.exp.aagaca+cps.obs.exp.aagacc+cps.obs.exp.aagacg+cps.obs.exp.aagact+cps.obs.exp.aagaga+cps.obs.exp.aagagc+cps.obs.exp.aagagg+cps.obs.exp.aagagt+cps.obs.exp.aagata+cps.obs.exp.aagatc+cps.obs.exp.aagatg+cps.obs.exp.aagatt+cps.obs.exp.aagcaa+cps.obs.exp.aagcac+cps.obs.exp.aagcag+cps.obs.exp.aagcat+cps.obs.exp.aagcca+cps.obs.exp.aagccc+cps.obs.exp.aagccg+cps.obs.exp.aagcct+cps.obs.exp.aagcga+cps.obs.exp.aagcgc+cps.obs.exp.aagcgg+cps.obs.exp.aagcgt+cps.obs.exp.aagcta+cps.obs.exp.aagctc+cps.obs.exp.aagctg+cps.obs.exp.aagctt+cps.obs.exp.aaggaa+cps.obs.exp.aaggac+cps.obs.exp.aaggag+cps.obs.exp.aaggat+cps.obs.exp.aaggca+cps.obs.exp.aaggcc+cps.obs.exp.aaggcg+cps.obs.exp.aaggct+cps.obs.exp.aaggga+cps.obs.exp.aagggc+cps.obs.exp.aagggg+cps.obs.exp.aagggt+cps.obs.exp.aaggta+cps.obs.exp.aaggtc+cps.obs.exp.aaggtg+cps.obs.exp.aaggtt+cps.obs.exp.aagtaa+cps.obs.exp.aagtac+cps.obs.exp.aagtag+cps.obs.exp.aagtat+cps.obs.exp.aagtca+cps.obs.exp.aagtcc+cps.obs.exp.aagtcg+cps.obs.exp.aagtct+cps.obs.exp.aagtga+cps.obs.exp.aagtgc+cps.obs.exp.aagtgg+cps.obs.exp.aagtgt+cps.obs.exp.aagtta+cps.obs.exp.aagttc+cps.obs.exp.aagttg+cps.obs.exp.aagttt+cps.obs.exp.aataaa+cps.obs.exp.aataac+cps.obs.exp.aataag+cps.obs.exp.aataat+cps.obs.exp.aataca+cps.obs.exp.aatacc+cps.obs.exp.aatacg+cps.obs.exp.aatact+cps.obs.exp.aataga+cps.obs.exp.aatagc+cps.obs.exp.aatagg+cps.obs.exp.aatagt+cps.obs.exp.aatata+cps.obs.exp.aatatc+cps.obs.exp.aatatg+cps.obs.exp.aatatt+cps.obs.exp.aatcaa+cps.obs.exp.aatcac+cps.obs.exp.aatcag+cps.obs.exp.aatcat+cps.obs.exp.aatcca+cps.obs.exp.aatccc+cps.obs.exp.aatccg+cps.obs.exp.aatcct+cps.obs.exp.aatcga+cps.obs.exp.aatcgc+cps.obs.exp.aatcgg+cps.obs.exp.aatcgt+cps.obs.exp.aatcta+cps.obs.exp.aatctc+cps.obs.exp.aatctg+cps.obs.exp.aatctt+cps.obs.exp.aatgaa+cps.obs.exp.aatgac+cps.obs.exp.aatgag+cps.obs.exp.aatgat+cps.obs.exp.aatgca+cps.obs.exp.aatgcc+cps.obs.exp.aatgcg+cps.obs.exp.aatgct+cps.obs.exp.aatgga+cps.obs.exp.aatggc+cps.obs.exp.aatggg+cps.obs.exp.aatggt+cps.obs.exp.aatgta+cps.obs.exp.aatgtc+cps.obs.exp.aatgtg+cps.obs.exp.aatgtt+cps.obs.exp.aattaa+cps.obs.exp.aattac+cps.obs.exp.aattag+cps.obs.exp.aattat+cps.obs.exp.aattca+cps.obs.exp.aattcc+cps.obs.exp.aattcg+cps.obs.exp.aattct+cps.obs.exp.aattga+cps.obs.exp.aattgc+cps.obs.exp.aattgg+cps.obs.exp.aattgt+cps.obs.exp.aattta+cps.obs.exp.aatttc+cps.obs.exp.aatttg+cps.obs.exp.aatttt+cps.obs.exp.acaaaa+cps.obs.exp.acaaac+cps.obs.exp.acaaag+cps.obs.exp.acaaat+cps.obs.exp.acaaca+cps.obs.exp.acaacc+cps.obs.exp.acaacg+cps.obs.exp.acaact+cps.obs.exp.acaaga+cps.obs.exp.acaagc+cps.obs.exp.acaagg+cps.obs.exp.acaagt+cps.obs.exp.acaata+cps.obs.exp.acaatc+cps.obs.exp.acaatg+cps.obs.exp.acaatt+cps.obs.exp.acacaa+cps.obs.exp.acacac+cps.obs.exp.acacag+cps.obs.exp.acacat+cps.obs.exp.acacca+cps.obs.exp.acaccc+cps.obs.exp.acaccg+cps.obs.exp.acacct+cps.obs.exp.acacga+cps.obs.exp.acacgc+cps.obs.exp.acacgg+cps.obs.exp.acacgt+cps.obs.exp.acacta+cps.obs.exp.acactc+cps.obs.exp.acactg+cps.obs.exp.acactt+cps.obs.exp.acagaa+cps.obs.exp.acagac+cps.obs.exp.acagag+cps.obs.exp.acagat+cps.obs.exp.acagca+cps.obs.exp.acagcc+cps.obs.exp.acagcg+cps.obs.exp.acagct+cps.obs.exp.acagga+cps.obs.exp.acaggc+cps.obs.exp.acaggg+cps.obs.exp.acaggt+cps.obs.exp.acagta+cps.obs.exp.acagtc+cps.obs.exp.acagtg+cps.obs.exp.acagtt+cps.obs.exp.acataa+cps.obs.exp.acatac+cps.obs.exp.acatag+cps.obs.exp.acatat+cps.obs.exp.acatca+cps.obs.exp.acatcc+cps.obs.exp.acatcg+cps.obs.exp.acatct+cps.obs.exp.acatga+cps.obs.exp.acatgc+cps.obs.exp.acatgg+cps.obs.exp.acatgt+cps.obs.exp.acatta+cps.obs.exp.acattc+cps.obs.exp.acattg+cps.obs.exp.acattt+cps.obs.exp.accaaa+cps.obs.exp.accaac+cps.obs.exp.accaag+cps.obs.exp.accaat+cps.obs.exp.accaca+cps.obs.exp.accacc+cps.obs.exp.accacg+cps.obs.exp.accact+cps.obs.exp.accaga+cps.obs.exp.accagc+cps.obs.exp.accagg+cps.obs.exp.accagt+cps.obs.exp.accata+cps.obs.exp.accatc+cps.obs.exp.accatg+cps.obs.exp.accatt+cps.obs.exp.acccaa+cps.obs.exp.acccac+cps.obs.exp.acccag+cps.obs.exp.acccat+cps.obs.exp.acccca+cps.obs.exp.accccc+cps.obs.exp.accccg+cps.obs.exp.acccct+cps.obs.exp.acccga+cps.obs.exp.acccgc+cps.obs.exp.acccgg+cps.obs.exp.acccgt+cps.obs.exp.acccta+cps.obs.exp.accctc+cps.obs.exp.accctg+cps.obs.exp.accctt+cps.obs.exp.accgaa+cps.obs.exp.accgac+cps.obs.exp.accgag+cps.obs.exp.accgat+cps.obs.exp.accgca+cps.obs.exp.accgcc+cps.obs.exp.accgcg+cps.obs.exp.accgct+cps.obs.exp.accgga+cps.obs.exp.accggc+cps.obs.exp.accggg+cps.obs.exp.accggt+cps.obs.exp.accgta+cps.obs.exp.accgtc+cps.obs.exp.accgtg+cps.obs.exp.accgtt+cps.obs.exp.acctaa+cps.obs.exp.acctac+cps.obs.exp.acctag+cps.obs.exp.acctat+cps.obs.exp.acctca+cps.obs.exp.acctcc+cps.obs.exp.acctcg+cps.obs.exp.acctct+cps.obs.exp.acctga+cps.obs.exp.acctgc+cps.obs.exp.acctgg+cps.obs.exp.acctgt+cps.obs.exp.acctta+cps.obs.exp.accttc+cps.obs.exp.accttg+cps.obs.exp.accttt+cps.obs.exp.acgaaa+cps.obs.exp.acgaac+cps.obs.exp.acgaag+cps.obs.exp.acgaat+cps.obs.exp.acgaca+cps.obs.exp.acgacc+cps.obs.exp.acgacg+cps.obs.exp.acgact+cps.obs.exp.acgaga+cps.obs.exp.acgagc+cps.obs.exp.acgagg+cps.obs.exp.acgagt+cps.obs.exp.acgata+cps.obs.exp.acgatc+cps.obs.exp.acgatg+cps.obs.exp.acgatt+cps.obs.exp.acgcaa+cps.obs.exp.acgcac+cps.obs.exp.acgcag+cps.obs.exp.acgcat+cps.obs.exp.acgcca+cps.obs.exp.acgccc+cps.obs.exp.acgccg+cps.obs.exp.acgcct+cps.obs.exp.acgcga+cps.obs.exp.acgcgc+cps.obs.exp.acgcgg+cps.obs.exp.acgcgt+cps.obs.exp.acgcta+cps.obs.exp.acgctc+cps.obs.exp.acgctg+cps.obs.exp.acgctt+cps.obs.exp.acggaa+cps.obs.exp.acggac+cps.obs.exp.acggag+cps.obs.exp.acggat+cps.obs.exp.acggca+cps.obs.exp.acggcc+cps.obs.exp.acggcg+cps.obs.exp.acggct+cps.obs.exp.acggga+cps.obs.exp.acgggc+cps.obs.exp.acgggg+cps.obs.exp.acgggt+cps.obs.exp.acggta+cps.obs.exp.acggtc+cps.obs.exp.acggtg+cps.obs.exp.acggtt+cps.obs.exp.acgtaa+cps.obs.exp.acgtac+cps.obs.exp.acgtag+cps.obs.exp.acgtat+cps.obs.exp.acgtca+cps.obs.exp.acgtcc+cps.obs.exp.acgtcg+cps.obs.exp.acgtct+cps.obs.exp.acgtga+cps.obs.exp.acgtgc+cps.obs.exp.acgtgg+cps.obs.exp.acgtgt+cps.obs.exp.acgtta+cps.obs.exp.acgttc+cps.obs.exp.acgttg+cps.obs.exp.acgttt+cps.obs.exp.actaaa+cps.obs.exp.actaac+cps.obs.exp.actaag+cps.obs.exp.actaat+cps.obs.exp.actaca+cps.obs.exp.actacc+cps.obs.exp.actacg+cps.obs.exp.actact+cps.obs.exp.actaga+cps.obs.exp.actagc+cps.obs.exp.actagg+cps.obs.exp.actagt+cps.obs.exp.actata+cps.obs.exp.actatc+cps.obs.exp.actatg+cps.obs.exp.actatt+cps.obs.exp.actcaa+cps.obs.exp.actcac+cps.obs.exp.actcag+cps.obs.exp.actcat+cps.obs.exp.actcca+cps.obs.exp.actccc+cps.obs.exp.actccg+cps.obs.exp.actcct+cps.obs.exp.actcga+cps.obs.exp.actcgc+cps.obs.exp.actcgg+cps.obs.exp.actcgt+cps.obs.exp.actcta+cps.obs.exp.actctc+cps.obs.exp.actctg+cps.obs.exp.actctt+cps.obs.exp.actgaa+cps.obs.exp.actgac+cps.obs.exp.actgag+cps.obs.exp.actgat+cps.obs.exp.actgca+cps.obs.exp.actgcc+cps.obs.exp.actgcg+cps.obs.exp.actgct+cps.obs.exp.actgga+cps.obs.exp.actggc+cps.obs.exp.actggg+cps.obs.exp.actggt+cps.obs.exp.actgta+cps.obs.exp.actgtc+cps.obs.exp.actgtg+cps.obs.exp.actgtt+cps.obs.exp.acttaa+cps.obs.exp.acttac+cps.obs.exp.acttag+cps.obs.exp.acttat+cps.obs.exp.acttca+cps.obs.exp.acttcc+cps.obs.exp.acttcg+cps.obs.exp.acttct+cps.obs.exp.acttga+cps.obs.exp.acttgc+cps.obs.exp.acttgg+cps.obs.exp.acttgt+cps.obs.exp.acttta+cps.obs.exp.actttc+cps.obs.exp.actttg+cps.obs.exp.actttt+cps.obs.exp.agaaaa+cps.obs.exp.agaaac+cps.obs.exp.agaaag+cps.obs.exp.agaaat+cps.obs.exp.agaaca+cps.obs.exp.agaacc+cps.obs.exp.agaacg+cps.obs.exp.agaact+cps.obs.exp.agaaga+cps.obs.exp.agaagc+cps.obs.exp.agaagg+cps.obs.exp.agaagt+cps.obs.exp.agaata+cps.obs.exp.agaatc+cps.obs.exp.agaatg+cps.obs.exp.agaatt+cps.obs.exp.agacaa+cps.obs.exp.agacac+cps.obs.exp.agacag+cps.obs.exp.agacat+cps.obs.exp.agacca+cps.obs.exp.agaccc+cps.obs.exp.agaccg+cps.obs.exp.agacct+cps.obs.exp.agacga+cps.obs.exp.agacgc+cps.obs.exp.agacgg+cps.obs.exp.agacgt+cps.obs.exp.agacta+cps.obs.exp.agactc+cps.obs.exp.agactg+cps.obs.exp.agactt+cps.obs.exp.agagaa+cps.obs.exp.agagac+cps.obs.exp.agagag+cps.obs.exp.agagat+cps.obs.exp.agagca+cps.obs.exp.agagcc+cps.obs.exp.agagcg+cps.obs.exp.agagct+cps.obs.exp.agagga+cps.obs.exp.agaggc+cps.obs.exp.agaggg+cps.obs.exp.agaggt+cps.obs.exp.agagta+cps.obs.exp.agagtc+cps.obs.exp.agagtg+cps.obs.exp.agagtt+cps.obs.exp.agataa+cps.obs.exp.agatac+cps.obs.exp.agatag+cps.obs.exp.agatat+cps.obs.exp.agatca+cps.obs.exp.agatcc+cps.obs.exp.agatcg+cps.obs.exp.agatct+cps.obs.exp.agatga+cps.obs.exp.agatgc+cps.obs.exp.agatgg+cps.obs.exp.agatgt+cps.obs.exp.agatta+cps.obs.exp.agattc+cps.obs.exp.agattg+cps.obs.exp.agattt+cps.obs.exp.agcaaa+cps.obs.exp.agcaac+cps.obs.exp.agcaag+cps.obs.exp.agcaat+cps.obs.exp.agcaca+cps.obs.exp.agcacc+cps.obs.exp.agcacg+cps.obs.exp.agcact+cps.obs.exp.agcaga+cps.obs.exp.agcagc+cps.obs.exp.agcagg+cps.obs.exp.agcagt+cps.obs.exp.agcata+cps.obs.exp.agcatc+cps.obs.exp.agcatg+cps.obs.exp.agcatt+cps.obs.exp.agccaa+cps.obs.exp.agccac+cps.obs.exp.agccag+cps.obs.exp.agccat+cps.obs.exp.agccca+cps.obs.exp.agcccc+cps.obs.exp.agcccg+cps.obs.exp.agccct+cps.obs.exp.agccga+cps.obs.exp.agccgc+cps.obs.exp.agccgg+cps.obs.exp.agccgt+cps.obs.exp.agccta+cps.obs.exp.agcctc+cps.obs.exp.agcctg+cps.obs.exp.agcctt+cps.obs.exp.agcgaa+cps.obs.exp.agcgac+cps.obs.exp.agcgag+cps.obs.exp.agcgat+cps.obs.exp.agcgca+cps.obs.exp.agcgcc+cps.obs.exp.agcgcg+cps.obs.exp.agcgct+cps.obs.exp.agcgga+cps.obs.exp.agcggc+cps.obs.exp.agcggg+cps.obs.exp.agcggt+cps.obs.exp.agcgta+cps.obs.exp.agcgtc+cps.obs.exp.agcgtg+cps.obs.exp.agcgtt+cps.obs.exp.agctaa+cps.obs.exp.agctac+cps.obs.exp.agctag+cps.obs.exp.agctat+cps.obs.exp.agctca+cps.obs.exp.agctcc+cps.obs.exp.agctcg+cps.obs.exp.agctct+cps.obs.exp.agctga+cps.obs.exp.agctgc+cps.obs.exp.agctgg+cps.obs.exp.agctgt+cps.obs.exp.agctta+cps.obs.exp.agcttc+cps.obs.exp.agcttg+cps.obs.exp.agcttt+cps.obs.exp.aggaaa+cps.obs.exp.aggaac+cps.obs.exp.aggaag+cps.obs.exp.aggaat+cps.obs.exp.aggaca+cps.obs.exp.aggacc+cps.obs.exp.aggacg+cps.obs.exp.aggact+cps.obs.exp.aggaga+cps.obs.exp.aggagc+cps.obs.exp.aggagg+cps.obs.exp.aggagt+cps.obs.exp.aggata+cps.obs.exp.aggatc+cps.obs.exp.aggatg+cps.obs.exp.aggatt+cps.obs.exp.aggcaa+cps.obs.exp.aggcac+cps.obs.exp.aggcag+cps.obs.exp.aggcat+cps.obs.exp.aggcca+cps.obs.exp.aggccc+cps.obs.exp.aggccg+cps.obs.exp.aggcct+cps.obs.exp.aggcga+cps.obs.exp.aggcgc+cps.obs.exp.aggcgg+cps.obs.exp.aggcgt+cps.obs.exp.aggcta+cps.obs.exp.aggctc+cps.obs.exp.aggctg+cps.obs.exp.aggctt+cps.obs.exp.agggaa+cps.obs.exp.agggac+cps.obs.exp.agggag+cps.obs.exp.agggat+cps.obs.exp.agggca+cps.obs.exp.agggcc+cps.obs.exp.agggcg+cps.obs.exp.agggct+cps.obs.exp.agggga+cps.obs.exp.aggggc+cps.obs.exp.aggggg+cps.obs.exp.aggggt+cps.obs.exp.agggta+cps.obs.exp.agggtc+cps.obs.exp.agggtg+cps.obs.exp.agggtt+cps.obs.exp.aggtaa+cps.obs.exp.aggtac+cps.obs.exp.aggtag+cps.obs.exp.aggtat+cps.obs.exp.aggtca+cps.obs.exp.aggtcc+cps.obs.exp.aggtcg+cps.obs.exp.aggtct+cps.obs.exp.aggtga+cps.obs.exp.aggtgc+cps.obs.exp.aggtgg+cps.obs.exp.aggtgt+cps.obs.exp.aggtta+cps.obs.exp.aggttc+cps.obs.exp.aggttg+cps.obs.exp.aggttt+cps.obs.exp.agtaaa+cps.obs.exp.agtaac+cps.obs.exp.agtaag+cps.obs.exp.agtaat+cps.obs.exp.agtaca+cps.obs.exp.agtacc+cps.obs.exp.agtacg+cps.obs.exp.agtact+cps.obs.exp.agtaga+cps.obs.exp.agtagc+cps.obs.exp.agtagg+cps.obs.exp.agtagt+cps.obs.exp.agtata+cps.obs.exp.agtatc+cps.obs.exp.agtatg+cps.obs.exp.agtatt+cps.obs.exp.agtcaa+cps.obs.exp.agtcac+cps.obs.exp.agtcag+cps.obs.exp.agtcat+cps.obs.exp.agtcca+cps.obs.exp.agtccc+cps.obs.exp.agtccg+cps.obs.exp.agtcct+cps.obs.exp.agtcga+cps.obs.exp.agtcgc+cps.obs.exp.agtcgg+cps.obs.exp.agtcgt+cps.obs.exp.agtcta+cps.obs.exp.agtctc+cps.obs.exp.agtctg+cps.obs.exp.agtctt+cps.obs.exp.agtgaa+cps.obs.exp.agtgac+cps.obs.exp.agtgag+cps.obs.exp.agtgat+cps.obs.exp.agtgca+cps.obs.exp.agtgcc+cps.obs.exp.agtgcg+cps.obs.exp.agtgct+cps.obs.exp.agtgga+cps.obs.exp.agtggc+cps.obs.exp.agtggg+cps.obs.exp.agtggt+cps.obs.exp.agtgta+cps.obs.exp.agtgtc+cps.obs.exp.agtgtg+cps.obs.exp.agtgtt+cps.obs.exp.agttaa+cps.obs.exp.agttac+cps.obs.exp.agttag+cps.obs.exp.agttat+cps.obs.exp.agttca+cps.obs.exp.agttcc+cps.obs.exp.agttcg+cps.obs.exp.agttct+cps.obs.exp.agttga+cps.obs.exp.agttgc+cps.obs.exp.agttgg+cps.obs.exp.agttgt+cps.obs.exp.agttta+cps.obs.exp.agtttc+cps.obs.exp.agtttg+cps.obs.exp.agtttt+cps.obs.exp.ataaaa+cps.obs.exp.ataaac+cps.obs.exp.ataaag+cps.obs.exp.ataaat+cps.obs.exp.ataaca+cps.obs.exp.ataacc+cps.obs.exp.ataacg+cps.obs.exp.ataact+cps.obs.exp.ataaga+cps.obs.exp.ataagc+cps.obs.exp.ataagg+cps.obs.exp.ataagt+cps.obs.exp.ataata+cps.obs.exp.ataatc+cps.obs.exp.ataatg+cps.obs.exp.ataatt+cps.obs.exp.atacaa+cps.obs.exp.atacac+cps.obs.exp.atacag+cps.obs.exp.atacat+cps.obs.exp.atacca+cps.obs.exp.ataccc+cps.obs.exp.ataccg+cps.obs.exp.atacct+cps.obs.exp.atacga+cps.obs.exp.atacgc+cps.obs.exp.atacgg+cps.obs.exp.atacgt+cps.obs.exp.atacta+cps.obs.exp.atactc+cps.obs.exp.atactg+cps.obs.exp.atactt+cps.obs.exp.atagaa+cps.obs.exp.atagac+cps.obs.exp.atagag+cps.obs.exp.atagat+cps.obs.exp.atagca+cps.obs.exp.atagcc+cps.obs.exp.atagcg+cps.obs.exp.atagct+cps.obs.exp.atagga+cps.obs.exp.ataggc+cps.obs.exp.ataggg+cps.obs.exp.ataggt+cps.obs.exp.atagta+cps.obs.exp.atagtc+cps.obs.exp.atagtg+cps.obs.exp.atagtt+cps.obs.exp.atataa+cps.obs.exp.atatac+cps.obs.exp.atatag+cps.obs.exp.atatat+cps.obs.exp.atatca+cps.obs.exp.atatcc+cps.obs.exp.atatcg+cps.obs.exp.atatct+cps.obs.exp.atatga+cps.obs.exp.atatgc+cps.obs.exp.atatgg+cps.obs.exp.atatgt+cps.obs.exp.atatta+cps.obs.exp.atattc+cps.obs.exp.atattg+cps.obs.exp.atattt+cps.obs.exp.atcaaa+cps.obs.exp.atcaac+cps.obs.exp.atcaag+cps.obs.exp.atcaat+cps.obs.exp.atcaca+cps.obs.exp.atcacc+cps.obs.exp.atcacg+cps.obs.exp.atcact+cps.obs.exp.atcaga+cps.obs.exp.atcagc+cps.obs.exp.atcagg+cps.obs.exp.atcagt+cps.obs.exp.atcata+cps.obs.exp.atcatc+cps.obs.exp.atcatg+cps.obs.exp.atcatt+cps.obs.exp.atccaa+cps.obs.exp.atccac+cps.obs.exp.atccag+cps.obs.exp.atccat+cps.obs.exp.atccca+cps.obs.exp.atcccc+cps.obs.exp.atcccg+cps.obs.exp.atccct+cps.obs.exp.atccga+cps.obs.exp.atccgc+cps.obs.exp.atccgg+cps.obs.exp.atccgt+cps.obs.exp.atccta+cps.obs.exp.atcctc+cps.obs.exp.atcctg+cps.obs.exp.atcctt+cps.obs.exp.atcgaa+cps.obs.exp.atcgac+cps.obs.exp.atcgag+cps.obs.exp.atcgat+cps.obs.exp.atcgca+cps.obs.exp.atcgcc+cps.obs.exp.atcgcg+cps.obs.exp.atcgct+cps.obs.exp.atcgga+cps.obs.exp.atcggc+cps.obs.exp.atcggg+cps.obs.exp.atcggt+cps.obs.exp.atcgta+cps.obs.exp.atcgtc+cps.obs.exp.atcgtg+cps.obs.exp.atcgtt+cps.obs.exp.atctaa+cps.obs.exp.atctac+cps.obs.exp.atctag+cps.obs.exp.atctat+cps.obs.exp.atctca+cps.obs.exp.atctcc+cps.obs.exp.atctcg+cps.obs.exp.atctct+cps.obs.exp.atctga+cps.obs.exp.atctgc+cps.obs.exp.atctgg+cps.obs.exp.atctgt+cps.obs.exp.atctta+cps.obs.exp.atcttc+cps.obs.exp.atcttg+cps.obs.exp.atcttt+cps.obs.exp.atgaaa+cps.obs.exp.atgaac+cps.obs.exp.atgaag+cps.obs.exp.atgaat+cps.obs.exp.atgaca+cps.obs.exp.atgacc+cps.obs.exp.atgacg+cps.obs.exp.atgact+cps.obs.exp.atgaga+cps.obs.exp.atgagc+cps.obs.exp.atgagg+cps.obs.exp.atgagt+cps.obs.exp.atgata+cps.obs.exp.atgatc+cps.obs.exp.atgatg+cps.obs.exp.atgatt+cps.obs.exp.atgcaa+cps.obs.exp.atgcac+cps.obs.exp.atgcag+cps.obs.exp.atgcat+cps.obs.exp.atgcca+cps.obs.exp.atgccc+cps.obs.exp.atgccg+cps.obs.exp.atgcct+cps.obs.exp.atgcga+cps.obs.exp.atgcgc+cps.obs.exp.atgcgg+cps.obs.exp.atgcgt+cps.obs.exp.atgcta+cps.obs.exp.atgctc+cps.obs.exp.atgctg+cps.obs.exp.atgctt+cps.obs.exp.atggaa+cps.obs.exp.atggac+cps.obs.exp.atggag+cps.obs.exp.atggat+cps.obs.exp.atggca+cps.obs.exp.atggcc+cps.obs.exp.atggcg+cps.obs.exp.atggct+cps.obs.exp.atggga+cps.obs.exp.atgggc+cps.obs.exp.atgggg+cps.obs.exp.atgggt+cps.obs.exp.atggta+cps.obs.exp.atggtc+cps.obs.exp.atggtg+cps.obs.exp.atggtt+cps.obs.exp.atgtaa+cps.obs.exp.atgtac+cps.obs.exp.atgtag+cps.obs.exp.atgtat+cps.obs.exp.atgtca+cps.obs.exp.atgtcc+cps.obs.exp.atgtcg+cps.obs.exp.atgtct+cps.obs.exp.atgtga+cps.obs.exp.atgtgc+cps.obs.exp.atgtgg+cps.obs.exp.atgtgt+cps.obs.exp.atgtta+cps.obs.exp.atgttc+cps.obs.exp.atgttg+cps.obs.exp.atgttt+cps.obs.exp.attaaa+cps.obs.exp.attaac+cps.obs.exp.attaag+cps.obs.exp.attaat+cps.obs.exp.attaca+cps.obs.exp.attacc+cps.obs.exp.attacg+cps.obs.exp.attact+cps.obs.exp.attaga+cps.obs.exp.attagc+cps.obs.exp.attagg+cps.obs.exp.attagt+cps.obs.exp.attata+cps.obs.exp.attatc+cps.obs.exp.attatg+cps.obs.exp.attatt+cps.obs.exp.attcaa+cps.obs.exp.attcac+cps.obs.exp.attcag+cps.obs.exp.attcat+cps.obs.exp.attcca+cps.obs.exp.attccc+cps.obs.exp.attccg+cps.obs.exp.attcct+cps.obs.exp.attcga+cps.obs.exp.attcgc+cps.obs.exp.attcgg+cps.obs.exp.attcgt+cps.obs.exp.attcta+cps.obs.exp.attctc+cps.obs.exp.attctg+cps.obs.exp.attctt+cps.obs.exp.attgaa+cps.obs.exp.attgac+cps.obs.exp.attgag+cps.obs.exp.attgat+cps.obs.exp.attgca+cps.obs.exp.attgcc+cps.obs.exp.attgcg+cps.obs.exp.attgct+cps.obs.exp.attgga+cps.obs.exp.attggc+cps.obs.exp.attggg+cps.obs.exp.attggt+cps.obs.exp.attgta+cps.obs.exp.attgtc+cps.obs.exp.attgtg+cps.obs.exp.attgtt+cps.obs.exp.atttaa+cps.obs.exp.atttac+cps.obs.exp.atttag+cps.obs.exp.atttat+cps.obs.exp.atttca+cps.obs.exp.atttcc+cps.obs.exp.atttcg+cps.obs.exp.atttct+cps.obs.exp.atttga+cps.obs.exp.atttgc+cps.obs.exp.atttgg+cps.obs.exp.atttgt+cps.obs.exp.atttta+cps.obs.exp.attttc+cps.obs.exp.attttg+cps.obs.exp.attttt+cps.obs.exp.caaaaa+cps.obs.exp.caaaac+cps.obs.exp.caaaag+cps.obs.exp.caaaat+cps.obs.exp.caaaca+cps.obs.exp.caaacc+cps.obs.exp.caaacg+cps.obs.exp.caaact+cps.obs.exp.caaaga+cps.obs.exp.caaagc+cps.obs.exp.caaagg+cps.obs.exp.caaagt+cps.obs.exp.caaata+cps.obs.exp.caaatc+cps.obs.exp.caaatg+cps.obs.exp.caaatt+cps.obs.exp.caacaa+cps.obs.exp.caacac+cps.obs.exp.caacag+cps.obs.exp.caacat+cps.obs.exp.caacca+cps.obs.exp.caaccc+cps.obs.exp.caaccg+cps.obs.exp.caacct+cps.obs.exp.caacga+cps.obs.exp.caacgc+cps.obs.exp.caacgg+cps.obs.exp.caacgt+cps.obs.exp.caacta+cps.obs.exp.caactc+cps.obs.exp.caactg+cps.obs.exp.caactt+cps.obs.exp.caagaa+cps.obs.exp.caagac+cps.obs.exp.caagag+cps.obs.exp.caagat+cps.obs.exp.caagca+cps.obs.exp.caagcc+cps.obs.exp.caagcg+cps.obs.exp.caagct+cps.obs.exp.caagga+cps.obs.exp.caaggc+cps.obs.exp.caaggg+cps.obs.exp.caaggt+cps.obs.exp.caagta+cps.obs.exp.caagtc+cps.obs.exp.caagtg+cps.obs.exp.caagtt+cps.obs.exp.caataa+cps.obs.exp.caatac+cps.obs.exp.caatag+cps.obs.exp.caatat+cps.obs.exp.caatca+cps.obs.exp.caatcc+cps.obs.exp.caatcg+cps.obs.exp.caatct+cps.obs.exp.caatga+cps.obs.exp.caatgc+cps.obs.exp.caatgg+cps.obs.exp.caatgt+cps.obs.exp.caatta+cps.obs.exp.caattc+cps.obs.exp.caattg+cps.obs.exp.caattt+cps.obs.exp.cacaaa+cps.obs.exp.cacaac+cps.obs.exp.cacaag+cps.obs.exp.cacaat+cps.obs.exp.cacaca+cps.obs.exp.cacacc+cps.obs.exp.cacacg+cps.obs.exp.cacact+cps.obs.exp.cacaga+cps.obs.exp.cacagc+cps.obs.exp.cacagg+cps.obs.exp.cacagt+cps.obs.exp.cacata+cps.obs.exp.cacatc+cps.obs.exp.cacatg+cps.obs.exp.cacatt+cps.obs.exp.caccaa+cps.obs.exp.caccac+cps.obs.exp.caccag+cps.obs.exp.caccat+cps.obs.exp.caccca+cps.obs.exp.cacccc+cps.obs.exp.cacccg+cps.obs.exp.caccct+cps.obs.exp.caccga+cps.obs.exp.caccgc+cps.obs.exp.caccgg+cps.obs.exp.caccgt+cps.obs.exp.caccta+cps.obs.exp.cacctc+cps.obs.exp.cacctg+cps.obs.exp.cacctt+cps.obs.exp.cacgaa+cps.obs.exp.cacgac+cps.obs.exp.cacgag+cps.obs.exp.cacgat+cps.obs.exp.cacgca+cps.obs.exp.cacgcc+cps.obs.exp.cacgcg+cps.obs.exp.cacgct+cps.obs.exp.cacgga+cps.obs.exp.cacggc+cps.obs.exp.cacggg+cps.obs.exp.cacggt+cps.obs.exp.cacgta+cps.obs.exp.cacgtc+cps.obs.exp.cacgtg+cps.obs.exp.cacgtt+cps.obs.exp.cactaa+cps.obs.exp.cactac+cps.obs.exp.cactag+cps.obs.exp.cactat+cps.obs.exp.cactca+cps.obs.exp.cactcc+cps.obs.exp.cactcg+cps.obs.exp.cactct+cps.obs.exp.cactga+cps.obs.exp.cactgc+cps.obs.exp.cactgg+cps.obs.exp.cactgt+cps.obs.exp.cactta+cps.obs.exp.cacttc+cps.obs.exp.cacttg+cps.obs.exp.cacttt+cps.obs.exp.cagaaa+cps.obs.exp.cagaac+cps.obs.exp.cagaag+cps.obs.exp.cagaat+cps.obs.exp.cagaca+cps.obs.exp.cagacc+cps.obs.exp.cagacg+cps.obs.exp.cagact+cps.obs.exp.cagaga+cps.obs.exp.cagagc+cps.obs.exp.cagagg+cps.obs.exp.cagagt+cps.obs.exp.cagata+cps.obs.exp.cagatc+cps.obs.exp.cagatg+cps.obs.exp.cagatt+cps.obs.exp.cagcaa+cps.obs.exp.cagcac+cps.obs.exp.cagcag+cps.obs.exp.cagcat+cps.obs.exp.cagcca+cps.obs.exp.cagccc+cps.obs.exp.cagccg+cps.obs.exp.cagcct+cps.obs.exp.cagcga+cps.obs.exp.cagcgc+cps.obs.exp.cagcgg+cps.obs.exp.cagcgt+cps.obs.exp.cagcta+cps.obs.exp.cagctc+cps.obs.exp.cagctg+cps.obs.exp.cagctt+cps.obs.exp.caggaa+cps.obs.exp.caggac+cps.obs.exp.caggag+cps.obs.exp.caggat+cps.obs.exp.caggca+cps.obs.exp.caggcc+cps.obs.exp.caggcg+cps.obs.exp.caggct+cps.obs.exp.caggga+cps.obs.exp.cagggc+cps.obs.exp.cagggg+cps.obs.exp.cagggt+cps.obs.exp.caggta+cps.obs.exp.caggtc+cps.obs.exp.caggtg+cps.obs.exp.caggtt+cps.obs.exp.cagtaa+cps.obs.exp.cagtac+cps.obs.exp.cagtag+cps.obs.exp.cagtat+cps.obs.exp.cagtca+cps.obs.exp.cagtcc+cps.obs.exp.cagtcg+cps.obs.exp.cagtct+cps.obs.exp.cagtga+cps.obs.exp.cagtgc+cps.obs.exp.cagtgg+cps.obs.exp.cagtgt+cps.obs.exp.cagtta+cps.obs.exp.cagttc+cps.obs.exp.cagttg+cps.obs.exp.cagttt+cps.obs.exp.cataaa+cps.obs.exp.cataac+cps.obs.exp.cataag+cps.obs.exp.cataat+cps.obs.exp.cataca+cps.obs.exp.catacc+cps.obs.exp.catacg+cps.obs.exp.catact+cps.obs.exp.cataga+cps.obs.exp.catagc+cps.obs.exp.catagg+cps.obs.exp.catagt+cps.obs.exp.catata+cps.obs.exp.catatc+cps.obs.exp.catatg+cps.obs.exp.catatt+cps.obs.exp.catcaa+cps.obs.exp.catcac+cps.obs.exp.catcag+cps.obs.exp.catcat+cps.obs.exp.catcca+cps.obs.exp.catccc+cps.obs.exp.catccg+cps.obs.exp.catcct+cps.obs.exp.catcga+cps.obs.exp.catcgc+cps.obs.exp.catcgg+cps.obs.exp.catcgt+cps.obs.exp.catcta+cps.obs.exp.catctc+cps.obs.exp.catctg+cps.obs.exp.catctt+cps.obs.exp.catgaa+cps.obs.exp.catgac+cps.obs.exp.catgag+cps.obs.exp.catgat+cps.obs.exp.catgca+cps.obs.exp.catgcc+cps.obs.exp.catgcg+cps.obs.exp.catgct+cps.obs.exp.catgga+cps.obs.exp.catggc+cps.obs.exp.catggg+cps.obs.exp.catggt+cps.obs.exp.catgta+cps.obs.exp.catgtc+cps.obs.exp.catgtg+cps.obs.exp.catgtt+cps.obs.exp.cattaa+cps.obs.exp.cattac+cps.obs.exp.cattag+cps.obs.exp.cattat+cps.obs.exp.cattca+cps.obs.exp.cattcc+cps.obs.exp.cattcg+cps.obs.exp.cattct+cps.obs.exp.cattga+cps.obs.exp.cattgc+cps.obs.exp.cattgg+cps.obs.exp.cattgt+cps.obs.exp.cattta+cps.obs.exp.catttc+cps.obs.exp.catttg+cps.obs.exp.catttt+cps.obs.exp.ccaaaa+cps.obs.exp.ccaaac+cps.obs.exp.ccaaag+cps.obs.exp.ccaaat+cps.obs.exp.ccaaca+cps.obs.exp.ccaacc+cps.obs.exp.ccaacg+cps.obs.exp.ccaact+cps.obs.exp.ccaaga+cps.obs.exp.ccaagc+cps.obs.exp.ccaagg+cps.obs.exp.ccaagt+cps.obs.exp.ccaata+cps.obs.exp.ccaatc+cps.obs.exp.ccaatg+cps.obs.exp.ccaatt+cps.obs.exp.ccacaa+cps.obs.exp.ccacac+cps.obs.exp.ccacag+cps.obs.exp.ccacat+cps.obs.exp.ccacca+cps.obs.exp.ccaccc+cps.obs.exp.ccaccg+cps.obs.exp.ccacct+cps.obs.exp.ccacga+cps.obs.exp.ccacgc+cps.obs.exp.ccacgg+cps.obs.exp.ccacgt+cps.obs.exp.ccacta+cps.obs.exp.ccactc+cps.obs.exp.ccactg+cps.obs.exp.ccactt+cps.obs.exp.ccagaa+cps.obs.exp.ccagac+cps.obs.exp.ccagag+cps.obs.exp.ccagat+cps.obs.exp.ccagca+cps.obs.exp.ccagcc+cps.obs.exp.ccagcg+cps.obs.exp.ccagct+cps.obs.exp.ccagga+cps.obs.exp.ccaggc+cps.obs.exp.ccaggg+cps.obs.exp.ccaggt+cps.obs.exp.ccagta+cps.obs.exp.ccagtc+cps.obs.exp.ccagtg+cps.obs.exp.ccagtt+cps.obs.exp.ccataa+cps.obs.exp.ccatac+cps.obs.exp.ccatag+cps.obs.exp.ccatat+cps.obs.exp.ccatca+cps.obs.exp.ccatcc+cps.obs.exp.ccatcg+cps.obs.exp.ccatct+cps.obs.exp.ccatga+cps.obs.exp.ccatgc+cps.obs.exp.ccatgg+cps.obs.exp.ccatgt+cps.obs.exp.ccatta+cps.obs.exp.ccattc+cps.obs.exp.ccattg+cps.obs.exp.ccattt+cps.obs.exp.cccaaa+cps.obs.exp.cccaac+cps.obs.exp.cccaag+cps.obs.exp.cccaat+cps.obs.exp.cccaca+cps.obs.exp.cccacc+cps.obs.exp.cccacg+cps.obs.exp.cccact+cps.obs.exp.cccaga+cps.obs.exp.cccagc+cps.obs.exp.cccagg+cps.obs.exp.cccagt+cps.obs.exp.cccata+cps.obs.exp.cccatc+cps.obs.exp.cccatg+cps.obs.exp.cccatt+cps.obs.exp.ccccaa+cps.obs.exp.ccccac+cps.obs.exp.ccccag+cps.obs.exp.ccccat+cps.obs.exp.ccccca+cps.obs.exp.cccccc+cps.obs.exp.cccccg+cps.obs.exp.ccccct+cps.obs.exp.ccccga+cps.obs.exp.ccccgc+cps.obs.exp.ccccgg+cps.obs.exp.ccccgt+cps.obs.exp.ccccta+cps.obs.exp.cccctc+cps.obs.exp.cccctg+cps.obs.exp.cccctt+cps.obs.exp.cccgaa+cps.obs.exp.cccgac+cps.obs.exp.cccgag+cps.obs.exp.cccgat+cps.obs.exp.cccgca+cps.obs.exp.cccgcc+cps.obs.exp.cccgcg+cps.obs.exp.cccgct+cps.obs.exp.cccgga+cps.obs.exp.cccggc+cps.obs.exp.cccggg+cps.obs.exp.cccggt+cps.obs.exp.cccgta+cps.obs.exp.cccgtc+cps.obs.exp.cccgtg+cps.obs.exp.cccgtt+cps.obs.exp.ccctaa+cps.obs.exp.ccctac+cps.obs.exp.ccctag+cps.obs.exp.ccctat+cps.obs.exp.ccctca+cps.obs.exp.ccctcc+cps.obs.exp.ccctcg+cps.obs.exp.ccctct+cps.obs.exp.ccctga+cps.obs.exp.ccctgc+cps.obs.exp.ccctgg+cps.obs.exp.ccctgt+cps.obs.exp.ccctta+cps.obs.exp.cccttc+cps.obs.exp.cccttg+cps.obs.exp.cccttt+cps.obs.exp.ccgaaa+cps.obs.exp.ccgaac+cps.obs.exp.ccgaag+cps.obs.exp.ccgaat+cps.obs.exp.ccgaca+cps.obs.exp.ccgacc+cps.obs.exp.ccgacg+cps.obs.exp.ccgact+cps.obs.exp.ccgaga+cps.obs.exp.ccgagc+cps.obs.exp.ccgagg+cps.obs.exp.ccgagt+cps.obs.exp.ccgata+cps.obs.exp.ccgatc+cps.obs.exp.ccgatg+cps.obs.exp.ccgatt+cps.obs.exp.ccgcaa+cps.obs.exp.ccgcac+cps.obs.exp.ccgcag+cps.obs.exp.ccgcat+cps.obs.exp.ccgcca+cps.obs.exp.ccgccc+cps.obs.exp.ccgccg+cps.obs.exp.ccgcct+cps.obs.exp.ccgcga+cps.obs.exp.ccgcgc+cps.obs.exp.ccgcgg+cps.obs.exp.ccgcgt+cps.obs.exp.ccgcta+cps.obs.exp.ccgctc+cps.obs.exp.ccgctg+cps.obs.exp.ccgctt+cps.obs.exp.ccggaa+cps.obs.exp.ccggac+cps.obs.exp.ccggag+cps.obs.exp.ccggat+cps.obs.exp.ccggca+cps.obs.exp.ccggcc+cps.obs.exp.ccggcg+cps.obs.exp.ccggct+cps.obs.exp.ccggga+cps.obs.exp.ccgggc+cps.obs.exp.ccgggg+cps.obs.exp.ccgggt+cps.obs.exp.ccggta+cps.obs.exp.ccggtc+cps.obs.exp.ccggtg+cps.obs.exp.ccggtt+cps.obs.exp.ccgtaa+cps.obs.exp.ccgtac+cps.obs.exp.ccgtag+cps.obs.exp.ccgtat+cps.obs.exp.ccgtca+cps.obs.exp.ccgtcc+cps.obs.exp.ccgtcg+cps.obs.exp.ccgtct+cps.obs.exp.ccgtga+cps.obs.exp.ccgtgc+cps.obs.exp.ccgtgg+cps.obs.exp.ccgtgt+cps.obs.exp.ccgtta+cps.obs.exp.ccgttc+cps.obs.exp.ccgttg+cps.obs.exp.ccgttt+cps.obs.exp.cctaaa+cps.obs.exp.cctaac+cps.obs.exp.cctaag+cps.obs.exp.cctaat+cps.obs.exp.cctaca+cps.obs.exp.cctacc+cps.obs.exp.cctacg+cps.obs.exp.cctact+cps.obs.exp.cctaga+cps.obs.exp.cctagc+cps.obs.exp.cctagg+cps.obs.exp.cctagt+cps.obs.exp.cctata+cps.obs.exp.cctatc+cps.obs.exp.cctatg+cps.obs.exp.cctatt+cps.obs.exp.cctcaa+cps.obs.exp.cctcac+cps.obs.exp.cctcag+cps.obs.exp.cctcat+cps.obs.exp.cctcca+cps.obs.exp.cctccc+cps.obs.exp.cctccg+cps.obs.exp.cctcct+cps.obs.exp.cctcga+cps.obs.exp.cctcgc+cps.obs.exp.cctcgg+cps.obs.exp.cctcgt+cps.obs.exp.cctcta+cps.obs.exp.cctctc+cps.obs.exp.cctctg+cps.obs.exp.cctctt+cps.obs.exp.cctgaa+cps.obs.exp.cctgac+cps.obs.exp.cctgag+cps.obs.exp.cctgat+cps.obs.exp.cctgca+cps.obs.exp.cctgcc+cps.obs.exp.cctgcg+cps.obs.exp.cctgct+cps.obs.exp.cctgga+cps.obs.exp.cctggc+cps.obs.exp.cctggg+cps.obs.exp.cctggt+cps.obs.exp.cctgta+cps.obs.exp.cctgtc+cps.obs.exp.cctgtg+cps.obs.exp.cctgtt+cps.obs.exp.ccttaa+cps.obs.exp.ccttac+cps.obs.exp.ccttag+cps.obs.exp.ccttat+cps.obs.exp.ccttca+cps.obs.exp.ccttcc+cps.obs.exp.ccttcg+cps.obs.exp.ccttct+cps.obs.exp.ccttga+cps.obs.exp.ccttgc+cps.obs.exp.ccttgg+cps.obs.exp.ccttgt+cps.obs.exp.ccttta+cps.obs.exp.cctttc+cps.obs.exp.cctttg+cps.obs.exp.cctttt+cps.obs.exp.cgaaaa+cps.obs.exp.cgaaac+cps.obs.exp.cgaaag+cps.obs.exp.cgaaat+cps.obs.exp.cgaaca+cps.obs.exp.cgaacc+cps.obs.exp.cgaacg+cps.obs.exp.cgaact+cps.obs.exp.cgaaga+cps.obs.exp.cgaagc+cps.obs.exp.cgaagg+cps.obs.exp.cgaagt+cps.obs.exp.cgaata+cps.obs.exp.cgaatc+cps.obs.exp.cgaatg+cps.obs.exp.cgaatt+cps.obs.exp.cgacaa+cps.obs.exp.cgacac+cps.obs.exp.cgacag+cps.obs.exp.cgacat+cps.obs.exp.cgacca+cps.obs.exp.cgaccc+cps.obs.exp.cgaccg+cps.obs.exp.cgacct+cps.obs.exp.cgacga+cps.obs.exp.cgacgc+cps.obs.exp.cgacgg+cps.obs.exp.cgacgt+cps.obs.exp.cgacta+cps.obs.exp.cgactc+cps.obs.exp.cgactg+cps.obs.exp.cgactt+cps.obs.exp.cgagaa+cps.obs.exp.cgagac+cps.obs.exp.cgagag+cps.obs.exp.cgagat+cps.obs.exp.cgagca+cps.obs.exp.cgagcc+cps.obs.exp.cgagcg+cps.obs.exp.cgagct+cps.obs.exp.cgagga+cps.obs.exp.cgaggc+cps.obs.exp.cgaggg+cps.obs.exp.cgaggt+cps.obs.exp.cgagta+cps.obs.exp.cgagtc+cps.obs.exp.cgagtg+cps.obs.exp.cgagtt+cps.obs.exp.cgataa+cps.obs.exp.cgatac+cps.obs.exp.cgatag+cps.obs.exp.cgatat+cps.obs.exp.cgatca+cps.obs.exp.cgatcc+cps.obs.exp.cgatcg+cps.obs.exp.cgatct+cps.obs.exp.cgatga+cps.obs.exp.cgatgc+cps.obs.exp.cgatgg+cps.obs.exp.cgatgt+cps.obs.exp.cgatta+cps.obs.exp.cgattc+cps.obs.exp.cgattg+cps.obs.exp.cgattt+cps.obs.exp.cgcaaa+cps.obs.exp.cgcaac+cps.obs.exp.cgcaag+cps.obs.exp.cgcaat+cps.obs.exp.cgcaca+cps.obs.exp.cgcacc+cps.obs.exp.cgcacg+cps.obs.exp.cgcact+cps.obs.exp.cgcaga+cps.obs.exp.cgcagc+cps.obs.exp.cgcagg+cps.obs.exp.cgcagt+cps.obs.exp.cgcata+cps.obs.exp.cgcatc+cps.obs.exp.cgcatg+cps.obs.exp.cgcatt+cps.obs.exp.cgccaa+cps.obs.exp.cgccac+cps.obs.exp.cgccag+cps.obs.exp.cgccat+cps.obs.exp.cgccca+cps.obs.exp.cgcccc+cps.obs.exp.cgcccg+cps.obs.exp.cgccct+cps.obs.exp.cgccga+cps.obs.exp.cgccgc+cps.obs.exp.cgccgg+cps.obs.exp.cgccgt+cps.obs.exp.cgccta+cps.obs.exp.cgcctc+cps.obs.exp.cgcctg+cps.obs.exp.cgcctt+cps.obs.exp.cgcgaa+cps.obs.exp.cgcgac+cps.obs.exp.cgcgag+cps.obs.exp.cgcgat+cps.obs.exp.cgcgca+cps.obs.exp.cgcgcc+cps.obs.exp.cgcgcg+cps.obs.exp.cgcgct+cps.obs.exp.cgcgga+cps.obs.exp.cgcggc+cps.obs.exp.cgcggg+cps.obs.exp.cgcggt+cps.obs.exp.cgcgta+cps.obs.exp.cgcgtc+cps.obs.exp.cgcgtg+cps.obs.exp.cgcgtt+cps.obs.exp.cgctaa+cps.obs.exp.cgctac+cps.obs.exp.cgctag+cps.obs.exp.cgctat+cps.obs.exp.cgctca+cps.obs.exp.cgctcc+cps.obs.exp.cgctcg+cps.obs.exp.cgctct+cps.obs.exp.cgctga+cps.obs.exp.cgctgc+cps.obs.exp.cgctgg+cps.obs.exp.cgctgt+cps.obs.exp.cgctta+cps.obs.exp.cgcttc+cps.obs.exp.cgcttg+cps.obs.exp.cgcttt+cps.obs.exp.cggaaa+cps.obs.exp.cggaac+cps.obs.exp.cggaag+cps.obs.exp.cggaat+cps.obs.exp.cggaca+cps.obs.exp.cggacc+cps.obs.exp.cggacg+cps.obs.exp.cggact+cps.obs.exp.cggaga+cps.obs.exp.cggagc+cps.obs.exp.cggagg+cps.obs.exp.cggagt+cps.obs.exp.cggata+cps.obs.exp.cggatc+cps.obs.exp.cggatg+cps.obs.exp.cggatt+cps.obs.exp.cggcaa+cps.obs.exp.cggcac+cps.obs.exp.cggcag+cps.obs.exp.cggcat+cps.obs.exp.cggcca+cps.obs.exp.cggccc+cps.obs.exp.cggccg+cps.obs.exp.cggcct+cps.obs.exp.cggcga+cps.obs.exp.cggcgc+cps.obs.exp.cggcgg+cps.obs.exp.cggcgt+cps.obs.exp.cggcta+cps.obs.exp.cggctc+cps.obs.exp.cggctg+cps.obs.exp.cggctt+cps.obs.exp.cgggaa+cps.obs.exp.cgggac+cps.obs.exp.cgggag+cps.obs.exp.cgggat+cps.obs.exp.cgggca+cps.obs.exp.cgggcc+cps.obs.exp.cgggcg+cps.obs.exp.cgggct+cps.obs.exp.cgggga+cps.obs.exp.cggggc+cps.obs.exp.cggggg+cps.obs.exp.cggggt+cps.obs.exp.cgggta+cps.obs.exp.cgggtc+cps.obs.exp.cgggtg+cps.obs.exp.cgggtt+cps.obs.exp.cggtaa+cps.obs.exp.cggtac+cps.obs.exp.cggtag+cps.obs.exp.cggtat+cps.obs.exp.cggtca+cps.obs.exp.cggtcc+cps.obs.exp.cggtcg+cps.obs.exp.cggtct+cps.obs.exp.cggtga+cps.obs.exp.cggtgc+cps.obs.exp.cggtgg+cps.obs.exp.cggtgt+cps.obs.exp.cggtta+cps.obs.exp.cggttc+cps.obs.exp.cggttg+cps.obs.exp.cggttt+cps.obs.exp.cgtaaa+cps.obs.exp.cgtaac+cps.obs.exp.cgtaag+cps.obs.exp.cgtaat+cps.obs.exp.cgtaca+cps.obs.exp.cgtacc+cps.obs.exp.cgtacg+cps.obs.exp.cgtact+cps.obs.exp.cgtaga+cps.obs.exp.cgtagc+cps.obs.exp.cgtagg+cps.obs.exp.cgtagt+cps.obs.exp.cgtata+cps.obs.exp.cgtatc+cps.obs.exp.cgtatg+cps.obs.exp.cgtatt+cps.obs.exp.cgtcaa+cps.obs.exp.cgtcac+cps.obs.exp.cgtcag+cps.obs.exp.cgtcat+cps.obs.exp.cgtcca+cps.obs.exp.cgtccc+cps.obs.exp.cgtccg+cps.obs.exp.cgtcct+cps.obs.exp.cgtcga+cps.obs.exp.cgtcgc+cps.obs.exp.cgtcgg+cps.obs.exp.cgtcgt+cps.obs.exp.cgtcta+cps.obs.exp.cgtctc+cps.obs.exp.cgtctg+cps.obs.exp.cgtctt+cps.obs.exp.cgtgaa+cps.obs.exp.cgtgac+cps.obs.exp.cgtgag+cps.obs.exp.cgtgat+cps.obs.exp.cgtgca+cps.obs.exp.cgtgcc+cps.obs.exp.cgtgcg+cps.obs.exp.cgtgct+cps.obs.exp.cgtgga+cps.obs.exp.cgtggc+cps.obs.exp.cgtggg+cps.obs.exp.cgtggt+cps.obs.exp.cgtgta+cps.obs.exp.cgtgtc+cps.obs.exp.cgtgtg+cps.obs.exp.cgtgtt+cps.obs.exp.cgttaa+cps.obs.exp.cgttac+cps.obs.exp.cgttag+cps.obs.exp.cgttat+cps.obs.exp.cgttca+cps.obs.exp.cgttcc+cps.obs.exp.cgttcg+cps.obs.exp.cgttct+cps.obs.exp.cgttga+cps.obs.exp.cgttgc+cps.obs.exp.cgttgg+cps.obs.exp.cgttgt+cps.obs.exp.cgttta+cps.obs.exp.cgtttc+cps.obs.exp.cgtttg+cps.obs.exp.cgtttt+cps.obs.exp.ctaaaa+cps.obs.exp.ctaaac+cps.obs.exp.ctaaag+cps.obs.exp.ctaaat+cps.obs.exp.ctaaca+cps.obs.exp.ctaacc+cps.obs.exp.ctaacg+cps.obs.exp.ctaact+cps.obs.exp.ctaaga+cps.obs.exp.ctaagc+cps.obs.exp.ctaagg+cps.obs.exp.ctaagt+cps.obs.exp.ctaata+cps.obs.exp.ctaatc+cps.obs.exp.ctaatg+cps.obs.exp.ctaatt+cps.obs.exp.ctacaa+cps.obs.exp.ctacac+cps.obs.exp.ctacag+cps.obs.exp.ctacat+cps.obs.exp.ctacca+cps.obs.exp.ctaccc+cps.obs.exp.ctaccg+cps.obs.exp.ctacct+cps.obs.exp.ctacga+cps.obs.exp.ctacgc+cps.obs.exp.ctacgg+cps.obs.exp.ctacgt+cps.obs.exp.ctacta+cps.obs.exp.ctactc+cps.obs.exp.ctactg+cps.obs.exp.ctactt+cps.obs.exp.ctagaa+cps.obs.exp.ctagac+cps.obs.exp.ctagag+cps.obs.exp.ctagat+cps.obs.exp.ctagca+cps.obs.exp.ctagcc+cps.obs.exp.ctagcg+cps.obs.exp.ctagct+cps.obs.exp.ctagga+cps.obs.exp.ctaggc+cps.obs.exp.ctaggg+cps.obs.exp.ctaggt+cps.obs.exp.ctagta+cps.obs.exp.ctagtc+cps.obs.exp.ctagtg+cps.obs.exp.ctagtt+cps.obs.exp.ctataa+cps.obs.exp.ctatac+cps.obs.exp.ctatag+cps.obs.exp.ctatat+cps.obs.exp.ctatca+cps.obs.exp.ctatcc+cps.obs.exp.ctatcg+cps.obs.exp.ctatct+cps.obs.exp.ctatga+cps.obs.exp.ctatgc+cps.obs.exp.ctatgg+cps.obs.exp.ctatgt+cps.obs.exp.ctatta+cps.obs.exp.ctattc+cps.obs.exp.ctattg+cps.obs.exp.ctattt+cps.obs.exp.ctcaaa+cps.obs.exp.ctcaac+cps.obs.exp.ctcaag+cps.obs.exp.ctcaat+cps.obs.exp.ctcaca+cps.obs.exp.ctcacc+cps.obs.exp.ctcacg+cps.obs.exp.ctcact+cps.obs.exp.ctcaga+cps.obs.exp.ctcagc+cps.obs.exp.ctcagg+cps.obs.exp.ctcagt+cps.obs.exp.ctcata+cps.obs.exp.ctcatc+cps.obs.exp.ctcatg+cps.obs.exp.ctcatt+cps.obs.exp.ctccaa+cps.obs.exp.ctccac+cps.obs.exp.ctccag+cps.obs.exp.ctccat+cps.obs.exp.ctccca+cps.obs.exp.ctcccc+cps.obs.exp.ctcccg+cps.obs.exp.ctccct+cps.obs.exp.ctccga+cps.obs.exp.ctccgc+cps.obs.exp.ctccgg+cps.obs.exp.ctccgt+cps.obs.exp.ctccta+cps.obs.exp.ctcctc+cps.obs.exp.ctcctg+cps.obs.exp.ctcctt+cps.obs.exp.ctcgaa+cps.obs.exp.ctcgac+cps.obs.exp.ctcgag+cps.obs.exp.ctcgat+cps.obs.exp.ctcgca+cps.obs.exp.ctcgcc+cps.obs.exp.ctcgcg+cps.obs.exp.ctcgct+cps.obs.exp.ctcgga+cps.obs.exp.ctcggc+cps.obs.exp.ctcggg+cps.obs.exp.ctcggt+cps.obs.exp.ctcgta+cps.obs.exp.ctcgtc+cps.obs.exp.ctcgtg+cps.obs.exp.ctcgtt+cps.obs.exp.ctctaa+cps.obs.exp.ctctac+cps.obs.exp.ctctag+cps.obs.exp.ctctat+cps.obs.exp.ctctca+cps.obs.exp.ctctcc+cps.obs.exp.ctctcg+cps.obs.exp.ctctct+cps.obs.exp.ctctga+cps.obs.exp.ctctgc+cps.obs.exp.ctctgg+cps.obs.exp.ctctgt+cps.obs.exp.ctctta+cps.obs.exp.ctcttc+cps.obs.exp.ctcttg+cps.obs.exp.ctcttt+cps.obs.exp.ctgaaa+cps.obs.exp.ctgaac+cps.obs.exp.ctgaag+cps.obs.exp.ctgaat+cps.obs.exp.ctgaca+cps.obs.exp.ctgacc+cps.obs.exp.ctgacg+cps.obs.exp.ctgact+cps.obs.exp.ctgaga+cps.obs.exp.ctgagc+cps.obs.exp.ctgagg+cps.obs.exp.ctgagt+cps.obs.exp.ctgata+cps.obs.exp.ctgatc+cps.obs.exp.ctgatg+cps.obs.exp.ctgatt+cps.obs.exp.ctgcaa+cps.obs.exp.ctgcac+cps.obs.exp.ctgcag+cps.obs.exp.ctgcat+cps.obs.exp.ctgcca+cps.obs.exp.ctgccc+cps.obs.exp.ctgccg+cps.obs.exp.ctgcct+cps.obs.exp.ctgcga+cps.obs.exp.ctgcgc+cps.obs.exp.ctgcgg+cps.obs.exp.ctgcgt+cps.obs.exp.ctgcta+cps.obs.exp.ctgctc+cps.obs.exp.ctgctg+cps.obs.exp.ctgctt+cps.obs.exp.ctggaa+cps.obs.exp.ctggac+cps.obs.exp.ctggag+cps.obs.exp.ctggat+cps.obs.exp.ctggca+cps.obs.exp.ctggcc+cps.obs.exp.ctggcg+cps.obs.exp.ctggct+cps.obs.exp.ctggga+cps.obs.exp.ctgggc+cps.obs.exp.ctgggg+cps.obs.exp.ctgggt+cps.obs.exp.ctggta+cps.obs.exp.ctggtc+cps.obs.exp.ctggtg+cps.obs.exp.ctggtt+cps.obs.exp.ctgtaa+cps.obs.exp.ctgtac+cps.obs.exp.ctgtag+cps.obs.exp.ctgtat+cps.obs.exp.ctgtca+cps.obs.exp.ctgtcc+cps.obs.exp.ctgtcg+cps.obs.exp.ctgtct+cps.obs.exp.ctgtga+cps.obs.exp.ctgtgc+cps.obs.exp.ctgtgg+cps.obs.exp.ctgtgt+cps.obs.exp.ctgtta+cps.obs.exp.ctgttc+cps.obs.exp.ctgttg+cps.obs.exp.ctgttt+cps.obs.exp.cttaaa+cps.obs.exp.cttaac+cps.obs.exp.cttaag+cps.obs.exp.cttaat+cps.obs.exp.cttaca+cps.obs.exp.cttacc+cps.obs.exp.cttacg+cps.obs.exp.cttact+cps.obs.exp.cttaga+cps.obs.exp.cttagc+cps.obs.exp.cttagg+cps.obs.exp.cttagt+cps.obs.exp.cttata+cps.obs.exp.cttatc+cps.obs.exp.cttatg+cps.obs.exp.cttatt+cps.obs.exp.cttcaa+cps.obs.exp.cttcac+cps.obs.exp.cttcag+cps.obs.exp.cttcat+cps.obs.exp.cttcca+cps.obs.exp.cttccc+cps.obs.exp.cttccg+cps.obs.exp.cttcct+cps.obs.exp.cttcga+cps.obs.exp.cttcgc+cps.obs.exp.cttcgg+cps.obs.exp.cttcgt+cps.obs.exp.cttcta+cps.obs.exp.cttctc+cps.obs.exp.cttctg+cps.obs.exp.cttctt+cps.obs.exp.cttgaa+cps.obs.exp.cttgac+cps.obs.exp.cttgag+cps.obs.exp.cttgat+cps.obs.exp.cttgca+cps.obs.exp.cttgcc+cps.obs.exp.cttgcg+cps.obs.exp.cttgct+cps.obs.exp.cttgga+cps.obs.exp.cttggc+cps.obs.exp.cttggg+cps.obs.exp.cttggt+cps.obs.exp.cttgta+cps.obs.exp.cttgtc+cps.obs.exp.cttgtg+cps.obs.exp.cttgtt+cps.obs.exp.ctttaa+cps.obs.exp.ctttac+cps.obs.exp.ctttag+cps.obs.exp.ctttat+cps.obs.exp.ctttca+cps.obs.exp.ctttcc+cps.obs.exp.ctttcg+cps.obs.exp.ctttct+cps.obs.exp.ctttga+cps.obs.exp.ctttgc+cps.obs.exp.ctttgg+cps.obs.exp.ctttgt+cps.obs.exp.ctttta+cps.obs.exp.cttttc+cps.obs.exp.cttttg+cps.obs.exp.cttttt+cps.obs.exp.gaaaaa+cps.obs.exp.gaaaac+cps.obs.exp.gaaaag+cps.obs.exp.gaaaat+cps.obs.exp.gaaaca+cps.obs.exp.gaaacc+cps.obs.exp.gaaacg+cps.obs.exp.gaaact+cps.obs.exp.gaaaga+cps.obs.exp.gaaagc+cps.obs.exp.gaaagg+cps.obs.exp.gaaagt+cps.obs.exp.gaaata+cps.obs.exp.gaaatc+cps.obs.exp.gaaatg+cps.obs.exp.gaaatt+cps.obs.exp.gaacaa+cps.obs.exp.gaacac+cps.obs.exp.gaacag+cps.obs.exp.gaacat+cps.obs.exp.gaacca+cps.obs.exp.gaaccc+cps.obs.exp.gaaccg+cps.obs.exp.gaacct+cps.obs.exp.gaacga+cps.obs.exp.gaacgc+cps.obs.exp.gaacgg+cps.obs.exp.gaacgt+cps.obs.exp.gaacta+cps.obs.exp.gaactc+cps.obs.exp.gaactg+cps.obs.exp.gaactt+cps.obs.exp.gaagaa+cps.obs.exp.gaagac+cps.obs.exp.gaagag+cps.obs.exp.gaagat+cps.obs.exp.gaagca+cps.obs.exp.gaagcc+cps.obs.exp.gaagcg+cps.obs.exp.gaagct+cps.obs.exp.gaagga+cps.obs.exp.gaaggc+cps.obs.exp.gaaggg+cps.obs.exp.gaaggt+cps.obs.exp.gaagta+cps.obs.exp.gaagtc+cps.obs.exp.gaagtg+cps.obs.exp.gaagtt+cps.obs.exp.gaataa+cps.obs.exp.gaatac+cps.obs.exp.gaatag+cps.obs.exp.gaatat+cps.obs.exp.gaatca+cps.obs.exp.gaatcc+cps.obs.exp.gaatcg+cps.obs.exp.gaatct+cps.obs.exp.gaatga+cps.obs.exp.gaatgc+cps.obs.exp.gaatgg+cps.obs.exp.gaatgt+cps.obs.exp.gaatta+cps.obs.exp.gaattc+cps.obs.exp.gaattg+cps.obs.exp.gaattt+cps.obs.exp.gacaaa+cps.obs.exp.gacaac+cps.obs.exp.gacaag+cps.obs.exp.gacaat+cps.obs.exp.gacaca+cps.obs.exp.gacacc+cps.obs.exp.gacacg+cps.obs.exp.gacact+cps.obs.exp.gacaga+cps.obs.exp.gacagc+cps.obs.exp.gacagg+cps.obs.exp.gacagt+cps.obs.exp.gacata+cps.obs.exp.gacatc+cps.obs.exp.gacatg+cps.obs.exp.gacatt+cps.obs.exp.gaccaa+cps.obs.exp.gaccac+cps.obs.exp.gaccag+cps.obs.exp.gaccat+cps.obs.exp.gaccca+cps.obs.exp.gacccc+cps.obs.exp.gacccg+cps.obs.exp.gaccct+cps.obs.exp.gaccga+cps.obs.exp.gaccgc+cps.obs.exp.gaccgg+cps.obs.exp.gaccgt+cps.obs.exp.gaccta+cps.obs.exp.gacctc+cps.obs.exp.gacctg+cps.obs.exp.gacctt+cps.obs.exp.gacgaa+cps.obs.exp.gacgac+cps.obs.exp.gacgag+cps.obs.exp.gacgat+cps.obs.exp.gacgca+cps.obs.exp.gacgcc+cps.obs.exp.gacgcg+cps.obs.exp.gacgct+cps.obs.exp.gacgga+cps.obs.exp.gacggc+cps.obs.exp.gacggg+cps.obs.exp.gacggt+cps.obs.exp.gacgta+cps.obs.exp.gacgtc+cps.obs.exp.gacgtg+cps.obs.exp.gacgtt+cps.obs.exp.gactaa+cps.obs.exp.gactac+cps.obs.exp.gactag+cps.obs.exp.gactat+cps.obs.exp.gactca+cps.obs.exp.gactcc+cps.obs.exp.gactcg+cps.obs.exp.gactct+cps.obs.exp.gactga+cps.obs.exp.gactgc+cps.obs.exp.gactgg+cps.obs.exp.gactgt+cps.obs.exp.gactta+cps.obs.exp.gacttc+cps.obs.exp.gacttg+cps.obs.exp.gacttt+cps.obs.exp.gagaaa+cps.obs.exp.gagaac+cps.obs.exp.gagaag+cps.obs.exp.gagaat+cps.obs.exp.gagaca+cps.obs.exp.gagacc+cps.obs.exp.gagacg+cps.obs.exp.gagact+cps.obs.exp.gagaga+cps.obs.exp.gagagc+cps.obs.exp.gagagg+cps.obs.exp.gagagt+cps.obs.exp.gagata+cps.obs.exp.gagatc+cps.obs.exp.gagatg+cps.obs.exp.gagatt+cps.obs.exp.gagcaa+cps.obs.exp.gagcac+cps.obs.exp.gagcag+cps.obs.exp.gagcat+cps.obs.exp.gagcca+cps.obs.exp.gagccc+cps.obs.exp.gagccg+cps.obs.exp.gagcct+cps.obs.exp.gagcga+cps.obs.exp.gagcgc+cps.obs.exp.gagcgg+cps.obs.exp.gagcgt+cps.obs.exp.gagcta+cps.obs.exp.gagctc+cps.obs.exp.gagctg+cps.obs.exp.gagctt+cps.obs.exp.gaggaa+cps.obs.exp.gaggac+cps.obs.exp.gaggag+cps.obs.exp.gaggat+cps.obs.exp.gaggca+cps.obs.exp.gaggcc+cps.obs.exp.gaggcg+cps.obs.exp.gaggct+cps.obs.exp.gaggga+cps.obs.exp.gagggc+cps.obs.exp.gagggg+cps.obs.exp.gagggt+cps.obs.exp.gaggta+cps.obs.exp.gaggtc+cps.obs.exp.gaggtg+cps.obs.exp.gaggtt+cps.obs.exp.gagtaa+cps.obs.exp.gagtac+cps.obs.exp.gagtag+cps.obs.exp.gagtat+cps.obs.exp.gagtca+cps.obs.exp.gagtcc+cps.obs.exp.gagtcg+cps.obs.exp.gagtct+cps.obs.exp.gagtga+cps.obs.exp.gagtgc+cps.obs.exp.gagtgg+cps.obs.exp.gagtgt+cps.obs.exp.gagtta+cps.obs.exp.gagttc+cps.obs.exp.gagttg+cps.obs.exp.gagttt+cps.obs.exp.gataaa+cps.obs.exp.gataac+cps.obs.exp.gataag+cps.obs.exp.gataat+cps.obs.exp.gataca+cps.obs.exp.gatacc+cps.obs.exp.gatacg+cps.obs.exp.gatact+cps.obs.exp.gataga+cps.obs.exp.gatagc+cps.obs.exp.gatagg+cps.obs.exp.gatagt+cps.obs.exp.gatata+cps.obs.exp.gatatc+cps.obs.exp.gatatg+cps.obs.exp.gatatt+cps.obs.exp.gatcaa+cps.obs.exp.gatcac+cps.obs.exp.gatcag+cps.obs.exp.gatcat+cps.obs.exp.gatcca+cps.obs.exp.gatccc+cps.obs.exp.gatccg+cps.obs.exp.gatcct+cps.obs.exp.gatcga+cps.obs.exp.gatcgc+cps.obs.exp.gatcgg+cps.obs.exp.gatcgt+cps.obs.exp.gatcta+cps.obs.exp.gatctc+cps.obs.exp.gatctg+cps.obs.exp.gatctt+cps.obs.exp.gatgaa+cps.obs.exp.gatgac+cps.obs.exp.gatgag+cps.obs.exp.gatgat+cps.obs.exp.gatgca+cps.obs.exp.gatgcc+cps.obs.exp.gatgcg+cps.obs.exp.gatgct+cps.obs.exp.gatgga+cps.obs.exp.gatggc+cps.obs.exp.gatggg+cps.obs.exp.gatggt+cps.obs.exp.gatgta+cps.obs.exp.gatgtc+cps.obs.exp.gatgtg+cps.obs.exp.gatgtt+cps.obs.exp.gattaa+cps.obs.exp.gattac+cps.obs.exp.gattag+cps.obs.exp.gattat+cps.obs.exp.gattca+cps.obs.exp.gattcc+cps.obs.exp.gattcg+cps.obs.exp.gattct+cps.obs.exp.gattga+cps.obs.exp.gattgc+cps.obs.exp.gattgg+cps.obs.exp.gattgt+cps.obs.exp.gattta+cps.obs.exp.gatttc+cps.obs.exp.gatttg+cps.obs.exp.gatttt+cps.obs.exp.gcaaaa+cps.obs.exp.gcaaac+cps.obs.exp.gcaaag+cps.obs.exp.gcaaat+cps.obs.exp.gcaaca+cps.obs.exp.gcaacc+cps.obs.exp.gcaacg+cps.obs.exp.gcaact+cps.obs.exp.gcaaga+cps.obs.exp.gcaagc+cps.obs.exp.gcaagg+cps.obs.exp.gcaagt+cps.obs.exp.gcaata+cps.obs.exp.gcaatc+cps.obs.exp.gcaatg+cps.obs.exp.gcaatt+cps.obs.exp.gcacaa+cps.obs.exp.gcacac+cps.obs.exp.gcacag+cps.obs.exp.gcacat+cps.obs.exp.gcacca+cps.obs.exp.gcaccc+cps.obs.exp.gcaccg+cps.obs.exp.gcacct+cps.obs.exp.gcacga+cps.obs.exp.gcacgc+cps.obs.exp.gcacgg+cps.obs.exp.gcacgt+cps.obs.exp.gcacta+cps.obs.exp.gcactc+cps.obs.exp.gcactg+cps.obs.exp.gcactt+cps.obs.exp.gcagaa+cps.obs.exp.gcagac+cps.obs.exp.gcagag+cps.obs.exp.gcagat+cps.obs.exp.gcagca+cps.obs.exp.gcagcc+cps.obs.exp.gcagcg+cps.obs.exp.gcagct+cps.obs.exp.gcagga+cps.obs.exp.gcaggc+cps.obs.exp.gcaggg+cps.obs.exp.gcaggt+cps.obs.exp.gcagta+cps.obs.exp.gcagtc+cps.obs.exp.gcagtg+cps.obs.exp.gcagtt+cps.obs.exp.gcataa+cps.obs.exp.gcatac+cps.obs.exp.gcatag+cps.obs.exp.gcatat+cps.obs.exp.gcatca+cps.obs.exp.gcatcc+cps.obs.exp.gcatcg+cps.obs.exp.gcatct+cps.obs.exp.gcatga+cps.obs.exp.gcatgc+cps.obs.exp.gcatgg+cps.obs.exp.gcatgt+cps.obs.exp.gcatta+cps.obs.exp.gcattc+cps.obs.exp.gcattg+cps.obs.exp.gcattt+cps.obs.exp.gccaaa+cps.obs.exp.gccaac+cps.obs.exp.gccaag+cps.obs.exp.gccaat+cps.obs.exp.gccaca+cps.obs.exp.gccacc+cps.obs.exp.gccacg+cps.obs.exp.gccact+cps.obs.exp.gccaga+cps.obs.exp.gccagc+cps.obs.exp.gccagg+cps.obs.exp.gccagt+cps.obs.exp.gccata+cps.obs.exp.gccatc+cps.obs.exp.gccatg+cps.obs.exp.gccatt+cps.obs.exp.gcccaa+cps.obs.exp.gcccac+cps.obs.exp.gcccag+cps.obs.exp.gcccat+cps.obs.exp.gcccca+cps.obs.exp.gccccc+cps.obs.exp.gccccg+cps.obs.exp.gcccct+cps.obs.exp.gcccga+cps.obs.exp.gcccgc+cps.obs.exp.gcccgg+cps.obs.exp.gcccgt+cps.obs.exp.gcccta+cps.obs.exp.gccctc+cps.obs.exp.gccctg+cps.obs.exp.gccctt+cps.obs.exp.gccgaa+cps.obs.exp.gccgac+cps.obs.exp.gccgag+cps.obs.exp.gccgat+cps.obs.exp.gccgca+cps.obs.exp.gccgcc+cps.obs.exp.gccgcg+cps.obs.exp.gccgct+cps.obs.exp.gccgga+cps.obs.exp.gccggc+cps.obs.exp.gccggg+cps.obs.exp.gccggt+cps.obs.exp.gccgta+cps.obs.exp.gccgtc+cps.obs.exp.gccgtg+cps.obs.exp.gccgtt+cps.obs.exp.gcctaa+cps.obs.exp.gcctac+cps.obs.exp.gcctag+cps.obs.exp.gcctat+cps.obs.exp.gcctca+cps.obs.exp.gcctcc+cps.obs.exp.gcctcg+cps.obs.exp.gcctct+cps.obs.exp.gcctga+cps.obs.exp.gcctgc+cps.obs.exp.gcctgg+cps.obs.exp.gcctgt+cps.obs.exp.gcctta+cps.obs.exp.gccttc+cps.obs.exp.gccttg+cps.obs.exp.gccttt+cps.obs.exp.gcgaaa+cps.obs.exp.gcgaac+cps.obs.exp.gcgaag+cps.obs.exp.gcgaat+cps.obs.exp.gcgaca+cps.obs.exp.gcgacc+cps.obs.exp.gcgacg+cps.obs.exp.gcgact+cps.obs.exp.gcgaga+cps.obs.exp.gcgagc+cps.obs.exp.gcgagg+cps.obs.exp.gcgagt+cps.obs.exp.gcgata+cps.obs.exp.gcgatc+cps.obs.exp.gcgatg+cps.obs.exp.gcgatt+cps.obs.exp.gcgcaa+cps.obs.exp.gcgcac+cps.obs.exp.gcgcag+cps.obs.exp.gcgcat+cps.obs.exp.gcgcca+cps.obs.exp.gcgccc+cps.obs.exp.gcgccg+cps.obs.exp.gcgcct+cps.obs.exp.gcgcga+cps.obs.exp.gcgcgc+cps.obs.exp.gcgcgg+cps.obs.exp.gcgcgt+cps.obs.exp.gcgcta+cps.obs.exp.gcgctc+cps.obs.exp.gcgctg+cps.obs.exp.gcgctt+cps.obs.exp.gcggaa+cps.obs.exp.gcggac+cps.obs.exp.gcggag+cps.obs.exp.gcggat+cps.obs.exp.gcggca+cps.obs.exp.gcggcc+cps.obs.exp.gcggcg+cps.obs.exp.gcggct+cps.obs.exp.gcggga+cps.obs.exp.gcgggc+cps.obs.exp.gcgggg+cps.obs.exp.gcgggt+cps.obs.exp.gcggta+cps.obs.exp.gcggtc+cps.obs.exp.gcggtg+cps.obs.exp.gcggtt+cps.obs.exp.gcgtaa+cps.obs.exp.gcgtac+cps.obs.exp.gcgtag+cps.obs.exp.gcgtat+cps.obs.exp.gcgtca+cps.obs.exp.gcgtcc+cps.obs.exp.gcgtcg+cps.obs.exp.gcgtct+cps.obs.exp.gcgtga+cps.obs.exp.gcgtgc+cps.obs.exp.gcgtgg+cps.obs.exp.gcgtgt+cps.obs.exp.gcgtta+cps.obs.exp.gcgttc+cps.obs.exp.gcgttg+cps.obs.exp.gcgttt+cps.obs.exp.gctaaa+cps.obs.exp.gctaac+cps.obs.exp.gctaag+cps.obs.exp.gctaat+cps.obs.exp.gctaca+cps.obs.exp.gctacc+cps.obs.exp.gctacg+cps.obs.exp.gctact+cps.obs.exp.gctaga+cps.obs.exp.gctagc+cps.obs.exp.gctagg+cps.obs.exp.gctagt+cps.obs.exp.gctata+cps.obs.exp.gctatc+cps.obs.exp.gctatg+cps.obs.exp.gctatt+cps.obs.exp.gctcaa+cps.obs.exp.gctcac+cps.obs.exp.gctcag+cps.obs.exp.gctcat+cps.obs.exp.gctcca+cps.obs.exp.gctccc+cps.obs.exp.gctccg+cps.obs.exp.gctcct+cps.obs.exp.gctcga+cps.obs.exp.gctcgc+cps.obs.exp.gctcgg+cps.obs.exp.gctcgt+cps.obs.exp.gctcta+cps.obs.exp.gctctc+cps.obs.exp.gctctg+cps.obs.exp.gctctt+cps.obs.exp.gctgaa+cps.obs.exp.gctgac+cps.obs.exp.gctgag+cps.obs.exp.gctgat+cps.obs.exp.gctgca+cps.obs.exp.gctgcc+cps.obs.exp.gctgcg+cps.obs.exp.gctgct+cps.obs.exp.gctgga+cps.obs.exp.gctggc+cps.obs.exp.gctggg+cps.obs.exp.gctggt+cps.obs.exp.gctgta+cps.obs.exp.gctgtc+cps.obs.exp.gctgtg+cps.obs.exp.gctgtt+cps.obs.exp.gcttaa+cps.obs.exp.gcttac+cps.obs.exp.gcttag+cps.obs.exp.gcttat+cps.obs.exp.gcttca+cps.obs.exp.gcttcc+cps.obs.exp.gcttcg+cps.obs.exp.gcttct+cps.obs.exp.gcttga+cps.obs.exp.gcttgc+cps.obs.exp.gcttgg+cps.obs.exp.gcttgt+cps.obs.exp.gcttta+cps.obs.exp.gctttc+cps.obs.exp.gctttg+cps.obs.exp.gctttt+cps.obs.exp.ggaaaa+cps.obs.exp.ggaaac+cps.obs.exp.ggaaag+cps.obs.exp.ggaaat+cps.obs.exp.ggaaca+cps.obs.exp.ggaacc+cps.obs.exp.ggaacg+cps.obs.exp.ggaact+cps.obs.exp.ggaaga+cps.obs.exp.ggaagc+cps.obs.exp.ggaagg+cps.obs.exp.ggaagt+cps.obs.exp.ggaata+cps.obs.exp.ggaatc+cps.obs.exp.ggaatg+cps.obs.exp.ggaatt+cps.obs.exp.ggacaa+cps.obs.exp.ggacac+cps.obs.exp.ggacag+cps.obs.exp.ggacat+cps.obs.exp.ggacca+cps.obs.exp.ggaccc+cps.obs.exp.ggaccg+cps.obs.exp.ggacct+cps.obs.exp.ggacga+cps.obs.exp.ggacgc+cps.obs.exp.ggacgg+cps.obs.exp.ggacgt+cps.obs.exp.ggacta+cps.obs.exp.ggactc+cps.obs.exp.ggactg+cps.obs.exp.ggactt+cps.obs.exp.ggagaa+cps.obs.exp.ggagac+cps.obs.exp.ggagag+cps.obs.exp.ggagat+cps.obs.exp.ggagca+cps.obs.exp.ggagcc+cps.obs.exp.ggagcg+cps.obs.exp.ggagct+cps.obs.exp.ggagga+cps.obs.exp.ggaggc+cps.obs.exp.ggaggg+cps.obs.exp.ggaggt+cps.obs.exp.ggagta+cps.obs.exp.ggagtc+cps.obs.exp.ggagtg+cps.obs.exp.ggagtt+cps.obs.exp.ggataa+cps.obs.exp.ggatac+cps.obs.exp.ggatag+cps.obs.exp.ggatat+cps.obs.exp.ggatca+cps.obs.exp.ggatcc+cps.obs.exp.ggatcg+cps.obs.exp.ggatct+cps.obs.exp.ggatga+cps.obs.exp.ggatgc+cps.obs.exp.ggatgg+cps.obs.exp.ggatgt+cps.obs.exp.ggatta+cps.obs.exp.ggattc+cps.obs.exp.ggattg+cps.obs.exp.ggattt+cps.obs.exp.ggcaaa+cps.obs.exp.ggcaac+cps.obs.exp.ggcaag+cps.obs.exp.ggcaat+cps.obs.exp.ggcaca+cps.obs.exp.ggcacc+cps.obs.exp.ggcacg+cps.obs.exp.ggcact+cps.obs.exp.ggcaga+cps.obs.exp.ggcagc+cps.obs.exp.ggcagg+cps.obs.exp.ggcagt+cps.obs.exp.ggcata+cps.obs.exp.ggcatc+cps.obs.exp.ggcatg+cps.obs.exp.ggcatt+cps.obs.exp.ggccaa+cps.obs.exp.ggccac+cps.obs.exp.ggccag+cps.obs.exp.ggccat+cps.obs.exp.ggccca+cps.obs.exp.ggcccc+cps.obs.exp.ggcccg+cps.obs.exp.ggccct+cps.obs.exp.ggccga+cps.obs.exp.ggccgc+cps.obs.exp.ggccgg+cps.obs.exp.ggccgt+cps.obs.exp.ggccta+cps.obs.exp.ggcctc+cps.obs.exp.ggcctg+cps.obs.exp.ggcctt+cps.obs.exp.ggcgaa+cps.obs.exp.ggcgac+cps.obs.exp.ggcgag+cps.obs.exp.ggcgat+cps.obs.exp.ggcgca+cps.obs.exp.ggcgcc+cps.obs.exp.ggcgcg+cps.obs.exp.ggcgct+cps.obs.exp.ggcgga+cps.obs.exp.ggcggc+cps.obs.exp.ggcggg+cps.obs.exp.ggcggt+cps.obs.exp.ggcgta+cps.obs.exp.ggcgtc+cps.obs.exp.ggcgtg+cps.obs.exp.ggcgtt+cps.obs.exp.ggctaa+cps.obs.exp.ggctac+cps.obs.exp.ggctag+cps.obs.exp.ggctat+cps.obs.exp.ggctca+cps.obs.exp.ggctcc+cps.obs.exp.ggctcg+cps.obs.exp.ggctct+cps.obs.exp.ggctga+cps.obs.exp.ggctgc+cps.obs.exp.ggctgg+cps.obs.exp.ggctgt+cps.obs.exp.ggctta+cps.obs.exp.ggcttc+cps.obs.exp.ggcttg+cps.obs.exp.ggcttt+cps.obs.exp.gggaaa+cps.obs.exp.gggaac+cps.obs.exp.gggaag+cps.obs.exp.gggaat+cps.obs.exp.gggaca+cps.obs.exp.gggacc+cps.obs.exp.gggacg+cps.obs.exp.gggact+cps.obs.exp.gggaga+cps.obs.exp.gggagc+cps.obs.exp.gggagg+cps.obs.exp.gggagt+cps.obs.exp.gggata+cps.obs.exp.gggatc+cps.obs.exp.gggatg+cps.obs.exp.gggatt+cps.obs.exp.gggcaa+cps.obs.exp.gggcac+cps.obs.exp.gggcag+cps.obs.exp.gggcat+cps.obs.exp.gggcca+cps.obs.exp.gggccc+cps.obs.exp.gggccg+cps.obs.exp.gggcct+cps.obs.exp.gggcga+cps.obs.exp.gggcgc+cps.obs.exp.gggcgg+cps.obs.exp.gggcgt+cps.obs.exp.gggcta+cps.obs.exp.gggctc+cps.obs.exp.gggctg+cps.obs.exp.gggctt+cps.obs.exp.ggggaa+cps.obs.exp.ggggac+cps.obs.exp.ggggag+cps.obs.exp.ggggat+cps.obs.exp.ggggca+cps.obs.exp.ggggcc+cps.obs.exp.ggggcg+cps.obs.exp.ggggct+cps.obs.exp.ggggga+cps.obs.exp.gggggc+cps.obs.exp.gggggg+cps.obs.exp.gggggt+cps.obs.exp.ggggta+cps.obs.exp.ggggtc+cps.obs.exp.ggggtg+cps.obs.exp.ggggtt+cps.obs.exp.gggtaa+cps.obs.exp.gggtac+cps.obs.exp.gggtag+cps.obs.exp.gggtat+cps.obs.exp.gggtca+cps.obs.exp.gggtcc+cps.obs.exp.gggtcg+cps.obs.exp.gggtct+cps.obs.exp.gggtga+cps.obs.exp.gggtgc+cps.obs.exp.gggtgg+cps.obs.exp.gggtgt+cps.obs.exp.gggtta+cps.obs.exp.gggttc+cps.obs.exp.gggttg+cps.obs.exp.gggttt+cps.obs.exp.ggtaaa+cps.obs.exp.ggtaac+cps.obs.exp.ggtaag+cps.obs.exp.ggtaat+cps.obs.exp.ggtaca+cps.obs.exp.ggtacc+cps.obs.exp.ggtacg+cps.obs.exp.ggtact+cps.obs.exp.ggtaga+cps.obs.exp.ggtagc+cps.obs.exp.ggtagg+cps.obs.exp.ggtagt+cps.obs.exp.ggtata+cps.obs.exp.ggtatc+cps.obs.exp.ggtatg+cps.obs.exp.ggtatt+cps.obs.exp.ggtcaa+cps.obs.exp.ggtcac+cps.obs.exp.ggtcag+cps.obs.exp.ggtcat+cps.obs.exp.ggtcca+cps.obs.exp.ggtccc+cps.obs.exp.ggtccg+cps.obs.exp.ggtcct+cps.obs.exp.ggtcga+cps.obs.exp.ggtcgc+cps.obs.exp.ggtcgg+cps.obs.exp.ggtcgt+cps.obs.exp.ggtcta+cps.obs.exp.ggtctc+cps.obs.exp.ggtctg+cps.obs.exp.ggtctt+cps.obs.exp.ggtgaa+cps.obs.exp.ggtgac+cps.obs.exp.ggtgag+cps.obs.exp.ggtgat+cps.obs.exp.ggtgca+cps.obs.exp.ggtgcc+cps.obs.exp.ggtgcg+cps.obs.exp.ggtgct+cps.obs.exp.ggtgga+cps.obs.exp.ggtggc+cps.obs.exp.ggtggg+cps.obs.exp.ggtggt+cps.obs.exp.ggtgta+cps.obs.exp.ggtgtc+cps.obs.exp.ggtgtg+cps.obs.exp.ggtgtt+cps.obs.exp.ggttaa+cps.obs.exp.ggttac+cps.obs.exp.ggttag+cps.obs.exp.ggttat+cps.obs.exp.ggttca+cps.obs.exp.ggttcc+cps.obs.exp.ggttcg+cps.obs.exp.ggttct+cps.obs.exp.ggttga+cps.obs.exp.ggttgc+cps.obs.exp.ggttgg+cps.obs.exp.ggttgt+cps.obs.exp.ggttta+cps.obs.exp.ggtttc+cps.obs.exp.ggtttg+cps.obs.exp.ggtttt+cps.obs.exp.gtaaaa+cps.obs.exp.gtaaac+cps.obs.exp.gtaaag+cps.obs.exp.gtaaat+cps.obs.exp.gtaaca+cps.obs.exp.gtaacc+cps.obs.exp.gtaacg+cps.obs.exp.gtaact+cps.obs.exp.gtaaga+cps.obs.exp.gtaagc+cps.obs.exp.gtaagg+cps.obs.exp.gtaagt+cps.obs.exp.gtaata+cps.obs.exp.gtaatc+cps.obs.exp.gtaatg+cps.obs.exp.gtaatt+cps.obs.exp.gtacaa+cps.obs.exp.gtacac+cps.obs.exp.gtacag+cps.obs.exp.gtacat+cps.obs.exp.gtacca+cps.obs.exp.gtaccc+cps.obs.exp.gtaccg+cps.obs.exp.gtacct+cps.obs.exp.gtacga+cps.obs.exp.gtacgc+cps.obs.exp.gtacgg+cps.obs.exp.gtacgt+cps.obs.exp.gtacta+cps.obs.exp.gtactc+cps.obs.exp.gtactg+cps.obs.exp.gtactt+cps.obs.exp.gtagaa+cps.obs.exp.gtagac+cps.obs.exp.gtagag+cps.obs.exp.gtagat+cps.obs.exp.gtagca+cps.obs.exp.gtagcc+cps.obs.exp.gtagcg+cps.obs.exp.gtagct+cps.obs.exp.gtagga+cps.obs.exp.gtaggc+cps.obs.exp.gtaggg+cps.obs.exp.gtaggt+cps.obs.exp.gtagta+cps.obs.exp.gtagtc+cps.obs.exp.gtagtg+cps.obs.exp.gtagtt+cps.obs.exp.gtataa+cps.obs.exp.gtatac+cps.obs.exp.gtatag+cps.obs.exp.gtatat+cps.obs.exp.gtatca+cps.obs.exp.gtatcc+cps.obs.exp.gtatcg+cps.obs.exp.gtatct+cps.obs.exp.gtatga+cps.obs.exp.gtatgc+cps.obs.exp.gtatgg+cps.obs.exp.gtatgt+cps.obs.exp.gtatta+cps.obs.exp.gtattc+cps.obs.exp.gtattg+cps.obs.exp.gtattt+cps.obs.exp.gtcaaa+cps.obs.exp.gtcaac+cps.obs.exp.gtcaag+cps.obs.exp.gtcaat+cps.obs.exp.gtcaca+cps.obs.exp.gtcacc+cps.obs.exp.gtcacg+cps.obs.exp.gtcact+cps.obs.exp.gtcaga+cps.obs.exp.gtcagc+cps.obs.exp.gtcagg+cps.obs.exp.gtcagt+cps.obs.exp.gtcata+cps.obs.exp.gtcatc+cps.obs.exp.gtcatg+cps.obs.exp.gtcatt+cps.obs.exp.gtccaa+cps.obs.exp.gtccac+cps.obs.exp.gtccag+cps.obs.exp.gtccat+cps.obs.exp.gtccca+cps.obs.exp.gtcccc+cps.obs.exp.gtcccg+cps.obs.exp.gtccct+cps.obs.exp.gtccga+cps.obs.exp.gtccgc+cps.obs.exp.gtccgg+cps.obs.exp.gtccgt+cps.obs.exp.gtccta+cps.obs.exp.gtcctc+cps.obs.exp.gtcctg+cps.obs.exp.gtcctt+cps.obs.exp.gtcgaa+cps.obs.exp.gtcgac+cps.obs.exp.gtcgag+cps.obs.exp.gtcgat+cps.obs.exp.gtcgca+cps.obs.exp.gtcgcc+cps.obs.exp.gtcgcg+cps.obs.exp.gtcgct+cps.obs.exp.gtcgga+cps.obs.exp.gtcggc+cps.obs.exp.gtcggg+cps.obs.exp.gtcggt+cps.obs.exp.gtcgta+cps.obs.exp.gtcgtc+cps.obs.exp.gtcgtg+cps.obs.exp.gtcgtt+cps.obs.exp.gtctaa+cps.obs.exp.gtctac+cps.obs.exp.gtctag+cps.obs.exp.gtctat+cps.obs.exp.gtctca+cps.obs.exp.gtctcc+cps.obs.exp.gtctcg+cps.obs.exp.gtctct+cps.obs.exp.gtctga+cps.obs.exp.gtctgc+cps.obs.exp.gtctgg+cps.obs.exp.gtctgt+cps.obs.exp.gtctta+cps.obs.exp.gtcttc+cps.obs.exp.gtcttg+cps.obs.exp.gtcttt+cps.obs.exp.gtgaaa+cps.obs.exp.gtgaac+cps.obs.exp.gtgaag+cps.obs.exp.gtgaat+cps.obs.exp.gtgaca+cps.obs.exp.gtgacc+cps.obs.exp.gtgacg+cps.obs.exp.gtgact+cps.obs.exp.gtgaga+cps.obs.exp.gtgagc+cps.obs.exp.gtgagg+cps.obs.exp.gtgagt+cps.obs.exp.gtgata+cps.obs.exp.gtgatc+cps.obs.exp.gtgatg+cps.obs.exp.gtgatt+cps.obs.exp.gtgcaa+cps.obs.exp.gtgcac+cps.obs.exp.gtgcag+cps.obs.exp.gtgcat+cps.obs.exp.gtgcca+cps.obs.exp.gtgccc+cps.obs.exp.gtgccg+cps.obs.exp.gtgcct+cps.obs.exp.gtgcga+cps.obs.exp.gtgcgc+cps.obs.exp.gtgcgg+cps.obs.exp.gtgcgt+cps.obs.exp.gtgcta+cps.obs.exp.gtgctc+cps.obs.exp.gtgctg+cps.obs.exp.gtgctt+cps.obs.exp.gtggaa+cps.obs.exp.gtggac+cps.obs.exp.gtggag+cps.obs.exp.gtggat+cps.obs.exp.gtggca+cps.obs.exp.gtggcc+cps.obs.exp.gtggcg+cps.obs.exp.gtggct+cps.obs.exp.gtggga+cps.obs.exp.gtgggc+cps.obs.exp.gtgggg+cps.obs.exp.gtgggt+cps.obs.exp.gtggta+cps.obs.exp.gtggtc+cps.obs.exp.gtggtg+cps.obs.exp.gtggtt+cps.obs.exp.gtgtaa+cps.obs.exp.gtgtac+cps.obs.exp.gtgtag+cps.obs.exp.gtgtat+cps.obs.exp.gtgtca+cps.obs.exp.gtgtcc+cps.obs.exp.gtgtcg+cps.obs.exp.gtgtct+cps.obs.exp.gtgtga+cps.obs.exp.gtgtgc+cps.obs.exp.gtgtgg+cps.obs.exp.gtgtgt+cps.obs.exp.gtgtta+cps.obs.exp.gtgttc+cps.obs.exp.gtgttg+cps.obs.exp.gtgttt+cps.obs.exp.gttaaa+cps.obs.exp.gttaac+cps.obs.exp.gttaag+cps.obs.exp.gttaat+cps.obs.exp.gttaca+cps.obs.exp.gttacc+cps.obs.exp.gttacg+cps.obs.exp.gttact+cps.obs.exp.gttaga+cps.obs.exp.gttagc+cps.obs.exp.gttagg+cps.obs.exp.gttagt+cps.obs.exp.gttata+cps.obs.exp.gttatc+cps.obs.exp.gttatg+cps.obs.exp.gttatt+cps.obs.exp.gttcaa+cps.obs.exp.gttcac+cps.obs.exp.gttcag+cps.obs.exp.gttcat+cps.obs.exp.gttcca+cps.obs.exp.gttccc+cps.obs.exp.gttccg+cps.obs.exp.gttcct+cps.obs.exp.gttcga+cps.obs.exp.gttcgc+cps.obs.exp.gttcgg+cps.obs.exp.gttcgt+cps.obs.exp.gttcta+cps.obs.exp.gttctc+cps.obs.exp.gttctg+cps.obs.exp.gttctt+cps.obs.exp.gttgaa+cps.obs.exp.gttgac+cps.obs.exp.gttgag+cps.obs.exp.gttgat+cps.obs.exp.gttgca+cps.obs.exp.gttgcc+cps.obs.exp.gttgcg+cps.obs.exp.gttgct+cps.obs.exp.gttgga+cps.obs.exp.gttggc+cps.obs.exp.gttggg+cps.obs.exp.gttggt+cps.obs.exp.gttgta+cps.obs.exp.gttgtc+cps.obs.exp.gttgtg+cps.obs.exp.gttgtt+cps.obs.exp.gtttaa+cps.obs.exp.gtttac+cps.obs.exp.gtttag+cps.obs.exp.gtttat+cps.obs.exp.gtttca+cps.obs.exp.gtttcc+cps.obs.exp.gtttcg+cps.obs.exp.gtttct+cps.obs.exp.gtttga+cps.obs.exp.gtttgc+cps.obs.exp.gtttgg+cps.obs.exp.gtttgt+cps.obs.exp.gtttta+cps.obs.exp.gttttc+cps.obs.exp.gttttg+cps.obs.exp.gttttt+cps.obs.exp.taaaaa+cps.obs.exp.taaaac+cps.obs.exp.taaaag+cps.obs.exp.taaaat+cps.obs.exp.taaaca+cps.obs.exp.taaacc+cps.obs.exp.taaacg+cps.obs.exp.taaact+cps.obs.exp.taaaga+cps.obs.exp.taaagc+cps.obs.exp.taaagg+cps.obs.exp.taaagt+cps.obs.exp.taaata+cps.obs.exp.taaatc+cps.obs.exp.taaatg+cps.obs.exp.taaatt+cps.obs.exp.taacaa+cps.obs.exp.taacac+cps.obs.exp.taacag+cps.obs.exp.taacat+cps.obs.exp.taacca+cps.obs.exp.taaccc+cps.obs.exp.taaccg+cps.obs.exp.taacct+cps.obs.exp.taacga+cps.obs.exp.taacgc+cps.obs.exp.taacgg+cps.obs.exp.taacgt+cps.obs.exp.taacta+cps.obs.exp.taactc+cps.obs.exp.taactg+cps.obs.exp.taactt+cps.obs.exp.taagaa+cps.obs.exp.taagac+cps.obs.exp.taagag+cps.obs.exp.taagat+cps.obs.exp.taagca+cps.obs.exp.taagcc+cps.obs.exp.taagcg+cps.obs.exp.taagct+cps.obs.exp.taagga+cps.obs.exp.taaggc+cps.obs.exp.taaggg+cps.obs.exp.taaggt+cps.obs.exp.taagta+cps.obs.exp.taagtc+cps.obs.exp.taagtg+cps.obs.exp.taagtt+cps.obs.exp.taataa+cps.obs.exp.taatac+cps.obs.exp.taatag+cps.obs.exp.taatat+cps.obs.exp.taatca+cps.obs.exp.taatcc+cps.obs.exp.taatcg+cps.obs.exp.taatct+cps.obs.exp.taatga+cps.obs.exp.taatgc+cps.obs.exp.taatgg+cps.obs.exp.taatgt+cps.obs.exp.taatta+cps.obs.exp.taattc+cps.obs.exp.taattg+cps.obs.exp.taattt+cps.obs.exp.tacaaa+cps.obs.exp.tacaac+cps.obs.exp.tacaag+cps.obs.exp.tacaat+cps.obs.exp.tacaca+cps.obs.exp.tacacc+cps.obs.exp.tacacg+cps.obs.exp.tacact+cps.obs.exp.tacaga+cps.obs.exp.tacagc+cps.obs.exp.tacagg+cps.obs.exp.tacagt+cps.obs.exp.tacata+cps.obs.exp.tacatc+cps.obs.exp.tacatg+cps.obs.exp.tacatt+cps.obs.exp.taccaa+cps.obs.exp.taccac+cps.obs.exp.taccag+cps.obs.exp.taccat+cps.obs.exp.taccca+cps.obs.exp.tacccc+cps.obs.exp.tacccg+cps.obs.exp.taccct+cps.obs.exp.taccga+cps.obs.exp.taccgc+cps.obs.exp.taccgg+cps.obs.exp.taccgt+cps.obs.exp.taccta+cps.obs.exp.tacctc+cps.obs.exp.tacctg+cps.obs.exp.tacctt+cps.obs.exp.tacgaa+cps.obs.exp.tacgac+cps.obs.exp.tacgag+cps.obs.exp.tacgat+cps.obs.exp.tacgca+cps.obs.exp.tacgcc+cps.obs.exp.tacgcg+cps.obs.exp.tacgct+cps.obs.exp.tacgga+cps.obs.exp.tacggc+cps.obs.exp.tacggg+cps.obs.exp.tacggt+cps.obs.exp.tacgta+cps.obs.exp.tacgtc+cps.obs.exp.tacgtg+cps.obs.exp.tacgtt+cps.obs.exp.tactaa+cps.obs.exp.tactac+cps.obs.exp.tactag+cps.obs.exp.tactat+cps.obs.exp.tactca+cps.obs.exp.tactcc+cps.obs.exp.tactcg+cps.obs.exp.tactct+cps.obs.exp.tactga+cps.obs.exp.tactgc+cps.obs.exp.tactgg+cps.obs.exp.tactgt+cps.obs.exp.tactta+cps.obs.exp.tacttc+cps.obs.exp.tacttg+cps.obs.exp.tacttt+cps.obs.exp.tagaaa+cps.obs.exp.tagaac+cps.obs.exp.tagaag+cps.obs.exp.tagaat+cps.obs.exp.tagaca+cps.obs.exp.tagacc+cps.obs.exp.tagacg+cps.obs.exp.tagact+cps.obs.exp.tagaga+cps.obs.exp.tagagc+cps.obs.exp.tagagg+cps.obs.exp.tagagt+cps.obs.exp.tagata+cps.obs.exp.tagatc+cps.obs.exp.tagatg+cps.obs.exp.tagatt+cps.obs.exp.tagcaa+cps.obs.exp.tagcac+cps.obs.exp.tagcag+cps.obs.exp.tagcat+cps.obs.exp.tagcca+cps.obs.exp.tagccc+cps.obs.exp.tagccg+cps.obs.exp.tagcct+cps.obs.exp.tagcga+cps.obs.exp.tagcgc+cps.obs.exp.tagcgg+cps.obs.exp.tagcgt+cps.obs.exp.tagcta+cps.obs.exp.tagctc+cps.obs.exp.tagctg+cps.obs.exp.tagctt+cps.obs.exp.taggaa+cps.obs.exp.taggac+cps.obs.exp.taggag+cps.obs.exp.taggat+cps.obs.exp.taggca+cps.obs.exp.taggcc+cps.obs.exp.taggcg+cps.obs.exp.taggct+cps.obs.exp.taggga+cps.obs.exp.tagggc+cps.obs.exp.tagggg+cps.obs.exp.tagggt+cps.obs.exp.taggta+cps.obs.exp.taggtc+cps.obs.exp.taggtg+cps.obs.exp.taggtt+cps.obs.exp.tagtaa+cps.obs.exp.tagtac+cps.obs.exp.tagtag+cps.obs.exp.tagtat+cps.obs.exp.tagtca+cps.obs.exp.tagtcc+cps.obs.exp.tagtcg+cps.obs.exp.tagtct+cps.obs.exp.tagtga+cps.obs.exp.tagtgc+cps.obs.exp.tagtgg+cps.obs.exp.tagtgt+cps.obs.exp.tagtta+cps.obs.exp.tagttc+cps.obs.exp.tagttg+cps.obs.exp.tagttt+cps.obs.exp.tataaa+cps.obs.exp.tataac+cps.obs.exp.tataag+cps.obs.exp.tataat+cps.obs.exp.tataca+cps.obs.exp.tatacc+cps.obs.exp.tatacg+cps.obs.exp.tatact+cps.obs.exp.tataga+cps.obs.exp.tatagc+cps.obs.exp.tatagg+cps.obs.exp.tatagt+cps.obs.exp.tatata+cps.obs.exp.tatatc+cps.obs.exp.tatatg+cps.obs.exp.tatatt+cps.obs.exp.tatcaa+cps.obs.exp.tatcac+cps.obs.exp.tatcag+cps.obs.exp.tatcat+cps.obs.exp.tatcca+cps.obs.exp.tatccc+cps.obs.exp.tatccg+cps.obs.exp.tatcct+cps.obs.exp.tatcga+cps.obs.exp.tatcgc+cps.obs.exp.tatcgg+cps.obs.exp.tatcgt+cps.obs.exp.tatcta+cps.obs.exp.tatctc+cps.obs.exp.tatctg+cps.obs.exp.tatctt+cps.obs.exp.tatgaa+cps.obs.exp.tatgac+cps.obs.exp.tatgag+cps.obs.exp.tatgat+cps.obs.exp.tatgca+cps.obs.exp.tatgcc+cps.obs.exp.tatgcg+cps.obs.exp.tatgct+cps.obs.exp.tatgga+cps.obs.exp.tatggc+cps.obs.exp.tatggg+cps.obs.exp.tatggt+cps.obs.exp.tatgta+cps.obs.exp.tatgtc+cps.obs.exp.tatgtg+cps.obs.exp.tatgtt+cps.obs.exp.tattaa+cps.obs.exp.tattac+cps.obs.exp.tattag+cps.obs.exp.tattat+cps.obs.exp.tattca+cps.obs.exp.tattcc+cps.obs.exp.tattcg+cps.obs.exp.tattct+cps.obs.exp.tattga+cps.obs.exp.tattgc+cps.obs.exp.tattgg+cps.obs.exp.tattgt+cps.obs.exp.tattta+cps.obs.exp.tatttc+cps.obs.exp.tatttg+cps.obs.exp.tatttt+cps.obs.exp.tcaaaa+cps.obs.exp.tcaaac+cps.obs.exp.tcaaag+cps.obs.exp.tcaaat+cps.obs.exp.tcaaca+cps.obs.exp.tcaacc+cps.obs.exp.tcaacg+cps.obs.exp.tcaact+cps.obs.exp.tcaaga+cps.obs.exp.tcaagc+cps.obs.exp.tcaagg+cps.obs.exp.tcaagt+cps.obs.exp.tcaata+cps.obs.exp.tcaatc+cps.obs.exp.tcaatg+cps.obs.exp.tcaatt+cps.obs.exp.tcacaa+cps.obs.exp.tcacac+cps.obs.exp.tcacag+cps.obs.exp.tcacat+cps.obs.exp.tcacca+cps.obs.exp.tcaccc+cps.obs.exp.tcaccg+cps.obs.exp.tcacct+cps.obs.exp.tcacga+cps.obs.exp.tcacgc+cps.obs.exp.tcacgg+cps.obs.exp.tcacgt+cps.obs.exp.tcacta+cps.obs.exp.tcactc+cps.obs.exp.tcactg+cps.obs.exp.tcactt+cps.obs.exp.tcagaa+cps.obs.exp.tcagac+cps.obs.exp.tcagag+cps.obs.exp.tcagat+cps.obs.exp.tcagca+cps.obs.exp.tcagcc+cps.obs.exp.tcagcg+cps.obs.exp.tcagct+cps.obs.exp.tcagga+cps.obs.exp.tcaggc+cps.obs.exp.tcaggg+cps.obs.exp.tcaggt+cps.obs.exp.tcagta+cps.obs.exp.tcagtc+cps.obs.exp.tcagtg+cps.obs.exp.tcagtt+cps.obs.exp.tcataa+cps.obs.exp.tcatac+cps.obs.exp.tcatag+cps.obs.exp.tcatat+cps.obs.exp.tcatca+cps.obs.exp.tcatcc+cps.obs.exp.tcatcg+cps.obs.exp.tcatct+cps.obs.exp.tcatga+cps.obs.exp.tcatgc+cps.obs.exp.tcatgg+cps.obs.exp.tcatgt+cps.obs.exp.tcatta+cps.obs.exp.tcattc+cps.obs.exp.tcattg+cps.obs.exp.tcattt+cps.obs.exp.tccaaa+cps.obs.exp.tccaac+cps.obs.exp.tccaag+cps.obs.exp.tccaat+cps.obs.exp.tccaca+cps.obs.exp.tccacc+cps.obs.exp.tccacg+cps.obs.exp.tccact+cps.obs.exp.tccaga+cps.obs.exp.tccagc+cps.obs.exp.tccagg+cps.obs.exp.tccagt+cps.obs.exp.tccata+cps.obs.exp.tccatc+cps.obs.exp.tccatg+cps.obs.exp.tccatt+cps.obs.exp.tcccaa+cps.obs.exp.tcccac+cps.obs.exp.tcccag+cps.obs.exp.tcccat+cps.obs.exp.tcccca+cps.obs.exp.tccccc+cps.obs.exp.tccccg+cps.obs.exp.tcccct+cps.obs.exp.tcccga+cps.obs.exp.tcccgc+cps.obs.exp.tcccgg+cps.obs.exp.tcccgt+cps.obs.exp.tcccta+cps.obs.exp.tccctc+cps.obs.exp.tccctg+cps.obs.exp.tccctt+cps.obs.exp.tccgaa+cps.obs.exp.tccgac+cps.obs.exp.tccgag+cps.obs.exp.tccgat+cps.obs.exp.tccgca+cps.obs.exp.tccgcc+cps.obs.exp.tccgcg+cps.obs.exp.tccgct+cps.obs.exp.tccgga+cps.obs.exp.tccggc+cps.obs.exp.tccggg+cps.obs.exp.tccggt+cps.obs.exp.tccgta+cps.obs.exp.tccgtc+cps.obs.exp.tccgtg+cps.obs.exp.tccgtt+cps.obs.exp.tcctaa+cps.obs.exp.tcctac+cps.obs.exp.tcctag+cps.obs.exp.tcctat+cps.obs.exp.tcctca+cps.obs.exp.tcctcc+cps.obs.exp.tcctcg+cps.obs.exp.tcctct+cps.obs.exp.tcctga+cps.obs.exp.tcctgc+cps.obs.exp.tcctgg+cps.obs.exp.tcctgt+cps.obs.exp.tcctta+cps.obs.exp.tccttc+cps.obs.exp.tccttg+cps.obs.exp.tccttt+cps.obs.exp.tcgaaa+cps.obs.exp.tcgaac+cps.obs.exp.tcgaag+cps.obs.exp.tcgaat+cps.obs.exp.tcgaca+cps.obs.exp.tcgacc+cps.obs.exp.tcgacg+cps.obs.exp.tcgact+cps.obs.exp.tcgaga+cps.obs.exp.tcgagc+cps.obs.exp.tcgagg+cps.obs.exp.tcgagt+cps.obs.exp.tcgata+cps.obs.exp.tcgatc+cps.obs.exp.tcgatg+cps.obs.exp.tcgatt+cps.obs.exp.tcgcaa+cps.obs.exp.tcgcac+cps.obs.exp.tcgcag+cps.obs.exp.tcgcat+cps.obs.exp.tcgcca+cps.obs.exp.tcgccc+cps.obs.exp.tcgccg+cps.obs.exp.tcgcct+cps.obs.exp.tcgcga+cps.obs.exp.tcgcgc+cps.obs.exp.tcgcgg+cps.obs.exp.tcgcgt+cps.obs.exp.tcgcta+cps.obs.exp.tcgctc+cps.obs.exp.tcgctg+cps.obs.exp.tcgctt+cps.obs.exp.tcggaa+cps.obs.exp.tcggac+cps.obs.exp.tcggag+cps.obs.exp.tcggat+cps.obs.exp.tcggca+cps.obs.exp.tcggcc+cps.obs.exp.tcggcg+cps.obs.exp.tcggct+cps.obs.exp.tcggga+cps.obs.exp.tcgggc+cps.obs.exp.tcgggg+cps.obs.exp.tcgggt+cps.obs.exp.tcggta+cps.obs.exp.tcggtc+cps.obs.exp.tcggtg+cps.obs.exp.tcggtt+cps.obs.exp.tcgtaa+cps.obs.exp.tcgtac+cps.obs.exp.tcgtag+cps.obs.exp.tcgtat+cps.obs.exp.tcgtca+cps.obs.exp.tcgtcc+cps.obs.exp.tcgtcg+cps.obs.exp.tcgtct+cps.obs.exp.tcgtga+cps.obs.exp.tcgtgc+cps.obs.exp.tcgtgg+cps.obs.exp.tcgtgt+cps.obs.exp.tcgtta+cps.obs.exp.tcgttc+cps.obs.exp.tcgttg+cps.obs.exp.tcgttt+cps.obs.exp.tctaaa+cps.obs.exp.tctaac+cps.obs.exp.tctaag+cps.obs.exp.tctaat+cps.obs.exp.tctaca+cps.obs.exp.tctacc+cps.obs.exp.tctacg+cps.obs.exp.tctact+cps.obs.exp.tctaga+cps.obs.exp.tctagc+cps.obs.exp.tctagg+cps.obs.exp.tctagt+cps.obs.exp.tctata+cps.obs.exp.tctatc+cps.obs.exp.tctatg+cps.obs.exp.tctatt+cps.obs.exp.tctcaa+cps.obs.exp.tctcac+cps.obs.exp.tctcag+cps.obs.exp.tctcat+cps.obs.exp.tctcca+cps.obs.exp.tctccc+cps.obs.exp.tctccg+cps.obs.exp.tctcct+cps.obs.exp.tctcga+cps.obs.exp.tctcgc+cps.obs.exp.tctcgg+cps.obs.exp.tctcgt+cps.obs.exp.tctcta+cps.obs.exp.tctctc+cps.obs.exp.tctctg+cps.obs.exp.tctctt+cps.obs.exp.tctgaa+cps.obs.exp.tctgac+cps.obs.exp.tctgag+cps.obs.exp.tctgat+cps.obs.exp.tctgca+cps.obs.exp.tctgcc+cps.obs.exp.tctgcg+cps.obs.exp.tctgct+cps.obs.exp.tctgga+cps.obs.exp.tctggc+cps.obs.exp.tctggg+cps.obs.exp.tctggt+cps.obs.exp.tctgta+cps.obs.exp.tctgtc+cps.obs.exp.tctgtg+cps.obs.exp.tctgtt+cps.obs.exp.tcttaa+cps.obs.exp.tcttac+cps.obs.exp.tcttag+cps.obs.exp.tcttat+cps.obs.exp.tcttca+cps.obs.exp.tcttcc+cps.obs.exp.tcttcg+cps.obs.exp.tcttct+cps.obs.exp.tcttga+cps.obs.exp.tcttgc+cps.obs.exp.tcttgg+cps.obs.exp.tcttgt+cps.obs.exp.tcttta+cps.obs.exp.tctttc+cps.obs.exp.tctttg+cps.obs.exp.tctttt+cps.obs.exp.tgaaaa+cps.obs.exp.tgaaac+cps.obs.exp.tgaaag+cps.obs.exp.tgaaat+cps.obs.exp.tgaaca+cps.obs.exp.tgaacc+cps.obs.exp.tgaacg+cps.obs.exp.tgaact+cps.obs.exp.tgaaga+cps.obs.exp.tgaagc+cps.obs.exp.tgaagg+cps.obs.exp.tgaagt+cps.obs.exp.tgaata+cps.obs.exp.tgaatc+cps.obs.exp.tgaatg+cps.obs.exp.tgaatt+cps.obs.exp.tgacaa+cps.obs.exp.tgacac+cps.obs.exp.tgacag+cps.obs.exp.tgacat+cps.obs.exp.tgacca+cps.obs.exp.tgaccc+cps.obs.exp.tgaccg+cps.obs.exp.tgacct+cps.obs.exp.tgacga+cps.obs.exp.tgacgc+cps.obs.exp.tgacgg+cps.obs.exp.tgacgt+cps.obs.exp.tgacta+cps.obs.exp.tgactc+cps.obs.exp.tgactg+cps.obs.exp.tgactt+cps.obs.exp.tgagaa+cps.obs.exp.tgagac+cps.obs.exp.tgagag+cps.obs.exp.tgagat+cps.obs.exp.tgagca+cps.obs.exp.tgagcc+cps.obs.exp.tgagcg+cps.obs.exp.tgagct+cps.obs.exp.tgagga+cps.obs.exp.tgaggc+cps.obs.exp.tgaggg+cps.obs.exp.tgaggt+cps.obs.exp.tgagta+cps.obs.exp.tgagtc+cps.obs.exp.tgagtg+cps.obs.exp.tgagtt+cps.obs.exp.tgataa+cps.obs.exp.tgatac+cps.obs.exp.tgatag+cps.obs.exp.tgatat+cps.obs.exp.tgatca+cps.obs.exp.tgatcc+cps.obs.exp.tgatcg+cps.obs.exp.tgatct+cps.obs.exp.tgatga+cps.obs.exp.tgatgc+cps.obs.exp.tgatgg+cps.obs.exp.tgatgt+cps.obs.exp.tgatta+cps.obs.exp.tgattc+cps.obs.exp.tgattg+cps.obs.exp.tgattt+cps.obs.exp.tgcaaa+cps.obs.exp.tgcaac+cps.obs.exp.tgcaag+cps.obs.exp.tgcaat+cps.obs.exp.tgcaca+cps.obs.exp.tgcacc+cps.obs.exp.tgcacg+cps.obs.exp.tgcact+cps.obs.exp.tgcaga+cps.obs.exp.tgcagc+cps.obs.exp.tgcagg+cps.obs.exp.tgcagt+cps.obs.exp.tgcata+cps.obs.exp.tgcatc+cps.obs.exp.tgcatg+cps.obs.exp.tgcatt+cps.obs.exp.tgccaa+cps.obs.exp.tgccac+cps.obs.exp.tgccag+cps.obs.exp.tgccat+cps.obs.exp.tgccca+cps.obs.exp.tgcccc+cps.obs.exp.tgcccg+cps.obs.exp.tgccct+cps.obs.exp.tgccga+cps.obs.exp.tgccgc+cps.obs.exp.tgccgg+cps.obs.exp.tgccgt+cps.obs.exp.tgccta+cps.obs.exp.tgcctc+cps.obs.exp.tgcctg+cps.obs.exp.tgcctt+cps.obs.exp.tgcgaa+cps.obs.exp.tgcgac+cps.obs.exp.tgcgag+cps.obs.exp.tgcgat+cps.obs.exp.tgcgca+cps.obs.exp.tgcgcc+cps.obs.exp.tgcgcg+cps.obs.exp.tgcgct+cps.obs.exp.tgcgga+cps.obs.exp.tgcggc+cps.obs.exp.tgcggg+cps.obs.exp.tgcggt+cps.obs.exp.tgcgta+cps.obs.exp.tgcgtc+cps.obs.exp.tgcgtg+cps.obs.exp.tgcgtt+cps.obs.exp.tgctaa+cps.obs.exp.tgctac+cps.obs.exp.tgctag+cps.obs.exp.tgctat+cps.obs.exp.tgctca+cps.obs.exp.tgctcc+cps.obs.exp.tgctcg+cps.obs.exp.tgctct+cps.obs.exp.tgctga+cps.obs.exp.tgctgc+cps.obs.exp.tgctgg+cps.obs.exp.tgctgt+cps.obs.exp.tgctta+cps.obs.exp.tgcttc+cps.obs.exp.tgcttg+cps.obs.exp.tgcttt+cps.obs.exp.tggaaa+cps.obs.exp.tggaac+cps.obs.exp.tggaag+cps.obs.exp.tggaat+cps.obs.exp.tggaca+cps.obs.exp.tggacc+cps.obs.exp.tggacg+cps.obs.exp.tggact+cps.obs.exp.tggaga+cps.obs.exp.tggagc+cps.obs.exp.tggagg+cps.obs.exp.tggagt+cps.obs.exp.tggata+cps.obs.exp.tggatc+cps.obs.exp.tggatg+cps.obs.exp.tggatt+cps.obs.exp.tggcaa+cps.obs.exp.tggcac+cps.obs.exp.tggcag+cps.obs.exp.tggcat+cps.obs.exp.tggcca+cps.obs.exp.tggccc+cps.obs.exp.tggccg+cps.obs.exp.tggcct+cps.obs.exp.tggcga+cps.obs.exp.tggcgc+cps.obs.exp.tggcgg+cps.obs.exp.tggcgt+cps.obs.exp.tggcta+cps.obs.exp.tggctc+cps.obs.exp.tggctg+cps.obs.exp.tggctt+cps.obs.exp.tgggaa+cps.obs.exp.tgggac+cps.obs.exp.tgggag+cps.obs.exp.tgggat+cps.obs.exp.tgggca+cps.obs.exp.tgggcc+cps.obs.exp.tgggcg+cps.obs.exp.tgggct+cps.obs.exp.tgggga+cps.obs.exp.tggggc+cps.obs.exp.tggggg+cps.obs.exp.tggggt+cps.obs.exp.tgggta+cps.obs.exp.tgggtc+cps.obs.exp.tgggtg+cps.obs.exp.tgggtt+cps.obs.exp.tggtaa+cps.obs.exp.tggtac+cps.obs.exp.tggtag+cps.obs.exp.tggtat+cps.obs.exp.tggtca+cps.obs.exp.tggtcc+cps.obs.exp.tggtcg+cps.obs.exp.tggtct+cps.obs.exp.tggtga+cps.obs.exp.tggtgc+cps.obs.exp.tggtgg+cps.obs.exp.tggtgt+cps.obs.exp.tggtta+cps.obs.exp.tggttc+cps.obs.exp.tggttg+cps.obs.exp.tggttt+cps.obs.exp.tgtaaa+cps.obs.exp.tgtaac+cps.obs.exp.tgtaag+cps.obs.exp.tgtaat+cps.obs.exp.tgtaca+cps.obs.exp.tgtacc+cps.obs.exp.tgtacg+cps.obs.exp.tgtact+cps.obs.exp.tgtaga+cps.obs.exp.tgtagc+cps.obs.exp.tgtagg+cps.obs.exp.tgtagt+cps.obs.exp.tgtata+cps.obs.exp.tgtatc+cps.obs.exp.tgtatg+cps.obs.exp.tgtatt+cps.obs.exp.tgtcaa+cps.obs.exp.tgtcac+cps.obs.exp.tgtcag+cps.obs.exp.tgtcat+cps.obs.exp.tgtcca+cps.obs.exp.tgtccc+cps.obs.exp.tgtccg+cps.obs.exp.tgtcct+cps.obs.exp.tgtcga+cps.obs.exp.tgtcgc+cps.obs.exp.tgtcgg+cps.obs.exp.tgtcgt+cps.obs.exp.tgtcta+cps.obs.exp.tgtctc+cps.obs.exp.tgtctg+cps.obs.exp.tgtctt+cps.obs.exp.tgtgaa+cps.obs.exp.tgtgac+cps.obs.exp.tgtgag+cps.obs.exp.tgtgat+cps.obs.exp.tgtgca+cps.obs.exp.tgtgcc+cps.obs.exp.tgtgcg+cps.obs.exp.tgtgct+cps.obs.exp.tgtgga+cps.obs.exp.tgtggc+cps.obs.exp.tgtggg+cps.obs.exp.tgtggt+cps.obs.exp.tgtgta+cps.obs.exp.tgtgtc+cps.obs.exp.tgtgtg+cps.obs.exp.tgtgtt+cps.obs.exp.tgttaa+cps.obs.exp.tgttac+cps.obs.exp.tgttag+cps.obs.exp.tgttat+cps.obs.exp.tgttca+cps.obs.exp.tgttcc+cps.obs.exp.tgttcg+cps.obs.exp.tgttct+cps.obs.exp.tgttga+cps.obs.exp.tgttgc+cps.obs.exp.tgttgg+cps.obs.exp.tgttgt+cps.obs.exp.tgttta+cps.obs.exp.tgtttc+cps.obs.exp.tgtttg+cps.obs.exp.tgtttt+cps.obs.exp.ttaaaa+cps.obs.exp.ttaaac+cps.obs.exp.ttaaag+cps.obs.exp.ttaaat+cps.obs.exp.ttaaca+cps.obs.exp.ttaacc+cps.obs.exp.ttaacg+cps.obs.exp.ttaact+cps.obs.exp.ttaaga+cps.obs.exp.ttaagc+cps.obs.exp.ttaagg+cps.obs.exp.ttaagt+cps.obs.exp.ttaata+cps.obs.exp.ttaatc+cps.obs.exp.ttaatg+cps.obs.exp.ttaatt+cps.obs.exp.ttacaa+cps.obs.exp.ttacac+cps.obs.exp.ttacag+cps.obs.exp.ttacat+cps.obs.exp.ttacca+cps.obs.exp.ttaccc+cps.obs.exp.ttaccg+cps.obs.exp.ttacct+cps.obs.exp.ttacga+cps.obs.exp.ttacgc+cps.obs.exp.ttacgg+cps.obs.exp.ttacgt+cps.obs.exp.ttacta+cps.obs.exp.ttactc+cps.obs.exp.ttactg+cps.obs.exp.ttactt+cps.obs.exp.ttagaa+cps.obs.exp.ttagac+cps.obs.exp.ttagag+cps.obs.exp.ttagat+cps.obs.exp.ttagca+cps.obs.exp.ttagcc+cps.obs.exp.ttagcg+cps.obs.exp.ttagct+cps.obs.exp.ttagga+cps.obs.exp.ttaggc+cps.obs.exp.ttaggg+cps.obs.exp.ttaggt+cps.obs.exp.ttagta+cps.obs.exp.ttagtc+cps.obs.exp.ttagtg+cps.obs.exp.ttagtt+cps.obs.exp.ttataa+cps.obs.exp.ttatac+cps.obs.exp.ttatag+cps.obs.exp.ttatat+cps.obs.exp.ttatca+cps.obs.exp.ttatcc+cps.obs.exp.ttatcg+cps.obs.exp.ttatct+cps.obs.exp.ttatga+cps.obs.exp.ttatgc+cps.obs.exp.ttatgg+cps.obs.exp.ttatgt+cps.obs.exp.ttatta+cps.obs.exp.ttattc+cps.obs.exp.ttattg+cps.obs.exp.ttattt+cps.obs.exp.ttcaaa+cps.obs.exp.ttcaac+cps.obs.exp.ttcaag+cps.obs.exp.ttcaat+cps.obs.exp.ttcaca+cps.obs.exp.ttcacc+cps.obs.exp.ttcacg+cps.obs.exp.ttcact+cps.obs.exp.ttcaga+cps.obs.exp.ttcagc+cps.obs.exp.ttcagg+cps.obs.exp.ttcagt+cps.obs.exp.ttcata+cps.obs.exp.ttcatc+cps.obs.exp.ttcatg+cps.obs.exp.ttcatt+cps.obs.exp.ttccaa+cps.obs.exp.ttccac+cps.obs.exp.ttccag+cps.obs.exp.ttccat+cps.obs.exp.ttccca+cps.obs.exp.ttcccc+cps.obs.exp.ttcccg+cps.obs.exp.ttccct+cps.obs.exp.ttccga+cps.obs.exp.ttccgc+cps.obs.exp.ttccgg+cps.obs.exp.ttccgt+cps.obs.exp.ttccta+cps.obs.exp.ttcctc+cps.obs.exp.ttcctg+cps.obs.exp.ttcctt+cps.obs.exp.ttcgaa+cps.obs.exp.ttcgac+cps.obs.exp.ttcgag+cps.obs.exp.ttcgat+cps.obs.exp.ttcgca+cps.obs.exp.ttcgcc+cps.obs.exp.ttcgcg+cps.obs.exp.ttcgct+cps.obs.exp.ttcgga+cps.obs.exp.ttcggc+cps.obs.exp.ttcggg+cps.obs.exp.ttcggt+cps.obs.exp.ttcgta+cps.obs.exp.ttcgtc+cps.obs.exp.ttcgtg+cps.obs.exp.ttcgtt+cps.obs.exp.ttctaa+cps.obs.exp.ttctac+cps.obs.exp.ttctag+cps.obs.exp.ttctat+cps.obs.exp.ttctca+cps.obs.exp.ttctcc+cps.obs.exp.ttctcg+cps.obs.exp.ttctct+cps.obs.exp.ttctga+cps.obs.exp.ttctgc+cps.obs.exp.ttctgg+cps.obs.exp.ttctgt+cps.obs.exp.ttctta+cps.obs.exp.ttcttc+cps.obs.exp.ttcttg+cps.obs.exp.ttcttt+cps.obs.exp.ttgaaa+cps.obs.exp.ttgaac+cps.obs.exp.ttgaag+cps.obs.exp.ttgaat+cps.obs.exp.ttgaca+cps.obs.exp.ttgacc+cps.obs.exp.ttgacg+cps.obs.exp.ttgact+cps.obs.exp.ttgaga+cps.obs.exp.ttgagc+cps.obs.exp.ttgagg+cps.obs.exp.ttgagt+cps.obs.exp.ttgata+cps.obs.exp.ttgatc+cps.obs.exp.ttgatg+cps.obs.exp.ttgatt+cps.obs.exp.ttgcaa+cps.obs.exp.ttgcac+cps.obs.exp.ttgcag+cps.obs.exp.ttgcat+cps.obs.exp.ttgcca+cps.obs.exp.ttgccc+cps.obs.exp.ttgccg+cps.obs.exp.ttgcct+cps.obs.exp.ttgcga+cps.obs.exp.ttgcgc+cps.obs.exp.ttgcgg+cps.obs.exp.ttgcgt+cps.obs.exp.ttgcta+cps.obs.exp.ttgctc+cps.obs.exp.ttgctg+cps.obs.exp.ttgctt+cps.obs.exp.ttggaa+cps.obs.exp.ttggac+cps.obs.exp.ttggag+cps.obs.exp.ttggat+cps.obs.exp.ttggca+cps.obs.exp.ttggcc+cps.obs.exp.ttggcg+cps.obs.exp.ttggct+cps.obs.exp.ttggga+cps.obs.exp.ttgggc+cps.obs.exp.ttgggg+cps.obs.exp.ttgggt+cps.obs.exp.ttggta+cps.obs.exp.ttggtc+cps.obs.exp.ttggtg+cps.obs.exp.ttggtt+cps.obs.exp.ttgtaa+cps.obs.exp.ttgtac+cps.obs.exp.ttgtag+cps.obs.exp.ttgtat+cps.obs.exp.ttgtca+cps.obs.exp.ttgtcc+cps.obs.exp.ttgtcg+cps.obs.exp.ttgtct+cps.obs.exp.ttgtga+cps.obs.exp.ttgtgc+cps.obs.exp.ttgtgg+cps.obs.exp.ttgtgt+cps.obs.exp.ttgtta+cps.obs.exp.ttgttc+cps.obs.exp.ttgttg+cps.obs.exp.ttgttt+cps.obs.exp.tttaaa+cps.obs.exp.tttaac+cps.obs.exp.tttaag+cps.obs.exp.tttaat+cps.obs.exp.tttaca+cps.obs.exp.tttacc+cps.obs.exp.tttacg+cps.obs.exp.tttact+cps.obs.exp.tttaga+cps.obs.exp.tttagc+cps.obs.exp.tttagg+cps.obs.exp.tttagt+cps.obs.exp.tttata+cps.obs.exp.tttatc+cps.obs.exp.tttatg+cps.obs.exp.tttatt+cps.obs.exp.tttcaa+cps.obs.exp.tttcac+cps.obs.exp.tttcag+cps.obs.exp.tttcat+cps.obs.exp.tttcca+cps.obs.exp.tttccc+cps.obs.exp.tttccg+cps.obs.exp.tttcct+cps.obs.exp.tttcga+cps.obs.exp.tttcgc+cps.obs.exp.tttcgg+cps.obs.exp.tttcgt+cps.obs.exp.tttcta+cps.obs.exp.tttctc+cps.obs.exp.tttctg+cps.obs.exp.tttctt+cps.obs.exp.tttgaa+cps.obs.exp.tttgac+cps.obs.exp.tttgag+cps.obs.exp.tttgat+cps.obs.exp.tttgca+cps.obs.exp.tttgcc+cps.obs.exp.tttgcg+cps.obs.exp.tttgct+cps.obs.exp.tttgga+cps.obs.exp.tttggc+cps.obs.exp.tttggg+cps.obs.exp.tttggt+cps.obs.exp.tttgta+cps.obs.exp.tttgtc+cps.obs.exp.tttgtg+cps.obs.exp.tttgtt+cps.obs.exp.ttttaa+cps.obs.exp.ttttac+cps.obs.exp.ttttag+cps.obs.exp.ttttat+cps.obs.exp.ttttca+cps.obs.exp.ttttcc+cps.obs.exp.ttttcg+cps.obs.exp.ttttct+cps.obs.exp.ttttga+cps.obs.exp.ttttgc+cps.obs.exp.ttttgg+cps.obs.exp.ttttgt+cps.obs.exp.ttttta+cps.obs.exp.tttttc+cps.obs.exp.tttttg+cps.obs.exp.tttttt
 
#cpb.n3<-cps.obs.exp.aaaaaa+cps.obs.exp.aaaaac+cps.obs.exp.aaaaag+cps.obs.exp.aaaaat+cps.obs.exp.aaaaca+cps.obs.exp.aaaacc+cps.obs.exp.aaaacg+cps.obs.exp.aaaact+cps.obs.exp.aaaaga+cps.obs.exp.aaaagc+cps.obs.exp.aaaagg+cps.obs.exp.aaaagt+cps.obs.exp.aaaata+cps.obs.exp.aaaatc+cps.obs.exp.aaaatg+cps.obs.exp.aaaatt+cps.obs.exp.aaacaa+cps.obs.exp.aaacac+cps.obs.exp.aaacag+cps.obs.exp.aaacat+cps.obs.exp.aaacca+cps.obs.exp.aaaccc+cps.obs.exp.aaaccg+cps.obs.exp.aaacct+cps.obs.exp.aaacga+cps.obs.exp.aaacgc+cps.obs.exp.aaacgg+cps.obs.exp.aaacgt+cps.obs.exp.aaacta+cps.obs.exp.aaactc+cps.obs.exp.aaactg+cps.obs.exp.aaactt+cps.obs.exp.aaagaa+cps.obs.exp.aaagac+cps.obs.exp.aaagag+cps.obs.exp.aaagat+cps.obs.exp.aaagca+cps.obs.exp.aaagcc+cps.obs.exp.aaagcg+cps.obs.exp.aaagct+cps.obs.exp.aaagga+cps.obs.exp.aaaggc+cps.obs.exp.aaaggg+cps.obs.exp.aaaggt+cps.obs.exp.aaagta+cps.obs.exp.aaagtc+cps.obs.exp.aaagtg+cps.obs.exp.aaagtt+cps.obs.exp.aaatac+cps.obs.exp.aaatat+cps.obs.exp.aaatca+cps.obs.exp.aaatcc+cps.obs.exp.aaatcg+cps.obs.exp.aaatct+cps.obs.exp.aaatgc+cps.obs.exp.aaatgg+cps.obs.exp.aaatgt+cps.obs.exp.aaatta+cps.obs.exp.aaattc+cps.obs.exp.aaattg+cps.obs.exp.aaattt+cps.obs.exp.aacaaa+cps.obs.exp.aacaac+cps.obs.exp.aacaag+cps.obs.exp.aacaat+cps.obs.exp.aacaca+cps.obs.exp.aacacc+cps.obs.exp.aacacg+cps.obs.exp.aacact+cps.obs.exp.aacaga+cps.obs.exp.aacagc+cps.obs.exp.aacagg+cps.obs.exp.aacagt+cps.obs.exp.aacata+cps.obs.exp.aacatc+cps.obs.exp.aacatg+cps.obs.exp.aacatt+cps.obs.exp.aaccaa+cps.obs.exp.aaccac+cps.obs.exp.aaccag+cps.obs.exp.aaccat+cps.obs.exp.aaccca+cps.obs.exp.aacccc+cps.obs.exp.aacccg+cps.obs.exp.aaccct+cps.obs.exp.aaccga+cps.obs.exp.aaccgc+cps.obs.exp.aaccgg+cps.obs.exp.aaccgt+cps.obs.exp.aaccta+cps.obs.exp.aacctc+cps.obs.exp.aacctg+cps.obs.exp.aacctt+cps.obs.exp.aacgaa+cps.obs.exp.aacgac+cps.obs.exp.aacgag+cps.obs.exp.aacgat+cps.obs.exp.aacgca+cps.obs.exp.aacgcc+cps.obs.exp.aacgcg+cps.obs.exp.aacgct+cps.obs.exp.aacgga+cps.obs.exp.aacggc+cps.obs.exp.aacggg+cps.obs.exp.aacggt+cps.obs.exp.aacgta+cps.obs.exp.aacgtc+cps.obs.exp.aacgtg+cps.obs.exp.aacgtt+cps.obs.exp.aactac+cps.obs.exp.aactat+cps.obs.exp.aactca+cps.obs.exp.aactcc+cps.obs.exp.aactcg+cps.obs.exp.aactct+cps.obs.exp.aactgc+cps.obs.exp.aactgg+cps.obs.exp.aactgt+cps.obs.exp.aactta+cps.obs.exp.aacttc+cps.obs.exp.aacttg+cps.obs.exp.aacttt+cps.obs.exp.aagaaa+cps.obs.exp.aagaac+cps.obs.exp.aagaag+cps.obs.exp.aagaat+cps.obs.exp.aagaca+cps.obs.exp.aagacc+cps.obs.exp.aagacg+cps.obs.exp.aagact+cps.obs.exp.aagaga+cps.obs.exp.aagagc+cps.obs.exp.aagagg+cps.obs.exp.aagagt+cps.obs.exp.aagata+cps.obs.exp.aagatc+cps.obs.exp.aagatg+cps.obs.exp.aagatt+cps.obs.exp.aagcaa+cps.obs.exp.aagcac+cps.obs.exp.aagcag+cps.obs.exp.aagcat+cps.obs.exp.aagcca+cps.obs.exp.aagccc+cps.obs.exp.aagccg+cps.obs.exp.aagcct+cps.obs.exp.aagcga+cps.obs.exp.aagcgc+cps.obs.exp.aagcgg+cps.obs.exp.aagcgt+cps.obs.exp.aagcta+cps.obs.exp.aagctc+cps.obs.exp.aagctg+cps.obs.exp.aagctt+cps.obs.exp.aaggaa+cps.obs.exp.aaggac+cps.obs.exp.aaggag+cps.obs.exp.aaggat+cps.obs.exp.aaggca+cps.obs.exp.aaggcc+cps.obs.exp.aaggcg+cps.obs.exp.aaggct+cps.obs.exp.aaggga+cps.obs.exp.aagggc+cps.obs.exp.aagggg+cps.obs.exp.aagggt+cps.obs.exp.aaggta+cps.obs.exp.aaggtc+cps.obs.exp.aaggtg+cps.obs.exp.aaggtt+cps.obs.exp.aagtac+cps.obs.exp.aagtat+cps.obs.exp.aagtca+cps.obs.exp.aagtcc+cps.obs.exp.aagtcg+cps.obs.exp.aagtct+cps.obs.exp.aagtgc+cps.obs.exp.aagtgg+cps.obs.exp.aagtgt+cps.obs.exp.aagtta+cps.obs.exp.aagttc+cps.obs.exp.aagttg+cps.obs.exp.aagttt+cps.obs.exp.aataaa+cps.obs.exp.aataac+cps.obs.exp.aataag+cps.obs.exp.aataat+cps.obs.exp.aataca+cps.obs.exp.aatacc+cps.obs.exp.aatacg+cps.obs.exp.aatact+cps.obs.exp.aataga+cps.obs.exp.aatagc+cps.obs.exp.aatagg+cps.obs.exp.aatagt+cps.obs.exp.aatata+cps.obs.exp.aatatc+cps.obs.exp.aatatg+cps.obs.exp.aatatt+cps.obs.exp.aatcaa+cps.obs.exp.aatcac+cps.obs.exp.aatcag+cps.obs.exp.aatcat+cps.obs.exp.aatcca+cps.obs.exp.aatccc+cps.obs.exp.aatccg+cps.obs.exp.aatcct+cps.obs.exp.aatcga+cps.obs.exp.aatcgc+cps.obs.exp.aatcgg+cps.obs.exp.aatcgt+cps.obs.exp.aatcta+cps.obs.exp.aatctc+cps.obs.exp.aatctg+cps.obs.exp.aatctt+cps.obs.exp.aatgaa+cps.obs.exp.aatgac+cps.obs.exp.aatgag+cps.obs.exp.aatgat+cps.obs.exp.aatgca+cps.obs.exp.aatgcc+cps.obs.exp.aatgcg+cps.obs.exp.aatgct+cps.obs.exp.aatgga+cps.obs.exp.aatggc+cps.obs.exp.aatggg+cps.obs.exp.aatggt+cps.obs.exp.aatgta+cps.obs.exp.aatgtc+cps.obs.exp.aatgtg+cps.obs.exp.aatgtt+cps.obs.exp.aattac+cps.obs.exp.aattat+cps.obs.exp.aattca+cps.obs.exp.aattcc+cps.obs.exp.aattcg+cps.obs.exp.aattct+cps.obs.exp.aattgc+cps.obs.exp.aattgg+cps.obs.exp.aattgt+cps.obs.exp.aattta+cps.obs.exp.aatttc+cps.obs.exp.aatttg+cps.obs.exp.aatttt+cps.obs.exp.acaaaa+cps.obs.exp.acaaac+cps.obs.exp.acaaag+cps.obs.exp.acaaat+cps.obs.exp.acaaca+cps.obs.exp.acaacc+cps.obs.exp.acaacg+cps.obs.exp.acaact+cps.obs.exp.acaaga+cps.obs.exp.acaagc+cps.obs.exp.acaagg+cps.obs.exp.acaagt+cps.obs.exp.acaata+cps.obs.exp.acaatc+cps.obs.exp.acaatg+cps.obs.exp.acaatt+cps.obs.exp.acacaa+cps.obs.exp.acacac+cps.obs.exp.acacag+cps.obs.exp.acacat+cps.obs.exp.acacca+cps.obs.exp.acaccc+cps.obs.exp.acaccg+cps.obs.exp.acacct+cps.obs.exp.acacga+cps.obs.exp.acacgc+cps.obs.exp.acacgg+cps.obs.exp.acacgt+cps.obs.exp.acacta+cps.obs.exp.acactc+cps.obs.exp.acactg+cps.obs.exp.acactt+cps.obs.exp.acagaa+cps.obs.exp.acagac+cps.obs.exp.acagag+cps.obs.exp.acagat+cps.obs.exp.acagca+cps.obs.exp.acagcc+cps.obs.exp.acagcg+cps.obs.exp.acagct+cps.obs.exp.acagga+cps.obs.exp.acaggc+cps.obs.exp.acaggg+cps.obs.exp.acaggt+cps.obs.exp.acagta+cps.obs.exp.acagtc+cps.obs.exp.acagtg+cps.obs.exp.acagtt+cps.obs.exp.acatac+cps.obs.exp.acatat+cps.obs.exp.acatca+cps.obs.exp.acatcc+cps.obs.exp.acatcg+cps.obs.exp.acatct+cps.obs.exp.acatgc+cps.obs.exp.acatgg+cps.obs.exp.acatgt+cps.obs.exp.acatta+cps.obs.exp.acattc+cps.obs.exp.acattg+cps.obs.exp.acattt+cps.obs.exp.accaaa+cps.obs.exp.accaac+cps.obs.exp.accaag+cps.obs.exp.accaat+cps.obs.exp.accaca+cps.obs.exp.accacc+cps.obs.exp.accacg+cps.obs.exp.accact+cps.obs.exp.accaga+cps.obs.exp.accagc+cps.obs.exp.accagg+cps.obs.exp.accagt+cps.obs.exp.accata+cps.obs.exp.accatc+cps.obs.exp.accatg+cps.obs.exp.accatt+cps.obs.exp.acccaa+cps.obs.exp.acccac+cps.obs.exp.acccag+cps.obs.exp.acccat+cps.obs.exp.acccca+cps.obs.exp.accccc+cps.obs.exp.accccg+cps.obs.exp.acccct+cps.obs.exp.acccga+cps.obs.exp.acccgc+cps.obs.exp.acccgg+cps.obs.exp.acccgt+cps.obs.exp.acccta+cps.obs.exp.accctc+cps.obs.exp.accctg+cps.obs.exp.accctt+cps.obs.exp.accgaa+cps.obs.exp.accgac+cps.obs.exp.accgag+cps.obs.exp.accgat+cps.obs.exp.accgca+cps.obs.exp.accgcc+cps.obs.exp.accgcg+cps.obs.exp.accgct+cps.obs.exp.accgga+cps.obs.exp.accggc+cps.obs.exp.accggg+cps.obs.exp.accggt+cps.obs.exp.accgta+cps.obs.exp.accgtc+cps.obs.exp.accgtg+cps.obs.exp.accgtt+cps.obs.exp.acctac+cps.obs.exp.acctat+cps.obs.exp.acctca+cps.obs.exp.acctcc+cps.obs.exp.acctcg+cps.obs.exp.acctct+cps.obs.exp.acctgc+cps.obs.exp.acctgg+cps.obs.exp.acctgt+cps.obs.exp.acctta+cps.obs.exp.accttc+cps.obs.exp.accttg+cps.obs.exp.accttt+cps.obs.exp.acgaaa+cps.obs.exp.acgaac+cps.obs.exp.acgaag+cps.obs.exp.acgaat+cps.obs.exp.acgaca+cps.obs.exp.acgacc+cps.obs.exp.acgacg+cps.obs.exp.acgact+cps.obs.exp.acgaga+cps.obs.exp.acgagc+cps.obs.exp.acgagg+cps.obs.exp.acgagt+cps.obs.exp.acgata+cps.obs.exp.acgatc+cps.obs.exp.acgatg+cps.obs.exp.acgatt+cps.obs.exp.acgcaa+cps.obs.exp.acgcac+cps.obs.exp.acgcag+cps.obs.exp.acgcat+cps.obs.exp.acgcca+cps.obs.exp.acgccc+cps.obs.exp.acgccg+cps.obs.exp.acgcct+cps.obs.exp.acgcga+cps.obs.exp.acgcgc+cps.obs.exp.acgcgg+cps.obs.exp.acgcgt+cps.obs.exp.acgcta+cps.obs.exp.acgctc+cps.obs.exp.acgctg+cps.obs.exp.acgctt+cps.obs.exp.acggaa+cps.obs.exp.acggac+cps.obs.exp.acggag+cps.obs.exp.acggat+cps.obs.exp.acggca+cps.obs.exp.acggcc+cps.obs.exp.acggcg+cps.obs.exp.acggct+cps.obs.exp.acggga+cps.obs.exp.acgggc+cps.obs.exp.acgggg+cps.obs.exp.acgggt+cps.obs.exp.acggta+cps.obs.exp.acggtc+cps.obs.exp.acggtg+cps.obs.exp.acggtt+cps.obs.exp.acgtac+cps.obs.exp.acgtat+cps.obs.exp.acgtca+cps.obs.exp.acgtcc+cps.obs.exp.acgtcg+cps.obs.exp.acgtct+cps.obs.exp.acgtgc+cps.obs.exp.acgtgg+cps.obs.exp.acgtgt+cps.obs.exp.acgtta+cps.obs.exp.acgttc+cps.obs.exp.acgttg+cps.obs.exp.acgttt+cps.obs.exp.actaaa+cps.obs.exp.actaac+cps.obs.exp.actaag+cps.obs.exp.actaat+cps.obs.exp.actaca+cps.obs.exp.actacc+cps.obs.exp.actacg+cps.obs.exp.actact+cps.obs.exp.actaga+cps.obs.exp.actagc+cps.obs.exp.actagg+cps.obs.exp.actagt+cps.obs.exp.actata+cps.obs.exp.actatc+cps.obs.exp.actatg+cps.obs.exp.actatt+cps.obs.exp.actcaa+cps.obs.exp.actcac+cps.obs.exp.actcag+cps.obs.exp.actcat+cps.obs.exp.actcca+cps.obs.exp.actccc+cps.obs.exp.actccg+cps.obs.exp.actcct+cps.obs.exp.actcga+cps.obs.exp.actcgc+cps.obs.exp.actcgg+cps.obs.exp.actcgt+cps.obs.exp.actcta+cps.obs.exp.actctc+cps.obs.exp.actctg+cps.obs.exp.actctt+cps.obs.exp.actgaa+cps.obs.exp.actgac+cps.obs.exp.actgag+cps.obs.exp.actgat+cps.obs.exp.actgca+cps.obs.exp.actgcc+cps.obs.exp.actgcg+cps.obs.exp.actgct+cps.obs.exp.actgga+cps.obs.exp.actggc+cps.obs.exp.actggg+cps.obs.exp.actggt+cps.obs.exp.actgta+cps.obs.exp.actgtc+cps.obs.exp.actgtg+cps.obs.exp.actgtt+cps.obs.exp.acttac+cps.obs.exp.acttat+cps.obs.exp.acttca+cps.obs.exp.acttcc+cps.obs.exp.acttcg+cps.obs.exp.acttct+cps.obs.exp.acttgc+cps.obs.exp.acttgg+cps.obs.exp.acttgt+cps.obs.exp.acttta+cps.obs.exp.actttc+cps.obs.exp.actttg+cps.obs.exp.actttt+cps.obs.exp.agaaaa+cps.obs.exp.agaaac+cps.obs.exp.agaaag+cps.obs.exp.agaaat+cps.obs.exp.agaaca+cps.obs.exp.agaacc+cps.obs.exp.agaacg+cps.obs.exp.agaact+cps.obs.exp.agaaga+cps.obs.exp.agaagc+cps.obs.exp.agaagg+cps.obs.exp.agaagt+cps.obs.exp.agaata+cps.obs.exp.agaatc+cps.obs.exp.agaatg+cps.obs.exp.agaatt+cps.obs.exp.agacaa+cps.obs.exp.agacac+cps.obs.exp.agacag+cps.obs.exp.agacat+cps.obs.exp.agacca+cps.obs.exp.agaccc+cps.obs.exp.agaccg+cps.obs.exp.agacct+cps.obs.exp.agacga+cps.obs.exp.agacgc+cps.obs.exp.agacgg+cps.obs.exp.agacgt+cps.obs.exp.agacta+cps.obs.exp.agactc+cps.obs.exp.agactg+cps.obs.exp.agactt+cps.obs.exp.agagaa+cps.obs.exp.agagac+cps.obs.exp.agagag+cps.obs.exp.agagat+cps.obs.exp.agagca+cps.obs.exp.agagcc+cps.obs.exp.agagcg+cps.obs.exp.agagct+cps.obs.exp.agagga+cps.obs.exp.agaggc+cps.obs.exp.agaggg+cps.obs.exp.agaggt+cps.obs.exp.agagta+cps.obs.exp.agagtc+cps.obs.exp.agagtg+cps.obs.exp.agagtt+cps.obs.exp.agatac+cps.obs.exp.agatat+cps.obs.exp.agatca+cps.obs.exp.agatcc+cps.obs.exp.agatcg+cps.obs.exp.agatct+cps.obs.exp.agatgc+cps.obs.exp.agatgg+cps.obs.exp.agatgt+cps.obs.exp.agatta+cps.obs.exp.agattc+cps.obs.exp.agattg+cps.obs.exp.agattt+cps.obs.exp.agcaaa+cps.obs.exp.agcaac+cps.obs.exp.agcaag+cps.obs.exp.agcaat+cps.obs.exp.agcaca+cps.obs.exp.agcacc+cps.obs.exp.agcacg+cps.obs.exp.agcact+cps.obs.exp.agcaga+cps.obs.exp.agcagc+cps.obs.exp.agcagg+cps.obs.exp.agcagt+cps.obs.exp.agcata+cps.obs.exp.agcatc+cps.obs.exp.agcatg+cps.obs.exp.agcatt+cps.obs.exp.agccaa+cps.obs.exp.agccac+cps.obs.exp.agccag+cps.obs.exp.agccat+cps.obs.exp.agccca+cps.obs.exp.agcccc+cps.obs.exp.agcccg+cps.obs.exp.agccct+cps.obs.exp.agccga+cps.obs.exp.agccgc+cps.obs.exp.agccgg+cps.obs.exp.agccgt+cps.obs.exp.agccta+cps.obs.exp.agcctc+cps.obs.exp.agcctg+cps.obs.exp.agcctt+cps.obs.exp.agcgaa+cps.obs.exp.agcgac+cps.obs.exp.agcgag+cps.obs.exp.agcgat+cps.obs.exp.agcgca+cps.obs.exp.agcgcc+cps.obs.exp.agcgcg+cps.obs.exp.agcgct+cps.obs.exp.agcgga+cps.obs.exp.agcggc+cps.obs.exp.agcggg+cps.obs.exp.agcggt+cps.obs.exp.agcgta+cps.obs.exp.agcgtc+cps.obs.exp.agcgtg+cps.obs.exp.agcgtt+cps.obs.exp.agctac+cps.obs.exp.agctat+cps.obs.exp.agctca+cps.obs.exp.agctcc+cps.obs.exp.agctcg+cps.obs.exp.agctct+cps.obs.exp.agctgc+cps.obs.exp.agctgg+cps.obs.exp.agctgt+cps.obs.exp.agctta+cps.obs.exp.agcttc+cps.obs.exp.agcttg+cps.obs.exp.agcttt+cps.obs.exp.aggaaa+cps.obs.exp.aggaac+cps.obs.exp.aggaag+cps.obs.exp.aggaat+cps.obs.exp.aggaca+cps.obs.exp.aggacc+cps.obs.exp.aggacg+cps.obs.exp.aggact+cps.obs.exp.aggaga+cps.obs.exp.aggagc+cps.obs.exp.aggagg+cps.obs.exp.aggagt+cps.obs.exp.aggata+cps.obs.exp.aggatc+cps.obs.exp.aggatg+cps.obs.exp.aggatt+cps.obs.exp.aggcaa+cps.obs.exp.aggcac+cps.obs.exp.aggcag+cps.obs.exp.aggcat+cps.obs.exp.aggcca+cps.obs.exp.aggccc+cps.obs.exp.aggccg+cps.obs.exp.aggcct+cps.obs.exp.aggcga+cps.obs.exp.aggcgc+cps.obs.exp.aggcgg+cps.obs.exp.aggcgt+cps.obs.exp.aggcta+cps.obs.exp.aggctc+cps.obs.exp.aggctg+cps.obs.exp.aggctt+cps.obs.exp.agggaa+cps.obs.exp.agggac+cps.obs.exp.agggag+cps.obs.exp.agggat+cps.obs.exp.agggca+cps.obs.exp.agggcc+cps.obs.exp.agggcg+cps.obs.exp.agggct+cps.obs.exp.agggga+cps.obs.exp.aggggc+cps.obs.exp.aggggg+cps.obs.exp.aggggt+cps.obs.exp.agggta+cps.obs.exp.agggtc+cps.obs.exp.agggtg+cps.obs.exp.agggtt+cps.obs.exp.aggtac+cps.obs.exp.aggtat+cps.obs.exp.aggtca+cps.obs.exp.aggtcc+cps.obs.exp.aggtcg+cps.obs.exp.aggtct+cps.obs.exp.aggtgc+cps.obs.exp.aggtgg+cps.obs.exp.aggtgt+cps.obs.exp.aggtta+cps.obs.exp.aggttc+cps.obs.exp.aggttg+cps.obs.exp.aggttt+cps.obs.exp.agtaaa+cps.obs.exp.agtaac+cps.obs.exp.agtaag+cps.obs.exp.agtaat+cps.obs.exp.agtaca+cps.obs.exp.agtacc+cps.obs.exp.agtacg+cps.obs.exp.agtact+cps.obs.exp.agtaga+cps.obs.exp.agtagc+cps.obs.exp.agtagg+cps.obs.exp.agtagt+cps.obs.exp.agtata+cps.obs.exp.agtatc+cps.obs.exp.agtatg+cps.obs.exp.agtatt+cps.obs.exp.agtcaa+cps.obs.exp.agtcac+cps.obs.exp.agtcag+cps.obs.exp.agtcat+cps.obs.exp.agtcca+cps.obs.exp.agtccc+cps.obs.exp.agtccg+cps.obs.exp.agtcct+cps.obs.exp.agtcga+cps.obs.exp.agtcgc+cps.obs.exp.agtcgg+cps.obs.exp.agtcgt+cps.obs.exp.agtcta+cps.obs.exp.agtctc+cps.obs.exp.agtctg+cps.obs.exp.agtctt+cps.obs.exp.agtgaa+cps.obs.exp.agtgac+cps.obs.exp.agtgag+cps.obs.exp.agtgat+cps.obs.exp.agtgca+cps.obs.exp.agtgcc+cps.obs.exp.agtgcg+cps.obs.exp.agtgct+cps.obs.exp.agtgga+cps.obs.exp.agtggc+cps.obs.exp.agtggg+cps.obs.exp.agtggt+cps.obs.exp.agtgta+cps.obs.exp.agtgtc+cps.obs.exp.agtgtg+cps.obs.exp.agtgtt+cps.obs.exp.agttac+cps.obs.exp.agttat+cps.obs.exp.agttca+cps.obs.exp.agttcc+cps.obs.exp.agttcg+cps.obs.exp.agttct+cps.obs.exp.agttgc+cps.obs.exp.agttgg+cps.obs.exp.agttgt+cps.obs.exp.agttta+cps.obs.exp.agtttc+cps.obs.exp.agtttg+cps.obs.exp.agtttt+cps.obs.exp.ataaaa+cps.obs.exp.ataaac+cps.obs.exp.ataaag+cps.obs.exp.ataaat+cps.obs.exp.ataaca+cps.obs.exp.ataacc+cps.obs.exp.ataacg+cps.obs.exp.ataact+cps.obs.exp.ataaga+cps.obs.exp.ataagc+cps.obs.exp.ataagg+cps.obs.exp.ataagt+cps.obs.exp.ataata+cps.obs.exp.ataatc+cps.obs.exp.ataatg+cps.obs.exp.ataatt+cps.obs.exp.atacaa+cps.obs.exp.atacac+cps.obs.exp.atacag+cps.obs.exp.atacat+cps.obs.exp.atacca+cps.obs.exp.ataccc+cps.obs.exp.ataccg+cps.obs.exp.atacct+cps.obs.exp.atacga+cps.obs.exp.atacgc+cps.obs.exp.atacgg+cps.obs.exp.atacgt+cps.obs.exp.atacta+cps.obs.exp.atactc+cps.obs.exp.atactg+cps.obs.exp.atactt+cps.obs.exp.atagaa+cps.obs.exp.atagac+cps.obs.exp.atagag+cps.obs.exp.atagat+cps.obs.exp.atagca+cps.obs.exp.atagcc+cps.obs.exp.atagcg+cps.obs.exp.atagct+cps.obs.exp.atagga+cps.obs.exp.ataggc+cps.obs.exp.ataggg+cps.obs.exp.ataggt+cps.obs.exp.atagta+cps.obs.exp.atagtc+cps.obs.exp.atagtg+cps.obs.exp.atagtt+cps.obs.exp.atatac+cps.obs.exp.atatat+cps.obs.exp.atatca+cps.obs.exp.atatcc+cps.obs.exp.atatcg+cps.obs.exp.atatct+cps.obs.exp.atatgc+cps.obs.exp.atatgg+cps.obs.exp.atatgt+cps.obs.exp.atatta+cps.obs.exp.atattc+cps.obs.exp.atattg+cps.obs.exp.atattt+cps.obs.exp.atcaaa+cps.obs.exp.atcaac+cps.obs.exp.atcaag+cps.obs.exp.atcaat+cps.obs.exp.atcaca+cps.obs.exp.atcacc+cps.obs.exp.atcacg+cps.obs.exp.atcact+cps.obs.exp.atcaga+cps.obs.exp.atcagc+cps.obs.exp.atcagg+cps.obs.exp.atcagt+cps.obs.exp.atcata+cps.obs.exp.atcatc+cps.obs.exp.atcatg+cps.obs.exp.atcatt+cps.obs.exp.atccaa+cps.obs.exp.atccac+cps.obs.exp.atccag+cps.obs.exp.atccat+cps.obs.exp.atccca+cps.obs.exp.atcccc+cps.obs.exp.atcccg+cps.obs.exp.atccct+cps.obs.exp.atccga+cps.obs.exp.atccgc+cps.obs.exp.atccgg+cps.obs.exp.atccgt+cps.obs.exp.atccta+cps.obs.exp.atcctc+cps.obs.exp.atcctg+cps.obs.exp.atcctt+cps.obs.exp.atcgaa+cps.obs.exp.atcgac+cps.obs.exp.atcgag+cps.obs.exp.atcgat+cps.obs.exp.atcgca+cps.obs.exp.atcgcc+cps.obs.exp.atcgcg+cps.obs.exp.atcgct+cps.obs.exp.atcgga+cps.obs.exp.atcggc+cps.obs.exp.atcggg+cps.obs.exp.atcggt+cps.obs.exp.atcgta+cps.obs.exp.atcgtc+cps.obs.exp.atcgtg+cps.obs.exp.atcgtt+cps.obs.exp.atctac+cps.obs.exp.atctat+cps.obs.exp.atctca+cps.obs.exp.atctcc+cps.obs.exp.atctcg+cps.obs.exp.atctct+cps.obs.exp.atctgc+cps.obs.exp.atctgg+cps.obs.exp.atctgt+cps.obs.exp.atctta+cps.obs.exp.atcttc+cps.obs.exp.atcttg+cps.obs.exp.atcttt+cps.obs.exp.atgaaa+cps.obs.exp.atgaac+cps.obs.exp.atgaag+cps.obs.exp.atgaat+cps.obs.exp.atgaca+cps.obs.exp.atgacc+cps.obs.exp.atgacg+cps.obs.exp.atgact+cps.obs.exp.atgaga+cps.obs.exp.atgagc+cps.obs.exp.atgagg+cps.obs.exp.atgagt+cps.obs.exp.atgata+cps.obs.exp.atgatc+cps.obs.exp.atgatg+cps.obs.exp.atgatt+cps.obs.exp.atgcaa+cps.obs.exp.atgcac+cps.obs.exp.atgcag+cps.obs.exp.atgcat+cps.obs.exp.atgcca+cps.obs.exp.atgccc+cps.obs.exp.atgccg+cps.obs.exp.atgcct+cps.obs.exp.atgcga+cps.obs.exp.atgcgc+cps.obs.exp.atgcgg+cps.obs.exp.atgcgt+cps.obs.exp.atgcta+cps.obs.exp.atgctc+cps.obs.exp.atgctg+cps.obs.exp.atgctt+cps.obs.exp.atggaa+cps.obs.exp.atggac+cps.obs.exp.atggag+cps.obs.exp.atggat+cps.obs.exp.atggca+cps.obs.exp.atggcc+cps.obs.exp.atggcg+cps.obs.exp.atggct+cps.obs.exp.atggga+cps.obs.exp.atgggc+cps.obs.exp.atgggg+cps.obs.exp.atgggt+cps.obs.exp.atggta+cps.obs.exp.atggtc+cps.obs.exp.atggtg+cps.obs.exp.atggtt+cps.obs.exp.atgtac+cps.obs.exp.atgtat+cps.obs.exp.atgtca+cps.obs.exp.atgtcc+cps.obs.exp.atgtcg+cps.obs.exp.atgtct+cps.obs.exp.atgtgc+cps.obs.exp.atgtgg+cps.obs.exp.atgtgt+cps.obs.exp.atgtta+cps.obs.exp.atgttc+cps.obs.exp.atgttg+cps.obs.exp.atgttt+cps.obs.exp.attaaa+cps.obs.exp.attaac+cps.obs.exp.attaag+cps.obs.exp.attaat+cps.obs.exp.attaca+cps.obs.exp.attacc+cps.obs.exp.attacg+cps.obs.exp.attact+cps.obs.exp.attaga+cps.obs.exp.attagc+cps.obs.exp.attagg+cps.obs.exp.attagt+cps.obs.exp.attata+cps.obs.exp.attatc+cps.obs.exp.attatg+cps.obs.exp.attatt+cps.obs.exp.attcaa+cps.obs.exp.attcac+cps.obs.exp.attcag+cps.obs.exp.attcat+cps.obs.exp.attcca+cps.obs.exp.attccc+cps.obs.exp.attccg+cps.obs.exp.attcct+cps.obs.exp.attcga+cps.obs.exp.attcgc+cps.obs.exp.attcgg+cps.obs.exp.attcgt+cps.obs.exp.attcta+cps.obs.exp.attctc+cps.obs.exp.attctg+cps.obs.exp.attctt+cps.obs.exp.attgaa+cps.obs.exp.attgac+cps.obs.exp.attgag+cps.obs.exp.attgat+cps.obs.exp.attgca+cps.obs.exp.attgcc+cps.obs.exp.attgcg+cps.obs.exp.attgct+cps.obs.exp.attgga+cps.obs.exp.attggc+cps.obs.exp.attggg+cps.obs.exp.attggt+cps.obs.exp.attgta+cps.obs.exp.attgtc+cps.obs.exp.attgtg+cps.obs.exp.attgtt+cps.obs.exp.atttac+cps.obs.exp.atttat+cps.obs.exp.atttca+cps.obs.exp.atttcc+cps.obs.exp.atttcg+cps.obs.exp.atttct+cps.obs.exp.atttgc+cps.obs.exp.atttgg+cps.obs.exp.atttgt+cps.obs.exp.atttta+cps.obs.exp.attttc+cps.obs.exp.attttg+cps.obs.exp.attttt+cps.obs.exp.caaaaa+cps.obs.exp.caaaac+cps.obs.exp.caaaag+cps.obs.exp.caaaat+cps.obs.exp.caaaca+cps.obs.exp.caaacc+cps.obs.exp.caaacg+cps.obs.exp.caaact+cps.obs.exp.caaaga+cps.obs.exp.caaagc+cps.obs.exp.caaagg+cps.obs.exp.caaagt+cps.obs.exp.caaata+cps.obs.exp.caaatc+cps.obs.exp.caaatg+cps.obs.exp.caaatt+cps.obs.exp.caacaa+cps.obs.exp.caacac+cps.obs.exp.caacag+cps.obs.exp.caacat+cps.obs.exp.caacca+cps.obs.exp.caaccc+cps.obs.exp.caaccg+cps.obs.exp.caacct+cps.obs.exp.caacga+cps.obs.exp.caacgc+cps.obs.exp.caacgg+cps.obs.exp.caacgt+cps.obs.exp.caacta+cps.obs.exp.caactc+cps.obs.exp.caactg+cps.obs.exp.caactt+cps.obs.exp.caagaa+cps.obs.exp.caagac+cps.obs.exp.caagag+cps.obs.exp.caagat+cps.obs.exp.caagca+cps.obs.exp.caagcc+cps.obs.exp.caagcg+cps.obs.exp.caagct+cps.obs.exp.caagga+cps.obs.exp.caaggc+cps.obs.exp.caaggg+cps.obs.exp.caaggt+cps.obs.exp.caagta+cps.obs.exp.caagtc+cps.obs.exp.caagtg+cps.obs.exp.caagtt+cps.obs.exp.caatac+cps.obs.exp.caatat+cps.obs.exp.caatca+cps.obs.exp.caatcc+cps.obs.exp.caatcg+cps.obs.exp.caatct+cps.obs.exp.caatgc+cps.obs.exp.caatgg+cps.obs.exp.caatgt+cps.obs.exp.caatta+cps.obs.exp.caattc+cps.obs.exp.caattg+cps.obs.exp.caattt+cps.obs.exp.cacaaa+cps.obs.exp.cacaac+cps.obs.exp.cacaag+cps.obs.exp.cacaat+cps.obs.exp.cacaca+cps.obs.exp.cacacc+cps.obs.exp.cacacg+cps.obs.exp.cacact+cps.obs.exp.cacaga+cps.obs.exp.cacagc+cps.obs.exp.cacagg+cps.obs.exp.cacagt+cps.obs.exp.cacata+cps.obs.exp.cacatc+cps.obs.exp.cacatg+cps.obs.exp.cacatt+cps.obs.exp.caccaa+cps.obs.exp.caccac+cps.obs.exp.caccag+cps.obs.exp.caccat+cps.obs.exp.caccca+cps.obs.exp.cacccc+cps.obs.exp.cacccg+cps.obs.exp.caccct+cps.obs.exp.caccga+cps.obs.exp.caccgc+cps.obs.exp.caccgg+cps.obs.exp.caccgt+cps.obs.exp.caccta+cps.obs.exp.cacctc+cps.obs.exp.cacctg+cps.obs.exp.cacctt+cps.obs.exp.cacgaa+cps.obs.exp.cacgac+cps.obs.exp.cacgag+cps.obs.exp.cacgat+cps.obs.exp.cacgca+cps.obs.exp.cacgcc+cps.obs.exp.cacgcg+cps.obs.exp.cacgct+cps.obs.exp.cacgga+cps.obs.exp.cacggc+cps.obs.exp.cacggg+cps.obs.exp.cacggt+cps.obs.exp.cacgta+cps.obs.exp.cacgtc+cps.obs.exp.cacgtg+cps.obs.exp.cacgtt+cps.obs.exp.cactac+cps.obs.exp.cactat+cps.obs.exp.cactca+cps.obs.exp.cactcc+cps.obs.exp.cactcg+cps.obs.exp.cactct+cps.obs.exp.cactgc+cps.obs.exp.cactgg+cps.obs.exp.cactgt+cps.obs.exp.cactta+cps.obs.exp.cacttc+cps.obs.exp.cacttg+cps.obs.exp.cacttt+cps.obs.exp.cagaaa+cps.obs.exp.cagaac+cps.obs.exp.cagaag+cps.obs.exp.cagaat+cps.obs.exp.cagaca+cps.obs.exp.cagacc+cps.obs.exp.cagacg+cps.obs.exp.cagact+cps.obs.exp.cagaga+cps.obs.exp.cagagc+cps.obs.exp.cagagg+cps.obs.exp.cagagt+cps.obs.exp.cagata+cps.obs.exp.cagatc+cps.obs.exp.cagatg+cps.obs.exp.cagatt+cps.obs.exp.cagcaa+cps.obs.exp.cagcac+cps.obs.exp.cagcag+cps.obs.exp.cagcat+cps.obs.exp.cagcca+cps.obs.exp.cagccc+cps.obs.exp.cagccg+cps.obs.exp.cagcct+cps.obs.exp.cagcga+cps.obs.exp.cagcgc+cps.obs.exp.cagcgg+cps.obs.exp.cagcgt+cps.obs.exp.cagcta+cps.obs.exp.cagctc+cps.obs.exp.cagctg+cps.obs.exp.cagctt+cps.obs.exp.caggaa+cps.obs.exp.caggac+cps.obs.exp.caggag+cps.obs.exp.caggat+cps.obs.exp.caggca+cps.obs.exp.caggcc+cps.obs.exp.caggcg+cps.obs.exp.caggct+cps.obs.exp.caggga+cps.obs.exp.cagggc+cps.obs.exp.cagggg+cps.obs.exp.cagggt+cps.obs.exp.caggta+cps.obs.exp.caggtc+cps.obs.exp.caggtg+cps.obs.exp.caggtt+cps.obs.exp.cagtac+cps.obs.exp.cagtat+cps.obs.exp.cagtca+cps.obs.exp.cagtcc+cps.obs.exp.cagtcg+cps.obs.exp.cagtct+cps.obs.exp.cagtgc+cps.obs.exp.cagtgg+cps.obs.exp.cagtgt+cps.obs.exp.cagtta+cps.obs.exp.cagttc+cps.obs.exp.cagttg+cps.obs.exp.cagttt+cps.obs.exp.cataaa+cps.obs.exp.cataac+cps.obs.exp.cataag+cps.obs.exp.cataat+cps.obs.exp.cataca+cps.obs.exp.catacc+cps.obs.exp.catacg+cps.obs.exp.catact+cps.obs.exp.cataga+cps.obs.exp.catagc+cps.obs.exp.catagg+cps.obs.exp.catagt+cps.obs.exp.catata+cps.obs.exp.catatc+cps.obs.exp.catatg+cps.obs.exp.catatt+cps.obs.exp.catcaa+cps.obs.exp.catcac+cps.obs.exp.catcag+cps.obs.exp.catcat+cps.obs.exp.catcca+cps.obs.exp.catccc+cps.obs.exp.catccg+cps.obs.exp.catcct+cps.obs.exp.catcga+cps.obs.exp.catcgc+cps.obs.exp.catcgg+cps.obs.exp.catcgt+cps.obs.exp.catcta+cps.obs.exp.catctc+cps.obs.exp.catctg+cps.obs.exp.catctt+cps.obs.exp.catgaa+cps.obs.exp.catgac+cps.obs.exp.catgag+cps.obs.exp.catgat+cps.obs.exp.catgca+cps.obs.exp.catgcc+cps.obs.exp.catgcg+cps.obs.exp.catgct+cps.obs.exp.catgga+cps.obs.exp.catggc+cps.obs.exp.catggg+cps.obs.exp.catggt+cps.obs.exp.catgta+cps.obs.exp.catgtc+cps.obs.exp.catgtg+cps.obs.exp.catgtt+cps.obs.exp.cattac+cps.obs.exp.cattat+cps.obs.exp.cattca+cps.obs.exp.cattcc+cps.obs.exp.cattcg+cps.obs.exp.cattct+cps.obs.exp.cattgc+cps.obs.exp.cattgg+cps.obs.exp.cattgt+cps.obs.exp.cattta+cps.obs.exp.catttc+cps.obs.exp.catttg+cps.obs.exp.catttt+cps.obs.exp.ccaaaa+cps.obs.exp.ccaaac+cps.obs.exp.ccaaag+cps.obs.exp.ccaaat+cps.obs.exp.ccaaca+cps.obs.exp.ccaacc+cps.obs.exp.ccaacg+cps.obs.exp.ccaact+cps.obs.exp.ccaaga+cps.obs.exp.ccaagc+cps.obs.exp.ccaagg+cps.obs.exp.ccaagt+cps.obs.exp.ccaata+cps.obs.exp.ccaatc+cps.obs.exp.ccaatg+cps.obs.exp.ccaatt+cps.obs.exp.ccacaa+cps.obs.exp.ccacac+cps.obs.exp.ccacag+cps.obs.exp.ccacat+cps.obs.exp.ccacca+cps.obs.exp.ccaccc+cps.obs.exp.ccaccg+cps.obs.exp.ccacct+cps.obs.exp.ccacga+cps.obs.exp.ccacgc+cps.obs.exp.ccacgg+cps.obs.exp.ccacgt+cps.obs.exp.ccacta+cps.obs.exp.ccactc+cps.obs.exp.ccactg+cps.obs.exp.ccactt+cps.obs.exp.ccagaa+cps.obs.exp.ccagac+cps.obs.exp.ccagag+cps.obs.exp.ccagat+cps.obs.exp.ccagca+cps.obs.exp.ccagcc+cps.obs.exp.ccagcg+cps.obs.exp.ccagct+cps.obs.exp.ccagga+cps.obs.exp.ccaggc+cps.obs.exp.ccaggg+cps.obs.exp.ccaggt+cps.obs.exp.ccagta+cps.obs.exp.ccagtc+cps.obs.exp.ccagtg+cps.obs.exp.ccagtt+cps.obs.exp.ccatac+cps.obs.exp.ccatat+cps.obs.exp.ccatca+cps.obs.exp.ccatcc+cps.obs.exp.ccatcg+cps.obs.exp.ccatct+cps.obs.exp.ccatgc+cps.obs.exp.ccatgg+cps.obs.exp.ccatgt+cps.obs.exp.ccatta+cps.obs.exp.ccattc+cps.obs.exp.ccattg+cps.obs.exp.ccattt+cps.obs.exp.cccaaa+cps.obs.exp.cccaac+cps.obs.exp.cccaag+cps.obs.exp.cccaat+cps.obs.exp.cccaca+cps.obs.exp.cccacc+cps.obs.exp.cccacg+cps.obs.exp.cccact+cps.obs.exp.cccaga+cps.obs.exp.cccagc+cps.obs.exp.cccagg+cps.obs.exp.cccagt+cps.obs.exp.cccata+cps.obs.exp.cccatc+cps.obs.exp.cccatg+cps.obs.exp.cccatt+cps.obs.exp.ccccaa+cps.obs.exp.ccccac+cps.obs.exp.ccccag+cps.obs.exp.ccccat+cps.obs.exp.ccccca+cps.obs.exp.cccccc+cps.obs.exp.cccccg+cps.obs.exp.ccccct+cps.obs.exp.ccccga+cps.obs.exp.ccccgc+cps.obs.exp.ccccgg+cps.obs.exp.ccccgt+cps.obs.exp.ccccta+cps.obs.exp.cccctc+cps.obs.exp.cccctg+cps.obs.exp.cccctt+cps.obs.exp.cccgaa+cps.obs.exp.cccgac+cps.obs.exp.cccgag+cps.obs.exp.cccgat+cps.obs.exp.cccgca+cps.obs.exp.cccgcc+cps.obs.exp.cccgcg+cps.obs.exp.cccgct+cps.obs.exp.cccgga+cps.obs.exp.cccggc+cps.obs.exp.cccggg+cps.obs.exp.cccggt+cps.obs.exp.cccgta+cps.obs.exp.cccgtc+cps.obs.exp.cccgtg+cps.obs.exp.cccgtt+cps.obs.exp.ccctac+cps.obs.exp.ccctat+cps.obs.exp.ccctca+cps.obs.exp.ccctcc+cps.obs.exp.ccctcg+cps.obs.exp.ccctct+cps.obs.exp.ccctgc+cps.obs.exp.ccctgg+cps.obs.exp.ccctgt+cps.obs.exp.ccctta+cps.obs.exp.cccttc+cps.obs.exp.cccttg+cps.obs.exp.cccttt+cps.obs.exp.ccgaaa+cps.obs.exp.ccgaac+cps.obs.exp.ccgaag+cps.obs.exp.ccgaat+cps.obs.exp.ccgaca+cps.obs.exp.ccgacc+cps.obs.exp.ccgacg+cps.obs.exp.ccgact+cps.obs.exp.ccgaga+cps.obs.exp.ccgagc+cps.obs.exp.ccgagg+cps.obs.exp.ccgagt+cps.obs.exp.ccgata+cps.obs.exp.ccgatc+cps.obs.exp.ccgatg+cps.obs.exp.ccgatt+cps.obs.exp.ccgcaa+cps.obs.exp.ccgcac+cps.obs.exp.ccgcag+cps.obs.exp.ccgcat+cps.obs.exp.ccgcca+cps.obs.exp.ccgccc+cps.obs.exp.ccgccg+cps.obs.exp.ccgcct+cps.obs.exp.ccgcga+cps.obs.exp.ccgcgc+cps.obs.exp.ccgcgg+cps.obs.exp.ccgcgt+cps.obs.exp.ccgcta+cps.obs.exp.ccgctc+cps.obs.exp.ccgctg+cps.obs.exp.ccgctt+cps.obs.exp.ccggaa+cps.obs.exp.ccggac+cps.obs.exp.ccggag+cps.obs.exp.ccggat+cps.obs.exp.ccggca+cps.obs.exp.ccggcc+cps.obs.exp.ccggcg+cps.obs.exp.ccggct+cps.obs.exp.ccggga+cps.obs.exp.ccgggc+cps.obs.exp.ccgggg+cps.obs.exp.ccgggt+cps.obs.exp.ccggta+cps.obs.exp.ccggtc+cps.obs.exp.ccggtg+cps.obs.exp.ccggtt+cps.obs.exp.ccgtac+cps.obs.exp.ccgtat+cps.obs.exp.ccgtca+cps.obs.exp.ccgtcc+cps.obs.exp.ccgtcg+cps.obs.exp.ccgtct+cps.obs.exp.ccgtgc+cps.obs.exp.ccgtgg+cps.obs.exp.ccgtgt+cps.obs.exp.ccgtta+cps.obs.exp.ccgttc+cps.obs.exp.ccgttg+cps.obs.exp.ccgttt+cps.obs.exp.cctaaa+cps.obs.exp.cctaac+cps.obs.exp.cctaag+cps.obs.exp.cctaat+cps.obs.exp.cctaca+cps.obs.exp.cctacc+cps.obs.exp.cctacg+cps.obs.exp.cctact+cps.obs.exp.cctaga+cps.obs.exp.cctagc+cps.obs.exp.cctagg+cps.obs.exp.cctagt+cps.obs.exp.cctata+cps.obs.exp.cctatc+cps.obs.exp.cctatg+cps.obs.exp.cctatt+cps.obs.exp.cctcaa+cps.obs.exp.cctcac+cps.obs.exp.cctcag+cps.obs.exp.cctcat+cps.obs.exp.cctcca+cps.obs.exp.cctccc+cps.obs.exp.cctccg+cps.obs.exp.cctcct+cps.obs.exp.cctcga+cps.obs.exp.cctcgc+cps.obs.exp.cctcgg+cps.obs.exp.cctcgt+cps.obs.exp.cctcta+cps.obs.exp.cctctc+cps.obs.exp.cctctg+cps.obs.exp.cctctt+cps.obs.exp.cctgaa+cps.obs.exp.cctgac+cps.obs.exp.cctgag+cps.obs.exp.cctgat+cps.obs.exp.cctgca+cps.obs.exp.cctgcc+cps.obs.exp.cctgcg+cps.obs.exp.cctgct+cps.obs.exp.cctgga+cps.obs.exp.cctggc+cps.obs.exp.cctggg+cps.obs.exp.cctggt+cps.obs.exp.cctgta+cps.obs.exp.cctgtc+cps.obs.exp.cctgtg+cps.obs.exp.cctgtt+cps.obs.exp.ccttac+cps.obs.exp.ccttat+cps.obs.exp.ccttca+cps.obs.exp.ccttcc+cps.obs.exp.ccttcg+cps.obs.exp.ccttct+cps.obs.exp.ccttgc+cps.obs.exp.ccttgg+cps.obs.exp.ccttgt+cps.obs.exp.ccttta+cps.obs.exp.cctttc+cps.obs.exp.cctttg+cps.obs.exp.cctttt+cps.obs.exp.cgaaaa+cps.obs.exp.cgaaac+cps.obs.exp.cgaaag+cps.obs.exp.cgaaat+cps.obs.exp.cgaaca+cps.obs.exp.cgaacc+cps.obs.exp.cgaacg+cps.obs.exp.cgaact+cps.obs.exp.cgaaga+cps.obs.exp.cgaagc+cps.obs.exp.cgaagg+cps.obs.exp.cgaagt+cps.obs.exp.cgaata+cps.obs.exp.cgaatc+cps.obs.exp.cgaatg+cps.obs.exp.cgaatt+cps.obs.exp.cgacaa+cps.obs.exp.cgacac+cps.obs.exp.cgacag+cps.obs.exp.cgacat+cps.obs.exp.cgacca+cps.obs.exp.cgaccc+cps.obs.exp.cgaccg+cps.obs.exp.cgacct+cps.obs.exp.cgacga+cps.obs.exp.cgacgc+cps.obs.exp.cgacgg+cps.obs.exp.cgacgt+cps.obs.exp.cgacta+cps.obs.exp.cgactc+cps.obs.exp.cgactg+cps.obs.exp.cgactt+cps.obs.exp.cgagaa+cps.obs.exp.cgagac+cps.obs.exp.cgagag+cps.obs.exp.cgagat+cps.obs.exp.cgagca+cps.obs.exp.cgagcc+cps.obs.exp.cgagcg+cps.obs.exp.cgagct+cps.obs.exp.cgagga+cps.obs.exp.cgaggc+cps.obs.exp.cgaggg+cps.obs.exp.cgaggt+cps.obs.exp.cgagta+cps.obs.exp.cgagtc+cps.obs.exp.cgagtg+cps.obs.exp.cgagtt+cps.obs.exp.cgatac+cps.obs.exp.cgatat+cps.obs.exp.cgatca+cps.obs.exp.cgatcc+cps.obs.exp.cgatcg+cps.obs.exp.cgatct+cps.obs.exp.cgatgc+cps.obs.exp.cgatgg+cps.obs.exp.cgatgt+cps.obs.exp.cgatta+cps.obs.exp.cgattc+cps.obs.exp.cgattg+cps.obs.exp.cgattt+cps.obs.exp.cgcaaa+cps.obs.exp.cgcaac+cps.obs.exp.cgcaag+cps.obs.exp.cgcaat+cps.obs.exp.cgcaca+cps.obs.exp.cgcacc+cps.obs.exp.cgcacg+cps.obs.exp.cgcact+cps.obs.exp.cgcaga+cps.obs.exp.cgcagc+cps.obs.exp.cgcagg+cps.obs.exp.cgcagt+cps.obs.exp.cgcata+cps.obs.exp.cgcatc+cps.obs.exp.cgcatg+cps.obs.exp.cgcatt+cps.obs.exp.cgccaa+cps.obs.exp.cgccac+cps.obs.exp.cgccag+cps.obs.exp.cgccat+cps.obs.exp.cgccca+cps.obs.exp.cgcccc+cps.obs.exp.cgcccg+cps.obs.exp.cgccct+cps.obs.exp.cgccga+cps.obs.exp.cgccgc+cps.obs.exp.cgccgg+cps.obs.exp.cgccgt+cps.obs.exp.cgccta+cps.obs.exp.cgcctc+cps.obs.exp.cgcctg+cps.obs.exp.cgcctt+cps.obs.exp.cgcgaa+cps.obs.exp.cgcgac+cps.obs.exp.cgcgag+cps.obs.exp.cgcgat+cps.obs.exp.cgcgca+cps.obs.exp.cgcgcc+cps.obs.exp.cgcgcg+cps.obs.exp.cgcgct+cps.obs.exp.cgcgga+cps.obs.exp.cgcggc+cps.obs.exp.cgcggg+cps.obs.exp.cgcggt+cps.obs.exp.cgcgta+cps.obs.exp.cgcgtc+cps.obs.exp.cgcgtg+cps.obs.exp.cgcgtt+cps.obs.exp.cgctac+cps.obs.exp.cgctat+cps.obs.exp.cgctca+cps.obs.exp.cgctcc+cps.obs.exp.cgctcg+cps.obs.exp.cgctct+cps.obs.exp.cgctgc+cps.obs.exp.cgctgg+cps.obs.exp.cgctgt+cps.obs.exp.cgctta+cps.obs.exp.cgcttc+cps.obs.exp.cgcttg+cps.obs.exp.cgcttt+cps.obs.exp.cggaaa+cps.obs.exp.cggaac+cps.obs.exp.cggaag+cps.obs.exp.cggaat+cps.obs.exp.cggaca+cps.obs.exp.cggacc+cps.obs.exp.cggacg+cps.obs.exp.cggact+cps.obs.exp.cggaga+cps.obs.exp.cggagc+cps.obs.exp.cggagg+cps.obs.exp.cggagt+cps.obs.exp.cggata+cps.obs.exp.cggatc+cps.obs.exp.cggatg+cps.obs.exp.cggatt+cps.obs.exp.cggcaa+cps.obs.exp.cggcac+cps.obs.exp.cggcag+cps.obs.exp.cggcat+cps.obs.exp.cggcca+cps.obs.exp.cggccc+cps.obs.exp.cggccg+cps.obs.exp.cggcct+cps.obs.exp.cggcga+cps.obs.exp.cggcgc+cps.obs.exp.cggcgg+cps.obs.exp.cggcgt+cps.obs.exp.cggcta+cps.obs.exp.cggctc+cps.obs.exp.cggctg+cps.obs.exp.cggctt+cps.obs.exp.cgggaa+cps.obs.exp.cgggac+cps.obs.exp.cgggag+cps.obs.exp.cgggat+cps.obs.exp.cgggca+cps.obs.exp.cgggcc+cps.obs.exp.cgggcg+cps.obs.exp.cgggct+cps.obs.exp.cgggga+cps.obs.exp.cggggc+cps.obs.exp.cggggg+cps.obs.exp.cggggt+cps.obs.exp.cgggta+cps.obs.exp.cgggtc+cps.obs.exp.cgggtg+cps.obs.exp.cgggtt+cps.obs.exp.cggtac+cps.obs.exp.cggtat+cps.obs.exp.cggtca+cps.obs.exp.cggtcc+cps.obs.exp.cggtcg+cps.obs.exp.cggtct+cps.obs.exp.cggtgc+cps.obs.exp.cggtgg+cps.obs.exp.cggtgt+cps.obs.exp.cggtta+cps.obs.exp.cggttc+cps.obs.exp.cggttg+cps.obs.exp.cggttt+cps.obs.exp.cgtaaa+cps.obs.exp.cgtaac+cps.obs.exp.cgtaag+cps.obs.exp.cgtaat+cps.obs.exp.cgtaca+cps.obs.exp.cgtacc+cps.obs.exp.cgtacg+cps.obs.exp.cgtact+cps.obs.exp.cgtaga+cps.obs.exp.cgtagc+cps.obs.exp.cgtagg+cps.obs.exp.cgtagt+cps.obs.exp.cgtata+cps.obs.exp.cgtatc+cps.obs.exp.cgtatg+cps.obs.exp.cgtatt+cps.obs.exp.cgtcaa+cps.obs.exp.cgtcac+cps.obs.exp.cgtcag+cps.obs.exp.cgtcat+cps.obs.exp.cgtcca+cps.obs.exp.cgtccc+cps.obs.exp.cgtccg+cps.obs.exp.cgtcct+cps.obs.exp.cgtcga+cps.obs.exp.cgtcgc+cps.obs.exp.cgtcgg+cps.obs.exp.cgtcgt+cps.obs.exp.cgtcta+cps.obs.exp.cgtctc+cps.obs.exp.cgtctg+cps.obs.exp.cgtctt+cps.obs.exp.cgtgaa+cps.obs.exp.cgtgac+cps.obs.exp.cgtgag+cps.obs.exp.cgtgat+cps.obs.exp.cgtgca+cps.obs.exp.cgtgcc+cps.obs.exp.cgtgcg+cps.obs.exp.cgtgct+cps.obs.exp.cgtgga+cps.obs.exp.cgtggc+cps.obs.exp.cgtggg+cps.obs.exp.cgtggt+cps.obs.exp.cgtgta+cps.obs.exp.cgtgtc+cps.obs.exp.cgtgtg+cps.obs.exp.cgtgtt+cps.obs.exp.cgttac+cps.obs.exp.cgttat+cps.obs.exp.cgttca+cps.obs.exp.cgttcc+cps.obs.exp.cgttcg+cps.obs.exp.cgttct+cps.obs.exp.cgttgc+cps.obs.exp.cgttgg+cps.obs.exp.cgttgt+cps.obs.exp.cgttta+cps.obs.exp.cgtttc+cps.obs.exp.cgtttg+cps.obs.exp.cgtttt+cps.obs.exp.ctaaaa+cps.obs.exp.ctaaac+cps.obs.exp.ctaaag+cps.obs.exp.ctaaat+cps.obs.exp.ctaaca+cps.obs.exp.ctaacc+cps.obs.exp.ctaacg+cps.obs.exp.ctaact+cps.obs.exp.ctaaga+cps.obs.exp.ctaagc+cps.obs.exp.ctaagg+cps.obs.exp.ctaagt+cps.obs.exp.ctaata+cps.obs.exp.ctaatc+cps.obs.exp.ctaatg+cps.obs.exp.ctaatt+cps.obs.exp.ctacaa+cps.obs.exp.ctacac+cps.obs.exp.ctacag+cps.obs.exp.ctacat+cps.obs.exp.ctacca+cps.obs.exp.ctaccc+cps.obs.exp.ctaccg+cps.obs.exp.ctacct+cps.obs.exp.ctacga+cps.obs.exp.ctacgc+cps.obs.exp.ctacgg+cps.obs.exp.ctacgt+cps.obs.exp.ctacta+cps.obs.exp.ctactc+cps.obs.exp.ctactg+cps.obs.exp.ctactt+cps.obs.exp.ctagaa+cps.obs.exp.ctagac+cps.obs.exp.ctagag+cps.obs.exp.ctagat+cps.obs.exp.ctagca+cps.obs.exp.ctagcc+cps.obs.exp.ctagcg+cps.obs.exp.ctagct+cps.obs.exp.ctagga+cps.obs.exp.ctaggc+cps.obs.exp.ctaggg+cps.obs.exp.ctaggt+cps.obs.exp.ctagta+cps.obs.exp.ctagtc+cps.obs.exp.ctagtg+cps.obs.exp.ctagtt+cps.obs.exp.ctatac+cps.obs.exp.ctatat+cps.obs.exp.ctatca+cps.obs.exp.ctatcc+cps.obs.exp.ctatcg+cps.obs.exp.ctatct+cps.obs.exp.ctatgc+cps.obs.exp.ctatgg+cps.obs.exp.ctatgt+cps.obs.exp.ctatta+cps.obs.exp.ctattc+cps.obs.exp.ctattg+cps.obs.exp.ctattt+cps.obs.exp.ctcaaa+cps.obs.exp.ctcaac+cps.obs.exp.ctcaag+cps.obs.exp.ctcaat+cps.obs.exp.ctcaca+cps.obs.exp.ctcacc+cps.obs.exp.ctcacg+cps.obs.exp.ctcact+cps.obs.exp.ctcaga+cps.obs.exp.ctcagc+cps.obs.exp.ctcagg+cps.obs.exp.ctcagt+cps.obs.exp.ctcata+cps.obs.exp.ctcatc+cps.obs.exp.ctcatg+cps.obs.exp.ctcatt+cps.obs.exp.ctccaa+cps.obs.exp.ctccac+cps.obs.exp.ctccag+cps.obs.exp.ctccat+cps.obs.exp.ctccca+cps.obs.exp.ctcccc+cps.obs.exp.ctcccg+cps.obs.exp.ctccct+cps.obs.exp.ctccga+cps.obs.exp.ctccgc+cps.obs.exp.ctccgg+cps.obs.exp.ctccgt+cps.obs.exp.ctccta+cps.obs.exp.ctcctc+cps.obs.exp.ctcctg+cps.obs.exp.ctcctt+cps.obs.exp.ctcgaa+cps.obs.exp.ctcgac+cps.obs.exp.ctcgag+cps.obs.exp.ctcgat+cps.obs.exp.ctcgca+cps.obs.exp.ctcgcc+cps.obs.exp.ctcgcg+cps.obs.exp.ctcgct+cps.obs.exp.ctcgga+cps.obs.exp.ctcggc+cps.obs.exp.ctcggg+cps.obs.exp.ctcggt+cps.obs.exp.ctcgta+cps.obs.exp.ctcgtc+cps.obs.exp.ctcgtg+cps.obs.exp.ctcgtt+cps.obs.exp.ctctac+cps.obs.exp.ctctat+cps.obs.exp.ctctca+cps.obs.exp.ctctcc+cps.obs.exp.ctctcg+cps.obs.exp.ctctct+cps.obs.exp.ctctgc+cps.obs.exp.ctctgg+cps.obs.exp.ctctgt+cps.obs.exp.ctctta+cps.obs.exp.ctcttc+cps.obs.exp.ctcttg+cps.obs.exp.ctcttt+cps.obs.exp.ctgaaa+cps.obs.exp.ctgaac+cps.obs.exp.ctgaag+cps.obs.exp.ctgaat+cps.obs.exp.ctgaca+cps.obs.exp.ctgacc+cps.obs.exp.ctgacg+cps.obs.exp.ctgact+cps.obs.exp.ctgaga+cps.obs.exp.ctgagc+cps.obs.exp.ctgagg+cps.obs.exp.ctgagt+cps.obs.exp.ctgata+cps.obs.exp.ctgatc+cps.obs.exp.ctgatg+cps.obs.exp.ctgatt+cps.obs.exp.ctgcaa+cps.obs.exp.ctgcac+cps.obs.exp.ctgcag+cps.obs.exp.ctgcat+cps.obs.exp.ctgcca+cps.obs.exp.ctgccc+cps.obs.exp.ctgccg+cps.obs.exp.ctgcct+cps.obs.exp.ctgcga+cps.obs.exp.ctgcgc+cps.obs.exp.ctgcgg+cps.obs.exp.ctgcgt+cps.obs.exp.ctgcta+cps.obs.exp.ctgctc+cps.obs.exp.ctgctg+cps.obs.exp.ctgctt+cps.obs.exp.ctggaa+cps.obs.exp.ctggac+cps.obs.exp.ctggag+cps.obs.exp.ctggat+cps.obs.exp.ctggca+cps.obs.exp.ctggcc+cps.obs.exp.ctggcg+cps.obs.exp.ctggct+cps.obs.exp.ctggga+cps.obs.exp.ctgggc+cps.obs.exp.ctgggg+cps.obs.exp.ctgggt+cps.obs.exp.ctggta+cps.obs.exp.ctggtc+cps.obs.exp.ctggtg+cps.obs.exp.ctggtt+cps.obs.exp.ctgtac+cps.obs.exp.ctgtat+cps.obs.exp.ctgtca+cps.obs.exp.ctgtcc+cps.obs.exp.ctgtcg+cps.obs.exp.ctgtct+cps.obs.exp.ctgtgc+cps.obs.exp.ctgtgg+cps.obs.exp.ctgtgt+cps.obs.exp.ctgtta+cps.obs.exp.ctgttc+cps.obs.exp.ctgttg+cps.obs.exp.ctgttt+cps.obs.exp.cttaaa+cps.obs.exp.cttaac+cps.obs.exp.cttaag+cps.obs.exp.cttaat+cps.obs.exp.cttaca+cps.obs.exp.cttacc+cps.obs.exp.cttacg+cps.obs.exp.cttact+cps.obs.exp.cttaga+cps.obs.exp.cttagc+cps.obs.exp.cttagg+cps.obs.exp.cttagt+cps.obs.exp.cttata+cps.obs.exp.cttatc+cps.obs.exp.cttatg+cps.obs.exp.cttatt+cps.obs.exp.cttcaa+cps.obs.exp.cttcac+cps.obs.exp.cttcag+cps.obs.exp.cttcat+cps.obs.exp.cttcca+cps.obs.exp.cttccc+cps.obs.exp.cttccg+cps.obs.exp.cttcct+cps.obs.exp.cttcga+cps.obs.exp.cttcgc+cps.obs.exp.cttcgg+cps.obs.exp.cttcgt+cps.obs.exp.cttcta+cps.obs.exp.cttctc+cps.obs.exp.cttctg+cps.obs.exp.cttctt+cps.obs.exp.cttgaa+cps.obs.exp.cttgac+cps.obs.exp.cttgag+cps.obs.exp.cttgat+cps.obs.exp.cttgca+cps.obs.exp.cttgcc+cps.obs.exp.cttgcg+cps.obs.exp.cttgct+cps.obs.exp.cttgga+cps.obs.exp.cttggc+cps.obs.exp.cttggg+cps.obs.exp.cttggt+cps.obs.exp.cttgta+cps.obs.exp.cttgtc+cps.obs.exp.cttgtg+cps.obs.exp.cttgtt+cps.obs.exp.ctttac+cps.obs.exp.ctttat+cps.obs.exp.ctttca+cps.obs.exp.ctttcc+cps.obs.exp.ctttcg+cps.obs.exp.ctttct+cps.obs.exp.ctttgc+cps.obs.exp.ctttgg+cps.obs.exp.ctttgt+cps.obs.exp.ctttta+cps.obs.exp.cttttc+cps.obs.exp.cttttg+cps.obs.exp.cttttt+cps.obs.exp.gaaaaa+cps.obs.exp.gaaaac+cps.obs.exp.gaaaag+cps.obs.exp.gaaaat+cps.obs.exp.gaaaca+cps.obs.exp.gaaacc+cps.obs.exp.gaaacg+cps.obs.exp.gaaact+cps.obs.exp.gaaaga+cps.obs.exp.gaaagc+cps.obs.exp.gaaagg+cps.obs.exp.gaaagt+cps.obs.exp.gaaata+cps.obs.exp.gaaatc+cps.obs.exp.gaaatg+cps.obs.exp.gaaatt+cps.obs.exp.gaacaa+cps.obs.exp.gaacac+cps.obs.exp.gaacag+cps.obs.exp.gaacat+cps.obs.exp.gaacca+cps.obs.exp.gaaccc+cps.obs.exp.gaaccg+cps.obs.exp.gaacct+cps.obs.exp.gaacga+cps.obs.exp.gaacgc+cps.obs.exp.gaacgg+cps.obs.exp.gaacgt+cps.obs.exp.gaacta+cps.obs.exp.gaactc+cps.obs.exp.gaactg+cps.obs.exp.gaactt+cps.obs.exp.gaagaa+cps.obs.exp.gaagac+cps.obs.exp.gaagag+cps.obs.exp.gaagat+cps.obs.exp.gaagca+cps.obs.exp.gaagcc+cps.obs.exp.gaagcg+cps.obs.exp.gaagct+cps.obs.exp.gaagga+cps.obs.exp.gaaggc+cps.obs.exp.gaaggg+cps.obs.exp.gaaggt+cps.obs.exp.gaagta+cps.obs.exp.gaagtc+cps.obs.exp.gaagtg+cps.obs.exp.gaagtt+cps.obs.exp.gaatac+cps.obs.exp.gaatat+cps.obs.exp.gaatca+cps.obs.exp.gaatcc+cps.obs.exp.gaatcg+cps.obs.exp.gaatct+cps.obs.exp.gaatgc+cps.obs.exp.gaatgg+cps.obs.exp.gaatgt+cps.obs.exp.gaatta+cps.obs.exp.gaattc+cps.obs.exp.gaattg+cps.obs.exp.gaattt+cps.obs.exp.gacaaa+cps.obs.exp.gacaac+cps.obs.exp.gacaag+cps.obs.exp.gacaat+cps.obs.exp.gacaca+cps.obs.exp.gacacc+cps.obs.exp.gacacg+cps.obs.exp.gacact+cps.obs.exp.gacaga+cps.obs.exp.gacagc+cps.obs.exp.gacagg+cps.obs.exp.gacagt+cps.obs.exp.gacata+cps.obs.exp.gacatc+cps.obs.exp.gacatg+cps.obs.exp.gacatt+cps.obs.exp.gaccaa+cps.obs.exp.gaccac+cps.obs.exp.gaccag+cps.obs.exp.gaccat+cps.obs.exp.gaccca+cps.obs.exp.gacccc+cps.obs.exp.gacccg+cps.obs.exp.gaccct+cps.obs.exp.gaccga+cps.obs.exp.gaccgc+cps.obs.exp.gaccgg+cps.obs.exp.gaccgt+cps.obs.exp.gaccta+cps.obs.exp.gacctc+cps.obs.exp.gacctg+cps.obs.exp.gacctt+cps.obs.exp.gacgaa+cps.obs.exp.gacgac+cps.obs.exp.gacgag+cps.obs.exp.gacgat+cps.obs.exp.gacgca+cps.obs.exp.gacgcc+cps.obs.exp.gacgcg+cps.obs.exp.gacgct+cps.obs.exp.gacgga+cps.obs.exp.gacggc+cps.obs.exp.gacggg+cps.obs.exp.gacggt+cps.obs.exp.gacgta+cps.obs.exp.gacgtc+cps.obs.exp.gacgtg+cps.obs.exp.gacgtt+cps.obs.exp.gactac+cps.obs.exp.gactat+cps.obs.exp.gactca+cps.obs.exp.gactcc+cps.obs.exp.gactcg+cps.obs.exp.gactct+cps.obs.exp.gactgc+cps.obs.exp.gactgg+cps.obs.exp.gactgt+cps.obs.exp.gactta+cps.obs.exp.gacttc+cps.obs.exp.gacttg+cps.obs.exp.gacttt+cps.obs.exp.gagaaa+cps.obs.exp.gagaac+cps.obs.exp.gagaag+cps.obs.exp.gagaat+cps.obs.exp.gagaca+cps.obs.exp.gagacc+cps.obs.exp.gagacg+cps.obs.exp.gagact+cps.obs.exp.gagaga+cps.obs.exp.gagagc+cps.obs.exp.gagagg+cps.obs.exp.gagagt+cps.obs.exp.gagata+cps.obs.exp.gagatc+cps.obs.exp.gagatg+cps.obs.exp.gagatt+cps.obs.exp.gagcaa+cps.obs.exp.gagcac+cps.obs.exp.gagcag+cps.obs.exp.gagcat+cps.obs.exp.gagcca+cps.obs.exp.gagccc+cps.obs.exp.gagccg+cps.obs.exp.gagcct+cps.obs.exp.gagcga+cps.obs.exp.gagcgc+cps.obs.exp.gagcgg+cps.obs.exp.gagcgt+cps.obs.exp.gagcta+cps.obs.exp.gagctc+cps.obs.exp.gagctg+cps.obs.exp.gagctt+cps.obs.exp.gaggaa+cps.obs.exp.gaggac+cps.obs.exp.gaggag+cps.obs.exp.gaggat+cps.obs.exp.gaggca+cps.obs.exp.gaggcc+cps.obs.exp.gaggcg+cps.obs.exp.gaggct+cps.obs.exp.gaggga+cps.obs.exp.gagggc+cps.obs.exp.gagggg+cps.obs.exp.gagggt+cps.obs.exp.gaggta+cps.obs.exp.gaggtc+cps.obs.exp.gaggtg+cps.obs.exp.gaggtt+cps.obs.exp.gagtac+cps.obs.exp.gagtat+cps.obs.exp.gagtca+cps.obs.exp.gagtcc+cps.obs.exp.gagtcg+cps.obs.exp.gagtct+cps.obs.exp.gagtgc+cps.obs.exp.gagtgg+cps.obs.exp.gagtgt+cps.obs.exp.gagtta+cps.obs.exp.gagttc+cps.obs.exp.gagttg+cps.obs.exp.gagttt+cps.obs.exp.gataaa+cps.obs.exp.gataac+cps.obs.exp.gataag+cps.obs.exp.gataat+cps.obs.exp.gataca+cps.obs.exp.gatacc+cps.obs.exp.gatacg+cps.obs.exp.gatact+cps.obs.exp.gataga+cps.obs.exp.gatagc+cps.obs.exp.gatagg+cps.obs.exp.gatagt+cps.obs.exp.gatata+cps.obs.exp.gatatc+cps.obs.exp.gatatg+cps.obs.exp.gatatt+cps.obs.exp.gatcaa+cps.obs.exp.gatcac+cps.obs.exp.gatcag+cps.obs.exp.gatcat+cps.obs.exp.gatcca+cps.obs.exp.gatccc+cps.obs.exp.gatccg+cps.obs.exp.gatcct+cps.obs.exp.gatcga+cps.obs.exp.gatcgc+cps.obs.exp.gatcgg+cps.obs.exp.gatcgt+cps.obs.exp.gatcta+cps.obs.exp.gatctc+cps.obs.exp.gatctg+cps.obs.exp.gatctt+cps.obs.exp.gatgaa+cps.obs.exp.gatgac+cps.obs.exp.gatgag+cps.obs.exp.gatgat+cps.obs.exp.gatgca+cps.obs.exp.gatgcc+cps.obs.exp.gatgcg+cps.obs.exp.gatgct+cps.obs.exp.gatgga+cps.obs.exp.gatggc+cps.obs.exp.gatggg+cps.obs.exp.gatggt+cps.obs.exp.gatgta+cps.obs.exp.gatgtc+cps.obs.exp.gatgtg+cps.obs.exp.gatgtt+cps.obs.exp.gattac+cps.obs.exp.gattat+cps.obs.exp.gattca+cps.obs.exp.gattcc+cps.obs.exp.gattcg+cps.obs.exp.gattct+cps.obs.exp.gattgc+cps.obs.exp.gattgg+cps.obs.exp.gattgt+cps.obs.exp.gattta+cps.obs.exp.gatttc+cps.obs.exp.gatttg+cps.obs.exp.gatttt+cps.obs.exp.gcaaaa+cps.obs.exp.gcaaac+cps.obs.exp.gcaaag+cps.obs.exp.gcaaat+cps.obs.exp.gcaaca+cps.obs.exp.gcaacc+cps.obs.exp.gcaacg+cps.obs.exp.gcaact+cps.obs.exp.gcaaga+cps.obs.exp.gcaagc+cps.obs.exp.gcaagg+cps.obs.exp.gcaagt+cps.obs.exp.gcaata+cps.obs.exp.gcaatc+cps.obs.exp.gcaatg+cps.obs.exp.gcaatt+cps.obs.exp.gcacaa+cps.obs.exp.gcacac+cps.obs.exp.gcacag+cps.obs.exp.gcacat+cps.obs.exp.gcacca+cps.obs.exp.gcaccc+cps.obs.exp.gcaccg+cps.obs.exp.gcacct+cps.obs.exp.gcacga+cps.obs.exp.gcacgc+cps.obs.exp.gcacgg+cps.obs.exp.gcacgt+cps.obs.exp.gcacta+cps.obs.exp.gcactc+cps.obs.exp.gcactg+cps.obs.exp.gcactt+cps.obs.exp.gcagaa+cps.obs.exp.gcagac+cps.obs.exp.gcagag+cps.obs.exp.gcagat+cps.obs.exp.gcagca+cps.obs.exp.gcagcc+cps.obs.exp.gcagcg+cps.obs.exp.gcagct+cps.obs.exp.gcagga+cps.obs.exp.gcaggc+cps.obs.exp.gcaggg+cps.obs.exp.gcaggt+cps.obs.exp.gcagta+cps.obs.exp.gcagtc+cps.obs.exp.gcagtg+cps.obs.exp.gcagtt+cps.obs.exp.gcatac+cps.obs.exp.gcatat+cps.obs.exp.gcatca+cps.obs.exp.gcatcc+cps.obs.exp.gcatcg+cps.obs.exp.gcatct+cps.obs.exp.gcatgc+cps.obs.exp.gcatgg+cps.obs.exp.gcatgt+cps.obs.exp.gcatta+cps.obs.exp.gcattc+cps.obs.exp.gcattg+cps.obs.exp.gcattt+cps.obs.exp.gccaaa+cps.obs.exp.gccaac+cps.obs.exp.gccaag+cps.obs.exp.gccaat+cps.obs.exp.gccaca+cps.obs.exp.gccacc+cps.obs.exp.gccacg+cps.obs.exp.gccact+cps.obs.exp.gccaga+cps.obs.exp.gccagc+cps.obs.exp.gccagg+cps.obs.exp.gccagt+cps.obs.exp.gccata+cps.obs.exp.gccatc+cps.obs.exp.gccatg+cps.obs.exp.gccatt+cps.obs.exp.gcccaa+cps.obs.exp.gcccac+cps.obs.exp.gcccag+cps.obs.exp.gcccat+cps.obs.exp.gcccca+cps.obs.exp.gccccc+cps.obs.exp.gccccg+cps.obs.exp.gcccct+cps.obs.exp.gcccga+cps.obs.exp.gcccgc+cps.obs.exp.gcccgg+cps.obs.exp.gcccgt+cps.obs.exp.gcccta+cps.obs.exp.gccctc+cps.obs.exp.gccctg+cps.obs.exp.gccctt+cps.obs.exp.gccgaa+cps.obs.exp.gccgac+cps.obs.exp.gccgag+cps.obs.exp.gccgat+cps.obs.exp.gccgca+cps.obs.exp.gccgcc+cps.obs.exp.gccgcg+cps.obs.exp.gccgct+cps.obs.exp.gccgga+cps.obs.exp.gccggc+cps.obs.exp.gccggg+cps.obs.exp.gccggt+cps.obs.exp.gccgta+cps.obs.exp.gccgtc+cps.obs.exp.gccgtg+cps.obs.exp.gccgtt+cps.obs.exp.gcctac+cps.obs.exp.gcctat+cps.obs.exp.gcctca+cps.obs.exp.gcctcc+cps.obs.exp.gcctcg+cps.obs.exp.gcctct+cps.obs.exp.gcctgc+cps.obs.exp.gcctgg+cps.obs.exp.gcctgt+cps.obs.exp.gcctta+cps.obs.exp.gccttc+cps.obs.exp.gccttg+cps.obs.exp.gccttt+cps.obs.exp.gcgaaa+cps.obs.exp.gcgaac+cps.obs.exp.gcgaag+cps.obs.exp.gcgaat+cps.obs.exp.gcgaca+cps.obs.exp.gcgacc+cps.obs.exp.gcgacg+cps.obs.exp.gcgact+cps.obs.exp.gcgaga+cps.obs.exp.gcgagc+cps.obs.exp.gcgagg+cps.obs.exp.gcgagt+cps.obs.exp.gcgata+cps.obs.exp.gcgatc+cps.obs.exp.gcgatg+cps.obs.exp.gcgatt+cps.obs.exp.gcgcaa+cps.obs.exp.gcgcac+cps.obs.exp.gcgcag+cps.obs.exp.gcgcat+cps.obs.exp.gcgcca+cps.obs.exp.gcgccc+cps.obs.exp.gcgccg+cps.obs.exp.gcgcct+cps.obs.exp.gcgcga+cps.obs.exp.gcgcgc+cps.obs.exp.gcgcgg+cps.obs.exp.gcgcgt+cps.obs.exp.gcgcta+cps.obs.exp.gcgctc+cps.obs.exp.gcgctg+cps.obs.exp.gcgctt+cps.obs.exp.gcggaa+cps.obs.exp.gcggac+cps.obs.exp.gcggag+cps.obs.exp.gcggat+cps.obs.exp.gcggca+cps.obs.exp.gcggcc+cps.obs.exp.gcggcg+cps.obs.exp.gcggct+cps.obs.exp.gcggga+cps.obs.exp.gcgggc+cps.obs.exp.gcgggg+cps.obs.exp.gcgggt+cps.obs.exp.gcggta+cps.obs.exp.gcggtc+cps.obs.exp.gcggtg+cps.obs.exp.gcggtt+cps.obs.exp.gcgtac+cps.obs.exp.gcgtat+cps.obs.exp.gcgtca+cps.obs.exp.gcgtcc+cps.obs.exp.gcgtcg+cps.obs.exp.gcgtct+cps.obs.exp.gcgtgc+cps.obs.exp.gcgtgg+cps.obs.exp.gcgtgt+cps.obs.exp.gcgtta+cps.obs.exp.gcgttc+cps.obs.exp.gcgttg+cps.obs.exp.gcgttt+cps.obs.exp.gctaaa+cps.obs.exp.gctaac+cps.obs.exp.gctaag+cps.obs.exp.gctaat+cps.obs.exp.gctaca+cps.obs.exp.gctacc+cps.obs.exp.gctacg+cps.obs.exp.gctact+cps.obs.exp.gctaga+cps.obs.exp.gctagc+cps.obs.exp.gctagg+cps.obs.exp.gctagt+cps.obs.exp.gctata+cps.obs.exp.gctatc+cps.obs.exp.gctatg+cps.obs.exp.gctatt+cps.obs.exp.gctcaa+cps.obs.exp.gctcac+cps.obs.exp.gctcag+cps.obs.exp.gctcat+cps.obs.exp.gctcca+cps.obs.exp.gctccc+cps.obs.exp.gctccg+cps.obs.exp.gctcct+cps.obs.exp.gctcga+cps.obs.exp.gctcgc+cps.obs.exp.gctcgg+cps.obs.exp.gctcgt+cps.obs.exp.gctcta+cps.obs.exp.gctctc+cps.obs.exp.gctctg+cps.obs.exp.gctctt+cps.obs.exp.gctgaa+cps.obs.exp.gctgac+cps.obs.exp.gctgag+cps.obs.exp.gctgat+cps.obs.exp.gctgca+cps.obs.exp.gctgcc+cps.obs.exp.gctgcg+cps.obs.exp.gctgct+cps.obs.exp.gctgga+cps.obs.exp.gctggc+cps.obs.exp.gctggg+cps.obs.exp.gctggt+cps.obs.exp.gctgta+cps.obs.exp.gctgtc+cps.obs.exp.gctgtg+cps.obs.exp.gctgtt+cps.obs.exp.gcttac+cps.obs.exp.gcttat+cps.obs.exp.gcttca+cps.obs.exp.gcttcc+cps.obs.exp.gcttcg+cps.obs.exp.gcttct+cps.obs.exp.gcttgc+cps.obs.exp.gcttgg+cps.obs.exp.gcttgt+cps.obs.exp.gcttta+cps.obs.exp.gctttc+cps.obs.exp.gctttg+cps.obs.exp.gctttt+cps.obs.exp.ggaaaa+cps.obs.exp.ggaaac+cps.obs.exp.ggaaag+cps.obs.exp.ggaaat+cps.obs.exp.ggaaca+cps.obs.exp.ggaacc+cps.obs.exp.ggaacg+cps.obs.exp.ggaact+cps.obs.exp.ggaaga+cps.obs.exp.ggaagc+cps.obs.exp.ggaagg+cps.obs.exp.ggaagt+cps.obs.exp.ggaata+cps.obs.exp.ggaatc+cps.obs.exp.ggaatg+cps.obs.exp.ggaatt+cps.obs.exp.ggacaa+cps.obs.exp.ggacac+cps.obs.exp.ggacag+cps.obs.exp.ggacat+cps.obs.exp.ggacca+cps.obs.exp.ggaccc+cps.obs.exp.ggaccg+cps.obs.exp.ggacct+cps.obs.exp.ggacga+cps.obs.exp.ggacgc+cps.obs.exp.ggacgg+cps.obs.exp.ggacgt+cps.obs.exp.ggacta+cps.obs.exp.ggactc+cps.obs.exp.ggactg+cps.obs.exp.ggactt+cps.obs.exp.ggagaa+cps.obs.exp.ggagac+cps.obs.exp.ggagag+cps.obs.exp.ggagat+cps.obs.exp.ggagca+cps.obs.exp.ggagcc+cps.obs.exp.ggagcg+cps.obs.exp.ggagct+cps.obs.exp.ggagga+cps.obs.exp.ggaggc+cps.obs.exp.ggaggg+cps.obs.exp.ggaggt+cps.obs.exp.ggagta+cps.obs.exp.ggagtc+cps.obs.exp.ggagtg+cps.obs.exp.ggagtt+cps.obs.exp.ggatac+cps.obs.exp.ggatat+cps.obs.exp.ggatca+cps.obs.exp.ggatcc+cps.obs.exp.ggatcg+cps.obs.exp.ggatct+cps.obs.exp.ggatgc+cps.obs.exp.ggatgg+cps.obs.exp.ggatgt+cps.obs.exp.ggatta+cps.obs.exp.ggattc+cps.obs.exp.ggattg+cps.obs.exp.ggattt+cps.obs.exp.ggcaaa+cps.obs.exp.ggcaac+cps.obs.exp.ggcaag+cps.obs.exp.ggcaat+cps.obs.exp.ggcaca+cps.obs.exp.ggcacc+cps.obs.exp.ggcacg+cps.obs.exp.ggcact+cps.obs.exp.ggcaga+cps.obs.exp.ggcagc+cps.obs.exp.ggcagg+cps.obs.exp.ggcagt+cps.obs.exp.ggcata+cps.obs.exp.ggcatc+cps.obs.exp.ggcatg+cps.obs.exp.ggcatt+cps.obs.exp.ggccaa+cps.obs.exp.ggccac+cps.obs.exp.ggccag+cps.obs.exp.ggccat+cps.obs.exp.ggccca+cps.obs.exp.ggcccc+cps.obs.exp.ggcccg+cps.obs.exp.ggccct+cps.obs.exp.ggccga+cps.obs.exp.ggccgc+cps.obs.exp.ggccgg+cps.obs.exp.ggccgt+cps.obs.exp.ggccta+cps.obs.exp.ggcctc+cps.obs.exp.ggcctg+cps.obs.exp.ggcctt+cps.obs.exp.ggcgaa+cps.obs.exp.ggcgac+cps.obs.exp.ggcgag+cps.obs.exp.ggcgat+cps.obs.exp.ggcgca+cps.obs.exp.ggcgcc+cps.obs.exp.ggcgcg+cps.obs.exp.ggcgct+cps.obs.exp.ggcgga+cps.obs.exp.ggcggc+cps.obs.exp.ggcggg+cps.obs.exp.ggcggt+cps.obs.exp.ggcgta+cps.obs.exp.ggcgtc+cps.obs.exp.ggcgtg+cps.obs.exp.ggcgtt+cps.obs.exp.ggctac+cps.obs.exp.ggctat+cps.obs.exp.ggctca+cps.obs.exp.ggctcc+cps.obs.exp.ggctcg+cps.obs.exp.ggctct+cps.obs.exp.ggctgc+cps.obs.exp.ggctgg+cps.obs.exp.ggctgt+cps.obs.exp.ggctta+cps.obs.exp.ggcttc+cps.obs.exp.ggcttg+cps.obs.exp.ggcttt+cps.obs.exp.gggaaa+cps.obs.exp.gggaac+cps.obs.exp.gggaag+cps.obs.exp.gggaat+cps.obs.exp.gggaca+cps.obs.exp.gggacc+cps.obs.exp.gggacg+cps.obs.exp.gggact+cps.obs.exp.gggaga+cps.obs.exp.gggagc+cps.obs.exp.gggagg+cps.obs.exp.gggagt+cps.obs.exp.gggata+cps.obs.exp.gggatc+cps.obs.exp.gggatg+cps.obs.exp.gggatt+cps.obs.exp.gggcaa+cps.obs.exp.gggcac+cps.obs.exp.gggcag+cps.obs.exp.gggcat+cps.obs.exp.gggcca+cps.obs.exp.gggccc+cps.obs.exp.gggccg+cps.obs.exp.gggcct+cps.obs.exp.gggcga+cps.obs.exp.gggcgc+cps.obs.exp.gggcgg+cps.obs.exp.gggcgt+cps.obs.exp.gggcta+cps.obs.exp.gggctc+cps.obs.exp.gggctg+cps.obs.exp.gggctt+cps.obs.exp.ggggaa+cps.obs.exp.ggggac+cps.obs.exp.ggggag+cps.obs.exp.ggggat+cps.obs.exp.ggggca+cps.obs.exp.ggggcc+cps.obs.exp.ggggcg+cps.obs.exp.ggggct+cps.obs.exp.ggggga+cps.obs.exp.gggggc+cps.obs.exp.gggggg+cps.obs.exp.gggggt+cps.obs.exp.ggggta+cps.obs.exp.ggggtc+cps.obs.exp.ggggtg+cps.obs.exp.ggggtt+cps.obs.exp.gggtac+cps.obs.exp.gggtat+cps.obs.exp.gggtca+cps.obs.exp.gggtcc+cps.obs.exp.gggtcg+cps.obs.exp.gggtct+cps.obs.exp.gggtgc+cps.obs.exp.gggtgg+cps.obs.exp.gggtgt+cps.obs.exp.gggtta+cps.obs.exp.gggttc+cps.obs.exp.gggttg+cps.obs.exp.gggttt+cps.obs.exp.ggtaaa+cps.obs.exp.ggtaac+cps.obs.exp.ggtaag+cps.obs.exp.ggtaat+cps.obs.exp.ggtaca+cps.obs.exp.ggtacc+cps.obs.exp.ggtacg+cps.obs.exp.ggtact+cps.obs.exp.ggtaga+cps.obs.exp.ggtagc+cps.obs.exp.ggtagg+cps.obs.exp.ggtagt+cps.obs.exp.ggtata+cps.obs.exp.ggtatc+cps.obs.exp.ggtatg+cps.obs.exp.ggtatt+cps.obs.exp.ggtcaa+cps.obs.exp.ggtcac+cps.obs.exp.ggtcag+cps.obs.exp.ggtcat+cps.obs.exp.ggtcca+cps.obs.exp.ggtccc+cps.obs.exp.ggtccg+cps.obs.exp.ggtcct+cps.obs.exp.ggtcga+cps.obs.exp.ggtcgc+cps.obs.exp.ggtcgg+cps.obs.exp.ggtcgt+cps.obs.exp.ggtcta+cps.obs.exp.ggtctc+cps.obs.exp.ggtctg+cps.obs.exp.ggtctt+cps.obs.exp.ggtgaa+cps.obs.exp.ggtgac+cps.obs.exp.ggtgag+cps.obs.exp.ggtgat+cps.obs.exp.ggtgca+cps.obs.exp.ggtgcc+cps.obs.exp.ggtgcg+cps.obs.exp.ggtgct+cps.obs.exp.ggtgga+cps.obs.exp.ggtggc+cps.obs.exp.ggtggg+cps.obs.exp.ggtggt+cps.obs.exp.ggtgta+cps.obs.exp.ggtgtc+cps.obs.exp.ggtgtg+cps.obs.exp.ggtgtt+cps.obs.exp.ggttac+cps.obs.exp.ggttat+cps.obs.exp.ggttca+cps.obs.exp.ggttcc+cps.obs.exp.ggttcg+cps.obs.exp.ggttct+cps.obs.exp.ggttgc+cps.obs.exp.ggttgg+cps.obs.exp.ggttgt+cps.obs.exp.ggttta+cps.obs.exp.ggtttc+cps.obs.exp.ggtttg+cps.obs.exp.ggtttt+cps.obs.exp.gtaaaa+cps.obs.exp.gtaaac+cps.obs.exp.gtaaag+cps.obs.exp.gtaaat+cps.obs.exp.gtaaca+cps.obs.exp.gtaacc+cps.obs.exp.gtaacg+cps.obs.exp.gtaact+cps.obs.exp.gtaaga+cps.obs.exp.gtaagc+cps.obs.exp.gtaagg+cps.obs.exp.gtaagt+cps.obs.exp.gtaata+cps.obs.exp.gtaatc+cps.obs.exp.gtaatg+cps.obs.exp.gtaatt+cps.obs.exp.gtacaa+cps.obs.exp.gtacac+cps.obs.exp.gtacag+cps.obs.exp.gtacat+cps.obs.exp.gtacca+cps.obs.exp.gtaccc+cps.obs.exp.gtaccg+cps.obs.exp.gtacct+cps.obs.exp.gtacga+cps.obs.exp.gtacgc+cps.obs.exp.gtacgg+cps.obs.exp.gtacgt+cps.obs.exp.gtacta+cps.obs.exp.gtactc+cps.obs.exp.gtactg+cps.obs.exp.gtactt+cps.obs.exp.gtagaa+cps.obs.exp.gtagac+cps.obs.exp.gtagag+cps.obs.exp.gtagat+cps.obs.exp.gtagca+cps.obs.exp.gtagcc+cps.obs.exp.gtagcg+cps.obs.exp.gtagct+cps.obs.exp.gtagga+cps.obs.exp.gtaggc+cps.obs.exp.gtaggg+cps.obs.exp.gtaggt+cps.obs.exp.gtagta+cps.obs.exp.gtagtc+cps.obs.exp.gtagtg+cps.obs.exp.gtagtt+cps.obs.exp.gtatac+cps.obs.exp.gtatat+cps.obs.exp.gtatca+cps.obs.exp.gtatcc+cps.obs.exp.gtatcg+cps.obs.exp.gtatct+cps.obs.exp.gtatgc+cps.obs.exp.gtatgg+cps.obs.exp.gtatgt+cps.obs.exp.gtatta+cps.obs.exp.gtattc+cps.obs.exp.gtattg+cps.obs.exp.gtattt+cps.obs.exp.gtcaaa+cps.obs.exp.gtcaac+cps.obs.exp.gtcaag+cps.obs.exp.gtcaat+cps.obs.exp.gtcaca+cps.obs.exp.gtcacc+cps.obs.exp.gtcacg+cps.obs.exp.gtcact+cps.obs.exp.gtcaga+cps.obs.exp.gtcagc+cps.obs.exp.gtcagg+cps.obs.exp.gtcagt+cps.obs.exp.gtcata+cps.obs.exp.gtcatc+cps.obs.exp.gtcatg+cps.obs.exp.gtcatt+cps.obs.exp.gtccaa+cps.obs.exp.gtccac+cps.obs.exp.gtccag+cps.obs.exp.gtccat+cps.obs.exp.gtccca+cps.obs.exp.gtcccc+cps.obs.exp.gtcccg+cps.obs.exp.gtccct+cps.obs.exp.gtccga+cps.obs.exp.gtccgc+cps.obs.exp.gtccgg+cps.obs.exp.gtccgt+cps.obs.exp.gtccta+cps.obs.exp.gtcctc+cps.obs.exp.gtcctg+cps.obs.exp.gtcctt+cps.obs.exp.gtcgaa+cps.obs.exp.gtcgac+cps.obs.exp.gtcgag+cps.obs.exp.gtcgat+cps.obs.exp.gtcgca+cps.obs.exp.gtcgcc+cps.obs.exp.gtcgcg+cps.obs.exp.gtcgct+cps.obs.exp.gtcgga+cps.obs.exp.gtcggc+cps.obs.exp.gtcggg+cps.obs.exp.gtcggt+cps.obs.exp.gtcgta+cps.obs.exp.gtcgtc+cps.obs.exp.gtcgtg+cps.obs.exp.gtcgtt+cps.obs.exp.gtctac+cps.obs.exp.gtctat+cps.obs.exp.gtctca+cps.obs.exp.gtctcc+cps.obs.exp.gtctcg+cps.obs.exp.gtctct+cps.obs.exp.gtctgc+cps.obs.exp.gtctgg+cps.obs.exp.gtctgt+cps.obs.exp.gtctta+cps.obs.exp.gtcttc+cps.obs.exp.gtcttg+cps.obs.exp.gtcttt+cps.obs.exp.gtgaaa+cps.obs.exp.gtgaac+cps.obs.exp.gtgaag+cps.obs.exp.gtgaat+cps.obs.exp.gtgaca+cps.obs.exp.gtgacc+cps.obs.exp.gtgacg+cps.obs.exp.gtgact+cps.obs.exp.gtgaga+cps.obs.exp.gtgagc+cps.obs.exp.gtgagg+cps.obs.exp.gtgagt+cps.obs.exp.gtgata+cps.obs.exp.gtgatc+cps.obs.exp.gtgatg+cps.obs.exp.gtgatt+cps.obs.exp.gtgcaa+cps.obs.exp.gtgcac+cps.obs.exp.gtgcag+cps.obs.exp.gtgcat+cps.obs.exp.gtgcca+cps.obs.exp.gtgccc+cps.obs.exp.gtgccg+cps.obs.exp.gtgcct+cps.obs.exp.gtgcga+cps.obs.exp.gtgcgc+cps.obs.exp.gtgcgg+cps.obs.exp.gtgcgt+cps.obs.exp.gtgcta+cps.obs.exp.gtgctc+cps.obs.exp.gtgctg+cps.obs.exp.gtgctt+cps.obs.exp.gtggaa+cps.obs.exp.gtggac+cps.obs.exp.gtggag+cps.obs.exp.gtggat+cps.obs.exp.gtggca+cps.obs.exp.gtggcc+cps.obs.exp.gtggcg+cps.obs.exp.gtggct+cps.obs.exp.gtggga+cps.obs.exp.gtgggc+cps.obs.exp.gtgggg+cps.obs.exp.gtgggt+cps.obs.exp.gtggta+cps.obs.exp.gtggtc+cps.obs.exp.gtggtg+cps.obs.exp.gtggtt+cps.obs.exp.gtgtac+cps.obs.exp.gtgtat+cps.obs.exp.gtgtca+cps.obs.exp.gtgtcc+cps.obs.exp.gtgtcg+cps.obs.exp.gtgtct+cps.obs.exp.gtgtgc+cps.obs.exp.gtgtgg+cps.obs.exp.gtgtgt+cps.obs.exp.gtgtta+cps.obs.exp.gtgttc+cps.obs.exp.gtgttg+cps.obs.exp.gtgttt+cps.obs.exp.gttaaa+cps.obs.exp.gttaac+cps.obs.exp.gttaag+cps.obs.exp.gttaat+cps.obs.exp.gttaca+cps.obs.exp.gttacc+cps.obs.exp.gttacg+cps.obs.exp.gttact+cps.obs.exp.gttaga+cps.obs.exp.gttagc+cps.obs.exp.gttagg+cps.obs.exp.gttagt+cps.obs.exp.gttata+cps.obs.exp.gttatc+cps.obs.exp.gttatg+cps.obs.exp.gttatt+cps.obs.exp.gttcaa+cps.obs.exp.gttcac+cps.obs.exp.gttcag+cps.obs.exp.gttcat+cps.obs.exp.gttcca+cps.obs.exp.gttccc+cps.obs.exp.gttccg+cps.obs.exp.gttcct+cps.obs.exp.gttcga+cps.obs.exp.gttcgc+cps.obs.exp.gttcgg+cps.obs.exp.gttcgt+cps.obs.exp.gttcta+cps.obs.exp.gttctc+cps.obs.exp.gttctg+cps.obs.exp.gttctt+cps.obs.exp.gttgaa+cps.obs.exp.gttgac+cps.obs.exp.gttgag+cps.obs.exp.gttgat+cps.obs.exp.gttgca+cps.obs.exp.gttgcc+cps.obs.exp.gttgcg+cps.obs.exp.gttgct+cps.obs.exp.gttgga+cps.obs.exp.gttggc+cps.obs.exp.gttggg+cps.obs.exp.gttggt+cps.obs.exp.gttgta+cps.obs.exp.gttgtc+cps.obs.exp.gttgtg+cps.obs.exp.gttgtt+cps.obs.exp.gtttac+cps.obs.exp.gtttat+cps.obs.exp.gtttca+cps.obs.exp.gtttcc+cps.obs.exp.gtttcg+cps.obs.exp.gtttct+cps.obs.exp.gtttgc+cps.obs.exp.gtttgg+cps.obs.exp.gtttgt+cps.obs.exp.gtttta+cps.obs.exp.gttttc+cps.obs.exp.gttttg+cps.obs.exp.gttttt+cps.obs.exp.tacaaa+cps.obs.exp.tacaac+cps.obs.exp.tacaag+cps.obs.exp.tacaat+cps.obs.exp.tacaca+cps.obs.exp.tacacc+cps.obs.exp.tacacg+cps.obs.exp.tacact+cps.obs.exp.tacaga+cps.obs.exp.tacagc+cps.obs.exp.tacagg+cps.obs.exp.tacagt+cps.obs.exp.tacata+cps.obs.exp.tacatc+cps.obs.exp.tacatg+cps.obs.exp.tacatt+cps.obs.exp.taccaa+cps.obs.exp.taccac+cps.obs.exp.taccag+cps.obs.exp.taccat+cps.obs.exp.taccca+cps.obs.exp.tacccc+cps.obs.exp.tacccg+cps.obs.exp.taccct+cps.obs.exp.taccga+cps.obs.exp.taccgc+cps.obs.exp.taccgg+cps.obs.exp.taccgt+cps.obs.exp.taccta+cps.obs.exp.tacctc+cps.obs.exp.tacctg+cps.obs.exp.tacctt+cps.obs.exp.tacgaa+cps.obs.exp.tacgac+cps.obs.exp.tacgag+cps.obs.exp.tacgat+cps.obs.exp.tacgca+cps.obs.exp.tacgcc+cps.obs.exp.tacgcg+cps.obs.exp.tacgct+cps.obs.exp.tacgga+cps.obs.exp.tacggc+cps.obs.exp.tacggg+cps.obs.exp.tacggt+cps.obs.exp.tacgta+cps.obs.exp.tacgtc+cps.obs.exp.tacgtg+cps.obs.exp.tacgtt+cps.obs.exp.tactac+cps.obs.exp.tactat+cps.obs.exp.tactca+cps.obs.exp.tactcc+cps.obs.exp.tactcg+cps.obs.exp.tactct+cps.obs.exp.tactgc+cps.obs.exp.tactgg+cps.obs.exp.tactgt+cps.obs.exp.tactta+cps.obs.exp.tacttc+cps.obs.exp.tacttg+cps.obs.exp.tacttt+cps.obs.exp.tataaa+cps.obs.exp.tataac+cps.obs.exp.tataag+cps.obs.exp.tataat+cps.obs.exp.tataca+cps.obs.exp.tatacc+cps.obs.exp.tatacg+cps.obs.exp.tatact+cps.obs.exp.tataga+cps.obs.exp.tatagc+cps.obs.exp.tatagg+cps.obs.exp.tatagt+cps.obs.exp.tatata+cps.obs.exp.tatatc+cps.obs.exp.tatatg+cps.obs.exp.tatatt+cps.obs.exp.tatcaa+cps.obs.exp.tatcac+cps.obs.exp.tatcag+cps.obs.exp.tatcat+cps.obs.exp.tatcca+cps.obs.exp.tatccc+cps.obs.exp.tatccg+cps.obs.exp.tatcct+cps.obs.exp.tatcga+cps.obs.exp.tatcgc+cps.obs.exp.tatcgg+cps.obs.exp.tatcgt+cps.obs.exp.tatcta+cps.obs.exp.tatctc+cps.obs.exp.tatctg+cps.obs.exp.tatctt+cps.obs.exp.tatgaa+cps.obs.exp.tatgac+cps.obs.exp.tatgag+cps.obs.exp.tatgat+cps.obs.exp.tatgca+cps.obs.exp.tatgcc+cps.obs.exp.tatgcg+cps.obs.exp.tatgct+cps.obs.exp.tatgga+cps.obs.exp.tatggc+cps.obs.exp.tatggg+cps.obs.exp.tatggt+cps.obs.exp.tatgta+cps.obs.exp.tatgtc+cps.obs.exp.tatgtg+cps.obs.exp.tatgtt+cps.obs.exp.tattac+cps.obs.exp.tattat+cps.obs.exp.tattca+cps.obs.exp.tattcc+cps.obs.exp.tattcg+cps.obs.exp.tattct+cps.obs.exp.tattgc+cps.obs.exp.tattgg+cps.obs.exp.tattgt+cps.obs.exp.tattta+cps.obs.exp.tatttc+cps.obs.exp.tatttg+cps.obs.exp.tatttt+cps.obs.exp.tcaaaa+cps.obs.exp.tcaaac+cps.obs.exp.tcaaag+cps.obs.exp.tcaaat+cps.obs.exp.tcaaca+cps.obs.exp.tcaacc+cps.obs.exp.tcaacg+cps.obs.exp.tcaact+cps.obs.exp.tcaaga+cps.obs.exp.tcaagc+cps.obs.exp.tcaagg+cps.obs.exp.tcaagt+cps.obs.exp.tcaata+cps.obs.exp.tcaatc+cps.obs.exp.tcaatg+cps.obs.exp.tcaatt+cps.obs.exp.tcacaa+cps.obs.exp.tcacac+cps.obs.exp.tcacag+cps.obs.exp.tcacat+cps.obs.exp.tcacca+cps.obs.exp.tcaccc+cps.obs.exp.tcaccg+cps.obs.exp.tcacct+cps.obs.exp.tcacga+cps.obs.exp.tcacgc+cps.obs.exp.tcacgg+cps.obs.exp.tcacgt+cps.obs.exp.tcacta+cps.obs.exp.tcactc+cps.obs.exp.tcactg+cps.obs.exp.tcactt+cps.obs.exp.tcagaa+cps.obs.exp.tcagac+cps.obs.exp.tcagag+cps.obs.exp.tcagat+cps.obs.exp.tcagca+cps.obs.exp.tcagcc+cps.obs.exp.tcagcg+cps.obs.exp.tcagct+cps.obs.exp.tcagga+cps.obs.exp.tcaggc+cps.obs.exp.tcaggg+cps.obs.exp.tcaggt+cps.obs.exp.tcagta+cps.obs.exp.tcagtc+cps.obs.exp.tcagtg+cps.obs.exp.tcagtt+cps.obs.exp.tcatac+cps.obs.exp.tcatat+cps.obs.exp.tcatca+cps.obs.exp.tcatcc+cps.obs.exp.tcatcg+cps.obs.exp.tcatct+cps.obs.exp.tcatgc+cps.obs.exp.tcatgg+cps.obs.exp.tcatgt+cps.obs.exp.tcatta+cps.obs.exp.tcattc+cps.obs.exp.tcattg+cps.obs.exp.tcattt+cps.obs.exp.tccaaa+cps.obs.exp.tccaac+cps.obs.exp.tccaag+cps.obs.exp.tccaat+cps.obs.exp.tccaca+cps.obs.exp.tccacc+cps.obs.exp.tccacg+cps.obs.exp.tccact+cps.obs.exp.tccaga+cps.obs.exp.tccagc+cps.obs.exp.tccagg+cps.obs.exp.tccagt+cps.obs.exp.tccata+cps.obs.exp.tccatc+cps.obs.exp.tccatg+cps.obs.exp.tccatt+cps.obs.exp.tcccaa+cps.obs.exp.tcccac+cps.obs.exp.tcccag+cps.obs.exp.tcccat+cps.obs.exp.tcccca+cps.obs.exp.tccccc+cps.obs.exp.tccccg+cps.obs.exp.tcccct+cps.obs.exp.tcccga+cps.obs.exp.tcccgc+cps.obs.exp.tcccgg+cps.obs.exp.tcccgt+cps.obs.exp.tcccta+cps.obs.exp.tccctc+cps.obs.exp.tccctg+cps.obs.exp.tccctt+cps.obs.exp.tccgaa+cps.obs.exp.tccgac+cps.obs.exp.tccgag+cps.obs.exp.tccgat+cps.obs.exp.tccgca+cps.obs.exp.tccgcc+cps.obs.exp.tccgcg+cps.obs.exp.tccgct+cps.obs.exp.tccgga+cps.obs.exp.tccggc+cps.obs.exp.tccggg+cps.obs.exp.tccggt+cps.obs.exp.tccgta+cps.obs.exp.tccgtc+cps.obs.exp.tccgtg+cps.obs.exp.tccgtt+cps.obs.exp.tcctac+cps.obs.exp.tcctat+cps.obs.exp.tcctca+cps.obs.exp.tcctcc+cps.obs.exp.tcctcg+cps.obs.exp.tcctct+cps.obs.exp.tcctgc+cps.obs.exp.tcctgg+cps.obs.exp.tcctgt+cps.obs.exp.tcctta+cps.obs.exp.tccttc+cps.obs.exp.tccttg+cps.obs.exp.tccttt+cps.obs.exp.tcgaaa+cps.obs.exp.tcgaac+cps.obs.exp.tcgaag+cps.obs.exp.tcgaat+cps.obs.exp.tcgaca+cps.obs.exp.tcgacc+cps.obs.exp.tcgacg+cps.obs.exp.tcgact+cps.obs.exp.tcgaga+cps.obs.exp.tcgagc+cps.obs.exp.tcgagg+cps.obs.exp.tcgagt+cps.obs.exp.tcgata+cps.obs.exp.tcgatc+cps.obs.exp.tcgatg+cps.obs.exp.tcgatt+cps.obs.exp.tcgcaa+cps.obs.exp.tcgcac+cps.obs.exp.tcgcag+cps.obs.exp.tcgcat+cps.obs.exp.tcgcca+cps.obs.exp.tcgccc+cps.obs.exp.tcgccg+cps.obs.exp.tcgcct+cps.obs.exp.tcgcga+cps.obs.exp.tcgcgc+cps.obs.exp.tcgcgg+cps.obs.exp.tcgcgt+cps.obs.exp.tcgcta+cps.obs.exp.tcgctc+cps.obs.exp.tcgctg+cps.obs.exp.tcgctt+cps.obs.exp.tcggaa+cps.obs.exp.tcggac+cps.obs.exp.tcggag+cps.obs.exp.tcggat+cps.obs.exp.tcggca+cps.obs.exp.tcggcc+cps.obs.exp.tcggcg+cps.obs.exp.tcggct+cps.obs.exp.tcggga+cps.obs.exp.tcgggc+cps.obs.exp.tcgggg+cps.obs.exp.tcgggt+cps.obs.exp.tcggta+cps.obs.exp.tcggtc+cps.obs.exp.tcggtg+cps.obs.exp.tcggtt+cps.obs.exp.tcgtac+cps.obs.exp.tcgtat+cps.obs.exp.tcgtca+cps.obs.exp.tcgtcc+cps.obs.exp.tcgtcg+cps.obs.exp.tcgtct+cps.obs.exp.tcgtgc+cps.obs.exp.tcgtgg+cps.obs.exp.tcgtgt+cps.obs.exp.tcgtta+cps.obs.exp.tcgttc+cps.obs.exp.tcgttg+cps.obs.exp.tcgttt+cps.obs.exp.tctaaa+cps.obs.exp.tctaac+cps.obs.exp.tctaag+cps.obs.exp.tctaat+cps.obs.exp.tctaca+cps.obs.exp.tctacc+cps.obs.exp.tctacg+cps.obs.exp.tctact+cps.obs.exp.tctaga+cps.obs.exp.tctagc+cps.obs.exp.tctagg+cps.obs.exp.tctagt+cps.obs.exp.tctata+cps.obs.exp.tctatc+cps.obs.exp.tctatg+cps.obs.exp.tctatt+cps.obs.exp.tctcaa+cps.obs.exp.tctcac+cps.obs.exp.tctcag+cps.obs.exp.tctcat+cps.obs.exp.tctcca+cps.obs.exp.tctccc+cps.obs.exp.tctccg+cps.obs.exp.tctcct+cps.obs.exp.tctcga+cps.obs.exp.tctcgc+cps.obs.exp.tctcgg+cps.obs.exp.tctcgt+cps.obs.exp.tctcta+cps.obs.exp.tctctc+cps.obs.exp.tctctg+cps.obs.exp.tctctt+cps.obs.exp.tctgaa+cps.obs.exp.tctgac+cps.obs.exp.tctgag+cps.obs.exp.tctgat+cps.obs.exp.tctgca+cps.obs.exp.tctgcc+cps.obs.exp.tctgcg+cps.obs.exp.tctgct+cps.obs.exp.tctgga+cps.obs.exp.tctggc+cps.obs.exp.tctggg+cps.obs.exp.tctggt+cps.obs.exp.tctgta+cps.obs.exp.tctgtc+cps.obs.exp.tctgtg+cps.obs.exp.tctgtt+cps.obs.exp.tcttac+cps.obs.exp.tcttat+cps.obs.exp.tcttca+cps.obs.exp.tcttcc+cps.obs.exp.tcttcg+cps.obs.exp.tcttct+cps.obs.exp.tcttgc+cps.obs.exp.tcttgg+cps.obs.exp.tcttgt+cps.obs.exp.tcttta+cps.obs.exp.tctttc+cps.obs.exp.tctttg+cps.obs.exp.tctttt+cps.obs.exp.tgcaaa+cps.obs.exp.tgcaac+cps.obs.exp.tgcaag+cps.obs.exp.tgcaat+cps.obs.exp.tgcaca+cps.obs.exp.tgcacc+cps.obs.exp.tgcacg+cps.obs.exp.tgcact+cps.obs.exp.tgcaga+cps.obs.exp.tgcagc+cps.obs.exp.tgcagg+cps.obs.exp.tgcagt+cps.obs.exp.tgcata+cps.obs.exp.tgcatc+cps.obs.exp.tgcatg+cps.obs.exp.tgcatt+cps.obs.exp.tgccaa+cps.obs.exp.tgccac+cps.obs.exp.tgccag+cps.obs.exp.tgccat+cps.obs.exp.tgccca+cps.obs.exp.tgcccc+cps.obs.exp.tgcccg+cps.obs.exp.tgccct+cps.obs.exp.tgccga+cps.obs.exp.tgccgc+cps.obs.exp.tgccgg+cps.obs.exp.tgccgt+cps.obs.exp.tgccta+cps.obs.exp.tgcctc+cps.obs.exp.tgcctg+cps.obs.exp.tgcctt+cps.obs.exp.tgcgaa+cps.obs.exp.tgcgac+cps.obs.exp.tgcgag+cps.obs.exp.tgcgat+cps.obs.exp.tgcgca+cps.obs.exp.tgcgcc+cps.obs.exp.tgcgcg+cps.obs.exp.tgcgct+cps.obs.exp.tgcgga+cps.obs.exp.tgcggc+cps.obs.exp.tgcggg+cps.obs.exp.tgcggt+cps.obs.exp.tgcgta+cps.obs.exp.tgcgtc+cps.obs.exp.tgcgtg+cps.obs.exp.tgcgtt+cps.obs.exp.tgctac+cps.obs.exp.tgctat+cps.obs.exp.tgctca+cps.obs.exp.tgctcc+cps.obs.exp.tgctcg+cps.obs.exp.tgctct+cps.obs.exp.tgctgc+cps.obs.exp.tgctgg+cps.obs.exp.tgctgt+cps.obs.exp.tgctta+cps.obs.exp.tgcttc+cps.obs.exp.tgcttg+cps.obs.exp.tgcttt+cps.obs.exp.tggaaa+cps.obs.exp.tggaac+cps.obs.exp.tggaag+cps.obs.exp.tggaat+cps.obs.exp.tggaca+cps.obs.exp.tggacc+cps.obs.exp.tggacg+cps.obs.exp.tggact+cps.obs.exp.tggaga+cps.obs.exp.tggagc+cps.obs.exp.tggagg+cps.obs.exp.tggagt+cps.obs.exp.tggata+cps.obs.exp.tggatc+cps.obs.exp.tggatg+cps.obs.exp.tggatt+cps.obs.exp.tggcaa+cps.obs.exp.tggcac+cps.obs.exp.tggcag+cps.obs.exp.tggcat+cps.obs.exp.tggcca+cps.obs.exp.tggccc+cps.obs.exp.tggccg+cps.obs.exp.tggcct+cps.obs.exp.tggcga+cps.obs.exp.tggcgc+cps.obs.exp.tggcgg+cps.obs.exp.tggcgt+cps.obs.exp.tggcta+cps.obs.exp.tggctc+cps.obs.exp.tggctg+cps.obs.exp.tggctt+cps.obs.exp.tgggaa+cps.obs.exp.tgggac+cps.obs.exp.tgggag+cps.obs.exp.tgggat+cps.obs.exp.tgggca+cps.obs.exp.tgggcc+cps.obs.exp.tgggcg+cps.obs.exp.tgggct+cps.obs.exp.tgggga+cps.obs.exp.tggggc+cps.obs.exp.tggggg+cps.obs.exp.tggggt+cps.obs.exp.tgggta+cps.obs.exp.tgggtc+cps.obs.exp.tgggtg+cps.obs.exp.tgggtt+cps.obs.exp.tggtac+cps.obs.exp.tggtat+cps.obs.exp.tggtca+cps.obs.exp.tggtcc+cps.obs.exp.tggtcg+cps.obs.exp.tggtct+cps.obs.exp.tggtgc+cps.obs.exp.tggtgg+cps.obs.exp.tggtgt+cps.obs.exp.tggtta+cps.obs.exp.tggttc+cps.obs.exp.tggttg+cps.obs.exp.tggttt+cps.obs.exp.tgtaaa+cps.obs.exp.tgtaac+cps.obs.exp.tgtaag+cps.obs.exp.tgtaat+cps.obs.exp.tgtaca+cps.obs.exp.tgtacc+cps.obs.exp.tgtacg+cps.obs.exp.tgtact+cps.obs.exp.tgtaga+cps.obs.exp.tgtagc+cps.obs.exp.tgtagg+cps.obs.exp.tgtagt+cps.obs.exp.tgtata+cps.obs.exp.tgtatc+cps.obs.exp.tgtatg+cps.obs.exp.tgtatt+cps.obs.exp.tgtcaa+cps.obs.exp.tgtcac+cps.obs.exp.tgtcag+cps.obs.exp.tgtcat+cps.obs.exp.tgtcca+cps.obs.exp.tgtccc+cps.obs.exp.tgtccg+cps.obs.exp.tgtcct+cps.obs.exp.tgtcga+cps.obs.exp.tgtcgc+cps.obs.exp.tgtcgg+cps.obs.exp.tgtcgt+cps.obs.exp.tgtcta+cps.obs.exp.tgtctc+cps.obs.exp.tgtctg+cps.obs.exp.tgtctt+cps.obs.exp.tgtgaa+cps.obs.exp.tgtgac+cps.obs.exp.tgtgag+cps.obs.exp.tgtgat+cps.obs.exp.tgtgca+cps.obs.exp.tgtgcc+cps.obs.exp.tgtgcg+cps.obs.exp.tgtgct+cps.obs.exp.tgtgga+cps.obs.exp.tgtggc+cps.obs.exp.tgtggg+cps.obs.exp.tgtggt+cps.obs.exp.tgtgta+cps.obs.exp.tgtgtc+cps.obs.exp.tgtgtg+cps.obs.exp.tgtgtt+cps.obs.exp.tgttac+cps.obs.exp.tgttat+cps.obs.exp.tgttca+cps.obs.exp.tgttcc+cps.obs.exp.tgttcg+cps.obs.exp.tgttct+cps.obs.exp.tgttgc+cps.obs.exp.tgttgg+cps.obs.exp.tgttgt+cps.obs.exp.tgttta+cps.obs.exp.tgtttc+cps.obs.exp.tgtttg+cps.obs.exp.tgtttt+cps.obs.exp.ttaaaa+cps.obs.exp.ttaaac+cps.obs.exp.ttaaag+cps.obs.exp.ttaaat+cps.obs.exp.ttaaca+cps.obs.exp.ttaacc+cps.obs.exp.ttaacg+cps.obs.exp.ttaact+cps.obs.exp.ttaaga+cps.obs.exp.ttaagc+cps.obs.exp.ttaagg+cps.obs.exp.ttaagt+cps.obs.exp.ttaata+cps.obs.exp.ttaatc+cps.obs.exp.ttaatg+cps.obs.exp.ttaatt+cps.obs.exp.ttacaa+cps.obs.exp.ttacac+cps.obs.exp.ttacag+cps.obs.exp.ttacat+cps.obs.exp.ttacca+cps.obs.exp.ttaccc+cps.obs.exp.ttaccg+cps.obs.exp.ttacct+cps.obs.exp.ttacga+cps.obs.exp.ttacgc+cps.obs.exp.ttacgg+cps.obs.exp.ttacgt+cps.obs.exp.ttacta+cps.obs.exp.ttactc+cps.obs.exp.ttactg+cps.obs.exp.ttactt+cps.obs.exp.ttagaa+cps.obs.exp.ttagac+cps.obs.exp.ttagag+cps.obs.exp.ttagat+cps.obs.exp.ttagca+cps.obs.exp.ttagcc+cps.obs.exp.ttagcg+cps.obs.exp.ttagct+cps.obs.exp.ttagga+cps.obs.exp.ttaggc+cps.obs.exp.ttaggg+cps.obs.exp.ttaggt+cps.obs.exp.ttagta+cps.obs.exp.ttagtc+cps.obs.exp.ttagtg+cps.obs.exp.ttagtt+cps.obs.exp.ttatac+cps.obs.exp.ttatat+cps.obs.exp.ttatca+cps.obs.exp.ttatcc+cps.obs.exp.ttatcg+cps.obs.exp.ttatct+cps.obs.exp.ttatgc+cps.obs.exp.ttatgg+cps.obs.exp.ttatgt+cps.obs.exp.ttatta+cps.obs.exp.ttattc+cps.obs.exp.ttattg+cps.obs.exp.ttattt+cps.obs.exp.ttcaaa+cps.obs.exp.ttcaac+cps.obs.exp.ttcaag+cps.obs.exp.ttcaat+cps.obs.exp.ttcaca+cps.obs.exp.ttcacc+cps.obs.exp.ttcacg+cps.obs.exp.ttcact+cps.obs.exp.ttcaga+cps.obs.exp.ttcagc+cps.obs.exp.ttcagg+cps.obs.exp.ttcagt+cps.obs.exp.ttcata+cps.obs.exp.ttcatc+cps.obs.exp.ttcatg+cps.obs.exp.ttcatt+cps.obs.exp.ttccaa+cps.obs.exp.ttccac+cps.obs.exp.ttccag+cps.obs.exp.ttccat+cps.obs.exp.ttccca+cps.obs.exp.ttcccc+cps.obs.exp.ttcccg+cps.obs.exp.ttccct+cps.obs.exp.ttccga+cps.obs.exp.ttccgc+cps.obs.exp.ttccgg+cps.obs.exp.ttccgt+cps.obs.exp.ttccta+cps.obs.exp.ttcctc+cps.obs.exp.ttcctg+cps.obs.exp.ttcctt+cps.obs.exp.ttcgaa+cps.obs.exp.ttcgac+cps.obs.exp.ttcgag+cps.obs.exp.ttcgat+cps.obs.exp.ttcgca+cps.obs.exp.ttcgcc+cps.obs.exp.ttcgcg+cps.obs.exp.ttcgct+cps.obs.exp.ttcgga+cps.obs.exp.ttcggc+cps.obs.exp.ttcggg+cps.obs.exp.ttcggt+cps.obs.exp.ttcgta+cps.obs.exp.ttcgtc+cps.obs.exp.ttcgtg+cps.obs.exp.ttcgtt+cps.obs.exp.ttctac+cps.obs.exp.ttctat+cps.obs.exp.ttctca+cps.obs.exp.ttctcc+cps.obs.exp.ttctcg+cps.obs.exp.ttctct+cps.obs.exp.ttctgc+cps.obs.exp.ttctgg+cps.obs.exp.ttctgt+cps.obs.exp.ttctta+cps.obs.exp.ttcttc+cps.obs.exp.ttcttg+cps.obs.exp.ttcttt+cps.obs.exp.ttgaaa+cps.obs.exp.ttgaac+cps.obs.exp.ttgaag+cps.obs.exp.ttgaat+cps.obs.exp.ttgaca+cps.obs.exp.ttgacc+cps.obs.exp.ttgacg+cps.obs.exp.ttgact+cps.obs.exp.ttgaga+cps.obs.exp.ttgagc+cps.obs.exp.ttgagg+cps.obs.exp.ttgagt+cps.obs.exp.ttgata+cps.obs.exp.ttgatc+cps.obs.exp.ttgatg+cps.obs.exp.ttgatt+cps.obs.exp.ttgcaa+cps.obs.exp.ttgcac+cps.obs.exp.ttgcag+cps.obs.exp.ttgcat+cps.obs.exp.ttgcca+cps.obs.exp.ttgccc+cps.obs.exp.ttgccg+cps.obs.exp.ttgcct+cps.obs.exp.ttgcga+cps.obs.exp.ttgcgc+cps.obs.exp.ttgcgg+cps.obs.exp.ttgcgt+cps.obs.exp.ttgcta+cps.obs.exp.ttgctc+cps.obs.exp.ttgctg+cps.obs.exp.ttgctt+cps.obs.exp.ttggaa+cps.obs.exp.ttggac+cps.obs.exp.ttggag+cps.obs.exp.ttggat+cps.obs.exp.ttggca+cps.obs.exp.ttggcc+cps.obs.exp.ttggcg+cps.obs.exp.ttggct+cps.obs.exp.ttggga+cps.obs.exp.ttgggc+cps.obs.exp.ttgggg+cps.obs.exp.ttgggt+cps.obs.exp.ttggta+cps.obs.exp.ttggtc+cps.obs.exp.ttggtg+cps.obs.exp.ttggtt+cps.obs.exp.ttgtac+cps.obs.exp.ttgtat+cps.obs.exp.ttgtca+cps.obs.exp.ttgtcc+cps.obs.exp.ttgtcg+cps.obs.exp.ttgtct+cps.obs.exp.ttgtgc+cps.obs.exp.ttgtgg+cps.obs.exp.ttgtgt+cps.obs.exp.ttgtta+cps.obs.exp.ttgttc+cps.obs.exp.ttgttg+cps.obs.exp.ttgttt+cps.obs.exp.tttaaa+cps.obs.exp.tttaac+cps.obs.exp.tttaag+cps.obs.exp.tttaat+cps.obs.exp.tttaca+cps.obs.exp.tttacc+cps.obs.exp.tttacg+cps.obs.exp.tttact+cps.obs.exp.tttaga+cps.obs.exp.tttagc+cps.obs.exp.tttagg+cps.obs.exp.tttagt+cps.obs.exp.tttata+cps.obs.exp.tttatc+cps.obs.exp.tttatg+cps.obs.exp.tttatt+cps.obs.exp.tttcaa+cps.obs.exp.tttcac+cps.obs.exp.tttcag+cps.obs.exp.tttcat+cps.obs.exp.tttcca+cps.obs.exp.tttccc+cps.obs.exp.tttccg+cps.obs.exp.tttcct+cps.obs.exp.tttcga+cps.obs.exp.tttcgc+cps.obs.exp.tttcgg+cps.obs.exp.tttcgt+cps.obs.exp.tttcta+cps.obs.exp.tttctc+cps.obs.exp.tttctg+cps.obs.exp.tttctt+cps.obs.exp.tttgaa+cps.obs.exp.tttgac+cps.obs.exp.tttgag+cps.obs.exp.tttgat+cps.obs.exp.tttgca+cps.obs.exp.tttgcc+cps.obs.exp.tttgcg+cps.obs.exp.tttgct+cps.obs.exp.tttgga+cps.obs.exp.tttggc+cps.obs.exp.tttggg+cps.obs.exp.tttggt+cps.obs.exp.tttgta+cps.obs.exp.tttgtc+cps.obs.exp.tttgtg+cps.obs.exp.tttgtt+cps.obs.exp.ttttac+cps.obs.exp.ttttat+cps.obs.exp.ttttca+cps.obs.exp.ttttcc+cps.obs.exp.ttttcg+cps.obs.exp.ttttct+cps.obs.exp.ttttgc+cps.obs.exp.ttttgg+cps.obs.exp.ttttgt+cps.obs.exp.ttttta+cps.obs.exp.tttttc+cps.obs.exp.tttttg+cps.obs.exp.tttttt

#Read CPS table without Inf
#cps.df.mod<-read.table(file="/home2/dmacedod/test/dn3_6deg_original_mod.txt", header=TRUE)

#Generate sum of all codon pair per sequence
cps.df.mod.sum<-apply(cps.df.final, 1, function(x) sum(x))

#Calculate CPB
cpb.n3<-cps.df.mod.sum/3693

#Write as a table (CPB)
write.table(cpb.n3, file="/home2/dmacedod/test/dn3_6deg_cpb_result.txt")

#Graph 
#png("CPB/n3_cpb.png")
#hist(cpb.n3, main="CPB (n3)", xlab=NULL)
#abline(v=cpb.n3[1], col=2,lty=2)
#text (0.65, 250, "WT", cex=0.5, col=2)
#dev.off()


               

         
               
               
