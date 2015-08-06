#Read fasta file
dn3.6deg <- read.fasta(file = "3dn_6deg/P1_3dn_6deg_with_original_seq.fas")                     
n3 <- read.fasta(file = "3N/P1_3n_with_original_seq.fas")                     
dn23 <- read.fasta(file = "dN23/P1_dn23_with_original_seq.fas")                     
dn31 <- read.fasta(file = "dN31/P1_dn31_with_original_seq.fas")                     

#Count dinucleotide frequency in the sequence (Package: plyr)
#lapply(n3, function(x) count(x,2))
#Generate a Table with dinucleotide frequency 
dn3.6deg.dinuc<-ldply(dn3.6deg, function(x) count(x,2)) 
n3.dinuc<-ldply(n3, function(x) count(x,2)) 
dn23.dinuc<-ldply(dn23, function(x) count(x,2)) 
dn31.dinuc<-ldply(dn31, function(x) count(x,2)) 

#ApA in the sequence
dn3.6deg.ApA<-dn3.6deg.dinuc$aa
n3.ApA<-n3.dinuc$aa
dn23.ApA<-dn23.dinuc$aa
dn31.ApA<-dn31.dinuc$aa

#ApC in the sequence
dn3.6deg.ApC<-dn3.6deg.dinuc$ac
n3.ApC<-n3.dinuc$ac
dn23.ApC<-dn23.dinuc$ac
dn31.ApC<-dn31.dinuc$ac

#ApG in the sequence
dn3.6deg.ApG<-dn3.6deg.dinuc$ag
n3.ApG<-n3.dinuc$ag
dn23.ApG<-dn23.dinuc$ag
dn31.ApG<-dn31.dinuc$ag

#ApT in the sequence
dn3.6deg.ApT<-dn3.6deg.dinuc$at
n3.ApT<-n3.dinuc$at
dn23.ApT<-dn23.dinuc$at
dn31.ApT<-dn31.dinuc$at

#CpA in the sequence
dn3.6deg.CpA<-dn3.6deg.dinuc$ca
n3.CpA<-n3.dinuc$ca
dn23.CpA<-dn23.dinuc$ca
dn31.CpA<-dn31.dinuc$ca

#CpC in the sequence
dn3.6deg.CpC<-dn3.6deg.dinuc$cc
n3.CpC<-n3.dinuc$cc
dn23.CpC<-dn23.dinuc$cc
dn31.CpC<-dn31.dinuc$cc

#CpG in the sequence
dn3.6deg.CpG<-dn3.6deg.dinuc$cg
n3.CpG<-n3.dinuc$cg
dn23.CpG<-dn23.dinuc$cg
dn31.CpG<-dn31.dinuc$cg

#CpT in the sequence
dn3.6deg.CpT<-dn3.6deg.dinuc$ct
n3.CpT<-n3.dinuc$ct
dn23.CpT<-dn23.dinuc$ct
dn31.CpT<-dn31.dinuc$ct


#GpA in the sequence
dn3.6deg.GpA<-dn3.6deg.dinuc$ga
n3.GpA<-n3.dinuc$ga
dn23.GpA<-dn23.dinuc$ga
dn31.GpA<-dn31.dinuc$ga

#GpC in the sequence
dn3.6deg.GpC<-dn3.6deg.dinuc$gc
n3.GpC<-n3.dinuc$gc
dn23.GpC<-dn23.dinuc$gc
dn31.GpC<-dn31.dinuc$gc

#GpG in the sequence
dn3.6deg.GpG<-dn3.6deg.dinuc$gg
n3.GpG<-n3.dinuc$gg
dn23.GpG<-dn23.dinuc$gg
dn31.GpG<-dn31.dinuc$gg

#GpT in the sequence
dn3.6deg.GpT<-dn3.6deg.dinuc$gt
n3.GpT<-n3.dinuc$gt
dn23.GpT<-dn23.dinuc$gt
dn31.GpT<-dn31.dinuc$gt


#TpA in the sequence
dn3.6deg.TpA<-dn3.6deg.dinuc$ta
n3.TpA<-n3.dinuc$ta
dn23.TpA<-dn23.dinuc$ta
dn31.TpA<-dn31.dinuc$ta

#TpC in the sequence
dn3.6deg.TpC<-dn3.6deg.dinuc$tc
n3.TpC<-n3.dinuc$tc
dn23.TpC<-dn23.dinuc$tc
dn31.TpC<-dn31.dinuc$tc

#TpG in the sequence
dn3.6deg.TpG<-dn3.6deg.dinuc$tg
n3.TpG<-n3.dinuc$tg
dn23.TpG<-dn23.dinuc$tg
dn31.TpG<-dn31.dinuc$tg

#TpT in the sequence
dn3.6deg.TpT<-dn3.6deg.dinuc$tt
n3.TpT<-n3.dinuc$tt
dn23.TpT<-dn23.dinuc$tt
dn31.TpT<-dn31.dinuc$tt


#Count dinucleotide frequency in the sequence
#lapply(n3, function(x) count(x,2))
#Nucleotide composition
dn3.6deg.comp<-ldply(dn3.6deg, function(x) count(x,1)) 
n3.comp<-ldply(n3, function(x) count(x,1)) 
dn23.comp<-ldply(dn23, function(x) count(x,1)) 
dn31.comp<-ldply(dn31, function(x) count(x,1)) 

#A composition in the sequence
dn3.6deg.A<-dn3.6deg.comp$a
n3.A<-n3.comp$a
dn23.A<-dn23.comp$a
dn31.A<-dn31.comp$a

#C composition in the sequence
dn3.6deg.C<-dn3.6deg.comp$c
n3.C<-n3.comp$c
dn23.C<-dn23.comp$c
dn31.C<-dn31.comp$c

#G composition in the sequence
dn3.6deg.G<-dn3.6deg.comp$g
n3.G<-n3.comp$g
dn23.G<-dn23.comp$g
dn31.G<-dn31.comp$g

#T composition in the sequence
dn3.6deg.T<-dn3.6deg.comp$t
n3.T<-n3.comp$t
dn23.T<-dn23.comp$t
dn31.T<-dn31.comp$t

#Ratio formula (number of Dinuc/(number of A * number of A))
#Ratio ApA
dn3.6deg.apa.ratio<-dn3.6deg.ApA/(dn3.6deg.A*dn3.6deg.A)
n3.apa.ratio<-n3.ApA/(n3.A*n3.A)
dn23.apa.ratio<-dn23.ApA/(dn23.A*dn23.A)
dn31.apa.ratio<-dn31.ApA/(dn31.A*dn31.A)

#Ratio ApC
dn3.6deg.apc.ratio<-dn3.6deg.ApC/(dn3.6deg.A*dn3.6deg.C)
n3.apc.ratio<-n3.ApC/(n3.A*n3.C)
dn23.apc.ratio<-dn23.ApC/(dn23.A*dn23.C)
dn31.apc.ratio<-dn31.ApC/(dn31.A*dn31.C)

#Ratio ApG
dn3.6deg.apg.ratio<-dn3.6deg.ApG/(dn3.6deg.A*dn3.6deg.G)
n3.apg.ratio<-n3.ApG/(n3.A*n3.G)
dn23.apg.ratio<-dn23.ApG/(dn23.A*dn23.G)
dn31.apg.ratio<-dn31.ApG/(dn31.A*dn31.G)

#Ratio ApT
dn3.6deg.apt.ratio<-dn3.6deg.ApT/(dn3.6deg.A*dn3.6deg.T)
n3.apt.ratio<-n3.ApT/(n3.A*n3.T)
dn23.apt.ratio<-dn23.ApT/(dn23.A*dn23.T)
dn31.apt.ratio<-dn31.ApT/(dn31.A*dn31.T)



#Ratio CpA
dn3.6deg.cpa.ratio<-dn3.6deg.CpA/(dn3.6deg.C*dn3.6deg.A)
n3.cpa.ratio<-n3.CpA/(n3.C*n3.A)
dn23.cpa.ratio<-dn23.CpA/(dn23.C*dn23.A)
dn31.cpa.ratio<-dn31.CpA/(dn31.C*dn31.A)

#Ratio CpC
dn3.6deg.cpc.ratio<-dn3.6deg.CpC/(dn3.6deg.C*dn3.6deg.C)
n3.cpc.ratio<-n3.CpC/(n3.C*n3.C)
dn23.cpc.ratio<-dn23.CpC/(dn23.C*dn23.C)
dn31.cpc.ratio<-dn31.CpC/(dn31.C*dn31.C)

#Ratio CpG
dn3.6deg.cpg.ratio<-dn3.6deg.CpG/(dn3.6deg.C*dn3.6deg.G)
n3.cpg.ratio<-n3.CpG/(n3.C*n3.G)
dn23.cpg.ratio<-dn23.CpG/(dn23.C*dn23.G)
dn31.cpg.ratio<-dn31.CpG/(dn31.C*dn31.G)

#Ratio CpT
dn3.6deg.cpt.ratio<-dn3.6deg.CpT/(dn3.6deg.C*dn3.6deg.T)
n3.cpt.ratio<-n3.CpT/(n3.C*n3.T)
dn23.cpt.ratio<-dn23.CpT/(dn23.C*dn23.T)
dn31.cpt.ratio<-dn31.CpT/(dn31.C*dn31.T)


#Ratio GpA
dn3.6deg.gpa.ratio<-dn3.6deg.GpA/(dn3.6deg.G*dn3.6deg.A)
n3.gpa.ratio<-n3.GpA/(n3.G*n3.A)
dn23.gpa.ratio<-dn23.GpA/(dn23.G*dn23.A)
dn31.gpa.ratio<-dn31.GpA/(dn31.G*dn31.A)

#Ratio GpC
dn3.6deg.gpc.ratio<-dn3.6deg.GpC/(dn3.6deg.G*dn3.6deg.C)
n3.gpc.ratio<-n3.GpC/(n3.G*n3.C)
dn23.gpc.ratio<-dn23.GpC/(dn23.G*dn23.C)
dn31.gpc.ratio<-dn31.GpC/(dn31.G*dn31.C)

#Ratio GpG
dn3.6deg.gpg.ratio<-dn3.6deg.GpG/(dn3.6deg.G*dn3.6deg.G)
n3.gpg.ratio<-n3.GpG/(n3.G*n3.G)
dn23.gpg.ratio<-dn23.GpG/(dn23.G*dn23.G)
dn31.gpg.ratio<-dn31.GpG/(dn31.G*dn31.G)

#Ratio GpT
dn3.6deg.gpt.ratio<-dn3.6deg.GpT/(dn3.6deg.G*dn3.6deg.T)
n3.gpt.ratio<-n3.GpT/(n3.G*n3.T)
dn23.gpt.ratio<-dn23.GpT/(dn23.G*dn23.T)
dn31.gpt.ratio<-dn31.GpT/(dn31.G*dn31.T)



#Ratio TpA
dn3.6deg.tpa.ratio<-dn3.6deg.TpA/(dn3.6deg.T*dn3.6deg.A)
n3.tpa.ratio<-n3.TpA/(n3.T*n3.A)
dn23.tpa.ratio<-dn23.TpA/(dn23.T*dn23.A)
dn31.tpa.ratio<-dn31.TpA/(dn31.T*dn31.A)

#Ratio TpC
dn3.6deg.tpc.ratio<-dn3.6deg.TpC/(dn3.6deg.T*dn3.6deg.C)
n3.tpc.ratio<-n3.TpC/(n3.T*n3.C)
dn23.tpc.ratio<-dn23.TpC/(dn23.T*dn23.C)
dn31.tpc.ratio<-dn31.TpC/(dn31.T*dn31.C)

#Ratio TpG
dn3.6deg.tpg.ratio<-dn3.6deg.TpG/(dn3.6deg.T*dn3.6deg.G)
n3.tpg.ratio<-n3.TpG/(n3.T*n3.G)
dn23.tpg.ratio<-dn23.TpG/(dn23.T*dn23.G)
dn31.tpg.ratio<-dn31.TpG/(dn31.T*dn31.G)

#Ratio TpT
dn3.6deg.tpt.ratio<-dn3.6deg.TpT/(dn3.6deg.T*dn3.6deg.T)
n3.tpt.ratio<-n3.TpT/(n3.T*n3.T)
dn23.tpt.ratio<-dn23.TpT/(dn23.T*dn23.T)
dn31.tpt.ratio<-dn31.TpT/(dn31.T*dn31.T)

#Write table
#ApA ratio 
write.table(dn3.6deg.apa.ratio, file="Dinucleotide table/dn3_6deg_apa_ratio_new.txt", sep="\t")
write.table(n3.apa.ratio, file="Dinucleotide table/n3_apa_ratio_new.txt", sep="\t")
write.table(dn23.apa.ratio, file="Dinucleotide table/dn23_apa_ratio_new.txt", sep="\t")
write.table(dn31.apa.ratio, file="Dinucleotide table/dn31_apa_ratio_new.txt", sep="\t")

#ApC ratio 
write.table(dn3.6deg.apc.ratio, file="Dinucleotide table/dn3_6deg_apc_ratio_new.txt", sep="\t")
write.table(n3.apc.ratio, file="Dinucleotide table/n3_apc_ratio_new.txt", sep="\t")
write.table(dn23.apc.ratio, file="Dinucleotide table/dn23_apc_ratio_new.txt", sep="\t")
write.table(dn31.apc.ratio, file="Dinucleotide table/dn31_apc_ratio_new.txt", sep="\t")

#ApG ratio 
write.table(dn3.6deg.apg.ratio, file="Dinucleotide table/dn3_6deg_apg_ratio_new.txt", sep="\t")
write.table(n3.apg.ratio, file="Dinucleotide table/n3_apg_ratio_new.txt", sep="\t")
write.table(dn23.apg.ratio, file="Dinucleotide table/dn23_apg_ratio_new.txt", sep="\t")
write.table(dn31.apg.ratio, file="Dinucleotide table/dn31_apg_ratio_new.txt", sep="\t")

#ApT ratio 
write.table(dn3.6deg.apt.ratio, file="Dinucleotide table/dn3_6deg_apt_ratio_new.txt", sep="\t")
write.table(n3.apt.ratio, file="Dinucleotide table/n3_apt_ratio_new.txt", sep="\t")
write.table(dn23.apt.ratio, file="Dinucleotide table/dn23_apt_ratio_new.txt", sep="\t")
write.table(dn31.apt.ratio, file="Dinucleotide table/dn31_apt_ratio_new.txt", sep="\t")



#CpA ratio 
write.table(dn3.6deg.cpa.ratio, file="Dinucleotide table/dn3_6deg_cpa_ratio_new.txt", sep="\t")
write.table(n3.cpa.ratio, file="Dinucleotide table/n3_cpa_ratio_new.txt", sep="\t")
write.table(dn23.cpa.ratio, file="Dinucleotide table/dn23_cpa_ratio_new.txt", sep="\t")
write.table(dn31.cpa.ratio, file="Dinucleotide table/dn31_cpa_ratio_new.txt", sep="\t")

#CpC ratio 
write.table(dn3.6deg.cpc.ratio, file="Dinucleotide table/dn3_6deg_cpc_ratio_new.txt", sep="\t")
write.table(n3.cpc.ratio, file="Dinucleotide table/n3_cpc_ratio_new.txt", sep="\t")
write.table(dn23.cpc.ratio, file="Dinucleotide table/dn23_cpc_ratio_new.txt", sep="\t")
write.table(dn31.cpc.ratio, file="Dinucleotide table/dn31_cpc_ratio_new.txt", sep="\t")

#CpG ratio 
write.table(dn3.6deg.cpg.ratio, file="Dinucleotide table/dn3_6deg_cpg_ratio_new.txt", sep="\t")
write.table(n3.cpg.ratio, file="Dinucleotide table/n3_cpg_ratio_new.txt", sep="\t")
write.table(dn23.cpg.ratio, file="Dinucleotide table/dn23_cpg_ratio_new.txt", sep="\t")
write.table(dn31.cpg.ratio, file="Dinucleotide table/dn31_cpg_ratio_new.txt", sep="\t")

#CpT ratio 
write.table(dn3.6deg.cpt.ratio, file="Dinucleotide table/dn3_6deg_cpt_ratio_new.txt", sep="\t")
write.table(n3.cpt.ratio, file="Dinucleotide table/n3_cpt_ratio_new.txt", sep="\t")
write.table(dn23.cpt.ratio, file="Dinucleotide table/dn23_cpt_ratio_new.txt", sep="\t")
write.table(dn31.cpt.ratio, file="Dinucleotide table/dn31_cpt_ratio_new.txt", sep="\t")



#GpA ratio 
write.table(dn3.6deg.gpa.ratio, file="Dinucleotide table/dn3_6deg_gpa_ratio_new.txt", sep="\t")
write.table(n3.gpa.ratio, file="Dinucleotide table/n3_gpa_ratio_new.txt", sep="\t")
write.table(dn23.gpa.ratio, file="Dinucleotide table/dn23_gpa_ratio_new.txt", sep="\t")
write.table(dn31.gpa.ratio, file="Dinucleotide table/dn31_gpa_ratio_new.txt", sep="\t")

#GpC ratio 
write.table(dn3.6deg.gpc.ratio, file="Dinucleotide table/dn3_6deg_gpc_ratio_new.txt", sep="\t")
write.table(n3.gpc.ratio, file="Dinucleotide table/n3_gpc_ratio_new.txt", sep="\t")
write.table(dn23.gpc.ratio, file="Dinucleotide table/dn23_gpc_ratio_new.txt", sep="\t")
write.table(dn31.gpc.ratio, file="Dinucleotide table/dn31_gpc_ratio_new.txt", sep="\t")

#GpG ratio 
write.table(dn3.6deg.gpg.ratio, file="Dinucleotide table/dn3_6deg_gpg_ratio_new.txt", sep="\t")
write.table(n3.gpg.ratio, file="Dinucleotide table/n3_gpg_ratio_new.txt", sep="\t")
write.table(dn23.gpg.ratio, file="Dinucleotide table/dn23_gpg_ratio_new.txt", sep="\t")
write.table(dn31.gpg.ratio, file="Dinucleotide table/dn31_gpg_ratio_new.txt", sep="\t")

#GpT ratio 
write.table(dn3.6deg.gpt.ratio, file="Dinucleotide table/dn3_6deg_gpt_ratio_new.txt", sep="\t")
write.table(n3.gpt.ratio, file="Dinucleotide table/n3_gpt_ratio_new.txt", sep="\t")
write.table(dn23.gpt.ratio, file="Dinucleotide table/dn23_gpt_ratio_new.txt", sep="\t")
write.table(dn31.gpt.ratio, file="Dinucleotide table/dn31_gpt_ratio_new.txt", sep="\t")



#TpA ratio 
write.table(dn3.6deg.tpa.ratio, file="Dinucleotide table/dn3_6deg_tpa_ratio_new.txt", sep="\t")
write.table(n3.tpa.ratio, file="Dinucleotide table/n3_tpa_ratio_new.txt", sep="\t")
write.table(dn23.tpa.ratio, file="Dinucleotide table/dn23_tpa_ratio_new.txt", sep="\t")
write.table(dn31.tpa.ratio, file="Dinucleotide table/dn31_tpa_ratio_new.txt", sep="\t")

#TpC ratio 
write.table(dn3.6deg.tpc.ratio, file="Dinucleotide table/dn3_6deg_tpc_ratio_new.txt", sep="\t")
write.table(n3.tpc.ratio, file="Dinucleotide table/n3_tpc_ratio_new.txt", sep="\t")
write.table(dn23.tpc.ratio, file="Dinucleotide table/dn23_tpc_ratio_new.txt", sep="\t")
write.table(dn31.tpc.ratio, file="Dinucleotide table/dn31_tpc_ratio_new.txt", sep="\t")

#TpG ratio 
write.table(dn3.6deg.tpg.ratio, file="Dinucleotide table/dn3_6deg_tpg_ratio_new.txt", sep="\t")
write.table(n3.tpg.ratio, file="Dinucleotide table/n3_tpg_ratio_new.txt", sep="\t")
write.table(dn23.tpg.ratio, file="Dinucleotide table/dn23_tpg_ratio_new.txt", sep="\t")
write.table(dn31.tpg.ratio, file="Dinucleotide table/dn31_tpg_ratio_new.txt", sep="\t")

#TpT ratio 
write.table(dn3.6deg.tpt.ratio, file="Dinucleotide table/dn3_6deg_tpt_ratio_new.txt", sep="\t")
write.table(n3.tpt.ratio, file="Dinucleotide table/n3_tpt_ratio_new.txt", sep="\t")
write.table(dn23.tpt.ratio, file="Dinucleotide table/dn23_tpt_ratio_new.txt", sep="\t")
write.table(dn31.tpt.ratio, file="Dinucleotide table/dn31_tpt_ratio_new.txt", sep="\t")


