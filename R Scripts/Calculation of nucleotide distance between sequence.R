#Nucleotide distance 

#Read fasta sequence
dn3.6deg.seq<-read.dna(file="3dn_6deg/P1_3dn_6deg_with_original_seq.fas", format="fasta")
n3.seq<-read.dna(file="3N/P1_3n_with_original_seq.fas", format="fasta")
dn23.seq<-read.dna(file="dN23/P1_dn23_with_original_seq.fas", format="fasta")
dn31.seq<-read.dna(file="dN31/P1_dn31_with_original_seq.fas", format="fasta")


#Compute the distance between the sequences
dn3.6deg.nuc.dist<-dist.dna(dn3.6deg, model="N", as.matrix=TRUE)
n3.nuc.dist<-dist.dna(n3, model="N", as.matrix=TRUE)
dn23.nuc.dist<-dist.dna(dn23, model="N", as.matrix=TRUE)
dn31.nuc.dist<-dist.dna(dn31, model="N", as.matrix=TRUE)

#Save table
write.table(dn3.6deg.nuc.dist, file="Nuc dist/dn3_6deg_nuc_dist.txt", row.names=TRUE, col.names=TRUE)
write.table(n3.nuc.dist, file="Nuc dist/n3_nuc_dist.txt", row.names=TRUE, col.names=TRUE)
write.table(dn23.nuc.dist, file="Nuc dist/dn23_nuc_dist.txt", row.names=TRUE, col.names=TRUE)
write.table(dn31.nuc.dist, file="Nuc dist/dn31_nuc_dist.txt", row.names=TRUE, col.names=TRUE)

#Showing distance between WT and replicate sequences
dn3.6deg.nuc.dist[1,]
n3.nuc.dist[1,]
dn23.nuc.dist[1,]
dn31.nuc.dist[1,]







