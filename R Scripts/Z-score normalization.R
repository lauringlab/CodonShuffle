#Nucleotide distance 
dn3.6deg.nuc.dist<-read.table(file="3dn_6deg/3dn_6deg_nucleotide difference_wild_type.txt", header=TRUE)
n3.nuc.dist<-read.table(file="3N/3n_nucleotide difference_wild_type.txt", header=TRUE)
dn23.nuc.dist<-read.table(file="dn23/dn23_nucleotide difference_wild_type.txt", header=TRUE)
dn31.nuc.dist<-read.table(file="dn31/dn31_nucleotide difference_wild_type.txt", header=TRUE)

#Mfold
dn3.6deg.fold <- read.table(file = "Unafold//3dn_6deg_RNA.fas.dG")                     
n3.fold <- read.table(file = "Unafold//3n_RNA.fas.dG")                     
dn23.fold <- read.table(file = "Unafold/dn23_RNA.fas.dG")                     
dn31.fold <- read.table(file = "Unafold/dn31_RNA.fas.dG")                     

#ENC
dn3.6deg.enc <- read.table(file = "ENC/dn3_6deg_enc.txt", header=TRUE)                     
n3.enc <- read.table(file = "ENC/n3_enc.txt", header=TRUE)
dn23.enc <- read.table(file = "ENC/dn23_enc.txt", header=TRUE)
dn31.enc <- read.table(file = "ENC/dn31_enc.txt", header=TRUE)


#CAI
dn3.6deg.cai <- read.table(file = "CAI/dn3_6deg.cai", header=TRUE)                     
n3.cai <- read.table(file = "CAI/n3.cai", header=TRUE)                     
dn23.cai <- read.table(file = "CAI/dn23.cai", header=TRUE)                     
dn31.cai <- read.table(file = "CAI/dn31.cai", header=TRUE)                     


#CPB
dn3.6deg.cpb <- read.table(file = "CPB/dn3_6deg_cpb.txt", header=TRUE)
n3.cpb <- read.table(file = "CPB/n3_cpb.txt", header=TRUE)
dn23.cpb <- read.table(file = "CPB/dn23_cpb.txt", header=TRUE)
dn31.cpb <- read.table(file = "CPB/dn31_cpb.txt", header=TRUE)

#Dinuc

#Trinuc



#Z-score
#Distance
dn3.6deg.nuc.dist.z<-(dn3.6deg.nuc.dist$Y-mean(dn3.6deg.nuc.dist$Y))/sd(dn3.6deg.nuc.dist$Y)
n3.nuc.dist.z<-(n3.nuc.dist$Y-mean(n3.nuc.dist$Y))/sd(n3.nuc.dist$Y)
dn23.nuc.dist.z<-(dn23.nuc.dist$Y-mean(dn23.nuc.dist$Y))/sd(dn23.nuc.dist$Y)
dn31.nuc.dist.z<-(dn31.nuc.dist$Y-mean(dn31.nuc.dist$Y))/sd(dn31.nuc.dist$Y)

#Mfold
dn3.6deg.fold.z<-(dn3.6deg.fold$V2-mean(dn3.6deg.fold$V2))/sd(dn3.6deg.fold$V2)
n3.fold.z<-(n3.fold$V2-mean(n3.fold$V2))/sd(n3.fold$V2)
dn23.fold.z<-(dn23.fold$V2-mean(dn23.fold$V2))/sd(dn23.fold$V2)
dn31.fold.z<-(dn31.fold$V2-mean(dn31.fold$V2))/sd(dn31.fold$V2)

#Cai
dn3.6deg.cai.z<-(dn3.6deg.cai$CAI-mean(dn3.6deg.cai$CAI))/sd(dn3.6deg.cai$CAI)
n3.cai.z<-(n3.cai$CAI-mean(n3.cai$CAI))/sd(n3.cai$CAI)
dn23.cai.z<-(dn23.cai$CAI-mean(dn23.cai$CAI))/sd(dn23.cai$CAI)
dn31.cai.z<-(dn31.cai$CAI-mean(dn31.cai$CAI))/sd(dn31.cai$CAI)

#ENC
dn3.6deg.enc.z<-(dn3.6deg.enc$Nc-mean(dn3.6deg.enc$Nc))/sd(dn3.6deg.enc$Nc)
n3.enc.z<-(n3.enc$Nc-mean(n3.enc$Nc))/sd(n3.enc$Nc)
dn23.enc.z<-(dn23.enc$Nc-mean(dn23.enc$Nc))/sd(dn23.enc$Nc)
dn31.enc.z<-(dn31.enc$Nc-mean(dn31.enc$Nc))/sd(dn31.enc$Nc)

#CPB
dn3.6deg.cpb.z<-(dn3.6deg.cpb-mean(dn3.6deg.cpb$x))/sd(dn3.6deg.cpb$x)
n3.cpb.z<-(n3.cpb-mean(n3.cpb$x))/sd(n3.cpb$x)
dn23.cpb.z<-(dn23.cpb-mean(dn23.cpb$x))/sd(dn23.cpb$x)
dn31.cpb.z<-(dn31.cpb-mean(dn31.cpb$x))/sd(dn31.cpb$x)

#Dinuc

#Trinuc





#Graph
png("Scripts_distance_Z.png")
par(mfrow=c(2,2))
hist(dn3.6deg.nuc.dist.z, main = "dn3 6deg nucleotide distance", xlab=NULL)
hist(n3.nuc.dist.z, main = "n3 nucleotide distance", xlab=NULL)
hist(dn23.nuc.dist.z, main = "dn23 nucleotide distance", xlab=NULL)
hist(dn31.nuc.dist.z, main = "dn31 nucleotide distance", xlab=NULL)
text (0.65, 250, "WT", cex=0.5, col=2)
dev.off()

png("Scripts_mfold_Z_new.png", width = 900, height = 900)
par(mfrow=c(2,2), mgp=c(2.3, 0.5, 0))
hist(dn3.6deg.fold.z, main = "Mfold (dn3 6deg)", xlab=NULL, cex.axis=2, cex.lab=2, cex.main=3)
abline(v=dn3.6deg.fold.z[1], col=2,lty=2, lwd=4)
hist(n3.fold.z, main = "Mfold (n3)", xlab=NULL, cex.axis=2, cex.lab=2, cex.main=3)
abline(v=n3.fold.z[1], col=2,lty=2, lwd=4)
hist(dn23.fold.z, main = "Mfold (dn23)", xlab=NULL, cex.axis=2, cex.lab=2, cex.main=3)
abline(v=dn23.fold.z[1], col=2,lty=2, lwd=4)
hist(dn31.fold.z, main = "Mfold (dn31)", xlab=NULL, cex.axis=2, cex.lab=2, cex.main=3)
abline(v=dn31.fold.z[1], col=2,lty=2, lwd=4)
text (0.65, 250, "WT", cex=0.5, col=2)
dev.off()

png("Scripts_cai_Z.png")
par(mfrow=c(2,2))
hist(dn3.6deg.cai.z, main = "dn3 6deg nucleotide cai", xlab=NULL)
abline(v=dn3.6deg.cai.z[1], col=2,lty=2)
hist(n3.cai.z, main = "n3 nucleotide cai", xlab=NULL)
abline(v=n3.cai.z[1], col=2,lty=2)
hist(dn23.cai.z, main = "dn23 nucleotide cai", xlab=NULL)
abline(v=dn23.cai.z[1], col=2,lty=2)
hist(dn31.cai.z, main = "dn31 nucleotide cai", xlab=NULL)
abline(v=dn31.cai.z[1], col=2,lty=2)
text (0.65, 250, "WT", cex=0.5, col=2)
dev.off()

png("Scripts_enc_Z.png")
par(mfrow=c(2,2))
hist(dn3.6deg.enc.z, main = "dn3 6deg nucleotide enc", xlab=NULL)
abline(v=dn3.6deg.enc.z[1], col=2,lty=2)
hist(n3.enc.z, main = "n3 nucleotide enc", xlab=NULL)
abline(v=n3.enc.z[1], col=2,lty=2)
hist(dn23.enc.z, main = "dn23 nucleotide enc", xlab=NULL)
abline(v=dn23.enc.z[1], col=2,lty=2)
hist(dn31.enc.z, main = "dn31 nucleotide enc", xlab=NULL)
abline(v=dn31.enc.z[1], col=2,lty=2)
text (0.65, 250, "WT", cex=0.5, col=2)
dev.off()


png("Scripts_cpb_Z.png")
par(mfrow=c(2,2))
hist(dn3.6deg.cpb.z$x, main = "dn3 6deg nucleotide cpb", xlab=NULL)
abline(v=dn3.6deg.cpb.z[1,], col=2,lty=2)
hist(n3.cpb.z$x, main = "n3 nucleotide cpb", xlab=NULL)
abline(v=n3.cpb.z[1,], col=2,lty=2)
hist(dn23.cpb.z$x, main = "dn23 nucleotide cpb", xlab=NULL)
abline(v=dn23.cpb.z[1,], col=2,lty=2)
hist(dn31.cpb.z$x, main = "dn31 nucleotide cpb", xlab=NULL)
abline(v=dn31.cpb.z[1,], col=2,lty=2)
text (0.65, 250, "WT", cex=0.5, col=2)
dev.off()

#Write table
#Distance
write.table(dn3.6deg.nuc.dist.z, file="Z_norm/dn3_6deg_nuc_dist_z.txt", sep="\t")
write.table(n3.nuc.dist.z, file="Z_norm/n3_nuc_dist_z.txt", sep="\t")
write.table(dn23.nuc.dist.z, file="Z_norm/dn23_nuc_dist_z.txt", sep="\t")
write.table(dn31.nuc.dist.z, file="Z_norm/dn31_nuc_dist_z.txt", sep="\t")

#Mfold
write.table(dn3.6deg.fold.z, file="Z_norm/dn3_6deg_fold_z.txt", sep="\t")
write.table(n3.fold.z, file="Z_norm/n3_fold_z.txt", sep="\t")
write.table(dn23.fold.z, file="Z_norm/dn23_fold_z.txt", sep="\t")
write.table(dn31.fold.z, file="Z_norm/dn31_fold_z.txt", sep="\t")

#Cai
write.table(dn3.6deg.cai.z, file="Z_norm/dn3_6deg_cai_z.txt", sep="\t")
write.table(n3.cai.z, file="Z_norm/n3_cai_z.txt", sep="\t")
write.table(dn23.cai.z, file="Z_norm/dn23_cai_z.txt", sep="\t")
write.table(dn31.cai.z, file="Z_norm/dn31_cai_z.txt", sep="\t")


#ENC
write.table(dn3.6deg.enc.z, file="Z_norm/dn3_6deg_enc_z.txt", sep="\t")
write.table(n3.enc.z, file="Z_norm/n3_enc_z.txt", sep="\t")
write.table(dn23.enc.z, file="Z_norm/dn23_enc_z.txt", sep="\t")
write.table(dn31.enc.z, file="Z_norm/dn31_enc_z.txt", sep="\t")


#CPB
write.table(dn3.6deg.cpb.z, file="Z_norm/dn3_6deg_cpb_z.txt", sep="\t")
write.table(n3.cpb.z, file="Z_norm/n3_cpb_z.txt", sep="\t")
write.table(dn23.cpb.z, file="Z_norm/dn23_cpb_z.txt", sep="\t")
write.table(dn31.cpb.z, file="Z_norm/dn31_cpb_z.txt", sep="\t")
