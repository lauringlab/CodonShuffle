#Read Enc table 
dn3.6deg.enc<-read.table(file="ENC/dn3_6deg_enc.txt", header=TRUE)
n3.enc<-read.table(file="ENC/n3_enc.txt", header=TRUE)
dn23.enc<-read.table(file="ENC/dn23_enc.txt", header=TRUE)
dn31.enc<-read.table(file="ENC/dn31_enc.txt", header=TRUE)


#Graph each script individually
#dn3 6deg
png("ENC/3dn_6deg_enc.png")
hist(dn3.6deg.enc$Nc, main="Effective number of codon random sequence (dn3 6deg)", xlab=NULL)
abline(v=dn3.6deg.enc[1,2], col=2,lty=2, lwd=4)
text (58.8, 180, "WT", cex=1.5, col=2)
dev.off()

#n3 
png("ENC/3n_enc.png")
hist(n3.enc$Nc, main="Effective number of codon random sequence (n3)", xlab=NULL)
abline(v=n3.enc[1,2], col=2,lty=2, lwd=4)
text (58.8, 180, "WT", cex=1.5, col=2)
dev.off()

#dn23 
png("ENC/dn23_enc.png")
hist(dn23.enc$Nc, main="Effective number of codon random sequence (dn23)", xlab=NULL)
abline(v=dn23.enc[1,2], col=2,lty=2, lwd=4)
text (58.8, 180, "WT", cex=1.5, col=2)
dev.off()

#dn31 
png("ENC/dn31_enc.png")
hist(dn31.enc$Nc, main="Effective number of codon random sequence (dn31)", xlab=NULL)
abline(v=dn31.enc[1,2], col=2,lty=2, lwd=4)
text (58.8, 180, "WT", cex=1.5, col=2)
dev.off()

