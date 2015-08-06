#Read Cai table 
dn3.6deg.cai<-read.table(file="CAI/dn3_6deg.cai", header=TRUE)
n3.cai<-read.table(file="CAI/n3.cai", header=TRUE)
dn23.cai<-read.table(file="CAI/dn23.cai", header=TRUE)
dn31.cai<-read.table(file="CAI/dn31.cai", header=TRUE)

#Remove first and third column of the table
dn3.6deg.cai<-dn3.6deg.cai[,-1]
dn3.6deg.cai<-dn3.6deg.cai[,-2]
n3.cai<-n3.cai[,-1]
n3.cai<-n3.cai[,-2]
dn23.cai<-dn23.cai[,-1]
dn23.cai<-dn23.cai[,-2]
dn31.cai<-dn31.cai[,-1]
dn31.cai<-dn31.cai[,-2]

#Graph each script individually
#dn3 6deg 
png("cai/3dn_6deg_cai.png")
hist(dn3.6deg.cai$CAI, main="Codon adptation index random sequence (dn3 6deg)", xlab=NULL)
abline(v=dn3.6deg.cai[1,2], col=2,lty=2, lwd=4)
text (0.7155, 180, "WT", cex=1.5, col=2)
dev.off()

#n3
png("cai/3n_cai.png")
hist(n3.cai$CAI, main="Codon adptation index random sequence (n3)", xlab=NULL)
abline(v=n3.cai[1,2], col=2,lty=2, lwd=4)
text (0.7155, 180, "WT", cex=1.5, col=2)
dev.off()

#dn23
png("cai/dn23_cai.png")
hist(dn23.cai$CAI, main="Codon adptation index random sequence (dn23)", xlab=NULL)
abline(v=dn3.6deg.enc[1,2], col=2,lty=2, lwd=4)
text (0.7155, 180, "WT", cex=1.5, col=2)
dev.off()

#dn31 
png("cai/dn31_cai.png")
hist(dn31.cai$CAI, main="Codon adptation index random sequence (dn31)", xlab=NULL)
abline(v=dn3.6deg.enc[1,2], col=2,lty=2, lwd=4)
text (0.7155, 180, "WT", cex=1.5, col=2)
dev.off()

