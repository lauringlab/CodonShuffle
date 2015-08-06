#Minimum free energy (Mfold)
dn3.6deg.fold <- read.table(file = "Unafold//3dn_6deg_RNA.fas.dG")                     
n3.fold <- read.table(file = "Unafold//3n_RNA.fas.dG")                     
dn23.fold <- read.table(file = "Unafold/dn23_RNA.fas.dG")                 
dn31.fold <- read.table(file = "Unafold/dn31_RNA.fas.dG")                     
entero.fold <- read.table(file = "Unafold/entero_final_seq_mod_RNA.fas.dG")                     

#par(mgp=c(2.3, 0.5, 0))
 
#Graph all scripts together 
png("mfold_new_all scripts.png", width = 900, height = 900)
par(mfrow=c(2,2), mgp=c(2.3, 0.5, 0))
hist(dn3.6deg.fold$V2, main = "Mfold (dn3_6deg)", xlab = "Energy", xlim=c(-600,-700), cex.axis=2.5, cex.lab=2.5, cex.main=3)
abline(v=dn3.6deg.fold[1,2], col=2,lty=2, lwd=4)
text (-668, 180, "WT", cex=1.5, col=2)
hist(n3.fold$V2, main = "Mfold (n3)", xlab = "Energy", xlim=c(-600,-700), cex.axis=2.5, cex.lab=2.5, cex.main=3)
abline(v=n3.fold[1,2], col=2,lty=2, lwd=4)
text (-668, 130, "WT", cex=1.5, col=2)
hist(dn23.fold$V2, main = "Mfold (dn23)", xlab = "Energy", xlim=c(-600,-700), cex.axis=2.5, cex.lab=2.5, cex.main=3)
abline(v=dn23.fold[1,2], col=2,lty=2, lwd=4)
text (-668, 220, "WT", cex=1.5, col=2)
hist(dn31.fold$V2, main = "Mfold (dn31)", xlab = "Energy", xlim=c(-600,-700), cex.axis=2.5, cex.lab=2.5, cex.main=3)
abline(v=dn31.fold[1,2], col=2,lty=2, lwd=4)
text (-670, 180, "WT", cex=1.5, col=2)
dev.off()



#Graph each script individually 
#dn3 6deg
png("dn3_6deg_mfold_new.png", width = 900, height = 900)
par(mfrow=c(2,2), mgp=c(2.3, 0.5, 0))
hist(dn3.6deg.fold$V2, main = "Mfold (dn3_6deg)", xlab = "Energy", xlim=c(-600,-700), cex.axis=2.5, cex.lab=2.5, cex.main=3)
abline(v=dn3.6deg.fold[1,2], col=2,lty=2, lwd=4)
text (-668, 180, "WT", cex=1.5, col=2)
dev.off()

#n3
png("n3_mfold_new.png", width = 900, height = 900)
par(mgp=c(2.3, 0.5, 0))
hist(n3.fold$V2, main = "Mfold (n3)", xlab = "Energy", xlim=c(-600,-700), cex.axis=2.5, cex.lab=2.5, cex.main=3)
abline(v=n3.fold[1,2], col=2,lty=2, lwd=4)
text (-668, 130, "WT", cex=1.5, col=2)
dev.off()

#dn23
png("dn23_mfold_new.png", width = 900, height = 900)
par(mgp=c(2.3, 0.5, 0))
hist(dn23.fold$V2, main = "Mfold (dn23)", xlab = "Energy", xlim=c(-600,-700), cex.axis=2.5, cex.lab=2.5, cex.main=3)
abline(v=dn23.fold[1,2], col=2,lty=2, lwd=4)
text (-668, 220, "WT", cex=1.5, col=2)
dev.off()

#dn31
png("dn31_mfold_new.png", width = 900, height = 900)
par(mgp=c(2.3, 0.5, 0))
hist(dn31.fold$V2, main = "Mfold (dn31)", xlab = "Energy", xlim=c(-600,-700), cex.axis=2.5, cex.lab=2.5, cex.main=3)
abline(v=dn31.fold[1,2], col=2,lty=2, lwd=4)
text (-670, 180, "WT", cex=1.5, col=2)
dev.off()

#Enterovirus (reference)
png("entero_mfold_new.png", width = 900, height = 900)
par(mgp=c(2.3, 0.5, 0))
hist(entero.fold$V2, main = "Mfold (enterovirus)", xlab = "Energy", cex.axis=2.5, cex.lab=2.5, cex.main=3.5)
abline(v=entero.fold[1,2], col=2,lty=2, lwd=4)
text (-750, 50, "WT", cex=1.5, col=2)
dev.off()



