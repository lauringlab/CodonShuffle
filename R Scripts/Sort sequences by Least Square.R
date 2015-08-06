#Select sequences 

#Read table
dn3.6deg.ham<-read.table(file="Haming distance/dn3_6deg_ham_new_correct.txt")
n3.ham<-read.table(file="Haming distance/n3_ham_new_correct.txt")
dn23.ham<-read.table(file="Haming distance/dn23_ham_new_correct.txt")
dn31.ham<-read.table(file="Haming distance/dn31_ham_new_correct.txt")



#Order table 
dn3.6deg.ham.sort<-dn3.6deg.ham[order(dn3.6deg.ham$dn3.6deg.least.mod),]
n3.ham.sort<-n3.ham[order(n3.ham$n3.least.mod),]
dn23.ham.sort<-dn23.ham[order(dn23.ham$dn23.least.mod),]
dn31.ham.sort<-dn31.ham[order(dn31.ham$dn31.least.mod),]

#Remove column X 
dn3.6deg.ham.sort.new<-dn3.6deg.ham.sort[,-2]
n3.ham.sort.new<-n3.ham.sort[,-2]
dn23.ham.sort.new<-dn23.ham.sort[,-2]
dn31.ham.sort.new<-dn31.ham.sort[,-2]



#Save table 
write.table(dn3.6deg.ham.sort.new, file="Haming distance/dn3_6deg_ham_new_correct_sort.txt", sep="\t")
write.table(n3.ham.sort.new, file="Haming distance/n3_ham_new_correct_sort.txt", sep="\t")
write.table(dn23.ham.sort.new, file="Haming distance/dn23_ham_new_correct_sort.txt", sep="\t")
write.table(dn31.ham.sort.new, file="Haming distance/dn31_ham_new_correct_sort.txt", sep="\t")



