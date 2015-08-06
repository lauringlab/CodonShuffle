#Read file 
dn231.ham<-read.table(file="Haming distance/dn3_6deg_ham_new_correct.txt", header=TRUE)
n3.ham<-read.table(file="Haming distance/n3_ham_new_correct.txt", header=TRUE)
dn23.ham<-read.table(file="Haming distance/dn23_ham_new_correct.txt", header=TRUE)
dn31.ham<-read.table(file="Haming distance/dn31_ham_new_correct.txt", header=TRUE)


#Eliminate one column
dn231.ham<-dn231.ham[,-2]
n3.ham<-n3.ham[,-2]
dn23.ham<-dn23.ham[,-2]
dn31.ham<-dn31.ham[,-2]

#Change the name of the column
colnames(dn231.ham)<-c("ls", "nd")
colnames(n3.ham)<-c("ls", "nd")
colnames(dn23.ham)<-c("ls", "nd")
colnames(dn31.ham)<-c("ls", "nd")

#Create the graph
dn231.ham.graph<-ggplot(dn231.ham,aes(ls, nd))+geom_point()+ggtitle("Hamming Distance (dn231)")+xlab("Least Square")+ylab("Nucleotide distance")+scale_y_continuous(limits=c(450, 650))+scale_x_continuous(limits=c(0, 11))
n3.ham.graph<-ggplot(n3.ham,aes(ls, nd))+geom_point()+ggtitle("Hamming Distance (n3)")+xlab("Least Square")+ylab("Nucleotide distance")+scale_y_continuous(limits=c(450, 650))+scale_x_continuous(limits=c(0, 11))
dn23.ham.graph<-ggplot(dn23.ham,aes(ls, nd))+geom_point()+ggtitle("Hamming Distance (dn23)")+xlab("Least Square")+ylab("Nucleotide distance")+scale_y_continuous(limits=c(450, 650))+scale_x_continuous(limits=c(0, 11))
dn31.ham.graph<-ggplot(dn31.ham,aes(ls, nd))+geom_point()+ggtitle("Hamming Distance (dn31)")+xlab("Least Square")+ylab("Nucleotide distance")+scale_y_continuous(limits=c(450, 650))+scale_x_continuous(limits=c(0, 11))


#Save the graphs

pdf(file="/Users/danieljorge/Desktop/Documents Lab/Scripts/Poliovirus_P1/Graphs/dn231_ggplot_hamming_distance.pdf")
ggplot(dn231.ham,aes(ls, nd))+geom_point()+ggtitle("Hamming Distance (dn231)")+xlab("Least Square")+ylab("Nucleotide distance")+scale_y_continuous(limits=c(450, 650))+scale_x_continuous(limits=c(0, 11))
dev.off()

pdf(file="/Users/danieljorge/Desktop/Documents Lab/Scripts/Poliovirus_P1/Graphs/n3_ggplot_hamming_distance.pdf")
ggplot(n3.ham,aes(ls, nd))+geom_point()+ggtitle("Hamming Distance (n3)")+xlab("Least Square")+ylab("Nucleotide distance")+scale_y_continuous(limits=c(450, 650))+scale_x_continuous(limits=c(0, 11))
dev.off()

pdf(file="/Users/danieljorge/Desktop/Documents Lab/Scripts/Poliovirus_P1/Graphs/dn23_ggplot_hamming_distance.pdf")
ggplot(dn23.ham,aes(ls, nd))+geom_point()+ggtitle("Hamming Distance (dn23)")+xlab("Least Square")+ylab("Nucleotide distance")+scale_y_continuous(limits=c(450, 650))+scale_x_continuous(limits=c(0, 11))
dev.off()

pdf(file="/Users/danieljorge/Desktop/Documents Lab/Scripts/Poliovirus_P1/Graphs/dn31_ggplot_hamming_distance.pdf")
ggplot(dn31.ham,aes(ls, nd))+geom_point()+ggtitle("Hamming Distance (dn31)")+xlab("Least Square")+ylab("Nucleotide distance")+scale_y_continuous(limits=c(450, 650))+scale_x_continuous(limits=c(0, 11))
dev.off()

#Run multiplot script before (see end of file)
multiplot(dn231.ham.graph, n3.ham.graph, dn23.ham.graph, dn31.ham.graph, cols=2)



