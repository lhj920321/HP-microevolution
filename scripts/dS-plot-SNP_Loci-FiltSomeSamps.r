##dS boxplot  and SNP Loci Num 
library("ggplot2")
library("cowplot")
path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-FiltSomeSamps/"
Data_File <- paste(path,"Group.E.dNdS.txt",sep = "")
Data <- read.table(Data_File,sep = "\t",header = T)
outFigureF <- paste( path,"dNdS_select-FiltSomeSamps.Stat.pdf",sep = "")

E <-ggplot(Data,aes(x=reorder(person, dS, mean),y=dS)) +
  geom_boxplot() + ylab(" dS of experience group") + xlab("Person")+
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  #+ ylim(0,1)



Data_File <- paste(path,"Group.C.dNdS.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
C <-ggplot(Data,aes(x=reorder(person, dS, mean),y=dS)) +
  geom_boxplot() + ylab(" dS of control group") + xlab("Person")+
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  #+ ylim(0,1)


SNPLociPath <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/"
Data_File <- paste(SNPLociPath,"person_SNP_Loci.0.95.stat-E.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
Data$person <- factor(Data$person, levels=c("P3","P2","P10", "P11", "P1", "P4","P12","P7","P9", "P5"))  #, ordered=TRUE
E_SNP  <- ggplot(Data, aes(x = person, y = SNPLoci, group = factor(1))) + 
  geom_bar(stat = "identity", width = 0.5) + ylab(" SNP Loci of experience group") + xlab("Person")+ ylim(0,6000) +  
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  #+ ylim(0,1)



SNPLociPath <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/"
Data_File <- paste(SNPLociPath,"person_SNP_Loci.0.95.stat-C.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)

Data$person <- factor(Data$person, levels=c("P15","P19", "P17","P22","P25", "P21","P18","P20", "P16","P23","P24"))
C_SNP <- ggplot(Data, aes(x = person, y = SNPLoci, group = factor(1))) + 
  geom_bar(stat = "identity", width = 0.5) + ylab(" SNP Loci of control group") + xlab("Person") +  ylim(0,6000) + 
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  

plot_grid(E,C,E_SNP,C_SNP,labels = c("A","B","C","D"),ncol = 2,nrow = 2)  #
ggsave(outFigureF, width=10, height=6)


