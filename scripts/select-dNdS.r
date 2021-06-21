##dS boxplot
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
Data$person <- factor(Data$person, levels=c("P2","P3", "P11", "P1","P10","P12", "P4","P6","P7", "P5","P9","P8"))  #, ordered=TRUE
E_SNP  <- ggplot(Data, aes(x = person, y = SNPLoci/10, group = factor(1))) + 
  geom_bar(stat = "identity", width = 0.5) + ylab(" SNP Loci of experience group") + xlab("Person")+ ylim(0,600) +  
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
   C_SNP <- ggplot(Data, aes(x = person, y = SNPLoci/10, group = factor(1))) + 
     geom_bar(stat = "identity", width = 0.5) + ylab(" SNP Loci of control group") + xlab("Person") + ylim(0,600) + 
     theme_bw() + 
     theme(panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),
           legend.position = "right")  +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))  #+ ylim(0,1)
   
plot_grid(E,C,E_SNP,C_SNP,labels = c("A","B","C","D"),ncol = 2,nrow = 2)  #
ggsave(outFigureF, width=10, height=10)



#####------------------



##dN
library("ggplot2")
  path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-FiltSomeSamps/"
Data_File <- paste(path,"allPerson.dNdS.txt",sep = "")
Data <- read.table(Data_File,sep = "\t",header = T)
outFigureF <- paste( Data_File,".boxplot.dN.pdf",sep = "")
ggplot(Data,aes(x=reorder(person, dN, mean),y=dN,fill=Samp_Group)) +   ###
  geom_boxplot() +
  #geom_jitter(width =0.2,shape = 21,size=0.5) + 
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  #+ ylim(0,1)
ggsave(outFigureF, width=10, height=4)

##dS 
library("ggplot2")
path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-FiltSomeSamps/"
Data_File <- paste(path,"allPerson.dNdS.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
outFigureF <- paste( Data_File,".boxplot.dS.pdf",sep = "")
ggplot(Data,aes(x=reorder(person, dS, mean),y=dS,fill=Samp_Group)) +
  geom_boxplot() +
  #geom_jitter(width =0.2,shape = 21,size=0.5) + 
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  #+ ylim(0,1)
ggsave(outFigureF, width=10, height=4)




##dN/dS 
library("ggplot2")
path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-FiltSomeSamps/"
Data_File <- paste(path,"allPerson.dNdS.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
outFigureF <- paste( Data_File,".boxplot.dNdS.pdf",sep = "")
ggplot(Data,aes(x=reorder(person, dNdS, mean),y=dNdS,fill=Samp_Group)) +  ###
  geom_boxplot() +
  #geom_point(alpha = 0.1, shape = 21,width =0.2,size=0.5) + 
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  #+ ylim(0,1)
ggsave(outFigureF, width=10, height=4)




  


#######-------------------------------------
#P value  ----  group test
#dN 
library("ggplot2")
path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-FiltSomeSamps/"
Data_File <- paste(path,"allPerson.dNdS.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
outFigureF <- paste( Data_File,".boxplot.dN.GroupTest.pdf",sep = "")

rankSumTest <- wilcox.test(dN~Samp_Group,data = Data,altermative="two.sides")
print(rankSumTest)
ggplot(Data,aes(x=Samp_Group,y=dN,fill=Samp_Group)) +
  geom_boxplot() +
  #stat_compare_means(label =  "p.format", label.x = 1.5)+
  #geom_jitter(aes(fill=type),width =0.2,shape = 21,size=1) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ ylim(0,5)
ggsave(outFigureF, width=10, height=4)


#dS 
library("ggplot2")
path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-FiltSomeSamps/"
Data_File <- paste(path,"allPerson.dNdS.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
outFigureF <- paste( Data_File,".boxplot.dS.GroupTest.pdf",sep = "")

rankSumTest <- wilcox.test(dS~Samp_Group,data = Data,altermative="two.sides")
print(rankSumTest)
ggplot(Data,aes(x=Samp_Group,y=dS,fill=Samp_Group)) +
  geom_boxplot() +
  #stat_compare_means(label =  "p.format", label.x = 1.5)+
  #geom_jitter(aes(fill=type),width =0.2,shape = 21,size=1) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ ylim(0,5)
ggsave(outFigureF, width=10, height=4)



#dN/dS 
library("ggplot2")
path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-FiltSomeSamps/"
Data_File <- paste(path,"allPerson.dNdS.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
outFigureF <- paste( Data_File,".boxplot.dNdS.GroupTest.pdf",sep = "")

rankSumTest <- wilcox.test(dNdS~Samp_Group,data = Data,altermative="two.sides")
print(rankSumTest)
ggplot(Data,aes(x=Samp_Group,y=dNdS,fill=Samp_Group)) +
  geom_boxplot() +
  #stat_compare_means(label =  "p.format", label.x = 1.5)+
  #geom_jitter(aes(fill=type),width =0.2,shape = 21,size=1) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ ylim(0,5)
ggsave(outFigureF, width=10, height=4)







  
  ##----------------------------------
#positive  genes
  personList = c("P1","P2","P3","P4","P5","P6","P7","P9","P10","P12","P15","P16","P20","P21","P22","P23","P24")
  #personList = c("P4")
  
  path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-FiltSomeSamps/"
  for (person in personList) {
  Data_File <- paste(path,person,".codeml.dNdS.over1.txt",sep = "")
  print(Data_File)
  Data <- read.table(Data_File,sep = "\t",header = T)
  #GeneGrp <- Data$GeneGrp
 
  outFigureF <- paste( Data_File,".barplot.pdf",sep = "")
  
  ggplot(Data, aes(x=Samp1__Samp2, y=dNdS)) + geom_bar(stat="identity") + ylim(0,5) + facet_grid(GeneGrp ~ .)
  
  ggsave(outFigureF, width=5, height=50,limitsize = FALSE)
  
  }
  
  
  

  
  
