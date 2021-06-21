##dS boxplot  and SNP Loci Num 
library("ggplot2")
library("cowplot")
path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-allSamp/"
Data_File <- paste(path,"Group.E.dNdS.txt",sep = "")
Data <- read.table(Data_File,sep = "\t",header = T)
outFigureF <- paste( path,"dNdS_select-allSamp.CoreGnm.pdf",sep = "")

E <-ggplot(Data,aes(x=reorder(person, dS, mean),y=dS)) +
 # geom_boxplot( outlier.shape = NA) + ylab(" dS ") + xlab("")+ ylim(0,0.4) +  
  geom_boxplot() + ylab(" dS ") + xlab("")+ ylim(0,0.4) +  
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  #+ ylim(0,1)



Data_File <- paste(path,"Group.C.dNdS.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
C <-ggplot(Data,aes(x=reorder(person, dS, mean),y=dS)) +
  geom_boxplot() + xlab("")+ ylab("") + ylim(0,0.4) + 
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

  
#E
SNPLociPath <- "/shared/liuhj/HP/process/Coregenome-SNP-Loci/"
Data_File <- paste(SNPLociPath,"coreGnm.Diff-Loci-E.stat.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
Data$person <- factor(Data$person, levels=c("P14","P13","P2","P3", "P11", "P1","P10","P12", "P4","P6","P7", "P5","P9","P8"))  #, ordered=TRUE
E_SNP_1  <- ggplot(Data, aes(x = person, y = SNPLoci/1000, group = factor(1))) + 
  geom_bar(stat = "identity") + ylab(" Number of SNP Loci ") + xlab("")+ ylim(0,4.0) +  
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


E_SNP_2 <- ggplot(Data,aes(x=person,y=SNPLoci/1000,group = factor(1))) + 
  geom_bar(stat='identity') +
  labs(x=NULL,y=NULL,fill=NULL) +   #不要标签
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +     #去掉X轴和X轴的文字
  coord_cartesian(ylim = c(24.5,25)) +  #设置上面一半的值域
  scale_y_continuous(breaks = c(24.5,25,0.5)) #以0.55为单位划分Y轴

E_SNP = ggarrange(E_SNP_2,E_SNP_1,heights=c(1/5, 4/5),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")


## core genome SNP Loci 
#C
SNPLociPath <-"/shared/liuhj/HP/process/Coregenome-SNP-Loci/"
Data_File <- paste(SNPLociPath,"coreGnm.Diff-Loci-C.stat.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)

Data$person <- factor(Data$person, levels=c("P15","P19", "P17","P22","P25", "P21","P18","P20", "P16","P23","P24"))
C_SNP_1 <- ggplot(Data, aes(x = person, y = SNPLoci/1000, group = factor(1))) + 
  geom_bar(stat = "identity") + xlab("")   + ylab("") + ylim(0,4.0) + 
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


C_SNP_2 <- ggplot(Data,aes(x=person,y=SNPLoci/1000,group = factor(1))) + 
  geom_bar(stat='identity') +
  labs(x=NULL,y=NULL,fill=NULL) +   #不要标签
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +     #去掉X轴和X轴的文字
  coord_cartesian(ylim = c(24.5,25)) +  #设置上面一半的值域
  scale_y_continuous(breaks = c(24.5,25,0.5)) #以0.55为单位划分Y轴

C_SNP = ggarrange(C_SNP_2,C_SNP_1,heights=c(1/5, 4/5),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")


## core genome SNP Loci Anno Stat
#C
path <- "/shared/liuhj/HP/process/Coregenome-SNP-Loci/"
Data_File <- paste(path,"coreGnm.Diff-Loci.stat.Cgroup.Anno.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
Data_mele <- melt(Data,id.vars="person",variable.name = "Anno",value.name= "AnnoNum")
Data_mele$person <- factor(Data_mele$person, levels=c("P15","P19", "P17","P22","P25", "P21","P18","P20", "P16","P23","P24"))

C_Anno <-  ggplot(Data_mele,aes(x=person,y=AnnoNum/1000,fill = Anno)) +
  geom_bar(stat = "identity") + ylim(0,4.0) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none')  +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

#outFigure <- paste( Data_File,".bar.pdf",sep = "")
#ggsave(outFigure, width=10, height=4)



##C
path <- "/shared/liuhj/HP/process/Coregenome-SNP-Loci/"
Data_File <- paste(path,"coreGnm.Diff-Loci.stat.Egroup.Anno.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
Data_mele <- melt(Data,id.vars="person",variable.name = "Anno",value.name= "AnnoNum")
Data_mele$person <- factor(Data_mele$person, levels=c("P14","P13","P2","P3", "P11", "P1","P10","P12", "P4","P6","P7", "P5","P9","P8"))  #, ordered=TRUE

E_Anno <-  ggplot(Data_mele,aes(x=person,y=AnnoNum/1000,fill = Anno)) +
  geom_bar(stat = "identity") + ylim(0,4.0) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none')  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

#outFigure <- paste( Data_File,".bar.pdf",sep = "")
#ggsave(outFigure, width=10, height=4)


plot_grid(E,C,E_SNP,C_SNP,ncol = 2,nrow = 2)  #,labels = c("A","B","C","D")
ggsave(outFigureF, width=6, height=6)

  

##########---------
#boxplot,and test (diff loci num : E and C)
##SNP Loci of Each person based on the core genome of the assembly genome
library("ggplot2")
path <- "/shared/liuhj/HP/process/Coregenome-SNP-Loci/"

Data_File <- paste(path,"coreGnm.Diff-Loci.stat_noP8.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)

rankSumTest <- wilcox.test(SNPLoci~group,data = Data,altermative="two.sides")
print(rankSumTest)
Data$group <- factor(Data$group, levels=c("E","C"))  #, ordered=TRUE

ggplot(data=Data,aes(x=group,y=SNPLoci)) +
  stat_boxplot(geom= 'errorbar', width = 0.3) +
  geom_boxplot() +
  #geom_jitter(aes(fill=flag),width =0.2,shape = 21,size=2) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right") 
outFigureF <- paste( Data_File,".boxplot.pdf",sep = "")
ggsave(outFigureF, width=4, height=5)
print(rankSumTest)




##bar plot of SNP Anno (all person)
path <- "/shared/liuhj/HP/process/Coregenome-SNP-Loci/"
Data_File <- paste(path,"coreGnm.Diff-Loci.stat.Anno.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
Data_mele <- melt(Data,id.vars="person",variable.names = "group",value.name= "AnnoNum")

ggplot(Data_mele,aes(x=person,y=AnnoNum,fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  
outFigureF <- paste( Data_File,".bar.pdf",sep = "")
ggsave(outFigureF, width=10, height=4)





##########---------
#boxplot,and test (dS : E and C)
##dS of Each person based on the core genome of the assembly genome
library("ggplot2")
path <- "/shared/liuhj/HP/process/dN-dS-caculate-CoreGnm/dNdS_select-allSamp/"

Data_File <- paste(path,"allPerson.dNdS-withoutP8.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
rankSumTest <- wilcox.test(dS~Samp_Group,data = Data,altermative="two.sides")
print(rankSumTest)
Data$Samp_Group <- factor(Data$Samp_Group, levels=c("E","C"))  #, ordered=TRUE  
ggplot(data=Data,aes(x=Samp_Group,y=dS)) +
  #stat_boxplot(geom= 'errorbar', width = 0.3) +
  geom_boxplot() +
 # geom_boxplot(outlier.colour = NA) +
  # ylim(0,0.00005) + 
  #geom_jitter(aes(fill=flag),width =0.2,shape = 21,size=2) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right") 
outFigureF <- paste( Data_File,".boxplot.pdf",sep = "")
ggsave(outFigureF, width=4, height=5)
print(rankSumTest)





