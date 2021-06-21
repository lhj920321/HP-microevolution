##SNP Loci of Each person 
library("ggplot2")
library("reshape2")
path <- "/shared/liuhj/HP/process/Coregenome-SNP-Loci/"

Data_File <- paste(path,"coreGnm.Diff-Loci.stat.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)
Num <- Data$SNPLoci
flag <- Data$group
rankSumTest <- wilcox.test(Num~flag,data = Data,altermative="two.sides")
print(rankSumTest)
outFigureF <- paste( Data_File,".boxplot.pdf",sep = "")

ggplot(data=Data,aes(x=flag,y=Num)) +
  geom_boxplot() +
  geom_jitter(aes(fill=flag),width =0.2,shape = 21,size=2) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  
ggsave(outFigureF, width=3, height=2)

print(rankSumTest)




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






