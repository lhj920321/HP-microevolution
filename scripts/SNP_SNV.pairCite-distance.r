##SNP Loci of Each person 
library("ggplot2")
path <- "/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/SNP_SNV_distance/"

Data_File <- paste(path,"SNPminFreq0.95.P8-E-d_distance.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)

sampID <- Data$sampID
type <- Data$type
Distance <- Data$Distance

outFigureF <- paste( Data_File,".boxplot.pdf",sep = "")
rankSumTest <- wilcox.test(Distance~type,data = Data,altermative="two.sides")
print(rankSumTest)

ggplot(data=Data,aes(x=sampID,y=Distance,fill=type)) +
  geom_boxplot() +
  #stat_compare_means(label =  "p.format", label.x = 1.5)+
  #geom_jitter(aes(fill=type),width =0.2,shape = 21,size=1) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(0,5000)
ggsave(outFigureF, width=10, height=4)

  
