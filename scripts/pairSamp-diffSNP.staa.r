##pair samps -diff SNP loci      
library("ggplot2")

path <- "/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/pairSamp-diffSNPLoci_noSNVLoci/"

Data_File <- paste(path,"allPersons.iSNV.pairSamps.diffSNPLoci.Stat.txt",sep = "")
print(Data_File)
Data <- read.table(Data_File,sep = "\t",header = T)


Num <- Data$Num_SNPLoci
flag <- Data$group

rankSumTest <- wilcox.test(Num~flag,data = Data,altermative="two.sides")
print(rankSumTest)

outFigureF <- paste( Data_File,".boxplot.pdf",sep = "")
  
ggplot(data=Data,aes(x=flag,y=Num)) +
  geom_boxplot(position = position_dodge(width = 0.6), outlier.colour="NA", width = 0.7, show.legend = TRUE) +
  stat_compare_means(label =  "p.format", label.x = 1.5)+
  geom_dotplot(binwidth = 1,binaxis = "y",stackdir = "center",fill="NA",dotsize=5) +  
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "right")  
ggsave(outFigureF, width=3, height=2)

  
    
        
        ##SNP Loci of Each person 
        library("ggplot2")
        path <- "/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/"
        
        Data_File <- paste(path,"person_SNP_Loci.0.9.stat.txt",sep = "")
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
        
