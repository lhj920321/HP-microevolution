          ## core genome SNP Loci 
          #
      library("ggplot2")
          Path <-"/shared/liuhj/HP/process/plot/SNV_Stat/"
          Data_File <- paste(Path,"allPerson.NGS.SNV_SNP.0.95.Stat.txt",sep = "")
          print(Data_File)
          Data <- read.table(Data_File,sep = "\t",header = T)
          #ggplot(Data,aes(x = sample, y = Filted_snvCount, fill=seqTech )) + 
            #geom_bar(stat = "identity",position = "dodge") + xlab("")   + ylab("") + 
          ggplot(Data,aes(x = sample, y = Filted_snvCount)) + 
            geom_bar(stat = "identity") + xlab("")   + ylab("") + ylim(0,10000) + 
            theme_bw() + #ylim(0,10000) + 
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  legend.position = "right")  +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
          
          outFigureF <- paste( Data_File,".bar.pdf",sep = "")
          ggsave(outFigureF, width=8, height=3)
        
          
        ####boxplot
    
        rankSumTest <- wilcox.test(Filted_snvCount~group,data = Data,altermative="two.sides")
        print(rankSumTest)
        Data$group <- factor(Data$group, levels=c("E","C"))  #, ordered=TRUE
        ggplot(data=Data,aes(x=group,y=Filted_snvCount)) +
          stat_boxplot(geom= 'errorbar', width = 0.3) +
          geom_boxplot(outlier.colour = NA) + ylim(0,6000) +
          #geom_jitter(aes(fill=group),width =0.2,shape = 21,size=1) +
          theme_bw() + 
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                legend.position = "right") 
        outFigureF <- paste( Data_File,".boxplot.pdf",sep = "")
        ggsave(outFigureF, width=3, height=2.5)
        print(rankSumTest)
      
    
 ###No.  CCS vs NGS    
    Path <-"/shared/liuhj/HP/process/plot/SNV_Stat/"
    Data_File <- paste(Path,"CCS-NGS.SNV_0.05-0.95.Stat.txt",sep = "")
    print(Data_File)
    Data <- read.table(Data_File,sep = "\t",header = T)
    
    Data$seqTech <- factor(Data$seqTech, levels=c("NGS","CCS"))  #, ordered=TRUE
    ggplot(Data,aes(x = sample, y = Filted_snvCount, fill=seqTech )) + 
    geom_bar(stat = "identity",position = "dodge") + xlab("")   + ylab("") + 
      theme_bw() + #ylim(0,10000) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            legend.position = "right")  +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    
    outFigureF <- paste( Path,"CCS-NGS.SNV_0.05-0.95.Stat.pdf",sep = "")
    ggsave(outFigureF, width=8, height=3)
    
    
    ###No.  CCS vs NGS   - 0.9 
    Path <-"/shared/liuhj/HP/process/plot/SNV_Stat/"
    Data_File <- paste(Path,"CCS-NGS.SNV_0.1-0.9.Stat.txt",sep = "")
    print(Data_File)
    Data <- read.table(Data_File,sep = "\t",header = T)
    
    Data$seqTech <- factor(Data$seqTech, levels=c("NGS","CCS"))  #, ordered=TRUE
    ggplot(Data,aes(x = sample, y = Filted_snvCount, fill=seqTech )) + 
      geom_bar(stat = "identity",position = "dodge") + xlab("")   + ylab("") + 
      theme_bw() + #ylim(0,10000) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            legend.position = "right")  +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    outFigureF <- paste( Path,"CCS-NGS.SNV_0.1-0.9.Stat.pdf",sep = "")
    ggsave(outFigureF, width=8, height=3)
    
    
  ## SNV venn fig plot 
    library("VennDiagram")
  CCS_Samps = c("P18-C_bR-j_2_P18-C_bR-x","P18-C_bR-x_2_P18-C_bR-x",
               "P19-C_bR-d_2_P19-C_bR-j","P19-C_bR-j_2_P19-C_bR-j",
               "P19-C_bR-t_2_P19-C_bR-j","P1-E-d_2_P1-E-j",
               "P1-E-j_2_P1-E-j","P1-E-t_2_P1-E-j","P1-E-x_2_P1-E-j",
               "P20-C_bR-d_2_P20-C_bR-t","P20-C_bR-j_2_P20-C_bR-t",
               "P20-C_bR-t_2_P20-C_bR-t","P20-C_bR-x_2_P20-C_bR-t",
               "P21-C_sL-d_2_P21-C_sL-j","P21-C_sL-j_2_P21-C_sL-j",
               "P21-C_sL-t_2_P21-C_sL-j","P21-C_sL-x_2_P21-C_sL-j",
               "P22-C_sL-d_2_P22-C_sL-d",
               "P23-C_sL-j_2_P23-C_sL-j","P23-C_sL-x_2_P23-C_sL-j",
               "P24-C_sC-j_2_P24-C_sC-t","P24-C_sC-t_2_P24-C_sC-t",
               "P24-C_sC-x_2_P24-C_sC-t","P25-C_sC-j_2_P25-C_sC-t",
               "P25-C_sC-t_2_P25-C_sC-t","P2-E-t_2_P2-E-t",
               "P4-E-d_2_P4-E-j","P4-E-j_2_P4-E-j","P4-E-x_2_P4-E-j",
               "P6-E-d_2_P6-E-t","P6-E-t_2_P6-E-t",
               "P7-E-d_2_P7-E-j","P7-E-j_2_P7-E-j")  #"P22-C_sL-t_2_P22-C_sL-d","P5-E-d_2_P5-E-j"
 
  #CCS_Samps = c("P7-E-d_2_P7-E-j","P7-E-j_2_P7-E-j")
   for (CCS_Samp in CCS_Samps) {
      print(CCS_Samp)
    CCSPath <-"/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/all_CCS_SNVPosi/"
    Data_File <- paste(CCSPath,CCS_Samp,".SNV.PosiAndFreq.txt",sep = "")
    print(Data_File)
    CCSData <- read.table(Data_File,sep = "\t",header = T)
    CCS_SNVPosi = CCSData$SNVPosi
    
    NGSPath <-"/shared/liuhj/HP/process/NGS_map2PersonResptGnm/all_NGS_SNVPosi/"
    Data_File <- paste(NGSPath,CCS_Samp,".SNV.PosiAndFreq.txt",sep = "")
    print(Data_File)
    NGSData <- read.table(Data_File,sep = "\t",header = T)
    NGS_SNVPosi =NGSData$SNVPosi
    
    input  <- list(NGS_SNVPosi,CCS_SNVPosi)
    OutPath <-"/shared/liuhj/HP/process/plot/SNV_CCS-NGS_VennFigs/"
      
    outF = paste(OutPath,CCS_Samp,".SNV.PosiAndFreq.VennFig.pdf",sep = "")
    vennFig <- venn.diagram(list(NGS=NGS_SNVPosi,CCS=CCS_SNVPosi), fill=c("red","green"), alpha=c(0.3,0.3),cex = 3, filename = NULL)
    pdf(file=outF)
    grid.draw(vennFig)
    dev.off()
    
   }  
  
    
############NGS Freq bar figs
      personLst <- c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24","P25")
    personLst <- c("P11")
    for (person in personLst) {
        if(person == "P1"){
        ymax=500
      }
      PersonPath <-paste("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/person_SNV_freqDistrib/",sep = "")
    Data_File <- paste(PersonPath,person,".iSNV_0.95.Freq_distrib.txt",sep = "")
    print(Data_File)
    
    Data <- read.table(Data_File,sep = "\t",header = F)
    ggplot(Data,aes(x = V2, y = V3,fill=V4)) + 
      geom_bar(stat = "identity") + xlab("")   + ylab("") + 
      #scale_fill_manual(values = c("red", "grey") )+
      theme_bw() + #ylim(0,10000) + 
      geom_hline(aes(yintercept=50),linetype="dashed") + 
      facet_grid(. ~ V1) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            legend.position = "right")  +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    outFigureF <- paste( Data_File,".bar.pdf",sep = "")
    ggsave(outFigureF, width=8, height=2)
  
    }
    
 
      
      ##test perosn
          PersonPath <-paste("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/person_SNV_freqDistrib/",sep = "")
          Data_File <- paste(PersonPath,"allPerson.Freq_distrub.txt",sep = "")
          print(Data_File)
          
          Data <- read.table(Data_File,sep = "\t",header = F)
          ggplot(Data,aes(x = V2, y = V3,fill=V5)) + 
            geom_bar(stat = "identity") + xlab("")   + ylab("") + 
            #scale_fill_manual(values = c("red", "grey") )+
            theme_bw() + #ylim(0,10000) + 
            geom_hline(aes(yintercept=50),linetype="dashed") + 
            facet_grid(V4 ~ V1) +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  legend.position = "right")  +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
          outFigureF <- paste( Data_File,".bar.pdf",sep = "")
          ggsave(outFigureF, width=8, height=50,limitsize = FALSE)

    
    
  ### CCS Freq bar figs
  personLst <- c("P1","P2","P4","P5","P6","P7","P18","P19","P20","P21","P22","P23","P24","P25")
  personLst <- c("P1","P2")
  for (person in personLst) {

    
    PersonPath <-paste("/shared/liuhj/HP/process/CCS_map2PersonRespectSamp/person_SNV_freqDistrib/",sep = "")
    #Files <- list.files(PersonPath,pattern="*Freq_distribu.txt") 
    #for (File in Files) {
      Data_File <- paste(PersonPath,person,".FreqDistrib.txt",sep = "")
      print(Data_File)
      Data <- read.table(Data_File,sep = "\t",header = F)
      ggplot(Data,aes(x = V2, y = V3)) + 
        geom_bar(stat = "identity") + xlab("")   + ylab("") +
         
        theme_bw() + #ylim(0,10000) + 
        geom_hline(aes(yintercept=50),linetype="dashed") + 
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              legend.position = "right")  +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
      outFigureF <- paste( Data_File,".bar.pdf",sep = "")
      ggsave(outFigureF, width=10, height=2)
      #ggplot(Data, aes(x=Samp1__Samp2, y=dNdS)) + geom_bar(stat="identity") + ylim(0,5)  
   # }
  }

  
    