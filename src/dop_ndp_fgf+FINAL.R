library("DESeq2")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")
library("openxlsx")
library("RColorBrewer")

data <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
rownames(data) <- data[,1]
#datano0 <-data[rowSums(data[,7:48]>0),,drop=FALSE]#Clean the rows that sum 0

dataclean <- subset(data[,28:48]) #take the sample columns 7:48 for all of them, 28 till the end for fgf8+


for (c in 1:ncol(dataclean)){ # rename the columns so its easier to read
  name <- strsplit(colnames(dataclean[c]),"_")[[1]][4:5] #split the whole line and take the interesting parts
  p1<-strsplit(name[1], "\\.")[[1]][3]#select and clean part one (FGF and day)
  #f1 <- paste(p1[1],p1[2], sep=".")
  p2<-strsplit(name[2],"\\.")[[1]][1:2]#Select and clean the second part
  
  if(p2[2]=="bam"){#check and delete the bam extension
    f2<- p2[1]
    fname <- paste(p1, f2, sep=".")#joint the whole name and rename
    colnames(dataclean)[c]<-fname
  }
  else{
    f2 <- paste(p2[1],p2[2], sep=".")
    fname <- paste(p1, f2, sep=".")#join the whole name and rename
    colnames(dataclean)[c]<-fname
  }
}
rm(f2,p2,p1,name,c) #clean


coldata <- data.frame(day=factor(c(rep("day_16",7),rep("day_30",7),rep("day_60",7)),
                                 labels = c("Day 16","Day 30","Day 60")), 
                      cell_type=factor(rep(c("Dopamine","Dopamine","Dopamine","No_Dopamine","No_Dopamine","No_Dopamine","No_Dopamine"),3),
                                       labels = c("Dopaminergic neurons","Floorplate")))

pval=c(0.01)
outputDir <- "./output/Dop_NdopFGF+09-10/pval_"


for(pval in pval){
  
  # Day DEsign#
  #############
  #create the object
  dds_d <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~day)
  levels(dds_d$day)
  #ddsc$Timepoint <- relevel(ddsc$Timepoint, ref = "Day_16") #Stablish day 16 as reference
  dds_d<-DESeq(dds_d)
  
  dir.create(paste0(outputDir,pval,"/day/"),recursive = TRUE)
  outdir <- paste0(outputDir,pval,"/day/")
  dir.create(paste0(outputDir,pval,"/Table"),recursive = TRUE)
  outdirT<-paste0(outputDir,pval,"/Table/")
  
  res_d <- results(dds_d, alpha = pval)
  res_d<-res_d[order(res_d$padj),]
  head(res_d,50)
  
  
  # sink(file = (paste0(outdir,"Summary_Day.txt")))
  # summary(res_d)
  # sink()
  
  #Plot MA
  res_d_df <- as.data.frame(res_d)
  res_d_df <- na.omit(res_d_df)
  res_d_df$state=rep(3,nrow(res_d_df))
  for (r in 1:nrow(res_d_df)){
    if(res_d_df$padj[r]<pval & res_d_df$log2FoldChange[r]>0){
      res_d_df$state[r] <- "Upregulated"
    }
    if(res_d_df$padj[r]<pval & res_d_df$log2FoldChange[r]<0){
      res_d_df$state[r] <- "Downregulated"
    }
    if(res_d_df$padj[r]>pval) {
      res_d_df$state[r] <- "Not significant"
    }
  }
  res_d_df$state <- as.factor(res_d_df$state)
  res_d_df$state <- factor(res_d_df$state, levels = c("Upregulated","Downregulated","Not significant"))
  expcounts = res_d_df %>%
    group_by(state) %>%
    summarise(Count = n()) %>%
    mutate(share = round(Count / sum(Count), digits = 2)) 
  
  ggplot()+
    geom_point(data=res_d_df,aes(x=log2(baseMean),
                               y=log2FoldChange, color=state),size=0.3)+
    theme_minimal()+
    scale_color_manual(name = "Expression",
                       labels=paste0(expcounts$state," [", expcounts$Count,"]"),
                       values = c("#F43E3E",
                                  "#3EA1F4",
                                  "#888686"))+
    labs(title="MA plot p<0.01")
  ggsave(filename =paste0("MAPlot_Day_",pval,".pdf"), device="pdf",path = outdir, width = 16,height = 9,units = "cm" )
  
  
  
  #Table for GO top 50
  write.xlsx(subset(res_d, res_d$padj<pval),file =paste0(outdirT,"Top_significant",pval,"_Day.xlsx"))
  
  
  #HEATMAP
  
  vsd <- varianceStabilizingTransformation(dds_d, blind = TRUE)
  top100 <- res_d[1:100,]
  counts_sorted<-counts(dds_d)
  counts_sorted<-counts_sorted[match(rownames(top100),rownames(counts_sorted)),]
  
  countTop50 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:50,]
  df<- as.data.frame(colData(dds_d)[,c("day","cell_type")])
  colnames(df) <- c("Timepoint", "Cell population")
  
  
  plotdata <- assay(vsd)[rownames(countTop50),]
  plotdf <- plotdata[order(plotdata[,18],decreasing = TRUE),]
  
  
  my_colour = list(
    "Timepoint" = c("Day 16" = "#f4d63e", "Day 30" = "#f4a53e","Day 60"="#F43E3E"),
    "Cell population" = c("Dopaminergic neurons" = "#14de3c", "Floorplate" = "#de149e"))
  

  
  p50NT<-pheatmap(plotdf,  
                  cluster_rows=FALSE,
                  show_rownames=TRUE, 
                  cluster_cols=FALSE, 
                  annotation_col=df, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late","VLMC"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Day_top50byPvalVSD_NoTree_day_cell.pdf",plot = p50NT, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )
  

  p50NT<-pheatmap(plotdf,  
                  cluster_rows=TRUE,
                  show_rownames=TRUE, 
                  cluster_cols=FALSE, 
                  annotation_col=df, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late","VLMC"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Day_top50byPvalVSD_ClusterRow.pdf",plot = p50NT, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )

  
  
  
  countTop100 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:100,]
  df<- as.data.frame(colData(dds_d)[,c("day","cell_type")])
  colnames(df) <- c("Timepoint", "Cell population")
  
  plotdata <- assay(vsd)[rownames(countTop100),]
  plotdf <- plotdata[order(plotdata[,18],decreasing = TRUE),]
  
  p100NT<-pheatmap(plotdf, 
                   cluster_rows=FALSE,
                   show_rownames=TRUE, 
                   cluster_cols=FALSE, 
                   annotation_col=df, 
                   annotation_colors = my_colour,
                   quotes=FALSE,
                   labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late","VLMC"),3),
                   angle_col = 45,
                   scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Day_top100byPvalVSD_NoTree_day_cell.pdf",plot = p100NT, device="pdf",path = outdir, width = 21,height = 32,units = "cm" )
  


  
  # PCA plot
  rld <- rlogTransformation(dds_d,blind=TRUE)
  pca<-plotPCA(rld, intgroup=c("day"))
  ggsave(filename ="PCA_INT_Day.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
  
  pca<-plotPCA(rld, intgroup=c("cell_type"))
  ggsave(filename ="PCA_INT_CellType.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
  
# celltype DEsign#
  #############
  #create the object
  dds_c <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~cell_type)
  levels(dds_c$day)
  #dds_c$Timepoint <- relevel(dds_c$Timepoint, ref = "Day_16") #Stablish day 16 as reference
  dds_c<-DESeq(dds_c)
  
  dir.create(paste0(outputDir,pval,"/Cell/"),recursive = TRUE)
  outdir <- paste0(outputDir,pval,"/Cell/")
  
  
  res_c <- results(dds_c, alpha = pval)
  res_c<-res_c[order(res_c$padj),]
  head(res_c,50)
  
  
  sink(file = (paste0(outdir,"Summary_cell.txt")))
  summary(res_c)
  sink()
  
  #Plot MA
  res_c_df <- as.data.frame(res_c)
  res_c_df <- as.data.frame(res_c)
  res_c_df <- na.omit(res_c_df)
  res_c_df$state=rep(3,nrow(res_c_df))
  for (r in 1:nrow(res_c_df)){
    if(res_c_df$padj[r]<pval & res_c_df$log2FoldChange[r]>0){
      res_c_df$state[r] <- "Upregulated"
    }
    if(res_c_df$padj[r]<pval & res_c_df$log2FoldChange[r]<0){
      res_c_df$state[r] <- "Downregulated"
    }
    if(res_c_df$padj[r]>pval) {
      res_c_df$state[r] <- "Not significant"
    }
  }
  res_c_df$state <- as.factor(res_c_df$state)
  res_c_df$state <- factor(res_c_df$state, levels = c("Upregulated","Downregulated","Not significant"))
  expcounts = res_c_df %>%
    group_by(state) %>%
    summarise(Count = n()) %>%
    mutate(share = round(Count / sum(Count), digits = 2)) 
  
  ggplot()+
    geom_point(data=res_c_df,aes(x=log2(baseMean),
                                 y=log2FoldChange, color=state),size=0.3)+
    theme_minimal()+
    scale_color_manual(name = "Expression",
                       labels=paste0(expcounts$state," [", expcounts$Count,"]"),
                       values = c("#F43E3E",
                                  "#3EA1F4",
                                  "#888686"))+
    labs(title="MA plot p<0.01")
  ggsave(filename =paste0("MAPlot_Cell_",pval,".pdf"), device="pdf",path = outdir, width = 16,height = 9,units = "cm" )
  
  
  
  #Table for GO top 50
  
  write.xlsx(subset(res_c, res_c$padj<pval),file =paste0(outdirT,"Top_significant",pval,"_Cell.xlsx"))
  
  
  #HEATMAP
  
  vsd <- varianceStabilizingTransformation(dds_c, blind = TRUE)
  top100 <- res_c[1:100,]
  
  counts_sorted<-counts(dds_c)
  counts_sorted<-counts_sorted[match(rownames(top100),rownames(counts_sorted)),]
  
  
  countTop50 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:50,]

  df<- as.data.frame(colData(dds_c)[,c("day","cell_type")])
  colnames(df) <- c("Timepoint", "Cell population")
  
  plotdata <- assay(vsd)[rownames(countTop50),]
  plotdf <- plotdata[order(plotdata[,2],decreasing = TRUE),]
  
  p50NT<-pheatmap(plotdf, 
                  cluster_rows=FALSE,
                  show_rownames=TRUE, 
                  cluster_cols=FALSE, 
                  annotation_col=df, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late","VLMC"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Cell_top50byPvalVSD_NoTree_day_cell.pdf",plot = p50NT, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )

  
  p50NT<-pheatmap(plotdf, 
                  cluster_rows=TRUE,
                  show_rownames=TRUE, 
                  cluster_cols=FALSE, 
                  annotation_col=df, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late","VLMC"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Cell_top50byPvalVSD_ClusterRow.pdf",plot = p50NT, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )

  
  countTop100 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:100,]
  df<- as.data.frame(colData(dds_c)[,c("day","cell_type")])
  colnames(df) <- c("Timepoint", "Cell population")
  
  plotdata <- assay(vsd)[rownames(countTop100),]
  plotdf <- plotdata[order(plotdata[,18],decreasing = TRUE),]
  
  p100NT<-pheatmap(plotdf, 
                   cluster_rows=FALSE,
                   show_rownames=TRUE, 
                   cluster_cols=FALSE, 
                   annotation_col=df, 
                   annotation_colors = my_colour,
                   quotes=FALSE,
                   labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late","VLMC"),3),
                   angle_col = 45,
                   scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Cell_top100byPvalVSD_NoTree_day_cell.pdf",plot = p100NT, device="pdf",path = outdir, width = 21,height = 32,units = "cm" )
  


  
  # PCA plot
  rld <- rlogTransformation(dds_c,blind=TRUE)
  pca<-plotPCA(rld, intgroup=c("day"))
  ggsave(filename ="PCA_INT_Day.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
  
  pca<-plotPCA(rld, intgroup=c("cell_type"))
  ggsave(filename ="PCA_INT_CellType.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  

}




#Violin
###########
outV <- paste0(outputDir,"0.01/ViolinCounts/")
dir.create(outV,recursive = TRUE)

LNCdata <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
rownames(LNCdata) <- data[,1]

countLNC <- subset(LNCdata[,28:48]) #take the sample columns 7:48 for all of them, 28 till the end for fgf8+


for (c in 1:ncol(countLNC)){ # rename the columns so its easier to read
  name <- strsplit(colnames(countLNC[c]),"_")[[1]][4:5] #split the whole line and take the interesting parts
  #p1<-strsplit(name[1], "\\.")[[1]][3]#select and clean part one (FGF and day)
  #f1 <- paste(p1[1],p1[2], sep=".")
  p2<-strsplit(name[2],"\\.")[[1]][1:2]#Select and clean the second part
  
  if(p2[2]=="bam"){#check and delete the bam extension
    f2<- p2[1]
    #fname <- paste(p1, f2, sep=".")#joint the whole name and rename
    colnames(countLNC)[c]<-f2
  }
  else{
    f2 <- paste(p2[1],p2[2], sep=".")
    #fname <- paste(p1, f2, sep=".")#join the whole name and rename
    colnames(countLNC)[c]<-f2
  }
}
rm(f2,p2,name,c,LNCdata) #clean

plotCountLNC <- data.frame(day=factor(c(rep("day_16",7),rep("day_30",7),rep("day_60",7)),labels = c("Day 16","Day 30","Day 60")),
                        cell=factor(rep(c("Dopamine","Dopamine","Dopamine","No_Dopamine","No_Dopamine","No_Dopamine","No_Dopamine"),3),labels = c("Dopaminergic neurons","Inmature cells")),
                        count=colSums(countLNC))


PCGdata <- read.table(file = "./data/Fcounts_2Ddiff_PcodingGenes_s1.txt", header = TRUE)
rownames(PCGdata) <- PCGdata[,1]
countPCG <- subset(PCGdata[,28:48])

for (c in 1:ncol(countPCG)){ # rename the columns so its easier to read
  name <- strsplit(colnames(countPCG[c]),"_")[[1]][4:5] #split the whole line and take the interesting parts
  #p1<-strsplit(name[1], "\\.")[[1]][3]#select and clean part one (FGF and day)
  #f1 <- paste(p1[1],p1[2], sep=".")
  p2<-strsplit(name[2],"\\.")[[1]][1:2]#Select and clean the second part
  
  if(p2[2]=="bam"){#check and delete the bam extension
    f2<- p2[1]
    #fname <- paste(p1, f2, sep=".")#joint the whole name and rename
    colnames(countPCG)[c]<-f2
  }
  else{
    f2 <- paste(p2[1],p2[2], sep=".")
    #fname <- paste(p1, f2, sep=".")#join the whole name and rename
    colnames(countPCG)[c]<-f2
  }
}
rm(f2,p2,name,c,PCGdata) #clean


plotCountPCG <- data.frame(day=factor(c(rep("day_16",7),rep("day_30",7),rep("day_60",7)),labels = c("Day 16","Day 30","Day 60")),
                           cell=factor(rep(c("Dopamine","Dopamine","Dopamine","No_Dopamine","No_Dopamine","No_Dopamine","No_Dopamine"),3),labels = c("Dopaminergic neurons","Inmature cells")),
                           count=colSums(countPCG))

totalCounts=data.frame(sample=factor(c(rep("lncRNA",21),rep("PCG",21))))
totalCounts <- cbind(totalCounts, rbind(plotCountLNC,plotCountPCG))

violin <- ggplot(totalCounts,aes(x=cell, y=count, fill=cell))+
  facet_grid(rows = vars(sample),cols = vars(day),scales = "free")+
  geom_violin() +
  scale_fill_manual(name="Cell population",
                    labels=c("Dopaminergic neurons","Inmature cells"),
                    values = c("#F43E3E","#3EA1F4"))+
  theme_minimal()+
  labs(y="Total lncRNA counts")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

ggsave("countVplot.pdf", plot = violin,path = outV, device = "pdf", width = 21,height = 9, units = "cm" )


violin <- ggplot(totalCounts,aes(x=sample, y=count, fill=sample))+
  facet_grid(cols = vars(day))+
  geom_violin() +
  scale_fill_manual(name="Sample",
                    labels=c("lncRNA","Protein coding genes"),
                    values = c("#F43E3E","#3EA1F4"))+
  theme_minimal()+
  labs(y="Total lncRNA counts")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

ggsave("countVplotlnc_pcg.pdf", plot = violin,path = outV, device = "pdf", width = 21,height = 9, units = "cm" )

