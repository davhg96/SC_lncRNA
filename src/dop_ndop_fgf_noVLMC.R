library("DESeq2")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")
library("openxlsx")


data <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
rownames(data) <- data[,1]
#datano0 <-data[rowSums(data[,7:48]>0),,drop=FALSE]#Clean the rows that sum 0

dataclean <- subset(data[,c(28:33,35:40,42:47)]) #take the sample columns 6 out of 7 from 28 onwards


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

ordered_samples <- c("day16.FP.Cycling","day16.FP.Early","day16.FP.Late","day16.DA.E1","day16.DA.1","day16.DA.2",
                     "day30.FP.Cycling","day30.FP.Early","day30.FP.Late","day30.DA.E1","day30.DA.1","day30.DA.2",
                     "day60.FP.Cycling","day60.FP.Early","day60.FP.Late","day60.DA.E1","day60.DA.1","day60.DA.2")
dataclean <-  dataclean[,ordered_samples ]

rm(f2,p2,p1,name,c, data, ordered_samples) #clean




coldata <- data.frame(day=factor(c(rep("day_16",6),rep("day_30",6),rep("day_60",6)),
                                 labels  = c("Day 16","Day 30","Day 60")), 
                      cell_type=factor(rep(c("No_Dopamine","No_Dopamine","No_Dopamine","Dopamine","Dopamine","Dopamine"),3),
                                  labels   = c("Dopaminergic neurons","Floorplate")))

pval=c(0.01)
outputDir <- "./output/FinalPlots/pval_"


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
  
  res_d_30_16 <- results(dds_d, alpha = pval, contrast = c("day", "Day 30","Day 16"))
  res_d_30_16<-res_d_30_16[order(res_d_30_16$padj),]
  
  head(res_d_30_16,50)
  
  res_d_60_16 <- results(dds_d, alpha = pval, contrast = c("day", "Day 60","Day 16"))
  res_d_60_16<-res_d_60_16[order(res_d_60_16$padj),]
  head(res_d_60_16,50)
  
  res_d_30_16_df <- as.data.frame(res_d_30_16)
  res_d_30_16_df <- na.omit(res_d_30_16_df)
  
  res_d_60_16_df <- as.data.frame(res_d_60_16)
  res_d_60_16_df <- na.omit(res_d_60_16_df)
  
  
  clasifyExp <- function(resdf){
    resdf$state=rep(3,nrow(resdf))
    for (r in 1:nrow(resdf)){
      if(resdf$padj[r]<pval & resdf$log2FoldChange[r]>0){
        resdf$state[r] <- "Upregulated"
      }
      if(resdf$padj[r]<pval & resdf$log2FoldChange[r]<0){
        resdf$state[r] <- "Downregulated"
      }
      if(resdf$padj[r]>pval) {
        resdf$state[r] <- "Not significant"
      }
    }
    resdf$state <- as.factor(resdf$state)
    resdf$state <- factor(resdf$state, levels = c("Upregulated","Downregulated","Not significant")) #Reorder the levels to the 
    #order we want in the plot
    return(resdf)
  }# anotate the expression
  
  res_d_30_16_df <- clasifyExp(res_d_30_16_df)
  res_d_60_16_df <- clasifyExp(res_d_60_16_df)
  
 

  #####
  #Plots day 30 and 16 (16 as reference)
  ######
  dir.create(paste0(outputDir,pval,"/day/day30vs16/"),recursive = TRUE)
  outdird1 <- paste0(outputDir,pval,"/day/day30vs16/")
  
  expcounts = res_d_30_16_df %>%
    group_by(state) %>%
    summarise(Count = n()) %>%
    mutate(share = round(Count / sum(Count), digits = 2)) 
  
  

  ggplot()+
    geom_point(data=res_d_30_16_df,aes(x=log2(baseMean),
                                 y=log2FoldChange, color=state),size=0.3)+
    theme_minimal()+
    scale_color_manual(name = "Expression",
                       labels=paste0(expcounts$state," [", expcounts$Count,"]"),
                       values = c("#F43E3E",
                                  "#3EA1F4",
                                  "#888686"))+
    labs(title="MA plot p<0.01")
  ggsave(filename =paste0("MAPlot_day30Vs16(reference)Pval",pval,".pdf"), device="pdf",path = outdird1, width = 16,height = 9,units = "cm" )
  
  
  #Table for GO top 50
  write.xlsx(rownames_to_column(as.data.frame(subset(res_d_30_16_df, res_d_30_16_df$padj<pval)), var = "ID"),file =paste0(outdirT,"Top_significant_day30Vs16(reference)Pval",pval,".xlsx"))
  
  
  ###### 
  #HEATMAP day30vs16
  #####
  vsd <- varianceStabilizingTransformation(dds_d, blind = TRUE)
  top100 <- as.data.frame(res_d_30_16[1:100,])
  rownames(top100)
  top100 <- top100[!grepl("MIR*",rownames(top100)),] #Delete lines that start with MIR (microRNAs)
  
  counts_sorted<-counts(dds_d)
  counts_sorted<-counts_sorted[match(rownames(top100),rownames(counts_sorted)),]
  
  countTop50 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:50,]


  plotdata <- assay(vsd)[rownames(countTop50),]
  #plotdf <- plotdata[order(plotdata[,18],decreasing = TRUE),]
  
  p50NT<-pheatmap(plotdata,  
                  cluster_rows=FALSE,
                  show_rownames=TRUE, 
                  show_colnames = TRUE,
                  cluster_cols=FALSE, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  #labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Top50_Day30Vs16(reference)NoTree.pdf",plot = p50NT, device="pdf",path = outdird1, width = 21,height = 21,units = "cm" )
  
  
  p50NT<-pheatmap(plotdata,  
                  cluster_rows=TRUE,
                  show_rownames=TRUE, 
                  show_colnames = TRUE,
                  cluster_cols=FALSE, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  #labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Top50_Day30Vs16(reference)_Tree.pdf",plot = p50NT, device="pdf",path = outdird1, width = 21,height = 21,units = "cm" )
  

  # Heatmap of similarity between replicates
  distVSD <- dist(t(assay(vsd)))
  matrix <- as.matrix(distVSD)
  rownames(matrix) <- paste(vsd$day,rep(c("FP.Cycling","FP.Early","FP.Late","DA.E1","DA.1","DA.2"),3), sep = "-")
  colnames(matrix) <- paste(vsd$day,rep(c("FP.Cycling","FP.Early","FP.Late","DA.E1","DA.1","DA.2"),3), sep = "-")
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE,show_colnames = TRUE, fontsize = 15,
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_day_NoTree.pdf",plot = dheatmap, device="pdf",path = outdird1, width = 25,height = 25,units = "cm" )
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = TRUE,show_colnames = TRUE, 
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_day_Tree.pdf",plot = dheatmap, device="pdf",path = outdird1, width = 25,height = 25,units = "cm" )
  
  #####
  #Plot MA day 60 and 16 (16 as reference)
  ######
  
  dir.create(paste0(outputDir,pval,"/day/day60vs16/"),recursive = TRUE)
  outdird1 <- paste0(outputDir,pval,"/day/day60vs16/")
  
  expcounts = res_d_60_16_df %>%
    group_by(state) %>%
    summarise(Count = n()) %>%
    mutate(share = round(Count / sum(Count), digits = 2)) 
  
  ggplot()+
    geom_point(data=res_d_60_16_df,aes(x=log2(baseMean),
                                       y=log2FoldChange, color=state),size=0.3)+
    theme_minimal()+
    scale_color_manual(name = "Expression",
                       labels=paste0(expcounts$state," [", expcounts$Count,"]"),
                       values = c("#F43E3E",
                                  "#3EA1F4",
                                  "#888686"))+
    labs(title="MA plot p<0.01")
  ggsave(filename =paste0("MAPlot_day60Vs16(reference)Pval",pval,".pdf"), device="pdf",path = outdird1, width = 16,height = 9,units = "cm" )
  
  
  #Table for GO top 50
  write.xlsx(rownames_to_column(as.data.frame(subset(res_d_60_16_df, res_d_60_16_df$padj<pval)), var = "ID"),file =paste0(outdirT,"Top_significant_MAPlot_day60Vs16(reference)Pva",pval,".xlsx"))

  ###### 
  #HEATMAP day60Vs16
  #####
  
  vsd <- varianceStabilizingTransformation(dds_d, blind = TRUE)
  top100 <- as.data.frame(res_d_60_16[1:100,])
  rownames(top100)
  top100 <- top100[!grepl("MIR*",rownames(top100)),] #Delete lines that start with MIR (microRNAs)
  
  counts_sorted<-counts(dds_d)
  counts_sorted<-counts_sorted[match(rownames(top100),rownames(counts_sorted)),]
  
  countTop50 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:50,]
  
  plotdata <- assay(vsd)[rownames(countTop50),]
  
  p50NT<-pheatmap(plotdata,  
                  cluster_rows=FALSE,
                  show_rownames=TRUE, 
                  show_colnames = TRUE,
                  cluster_cols=FALSE, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  #labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Top50_Day60Vs16(reference)NoTree.pdf",plot = p50NT, device="pdf",path = outdird1, width = 21,height = 21,units = "cm" )
  
  
  p50NT<-pheatmap(plotdata,  
                  cluster_rows=TRUE,
                  show_rownames=TRUE, 
                  show_colnames = TRUE,
                  cluster_cols=FALSE, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  #labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Top50_Day60Vs16(reference)_Tree.pdf",plot = p50NT, device="pdf",path = outdird1, width = 21,height = 21,units = "cm" )
  
  
  # Heatmap of similarity between replicates
  distVSD <- dist(t(assay(vsd)))
  matrix <- as.matrix(distVSD)
  rownames(matrix) <- paste(vsd$day,rep(c("FP.Cycling","FP.Early","FP.Late","DA.E1","DA.1","DA.2"),3), sep = "-")
  colnames(matrix) <- paste(vsd$day,rep(c("FP.Cycling","FP.Early","FP.Late","DA.E1","DA.1","DA.2"),3), sep = "-")
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE,show_colnames = TRUE, fontsize = 15,
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_day_NoTree.pdf",plot = dheatmap, device="pdf",path = outdird1, width = 25,height = 25,units = "cm" )
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = TRUE,show_colnames = TRUE, 
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_day_Tree.pdf",plot = dheatmap, device="pdf",path = outdird1, width = 25,height = 25,units = "cm" )
  
  #####
 
  #####
  # PCA plot
  #####
  rld <- rlogTransformation(dds_d,blind=TRUE)
  pca<-plotPCA(rld, intgroup=c("day"))
  ggsave(filename ="PCA_INT_Day.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
  
  pca<-plotPCA(rld, intgroup=c("cell_type"))
  ggsave(filename ="PCA_INT_CellType.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
  #####
  # celltype DEsign#
  #############
  #create the object
  dds_c <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~cell_type)
  levels(dds_c$day)
  #dds_c$Timepoint <- relevel(dds_c$Timepoint, ref = "Day_16") #Stablish day 16 as reference
  dds_c<-DESeq(dds_c)
  
  dir.create(paste0(outputDir,pval,"/Cell/"),recursive = TRUE)
  outdir <- paste0(outputDir,pval,"/Cell/")
  
  
  res_c <- results(dds_c, alpha = pval,contrast = c("cell_type","Dopaminergic neurons", "Floorplate"))
  res_c<-res_c[order(res_c$padj),]
  head(res_c,50)
  
  
  
  #####
  #Plots Celldesign Floorplate as reference
  #####
  res_c_df <- as.data.frame(res_c)
  res_c_df <- na.omit(res_c_df)
  
  res_c_df <- clasifyExp(res_c_df)
  
  
  
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
  ggsave(filename =paste0("MAPlot_Cell_DopVs Floorplate(reference)_",pval,".pdf"), device="pdf",path = outdir, width = 16,height = 9,units = "cm" )
  
  
  
  #Table for GO top 50
  
  write.xlsx(rownames_to_column(as.data.frame(subset(res_c, res_c$padj<pval)),var = "ID"),file =paste0(outdirT,"Top_significant",pval,"_Cell_dopvsFloorplate(reference).xlsx"))
  
  #####
  #HEATMAP cell types
  #####
  
  vsd <- varianceStabilizingTransformation(dds_c, blind = TRUE)
  top100 <- as.data.frame(res_c_df[1:100,])
  rownames(top100)
  top100 <- top100[!grepl("MIR*",rownames(top100)),] #Delete lines that start with MIR (microRNAs)
  
  counts_sorted<-counts(dds_c)
  counts_sorted<-counts_sorted[match(rownames(top100),rownames(counts_sorted)),]
  
  countTop50 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:50,]
  
  plotdata <- assay(vsd)[rownames(countTop50),]
  
  p50NT<-pheatmap(plotdata,  
                  cluster_rows=FALSE,
                  show_rownames=TRUE, 
                  show_colnames = TRUE,
                  cluster_cols=FALSE, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  #labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Top50_Cell_dopvsFloorplate(reference)NoTree.pdf",plot = p50NT, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )
  
  
  p50NT<-pheatmap(plotdata,  
                  cluster_rows=TRUE,
                  show_rownames=TRUE,
                  show_colnames = TRUE,
                  cluster_cols=FALSE, 
                  annotation_colors = my_colour,
                  quotes=FALSE,
                  #labels_col = rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late"),3),
                  angle_col = 45,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Top50_Cell_dopvsFloorplate(reference)Tree.pdf",plot = p50NT, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )
  
  
 
  
  
  # Heatmap of similarity between replicates
  distVSD <- dist(t(assay(vsd)))
  matrix <- as.matrix(distVSD)
  rownames(matrix) <- paste(vsd$day,rep(c("FP.Cycling","FP.Early","FP.Late","DA.E1","DA.1","DA.2"),3), sep = "-")
  colnames(matrix) <- paste(vsd$day,rep(c("FP.Cycling","FP.Early","FP.Late","DA.E1","DA.1","DA.2"),3), sep = "-")
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE,show_colnames = TRUE, fontsize = 15,
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_cell_NoTree.pdf",plot = dheatmap, device="pdf",path = outdir, width = 25,height = 25,units = "cm" )
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = TRUE,show_colnames = TRUE, 
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_cell_Tree.pdf",plot = dheatmap, device="pdf",path = outdir, width = 25,height = 25,units = "cm" )
  
  
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

violin <-  ggplot(totalCounts,aes(x=day, y=count, fill=day))+
  facet_grid(rows = vars(sample))+
  geom_violin() +
  scale_fill_manual(name="Timepoint",
                    labels=c("Day 16","Day 30","Day 60"),
                    values = c("#f4d63e","#f4a53e","#F43E3E"))+
  theme_minimal()+
  labs(y="Total lncRNA counts")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

ggsave("countVplotSameScale.pdf", plot = violin,path = outV, device = "pdf", width = 21,height = 9, units = "cm" )

violin <-  ggplot(totalCounts,aes(x=day, y=count, fill=day))+
  facet_grid(rows = vars(sample),scales = "free")+
  geom_violin() +
  scale_fill_manual(name="Timepoint",
                    labels=c("Day 16","Day 30","Day 60"),
                    values = c("#f4d63e","#f4a53e","#F43E3E"))+
  theme_minimal()+
  labs(y="Total lncRNA counts")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

ggsave("countVplotFreeScale.pdf", plot = violin,path = outV, device = "pdf", width = 21,height = 9, units = "cm" )


violin <- ggplot(totalCounts,aes(x=sample, y=count, fill=sample))+
  facet_grid(cols = vars(day))+
  geom_violin() +
  scale_fill_manual(name="Sample",
                    labels=c("lncRNA","Protein coding genes"),
                    values = c("#F43E3E","#3EA1F4"))+
  theme_minimal()+
  labs(y="Total counts")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

ggsave("countVplotlnc_pcg.pdf", plot = violin,path = outV, device = "pdf", width = 21,height = 9, units = "cm" )


violin <- ggplot(plotCountLNC,aes(x=day, y=count, fill=day))+
  geom_violin() +
  scale_fill_manual(name="Timepoint",
                    labels=c("Day 16","Day 30","Day 60"),
                    values = c("#f4d63e","#f4a53e","#F43E3E"))+
  theme_minimal()+
  labs(y="Total counts")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12))

ggsave("countVplotLNC.pdf", plot = violin,path = outV, device = "pdf", width = 16 ,height = 9, units = "cm" )


# violin <-  
ggplot(plotCountPCG,aes(x=day, y=count, fill=day))+
  geom_violin() +
  scale_fill_manual(name="Timepoint",
                    labels=c("Day 16","Day 30","Day 60"),
                    values = c("#f4d63e","#f4a53e","#F43E3E"))+
  theme_minimal()+
  labs(y="Total counts")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12))

ggsave("countVplotPCG.pdf", plot = violin,path = outV, device = "pdf", width = 16 ,height = 9, units = "cm" )

