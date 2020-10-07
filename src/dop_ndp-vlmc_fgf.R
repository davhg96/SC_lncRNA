library("DESeq2")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")
library("openxlsx")
library("RColorBrewer")

data <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
rownames(data) <- data[,1]


dataclean <- subset(data[,28:48]) #take the sample columns 7:48 for all of them, 28 till the end for fgf8+


for (c in 1:ncol(dataclean)){ # rename the columns so its easier to read
  name <- strsplit(colnames(dataclean[c]),"_")[[1]][4:5] #split the whole line and take the interesting parts
  p1<-strsplit(name[1], "\\.")[[1]][2:3]#select and clean part one (FGF and day)
  f1 <- paste(p1[1],p1[2], sep=".")
  p2<-strsplit(name[2],"\\.")[[1]][1:2]#Select and clean the second part
  
  if(p2[2]=="bam"){#check and delete the bam extension
    f2<- p2[1]
    fname <- paste(f1, f2, sep=".")#joint the whole name and rename
    colnames(dataclean)[c]<-fname
  }
  else{
    f2 <- paste(p2[1],p2[2], sep=".")
    fname <- paste(f1, f2, sep=".")#join the whole name and rename
    colnames(dataclean)[c]<-fname
  }
}
rm(f1,f2,p1,p2,name,fname,c,data) #clean


coldata <- data.frame(day=factor(c(rep("day_16",7),rep("day_30",7),rep("day_60",7))), 
                      cell_type=factor(rep(c("Dopamine","Dopamine","Dopamine","No_Dopamine","No_Dopamine","No_Dopamine","VLMC"),3)))

pval=c(0.05,0.01)
outputDir <- "./output/Dop_Ndop-VLMC_FGF+/pval_"

for(pval in pval){
  
  #############
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
  NotSig <-subset(res_d_df, padj>pval)
  Sig <-subset(res_d_df, padj<pval)
  Sig<-rownames_to_column(Sig,var = "ID")
  
  over<-subset(Sig,Sig$log2FoldChange>0)
  under<-subset(Sig,Sig$log2FoldChange<0)
  
  ggplot()+
    geom_point(data=NotSig,aes(x=log2(baseMean),
                               y=log2FoldChange, colour="Not Significant"))+
    geom_point(data=over,aes(x=log2(baseMean),
                             y=log2FoldChange,colour="Upregulated"))+
    geom_point(data=under,aes(x=log2(baseMean),
                              y=log2FoldChange,colour="Downregulated"))+
    scale_color_manual(name = "",
                       values = c("#0000ff",
                                  "#000000",
                                  "#ff0000"))
  ggsave(filename =paste0("MAPlot_Day_",pval,".pdf"), device="pdf",path = outdir, width = 16,height = 9,units = "cm" )
  
  
  
  #Table for GO top 50
  write.xlsx(Sig,file =paste0(outdirT,"Top_significant_Day.xlsx"))
  
  
  #HEATMAP
  
  vsd <- varianceStabilizingTransformation(dds_d, blind = TRUE)
  top100 <- res_d[1:100,]
  
  counts_sorted<-counts(dds_d)
  counts_sorted<-counts_sorted[match(rownames(top100),rownames(counts_sorted)),]
  
  
  countTop50 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:50,]
  
  
  df<- as.data.frame(colData(dds_d)[,c("day","cell_type")])
  p50NT<-pheatmap(assay(vsd)[rownames(countTop50),], 
                  cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df, quotes=FALSE,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Day_top50byPvalVSD_NoTree_day_cell.pdf",plot = p50NT, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )
  
  p50T<-pheatmap(assay(vsd)[rownames(countTop50),], 
                 cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, quotes=FALSE,
                 scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Day_top50byPvalVSD_Tree_day_cell.pdf",plot = p50T, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )
  
  
  countTop100 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:100,]
  
  df<- as.data.frame(colData(dds_d)[,c("day","cell_type")])
  p100NT<-pheatmap(assay(vsd)[rownames(countTop100),], 
                   cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df, quotes=FALSE,
                   scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Day_top100byPvalVSD_NoTree_day_cell.pdf",plot = p100NT, device="pdf",path = outdir, width = 21,height = 32,units = "cm" )
  
  p100T<-pheatmap(assay(vsd)[rownames(countTop100),], 
                  cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, quotes=FALSE,
                  scale = "row",col = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Day_top100byPvalVSD_Tree_day_cell.pdf",plot = p100T, device="pdf",path = outdir, width = 21,height = 32,units = "cm" )
  
  
  
  # Heatmap of similarity between replicates
  distVSD <- dist(t(assay(vsd)))
  matrix <- as.matrix(distVSD)
  rownames(matrix) <- paste(vsd$day,vsd$cell_type, sep = "-")
  colnames(matrix) <- paste(vsd$day,vsd$cell_type, sep = "-")
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE,show_colnames = TRUE, fontsize = 15,
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_day_NoTree.pdf",plot = dheatmap, device="pdf",path = outdir, width = 25,height = 25,units = "cm" )
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = TRUE,show_colnames = TRUE, 
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_day_Tree.pdf",plot = dheatmap, device="pdf",path = outdir, width = 25,height = 25,units = "cm" )
  
  # PCA plot
  rld <- rlogTransformation(dds_d,blind=TRUE)
  pca<-plotPCA(rld, intgroup=c("day"))
  ggsave(filename ="PCA_INT_Day.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
  
  pca<-plotPCA(rld, intgroup=c("cell_type"))
  ggsave(filename ="PCA_INT_CellType.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
  
  #############
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
  NotSig <-subset(res_c_df, padj>pval)
  Sig <-subset(res_c_df, padj<pval)
  Sig<-rownames_to_column(Sig,var = "ID")
  
  over<-subset(Sig,Sig$log2FoldChange>0)
  under<-subset(Sig,Sig$log2FoldChange<0)
  
  ggplot()+
    geom_point(data=NotSig,aes(x=log2(baseMean),
                               y=log2FoldChange, colour="Not Significant"))+
    geom_point(data=over,aes(x=log2(baseMean),
                             y=log2FoldChange,colour="Upregulated"))+
    geom_point(data=under,aes(x=log2(baseMean),
                              y=log2FoldChange,colour="Downregulated"))+
    scale_color_manual(name = "",
                       values = c("#0000ff",
                                  "#000000",
                                  "#ff0000"))
  ggsave(filename =paste0("MAPlot_Cell_",pval,".pdf"), device="pdf",path = outdir, width = 16,height = 9,units = "cm" )
  
  
  
  #Table for GO top 50
  
  write.xlsx(Sig,file =paste0(outdirT,"Top_significant_Cell.xlsx"))
  
  
  #HEATMAP
  
  vsd <- varianceStabilizingTransformation(dds_c, blind = TRUE)
  top100 <- res_c[1:100,]
  
  counts_sorted<-counts(dds_c)
  counts_sorted<-counts_sorted[match(rownames(top100),rownames(counts_sorted)),]
  
  
  countTop50 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:50,]
  
  
  df<- as.data.frame(colData(dds_c)[,c("day","cell_type")])
  p50NT<-pheatmap(assay(vsd)[rownames(countTop50),], 
                  cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df, quotes=FALSE,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Cell_top50byPvalVSD_NoTree_day_cell.pdf",plot = p50NT, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )
  
  p50T<-pheatmap(assay(vsd)[rownames(countTop50),], 
                 cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, quotes=FALSE,
                 scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Cell_top50byPvalVSD_Tree_day_cell.pdf",plot = p50T, device="pdf",path = outdir, width = 21,height = 21,units = "cm" )
  
  
  
  countTop100 <- subset(counts_sorted,  rownames(counts_sorted) %in% rownames(top100))[1:100,]
  
  df<- as.data.frame(colData(dds_c)[,c("day","cell_type")])
  p100NT<-pheatmap(assay(vsd)[rownames(countTop100),], 
                   cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df, quotes=FALSE,
                   scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Cell_top100byPvalVSD_NoTree_day_cell.pdf",plot = p100NT, device="pdf",path = outdir, width = 21,height = 32,units = "cm" )
  
  p100T<-pheatmap(assay(vsd)[rownames(countTop100),], 
                  cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, quotes=FALSE,
                  scale = "row",color = colorRampPalette(c("blue", "white",  "red"))(50))
  ggsave(filename ="Heatmap_Cell_top100byPvalVSD_Tree_day_cell.pdf",plot = p100T, device="pdf",path = outdir, width = 21,height = 32,units = "cm" )
  
  
  
  # Heatmap of similarity between replicates
  distVSD <- dist(t(assay(vsd)))
  matrix <- as.matrix(distVSD)
  rownames(matrix) <- paste(vsd$day,vsd$cell_type, sep = "-")
  colnames(matrix) <- paste(vsd$day,vsd$cell_type, sep = "-")
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE,show_colnames = TRUE, fontsize = 15,
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_Cell_NoTree.pdf",plot = dheatmap, device="pdf",path = outdir, width = 25,height = 25,units = "cm" )
  
  
  dheatmap<-pheatmap(matrix,clustering_distance_rows = distVSD, 
                     clustering_distance_cols = distVSD,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = TRUE,show_colnames = TRUE, 
                     color = hmcol, main = "Distance Matrix")
  ggsave(filename ="Heatmap_Distances_Cell_Tree.pdf",plot = dheatmap, device="pdf",path = outdir, width = 25,height = 25,units = "cm" )
  
  # PCA plot
  rld <- rlogTransformation(dds_c,blind=TRUE)
  pca<-plotPCA(rld, intgroup=c("day"))
  ggsave(filename ="PCA_INT_Day.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
  
  pca<-plotPCA(rld, intgroup=c("cell_type"))
  ggsave(filename ="PCA_INT_CellType.pdf",plot = pca, device="pdf",path = outdir, width = 15,height = 15,units = "cm" )
  
}
