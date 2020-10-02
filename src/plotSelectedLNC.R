library(Seurat)
library(tidyverse)
library(patchwork)
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("openxlsx")


suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))
suppressMessages(require(pheatmap))
suppressMessages(require(stringr))
suppressMessages(require(knitr))
suppressMessages(require(data.table))
suppressMessages(require(reticulate))
suppressMessages(require(futile.logger))

seurat.data <- readRDS("./data/2ddiff.fgf8.d16d30d60.seurat.rds")

count.data <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
rownames(count.data) <- count.data[,1]
#datano0 <-data[rowSums(data[,7:48]>0),,drop=FALSE]#Clean the rows that sum 0
dataclean <- subset(count.data[,7:48]) #take the sample columns

for (c in 1:42){ # rename the columns so its easier to read
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
rm(f1,f2,p1,p2,name,fname,c,count.data) #clean

coldata <- data.frame(FGF=factor(c(rep("FGF8_minus",21),rep("FGF8_plus",21))),
                      day=factor(rep(c(rep("day_16",7),rep("day_30",7),rep("day_60",7)),2)), 
                      cell_type=factor(rep(c("Dopamine","Dopamine","Dopamine","No_Dopamine","No_Dopamine","No_Dopamine","No_Dopamine"),6)))

pval=c(0.01)

markers <- read.xlsx("./data/LncRNAToCheck.xlsx")

outdir <- paste0("./output/selectedLNC")
dir.create(paste0(outdir),recursive = TRUE)

####
#Cell design
######
dds_c <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~cell_type)
levels(dds_c$day)
#dds_c$Timepoint <- relevel(dds_c$Timepoint, ref = "Day_16") #Stablish day 16 as reference
dds_c<-DESeq(dds_c)

res_c <- results(dds_c, alpha = pval)
res_c<-res_c[order(res_c$padj),]
head(res_c,50)

allmarkers<-unique(c(markers$Corrected_Ids,markers$Ale_List))

#Plot MA
res_c_df <- as.data.frame(res_c)
NotSig <-subset(res_c_df, padj>pval)
Sig <-subset(res_c_df, padj<pval)
Sig<-rownames_to_column(Sig,var = "ID")

over<-subset(Sig,Sig$log2FoldChange>0)
under<-subset(Sig,Sig$log2FoldChange<0)
res_label<- subset(Sig, Sig$ID %in% allmarkers )


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
                                "#ff0000"))+
  ggrepel::geom_label_repel(data=res_label, aes(label=ID,x=log2(baseMean),
                                                y=log2FoldChange))


ggsave(filename =paste0("MAPlot_Cell_Markers.pdf"), device="pdf",path = outdir )


vsd.c <- varianceStabilizingTransformation(dds_c, blind = TRUE)
markers.heat<- subset(res_c, rownames(res_c) %in% markers$Corrected_Ids)

count.markers.c <- subset(counts(dds_c), rownames(counts(dds_c)) %in% rownames(markers.heat))

df<- as.data.frame(colData(dds_c)[,c("day","cell_type")])
Heat.markers.C<-pheatmap(assay(vsd.c)[rownames(count.markers.c),], cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df, quotes=FALSE,scale = "row")
ggsave(filename ="Heat_Markers_cellDesign.pdf",plot = Heat.markers.C, device="pdf",path = outdir, width = 25, height = 15, units="cm")

Heat.markers.C<-pheatmap(assay(vsd.c)[rownames(count.markers.c),], cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, quotes=FALSE,scale = "row")
ggsave(filename ="Heat_Markers_cellDesign_Tree.pdf",plot = Heat.markers.C, device="pdf",path = outdir,width = 25, height = 25, units="cm")





####⌂
#day-cell
######
dds_dc <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~day+cell_type)
#dds_c$Timepoint <- relevel(dds_c$Timepoint, ref = "Day_16") #Stablish day 16 as reference
dds_dc<-DESeq(dds_dc)

res_dc <- results(dds_dc, alpha = pval)
res_dc<-res_dc[order(res_dc$padj),]
head(res_dc,50)



#Plot MA
res_dc_df <- as.data.frame(res_dc)
NotSig <-subset(res_dc_df, padj>pval)
Sig <-subset(res_dc_df, padj<pval)
Sig<-rownames_to_column(Sig,var = "ID")

over<-subset(Sig,Sig$log2FoldChange>0)
under<-subset(Sig,Sig$log2FoldChange<0)
res_label<- subset(Sig, Sig$ID %in% allmarkers )


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
                                "#ff0000"))+
  ggrepel::geom_label_repel(data=res_label, aes(label=ID,x=log2(baseMean),
                                                y=log2FoldChange))


ggsave(filename =paste0("MAPlot_DAY-Cell_Markers.pdf"), device="pdf",path = outdir )


vsd.dc <- varianceStabilizingTransformation(dds_dc, blind = TRUE)
markers.heat<- subset(res_dc, rownames(res_dc) %in% markers$Corrected_Ids)

count.markers.dc <- subset(counts(dds_dc), rownames(counts(dds_dc)) %in% rownames(markers.heat))

df<- as.data.frame(colData(dds_dc)[,c("day","cell_type")])
Heat.markers.DC<-pheatmap(assay(vsd)[rownames(count.markers.dc),], cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df, quotes=FALSE,scale = "row")
ggsave(filename ="Heat_Markers_day-cellDesign.pdf",plot = Heat.markers.C, device="pdf",path = outdir, width = 25, height = 15, units="cm")

Heat.markers.DC<-pheatmap(assay(vsd)[rownames(count.markers.c),], cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, quotes=FALSE,scale = "row")
ggsave(filename ="Heat_Markers_daycellDesign_Tree.pdf",plot = Heat.markers.C, device="pdf",path = outdir,width = 25, height = 25, units="cm")


#####
#Seurat plot
####
seurat.data@active.assay <-"RNA"
plot_z <- FeaturePlot(seurat.data, features = unique(c(markers$Corrected_Ids,markers$Ale_List)),
            reduction = "umap", col=c("yellow","red","black"),pt.size = 0.5)
ggsave(filename = "UMAP_Feature_plots.pdf",plot=plot_z,path = outdir, device = "pdf", width = 70, height = 50, units = "cm")
