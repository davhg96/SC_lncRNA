library("DESeq2")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")
library("openxlsx")
library("RColorBrewer")
library("GenomicRanges")


D_cell <- read.xlsx("./output/Dop_NdopFGF+/pval_0.01/table/Top_significant_Cell.xlsx",rowNames = TRUE)
D_day <- read.xlsx("./output/Dop_NdopFGF+/pval_0.01/table/Top_significant_Day.xlsx",rowNames = TRUE)


common <- subset(D_day, rownames(D_day) %in% rownames(D_cell))                 
common <- cbind(common, D_cell[rownames(common),])#First 6 cols are for Day and last 6 for cell DF

PCGdata <- read.table(file = "./data/Fcounts_2Ddiff_PcodingGenes_s1.txt", header = TRUE)
PCGdata <- column_to_rownames(PCGdata, var = "Geneid")
Lncdata <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
Lncdata <- column_to_rownames(Lncdata, var = "Geneid")


takeCoordinates <- function(df){
  for (c in 1: nrow(df)){
    df[c,1] <- strsplit(df[c,1], ";" )[[1]][1]
    df[c,2] <- strsplit(df[c,2], ";" )[[1]][1]
    df[c,3] <- tail(strsplit(df[c,2], ";" )[[1]],1)
  }
  return(df)
}

coord_PCG <- PCGdata[,1:3]
coord_D_cell <- subset(Lncdata[,1:3], rownames(Lncdata) %in% rownames(D_cell))
coord_D_day <- subset(Lncdata[,1:3], rownames(Lncdata) %in% rownames(D_day))

coord_PCG <- takeCoordinates(coord_PCG)
coord_D_cell <- takeCoordinates(coord_D_cell)
coord_D_day <- takeCoordinates(coord_D_day)


###################
data <- read.table(file = "./data/Fcounts_2Ddiff_PcodingGenes_s1.txt", header = TRUE)
rownames(data) <- data[,1]
dataclean <- subset(data[,28:48]) #take the sample columns


for (c in 1:21){ # rename the columns so its easier to read
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
                      cell_type=factor(rep(c("Dopamine","Dopamine","Dopamine","No_Dopamine","No_Dopamine","No_Dopamine","No_Dopamine"),3)))

pval=c(0.01)

dds_PCG_D <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~day)
dds_PCG_D<-DESeq(dds_PCG_D)
res_PCG_D <- results(dds_PCG_D, alpha = pval)
res_PCG_D<-res_PCG_D[order(res_PCG_D$padj),]
sig_PCG_D <- subset(res_PCG_D,padj<pval)
coord_sig_PCG_D <- subset(coord_PCG, rownames(coord_PCG) %in% rownames(sig_PCG_D))

dds_PCG_C <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~cell_type)
dds_PCG_C<-DESeq(dds_PCG_C)
res_PCG_C <- results(dds_PCG_C, alpha = pval)
res_PCG_C<-res_PCG_C[order(res_PCG_C$padj),]
sig_PCG_C <- subset(res_PCG_C,padj<pval)
coord_sig_PCG_C <- subset(coord_PCG, rownames(coord_PCG) %in% rownames(sig_PCG_C))


#Need to create the objects first
findOverlapPairs(coord_D_cell, coord_sig_PCG_C)
