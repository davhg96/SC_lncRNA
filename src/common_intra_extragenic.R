library("DESeq2")
library("tidyverse")
library("openxlsx")
library("RColorBrewer")
library("GenomicRanges")

outdir <- "./output/overlap"

D_cell <- read.xlsx("./output/Dop_NdopFGF+/pval_0.01/table/Top_significant_Cell.xlsx",rowNames = TRUE) #Diff expr
D_day <- read.xlsx("./output/Dop_NdopFGF+/pval_0.01/table/Top_significant_Day.xlsx",rowNames = TRUE) #Diff expr


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

coord_PCG <- PCGdata[,1:3]#All genes
coord_D_cell <- subset(Lncdata[,1:3], rownames(Lncdata) %in% rownames(D_cell)) #diff cell design
coord_D_day <- subset(Lncdata[,1:3], rownames(Lncdata) %in% rownames(D_day)) #diff day design

coord_PCG <- takeCoordinates(coord_PCG) #all pcg coords
coord_D_cell <- takeCoordinates(coord_D_cell) # coord diff cell
coord_D_cell <- rownames_to_column(coord_D_cell, var = "ID") # add the ids as a column
coord_D_day <- takeCoordinates(coord_D_day) # coord diff day
coord_D_day <- rownames_to_column(coord_D_day, var = "ID")# add the ids in a column

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

dds_PCG_D <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~day) # diff exp on genes day
dds_PCG_D<-DESeq(dds_PCG_D)
res_PCG_D <- results(dds_PCG_D, alpha = pval)
res_PCG_D<-res_PCG_D[order(res_PCG_D$padj),]
sig_PCG_D <- subset(res_PCG_D,padj<pval)
coord_sig_PCG_D <- subset(coord_PCG, rownames(coord_PCG) %in% rownames(sig_PCG_D)) #coord of sig exp by day
for(c in 1:nrow(coord_sig_PCG_D)){ coord_sig_PCG_D[c,1] <- paste0("chr",coord_sig_PCG_D[c,1]);rm(c)}
coord_sig_PCG_D <- rownames_to_column(coord_sig_PCG_D, var = "ID") # add the ids to a

dds_PCG_C <- DESeqDataSetFromMatrix(countData = dataclean, colData = coldata, design = ~cell_type) #diff exp genes by cell type
dds_PCG_C<-DESeq(dds_PCG_C)
res_PCG_C <- results(dds_PCG_C, alpha = pval)
res_PCG_C<-res_PCG_C[order(res_PCG_C$padj),]
sig_PCG_C <- subset(res_PCG_C,padj<pval)
coord_sig_PCG_C <- subset(coord_PCG, rownames(coord_PCG) %in% rownames(sig_PCG_C)) #coord of sig exp bu celltype
for(c in 1:nrow(coord_sig_PCG_C)){ coord_sig_PCG_C[c,1] <- paste0("chr",coord_sig_PCG_C[c,1] );rm(c)}
coord_sig_PCG_C <- rownames_to_column(coord_sig_PCG_C, var = "ID") #add ids to a column


get_insert_info <- function(query_pos_df, subject_pos_df, query_dif_exp, subject_dif_exp, outputdir, assay){
  outdir <- paste0(outputdir,"/insertionAnalysis/")
  dir.create(outdir, recursive = TRUE,showWarnings = FALSE)
  #create the overlap objects
  gr_query <- makeGRangesFromDataFrame(query_pos_df,seqnames.field = "Chr", start.field = "Start", end.field = "End",keep.extra.columns = TRUE)
  gr_subject <- makeGRangesFromDataFrame(subject_pos_df,seqnames.field = "Chr", start.field = "Start", end.field = "End",keep.extra.columns = TRUE)
  overlap <- findOverlaps(gr_query,gr_subject,type = "within") #save the overlap
  #Info for the report
  total <- length(gr_query)
  intra_counts <- sum(countOverlaps(gr_query,gr_subject))
  extra_counts <- total - intra_counts
  #extract the Ids from the hits
  query_hit_Id <- gr_query$ID[queryHits(overlap)]
  subject_hit_Id <- gr_subject$ID[subjectHits(overlap)]
  #extract the IDs from the expresssion data
  query_exp <- subset(query_dif_exp, rownames(query_dif_exp) %in% query_hit_Id)
  print(paste0("query Nrow", nrow(query_exp)))
  subject_exp <-  subset (subject_dif_exp, rownames(subject_dif_exp) %in% subject_hit_Id)
  print(paste0("subject Nrow",nrow(subject_exp)))
  #Create the result DF
  result <- data.frame(query_ID=rownames(query_exp),
                       queryLFC=query_exp$log2FoldChange,
                       subject_ID=rownames(subject_exp),
                       subjectLFC=subject_exp$log2FoldChange,
                       correlation=rep(3,nrow(query_exp)),
                       correlation_text=rep(4,nrow(query_exp)))
  for (c in 1:nrow(result)){
    if(result$queryLFC[c] < 0 & result$subjectLFC[c] < 0 | result$queryLFC[c] > 0 & result$subjectLFC[c] > 0 ){
      result$correlation[c]=1
      result$correlation_text[c]="Correlation"
    }
    else{
      result$correlation[c]=0
      result$correlation_text[c]="No Correlation"
    }
  }
  result$correlation_text <- as.factor(result$correlation_text)
  sink(paste0(outdir,"summary",assay,".txt"))
  cat("Total lncRNA", "\t", total,"\n",
  "extragenic counts", "\t", extra_counts, "\n",
  "Intragenic Counts", "\t", intra_counts, "\n",
  "Number of inseritions with correlated counts","\t", sum(result$correlation), "\n",
  "Out of:","\t",nrow(result))
  sink()
  
  ggplot(result)+
    geom_bar(aes(x=correlation_text, y = ((..count..)/sum(..count..)*100),fill=correlation_text))+
    xlab("% of correlation of lncRNAs and genes they are inserted in")+
    ylab("%")+ylim(0,100)+
    labs(fill="Legend")+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  ggsave(paste0("barplot",assay,".pdf"),device = "pdf",path = outdir, width = 9, height = 16, units = "cm")  
  
  return(result)
}



####overlaps cell design
celloverlap <- get_insert_info(query_pos_df = coord_D_cell, subject_pos_df = coord_sig_PCG_C, query_dif_exp = D_cell,subject_dif_exp = sig_PCG_C, outputdir = outdir, assay = "cell")



  
# gr_D_Cell <- makeGRangesFromDataFrame(coord_D_cell,seqnames.field = "Chr", start.field = "Start", end.field = "End",keep.extra.columns = TRUE) #create objects
# gr_PCG_C <- makeGRangesFromDataFrame(coord_sig_PCG_C,seqnames.field = "Chr", start.field = "Start", end.field = "End",keep.extra.columns = TRUE) #create objects
# 
# Overlap_cell <- findOverlaps(gr_D_Cell,gr_PCG_C,type = "within") #save the overlap
# sum(countOverlaps(gr_D_Cell,gr_PCG_C)) #90 hits



#Overlaps day design
dayoverlap <- get_insert_info(query_pos_df = coord_D_day, subject_pos_df = coord_sig_PCG_D, query_dif_exp = D_day,subject_dif_exp = sig_PCG_D, outputdir = outdir, assay = "day")

# gr_D_day <- makeGRangesFromDataFrame(coord_D_day,seqnames.field = "Chr", start.field = "Start", end.field = "End",keep.extra.columns = TRUE) #create objects
# gr_PCG_D <- makeGRangesFromDataFrame(coord_sig_PCG_D,seqnames.field = "Chr", start.field = "Start", end.field = "End",keep.extra.columns = TRUE) #create objects
# 
# Overlap_day <- findOverlaps(gr_D_day,gr_PCG_D,type = "within") #save the overlapping
# sum(countOverlaps(gr_D_day,gr_PCG_D)) #260 hits
# 

#Next: extract the intragenic lncRNA and their respective PCG, and barplot % where the correlation is true or false