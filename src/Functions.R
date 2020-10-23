library("DESeq2")
library("pheatmap")
library("tidyverse")
library("patchwork")
library("RColorBrewer")
library("openxlsx")
library("GenomicRanges")

cleanpcg <- function(LNCdata, PCGdata){
  PCGdata <-  PCGdata[! PCGdata$Geneid %in% LNCdata$Geneid,]#Clean all lncRNA from the list
  PCGdata$ic <- PCGdata$Geneid
  LNnames <- c()
  for(c in 1:nrow(PCGdata)){
    PCGdata$ic[c] <- strsplit(PCGdata$Geneid[c], "\\.")[[1]][1]
  }
  for(n in 1:nrow(LNCdata)){
    LNnames <- c(LNnames,strsplit(LNCdata$Geneid[n], "\\.")[[1]][1])
  }
  result <- PCGdata[! PCGdata$ic %in% LNnames,]
  result$ic <- NULL
  return(result)
}

cleanSampleNAmes <- function(datadf){
  for (c in 1:ncol(datadf)){ # rename the columns so its easier to read
    name <- strsplit(colnames(datadf[c]),"_")[[1]][4:5] #split the whole line and take the interesting parts
    p1<-strsplit(name[1], "\\.")[[1]][3]#select and clean part one (FGF and day)
    #f1 <- paste(p1[1],p1[2], sep=".")
    p2<-strsplit(name[2],"\\.")[[1]][1:2]#Select and clean the second part
    
    if(p2[2]=="bam"){#check and delete the bam extension
      f2<- p2[1]
      fname <- paste(p1, f2, sep=".")#joint the whole name and rename
      colnames(datadf)[c]<-fname
    }
    else{
      f2 <- paste(p2[1],p2[2], sep=".")
      fname <- paste(p1, f2, sep=".")#join the whole name and rename
      colnames(datadf)[c]<-fname
    }
  }
  return(datadf)
}#Used to clean the sample names

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
}# anotate the expression for MA plots


takeCoordinates <- function(df){
  for (c in 1: nrow(df)){
    df[c,1] <- strsplit(df[c,1], ";" )[[1]][1]#Chr
    df[c,2] <- strsplit(df[c,2], ";" )[[1]][1] #Start
    df[c,3] <- tail(strsplit(df[c,3], ";" )[[1]],1) #end
    df[c,4] <- strsplit(df[c,4], ";" )[[1]][1] #strand
  }
  return(df)
} 


get_insert_info <- function(query_pos_df, subject_pos_df, query_dif_exp, subject_dif_exp){
  
  #create the overlap objects
  gr_query <- makeGRangesFromDataFrame(query_pos_df,seqnames.field = "Chr", 
                                       start.field = "Start", 
                                       end.field = "End",
                                       strand.field= "Strand",
                                       keep.extra.columns = TRUE)
  gr_subject <- makeGRangesFromDataFrame(subject_pos_df,seqnames.field = "Chr", 
                                         start.field = "Start", 
                                         end.field = "End",
                                         strand.field= "Strand",
                                         keep.extra.columns = TRUE)
  overlap <- findOverlaps(gr_query,gr_subject,type = "within") #save the overlap only intragenic
  
  #Info for the report
  total <- length(gr_query)
  intra_counts <- sum(countOverlaps(gr_query,gr_subject))
  extra_counts <- total - intra_counts
  #extract the Ids from the hits
  query_hit_Id <- gr_query$ID[queryHits(overlap)]
  
  subject_hit_Id <- gr_subject$ID[subjectHits(overlap)]
  
  
  #Create the result placeholder
  result <- data.frame(query_ID=rep("ID", length(query_hit_Id)),
                       queryLFC=rep("qLFC", length(query_hit_Id)),
                       queryStrand=rep("qLFC", length(query_hit_Id)),
                       subject_ID=rep("ID", length(query_hit_Id)),
                       subjectLFC=rep("sLFC", length(query_hit_Id)))
  #populate the DF
  for (c in 1:length(query_hit_Id)){
    
    result$query_ID[c] <- query_hit_Id[c]
    result$queryLFC[c] <- query_dif_exp[query_hit_Id[c],"log2FoldChange"]
    result$queryStrand[c] <- query_pos_df[query_hit_Id[c],"Strand"]
    result$subject_ID[c] <- subject_hit_Id[c] 
    result$subjectLFC[c] <- subject_dif_exp[subject_hit_Id[c],"log2FoldChange"]
    
  }
  
  result <- na.omit(result)#Clean NA rows
  result$queryLFC <- as.numeric(result$queryLFC)
  result$subjectLFC <- as.numeric(result$subjectLFC)
  return(result)
}