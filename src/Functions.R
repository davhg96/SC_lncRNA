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

clasifyExp <- function(resdf, pval){
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


get_insert_info <- function(query_pos_df, subject_pos_df, query_dif_exp, subject_dif_exp, outputdir, assay,maxgap=5000, ignore_strand=TRUE){
  
  
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
  overlap <- findOverlaps(gr_query,gr_subject,type = "any", maxgap = maxgap, ignore.strand=ignore_strand ) #save the overlap

  #Info for the report
  total <- length(gr_query)
  intra_counts <- sum(countOverlaps(gr_query,gr_subject))
  extra_counts <- total - intra_counts
  
  #extract the Ids from the hits
  hit_ID <- data.frame(query_ID = gr_query$ID[queryHits(overlap)],
                       subject_ID = gr_subject$ID[subjectHits(overlap)])
  

  
  #Create the result DF
  result <- data.frame(query_ID=hit_ID$query_ID,
                       queryLFC=rep(1,nrow(hit_ID)),
                       queryBM=rep(1,nrow(hit_ID)),
                       queryStrand=rep(3,nrow(hit_ID)),
                       subject_ID=hit_ID$subject_ID,
                       subjectLFC=rep(3,nrow(hit_ID)),
                       subjectBM=rep(3,nrow(hit_ID)),
                       subjectStrand=rep(3,nrow(hit_ID)),
                       correlation=rep(3,nrow(hit_ID)),
                       correlation_text=rep(3,nrow(hit_ID)),
                       sense_rel =rep(3,nrow(hit_ID)))


  for (c in 1:nrow(result)){
    result$queryLFC[c] <- query_dif_exp[result$query_ID[c], "log2FoldChange"]
    result$queryBM[c] <- query_dif_exp[result$query_ID[c], "baseMean"]
    result$queryStrand[c] <- query_pos_df[result$query_ID[c],"Strand"]
    
    result$subjectLFC[c] <- subject_dif_exp[result$subject_ID[c], "log2FoldChange"]
    result$subjectBM[c] <- subject_dif_exp[result$subject_ID[c],"baseMean"]
    result$subjectStrand[c] <- subject_pos_df[result$subject_ID[c],"Strand"]
  }
  result <- na.omit(result)#Clean NA rows

  
  for (c in 1:nrow(result)){
    if(result$queryLFC[c] < 0 & result$subjectLFC[c] < 0 | result$queryLFC[c] > 0 & result$subjectLFC[c] > 0 ){
      result$correlation[c] = 1
      result$correlation_text[c]="Correlation"
    }
    else{
      result$correlation[c]=0
      result$correlation_text[c]="No Correlation"
    }
    
    if(result$queryStrand[c]==result$subjectStrand[c]){
      result$sense_rel[c] <- "Sense"
    }
    else{
      result$sense_rel[c] <- "Antisense"
    }
      
  }
  result$correlation_text <- as.factor(result$correlation_text)
  result$sense_rel <- as.factor(result$sense_rel)
 
  result$queryLFC <- as.numeric(result$queryLFC)
  result$subjectLFC <- as.numeric(result$subjectLFC)
return(result)
}
