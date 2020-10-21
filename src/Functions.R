library("DESeq2")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")
library("openxlsx")
library("GenomicRanges")

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