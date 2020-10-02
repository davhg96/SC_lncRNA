library(Seurat)
library(tidyverse)
library(patchwork)

data <- readRDS("./data/2ddiff.fgf8.d16d30d60.seurat.rds")

GetGeneList <- function(directory, header=FALSE){
  FileList <- list.files(directory, pattern = ".txt") # ls the directory
  ResList <- list() #placeholder
  for (file in FileList){
    name<-str_remove(file,pattern = ".txt")
    path <- paste0(directory,file)
    Filedata <- as.list(read.delim(path, sep = "\n", header = header)) #Extract the info
    ResList[[name]] <- Filedata[[1]] #Add the info to the result object
    
    
  }
  
  return(ResList)
}

featuresRNA<-GetGeneList("./data/Features/lncRNA/")

dir.create("./output/FeaturePlots", recursive = TRUE)
outdir<-"./output/FeaturePlots"

for(name in names(featuresRNA)){
  img1<-DotPlot(data, assay = "RNA", features = featuresRNA[[name]]) + coord_flip()
  ggsave(filename=paste0("DotPlot_",name,".pdf"),plot = img1,path = outdir, device = "pdf" )
  

  img2<-VlnPlot(data,assay = "RNA", features =  featuresRNA[[name]] ,pt.size = 0, ncol = 2)
  ggsave(filename=paste0("VLNPlot_",name,"_","gene",".pdf"),plot = img2,path = outdir, device = "pdf", height = 40, width = 20, units = "cm" )
}



featuresCellTypes<-GetGeneList("./data/Features/CellTypes/", header = TRUE)
