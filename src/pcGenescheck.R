library("DESeq2")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")
library("openxlsx")


dataPCG <- read.table(file = "./data/Fcounts_2Ddiff_PcodingGenes_s1.txt", header = TRUE)
rownames(dataPCG) <- dataPCG[,1]
dataRNA <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
rownames(dataRNA) <- dataRNA[,1]


PCGclean <- subset(dataPCG[,7:48]) #take the sample columns
RNAclean <- subset(dataRNA[,7:48]) #take the sample columns
rm(dataPCG,dataRNA)

cleanDF <- function(DF){
for (c in 1:42){ # rename the columns so its easier to read
  name <- strsplit(colnames(DF[c]),"_")[[1]][4:5] #split the whole line and take the interesting parts
  p1<-strsplit(name[1], "\\.")[[1]][2:3]#select and clean part one (FGF and day)
  f1 <- paste(p1[1],p1[2], sep=".")
  p2<-strsplit(name[2],"\\.")[[1]][1:2]#Select and clean the second part
  
  if(p2[2]=="bam"){#check and delete the bam extension
    f2<- p2[1]
    fname <- paste(f1, f2, sep=".")#joint the whole name and rename
    colnames(DF)[c]<-fname
  }
  else{
    f2 <- paste(p2[1],p2[2], sep=".")
    fname <- paste(f1, f2, sep=".")#join the whole name and rename
    colnames(DF)[c]<-fname
  }
}
rm(f1,f2,p1,p2,name,fname,c) #clean
return(DF)
}

PCGclean <- cleanDF(PCGclean)
RNAclean <- cleanDF(RNAclean)

common<-merge(RNAclean, PCGclean, by="row.names")
common<-column_to_rownames(common,var = "Row.names")

#MArkerBarplot
####### 
markers<- c("MEG3","MKI67","TH")
markerDF <- subset(PCGclean, rownames(PCGclean)%in% markers)
meanmarkers<- data.frame(row.names = rownames(markerDF))
for (c in 1:21){
  name<-paste(strsplit(colnames(markerDF[c]),"\\.")[[1]][2:4],sep = ".")
  cname<-paste(name[1],name[2],name[3], sep = ".")
  columns<- c(c,c+21)
  meanmarkers[[cname]] <- rowMeans(markerDF[,columns])
}

meanmarkers<- rownames_to_column(meanmarkers, "marker")
meanmarkers<- gather(meanmarkers,cell, counts, colnames(meanmarkers[,2:22] ))

plot<- ggplot(meanmarkers,aes(x=cell, y=counts, fill=cell)) +
  geom_bar(position = "dodge", stat="identity") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 10)) +
  xlab("Sample") +
  ylab("Counts") +
  facet_wrap(~marker)

ggsave("Barplot_Sample.pdf", plot = plot, device = "pdf", path = "./output")



#Boxplot

#####

CType <- c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late","VLMC")

RNAmean<-data.frame(row.names = rownames(RNAclean))

for (c in 1:21){
  name<-paste(strsplit(colnames(RNAclean[c]),"\\.")[[1]][2:4],sep = ".")
  cname<-paste(name[1],name[2],name[3], sep = ".")
  columns<- c(c,c+21)
  RNAmean[[cname]] <- rowMeans(RNAclean[,columns])
}

RNABox<-gather(RNAmean,cell,counts,c(1:21))
RNABox$day<-RNABox$cell
RNABox<-filter(RNABox,counts!=0 &counts<1000)
for(c in 1:nrow(RNABox)){#This is pure bruteforce
  cell<-paste(strsplit(RNABox$cell[c],"\\.")[[1]][2],strsplit(RNABox$cell[c],"\\.")[[1]][3],sep = ".")
  day<-strsplit(RNABox$cell[c],"\\.")[[1]][1]
  RNABox$cell[c]<-cell
  RNABox$day[c]<-day
}

BX<-ggplot(RNABox, aes(x=cell, y=counts, fill=day)) + 
  geom_boxplot() +
  facet_wrap(~cell, scale="free")
ggsave("boxplot.pdf", plot = BX, device = "pdf", path = "./output")




RNAsum<-data.frame(day=factor(c(rep("day_16",7),rep("day_30",7),rep("day_60",7))), 
                  cell_type=factor(rep(c("DA.1","DA.2","DA.E1","FP.cycling","FP.early","FP.late","VLMC"),3)),
                  counts=colSums(RNAmean))

BP_C<-ggplot(RNAsum,aes(x=day, y=counts, fill=day)) +
  geom_bar(position = "dodge", stat="identity") +
  xlab("Sample") +
  ylab("Counts") +
  facet_wrap(~cell_type)
ggsave("Barplot_cell.pdf", plot = BP_C, device = "pdf", path = "./output")
