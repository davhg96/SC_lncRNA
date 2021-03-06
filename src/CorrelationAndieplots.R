{
source("./src/Functions.R") #Useful functions + packages
outputdir <- ("./output/reproduce/pval_0.01/")
dir.create(outputdir, recursive = TRUE,showWarnings = FALSE)
pval=c(0.01)
}

LNCdata <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
rownames(LNCdata) <- LNCdata[,1]
LNCclean <- subset(LNCdata[,c(28:33,35:40,42:47)]) #take everything but fgf- and vlmc
LNCclean <- cleanSampleNAmes(LNCclean)

PCGdata <- read.table(file = "./data/Fcounts_2Ddiff_PcodingGenes_s1.txt", header = TRUE)
PCGdata <- cleanpcg(LNCdata,PCGdata)
rownames(PCGdata) <- PCGdata[,1]
PCGclean <- subset(PCGdata[,c(28:33,35:40,42:47)]) #take everything but fgf- and vlmc
PCGclean <- cleanSampleNAmes(PCGclean)


problematic <- c(16114,16115)


coordLNC <- takeCoordinates(LNCdata[,2:5])
coordLNC$ID <- rownames(coordLNC)
coordLNC <- coordLNC[-problematic,]

coordPCG <- takeCoordinates(PCGdata[,2:5])

for(c in 1:nrow(coordPCG)){ coordPCG[c,1] <- paste0("chr",coordPCG[c,1] );rm(c)} #make the chr names match
coordPCG$ID <- rownames(coordPCG)

coldata <- data.frame(day=factor(c(rep("day_16",6),rep("day_30",6),rep("day_60",6)),
                                 labels  = c("Day 16","Day 30","Day 60")), 
                      cell_type=factor(rep(c("No_Dopamine","No_Dopamine","No_Dopamine","Dopamine","Dopamine","Dopamine"),3),
                                       labels   = c("Dopaminergic neurons","Floorplate")))


# day design ----
dds_LNCd <- DESeqDataSetFromMatrix(countData = LNCclean, colData = coldata, design = ~day)
dds_LNCd<-DESeq(dds_LNCd)
res_LNCd_60_16 <- results(dds_LNCd, alpha = pval, contrast = c("day", "Day 60","Day 16"))
res_LNCd_60_16<-res_LNCd_60_16[order(res_LNCd_60_16$padj),]
res_LNCd_60_16 <- as.data.frame(res_LNCd_60_16)
expressedLNCd <- subset(res_LNCd_60_16, res_LNCd_60_16$baseMean>1)

dds_PCGd <- DESeqDataSetFromMatrix(countData = PCGclean, colData = coldata, design = ~day)
dds_PCGd<-DESeq(dds_PCGd)
res_PCGd <- results(dds_PCGd, alpha = pval,contrast = c("day", "Day 60","Day 16"))
res_PCGd<-res_PCGd[order(res_PCGd$padj),]
res_PCGd <- as.data.frame(res_PCGd)
expressedPCGd <- subset(res_PCGd, res_PCGd$baseMean>1)


# cell design ----

dds_LNCc <- DESeqDataSetFromMatrix(countData = LNCclean, colData = coldata, design = ~cell_type)
dds_LNCc<-DESeq(dds_LNCc)
res_LNCc <- results(dds_LNCc, alpha = pval,contrast = c("cell_type","Dopaminergic neurons", "Floorplate"))
res_LNCc <- as.data.frame(res_LNCc)


dds_PCGc <- DESeqDataSetFromMatrix(countData = PCGclean, colData = coldata, design = ~cell_type)
dds_PCGc<-DESeq(dds_PCGc)
res_PCGc <- results(dds_PCGc, alpha = pval,contrast = c("cell_type","Dopaminergic neurons", "Floorplate"))
res_PCGc <- as.data.frame(res_PCGc)





#Nearby correlation ----
outdir <- paste0(outputdir,"Nearby_correlation/")
dir.create(outdir, recursive = TRUE,showWarnings = FALSE)
res_LNCc <- subset(res_LNCc, res_LNCc$pvalue<pval)

celloverlap <- get_insert_info(query_pos_df = coordLNC, 
                               subject_pos_df = coordPCG, 
                               query_dif_exp = res_LNCc,
                               subject_dif_exp = res_PCGc,
                               outputdir = outputdir,
                               assay = "cell",
                               ignore_strand = TRUE)




corr <- round(cor(celloverlap$queryLFC, celloverlap$subjectLFC, method = "pearson"), 2)

ggplot()+
  geom_point(data=celloverlap,aes(x=queryLFC,
                                  y=subjectLFC),size=0.3, color="#000000")+
  theme_minimal()+
  scale_x_continuous(name="lncRNA log fold change", limits=c(-5, 5)) +
  scale_y_continuous(name="PC gene log fold change", limits=c(-5, 5)) +
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  labs(title=paste0("Nearby expression Pearson corr=", corr))
ggsave(filename =paste0("Nearby expression cell.pdf"), device="pdf",path = outdir, width = 16,height = 9,units = "cm" )




# intra intergenic----

dayoverlap <- get_insert_info(query_pos_df = coordLNC, 
                               subject_pos_df = coordPCG, 
                               query_dif_exp = expressedLNCd,
                               subject_dif_exp = res_PCGd,
                               outputdir = outputdir,
                               assay = "cell",
                               ignore_strand = TRUE)


intragenic <- subset(expressedLNCd, rownames(expressedLNCd) %in% unique(dayoverlap$query_ID))#get the expressed ones
intragenic <- nrow(intragenic)
intergenic <- coordLNC[! rownames(coordLNC) %in% dayoverlap$query_ID,]#take em all
intergenic <- subset(intergenic, rownames(intergenic) %in% rownames(expressedLNCd) )#Get the expressed ones
intergenic <- nrow(intergenic)

plotdf <- data.frame(row.names =c("Intragenic","Intergenic"),
                     Count=c(intragenic, intergenic))

plotdf <- cbind(plotdf,prop.table(plotdf))
plotdf <- rownames_to_column(plotdf, var = "Sample")
colnames(plotdf) <- c("Sample","Counts","Perc")


ggplot(data = plotdf, aes(x="", y=Counts, fill=Sample))+
  geom_col(width = 1,color="white", position = "stack")  +
  coord_polar(theta = "y") +
  scale_fill_manual(name="",
                    breaks = c("Intragenic", "Intergenic"),
                    labels = paste(plotdf$Sample,"[",round(plotdf$Perc, 4)*100, "%]"),
                    values=c("#83a9a7","#cbcbcb"))+
  # geom_text(aes(label =Counts), 
  #           position = position_stack(vjust = 0.5)) +
  labs(title = "Count distribution") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, color = "#000000"),
        legend.title = element_text(size = rel(1.1)),
        legend.text=element_text(size=rel(1)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())# remove background, grid, numeric labels


ggsave("Intragenic_Intergenic.pdf",device = "pdf", dpi = 600, path = outdir )



# piechart distribution, Strand correlation ----
outdir <- outputdir
totalLNCcounts <- nrow(expressedLNCd)
totalPCGcounts <-  nrow(expressedPCGd)
totalcounts <- sum(totalLNCcounts,totalPCGcounts)

count_df <- data.frame(sample=c("lncRNA","Protein coding genes"),
                      counts=c(totalLNCcounts,totalPCGcounts),
                      perc=c((totalLNCcounts/totalcounts)*100,(totalPCGcounts/totalcounts)*100 ))

coordLNC <- LNCdata[,2:5]

coordLNC <- takeCoordinates(coordLNC)
coordLNC <- subset(coordLNC,rownames(coordLNC) %in% rownames(expressedLNCd))
count_df1 <- data.frame(strand=c("Antisense","Sense"),
                     counts=table(coordLNC$Strand),
                     perc=prop.table(table(coordLNC$Strand))*100)
count_df1 <- count_df1[2:1,c(1,3,5)]
colnames(count_df1) <- colnames(count_df)



plotdf <- data.frame(side=factor(c(1,1,2,2), 
                                 labels = c("Total count distribution","lncRNA strand distribution")),
                     sample=factor(c(1,2,3,4),
                     labels=c("Protein Coding genes","lncRNA","sense lncRNA","Antisense lncRNA")))

plotdf <- cbind(plotdf,rbind(count_df[2:1,2:3], count_df1[,2:3]))



ggplot(data = plotdf, aes(x="", y=perc, fill=sample))+
  geom_col(width = 1,color="white", position = "stack")  +
  coord_polar(theta = "y") +
  scale_fill_manual(name="",
                   values=c("#f4a53e","#09b328","#06ccd6","#074ca6"))+
  facet_wrap(~side)+
  geom_text(aes(label =counts), 
            position = position_stack(vjust = 0.5)) +
  labs(title = "Count distribution") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, color = "#000000"),
        legend.title = element_text(size = rel(1.1)),
        legend.text=element_text(size=rel(1)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())# remove background, grid, numeric labels


ggsave("Count_distributionPCG-LNCRNA.pdf",device = "pdf", dpi = 600, path = outdir )

rm(count_df, coordLNC, plotdf, count_df1,totalLNCcounts,totalcounts, totalPCGcounts)


# length distribution (detected) ----
lengths <- subset(LNCdata, rownames(LNCdata) %in%  rownames(expressedLNCd))
lengths <- lengths[,c("Geneid","Length")]
lengths <- lengths[order(lengths$Length),]

lengths %>% 
  filter(Length<25000) %>% 
ggplot(  aes(x=Length)) +
  geom_histogram( binwidth=100, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("lncRNA length distribution lncRNAs < 25Kb") +
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  theme_minimal() +
  theme(
    plot.title = element_text(size=15)
  )

ggsave("Length distribution-DETECTED.pdf",device = "pdf", dpi = 600, path = outputdir )




rm(list=setdiff(ls(), c("pval", "outputdir")))


