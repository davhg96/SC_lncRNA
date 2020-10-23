
source("./src/Functions.R") #Useful functions + packages

outputdir <- ("./output/FinalPlots/pval_0.01/")
dir.create(outputdir, recursive = TRUE,showWarnings = FALSE)
pval=c(0.01)

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

expressedLNC <- subset(res_LNCd_60_16, res_LNCd_60_16$baseMean>1)


dds_PCGd <- DESeqDataSetFromMatrix(countData = PCGclean, colData = coldata, design = ~day)
dds_PCGd<-DESeq(dds_PCGd)
res_PCGd <- results(dds_PCGd, alpha = pval,contrast = c("day", "Day 60","Day 16"))
res_PCGd<-res_PCGd[order(res_PCGd$padj),]
res_PCGd <- as.data.frame(res_PCGd)

expressedPCG <- subset(res_PCGd, res_PCGd$baseMean>1)


# cell design ----

dds_LNCc <- DESeqDataSetFromMatrix(countData = LNCclean, colData = coldata, design = ~cell_type)
dds_LNCc<-DESeq(dds_LNCc)
res_LNCc <- results(dds_LNCc, alpha = pval,contrast = c("cell_type","Dopaminergic neurons", "Floorplate"))
res_LNCc<-res_c[order(res_LNCc$padj),]
res_LNCc <- as.data.frame(res_LNCc)


dds_PCGc <- DESeqDataSetFromMatrix(countData = PCGclean, colData = coldata, design = ~cell_type)
dds_PCGc<-DESeq(dds_PCGc)
res_PCGc <- results(dds_PCGc, alpha = pval,contrast = c("cell_type","Dopaminergic neurons", "Floorplate"))
res_PCGc<-res_PCGc[order(res_PCGc$padj),]
res_PCGc <- as.data.frame(res_PCGc)





#Nearby correlation ----
outdir <- paste0(outputdir,"Nearby_correlation/")
dir.create(outdir, recursive = TRUE,showWarnings = FALSE)



celloverlap <- get_insert_info(query_pos_df = coordLNC, 
                               subject_pos_df = coordPCG, 
                               query_dif_exp = res_LNCc,
                               subject_dif_exp = res_PCGc)

celloverlap <- subset(celloverlap, celloverlap$query_ID %in% rownames(expressedLNC))

corr <- cor(celloverlap$queryLFC, celloverlap$subjectLFC, method = "pearson")

ggplot()+
  geom_point(data=celloverlap,aes(x=queryLFC,
                                  y=subjectLFC),size=0.3, color="#0914e0")+
  theme_minimal()+
  scale_x_continuous(name="lncRNA log fold change", limits=c(-5, 5)) +
  scale_y_continuous(name="PC gene log fold change", limits=c(-5, 5)) +
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  labs(title=paste0("Nearby expression Pearson corr=", corr))
ggsave(filename =paste0("Nearby expression cell.pdf"), device="pdf",path = outdir, width = 16,height = 9,units = "cm" )


dayoverlap <- get_insert_info(query_pos_df = coordLNC, 
                              subject_pos_df = coordPCG, 
                              query_dif_exp = res_LNCd_60_16,
                              subject_dif_exp = res_PCGd)

dayoverlap <- subset(dayoverlap, dayoverlap$query_ID %in% rownames(expressedLNC))
corr <- cor(dayoverlap$queryLFC, dayoverlap$subjectLFC, method = "pearson")

ggplot()+
  geom_point(data=dayoverlap,aes(x=queryLFC,
                                  y=subjectLFC),size=0.3,color="#0914e0")+
  theme_minimal()+
  scale_x_continuous(name="lncRNA log fold change", limits=c(-5, 5)) +
  scale_y_continuous(name="PC gene log fold change", limits=c(-5, 5)) +
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  labs(title=paste0("Nearby expression Pearson corr=", corr))
ggsave(filename =paste0("Nearby expression day.pdf"), device="pdf",path = outdir, width = 16,height = 9,units = "cm" )

# strand intra inter gen analysis ----
#Detected
intragenic <- subset(dayoverlap, dayoverlap$query_ID %in% rownames(expressedLNC))#get the expressed ones
intragenic <- cbind(table(intragenic$queryStrand), prop.table(table(intragenic$queryStrand)))
intergenic <- coordLNC[! rownames(coordLNC) %in% dayoverlap$query_ID,]#take em all
intergenic <- subset(intergenic, rownames(intergenic) %in% rownames(expressedLNC) )#Get the expressed ones
intergenic <- cbind(table(intergenic$Strand), prop.table(table(intergenic$Strand)))

plotdf <- data.frame(sample=factor(c("Intragenic Sense","Intragenic Antisense","Intergenic Sense"  ,"Intergenic Antisense"), levels =c("Intragenic Sense","Intragenic Antisense","Intergenic Sense"  ,"Intergenic Antisense") ))
plotdf <- cbind(plotdf, rbind(intragenic[2:1,], intergenic[2:1,]))
colnames(plotdf) <- c("Sample","Counts","Perc")


ggplot(data = plotdf, aes(x="", y=Counts, fill=Sample))+
  geom_col(width = 1,color="white", position = "stack")  +
  coord_polar(theta = "y") +
  scale_fill_manual(name="",
                    values=c("#35de62","#009c29","#06ccd6","#074ca6"))+
  geom_text(aes(label =Counts), 
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


ggsave("Count_distribution-Detected-PCG-LNCRNAInter-Intragenic.pdf",device = "pdf", dpi = 600, path = outdir )

#All
intragenic <- cbind(table(dayoverlap$queryStrand), prop.table(table(dayoverlap$queryStrand)))
intergenic <- coordLNC[! rownames(coordLNC) %in% dayoverlap$query_ID,]
intergenic <- cbind(table(intergenic$Strand), prop.table(table(intergenic$Strand)))


plotdf <- data.frame(sample=factor(c("Intragenic Sense","Intragenic Antisense","Intergenic Sense"  ,"Intergenic Antisense"), levels =c("Intragenic Sense","Intragenic Antisense","Intergenic Sense"  ,"Intergenic Antisense") ))
plotdf <- cbind(plotdf, rbind(intragenic[2:1,], intergenic[2:1,]))
colnames(plotdf) <- c("Sample","Counts","Perc")


ggplot(data = plotdf, aes(x="", y=Counts, fill=Sample))+
  geom_col(width = 1,color="white", position = "stack")  +
  coord_polar(theta = "y") +
  scale_fill_manual(name="",
                    values=c("#35de62","#009c29","#06ccd6","#074ca6"))+
  geom_text(aes(label =Counts), 
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


ggsave("Count_distributionPCG-LNCRNAInter-Intragenic.pdf",device = "pdf", dpi = 600, path = outdir )

# piechart distribution, Strand correlation ----
outdir <- outputdir
totalLNCcounts <- nrow(expressedLNC)
totalPCGcounts <-  nrow(expressedPCG)
totalcounts <- sum(totalLNCcounts,totalPCGcounts)

count_df <- data.frame(sample=c("lncRNA","Protein coding genes"),
                      counts=c(totalLNCcounts,totalPCGcounts),
                      perc=c((totalLNCcounts/totalcounts)*100,(totalPCGcounts/totalcounts)*100 ))

coordLNC <- LNCdata[,2:5]

coordLNC <- takeCoordinates(coordLNC)
coordLNC <- subset(coordLNC,rownames(coordLNC) %in% rownames(expressedLNC))
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
lengths <- subset(LNCdata, rownames(LNCdata) %in%  rownames(expressedLNC))
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

# Common expression -------------------------------------------------------
outdirT <- paste0(outputdir,"Table/")
dir.create(outdirT, recursive = TRUE, showWarnings = FALSE)

common <- data.frame(ID=rownames(res_LNCd_60_16),
                     dayLFC=subset(res_LNCd_60_16$log2FoldChange, rownames(res_LNCd_60_16) %in% rownames(res_LNCc)),
                     dayPADJ=subset(res_LNCd_60_16$padj, rownames(res_LNCd_60_16) %in% rownames(res_LNCc)),
                     cellLFC=subset(res_LNCc$log2FoldChange, rownames(res_LNCd_60_16) %in% rownames(res_LNCc)),
                     cellPADJ=subset(res_LNCc$padj, rownames(res_LNCd_60_16) %in% rownames(res_LNCc)))
common <- na.omit(common)
write.xlsx(common,file =paste0(outdirT, "Common_NoFilter_DayVsCell.xlsx"))

res_c_sig <- subset(res_LNCc, res_LNCc$padj<pval)
res_d_60_16_sig <- subset(res_LNCd_60_16, res_LNCd_60_16$padj<pval)
res_d_60_16_sig <- res_d_60_16_sig[rownames(res_c_sig),]
rownames(subset(res_d_60_16_sig, rownames(res_d_60_16_sig)%in%rownames(res_c_sig)))

common_sig <- data.frame(ID=rownames(subset(res_d_60_16_sig, rownames(res_d_60_16_sig)%in%rownames(res_c_sig))),
                         dayLFC=subset(res_d_60_16_sig$log2FoldChange, rownames(res_d_60_16_sig) %in% rownames(res_c_sig)),
                         dayPADJ=subset(res_d_60_16_sig$padj, rownames(res_d_60_16_sig) %in% rownames(res_c_sig)),
                         cellLFC=subset(res_c_sig$log2FoldChange, rownames(res_d_60_16_sig) %in% rownames(res_c_sig)),
                         cellPADJ=subset(res_c_sig$padj, rownames(res_d_60_16_sig) %in% rownames(res_c_sig)))
common_sig <- na.omit(common_sig)
write.xlsx(common_sig,file =paste0(outdirT, "Common_Significant_DayVsCell.xlsx"))


rm(coldata,pval, outdirT,res_LNCc,res_c_sig ,res_LNCd_60_16,res_d_60_16_sig, common, common_sig)


