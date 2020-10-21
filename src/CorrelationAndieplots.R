
source("./src/Functions.R") #Useful functions + packages

outdir <- ("./output/FinalPlots/pval_0.01/")
dir.create(outdir, recursive = TRUE,showWarnings = FALSE)

LNCdata <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)
rownames(LNCdata) <- LNCdata[,1]
LNCclean <- subset(LNCdata[,c(28:33,35:40,42:47)]) #take everything but fgf- and vlmc
LNCclean <- cleanSampleNAmes(LNCclean)

PCGdata <- read.table(file = "./data/Fcounts_2Ddiff_PcodingGenes_s1.txt", header = TRUE)
PCGdata <-  PCGdata[! PCGdata$Geneid %in% LNCdata$Geneid,]#Clean all lncRNA from the list
rownames(PCGdata) <- PCGdata[,1]
PCGclean <- subset(PCGdata[,c(28:33,35:40,42:47)]) #take everything but fgf- and vlmc
PCGclean <- cleanSampleNAmes(PCGclean)




# piechart distribution ---------------------------------------------------


totalLNCcounts <- sum(colSums(LNCclean))
totalPCGcounts <- sum(colSums(PCGclean))
totalcounts <- sum(totalLNCcounts,totalPCGcounts)

count_df <- data.frame(sample=c("lncRNA","Protein coding genes"),
                      counts=c(totalLNCcounts,totalPCGcounts),
                      perc=c((totalLNCcounts/totalcounts)*100,(totalPCGcounts/totalcounts)*100 ))

rm(totalLNCcounts,totalcounts, totalPCGcounts)

ggplot(data = count_df, aes(x="", y=perc, fill=sample))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y") +
  scale_fill_manual(name="Legend", 
                    values=c("#F4D03F","#82E0AA"))+
                    # labels=paste0(count_df$sample," [", count_df$counts,"]"))+
  geom_text(aes(label = paste0(round(perc), "%")), 
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
rm(count_df)



# Common expression -------------------------------------------------------
coldata <- data.frame(day=factor(c(rep("day_16",6),rep("day_30",6),rep("day_60",6)),
                                 labels  = c("Day 16","Day 30","Day 60")), 
                      cell_type=factor(rep(c("No_Dopamine","No_Dopamine","No_Dopamine","Dopamine","Dopamine","Dopamine"),3),
                                       labels   = c("Dopaminergic neurons","Floorplate")))

pval=c(0.01)

outdirT<-paste0(outdir,"Table/")
dir.create(outdirT,recursive = TRUE,showWarnings = FALSE )

dds_d <- DESeqDataSetFromMatrix(countData = LNCclean, colData = coldata, design = ~day)
dds_d<-DESeq(dds_d)
res_d_60_16 <- results(dds_d, alpha = pval, contrast = c("day", "Day 60","Day 16"))
res_d_60_16<-res_d_60_16[order(res_d_60_16$padj),]
res_d_60_16 <- as.data.frame(res_d_60_16)

dds_c <- DESeqDataSetFromMatrix(countData = LNCclean, colData = coldata, design = ~cell_type)
dds_c<-DESeq(dds_c)
res_c <- results(dds_c, alpha = pval,contrast = c("cell_type","Dopaminergic neurons", "Floorplate"))
res_c<-res_c[order(res_c$padj),]
res_c <- as.data.frame(res_c)
res_c <- res_c[rownames(res_d_60_16),]#make the df have the same order

common <- data.frame(ID=rownames(res_d_60_16),
                     dayLFC=subset(res_d_60_16$log2FoldChange, rownames(res_d_60_16) %in% rownames(res_c)),
                     dayPADJ=subset(res_d_60_16$padj, rownames(res_d_60_16) %in% rownames(res_c)),
                     cellLFC=subset(res_c$log2FoldChange, rownames(res_d_60_16) %in% rownames(res_c)),
                     cellPADJ=subset(res_c$padj, rownames(res_d_60_16) %in% rownames(res_c)))
write.xlsx(common,file =paste0(outdirT, "Common_NoFilter_DayVsCell.xlsx"))

res_c_sig <- subset(res_c, res_c$padj<pval)
res_d_60_16_sig <- subset(res_d_60_16, res_d_60_16$padj<pval)
res_d_60_16_sig <- res_d_60_16_sig[rownames(res_c_sig),]
rownames(subset(res_d_60_16_sig, rownames(res_d_60_16_sig)%in%rownames(res_c_sig)))

common_sig <- data.frame(ID=rownames(subset(res_d_60_16_sig, rownames(res_d_60_16_sig)%in%rownames(res_c_sig))),
                         dayLFC=subset(res_d_60_16_sig$log2FoldChange, rownames(res_d_60_16_sig) %in% rownames(res_c_sig)),
                         dayPADJ=subset(res_d_60_16_sig$padj, rownames(res_d_60_16_sig) %in% rownames(res_c_sig)),
                         cellLFC=subset(res_c_sig$log2FoldChange, rownames(res_d_60_16_sig) %in% rownames(res_c_sig)),
                         cellPADJ=subset(res_c_sig$padj, rownames(res_d_60_16_sig) %in% rownames(res_c_sig)))

write.xlsx(common_sig,file =paste0(outdirT, "Common_Significant_DayVsCell.xlsx"))


rm(coldata,pval, outdirT, dds_d,dds_c,res_c,res_c_sig,res_d_60_16,res_d_60_16_sig, common, common_sig)


# Strand correlation ----


#Nearby correlation ----
