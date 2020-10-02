plotMAwithLabels <- function(res, threshold, maxNum=20) {
  res_df <- res %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj < threshold) %>%
    dplyr::mutate(regulation = as.factor( - (2*(log2FoldChange > 0) - 1)*significant))
  res_df$regulation[is.na(res_df$regulation)] <- “0”
  res_labelled <- res_df %>%
    dplyr::filter(significant)
  flog.info(“res_labelled: %s”, paste(dim(res_labelled), collapse=‘x’))
  p <- ggplot(res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=regulation)) +
    geom_point()
  p <- p + scale_colour_manual(values=c(“1”=“blue”, “-1”=“red”, “0”=“black”, “NA”=“lightgray”),
                               labels=c(“1”=“DOWN”, “-1”=“UP”, “0”=“Not significant”))
  if(dim(res_labelled)[1] < maxNum) {
    p <- p + ggrepel::geom_label_repel(data=res_labelled, aes(label=Gene))
  }
  p
}



suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))
suppressMessages(require(pheatmap))
suppressMessages(require(stringr))
suppressMessages(require(knitr))
suppressMessages(require(data.table))
suppressMessages(require(reticulate))
suppressMessages(require(futile.logger))
# Set up outDir
dataDir <- ‘/Users/yo4230sh/DropboxDNLP/2D_Diff_scRNAseq/scRNASeq_gene_level/raw_data/’
out_dir <- file.path(‘/Users/yo4230sh/Desktop/’, paste0(“LncRNA_Featureplot”, Sys.Date()))
dir.create(out_dir, recursive = T, showWarnings = F)
inputRDSfile <- paste0(dataDir, ‘2ddiff.fgf8.d16d30d60.seurat.rds’)
rdsData <- readRDS(inputRDSfile)
rdsData@active.assay <-“RNA”
plot_z <- FeaturePlot(rdsData, features = c(“MEG3”,“DLK1"),
            reduction = “umap”, cols = c(“yellow”,“red”,“black”), pt.size = 0.5)
ggsave(plot_z, filename = file.path(out_dir, ‘FeaturePlot.png’), width = 12, height = 6)