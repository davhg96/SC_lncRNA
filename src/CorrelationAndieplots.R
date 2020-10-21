library("DESeq2")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")
library("openxlsx")
library("GenomicRanges")

PCGdata <- read.table(file = "./data/Fcounts_2Ddiff_PcodingGenes_s1.txt", header = TRUE)
LNCdata <- read.table(file = "./data/Fcounts_2Ddiff_LncRNA_s1.txt", header = TRUE)