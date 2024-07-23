setwd("/home/micah/TWU/workshop_2024/twupa_bioinformaticsws_2024/data/genomefiles/de")
df <- read.table("GSE177494_Raw_gene_counts_matrix.txt",header=TRUE)

coldf <- data.frame(names = colnames(df)[2:85])

get_reagent <- function(x){
  strsplit(x,"_")[[1]][2]
}

get_dosage <- function(x){
  reagent <- strsplit(x,"_")[[1]][2]
  if (reagent == "DMSO"){
    return("NA")
  } else if (reagent %in% c("RBN012397","RBN012337")){
    return(strsplit(x,"_")[[1]][3])
  }
}

get_time <- function(x){
  reagent <- strsplit(x,"_")[[1]][2]
  if (reagent == "DMSO"){
    return(strsplit(x,"_")[[1]][3])
  } else if (reagent %in% c("RBN012397","RBN012337")){
    return(strsplit(x,"_")[[1]][4])
  }
}

get_replicate <- function(x){
  reagent <- strsplit(x,"_")[[1]][2]
  if (reagent == "DMSO"){
    return(strsplit(x,"_")[[1]][4])
  } else if (reagent %in% c("RBN012397","RBN012337")){
    return(strsplit(x,"_")[[1]][5])
  }
}

coldf$reagent <- factor(unlist(lapply(coldf$names, get_reagent)))
coldf$dosage <- factor(unlist(lapply(coldf$names, get_dosage)))
coldf$time <- factor(unlist(lapply(coldf$names, get_time)))
coldf$replicate <- factor(unlist(lapply(coldf$names, get_replicate)))
coldf

library(DESeq2)

rownames(df)<-df[,1]
df <- df[,2:85]

dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("10nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012397","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("10nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012397","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012397_dose_10nM_time_6hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012397_dose_10nM_time_6hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])



dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("10nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012397","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("10nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012397","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012397_dose_10nM_time_24hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012397_dose_10nM_time_24hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])


dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("30nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012397","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("30nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012397","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012397_dose_30nM_time_6hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012397_dose_30nM_time_6hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])


dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("30nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012397","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("30nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012397","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012397_dose_30nM_time_24hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012397_dose_30nM_time_24hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])


dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("300nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012397","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("300nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012397","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012397_dose_300nM_time_6hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012397_dose_300nM_time_6hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])


dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("300nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012397","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("300nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012397","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012397_dose_300nM_time_24hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012397_dose_300nM_time_24hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])

length(intersect(RBN012397_dose_300nM_time_24hr_up,RBN012397_dose_10nM_time_24hr_up ))

library("ggvenn")
A <- list(RBN012397_300nM=RBN012397_dose_300nM_time_24hr_up,RBN012397_10nM=RBN012397_dose_10nM_time_24hr_up,RBN012397_30nM=RBN012397_dose_30nM_time_24hr_up)
ggvenn(A)

library("VennDiagram")
draw.pairwise.venn(area1 = 217+189, area2 = 209, cross.area= 189)

library("eulerr")
fit1 <- euler(c("300nM"=201,"10nM"=19,"30nM"=5,"300nM&10nM"=82, "300nM&30nM"=16,"10nM&30nM"=1,"300nM&10nM&30nM"=107))
plot(fit1)

dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("10nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012337","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("10nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012337","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012337_dose_10nM_time_6hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012337_dose_10nM_time_6hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])



dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("10nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012337","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("10nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012337","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012337_dose_10nM_time_24hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012337_dose_10nM_time_24hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])


dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("30nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012337","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("30nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012337","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012337_dose_30nM_time_6hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012337_dose_30nM_time_6hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])


dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("30nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012337","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("30nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012337","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012337_dose_30nM_time_24hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012337_dose_30nM_time_24hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])


dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("300nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012337","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("300nM","NA") & coldf$time == "6hr" & coldf$reagent %in% c("RBN012337","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012337_dose_300nM_time_6hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012337_dose_300nM_time_6hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])


dds <- DESeqDataSetFromMatrix(countData=round(df[,coldf$dosage %in% c("300nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012337","DMSO")]), 
                              colData=coldf[coldf$dosage %in% c("300nM","NA") & coldf$time == "24hr" & coldf$reagent %in% c("RBN012337","DMSO"),],
                              design = ~ reagent)

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

RBN012337_dose_300nM_time_24hr_up <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange > 0,])
RBN012337_dose_300nM_time_24hr_down <- rownames(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.1 & resOrdered$log2FoldChange < 0,])
