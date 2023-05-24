# This scritps aims to identify DE lists of RNAseq samples (second run)


# Laod R libraries
library(DT)
library(data.table)
library(stringr)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(reshape2)
library(pheatmap)

library(purrr)     # map function
library(readr)     # read_tsv
library(dplyr)
library(plyr)
library(plotly)

library(RColorBrewer)
library(affy)
library(viridis)


# Load count data (raw counts, STAR output)
counts<-read.table("/Users/ieo4890/derenzini/datasets/RNAseq/total_counts_BCL2_massa_2run.txt")

# prepare the count matrix
countdata<-counts

# remove the gene ensembl version
# rownames(countdata)<-gsub("\\..*$","",rownames(countdata))
#pseudoCount = log2(countdata + 1)

# prepare the condition info
sample_condition<-c("E_NT","E_DOXO","E_ACTD","E_CX","BCL2_NT","BCL2_DOXO","BCL2_ACTD","BCL2_CX",
                    "E_NT","E_DOXO","E_ACTD","E_CX","BCL2_NT","BCL2_DOXO","BCL2_ACTD","BCL2_CX",
                    "E_NT","E_DOXO","E_ACTD","E_CX","BCL2_NT","BCL2_DOXO","BCL2_ACTD","BCL2_CX")
names(sample_condition)<-colnames(countdata)


### Preparation dds object
  # prepare table of sample information
coldata<-data.frame(condition=c("E_NT","E_DOXO","E_ACTD","E_CX","BCL2_NT","BCL2_DOXO","BCL2_ACTD","BCL2_CX","E_NT","E_DOXO","E_ACTD","E_CX","BCL2_NT","BCL2_DOXO","BCL2_ACTD","BCL2_CX","E_NT","E_DOXO","E_ACTD","E_CX","BCL2_NT","BCL2_DOXO","BCL2_ACTD","BCL2_CX"),Type=rep("paired-end",24) )
rownames(coldata)<-c( "S31137","S31138","S31139","S31140","S31141","S31142","S31143","S31144","S31145", "S31146","S31147","S31148","S31149","S31150","S31151","S31152","S31153","S31154","S31155","S31156", "S31157","S31158","S31159", "S31160")

head(countdata)
head(coldata)

# the columns of count matrix and the rows of column data must be in the same order
if(!all(rownames(coldata)==colnames(countdata))){
  countdata<-countdata[,rownames(coldata)]
}

all(rownames(coldata)==colnames(countdata))


# global DESeqData object
dds<-DESeqDataSetFromMatrix(countData=countdata,
                            colData=coldata,
                            design= ~ condition)

# track memory of original dds object
dds_original<-dds

# Gene Biotype
gene_biotype<-read.table("/Users/ieo4890/derenzini/datasets/RNAseq/gene_biotype.txt",skip=5)
# filter only gene_type
gene_biotype<-filter(gene_biotype,V3 == "gene_type")
# drop levels (e.g. column V4 has inside levels)
gene_biotype<-droplevels(gene_biotype)

### Exploratory Analysis: PCA 

  # Load object useful for PCA
rld<-readRDS("/Users/ieo4890/derenzini/datasets/RNAseq/rlog_bcl2_massa_global_2run.rds")

  # PCA computation
data <- plotPCA(rld, intgroup =c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

 # Extract the 3D PCA coordinates
rv <- rowVars(assay(rld))
ntop = 500
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))
percentVar_self <- pca$sdev^2/sum(pca$sdev^2)

pca_bcl2<-as.data.frame(pca$x[,1:3])
saveRDS(pca_bcl2,file="3d_pca_bcl2.rds")

 # PCA visualization

ggplot(data,aes(PC1,PC2,colour=condition))+geom_point(size=3)+
  ggtitle("All samples")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))+
  geom_text(aes(label=gsub("SQ_","",rownames(data))),hjust=0, vjust=2,size=3.5)+
  xlim(c(-20,30))

### Global DE Analysis

  # filtering: at least 20 counts for gene, is is not necesary this pre filtering
dds<-dds[rowSums(counts(dds))>20,]

dds.norm <-  estimateSizeFactors(dds)

sizeFactors(dds.norm)

  # Checking the normalization
epsilon<-1

 # Performing estimation of dispersion parameter
dds.disp <- estimateDispersions(dds.norm)

alpha <- 0.05
waldTestResult <- nbinomWaldTest(dds.disp)

####### Empty NT vs E doxo
   ## Comparison Empty NT vs EMPTY DOXO

cmp<-c("E_NT","E_DOXO")
pos<-which(coldata[,1] %in% cmp) 

res_Ent_vs_Edoxo <- results(waldTestResult, alpha=alpha,contrast=c("condition", "E_DOXO", "E_NT"), pAdjustMethod="BH") 
# Order the table by decreasing p-value 
res_Ent_vs_Edoxo <- res_Ent_vs_Edoxo[order(res_Ent_vs_Edoxo$padj),]

summary(res_Ent_vs_Edoxo) 
table(res_Ent_vs_Edoxo$padj < 0.05) 


# recover gene annotation 
normalized.counts <- as.data.frame(counts(dds.norm, normalized=TRUE )) 
gene_ids = rownames(normalized.counts) 

# discard NA values 
res_sig_Ent_vs_Edoxo <- res_Ent_vs_Edoxo[!is.na(res_Ent_vs_Edoxo$padj),] 

#recover gene annotation for DE genes 
res_sig_Ent_vs_Edoxo$ensembl_gene_id = rownames(res_sig_Ent_vs_Edoxo) 
res_sig_Ent_vs_Edoxo <- as.data.frame(res_sig_Ent_vs_Edoxo) 

# annotation directly from transcriptome (gene_biotype.txt) 
genes_info=gene_biotype[match(rownames(res_sig_Ent_vs_Edoxo),gene_biotype$V2),] 

# mantain only useful information 
genes_info=genes_info[,c(2,4,6)] 

colnames(genes_info)<-c("ensembl_gene_id","gene_type","gene_name")


res.anno_Ent_vs_Edoxo <- merge (res_sig_Ent_vs_Edoxo, genes_info, by ="ensembl_gene_id", all.x = TRUE) 
res.anno_Ent_vs_Edoxo<- res.anno_Ent_vs_Edoxo[order(res.anno_Ent_vs_Edoxo$padj),] 

 # !!!Note: this file is saved /scripts/BCL2/massa/
#write.csv(res.anno_Ent_vs_Edoxo,file = "DE_Ent_vs_Edoxo.csv",row.names = F) 

 # heatmap on DE genes (heatmap on DE genes with padj<0.05)

#rm(temp,temp_extended)     
# all DE genes 
alpha=0.05 
# recover gene annotation 
temp<-res_Ent_vs_Edoxo[res_Ent_vs_Edoxo$padj <= alpha & !is.na(res_Ent_vs_Edoxo$padj),] 
temp$ensembl_gene_id<-rownames(temp) 

temp_extended<-merge(as.data.frame(temp),genes_info,by.y="ensembl_gene_id",all.x=TRUE) 
# reorder according the pvalue 
temp_extended<-temp_extended[order(temp_extended$padj),] 

# protein_coding_filter
temp_extended_protein<-filter(temp_extended,gene_type=="protein_coding")

# !!!Note: this file is saved /scripts/BCL2/massa/
#write.csv(temp_extended,file="DE_analysis_filtered_Ent_vs_Edoxo.csv",row.names = F) 

gene.kept_Ent_vs_Edoxo <- temp_extended$ensembl_gene_id 
gene.kept.symbol_Ent_vs_Edoxo<-as.character(temp_extended$gene_name) 

countTable.kept_Ent_vs_Edoxo <- log2(counts(dds.norm,normalized=T) +1)[gene.kept_Ent_vs_Edoxo, ] 
rownames(countTable.kept_Ent_vs_Edoxo)<-gene.kept.symbol_Ent_vs_Edoxo 

# select only samples of investigation 
countTable.kept_Ent_vs_Edoxo<-countTable.kept_Ent_vs_Edoxo[,pos] 


# insert the annotation column 
annotation_col = data.frame(Condition = factor(rep(c("E_NT","E_DOXO"),times=3))) 
rownames(annotation_col)<-colnames(countTable.kept_Ent_vs_Edoxo) 
color_annotation_col = list(Condition = c(E_DOXO ="#3366ff", E_NT="#ff9900")) 


pheatmap(as.matrix((countTable.kept_Ent_vs_Edoxo)),scale="row", 
         cluster_rows = T,cluster_cols=T, 
         annotation_col = annotation_col, 
         annotation_colors=color_annotation_col[1], 
         #cellheight=1.8, 
         show_rownames  = F, 
         main=" Empty NT vs Empty DOXO", 
         fontsize_row = 4 
)


## Empty NT vs EMPTY ACTD 
cmp<-c("E_NT","E_ACTD") 
pos<-which(coldata[,1] %in% cmp) 


# DE  
res_Ent_vs_Eactd <- results(waldTestResult, alpha=alpha,contrast=c("condition", "E_ACTD", "E_NT"), pAdjustMethod="BH") 
# Order the table by decreasing p-valuer 
res_Ent_vs_Eactd <- res_Ent_vs_Eactd[order(res_Ent_vs_Eactd$padj),] 

summary(res_Ent_vs_Eactd) 
table(res_Ent_vs_Eactd$padj < 0.05) 


# recover gene annotation 
normalized.counts <- as.data.frame(counts(dds.norm, normalized=TRUE )) 
gene_ids = rownames(normalized.counts) 

# discard NA values 
res_sig_Ent_vs_Eactd <- res_Ent_vs_Eactd[!is.na(res_Ent_vs_Eactd$padj),] 

#recover gene annotation for DE genes 
res_sig_Ent_vs_Eactd$ensembl_gene_id = rownames(res_sig_Ent_vs_Eactd) 
res_sig_Ent_vs_Eactd <- as.data.frame(res_sig_Ent_vs_Eactd) 

# annotation directly from transcriptome (gene_biotype.txt) 
genes_info=gene_biotype[match(rownames(res_sig_Ent_vs_Eactd),gene_biotype$V2),] 

# mantain only useful information 
genes_info=genes_info[,c(2,4,6)] 
colnames(genes_info)<-c("ensembl_gene_id","gene_type","gene_name") 


res.anno_Ent_vs_Eactd <- merge (res_sig_Ent_vs_Eactd, genes_info, by ="ensembl_gene_id", all.x = TRUE) 
res.anno_Ent_vs_Eactd<- res.anno_Ent_vs_Eactd[order(res.anno_Ent_vs_Eactd$padj),] 

write.csv(res.anno_Ent_vs_Eactd,file = "DE_Ent_vs_Eactd.csv",row.names = F) 


# heatmap on DE genes 
rm(temp,temp_extended) 
# all DE genes 
alpha=0.05 
# recover gene annotation 
temp<-res_Ent_vs_Eactd[res_Ent_vs_Eactd$padj <= alpha & !is.na(res_Ent_vs_Eactd$padj),] 
temp$ensembl_gene_id<-rownames(temp) 

temp_extended<-merge(as.data.frame(temp),genes_info,by.y="ensembl_gene_id",all.x=TRUE) 
# reorder according the pvalue 
temp_extended<-temp_extended[order(temp_extended$padj),] 

# protein_coding_filter
temp_extended_protein<-filter(temp_extended,gene_type=="protein_coding")

#write.csv(temp_extended,file = "DE_analysis_filtered_Ent_vs_Eactd.csv",row.names = F) 

gene.kept_Ent_vs_Eactd <- temp_extended$ensembl_gene_id 
gene.kept.symbol_Ent_vs_Eactd<-as.character(temp_extended$gene_name) 

countTable.kept_Ent_vs_Eactd <- log2(counts(dds.norm,normalized=T) +1)[gene.kept_Ent_vs_Eactd, ] 
rownames(countTable.kept_Ent_vs_Eactd)<-gene.kept.symbol_Ent_vs_Eactd 

# select only samples of investigation 
countTable.kept_Ent_vs_Eactd<-countTable.kept_Ent_vs_Eactd[,pos] 


# insert the annotation column 
annotation_col = data.frame(Condition = factor(rep(c("E_NT","E_ACTD"),times=3))) 
rownames(annotation_col)<-colnames(countTable.kept_Ent_vs_Eactd) 
color_annotation_col = list(Condition = c(E_ACTD ="#3366ff", E_NT="#ff9900")) 


pheatmap(as.matrix((countTable.kept_Ent_vs_Eactd)),scale="row", 
         cluster_rows = T,cluster_cols=T, 
         annotation_col = annotation_col, 
         annotation_colors=color_annotation_col[1], 
         #clustering_distance_rows="correlation", 
         #clustering_distance_cols="euclidean", 
         #clustering_method="ward.D2", 
         #cellheight=1.8, 
         show_rownames  = F, 
         main=" Empty NT vs Empty ACTD", 
         fontsize_row = 4 
) 


### Empty NT vs EMPTY CX 
cmp<-c("E_NT","E_CX") 
pos<-which(coldata[,1] %in% cmp) 

# DE 
res_Ent_vs_Ecx <- results(waldTestResult, alpha=alpha,contrast=c("condition", "E_CX", "E_NT"), pAdjustMethod="BH") 
# Order the table by decreasing p-valuer 
res_Ent_vs_Ecx <- res_Ent_vs_Ecx[order(res_Ent_vs_Ecx$padj),] 

summary(res_Ent_vs_Ecx) 
table(res_Ent_vs_Ecx$padj < 0.05) 


# recover gene annotation 
normalized.counts <- as.data.frame(counts(dds.norm, normalized=TRUE )) 
gene_ids = rownames(normalized.counts) 

# discard NA values 
res_sig_Ent_vs_Ecx <- res_Ent_vs_Ecx[!is.na(res_Ent_vs_Ecx$padj),] 

#recover gene annotation for DE genes 
res_sig_Ent_vs_Ecx$ensembl_gene_id = rownames(res_sig_Ent_vs_Ecx) 
res_sig_Ent_vs_Ecx <- as.data.frame(res_sig_Ent_vs_Ecx) 

# annotation directly from transcriptome (gene_biotype.txt) 
genes_info=gene_biotype[match(rownames(res_sig_Ent_vs_Ecx),gene_biotype$V2),] 

# mantain only useful information 
genes_info=genes_info[,c(2,4,6)] 

colnames(genes_info)<-c("ensembl_gene_id","gene_type","gene_name") 


res.anno_Ent_vs_Ecx <- merge (res_sig_Ent_vs_Ecx, genes_info, by ="ensembl_gene_id", all.x = TRUE) 
res.anno_Ent_vs_Ecx<- res.anno_Ent_vs_Ecx[order(res.anno_Ent_vs_Ecx$padj),] 

#write.csv(res.anno_Ent_vs_Ecx,file ="DE_Ent_vs_Ecx.csv",row.names = F) 


# heatmap on DE genes 

# all DE genes 
rm(temp,temp_extended) 
alpha=0.05 
# recover gene annotation 
temp<-res_Ent_vs_Ecx[res_Ent_vs_Ecx$padj <= alpha & !is.na(res_Ent_vs_Ecx$padj),] 
temp$ensembl_gene_id<-rownames(temp) 

temp_extended<-merge(as.data.frame(temp),genes_info,by.y="ensembl_gene_id",all.x=TRUE) 
# reorder according the pvalue 
temp_extended<-temp_extended[order(temp_extended$padj),] 

# protein_coding_filter
temp_extended_protein<-filter(temp_extended,gene_type=="protein_coding")


#write.csv(temp_extended,file = "DE_analysis_filtered_Ent_vs_Ecx.csv",row.names = F) 

gene.kept_Ent_vs_Ecx <- temp_extended$ensembl_gene_id 
gene.kept.symbol_Ent_vs_Ecx<-as.character(temp_extended$gene_name) 

countTable.kept_Ent_vs_Ecx <- log2(counts(dds.norm,normalized=T) +1)[gene.kept_Ent_vs_Ecx, ] 
rownames(countTable.kept_Ent_vs_Ecx)<-gene.kept.symbol_Ent_vs_Ecx 

# select only samples of investigation 
countTable.kept_Ent_vs_Ecx<-countTable.kept_Ent_vs_Ecx[,pos] 


# insert the annotation column 
annotation_col = data.frame(Condition = factor(rep(c("E_CNT","E_CX"),times=3))) 
rownames(annotation_col)<-colnames(countTable.kept_Ent_vs_Ecx) 
color_annotation_col = list(Condition = c(E_CX ="#3366ff", E_CNT="#ff9900")) 


pheatmap(as.matrix((countTable.kept_Ent_vs_Ecx)[,c(1,3,5,2,4,6)]),scale="row", 
         cluster_rows = T,cluster_cols=F, 
         annotation_col = annotation_col, 
         annotation_colors=color_annotation_col[1], 
         #cellheight=1.8, 
         show_rownames  = F, 
         main=" Empty NT vs Empty CX", 
         fontsize_row = 4 
) 




cmp<-c("BCL2_NT","BCL2_DOXO") 
pos<-which(coldata[,1] %in% cmp) 

# DE 
res_BCL2nt_vs_BCL2doxo <- results(waldTestResult, alpha=alpha,contrast=c("condition", "BCL2_DOXO", "BCL2_NT"), pAdjustMethod="BH") 
# Order the table by decreasing p-valuer 
res_BCL2nt_vs_BCL2doxo <- res_BCL2nt_vs_BCL2doxo[order(res_BCL2nt_vs_BCL2doxo$padj),] 

summary(res_BCL2nt_vs_BCL2doxo) 
table(res_BCL2nt_vs_BCL2doxo$padj < 0.05) 


# recover gene annotation 
normalized.counts <- as.data.frame(counts(dds.norm, normalized=TRUE )) 
gene_ids = rownames(normalized.counts) 

# discard NA values 
res_sig_BCL2nt_vs_BCL2doxo <- res_BCL2nt_vs_BCL2doxo[!is.na(res_BCL2nt_vs_BCL2doxo$padj),] 

#recover gene annotation for DE genes 
res_sig_BCL2nt_vs_BCL2doxo$ensembl_gene_id = rownames(res_sig_BCL2nt_vs_BCL2doxo) 
res_sig_BCL2nt_vs_BCL2doxo <- as.data.frame(res_sig_BCL2nt_vs_BCL2doxo) 

# annotation directly from transcriptome (gene_biotype.txt) 
genes_info=gene_biotype[match(rownames(res_sig_BCL2nt_vs_BCL2doxo),gene_biotype$V2),] 

# mantain only useful information 
genes_info=genes_info[,c(2,4,6)] 

colnames(genes_info)<-c("ensembl_gene_id","gene_type","gene_name") 


res.anno_BCL2nt_vs_BCL2doxo <- merge (res_sig_BCL2nt_vs_BCL2doxo, genes_info, by ="ensembl_gene_id", all.x = TRUE) 
res.anno_BCL2nt_vs_BCL2doxo<- res.anno_BCL2nt_vs_BCL2doxo[order(res.anno_BCL2nt_vs_BCL2doxo$padj),] 

#write.csv(res.anno_BCL2nt_vs_BCL2doxo,file ="DE_BCL2nt_vs_BCL2doxo.csv",row.names = F) 


# heatmap on DE genes 

# all DE genes 
rm(temp,temp_extended) 
alpha=0.05 
# recover gene annotation 
temp<-res_BCL2nt_vs_BCL2doxo[res_BCL2nt_vs_BCL2doxo$padj <= alpha & !is.na(res_BCL2nt_vs_BCL2doxo$padj),] 
temp$ensembl_gene_id<-rownames(temp) 

temp_extended<-merge(as.data.frame(temp),genes_info,by.y="ensembl_gene_id",all.x=TRUE) 
# reorder according the pvalue 
temp_extended<-temp_extended[order(temp_extended$padj),] 

# protein_coding_filter
temp_extended_protein<-filter(temp_extended,gene_type=="protein_coding")


#write.csv(temp_extended,file = "DE_analysis_filtered_BCL2nt_vs_BCL2doxo.csv",row.names = F) 

gene.kept_BCL2nt_vs_BCL2doxo <- temp_extended$ensembl_gene_id 
gene.kept.symbol_BCL2nt_vs_BCL2doxo<-as.character(temp_extended$gene_name) 

countTable.kept_BCL2nt_vs_BCL2doxo <- log2(counts(dds.norm,normalized=T) +1)[gene.kept_BCL2nt_vs_BCL2doxo, ] 
rownames(countTable.kept_BCL2nt_vs_BCL2doxo)<-gene.kept.symbol_BCL2nt_vs_BCL2doxo 

# select only samples of investigation 
countTable.kept_BCL2nt_vs_BCL2doxo<-countTable.kept_BCL2nt_vs_BCL2doxo[,pos] 


# insert the annotation column 
annotation_col = data.frame(Condition = factor(rep(c("BCL2_NT","BCL2_DOXO"),times=3))) 
rownames(annotation_col)<-colnames(countTable.kept_BCL2nt_vs_BCL2doxo) 
color_annotation_col = list(Condition = c(BCL2_NT ="#3366ff", BCL2_DOXO="#ff9900")) 


out<-pheatmap(as.matrix((countTable.kept_BCL2nt_vs_BCL2doxo)),scale="row", 
         cluster_rows = T,cluster_cols=T, 
         annotation_col = annotation_col, 
         annotation_colors=color_annotation_col[1], 
         #cellheight=1.8, 
         show_rownames  = F, 
         main=" BCL2 NT vs BCL2 DOXO", 
         fontsize_row = 4
) 

res <- as.matrix((countTable.kept_BCL2nt_vs_BCL2doxo))[c(out$tree_row[["order"]]),out$tree_col[["order"]]]
res<-as.data.frame(res)

temp_out_cluster<-(sort(cutree(out$tree_row, k=2)))
out_cluster<-data.frame(order=sort(cutree(out$tree_row, k=2)))
out_cluster$gene<-names(temp_out_cluster)

res$cluster<-out_cluster$order[ order(match(out_cluster$gene,rownames(res)))]

### manual division
row_clust2<-res[(2080:3641),]
row_clust1<-res[(1:2079),]

new_clust<-rbind(row_clust2,row_clust1)

pheatmap(as.matrix((new_clust)),scale="row", 
         cluster_rows = F,cluster_cols=F, 
         annotation_col = annotation_col, 
         annotation_colors=color_annotation_col[1], 
         #cellheight=1.8, 
         show_rownames  = F, 
         main=" BCL2 NT vs BCL2 DOXO", 
         fontsize_row = 4 
     
) 


# prepare matrix for pheatmap

### BCL2 NT vs BCL2 ACTD
cmp<-c("BCL2_NT","BCL2_ACTD") 
pos<-which(coldata[,1] %in% cmp) 


res_BCL2nt_vs_BCL2actd <- results(waldTestResult, alpha=alpha,contrast=c("condition", "BCL2_ACTD", "BCL2_NT"), pAdjustMethod="BH") 
# Order the table by decreasing p-valuer 
res_BCL2nt_vs_BCL2actd <- res_BCL2nt_vs_BCL2actd[order(res_BCL2nt_vs_BCL2actd$padj),] 

summary(res_BCL2nt_vs_BCL2actd) 
table(res_BCL2nt_vs_BCL2actd$padj < 0.05) 


# recover gene annotation 
normalized.counts <- as.data.frame(counts(dds.norm, normalized=TRUE )) 
gene_ids = rownames(normalized.counts) 

# discard NA values 
res_sig_BCL2nt_vs_BCL2actd <- res_BCL2nt_vs_BCL2actd[!is.na(res_BCL2nt_vs_BCL2actd$padj),] 

#recover gene annotation for DE genes 
res_sig_BCL2nt_vs_BCL2actd$ensembl_gene_id = rownames(res_sig_BCL2nt_vs_BCL2actd) 
res_sig_BCL2nt_vs_BCL2actd <- as.data.frame(res_sig_BCL2nt_vs_BCL2actd) 

# annotation directly from transcriptome (gene_biotype.txt) 
genes_info=gene_biotype[match(rownames(res_sig_BCL2nt_vs_BCL2actd),gene_biotype$V2),] 

# mantain only useful information 
genes_info=genes_info[,c(2,4,6)] 

colnames(genes_info)<-c("ensembl_gene_id","gene_type","gene_name") 


res.anno_BCL2nt_vs_BCL2actd <- merge (res_sig_BCL2nt_vs_BCL2actd, genes_info, by ="ensembl_gene_id", all.x = TRUE) 
res.anno_BCL2nt_vs_BCL2actd<- res.anno_BCL2nt_vs_BCL2actd[order(res.anno_BCL2nt_vs_BCL2actd$padj),] 

#write.csv(res.anno_BCL2nt_vs_BCL2actd,file ="DE_BCL2nt_vs_BCL2actd.csv",row.names = F) 

# heatmap on DE genes 

# all DE genes 
rm(temp,temp_extended) 
alpha=0.05 
# recover gene annotation 
temp<-res_BCL2nt_vs_BCL2actd[res_BCL2nt_vs_BCL2actd$padj <= alpha & !is.na(res_BCL2nt_vs_BCL2actd$padj),] 
temp$ensembl_gene_id<-rownames(temp) 

temp_extended<-merge(as.data.frame(temp),genes_info,by.y="ensembl_gene_id",all.x=TRUE) 
# reorder according the pvalue 
temp_extended<-temp_extended[order(temp_extended$padj),] 

# protein_coding_filter
temp_extended_protein<-filter(temp_extended,gene_type=="protein_coding")


#write.csv(temp_extended,file = "DE_analysis_filtered_BCL2nt_vs_BCL2actd.csv",row.names = F) 

gene.kept_BCL2nt_vs_BCL2actd <- temp_extended$ensembl_gene_id 
gene.kept.symbol_BCL2nt_vs_BCL2actd<-as.character(temp_extended$gene_name) 

countTable.kept_BCL2nt_vs_BCL2actd <- log2(counts(dds.norm,normalized=T) +1)[gene.kept_BCL2nt_vs_BCL2actd, ] 
rownames(countTable.kept_BCL2nt_vs_BCL2actd)<-gene.kept.symbol_BCL2nt_vs_BCL2actd 

# select only samples of investigation 
countTable.kept_BCL2nt_vs_BCL2actd<-countTable.kept_BCL2nt_vs_BCL2actd[,pos] 


# insert the annotation column 
annotation_col = data.frame(Condition = factor(rep(c("BCL2_NT","BCL2_ACTD"),times=3))) 
rownames(annotation_col)<-colnames(countTable.kept_BCL2nt_vs_BCL2actd) 
color_annotation_col = list(Condition = c(BCL2_NT ="#3366ff", BCL2_ACTD="#ff9900")) 


pheatmap(as.matrix((countTable.kept_BCL2nt_vs_BCL2actd)[,c(5,1,3,6,2,4)]),scale="row", 
         cluster_rows = T,cluster_cols=F, 
         annotation_col = annotation_col, 
         annotation_colors=color_annotation_col[1], 
         #cellheight=1.8, 
         show_rownames  = F, 
         main=" BCL2 NT vs BCL2 ACTD", 
         fontsize_row = 4 
) 




### BCL2 NT vs BCL2 CX
cmp<-c("BCL2_NT","BCL2_CX") 
pos<-which(coldata[,1] %in% cmp) 


res_BCL2nt_vs_BCL2cx <- results(waldTestResult, alpha=alpha,contrast=c("condition", "BCL2_CX", "BCL2_NT"), pAdjustMethod="BH") 
# Order the table by decreasing p-valuer 
res_BCL2nt_vs_BCL2cx <- res_BCL2nt_vs_BCL2cx[order(res_BCL2nt_vs_BCL2cx$padj),] 

summary(res_BCL2nt_vs_BCL2cx) 
table(res_BCL2nt_vs_BCL2cx$padj < 0.05) 


# recover gene annotation 
normalized.counts <- as.data.frame(counts(dds.norm, normalized=TRUE )) 
gene_ids = rownames(normalized.counts) 

# discard NA values 
res_sig_BCL2nt_vs_BCL2cx <- res_BCL2nt_vs_BCL2cx[!is.na(res_BCL2nt_vs_BCL2cx$padj),] 

#recover gene annotation for DE genes 
res_sig_BCL2nt_vs_BCL2cx$ensembl_gene_id = rownames(res_sig_BCL2nt_vs_BCL2cx) 
res_sig_BCL2nt_vs_BCL2cx <- as.data.frame(res_sig_BCL2nt_vs_BCL2cx) 

# annotation directly from transcriptome (gene_biotype.txt) 
genes_info=gene_biotype[match(rownames(res_sig_BCL2nt_vs_BCL2cx),gene_biotype$V2),] 

# mantain only useful information 
genes_info=genes_info[,c(2,4,6)] 

colnames(genes_info)<-c("ensembl_gene_id","gene_type","gene_name") 


res.anno_BCL2nt_vs_BCL2cx <- merge (res_sig_BCL2nt_vs_BCL2cx, genes_info, by ="ensembl_gene_id", all.x = TRUE) 
res.anno_BCL2nt_vs_BCL2cx<- res.anno_BCL2nt_vs_BCL2cx[order(res.anno_BCL2nt_vs_BCL2cx$padj),] 

#write.csv(res.anno_BCL2nt_vs_BCL2cx,file ="DE_BCL2nt_vs_BCL2cx.csv",row.names = F) 


# heatmap on DE genes 

# all DE genes 
rm(temp,temp_extended) 
alpha=0.05 
# recover gene annotation 
temp<-res_BCL2nt_vs_BCL2cx[res_BCL2nt_vs_BCL2cx$padj <= alpha & !is.na(res_BCL2nt_vs_BCL2cx$padj),] 
temp$ensembl_gene_id<-rownames(temp) 

temp_extended<-merge(as.data.frame(temp),genes_info,by.y="ensembl_gene_id",all.x=TRUE) 
# reorder according the pvalue 
temp_extended<-temp_extended[order(temp_extended$padj),] 

# protein_coding_filter
temp_extended_protein<-filter(temp_extended,gene_type=="protein_coding")


#write.csv(temp_extended,file = "DE_analysis_filtered_BCL2nt_vs_BCL2cx.csv",row.names = F) 

gene.kept_BCL2nt_vs_BCL2cx <- temp_extended$ensembl_gene_id 
gene.kept.symbol_BCL2nt_vs_BCL2cx<-as.character(temp_extended$gene_name) 

countTable.kept_BCL2nt_vs_BCL2cx <- log2(counts(dds.norm,normalized=T) +1)[gene.kept_BCL2nt_vs_BCL2cx, ] 
rownames(countTable.kept_BCL2nt_vs_BCL2cx)<-gene.kept.symbol_BCL2nt_vs_BCL2cx 

# select only samples of investigation 
countTable.kept_BCL2nt_vs_BCL2cx<-countTable.kept_BCL2nt_vs_BCL2cx[,pos] 


# insert the annotation column 
annotation_col = data.frame(Condition = factor(rep(c("BCL2_NT","BCL2_CX"),times=3))) 
rownames(annotation_col)<-colnames(countTable.kept_BCL2nt_vs_BCL2cx) 
color_annotation_col = list(Condition = c(BCL2_NT ="#3366ff", BCL2_CX="#ff9900")) 


pheatmap(as.matrix((countTable.kept_BCL2nt_vs_BCL2cx)[,c(5,1,3,6,2,4)]),scale="row", 
         cluster_rows = T,cluster_cols=F, 
         annotation_col = annotation_col, 
         annotation_colors=color_annotation_col[1], 
         #cellheight=1.8, 
         show_rownames  = F, 
         main=" BCL2 NT vs BCL2 CX", 
         fontsize_row = 4 
) 






