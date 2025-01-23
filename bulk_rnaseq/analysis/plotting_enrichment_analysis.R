# ## This script aims to compute Enrichment Analysis for DE genes (RNAseq)
# Author: Saveria Mazzara

# load libraries 
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA)
library(venn)
library(stringr)


# load DE list
len<-10
DEG_list <- vector(mode = "list", length = len)


#1. 
temp_Ent_vs_BCL2nt <- read.csv("DE_analysis_filtered_Ent_vs_BCL2nt.csv")

#2. 
temp_Ent_vs_Edoxo <- read.csv("DE_analysis_filtered_Ent_vs_Edoxo.csv")

#3. 
temp_Ent_vs_Eactd <- read.csv("DE_analysis_filtered_Ent_vs_Eactd.csv")

#4. 
temp_Ent_vs_Ecx <- read.csv("DE_analysis_filtered_Ent_vs_Ecx.csv")

#5.
temp_BCL2nt_vs_BCL2doxo <- read.csv("DE_analysis_filtered_BCL2nt_vs_BCL2doxo.csv")

#6. 
temp_BCL2nt_vs_BCL2actd <- read.csv("DE_analysis_filtered_BCL2nt_vs_BCL2actd.csv")

#7. 
temp_BCL2nt_vs_BCL2cx <- read.csv("DE_analysis_filtered_BCL2nt_vs_BCL2cx.csv")

#8. 
temp_Edoxo_vs_BCL2doxo <- read.csv("DE_analysis_filtered_Edoxo_vs_BCL2doxo.csv")

#9. 
temp_Eactd_vs_BCL2actd <- read.csv("DE_analysis_filtered_Eactd_vs_BCL2actd.csv")

#10. 
temp_Ecx_vs_BCL2cx <- read.csv("DE_analysis_filtered_Ecx_vs_BCL2cx.csv")


# List: DE genes
DEG_list[[1]]<-temp_Ent_vs_BCL2nt$gene_name
DEG_list[[2]]<-temp_Ent_vs_Edoxo$gene_name
DEG_list[[3]]<-temp_Ent_vs_Eactd$gene_name
DEG_list[[4]]<-temp_Ent_vs_Ecx$gene_name
DEG_list[[5]]<-temp_BCL2nt_vs_BCL2doxo$gene_name
DEG_list[[6]]<-temp_BCL2nt_vs_BCL2actd$gene_name
DEG_list[[7]]<-temp_BCL2nt_vs_BCL2cx$gene_name
DEG_list[[8]]<-temp_Edoxo_vs_BCL2doxo$gene_name
DEG_list[[9]]<-temp_Eactd_vs_BCL2actd$gene_name
DEG_list[[10]]<-temp_Ecx_vs_BCL2cx$gene_name


names(DEG_list)<-c("Ent_vs_BCL2nt","Ent_vs_Edoxo","Ent_vs_Eactd","Ent_vs_Ecx","BCL2nt_vs_BCL2doxo",
                   "BCL2nt_vs_BCL2actd","BCL2nt_vs_BCL2cx","Edoxo_vs_BCL2doxo","Eactd_vs_BCL2actd","Ecx_vs_BCL2cx")


# COnversion to EntrezID
# 
list1<-bitr((DEG_list[[1]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list2<-bitr((DEG_list[[2]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list3<-bitr((DEG_list[[3]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list4<-bitr((DEG_list[[4]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list5<-bitr((DEG_list[[5]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list6<-bitr((DEG_list[[6]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list7<-bitr((DEG_list[[7]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list8<-bitr((DEG_list[[8]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list9<-bitr((DEG_list[[9]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
list10<-bitr((DEG_list[[10]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")


my_gcSample<-list()
my_gcSample[[1]]<-list1$ENTREZID
my_gcSample[[2]]<-list2$ENTREZID
my_gcSample[[3]]<-list3$ENTREZID
my_gcSample[[4]]<-list4$ENTREZID
my_gcSample[[5]]<-list5$ENTREZID
my_gcSample[[6]]<-list6$ENTREZID
my_gcSample[[7]]<-list7$ENTREZID
my_gcSample[[8]]<-list8$ENTREZID
my_gcSample[[9]]<-list9$ENTREZID
my_gcSample[[10]]<-list10$ENTREZID


names(my_gcSample)<-c("Ent_vs_BCL2nt","Ent_vs_Edoxo","Ent_vs_Eactd","Ent_vs_Ecx","BCL2nt_vs_BCL2doxo",
                      "BCL2nt_vs_BCL2actd","BCL2nt_vs_BCL2cx","Edoxo_vs_BCL2doxo","Eactd_vs_BCL2actd","Ecx_vs_BCL2cx")


# GO (optional)
#ck<-compareCluster(geneClusters = my_gcSample, fun="enrichGO",OrgDb="org.Hs.eg.db", pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05)

#p<-dotplot(ck)
# change the angle x label
#p + theme(axis.text.x=element_text(angle=90, hjust=1))

# KEGG pathways
ck_kegg<-compareCluster(geneClusters = my_gcSample, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)
ck_kegg_partial<-compareCluster(geneClusters = my_gcSample[2:7], fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)


ck_kegg_partial <- setReadable(ck_kegg_partial, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck_kegg_partial) 

 # Visualization 
p_kegg<-dotplot(ck_kegg)
# change the angle x label
p_kegg + theme(axis.text.x=element_text(angle=90, hjust=1))

dotplot(ck_kegg_partial,showCategory=6)+ theme(axis.text.x=element_text(angle=90, hjust=1))
cnetplot(ck_kegg_partial,showCategory=c("Apoptosis","Cell cycle","p53 signaling pathway"))
