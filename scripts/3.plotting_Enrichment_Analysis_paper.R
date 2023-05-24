# This scritp aims to compute Pathway Enrichment Analysis for DE genes (RNAseq)

# Load libraries
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


# GO
ck<-compareCluster(geneClusters = my_gcSample, fun="enrichGO",OrgDb="org.Hs.eg.db", pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05)

p<-dotplot(ck)
# change the angle x label
p + theme(axis.text.x=element_text(angle=90, hjust=1))

# KEGG pathways
ck_kegg<-compareCluster(geneClusters = my_gcSample, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)
ck_kegg_partial<-compareCluster(geneClusters = my_gcSample[2:7], fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)


ck_kegg_partial <- setReadable(ck_kegg_partial, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck_kegg_partial) 

# to have results in readable way
#prova<-ck_kegg_partial@compareClusterResult

 # Visualization 

#p_kegg<-dotplot(ck_kegg,showCategory=21)
p_kegg<-dotplot(ck_kegg)
# change the angle x label
p_kegg + theme(axis.text.x=element_text(angle=90, hjust=1))

dotplot(ck_kegg_partial,showCategory=6)+ theme(axis.text.x=element_text(angle=90, hjust=1))
cnetplot(ck_kegg_partial,showCategory=c("Apoptosis","Cell cycle","p53 signaling pathway"))


#### Apoptosis Genes
p_kegg_apopt<-ck_kegg_partial@compareClusterResult
p_kegg_apopt<-filter(p_kegg_apopt,ID=="hsa04210")

# older version
# Ent vs BCL2 bt
# genes_apop_ent_vs_bcl2nt<-list1$SYMBOL[match(unlist( str_split(p_kegg_apopt$geneID[1], "/", n = Inf, simplify = FALSE)),list1$ENTREZID)]
# 
# # Edoxo vs BCL2 doxo
# genes_apop_edoxo_vs_bcl2doxo<-list8$SYMBOL[match(unlist( str_split(p_kegg_apopt$geneID[5], "/", n = Inf, simplify = FALSE)),list8$ENTREZID)]
# 
# 
# # Ent vs Edoxo
# genes_apop_ent_vs_edoxo<-list2$SYMBOL[match(unlist( str_split(p_kegg_apopt$geneID[1], "/", n = Inf, simplify = FALSE)),list2$ENTREZID)]
# 
# # Ent vs Eactd
# genes_apop_ent_vs_eactd<-list3$SYMBOL[match(unlist( str_split(p_kegg_apopt$geneID[2], "/", n = Inf, simplify = FALSE)),list3$ENTREZID)]
# 
# # Ent vs Ecx
# genes_apop_ent_vs_ecx<-list4$SYMBOL[match(unlist( str_split(p_kegg_apopt$geneID[3], "/", n = Inf, simplify = FALSE)),list4$ENTREZID)]


# Read TP53 target database: TRANSFAC, Fischer, iRegulon

TP53_targets<- as.data.frame(read_excel("~/derenzini/datasets/TP53_targets_transfac.xlsx",  col_names = FALSE))
colnames(TP53_targets)<-"gene_name"
# Ent vs Edoxo
genes_apop_ent_vs_edoxo<-unlist( str_split(p_kegg_apopt$geneID[1], "/", n = Inf, simplify = FALSE))
length(intersect(genes_apop_ent_vs_edoxo,TP53_targets$gene_name))
# Ent vs Eactd
genes_apop_ent_vs_eactd<-unlist( str_split(p_kegg_apopt$geneID[2], "/", n = Inf, simplify = FALSE))
length(intersect(genes_apop_ent_vs_eactd,TP53_targets$gene_name))
# Ent vs Ecx
  genes_apop_ent_vs_ecx<-unlist( str_split(p_kegg_apopt$geneID[3], "/", n = Inf, simplify = FALSE))
  length(intersect(genes_apop_ent_vs_ecx,TP53_targets$gene_name))
# BCL2nt vs BCL2cx
  genes_apop_bcl2nt_vs_bcl2cx<-unlist( str_split(p_kegg_apopt$geneID[4], "/", n = Inf, simplify = FALSE))
  length(intersect(genes_apop_bcl2nt_vs_bcl2cx,TP53_targets$gene_name))

# list apoptosis genes
apop_list <- vector(mode = "list", length = 4)
apop_list[[1]]<-genes_apop_ent_vs_edoxo
apop_list[[2]]<-genes_apop_ent_vs_ecx
apop_list[[3]]<-genes_apop_ent_vs_eactd
names(apop_list)<-c("Ent_vs_Edoxo","Ent_vs_Ecx","Ent_vs_Eactd")
# plotting Venn diagram
venn(apop_list[1:3],col="black",zcolor=viridis(3),ilcs = 1.3,sncs=0.9)





# Reactome pathways
ck_reac<-compareCluster(geneClusters = my_gcSample, fun="enrichPathway", pvalueCutoff=0.05)
p_reac<-dotplot(ck_reac)
# change the angle x label
p_reac + theme(axis.text.x=element_text(angle=90, hjust=1))


