# ## This script aims to highlight the attenuation behaviour of p53 targets in BCL-2 over-expression
# Plotting Eulero Venn Diagram

# Author: Saveria Mazzara

# load libraries
library(venn)
library(readr)
library(dplyr)
library(viridis)
library(ggplot2)
library(pheatmap)
library(readxl)
#library(ggalt)


len<-10
DE_list <- vector(mode = "list", length = len)

DE_list_filt<-vector(mode="list",length= len)


#1. 
temp_Ent_vs_BCL2nt <- read.csv("DE_analysis_filtered_Ent_vs_BCL2nt.csv")
temp_Ent_vs_BCL2nt<-filter(temp_Ent_vs_BCL2nt, abs(log2FoldChange)>log2(1.5))
temp_DE_Ent_vs_BCL2nt <- read.csv("DE_Ent_vs_BCL2nt.csv")
#2. 
temp_Ent_vs_Edoxo <- read.csv("DE_analysis_filtered_Ent_vs_Edoxo.csv")
temp_Ent_vs_Edoxo<-filter(temp_Ent_vs_Edoxo, abs(log2FoldChange)>log2(1.5))
temp_DE_Ent_vs_Edoxo <- read.csv("DE_Ent_vs_Edoxo.csv")
#3. 
temp_Ent_vs_Eactd <- read.csv("DE_analysis_filtered_Ent_vs_Eactd.csv")
temp_Ent_vs_Eactd <-filter(temp_Ent_vs_Eactd, abs(log2FoldChange)>log2(1.5))
temp_DE_Ent_vs_Eactd <- read.csv("DE_Ent_vs_Eactd.csv")
#4. 
temp_Ent_vs_Ecx <- read.csv("DE_analysis_filtered_Ent_vs_Ecx.csv")
temp_Ent_vs_Ecx <-filter(temp_Ent_vs_Ecx, abs(log2FoldChange)>log2(1.5))
temp_DE_Ent_vs_Ecx <- read.csv("DE_Ent_vs_Ecx.csv")
#5.
temp_BCL2nt_vs_BCL2doxo <- read.csv("DE_analysis_filtered_BCL2nt_vs_BCL2doxo.csv")
temp_BCL2nt_vs_BCL2doxo <-filter(temp_BCL2nt_vs_BCL2doxo, abs(log2FoldChange)>log2(1.5))
temp_DE_BCL2nt_vs_BCL2doxo <- read.csv("DE_BCL2nt_vs_BCL2doxo.csv")
#6. 
temp_BCL2nt_vs_BCL2actd <- read.csv("DE_analysis_filtered_BCL2nt_vs_BCL2actd.csv")
temp_BCL2nt_vs_BCL2actd <-filter(temp_BCL2nt_vs_BCL2actd, abs(log2FoldChange)>log2(1.5))
temp_DE_BCL2nt_vs_BCL2actd <- read.csv("DE_BCL2nt_vs_BCL2actd.csv")
#7. 
temp_BCL2nt_vs_BCL2cx <- read.csv("DE_analysis_filtered_BCL2nt_vs_BCL2cx.csv")
temp_BCL2nt_vs_BCL2cx <-filter(temp_BCL2nt_vs_BCL2cx, abs(log2FoldChange)>log2(1.5))
temp_DE_BCL2nt_vs_BCL2cx <- read.csv("DE_BCL2nt_vs_BCL2cx.csv")
#8. 
temp_Edoxo_vs_BCL2doxo <- read.csv("DE_analysis_filtered_Edoxo_vs_BCL2doxo.csv")
temp_Edoxo_vs_BCL2doxo <-filter(temp_Edoxo_vs_BCL2doxo, abs(log2FoldChange)>log2(1.5))
temp_DE_Edoxo_vs_BCL2doxo <- read.csv("DE_Edoxo_vs_BCL2doxo.csv")
#9. 
temp_Eactd_vs_BCL2actd <- read.csv("DE_analysis_filtered_Eactd_vs_BCL2actd.csv")
temp_Eactd_vs_BCL2actd <-filter(temp_Eactd_vs_BCL2actd, abs(log2FoldChange)>log2(1.5))
temp_DE_Eactd_vs_BCL2actd <- read.csv("DE_Eactd_vs_BCL2actd.csv")
#10. 
temp_Ecx_vs_BCL2cx <- read.csv("DE_analysis_filtered_Ecx_vs_BCL2cx.csv")
temp_Ecx_vs_BCL2cx <-filter(temp_Ecx_vs_BCL2cx, abs(log2FoldChange)>log2(1.5))
temp_DE_Ecx_vs_BCL2cx <- read.csv("DE_Ecx_vs_BCL2cx.csv")

# List: DE genes
DE_list_filt[[1]]<-temp_Ent_vs_BCL2nt$gene_name
DE_list_filt[[2]]<-temp_Ent_vs_Edoxo$gene_name
DE_list_filt[[3]]<-temp_Ent_vs_Eactd$gene_name
DE_list_filt[[4]]<-temp_Ent_vs_Ecx$gene_name
DE_list_filt[[5]]<-temp_BCL2nt_vs_BCL2doxo$gene_name
DE_list_filt[[6]]<-temp_BCL2nt_vs_BCL2actd$gene_name
DE_list_filt[[7]]<-temp_BCL2nt_vs_BCL2cx$gene_name
DE_list_filt[[8]]<-temp_Edoxo_vs_BCL2doxo$gene_name
DE_list_filt[[9]]<-temp_Eactd_vs_BCL2actd$gene_name
DE_list_filt[[10]]<-temp_Ecx_vs_BCL2cx$gene_name

names(DE_list_filt)<-c("Ent_vs_BCL2nt","Ent_vs_Edoxo","Ent_vs_Eactd","Ent_vs_Ecx","BCL2nt_vs_BCL2doxo",
                       "BCL2nt_vs_BCL2actd","BCL2nt_vs_BCL2cx","Edoxo_vs_BCL2doxo","Eactd_vs_BCL2actd","Ecx_vs_BCL2cx")


# Venn:plotting Venn diagram
venn(DE_list_filt[2:7],col="black",zcolor=viridis(6),ilcs = 0.8)

# Empty
venn(DE_list_filt[2:4],col="black",zcolor=viridis(6),ilcs = 1.4,sncs=1.4,borders=FALSE)

# BCL2
venn(DE_list_filt[5:7],col="black",zcolor=viridis(6),ilcs = 1.4,sncs=1.4,borders=FALSE)

#### Modulation 
# Merging Lists

# Modulation BCL2
# common between the 3 groups
common_genes<-c(intersect(DE_list_filt[[7]],intersect(DE_list_filt[[5]],DE_list_filt[[6]])))
# common between doxo and actd
common_doxo_actd<-setdiff(intersect(DE_list_filt[[5]],DE_list_filt[[6]]),common_genes)
# common doxo and cx
common_doxo_cx<-setdiff(intersect(DE_list_filt[[5]],DE_list_filt[[7]]),common_genes)
# common actd and cx
common_actd_cx<-setdiff(intersect(DE_list_filt[[6]],DE_list_filt[[7]]),common_genes)

common<-c(common_genes,common_doxo_actd,common_doxo_cx,common_actd_cx)
genes_4_heat<-c(setdiff(DE_list_filt[[5]],common),
                setdiff(DE_list_filt[[6]],common),
                setdiff(DE_list_filt[[7]],common),
                common)

# Modulation Empty
# common between the 3 groups
common_genes_em<-c(intersect(DE_list_filt[[4]],intersect(DE_list_filt[[2]],DE_list_filt[[3]])))
# common between doxo and actd
common_doxo_actd_em<-setdiff(intersect(DE_list_filt[[2]],DE_list_filt[[3]]),common_genes_em)
# common doxo and cx
common_doxo_cx_em<-setdiff(intersect(DE_list_filt[[2]],DE_list_filt[[4]]),common_genes_em)
# common actd and cx
common_actd_cx_em<-setdiff(intersect(DE_list_filt[[3]],DE_list_filt[[4]]),common_genes_em)

common_em<-c(common_genes_em,common_doxo_actd_em,common_doxo_cx_em,common_actd_cx_em)
genes_4_heat_em<-c(setdiff(DE_list_filt[[2]],common_em),
                   setdiff(DE_list_filt[[3]],common_em),
                   setdiff(DE_list_filt[[4]],common_em),
                   common_em)


# 2523 genes in common
length(intersect(genes_4_heat,genes_4_heat_em))


list_bcl2<-genes_4_heat      # 3895
list_empty<-genes_4_heat_em  # 4083

# 1560 Empty, 2523 Common, 1372 BCL2
common_overall<-intersect(list_bcl2,list_empty)
spec_empty<-setdiff(list_empty,list_bcl2)
spec_bcl2<-setdiff(list_bcl2,list_empty)

# all DE genes: 5455 genes (aka Consensus DEGs list)
# genes_4_plot<-c(spec_empty,spec_bcl2,common_overall)
# prepare format for excel
#genes_4_plot<-data.frame(Gene=genes_4_plot)
#write.xlsx(genes_4_plot,file="Consensus_DEGs_list.xlsx")

# TRANSFAC: transcription target DE genes 
# Load the TRANSFAC list
TP53_targets<- as.data.frame(read_excel("local_data/TP53_targets_transfac.xlsx",  col_names = FALSE))
colnames(TP53_targets)<-"gene_name"

tp53_transfac_spec_empty<-(intersect(spec_empty,TP53_targets$gene_name))
tp53_transfac_spec_bcl2<-(intersect(spec_bcl2,TP53_targets$gene_name))
tp53_transfac_common<-(intersect(common_overall,TP53_targets$gene_name))




