# This script aims to highlight the attenuation behaviour of p53 targets in BCL-2 over-expression
# (this script produce the Fig5)


# Ploting Eulero Venn Diagram

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

# DE genes: barplot

df_ngenes<-data.frame(comp=rep(names(DE_list_filt),each=2),group=c(rep(c("up","down"),10)),
                      value=c(145,16,2075,863,1392,875,378,111,1649,1085,1173,1036,270,45,639,891,110,59,172,102),
                      label_y=c(145,250,2075,2938,1392,2267,378,550,1649,2734,1173,2209,270,370,639,1530,110,220,172,300))
# fix order
# df_ngenes$comp <- factor(df_ngenes$comp,levels = c(df_ngenes$comp))

ggplot(df_ngenes, aes(fill=group, y=value, x=comp)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip()+
  #scale_fill_viridis(discrete = T) +
  scale_fill_manual(values=c("#B4D6C1","#358873"))+
  ggtitle("# DE genes") +
  xlab("comparisons")+
  ylab("n DE genes")+
  theme(text = element_text(size=20))+
  geom_text(aes(y = label_y, label = value), vjust = 1.5, colour = "black")


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
TP53_targets<- as.data.frame(read_excel("~/derenzini/datasets/TP53_targets_transfac.xlsx",  col_names = FALSE))
colnames(TP53_targets)<-"gene_name"

tp53_transfac_spec_empty<-(intersect(spec_empty,TP53_targets$gene_name))
tp53_transfac_spec_bcl2<-(intersect(spec_bcl2,TP53_targets$gene_name))
tp53_transfac_common<-(intersect(common_overall,TP53_targets$gene_name))

# FISHCER: transcription target DE genes
#tp53_fischer_spec_empty<-(intersect(spec_empty,TP53_targets_f$gene_name))
#tp53_fischer_spec_bcl2<-(intersect(spec_bcl2,TP53_targets_f$gene_name))
#tp53_fischer_common<-(intersect(common_overall,TP53_targets_f$gene_name))

genes_4_plot<-c( tp53_transfac_spec_empty,tp53_transfac_spec_bcl2,tp53_transfac_common)
# genes_4_plot<-c( tp53_fischer_spec_empty,tp53_fischer_spec_bcl2,tp53_fischer_common)

df_plot<-data.frame(gene_name=c(tp53_transfac_spec_empty,tp53_transfac_spec_bcl2, tp53_transfac_common),
                    comparison=c(rep("empty",length(tp53_transfac_spec_empty)),rep("bcl2",length(tp53_transfac_spec_bcl2)),rep("common",length(tp53_transfac_common))),
                    fc_doxo_empty=temp_DE_Ent_vs_Edoxo[match(genes_4_plot,temp_DE_Ent_vs_Edoxo$gene_name),3],
                    fc_actd_empty=temp_DE_Ent_vs_Eactd[match(genes_4_plot,temp_DE_Ent_vs_Eactd$gene_name),3],
                    fc_cx_empty=temp_DE_Ent_vs_Ecx[match(genes_4_plot,temp_DE_Ent_vs_Ecx$gene_name),3],
                    fc_doxo_bcl2=temp_DE_BCL2nt_vs_BCL2doxo[match(genes_4_plot,temp_DE_BCL2nt_vs_BCL2doxo$gene_name),3],
                    fc_actd_bcl2=temp_DE_BCL2nt_vs_BCL2actd[match(genes_4_plot,temp_DE_BCL2nt_vs_BCL2actd$gene_name),3],
                    fc_cx_bcl2=temp_DE_BCL2nt_vs_BCL2cx[match(genes_4_plot,temp_DE_BCL2nt_vs_BCL2cx$gene_name),3])

# rename column
colnames(df_plot)<-c("gene_name", "comparison" ,"fc_doxo_empty", "fc_actd_empty", "fc_cx_empty","fc_doxo_bcl2","fc_actd_bcl2","fc_cx_bcl2")  

# change NaN in NA
df_plot[df_plot == "NaN" ] <- NA

rownames(df_plot)<-df_plot$gene_name
# save data
#write.xlsx(df_plot,file="modulation_heatmap_targets_fischer_rnaseq_fc.xlsx")
#write.xlsx(df_plot,file="modulation_heatmap_targets_fischer_rnaseq_fc.xlsx")

### Heatmap

# annotation rows
n_empty<-length(tp53_transfac_spec_empty)
n_bcl2<-length(tp53_transfac_spec_bcl2)
n_overall<-length(tp53_transfac_common)


annotation_row_class=data.frame(GeneClass=c(rep("Empty",n_empty),rep("BCL2",n_bcl2),rep("Common",n_overall)  )) 
rownames(annotation_row_class)<-rownames(df_plot)

# color code
annotation_colors = list(
  GeneClass=c(Empty="#F48024",BCL2="#18A3AC",Common="#B4B9BF"))  # DMSO="#218B82",

my_palette <- colorRampPalette(c("#D4070F", "red","ivory","blue" ,"#1C1AAF"))(n = 99)

paletteLength <- 50
my_palette <- colorRampPalette(c("#D4070F", "red","ivory","blue" ,"#1C1AAF"))(paletteLength)
myBreaks <- c(seq(min(na.omit(df_plot[,-c(1:2)])),0,  length.out=ceiling(paletteLength/2) + 1),
              seq(max(max(na.omit(df_plot[,-c(1:2)])))/paletteLength, max(na.omit(df_plot[,-c(1:2)])), length.out=floor(paletteLength/2)))


# Overall heatmap

pheatmap(((df_plot[,-c(1:2)])),scale="none",
         cluster_rows = F,
         cluster_cols=F,
         row.names=T,
         annotation_row = annotation_row_class,
         annotation_colors = annotation_colors,
         gaps_col = c(3),
         color=rev(my_palette),
         breaks=myBreaks,
         fontsize_col = 8,
         fontsize_row = 3,
         #cellheight = 1,
         cellwidth = 10,
         border_color=NA,
         na_col="grey54"
         #clustering_distance_rows="correlation",
         #clustering_distance_cols="correlation",
         #clustering_method = "average"
)



##### ACTD
df_actd<-df_plot


ggplot(df_actd,aes(x = fc_actd_bcl2, y =fc_actd_empty)) +
  geom_point(color="#8EA4C8", size=4) + 
  # geom_point(aes(color=factor(df_actd$comparison)),  size=4) +    
  geom_abline(slope = 1, intercept = 0, alpha = 0.6) +
  ggtitle("ACTD treatment") +
  geom_text(
    size=2, hjust=0.001,
    data= filter(df_actd, ( (abs(fc_actd_empty) > log2(1.5) | abs(fc_actd_bcl2) > log2(1.5)) & ((fc_actd_bcl2-fc_actd_empty)>=1 | (fc_actd_empty-fc_actd_bcl2)>=1 ))), # Filter data first
    aes(label=gene_name),
    nudge_x = 0.05, nudge_y = 0.05,
  ) +
  # scale_color_manual(values=c("#18A3AC", "#B4B9BF", "#F48024"))+
  # scale_colour_manual(values = c("#18A3AC", "#B4B9BF", "#F48024","#ffb577","#ffd9aa","#90BD31","#D3E4AD"))+
  xlab(" log2 FC bcl2")+
  ylab(" log2 FC empty")+
  labs(color='Gene Type')+
  theme(plot.title = element_text(size=20,face = "bold.italic"),
        axis.text.x = element_text(size=20,color = "black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.title.x =element_text(size=20,color="black"),
        axis.title.y=element_text(size=20,color="black"),
        legend.text = element_text(size=20),
        legend.title= element_text(size=20))

# Delta FC
delta_fc_actd<-data.frame(delta=df_actd$fc_actd_bcl2-df_actd$fc_actd_empty,gene_name=df_actd$gene_name,
                          fc_empty=df_actd$fc_actd_empty,
                          fc_bcl2=df_actd$fc_actd_bcl2,
                          comparison=df_actd$comparison,treatment=rep("ACTD",dim(df_actd)[1]))


ggplot(delta_fc_actd,aes(y=delta,x=treatment))+
  #geom_jitter(aes(color=factor(delta_fc_actd$comparison)),  size=3,position = position_jitter(seed = 1))+
  geom_jitter(color="#8EA4C8", size=2.5, position = position_jitter(seed = 1))+  
  geom_hline(yintercept= -log2(1.5), linetype="dashed", color = "red")+
  geom_hline(yintercept= log2(1.5), linetype="dashed", color = "red")+
  ggtitle("ACTD treatment") +
  geom_text(size=3,
            position = position_jitter(seed = 1),
            label=delta_fc_actd$gene_name) 
#   scale_color_manual(values=c("#18A3AC", "#B4B9BF", "#F48024"))

##### DOXO

df_doxo<-df_plot

ggplot(df_doxo,aes(x = fc_doxo_bcl2, y =fc_doxo_empty)) +
  geom_point(color="#8EA4C8", size=4) +       #aes(color=factor(df_doxo$comparison))
  geom_abline(slope = 1, intercept = 0, alpha = 0.6) +
  ggtitle("DOXO treatment") +
  geom_text(
    size=5, hjust=0.001,
    data= filter(df_doxo, ( (abs(fc_doxo_empty) > log2(2) | abs(fc_doxo_bcl2) > log2(2)) & ((fc_doxo_bcl2-fc_doxo_empty)>=1 | (fc_doxo_empty-fc_doxo_bcl2)>=1 ))), # Filter data first
    aes(label=gene_name),
    nudge_x = 0.05, nudge_y = 0.05,
  ) +
  # scale_color_manual(values=c("#18A3AC", "#B4B9BF", "#F48024"))+
  # scale_colour_manual(values = c("#18A3AC", "#B4B9BF", "#F48024","#ffb577","#ffd9aa","#90BD31","#D3E4AD"))+
  xlab(" log2 FC bcl2")+
  ylab(" log2 FC empty")+
  labs(color='Gene Type')+
  theme(plot.title = element_text(size=20,face = "bold.italic"),
        axis.text.x = element_text(size=20,color = "black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.title.x =element_text(size=20,color="black"),
        axis.title.y=element_text(size=20,color="black"),
        legend.text = element_text(size=20),
        legend.title= element_text(size=20))

# Delta FC
delta_fc_doxo<-data.frame(delta=df_doxo$fc_doxo_bcl2-df_doxo$fc_doxo_empty,gene_name=df_doxo$gene_name,
                          fc_empty=df_doxo$fc_doxo_empty,                          fc_bcl2=df_doxo$fc_doxo_bcl2,
                          comparison=df_doxo$comparison,treatment=rep("DOXO",dim(df_doxo)[1]))


ggplot(delta_fc_doxo,aes(y=delta,x=treatment))+
  geom_jitter(color="#8EA4C8", size=2.5, position = position_jitter(seed = 1))+  
  geom_hline(yintercept= -log2(1.5), linetype="dashed", color = "red")+
  geom_hline(yintercept= log2(1.5), linetype="dashed", color = "red")+
  ggtitle("DOXO treatment") +
  geom_text(size=3,
            position = position_jitter(seed = 1),
            label=delta_fc_doxo$gene_name) 



##### CX

df_cx<-df_plot

ggplot(df_cx,aes(x = fc_cx_bcl2, y =fc_cx_empty)) +
  geom_point(color="#8EA4C8", size=4) +       #aes(color=factor(df_doxo$comparison))
  geom_abline(slope = 1, intercept = 0, alpha = 0.6) +
  ggtitle("CX treatment") +
  xlim(c(-1, 2)) + 
  geom_text(
    size=5, hjust=0.001,
    data= filter(df_doxo, ( (abs(fc_cx_empty) > log2(1.5) | abs(fc_cx_bcl2) > log2(1.5)) & ((fc_cx_bcl2-fc_cx_empty)>=1 | (fc_cx_empty-fc_cx_bcl2)>=1 ))), # Filter data first
    aes(label=gene_name),
    nudge_x = 0.05, nudge_y = 0.05,
  ) +
  # scale_color_manual(values=c("#18A3AC", "#B4B9BF", "#F48024"))+
  # scale_colour_manual(values = c("#18A3AC", "#B4B9BF", "#F48024","#ffb577","#ffd9aa","#90BD31","#D3E4AD"))+
  xlab(" log2 FC bcl2")+
  ylab(" log2 FC empty")+
  labs(color='Gene Type')+
  theme(plot.title = element_text(size=20,face = "bold.italic"),
        axis.text.x = element_text(size=20,color = "black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.title.x =element_text(size=20,color="black"),
        axis.title.y=element_text(size=20,color="black"),
        legend.text = element_text(size=20),
        legend.title= element_text(size=20))

# Delta FC
delta_fc_cx<-data.frame(delta=df_cx$fc_cx_bcl2-df_cx$fc_cx_empty,gene_name=df_cx$gene_name,
                        fc_empty=df_cx$fc_cx_empty,
                        fc_bcl2=df_cx$fc_cx_bcl2,
                        comparison=df_cx$comparison,treatment=rep("CX",dim(df_cx)[1]))


ggplot(delta_fc_cx,aes(y=delta,x=treatment))+
  geom_jitter(color="#8EA4C8", size=2.5, position = position_jitter(seed = 1))+  
  geom_hline(yintercept= -log2(1.5), linetype="dashed", color = "red")+
  geom_hline(yintercept= log2(1.5), linetype="dashed", color = "red")+
  ggtitle("CX treatment") +
  geom_text(size=3,
            position = position_jitter(seed = 1),
            label=delta_fc_cx$gene_name) 


target_actd<-delta_fc_actd[(which(abs(delta_fc_actd$delta) >= log2(1.5))),]
target_actd<-target_actd[order(target_actd$delta),]

target_doxo<-delta_fc_doxo[(which(abs(delta_fc_doxo$delta) >= log2(1.5))),]
target_doxo<-target_doxo[order(target_doxo$delta),]

target_cx<-delta_fc_cx[(which(abs(delta_fc_cx$delta) >= log2(1.5))),]
target_cx<-target_cx[order(target_cx$delta),]



###### Dumbbell Plot
# ref: https://towardsdatascience.com/create-dumbbell-plots-to-visualize-group-differences-in-r-3536b7d0a19a
df_dumb_actd<-delta_fc_actd[match(target_actd$gene_name,delta_fc_actd$gene_name),-c(5,6)]
df_dumb_actd<-df_dumb_actd[order(abs(df_dumb_actd$delta)),]
df_dumb_actd<-df_dumb_actd[order((df_dumb_actd$delta)),]
#df_dumb_actd<-df_dumb_actd[order((df_dumb_actd$delta)),]

write.xlsx(df_dumb_actd,file="p53targets_actd_rnaseq_fc.xlsx")


specie <- c(rep("LIF" , 2) , rep("TTLL6",2) , rep("LRRC39",2) , rep("TNFRSF9",2),rep( "RNF224",2),
            rep( "DHDH",2),rep( "CDKN1A",2),rep("RNF151",2),rep("KIFC3",2),rep("SPOCK1",2),rep("RBM47",2),
            rep("CUEDC1",2),rep("NAV1",2),rep("IER3",2),rep("SLC25A32",2),rep("IRF7",2),rep("CCDC120",2),
            rep("CKB",2),rep("MICALL2",2),rep("GALR2",2),rep("SDSL",2),rep("IL17RD",2),rep("NCKAP1",2),   
            rep("PRX",2),rep("KSR1",2),rep("BCAS3",2),rep("CCDC62",2),rep("HRK",2),rep("DNAJC5G",2)) 
condition <- rep(c("bbempty" , "bcl2" ) , 29)
value<-unlist(c(df_dumb_actd[1,3:4],df_dumb_actd[2,3:4],df_dumb_actd[3,3:4],df_dumb_actd[4,3:4],df_dumb_actd[5,3:4],
                df_dumb_actd[6,3:4],df_dumb_actd[7,3:4],df_dumb_actd[8,3:4],df_dumb_actd[9,3:4],df_dumb_actd[10,3:4],
                df_dumb_actd[11,3:4],df_dumb_actd[12,3:4],df_dumb_actd[13,3:4],df_dumb_actd[14,3:4],df_dumb_actd[15,3:4],
                df_dumb_actd[16,3:4],df_dumb_actd[17,3:4],df_dumb_actd[18,3:4],df_dumb_actd[19,3:4],df_dumb_actd[20,3:4],
                df_dumb_actd[21,3:4],df_dumb_actd[22,3:4],df_dumb_actd[23,3:4],df_dumb_actd[24,3:4],df_dumb_actd[25,3:4],
                df_dumb_actd[26,3:4],df_dumb_actd[27,3:4],df_dumb_actd[28,3:4],df_dumb_actd[29,3:4]
                ))

dataprova <- data.frame(specie,condition,value)
dataprova$specie <- factor(dataprova$specie,
                                 levels=c('LIF',  'TTLL6', 'LRRC39','TNFRSF9','RNF224','DHDH','CDKN1A','RNF151','KIFC3','SPOCK1',
                                          'RBM47','CUEDC1','NAV1','IER3','SLC25A32','IRF7','CCDC120','CKB','MICALL2','GALR2',
                                          'SDSL','IL17RD','NCKAP1','PRX', 'KSR1','BCAS3','CCDC62','HRK','DNAJC5G'
                                 ),
                                 ordered = TRUE)
#dataprova$condition<-factor(dataprova$condition,
#                            levels=c('bbempty','bcl2'),ordered=TRUE)


ggplot(dataprova, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="dodge", stat="identity")+theme(axis.text.x=element_text(angle=45, hjust=0.9))+ylab("fc_empty/fc_bcl2")

paletteLength <- 50
my_palette <- colorRampPalette(c("#D4070F", "red","ivory","blue" ,"#1C1AAF"))(paletteLength)
myBreaks <- c(seq(min(na.omit(df_dumb_actd[,-c(1:2)])),0,  length.out=ceiling(paletteLength/2) + 1),
              seq(max(max(na.omit(df_dumb_actd[,-c(1:2)])))/paletteLength, max(na.omit(df_dumb_actd[,-c(1:2)])), length.out=floor(paletteLength/2)))


# Overall heatmap

pheatmap(((df_dumb_actd[,-c(1:2)])),scale="none",
         cluster_rows = F,
         cluster_cols=F,
         row.names=T,
         #annotation_row = annotation_row_class,
         #annotation_colors = annotation_colors,
         #gaps_col = c(3),
         color=rev(my_palette),
         breaks=myBreaks,
         fontsize_col = 8,
         fontsize_row = 3,
         #cellheight = 1,
         cellwidth = 10,
         border_color=NA,
         na_col="grey54"
         #clustering_distance_rows="correlation",
         #clustering_distance_cols="correlation",
         #clustering_method = "average"
)




# set the order 
# df_dumb_actd$gene_name<-factor(df_dumb_actd$gene_name,
#                                levels=(df_dumb_actd$gene_name))
# df_dumb_actd$gene_name <- factor(df_dumb_actd$gene_name,
#                                 levels=c('HRK','CCDC62','DNAJC5G',
#                                          'SLC25A32','SPOCK1','BCAS3','SDSL','RNF151','NAV1', 'KSR1','IL17RD', 'PRX',  'IRF7','TNFRSF9','DHDH',
#                                          'CCDC120','CKB','CUEDC1','KIFC3','MICALL2','IER3', 'NCKAP1','RBM47','GALR2','LRRC39',
#                                        'CDKN1A','RNF224','TTLL6','LIF'),
#                                 ordered = TRUE)
df_dumb_actd$gene_name <- factor(df_dumb_actd$gene_name,
                                 levels=c('BCAS3','KSR1','PRX', 'NCKAP1','IL17RD','SDSL','GALR2','MICALL2','CKB','CCDC120','IRF7',
                                          'IER3','NAV1','CUEDC1','HRK','RBM47','KIFC3','RNF151','CDKN1A','DHDH','RNF224','TNFRSF9','LRRC39',
                                          'TTLL6','LIF',
                                          'CCDC62', 'SLC25A32','SPOCK1', 'DNAJC5G' 
                                 ),
                                 ordered = TRUE)


ggplot(df_dumb_actd,aes(y=gene_name,x=fc_empty,xend=fc_bcl2))+
  geom_dumbbell(size=1,color="gray24",size_x=4,size_xend=4,colour_x="seagreen",colour_xend="purple",alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()+
  labs(x="delta fold change ")+
  geom_hline(yintercept=c(11.5,25.5,27.5), linetype="dashed", 
             color = "#84A6D6", size=1)+
  geom_hline(yintercept=c(25.5),
             color = "#4382BB", size=1.3)


geom_hline(yintercept=c(6.5), linetype="3313", color = "#84A6D6", size=1)+
  geom_hline(yintercept=c(12.5), color = "#4382BB", size=1.3)

# geom_text(data=filter(df_dumb_actd, gene_name=="TTLL6"),
#           aes(x=fc_bcl2, y=gene_name, label="fc bcl2"),
#           color="#b8627dff", size=3, vjust=-1.5, fontface="bold") + coord_flip()
#    geom_text(data=filter(df_dumb_actd, gene_name=="TTLL6"),
#              aes(x=fc_empty, y=gene_name, label="fc empty"),
#              color="#eb8055ff", size=3, vjust=-1.5, fontface="bold")   

#geom_dumbbell(data=infected, aes(y=concerned, x=fc_empty, xend=fc_bcl2), size=1.5, color="#b2b2b2", size_x=3, size_xend = 3, colour_x = red, colour_xend = blue)

# DOXO
df_dumb_doxo<-delta_fc_doxo[match(target_doxo$gene_name,delta_fc_doxo$gene_name),-c(5,6)]
df_dumb_doxo<-df_dumb_doxo[order(abs(df_dumb_doxo$delta)),]  
df_dumb_doxo<-df_dumb_doxo[order((df_dumb_doxo$delta)),] 

write.xlsx(df_dumb_doxo,file="p53targets_doxo_rnaseq_fc.xlsx")


# set the order 
# df_dumb_doxo$gene_name <- factor(df_dumb_doxo$gene_name,
#                                  levels=c('HRK','YWHAH','VPS72','GTPBP4','SNRPB','HNRNPF','SLC3A2','CENPM','ZNF408','RFFL',
#                                           'ABHD14A','ZSCAN21','DDX58','OSER1','PSMC3IP','CTNNBIP1','AEN','LZTS3','DNAJC5G',
#                                           'CD80','GALR2',
#                                     'SLC23A1','MEF2B','BNIPL','CPT1B','SYVN1','RFXAP','SPOCK1','SLC35C1','CHRNA10','NDUFV2',
#                                           'RAI1','SIPA1L3','SDSL','NAV1', 'IRF7','CCDC120','CKB','GADD45B','CUEDC1','MICALL2','FOS',
#                                           'KIFC3','DHDH','RBM47','LRRC39','RNF151','SLC35G2','RPL3L','RNF224','CDKN1A','TTLL6','LIF'
#                                           ),
#                                  ordered = TRUE)


df_dumb_doxo$gene_name<-factor(df_dumb_doxo$gene_name,
                               levels=c('RAI1','SDSL','YWHAH','GTPBP4','HNRNPF','VPS72','CCDC120','IRF7','FOS','SNRPB','MICALL2',
                                        'DHDH','HRK','KIFC3','CUEDC1','RBM47','SIPA1L3','CDKN1A','CKB','NAV1','RNF224','GADD45B',
                                        'LIF','SLC35G2','RPL3L','RNF151', 'TTLL6','LRRC39',   
                                        'CTNNBIP1','CENPM','RFXAP','SLC23A1','OSER1','DDX58','CPT1B','SLC3A2','ZSCAN21','SYVN1','RFFL',
                                        'PSMC3IP','AEN','NDUFV2','ABHD14A','DNAJC5G','MEF2B','GALR2','ZNF408','SLC35C1','CHRNA10','LZTS3',
                                        'CD80','BNIPL','SPOCK1'
                               ),
                               ordered = TRUE)

#df_dumb_doxo$gene_name<-factor(df_dumb_doxo$gene_name,
#                               levels=(df_dumb_doxo$gene_name))



ggplot(df_dumb_doxo,aes(y=gene_name,x=fc_empty,xend=fc_bcl2))+
  geom_dumbbell(size=1,color="gray24",size_x=4,size_xend=4,colour_x="seagreen",colour_xend="purple",alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()+
  labs(x="delta fold change ")+
  geom_hline(yintercept=c(16.5,49.5), linetype="3313", 
             color = "#84A6D6", size=1)+
  geom_hline(yintercept=c(28.5),
             color = "#4382BB", size=1.3)



# CX
df_dumb_cx<-delta_fc_cx[match(target_cx$gene_name,delta_fc_cx$gene_name),-c(5,6)]
df_dumb_cx<-df_dumb_cx[order(abs(df_dumb_cx$delta)),]
df_dumb_cx<-df_dumb_cx[order((df_dumb_cx$delta)),]


write.xlsx(df_dumb_cx,file="p53targets_cx_rnaseq_fc.xlsx")


# set the order 
# df_dumb_cx$gene_name<-factor(df_dumb_cx$gene_name,
#                              levels=(df_dumb_cx$gene_name))
# df_dumb_cx$gene_name <- factor(df_dumb_cx$gene_name,
#                                  levels=c('DDX58',
#                                           'CMTM8','SPOCK1','NAV1','MICALL2','CUEDC1','IER3','NCKAP1','RBM47','MYRF','GALR2','RNF224',
#                                           'CDKN1A','LIF'),
#                                  ordered = TRUE)
df_dumb_cx$gene_name <- factor(df_dumb_cx$gene_name,
                               levels=c('GALR2','RBM47','SPOCK1','CUEDC1','NAV1','NCKAP1','IER3','MYRF',
                                        'MICALL2','CDKN1A','RNF224','LIF',
                                        'CMTM8','DDX58'),
                               ordered = TRUE)


ggplot(df_dumb_cx,aes(y=gene_name,x=fc_empty,xend=fc_bcl2))+
  geom_dumbbell(size=1,color="gray24",size_x=4,size_xend=4,colour_x="seagreen",colour_xend="purple",alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip() +
  labs(x="delta fold change ")  +
  geom_hline(yintercept=c(6.5), linetype="3313", color = "#84A6D6", size=1)+
  geom_hline(yintercept=c(12.5), color = "#4382BB", size=1.3)





#############################################################################################
####### Heatmaps on single p53 targets for each treatment with expression values (normcounts)
##### ACTD
# Heatmaps with replicates
target_actd<-delta_fc_actd[(which(abs(delta_fc_actd$delta) >= log2(1.5))),]

target_actd<-target_actd[order(target_actd$delta),]


# recover gene annotation
normalized.counts <- as.data.frame(counts(dds.norm, normalized=TRUE ))
gene_ids = rownames(normalized.counts)

#recover gene annotation for DE genes
normalized.counts$ensembl_gene_id = rownames(normalized.counts)

# annotation directly from transcriptome (gene_biotype.txt)
genes_info=gene_biotype[match(rownames(normalized.counts),gene_biotype$V2),]

# mantain only useful information
genes_info=genes_info[,c(2,4,6)]

colnames(genes_info)<-c("ensembl_gene_id","gene_type","gene_name")

# create data frame with gene name
norm.counts.anno <- merge (normalized.counts, genes_info, by ="ensembl_gene_id", all.x = TRUE)

# actd
df_targets_actd<-norm.counts.anno[match(target_actd$gene_name,norm.counts.anno$gene_name),]

df_targets_actd<-df_targets_actd[,c(2,10,18,4,12,20,6,14,22,8,16,24)]


rownames(df_targets_actd)<-norm.counts.anno$gene_name[match(target_actd$gene_name,norm.counts.anno$gene_name)]

# insert the annotation column
annotation_col = data.frame(Condition = factor(c(rep("empty_nt",3),rep("empty_actd",3),rep("bcl2_nt",3),rep("bcl2_actd",3))))
rownames(annotation_col)<-colnames(df_targets_actd)

# define the annotation
annotation_row = data.frame(
  delta_fc = factor(rep(c("neg", "pos"), c(26, 3))))

rownames(annotation_row) = rownames(df_targets_actd)
#color_annotation_row=c(neg="#C2D9E1",pos="#E5B3BB")

color_annotation_col = list(Condition = c(empty_nt ="#DCDBD9", empty_actd="#8FA2A6",bcl2_nt="#D5E4C3",bcl2_actd="#98D4BB"),
                            delta_fc=c(neg="#84A6D6",pos="#DC828F"))

#my_palette <- colorRampPalette(c("purple", "gray54" ,"yellow"))(n = 29)
my_palette <- colorRampPalette(c("white", "ivory","gold" ,"indianred","purple"))(n = 29)

pheatmap(df_targets_actd,scale="row",
         cluster_rows =F,
         cluster_cols=F,
         row.names=T,
         gaps_col = c(3,6,9),
         gaps_row = c(26),
         annotation_col = annotation_col,
         annotation_colors=color_annotation_col,
         annotation_row=annotation_row,
         # color=(my_palette),
         border_color=NA,
         na_col="grey54")

# check the expression for the 29 genes
pdf(file="targets_genes_actd.pdf",onefile=TRUE)  

# boxplots for the 29 actd targets
for (i in 1:dim(df_targets_actd)[1]){
  z<-data.frame(value=t(df_targets_actd[i,]),sample=colnames(df_targets_actd))
  colnames(z)[1]<-"value"
  z$sample<-c(rep("empty_nt",3),rep("empty_actd",3),rep("bcl2_nt",3),rep("bcl2_actd",3))
  z$sample <- factor(z$sample,levels =  c('empty_nt','empty_actd','bcl2_nt','bcl2_actd'),ordered = TRUE)
  p<-ggplot(z,aes(x=sample,y=value))+geom_boxplot(aes(col=sample))+
    geom_jitter(position=position_jitter(0.2))+
    ggtitle(rownames(df_targets_actd)[i])
  print(p)
}

dev.off() 



#### DOXO
# Heatmaps with replicates
target_doxo<-delta_fc_doxo[(which(abs(delta_fc_doxo$delta) >= log2(1.5))),]

target_doxo<-target_doxo[order(target_doxo$delta),]


# recover gene annotation
normalized.counts <- as.data.frame(counts(dds.norm, normalized=TRUE ))
gene_ids = rownames(normalized.counts)

#recover gene annotation for DE genes
normalized.counts$ensembl_gene_id = rownames(normalized.counts)

# annotation directly from transcriptome (gene_biotype.txt)
genes_info=gene_biotype[match(rownames(normalized.counts),gene_biotype$V2),]

# mantain only useful information
genes_info=genes_info[,c(2,4,6)]

colnames(genes_info)<-c("ensembl_gene_id","gene_type","gene_name")

# create data frame with gene name
norm.counts.anno <- merge (normalized.counts, genes_info, by ="ensembl_gene_id", all.x = TRUE)

# doxo
df_targets_doxo<-norm.counts.anno[match(target_doxo$gene_name,norm.counts.anno$gene_name),]

df_targets_doxo<-df_targets_doxo[,c(2,10,18,3,11,19,6,14,22,7,15,23)]


rownames(df_targets_doxo)<-norm.counts.anno$gene_name[match(target_doxo$gene_name,norm.counts.anno$gene_name)]

# insert the annotation column
annotation_col = data.frame(Condition = factor(c(rep("empty_nt",3),rep("empty_doxo",3),rep("bcl2_nt",3),rep("bcl2_doxo",3))))
rownames(annotation_col)<-colnames(df_targets_doxo)

# define the annotation
annotation_row = data.frame(
  delta_fc = factor(rep(c("neg", "pos"), c(32, 21))))

rownames(annotation_row) = rownames(df_targets_doxo)
#color_annotation_row=c(neg="#C2D9E1",pos="#E5B3BB")

color_annotation_col = list(Condition = c(empty_nt ="#DCDBD9", empty_doxo="#8FA2A6",bcl2_nt="#D5E4C3",bcl2_doxo="#98D4BB"),
                            delta_fc=c(neg="#84A6D6",pos="#DC828F"))


#neg="#C2D9E1",pos="#E5B3BB")
#84A6D6
#A2C4C6

#my_palette <- colorRampPalette(c("purple", "gray54" ,"yellow"))(n = 29)
# my_palette <- colorRampPalette(c("white", "ivory","gold" ,"indianred","purple"))(n = 29)

pheatmap(df_targets_doxo,scale="row",
         cluster_rows =F,
         cluster_cols=F,
         row.names=T,
         gaps_col = c(3,6,9),
         gaps_row = c(32),
         annotation_col = annotation_col,
         annotation_colors=color_annotation_col,
         annotation_row=annotation_row,
         # color=(my_palette),
         # breaks=myBreaks,
         # fontsize_col = 8,
         #  fontsize_row = 3,
         #cellheight = 1,
         #  cellwidth = 10,
         border_color=NA,
         na_col="grey54")

# check the expression for the 53 genes
pdf(file="targets_genes_doxo.pdf",onefile=TRUE)  

# boxplots for the 53 doxo targets
for (i in 1:dim(df_targets_doxo)[1]){
  z<-data.frame(value=t(df_targets_doxo[i,]),sample=colnames(df_targets_doxo))
  colnames(z)[1]<-"value"
  z$sample<-c(rep("empty_nt",3),rep("empty_doxo",3),rep("bcl2_nt",3),rep("bcl2_doxo",3))
  z$sample <- factor(z$sample,levels =  c('empty_nt','empty_doxo','bcl2_nt','bcl2_doxo'),ordered = TRUE)
  p<-ggplot(z,aes(x=sample,y=value))+geom_boxplot(aes(col=sample))+
    geom_jitter(position=position_jitter(0.2))+
    ggtitle(rownames(df_targets_doxo)[i])
  print(p)
}

dev.off() 


#### CX
# Heatmaps with replicates
target_cx<-delta_fc_cx[(which(abs(delta_fc_cx$delta) >= log2(1.5))),]

target_cx<-target_cx[order(target_cx$delta),]


# cx
df_targets_cx<-norm.counts.anno[match(target_cx$gene_name,norm.counts.anno$gene_name),]

df_targets_cx<-df_targets_cx[,c(2,10,18,5,13,21,6,14,22,9,17,25)]


rownames(df_targets_cx)<-norm.counts.anno$gene_name[match(target_cx$gene_name,norm.counts.anno$gene_name)]

# insert the annotation column
annotation_col = data.frame(Condition = factor(c(rep("empty_nt",3),rep("empty_cx",3),rep("bcl2_nt",3),rep("bcl2_cx",3))))
rownames(annotation_col)<-colnames(df_targets_cx)

# define the annotation
annotation_row = data.frame(
  delta_fc = factor(rep(c("neg", "pos"), c(13, 1))))

rownames(annotation_row) = rownames(df_targets_cx)

color_annotation_col = list(Condition = c(empty_nt ="#DCDBD9", empty_cx="#8FA2A6",bcl2_nt="#D5E4C3",bcl2_cx="#98D4BB"),
                            delta_fc=c(neg="#84A6D6",pos="#DC828F"))


#neg="#C2D9E1",pos="#E5B3BB")
#84A6D6
#A2C4C6

#my_palette <- colorRampPalette(c("purple", "gray54" ,"yellow"))(n = 29)
# my_palette <- colorRampPalette(c("white", "ivory","gold" ,"indianred","purple"))(n = 29)

pheatmap(df_targets_cx,scale="row",
         cluster_rows =F,
         cluster_cols=F,
         row.names=T,
         gaps_col = c(3,6,9),
         gaps_row = c(13),
         annotation_col = annotation_col,
         annotation_colors=color_annotation_col,
         annotation_row=annotation_row,
         # color=(my_palette),
         # breaks=myBreaks,
         # fontsize_col = 8,
         #  fontsize_row = 3,
         #cellheight = 1,
         #  cellwidth = 10,
         border_color=NA,
         na_col="grey54")

# check the expression for the 14 genes
pdf(file="targets_genes_cx.pdf",onefile=TRUE)  

# boxplots for the 53 doxo targets
for (i in 1:dim(df_targets_cx)[1]){
  z<-data.frame(value=t(df_targets_cx[i,]),sample=colnames(df_targets_cx))
  colnames(z)[1]<-"value"
  z$sample<-c(rep("empty_nt",3),rep("empty_cx",3),rep("bcl2_nt",3),rep("bcl2_cx",3))
  z$sample <- factor(z$sample,levels =  c('empty_nt','empty_cx','bcl2_nt','bcl2_cx'),ordered = TRUE)
  p<-ggplot(z,aes(x=sample,y=value))+geom_boxplot(aes(col=sample))+
    geom_jitter(position=position_jitter(0.2))+
    ggtitle(rownames(df_targets_cx)[i])
  print(p)
}

dev.off() 

#############################################################################################

###### Dumbbell Plot
# ref: https://towardsdatascience.com/create-dumbbell-plots-to-visualize-group-differences-in-r-3536b7d0a19a
df_dumb_actd<-delta_fc_actd[match(target_actd$gene_name,delta_fc_actd$gene_name),-c(5,6)]
df_dumb_actd<-df_dumb_actd[order(abs(df_dumb_actd$delta)),]
#df_dumb_actd<-df_dumb_actd[order((df_dumb_actd$delta)),]

# set the order 
# df_dumb_actd$gene_name<-factor(df_dumb_actd$gene_name,
#                                levels=(df_dumb_actd$gene_name))
# df_dumb_actd$gene_name <- factor(df_dumb_actd$gene_name,
#                                 levels=c('HRK','CCDC62','DNAJC5G',
#                                          'SLC25A32','SPOCK1','BCAS3','SDSL','RNF151','NAV1', 'KSR1','IL17RD', 'PRX',  'IRF7','TNFRSF9','DHDH',
#                                          'CCDC120','CKB','CUEDC1','KIFC3','MICALL2','IER3', 'NCKAP1','RBM47','GALR2','LRRC39',
#                                        'CDKN1A','RNF224','TTLL6','LIF'),
#                                 ordered = TRUE)
df_dumb_actd$gene_name <- factor(df_dumb_actd$gene_name,
                                 levels=c('BCAS3','KSR1','PRX', 'NCKAP1','IL17RD','SDSL','GALR2','MICALL2','CKB','CCDC120','IRF7',
                                          'IER3','NAV1','CUEDC1','HRK','RBM47','KIFC3','RNF151','CDKN1A','DHDH','RNF224','TNFRSF9','LRRC39',
                                          'TTLL6','LIF',
                                          'CCDC62', 'SLC25A32','SPOCK1', 'DNAJC5G' 
                                 ),
                                 ordered = TRUE)


ggplot(df_dumb_actd,aes(y=gene_name,x=fc_empty,xend=fc_bcl2))+
  geom_dumbbell(size=1,color="gray24",size_x=4,size_xend=4,colour_x="seagreen",colour_xend="purple",alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()+
  labs(x="delta fold change ")+
  geom_hline(yintercept=c(11.5,25.5,27.5), linetype="dashed", 
             color = "#84A6D6", size=1)+
  geom_hline(yintercept=c(25.5),
             color = "#4382BB", size=1.3)


geom_hline(yintercept=c(6.5), linetype="3313", color = "#84A6D6", size=1)+
  geom_hline(yintercept=c(12.5), color = "#4382BB", size=1.3)

# geom_text(data=filter(df_dumb_actd, gene_name=="TTLL6"),
#           aes(x=fc_bcl2, y=gene_name, label="fc bcl2"),
#           color="#b8627dff", size=3, vjust=-1.5, fontface="bold") + coord_flip()
#    geom_text(data=filter(df_dumb_actd, gene_name=="TTLL6"),
#              aes(x=fc_empty, y=gene_name, label="fc empty"),
#              color="#eb8055ff", size=3, vjust=-1.5, fontface="bold")   

#geom_dumbbell(data=infected, aes(y=concerned, x=fc_empty, xend=fc_bcl2), size=1.5, color="#b2b2b2", size_x=3, size_xend = 3, colour_x = red, colour_xend = blue)

# DOXO
df_dumb_doxo<-delta_fc_doxo[match(target_doxo$gene_name,delta_fc_doxo$gene_name),-c(5,6)]
df_dumb_doxo<-df_dumb_doxo[order(abs(df_dumb_doxo$delta)),]  
#df_dumb_doxo<-df_dumb_doxo[order((df_dumb_doxo$delta)),] 

# set the order 
# df_dumb_doxo$gene_name <- factor(df_dumb_doxo$gene_name,
#                                  levels=c('HRK','YWHAH','VPS72','GTPBP4','SNRPB','HNRNPF','SLC3A2','CENPM','ZNF408','RFFL',
#                                           'ABHD14A','ZSCAN21','DDX58','OSER1','PSMC3IP','CTNNBIP1','AEN','LZTS3','DNAJC5G',
#                                           'CD80','GALR2',
#                                     'SLC23A1','MEF2B','BNIPL','CPT1B','SYVN1','RFXAP','SPOCK1','SLC35C1','CHRNA10','NDUFV2',
#                                           'RAI1','SIPA1L3','SDSL','NAV1', 'IRF7','CCDC120','CKB','GADD45B','CUEDC1','MICALL2','FOS',
#                                           'KIFC3','DHDH','RBM47','LRRC39','RNF151','SLC35G2','RPL3L','RNF224','CDKN1A','TTLL6','LIF'
#                                           ),
#                                  ordered = TRUE)


df_dumb_doxo$gene_name<-factor(df_dumb_doxo$gene_name,
                               levels=c('RAI1','SDSL','YWHAH','GTPBP4','HNRNPF','VPS72','CCDC120','IRF7','FOS','SNRPB','MICALL2',
                                        'DHDH','HRK','KIFC3','CUEDC1','RBM47','SIPA1L3','CDKN1A','CKB','NAV1','RNF224','GADD45B',
                                        'LIF','SLC35G2','RPL3L','RNF151', 'TTLL6','LRRC39',   
                                        'CTNNBIP1','CENPM','RFXAP','SLC23A1','OSER1','DDX58','CPT1B','SLC3A2','ZSCAN21','SYVN1','RFFL',
                                        'PSMC3IP','AEN','NDUFV2','ABHD14A','DNAJC5G','MEF2B','GALR2','ZNF408','SLC35C1','CHRNA10','LZTS3',
                                        'CD80','BNIPL','SPOCK1'
                               ),
                               ordered = TRUE)

#df_dumb_doxo$gene_name<-factor(df_dumb_doxo$gene_name,
#                               levels=(df_dumb_doxo$gene_name))



ggplot(df_dumb_doxo,aes(y=gene_name,x=fc_empty,xend=fc_bcl2))+
  geom_dumbbell(size=1,color="gray24",size_x=4,size_xend=4,colour_x="seagreen",colour_xend="purple",alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()+
  labs(x="delta fold change ")+
  geom_hline(yintercept=c(16.5,49.5), linetype="3313", 
             color = "#84A6D6", size=1)+
  geom_hline(yintercept=c(28.5),
             color = "#4382BB", size=1.3)



# CX
df_dumb_cx<-delta_fc_cx[match(target_cx$gene_name,delta_fc_cx$gene_name),-c(5,6)]
df_dumb_cx<-df_dumb_cx[order(abs(df_dumb_cx$delta)),]
# df_dumb_cx<-df_dumb_cx[order((df_dumb_cx$delta)),]


# set the order 
# df_dumb_cx$gene_name<-factor(df_dumb_cx$gene_name,
#                              levels=(df_dumb_cx$gene_name))
# df_dumb_cx$gene_name <- factor(df_dumb_cx$gene_name,
#                                  levels=c('DDX58',
#                                           'CMTM8','SPOCK1','NAV1','MICALL2','CUEDC1','IER3','NCKAP1','RBM47','MYRF','GALR2','RNF224',
#                                           'CDKN1A','LIF'),
#                                  ordered = TRUE)
df_dumb_cx$gene_name <- factor(df_dumb_cx$gene_name,
                               levels=c('GALR2','RBM47','SPOCK1','CUEDC1','NAV1','NCKAP1','IER3','MYRF',
                                        'MICALL2','CDKN1A','RNF224','LIF',
                                        'CMTM8','DDX58'),
                               ordered = TRUE)


ggplot(df_dumb_cx,aes(y=gene_name,x=fc_empty,xend=fc_bcl2))+
  geom_dumbbell(size=1,color="gray24",size_x=4,size_xend=4,colour_x="seagreen",colour_xend="purple",alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip() +
  labs(x="delta fold change ")  +
  geom_hline(yintercept=c(6.5), linetype="3313", color = "#84A6D6", size=1)+
  geom_hline(yintercept=c(12.5), color = "#4382BB", size=1.3)



data_raw<-data.frame(sampleID=c(paste0(colnames(dds),"_",sample_condition)),
                     gene=c(rep("LIF",24)),
                     value=c(counts(dds.norm)[which(rownames(dds.norm)=="ENSG00000128342.5"),]))



data_raw<-data_raw[c(1,9,17,3,11,19,5,13,21,7,15,23),]

# fix order
#data_raw$sampleID <- factor(data_raw$sampleID,levels = c(data_raw$sampleID))

data_raw$sampleID<-paste0(seq(1,12,1),"_",data_raw$sampleID)
data_raw$sampleID <- factor(data_raw$sampleID,levels = c(data_raw$sampleID))

ggplot(data_raw, aes(fill=gene, y=value, x=sampleID)) + 
  geom_bar(position="dodge", stat="identity")+coord_flip()+
  ggtitle("LIF (raw)")



# Read TP53 target database: TRANSFAC, Fischer, iRegulon

TP53_targets<- as.data.frame(read_excel("~/derenzini/datasets/TP53_targets_transfac.xlsx",  col_names = FALSE))
colnames(TP53_targets)<-"gene_name"

# Read Fischer list
TP53_targets_f<- as.data.frame(read_excel("~/derenzini/datasets/TP53_116_targets_fischer.xlsx",  col_names = TRUE))
colnames(TP53_targets_f)<-"gene_name"

# Absolute Targets
DE_targets<-unlist(lapply(DE_list_filt,function(x) {length(intersect(x,TP53_targets$gene_name))}))

DE_targets_fischer<-unlist(lapply(DE_list_filt,function(x) {length(intersect(x,TP53_targets_f$TP53_targets))}))

# N Targets: barplot
# TRANSFAC
df_ntargets<-data.frame(comp=names(DE_list_filt),
                        value=c(12,119,109,18,127,110,12,56,7,13))

# Fisher
df_ntargets<-data.frame(comp=names(DE_list_filt),
                        value=c(15,85,79,49,91,76,28,34,3,7))


# fix order
#df_ntargets$comp <- factor(df_ntargets$comp,levels = c(df_ntargets$comp))

ggplot(df_ntargets, aes(y=value, x=comp)) + 
  geom_bar(position="stack", stat="identity",fill="#7B92AA") +
  coord_flip()+
  #ggtitle("# TP53 targets") +
  xlab("comparisons")+
  ylab("n targets")+
  theme(text = element_text(size=20))+
  geom_text(aes(y = value-4, label = value), vjust = 1, colour = "white",size=6)


# Relative Targets: targets in the merging list

tp53_transfac_spec_empty<-(intersect(spec_empty,TP53_targets$gene_name))
tp53_transfac_spec_bcl2<-(intersect(spec_bcl2,TP53_targets$gene_name))
tp53_transfac_common<-(intersect(common_overall,TP53_targets$gene_name))



length(intersect(spec_empty,TP53_targets_f$gene_name))
length(intersect(spec_bcl2,TP53_targets_f$gene_name))
length(intersect(common_overall,TP53_targets_f$gene_name))
