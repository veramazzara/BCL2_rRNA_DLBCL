# This script aims to prepare BCL2 and MYC KM in the discovery and validation cohort
# (ex FigS2)

# Load libraries
library(openxlsx)
library(readxl)
library(survival)
library(survminer)
library(dplyr)

################# Trial Data: load data
new_data_dlcl04 <- as.data.frame(read_excel("Casi Bologna per Vera COMBO.xlsx",sheet = "DLCL04"))

# extract only useful data
data_dlcl04<-new_data_dlcl04[,c(1,63,38,41,52,78,79,71,72,74,75)]

# change colnames
colnames(data_dlcl04)<-c("CodicePaz","AgNOR", "COO","BCL2","MYC","BCL2_log2","MYC_log2","Mesi_PFS","Evento_PFS","Mesi_OS","Evento_OS")

# COO code: harmonize the code in COO
data_dlcl04$COO<-gsub("UNCLASSIFIED","Unclassified",data_dlcl04$COO)

# AgNOR
data_dlcl04$AgNOR<-as.numeric(data_dlcl04$AgNOR)

# IPI: maunal insertion of IPI
data_dlcl04$IPI<-c("High",rep("Intermediate-High",11),"High",rep("Intermediate-High",7),
                   "High",rep("Intermediate-High",3),"High",rep("Intermediate-High",5),
                   rep("High",2),"Intermediate-High","High",rep("Intermediate-High",3),
                   "High",rep("Intermediate-High",6),"High","Intermediate-High")

#### Median Threshold
data_dlcl04$BCL2_median_class<- data_dlcl04$BCL2_log2 > median(data_dlcl04$BCL2_log2)
data_dlcl04$MYC_median_class<- data_dlcl04$MYC_log2 > median(data_dlcl04$MYC_log2)

# change TRUE with High and FALSE with Low
data_dlcl04$MYC_median_class<-gsub("TRUE","High",data_dlcl04$MYC_median_class)
data_dlcl04$MYC_median_class<-gsub("FALSE","Low",data_dlcl04$MYC_median_class)

data_dlcl04$BCL2_median_class<-gsub("TRUE","High",data_dlcl04$BCL2_median_class)
data_dlcl04$BCL2_median_class<-gsub("FALSE","Low",data_dlcl04$BCL2_median_class)

# KM
fit_sig<-survfit(Surv(Mesi_PFS,Evento_PFS) ~ BCL2_median_class, data = data_dlcl04)

ggsurvplot(fit_sig,pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.title="",
           xlab="Months",
           ylab="Progression Free Survival",
           risk.table.col = "strata", # Change risk table color by groups
           #legend.labs = c( "High","Low"),
           legend.title="Groups",
          # palette=c("#c30101","#0a00ff"),
           ggtheme = theme_bw(),
           title="PFS BCL2 (DLCL04)")




################# Bologna Data: load data
# A. Complete dataset
# Bologna Data: load data
new_data_bo <- as.data.frame(read_excel("Casi Bologna per Vera.xlsx"))

# modify the code related to Evento_OS, Evento_PFS
new_data_bo$`ricaduta/pd_corretto`[which(new_data_bo$`ricaduta/pd_corretto`=="non trovato ID")]<-new_data_bo$`ricaduta/pd`[which(new_data_bo$`ricaduta/pd_corretto`=="non trovato ID")]

# add the correct code for OS
new_data_bo$decesso_corretto<-0
new_data_bo$decesso_corretto[which(new_data_bo$decesso==0)]<-1

# extract only useful data
data_bo<-new_data_bo[,c(1,104,10,41,52,28,25,24,31,105,34,17)]
              
# change colnames
colnames(data_bo)<-c("CodicePaz","AgNOR", "COO","BCL2","MYC", "Mesi_PFS","Evento_PFS","Evento_PFS_original","Mesi_OS","Evento_OS","Evento_OS_original","IPI")

# BCL2 and MYC: add log2 value
data_bo$BCL2_log2<-log2(data_bo$BCL2)
data_bo$MYC_log2<-log2(data_bo$MYC)

# PFS status: filter pts with missing value 
data_bo<-filter(data_bo,Evento_PFS!="missing")

# AgnOR: filter pts with NA 
data_bo<-filter(data_bo,AgNOR!="NA")

# PFS/OS: change the class type
data_bo$Evento_PFS<-as.numeric(data_bo$Evento_PFS)
data_bo$Evento_OS<-as.numeric(data_bo$Evento_OS)

# IPI: conversion from numeric to letteral code
data_bo$IPI<-gsub("0","Low-Low Intermediate",data_bo$IPI)
data_bo$IPI<-gsub("1","Low-Low Intermediate",data_bo$IPI)
data_bo$IPI<-gsub("2","Intermediate-High",data_bo$IPI)
data_bo$IPI<-gsub("3","High",data_bo$IPI)
data_bo$IPI<-gsub("4","High",data_bo$IPI)


###### Median Threshold
# define the High and Low Class
data_bo$BCL2_median_class<- data_bo$BCL2_log2 > median(data_bo$BCL2_log2)
data_bo$MYC_median_class<- data_bo$MYC_log2 > median(data_bo$MYC_log2)
data_bo$AgNOR_median_class<- data_bo$AgNOR > median(data_bo$AgNOR)

# change TRUE with High and FALSE with Low
data_bo$MYC_median_class<-gsub("TRUE","High",data_bo$MYC_median_class)
data_bo$MYC_median_class<-gsub("FALSE","Low",data_bo$MYC_median_class)

data_bo$BCL2_median_class<-gsub("TRUE","High",data_bo$BCL2_median_class)
data_bo$BCL2_median_class<-gsub("FALSE","Low",data_bo$BCL2_median_class)

data_bo$AgNOR_median_class<-gsub("TRUE","High",data_bo$AgNOR_median_class)
data_bo$AgNOR_median_class<-gsub("FALSE","Low",data_bo$AgNOR_median_class)

# KM
fit_sig<-survfit(Surv(Mesi_PFS,Evento_PFS_original) ~ MYC_median_class, data = data_bo)

ggsurvplot(fit_sig,pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.title="",
           title="PFS")


ggsurvplot(fit_sig,pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.title="",
           xlab="Months",
           ylab="Progression Free Survival",
           risk.table.col = "strata", # Change risk table color by groups
           legend.labs = c( "High","Low"),
           legend.title="Groups",
            palette=c("#c30101","#0a00ff"),
           ggtheme = theme_bw(),
           title="PFS MYC (BO)")



# B. Filtered dataset

new_data_bo_filt <- as.data.frame(read_excel("Casi Bologna per Vera COMBO.xlsx",sheet = "COMBO 101 (CASI BO FILTRATI)"))
colnames(new_data_bo_filt)[2]<-"CodicePaz"

temp<-semi_join(data_bo[,c(1,6)],new_data_bo_filt[,c(2,8)])

data_bo_filt<-semi_join(data_bo,temp)

######## Median threshold

# define the High and Low Class
data_bo_filt$BCL2_median_class<- data_bo_filt$BCL2_log2 > median(data_bo_filt$BCL2_log2)
data_bo_filt$MYC_median_class<- data_bo_filt$MYC_log2 > median(data_bo_filt$MYC_log2)
data_bo_filt$AgNOR_median_class<- data_bo_filt$AgNOR > median(data_bo_filt$AgNOR)

# change TRUE with High and FALSE with Low
data_bo_filt$MYC_median_class<-gsub("TRUE","High",data_bo_filt$MYC_median_class)
data_bo_filt$MYC_median_class<-gsub("FALSE","Low",data_bo_filt$MYC_median_class)

data_bo_filt$BCL2_median_class<-gsub("TRUE","High",data_bo_filt$BCL2_median_class)
data_bo_filt$BCL2_median_class<-gsub("FALSE","Low",data_bo_filt$BCL2_median_class)

data_bo_filt$AgNOR_median_class<-gsub("TRUE","High",data_bo_filt$AgNOR_median_class)
data_bo_filt$AgNOR_median_class<-gsub("FALSE","Low",data_bo_filt$AgNOR_median_class)

# DEXP
data_bo_filt$DEXP_median_class<-paste(data_bo_filt$MYC_median_class,data_bo_filt$BCL2_median_class,sep="_")

data_bo_filt$DEXP_median_class<-gsub("High_High","yes",data_bo_filt$DEXP_median_class)
data_bo_filt$DEXP_median_class<-gsub("High_Low","no",data_bo_filt$DEXP_median_class)
data_bo_filt$DEXP_median_class<-gsub("Low_High","no",data_bo_filt$DEXP_median_class)
data_bo_filt$DEXP_median_class<-gsub("Low_Low","no",data_bo_filt$DEXP_median_class)


# KM
fit_sig<-survfit(Surv(Mesi_PFS,Evento_OS_original) ~ AgNOR_median_class, data = data_bo_filt)

ggsurvplot(fit_sig,pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.title="",
           title="OS AgNOR (BO filtered)")




################# Combo Data: load data
##### Maxstat Threshold
# A1. COMBO 135
# Note: COMBO 135 sheet has wrong code for PFS and OS 

#new_data_combo <- as.data.frame(read_excel("Casi Bologna per Vera COMBO.xlsx",sheet = "COMBO 135"))
# Patient selection: extratct only Sample_id and match with the separate datasets
#patient_id<-new_data_combo$`codice paziente...2`

# Data Merging: merging bologna and dlc04 datasets
data_combo_135<-(data_bo[,c(1,2,3,4,5,13,14,6,8,9,11,12)])

colnames(data_combo_135)[c(9,11)]<-c("Evento_PFS","Evento_OS")
                      
data_combo_135<-rbind(data_combo_135,data_dlcl04[,c(1,2,3,4,5,6,7,8,9,10,11,12)])

 # Maxstat Threshold: select OS or PFS in the code
 # BCL2
  res.cut_bcl2 <- surv_cutpoint(data_combo_135, time = "Mesi_PFS", event = "Evento_PFS",
                              variables = c("BCL2_log2"))
  res.cat_bcl2 <- surv_categorize(res.cut_bcl2)

  # MYC
  res.cut_myc <- surv_cutpoint(data_combo_135, time = "Mesi_PFS", event = "Evento_PFS",
                             variables = c("MYC_log2"))
  res.cat_myc <- surv_categorize(res.cut_myc)

  # AgNOR
  res.cut_agnor <- surv_cutpoint(data_combo_135, time = "Mesi_PFS", event = "Evento_PFS",
                               variables = c("AgNOR"))
  res.cat_agnor <- surv_categorize(res.cut_agnor)

# Assign the class
data_combo_135$BCL2_class<-res.cat_bcl2$BCL2_log2
data_combo_135$MYC_class<-res.cat_myc$MYC_log2
data_combo_135$AgNOR_class<-res.cat_agnor$AgNOR

# DEXP
data_combo_135$DEXP_class<-paste(data_combo_135$MYC_class,data_combo_135$BCL2_class,sep="_")

data_combo_135$DEXP_class<-gsub("high_high","yes",data_combo_135$DEXP_class)
data_combo_135$DEXP_class<-gsub("high_low","no",data_combo_135$DEXP_class)
data_combo_135$DEXP_class<-gsub("low_high","no",data_combo_135$DEXP_class)
data_combo_135$DEXP_class<-gsub("low_low","no",data_combo_135$DEXP_class)

# KM
fit_sig<-survfit(Surv(Mesi_PFS,Evento_PFS) ~ AgNOR_class, data = data_combo_135)

ggsurvplot(fit_sig,pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.title="",
           title="PFS AgNOR (COMBO 135)")

# Multivariate analysis
#multi_cox_os<-coxph((Surv(Mesi_OS, Evento_OS) ~  BCL2_class + MYC_class + AgNOR_class + COO ), data=data_combo_135)
#ggforest(multi_cox_os,data=data_combo_135,fontsize = 1.2)

multi_cox_pfs<-coxph((Surv(Mesi_PFS, Evento_PFS) ~ BCL2_class + MYC_class + DEXP_class  + AgNOR_class + COO  ), data=data_combo_135)
ggforest(multi_cox_pfs,data=data_combo_135,fontsize = 1.2)

# A2. COMBO 101

# Data Merging: merging bologna and dlc04 datasets
data_combo_101<-(data_bo_filt[,c(1,2,3,4,5,13,14,6,8,9,11,12)])

colnames(data_combo_101)[c(9,11)]<-c("Evento_PFS","Evento_OS")

data_combo_101<-rbind(data_combo_101,data_dlcl04[,c(1,2,3,4,5,6,7,8,9,10,11,12)])
# Maxstat Threshold: select OS or PFS in the code
# BCL2
res.cut_bcl2 <- surv_cutpoint(data_combo_101, time = "Mesi_PFS", event = "Evento_PFS",
                              variables = c("BCL2_log2"))
res.cat_bcl2 <- surv_categorize(res.cut_bcl2)

# MYC
res.cut_myc <- surv_cutpoint(data_combo_101, time = "Mesi_PFS", event = "Evento_PFS",
                             variables = c("MYC_log2"))
res.cat_myc <- surv_categorize(res.cut_myc)

# AgNOR
res.cut_agnor <- surv_cutpoint(data_combo_101, time = "Mesi_PFS", event = "Evento_PFS",
                               variables = c("AgNOR"))
res.cat_agnor <- surv_categorize(res.cut_agnor)

# Assign the class
data_combo_101$BCL2_class<-res.cat_bcl2$BCL2_log2
data_combo_101$MYC_class<-res.cat_myc$MYC_log2
data_combo_101$AgNOR_class<-res.cat_agnor$AgNOR

# DEXP
data_combo_101$DEXP_class<-paste(data_combo_101$MYC_class,data_combo_101$BCL2_class,sep="_")

data_combo_101$DEXP_class<-gsub("high_high","yes",data_combo_101$DEXP_class)
data_combo_101$DEXP_class<-gsub("high_low","no",data_combo_101$DEXP_class)
data_combo_101$DEXP_class<-gsub("low_high","no",data_combo_101$DEXP_class)
data_combo_101$DEXP_class<-gsub("low_low","no",data_combo_101$DEXP_class)

# Set the reference
data_combo_101$DEXP_class <- relevel(factor(data_combo_101$DEXP_class), ref="no")
data_combo_101$BCL2_class <- relevel(factor(data_combo_101$BCL2_class), ref="low")

# KM
fit_sig<-survfit(Surv(Mesi_PFS,Evento_PFS) ~ BCL2_class, data = data_combo_101)

ggsurvplot(fit_sig,pval = TRUE,
           risk.table = TRUE, # Add risk tablenextflow youtube
           
           risk.table.title="",
           title="PFS BCL2 (COMBO 101)")

# Multivariate analysis
#multi_cox_os<-coxph((Surv(Mesi_OS, Evento_OS) ~ BCL2_class + MYC_class + DEXP_class+ AgNOR_class + COO ), data=data_combo_101)
#ggforest(multi_cox_os,data=data_combo_101,fontsize = 1.2)

multi_cox_pfs<-coxph((Surv(Mesi_PFS, Evento_PFS) ~ BCL2_class + MYC_class + DEXP_class  + AgNOR_class + COO  ), data=data_combo_101)
ggforest(multi_cox_pfs,data=data_combo_101,fontsize = 1.2)




###### Median Threshold: select OS or PFS in the code

# B1: COMBO 135

# define the High and Low Class

data_combo_135$BCL2_median_class<- data_combo_135$BCL2_log2 > median(data_combo_135$BCL2_log2)
data_combo_135$MYC_median_class<- data_combo_135$MYC_log2 > median(data_combo_135$MYC_log2)
data_combo_135$AgNOR_median_class<- data_combo_135$AgNOR > median(data_combo_135$AgNOR)

# change TRUE with High and FALSE with Low
data_combo_135$MYC_median_class<-gsub("TRUE","High",data_combo_135$MYC_median_class)
data_combo_135$MYC_median_class<-gsub("FALSE","Low",data_combo_135$MYC_median_class)

data_combo_135$BCL2_median_class<-gsub("TRUE","High",data_combo_135$BCL2_median_class)
data_combo_135$BCL2_median_class<-gsub("FALSE","Low",data_combo_135$BCL2_median_class)

data_combo_135$AgNOR_median_class<-gsub("TRUE","High",data_combo_135$AgNOR_median_class)
data_combo_135$AgNOR_median_class<-gsub("FALSE","Low",data_combo_135$AgNOR_median_class)

 # DEXP
data_combo_135$DEXP_median_class<-paste(data_combo_135$MYC_median_class,data_combo_135$BCL2_median_class,sep="_")

data_combo_135$DEXP_median_class<-gsub("High_High","yes",data_combo_135$DEXP_median_class)
data_combo_135$DEXP_median_class<-gsub("High_Low","no",data_combo_135$DEXP_median_class)
data_combo_135$DEXP_median_class<-gsub("Low_High","no",data_combo_135$DEXP_median_class)
data_combo_135$DEXP_median_class<-gsub("Low_Low","no",data_combo_135$DEXP_median_class)

# Set the reference
data_combo_135$BCL2_median_class <- relevel(factor(data_combo_135$BCL2_median_class), ref="Low")
data_combo_135$MYC_median_class <- relevel(factor(data_combo_135$MYC_median_class), ref="Low")
data_combo_135$AgNOR_median_class <- relevel(factor(data_combo_135$AgNOR_median_class), ref="Low")

data_combo_135$DEXP_median_class <- relevel(factor(data_combo_135$DEXP_median_class), ref="no")
data_combo_135$IPI <- relevel(factor(data_combo_135$IPI), ref="Low-Low Intermediate")
data_combo_135$COO <- relevel(factor(data_combo_135$COO), ref="GCB")

# KM
fit_sig<-survfit(Surv(Mesi_OS,Evento_OS) ~ AgNOR_median_class, data = data_combo_135)
fit_sig<-survfit(Surv(Mesi_PFS,Evento_PFS) ~ AgNOR_median_class, data = data_combo_135)



fit_sig<-survfit(Surv(Mesi_PFS,Evento_PFS) ~ BCL2_median_class, data = data_combo_135)
fit_sig<-survfit(Surv(Mesi_PFS,Evento_PFS) ~ MYC_median_class, data = data_combo_135)


ggsurvplot(fit_sig,pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.title="",
           title="PFS AgNOR (COMBO 135)")

ggsurvplot(fit_sig,
           pval = TRUE, conf.int = T,
           risk.table = TRUE, # Add risk table
           risk.table.title="",
           risk.table.height = 0.15,
           tables.theme = theme_cleantable(),
           xlab="Time (months)",
           ylab="OS",
           title="",
           legend.title="Groups",
           #palette=c("#c30101","#0a00ff"),
           palette=c("#0a00ff","#c30101"),   # OS
           #legend.labs = c("high", "low"),
           legend.labs = c("Low", "High"),   #OS
           risk.table.fontsize = 12,
           pval.size = 12,
           font.x = c(30, "bold.italic"),
           font.y = c(30, "bold.italic"),
           font.tickslab = c(22, "plain"),
           risk.table.y.text.col = T, # colour risk table text annotations.
           risk.table.y.text = FALSE,
           ggtheme = theme_bw()
           
)



# Multivariate analysis
multi_cox_os<-coxph((Surv(Mesi_OS, Evento_OS) ~ BCL2_median_class  + DEXP_median_class + AgNOR_median_class + COO ), data=data_combo_135)
ggforest(multi_cox_os,data=data_combo_135,fontsize = 1.2)

multi_cox_pfs<-coxph((Surv(Mesi_PFS, Evento_PFS) ~ BCL2_median_class  + DEXP_median_class  + AgNOR_median_class + COO  ), data=data_combo_135)
ggforest(multi_cox_pfs,data=data_combo_135,fontsize = 1.2)


# B2: COMBO 101

# define the High and Low Class
data_combo_101$BCL2_median_class<- data_combo_101$BCL2_log2 > median(data_combo_101$BCL2_log2)
data_combo_101$MYC_median_class<- data_combo_101$MYC_log2 > median(data_combo_101$MYC_log2)
data_combo_101$AgNOR_median_class<- data_combo_101$AgNOR > median(data_combo_101$AgNOR)

# change TRUE with High and FALSE with Low
data_combo_101$MYC_median_class<-gsub("TRUE","High",data_combo_101$MYC_median_class)
data_combo_101$MYC_median_class<-gsub("FALSE","Low",data_combo_101$MYC_median_class)

data_combo_101$BCL2_median_class<-gsub("TRUE","High",data_combo_101$BCL2_median_class)
data_combo_101$BCL2_median_class<-gsub("FALSE","Low",data_combo_101$BCL2_median_class)

data_combo_101$AgNOR_median_class<-gsub("TRUE","High",data_combo_101$AgNOR_median_class)
data_combo_101$AgNOR_median_class<-gsub("FALSE","Low",data_combo_101$AgNOR_median_class)

# DEXP
data_combo_101$DEXP_median_class<-paste(data_combo_101$MYC_median_class,data_combo_101$BCL2_median_class,sep="_")

data_combo_101$DEXP_median_class<-gsub("High_High","yes",data_combo_101$DEXP_median_class)
data_combo_101$DEXP_median_class<-gsub("High_Low","no",data_combo_101$DEXP_median_class)
data_combo_101$DEXP_median_class<-gsub("Low_High","no",data_combo_101$DEXP_median_class)
data_combo_101$DEXP_median_class<-gsub("Low_Low","no",data_combo_101$DEXP_median_class)


# Set the reference
data_combo_101$BCL2_median_class <- relevel(factor(data_combo_101$BCL2_median_class), ref="Low")
data_combo_101$MYC_median_class <- relevel(factor(data_combo_101$MYC_median_class), ref="Low")
data_combo_101$AgNOR_median_class <- relevel(factor(data_combo_101$AgNOR_median_class), ref="Low")

data_combo_101$DEXP_median_class <- relevel(factor(data_combo_101$DEXP_median_class), ref="no")
data_combo_101$IPI <- relevel(factor(data_combo_101$IPI), ref="Low-Low Intermediate")
data_combo_101$COO <- relevel(factor(data_combo_101$COO), ref="GCB")


# KM
fit_sig<-survfit(Surv(Mesi_OS,Evento_OS) ~ AgNOR_median_class, data = data_combo_101)

ggsurvplot(fit_sig,pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.title="",
           title="OS AgNOR (COMBO 101)")

# Multivariate analysis
multi_cox_os<-coxph((Surv(Mesi_OS, Evento_OS) ~ BCL2_median_class + DEXP_median_class + AgNOR_median_class + COO  ), data=data_combo_101)
ggforest(multi_cox_os,data=data_combo_101,fontsize = 1.2)

multi_cox_pfs<-coxph((Surv(Mesi_PFS, Evento_PFS) ~ BCL2_median_class + DEXP_median_class  + AgNOR_median_class + COO), data=data_combo_101)
ggforest(multi_cox_pfs,data=data_combo_101,fontsize = 1.2)



