#' ---
#' title: "SpatialDE GSEA"
#' author: "Jiaqiang Zhu"
#' date: "March 26th, 2020"
#' ---



#+ label=hdst, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
library(ggpubr)
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "CN24_D1"

iblock = 10

load(paste0(workdir,"/data/hdst/",isample,"_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
speres 		<- read.table(paste0(workdir,"/output/spatialDE/clean_",isample,"_idx_",iblock,"X_spe_nomt.csv"),head=T)
# speres[which(speres$g=="log_total_count"),"qval"]

speres 		<- speres[-which(speres$g=="log_total_count"),]

# it is weird, that the log total counts are also significant
genename 	<- sapply(strsplit(as.character(speres$g),split="'"),"[[",2)
pvals 		<- speres$pval
LLR_stat 	<- speres$LLR
names(pvals) <- names(LLR_stat)	<- genename

deg <- cbind.data.frame(symbol=genename,stat=LLR_stat,pval=pvals)

DEG <- deg
df 	<- bitr(names(LLR_stat), fromType = "SYMBOL",
	   toType = c( "ENTREZID"),
	   OrgDb = org.Mm.eg.db)

DEG <- merge(DEG,df,by.y='SYMBOL',by.x='symbol')

geneList  = DEG$stat
names(geneList) = DEG$ENTREZID
geneList = sort(geneList,decreasing = T)

#' HDST 10X resolution, where spatialDE identify 267 signal genes (including log_total_count)
#+ label=hdst_gsea, echo=T, warnings=F, message=F,eval=T
kk_gse <- gseKEGG(geneList     = geneList,
				organism     = 'mmu',
				nPerm        = 1000,
				minGSSize    = 50,
				pvalueCutoff = 0.9,
				verbose      = FALSE)
head(kk_gse,50)[,2:6]
nrow(kk_gse)





#+ label=slideseq, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
library(ggpubr)
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"

iblock = 200

speres 		<- read.table(paste0(workdir,"/output/spatialDE/clean_",isample,"_idx_",iblock,"X_spe_nomt.csv"),head=T)
# speres[which(speres$g=="log_total_count"),"qval"]

speres 		<- speres[-which(speres$g=="log_total_count"),]
if(sum(duplicated(speres$g))>0){
	# the LLR for duplicated results are the same
	speres 		<- speres[-which(duplicated(speres$g)),]
}

# it is weird, that the log total counts are also significant
genename 	<- sapply(strsplit(as.character(speres$g),split="'"),"[[",2)
pvals 		<- speres$pval
LLR_stat 	<- speres$LLR
names(pvals) <- names(LLR_stat)	<- genename

deg <- cbind.data.frame(symbol=genename,stat=LLR_stat,pval=pvals)

DEG <- deg
df 	<- bitr(names(LLR_stat), fromType = "SYMBOL",
	   toType = c( "ENTREZID"),
	   OrgDb = org.Mm.eg.db)

DEG <- merge(DEG,df,by.y='SYMBOL',by.x='symbol')

geneList  = DEG$stat
names(geneList) = DEG$ENTREZID
geneList = sort(geneList,decreasing = T)


#' Slideseq 50X resolution, where spatialDE identify 6915 signal genes (including log_total_count)
#+ label=slideseq_gsea, echo=T, warnings=F, message=F,eval=T
kk_gse <- gseKEGG(geneList     = geneList,
				organism     = 'mmu',
				nPerm        = 1000,
				minGSSize    = 50,
				pvalueCutoff = 0.9,
				verbose      = FALSE)
head(kk_gse,50)[,2:6]
nrow(kk_gse)

