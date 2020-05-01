#' ---
#' title: "Mouse Olfactory Bulb Analysis"
#' author: "Jiaqiang Zhu"
#' date: "April 30th, 2020"
#' ---


#+ label=analysis_code, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)

source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")

datapath = "/net/mulan/disk2/jiaqiang/SDE/Realdata/rawdata/"
countdata   <- read.table(paste0(datapath,"/Rep11_MOB_count_matrix-1.tsv"),check.names=F)
rn          <- rownames(countdata)
info        <- cbind.data.frame(x=as.numeric(sapply(strsplit(rn,split="x"),"[",1)),y=as.numeric(sapply(strsplit(rn,split="x"),"[",2)))
rownames(info) <- rn

sp_count    <- as(t(countdata),"sparseMatrix")

removed_cell <- which(as.vector(sp_sums_Rcpp(sp_count))<10)

if(length(removed_cell)>0){
	sp_count <- sp_count[,-removed_cell]
	info <- info[-removed_cell,]
}
numGene = nrow(sp_count)
numSample = nrow(info)

raw_loc  	<- as.matrix(info[,1:2])
t1 			<- proc.time()
KDC 		<- kdc_mk_spark(sp_count,raw_loc,numCores=20,filter=F,option="mixture",verbose=FALSE)
t2 			<- proc.time() - t1
runTime_kdc = t2

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v2/"
save(KDC,runTime_kdc,file=paste0(workdir,"/output/kdc/Rep11_MOB_ten_removed_kdc.rds"))




#+ label=comparison, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v2/"
load(paste0(workdir,"/output/kdc/Rep11_MOB_ten_removed_kdc.rds"))

kdc_gene <- rownames(KDC$res_mtest)[which(KDC$res_mtest$adjustedPval<0.05)]


top_KDC <- head(rownames(KDC$res_mtest[order(KDC$res_mtest$combinedPval),]),10)


data_path <- "/net/mulan/disk2/jiaqiang/SDE/revision/v1/realdata_summary"
spe <- read.csv(paste0(data_path,"/Rep11_MOB_all_genes_spatialDE_pval.csv"),row.names=1)

spe_gene 	<- rownames(spe)[which(spe$qval<0.05)]
spark 		<- read.csv(paste0(data_path,"/Rep11_MOB_all_genes_spark_pval.csv"),row.names=1)
spark_gene 	<- rownames(spark)[which(spark$adjusted_pvalue<0.05)]

sq_color 	<- c("hotpink","mediumorchid2","lightskyblue")
signal_list <- list(KDC=kdc_gene,SpatialDE=spe_gene,SPARK=spark_gene)
allgene 	<- unique(unlist(signal_list))


pd <- cbind.data.frame(KDC=allgene%in%signal_list[[1]],
						SpatialDE=allgene%in%signal_list[[2]],
						SPARK=allgene%in%signal_list[[3]])
rownames(pd) <- allgene
library(eulerr)
library(ggsci)

VennDiag 	<- euler(pd,shape="ellipse") ## other wise some number if missing 

venn_list 	<- plot(VennDiag, quantities = T, lty=rep(0,3),fill=sq_color,
    alpha=1,legend=T,labels=F)

#' ***
#' > Figure 1: Venn Diagram
#+ label=venn1,fig.width=8, fig.height=8, echo=F,fig.align="center",eval=T
venn_list


#' Top Genes By KDC only
#+ label=topgene, echo=F, warnings=F, message=F,eval=T
setdiff(top_KDC,union(spe_gene,spark_gene))


#+ label=gene_pattern, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Matrix)
library(Rcpp)

source("/net/fantasia/home/jiaqiang/mulan_temp/shortcut_function/simple.qqplot.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/freq_func.R")

countdata   <- read.table("/net/fantasia/home/jiaqiang/mulan_temp/RealData/rawdata/Rep11_MOB_count_matrix-1.tsv",check.names=F)
rn 			<- rownames(countdata)
info 		<- cbind.data.frame(x=as.numeric(sapply(strsplit(rn,split="x"),"[",1)),y=as.numeric(sapply(strsplit(rn,split="x"),"[",2)))
rownames(info) <- rn
countdata <- t(countdata)

mainGene 		<- c("Celf4","Mll2","Mrfap1","Rpl6","Rplp1")
sig_ct      	<- countdata[mainGene,]
vst_ct         	<- var_stabilize(countdata)
sig_vst_ct     	<- vst_ct[mainGene,]
rel_vst_ct  	<- apply(sig_vst_ct,1,relative_func)
pltdat      	<- cbind.data.frame(info[,1:2],rel_vst_ct)
pp      		<- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot2(pltdat,x,main=T,titlesize=1.9)})

library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],
			ncol = 3, nrow = 1)

pltdat2      	<- cbind.data.frame(info[,1:2],t(sig_ct))
pp2      		<- lapply(1:(ncol(pltdat2)-2),function(x){pattern_plot2(pltdat2,x,main=T,titlesize=1.9)})

library(ggpubr)
pp_list2 <- ggarrange(pp2[[1]], pp2[[2]], pp2[[3]],
			ncol = 3, nrow = 1)

#' ***
#' > Figure 2: Gene Expression Pattern (Raw Count)
#+ label=gep1,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
pp_list2



#' ***
#' > Figure 3: Gene Expression Pattern (VST Residual)
#+ label=gep2,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
pp_list





