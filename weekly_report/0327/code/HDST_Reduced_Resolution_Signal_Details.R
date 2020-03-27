#' ---
#' title: "HDST Reduced Resolution Signal Comparison Details"
#' author: "Jiaqiang Zhu"
#' date: "March 22rd, 2020"
#' ---


#' ***
#' > Comparison across methods 
#+ label=threemeth, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(SPARK)

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "CN24_D1"

venn_list <- list()
icount = 0
for(iblock in c(50,25,20,15,10)){
	icount = icount + 1
	# spark
	load(paste0(workdir,"/output/spark/clean_",isample,"_idx_",iblock,"X_spark_default_nomt.rds"))
	spark_conv 	<- spark@res_mtest[sapply(spark@res_vc,function(x){x$converged}),]
	spark_gene 	<- rownames(spark_conv)[spark_conv$adjusted_pvalue<0.05]

	# spatialDE
	spe 		<- read.table(paste0(workdir,"/output/spatialDE/clean_",isample,"_idx_",iblock,"X_spe_nomt.csv"),head=T)
	spe_gene_raw <- as.character(spe$g[spe$qval<0.05])
	if("log_total_count" %in% spe_gene_raw){
		spe_gene_raw <- spe_gene_raw[-which(spe_gene_raw=="log_total_count")]
	}
	spe_gene <- sapply(strsplit(spe_gene_raw,split="'"),"[[",2)

	# KDC
	load(paste0(workdir,"/output/kdc/clean_",isample,"_idx_",iblock,"X_Joint_kdc.rds"))
	mtgene <- which(substr(rownames(KDC$res_mtest),start=1,stop=3)=="mt-")
	if(length(mtgene)>0){
		kdc_nomt <- KDC$res_mtest[-mtgene,]
	} 
	kdc_gene <- rownames(kdc_nomt)[p.adjust(kdc_nomt$combinedPval,method="BY")<0.05]

	# venn
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
	if(icount==5){
		venn_list[[icount]] 	<- plot(VennDiag, quantities = T, lty=rep(0,3),fill=sq_color,
	    alpha=1,legend=T,labels=F)
	}else{
		venn_list[[icount]] 	<- plot(VennDiag, quantities = T, lty=rep(0,3),fill=sq_color,
	    alpha=1,legend=F,labels=F)
	}	
	cat(iblock,"X: ## Signal Genes in KDC and Filtered by SPARK:",length(setdiff(kdc_gene,rownames(spark_conv))),"\n")
	rm(list=setdiff(ls(),c("venn_list","isample","workdir","icount")))
}

# length(intersect(kdc_gene,spark_gene))
# length(intersect(kdc_gene,spe_gene))
# setdiff(spark_gene,kdc_gene)
# 515 genes that significant by KDC is also in SPARK dataset, just not significant
# length(intersect(setdiff(kdc_gene,spark_gene),rownames(spark_conv)))

library(ggpubr)
fig0 <- ggarrange(venn_list[[1]], venn_list[[2]], venn_list[[3]],
					venn_list[[4]],venn_list[[5]],
			labels = c("A", "B", "C","D","E"),
			font.label=list(size=20),
			ncol = 3, nrow = 2)


#' ***
#' > Figure 0: Venn Across Methods at Different Resoultions: 50, 25, 20, 15, 10
#+ label=venn1,fig.width=12, fig.height=8, fig.cap="Venn Diagrams",echo=F,fig.align="center"
fig0




#' ***
#' > Comparison across resolutions: SPARK
#+ label=spark, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(SPARK)

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "CN24_D1"


signal_list <- c()
icount = 0
for(iblock in c(50,25,20,15,10)){
	icount = icount + 1
	load(paste0(workdir,"/output/spark/clean_",isample,"_idx_",iblock,"X_spark_default_nomt.rds"))
	spark_conv 	<- spark@res_mtest[sapply(spark@res_vc,function(x){x$converged}),]
	spark_gene 	<- rownames(spark_conv)[spark_conv$adjusted_pvalue<0.05]
	signal_list[[icount]] <- spark_gene
	rm(spark,spark_gene,spark_conv)
}

names(signal_list) <- paste0("X",c(50,25,20,15,10))

#' Signal Genes at different resolutions
#+ echo=F, warnings=F, message=F,eval=T
sapply(signal_list,length)


#+ echo=F, warnings=F, message=F,eval=F
# example code
# a1 <- kdc_gene
# a2 <- spe_gene
# a11 <- setdiff(a1,union(a2,a3))
# a22 <- setdiff(a2,a1)
# a12 <- intersect(a1,a2)

# pd <- sapply(list(a11,a22,a12),length)
# library(eulerr)
# VennDiag <- euler(c("SPARK" = pd[1],"Moran"= pd[2],
#     "SPARK&Moran" = pd[3]),
# 	shape="ellipse")

# sq_color <- c("hotpink","mediumorchid2","lightskyblue","lightgreen")
# plot(VennDiag, quantities = list(cex = 2), lty=rep(0,2),
#      fill=sq_color[c(1,4)],alpha=0.5,legend=F,labels=F)


#+ echo=F, warnings=F, message=F,eval=T
allgene <- unique(unlist(signal_list))
pd <- cbind.data.frame(X50=allgene%in%signal_list[[1]],
						X25=allgene%in%signal_list[[2]],
						X20=allgene%in%signal_list[[3]],
						X15=allgene%in%signal_list[[4]],
						X10=allgene%in%signal_list[[5]])
rownames(pd) <- allgene
library(eulerr)
VennDiag <- euler(pd,shape="ellipse") ## other wise some number if missing 
library(ggsci)
fig1 <- plot(VennDiag, quantities =T, lty=rep(0,ncol(pd)),fill=pal_npg("nrc")(ncol(pd)),
    alpha=1,legend=F,labels=T)


#' ***
#' > Figure 1: Venn Diag within SPARK
#+ label=fig1,fig.width=8, fig.height=8,echo=F,fig.align="center"
fig1




#' ***
#' > Comparison across resolutions: KDC
#+ label=kdc, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "CN24_D1"

signal_list <- c()
icount = 0
for(iblock in c(50,25,20,15,10,5,1)){
	icount = icount + 1
	load(paste0(workdir,"/output/kdc/clean_",isample,"_idx_",iblock,"X_Joint_kdc.rds"))
	mtgene <- which(substr(rownames(KDC$res_mtest),start=1,stop=3)=="mt-")
	if(length(mtgene)>0){
		kdc_nomt <- KDC$res_mtest[-mtgene,]
	} 
	kdc_gene <- rownames(kdc_nomt)[p.adjust(kdc_nomt$combinedPval,method="BY")<0.05]

	signal_list[[icount]] <- kdc_gene
	rm(kdc_nomt,mtgene,KDC)
}

names(signal_list) <- paste0("X",c(50,25,20,15,10,5,1))

allgene <- unique(unlist(signal_list))
pd <- cbind.data.frame(X50=allgene%in%signal_list[[1]],
						X25=allgene%in%signal_list[[2]],
						X20=allgene%in%signal_list[[3]],
						X15=allgene%in%signal_list[[4]],
						X10=allgene%in%signal_list[[5]],
						X5=allgene%in%signal_list[[6]],
						X1=allgene%in%signal_list[[7]])
rownames(pd) <- allgene
library(eulerr)
VennDiag <- euler(pd,shape="ellipse") ## other wise some number if missing 
library(ggsci)
fig1 <- plot(VennDiag, quantities =T, lty=rep(0,ncol(pd)),fill=pal_npg("nrc")(ncol(pd)),
    alpha=1,legend=F,labels=T)

#' Signal Genes at different resolutions
#+ echo=F, warnings=F, message=F,eval=T
sapply(signal_list,length)

#' ***
#' > Figure 2: Venn Diag within KDC
#+ label=fig2,fig.width=8, fig.height=8,echo=F,fig.align="center"
fig1



#' ***
#' > Comparison across resolutions: KDC 1X vs. SPARK
#+ label=kdcspark, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "CN24_D1"

iblock = 50
load(paste0(workdir,"/data/hdst/clean_",isample,"_idx_",iblock,"X_sparseCount.rds"))

load(paste0(workdir,"/output/kdc/clean_",isample,"_idx_1X_Joint_kdc.rds"))
mtgene <- which(substr(rownames(KDC$res_mtest),start=1,stop=3)=="mt-")
if(length(mtgene)>0){
	kdc_nomt <- KDC$res_mtest[-mtgene,]
} 
kdc_gene <- rownames(kdc_nomt)[p.adjust(kdc_nomt$combinedPval,method="BY")<0.05]

signal_list <- c()
icount = 0
for(iblock in c(50,25,20,15,10)){
	icount = icount + 1
	load(paste0(workdir,"/output/spark/clean_",isample,"_idx_",iblock,"X_spark_default_nomt.rds"))
	spark_conv 	<- spark@res_mtest[sapply(spark@res_vc,function(x){x$converged}),]
	spark_gene 	<- rownames(spark_conv)[spark_conv$adjusted_pvalue<0.05]
	signal_list[[icount]] <- spark_gene
	rm(spark,spark_gene,spark_conv)
}

names(signal_list) <- paste0("X",c(50,25,20,15,10))
kdc_only <- Reduce(intersect,lapply(signal_list,function(x){setdiff(kdc_gene,x)}))

#' Number of Signal Genes in KDC not in SPARK 
#+ echo=F, warnings=F, message=F,eval=T
sapply(signal_list,function(x){length(intersect(kdc_gene,x))})

#' KDC only genes total expression
#+ echo=F, warnings=F, message=F,eval=T
apply(sp_count[kdc_only,],1,sum)

#' Number of Spots that KDC only genes having expression at 50X resolution
#+ echo=F, warnings=F, message=F,eval=T
apply(sp_count[kdc_only,],1,function(x){sum(x!=0)})

