#' ---
#' title: "Puck 180430_6 Reduced Resolution"
#' author: "Jiaqiang Zhu"
#' date: "March 22rd, 2020"
#' ---


#+ label=default, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
library(SPARK)

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"

isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	load(paste0(workdir,"/data/slideseq/clean_",isample,"_idx_",iblock,"X_sparseCount.rds"))
	raw_loc  <- as.data.frame(location[,1:2])
	spark <- CreateSPARKObject(counts = sp_count, location = raw_loc, percentage = 0.1, min_total_counts = 10)
	spark@lib_size <- apply(spark@counts, 2, sum)
	# spark@counts <- spark@counts[1:10,]
	t1 		<- proc.time()
	spark 	<- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
	    num_core = 10, verbose = T)
	t2 		<- proc.time()- t1

	spark <- spark.test(spark,  check_positive=T,verbose = T)
	t3 		<- proc.time()- t1

	runTimeAll = t3 
	runTimeTest = t3 - t2
	save(spark,runTimeAll,runTimeTest,file=paste0(workdir,"/output/spark/clean_",isample,"_idx_",iblock,"X_spark.rds"))
	rm(list=setdiff(ls(),c("workdir","isample")))
}

#+ label=default_nomt, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
library(SPARK)

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"

isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	load(paste0(workdir,"/data/slideseq/clean_",isample,"_idx_",iblock,"X_sparseCount.rds"))
	raw_loc  <- as.data.frame(location[,1:2])
	mtgene <- which(substr(rownames(sp_count),start=1,stop=3)=="mt-")
	if(length(mtgene)>0){
		sp_count <- sp_count[-mtgene,]
	} 

	spark <- CreateSPARKObject(counts = sp_count, location = raw_loc, percentage = 0.1, min_total_counts = 10)
	spark@lib_size <- apply(spark@counts, 2, sum)
	# spark@counts <- spark@counts[1:10,]
	t1 		<- proc.time()
	spark 	<- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
	    num_core = 10, verbose = T)
	t2 		<- proc.time()- t1

	spark <- spark.test(spark,  check_positive=T,verbose = T)
	t3 		<- proc.time()- t1

	runTimeAll = t3 
	runTimeTest = t3 - t2
	save(spark,runTimeAll,runTimeTest,file=paste0(workdir,"/output/spark/clean_",isample,"_idx_",iblock,"X_spark_nomt.rds"))
	rm(list=setdiff(ls(),c("workdir","isample")))
}


#+ label=default_nomt_per, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
library(SPARK)
workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
args    	<- as.numeric(commandArgs(TRUE))
iblock    	<- args[1]
isid    	<- args[2]

isample = "Puck_180430_6"
load(paste0(workdir,"/data/slideseq/clean_",isample,"_idx_",iblock,"X_sparseCount.rds"))
raw_loc  <- as.data.frame(location[,1:2])
mtgene <- which(substr(rownames(sp_count),start=1,stop=3)=="mt-")
if(length(mtgene)>0){
	sp_count <- sp_count[-mtgene,]
} 
set.seed(isid)
idx 			<- sample(1:ncol(sp_count))
perm_count 		<- sp_count[,idx]
colnames(perm_count) <- colnames(sp_count)
spark <- CreateSPARKObject(counts = perm_count, location = raw_loc, percentage = 0.1, min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
t1 		<- proc.time()
spark 	<- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
    num_core = 1, verbose = T)
t2 		<- proc.time()- t1
spark <- spark.test(spark,  check_positive=T,verbose = T)
t3 		<- proc.time()- t1

runTimeAll = t3 
runTimeTest = t3 - t2
save(spark,runTimeAll,runTimeTest,file=paste0(workdir,"/output/spark/clean_",isample,"_idx_",iblock,"X_spark_default_nomt_perm",isid,".rds"))


#' ***
#' **KDC Signals **
#+ label=KDC_signals, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"
for(iblock in c(200,100,50,1)){
	load(paste0(workdir,"/output/kdc/clean_",isample,"_idx_",iblock,"X_Joint_kdc.rds"))
	mtgene <- which(substr(rownames(KDC$res_mtest),start=1,stop=3)=="mt-")
	if(length(mtgene)>0){
		pvals <- KDC$res_mtest$combinedPval[-mtgene]
	} 
	cat("------------------------------------\n")
	cat(iblock,"X: numGene:",nrow(KDC$res_mtest),"\n")
	cat("signal genes:",sum(KDC$res_mtest$adjustedPval<0.05),"\n")
	cat("signal genes nomt:",sum(p.adjust(pvals,method="BY")<0.05),"\n")
	rm(pvals)
}


#' **SPARK with Default Setting: genes < 10% removed, cells < 10 totalcount removed**
#+ label=Default_res, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(SPARK)

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	load(paste0(workdir,"/output/spark/clean_",isample,"_idx_",iblock,"X_spark.rds"))
	cat("------------------------------------\n")
	cat(iblock,"X: sampleSize",ncol(spark@counts),"numGene",nrow(spark@counts),"\n")
	cat("converged genes:",sum(sapply(spark@res_vc,function(x){x$converged})),"\n")
	cat("signal genes:",sum(spark@res_mtest$adjusted_pvalue<0.05),"\n")
	rm(spark)
}

#' ***
#' **SPARK with Default Setting (No Mitochondrial): genes < 10% removed, cells < 10 totalcount removed**
#+ label=Default_res_nomt, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(SPARK)

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	load(paste0(workdir,"/output/spark/clean_",isample,"_idx_",iblock,"X_spark_nomt.rds"))
	cat("------------------------------------\n")
	cat(iblock,"X: sampleSize",ncol(spark@counts),"numGene",nrow(spark@counts),"\n")
	cat("converged genes:",sum(sapply(spark@res_vc,function(x){x$converged})),"\n")
	cat("signal genes:",sum(spark@res_mtest$adjusted_pvalue<0.05),"\n")
	rm(spark)
}


#' ***
#' **SpatialDE with Default Setting: genesSum < 3 removed, cells < 10 totalcount removed**
#+ label=spe_res, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	res <- read.table(paste0(workdir,"/output/spatialDE/clean_",isample,"_idx_",iblock,"X_spe.csv"),head=T)
	cat("------------------------------------\n")
	cat(iblock,"X: sampleSize",res$n[1],"numGene",nrow(res)-1,"\n")
	cat("signal genes:",sum(res$qval<0.05),"\n")
	rm(res)
}


#' ***
#' **SpatialDE with Default Setting (No Mitochondrial): genesSum < 3 removed, cells < 10 totalcount removed**
#+ label=spe_res_nomt, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	res <- read.table(paste0(workdir,"/output/spatialDE/clean_",isample,"_idx_",iblock,"X_spe_nomt.csv"),head=T)
	cat("------------------------------------\n")
	cat(iblock,"X: sampleSize",res$n[1],"numGene",nrow(res)-1,"\n")
	cat("signal genes:",sum(res$qval<0.05),"\n")
	rm(res)
}



#+ label = kdcper, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

sparkpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spark/"
spepath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spatialDE/"
kdcpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/kdc/"

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/result/"
isample = "Puck_180430_6"
for(iblock in c(200,100,50,1)){
	kdc_pval <- c()
	for(isid in 1:10){
		load(paste0(kdcpath,"/clean_",isample,"_idx_",iblock,"X_Joint_kdc_perm",isid,".rds"))
		kdc_pval 	<- cbind(kdc_pval,KDC$res_mtest$combinedPval)
		gn 			<- rownames(KDC$res_mtest)
		rm(KDC)
	}
	rownames(kdc_pval) <- gn
	save(kdc_pval,file=paste0(respath,"/clean_",isample,"_idx_",iblock,"X_KDC_permutation_pval.rds"))
	rm(kdc_pval)
}


#+ label=speper, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

sparkpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spark/"
spepath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spatialDE/"

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/result/"
qq_fig <- list()
isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	spe_pval <- c()
	for(isid in 1:10){
		spe 	<- read.table(paste0(spepath,"/clean_",isample,"_idx_",iblock,"X_spe_perm",isid,"_nomt.csv"),head=T)
		res1 	<- spe[-which(spe$g=="log_total_count"),]
		if(sum(duplicated(res1$g))>0){
			clean_res <- res1[-which(duplicated(res1$g)),]
		}else{
			clean_res <- res1
		}
		spe_pval 	<- cbind(spe_pval,clean_res$pval)
		rm(spe)
	}
	save(spe_pval,file=paste0(respath,"/clean_",isample,"_idx_",iblock,"X_spe_permutation_pval_nomt.rds"))
	rm(spe_pval)
}

# make sure the results are not the same, though the pval above are almost same, but the gene order is different 
# isid = 1
# spe1 <- read.table(paste0(spepath,"/clean_CN24_D1_idx_",iblock,"X_spe_perm",isid,"_nomt.csv"),head=T)
# isid = 2
# spe2 <- read.table(paste0(spepath,"/clean_CN24_D1_idx_",iblock,"X_spe_perm",isid,"_nomt.csv"),head=T)



#+ label=sparkper, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
library(SPARK)
sparkpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spark/"
spepath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spatialDE/"

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/result/"
isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	spark_pval <- c()
	for(isid in 1:10){
		load(paste0(sparkpath,"/clean_",isample,"_idx_",iblock,"X_spark_default_nomt_perm",isid,".rds"))
		spark_pval 	<- cbind(spark_pval,spark@res_mtest$combined_pvalue)
		rm(spark)
	}
	save(spark_pval,file=paste0(respath,"/clean_",isample,"_idx_",iblock,"X_spark_permutation_pval_nomt.rds"))
}


#+ label=sparkperC, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
library(SPARK)
sparkpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spark/"
spepath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spatialDE/"

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/result/"
isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	spark_pval <- c()
	for(isid in 1:10){
		load(paste0(sparkpath,"/clean_",isample,"_idx_",iblock,"X_spark_default_nomt_perm",isid,".rds"))
		conver_idx <- sapply(spark@res_vc,function(x){x$converged})
		# can't combined, as the converges changes
		spark_pval 	<- c(spark_pval,spark@res_mtest$combined_pvalue[conver_idx])
		rm(spark)
	}
	save(spark_pval,file=paste0(respath,"/clean_",isample,"_idx_",iblock,"X_spark_permutation_pval_nomt_converged.rds"))
}



#+ label=KDC_QQ, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
sparkpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spark/"
spepath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spatialDE/"

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/result/"

permlist <- list()
icount = 0
isample = "Puck_180430_6"
for(iblock in c(200,100,50,1)){
	load(paste0(respath,"/clean_",isample,"_idx_",iblock,"X_KDC_permutation_pval.rds"))
	icount = icount + 1
	permlist[[icount]] <- kdc_pval
	rm(kdc_pval)
}

names(permlist) <- paste0("X",c(200,100,50,1))
KDC_qq 			<- qplot_gg(permlist,cl=0,legend.position=c(0.15,0.8),pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=paste0("X",c(200,100,50,1)),
					        self_label=paste0("X",c(200,100,50,1)))



#+ label=SPE_QQ, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"KDC_qq"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
sparkpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spark/"
spepath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spatialDE/"

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/result/"

permlist <- list()
icount = 0
isample = "Puck_180430_6"
for(iblock in c(200,100,50)){
	icount = icount + 1
	load(paste0(respath,"/clean_",isample,"_idx_",iblock,"X_spe_permutation_pval_nomt.rds"))
	permlist[[icount]] <- spe_pval
	rm(spe_pval)
}
names(permlist) <- paste0("X",c(200,100,50))
spe_qq 			<- qplot_gg(permlist,cl=0,legend.position=c(0.15,0.8),pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=paste0("X",c(200,100,50)),
					        self_label=paste0("X",c(200,100,50)))


#+ label=spark_QQ, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("KDC_qq","spe_qq")))
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
sparkpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spark/"
spepath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spatialDE/"

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/result/"

permlist <- list()
icount = 0
isample = "Puck_180430_6"
blocklist <- c(200,100,50)
for(iblock in blocklist){
	icount = icount + 1
	load(paste0(respath,"/clean_",isample,"_idx_",iblock,"X_spark_permutation_pval_nomt.rds"))
	permlist[[icount]] <- spark_pval
	rm(spark_pval)
}
names(permlist) <- paste0("X",blocklist)
spark_qq 		<- qplot_gg(permlist,cl=0,legend.position=c(0.15,0.8),pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=paste0("X",blocklist),
					        self_label=paste0("X",blocklist))

#+ label=spark_con_QQ, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("KDC_qq","spe_qq","spark_qq")))
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
sparkpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spark/"
spepath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/output/spatialDE/"

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/result/"

permlist <- list()
icount = 0
isample = "Puck_180430_6"
blocklist <- c(200,100,50)
for(iblock in blocklist){
	icount = icount + 1
	load(paste0(respath,"/clean_",isample,"_idx_",iblock,"X_spark_permutation_pval_nomt_converged.rds"))
	permlist[[icount]] <- spark_pval
	rm(spark_pval)
}
names(permlist) <- paste0("X",blocklist)
spark_con_qq 		<- qplot_gg(permlist,cl=0,legend.position=c(0.15,0.8),pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=paste0("X",blocklist),
					        self_label=paste0("X",blocklist))


library(ggpubr)
fig1 <- ggarrange(KDC_qq, spe_qq, spark_qq,spark_con_qq,
			labels = c("A", "B", "C","D"),
			font.label=list(size=20),
			ncol = 2, nrow = 2)


#' ***
#' > Figure 1: Permutation QQ, A: KDC; B: spatialDE; C: SPARK D: SPARK Converged
#+ label=fig1,fig.width=8, fig.height=8, fig.cap="QQ Plot",echo=F,fig.align="center"
fig1


