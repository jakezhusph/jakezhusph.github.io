#' ---
#' title: "QQ Comparison"
#' author: "Jiaqiang Zhu"
#' date: "April 30th, 2020"
#' ---


#+ label=moderate1, echo=F, warnings=F, message=F,eval=T

rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/"


itheta = 5
mu0 = 0.5

col_base <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
kdc_pval <- pmm_pval <- lmm_sw_pval_vst <- spe_pval <- list()

isample = 300

p_kdc <- p_pmm <- p_lmm_sw_vst <- p_spe <-  c()

for(irpt in 1:10){
	# KDC
	load(paste0(outpath,"/kdc/sim_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))
	p_kdc <- cbind(p_kdc,KDC$res_mtest$combinedPval)
	rm(KDC)


	load(paste0(outpath,"/spark/sim_spark_pmm_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))
	p_pmm <- cbind(p_pmm,spark_pval$combined_pvalue)
	rm(spark_pval)

	load(paste0(outpath,"/spark/sim_spark_lmm_null_SS",isample,"_mu",mu0,"_theta",itheta,"_vst_rpt",irpt,".rds"))
	p_lmm_sw_vst <- cbind(p_lmm_sw_vst,CombinePValues(SW_pval))
	rm(SW_pval)

	# spatialDE
	res <- read.table(paste0(outpath,"/spatialde/sim_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
	res1 <- res[-which(res$g=="log_total_count"),]

	if(sum(duplicated(res1$g))>0){
		clean_res <- res1[-which(duplicated(res1$g)),]
	}else{
		clean_res <- res1
	}

	order_res <- clean_res[order(as.numeric(sapply(strsplit(as.character(clean_res$g),split="gene"),"[[",2))),]
	p_spe <- cbind(p_spe,order_res$pval)
	rm(res,res1,clean_res,order_res)
}


source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")
fig1     <- qplot_gg(list(KDC=p_kdc,
							LMM_SW_VST=p_lmm_sw_vst,
							SpatialDE=p_spe,
							SPARK=p_pmm),
							cl=0,
							col.base=col_base,
							legend.position=c(0.25,0.8),pt.size=2,
	                        ax.txt.size=15,ax.title.size=15,
	                        len.txt.size=1,
	                        factor_level=c("KDC","LMM_SW_VST","SpatialDE","SPARK"),
					        self_label=c("SPARK-X","SPARK-G","SpatialDE","SPARK"))
fig1 <- fig1+ ggtitle("SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))

#+ label=moderate2, echo=F, warnings=F, message=F,eval=T

rm(list=setdiff(ls(),"fig1"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/"


itheta = 5
mu0 = 0.5

col_base <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
kdc_pval <- pmm_pval <- lmm_sw_pval_vst <- spe_pval <- list()

isample = 1000

p_kdc <- p_pmm <- p_lmm_sw_vst <- p_spe <-  c()

for(irpt in 1:10){
	# KDC
	load(paste0(outpath,"/kdc/sim_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))
	p_kdc <- cbind(p_kdc,KDC$res_mtest$combinedPval)
	rm(KDC)


	load(paste0(outpath,"/spark/sim_spark_pmm_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))
	p_pmm <- cbind(p_pmm,spark_pval$combined_pvalue)
	rm(spark_pval)

	load(paste0(outpath,"/spark/sim_spark_lmm_null_SS",isample,"_mu",mu0,"_theta",itheta,"_vst_rpt",irpt,".rds"))
	p_lmm_sw_vst <- cbind(p_lmm_sw_vst,CombinePValues(SW_pval))
	rm(SW_pval)

	# spatialDE
	res <- read.table(paste0(outpath,"/spatialde/sim_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
	res1 <- res[-which(res$g=="log_total_count"),]

	if(sum(duplicated(res1$g))>0){
		clean_res <- res1[-which(duplicated(res1$g)),]
	}else{
		clean_res <- res1
	}

	order_res <- clean_res[order(as.numeric(sapply(strsplit(as.character(clean_res$g),split="gene"),"[[",2))),]
	p_spe <- cbind(p_spe,order_res$pval)
	rm(res,res1,clean_res,order_res)
}


source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")
fig2     <- qplot_gg(list(KDC=p_kdc,
							LMM_SW_VST=p_lmm_sw_vst,
							SpatialDE=p_spe,
							SPARK=p_pmm),
							cl=0,
							col.base=col_base,
							legend.position=c(0.25,0.8),pt.size=2,
	                        ax.txt.size=15,ax.title.size=15,
	                        len.txt.size=1,
	                        factor_level=c("KDC","LMM_SW_VST","SpatialDE","SPARK"),
					        self_label=c("SPARK-X","SPARK-G","SpatialDE","SPARK"))

fig2 <- fig2+ ggtitle("SS1000") + theme(plot.title = element_text(hjust = 0.5,size=15))

library(ggpubr)
pp1 <- ggarrange(fig1,fig2+rremove("legend"),
	labels = c("A","B"),font.label=list(size=15),ncol = 2, nrow = 1)
	
#' ***
#' > Figure 1: QQ plots, moderate expression (mu=0.5). (A) SS300  (B) SS1000
#+ label=pp1,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
pp1




#+ label=sparse1, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/"


itheta = 0.4
mu0 = 0.005

col_base <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
kdc_pval <- pmm_pval <- lmm_sw_pval_vst <- spe_pval <- list()

isample = 3000

p_kdc <- p_pmm <- p_lmm_sw_vst <- p_spe <-  c()

for(irpt in 1:10){
	# KDC
	load(paste0(outpath,"/kdc/sim_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))
	p_kdc <- cbind(p_kdc,KDC$res_mtest$combinedPval)
	rm(KDC)


	# load(paste0(outpath,"/spark/sim_spark_pmm_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))
	# p_pmm <- cbind(p_pmm,spark_pval$combined_pvalue)
	# rm(spark_pval)

	load(paste0(outpath,"/spark/sim_spark_lmm_null_SS",isample,"_mu",mu0,"_theta",itheta,"_vst_rpt",irpt,".rds"))
	p_lmm_sw_vst <- cbind(p_lmm_sw_vst,CombinePValues(SW_pval))
	rm(SW_pval)

	# spatialDE
	res <- read.table(paste0(outpath,"/spatialde/sim_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
	res1 <- res[-which(res$g=="log_total_count"),]

	if(sum(duplicated(res1$g))>0){
		clean_res <- res1[-which(duplicated(res1$g)),]
	}else{
		clean_res <- res1
	}

	order_res <- clean_res[order(as.numeric(sapply(strsplit(as.character(clean_res$g),split="gene"),"[[",2))),]
	p_spe <- cbind(p_spe,order_res$pval)
	rm(res,res1,clean_res,order_res)
}



source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")
fig1     <- qplot_gg(list(KDC=p_kdc,
							LMM_SW_VST=p_lmm_sw_vst,
							SpatialDE=p_spe),
							cl=0,
							col.base=col_base,
							legend.position=c(0.25,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1,
                            factor_level=c("KDC","LMM_SW_VST","SpatialDE"),
					        self_label=c("SPARK-X","SPARK-G","SpatialDE"))

fig1 <- fig1+ ggtitle("SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))


#+ label=sparse2, echo=F, warnings=F, message=F,eval=T

rm(list=setdiff(ls(),c("fig1")))
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/"


itheta = 0.4
mu0 = 0.005

col_base <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
kdc_pval <- pmm_pval <- lmm_sw_pval_vst <- spe_pval <- list()

isample = 10000

p_kdc <- p_pmm <- p_lmm_sw_vst <- p_spe <-  c()

for(irpt in 1:10){
	# KDC
	load(paste0(outpath,"/kdc/sim_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))
	p_kdc <- cbind(p_kdc,KDC$res_mtest$combinedPval)
	rm(KDC)


	# load(paste0(outpath,"/spark/sim_spark_pmm_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))
	# p_pmm <- cbind(p_pmm,spark_pval$combined_pvalue)
	# rm(spark_pval)

	load(paste0(outpath,"/spark/sim_spark_lmm_null_SS",isample,"_mu",mu0,"_theta",itheta,"_vst_rpt",irpt,".rds"))
	p_lmm_sw_vst <- cbind(p_lmm_sw_vst,CombinePValues(SW_pval))
	rm(SW_pval)

	# spatialDE
	res <- read.table(paste0(outpath,"/spatialde/sim_null_SS",isample,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
	res1 <- res[-which(res$g=="log_total_count"),]

	if(sum(duplicated(res1$g))>0){
		clean_res <- res1[-which(duplicated(res1$g)),]
	}else{
		clean_res <- res1
	}

	order_res <- clean_res[order(as.numeric(sapply(strsplit(as.character(clean_res$g),split="gene"),"[[",2))),]
	p_spe <- cbind(p_spe,order_res$pval)
	rm(res,res1,clean_res,order_res)
}



source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")
fig2     <- qplot_gg(list(KDC=p_kdc,
							LMM_SW_VST=p_lmm_sw_vst,
							SpatialDE=p_spe),
							cl=0,
							col.base=col_base,
							legend.position=c(0.25,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1,
                            factor_level=c("KDC","LMM_SW_VST","SpatialDE"),
					        self_label=c("SPARK-X","SPARK-G","SpatialDE"))
fig2 <- fig2+ ggtitle("SS10000") + theme(plot.title = element_text(hjust = 0.5,size=15))



library(ggpubr)
pp2 <- ggarrange(fig1+rremove("legend"),fig2+rremove("legend"),labels = c("A","B"),font.label=list(size=15),ncol = 2, nrow = 1)

#' ***
#' > Figure 2: QQ plots, sparse expression (mu=0.005). (A) SS3000  (B) SS10000
#+ label=pp2,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
pp2

