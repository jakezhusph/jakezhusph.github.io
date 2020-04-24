#' ---
#' title: "QQ Comparison for Simulation II"
#' author: "Jiaqiang Zhu"
#' date: "April 20th, 2020"
#' ---


#+ label=high, echo=F, warnings=F, message=F,eval=T

rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/"
cellProp = 0.2
ieffect = 2
itheta = 10
ipat = 1
mu0 = 1.5

col_base <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
pattern_list 	<- c("streak","hotspot")
kdc_pval 		<- pmm_pval <- lmm_sw_pval_qn <- lmm_RL_pval_qn <- lmm_sw_pval_vst <- lmm_RL_pval_vst <- spe_pval <- list()

icount 			<- 0 
for(isample in c(300,1000,3000)){
	icount = icount + 1
	p_kdc <- p_pmm <- p_lmm_sw_qn <- p_lmm_RL_qn <- p_lmm_sw_vst <- p_lmm_RL_vst <- p_spe <-  c()
	
	# cat("SS",isample,", mu:",mu0,",",pattern_list[ipat],"\n")
	for(irpt in 1:10){
		# KDC
		load(paste0(outpath,"/kdc/sim_kdc_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))
		p_kdc <- cbind(p_kdc,KDC$res_mtest$combinedPval)
		rm(KDC)

		# PMM
		load(paste0(outpath,"/spark/sim_spark_pmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))
		p_pmm <- cbind(p_pmm,spark_pval$combined_pvalue)
		rm(spark_pval)

		# LMM--QN
		load(paste0(outpath,"/spark/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_rpt",irpt,".rds"))
		p_lmm_sw_qn <- cbind(p_lmm_sw_qn,CombinePValues(SW_pval))
		p_lmm_RL_qn <- cbind(p_lmm_RL_qn,CombinePValues(RL_pval))

		rm(RL_pval,SW_pval)

		# LMM--VST
		load(paste0(outpath,"/spark/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_rpt",irpt,".rds"))
		p_lmm_sw_vst <- cbind(p_lmm_sw_vst,CombinePValues(SW_pval))
		p_lmm_RL_vst <- cbind(p_lmm_RL_vst,CombinePValues(RL_pval))
		rm(RL_pval,SW_pval)


		# spatialDE
		res <- read.table(paste0(outpath,"/spatialde/sim_spe_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
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


	kdc_pval[[icount]] <- p_kdc[101:1000,]
	pmm_pval[[icount]] <- p_pmm[101:1000,]
	lmm_sw_pval_qn[[icount]] <- p_lmm_sw_qn[101:1000,]
	lmm_RL_pval_qn[[icount]] <- p_lmm_RL_qn[101:1000,]

	lmm_sw_pval_vst[[icount]] <- p_lmm_sw_vst[101:1000,]
	lmm_RL_pval_vst[[icount]] <- p_lmm_RL_vst[101:1000,]

	spe_pval[[icount]] <- p_spe[101:1000,]
	rm(p_kdc,p_pmm,p_lmm_sw_qn,p_lmm_RL_qn,p_lmm_sw_vst,p_lmm_RL_vst,p_spe)
}


source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq1     <- qplot_gg(list(KDC=kdc_pval[[1]],
							PMM=pmm_pval[[1]],
							LMM_SW_QN=lmm_sw_pval_qn[[1]],
							LMM_RL_QN=lmm_RL_pval_qn[[1]],
							LMM_SW_VST=lmm_sw_pval_vst[[1]],
							LMM_RL_VST=lmm_RL_pval_vst[[1]],
							SpatialDE=spe_pval[[1]]),
							cl=0,
							col.base=col_base,
							legend.position=c(0.25,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1,
                            factor_level=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"),
					        self_label=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"))

comp_qq2     <- qplot_gg(list(KDC=kdc_pval[[2]],
							PMM=pmm_pval[[2]],
							LMM_SW_QN=lmm_sw_pval_qn[[2]],
							LMM_RL_QN=lmm_RL_pval_qn[[2]],
							LMM_SW_VST=lmm_sw_pval_vst[[2]],
							LMM_RL_VST=lmm_RL_pval_vst[[2]],
							SpatialDE=spe_pval[[2]]),
							cl=0,
							col.base=col_base,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"),
					        self_label=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"))


comp_qq3     <- qplot_gg(list(KDC=kdc_pval[[3]],
							PMM=pmm_pval[[3]],
							LMM_SW_QN=lmm_sw_pval_qn[[3]],
							LMM_RL_QN=lmm_RL_pval_qn[[3]],
							LMM_SW_VST=lmm_sw_pval_vst[[3]],
							LMM_RL_VST=lmm_RL_pval_vst[[3]],
							SpatialDE=spe_pval[[3]]),
							cl=0,
							col.base=col_base,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"),
					        self_label=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"))


fig1 <- comp_qq1 + ggtitle("SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))
fig2 <- comp_qq2 + ggtitle("SS1000") + theme(plot.title = element_text(hjust = 0.5,size=15))
fig3 <- comp_qq3 + ggtitle("SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))

library(ggpubr)
pp1  <- ggarrange(fig1,fig2,fig3,labels = c("A", "B","C"),font.label=list(size=15),ncol = 3, nrow = 1)

#' ***
#' > Figure 1: mu = 1.5
#+ label=fig1,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
pp1



#+ label=moderate, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/"
cellProp = 0.2
ieffect = 2
itheta = 10
ipat = 1
mu0 = 0.5

col_base <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
pattern_list 	<- c("streak","hotspot")
kdc_pval 		<- pmm_pval <- lmm_sw_pval_qn <- lmm_RL_pval_qn <- lmm_sw_pval_vst <- lmm_RL_pval_vst <- spe_pval <- list()

icount 			<- 0 
for(isample in c(300,1000,3000)){
	icount = icount + 1
	p_kdc <- p_pmm <- p_lmm_sw_qn <- p_lmm_RL_qn <- p_lmm_sw_vst <- p_lmm_RL_vst <- p_spe <-  c()
	
	# cat("SS",isample,", mu:",mu0,",",pattern_list[ipat],"\n")
	for(irpt in 1:10){
		# KDC
		load(paste0(outpath,"/kdc/sim_kdc_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))
		p_kdc <- cbind(p_kdc,KDC$res_mtest$combinedPval)
		rm(KDC)

		# PMM
		load(paste0(outpath,"/spark/sim_spark_pmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))
		p_pmm <- cbind(p_pmm,spark_pval$combined_pvalue)
		rm(spark_pval)

		# LMM--QN
		load(paste0(outpath,"/spark/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_rpt",irpt,".rds"))
		p_lmm_sw_qn <- cbind(p_lmm_sw_qn,CombinePValues(SW_pval))
		p_lmm_RL_qn <- cbind(p_lmm_RL_qn,CombinePValues(RL_pval))

		rm(RL_pval,SW_pval)

		# LMM--VST
		load(paste0(outpath,"/spark/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_rpt",irpt,".rds"))
		p_lmm_sw_vst <- cbind(p_lmm_sw_vst,CombinePValues(SW_pval))
		p_lmm_RL_vst <- cbind(p_lmm_RL_vst,CombinePValues(RL_pval))
		rm(RL_pval,SW_pval)


		# spatialDE
		res <- read.table(paste0(outpath,"/spatialde/sim_spe_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
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


	kdc_pval[[icount]] <- p_kdc[101:1000,]
	pmm_pval[[icount]] <- p_pmm[101:1000,]
	lmm_sw_pval_qn[[icount]] <- p_lmm_sw_qn[101:1000,]
	lmm_RL_pval_qn[[icount]] <- p_lmm_RL_qn[101:1000,]

	lmm_sw_pval_vst[[icount]] <- p_lmm_sw_vst[101:1000,]
	lmm_RL_pval_vst[[icount]] <- p_lmm_RL_vst[101:1000,]

	spe_pval[[icount]] <- p_spe[101:1000,]
	rm(p_kdc,p_pmm,p_lmm_sw_qn,p_lmm_RL_qn,p_lmm_sw_vst,p_lmm_RL_vst,p_spe)
}


source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq1     <- qplot_gg(list(KDC=kdc_pval[[1]],
							PMM=pmm_pval[[1]],
							LMM_SW_QN=lmm_sw_pval_qn[[1]],
							LMM_RL_QN=lmm_RL_pval_qn[[1]],
							LMM_SW_VST=lmm_sw_pval_vst[[1]],
							LMM_RL_VST=lmm_RL_pval_vst[[1]],
							SpatialDE=spe_pval[[1]]),
							cl=0,
							col.base=col_base,
							legend.position=c(0.25,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1,
                            factor_level=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"),
					        self_label=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"))

comp_qq2     <- qplot_gg(list(KDC=kdc_pval[[2]],
							PMM=pmm_pval[[2]],
							LMM_SW_QN=lmm_sw_pval_qn[[2]],
							LMM_RL_QN=lmm_RL_pval_qn[[2]],
							LMM_SW_VST=lmm_sw_pval_vst[[2]],
							LMM_RL_VST=lmm_RL_pval_vst[[2]],
							SpatialDE=spe_pval[[2]]),
							cl=0,
							col.base=col_base,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"),
					        self_label=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"))


comp_qq3     <- qplot_gg(list(KDC=kdc_pval[[3]],
							PMM=pmm_pval[[3]],
							LMM_SW_QN=lmm_sw_pval_qn[[3]],
							LMM_RL_QN=lmm_RL_pval_qn[[3]],
							LMM_SW_VST=lmm_sw_pval_vst[[3]],
							LMM_RL_VST=lmm_RL_pval_vst[[3]],
							SpatialDE=spe_pval[[3]]),
							cl=0,
							col.base=col_base,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"),
					        self_label=c("KDC","PMM","LMM_SW_QN","LMM_RL_QN","LMM_SW_VST","LMM_RL_VST","SpatialDE"))

fig1 <- comp_qq1 + ggtitle("SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))
fig2 <- comp_qq2 + ggtitle("SS1000") + theme(plot.title = element_text(hjust = 0.5,size=15))
fig3 <- comp_qq3 + ggtitle("SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))

library(ggpubr)
pp2  <- ggarrange(fig1,fig2,fig3,labels = c("A", "B","C"),font.label=list(size=15),ncol = 3, nrow = 1)

#' ***
#' > Figure 2: mu = 0.5
#+ label=fig2,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
pp2

