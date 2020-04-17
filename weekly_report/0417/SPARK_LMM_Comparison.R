#' ---
#' title: "SPARK LMM Comparison"
#' author: "Jiaqiang Zhu"
#' date: "April 17th, 2020"
#' ---

#' ALL Gaussian Kernels
#+ label=total_qn, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 3
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
# isid    	<- 2

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

source("/net/mulan/disk2/jiaqiang/kdc/code/spark/lmm/sw_lmm/sw_lmm.R")
for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_qn_rpt",isid,".rds"))	
	RL_pval_combined <- CombinePValues(RL_pval)
	sw_pval_combined <- CombinePValues(SW_pval)
	p1 <- cbind(p1,RL_pval_combined)	
	p2 <- cbind(p2,sw_pval_combined)
	# load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_qn_rpt",isid,".rds"))	
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(RL=as.vector(p1[101:1000,]),SW=as.vector(p2[101:1000,])),cl=0,
							legend.position=c(0.2,0.8),pt.size=2,
                            ax.txt.size=12,ax.title.size=12,
                            len.txt.size=1,
                            factor_level=c("RL","SW"),
					        self_label=c("RL","SW"))

fig1 		<- comp_qq + ggtitle("Total_Adjusted_QN") + theme(plot.title = element_text(hjust = 0.5,size=12))

#+ label=total_vst, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 3
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

source("/net/mulan/disk2/jiaqiang/kdc/code/spark/lmm/sw_lmm/sw_lmm.R")
for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_vst_rpt",isid,".rds"))	
	RL_pval_combined <- CombinePValues(RL_pval)
	sw_pval_combined <- CombinePValues(SW_pval)
	p1 <- cbind(p1,RL_pval_combined)	
	p2 <- cbind(p2,sw_pval_combined)
	# load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_qn_rpt",isid,".rds"))	
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(RL=as.vector(p1[101:1000,]),SW=as.vector(p2[101:1000,])),cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=12,ax.title.size=12,
                            len.txt.size=1,
                            factor_level=c("RL","SW"),
					        self_label=c("RL","SW"))

fig2 		<- comp_qq + ggtitle("Total_Adjusted_VST") + theme(plot.title = element_text(hjust = 0.5,size=12))


#+ label=qn, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 3
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
# isid    	<- 2

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

source("/net/mulan/disk2/jiaqiang/kdc/code/spark/lmm/sw_lmm/sw_lmm.R")
for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_qn_rpt",isid,".rds"))	
	RL_pval_combined <- CombinePValues(RL_pval)
	sw_pval_combined <- CombinePValues(SW_pval)
	p1 <- cbind(p1,RL_pval_combined)	
	p2 <- cbind(p2,sw_pval_combined)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(RL=as.vector(p1[101:1000,]),SW=as.vector(p2[101:1000,])),cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=12,ax.title.size=12,
                            len.txt.size=1,
                            factor_level=c("RL","SW"),
					        self_label=c("RL","SW"))

fig3 		<- comp_qq + ggtitle("QN Directly") + theme(plot.title = element_text(hjust = 0.5,size=12))


library(ggpubr)
pp_list <- ggarrange(fig1,fig2, fig3,
			labels = c("A", "B", "C"),font.label=list(size=12),ncol = 3, nrow = 1)

#' ***
#' > Figure 1: LMM With Different Normalization (Ten Gaussian Kernels)

#+ label=fig1,fig.width=9, fig.height=3, fig.cap="QQ Plots",echo=F,fig.align="center",eval=T
pp_list


#' ***
#' Five Gaussian Kernels
#+ label=total_qn_five, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 3
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
# isid    	<- 2

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

source("/net/mulan/disk2/jiaqiang/kdc/code/spark/lmm/sw_lmm/sw_lmm.R")
for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_qn_rpt",isid,".rds"))	
	RL_pval_combined <- CombinePValues(RL_pval[,3:7])
	sw_pval_combined <- CombinePValues(SW_pval[,3:7])
	p1 <- cbind(p1,RL_pval_combined)	
	p2 <- cbind(p2,sw_pval_combined)
	# load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_qn_rpt",isid,".rds"))	
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(RL=as.vector(p1[101:1000,]),SW=as.vector(p2[101:1000,])),cl=0,
							legend.position=c(0.2,0.8),pt.size=2,
                            ax.txt.size=12,ax.title.size=12,
                            len.txt.size=1,
                            factor_level=c("RL","SW"),
					        self_label=c("RL","SW"))

fig1 		<- comp_qq + ggtitle("Total_Adjusted_QN") + theme(plot.title = element_text(hjust = 0.5,size=12))

#+ label=total_vst_five, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 3
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

source("/net/mulan/disk2/jiaqiang/kdc/code/spark/lmm/sw_lmm/sw_lmm.R")
for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_vst_rpt",isid,".rds"))	
	RL_pval_combined <- CombinePValues(RL_pval[,3:7])
	sw_pval_combined <- CombinePValues(SW_pval[,3:7])
	p1 <- cbind(p1,RL_pval_combined)	
	p2 <- cbind(p2,sw_pval_combined)
	# load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_qn_rpt",isid,".rds"))	
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(RL=as.vector(p1[101:1000,]),SW=as.vector(p2[101:1000,])),cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=12,ax.title.size=12,
                            len.txt.size=1,
                            factor_level=c("RL","SW"),
					        self_label=c("RL","SW"))

fig2 		<- comp_qq + ggtitle("Total_Adjusted_VST") + theme(plot.title = element_text(hjust = 0.5,size=12))


#+ label=qn_five, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 3
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
# isid    	<- 2

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

source("/net/mulan/disk2/jiaqiang/kdc/code/spark/lmm/sw_lmm/sw_lmm.R")
for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_qn_rpt",isid,".rds"))	
	RL_pval_combined <- CombinePValues(RL_pval[,3:7])
	sw_pval_combined <- CombinePValues(SW_pval[,3:7])
	p1 <- cbind(p1,RL_pval_combined)	
	p2 <- cbind(p2,sw_pval_combined)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(RL=as.vector(p1[101:1000,]),SW=as.vector(p2[101:1000,])),cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=12,ax.title.size=12,
                            len.txt.size=1,
                            factor_level=c("RL","SW"),
					        self_label=c("RL","SW"))

fig3 		<- comp_qq + ggtitle("QN Directly") + theme(plot.title = element_text(hjust = 0.5,size=12))


library(ggpubr)
pp_list <- ggarrange(fig1,fig2, fig3,
			labels = c("A", "B", "C"),font.label=list(size=12),ncol = 3, nrow = 1)

#' ***
#' > Figure 2: LMM With Different Normalization (Five Kernels)

#+ label=fig2,fig.width=9, fig.height=3, fig.cap="QQ Plots",echo=F,fig.align="center",eval=T
pp_list


