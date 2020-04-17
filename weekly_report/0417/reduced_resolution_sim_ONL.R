#' ---
#' title: "HDST vs. 38X (Layer ONL)"
#' author: "Jiaqiang Zhu"
#' date: "April 17th, 2020"
#' ---


#' Three Patterns Examined for Power Comparison
#+ label=Pattern, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")


realpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/data/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

load(paste0(realpath,"/hdst/",isample,"_barcodes_clean.rds"))
xys <- cbind(
    x=as.numeric(sapply(strsplit(as.character(barcode_with_layers$bc),split="x"),"[[",1)),
    y=as.numeric(sapply(strsplit(as.character(barcode_with_layers$bc),split="x"),"[[",2))
    )

new_barcode_layer <- cbind.data.frame(x=xys[,1],y=xys[,2],layer=barcode_with_layers$abbr)

pltdat <- cbind.data.frame(x=xys[,1],y=xys[,2],
	E=as.numeric(barcode_with_layers$abbr=="E"),
	ONL=as.numeric(barcode_with_layers$abbr=="ONL"),
	GL=as.numeric(barcode_with_layers$abbr=="GL"))

pp    <- lapply(1:3,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})


library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],
			labels = c("A", "B", "C"),font.label=list(size=15),ncol = 3, nrow = 1)

#' ***
#' > Figure 0: Three Patterns (A): Layer E (B): Layer ONL  (C): Layer GL

#+ label=fig0,fig.width=9, fig.height=3, fig.cap="Pattern Plots",echo=F,fig.align="center",eval=T
pp_list


#+ label=kdc, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"
ilayer    	<- 2
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
layerlist 	<- c("E","ONL","GL")
p1 <- p2 <- p3 <- c()
for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_mtest[,1])	
	rm(KDC)
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_zerorm_binR_rpt",isid,".rds"))	
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(total_adj=as.vector(p1[101:1000,]),
							no_adj=as.vector(p2[101:1000,])),
							cl=0,
							legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("total_adj","no_adj"),
					        self_label=c("total_adj","no_adj"))
fig1 		<- comp_qq + ggtitle("KDC_NULL") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(list(total_adj=as.vector(p1),
							no_adj=as.vector(p2)),
							cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("total_adj","no_adj"),
					        self_label=c("total_adj","no_adj"))
fig2 		<- comp_qq2 + ggtitle("KDC_ALL") + theme(plot.title = element_text(hjust = 0.5,size=15))


#+ label=sparklm, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 2
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
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(RL=as.vector(p1[101:1000,]),SW=as.vector(p2[101:1000,])),cl=0,
							legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("RL","SW"),
					        self_label=c("RL","SW"))

fig3 		<- comp_qq + ggtitle("SPARK_LMM_NULL") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(list(RL=as.vector(p1),SW=as.vector(p2)),
							cl=0,legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("RL","SW"),
					        self_label=c("RL","SW"))

fig4 		<- comp_qq2 + ggtitle("SPARK_LMM_ALL") + theme(plot.title = element_text(hjust = 0.5,size=15))


#+ label=sparkpoi, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2","fig3","fig4")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 2
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
# isid    	<- 2

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_count_rpt",isid,".rds"))	
	p2 <- cbind(p2,spark_pval$combined_pvalue)	
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig5 		<- comp_qq + ggtitle("SPARK_PMM_NULL") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(as.vector(p2[1:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig6 		<- comp_qq2 + ggtitle("SPARK_PMM_ALL") + theme(plot.title = element_text(hjust = 0.5,size=15))


#+ label=spatialDE, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2","fig3","fig4","fig5","fig6")))
library(data.table)
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spatialDE/"

ilayer    	<- 2
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")
p1 <- p2 <- c()
for(irpt in 1:10){
	res <- read.table(paste0(outpath,"/sim_spe_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_rpt",irpt,".csv"),head=T)
	res1 <- res[-which(res$g=="log_total_count"),]
	if(sum(duplicated(res1$g))>0){
		clean_res <- res1[-which(duplicated(res1$g)),]
	}else{
		clean_res <- res1
	}
	order_res <- clean_res[order(as.numeric(sapply(strsplit(as.character(clean_res$g),split="gene"),"[[",2))),]
	p2 <- cbind(p2,order_res$pval)	
	rm(res,res1,order_res)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq        <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=1)
fig7 			<- comp_qq + ggtitle("SpatialDE_NULL") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     	<- qplot_gg(as.vector(p2[1:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig8 			<- comp_qq2 + ggtitle("SpatialDE_ALL") + theme(plot.title = element_text(hjust = 0.5,size=15))

#' ***
#' > Figure 1: 38X Resolution QQ Plots, ONL Layer
#+ label=figs,fig.width=12, fig.height=24, echo=F,fig.align="center"
all_figs <- ggarrange(fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8,
							ncol=2,nrow=4)
all_figs


#+ label=kdc1x, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"
ilayer    	<- 2
iblock    	<- 1
numZero 	<- 2000
sigZero 	<- 200
layerlist 	<- c("E","ONL","GL")
p1 <- p2 <- p3 <- c()
for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_mtest[,1])	
	rm(KDC)
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_zerorm_binR_rpt",isid,".rds"))	
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(total_adj=as.vector(p1[101:1000,]),
							no_adj=as.vector(p2[101:1000,])),
							cl=0,
							legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("total_adj","no_adj"),
					        self_label=c("total_adj","no_adj"))
fig1 		<- comp_qq + ggtitle("KDC_NULL") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(list(total_adj=as.vector(p1),
							no_adj=as.vector(p2)),
							cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("total_adj","no_adj"),
					        self_label=c("total_adj","no_adj"))
fig2 		<- comp_qq2 + ggtitle("KDC_ALL") + theme(plot.title = element_text(hjust = 0.5,size=15))

library(ggpubr)
pp_list <- ggarrange(fig1,fig2,ncol = 2, nrow = 1)

#' ***
#' > Figure 2: 1X Resolution QQ Plots, ONL Layer
#+ label=x1qq,fig.width=12, fig.height=6, echo=F,fig.align="center"
pp_list





#+ label=patHDST, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")

source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 2
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_HDST_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))


idx <- which(sp_sums_Rcpp(sp_sim_count)==0)
raw_loc  	<- as.matrix(location[,1:2])


if(length(idx)!=0){
	sp_sim_count_rm <- sp_sim_count[,-which(sp_sums_Rcpp(sp_sim_count)==0)]
	location_rm <- raw_loc[-which(sp_sums_Rcpp(sp_sim_count)==0),]
}else{
	sp_sim_count_rm <- sp_sim_count
	location_rm <- raw_loc
}
sp_sim_count_adj <- t(t(sp_sim_count_rm)/as.vector(sp_sums_Rcpp(sp_sim_count_rm)))

pltdat <- cbind.data.frame(x=location_rm[,1],y=location_rm[,2],
			t(as.matrix(sp_sim_count_adj[c(1:3,101:103),])))

pp    <- lapply(1:6,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})

library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],
					pp[[4]], pp[[5]], pp[[6]],ncol = 3, nrow = 2)

#' ***
#' > Figure 3: HDST Resolution Gene Expression Pattern
#+ label=fig3,fig.width=9, fig.height=6, echo=F,fig.align="center",eval=T
pp_list



#' ***
#' KDC HDST *p*-values
#+ label=patkdc1x, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

ilayer    	<- 2
iblock    	<- 1
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_stest[,1])
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}
rownames(p2) <- paste0("gene",1:1000)
p2[c(1:3,101:103),10]


#+ label=X38, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")
library(Rcpp)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 2
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_binR_rpt",isid,".rds"))

idx <- which(sp_sums_Rcpp(sp_sim_count)==0)
raw_loc  	<- as.matrix(location[,1:2])


if(length(idx)!=0){
	sp_sim_count_rm <- sp_sim_count[,-which(sp_sums_Rcpp(sp_sim_count)==0)]
	location_rm <- raw_loc[-which(sp_sums_Rcpp(sp_sim_count)==0),]
}else{
	sp_sim_count_rm <- sp_sim_count
	location_rm <- raw_loc
}
sp_sim_count_adj <- t(t(sp_sim_count_rm)/as.vector(sp_sums_Rcpp(sp_sim_count_rm)))

pltdat <- cbind.data.frame(x=location_rm[,1],y=location_rm[,2],
			t(as.matrix(sp_sim_count_adj[c(1:3,101:103),])))

pp    <- lapply(1:6,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=4,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})


library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],
					pp[[4]], pp[[5]], pp[[6]],ncol = 3, nrow = 2)


#' ***
#' > Figure 4: 38X Resolution Gene Expression Pattern
#+ label=fig4,fig.width=9, fig.height=6, echo=F,fig.align="center",eval=T
pp_list



#' ***
#' KDC 38X *p*-values 
#+ label=patkdc38x, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

ilayer    	<- 2
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_stest[,1])
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}
rownames(p2) <- paste0("gene",1:1000)
p2[c(1:3,101:103),10]

