#' ---
#' title: "KDC QQ Plots (TotalCount Adjusted)"
#' author: "Jiaqiang Zhu"
#' date: "April 16th, 2020"
#' ---


#+ label=one, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"
ilayer    	<- 1
iblock    	<- 38
numZero 	<- 1000
sigZero 	<- 100
layerlist 	<- c("E","ONL","GL")
p1 <- p2 <- p3 <- c()
for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_mtest[,1])	
	rm(KDC)
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,"_perm.rds"))	
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(NULL=as.vector(p1[101:1000,]),
							Permuted_NULL=as.vector(p2[101:1000,])),
							cl=0,
							legend.position=c(0.25,0.85),pt.size=1,
                            ax.txt.size=10,ax.title.size=12,
                            len.txt.size=0.8,
                            factor_level=c("NULL","Permuted_NULL"),
					        self_label=c("NULL","PermNULL"))
fig1 		<- comp_qq + ggtitle(paste0(numZero,"/",sigZero)) + theme(plot.title = element_text(hjust = 0.5,size=12))

#+ label=two, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"fig1"))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"
ilayer    	<- 1
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
layerlist 	<- c("E","ONL","GL")
p1 <- p2 <- p3 <- c()
for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_mtest[,1])	
	rm(KDC)
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,"_perm.rds"))	
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(NULL=as.vector(p1[101:1000,]),
							Permuted_NULL=as.vector(p2[101:1000,])),
							cl=0,
							legend.position="none",pt.size=1,
                            ax.txt.size=10,ax.title.size=12,
                            len.txt.size=1,
                            factor_level=c("NULL","Permuted_NULL"),
					        self_label=c("NULL","Permuted_NULL"))
fig2 		<- comp_qq + ggtitle(paste0(numZero,"/",sigZero)) + theme(plot.title = element_text(hjust = 0.5,size=12))

#+ label=five, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"
ilayer    	<- 1
iblock    	<- 38
numZero 	<- 5000
sigZero 	<- 500
layerlist 	<- c("E","ONL","GL")
p1 <- p2 <- p3 <- c()
for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_mtest[,1])	
	rm(KDC)
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,"_perm.rds"))	
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(list(NULL=as.vector(p1[101:1000,]),
							Permuted_NULL=as.vector(p2[101:1000,])),
							cl=0,
							legend.position="none",pt.size=1,
                            ax.txt.size=10,ax.title.size=12,
                            len.txt.size=1,
                            factor_level=c("NULL","Permuted_NULL"),
					        self_label=c("NULL","Permuted_NULL"))
fig3 		<- comp_qq + ggtitle(paste0(numZero,"/",sigZero)) + theme(plot.title = element_text(hjust = 0.5,size=12))

library(ggpubr)
pp_list <- ggarrange(fig1,fig2,fig3,ncol = 3, nrow = 1)

#' ***
#' > Figure 1: 38X, NULL vs. Permuted NULL
#+ label=fig1,fig.width=9, fig.height=3, echo=F,fig.align="center",eval=T
pp_list







