#' ---
#' title: "QQ Plots: HDST vs. 50X"
#' author: "Jiaqiang Zhu"
#' date: "March 29th, 2020"
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

ilayer    	<- 3
iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	
	# load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_perm_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_stest[,1])
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig1 		<- comp_qq + ggtitle("KDC") + theme(plot.title = element_text(hjust = 0.5,size=15))


#+ label=kdc_perm, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"fig1"))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

ilayer    	<- 3
iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	# load(paste0(outpath,"/sim_kdc_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_perm_rpt",isid,".rds"))	
	
	p1 <- cbind(p1,KDC$res_stest[,1])
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}


source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig2 		<- comp_qq + ggtitle("KDC_perm") + theme(plot.title = element_text(hjust = 0.5,size=15))


#+ label=spark, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer    	<- 3
iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200
# isid    	<- 2

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	
	p2 <- cbind(p2,spark_pval$combined_pvalue)	
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig3 		<- comp_qq + ggtitle("spark") + theme(plot.title = element_text(hjust = 0.5,size=15))


#+ label=spark_perm, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2","fig3")))
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"

ilayer    	<- 3
iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200
# isid    	<- 2

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_perm_rpt",isid,".rds"))	
	p2 <- cbind(p2,spark_pval$combined_pvalue)	
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig4 		<- comp_qq + ggtitle("spark_perm") + theme(plot.title = element_text(hjust = 0.5,size=15))

#' ***
#' 50X Resolution QQ Plots, GL
#+ label=figs,fig.width=12, fig.height=12, echo=F,fig.align="center"
all_figs <- ggarrange(fig1,fig2,fig3,fig4,
							ncol=2,nrow=2)
all_figs


#+ label=kdc1x, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

ilayer    	<- 3
iblock    	<- 1
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	
	
	p1 <- cbind(p1,KDC$res_stest[,1])
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig2 		<- comp_qq + ggtitle("KDC") + theme(plot.title = element_text(hjust = 0.5,size=15))

#' ***
#' 1X Resolution QQ Plots, GL
#+ label=x1qq,fig.width=6, fig.height=6, echo=F,fig.align="center"
fig2

