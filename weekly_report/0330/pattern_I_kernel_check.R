#' ---
#' title: "Pattern I with SPARK Power Detail Check"
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

#' ***
#' **Check if it is the Gaussian kernel playing a role in the pattern 1**

#' *50X with numZero=2000*

#+ label=X50, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"sim_power_func"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result/"
numZero = 2000
sigZero = 200

layerlist <- c("E","ONL","GL")
iblock = 50
ilayer = 1
isid = 1
load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,"_nofilter.rds"))

spark_pval[95:104,]



#' *50X with numZero=2000*

#+ label=X25, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"sim_power_func"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result/"
numZero = 2000
sigZero = 200

layerlist <- c("E","ONL","GL")
iblock = 25
ilayer = 1
isid = 1
load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,"_nofilter.rds"))

spark_pval[95:104,]

