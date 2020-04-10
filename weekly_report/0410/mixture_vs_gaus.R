#' ---
#' title: "All Gaussian vs. Mixture"
#' author: "Jiaqiang Zhu"
#' date: "April 10th, 2020"
#' ---


#+ label=sparklm, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
simpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer      <- 1
iblock      <- 38
numZero     <- 2000
sigZero     <- 200
# isid      <- 2

layerlist   <- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
    load(paste0(outpath,"/sim_spark_mixture_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_means_blocked_binR_lmm_rpt",isid,".rds"))  
    p2 <- cbind(p2,spark_pval$combined_pval)    
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig1        <- comp_qq + ggtitle("SPARK_LMM_NULL") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(as.vector(p2[1:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig2        <- comp_qq2 + ggtitle("SPARK_LMM_ALL") + theme(plot.title = element_text(hjust = 0.5,size=15))



#+ label=sparklmgg, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
simpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
ilayer      <- 1
iblock      <- 38
numZero     <- 2000
sigZero     <- 200
# isid      <- 2

layerlist   <- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
    load(paste0(outpath,"/sim_spark_mixture_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_means_blocked_gaussian_binR_lmm_rpt",isid,".rds")) 
    p2 <- cbind(p2,spark_pval$combined_pval)    
}

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig3        <- comp_qq + ggtitle("SPARK_LMM_Gaussian_NULL") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(as.vector(p2[1:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)
fig4        <- comp_qq2 + ggtitle("SPARK_LMM_Gaussian_ALL") + theme(plot.title = element_text(hjust = 0.5,size=15))

library(ggpubr)
pp_list <- ggarrange(fig1,fig2,fig3,fig4,ncol = 2, nrow = 2)


#' ***
#+ label=fig0,fig.width=12, fig.height=12, fig.cap="QQ Plots",echo=F,fig.align="center",eval=T
pp_list

