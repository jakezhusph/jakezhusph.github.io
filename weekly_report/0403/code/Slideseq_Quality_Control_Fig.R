#' ---
#' title: "Slideseq Quality Control Figures (Puck_180430_6)"
#' author: "Jiaqiang Zhu"
#' date: "April 3rd, 2020"
#' ---

#'***
#' > Variance vs. Mean Figures
#+ label=plot, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
library(fitdistrplus)
library(ggplot2)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/realdata/slideseq/"
load(paste0(workdir,"/data/Puck_180430_6_zero_removed.rds"))

numZero_gene    <- as.vector(ncol(sp_count)- sp_nz_count_Rcpp(sp_count,rowSums=T))
num_nonZero_gene <- as.vector(sp_nz_count_Rcpp(sp_count,rowSums=T))
nonZeroPercent  <- num_nonZero_gene/ncol(sp_count)
zeroPercent     <- 1-nonZeroPercent

mean_gene       <- as.vector(sp_means_Rcpp(sp_count,rowMeans=T))
var_gene        <- as.vector(sp_vars_Rcpp(sp_count,rowVars=T))
mt_gene_idx     <- which(substr(rownames(sp_count),start=1,stop=3)=="mt-")


# based on NaiveDE code: results should be similar
# https://github.com/Teichlab/NaiveDE/blob/master/NaiveDE/base.py
# the start value is not a concern here, 0.1 or 1 leads to the same estimates of phi
phi1 = coef(nls(var_gene ~ mean_gene + phi * mean_gene^2, start = list(phi = 0.1)))
px1 <- seq(min(mean_gene),max(mean_gene),length.out=100)
py1 <- px1 + phi1*px1^2
phi1_inv <- 1/phi1
expZero <- (phi1_inv/(phi1_inv+mean_gene))^phi1_inv



df  <- cbind.data.frame(meanx=mean_gene,varx=var_gene,zeroPercent=zeroPercent,expZ=expZero)
df2 <- cbind.data.frame(px=px1,py=py1)


var_gene_nomt   <- var_gene[-mt_gene_idx] 
mean_gene_nomt  <- mean_gene[-mt_gene_idx]
zeroPercent_nomt    <- zeroPercent[-mt_gene_idx]
phi2 = coef(nls(var_gene_nomt ~ mean_gene_nomt + phi * mean_gene_nomt^2, start = list(phi = 0.1)))
px2 <- seq(min(mean_gene_nomt),max(mean_gene_nomt),length.out=100)
py2 <- px2 + phi2*px2^2
phi2_inv <- 1/phi2
expZero_nomt <- (phi2_inv/(phi2_inv+mean_gene_nomt))^phi2_inv



df_nomt  <- cbind.data.frame(meanx=mean_gene_nomt,varx=var_gene_nomt,
                            zeroPercent=zeroPercent_nomt,expZ=expZero_nomt)
df_nomt2 <- cbind.data.frame(px=px2,py=py2)


p1  <- ggplot() +
        geom_point(data=df, aes(x=meanx, y=varx))+
        geom_line(data=df2, aes(x=px, y=py,linetype="solid"),color="red")+
        geom_abline(aes(intercept=0,slope=1,linetype="dashed"),color="green")+
        labs(title="With Mitochondrial Genes",
            x="Mean", y = "Variance")+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position = "none")


# log
p2  <- ggplot() +
        geom_point(data=df, aes(x=log(meanx), y=log(varx)))+
        geom_line(data=df2, aes(x=log(px), y=log(py),linetype="solid"),color="red")+
        geom_abline(aes(intercept=0,slope=1,linetype="dashed"),color="green")+
        labs(title="With Mitochondrial Genes",
            x="Mean (log)", y = "Variance (log)")+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position = "none")




p3  <- ggplot() +
        geom_point(data=df_nomt, aes(x=meanx, y=varx))+
        geom_line(data=df_nomt2, aes(x=px, y=py,linetype="solid"),color="red")+
        geom_abline(aes(intercept=0,slope=1,linetype="dashed"),color="green")+
        labs(title="Without Mitochondrial Genes",
            x="Mean", y = "Variance")+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position = "none")


p4  <- ggplot() +
        geom_point(data=df_nomt, aes(x=log(meanx), y=log(varx)))+
        geom_line(data=df_nomt2, aes(x=log(px), y=log(py),linetype="solid"),color="red")+
        geom_abline(aes(intercept=0,slope=1,linetype="dashed"),color="green")+
        labs(title="Without Mitochondrial Genes",
            x="Mean (log)", y = "Variance (log)")+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position = "none")


library(ggpubr)
# var-mean figures
# var_mean_fig  <- ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
var_mean_fig    <- ggarrange(p1,p3,p2,p4,ncol=2,nrow=2)



library(tidyr)
df_long <- gather(df, condition, measurement, zeroPercent:expZ, factor_key=TRUE)
df_nomt_long <- gather(df_nomt, condition, measurement, zeroPercent:expZ, factor_key=TRUE)

p5  <- ggplot(data=df_long, aes(x=meanx, y=measurement,color=condition)) +
        geom_point()+
        scale_color_manual(values=c("black", "grey"), 
                       name="Genes",
                       labels=c("Observed", "Expected"))+
        labs(title="With Mitochondrial Genes",
            x="Mean", y = "Fraction Zeros")+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position=c(0.8,0.8),
                legend.title=element_text(size=18), 
                legend.text=element_text(size=15))


p6  <- ggplot(data=df_nomt_long, aes(x=meanx, y=measurement,color=condition)) +
        geom_point()+
        scale_color_manual(values=c("black", "grey"), 
                       name="Genes",
                       labels=c("Observed", "Expected"))+
        labs(title="Without Mitochondrial Genes",
            x="Mean", y = "Fraction Zeros")+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position="none")

zero_mean_fig   <- ggarrange(p5,p6,ncol=2,nrow=1)


cat("phi estimate with Mitochondrial:",phi1)
cat("phi estimate without Mitochondrial:",phi2)


#+ label=fig1,fig.width=12, fig.height=12, echo=F,fig.align="center"
var_mean_fig


#'***
#' > Comparing observed with expected zeros 
#+ label=fig2,fig.width=12, fig.height=6, echo=F,fig.align="center"
zero_mean_fig



#+ label=savenpz, echo=F, warnings=F, message=F,eval=F
# save as npz and fit the nature biotech paper, 
# not yet know how to run that gene wise thing
rm(list=ls())
datapath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
npzpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/data/npz/"
library(reticulate) 
library(Matrix)

np <- import("numpy")
sample      = "CN24_D1"
load(paste0(datapath,"/data/hdst/",sample,"_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
gn <- rownames(sp_count)
np$savez(paste0(npzpath,"/clean_",sample,"_idx_HDST_sparseCount.npz"), sp_count=sp_count,info=location,gn=gn)


