#' ---
#' title: "Fit Negative Binomial Curve For Four Real Data"
#' author: "Jiaqiang Zhu"
#' date: "April 19th, 2020"
#' ---

#'***
#' > Mouse Olfactory Bulb ST Data
#+ label=mobst, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
library(fitdistrplus)
library(ggplot2)
library(Matrix)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")

datapath = "/net/mulan/disk2/jiaqiang/SDE/Realdata/rawdata/"
countdata   <- read.table(paste0(datapath,"/Rep11_MOB_count_matrix-1.tsv"),check.names=F)
rn          <- rownames(countdata)
info        <- cbind.data.frame(x=as.numeric(sapply(strsplit(rn,split="x"),"[",1)),y=as.numeric(sapply(strsplit(rn,split="x"),"[",2)))
rownames(info) <- rn

sp_count    <- as(t(countdata),"sparseMatrix")

numZero_gene    <- as.vector(ncol(sp_count)- sp_nz_count_Rcpp(sp_count,rowSums=T))
num_nonZero_gene <- as.vector(sp_nz_count_Rcpp(sp_count,rowSums=T))
nonZeroPercent  <- num_nonZero_gene/ncol(sp_count)
zeroPercent     <- 1-nonZeroPercent

mean_gene       <- as.vector(sp_means_Rcpp(sp_count,rowMeans=T))
var_gene        <- as.vector(sp_vars_Rcpp(sp_count,rowVars=T))
mt_gene_idx     <- which(substr(rownames(sp_count),start=1,stop=3)=="mt-")

phi1 = coef(nls(var_gene ~ mean_gene + phi * mean_gene^2, start = list(phi = 0.1)))
px1 <- seq(min(mean_gene),max(mean_gene),length.out=100)
py1 <- px1 + phi1*px1^2
phi1_inv <- 1/phi1
expZero <- (phi1_inv/(phi1_inv+mean_gene))^phi1_inv

df  <- cbind.data.frame(meanx=mean_gene,varx=var_gene,zeroPercent=zeroPercent,expZ=expZero)
df2 <- cbind.data.frame(px=px1,py=py1)

p1  <- ggplot() +
        geom_point(data=df, aes(x=log10(meanx), y=log10(varx)))+
        geom_line(data=df2, aes(x=log10(px), y=log10(py),linetype="solid"),color="red")+
        geom_abline(aes(intercept=0,slope=1,linetype="dashed"),color="green")+
        labs(title="Mouse Olfactory Bulb ST Data",
            x="Mean", y = "Variance")+
        scale_x_continuous(breaks=seq(-2,2,by=1),
                    labels=c(expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2)))+
        scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                    labels=c(expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2),
                            expression(10^3)))+
        # geom_text(aes(x=-2, y=3), label=expression(paste(phi,"=0.2164729")),parse = T,size=5)+
        geom_text(data = data.frame(),aes(x=-1, y=3,label="phi==0.2164729 *','~~ tilde(mu)==0.52672"), parse = T,size=5)+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position = "none")

cat("numGene:",nrow(sp_count),"; numSample:",ncol(sp_count)," phi:",phi1,"\n")
#+ echo=F, warnings=F, message=F,eval=T
print(summary(mean_gene))




#'***
#' > Breast Cancer ST Data
#+ label=BCst, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"p1"))
library(Rcpp)
library(fitdistrplus)
library(ggplot2)
library(Matrix)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")

datapath = "/net/mulan/disk2/jiaqiang/SDE/Realdata/rawdata/"
countdata   <- read.table(paste0(datapath,"/Layer2_BC_count_matrix-1.tsv"),check.names=F)
rn          <- rownames(countdata)
info        <- cbind.data.frame(x=as.numeric(sapply(strsplit(rn,split="x"),"[",1)),y=as.numeric(sapply(strsplit(rn,split="x"),"[",2)))
rownames(info) <- rn
sp_count    <- as(t(countdata),"sparseMatrix")

numZero_gene    <- as.vector(ncol(sp_count)- sp_nz_count_Rcpp(sp_count,rowSums=T))
num_nonZero_gene <- as.vector(sp_nz_count_Rcpp(sp_count,rowSums=T))
nonZeroPercent  <- num_nonZero_gene/ncol(sp_count)
zeroPercent     <- 1-nonZeroPercent

mean_gene       <- as.vector(sp_means_Rcpp(sp_count,rowMeans=T))
var_gene        <- as.vector(sp_vars_Rcpp(sp_count,rowVars=T))
mt_gene_idx     <- which(substr(rownames(sp_count),start=1,stop=3)=="mt-")


phi1 = coef(nls(var_gene ~ mean_gene + phi * mean_gene^2, start = list(phi = 0.1)))


px1 <- seq(min(mean_gene),max(mean_gene),length.out=100)
py1 <- px1 + phi1*px1^2
phi1_inv <- 1/phi1
expZero <- (phi1_inv/(phi1_inv+mean_gene))^phi1_inv

df  <- cbind.data.frame(meanx=mean_gene,varx=var_gene,zeroPercent=zeroPercent,expZ=expZero)
df2 <- cbind.data.frame(px=px1,py=py1)

p2  <- ggplot() +
        geom_point(data=df, aes(x=log10(meanx), y=log10(varx)))+
        geom_line(data=df2, aes(x=log10(px), y=log10(py),linetype="solid"),color="red")+
        geom_abline(aes(intercept=0,slope=1,linetype="dashed"),color="green")+
        labs(title="Breast Cancer ST Data",
            x="Mean", y = "Variance")+
        scale_x_continuous(breaks=seq(-2,2,by=1),
                    labels=c(expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2)))+
        scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                    labels=c(expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2),
                            expression(10^3)))+
        # geom_text(aes(x=-2, y=3), label=expression(paste(phi,"=0.2164729")),parse = T,size=5)+
        geom_text(data = data.frame(),aes(x=-1.5, y=3,label="phi==1.730921 *','~~ tilde(mu)==0.075697"), parse = T,size=5)+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position = "none")

cat("numGene:",nrow(sp_count),"; numSample:",ncol(sp_count)," phi:",phi1,"\n")
#+ echo=F, warnings=F, message=F,eval=T
print(summary(mean_gene))


#'***
#' > Mouse Olfactory HDST Data
#+ label=HDst, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("p1","p2")))
library(Rcpp)
library(fitdistrplus)
library(ggplot2)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")
datapath    = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
sample      = "CN24_D1"
load(paste0(datapath,"/data/hdst/",sample,"_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))

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


# log
p3  <- ggplot() +
        geom_point(data=df, aes(x=log10(meanx), y=log10(varx)))+
        geom_line(data=df2, aes(x=log10(px), y=log10(py),linetype="solid"),color="red")+
        geom_abline(aes(intercept=0,slope=1,linetype="dashed"),color="green")+
        labs(title="Mouse Olfactory HDST Data",
            x="Mean", y = "Variance")+
        scale_x_continuous(breaks=-6:2,
                    labels=c(expression(10^-6),
                            expression(10^-5),
                            expression(10^-4),
                            expression(10^-3),
                            expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2)))+
        scale_y_continuous(breaks=-6:3,
                    labels=c(expression(10^-6),
                            expression(10^-5),
                            expression(10^-4),
                            expression(10^-3),
                            expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2),
                            expression(10^3)))+
        # geom_text(aes(x=-2, y=3), label=expression(paste(phi,"=0.2164729")),parse = T,size=5)+
        geom_text(data = data.frame(),aes(x=-4, y=3,label="phi==1.009104 *','~~ tilde(mu)==0.0000662"), parse = T,size=5)+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position = "none")

cat("numGene:",nrow(sp_count),"; numSample:",ncol(sp_count)," phi:",phi1,"\n")
#+ echo=F, warnings=F, message=F,eval=T
print(summary(mean_gene))



#'***
#' > Slideseq Data
#+ label=slideseq, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("p1","p2","p3")))
library(Rcpp)
library(fitdistrplus)
library(ggplot2)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")

datapath    = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
load(paste0(datapath,"/data/slideseq/Puck_180430_6_zero_removed.rds"))

numZero_gene    <- as.vector(ncol(sp_count)- sp_nz_count_Rcpp(sp_count,rowSums=T))
num_nonZero_gene <- as.vector(sp_nz_count_Rcpp(sp_count,rowSums=T))
nonZeroPercent  <- num_nonZero_gene/ncol(sp_count)
zeroPercent     <- 1-nonZeroPercent

mean_gene       <- as.vector(sp_means_Rcpp(sp_count,rowMeans=T))
var_gene        <- as.vector(sp_vars_Rcpp(sp_count,rowVars=T))
mt_gene_idx     <- which(substr(rownames(sp_count),start=1,stop=3)=="mt-")

phi1 = coef(nls(var_gene ~ mean_gene + phi * mean_gene^2, start = list(phi = 0.1)))

px1 <- seq(min(mean_gene),max(mean_gene),length.out=100)
py1 <- px1 + phi1*px1^2
phi1_inv <- 1/phi1
expZero <- (phi1_inv/(phi1_inv+mean_gene))^phi1_inv

df  <- cbind.data.frame(meanx=mean_gene,varx=var_gene,zeroPercent=zeroPercent,expZ=expZero)
df2 <- cbind.data.frame(px=px1,py=py1)


# log
p4  <- ggplot() +
        geom_point(data=df, aes(x=log10(meanx), y=log10(varx)))+
        geom_line(data=df2, aes(x=log10(px), y=log10(py),linetype="solid"),color="red")+
        geom_abline(aes(intercept=0,slope=1,linetype="dashed"),color="green")+
        labs(title="Slideseq Data (Puck_180430_6)",
            x="Mean", y = "Variance")+
        scale_x_continuous(breaks=-6:2,
                    labels=c(expression(10^-6),
                            expression(10^-5),
                            expression(10^-4),
                            expression(10^-3),
                            expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2)))+
        scale_y_continuous(breaks=-6:3,
                    labels=c(expression(10^-6),
                            expression(10^-5),
                            expression(10^-4),
                            expression(10^-3),
                            expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2),
                            expression(10^3)))+
        # geom_text(aes(x=-2, y=3), label=expression(paste(phi,"=0.2164729")),parse = T,size=5)+
        geom_text(data = data.frame(),aes(x=-3, y=3,label="phi==3.545992 *','~~ tilde(mu)==0.001448"), parse = T,size=5)+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5,size=20),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=18),
                axis.title.x=element_text(size=18),
                legend.position = "none")


cat("numGene:",nrow(sp_count),"; numSample:",ncol(sp_count)," phi:",phi1,"\n")
#+ echo=F, warnings=F, message=F,eval=T
print(summary(mean_gene))


#'***
#' > Variance vs. Mean Figures
#+ label=fig1,fig.width=12, fig.height=12, echo=F,fig.align="center"
library(ggpubr)
var_mean_fig    <- ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
var_mean_fig
