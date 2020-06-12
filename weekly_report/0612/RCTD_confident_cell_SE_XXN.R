#' ---
#' title: "RCTD Confident Cell Analysis III: XXN"
#' author: "Jiaqiang Zhu"
#' date: "June 12th, 2020"
#' ---


#' **A dataset of 7177 cells and 3240 genes is analyzed **  
#+ label=rawcount, echo=F, warnings=F, message=F,eval=T

rm(list=ls())
library(RCTD)
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)
# source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
# sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")
workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"
load(paste0(workdir,"/kdc/fdr_kdc_RCTD_partgenes_decomposed.rds"))
load(paste0(workdir,"/puck_4kdc.rds"))

kdc_genes       <- as.character(kdc_fdr_res$rn[kdc_fdr_res$fdr<0.01])
kdc_nocov_genes <- as.character(kdc_fdr_res_nocov$rn[kdc_fdr_res_nocov$fdr<0.01])

both_genes      <- intersect(kdc_genes,kdc_nocov_genes)
no_cell_genes   <- setdiff(kdc_nocov_genes,kdc_genes)

figpath = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/figure/"

source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

ct_rel  <- apply(sp_count[both_genes,],1,relative_func)
pltdat  <- cbind.data.frame(location[,1:2],ct_rel)

colnames(pltdat) <- c("x","y",both_genes)

pp     <- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc2(pltdat,x,main=T,
                                pointsize=0.5,
                                titlesize=1.5,
                                min.pand=0.8,
                                max.pand=1.01,opt="D",
                                direct=-1,legend.position="none",
                                barheight=10 )})

library(ggpubr) 
ifig = 1
fig1 <- ggarrange(pp[[9*(ifig-1)+1]],pp[[9*(ifig-1)+2]],pp[[9*(ifig-1)+3]],
                    pp[[9*(ifig-1)+4]],pp[[9*(ifig-1)+5]],pp[[9*(ifig-1)+6]],
                    pp[[9*(ifig-1)+7]],pp[[9*(ifig-1)+8]],pp[[9*(ifig-1)+9]],
                    pp[[9*(ifig-1)+10]],
                    ncol = 5, nrow = 2)

#' *** 
#' > Figure 1: Expression Pattern for SE Genes Identified in both scenarios (Count)
#+ label=fig1,fig.width=20, fig.height=8, echo=F,fig.align="center",eval=T
fig1



#+ label=residuals, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(RCTD)
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)
# source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
# sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")
workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"
load(paste0(workdir,"/kdc/fdr_kdc_RCTD_partgenes_decomposed.rds"))
load(paste0(workdir,"/puck_4kdc.rds"))


RCTD_object <- readRDS(paste0(workdir,"/SplitPuckResults/results_decomposed.RDS"))
cate_x <- c()
for(itype in RCTD_object@cell_type_names){
    cate_x <- cbind(cate_x,as.numeric(RCTD_object@cell_labels==itype))
}
colnames(cate_x) <- RCTD_object@cell_type_names

kdc_genes       <- as.character(kdc_fdr_res$rn[kdc_fdr_res$fdr<0.01])
kdc_nocov_genes <- as.character(kdc_fdr_res_nocov$rn[kdc_fdr_res_nocov$fdr<0.01])

both_genes      <- intersect(kdc_genes,kdc_nocov_genes)
no_cell_genes   <- setdiff(kdc_nocov_genes,kdc_genes)

figpath = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/figure/"

source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")


# XTX_inv         <- solve(crossprod(cate_x,cate_x))
selected_sp_count <- sp_count[both_genes,]-sp_count[both_genes,]%*%cate_x%*%t(cate_x)/ncol(sp_count)



ct_rel  <- apply(selected_sp_count,1,relative_func)
pltdat  <- cbind.data.frame(location[,1:2],ct_rel)


colnames(pltdat) <- c("x","y",both_genes)

pp     <- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc2(pltdat,x,main=T,
                                pointsize=0.5,
                                titlesize=1.5,
                                min.pand=0.8,
                                max.pand=1.01,opt="D",
                                direct=-1,legend.position="none",
                                barheight=10 )})
library(ggpubr) 
ifig = 1
fig2 <- ggarrange(pp[[9*(ifig-1)+1]],pp[[9*(ifig-1)+2]],pp[[9*(ifig-1)+3]],
                    pp[[9*(ifig-1)+4]],pp[[9*(ifig-1)+5]],pp[[9*(ifig-1)+6]],
                    pp[[9*(ifig-1)+7]],pp[[9*(ifig-1)+8]],pp[[9*(ifig-1)+9]],
                    pp[[9*(ifig-1)+10]],
                    ncol = 5, nrow = 2)

#' *** 
#' > Figure 2: Expression Pattern for SE Genes Identified in both scenarios (Residuals)
#+ label=fig2,fig.width=20, fig.height=8, echo=F,fig.align="center",eval=T
fig2





#+ label=rawcount2, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(RCTD)
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)
workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"
load(paste0(workdir,"/kdc/fdr_kdc_RCTD_partgenes_decomposed.rds"))
load(paste0(workdir,"/puck_4kdc.rds"))

kdc_genes       <- as.character(kdc_fdr_res$rn[kdc_fdr_res$fdr<0.01])
kdc_nocov_genes <- as.character(kdc_fdr_res_nocov$rn[kdc_fdr_res_nocov$fdr<0.01])

both_genes      <- intersect(kdc_genes,kdc_nocov_genes)
no_cell_genes   <- setdiff(kdc_nocov_genes,kdc_genes)

figpath = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/figure/"

source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

ct_rel  <- apply(sp_count[no_cell_genes,],1,relative_func)
pltdat  <- cbind.data.frame(location[,1:2],ct_rel)

colnames(pltdat) <- c("x","y",no_cell_genes)

pp     <- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc2(pltdat,x,main=T,
                                pointsize=0.5,
                                titlesize=1.5,
                                min.pand=0.8,
                                max.pand=1.01,opt="D",
                                direct=-1,legend.position="none",
                                barheight=10 )})

library(ggpubr) 
ifig = 1
fig3 <- ggarrange(pp[[9*(ifig-1)+1]],pp[[9*(ifig-1)+2]],pp[[9*(ifig-1)+3]],
                    pp[[9*(ifig-1)+4]],pp[[9*(ifig-1)+5]],pp[[9*(ifig-1)+6]],
                    pp[[9*(ifig-1)+7]],pp[[9*(ifig-1)+8]],pp[[9*(ifig-1)+9]],
                    pp[[9*(ifig-1)+10]],
                    ncol = 5, nrow = 2)

#' *** 
#' > Figure 3: Expression Pattern for SE Genes Identified Only in NCT (Count)
#+ label=fig3,fig.width=20, fig.height=8, echo=F,fig.align="center",eval=T
fig3



#+ label=residuals2, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(RCTD)
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)
# source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
# sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")
workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"
load(paste0(workdir,"/kdc/fdr_kdc_RCTD_partgenes_decomposed.rds"))
load(paste0(workdir,"/puck_4kdc.rds"))


RCTD_object <- readRDS(paste0(workdir,"/SplitPuckResults/results_decomposed.RDS"))
cate_x <- c()
for(itype in RCTD_object@cell_type_names){
    cate_x <- cbind(cate_x,as.numeric(RCTD_object@cell_labels==itype))
}
colnames(cate_x) <- RCTD_object@cell_type_names

kdc_genes       <- as.character(kdc_fdr_res$rn[kdc_fdr_res$fdr<0.01])
kdc_nocov_genes <- as.character(kdc_fdr_res_nocov$rn[kdc_fdr_res_nocov$fdr<0.01])

both_genes      <- intersect(kdc_genes,kdc_nocov_genes)
no_cell_genes   <- setdiff(kdc_nocov_genes,kdc_genes)

figpath = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/figure/"

source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")


# XTX_inv         <- solve(crossprod(cate_x,cate_x))
selected_sp_count <- sp_count[no_cell_genes,]-sp_count[no_cell_genes,]%*%cate_x%*%t(cate_x)/ncol(sp_count)


ct_rel  <- apply(selected_sp_count,1,relative_func)
pltdat  <- cbind.data.frame(location[,1:2],ct_rel)


colnames(pltdat) <- c("x","y",no_cell_genes)

pp     <- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc2(pltdat,x,main=T,
                                pointsize=0.5,
                                titlesize=1.5,
                                min.pand=0.8,
                                max.pand=1.01,opt="D",
                                direct=-1,legend.position="none",
                                barheight=10 )})
library(ggpubr) 
ifig = 1
fig4 <- ggarrange(pp[[9*(ifig-1)+1]],pp[[9*(ifig-1)+2]],pp[[9*(ifig-1)+3]],
                    pp[[9*(ifig-1)+4]],pp[[9*(ifig-1)+5]],pp[[9*(ifig-1)+6]],
                    pp[[9*(ifig-1)+7]],pp[[9*(ifig-1)+8]],pp[[9*(ifig-1)+9]],
                    pp[[9*(ifig-1)+10]],
                    ncol = 5, nrow = 2)

#' *** 
#' > Figure 4: Expression Pattern for SE Genes Identified Only in NCT (Residuals)
#+ label=fig4,fig.width=20, fig.height=8, echo=F,fig.align="center",eval=T
fig4






