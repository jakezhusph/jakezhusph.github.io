#' ---
#' title: "RCTD All Cell Analysis"
#' author: "Jiaqiang Zhu"
#' date: "June 12th, 2020"
#' ---


#' **22207 cells and 17754 genes**  
#+ label=QQ_comparison, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(RCTD)
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)

ACAT <- function(Pvals,Weights=NULL){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
        stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
        stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(sum(Pvals==0)>=1)
    is.one<-(sum(Pvals==1)>=1)
    if (is.zero && is.one){
        stop("Cannot have both 0 and 1 p-values!")
    }
    if (is.zero){
        return(0)
    }
    if (is.one){
        # warning("There are p-values that are exactly 1!")
        return(1)
    }

    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(Weights)){
        Weights<-rep(1/length(Pvals),length(Pvals))
    }else if (length(Weights)!=length(Pvals)){
        stop("The length of weights should be the same as that of the p-values")
    }else if (sum(Weights<0)>0){
        stop("All the weights must be positive!")
    }else{
        Weights<-Weights/sum(Weights)
    }


    #### check if there are very small non-zero p values
    is.small<-(Pvals<1e-16)
    if (sum(is.small)==0){
        cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
    }else{
        cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
        cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
    }
    #### check if the test statistic is very large.
    if (cct.stat>1e+15){
        pval<-(1/cct.stat)/pi
    }else{
        pval<-1-pcauchy(cct.stat)
    }
    return(pval)
}


workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"
load(paste0(workdir,"/kdc/fdr_kdc_RCTD_partgenes_decomposed_XTX.rds"))
load(paste0(workdir,"/puck_4kdc.rds"))


figpath = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/figure/"
source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")


workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"

res1 <- res2 <- c()
for(iper in 1:10){
    load(paste0(workdir,"/kdc/kdc_RCTD_allcells_cut20_decomposed_per",iper,"_XTX.rds"))
    res1 <- cbind(res1,KDC$res_mtest$combinedPval)
    res2 <- cbind(res2,KDC_nocov$res_mtest$combinedPval)
    rm(KDC,KDC_nocov)
}

res3 <- c()
for(iper in 1:10){
    load(paste0(workdir,"/kdc/kdc_RCTD_allcells_cut20_decomposed_per",iper,".rds"))
    res3 <- cbind(res3,KDC$res_mtest$combinedPval)
    rm(KDC,KDC_nocov)
}



source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

comp_qq     <- qplot_gg(list(Adjusted_Cell_Type_XTX=as.vector(res1),
                            Adjusted_Cell_Type_XXN=as.vector(res3),
                            No_Cell_Type=as.vector(res2)),
                            cl=0,
                            # col.base=col_base,
                            legend.position=c(0.35,0.9),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1,
                            factor_level=c("Adjusted_Cell_Type_XTX","Adjusted_Cell_Type_XXN","No_Cell_Type"),
                            self_label=c("Adjusted_Cell_Type_XTX","Adjusted_Cell_Type_XXN","No_Cell_Type"))                      



#' *** 
#' > Figure 1: QQ Plots of KDC. 
#+ label=fig1,fig.width=4, fig.height=4, echo=F,fig.align="center",eval=T
comp_qq





#+ label=Signal_comparison, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"ACAT"))
library(RCTD)
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)

workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"


load(paste0(workdir,"/kdc/fdr_kdc_RCTD_allcells_cut20_decomposed.rds"))

NX_celltype         <- as.character(kdc_fdr_res$rn[kdc_fdr_res$fdr<0.01])
NX_without_celltype <- as.character(kdc_fdr_res_nocov$rn[kdc_fdr_res_nocov$fdr<0.01])


rm(kdc_fdr_res,kdc_fdr_res_nocov)
load(paste0(workdir,"/kdc/fdr_kdc_RCTD_allcells_cut20_decomposed_XTX.rds"))

XTX_celltype         <- as.character(kdc_fdr_res$rn[kdc_fdr_res$fdr<0.01])
XTX_without_celltype <- as.character(kdc_fdr_res_nocov$rn[kdc_fdr_res_nocov$fdr<0.01])



NX_allgene <- union(NX_celltype,NX_without_celltype)
NX_pd <- cbind.data.frame(Adjusted_Cell_Type=NX_allgene%in%NX_celltype,
                      No_Cell_Type=NX_allgene%in%NX_without_celltype
                      )
rownames(NX_pd) <- NX_allgene
library(eulerr)
VennDiag <- euler(NX_pd,shape="ellipse") ## other wise some number if missing 
library(ggsci)
venn1 <- plot(VennDiag, quantities =T, lty=rep(0,ncol(NX_pd)),fill=pal_npg("nrc")(ncol(NX_pd)),alpha=1,legend=F, labels = c("ACT", "NCT"))




XTX_allgene <- union(XTX_celltype,XTX_without_celltype)
XTX_pd <- cbind.data.frame(Adjusted_Cell_Type=XTX_allgene%in%XTX_celltype,
                      No_Cell_Type=XTX_allgene%in%XTX_without_celltype
                      )
rownames(XTX_pd) <- XTX_allgene
library(eulerr)
VennDiag <- euler(XTX_pd,shape="ellipse") ## other wise some number if missing 
library(ggsci)
venn2 <- plot(VennDiag, quantities =T, lty=rep(0,ncol(XTX_pd)),fill=pal_npg("nrc")(ncol(XTX_pd)),alpha=1,legend=F, labels = c("ACT", "NCT"))




library(ggpubr)
fig2 <- ggarrange(venn1,venn2,labels=c("A","B"),
                ncol = 2, nrow = 1,font.label=list(size=15))

#' *** 
#' > Figure 2: Venn Diagram of SE Genes By KDC. **(A)** XX^T^/N ; **(B)** X(X^T^X)^-1^X^T^. ACT: Adjusted for Cell Type. NCT: No Cell Type Information Used
#+ label=fig2,fig.width=14, fig.height=6, echo=F,fig.align="center",eval=T
fig2



#' *** 
#' > Enrichment results based on the SE genes after adjusting for cell type
#+ label=Enrichment_comparison, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"

for(iter in c(".",".XTX.")){
    load(paste0(workdir,"/enrichment/res.RCTD.allcells.cut20.enrichGO.diff2all.kdc",iter,"rds"))

    load(paste0(workdir,"/enrichment/res.RCTD.allcells.cut20.enrichKEGG.diff2all.kdc",iter,"rds"))

    if(iter==".XTX."){
        cat("==========================================\n")
        cat("XTX:\n")
        cat("==========================================\n")
    }else{
        cat("==========================================\n")
        cat("XXN:\n")
        cat("==========================================\n")
    }
    cat("Go Term Enrichment:",sum(ego$p.adjust<0.05),"\n")
    print(ego[1:10,3:7])
    cat("KEGG Enrichment:",sum(kk.diff$p.adjust<0.05),"\n")
    print(kk.diff[1:10,2:6])
}


