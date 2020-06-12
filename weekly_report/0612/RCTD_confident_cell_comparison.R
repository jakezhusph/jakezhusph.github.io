#' ---
#' title: "RCTD Confident Cell Comparison Across Methods"
#' author: "Jiaqiang Zhu"
#' date: "June 12th, 2020"
#' ---


#' **7177 Cells, 3240 Genes, Permutation Ten Times**  
#+ label=permutation, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/shortcut_function/fdr.cal.R")

ACAT <- function(Pvals,Weights=NULL){
    # check if there is NA
    if (sum(is.na(Pvals))>0){
        stop("Cannot have NAs in the p-values!")
    }
    # check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
        stop("P-values must be between 0 and 1!")
    }
    # check if there are pvals that are either exactly 0 or 1.
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

    # Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(Weights)){
        Weights<-rep(1/length(Pvals),length(Pvals))
    }else if (length(Weights)!=length(Pvals)){
        stop("The length of weights should be the same as that of the p-values")
    }else if (sum(Weights<0)>0){
        stop("All the weights must be positive!")
    }else{
        Weights<-Weights/sum(Weights)
    }


    # check if there are very small non-zero p values
    is.small<-(Pvals<1e-16)
    if (sum(is.small)==0){
        cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
    }else{
        cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
        cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
    }
    # check if the test statistic is very large.
    if (cct.stat>1e+15){
        pval<-(1/cct.stat)/pi
    }else{
        pval<-1-pcauchy(cct.stat)
    }
    return(pval)
}



workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"


vst_per <-  c()
for(iper in 1:10){
    tmp_vst <- c()
    for(iker in 1:10){
        load(paste0(workdir,"/spark/each_kernel/RCTD_spark_lmm_partgene_kernel",iker,"_per",iper,".rds"))
        tmp_vst <- cbind(tmp_vst,out2)
        rm(out2)
    }
    vst_per <- cbind(vst_per,apply(tmp_vst,1,ACAT))
}



spe_per <-  c()
for(iper in 1:10){
    raw_spe <- read.table(paste0(workdir,"/spatialde/puck_spe_partgene_per",iper,".csv"),head=T)
    tmp_spe <- raw_spe[-which(raw_spe$g=="log_total_count"),]
    if(sum(duplicated(tmp_spe$g))>0){
        spe_res <- tmp_spe[-which(duplicated(tmp_spe$g)),]
    }else{
        spe_res <- tmp_spe
    }
    suppressWarnings(spe_per <- cbind(spe_per,spe_res$pval))
    rm(spe_res)
}


kdc_per_nocov <-  c()
for(iper in 1:10){
    load(paste0(workdir,"/kdc/kdc_RCTD_partgenes_decomposed_per",iper,".rds"))
    kdc_per_nocov <- cbind(kdc_per_nocov,apply(KDC_nocov$res_stest,1,ACAT))
    rm(KDC,KDC_nocov)
}




source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")


col_base  <- c("#E64B35FF","#00A087FF","#3C5488FF","#F39B7FFF")

per_qq     <- qplot_gg(list(KDC=kdc_per_nocov,
                            LMM_SW_VST=vst_per,
                            SpatialDE=spe_per),
                            cl=0,
                            col.base=col_base,
                            legend.position=c(0.2,0.9),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("KDC","LMM_SW_VST","SpatialDE"),
                            self_label=c("SPARK-X","SPARK-G","SpatialDE"))



# kdc
load(paste0(workdir,"/kdc/kdc_RCTD_partgenes_decomposed.rds"))
p_kdc <- apply(KDC$res_stest,1,ACAT)
p_kdc_nocov <- apply(KDC_nocov$res_stest,1,ACAT)
col_base1  <- c("#E64B35FF","grey")
kdc_qq     <- qplot_gg(list(Real=p_kdc_nocov,
                            Permutation=kdc_per_nocov),
                            cl=0,
                            col.base=col_base1,
                            legend.position=c(0.2,0.9),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            self_label=c("Real","Permutation"))

# spark
res_vst <- c()
for(iker in 1:10){
    load(paste0(workdir,"/spark/each_kernel/RCTD_spark_lmm_partgene_kernel",iker,".rds"))
    res_vst <- cbind(res_vst,out2)
    rm(out2)
}

p_vst   <- apply(res_vst,1,ACAT)
col_base2  <- c("#00A087FF","grey")
spark_qq   <- qplot_gg(list(Real=p_vst,
                            Permutation=vst_per),
                            cl=0,
                            col.base=col_base2,
                            legend.position=c(0.2,0.9),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            self_label=c("Real","Permutation"))


# spatialde
raw_spe <- read.table(paste0(workdir,"/spatialde/puck_spe_partgene.csv"),head=T)
tmp_spe <- raw_spe[-which(raw_spe$g=="log_total_count"),]
if(sum(duplicated(tmp_spe$g))>0){
    spe_res <- tmp_spe[-which(duplicated(tmp_spe$g)),]
}else{
    spe_res <- tmp_spe
}
p_spe <- spe_res$pval

col_base3  <- c("#3C5488FF","grey")
spe_qq   <- qplot_gg(list(Real=p_spe,
                            Permutation=spe_per),
                            cl=0,
                            col.base=col_base3,
                            legend.position=c(0.2,0.9),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            self_label=c("Real","Permutation"))



library(ggpubr)
all_qq <- ggarrange(per_qq,kdc_qq,spark_qq,spe_qq,labels=c("A","B","C","D"),
                ncol = 4, nrow = 1,font.label=list(size=15))


#' *** 
#' > Figure 1: QQ plots  (A) Permutations  (B) SPARK-X  (C) SPARK-G  (D) SpatialDE
#+ label=per_fig,fig.width=20, fig.height=5, echo=F,fig.align="center",eval=T
all_qq



#' *** 
#' **Without Considering the Cell Type**  
#+ label=fdr_se, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"

load(paste0(workdir,"/spark/fdr_spark_lmm_RCTD_partgenes_decomposed.rds"))
spark_gene <- as.character(vst_fdr_res$rn)[vst_fdr_res$fdr<0.01]


load(paste0(workdir,"/kdc/fdr_kdc_RCTD_partgenes_decomposed.rds"))
kdc_gene <- as.character(kdc_fdr_res$rn)[kdc_fdr_res$fdr<0.01]

allgene <- union(kdc_gene,spark_gene)
se_pd <- cbind.data.frame(KDC=allgene%in%kdc_gene,
                      SPARK=allgene%in%spark_gene)
rownames(se_pd) <- allgene
library(eulerr)
VennDiag <- euler(se_pd,shape="ellipse") ## other wise some number if missing 
library(ggsci)
venn2 <- plot(VennDiag, quantities =T, lty=rep(0,ncol(se_pd)),fill=pal_npg("nrc")(ncol(se_pd)),alpha=1,legend=F, labels = c("SPARK-X", "SPARK-G"))


#' *** 
#' > Figure 2: Venn Diagram of SE Genes by Different Methods.
#+ label=fig2,fig.width=12, fig.height=8, echo=F,fig.align="center",eval=T
venn2



