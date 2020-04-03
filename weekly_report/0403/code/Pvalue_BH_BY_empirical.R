#' ---
#' title: "BH & BY & Empirical"
#' author: "Jiaqiang Zhu"
#' date: "April 3rd, 2020"
#' ---


#+ label=empirical, echo=F, warnings=F, message=F,eval=F
rm(list=setdiff(ls(),"res"))
source("/net/fantasia/home/jiaqiang/mulan_temp/shortcut_function/fdr.cal.R")
getpval <- function(object){
    mtgene <- which(substr(rownames(object$res_mtest),start=1,stop=3)=="mt-")
    if(length(mtgene)>0){
        pvals <- object$res_mtest$combinedPval[-mtgene]
    } 
    names(pvals) <- rownames(object$res_mtest)[-mtgene]
    return(pvals)
}

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "CN24_D1"

out <- c()
for(iblock in c(50,25,10,5,1)){
    if(iblock !=1){
        load(paste0(workdir,"/output/kdc/clean_CN24_D1_idx_",iblock,"X_mean_Joint_kdc.rds"))
    }else{
        load(paste0(workdir,"/output/kdc/clean_CN24_D1_idx_",iblock,"X_Joint_kdc.rds"))
    }

    realpval    <- getpval(KDC)
    rm(KDC)
    perm_pval   <- c()

    for(iper in 1:10){
        if(iblock !=1){
            load(paste0(workdir,"/output/kdc/clean_CN24_D1_idx_",iblock,"X_mean_Joint_kdc_perm",iper,".rds"))
        }else{
            load(paste0(workdir,"/output/kdc/clean_CN24_D1_idx_",iblock,"X_Joint_kdc_perm",iper,".rds"))
        }
        perm_pval <- cbind(perm_pval,getpval(KDC))
        rm(KDC)
    }

    # a <- fdr.cal(realpval,perm_pval,rn=names(realpval))
    fdr_res <- fdr.cal2(realpval,perm_pval,rn=names(realpval))
    # out   <- c(out,sum(fdr_res$fdr<0.05))

    save(fdr_res,file=paste0(workdir,"/result/clean_CN24_D1_idx_",iblock,"X_mean_kdc_fdr.rds"))
    cat(iblock,"done\n")
    rm(fdr_res,realpval,perm_pval)
}


#+ label=empirical_slideseq, echo=F, warnings=F, message=F,eval=F
rm(list=setdiff(ls(),"res"))
source("/net/fantasia/home/jiaqiang/mulan_temp/shortcut_function/fdr.cal.R")
getpval <- function(object){
    mtgene <- which(substr(rownames(object$res_mtest),start=1,stop=3)=="mt-")
    if(length(mtgene)>0){
        pvals <- object$res_mtest$combinedPval[-mtgene]
    } 
    names(pvals) <- rownames(object$res_mtest)[-mtgene]
    return(pvals)
}

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"

iblock = 1
load(paste0(workdir,"/output/kdc/clean_Puck_180430_6_idx_",iblock,"X_Joint_kdc.rds"))

realpval    <- getpval(KDC)
rm(KDC)
perm_pval   <- c()
for(iper in 1:10){
    load(paste0(workdir,"/output/kdc/clean_Puck_180430_6_idx_",iblock,"X_Joint_kdc_perm",iper,".rds"))
    perm_pval <- cbind(perm_pval,getpval(KDC))
    rm(KDC)
}
fdr_res <- fdr.cal2(realpval,perm_pval,rn=names(realpval))
save(fdr_res,file=paste0(workdir,"/result/clean_Puck_180430_6_idx_",iblock,"X_mean_kdc_fdr.rds"))



#' ***
#' **HDST**
#+ label=HDST, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "CN24_D1"
res <- c()
for(iblock in c(50,25,10,5,1)){
    if(iblock !=1){
        load(paste0(workdir,"/output/kdc/clean_CN24_D1_idx_",iblock,"X_mean_Joint_kdc.rds"))
    }else{
        load(paste0(workdir,"/output/kdc/clean_CN24_D1_idx_",iblock,"X_Joint_kdc.rds"))
    }
    mtgene <- which(substr(rownames(KDC$res_mtest),start=1,stop=3)=="mt-")
    if(length(mtgene)>0){
        pvals <- KDC$res_mtest$combinedPval[-mtgene]
    } 

    load(paste0(workdir,"/result/clean_CN24_D1_idx_",iblock,"X_mean_kdc_fdr.rds"))


    by <- sum(p.adjust(pvals,method="BY")<0.05)
    bh <- sum(p.adjust(pvals,method="BH")<0.05)
    em <- sum(fdr_res$fdr<0.05)
    res <- rbind(res,c(by,bh,em))

    # cat("------------------------------------\n")
    # cat(iblock,"X: numGene:",nrow(KDC$res_mtest),"\n")
    # cat("signal genes:",sum(KDC$res_mtest$adjustedPval<0.05),"\n")
    # cat("signal genes nomt:",sum(p.adjust(pvals,method="BY")<0.05),"\n")
    rm(pvals,by,bh,KDC,em)
}

rownames(res) <- paste0("X",c(50,25,10,5,1))
colnames(res) <- c("BY","BH","Empirical")
print(res)


#' ***
#' **Slideseq**
#+ label=Slideseq, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"
res <- c()
iblock = 1
if(iblock !=1){
    load(paste0(workdir,"/output/kdc/clean_",isample,"_idx_",iblock,"X_mean_Joint_kdc.rds"))
}else{
    load(paste0(workdir,"/output/kdc/clean_",isample,"_idx_",iblock,"X_Joint_kdc.rds"))
}
mtgene <- which(substr(rownames(KDC$res_mtest),start=1,stop=3)=="mt-")
if(length(mtgene)>0){
    pvals <- KDC$res_mtest$combinedPval[-mtgene]
} 

load(paste0(workdir,"/result/clean_",isample,"_idx_",iblock,"X_mean_kdc_fdr.rds"))
by <- sum(p.adjust(pvals,method="BY")<0.05)
bh <- sum(p.adjust(pvals,method="BH")<0.05)
em <- sum(fdr_res$fdr<0.05)
res <- rbind(res,c(by,bh,em))
colnames(res) <- c("BY","BH","Empirical")
rownames(res) <- "X1"
print(res)



#' ***
#' **Empirical FDR Calculation Function**
#+ label=fdrcal, echo=T, warnings=F, message=F,eval=T
rm(list=ls())
fdr.cal2 <- function(altpval,nullpval,rn=NULL){ 
    if(!is.null(rn)){
            order.rn <- rn[order(altpval)]
        }else{
            order.rn <- NULL
        }
    altpval.order   <- altpval[order(altpval)]
    numperm         <- ncol(nullpval)
    nullpval.order  <- sapply(1:numperm,function(x){nullpval[,x][order(nullpval[,x])]})
    out             <- c()
    for(i in 1:nrow(nullpval)){
        s   <- sum(nullpval.order[1:i,]<=altpval.order[i])
        e   <- s/(i*numperm)
        out <- c(out,e)
    }
    if(is.null(order.rn)){
        return(out)
    }else{
        res <- cbind.data.frame(rn=order.rn,fdr=out)
        return(res)
    } 
}













