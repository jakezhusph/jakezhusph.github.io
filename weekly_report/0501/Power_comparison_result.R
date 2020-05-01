#' ---
#' title: "Power Comparison"
#' author: "Jiaqiang Zhu"
#' date: "April 30th, 2020"
#' ---


#+ label=KDC, echo=F, warnings=F, message=F,eval=F
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sim_power_func <- function(pval,upfdr=0.2,numTrue=1000,totalGene=10000){
        out         <- cbind(pval,c(rep(1,numTrue),rep(0,totalGene-numTrue)))
        order_out   <- out[order(out[,1]),]
        all_fdr     <- sapply(1:totalGene,function(x){sum(order_out[1:x,2]==0)/x})
        power       <- sapply(seq(0,upfdr,by=0.01),function(y){sum(order_out[,2][which(all_fdr<=y)])/numTrue})
        output      <- cbind.data.frame(fdr=seq(0,upfdr,by=0.01),power=power)
        return(output)
}

outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/kdc/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"


pattern_list <- c("streak","hotspot")

for(cellProp in c(0.2,0.3)){
	for(ieffect in c(2,3)){
		for(itheta in c(0.2,0.4)){
			for(mu0 in c(0.005)){
				for(isample in c(10000,30000)){
					cat("SS",isample,", mu:",mu0,"\n")
					all_power <- list()
					for(ipat in 1:2){
						single_power <- c()
						for(irpt in 1:10){
							allpval <- c()
							load(paste0(outpath,"/sim_kdc_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))

							single_power <- cbind(single_power,sim_power_func(KDC$res_mtest$combinedPval,upfdr=0.1,numTrue=100,totalGene=1000)$power)
							rm(allpval)
						}
						all_power[[ipat]] <- single_power
					}


					kdc_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",pattern_list))

					save(kdc_power,file=paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
					print(kdc_power)
					cat("----------------------------------------------------------------\n")
				    cat("----------------------------------------------------------------\n")
				    cat("\n")
				}
			}
		}
	}
}


#+ label=lmm_vst, echo=F, warnings=F, message=F,eval=F
rm(list=setdiff(ls(),"sim_power_func"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

# cellProp = 0.2
# ieffect = 1.5
# itheta = 5

pattern_list <- c("streak","hotspot")
for(cellProp in c(0.2,0.3)){
	cat("cellProp",cellProp,"\n")
	cat("--------------------------------\n")
	for(ieffect in c(2,3)){
		cat("ieffect",ieffect,"\n")
		cat("-------------\n")
		for(itheta in c(0.2,0.4)){
			cat("itheta",itheta,"\n")
			cat("-----\n")
			for(mu0 in c(0.005)){
				for(isample in c(10000)){
					cat("SS",isample,", mu:",mu0,"\n")
					all_power <- list()
					for(ipat in 1:2){
						single_power <- c()
						for(irpt in 1:10){
							allpval <- c()
							if(file.exists(paste0(outpath,"/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_rpt",irpt,".rds"))){
								load(paste0(outpath,"/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_rpt",irpt,".rds"))
		
								single_power <- cbind(single_power,sim_power_func(CombinePValues(SW_pval),upfdr=0.1,numTrue=100,totalGene=1000)$power)
								rm(allpval)
							}

						}
						all_power[[ipat]] <- single_power
					}

					lmm_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",pattern_list))
		
					save(lmm_power,file=paste0(respath,"/sim_spark_lmm_SW_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
					print(lmm_power)
					cat("----------------------------------------------------------------\n")
				    cat("----------------------------------------------------------------\n")
				    cat("\n")
				}
			}
		}
	}
}





#+ label=spatialDE, echo=F, warnings=F, message=F,eval=F
rm(list=setdiff(ls(),"sim_power_func"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spatialde/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"


sim_power_func <- function(pval,upfdr=0.2,numTrue=1000,totalGene=10000){
        out         <- cbind(pval,c(rep(1,numTrue),rep(0,totalGene-numTrue)))
        order_out   <- out[order(out[,1]),]
        all_fdr     <- sapply(1:totalGene,function(x){sum(order_out[1:x,2]==0)/x})
        power       <- sapply(seq(0,upfdr,by=0.01),function(y){sum(order_out[,2][which(all_fdr<=y)])/numTrue})
        output      <- cbind.data.frame(fdr=seq(0,upfdr,by=0.01),power=power)
        return(output)
}

# cellProp = 0.2
# ieffect = 2
# itheta = 5
pattern_list <- c("streak","hotspot")

for(cellProp in c(0.2,0.3)){
	cat("cellProp",cellProp,"\n")
	cat("--------------------------------\n")
	for(ieffect in c(2,3)){
		cat("ieffect",ieffect,"\n")
		cat("-------------\n")
		for(itheta in c(0.2,0.4)){
			cat("itheta",itheta,"\n")
			cat("-----\n")
			for(mu0 in c(0.005)){
				for(isample in c(10000)){
					cat("SS",isample,", mu:",mu0,"\n")
					all_power <- list()
					for(ipat in 1:2){
						single_power <- c()
						for(irpt in 1:10){
							allpval <- c()

							# res <- read.table(paste0(outpath,"/sim_spe_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
							res <- read.table(paste0(outpath,"/sim_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
							res1 <- res[-which(res$g=="log_total_count"),]

							if(sum(duplicated(res1$g))>0){
								clean_res <- res1[-which(duplicated(res1$g)),]
							}else{
								clean_res <- res1
							}

							order_res <- clean_res[order(as.numeric(sapply(strsplit(as.character(clean_res$g),split="gene"),"[[",2))),]
							if(nrow(clean_res)<1000){
								if(nrow(clean_res)==0){
									single_power <- cbind(single_power,rep(0,11))
								}else{
									numsig <- length(intersect(as.numeric(sapply(strsplit(as.character(order_res$g),split="gene"),"[[",2)),1:100))
									single_power <- cbind(single_power,sim_power_func(order_res$pval,upfdr=0.1,numTrue=numsig,totalGene=nrow(clean_res))$power)
								}
							}else{
								single_power <- cbind(single_power,sim_power_func(order_res$pval,upfdr=0.1,numTrue=100,totalGene=1000)$power)
						
							}
							rm(res,res1,order_res,clean_res)
						}

						all_power[[ipat]] <- single_power
					}

					spe_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",pattern_list))
					save(spe_power,file=paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
					print(spe_power)
					cat("----------------------------------------------------------------\n")
					cat("----------------------------------------------------------------\n")
				}
			}
		}
	}
}


