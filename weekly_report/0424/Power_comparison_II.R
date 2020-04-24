#' ---
#' title: "Power Comparison for Simulation II (Gradient) "
#' author: "Jiaqiang Zhu"
#' date: "April 20th, 2020"
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


# ipat = 1
# isample = 300
# mu0 = 1.5
# cellProp = 0.2
# ieffect = 2
# itheta = 5

pattern_list <- c("gradient")

for(orderProb in c(0.3,0.5,0.7)){
	for(itheta in c(1)){
		for(mu0 in c(0.5,1.5)){
			for(isample in c(300,1000,3000)){
				cat("SS",isample,", mu:",mu0,"\n")
				all_power <- list()
				for(ipat in 1:1){
					single_power <- c()
					for(irpt in 1:10){
						allpval <- c()
						load(paste0(outpath,"/sim_kdc_",pattern_list[ipat],"_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))

						single_power <- cbind(single_power,sim_power_func(KDC$res_mtest$combinedPval,upfdr=0.1,numTrue=100,totalGene=1000)$power)
						rm(allpval)
					}
					all_power[[ipat]] <- single_power
				}


				kdc_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",pattern_list))

				save(kdc_power,file=paste0(respath,"/sim_kdc_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
				print(kdc_power)
				cat("----------------------------------------------------------------\n")
			    cat("----------------------------------------------------------------\n")
			    cat("\n")
			
		}
	}
}
}



#+ label=pmm, echo=F, warnings=F, message=F,eval=F
rm(list=setdiff(ls(),"sim_power_func"))
sim_power_func <- function(pval,upfdr=0.2,numTrue=1000,totalGene=10000){
        out         <- cbind(pval,c(rep(1,numTrue),rep(0,totalGene-numTrue)))
        order_out   <- out[order(out[,1]),]
        all_fdr     <- sapply(1:totalGene,function(x){sum(order_out[1:x,2]==0)/x})
        power       <- sapply(seq(0,upfdr,by=0.01),function(y){sum(order_out[,2][which(all_fdr<=y)])/numTrue})
        output      <- cbind.data.frame(fdr=seq(0,upfdr,by=0.01),power=power)
        return(output)
}
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

# cellProp = 0.2
# ieffect = 2
# itheta = 10

pattern_list <- c("gradient")
for(orderProb in c(0.3,0.5,0.7)){
	for(itheta in c(1)){
		for(mu0 in c(0.5,1.5)){
			for(isample in c(300,1000,3000)){
				cat("SS",isample,", mu:",mu0,"\n")
				all_power <- list()
				for(ipat in 1:1){
					single_power <- c()
					for(irpt in 1:10){
						allpval <- c()
						load(paste0(outpath,"/sim_spark_pmm_",pattern_list[ipat],"_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))

						single_power <- cbind(single_power,sim_power_func(spark_pval$combined_pvalue,upfdr=0.1,numTrue=100,totalGene=1000)$power)
						rm(allpval)
					}
					all_power[[ipat]] <- single_power
				}

				pmm_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",pattern_list))
				save(pmm_power,file=paste0(respath,"/sim_spark_pmm_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
				print(pmm_power)
				cat("----------------------------------------------------------------\n")
			    cat("----------------------------------------------------------------\n")
			    cat("\n")
			}
		}
	}
}



#+ label=lmm, echo=F, warnings=F, message=F,eval=F
rm(list=setdiff(ls(),"sim_power_func"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

# cellProp = 0.2
# ieffect = 1.5
# itheta = 5

pattern_list <- c("gradient")
for(orderProb in c(0.3,0.5,0.7)){
	for(itheta in c(5,10)){
		for(mu0 in c(0.5,1.5)){
			for(isample in c(300,1000,3000)){
				cat("SS",isample,", mu:",mu0,"\n")
				all_power <- list()
				for(ipat in 1:1){
					single_power <- c()
					for(irpt in 1:10){
						allpval <- c()
						if(file.exists(paste0(outpath,"/sim_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_qn_rpt",irpt,".rds"))){
							load(paste0(outpath,"/sim_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_qn_rpt",irpt,".rds"))

							# single_power <- cbind(single_power,sim_power_func(CombinePValues(RL_pval),upfdr=0.1,numTrue=100,totalGene=1000)$power)

							single_power <- cbind(single_power,sim_power_func(CombinePValues(SW_pval),upfdr=0.1,numTrue=100,totalGene=1000)$power)
							rm(allpval)
						}
						
					}
					all_power[[ipat]] <- single_power
				}


				lmm_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",pattern_list))

				# save(lmm_power,file=paste0(respath,"/sim_spark_lmm_RL_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_qn_power.rds"))
				save(lmm_power,file=paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_qn_power.rds"))
				print(lmm_power)
				cat("----------------------------------------------------------------\n")
			    cat("----------------------------------------------------------------\n")
			    cat("\n")
			}
		}
	}
}



#+ label=lmm_vst, echo=F, warnings=F, message=F,eval=F
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

# cellProp = 0.2
# ieffect = 1.5
# itheta = 5

pattern_list <- c("gradient")
for(orderProb in c(0.3,0.5,0.7)){
	for(itheta in c(5,10)){
		for(mu0 in c(0.5,1.5)){
			for(isample in c(300,1000,3000)){
				cat("SS",isample,", mu:",mu0,"\n")
				all_power <- list()
				for(ipat in 1:1){
					single_power <- c()
					for(irpt in 1:10){
						allpval <- c()
						if(file.exists(paste0(outpath,"/sim_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_rpt",irpt,".rds"))){
							load(paste0(outpath,"/sim_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_rpt",irpt,".rds"))

							# single_power <- cbind(single_power,sim_power_func(CombinePValues(RL_pval),upfdr=0.1,numTrue=100,totalGene=1000)$power)

							single_power <- cbind(single_power,sim_power_func(CombinePValues(SW_pval),upfdr=0.1,numTrue=100,totalGene=1000)$power)
							rm(allpval)
						}
		
					}
					all_power[[ipat]] <- single_power
				}

				lmm_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",pattern_list))

				# save(lmm_power,file=paste0(respath,"/sim_spark_lmm_RL_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
				save(lmm_power,file=paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
				print(lmm_power)
				cat("----------------------------------------------------------------\n")
			    cat("----------------------------------------------------------------\n")
			    cat("\n")
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

pattern_list <- c("gradient")
for(orderProb in c(0.3,0.5,0.7)){
	for(itheta in c(5,10)){
		for(mu0 in c(0.5,1.5)){
			for(isample in c(300,1000,3000)){
				cat("SS",isample,", mu:",mu0,"\n")
				all_power <- list()
				for(ipat in 1:1){
					single_power <- c()
					for(irpt in 1:10){
						allpval <- c()

						res <- read.table(paste0(outpath,"/sim_",pattern_list[ipat],"_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".csv"),head=T)

						res1 <- res[-which(res$g=="log_total_count"),]

						if(sum(duplicated(res1$g))>0){
							clean_res <- res1[-which(duplicated(res1$g)),]
						}else{
							clean_res <- res1
						}

						order_res <- clean_res[order(as.numeric(sapply(strsplit(as.character(clean_res$g),split="gene"),"[[",2))),]

						single_power <- cbind(single_power,sim_power_func(order_res$pval,upfdr=0.1,numTrue=100,totalGene=1000)$power)
						rm(res,res1,order_res,clean_res)
					}

					all_power[[ipat]] <- single_power
				}

				spe_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",pattern_list))
				save(spe_power,file=paste0(respath,"/sim_spe_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
				print(spe_power)
				cat("----------------------------------------------------------------\n")
				cat("----------------------------------------------------------------\n")
			}
		}
	}
}





#' **Power Figures (Fig.1-2): Sample Size Effect**  
#' orderProb = 0.2, phi = 0.1. SampleSize varies at 300,1000 and 3000. 
#+ label=powerfig1, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")
orderProb = 0.3
itheta = 10
mu0 = 1.5
for(isample in c(300,1000,3000)){
	icount = icount + 1
	load(paste0(respath,"/sim_spark_pmm_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
	load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
	lmm_vst_power = lmm_power
	rm(lmm_power)
	load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_qn_power.rds"))
	lmm_qn_power = lmm_power
	rm(lmm_power)
	load(paste0(respath,"/sim_kdc_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
	load(paste0(respath,"/sim_spe_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))

	single_sample_power_fig <- list()

	for(ipat in 1:1){
		pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
		if(ipat==1){
			single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.75,0.25),
								len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
		}else{
			single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
								len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
		}

		rm(pltdat)
	}
	
	power_fig[[icount]] <- single_sample_power_fig
}

library(ggpubr)
fig1 <- ggarrange(power_fig[[1]][[1]], power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)



#' ***
#' > Figure 1: Gradient (mu = 1.5)
#+ label=fig1,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig1



#+ label=powerfig2, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")
orderProb = 0.3
itheta = 10
mu0 = 0.5
for(isample in c(300,1000,3000)){
	icount = icount + 1
	load(paste0(respath,"/sim_spark_pmm_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
	load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
	lmm_vst_power = lmm_power
	rm(lmm_power)
	load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_qn_power.rds"))
	lmm_qn_power = lmm_power
	rm(lmm_power)
	load(paste0(respath,"/sim_kdc_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
	load(paste0(respath,"/sim_spe_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))

	single_sample_power_fig <- list()

	for(ipat in 1:1){
		pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
		if(ipat==1){
			single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.8,0.15),
								len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
		}else{
			single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
								len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
		}

		rm(pltdat)
	}
	
	power_fig[[icount]] <- single_sample_power_fig
}

library(ggpubr)
fig2 <- ggarrange(power_fig[[1]][[1]]+rremove("legend"), power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)


#' ***
#' > Figure 2: Gradient (mu = 0.5)

#+ label=fig2,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig2



#' **Power Figure : Fraction of Marked Cells Effect**  
#' SS = 300, mu = 0.5, phi = 0.1. Fraction of Marked Cells Varies at 0.3, 0.5 and 0.7

#+ label=powerfig3, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

itheta = 10
mu0 = 0.5
for(isample in c(300)){
	for(orderProb in c(0.3,0.5,0.7)){
		icount = icount + 1
		load(paste0(respath,"/sim_spark_pmm_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
		load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
		lmm_vst_power = lmm_power
		rm(lmm_power)
		load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_qn_power.rds"))
		lmm_qn_power = lmm_power
		rm(lmm_power)
		load(paste0(respath,"/sim_kdc_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
		load(paste0(respath,"/sim_spe_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))

		single_sample_power_fig <- list()

		for(ipat in 1:1){
			pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
			if(ipat==1){
				single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.8,0.15),
									len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
			                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
			}else{
				single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
									len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
			                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
			}

			rm(pltdat)
		}
		
		power_fig[[icount]] <- single_sample_power_fig
	}
}



library(ggpubr)
fig3 <- ggarrange(power_fig[[1]][[1]]+rremove("legend"), power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)


#' ***
#' > Figure 3: Gradient (mu=0.5)

#+ label=fig3,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig3





#' **Power Figure : Dispersion Effect**  
#' SS = 300, mu = 0.5, orderProb=0.3. Phi Varies at 0.1, 0.2 and 1.

#+ label=powerfig4, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

mu0 = 0.5
orderProb = 0.3
for(isample in c(300)){
	for(itheta in c(10,5,1)){
		icount = icount + 1
		load(paste0(respath,"/sim_spark_pmm_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
		load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
		lmm_vst_power = lmm_power
		rm(lmm_power)
		load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_qn_power.rds"))
		lmm_qn_power = lmm_power
		rm(lmm_power)
		load(paste0(respath,"/sim_kdc_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
		load(paste0(respath,"/sim_spe_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))

		single_sample_power_fig <- list()

		for(ipat in 1:1){
			pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
			if(ipat==1){
				single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.8,0.15),
									len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
			                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
			}else{
				single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
									len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
			                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
			}

			rm(pltdat)
		}
		
		power_fig[[icount]] <- single_sample_power_fig
	}
}



library(ggpubr)
fig4 <- ggarrange(power_fig[[1]][[1]]+rremove("legend"), power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"), labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)


#' ***
#' > Figure 4: Gradient 

#+ label=fig4,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig4



