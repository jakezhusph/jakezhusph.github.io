#' ---
#' title: "Power Comparison for Simulation Based on HDST (TotalCount Adjusted)"
#' author: "Jiaqiang Zhu"
#' date: "April 16th, 2020"
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

outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/result/"

numZero = 5000 
sigZero = 500
layerlist <- c("E","ONL","GL")
for(iblock in c(38,5,1)){
	if(iblock != 1){
		load(paste0("/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/sim_data_CN24_D1_",iblock,"X_ONL_zero",numZero,"_sigzero",sigZero,"_binR_rpt1.rds"))
	}else{
		load(paste0("/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/sim_data_CN24_D1_HDST_ONL_zero",numZero,"_sigzero",sigZero,"_rpt1.rds"))
	}
	
	cat(paste0(iblock,"X: number of Spots ",nrow(location),"\n"))
	all_power <- list()
	for(ilayer in 1:3){
		single_power <- c()
		for(isid in 1:10){
			allpval <- c()
			if(iblock!=1){
				load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))
			}else{
				load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_rpt",isid,".rds"))
			}
			single_power <- cbind(single_power,sim_power_func(KDC$res_mtest$combinedPval,upfdr=0.1,numTrue=100,totalGene=1000)$power)
			rm(allpval)
		}
		all_power[[ilayer]] <- single_power
	}
	kdc_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",layerlist))
	save(kdc_power,file=paste0(respath,"/sim_kdc_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_power.rds"))
	print(kdc_power)
	cat("----------------------------------------------------------------\n")
    cat("----------------------------------------------------------------\n")
    cat("\n")
}



#+ label=SPARK, echo=F, warnings=F, message=F,eval=F
rm(list=setdiff(ls(),"sim_power_func"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
sim_power_func <- function(pval,upfdr=0.2,numTrue=1000,totalGene=10000){
        out         <- cbind(pval,c(rep(1,numTrue),rep(0,totalGene-numTrue)))
        order_out   <- out[order(out[,1]),]
        all_fdr     <- sapply(1:totalGene,function(x){sum(order_out[1:x,2]==0)/x})
        power       <- sapply(seq(0,upfdr,by=0.01),function(y){sum(order_out[,2][which(all_fdr<=y)])/numTrue})
        output      <- cbind.data.frame(fdr=seq(0,upfdr,by=0.01),power=power)
        return(output)
}

outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/result/"

numZero = 2000
sigZero = 200

layerlist <- c("E","ONL","GL")

for(iblock in c(38)){
	load(paste0("/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/sim_data_CN24_D1_",iblock,"X_ONL_zero",numZero,"_sigzero",sigZero,"_means_binR_rpt1.rds"))

	cat(paste0(iblock,"X: number of Spots ",nrow(location),"\n"))
	all_power <- list()
	for(ilayer in 1:3){
		single_power <- c()
		for(isid in 1:10){
			# load(paste0(outpath,"/sim_spark_mixture_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_means_blocked_binR_lmm_rpt",isid,".rds"))
			load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_qn_rpt",isid,".rds"))	
			# RL_pval_combined <- CombinePValues(RL_pval)
			sw_pval_combined <- CombinePValues(SW_pval)

			single_power <- cbind(single_power,sim_power_func(sw_pval_combined,upfdr=0.1,numTrue=100,totalGene=1000)$power)
		}
		all_power[[ilayer]] <- single_power
	}

	spark_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",layerlist))
	save(spark_power,file=paste0(respath,"/sim_spark_mixture_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_qn_lmm_power.rds"))
	print(spark_power)
	cat("----------------------------------------------------------------\n")
    cat("----------------------------------------------------------------\n")
    cat("\n")
}



#+ label=spatialDE, echo=F, warnings=F, message=F,eval=F
# not run yet
rm(list=setdiff(ls(),"sim_power_func"))
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spatialDE/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/result/"

sim_power_func <- function(pval,upfdr=0.2,numTrue=1000,totalGene=10000){
        out         <- cbind(pval,c(rep(1,numTrue),rep(0,totalGene-numTrue)))
        order_out   <- out[order(out[,1]),]
        all_fdr     <- sapply(1:totalGene,function(x){sum(order_out[1:x,2]==0)/x})
        power       <- sapply(seq(0,upfdr,by=0.01),function(y){sum(order_out[,2][which(all_fdr<=y)])/numTrue})
        output      <- cbind.data.frame(fdr=seq(0,upfdr,by=0.01),power=power)
        return(output)
}

numZero = 5000
sigZero = 500

layerlist <- c("E","ONL","GL")
for(iblock in c(38)){
	load(paste0("/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/sim_data_CN24_D1_",iblock,"X_ONL_zero",numZero,"_sigzero",sigZero,"_means_binR_rpt1.rds"))

	cat(paste0(iblock,"X: number of Spots ",nrow(location),"\n"))
	all_power <- list()
	for(ilayer in 1:3){
		single_power <- c()
		for(irpt in 1:10){
			# res <- read.table(paste0(outpath,"/sim_spe_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_means_blocked_binR_rpt",isid,".csv"),head=T)
			res <- read.table(paste0(outpath,"/sim_spe_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_sum_binR_rpt",irpt,".csv"),head=T)
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
		all_power[[ilayer]] <- single_power
	}

	spe_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",layerlist))
	save(spe_power,file=paste0(respath,"/sim_spe_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_sum_binR_power.rds"))
	print(spe_power)
	cat("----------------------------------------------------------------\n")
	cat("----------------------------------------------------------------\n")
}




#' **Power Figures (Fig.1-3): numZero = 2000, sigZero = 200.**   
#+ label=powerfig, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/result/"

power_fig <- list()
icount = 0
for(sigZero in c(100,200,500)){

	icount = icount + 1
	numZero = sigZero*10
	# sigZero = 200
	for(iblock in c(38)){
		
		load(paste0(respath,"/sim_spark_mixture_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_qn_lmm_power.rds"))

		load(paste0(respath,"/sim_kdc_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_power.rds"))

		load(paste0(respath,"/sim_spe_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_sum_binR_power.rds"))

		single_sample_power_fig <- list()

		for(ilayer in 1:3){
			pltdat 		<- cbind.data.frame(KDC=kdc_power[,ilayer+1],SPARK=spark_power[,ilayer+1],SpatialDE=spe_power[,ilayer+1])
			if(ilayer==1){
				single_sample_power_fig[[ilayer]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.7,0.2),
									len.txt.size=1,ax.txt.size=10,ax.title.size=10)+
			                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
			}else{
				single_sample_power_fig[[ilayer]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
									len.txt.size=1,ax.txt.size=10,ax.title.size=10)+
			                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
			}

			rm(pltdat)
		}
	}
	power_fig[[icount]] <- single_sample_power_fig
}

library(ggpubr)
fig1 <- ggarrange(power_fig[[1]][[1]], power_fig[[1]][[2]],power_fig[[1]][[3]], 
			labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)
fig2 <- ggarrange(power_fig[[2]][[1]]+rremove("legend"), power_fig[[2]][[2]],power_fig[[2]][[3]], labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)
fig3 <- ggarrange(power_fig[[3]][[1]]+rremove("legend"), power_fig[[3]][[2]],power_fig[[3]][[3]], labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)


#' ***
#' > Figure 1: 38X, numZero = 1000, sigZero = 100

#+ label=fig1,fig.width=9, fig.height=3, echo=F,fig.align="center",eval=T
fig1

#' ***
#' > Figure 2: 38X, numZero = 2000, sigZero = 200

#+ label=fig2,fig.width=9, fig.height=3, echo=F,fig.align="center",eval=T
fig2


#' ***
#' > Figure 3: 38X, numZero = 5000, sigZero = 500

#+ label=fig3,fig.width=9, fig.height=3, echo=F,fig.align="center",eval=T
fig3


#'***
#' > **KDC with different sample size**
#+ label=power, echo=F, warnings=F, message=F,eval=T

rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/result//"

library(ggsci)
# pal_npg("nrc")(10)

power_fig <- list()

icount = 0
for(sigZero in c(100,200,500)){
	numZero = sigZero*10
	icount = icount + 1
	load(paste0(respath,"/sim_kdc_CN24_D1_38X_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_power.rds"))
	x38 <- kdc_power
	rm(kdc_power)

	load(paste0(respath,"/sim_kdc_CN24_D1_5X_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_power.rds"))
	x5 <- kdc_power
	rm(kdc_power)

	load(paste0(respath,"/sim_kdc_CN24_D1_1X_zero",numZero,"_sigzero",sigZero,"_sum_totaladj_binR_power.rds"))
	xhdst <- kdc_power
	rm(kdc_power)


	single_sample_power_fig <- list()
	ilaycount <- 0
	for(ilayer in c(1,2,3)){
		ilaycount 	<- ilaycount + 1
		pltdat 		<- cbind.data.frame(X38=x38[,ilayer+1],
										X5=x5[,ilayer+1],
										HDST=xhdst[,ilayer+1])
		if(ilaycount==1){
			single_sample_power_fig[[ilaycount]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position=c(0.7,0.2),
								len.txt.size=0.8,ax.txt.size=10,ax.title.size=10,col.base=pal_npg("nrc")(10))+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 2)))+
		                		theme(legend.key.size = unit(0.4, "cm"))
		}else{
			single_sample_power_fig[[ilaycount]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
								len.txt.size=0.8,ax.txt.size=10,ax.title.size=10,col.base=pal_npg("nrc")(10))+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 2)))+
		                		theme(legend.key.size = unit(0.4, "cm"))
		}

		rm(pltdat)
	}
	power_fig[[icount]] <- single_sample_power_fig
	rm(single_sample_power_fig,x38,x5,xhdst)
}


library(ggpubr)
fig1 <- ggarrange(power_fig[[1]][[1]], power_fig[[1]][[2]], power_fig[[1]][[3]],labels = c("A", "B", "C"),font.label=list(size=10),
			ncol = 3, nrow = 1)

fig2 <- ggarrange(power_fig[[2]][[1]]+rremove("legend"), power_fig[[2]][[2]], power_fig[[2]][[3]],labels = c("A", "B", "C"),font.label=list(size=10),
			ncol = 3, nrow = 1)

fig3 <- ggarrange(power_fig[[3]][[1]]+rremove("legend"), power_fig[[3]][[2]], power_fig[[3]][[3]],labels = c("A", "B", "C"),font.label=list(size=10),
			ncol = 3, nrow = 1)


#' ***
#' > Figure 4: numZero = 1000, sigZero = 100 

#+ label=fig4,fig.width=9, fig.height=3, echo=F,fig.align="center"
fig1



#' ***
#' > Figure 5: numZero = 2000, sigZero = 200

#+ label=fig5,fig.width=9, fig.height=3,echo=F,fig.align="center"
fig2


#' ***
#' > Figure 6: numZero = 5000, sigZero = 500

#+ label=fig6,fig.width=9, fig.height=3, echo=F,fig.align="center"
fig3

