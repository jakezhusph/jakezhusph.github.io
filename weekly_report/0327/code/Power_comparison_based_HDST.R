#' ---
#' title: "Power Comparison for Simulation Based on HDST"
#' author: "Jiaqiang Zhu"
#' date: "March 13th, 2020"
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

outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/output/kdc/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result/"


numZero = 5000 
sigZero = 500
layerlist <- c("E","ONL","GL")
for(iblock in c(50,25,20,15,10,5,1)){
	if(iblock != 1){
		load(paste0("/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/data/rds/sim_data_CN24_D1_",iblock,"X_ONL_zero",numZero,"_sigzero",sigZero,"_rpt1.rds"))
	}else{
		load(paste0("/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/data/rds/sim_data_CN24_D1_HDST_ONL_zero",numZero,"_sigzero",sigZero,"_rpt1.rds"))
	}
	
	cat(paste0(iblock,"X: number of Spots ",nrow(location),"\n"))
	all_power <- list()
	for(ilayer in 1:3){
		single_power <- c()
		for(isid in 1:10){
			allpval <- c()
			# if(iblock!=1){
				load(paste0(outpath,"sim_kdc_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))
		
			single_power <- cbind(single_power,sim_power_func(KDC$res_mtest$combinedPval,upfdr=0.1,numTrue=100,totalGene=1000)$power)
			rm(allpval)
		}
		all_power[[ilayer]] <- single_power
	}
	kdc_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",layerlist))
	save(kdc_power,file=paste0(respath,"/sim_kdc_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))
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

outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result/"

numZero = 5000
sigZero = 500

layerlist <- c("E","ONL","GL")

for(iblock in c(50,25,20,15,10)){
	load(paste0("/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/data/rds/sim_data_CN24_D1_",iblock,"X_ONL_zero",numZero,"_sigzero",sigZero,"_rpt1.rds"))
	cat(paste0(iblock,"X: number of Spots ",nrow(location),"\n"))
	all_power <- list()
	for(ilayer in 1:3){
		single_power <- c()
		for(isid in 1:10){
			load(paste0(outpath,"/sim_spark_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,"_nofilter.rds"))
			single_power <- cbind(single_power,sim_power_func(spark_pval$combined_pvalue,upfdr=0.1,numTrue=100,totalGene=1000)$power)
		}
		all_power[[ilayer]] <- single_power
	}

	spark_power <- setNames(data.frame(fdr=seq(0,0.1,by=0.01),sapply(all_power,function(x){apply(x,1,mean)})),c("FDR",layerlist))
	# save(spark_power,file=paste0(respath,"/sim_spark_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_power.rds"))
	print(spark_power)
	cat("----------------------------------------------------------------\n")
    cat("----------------------------------------------------------------\n")
    cat("\n")
}

#+ label=spatialDE, echo=F, warnings=F, message=F,eval=F
rm(list=setdiff(ls(),"sim_power_func"))
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/output/spatialDE/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result/"

sim_power_func <- function(pval,upfdr=0.2,numTrue=1000,totalGene=10000){
        out         <- cbind(pval,c(rep(1,numTrue),rep(0,totalGene-numTrue)))
        order_out   <- out[order(out[,1]),]
        all_fdr     <- sapply(1:totalGene,function(x){sum(order_out[1:x,2]==0)/x})
        power       <- sapply(seq(0,upfdr,by=0.01),function(y){sum(order_out[,2][which(all_fdr<=y)])/numTrue})
        output      <- cbind.data.frame(fdr=seq(0,upfdr,by=0.01),power=power)
        return(output)
}

numZero = 2000
sigZero = 200

layerlist <- c("E","ONL","GL")
for(iblock in c(50,25,20,15,10,5)){
	load(paste0("/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/data/rds/sim_data_CN24_D1_",iblock,"X_ONL_zero",numZero,"_sigzero",sigZero,"_rpt1.rds"))
	cat(paste0(iblock,"X: number of Spots ",nrow(location),"\n"))
	all_power <- list()
	for(ilayer in 1:3){
		single_power <- c()
		for(isid in 1:10){
			res 	<- read.table(paste0(outpath,"/sim_spe_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".csv"),head=T)
			res1 	<- res[-which(res$g=="log_total_count"),]

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
	save(spe_power,file=paste0(respath,"/sim_spe_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_power.rds"))
	print(spe_power)
	cat("----------------------------------------------------------------\n")
	cat("----------------------------------------------------------------\n")
}




#' **Power Figures (Fig.1-4): numZero = 2000, sigZero = 200.**   
#+ label=powerfig, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result/"

power_fig <- list()
icount = 0
numZero = 2000
sigZero = 200
for(iblock in c(50,25,20,15,10)){
	icount = icount + 1
	load(paste0(respath,"/sim_spark_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_power.rds"))
	load(paste0(respath,"/sim_kdc_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))

	# kdc_power <- kdc_power[,-3]

	single_sample_power_fig <- list()

	for(ilayer in 1:3){
		pltdat 		<- cbind.data.frame(KDC=kdc_power[,ilayer+1],SPARK=spark_power[,ilayer+1])
		if(ilayer==1){
			single_sample_power_fig[[ilayer]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.7,0.5),
								len.txt.size=1,ax.txt.size=10,ax.title.size=10)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
		}else{
			single_sample_power_fig[[ilayer]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
								len.txt.size=1,ax.txt.size=10,ax.title.size=10)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
		}

		rm(pltdat)
	}
	power_fig[[icount]] <- single_sample_power_fig

}

library(ggpubr)
fig1 <- ggarrange(power_fig[[1]][[1]], power_fig[[1]][[2]], 
			labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig2 <- ggarrange(power_fig[[2]][[1]]+rremove("legend"), power_fig[[2]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig3 <- ggarrange(power_fig[[3]][[1]]+rremove("legend"), power_fig[[3]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig4 <- ggarrange(power_fig[[4]][[1]]+rremove("legend"), power_fig[[4]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig5 <- ggarrange(power_fig[[5]][[1]]+rremove("legend"), power_fig[[5]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

#' ***
#' > Figure 1: 50X, numZero = 2000, sigZero = 200

#+ label=fig1,fig.width=8, fig.height=4, fig.cap="50X Power Comparison",echo=F,fig.align="center",eval=T
fig1


#' ***
#' > Figure 2: 25X, numZero = 2000, sigZero = 200

#+ label=fig2,fig.width=8, fig.height=4, fig.cap="25X Power Comparison",echo=F,fig.align="center",eval=T
fig2



#' ***
#' > Figure 3: 20X, numZero = 2000, sigZero = 200

#+ label=fig3,fig.width=8, fig.height=4,  fig.cap="20X Power Comparison",echo=F,fig.align="center",eval=T
fig3



#' ***
#' > Figure 4: 15X, numZero = 2000, sigZero = 200

#+ label=fig4,fig.width=8, fig.height=4, fig.cap="15X Power Comparison",echo=F,fig.align="center",eval=T
fig4

#' ***
#' > Figure 5: 10X, numZero = 2000, sigZero = 200

#+ label=fig5,fig.width=8, fig.height=4, fig.cap="10X Power Comparison",echo=F,fig.align="center",eval=T
fig5




#' **Power Figures II (Fig.5-8), numZero = 1000, sigZero = 100.**    
#+ label=powerfig2, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result/"

power_fig <- list()
icount = 0
numZero = 1000
sigZero = 100
for(iblock in c(50,25,20,15,10)){
	icount = icount + 1
	load(paste0(respath,"/sim_spark_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_power.rds"))
	load(paste0(respath,"/sim_kdc_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))

	# kdc_power <- kdc_power[,-3]

	single_sample_power_fig <- list()

	for(ilayer in 1:3){
		pltdat 		<- cbind.data.frame(KDC=kdc_power[,ilayer+1],SPARK=spark_power[,ilayer+1])
		if(ilayer==1){
			single_sample_power_fig[[ilayer]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.7,0.5),
								len.txt.size=1,ax.txt.size=10,ax.title.size=10)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
		}else{
			single_sample_power_fig[[ilayer]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
								len.txt.size=1,ax.txt.size=10,ax.title.size=10)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
		}

		rm(pltdat)
	}
	power_fig[[icount]] <- single_sample_power_fig

}

library(ggpubr)
fig1 <- ggarrange(power_fig[[1]][[1]], power_fig[[1]][[2]], 
			labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig2 <- ggarrange(power_fig[[2]][[1]]+rremove("legend"), power_fig[[2]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig3 <- ggarrange(power_fig[[3]][[1]]+rremove("legend"), power_fig[[3]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig4 <- ggarrange(power_fig[[4]][[1]]+rremove("legend"), power_fig[[4]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig5 <- ggarrange(power_fig[[5]][[1]]+rremove("legend"), power_fig[[5]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

#' ***
#' > Figure 6: 50X, numZero = 1000, sigZero = 100

#+ label=fig6,fig.width=8, fig.height=4, fig.cap="50X Power Comparison",echo=F,fig.align="center",eval=T
fig1


#' ***
#' > Figure 7: 25X, numZero = 1000, sigZero = 100

#+ label=fig7,fig.width=8, fig.height=4, fig.cap="25X Power Comparison",echo=F,fig.align="center",eval=T
fig2



#' ***
#' > Figure 8: 20X, numZero = 1000, sigZero = 100

#+ label=fig8,fig.width=8, fig.height=4,  fig.cap="20X Power Comparison",echo=F,fig.align="center",eval=T
fig3


#' ***
#' > Figure 9: 15X, numZero = 1000, sigZero = 100

#+ label=fig9,fig.width=8, fig.height=4, fig.cap="15X Power Comparison",echo=F,fig.align="center",eval=T
fig4

#' ***
#' > Figure 10: 10X, numZero = 1000, sigZero = 100

#+ label=fig10,fig.width=8, fig.height=4, fig.cap="10X Power Comparison",echo=F,fig.align="center",eval=T
fig5



#' **Power Figures III (Fig.9-Fig.13), numZero = 5000, sigZero = 500.**  
#+ label=powerfig3, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result/"

power_fig <- list()
icount = 0
numZero = 5000
sigZero = 500
for(iblock in c(50,25,20,15,10)){
	icount = icount + 1
	load(paste0(respath,"/sim_spark_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_power.rds"))
	load(paste0(respath,"/sim_kdc_CN24_D1_",iblock,"X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))

	# kdc_power <- kdc_power[,-3]

	single_sample_power_fig <- list()

	for(ilayer in 1:3){
		pltdat 		<- cbind.data.frame(KDC=kdc_power[,ilayer+1],SPARK=spark_power[,ilayer+1])
		if(ilayer==1){
			single_sample_power_fig[[ilayer]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.7,0.5),
								len.txt.size=1,ax.txt.size=10,ax.title.size=10)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
		}else{
			single_sample_power_fig[[ilayer]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
								len.txt.size=1,ax.txt.size=10,ax.title.size=10)+
		                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
		}

		rm(pltdat)
	}
	power_fig[[icount]] <- single_sample_power_fig

}

library(ggpubr)
fig1 <- ggarrange(power_fig[[1]][[1]], power_fig[[1]][[2]], 
			labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig2 <- ggarrange(power_fig[[2]][[1]]+rremove("legend"), power_fig[[2]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig3 <- ggarrange(power_fig[[3]][[1]]+rremove("legend"), power_fig[[3]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig4 <- ggarrange(power_fig[[4]][[1]]+rremove("legend"), power_fig[[4]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

fig5 <- ggarrange(power_fig[[5]][[1]]+rremove("legend"), power_fig[[5]][[2]],labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

#' ***
#' > Figure 11: 50X, numZero = 5000, sigZero = 500

#+ label=fig11,fig.width=8, fig.height=4, fig.cap="50X Power Comparison",echo=F,fig.align="center",eval=T
fig1


#' ***
#' > Figure 12: 25X, numZero = 5000, sigZero = 500

#+ label=fig12,fig.width=8, fig.height=4, fig.cap="25X Power Comparison",echo=F,fig.align="center",eval=T
fig2



#' ***
#' > Figure 13: 20X, numZero = 5000, sigZero = 500

#+ label=fig13,fig.width=8, fig.height=4,  fig.cap="20X Power Comparison",echo=F,fig.align="center",eval=T
fig3


#' ***
#' > Figure 14: 15X, numZero = 5000, sigZero = 500

#+ label=fig14,fig.width=8, fig.height=4, fig.cap="15X Power Comparison",echo=F,fig.align="center",eval=T
fig4

#' ***
#' > Figure 15: 10X, numZero = 5000, sigZero = 500

#+ label=fig15,fig.width=8, fig.height=4, fig.cap="10X Power Comparison",echo=F,fig.align="center",eval=T
fig5
