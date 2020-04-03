#' ---
#' title: "QQ Plots for Simulation Based on HDST (Average Count)"
#' author: "Jiaqiang Zhu"
#' date: "April 3rd, 2020"
#' ---

#' ***
#' **numZero = 2000, sigZero = 200**
#+ label=KDC2, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/"

ilayer    	<- 1
numZero 	<- 2000
sigZero 	<- 200
layerlist 	<- c("E","ONL","GL")

nullgene_qq <- allgene_qq <- list()
icount = 0
for(iblock in c(50,25,10)){
	icount = icount + 1
	p1 <- p2 <- c()
	for(isid in 1:10){
		load(paste0(outpath,"/kdc/sim_kdc_joint_mean_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	

		load(paste0(outpath,"/spark/sim_spark_mean_mixture_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))


		p1 <- cbind(p1,spark_pval$combined_pval)
		p2 <- cbind(p2,KDC$res_mtest[,1])	
		rm(KDC)
	}		

	if(icount == 1){
		null_qq     <- qplot_gg(list(KDC=as.vector(p2[101:1000,]),SPARK=as.vector(p1[101:1000,])),
							cl=0,legend.position=c(0.25,0.85),pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig1 		<- null_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))


		all_qq 		<- qplot_gg(list(KDC=as.vector(p2),SPARK=as.vector(p1)),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig2 		<- all_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))

	}else{
		null_qq     <- qplot_gg(list(KDC=as.vector(p2[101:1000,]),SPARK=as.vector(p1[101:1000,])),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig1 		<- null_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))

		all_qq 		<- qplot_gg(list(KDC=as.vector(p2),SPARK=as.vector(p1)),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig2 		<- all_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))


	}

	
	nullgene_qq[[icount]] <- fig1
	allgene_qq[[icount]] <- fig2
	rm(fig1,fig2)
}


library(ggpubr)

null_qq_fig 	<- ggarrange(nullgene_qq[[1]],nullgene_qq[[2]],nullgene_qq[[3]],
						ncol=3,nrow=1)


all_qq_fig 		<- ggarrange(allgene_qq[[1]],allgene_qq[[2]],allgene_qq[[3]],
						ncol=3,nrow=1)


#' ***
#' > Figure 1: Null Genes QQ plots, numZero = 2000, sigZero = 200

#+ label=null_qq_fig,fig.width=9, fig.height=3,echo=F,fig.align="center",eval=T
null_qq_fig


#' ***
#' > Figure 2: All Genes QQ plots, numZero = 2000, sigZero = 200

#+ label=all_qq_fig,fig.width=9, fig.height=3, echo=F,fig.align="center",eval=T
all_qq_fig




#' ***
#' **numZero = 1000, sigZero = 100**
#+ label=KDC1, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/"

ilayer    	<- 1
numZero 	<- 1000
sigZero 	<- 100
layerlist 	<- c("E","ONL","GL")

nullgene_qq <- allgene_qq <- list()
icount = 0
for(iblock in c(50,25,10)){
	icount = icount + 1
	p1 <- p2 <- c()
	for(isid in 1:10){
		load(paste0(outpath,"/kdc/sim_kdc_joint_mean_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	

		load(paste0(outpath,"/spark/sim_spark_mean_mixture_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))


		p1 <- cbind(p1,spark_pval$combined_pval)
		p2 <- cbind(p2,KDC$res_mtest[,1])	
		rm(KDC)
	}		

	if(icount == 1){
		null_qq     <- qplot_gg(list(KDC=as.vector(p2[101:1000,]),SPARK=as.vector(p1[101:1000,])),
							cl=0,legend.position=c(0.25,0.85),pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig1 		<- null_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))


		all_qq 		<- qplot_gg(list(KDC=as.vector(p2),SPARK=as.vector(p1)),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig2 		<- all_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))

	}else{
		null_qq     <- qplot_gg(list(KDC=as.vector(p2[101:1000,]),SPARK=as.vector(p1[101:1000,])),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig1 		<- null_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))

		all_qq 		<- qplot_gg(list(KDC=as.vector(p2),SPARK=as.vector(p1)),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig2 		<- all_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))


	}

	
	nullgene_qq[[icount]] <- fig1
	allgene_qq[[icount]] <- fig2
	rm(fig1,fig2)
}


library(ggpubr)

null_qq_fig 	<- ggarrange(nullgene_qq[[1]],nullgene_qq[[2]],nullgene_qq[[3]],
						ncol=3,nrow=1)


all_qq_fig 		<- ggarrange(allgene_qq[[1]],allgene_qq[[2]],allgene_qq[[3]],
						ncol=3,nrow=1)


#' ***
#' > Figure 3: Null Genes QQ plots, numZero = 1000, sigZero = 100

#+ label=null_qq_fig2,fig.width=9, fig.height=3,echo=F,fig.align="center",eval=T
null_qq_fig


#' ***
#' > Figure 4: All Genes QQ plots, numZero = 1000, sigZero = 100

#+ label=all_qq_fig2,fig.width=9, fig.height=3, echo=F,fig.align="center",eval=T
all_qq_fig






#' ***
#' **numZero = 5000, sigZero = 500**
#+ label=KDC5, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/"

ilayer    	<- 1
numZero 	<- 1000
sigZero 	<- 100
layerlist 	<- c("E","ONL","GL")

nullgene_qq <- allgene_qq <- list()
icount = 0
for(iblock in c(50,25,10)){
	icount = icount + 1
	p1 <- p2 <- c()
	for(isid in 1:10){
		load(paste0(outpath,"/kdc/sim_kdc_joint_mean_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	

		load(paste0(outpath,"/spark/sim_spark_mean_mixture_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))


		p1 <- cbind(p1,spark_pval$combined_pval)
		p2 <- cbind(p2,KDC$res_mtest[,1])	
		rm(KDC)
	}		

	if(icount == 1){
		null_qq     <- qplot_gg(list(KDC=as.vector(p2[101:1000,]),SPARK=as.vector(p1[101:1000,])),
							cl=0,legend.position=c(0.25,0.85),pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig1 		<- null_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))


		all_qq 		<- qplot_gg(list(KDC=as.vector(p2),SPARK=as.vector(p1)),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig2 		<- all_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))

	}else{
		null_qq     <- qplot_gg(list(KDC=as.vector(p2[101:1000,]),SPARK=as.vector(p1[101:1000,])),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig1 		<- null_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))

		all_qq 		<- qplot_gg(list(KDC=as.vector(p2),SPARK=as.vector(p1)),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("KDC","SPARK"),
					        self_label=c("KDC","SPARK"))

		fig2 		<- all_qq + ggtitle(paste0(iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))


	}

	
	nullgene_qq[[icount]] <- fig1
	allgene_qq[[icount]] <- fig2
	rm(fig1,fig2)
}


library(ggpubr)

null_qq_fig 	<- ggarrange(nullgene_qq[[1]],nullgene_qq[[2]],nullgene_qq[[3]],
						ncol=3,nrow=1)


all_qq_fig 		<- ggarrange(allgene_qq[[1]],allgene_qq[[2]],allgene_qq[[3]],
						ncol=3,nrow=1)


#' ***
#' > Figure 5: Null Genes QQ plots, numZero = 5000, sigZero = 500

#+ label=null_qq_fig3,fig.width=9, fig.height=3,echo=F,fig.align="center",eval=T
null_qq_fig


#' ***
#' > Figure 6: All Genes QQ plots, numZero = 5000, sigZero = 500

#+ label=all_qq_fig3,fig.width=9, fig.height=3, echo=F,fig.align="center",eval=T
all_qq_fig






#' ***
#' **KDC with 5X and HDST**
#+ label=KDC_large, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/"

ilayer    	<- 1
layerlist 	<- c("E","ONL","GL")

nullgene_qq <- allgene_qq <- list()
icount = 0
for(sigZero in c(100,200,500)){
	numZero = sigZero*10
	icount = icount + 1
	p1 <- p2 <- c()
	for(isid in 1:10){
		load(paste0(outpath,"/kdc/sim_kdc_joint_mean_CN24_D1_5X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	
		p1 <- cbind(p1,KDC$res_mtest[,1])
		rm(KDC)

		load(paste0(outpath,"/kdc/sim_kdc_joint_mean_CN24_D1_1X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))

		p2 <- cbind(p2,KDC$res_mtest[,1])	
		rm(KDC)

	}		

	if(icount == 1){
		null_qq     <- qplot_gg(list(X5=as.vector(p1[101:1000,]),HDST=as.vector(p2[101:1000,])),
							cl=0,legend.position=c(0.25,0.85),pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("X5","HDST"),
					        self_label=c("X5","HDST"))

		fig1 		<- null_qq + ggtitle(paste0(sigZero*10,"|",sigZero)) + theme(plot.title = element_text(hjust = 0.5,size=15))


		all_qq     <- qplot_gg(list(X5=as.vector(p1),HDST=as.vector(p2)),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("X5","HDST"),
					        self_label=c("X5","HDST"))

		fig2 		<- all_qq + ggtitle(paste0(sigZero*10,"|",sigZero))+ theme(plot.title = element_text(hjust = 0.5,size=15))

	}else{
		null_qq     <- qplot_gg(list(X5=as.vector(p1[101:1000,]),HDST=as.vector(p2[101:1000,])),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("X5","HDST"),
					        self_label=c("X5","HDST"))

		fig1 		<- null_qq + ggtitle(paste0(sigZero*10,"|",sigZero)) + theme(plot.title = element_text(hjust = 0.5,size=15))


		all_qq     <- qplot_gg(list(X5=as.vector(p1),HDST=as.vector(p2)),
							cl=0,legend.position="none",pt.size=2,
					        ax.txt.size=10,ax.title.size=10,len.txt.size=1,
					        factor_level=c("X5","HDST"),
					        self_label=c("X5","HDST"))

		fig2 		<- all_qq + ggtitle(paste0(sigZero*10,"|",sigZero)) + theme(plot.title = element_text(hjust = 0.5,size=15))


	}

	
	nullgene_qq[[icount]] <- fig1
	allgene_qq[[icount]] <- fig2
	rm(fig1,fig2)
}


library(ggpubr)

null_qq_fig 	<- ggarrange(nullgene_qq[[1]],nullgene_qq[[2]],nullgene_qq[[3]],
						ncol=3,nrow=1)


all_qq_fig 		<- ggarrange(allgene_qq[[1]],allgene_qq[[2]],allgene_qq[[3]],
						ncol=3,nrow=1)


#' ***
#' > Figure 7: Null Genes QQ plots, KDC 5X and HDST

#+ label=null_qq_fig4,fig.width=9, fig.height=3,echo=F,fig.align="center",eval=T
null_qq_fig


#' ***
#' > Figure 8: All Genes QQ plots,  KDC 5X and HDST

#+ label=all_qq_fig4,fig.width=9, fig.height=3, echo=F,fig.align="center",eval=T
all_qq_fig



