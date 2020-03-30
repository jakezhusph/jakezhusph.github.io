#' ---
#' title: "HDST vs. Reduced Resolution"
#' author: "Jiaqiang Zhu"
#' date: "March 29th, 2020"
#' ---


#' ***
#+ label = patHDST, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 1
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_HDST_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))

pltdat <- cbind.data.frame(x=location[,1],y=location[,2],
			t(as.matrix(sp_sim_count[c(1:2,101:102),])))

pp    <- lapply(1:6,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})

library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]],
					pp[[3]], pp[[4]],ncol = 4, nrow = 1)

#' ***
#' > Figure 1: HDST Resolution Gene Expression Pattern (E)
#+ label=fig1,fig.width=12, fig.height=3, echo=F,fig.align="center",eval=T
pp_list



#+ label=X5, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 1
iblock    	<- 5
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))

pltdat <- cbind.data.frame(x=location[,1],y=location[,2],
			t(as.matrix(sp_sim_count[c(1:2,101:102),])))


pp    <- lapply(1:4,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=4,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})


library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],pp[[4]],ncol = 4, nrow = 1)


#' ***
#' > Figure 2: 5X Resolution Gene Expression Pattern (E)
#+ label=fig2,fig.width=12, fig.height=3, echo=F,fig.align="center",eval=T
pp_list


#+ label=X10, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 1
iblock    	<- 10
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))

pltdat <- cbind.data.frame(x=location[,1],y=location[,2],
			t(as.matrix(sp_sim_count[c(1:2,101:102),])))


pp    <- lapply(1:4,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=4,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})


library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],pp[[4]],ncol = 4, nrow = 1)


#' ***
#' > Figure 3: 10X Resolution Gene Expression Pattern (E)
#+ label=fig3,fig.width=12, fig.height=3, echo=F,fig.align="center",eval=T
pp_list


#+ label=X25, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 1
iblock    	<- 25
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))

pltdat <- cbind.data.frame(x=location[,1],y=location[,2],
			t(as.matrix(sp_sim_count[c(1:2,101:102),])))


pp    <- lapply(1:4,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=4,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})


library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],pp[[4]],ncol = 4, nrow = 1)



#' ***
#' > Figure 4: 25X Resolution Gene Expression Pattern (E)
#+ label=fig4,fig.width=12, fig.height=3, echo=F,fig.align="center",eval=T
pp_list



#+ label=X50, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 1
iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))

pltdat <- cbind.data.frame(x=location[,1],y=location[,2],
			t(as.matrix(sp_sim_count[c(1:2,101:102),])))


pp    <- lapply(1:4,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=4,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})


library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],pp[[4]],ncol = 4, nrow = 1)


#' ***
#' > Figure 5: 50X Resolution Gene Expression Pattern
#+ label=fig5,fig.width=12, fig.height=3, echo=F,fig.align="center",eval=T
pp_list




#+ label=sim_qq, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

ilayer    	<- 1
# iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

sim_qq <- list()

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
icount <- 0
for(iblock in c(50,25,10,5)){
	icount = icount + 1
	p1 <- p2 <- c()

	for(isid in 1:10){
		load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	
		p1 <- cbind(p1,KDC$res_stest[,1])
		p2 <- cbind(p2,KDC$res_mtest[,1])	
		rm(KDC)
	}

	comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
	                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)

	fig1 		<- comp_qq + ggtitle(paste0("KDC_",iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))

	sim_qq[[icount]]<- fig1
}


#' ***
#' > Figure 6: QQ Plots of Null Gene (50X,25X,10X,5X)
#+ label=sim_qq_figs,fig.width=12, fig.height=3, echo=F,fig.align="center"
sim_qq_figs <- ggarrange(sim_qq[[1]],sim_qq[[2]],sim_qq[[3]],sim_qq[[4]],
							ncol=4,nrow=1)
sim_qq_figs



#+ label=perm_qq, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

ilayer    	<- 1
# iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

perm_qq <- list()

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
icount <- 0
for(iblock in c(50,25,10,5)){
	icount = icount + 1
	p1 <- p2 <- c()

	for(isid in 1:10){
		load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_perm_rpt",isid,".rds"))	
		p1 <- cbind(p1,KDC$res_stest[,1])
		p2 <- cbind(p2,KDC$res_mtest[,1])	
		rm(KDC)
	}

	comp_qq     <- qplot_gg(as.vector(p2[101:1000,]),cl=0,legend.position=c(0.15,0.8),pt.size=2,
	                            ax.txt.size=15,ax.title.size=15,len.txt.size=2)

	fig1 		<- comp_qq + ggtitle(paste0("KDC_",iblock,"X")) + theme(plot.title = element_text(hjust = 0.5,size=15))

	perm_qq[[icount]]<- fig1
}


#' ***
#' > Figure 7: QQ Plots of Permuted Null Gene (50X,25X,10X,5X)
#+ label=perm_qq_figs,fig.width=12, fig.height=3, echo=F,fig.align="center"
perm_qq_figs <- ggarrange(perm_qq[[1]],perm_qq[[2]],perm_qq[[3]],perm_qq[[4]],
							ncol=4,nrow=1)
perm_qq_figs


