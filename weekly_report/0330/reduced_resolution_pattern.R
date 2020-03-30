#' ---
#' title: "HDST vs. 50X Gene Pattern"
#' author: "Jiaqiang Zhu"
#' date: "March 30th, 2020"
#' ---




#+ label=HDST, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 3
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_HDST_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))

pltdat <- cbind.data.frame(x=location[,1],y=location[,2],
			t(as.matrix(sp_sim_count[c(1:3,101:103),])))

pp    <- lapply(1:6,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})

library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],
					pp[[4]], pp[[5]], pp[[6]],ncol = 3, nrow = 2)

#' ***
#' > Figure 1: HDST Resolution 
#+ label=fig1,fig.width=9, fig.height=6, echo=F,fig.align="center",eval=T
pp_list



#' ***
#' KDC HDST *p*-values
#+ label=kdc1x, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

ilayer    	<- 3
iblock    	<- 1
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_stest[,1])
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}
rownames(p2) <- paste0("gene",1:1000)
p2[c(1:3,101:103),10]


#+ label=X50, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

layerlist   <- c("E","ONL","GL")
isample     <- "CN24_D1"

ilayer    	<- 3
iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200
isid 		<- 10

load(paste0(simpath,"/sim_data_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))

pltdat <- cbind.data.frame(x=location[,1],y=location[,2],
			t(as.matrix(sp_sim_count[c(1:3,101:103),])))


pp    <- lapply(1:6,function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=4,titlesize=2,min.pand=0.9,max.pand=1.01,opt="C")})


library(ggpubr)
pp_list <- ggarrange(pp[[1]], pp[[2]], pp[[3]],
					pp[[4]], pp[[5]], pp[[6]],ncol = 3, nrow = 2)


#' ***
#' > Figure 2: 50X Resolution
#+ label=fig2,fig.width=9, fig.height=6, echo=F,fig.align="center",eval=T
pp_list



#' ***
#' KDC 50X *p*-values 
#+ label=kdc50x, echo=F, warnings=F, message=F,eval=T
rm(list=ls())

simpath	= "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/kdc/"

ilayer    	<- 3
iblock    	<- 50
numZero 	<- 2000
sigZero 	<- 200

layerlist 	<- c("E","ONL","GL")

p1 <- p2 <- c()

for(isid in 1:10){
	load(paste0(outpath,"/sim_kdc_joint_CN24_D1_",iblock,"X_",layerlist[ilayer],"_zero",numZero,"_sigzero",sigZero,"_rpt",isid,".rds"))	
	p1 <- cbind(p1,KDC$res_stest[,1])
	p2 <- cbind(p2,KDC$res_mtest[,1])	
	rm(KDC)
}
rownames(p2) <- paste0("gene",1:1000)
p2[c(1:3,101:103),10]

