#' ---
#' title: "Power Comparison for Sparse Simulation (HotSpot and Streak) "
#' author: "Jiaqiang Zhu"
#' date: "April 20th, 2020"
#' ---


#' **Power Figures: Sparse Setting in KDC (Mu = 0.005)**  
#+ label=powerfig1, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

itheta = 0.2
ieffect = 3
mu0 = 0.005
cellProp = 0.3

load(paste0(respath,"/sim_kdc_SS10000_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
p1 <- kdc_power
rm(kdc_power)
load(paste0(respath,"/sim_kdc_SS30000_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
p2 <- kdc_power
rm(kdc_power)


single_sample_power_fig <- list()

for(ipat in 1:2){
	pltdat 		<- cbind.data.frame(SS10000=p1[,ipat+1],
									SS30000=p2[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.2,0.75),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
	}

	rm(pltdat)
}


library(ggpubr)
fig1 <- ggarrange(single_sample_power_fig[[1]], single_sample_power_fig[[2]]+rremove("legend"),labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

#' ***
#' > Figure 1: SampleSize Effect (Phi = 5, Strength = 3, CellProp = 0.3)
#+ label=fig1,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
fig1




#+ label=powerfig2, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

itheta = 0.2
ieffect = 3
mu0 = 0.005
isample = 30000
load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp0.2_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
p1 <- kdc_power
rm(kdc_power)
load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp0.3_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
p2 <- kdc_power
rm(kdc_power)


single_sample_power_fig <- list()

for(ipat in 1:2){
	pltdat 		<- cbind.data.frame(F0.2=p1[,ipat+1],
									F0.3=p2[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.15,0.75),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
	}

	rm(pltdat)
}


library(ggpubr)
fig2 <- ggarrange(single_sample_power_fig[[1]], single_sample_power_fig[[2]]+rremove("legend"),labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

#' ***
#' > Figure 2: Fraction of Marked Cells Effect (Phi = 5, Strength = 3, SampleSize = 30000)
#+ label=fig2,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
fig2




#+ label=powerfig3, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

# itheta = 0.2
ieffect = 3
mu0 = 0.005
cellProp = 0.3
isample = 30000

load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta0.2_power.rds"))
p1 <- kdc_power
rm(kdc_power)
load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta0.4_power.rds"))
p2 <- kdc_power
rm(kdc_power)


single_sample_power_fig <- list()

for(ipat in 1:2){
	pltdat 		<- cbind.data.frame(Phi5=p1[,ipat+1],
									Phi2.5=p2[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.15,0.75),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
	}

	rm(pltdat)
}


library(ggpubr)
fig3 <- ggarrange(single_sample_power_fig[[1]], single_sample_power_fig[[2]]+rremove("legend"),labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

#' ***
#' > Figure 3: Dispersion Effect (SampleSize = 30000, Strength = 3, CellProp = 0.3)
#+ label=fig3,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
fig3




#+ label=powerfig4, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

itheta = 0.2
# ieffect = 3
mu0 = 0.005
cellProp = 0.3
isample = 30000

load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect2_theta",itheta,"_power.rds"))
p1 <- kdc_power
rm(kdc_power)
load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect3_theta",itheta,"_power.rds"))
p2 <- kdc_power
rm(kdc_power)


single_sample_power_fig <- list()

for(ipat in 1:2){
	pltdat 		<- cbind.data.frame(SE2=p1[,ipat+1],
									SE3=p2[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.15,0.75),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
	}

	rm(pltdat)
}


library(ggpubr)
fig4 <- ggarrange(single_sample_power_fig[[1]], single_sample_power_fig[[2]]+rremove("legend"),labels = c("A", "B"),font.label=list(size=10),ncol = 2, nrow = 1)

#' ***
#' > Figure 4: SE Strength Effect (SampleSize = 30000, Phi=5, CellProp = 0.3)
#+ label=fig4,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
fig4



