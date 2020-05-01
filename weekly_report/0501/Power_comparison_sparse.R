#' ---
#' title: "Power Comparison"
#' author: "Jiaqiang Zhu"
#' date: "April 30th, 2020"
#' ---


#' **Pattern Visualization**  
#+ label=pattern, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

simpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/data/rds/"

numCell = 1000
cellProp = 0.2 ## 0.1, 0.2, 0.3
irpt = 1
mu0 = 1.5
ieffect = 3
itheta = 5

load(paste0(simpath,"/sim_hotspot_SS",numCell,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))

pltdat1  	<- cbind.data.frame(location[,1:2],relative_func(sp_sim_count[1,]))
pltdat2  	<- cbind.data.frame(location[,1:2],relative_func(sp_sim_count[51,]))
pp1      	<- pattern_plot_kdc2(pltdat1,1,main=T,pointsize=2,
								min.pand=0.8,max.pand=1.01,opt="C",title="Hotspot")
pp2      	<- pattern_plot_kdc2(pltdat2,1,main=F,pointsize=2,
								min.pand=0.8,max.pand=1.01,opt="C",title="Hotspot")

rm(pltdat1,pltdat2,location,sp_sim_count)

load(paste0(simpath,"/sim_streak_SS",numCell,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))

pltdat1  	<- cbind.data.frame(location[,1:2],relative_func(sp_sim_count[1,]))
pltdat2  	<- cbind.data.frame(location[,1:2],relative_func(sp_sim_count[51,]))
pp3      	<- pattern_plot_kdc2(pltdat1,1,main=T,pointsize=2,
								min.pand=0.8,max.pand=1.01,opt="C",title="Streak")
pp4      	<- pattern_plot_kdc2(pltdat2,1,main=F,pointsize=2,
								min.pand=0.8,max.pand=1.01,opt="C",title="Streak")
rm(pltdat1,pltdat2,location,sp_sim_count)


orderProb = 0.5
load(paste0(simpath,"/sim_gradient_SS",numCell,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_rpt",irpt,".rds"))

pltdat1  	<- cbind.data.frame(location[,1:2],relative_func(sp_sim_count[1,]))
pltdat2  	<- cbind.data.frame(location[,1:2],relative_func(sp_sim_count[51,]))
pp5      	<- pattern_plot_kdc2(pltdat1,1,main=T,pointsize=2,
								min.pand=0.8,max.pand=1.01,opt="C",title="Gradient")
pp6      	<- pattern_plot_kdc2(pltdat2,1,main=F,pointsize=2,
								min.pand=0.8,max.pand=1.01,opt="C",title="Gradient")
rm(pltdat1,pltdat2,location,sp_sim_count)

library(ggpubr)
fig0 <- ggarrange(pp1,pp3,pp5,pp2,pp4,pp6,font.label=list(size=10),ncol = 3, nrow = 2)

#+ label=fig0,fig.width=12, fig.height=8, echo=F,fig.align="center",eval=T
fig0



#' **Power Figures: Moderate Expression**   
#+ label=powerfig1, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")
orderProb = 0.4
itheta = 5
mu0 = 0.5
isample = 300


load(paste0(respath,"/sim_spark_pmm_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
lmm_vst_power = lmm_power

load(paste0(respath,"/sim_kdc_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
load(paste0(respath,"/sim_spe_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))

single_sample_power_fig <- list()

for(ipat in 1:1){
	pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.75,0.25),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base,
							self_label=c("SPARK-X","SPARK","SPARK-G","SpatialDE"))+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
	}

	rm(pltdat)
}


gradient_power_fig <- single_sample_power_fig[[1]]



#+ label=powerfig2, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"gradient_power_fig"))
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

ieffect = 2
mu0 = 0.5
cellProp = 0.2
isample = 300
itheta = 5

load(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

load(paste0(respath,"/sim_spark_lmm_SW_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
lmm_vst_power = lmm_power

load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

load(paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

single_sample_power_fig <- list()

for(ipat in 1:2){
	pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.2,0.8),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base,
							self_label=c("SPARK-X","SPARK","SPARK-G","SpatialDE"))+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base,
							self_label=c("SPARK-X","SPARK","SPARK-G","SpatialDE"))+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}

	rm(pltdat)
}



streak_power_fig <- single_sample_power_fig[[1]]
hotspot_power_fig <- single_sample_power_fig[[2]]



library(ggpubr)
fig1 <- ggarrange(streak_power_fig, hotspot_power_fig+rremove("legend"),gradient_power_fig+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)

#' ***
#' SampleSize = 300, Mu=0.5, Moderate Strength, Theta = 0.2, CellProp = 0.2
#+ label=fig1,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig1








#+ label=powerfig3, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")
orderProb = 0.4
itheta = 5
mu0 = 0.5
isample = 1000


load(paste0(respath,"/sim_spark_pmm_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
lmm_vst_power = lmm_power

load(paste0(respath,"/sim_kdc_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
load(paste0(respath,"/sim_spe_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))

single_sample_power_fig <- list()

for(ipat in 1:1){
	pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.75,0.25),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base,
							self_label=c("SPARK-X","SPARK","SPARK-G","SpatialDE"))+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
	}

	rm(pltdat)
}


gradient_power_fig <- single_sample_power_fig[[1]]



#+ label=powerfig4, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),"gradient_power_fig"))
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

ieffect = 2
mu0 = 0.5
cellProp = 0.2
isample = 1000
itheta = 5

load(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

load(paste0(respath,"/sim_spark_lmm_SW_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
lmm_vst_power = lmm_power

load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

load(paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

single_sample_power_fig <- list()

for(ipat in 1:2){
	pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.2,0.8),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base,
							self_label=c("SPARK-X","SPARK","SPARK-G","SpatialDE"))+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base,
							self_label=c("SPARK-X","SPARK","SPARK-G","SpatialDE"))+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}

	rm(pltdat)
}



streak_power_fig <- single_sample_power_fig[[1]]
hotspot_power_fig <- single_sample_power_fig[[2]]



library(ggpubr)
fig2 <- ggarrange(streak_power_fig, hotspot_power_fig+rremove("legend"),gradient_power_fig+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)

#' ***
#' SampleSize = 1000, Mu=0.5, Moderate Strength, Theta = 0.2, CellProp = 0.2
#+ label=fig2,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig2











#+ label=powerfig3, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_gg_general.R")

respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")
orderProb = 0.4
itheta = 5
mu0 = 0.005
isample = 1000

load(paste0(respath,"/sim_spark_lmm_SW_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_vst_power.rds"))
lmm_vst_power = lmm_power

load(paste0(respath,"/sim_kdc_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))
load(paste0(respath,"/sim_spe_gradient_SS",isample,"_orderProb",orderProb,"_mu",mu0,"_theta",itheta,"_power.rds"))

single_sample_power_fig <- list()

for(ipat in 1:1){
	pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
	if(ipat==1){
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.75,0.25),
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base,
							self_label=c("SPARK-X","SPARK","SPARK-G","SpatialDE"))+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
	}else{
		single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position="none",
							len.txt.size=1,ax.txt.size=10,ax.title.size=10,col.base=col_base)+
	                		guides(colour = guide_legend(override.aes = list(shape = 15,size = 10)))
	}

	rm(pltdat)
}


gradient_power_fig <- single_sample_power_fig[[1]]




