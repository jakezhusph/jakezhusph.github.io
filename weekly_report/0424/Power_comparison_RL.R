#' ---
#' title: "Power Comparison for Simulation I -- RL (HotSpot and Streak) "
#' author: "Jiaqiang Zhu"
#' date: "April 20th, 2020"
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


orderProb = 0.7
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


#' **Power Figures (Fig.1-4): Sample Size Effect**  
#' Strength = 2, CellProp = 0.2, phi = 0.1. SampleSize varies at 300,1000 and 3000. 
#+ label=powerfig1, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")
cellProp = 0.2
ieffect = 2
itheta = 10
mu0 = 1.5
for(isample in c(300,1000,3000)){
	icount = icount + 1
	load(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
	load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
	lmm_vst_power = lmm_power
	rm(lmm_power)
	load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_power.rds"))
	lmm_qn_power = lmm_power
	rm(lmm_power)
	load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
	load(paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

	single_sample_power_fig <- list()

	for(ipat in 1:2){
		pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
		if(ipat==1){
			single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.3,0.85),
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
fig2 <- ggarrange(power_fig[[1]][[2]], power_fig[[2]][[2]],power_fig[[3]][[2]],labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)

#' ***
#' > Figure 1: Streak (mu = 1.5)

#+ label=fig1,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig1

#' ***
#' > Figure 2: Hotspot (mu = 1.5)
#+ label=fig2,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig2



#+ label=powerfig2, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")
cellProp = 0.2
ieffect = 2
itheta = 10
mu0 = 0.5
for(isample in c(300,1000,3000)){
	icount = icount + 1
	load(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
	load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
	lmm_vst_power = lmm_power
	rm(lmm_power)
	load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_power.rds"))
	lmm_qn_power = lmm_power
	rm(lmm_power)
	load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
	load(paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

	single_sample_power_fig <- list()

	for(ipat in 1:2){
		pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
		if(ipat==1){
			single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.3,0.85),
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
fig3 <- ggarrange(power_fig[[1]][[1]], power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)
fig4 <- ggarrange(power_fig[[1]][[2]], power_fig[[2]][[2]],power_fig[[3]][[2]],labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)

#' ***
#' > Figure 3: Streak (mu = 0.5)

#+ label=fig3,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig3

#' ***
#' > Figure 4: Hotspot (mu = 0.5)
#+ label=fig4,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig4





#' **Power Figures (Fig.5-6): Fraction of Marked Cells Effect**  
#' SS = 1000, strength = 2,  mu = 0.5, phi = 0.1. Fraction of Marked Cells Varies at 0.1, 0.2 and 0.3
#+ label=powerfig3, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

ieffect = 2
itheta = 10
mu0 = 0.5
for(isample in c(1000)){
	for(cellProp in c(0.1,0.2,0.3)){
		icount = icount + 1
		load(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
		lmm_vst_power = lmm_power
		rm(lmm_power)

		load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_power.rds"))
		lmm_qn_power = lmm_power
		rm(lmm_power)

		load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		load(paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		single_sample_power_fig <- list()

		for(ipat in 1:2){
			pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
			if(ipat==1){
				single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.3,0.85),
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
fig5 <- ggarrange(power_fig[[1]][[1]], power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)
fig6 <- ggarrange(power_fig[[1]][[2]], power_fig[[2]][[2]],power_fig[[3]][[2]],labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)

#' ***
#' > Figure 5: Streak

#+ label=fig5,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig5

#' ***
#' > Figure 6: Hotspot
#+ label=fig6,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig6





#' **Power Figures (Fig.7-8): SE Strength Effect**  
#' SS = 1000, cellProp = 0.2,  mu = 0.5, phi = 0.1. SE strength Varies at 1.5, 2 and 3
#+ label=powerfig4, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

itheta = 10
mu0 = 0.5
cellProp = 0.2
for(isample in c(1000)){
	for(ieffect in c(1.5,2,3)){
		icount = icount + 1
		load(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
		lmm_vst_power = lmm_power
		rm(lmm_power)

		load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_power.rds"))
		lmm_qn_power = lmm_power
		rm(lmm_power)

		load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		load(paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		single_sample_power_fig <- list()

		for(ipat in 1:2){
			pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
			if(ipat==1){
				single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.3,0.85),
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
fig7 <- ggarrange(power_fig[[1]][[1]], power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)
fig8 <- ggarrange(power_fig[[1]][[2]]+rremove("legend"), power_fig[[2]][[2]]+rremove("legend"),power_fig[[3]][[2]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)

#' ***
#' > Figure 7: Streak

#+ label=fig7,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig7

#' ***
#' > Figure 8: Hotspot
#+ label=fig8,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig8



#' **Power Figures (Fig.9-10): Dispersion Effect**  
#' SS = 1000, cellProp = 0.2,  mu = 0.5, Strength = 2. Phi Varies at 0.1,0.2,1
#+ label=powerfig5, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

ieffect = 2
mu0 = 0.5
cellProp = 0.2
for(isample in c(1000)){
	for(itheta in c(10,5,1)){
		icount = icount + 1
		load(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
		lmm_vst_power = lmm_power
		rm(lmm_power)

		load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_power.rds"))
		lmm_qn_power = lmm_power
		rm(lmm_power)

		load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		load(paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		single_sample_power_fig <- list()

		for(ipat in 1:2){
			pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
			if(ipat==1){
				single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.3,0.85),
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
fig9 <- ggarrange(power_fig[[1]][[1]], power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)
fig10 <- ggarrange(power_fig[[1]][[2]]+rremove("legend"), power_fig[[2]][[2]]+rremove("legend"),power_fig[[3]][[2]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)

#' ***
#' > Figure 9: Streak

#+ label=fig9,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig9

#' ***
#' > Figure 10: Hotspot
#+ label=fig10,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig10




#' **Power Figures (Fig.11-12): Mu Effect**  
#' SS = 3000, cellProp = 0.2,  itheta = 10, Strength = 3. Phi Varies at 0.1,0.2,1
#+ label=powerfig6, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

power_fig <- list()
icount = 0
pattern_list <- c("streak","hotspot")
col_base  <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
# col_base <- c("hotpink","mediumorchid2","darkred","lightskyblue")

itheta = 10
ieffect = 2
# mu0 = 0.5
cellProp = 0.2

for(isample in c(3000)){
	for(mu0 in c(0.05,0.5,1.5)){
		icount = icount + 1
		
		if(file.exists(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))){
			load(paste0(respath,"/sim_spark_pmm_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))
		}else{
			pmm_power <- cbind(FDR=seq(0,0.1,by=0.01),streak=0,hotspot=0)
		}
	
		load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_power.rds"))
		lmm_vst_power = lmm_power
		rm(lmm_power)

		load(paste0(respath,"/sim_spark_lmm_RL_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_power.rds"))
		lmm_qn_power = lmm_power
		rm(lmm_power)

		load(paste0(respath,"/sim_kdc_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		load(paste0(respath,"/sim_spe_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_power.rds"))

		single_sample_power_fig <- list()

		for(ipat in 1:2){
			pltdat 		<- cbind.data.frame(KDC=kdc_power[,ipat+1],SPARK=pmm_power[,ipat+1],SPARK_LMM_QN=lmm_qn_power[,ipat+1],SPARK_LMM_VST=lmm_vst_power[,ipat+1],SpatialDE=spe_power[,ipat+1])
			if(ipat==1){
				single_sample_power_fig[[ipat]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.3,0.85),
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
fig11 <- ggarrange(power_fig[[1]][[1]], power_fig[[2]][[1]]+rremove("legend"),power_fig[[3]][[1]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)
fig12 <- ggarrange(power_fig[[1]][[2]]+rremove("legend"), power_fig[[2]][[2]]+rremove("legend"),power_fig[[3]][[2]]+rremove("legend"),labels = c("A", "B","C"),font.label=list(size=10),ncol = 3, nrow = 1)

#' ***
#' > Figure 11: Streak

#+ label=fig11,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig11

#' ***
#' > Figure 12: Hotspot
#+ label=fig12,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig12

