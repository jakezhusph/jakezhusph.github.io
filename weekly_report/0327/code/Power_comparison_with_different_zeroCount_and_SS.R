#' ---
#' title: "Power Comparison with different ZeroCount and SampleSize"
#' author: "Jiaqiang Zhu"
#' date: "March 27th, 2020"
#' ---



#' Power Figures
#+ label=powerfig, echo=F, warnings=F, message=F,eval=T

rm(list=ls())
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v1/result//"

library(ggsci)
# pal_npg("nrc")(10)

power_fig <- list()

icount = 0
for(sigZero in c(100,200,500)){
	numZero = sigZero*10
	icount = icount + 1
	load(paste0(respath,"/sim_kdc_CN24_D1_50X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))
	x50 <- kdc_power
	rm(kdc_power)

	load(paste0(respath,"/sim_kdc_CN24_D1_25X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))
	x25 <- kdc_power
	rm(kdc_power)

	load(paste0(respath,"/sim_kdc_CN24_D1_20X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))
	x20 <- kdc_power
	rm(kdc_power)

	load(paste0(respath,"/sim_kdc_CN24_D1_15X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))
	x15 <- kdc_power
	rm(kdc_power)

	load(paste0(respath,"/sim_kdc_CN24_D1_10X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))
	x10 <- kdc_power
	rm(kdc_power)

	load(paste0(respath,"/sim_kdc_CN24_D1_5X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))
	x5 <- kdc_power
	rm(kdc_power)

	load(paste0(respath,"/sim_kdc_CN24_D1_1X_zero",numZero,"_sigzero",sigZero,"_everyother_power.rds"))
	xhdst <- kdc_power
	rm(kdc_power)


	single_sample_power_fig <- list()
	ilaycount <- 0
	for(ilayer in c(1,2,3)){
		ilaycount 	<- ilaycount + 1
		pltdat 		<- cbind.data.frame(X50=x50[,ilayer+1],
										X25=x25[,ilayer+1],X20 =x20[,ilayer+1],
										X15=x15[,ilayer+1],X10=x10[,ilayer+1],
										X5=x5[,ilayer+1],HDST=xhdst[,ilayer+1])
		if(ilaycount==1){
			single_sample_power_fig[[ilaycount]] <- powerplot_gg(x=seq(0,0.1,by=0.01),Power_dat=pltdat,linesize=2,
								ylim=c(0,1),xlim=c(0,0.105),legend.position=c(.2,0.8),
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
	rm(single_sample_power_fig,x50,x25,x20,x15,x10,x5,xhdst)
}


library(ggpubr)
fig1 <- ggarrange(power_fig[[1]][[1]]+rremove("legend"), power_fig[[1]][[2]], power_fig[[1]][[3]],
			labels = c("A", "B", "C"),
			font.label=list(size=10),
			ncol = 3, nrow = 1)

fig2 <- ggarrange(power_fig[[2]][[1]]+rremove("legend"), power_fig[[2]][[2]], power_fig[[2]][[3]],
			labels = c("A", "B", "C"),
			font.label=list(size=10),
			ncol = 3, nrow = 1)

fig3 <- ggarrange(power_fig[[3]][[1]], power_fig[[3]][[2]], power_fig[[3]][[3]],
			labels = c("A", "B", "C"),
			font.label=list(size=10),
			ncol = 3, nrow = 1)



#' ***
#' > Figure 1: numZero = 5000, sigZero = 500

#+ label=fig1,fig.width=9, fig.height=3, echo=F,fig.align="center"
fig3



#' ***
#' > Figure 2: numZero = 2000, sigZero = 200

#+ label=fig2,fig.width=9, fig.height=3,echo=F,fig.align="center"
fig2


#' ***
#' > Figure 3: numZero = 1000, sigZero = 100

#+ label=fig3,fig.width=9, fig.height=3, echo=F,fig.align="center"
fig1

