#' ---
#' title: "QQ Comparison for Simulation "
#' author: "Jiaqiang Zhu"
#' date: "April 20th, 2020"
#' ---


#+ label=KDC, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/kdc/"
cellProp = 0.2
ieffect = 2
itheta = 10

pattern_list <- c("streak","hotspot")
all_pval <- list()
icount <- 0 
for(isample in c(300,3000)){
	for(mu0 in c(0.5,1.5)){
		for(ipat in 1:2){
			icount = icount + 1
			p1 <- c()
			# cat("SS",isample,", mu:",mu0,",",pattern_list[ipat],"\n")
			for(irpt in 1:10){
				load(paste0(outpath,"/sim_kdc_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))
				p1 <- cbind(p1,KDC$res_mtest$combinedPval)
			}
			all_pval[[icount]] <- p1[101:1000,]
		}

	}
}

names(all_pval) <- paste0(rep(c("SS300","SS3000"),each=4),"_",rep(rep(c("mu0.5","mu1.5"),each=2),2),"_",rep(pattern_list,4))

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(all_pval[1:4],
							cl=0,
							legend.position=c(0.2,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS300_mu0.5_streak","SS300_mu0.5_hotspot","SS300_mu1.5_streak","SS300_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))


fig1 		<- comp_qq + ggtitle("KDC_NULL_SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))


comp_qq2     <- qplot_gg(all_pval[5:8],
							cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS3000_mu0.5_streak","SS3000_mu0.5_hotspot","SS3000_mu1.5_streak","SS3000_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))


fig2 		<- comp_qq2 + ggtitle("KDC_NULL_SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))






#+ label=pmm, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2")))
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"

cellProp = 0.2
ieffect = 2
itheta = 10
icount = 0
pattern_list <- c("streak","hotspot")
all_pval <- list()
for(isample in c(300,3000)){
	for(mu0 in c(0.5,1.5)){
		for(ipat in 1:2){
			icount = icount + 1
			p1 <- c()
			# cat("SS",isample,", mu:",mu0,",",pattern_list[ipat],"\n")
			for(irpt in 1:10){
				load(paste0(outpath,"/sim_spark_pmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))
				p1 <- cbind(p1,spark_pval$combined_pvalue)
			}
			all_pval[[icount]] <- p1[101:1000,]
		}
	}
}

names(all_pval) <- paste0(rep(c("SS300","SS3000"),each=4),"_",rep(rep(c("mu0.5","mu1.5"),each=2),2),"_",rep(pattern_list,4))

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(all_pval[1:4],
							cl=0,
							legend.position=c(0.2,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS300_mu0.5_streak","SS300_mu0.5_hotspot","SS300_mu1.5_streak","SS300_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))


fig3 		<- comp_qq + ggtitle("PMM_NULL_SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(all_pval[5:8],
							cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS3000_mu0.5_streak","SS3000_mu0.5_hotspot","SS3000_mu1.5_streak","SS3000_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))

fig4 		<- comp_qq2 + ggtitle("PMM_NULL_SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))





#+ label=lmm_sw, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2","fig3","fig4")))
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")

cellProp = 0.2
ieffect = 2
itheta = 10
icount = 0
pattern_list <- c("streak","hotspot")
all_pval <- list()

for(isample in c(300,3000)){
	for(mu0 in c(0.5,1.5)){
		for(ipat in 1:2){
			icount = icount + 1
			p1 <- c()
			# cat("SS",isample,", mu:",mu0,",",pattern_list[ipat],"\n")
			for(irpt in 1:10){
				load(paste0(outpath,"/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_rpt",irpt,".rds"))
				p1 <- cbind(p1,CombinePValues(SW_pval))
			}
			all_pval[[icount]] <- p1[101:1000,]
		}
	}
}

names(all_pval) <- paste0(rep(c("SS300","SS3000"),each=4),"_",rep(rep(c("mu0.5","mu1.5"),each=2),2),"_",rep(pattern_list,4))

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(all_pval[1:4],
							cl=0,
							legend.position=c(0.2,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS300_mu0.5_streak","SS300_mu0.5_hotspot","SS300_mu1.5_streak","SS300_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))


fig5 		<- comp_qq + ggtitle("LMM_SW_NULL_SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(all_pval[5:8],
							cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS3000_mu0.5_streak","SS3000_mu0.5_hotspot","SS3000_mu1.5_streak","SS3000_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))

fig6 		<- comp_qq2 + ggtitle("LMM_SW_NULL_SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))





#+ label=lmm_RL, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2","fig3","fig4","fig5","fig6")))
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spark/"
respath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/result/"
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")

cellProp = 0.2
ieffect = 2
itheta = 10
icount = 0
pattern_list <- c("streak","hotspot")
all_pval <- list()
for(isample in c(300,3000)){
	for(mu0 in c(0.5,1.5)){
		for(ipat in 1:2){
			icount = icount + 1
			p1 <- c()
			# cat("SS",isample,", mu:",mu0,",",pattern_list[ipat],"\n")
			for(irpt in 1:10){
				load(paste0(outpath,"/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_rpt",irpt,".rds"))
				p1 <- cbind(p1,CombinePValues(RL_pval))
			}
			all_pval[[icount]] <- p1[101:1000,]
		}
	}
}

names(all_pval) <- paste0(rep(c("SS300","SS3000"),each=4),"_",rep(rep(c("mu0.5","mu1.5"),each=2),2),"_",rep(pattern_list,4))

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(all_pval[1:4],
							cl=0,
							legend.position=c(0.2,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS300_mu0.5_streak","SS300_mu0.5_hotspot","SS300_mu1.5_streak","SS300_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))


fig7		<- comp_qq + ggtitle("LMM_RL_NULL_SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(all_pval[5:8],
							cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS3000_mu0.5_streak","SS3000_mu0.5_hotspot","SS3000_mu1.5_streak","SS3000_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))

fig8 		<- comp_qq2 + ggtitle("LMM_RL_NULL_SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))




#+ label=spatialDE, echo=F, warnings=F, message=F,eval=T
rm(list=setdiff(ls(),c("fig1","fig2","fig3","fig4","fig5","fig6","fig7","fig8")))
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/spatialde/"

cellProp = 0.2
ieffect = 2
itheta = 10
icount = 0
pattern_list <- c("streak","hotspot")
all_pval <- list()
for(isample in c(300,3000)){
	for(mu0 in c(0.5,1.5)){
		for(ipat in 1:2){
			icount = icount + 1
			p1 <- c()
			# cat("SS",isample,", mu:",mu0,",",pattern_list[ipat],"\n")
			for(irpt in 1:10){
				allpval <- c()

				res <- read.table(paste0(outpath,"/sim_spe_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".csv"),head=T)
				res1 <- res[-which(res$g=="log_total_count"),]

				if(sum(duplicated(res1$g))>0){
					clean_res <- res1[-which(duplicated(res1$g)),]
				}else{
					clean_res <- res1
				}

				order_res <- clean_res[order(as.numeric(sapply(strsplit(as.character(clean_res$g),split="gene"),"[[",2))),]

				p1 <- cbind(p1,order_res$pval)
				rm(res,res1,order_res,clean_res)
			}

			all_pval[[icount]] <- p1[101:1000,]
		}
	}
}

names(all_pval) <- paste0(rep(c("SS300","SS3000"),each=4),"_",rep(rep(c("mu0.5","mu1.5"),each=2),2),"_",rep(pattern_list,4))

source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq     <- qplot_gg(all_pval[1:4],
							cl=0,
							legend.position=c(0.2,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS300_mu0.5_streak","SS300_mu0.5_hotspot","SS300_mu1.5_streak","SS300_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))


fig9		<- comp_qq + ggtitle("SpatialDE_NULL_SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))

comp_qq2     <- qplot_gg(all_pval[5:8],
							cl=0,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1.2,
                            factor_level=c("SS3000_mu0.5_streak","SS3000_mu0.5_hotspot","SS3000_mu1.5_streak","SS3000_mu1.5_hotspot"),
					        self_label=c("streak0.5","hotspot0.5","streak1.5","hotspot1.5"))

fig10 		<- comp_qq2 + ggtitle("SpatialDE_NULL_SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))




library(ggpubr)
pp1  <- ggarrange(fig1,fig2,labels = c("A", "B"),font.label=list(size=15),ncol = 2, nrow = 1)
pp2  <- ggarrange(fig3,fig4,labels = c("A", "B"),font.label=list(size=15),ncol = 2, nrow = 1)
pp3  <- ggarrange(fig5,fig6,labels = c("A", "B"),font.label=list(size=15),ncol = 2, nrow = 1)
pp4  <- ggarrange(fig7,fig8,labels = c("A", "B"),font.label=list(size=15),ncol = 2, nrow = 1)
pp5  <- ggarrange(fig9,fig10,labels = c("A", "B"),font.label=list(size=15),ncol = 2, nrow = 1)



#' ***
#' > Figure 1: KDC
#+ label=fig1,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
pp1

#' ***
#' > Figure 2: PMM
#+ label=fig2,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
pp2


#' ***
#' > Figure 3: LMM-SW
#+ label=fig3,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
pp3


#' ***
#' > Figure 4: LMM-RL
#+ label=fig4,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
pp4


#' ***
#' > Figure 5: SpatialDE
#+ label=fig5,fig.width=8, fig.height=4, echo=F,fig.align="center",eval=T
pp5


