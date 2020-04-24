#' ---
#' title: "QQ Examination For RL-Pval"
#' author: "Jiaqiang Zhu"
#' date: "April 20th, 2020"
#' ---


#+ label=simple_simulation, echo=T, warnings=F, message=F,eval=F
rm(list=ls())
library(Matrix)
simpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/data/rds/"
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v3/output/spark/"

ilayer    	<- 1
iblock    	<- 38
numZero 	<- 2000
sigZero 	<- 200
isid    	<- 2
load(paste0(simpath,"/sim_data_CN24_D1_",iblock,"X_E_zero",numZero,"_sigzero",sigZero,"_means_binR_rpt",isid,".rds"))

library(CompQuadForm)
library(pracma)
library(mvtnorm)

numGene = 1000
numSample = nrow(location)

set.seed(1)
tau2 = 0.5 # tau2=1 works, but tau2=0.5 is not working here
e 	<- rmvnorm(numGene, mean = rep(0, numSample), sigma = tau2*diag(numSample))
y  	<- 10 + e

colnames(y) <- paste0("C",1:ncol(y))
rownames(y) <- paste0("gene",1:nrow(y))

library(Rcpp)
source("/net/mulan/disk2/jiaqiang/kdc/code/spark/lmm/sw_lmm/sw_lmm.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/spark/lmm/sw_lmm/sw_lmm.cpp")

ED <- as.matrix(dist(location))
lrang <- ComputeGaussianPL(ED)

res <- list()

iker = 5
cat(iker," ")
kernel_mat <- exp(-ED^2/(2*lrang[iker]^2))

res <- res2 <- list()
for(iker in 1:10){
	cat(iker," ")
	kernel_mat <- exp(-ED^2/(2*lrang[iker]^2))
	# kernel_mat_sq <- kernel_mat%*%kernel_mat
	out 	<- RLskat_lmm(y,kernel_mat)
	out2 	<- sw_test_fast(y,kernel_mat)
	res[[iker]] <- out
	res2[[iker]] <- out2
}

names(res) <- names(res2) <- paste0("p",1:10)
a <- sapply(res,function(x){x$sw})
b <- do.call(cbind,res2)

RL_pval <- a
SW_pval <- b
save(RL_pval,SW_pval,
	file=paste0("/net/mulan/disk2/jiaqiang/kdc/weekReports/04242020/temp_output/normal_simulation_RL_SW.rds"))


#+ label=simple_qq, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
load(paste0("/net/mulan/disk2/jiaqiang/kdc/weekReports/04242020/temp_output/normal_simulation_RL_SW.rds"))
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")

fig1 <- qplot_gg(list(RL=CombinePValues(RL_pval),SW=CombinePValues(SW_pval)),
					legend.position=c(0.25,0.8),
					pt.size=2,
                    ax.txt.size=15,
                    ax.title.size=15,
                    len.txt.size=1.2,
					factor_level=c("RL","SW"),
					 self_label=c("RL","SW"))

#' ***
#' > Figure 1: QQ Plot
#+ label=fig1,fig.width=6, fig.height=6, echo=F,fig.align="center",eval=T
fig1



#+ label=high, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
outpath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/simulation/v4/output/"
cellProp = 0.2
ieffect = 2
itheta = 10
ipat = 1
mu0 = 1.5

col_base <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
pattern_list 	<- c("streak","hotspot")
sw_all_pval_qn <- RL_all_pval_qn <- sw_all_pval_vst <- RL_all_pval_vst <- RL_five_pval_qn <- RL_five_pval_vst <- list()
icount 			<- 0 
for(isample in c(300,1000,3000)){
	icount = icount + 1
	sw_qn_all <- RL_qn_all <- sw_vst_all <- RL_vst_all <- RL_qn_five <- RL_vst_five <- c()
	# cat("SS",isample,", mu:",mu0,",",pattern_list[ipat],"\n")
	for(irpt in 1:10){

		# LMM--QN
		load(paste0(outpath,"/spark/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_qn_rpt",irpt,".rds"))
		sw_qn_all <- cbind(sw_qn_all,CombinePValues(SW_pval))
		RL_qn_all <- cbind(RL_qn_all,CombinePValues(RL_pval))

		RL_qn_five <- cbind(RL_qn_five,CombinePValues(RL_pval[,3:7]))

		rm(RL_pval,SW_pval)

		# LMM--VST
		load(paste0(outpath,"/spark/sim_spark_lmm_",pattern_list[ipat],"_SS",isample,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_vst_rpt",irpt,".rds"))
		sw_vst_all <- cbind(sw_vst_all,CombinePValues(SW_pval))
		RL_vst_all <- cbind(RL_vst_all,CombinePValues(RL_pval))
		RL_vst_five <- cbind(RL_vst_five,CombinePValues(RL_pval[,3:7]))

		rm(RL_pval,SW_pval)
	}

	sw_all_pval_qn[[icount]] <- sw_qn_all[101:1000,]
	RL_all_pval_qn[[icount]] <- RL_qn_all[101:1000,]

	sw_all_pval_vst[[icount]] <- sw_vst_all[101:1000,]
	RL_all_pval_vst[[icount]] <- RL_vst_all[101:1000,]

	RL_five_pval_qn[[icount]] <- RL_qn_five[101:1000,]
	RL_five_pval_vst[[icount]] <- RL_vst_five[101:1000,]

	rm(sw_qn_all,RL_qn_all,sw_vst_all,RL_vst_all,RL_qn_five,RL_vst_five)
}



source("/net/fantasia/home/jiaqiang/mulan_temp/common_function/simple_gg_general.R")
comp_qq1     <- qplot_gg(list(
							SW_QN_ALL=sw_all_pval_qn[[1]],
							RL_QN_ALL=RL_all_pval_qn[[1]],
							SW_VST_ALL=sw_all_pval_vst[[1]],
							RL_VST_ALL=RL_all_pval_vst[[1]],
							RL_QN_FIVE=RL_five_pval_qn[[1]],
							RL_VST_FIVE=RL_five_pval_vst[[1]]),
							cl=0,
							col.base=col_base,
							legend.position=c(0.25,0.8),pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1,
                            factor_level=c("SW_QN_ALL","RL_QN_ALL","SW_VST_ALL","RL_VST_ALL","RL_QN_FIVE","RL_VST_FIVE"),
					        self_label=c("SW_QN_ALL","RL_QN_ALL","SW_VST_ALL","RL_VST_ALL","RL_QN_FIVE","RL_VST_FIVE"))

comp_qq2     <- qplot_gg(list(
							SW_QN_ALL=sw_all_pval_qn[[2]],
							RL_QN_ALL=RL_all_pval_qn[[2]],
							SW_VST_ALL=sw_all_pval_vst[[2]],
							RL_VST_ALL=RL_all_pval_vst[[2]],
							RL_QN_FIVE=RL_five_pval_qn[[2]],
							RL_VST_FIVE=RL_five_pval_vst[[2]]),
							cl=0,
							col.base=col_base,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1,
                            factor_level=c("SW_QN_ALL","RL_QN_ALL","SW_VST_ALL","RL_VST_ALL","RL_QN_FIVE","RL_VST_FIVE"),
					        self_label=c("SW_QN_ALL","RL_QN_ALL","SW_VST_ALL","RL_VST_ALL","RL_QN_FIVE","RL_VST_FIVE"))

comp_qq3     <- qplot_gg(list(
							SW_QN_ALL=sw_all_pval_qn[[3]],
							RL_QN_ALL=RL_all_pval_qn[[3]],
							SW_VST_ALL=sw_all_pval_vst[[3]],
							RL_VST_ALL=RL_all_pval_vst[[3]],
							RL_QN_FIVE=RL_five_pval_qn[[3]],
							RL_VST_FIVE=RL_five_pval_vst[[3]]),
							cl=0,
							col.base=col_base,
							legend.position="none",pt.size=2,
                            ax.txt.size=15,ax.title.size=15,
                            len.txt.size=1,
                            factor_level=c("SW_QN_ALL","RL_QN_ALL","SW_VST_ALL","RL_VST_ALL","RL_QN_FIVE","RL_VST_FIVE"),
					        self_label=c("SW_QN_ALL","RL_QN_ALL","SW_VST_ALL","RL_VST_ALL","RL_QN_FIVE","RL_VST_FIVE"))

fig1 <- comp_qq1 + ggtitle("SS300") + theme(plot.title = element_text(hjust = 0.5,size=15))
fig2 <- comp_qq2 + ggtitle("SS1000") + theme(plot.title = element_text(hjust = 0.5,size=15))
fig3 <- comp_qq3 + ggtitle("SS3000") + theme(plot.title = element_text(hjust = 0.5,size=15))

library(ggpubr)
pp1  <- ggarrange(fig1,fig2,fig3,labels = c("A", "B","C"),font.label=list(size=15),ncol = 3, nrow = 1)

#' ***
#' > Figure 2: mu = 1.5
#+ label=fig2,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
pp1

