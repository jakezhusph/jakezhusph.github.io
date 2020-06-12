#' ---
#' title: "NMFreg on Slideseq"
#' author: "Jiaqiang Zhu"
#' date: "June 5th, 2020"
#' ---


#' **Repeat 1**  
#+ label=pattern1, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
outpath = "/net/mulan/disk2/jiaqiang/kdc/NMFreg_tutorial/jiaqiang/output/"

source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

irpt = 1

df_clust <- read.csv(paste0(outpath,"/df_clust_seed1_rpt",irpt,".csv"))
bead_deconv_df_norm <- read.csv(paste0(outpath,"/bead_deconv_df_norm_seed1_rpt",irpt,".csv"))


# non_zero <- apply(bead_deconv_df_norm[,-ncol(bead_deconv_df_norm)],2,function(x){sum(x!=0)})
# apply(bead_deconv_df_norm[,-ncol(bead_deconv_df_norm)],2,mean)


pltdat  <- cbind.data.frame(df_clust[,1:2],bead_deconv_df_norm[,-ncol(bead_deconv_df_norm)])

colnames(pltdat) <- c("x","y",colnames(bead_deconv_df_norm)[-ncol(bead_deconv_df_norm)])

pp     <- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc2(pltdat,x,main=T,
								pointsize=0.5,
								titlesize=2,
								min.pand=0.8,
								max.pand=1.01,opt="D",
								direct=-1,legend.position="none",
								barheight=10 )})
library(ggpubr)
fig0 <- ggarrange(pp[[1]],pp[[2]],pp[[5]],pp[[6]],
					pp[[7]],pp[[8]],pp[[9]],pp[[10]],
					font.label=list(size=10),ncol = 4, nrow = 2)
			
#+ label=fig0,fig.width=12, fig.height=8, echo=F,fig.align="center",eval=T
fig0




#' **Repeat 2**  
#+ label=pattern2, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
outpath = "/net/mulan/disk2/jiaqiang/kdc/NMFreg_tutorial/jiaqiang/output/"

source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

irpt = 2

df_clust <- read.csv(paste0(outpath,"/df_clust_seed1_rpt",irpt,".csv"))
bead_deconv_df_norm <- read.csv(paste0(outpath,"/bead_deconv_df_norm_seed1_rpt",irpt,".csv"))

pltdat  <- cbind.data.frame(df_clust[,1:2],bead_deconv_df_norm[,-ncol(bead_deconv_df_norm)])

colnames(pltdat) <- c("x","y",colnames(bead_deconv_df_norm)[-ncol(bead_deconv_df_norm)])

pp     <- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc2(pltdat,x,main=T,
								pointsize=0.5,
								titlesize=2,
								min.pand=0.8,
								max.pand=1.01,opt="D",
								direct=-1,legend.position="none",
								barheight=10 )})
library(ggpubr)
fig0 <- ggarrange(pp[[1]],pp[[2]],pp[[5]],pp[[6]],
					pp[[7]],pp[[8]],pp[[9]],pp[[10]],
					font.label=list(size=10),ncol = 4, nrow = 2)
			
#+ label=fig1,fig.width=12, fig.height=8, echo=F,fig.align="center",eval=T
fig0





#' **Repeat 3**  
#+ label=pattern3, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
outpath = "/net/mulan/disk2/jiaqiang/kdc/NMFreg_tutorial/jiaqiang/output/"

source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

irpt = 3

df_clust <- read.csv(paste0(outpath,"/df_clust_seed1_rpt",irpt,".csv"))
bead_deconv_df_norm <- read.csv(paste0(outpath,"/bead_deconv_df_norm_seed1_rpt",irpt,".csv"))

pltdat  <- cbind.data.frame(df_clust[,1:2],bead_deconv_df_norm[,-ncol(bead_deconv_df_norm)])

colnames(pltdat) <- c("x","y",colnames(bead_deconv_df_norm)[-ncol(bead_deconv_df_norm)])

pp     <- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc2(pltdat,x,main=T,
								pointsize=0.5,
								titlesize=2,
								min.pand=0.8,
								max.pand=1.01,opt="D",
								direct=-1,legend.position="none",
								barheight=10 )})
library(ggpubr)
fig0 <- ggarrange(pp[[1]],pp[[2]],pp[[5]],pp[[6]],
					pp[[7]],pp[[8]],pp[[9]],pp[[10]],
					font.label=list(size=10),ncol = 4, nrow = 2)
			
#+ label=fig2,fig.width=12, fig.height=8, echo=F,fig.align="center",eval=T
fig0



