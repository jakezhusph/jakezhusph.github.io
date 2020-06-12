#' ---
#' title: "RCTD Confident Cell Distribution"
#' author: "Jiaqiang Zhu"
#' date: "June 12th, 2020"
#' ---


#' **RCTD inferred cell type**  
#' In total, 22207 out of 25551 beads passed the filtering procedures; 14624 beads are rejected, 6717 beads are singlets; 230 beads are certain doublets and 636 beads are uncertain doublets. Singlets and certain doublets constitute the final 7177 cells. 3240 genes are used to determine the cell types. 
#+ label=pattern1, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(RCTD)
library(Matrix)
library(CompQuadForm)
library(parallel)
library(Rcpp)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")
workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"
load(paste0(workdir,"/kdc/fdr_kdc_RCTD_partgenes_decomposed_XTX.rds"))
load(paste0(workdir,"/puck_4kdc.rds"))

RCTD_object <- readRDS(paste0(workdir,"/SplitPuckResults/results_decomposed.RDS"))
cate_x <- c()
for(itype in RCTD_object@cell_type_names){
	cate_x <- cbind(cate_x,as.numeric(RCTD_object@cell_labels==itype))
}
colnames(cate_x) <- RCTD_object@cell_type_names


figpath = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/figure/"
source("/net/mulan/disk2/jiaqiang/kdc/shortcut_functions/freq_func_kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")


pltdat  	<- cbind.data.frame(location[,1:2],cate_x)


colnames(pltdat) <- c("x","y",colnames(cate_x))

pp     <- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc4(pltdat,x,main=T,
								pointsize=0.5,
								titlesize=1.5,
								min.pand=0.8,
								max.pand=1.01,
								pal=colorRampPalette(c("grey","darkblue"))
								)})
library(ggpubr)

fig0 <- ggarrange(pp[[1]],pp[[2]],pp[[3]],pp[[4]],
					pp[[5]],pp[[6]],pp[[7]],pp[[8]],
					pp[[9]],pp[[10]],pp[[11]],pp[[12]],
					pp[[13]],pp[[14]],pp[[15]],pp[[16]],
					pp[[17]],pp[[18]],pp[[19]],
							ncol = 4, nrow = 5)
			
#+ label=fig0,fig.width=12, fig.height=15, echo=F,fig.align="center",eval=T
fig0






