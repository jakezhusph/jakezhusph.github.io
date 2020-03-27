#' ---
#' title: "Slideseq Quality Control Figures (180430_6)"
#' author: "Jiaqiang Zhu"
#' date: "March 26th, 2020"
#' ---



#+ label=plot, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")
datapath = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
sample = "Puck_180430_6"
load(paste0(datapath,"/data/slideseq/Puck_180430_6_zero_removed.rds"))


numZero_gene 	<- as.vector(ncol(sp_count)- sp_nz_count_Rcpp(sp_count,rowSums=T))
num_nonZero_gene <- as.vector(sp_nz_count_Rcpp(sp_count,rowSums=T))
nonZeroPercent 	<- num_nonZero_gene/ncol(sp_count)

mean_gene 		<- as.vector(sp_means_Rcpp(sp_count,rowMeans=T))
var_gene 		<- as.vector(sp_vars_Rcpp(sp_count,rowVars=T))


nFeature_cell 	<- as.vector(sp_nz_count_Rcpp(sp_count,rowSums=F))
totalcount_cell <- as.vector(sp_sums_Rcpp(sp_count))

mt_gene_idx 	<- which(substr(rownames(sp_count),start=1,stop=3)=="mt-")
mt_gene_dat 	<- sp_count[mt_gene_idx,]
mtcount_cell 	<- as.vector(sp_sums_Rcpp(mt_gene_dat))
nFeature_nomt_cell 		<- as.vector(sp_nz_count_Rcpp(sp_count[-mt_gene_idx,],rowSums=F))
totalcount_nomt_cell 	<- as.vector(sp_sums_Rcpp(sp_count[-mt_gene_idx,]))

# gene level stats and figs
genedat 		<- data.frame(numZero=numZero_gene,
					nonZeroCount= num_nonZero_gene,
					nonZeroPercent= nonZeroPercent*100,
					mean_expr=mean_gene,
					var_expr=var_gene)

genedat_nomt 	<- genedat[-mt_gene_idx,]


# cell level stats and figures
celldat 		<- data.frame(totalcount = totalcount_cell,
						nFeature = nFeature_cell,
						percentMT = 100*mtcount_cell/totalcount_cell,
						total_MT_free = totalcount_nomt_cell,
						nFeature_MT_free = nFeature_nomt_cell)


#' Raw Data after Removing Genes Without Any Counts
#+ label=raw,echo=F, warnings=F, message=F,eval=T
cat("Number of Genes:",nrow(sp_count),"Number of Samples:",ncol(sp_count))


#' Non-Zero Counts For Each Gene
#+ label=num_nonZero_gene,echo=F, warnings=F, message=F,eval=T
quantile(num_nonZero_gene,probs=seq(0,1,by=0.1))


#' Percentage of Non-Zero Counts For Each Gene
#+ label=nonZeroPercent,echo=F, warnings=F, message=F,eval=T
quantile(nonZeroPercent,probs=seq(0,1,by=0.1))

#' Total Counts For Each Cell
#+ label=totalcount_cell,echo=F, warnings=F, message=F,eval=T
quantile(totalcount_cell,probs=seq(0,1,by=0.1))


#' Total Counts without Mitochondrial For Each Cell
#+ label=totalcount_nomt_cell,echo=F, warnings=F, message=F,eval=T
quantile(totalcount_nomt_cell,probs=seq(0,1,by=0.1))


#' Number of Genes Expressed in Each Cell
#+ label=nFeature_cell,echo=F, warnings=F, message=F,eval=T
quantile(nFeature_cell,probs=seq(0,1,by=0.1))





library(ggplot2)
library(scales)
library(ggpubr)

#' ***
#' > Figure 0: Variance and Mean Plots For Each Gene
#+ label=fig0,fig.width=12, fig.height=6, echo=F,fig.align="center"
mean_fig1 		<- ggplot(genedat, aes(x=mean_expr, y=var_expr)) + 
					geom_point() + 
					geom_abline(intercept=0,slope=1)+
					labs(title="All Genes",
       					x="Mean", y = "Variance")+
					theme(plot.title = element_text(hjust = 0.5,size=20),
							axis.text.y=element_text(size=15),
							axis.text.x=element_text(size=15),
							axis.title.y=element_text(size=18),
							axis.title.x=element_text(size=18))

mean_fig2 		<- ggplot(genedat_nomt, aes(x=mean_expr, y=var_expr)) + 
					geom_point() + 
					geom_abline(intercept=0,slope=1)+
					labs(title="Without Mitochondrial Genes",
       					x="Mean", y = "Variance")+
					theme(plot.title = element_text(hjust = 0.5,size=20),
							axis.text.y=element_text(size=15),
							axis.text.x=element_text(size=15),
							axis.title.y=element_text(size=18),
							axis.title.x=element_text(size=18))

# var-mean figures
var_mean_fig 	<- ggarrange(mean_fig1,mean_fig2,
							labels = c("A", "B"),
							font.label=list(size=20),
							ncol=2)


var_mean_fig


#' ***
#' > Figure 1: Percentage of Non-Zero Count For Each Gene, B and D are Zoom-In version of A and C at 0%-10% range
#+ label=fig1,fig.width=12, fig.height=12,echo=F,fig.align="center"

zero_fig1 <- ggplot(genedat, aes(y = nonZeroPercent,x="nonZeroPercent",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	theme_bw()+
	labs(title="All Genes")+
	theme(legend.position="none",
			plot.title = element_text(hjust = 0.5,size=20),
			axis.title.y = element_blank(),
			axis.title.x = element_blank(),
			axis.text.y = element_text(size=15),
			axis.text.x = element_text(size=15)) 

zero_fig2 <- ggplot(genedat, aes(y = nonZeroPercent,x="nonZeroPercent",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	# ylim(0, 0.1)+
	ylim(0, 10)+
	theme_bw()+
	labs(title="All Genes")+
	theme(legend.position="none",
			plot.title = element_text(hjust = 0.5,size=20),
			axis.title.y = element_blank(),
			axis.title.x = element_blank(),
			axis.text.y = element_text(size=15),
			axis.text.x = element_text(size=15)) 

zero_fig3 <- ggplot(genedat_nomt, aes(y = nonZeroPercent,x="nonZeroPercent",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	theme_bw()+
	labs(title="Without Mitochondrial")+
	theme(legend.position="none",
			plot.title = element_text(hjust = 0.5,size=20),
			axis.title.y = element_blank(),
			axis.title.x = element_blank(),
			axis.text.y = element_text(size=15),
			axis.text.x = element_text(size=15)) 

zero_fig4 <- ggplot(genedat_nomt, aes(y = nonZeroPercent,x="nonZeroPercent",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	# ylim(0, 0.1)+
	ylim(0, 10)+
	theme_bw()+
	labs(title="Without Mitochondrial")+
	theme(legend.position="none",
			plot.title = element_text(hjust = 0.5,size=20),
			axis.title.y = element_blank(),
			axis.title.x = element_blank(),
			axis.text.y = element_text(size=15),
			axis.text.x = element_text(size=15)) 

# zero percentage figures
zero_all_fig <- ggarrange(zero_fig1,zero_fig2,zero_fig3,zero_fig4,
						labels = c("A", "B", "C","D"),
						font.label=list(size=20),
						ncol=2,nrow=2)

zero_all_fig


#' ***
#' > Figure 2: Total Counts For Each Cell. In the Seurat default setting,cells that have >5% mitochondrial counts are also removed.

#+ label=fig2,fig.width=12, fig.height=4, echo=F,fig.align="center"
tc_fig1 <- ggplot(celldat, aes(y = totalcount,x="Total Counts",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	labs(title="All Genes")+
	theme_bw()+
	theme(legend.position="none",
		plot.title = element_text(hjust = 0.5,size=20),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_text(size=15),
		axis.text.x = element_text(size=15)) 


tc_fig2 <- ggplot(celldat, aes(y = total_MT_free,x="Total Counts",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	labs(title="Without Mitochondrial")+
	theme_bw()+
	theme(legend.position="none",
		plot.title = element_text(hjust = 0.5,size=20),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_text(size=15),
		axis.text.x = element_text(size=15)) 

MT_fig1 <- ggplot(celldat, aes(y = percentMT,x="MT/TotalCounts",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	labs(title="Percentage of MT Gene")+
	theme_bw()+
	theme(legend.position="none",
		plot.title = element_text(hjust = 0.5,size=20),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_text(size=15),
		axis.text.x = element_text(size=15)) 

# total count figures
totalcount_fig <- ggarrange(tc_fig1,tc_fig2,MT_fig1,ncol=3,
							labels = c("A", "B", "C"),
							font.label=list(size=20))


totalcount_fig

#' ***
#' > Figure 3: Number of Genes Expressed in Each Cell. In the Seurat default setting, cells that have unique feature counts over 2,500 or less than 200 are removed.
#+ label=fig3,fig.width=12, fig.height=6, echo=F,fig.align="center"
nf_fig1 <- ggplot(celldat, aes(y = nFeature,x="Number of Expressed Genes",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	labs(title="All Genes")+
	theme_bw()+
	theme(legend.position="none",
		plot.title = element_text(hjust = 0.5,size=20),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_text(size=15),
		axis.text.x = element_text(size=15)) 


nf_fig2 <- ggplot(celldat, aes(y = nFeature_MT_free,x="Number of Expressed Genes",fill=hue_pal()(1))) +
	geom_violin() +   
	geom_jitter(height=0,size=1) +
	labs(title="Without Mitochondrial")+
	theme_bw()+
	theme(legend.position="none",
		plot.title = element_text(hjust = 0.5,size=20),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_text(size=15),
		axis.text.x = element_text(size=15)) 

# numFeature figures
numFeature_fig <- ggarrange(nf_fig1,nf_fig2,
							labels = c("A", "B"),
							font.label=list(size=20),
							ncol=2)

numFeature_fig





