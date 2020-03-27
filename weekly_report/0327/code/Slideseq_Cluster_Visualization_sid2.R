#' ---
#' title: "Slideseq Cluster Visualization (180430_6,MT FREE)"
#' author: "Jiaqiang Zhu"
#' date: "March 24th, 2020"
#' ---


#' Numbers of Signal Genes within Each Cluster
#+ label=plot, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
library(ggpubr)
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "Puck_180430_6"
iblock = 1

load(paste0(workdir,"/data/slideseq/",isample,"_zero_removed.rds"))
load(paste0(workdir,"/output/kdc/clean_",isample,"_idx_",iblock,"X_Joint_kdc.rds"))
mtgene <- which(substr(rownames(KDC$res_mtest),start=1,stop=3)=="mt-")
if(length(mtgene)>0){
	pvals <- KDC$res_mtest$combinedPval[-mtgene]
} 
names(pvals)	<- rownames(KDC$res_mtest)[-mtgene]
pBY 			<- p.adjust(pvals,method="BY")

sigGenes 		<- sp_count[names(pBY)[which(pBY<0.05)],]

# log normalize
norm_count <- log(sigGenes+1)
nPC = 50

pc_type = "single"
if(pc_type=="both"){
	# calculate cell, gene simultaneously, slower than prcomp
	library(irlba)
	pca.results <- irlba(A =t(norm_count), nv = nPC)
	rd <- pca.results$v
}else{
	pca.results <- prcomp(t(norm_count), scale = TRUE)
	rd <- pca.results$rotation[,1:nPC]
}


# set.seed(20200324) # kmeans there is some randomness

# Three Clusters
numClust = 3
allfig_list <- list()
for(isid in 1:10){
	set.seed(isid)
	cl <- kmeans(rd, centers = numClust)
	centers_kmeans <- c()
	for(iclust in 1:length(unique(cl$cluster))){
		m <- apply(norm_count[cl$cluster==iclust,],2,mean)
		centers_kmeans <- cbind(centers_kmeans,m)

		ct_rel 	<- apply(centers_kmeans,2,relative_func)
		pltdat  <- cbind.data.frame(location[,1:2],ct_rel)
		colnames(pltdat) <- c("x","y",paste0("clust",1:ncol(ct_rel)))

		pp     	<- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.8,max.pand=1.01,opt="C")})	
	}
	fig1 	<- ggarrange(pp[[1]],pp[[2]],pp[[3]],ncol=3)
	cat("---------------------------------\n")
	cat("Seed",isid,":\n")
	print(table(cl$cluster))
	rm(cl,centers_kmeans,pp)
	allfig_list[[isid]] <- fig1
}



#' ***
#' > Figure 1: Seed 1
#+ label=fig1,fig.width=12, fig.height=4, echo=F,fig.align="center"
allfig_list[[1]]


#' ***
#' > Figure 2: Seed 2
#+ label=fig2,fig.width=12, fig.height=4, echo=F,fig.align="center"
allfig_list[[2]]


#' ***
#' > Figure 3: Seed 3
#+ label=fig3,fig.width=12, fig.height=4, echo=F,fig.align="center"
allfig_list[[3]]


#' ***
#' > Figure 4: Seed 4
#+ label=fig4,fig.width=12, fig.height=4,echo=F,fig.align="center"
allfig_list[[4]]


#' ***
#' > Figure 5: Seed 5
#+ label=fig5,fig.width=12, fig.height=4,echo=F,fig.align="center"
allfig_list[[5]]



#' ***
#' > Figure 6: Seed 6
#+ label=fig6,fig.width=12, fig.height=4,echo=F,fig.align="center"
allfig_list[[6]]




#' ***
#' > Figure 7: Seed 7
#+ label=fig7,fig.width=12, fig.height=4,echo=F,fig.align="center"
allfig_list[[7]]




#' ***
#' > Figure 8: Seed 8
#+ label=fig8,fig.width=12, fig.height=4,echo=F,fig.align="center"
allfig_list[[8]]




#' ***
#' > Figure 9: Seed 9
#+ label=fig9,fig.width=12, fig.height=4,echo=F,fig.align="center"
allfig_list[[9]]




#' ***
#' > Figure 10: Seed 10
#+ label=fig10,fig.width=12, fig.height=4,echo=F,fig.align="center"
allfig_list[[10]]






