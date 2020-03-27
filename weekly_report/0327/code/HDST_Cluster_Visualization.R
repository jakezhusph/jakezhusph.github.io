#' ---
#' title: "HDST Cluster Visualization (CN24_D1,MT FREE)"
#' author: "Jiaqiang Zhu"
#' date: "March 26th, 2020"
#' ---


#' Plot Code
#+ label=plot, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
library(ggpubr)
source("/net/mulan/disk2/jiaqiang/kdc/code/v5/kdc.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/freq_func.R")
source("/net/mulan/disk2/jiaqiang/SDE/mulan_temp/common_function/simple_plot.R")

workdir = "/net/mulan/disk2/jiaqiang/kdc/experiment/manuscript/realdata/v1/"
isample = "CN24_D1"
iblock = 1

load(paste0(workdir,"/data/hdst/CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
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


set.seed(1) # kmeans there is some randomness
# Two Clusters
numClust = 2
cl <- kmeans(rd, centers = numClust)
centers_kmeans <- c()
for(iclust in 1:length(unique(cl$cluster))){
	m <- apply(norm_count[cl$cluster==iclust,],2,mean)
	centers_kmeans <- cbind(centers_kmeans,m)
}

ct_rel 	<- apply(centers_kmeans,2,relative_func)
pltdat  <- cbind.data.frame(location[,1:2],ct_rel)
colnames(pltdat) <- c("x","y",paste0("clust",1:ncol(ct_rel)))

pp     	<- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.8,max.pand=1.01,opt="C")})

fig0 	<- ggarrange(pp[[1]],pp[[2]],ncol=3)
print(table(cl$cluster))
rm(cl,centers_kmeans,pp)

# Three Clusters
numClust = 3
cl <- kmeans(rd, centers = numClust)
centers_kmeans <- c()
for(iclust in 1:length(unique(cl$cluster))){
	m <- apply(norm_count[cl$cluster==iclust,],2,mean)
	centers_kmeans <- cbind(centers_kmeans,m)
}

ct_rel 	<- apply(centers_kmeans,2,relative_func)
pltdat  <- cbind.data.frame(location[,1:2],ct_rel)
colnames(pltdat) <- c("x","y",paste0("clust",1:ncol(ct_rel)))

pp     	<- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.8,max.pand=1.01,opt="C")})
fig1 	<- ggarrange(pp[[1]],pp[[2]],pp[[3]],ncol=3)
print(table(cl$cluster))
rm(cl,centers_kmeans,pp)

# Four Clusters
numClust = 4
cl <- kmeans(rd, centers = numClust)
centers_kmeans <- c()
for(iclust in 1:length(unique(cl$cluster))){
	m <- apply(norm_count[cl$cluster==iclust,],2,mean)
	centers_kmeans <- cbind(centers_kmeans,m)
}

ct_rel 	<- apply(centers_kmeans,2,relative_func)
pltdat  <- cbind.data.frame(location[,1:2],ct_rel)
colnames(pltdat) <- c("x","y",paste0("clust",1:ncol(ct_rel)))

pp     	<- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.8,max.pand=1.01,opt="C")})
fig2 	<- ggarrange(pp[[1]],pp[[2]],pp[[3]],pp[[4]],ncol=3,nrow=2)
print(table(cl$cluster))
rm(cl,centers_kmeans,pp)


# Five Clusters
numClust = 5
cl <- kmeans(rd, centers = numClust)
centers_kmeans <- c()
for(iclust in 1:length(unique(cl$cluster))){
	m <- apply(norm_count[cl$cluster==iclust,],2,mean)
	centers_kmeans <- cbind(centers_kmeans,m)
}

ct_rel 	<- apply(centers_kmeans,2,relative_func)
pltdat  <- cbind.data.frame(location[,1:2],ct_rel)
colnames(pltdat) <- c("x","y",paste0("clust",1:ncol(ct_rel)))

pp     	<- lapply(1:(ncol(pltdat)-2),function(x){pattern_plot_kdc(pltdat,x,main=T,pointsize=1,titlesize=2,min.pand=0.8,max.pand=1.01,opt="C")})
fig3 	<- ggarrange(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],ncol=3,nrow=2)
print(table(cl$cluster))


#' ***
#' > Figure 0: Two Clusters
#+ label=fig0,fig.width=12, fig.height=4, fig.cap="Two Clusters",echo=F,fig.align="center"
fig0


#' ***
#' > Figure 1: Three Clusters
#+ label=fig1,fig.width=12, fig.height=4, fig.cap="Three Clusters",echo=F,fig.align="center"
fig1


#' ***
#' > Figure 2: Four Clusters
#+ label=fig2,fig.width=12, fig.height=8, fig.cap="Four Clusters",echo=F,fig.align="center"
fig2


#' ***
#' > Figure 3: Five Clusters
#+ label=fig3,fig.width=12, fig.height=8, fig.cap="Five Clusters",echo=F,fig.align="center"
fig3





