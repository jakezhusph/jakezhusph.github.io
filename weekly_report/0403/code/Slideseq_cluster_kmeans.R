#' ---
#' title: "Kmeans Cluster for Slideseq"
#' author: "Jiaqiang Zhu"
#' date: "April 3rd, 2020"
#' ---

#' ***
#' Determine the optimal number of the clusters
#+ label=runcode, echo=F, warnings=F, message=F,eval=T
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


ep_search <- function(x){
	df 			<- cbind(1:length(x), x)
	line 		<- df[c(1, nrow(df)),]
	proj 		<- princurve::project_to_curve(df, line)
	optpoint 	<- which.max(proj$dist_ind)-1
	return(optpoint)
}


# this is a simple method for finding the "elbow" of a curve, from
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
opt_dim 	<- ep_search(pca.results$sdev)
rd 			<- pca.results$rotation[,1:opt_dim]

library(cluster)
# https://uc-r.github.io/kmeans_clustering#silo
avg_sil <- function(df,k) {
  km.res <- kmeans(df, centers = k, nstart = 25,iter.max=50)
  ss <- silhouette(km.res$cluster, dist(df))
  mean(ss[, 3])
}


set.seed(1)
library(tidyverse)
k.values 		<- 2:10

# avg_sil_values 	<- map_dbl(k.values, function(x)avg_sil(rd,x))
# plot(k.values, avg_sil_values,
#        type = "b", pch = 19, frame = FALSE, 
#        xlab = "Number of clusters K",
#        ylab = "Average Silhouettes")


# the above code can be done with fviz_nbclust
library(factoextra)
# a gg plot object
sil_plot <- fviz_nbclust(rd, kmeans,iter.max=50, method = "silhouette")
wss_plot <- fviz_nbclust(rd, kmeans,iter.max=50, method = "wss")


# Gap Stat
set.seed(123)
gap_stat <- clusGap(rd, FUN = kmeans, nstart = 25,
                    K.max = 10,iter.max=50, B = 50)

# print(gap_stat, method = "firstmax")
gap_plot <- fviz_gap_stat(gap_stat)


#' ***
#' > Figure 1: Elbow Method
#+ label=Elbow,fig.width=6, fig.height=6,echo=F,fig.align="center",eval=T
wss_plot


#' ***
#' > Figure 2: Silhouette method
#+ label=Silhouette,fig.width=6, fig.height=6,echo=F,fig.align="center",eval=T
sil_plot


#' ***
#' > Figure 3: Gap statistic
#+ label=Gap,fig.width=6, fig.height=6,echo=F,fig.align="center",eval=T
gap_plot






#' ***
#' **Clustergram**  
#' https://www.r-statistics.com/2010/06/clustergram-visualization-and-diagnostics-for-cluster-analysis-r-code/
#+ label=Clustergram, echo=F, warnings=F, message=F,eval=T
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


ep_search <- function(x){
	df 			<- cbind(1:length(x), x)
	line 		<- df[c(1, nrow(df)),]
	proj 		<- princurve::project_to_curve(df, line)
	optpoint 	<- which.max(proj$dist_ind)-1
	return(optpoint)
}


# this is a simple method for finding the "elbow" of a curve, from
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
opt_dim 	<- ep_search(pca.results$sdev)
rd 			<- pca.results$rotation[,1:opt_dim]

source_https <- function(url, ...) {
  # load package
  require(RCurl)
 
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}
source_https("https://raw.github.com/talgalili/R-code-snippets/master/clustergram.r")
iris_dat <- scale(iris[,-5]) # notice I am scaling the vectors)
set.seed(1)


#' ***
#' > Figure 4: Clustergram of the realdata
#+ label=cgslide,fig.width=6, fig.height=6,echo=F,fig.align="center",eval=T
par(cex.lab = 1.1, cex.main = 1.2)
clustergram(rd, k.range = 2:8, line.width = 0.004)


#' ***
#' > Figure 5: Clustergram example
#+ label=cgiris,fig.width=6, fig.height=6,echo=F,fig.align="center",eval=T
clustergram(iris_dat, k.range = 2:8, line.width = 0.004)




