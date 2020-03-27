#' ---
#' title: "Slideseq Enrichment (180430_6,MT FREE)"
#' author: "Jiaqiang Zhu"
#' date: "March 24th, 2020"
#' ---



#+ label=analysis, echo=F, warnings=F, message=F,eval=T
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
# pBY 			<- p.adjust(pvals,method="BY")
# sigGenes 		<- sp_count[names(pBY)[which(pBY<0.05)],]

combined_pvals  <- pvals
combined_stat 	<- abs(tan((0.5 - combined_pvals)*pi))
pBY <- p.adjust(na.omit(combined_pvals),method="BY")

deg <- data.frame(symbol=names(pvals), stats=combined_stat)
deg$g <- ifelse( pBY< 0.05, 'UP', 'stable' )
DEG <- deg


library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
df 	<- bitr(unique(deg$symbol), fromType = "SYMBOL",
	   toType = c( "ENTREZID"),
	   OrgDb = org.Mm.eg.db)

DEG <- merge(DEG,df,by.y='SYMBOL',by.x='symbol')

missing_gene <- setdiff(deg$symbol,DEG$symbol)
missing_signal_gene <- intersect(missing_gene,names(pBY)[pBY< 0.05])
gene_diff= DEG[DEG$g == 'UP','ENTREZID'] 
gene_all = as.character(DEG[ ,'ENTREZID'] )
geneList  = DEG$stats
names(geneList) = DEG$ENTREZID
geneList = sort(geneList,decreasing = T)



#' ***
#' Signal Genes that not mapped to ENTREZID
#+ label=missing_signal_gene, echo=F, warnings=F, message=F,eval=T
missing_signal_gene



#' ***
#' GO: BP MF CC togther
#+ label=GOterm, echo=T, warnings=F, message=F,eval=T
ego <- enrichGO(gene = gene_diff,
				universe      = gene_all,
				OrgDb         = org.Mm.eg.db,
				ont           = "All" ,
				pAdjustMethod = "BH",
				pvalueCutoff  = 0.99,
				qvalueCutoff  = 0.99,
				readable      = TRUE,
				pool 		  = TRUE)

sum(ego$p.adjust<0.05)
print(ego[1:10,3:7])

#' ***
#' > Figure 1: Dot Plots of GO
#+ label=fig1,fig.width=8, fig.height=8, fig.cap="Go Term Dot Plot",echo=F,fig.align="center"
dotplot(ego)

#' ***
#' KEGG
#+ label=KEGG, echo=T, warnings=F, message=F,eval=T
#over-representation test
kk.diff <- enrichKEGG(gene         = gene_diff,
				  organism     = 'mmu',
				  universe     = gene_all,
				  pvalueCutoff = 0.9,
				  qvalueCutoff =0.9)
sum(kk.diff$p.adjust<0.05)
head(kk.diff,10)[,2:6]

#' ***
#' > Figure 2: Dot Plots of KEGG
#+ label=fig2,fig.width=8, fig.height=8, fig.cap="KEGG Dot Plot",echo=F,fig.align="center"
dotplot(kk.diff)

#' ***
#' GSEA
#+ label=GSEA, echo=T, warnings=F, message=F,eval=T
kk_gse <- gseKEGG(geneList     = geneList,
				organism     = 'mmu',
				nPerm        = 1000,
				minGSSize    = 50,
				pvalueCutoff = 0.9,
				verbose      = FALSE)
head(kk_gse,20)[,2:6]
nrow(kk_gse)

#' ***
#' Reactome pathway 
#+ label=Reactome, echo=T, warnings=F, message=F,eval=T
library(ReactomePA)
pp <- enrichPathway(gene = gene_diff,
				 organism = 'mouse',
				 pvalueCutoff = 0.05)
head(pp)[,2:6]
sum(pp$p.adjust<0.05)

