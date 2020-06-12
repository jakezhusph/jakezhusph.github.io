#' ---
#' title: "RCTD Confident Set Data Exploration"
#' author: "Jiaqiang Zhu"
#' date: "June 12th, 2020"
#' ---


#'***
#' > **7177 cells and 3240 genes**  
#+ label=slideseq, echo=F, warnings=F, message=F,eval=T
rm(list=ls())
library(Rcpp)
library(fitdistrplus)
library(ggplot2)
source("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.R")
sourceCpp("/net/mulan/disk2/jiaqiang/kdc/code/v6/kdc.cpp")

workdir = "/net/mulan/disk2/jiaqiang/kdc/RCTD/RCTD_base_cut20_small/data/SpatialRNA/mydata/"

load(paste0(workdir,"/puck_4kdc.rds"))

numZero_gene    <- as.vector(ncol(sp_count)- sp_nz_count_Rcpp(sp_count,rowSums=T))
num_nonZero_gene <- as.vector(sp_nz_count_Rcpp(sp_count,rowSums=T))
nonZeroPercent  <- num_nonZero_gene/ncol(sp_count)
zeroPercent     <- 1-nonZeroPercent

mean_gene       <- as.vector(sp_means_Rcpp(sp_count,rowMeans=T))
var_gene        <- as.vector(sp_vars_Rcpp(sp_count,rowVars=T))
mt_gene_idx     <- which(substr(rownames(sp_count),start=1,stop=3)=="mt-")

phi1 = coef(nls(var_gene ~ mean_gene + phi * mean_gene^2, start = list(phi = 0.1)))

px1 <- seq(min(mean_gene),max(mean_gene),length.out=100)
py1 <- px1 + phi1*px1^2
phi1_inv <- 1/phi1
expZero <- (phi1_inv/(phi1_inv+mean_gene))^phi1_inv

df  <- cbind.data.frame(meanx=mean_gene,varx=var_gene,zeroPercent=zeroPercent,expZ=expZero)
df2 <- cbind.data.frame(px=px1,py=py1)



x = log10(mean_gene)
y = log10(var_gene)
loessMod10 <- loess(y ~ x,  span=0.15) 
smoothed10 <- predict(loessMod10) 
# lines(x=x[order(x)],y=smoothed10[order(x)],col="red")

df3 <- cbind.data.frame(px=x[order(x)],py=smoothed10[order(x)])


# log
p4  <- ggplot() +
        geom_point(data=df, aes(x=log10(meanx), y=log10(varx)))+
        geom_line(data=df2, aes(x=log10(px), y=log10(py),linetype="dashed"),color="red",size=3)+
        geom_line(data=df3, aes(x=px, y=py,linetype="dashed"),color="green",size=2)+
        geom_abline(aes(intercept=0,slope=1,linetype="solid"),size=2,color="grey")+
        # labs(title="Slideseq Data (Puck_180430_6)",
          labs(title="Mouse Cerebellum Slideseq Data (RCTD)",
            x="Mean", y = "Variance")+
        scale_x_continuous(breaks=-6:2,
                    labels=c(expression(10^-6),
                            expression(10^-5),
                            expression(10^-4),
                            expression(10^-3),
                            expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2)))+
        scale_y_continuous(breaks=-6:3,
                    labels=c(expression(10^-6),
                            expression(10^-5),
                            expression(10^-4),
                            expression(10^-3),
                            expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2),
                            expression(10^3)))+
        # geom_text(aes(x=-2, y=3), label=expression(paste(phi,"=0.2164729")),parse = T,size=5)+
        geom_text(data = data.frame(),aes(x=-1.8, y=3,label="phi==2.979427 *','~~ tilde(mu)==0.004876689"), parse = T,size= 5)+
        theme_bw()+
        theme(plot.title = element_blank(),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=15),
                axis.title.x=element_text(size=15),    
                 panel.border = element_blank(), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"),
                legend.position = "none")


cat("numGene:",nrow(sp_count),"; numSample:",ncol(sp_count)," phi:",phi1,"\n")
#+ echo=F, warnings=F, message=F,eval=T
print(summary(mean_gene))




library(tidyr)
df_long <- gather(df, condition, measurement, zeroPercent:expZ, factor_key=TRUE)
p6  <- ggplot(data=df_long, aes(x=log10(meanx), y=measurement,color=condition)) +
        geom_point(size=2)+
          scale_x_continuous(breaks=-6:2,
                    labels=c(expression(10^-6),
                            expression(10^-5),
                            expression(10^-4),
                            expression(10^-3),
                            expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2)))+
        scale_color_manual(values=c("black", "grey"), 
                       name="Genes",
                       labels=c("Observed", "Expected"))+
       labs(title="Mouse Cerebellum Slideseq Data (RCTD)",
            x="Mean", y = "Fraction Zeros")+
        theme_bw()+
       theme(plot.title = element_blank(),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=15),
                axis.title.x=element_text(size=15),    
                 panel.border = element_blank(), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"),
                 legend.position="none",
                 legend.text= element_text(size=rel(2)),
                  legend.title= element_text(size=rel(2)))







library(tidyr)
df_diff <- df
df_diff$diff <- df$zeroPercent-df$expZ
p10  <- ggplot(data=df_diff, aes(x=log10(meanx), y=diff)) +
        geom_point(size=2 )+
        scale_x_continuous(breaks=-6:2,
                    labels=c(expression(10^-6),
                            expression(10^-5),
                            expression(10^-4),
                            expression(10^-3),
                            expression(10^-2),
                            expression(10^-1),
                            expression(10^0),
                            expression(10^1),
                            expression(10^2)))+
        scale_y_continuous(lim=c(-0.5,0.5))+
       labs(title="Mouse Cerebellum Slideseq Data (RCTD)",
            x="Mean", y = paste0("Difference \n (Observed - Expected)"))+
        theme_bw()+
       theme(plot.title = element_blank(),
                axis.text.y=element_text(size=15),
                axis.text.x=element_text(size=15),
                axis.title.y=element_text(size=15),
                axis.title.x=element_text(size=15),    
                 panel.border = element_blank(), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"),
                 legend.position="none",
                 legend.text= element_text(size=rel(2)),
                  legend.title= element_text(size=rel(2)))


library(ggpubr)

fig0 <- ggarrange(p4,p6,p10,ncol = 3, nrow = 1)
			
#+ label=fig0,fig.width=12, fig.height=4, echo=F,fig.align="center",eval=T
fig0






