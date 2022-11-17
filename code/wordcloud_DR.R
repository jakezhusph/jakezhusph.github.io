rm(list=ls())
library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")

wordnames <- c("PCA","NMF","Deep Learning",
				"MDS","MF","ICA","Diffusion Map",
				"GPLVM")
d 		<- cbind.data.frame(word=wordnames,freq=c(200,150,130,100,100,110,80,60))

jpeg("/net/mulan/disk2/jiaqiang/PersonSite/DR.jpeg",height=1280,width=1920,res=300)
par(mar=c(0.1,0.1,0.1,0.1))
set.seed(2)
wordcloud(words = d$word, freq = d$freq, min.freq = 30,
          max.words=200, random.order=FALSE, rot.per=0.5, 
          colors=brewer.pal(8, "Dark2"))
dev.off()

