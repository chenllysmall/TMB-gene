#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")



library(limma)
library(sva)

tcgaExpFile="mRNAmatrix.txt"       
geoExpFile="geoMatrix.txt"     
geneFile="tcga.diff.txt"       
setwd("C:\\Users\\dell\\Desktop\\COAD\\4. intersect")      

rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)
tcga=log2(tcga+1)


group=sapply(strsplit(colnames(tcga),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tcga=tcga[,group==0]
tcga=t(tcga)
rownames(tcga)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tcga))
tcga=t(avereps(tcga))


rt=read.table(geoExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)


qx=as.numeric(quantile(geo, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
    geo[geo<0]=0
    geo=log2(geo+1)}
geo=normalizeBetweenArrays(geo)

sameGene=intersect(row.names(tcga),row.names(geo))
tcgaOut=tcga[sameGene,]
geoOut=geo[sameGene,]

all=cbind(tcgaOut,geoOut)
batchType=c(rep(1,ncol(tcgaOut)),rep(2,ncol(geoOut)))
outTab=ComBat(all, batchType, par.prior=TRUE)
tcgaOut=outTab[,colnames(tcgaOut)]
tcgaOut[tcgaOut<0]=0
geoOut=outTab[,colnames(geoOut)]
geoOut[geoOut<0]=0

tcgaTab=rbind(ID=colnames(tcgaOut), tcgaOut)
write.table(tcgaTab, file="TCGA.normalize.txt", sep="\t", quote=F, col.names=F)
geoTab=rbind(ID=colnames(geoOut), geoOut)
write.table(geoTab,file="GEO.normalize.txt",sep="\t",quote=F,col.names=F)

gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(tcgaOut))
tcgaShareExp=tcgaOut[sameGene,]
geoShareExp=geoOut[sameGene,]

tcgaShareExp=rbind(ID=colnames(tcgaShareExp),tcgaShareExp)
write.table(tcgaShareExp,file="TCGA.share.txt",sep="\t",quote=F,col.names=F)
geoShareExp=rbind(ID=colnames(geoShareExp),geoShareExp)
write.table(geoShareExp,file="GEO.share.txt",sep="\t",quote=F,col.names=F)


