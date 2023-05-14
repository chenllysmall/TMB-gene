#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)                
expFile="tcga.share.txt"     
cliFile="time.txt"            
setwd("C:\\Users\\dell\\Desktop\\COAD\\5. merge time")     


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)


sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli, data)
out=cbind(id=row.names(out), out)
write.table(out, file="TCGA.expTime.txt", sep="\t", row.names=F, quote=F)



