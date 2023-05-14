#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)          
setwd("C:\\Users\\dell\\Desktop\\COAD\\20. TMB\\2. TMB show")      

#??ȡע???ļ?
score=read.table("Trainrisk.txt", header=T, sep="\t", check.names=F, row.names=1)
clu=read.table("geneCluster.txt", header=F, sep="\t", check.names=F, row.names=1)
rownames(clu)=gsub("(.*?)\\_(.*?)", "\\2", rownames(clu))


sameSample=intersect(rownames(clu), row.names(score))
score=score[sameSample,]
clu=clu[sameSample,]
data=cbind(score[,c("risk", "fustat")], geneCluster=clu)
data$fustat=ifelse(data$fustat==1, "Dead", "Alive")
letter=c("Low","High","C","D","E","F","G")
uniqClu=levels(factor(data$geneCluster))
data$geneCluster=letter[match(data$geneCluster,uniqClu)]
colnames(data)=c("Risk score", "Fustat", "Risk")
outTab=rbind(Tumor_Sample_Barcode=colnames(data), data)
write.table(outTab, file="ann.txt", sep="\t", quote=F, col.names=F)


geneNum=10
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]


ann_colors=list()
col=c("blue", "red")
names(col)=c("Low", "High")
ann_colors[["Risk"]]=col
col=c("#088247", "#FFD121")
names(col)=levels(factor(data$Fustat))
ann_colors[["Fustat"]]=col
bioCol=c("red","blue","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
col=bioCol[1:length(levels(factor(data$Risk)))]
names(col)=levels(factor(data$Risk))
ann_colors[["Risk"]]=col


pdf(file="low.pdf", width=7, height=7)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures=colnames(data), genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()


pdf(file="high.pdf", width=7, height=7)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures=colnames(data), genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

