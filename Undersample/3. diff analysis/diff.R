#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("ggplot2")


#???ð?
library(limma)
library(ggplot2)
library(pheatmap)

logFCfilter=1         
fdrFilter=0.05       
expFile="TCGAExp.txt"      
setwd("C:\\Users\\dell\\Desktop\\COAD\\revised contents\\Undersample\\3. diff analysis")    


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]


group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])      
treatNum=length(group[group==0])     
Type=c(rep(1,conNum), rep(2,treatNum))


outTab=data.frame()
for(i in row.names(data)){
	rt=data.frame(expression=data[i,], Type=Type)
	wilcoxTest=wilcox.test(expression ~ Type, data=rt)
	pvalue=wilcoxTest$p.value
	conGeneMeans=mean(data[i,1:conNum])
	treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	logFC=log2(treatGeneMeans)-log2(conGeneMeans)
	conMed=median(data[i,1:conNum])
	treatMed=median(data[i,(conNum+1):ncol(data)])
	diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)


write.table(outTab,file="tcga.all.txt",sep="\t",row.names=F,quote=F)


outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="tcga.diff.txt",sep="\t",row.names=F,quote=F)


diffExp=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(diffExp,file="tcga.diffExp.txt",sep="\t",col.names=F,quote=F)


hmExp=log2(data[as.vector(outDiff[,1]),]+0.01)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",width=10,height=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
         cluster_cols =F,
         scale="row",
         show_colnames = F,
         show_rownames = T,
         fontsize = 9,
         fontsize_row=1.5,
         fontsize_col=8)
dev.off()


outTab$fdr=as.numeric(outTab$fdr)
outTab$logFC=as.numeric(outTab$logFC)
Significant=ifelse((outTab$fdr<fdrFilter & abs(outTab$logFC)>logFCfilter), ifelse(outTab$logFC>logFCfilter,"Up","Down"), "Not")

p = ggplot(outTab, aes(logFC, -log10(fdr)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("skyblue", "black", "red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

pdf("vol.pdf", width=6.2, height=5.5)
print(p)
dev.off()



