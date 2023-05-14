#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")



library(limma)
library(ggpubr)
riskFile="Trainrisk.txt"         
cliFile="tcgaClinical.txt"    
setwd("C:\\Users\\dell\\Desktop\\COAD\\11. corelation between model and clinical")    


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)


cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F] 


samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,"riskScore",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk, cli)


for(clinical in colnames(rt[,2:ncol(rt)])){
	data=rt[c("riskScore", clinical)]
	colnames(data)=c("riskScore", "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

	boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
		          xlab=clinical,
		          ylab="Risk score",
		          legend.title=clinical,
		          add = "jitter")+ 
	    stat_compare_means(comparisons = my_comparisons)

	pdf(file=paste0(clinical, ".pdf"), width=4.5, height=4.3)
	print(boxplot)
	dev.off()
}


