#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


library(limma)
library(ggpubr)
riskFile="Trainrisk.txt"      
immFile="infiltration_estimation_for_tcga.csv"     
setwd("C:\\Users\\dell\\Desktop\\COAD\\18. immune cell Diff")    


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)","\\1\\-\\2\\-\\3",rownames(immune))
immune=avereps(immune)


sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "risk", drop=F]
immune=immune[sameSample,]
data=cbind(risk, immune)


data$risk=factor(data$risk, levels=c("low", "high"))
type=levels(factor(data[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


for(i in colnames(data)[2:ncol(data)]){

	boxplot=ggboxplot(data, x="risk", y=i, fill="risk",
			          xlab="Risk",
			          ylab=i,
			          legend.title="Risk",
			          palette=c("skyblue", "#FF0000")
			          )+ 
		    stat_compare_means(comparisons=my_comparisons)
	wilcoxTest=wilcox.test(data[,i] ~ data[,"risk"])
	if(wilcoxTest$p.value<0.05){
		j=gsub("/", "-", i)
		pdf(file=paste0(j, ".pdf"), width=4, height=4.5)
		print(boxplot)
		dev.off()
	}
}


