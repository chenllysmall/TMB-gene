#install.packages("ggpubr")



library(ggpubr)
library(reshape2)
tmbFile="TMB.txt"     
scoreFile="totalrisk.txt"     
cluFile="geneCluster.txt"         
setwd("C:\\Users\\dell\\Desktop\\COAD\\20. TMB\\1. score TMB")      


tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)       
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)  
clu=read.table(cluFile, header=F, sep="\t", check.names=F, row.names=1)       


tmb=as.matrix(tmb)
tmb[tmb>quantile(tmb,0.99)]=quantile(tmb,0.99)
sameSample=intersect(row.names(tmb), row.names(score))
tmb=tmb[sameSample,]
score=score[sameSample,]
rownames(clu)=gsub("(.*?)\\_(.*?)", "\\2", rownames(clu))
clu=clu[sameSample,]
data=cbind(score, TMB=tmb, geneCluster=clu)

data=data[,c("riskScore", "risk","geneCluster", "TMB")]


data

letter=c("Low Risk","High risk","C","D","E","F","G")

uniqClu=levels(factor(data$geneCluster))

data$geneCluster=letter[match(data$geneCluster, uniqClu)]



group=levels(factor(data$risk))

data$risk=factor(data$risk, levels=c("low", "high"))
data$risk
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


boxplot=ggboxplot(data, x="risk", y="TMB", fill="risk",
		          xlab="",
		          ylab="Tumor Burden Mutation",
		          legend.title="Risk",
		          palette =  c("blue","red") )+ 
  stat_compare_means(comparisons = my_comparisons)
pdf(file="boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()


length=length(levels(factor(data$geneCluster)))
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(riskScore, TMB)) + 
		  xlab("Risk score")+ylab("Tumor Burden Mutation")+
		  geom_point(aes(colour=geneCluster))+
		  scale_color_manual(values=bioCol[1:length])+ 
		  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =riskScore, y =TMB))

pdf(file="cor.pdf", width=6, height=4.5)
print(p1)
dev.off()

