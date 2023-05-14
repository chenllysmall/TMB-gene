#install.packages("ggplot2")
#install.packages("ggpubr")


library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="Trainrisk.txt"    
cliFile="clinical.txt"            
trait="Fustat"                   
setwd("C:\\Users\\dell\\Desktop\\COAD\\19.scoreCli")     


score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,,drop=F], cli[sameSample,,drop=F])


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]

=rt[,c(trait, "risk")]
colnames(rt1)=c("trait", "risk")
df=as.data.frame(table(rt1))
#???y(df, .(risk), transform, percent = Freq/sum(Freq) * 100)
#?ٷֱ?�ly(df, .(risk), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$risk=factor(df$risk, levels=c("low", "high"))
#???ưٷt(df, aes(x = factor(risk), y = percent, fill = trait)) +
	   geom_bar(position = position_stack(), stat = "identity", width = .7) +
	   scale_fill_manual(values=bioCol)+
	   xlab("riskScore")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
	   geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
	   #coord_flip()+
	   theme_bw()
pdf(file="barplot.pdf", width=4, height=5)
print(p)
dev.off()

#???ñȽ,c(trait, "riskScore")]
colnames(rt2)=c("trait", "riskScore")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#???????=ggboxplot(rt2, x="trait", y="riskScore", fill="trait",
		          xlab=trait,
		          ylab="riskScore",
		          legend.title=trait,
		          palette=bioCol
		          )+ 
	    stat_compare_means(comparisons=my_comparisons)
pdf(file="boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()


######Vi