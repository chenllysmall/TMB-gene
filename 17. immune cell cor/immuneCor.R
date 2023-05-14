(!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")



library(limma)
library(pheatmap)
riskFile="Trainrisk.txt"      
immFile="infiltration_estimation_for_tcga.csv"     
setwd("C:\\Users\\dell\\Desktop\\COAD\\26. 17mune cell cor")    


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)


immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)



sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, c("risk", "riskScore")]
immune=immune[sameSample,]
data=cbind(risk, immune)


outTab=data.frame()
sigCell=c("risk","riskScore")
for(i in colnames(data)[3:ncol(data)]){
  if(sd(data[,i])<0.001){next}
  wilcoxTest=wilcox.test(data[,i] ~ data[,"risk"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.05){
    outTab=rbind(outTab,cbind(immune=i, pvalue))
    sigCell=c(sigCell, i)
  }
}
write.table(file="immuneCor.txt", outTab, sep="\t", quote=F, row.names=F)


data=data[,sigCell]
data=data[order(data[,"riskScore"]),]
annCol=data[,1:2]
annCol[,"risk"]=factor(annCol[,"risk"], unique(annCol[,"risk"]))
data=t(data[,(3:ncol(data))])
annRow=sapply(strsplit(rownames(data),"_"), '[', 2)
annRow=as.data.frame(annRow)
row.names(annRow)=row.names(data)
colnames(annRow)=c("Methods")
annRow[,"Methods"]=factor(annRow[,"Methods"], unique(annRow[,"Methods"]))
gapCol=as.vector(cumsum(table(annCol[,"risk"])))
gapRow=as.vector(cumsum(table(annRow[,"Methods"])))


risk=c("blue", "red")
names(risk)=c("low", "high")
ann_colors=list(risk=risk)


pdf("immHeatmap.pdf", width=9, height=6)
pheatmap(data,
         annotation=annCol,
         annotation_row=annRow,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         gaps_row=gapRow,
         gaps_col=gapCol,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=5,
         fontsize_col=6)
dev.off()




#### (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("scales")
#install.packages("ggplot2")
#install.packages("ggtext")


#???rary(limma)
library(scales)
library(ggplot2)
library(ggtext)
riskFile="Trainrisk.txt"      #???File="infiltration_estimation_for_tcga.csv"     #????Èk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#??Èune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)

#?Ô·eSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "riskScore"]
immune=immune[sameSample,]

#?Ô·s.numeric(risk)
outTab=data.frame()
for(i in colnames(immune)){
	y=as.numeric(immune[,i])
	corT=cor.test(x, y, method="spearman")
	cor=corT$estimate
	pvalue=corT$p.value
	if(pvalue<0.05){
		outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
	}
}
#???te.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

#???Result=read.table("corResult.txt", head=T, sep="\t")
corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))     #??(file="cor.pdf", width=10, height=7.5)       #??lot(data=b, aes(x=cor, y=immune, color=Software))+
	labs(x="Correlation coefficient",y="Immune cell")+
	geom_point(size=4.1)+
	theme(panel.background=element_rect(fill="white",size=1,color="black"),
	      panel.grid=element_line(color="grey75",size=0.5),
	      axis.ticks = element_line(size=0.5),
	      axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))
dev.off()

