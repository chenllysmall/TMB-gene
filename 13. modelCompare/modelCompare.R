#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")



library(limma)
library(survival)
library(survminer)
library(timeROC)

expFile="mRNAmatrix.txt"           
riskFile="Trainrisk.txt"    
geneFiles=c("Cao signature.txt", "Chang signature.txt", "Liang signature.txt", "Wang signature.txt")     
setwd("C:\\Users\\dell\\Desktop\\COAD\\13. modelCompare")             


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
data=data[,group==0]


riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
riskRT=riskRT[,c("futime","fustat","riskScore")]
colnames(riskRT)=c("futime","fustat","TMB-gene signature")

for(i in geneFiles){

	header=unlist(strsplit(i, "\\."))
	gene=read.table(i, header=F, sep="\t", check.names=F)
	sameGene=intersect(as.vector(gene[,1]), row.names(data))
	data1=data[sameGene,]
	colnames(data1)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data1))
	data1=t(data1)
	data1=avereps(data1)
	

	cli=riskRT[,c("futime", "fustat")]
	sameSample=intersect(row.names(data1), row.names(cli))
	data1=data1[sameSample,]
	cli=cli[sameSample,]
	data1=cbind(cli,data1)
	

	multiCox=coxph(Surv(futime, fustat) ~ ., data = data1)
	riskScore=predict(multiCox,type="risk", newdata=data1)
	data1=cbind(data1, riskScore)
	data1=data1[row.names(riskRT),]
	riskRT=cbind(riskRT, data1[,"riskScore"])
	colnames(riskRT)[ncol(riskRT)]=header[[1]]
}


riskOut=rbind(ID=colnames(riskRT), riskRT)
write.table(riskOut, file="risk.models.txt", sep="\t", col.names=F, quote=F)



bioSurvival=function(inputFile=null, outFile=null, varName=null){
	
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	
	rt$Type=ifelse(rt[,varName]>median(rt[,varName]), "high", "low")
	diff=survdiff(Surv(futime, fustat) ~ Type,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ Type, data = rt)
		

	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           title=varName,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=F,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)
	pdf(file=outFile, onefile = FALSE, width=6, height=5)
	print(surPlot)
	dev.off()
}


bioROC=function(inputFile=null, outFile=null, varName=null){

	rt=read.table(inputFile, header=T, sep="\t", check.names=F)

	ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
	               marker=rt[,varName], cause=1,
	               weighting='aalen',
	               times=c(1,3,5), ROC=TRUE)
	pdf(file=outFile, width=5, height=5)
	plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
	text(0.75, 0.24, varName, cex=1.2)
	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("green",'blue','red'),lwd=2,bty = 'n')
	dev.off()
}

for(varName in colnames(riskRT)[3:ncol(riskRT)]){
	bioSurvival(inputFile="risk.models.txt", outFile=paste0("sur.",varName,".pdf"), varName=varName)
	bioROC(inputFile="risk.models.txt", outFile=paste0("ROC.",varName,".pdf"), varName=varName)
}

