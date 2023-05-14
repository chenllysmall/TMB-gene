#install.packages("survival")
#install.packages("survminer")



library(survival)
library(survminer)
scoreFile="Trainrisk.txt"    
cliFile="tcgaClinical.txt"             
trait="N"                     
setwd("C:\\Users\\dell\\Desktop\\COAD\\11-2. cliSurvival")        


score=read.table(scoreFile, header=T, sep="\t",check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(cli), row.names(score))
score=score[sameSample,]
cli=cli[sameSample,]
data=cbind(futime=score[,1], fustat=score[,2], cli, risk=score[,"risk"])


rt=data[,c("futime", "fustat", trait, "risk")]
rt=rt[(rt[,trait]!="unknow"),]
colnames(rt)=c("futime", "fustat", "clinical", "risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]


for(j in names(tab)){
	rt1=rt[(rt[,"clinical"]==j),]
	tab1=table(rt1[,"risk"])
	tab1=tab1[tab1!=0]
	labels=names(tab1)
	if(length(labels)==2){
		titleName=j
		if((trait=="age") | (trait=="Age") | (trait=="AGE")){
			titleName=paste0("age", j)
		}
		diff=survdiff(Surv(futime, fustat) ~risk,data = rt1)
		pValue=1-pchisq(diff$chisq, df=1)
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f", pValue))
		}
		fit <- survfit(Surv(futime, fustat) ~ risk, data = rt1)
		bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
		bioCol=bioCol[1:length(levels(factor(rt1[,"risk"])))]

		surPlot=ggsurvplot(fit, 
			           data=rt1,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           title=paste0("Patients with ",titleName),
			           legend.title="risk",
			           legend.labs=labels,
			           font.legend=12,
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette=bioCol,
			           risk.table=TRUE,
			       	   risk.table.title="",
			           risk.table.col = "strata",
			           risk.table.height=.25)

		j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
		pdf(file=paste0(trait,"_",j,".pdf"), onefile = FALSE, width = 5.5, height =5)
		print(surPlot)
		dev.off()
	}
}

