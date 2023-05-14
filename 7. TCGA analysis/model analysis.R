
#install.packages("survivalROC")


library(survivalROC)      
setwd("C:\\Users\\dell\\Desktop\\COAD\\7. TCGA analysis")    
rt=read.table("Trainrisk.txt", header=T, sep="\t", check.names=F, row.names=1)   

#ROC
predictTime=5      
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time =predictTime, method="KM")
pdf(file="ROC.pdf", width=5, height=5)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="black", 
	xlab="False positive rate", ylab="True positive rate",
	lwd = 2, cex.main=1.2, cex.lab=1.2, cex.axis=1.2, font=1.2)
polygon(x=c(0,roc$FP,1,0),y=c(0,roc$TP,1,0),col="skyblue",border=NA)
text(0.85, 0.1, paste0("AUC=",sprintf("%.3f",roc$AUC)), cex=1.2)
segments(0,0,1,1,lty=2)
dev.off()

#
predictTime=5      
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time =predictTime, method="KM")
sum=roc$TP-roc$FP
cutOp=roc$cut.values[which.max(sum)]
cutTP=roc$TP[which.max(sum)]
cutFP=roc$FP[which.max(sum)]
pdf(file="ROC.cutoff.pdf",width=5,height=5)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="black", 
	xlab="False positive rate", ylab="True positive rate",
	lwd = 2, cex.main=1.2, cex.lab=1.2, cex.axis=1.2, font=1.2)
polygon(x=c(0,roc$FP,1,0),y=c(0,roc$TP,1,0),col="skyblue",border=NA)
segments(0,0,1,1,lty=2)
points(cutFP,cutTP, pch=20, col="red",cex=1.5)
text(cutFP+0.15,cutTP-0.05,paste0("Cutoff:",sprintf("%0.3f",cutOp)))
text(0.85, 0.1, paste0("AUC=",sprintf("%.3f",roc$AUC)), cex=1.2)
dev.off()


######
rocCol=c("red", "green", "blue")
aucText=c()
pdf(file="ROC.multiTime.pdf",width=5,height=5)

predictTime=1
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time=predictTime, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
  xlab="False positive rate", ylab="True positive rate",
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("One year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)
#
predictTime=3
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time =predictTime, method="KM")
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
aucText=c(aucText,paste0("Three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
#
predictTime=5
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker=rt$riskScore, predict.time =predictTime, method="KM")
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
aucText=c(aucText,paste0("Five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()


#install.packages("survival")
#install.packages("survminer")



library(survival)
library(survminer)


bioSurvival=function(inputFile=null,outFile=null){

  rt=read.table(inputFile, header=T, sep="\t")
ֵ
  diff=survdiff(Surv(futime, fustat) ~risk, data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%0.3f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  

  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=TRUE,
                     pval=pValue,
                     pval.size=6,
                     palette=c("red", "blue"),
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table=TRUE,
                     risk.table.title="",	  
                     risk.table.height=.25)
 
  pdf(file=outFile,onefile = FALSE,width = 5.5,height =5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="Trainrisk.txt",outFile="survival.pdf")





##Risk plot


inputFile="Trainrisk.txt"             
riskScoreFile="riskScore.pdf"    
survStatFile="survStat.pdf"      


rt=read.table(inputFile, header=T, sep="\t", row.names=1, check.names=F)   
rt=rt[order(rt$riskScore),]       


riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
lowMax=max(rt$riskScore[riskClass=="low"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file=riskScoreFile, width=8, height=3.5)
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("skyblue",lowLength),rep("red",highLength)) )
abline(h=lowMax,v=lowLength,lty=2)
legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","skyblue"),cex=1.1)
dev.off()


color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="skyblue"
pdf(file=survStatFile, width=8, height=3.5)
plot(rt$futime, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","skyblue"),cex=1.1)
abline(v=lowLength,lty=2)
dev.off()





library(pheatmap)

bioRiskPlot=function(inputFile=null, riskScoreFile=null, survStatFile=null, heatmapFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)   
  rt=rt[order(rt$riskScore),]     

rt1=rt[c(3:(ncol(rt)-2))]
rt1=log2(rt1+1)
rt1=t(rt1)

annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file=heatmapFile, width=6, height=2.6)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = F,
         cluster_rows = F,
         show_colnames = F,
         scale="row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         fontsize_col=3,
         fontsize=7,
         fontsize_row=8)
dev.off()
}


bioRiskPlot(inputFile="Trainrisk.txt",heatmapFile="train.heatmap.pdf")






