###

library("caret")
setwd("C:\\Users\\dell\\Desktop\\COAD\\8. TCGA inner validation\\1 divide")    
rt=read.table("Trainrisk.txt",sep="\t",header=T,check.names=F)
 

inTrain<-createDataPartition(y=rt[,3],p=0.5,list=F)           #p=0.5 
train<-rt[inTrain,]
test<-rt[-inTrain,]



write.table(train,file="riskTrain.txt",sep="\t",quote=F,row.names=F)
write.table(test,file="riskTest.txt",sep="\t",quote=F,row.names=F)
