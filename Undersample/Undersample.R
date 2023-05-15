
setwd("C:\\Users\\dell\\Desktop\\COAD\\revised contents\\Undersample") 
# 随机欠采样
random_undersample <- function(data, minority_class_label, desired_ratio) {
  # 提取多数类和少数类样本
  majority_class <- data[data$label != minority_class_label, ]
  minority_class <- data[data$label == minority_class_label, ]
  
  # 计算需要欠采样的多数类样本数量
  majority_samples <- nrow(minority_class) * desired_ratio
  
  # 随机选择多数类样本
  random_indices <- sample(1:nrow(majority_class), majority_samples)
  undersampled_majority <- majority_class[random_indices, ]
  
  # 合并欠采样后的样本
  undersampled_data <- rbind(minority_class, undersampled_majority)
  
  return(undersampled_data)
}

#1_读入筛选数据
exp <- read.table("mRNAmatrix.txt",header = T,check.names = F,sep = "\t",row.names = 1)

group=sapply(strsplit(colnames(exp),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       
treatNum=length(group[group==0])    
sampleType=c(rep("Normal",conNum), rep("Tumor",treatNum))

exp <- as.data.frame(t(exp))
exp$label <- sampleType

#2_进行随机欠采样提取数据
undersampled_data <- random_undersample(exp, minority_class_label = "Normal", desired_ratio = 1)

exp_filter <- as.data.frame(t(undersampled_data))[1:nrow(t(undersampled_data))-1,]
write.table(rbind(id=colnames(exp_filter),exp_filter),
            "mRNA_filter_matrix.txt",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = F)
