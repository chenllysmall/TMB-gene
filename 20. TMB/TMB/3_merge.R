

file_names <- dir("C:\\Users\\dell\\Desktop\\COAD\\1. download\\TMB\\mutation", pattern = "*masked", recursive = F, full.names = T)
df <- read.delim(file_names[1],comment.char = "#")
for (i in 2:length(file_names)) {
  df <- rbind(df, read.delim(file_names[i],comment.char = "#"))}
write.table(df,"df.txt",sep = "\t",row.names = F,quote = F)








