


dir.create('files')
file1<-list.files(list.files(getwd()),full.name=T)
file2<-list.files(getwd())
for (i in 1:length(file1)) file.copy(from=file1[i],to='files')





