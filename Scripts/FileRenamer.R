library(readxl)

input.folder<-"/home/docker/Fastq"

xl.to.read<-list.files(input.folder, pattern = "\\.xlsx", full.names = TRUE)
fastq<-list.files(input.folder, pattern = "\\.fastq", full.names = TRUE)
df<-read_xlsx(xl.to.read)
colnames(df)<-df[1,]
df<-df[-1,]
if(length(which(is.na (df$SequenceID)))>0) df<-df[-which(is.na (df$SequenceID)),]

for (i in 1:nrow(df)) {
  if(length(grep(df$Barcode[i], fastq))>0){
    file.rename(fastq[i], paste(input.folder, df$SequenceID[i], ".fastq", sep = "")) 
  }
}