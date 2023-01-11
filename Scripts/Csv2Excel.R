library(writexl)

files<-list.files()

if(length(grep(".csv$",files))==1){

  if(length(grep("xlsx$", files))>0){
    file.rename(files[grep("xlsx$",files)], gsub("xlsx$","xlsx_bk",files[grep("xlsx$",files)])) 
  }
  
  file.toget<-files[grep(".csv$",files)]
  if(length(file.toget)==1){
    
    opname<-gsub(".*Exp_","",file.toget)
    opname<-gsub("_[0-9]+.csv","",opname)
    
    df<-read.csv(file.toget, stringsAsFactors = TRUE)
    df$dummy<-NA
    df<-df[,c(1,2,3,5,4)]
    df<-apply(df,2,as.character)
    
    colnames(df)<-c("Posisjon  på PCR-plate",	"SequenceID",	"Barcode",	"Stamme-LabWare No.",	"Kons. (Ct.)")
    df$Barcode<-gsub("BC","barcode",df$Barcode)
    dummy<-as.data.frame(t(colnames(df)))
    colnames(dummy)<-colnames(df)
    df<-rbind(dummy, df)
  
    colnames(df)<-rep(" ",ncol(df))
    colnames(df)[2]<-gsub("-","",Sys.Date())
    write_xlsx(df, paste(opname,".xlsx",sep = ""))

  }
}