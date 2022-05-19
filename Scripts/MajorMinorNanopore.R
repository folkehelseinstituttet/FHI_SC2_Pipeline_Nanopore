library(seqinr)
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(Biostrings)
library(msa)
library("ggplot2")
library(nnet)
library(ggpubr)
library(writexl)
library(seqinr)

 cores<-10

results <-"/Noise/"   #Docker
temp <-paste(results,"rawnoise/",sep = "")   #Docker

bamfiles<-list.files(results, pattern = ".bam$", full.names = TRUE, recursive = TRUE)

if(length(bamfiles)>0){
  
  samp<-c(1:length(bamfiles))
  
  pb <- progress_bar$new(
    format = "Index: :samp.pb [:bar] :elapsed | eta: :eta",
    total = length(bamfiles),    # 100 
    width = 60)
  
  
  progress <- function(n){
    pb$tick(tokens = list(samp.pb = samp[n]))
  } 
  
  opts <- list(progress = progress)
  
  if(!dir.exists(temp)) dir.create(temp)
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  out.par<-foreach(i=1:length(bamfiles), .verbose=FALSE, .options.snow = opts) %dopar%{
    
    try(system(paste("cd /Noise && FINex -f ",bamfiles[i], " > ",
                     gsub(".*/",temp,gsub("\\.bam", "_NoisExtractorResult.tsv",bamfiles[i])), sep = "")))
    
  }
  stopCluster(cluster.cores)
  
  results.files<-list.files(temp, pattern = "_NoisExtractorResult\\.tsv$", full.names = TRUE)


 
  co<-0.1
  pb<-txtProgressBar(min = 0, max = length(results.files),initial = 0) #01072021 Nacho Garcia
  try(rm(extended.out))

  for (i in 1:length(results.files)) {
    
    try(rm(df), silent = TRUE)
    try(df<-read.csv(results.files[i],sep = "\t", header = FALSE), silent = TRUE)
   
    if(exists("df")){

      df<-df[,c(1,2,3,4,5,6,7)]
      
      colnames(df)<-c("Base","Noise","Reads","S1","S2","FreqMinor","NoiseMinor")
      df$NonMR<-df$Reads*df$Noise
      df$MR<-df$Reads-df$NonMR
      
      #stable.pos<-stable.pos[which(stable.pos %in% df$Base)]
      
      
      genome.position<-as.data.frame(c(1:29903))
      colnames(genome.position)<-"Base"
      df<-merge(genome.position, df, by="Base", all=TRUE)
      df$Outlier<-"NO"
      
      lowerq = quantile(df$FreqMinor,na.rm=TRUE)[2]
      upperq = quantile(df$FreqMinor,na.rm=TRUE)[4]
      iqr = upperq - lowerq 
      co<-(iqr * 20) + upperq
      
      df$Outlier[which(df$FreqMinor>co)]<-"YES"
      
      df$Reads[which(is.na(df$Reads))]<-0
      
      df$NoiseNP<-0
      df$NoiseNP[which(df$Outlier=="YES")]<-df$FreqMinor[which(df$Outlier=="YES")]
      
      if(length(which(df$Outlier=="YES"))>10 ){
        dummy<-df
        if(!dir.exists(paste(results,"fasta/",sep = ""))) try(dir.create(paste(results,"fasta/",sep = "")))
        
        genome.position<-as.data.frame(c(1:29903))
        colnames(genome.position)<-"Base"
        dummy<-merge(dummy,genome.position,by="Base",all.y=TRUE)
        dummy$S1[which(is.na(dummy$S1))]<-"N"
        dummy$S2[which(is.na(dummy$S2))]<-"N"
        dummy$S2[which(dummy$Reads<20)]<-"N"
        dummy$S1[which(dummy$Reads<20)]<-"N"
        
        dummy$S2[which(df$Outlier!="YES")]<-dummy$S1[which(df$Outlier!="YES")]
        
        if(length(which(dummy$S1!=dummy$S2))>5 & length(which(dummy$S1=="N"))<1000){
          
          write.fasta(paste(dummy$S1[which(dummy$S1 %in% c("A","T","C","G","N"))], collapse = ""),
                      file.out = gsub("_NoisExtractorResult.tsv","_S1.fa", gsub(".*/",paste(results,"fasta/",sep = ""),results.files[i])), 
                      names = gsub("_NoisExtractorResult.tsv","_S1.fa", gsub(".*/","",results.files[i])))
          
          write.fasta(paste(dummy$S2[which(dummy$S2 %in% c("A","T","C","G","N"))], collapse = ""),
                      file.out = gsub("_NoisExtractorResult.tsv","_S2.fa", gsub(".*/",paste(results,"fasta/",sep = ""),results.files[i])), 
                      names = gsub("_NoisExtractorResult.tsv","_S2.fa", gsub(".*/","",results.files[i])))   
        }
      

      

    }
    }
  }
    
    
    if(length(list.files("/Noise/fasta/"))>0){
      system("cat /Noise/fasta/*.fa > /Noise/fasta/Coinfections_total.fa")
      system(paste("pangorunner.sh", "/Noise/fasta/Coinfections_total.fa"))
      df2<- read.csv("/Noise/fasta/Coinfections_total.fa_pango.csv")
      
      df2$Sample<-"S1"
      df2$Sample[grep("_S2",df2$taxon)]<-"S2"
      df2$ID<-gsub( "_S.\\.fa","",df2$taxon)
      
      
      S1.df<-df2[which(df$Sample=="S1"),c("lineage","conflict","note","ID")]
      S2.df<-df2[which(df$Sample=="S2"),c("lineage","conflict","note","ID")]
      
      colnames(S1.df)[-4]<-paste("S1_",colnames(S1.df)[-4], sep = "")
      colnames(S2.df)[-4]<-paste("S2_",colnames(S2.df)[-4], sep = "")
      
      total<-merge(S1.df, S2.df, by="ID")
      total<-total[,c("ID", "S1_lineage", "S2_lineage")]
      colnames(total)<-c("Sample", "Major", "Minor")
      
      if(length(which(total$Major==total$Minor))>0){
        total<-total[-which(total$Major==total$Minor),]
      }
      date<-gsub("-","",Sys.Date())
      if(nrow(total)>0) write.csv(total, paste("/Noise/Coinfection_Results_", date,".csv",sep=""), row.names = FALSE)
      if(nrow(total)==0) write.table("No Coinfections Found!",paste("/Noise/Coinfection_Results_", date,".csv",sep=""), row.names = FALSE, quote = FALSE, col.names = FALSE)
      try(system("rm -rf /Noise/fasta"))
    }else{
      date<-gsub("-","",Sys.Date())
      write.table("No Coinfections Found!",paste("/Noise/Coinfection_Results_", date,".csv",sep=""), row.names = FALSE, quote = FALSE, col.names = FALSE)
    }
}
    