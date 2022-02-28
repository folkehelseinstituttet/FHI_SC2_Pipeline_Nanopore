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

# args=commandArgs(TRUE)
#  cores<-10
# if(length(args)==0){ 
#   cores<-10
# }else{
#   cores<-as.numeric(gsub("c","",args[1]))
#   }

cores<-10
# results <-"/home/docker/Fastq/"   #Docker
# results<-list.dirs(results)
# results<-results[grep("summaries$",results)]
# results<-paste(results,"/",sep = "")
# 
# common.folder <-"/home/docker/CommonFiles/"
# 
# 
results<-"/media/nacho/Data/contamintant_detection/Nanopore/"
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
  
  dir.create(temp)
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  out.par<-foreach(i=1:length(bamfiles), .verbose=FALSE, .options.snow = opts) %dopar%{
    
    try(system(paste("FINex -f ",bamfiles[i], " > ",
                     gsub(".*/",temp,gsub("\\.bam", "_NoisExtractorResult.tsv",bamfiles[i])), sep = "")))
    
  }
  stopCluster(cluster.cores)
  
  

# Plotting ----------------------------------------------------------------

  
  results.files<-list.files(temp, pattern = "_NoisExtractorResult\\.tsv$", full.names = TRUE)

  samp<-c(1:length(results.files))
  
  pb <- progress_bar$new(
    format = "Index: :samp.pb [:bar] :elapsed | eta: :eta",
    total = length(results.files),    # 100 
    width = 60)
  
  
  progress <- function(n){
    pb$tick(tokens = list(samp.pb = samp[n]))
  } 
  
  opts <- list(progress = progress)
  
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  stable.pos<-read.csv("/media/nacho/Data/contamintant_detection/stable_positions_feb2022.csv")
  stable.pos<-stable.pos$x
  
  out.par<-foreach(i=1:length(results.files), .verbose=FALSE, .options.snow = opts, .packages = c("ggplot2")) %dopar%{
    p.noise<-function(x){
      return(length(which(df$Noise[which(df$Base %in% stable.pos)]> x))/length(stable.pos))
    }
    

    
    df<-read.csv(results.files[i],sep = "\t",header = FALSE)
    
    df<-df[,c(1,2,3,6,7)]
    
    colnames(df)<-c("Base","Noise","Reads","FreqMinor","NoiseMinor")
    df$NonMR<-df$Reads*df$Noise
    df$MR<-df$Reads-df$NonMR
    
    stable.pos<-stable.pos[which(stable.pos %in% df$Base)]
    
    df$p.value<-unlist(lapply(df$Noise, p.noise)) 
    df$p.adj<-p.adjust(df$p.value, method = "fdr")  
    
    genome.position<-as.data.frame(c(1:29903))
    colnames(genome.position)<-"Base"
    df<-merge(genome.position, df, by="Base", all=TRUE)
    df$Outlier<-"NO"
    df$Outlier[which(df$NoiseNP>0.95)]<-"YES"
    
    df$Reads[which(is.na(df$Reads))]<-0
    df$NoiseNP<-1-df$p.adj
    df$NoiseNP[which(is.na(df$NoiseNP))]<-0
    names<-gsub("\\.sorted.*","",gsub("_S[0-9].*","",gsub("R[0-9].*","",gsub(".*/","",results.files[i]))))
    p<-ggplot(df)+
      geom_line(aes(Base, NoiseNP))+
      geom_point(data=subset(df, Outlier=="YES"),aes(Base, NoiseNP),col="red", alpha=0.3)+
      geom_point(data=subset(df, Reads<20),aes(Base, 0),col="blue", alpha=0.1)+
      ylim(0,1)+
      theme_minimal()+
      ggtitle(names)
    return(list(p,df))
  }
  
 stopCluster(cluster.cores)
  
  out.plots<-out.par
  
 
date<-gsub("-","",Sys.Date())
  
  if(length(out.plots)<=40 ){
    if(length(out.plots)>0){
      ggarrange(plotlist =  out.plots[1:length(out.plots)], ncol = 4, nrow = 10)
      ggsave(paste(results,"ResultsNoisExtractor_",date,"_",".pdf", sep=""), width = 420, height = 600, units = "mm") #A4
    }  
  }else{
    plotting<-TRUE
    start<-1
    end<-40
    counter<-0
    
    while(plotting){
      if(end==length(out.plots)) plotting<-FALSE
      ggarrange(plotlist =  out.plots[start:end], ncol = 4, nrow = 10)
      start<-end+1
      end<-end+40
      if(end>=length(out.plots)) end<-length(out.plots)
      counter<-counter+1
      
      ggsave(paste(results,"ResultsNoisExtractor_",date,"_","_",counter,".pdf", sep=""), width = 420, height = 600, units = "mm") #A4
      
    }
  }

 pdf.list<-list.files(results, full.names = TRUE, pattern = ".*NoisExtractor.*\\.pdf")
 if(length(pdf.list)>1){
 pdf_combine(pdf.list, output = gsub("_.\\.pdf","_Merged.pdf",pdf.list[1]))
 file.remove(pdf.list)}
  }
  
