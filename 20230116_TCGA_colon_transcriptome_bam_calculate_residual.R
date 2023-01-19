# This script is to calculate miRNA biogenesis efficiency (residual) and draw plot about correlation between expression levels of miRNA and transcirpt
# 2023/1/16 made

# activate packages
library(stringr)
library(ggplot2)

# function for caluculating outlier
cal.outlier <-function(x){
  q <-as.numeric(quantile(x))
  iqr <-IQR(x)
  outlier1 <-q[2]-iqr*1.5
  outlier2 <-q[4]+iqr*1.5
  outliers <-append(outlier1,outlier2)
  return(outliers)
}

# function for calculating minimum value
cal.min <-function(x){
  for (i in 1:ncol(x)) {
    x[,i] <-ifelse(x[,i]==0,1,x[,i])
  }
  m.table <-apply(x, 1,min)
  m <-min(m.table)
  return(m)
}

# import correspondence table between TCGA colon transcriptome bam and TCGA colon miRNA quantification file
# this table is located at "https://github.com/Ryosuke-Hirota/20221124_make_table_of_TCGA_transcriptome_bam_correspond_to_TCGA_miRNA_expression_data"
setwd("C:/Rdata/20221124_make_table_of_TCGA_transcriptome_bam_correspond_to_TCGA_miRNA_expression_data")
cor.table <-read.table("correspondence_table_between_TCGA_colon_transcriptome_bam_and_miRNA_qunatification_file.txt",sep="\t",header = T,stringsAsFactors = F)

# import list of transcripts that intersect with miRNAs in gencode v36
# this list is located at "https://github.com/Ryosuke-Hirota/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA"
setwd("C:/Rdata")
primir.list <-read.table("TCGA_hg38_transcript_intersect_with_miRNA.txt",sep="\t",header = F,stringsAsFactors = F)
primir.list[,2] <-primir.list[,2]-1
primir.list <-primir.list[primir.list[,2]!=primir.list[,8]&primir.list[,3]!=primir.list[,9],]

# remove duplicated gene names and make list regarding miRNA name and gene name
mir_transcipt <-primir.list[,c(4,11)]
mir_transcipt <-subset(mir_transcipt,!duplicated(mir_transcipt))
mir_transcipt <-mir_transcipt[order(mir_transcipt[,1]),]
primir.list <-primir.list[,c(4,11,10)]
transcripts <-unique(primir.list[,3])

# list TCGA colon transcript quantification files
# these files are located at "\\fsw-q02\okamura-lab\20221006_TCGA_colon_salmon_quant_transcriptome"
setwd("C:/Rdata/20221006_TCGA_colon_salmon_transcriptome_quantification")
transcript.quant <-list.files(path = "C:/Rdata/20221006_TCGA_colon_salmon_transcriptome_quantification",pattern = ".txt")
t.file.id <-gsub("_quant.txt","",transcript.quant)

# make table summarized TCGA colon transcript quantification
for (i in 1:length(transcript.quant)){
  transcript.file <-read.table(transcript.quant[i],sep = "\t",header = T,stringsAsFactors = F)
  t <-match(transcripts,transcript.file[,1])
  transcript.file <-transcript.file[t,]
  transcript.file <-transcript.file[,c(1,4)]
  colnames(transcript.file)[2] <-transcript.quant[i]
  if(i==1){
    transcript.quant.table <-transcript.file
  }else{
    transcript.quant.table <-merge(transcript.quant.table,transcript.file,by="Name")
  }}

# import miRNA quantification table
# this table is located at ""
setwd("C:/Rdata/20230105_TCGA_colon_miRNA_quantification")
miRNA.quant.table <-read.table("table_of_TCGA_colon_miRNA_quantifications.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# calculate minimum value of transcript expression level
transcript.table <-transcript.quant.table[,-1]
min.transcript <-cal.min(transcript.table)

# calculate minimum value of miRNA expression level
mirna.table <-miRNA.quant.table[,-1]
min.miRNA <-cal.min(mirna.table)

# make new directory
setwd("C:/Rdata")
dir.create("20230116_TCGA_colon_transcriptome_bam_calculate_residual")
setwd("C:/Rdata/20230116_TCGA_colon_transcriptome_bam_calculate_residual")
dir.create("plot")
dir.create("table")

# investigate correlation between transcript and miRNA, draw plot and make summary
for (i in 1:nrow(mir_transcipt)){
  # extract expression level of a certain miRNA 
  miRNA.df <-miRNA.quant.table[miRNA.quant.table[,1]==mir_transcipt[i,1],]
  miRNA.df <-as.data.frame(t(miRNA.df),stringsAsFactors = F)
  m.cor <-match(rownames(miRNA.df),cor.table[,5])
  miRNA.df[,2] <-rownames(miRNA.df)
  miRNA.df[,3] <-cor.table[m.cor,2]
  colnames(miRNA.df) <-c(mir_transcipt[i,1],"miRNA_file_name","transcript_file_name")
  rownames(miRNA.df) <-NULL
  miRNA.df <-miRNA.df[-1,]
  
  # extract expression level of a certain transcript
  transcript <-primir.list[primir.list[,1]==mir_transcipt[i,1]&primir.list[,2]==mir_transcipt[i,2],3]
  t <-match(transcript,transcript.quant.table[,1])
  transcript.df <-transcript.quant.table[t,]
  transcript.df <-transcript.df[,-1]
  transcript.df[1,] <-apply(transcript.df, 2, sum)
  transcript.df <-transcript.df[-2,] 
  transcript.df <-as.data.frame(t(transcript.df),stringsAsFactors = F)
  t.cor <-match(t.file.id,cor.table[,1])
  transcript.df[,2] <-cor.table[t.cor,2]
  colnames(transcript.df) <-c(mir_transcipt[i,2],"transcript_file_name")
  rownames(transcript.df) <-NULL
  
  # merge expression levels of miRNA and transcript
  mt.df <-merge(miRNA.df,transcript.df,by="transcript_file_name")
  mt.df <-subset(mt.df,!is.na(mt.df[,1]))
  mt.df <-mt.df[,c(4,2,1,3)]
  mt.df[,2] <-as.numeric(mt.df[,2])
  
  # remove zero expression
  ex.mt.df <-mt.df[mt.df[,1]!=0&mt.df[,2]!=0,]
  
  # normalize with log2 (after excluding cell lines without expression) 
  ex.mt.df[,1] <-log2(ex.mt.df[,1])
  ex.mt.df[,2] <-log2(ex.mt.df[,2])
  
  # calculate outliers
  t.outlier <-cal.outlier(ex.mt.df[,1])
  m.outlier <-cal.outlier(ex.mt.df[,2])
  
  # remove outliers
  ex.mt.df <-ex.mt.df[ex.mt.df[,1]>t.outlier[1]&ex.mt.df[,1]<t.outlier[2],]
  ex.mt.df <-ex.mt.df[ex.mt.df[,2]>m.outlier[1]&ex.mt.df[,2]<m.outlier[2],]
  
  # make empty summary
  sm <-as.data.frame(matrix(nrow = nrow(ex.mt.df),ncol = 3))
  colnames(sm) <-c("residual","transcriptome_bam","miRNA_quantification")
  
  # calculate correlation coefficient
  r <-try(cor.test(ex.mt.df[,1],ex.mt.df[,2],method="pearson"),silent = T)
  
  if(class(r)!="try-error"){
    # perform linear regression analysis and calculate residual
    l <-lm(ex.mt.df[,2]~ex.mt.df[,1],data = ex.mt.df)
    res <-residuals(l)
    
    # draw plot about correlation between expression levels of transcript and miRNA
    p <-ggplot()+
      geom_point(data=ex.mt.df,aes(x=ex.mt.df[,1],y=ex.mt.df[,2]))+
      geom_smooth(data=ex.mt.df,mapping = aes(x=ex.mt.df[,1],y=ex.mt.df[,2]),method="lm",formula='y~x',se=FALSE,colour="black",size=0.5)+
      labs(title=paste0("R =",signif(r$estimate,3),", p = ",signif(r$p.value,3),", n = ",nrow(ex.mt.df)),x=mir_transcipt[i,2],
           y=mir_transcipt[i,1])+ 
      theme_bw()+
      theme(legend.background = element_rect(fill = "white", colour = "black"))
    
    setwd("C:/Rdata/20230116_TCGA_colon_transcriptome_bam_calculate_residual/plot")
    ggsave(filename=paste0("plot_of_correlation_between_",mir_transcipt[i,1],"_and_",mir_transcipt[i,2],".pdf"),plot = p)
    
    # write summary
    sm[,1] <-res
    sm[,2] <-ex.mt.df[,3]
    sm[,3] <-ex.mt.df[,4]
  }
  # output summary
  setwd("C:/Rdata/20230116_TCGA_colon_transcriptome_bam_calculate_residual/table")
  write.table(sm,paste0("table_of_residual_about_",mir_transcipt[i,1],"_vs_",mir_transcipt[i,2],".txt"),sep="\t",row.names = F,quote = F)
}

