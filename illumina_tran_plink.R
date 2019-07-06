illumina_plink_tran<-function(final_report_name,map_name,out_name,coding_type){
  if(!require(readr)) install.packages("readr")
  if(!require(tidyverse)) install.packages("tidyverse")
  final_report<-read_tsv(final_report_name,col_names = F,skip = 10)
  final_report<-as.data.frame(final_report)
  names(final_report)<-c("SNP_Name","Sample_ID","Allele1_Forward","Allele2_Forward","Allele1_Top","Allele2_Top","Allele1_AB","Allele2_AB","GC_Score","X","Y")
  id_freq<-as.data.frame(table(final_report$Sample_ID))
  id_freq<-id_freq[which(id_freq$Freq>0),]
  if(max(id_freq$Freq)!= min(id_freq$Freq)){
    print("error in snp of sample")}else{
      sample_id<-sort(as.numeric(as.character(id_freq$Var1)))
      sample_n<-length(sample_id)
      allele_matrix<-matrix(nrow = sample_n,ncol = (6+2*max(id_freq$Freq)))
      odd_n<-seq(7,(5+2*max(id_freq$Freq)),by=2)
      evev_n<-seq(8,(6+2*max(id_freq$Freq)),by=2)
      if(coding_type=="Forward"){
        for(i in 1:sample_n){
          allele_matrix[i,2]<-sample_id[i]
          print(paste("====>individual",sample_id[i],"is transforming<===="))
          sample_position<-which(final_report$Sample_ID==sample_id[i])
          allele_matrix[i,odd_n]<-final_report[sample_position,3]
          allele_matrix[i,evev_n]<-final_report[sample_position,4]
        }
      }else if(coding_type=="Top"){
        for(i in 1:sample_n){
          allele_matrix[i,2]<-sample_id[i]
          print(paste("====>individual",sample_id[i],"is transforming<===="))
          sample_position<-which(final_report$Sample_ID==sample_id[i])
          allele_matrix[i,odd_n]<-final_report[sample_position,5]
          allele_matrix[i,evev_n]<-final_report[sample_position,6]
        }
      }else{
        for(i in 1:sample_n){
          allele_matrix[i,2]<-sample_id[i]
          print(paste("====>individual",sample_id[i],"is transforming<===="))
          sample_position<-which(final_report$Sample_ID==sample_id[i])
          allele_matrix[i,odd_n]<-final_report[sample_position,7]
          allele_matrix[i,evev_n]<-final_report[sample_position,8]
        }
      }
      print(paste("number of individual",sample_n,"done",sep=" "))
    }
  allele_matrix[allele_matrix=="-"]<-"0"
  allele_matrix[,1]<-"TEST"
  allele_matrix[,3:5]<-0
  allele_matrix[,6]<--9
  snp_list<-as.data.frame(final_report[sample_position,1])
  names(snp_list)<-"snp"
  map<-read.table("SNP_Map.txt",header = F,col.names = c(paste("V",c(1:9),sep = "")),fill = T)
  map<-map[-c(1,2),2:4]
  names(map)<-c("snp","chr","position")
  map<-merge(snp_list,map,by="snp",sort = F)
  map[,4]<-0
  write.table(allele_matrix,paste(out_name,".ped",sep = ""),quote = F,col.names = F,row.names = F)
  write.table(map[,c(2,1,4,3)],paste(out_name,".map",sep = ""),quote = F,col.names = F,row.names = F)
  print("done!")
}
