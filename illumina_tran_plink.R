illumina_plink_tran<-function(final_report_name,map_name,out_name,coding_type,callrate_threshold){
  if(!require(data.table)) install.packages("data.table")
  if(!require(ggplot2)) install.packages("ggplot2")
  list.files(pattern = 'DNAReport.csv')
  DNA_report_name <- gsub(pattern = 'FinalReport.txt','DNAReport.csv',final_report_name)
  DNA_report <- read.csv(DNA_report_name,skip = 2)
  names(DNA_report)
  DNA_report_fail <- DNA_report[which(DNA_report$Call_Rate<callrate_threshold),]
  ggplot(data = DNA_report,aes(x=Call_Rate))+
    geom_density(fill='gold3',color='gold3')+
    geom_vline(xintercept =callrate_threshold,color='red')+
    labs(title = paste0('Total samples:',nrow(DNA_report),'; ','Samples lower than threshold:',nrow(DNA_report_fail)))+
    theme_light()
  ggsave('CallRate_distribution.pdf')
  write.csv(DNA_report_fail,gsub(pattern = 'FinalReport.txt','DNAReport_failed.csv',final_report_name),row.names = F)
  final_report<-fread(final_report_name,header = F,skip = 10)
  final_report<-as.data.frame(final_report)
  names(final_report)<-c("SNP_Name","Sample_ID","Allele1_Forward","Allele2_Forward","Allele1_Top","Allele2_Top","Allele1_AB","Allele2_AB","GC_Score","X","Y")
  final_report$Sample_ID<-as.character(final_report$Sample_ID)
  id_freq<-as.data.frame(table(final_report$Sample_ID))
  id_freq<-id_freq[which(id_freq$Freq>0),]
  if(max(id_freq$Freq)!= min(id_freq$Freq)){
    print("error in snp of sample")}else{
      sample_id<-sort(as.character(id_freq$Var1))
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
  map<-fread(map_name,fill = T)
  map<-map[,2:4]
  names(map)<-c("snp","chr","position")
  map<-merge(snp_list,map,by="snp",sort = F)
  map$mendel<-0
  #allele_matrix<-as.data.frame(allele_matrix)
  fwrite(allele_matrix,paste(out_name,".ped",sep = ""),quote = F,col.names = F,row.names = F,sep = " ")
  fwrite(map[,c(2,1,4,3)],paste(out_name,".map",sep = ""),quote = F,col.names = F,row.names = F,sep = " ")
  print("done!")
}
