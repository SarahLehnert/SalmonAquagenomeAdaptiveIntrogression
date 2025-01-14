
#library 
library(qqman)

#Set directory
setwd("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/Results/")


#######################################################################################
#This script will save all the results for each chromosome together to use for analyses

data_table_to_save=NULL
for(i in 1:9){
#read in data (0s and 1s for ancestry assignment)
  data1_conf=read.table(paste0("ssa0",i,"_WGS_FinnmarkIntrogress_jan_2025.vit.txt"), header=F)
#remove first row of individual names  
  data_chr1_toplot=data1_conf[,-1]
  data_chr1_toplot=as.matrix(data_chr1_toplot)
#Get marker positions/windows  
  windows=read.table(paste0("ssa0",i,"_WGS_FinnmarkIntrogress_jan_2025.markers.txt"), header=F, fill=T)
#start position  
  start_poistion=windows$V2
  windows$Start_Position=gsub(start_poistion, pattern = paste0("ssa0",i,":"), replacement = "")
#end position  
  end_position=windows$V101
  windows$End_Position=gsub(end_position, pattern = paste0("ssa0",i,":"), replacement = "")
#Combine results for chromosme  
  data_chr1=as.data.frame(cbind(colMeans(data_chr1_toplot, na.rm = T), windows$Start_Position, windows$End_Position))

#ensure values are numeric
  data_chr1$V1=as.numeric(as.character(data_chr1$V1))
  data_chr1$V2=as.numeric(as.character(data_chr1$V2))
  data_chr1$V3=as.numeric(as.character(data_chr1$V3))
  
#combine data and add chr name  
  data_table_to_save=rbind(data_table_to_save, cbind(data_chr1, rep(paste0("ssa0",i))))
  
}

#Combined results for Ssa01 to Ssa09
Results_Ssa01_09=as.data.frame(data_table_to_save)
colnames(Results_Ssa01_09)=c("CanadaAncestry", "Start_Position","End_Position", "Chr")

#Repeat same script as above to get results for Ssa10 to Ssa29
data_table_to_save2=NULL
for(i in 10:29){
  
  data1_conf=read.table(paste0("ssa",i,"_WGS_FinnmarkIntrogress_jan_2025.vit.txt"), header=F)
  
  data_chr1_toplot=data1_conf[,-1]
  data_chr1_toplot=as.matrix(data_chr1_toplot)
  
  windows=read.table(paste0("ssa",i,"_WGS_FinnmarkIntrogress_jan_2025.markers.txt"), header=F, fill=T)
  start_poistion=windows$V2
  windows$Start_Position=gsub(start_poistion, pattern = paste0("ssa",i,":"), replacement = "")
  
  end_position=windows$V101
  windows$End_Position=gsub(end_position, pattern = paste0("ssa",i,":"), replacement = "")
  
  data_chr1=as.data.frame(cbind(colMeans(data_chr1_toplot, na.rm = T), windows$Start_Position, windows$End_Position))
  
  data_chr1$V1=as.numeric(as.character(data_chr1$V1))
  data_chr1$V2=as.numeric(as.character(data_chr1$V2))
  data_chr1$V3=as.numeric(as.character(data_chr1$V3))
  
  data_table_to_save2=rbind(data_table_to_save2, cbind(data_chr1, rep(paste0("ssa",i))))
  
}

Results_Ssa010_29=as.data.frame(data_table_to_save2)
colnames(Results_Ssa010_29)=c("CanadaAncestry", "Start_Position","End_Position", "Chr")

#combine results
all_PCAdmix_results=as.data.frame(rbind(Results_Ssa01_09,Results_Ssa010_29))

all_PCAdmix_results$Chr2=as.numeric(all_PCAdmix_results$Chr)

#Save all results together
write.table(all_PCAdmix_results, "All_PCAdmix_results_ProportionCanadian_jan2025.txt", quote = F, row.names = F, col.names = T, sep="\t")



#Plot manhattan plot for comparisons
all_PCAdmix_results$SNP=interaction(all_PCAdmix_results$Chr, all_PCAdmix_results$Start_Position, sep = "_")

##Selecting outlier regions

#Top 5 % (95%) were chosen as outlier regions
quantile(all_PCAdmix_results$CanadaAncestry, 0.95)

#Save outliers in seperate file
write.table(all_PCAdmix_results[which(all_PCAdmix_results$CanadaAncestry >= quantile(all_PCAdmix_results$CanadaAncestry, 0.95) ),], "Outliers_PCAdmix_results_ProportionCanadian_jan2025.txt", quote = F, row.names = F, col.names = T, sep="\t")

#Max position
all_PCAdmix_results[which.max(all_PCAdmix_results$CanadaAncestry),]

#Check results with quick manhattan plot
manhattan(x=all_PCAdmix_results, 
          snp = "SNP", chr = "Chr2",bp = "Start_Position", 
          p="CanadaAncestry", ylim=c(0, 0.8), logp = F, 
          ylab="Proportion Canada Ancestry" ,
          highlight =all_PCAdmix_results$SNP[which(all_PCAdmix_results$CanadaAncestry >= quantile(all_PCAdmix_results$CanadaAncestry, 0.95) )] )
