setwd("~/Desktop/Sarah/Salmon/WGS_Aquagenome/SweeD/Results/finnmark/")

#Combine results from SweeD reports into single file

file_list <- list.files()
rm(dataset)

dataset<-NULL

for (file in file_list){
  
  # if the merged dataset does exist, append to it
  
  temp_dataset <-read.table(file, header=T, sep="\t", fill=T, skip=2)
  temp_dataset$Chromosome = rep(file)  
  
  
  dataset<-rbind(dataset, temp_dataset)
  rm(temp_dataset)
  
}


dataset_final <- as.data.frame(sapply(dataset, gsub, pattern = "SweeD_Report.ssa", replacement=""), stringsAsFactors=FALSE)
dataset_final <- as.data.frame(sapply(dataset_final, gsub, pattern = "_wg_sweed_finnmark", replacement=""), stringsAsFactors=FALSE)

dataset_final$Chr2=as.numeric(dataset_final$Chromosome)
tail(dataset_final)
write.table(dataset_final, "Finnmark_Results_SweeD.txt", quote = F, row.names = F, col.names = T, sep="\t")

