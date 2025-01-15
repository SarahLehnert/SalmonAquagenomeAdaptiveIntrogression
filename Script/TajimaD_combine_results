
#Rscript to combine results for Tajima D for Finnmark populations

setwd("~/Desktop/Sarah/Salmon/WGS_Aquagenome/TajimaD/Results/finnmark/")
file_list <- list.files()
rm(dataset)

dataset<-NULL

for (file in file_list){
  
  # if the merged dataset does exist, append to it
  
  temp_dataset <-read.table(file, header=T, sep="\t", fill=T)
  temp_dataset$Chromosome = rep(file)  
  
  
  dataset<-rbind(dataset, temp_dataset)
  rm(temp_dataset)
  
}


dataset_final <- as.data.frame(sapply(dataset, gsub, pattern = "_wg.phased_finnmark_TajimaD.Tajima.D", replacement=""), stringsAsFactors=FALSE)
dataset_final <- as.data.frame(sapply(dataset_final, gsub, pattern = "ssa", replacement=""), stringsAsFactors=FALSE)

dataset_final$Chr2=as.numeric(dataset_final$Chromosome)
tail(dataset_final)
write.table(dataset_final, "~/Desktop/Sarah/Salmon/WGS_Aquagenome/TajimaD/Finnmark_Results_TajimaD.txt", quote = F, row.names = F, col.names = T, sep="\t")

