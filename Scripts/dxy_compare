#Libraries
library(ggplot2)
library(qqman)
library(windowscanr)

#Working directory
setwd("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/")

data_file<-NULL
for(i in 1:9){
  fin=read.table(paste0("Ssa0",i,"_dxy_finnmark_canada.txt"), header=T)
  data_file=rbind(data_file, fin)
  rm(fin)
}
for(i in 10:29){
  fin=read.table(paste0("Ssa",i,"_dxy_finnmark_canada.txt"), header=T)
  data_file=rbind(data_file, fin)
  rm(fin)
}
#Generate Z score
data_file$StandardDxy_Fin=scale(data_file$dxy, center = TRUE, scale = TRUE)
hist(data_file$StandardDxy_Fin)

data_file2<-NULL
for(i in 1:9){
  fin=read.table(paste0("Ssa0",i,"_dxy_southernnorway_canada.txt"), header=T)
  data_file2=rbind(data_file2, fin)
  rm(fin)
}

for(i in 10:29){
  fin=read.table(paste0("Ssa",i,"_dxy_southernnorway_canada.txt"), header=T)
  data_file2=rbind(data_file2, fin)
  rm(fin)
}

#Generate Z scores
data_file2$StandardDxy_Snor=scale(data_file2$dxy, center = TRUE, scale = TRUE)
head(data_file2)
hist(data_file2$dxy)

#all data combine
all_dxy_ssa1_29=as.data.frame(cbind(data_file2$StandardDxy_Snor,data_file2$dxy, data_file))
head(all_dxy_ssa1_29)
colnames(all_dxy_ssa1_29)=c("StanDxy_Snor", "dxy_southnor", "dxy_fin", "Position", "Chr", "StanDxy_Fin")

#Calculate difference in standardized dxy
all_dxy_ssa1_29$diff_dxy=all_dxy_ssa1_29$StanDxy_Fin-all_dxy_ssa1_29$StanDxy_Snor
#Values that are low indicate that Finnmark is more similar to Canada in this region (compared to S Norway and Canada)

#Add extra identifier (chr_Position)
all_dxy_ssa1_29$SNP=interaction(all_dxy_ssa1_29$Chr, all_dxy_ssa1_29$Position, sep="_")
head(all_dxy_ssa1_29)

#Save file -- Z scores used now -- note
write.table(all_dxy_ssa1_29, "Chr1_29_StandardizeDxy_SNor_Finnmark.txt", quote=F, row.names = F, col.names = T, sep="\t")

#Top 
order=all_dxy_ssa1_29[order(all_dxy_ssa1_29$diff_dxy, decreasing = F),]
order[1:10,]

all_dxy_ssa1_29$Chr2=as.numeric(all_dxy_ssa1_29$Chr)

min(all_dxy_ssa1_29$diff_dxy)
max(all_dxy_ssa1_29$diff_dxy)

#plot values
manhattan(all_dxy_ssa1_29, chr = "Chr2", bp = "Position", snp = "SNP", p = "diff_dxy", logp = F,
          ylim=c(min(all_dxy_ssa1_29$diff_dxy)+0.1,max(all_dxy_ssa1_29$diff_dxy)+0.1))


#sliding window
  
#Add fake chromosome and do window scan 
head(all_dxy_ssa1_29)

fake_data = list()

#add fake position at 1MBP for the max. --- for loop gets max value for each chromosome
for(i in 1:9){
  dat_EU=as.data.frame((max(all_dxy_ssa1_29$Position[which(all_dxy_ssa1_29$Chr==paste0("ssa0",i))])))
  dat_EU$i <- i  # maybe you want to keep track of which iteration produced it?
  fake_data[[i]] <- dat_EU # add it to your list  }
}

fake_data_chr = do.call(rbind, fake_data)

fake_data2 = list()

#add fake position at 1MBP for the max. --- for loop gets max value for each chromosome
for(i in 10:29){
  dat_EU=as.data.frame((max(all_dxy_ssa1_29$Position[which(all_dxy_ssa1_29$Chr==paste0("ssa",i))])))
  dat_EU$i <- i  # maybe you want to keep track of which iteration produced it?
  fake_data2[[i]] <- dat_EU # add it to your list  }
}

fake_data_chr2 = do.call(rbind, fake_data2)
colnames(fake_data_chr2)=c("Max_Pos", "Chr")
colnames(fake_data_chr)=c("Max_Pos", "Chr")


fake_data2_for_file=rbind(fake_data_chr,fake_data_chr2)

colnames(fake_data2_for_file)=c("Max_Pos", "Chr")

fake_data2_for_file$Fake_Pos=fake_data2_for_file$Max_Pos+1e6

#Create fake data frame to add to include extra windows at the end of chromosome ... Fake position have extra 1MBP added to the chromosome
#Make dataframe the same as the one you want to add to
head(all_dxy_ssa1_29)

fake_data_to_add_dxy=as.data.frame(cbind("StanDxy_Snor"=rep(NA),
                                         "dxy_southnor"=rep(NA),
                                         "dxy_fin"=rep(NA),
                                         "Position"= fake_data2_for_file$Fake_Pos,
                                         "Chr"=paste0("ssa", fake_data2_for_file$Chr), 
                                         "StanDxy_Fin"=rep(NA),
                                         "diff_dxy"=rep(NA),
                                         "SNP"=paste0("AX-fake",1:29),
                                         "Chr2"=c(1:29) 
))

#remove Chr that have NA from original dataframe 
all_dxy_ssa1_29$Chr2=as.numeric(all_dxy_ssa1_29$Chr)

#rbind both datasets (fake and real)
all_data_plus_extra_dxy=rbind(as.data.frame(all_dxy_ssa1_29), fake_data_to_add_dxy)

tail(all_data_plus_extra_dxy)

#Change anything to numeric if needed
all_data_plus_extra_dxy$Position=as.numeric(as.character(all_data_plus_extra_dxy$Position))
all_data_plus_extra_dxy$diff_dxy=as.numeric(as.character(all_data_plus_extra_dxy$diff_dxy))

head(all_data_plus_extra_dxy)
all_data_plus_extra_dxy$Chr2=as.factor(all_data_plus_extra_dxy$Chr2)

#Run window scaner to get mean dxy per window (100 Kbp windows)
dxy_sliding=winScan(x = all_data_plus_extra_dxy, groups = "Chr2", position = "Position",
                    win_size = 100000, win_step = 50000,  values = "diff_dxy", funs = "mean")

#make sure chromosme is numeric for plotting
dxy_sliding$Chr2=as.numeric(as.character(dxy_sliding$Chr2))

#Get lower 5% quantile for outliers
quantile(dxy_sliding$diff_dxy_mean, 0.05, na.rm=T)
min_y <-min(dxy_sliding$diff_dxy_mean, na.rm=T)-0.1
max_y <-max(dxy_sliding$diff_dxy_mean, na.rm=T)+0.1
  
dxy_sliding$SNP=interaction(dxy_sliding$Chr2, dxy_sliding$win_start, sep="_")
manhattan(dxy_sliding[!is.na(dxy_sliding$diff_dxy_mean),], chr = "Chr2", bp = "win_start", snp = "SNP",
          p = "diff_dxy_mean", logp = F, ylim=c(min_y, max_y), ylab="Dxy finnmark:canada - Dxy S.Norway:canada",
          genomewideline =  quantile(dxy_sliding$diff_dxy_mean, 0.05, na.rm=T)  )

outliers_windows_dxy<- (dxy_sliding[which(dxy_sliding$diff_dxy_mean <= quantile(dxy_sliding$diff_dxy_mean, 0.05, na.rm=T)),])
barplot(table(outliers_windows_dxy$Chr2))

#Save outlier results
write.table(dxy_sliding, "AllDyx_windows100kbp_Chr1_29_StandardizeDxy_SNor_Finnmark.txt", quote=F, row.names = F, col.names = T, sep="\t")

write.table(outliers_windows_dxy, "DxyOutliers_windows_Chr1_29_StandardizeDxy_SNor_Finnmark.txt", quote=F, row.names = F, col.names = T, sep="\t")

manhattan(dxy_sliding[which(dxy_sliding$Chr2=="23"),], chr = "Chr2", bp = "win_start", snp = "SNP",
          p = "diff_dxy_mean", logp = F, ylim=c(min_y, max_y), ylab="Dxy finnmark:canada - Dxy S.Norway:canada",
          genomewideline =  quantile(dxy_sliding$diff_dxy_mean, 0.05, na.rm=T)  )

