library(genepopedit)
library(data.table)
library(ggplot2)
library(dplyr)
library(windowscanr)

#Script to be updated for new data Jan 2025

#Read in dxy and Pcadmix results
#consider outliers the top 5% across the genome
#Dxy -> lowest 5% indicate Canada and Barents-White Sea are similar (and more so than Canada-Southern Norway)

#PCAdmix -> top 5% indicates regions where high proportion of canadian ancestry in Finnmark



#Read Dxy data - Remove NAs if present
dxy=read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/AllDyx_windows100kbp_Chr1_29_StandardizeDxy_SNor_Finnmark.txt", header=T)
head(dxy)
sum(is.na(dxy$diff_dxy_mean))

#order dxy values
dxy_order=dxy[order(dxy$diff_dxy_mean),]
dxy_order[is.na(dxy_order$diff_dxy_n),]
dxy_order[is.na(dxy_order$diff_dxy_mean),]

(nrow(dxy_order))*0.05

#Plot outliers -
ggplot()+geom_point(data=dxy_order, aes(x=win_start, y=diff_dxy_mean))+
  geom_point(data= dxy_order[which(dxy_order$diff_dxy_mean< quantile(dxy_order$diff_dxy_mean, 0.05, na.rm=T)),], aes(x=win_start, y=diff_dxy_mean), col="red")+
  facet_wrap(.~Chr2, scales = "free_x")+theme_classic()

#Read PCAdmix data - remove NAs if present
pcadmix=read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/Results/All_PCAdmix_results_ProportionCanadian_jan2025.txt", header=T)
head(pcadmix)
nrow(pcadmix)
pcadmix[is.na(pcadmix$CanadaAncestry),]

sum(is.na(pcadmix$CanadaAncestry))

pcadmix_order=pcadmix[order(pcadmix$CanadaAncestry, decreasing = T),]
(nrow(pcadmix_order))*0.05

#Plot top 5% outliers
ggplot()+geom_smooth(data=pcadmix_order, aes(x=Start_Position, y=CanadaAncestry,  col="PCAdmix"),method="loess", pch=21, span=0.1)+
 # geom_point(data=pcadmix_order[which(pcadmix_order$CanadaAncestry> quantile(pcadmix_order$CanadaAncestry, 0.95, na.rm=T)),],
       #      aes(x=Start_Position, y=CanadaAncestry), col="red")+
 theme(legend.position = "top")+
  facet_wrap(.~Chr2, scales = "free_x")+theme_classic()+
  geom_smooth(data=dxy_order, aes(x=win_start, y=diff_dxy_mean, col="Diff dxy"),method="loess", pch=21, span=0.1)
 # geom_point(data= dxy_order[which(dxy_order$diff_dxy_mean< quantile(dxy_order$diff_dxy_mean, 0.05, na.rm=T)),], aes(x=win_start, y=diff_dxy_mean), col="blue")


#Oultiers (top 5%) from each analysis
padmix_outliers <-  pcadmix_order[which(pcadmix_order$CanadaAncestry >= quantile(pcadmix_order$CanadaAncestry, 0.95)),]

dxy_outliers <- dxy_order[which(dxy_order$diff_dxy_mean <= quantile(dxy_order$diff_dxy_mean, 0.05, na.rm=T)),]
nrow(dxy_outliers)

#Save outlier results
write.table(dxy_outliers, "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/Dxy_Outliers_Canada_finnmark_jan2025.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(padmix_outliers, "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/Pcadmix_Outliers_Canada_finnmark_jan2025.txt", quote = F, row.names = F, col.names = T, sep="\t")

#Find overlap

#To look for overlap - run window scanner with NAs for all points that aren't outliers
#This will provide the number of outliers in each 1Mbp Window

dxy_order2 = dxy_order

dxy_order2$diff_dxy_mean[which(dxy_order2$diff_dxy_mean > quantile(dxy_order2$diff_dxy_mean, 0.05, na.rm=T))] <- NA
View(dxy_order2)


#Inlcude extra data in case it gets removed due to missing positions at end of chromosome windows
#this does not influence results but includes an extra window at end of chromosome so that all results are read in - otherwise the windowscanner excludes them.

fake_data = list()
#add fake position at 1MBP for the max. --- for loop gets max value for each chromosome
for(i in 1:29){
  dat_EU=as.data.frame((max(dxy_order2$win_start[which(dxy_order2$Chr2==i)])))
  dat_EU$i <- i  # maybe you want to keep track of which iteration produced it?
  fake_data[[i]] <- dat_EU # add it to your list  }
}

fake_data_chr = do.call(rbind, fake_data)
colnames(fake_data_chr)=c("Max_Pos", "Chr")
fake_data_chr$Fake_Pos=fake_data_chr$Max_Pos+1e6

#Create tentative data frame to add to include extra windows at the end of chromosome ... Extra position have extra 1MBP added to the chromosome
#Make dataframe the same as the one you want to add to
head(dxy_order2)

fake_data_to_add_FST=as.data.frame(cbind("Chr2"=fake_data_chr$Chr, 
                                         "win_start"= fake_data_chr$Fake_Pos,
                                         "win_end"=rep(NA),
                                         "win_mid"=rep(NA),
                                         "diff_dxy_n"=rep(NA),
                                         "diff_dxy_mean"=rep(NA), "SNP"=rep(NA)
))


head(fake_data_to_add_FST)
head(dxy_order2)

dxy_order2_extra=rbind(as.data.frame(dxy_order2), fake_data_to_add_FST)

dxy_sliding=winScan(x = dxy_order2_extra, groups = "Chr2", position = "win_start",
                    win_size = 1000000, win_step = 500000,  values = "diff_dxy_mean", funs = "mean")

#Remove all windows that have 0 outliers
outlier_window_dxy_1MBP=dxy_sliding[which(dxy_sliding$diff_dxy_mean_n>0),]


#This will provide the number of outliers in each 1Mbp Window
##Pcadmix outliers
pcadmix_order2 <- pcadmix_order
pcadmix_order2[which(pcadmix_order2$Chr=="ssa08"),]

pcadmix_order2$CanadaAncestry[which(pcadmix_order2$CanadaAncestry < quantile(pcadmix_order2$CanadaAncestry, 0.95))] <- NA

pcadmix_order2[which(pcadmix_order2$Chr=="ssa08"),]


#Inlcude extra data in case it gets removed due to missing positions at end of chromosome windows
fake_data = list()
#add fake position at 1MBP for the max. --- for loop gets max value for each chromosome
for(i in 1:29){
  dat_EU=as.data.frame((max(pcadmix_order2$Start_Position[which(pcadmix_order2$Chr2==i)])))
  dat_EU$i <- i  # maybe you want to keep track of which iteration produced it?
  fake_data[[i]] <- dat_EU # add it to your list  }
}

fake_data_chr = do.call(rbind, fake_data)
colnames(fake_data_chr)=c("Max_Pos", "Chr")
fake_data_chr$Fake_Pos=fake_data_chr$Max_Pos+1e6

#Create fake data frame to add to include extra windows at the end of chromosome ... Fake position have extra 1MBP added to the chromosome
#Make dataframe the same as the one you want to add to
head(pcadmix_order2)




fake_data_to_add_FST=as.data.frame(cbind(   "CanadaAncestry"=rep(NA),
                             "Start_Position"= fake_data_chr$Fake_Pos,
                             "End_Position"= fake_data_chr$Fake_Pos,
                                          "Chr"=rep(NA),
                                          "Chr2"=fake_data_chr$Chr
                                          
))
head(fake_data_to_add_FST)
head(pcadmix_order2)

pcadmix_order2_extra=rbind(as.data.frame(pcadmix_order2), fake_data_to_add_FST)

#View(pcadmix_order2_extra)
Pcadmix_sliding=winScan(x = pcadmix_order2_extra, groups = "Chr2", position = "Start_Position",
                       win_size = 1000000, win_step = 500000,  values = "CanadaAncestry", funs = "mean")

#remove windows without outliers
outlier_window_pcadmix_1MBP=Pcadmix_sliding[which(Pcadmix_sliding$CanadaAncestry_n>0),]

#Now create "SNP" for each dataset (SNP indicates Chr and Start Position)
#Necessary to merge datasets together to look for overlap
outlier_window_pcadmix_1MBP$SNP=interaction(outlier_window_pcadmix_1MBP$Chr2, outlier_window_pcadmix_1MBP$win_start, sep="_")
View(outlier_window_pcadmix_1MBP)
outlier_window_dxy_1MBP$SNP=interaction(outlier_window_dxy_1MBP$Chr2, outlier_window_dxy_1MBP$win_start, sep="_")
View(outlier_window_dxy_1MBP)


#Overlap in outliers between ALL DATASETS
all_overlap_1MBP=merge(outlier_window_pcadmix_1MBP, outlier_window_dxy_1MBP, by="SNP")
nrow(all_overlap_1MBP)


#Results of overlapping outlier windows
write.table(all_overlap_1MBP,
            "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Compare_Introgression_Outliers/Overlap_2metrics_introgression_jan2025.txt",
         quote = F, row.names = F, col.names = T, sep="\t")

#Data saved above was manually evaluated to determine number of blocks. For each window a block # was added, and then these data were used to summarize outlier regions:
block_counts<-read.csv("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Compare_Introgression_Outliers/Overlap_2metrics_introgression_Jan2025_count_blocks_in_ORDER.csv", header=T)
head(block_counts)

#Summary table with size of outlier blocks and positions
#Note that the number out outliers within the blocks for each metric represents sum of outliers detected for sliding window approach - therefore this should be adjusted for duplicates 
summary_table<- block_counts %>%
  group_by(Genomic_block_in_order) %>%
  summarise(Chr= min(Chr), Start = min(Window.start.position), End = max(Window.end.position), Pcadmix=sum(Pcadmix.outlier),
            Dyx=sum(Dxy.Outlier), Block_size = ((max(Window.end.position)-min(Window.start.position))+1) )

#Saved data used for Supplemental figure to show overlapping outlier blocks between the two metrics
write.table(summary_table, "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Compare_Introgression_Outliers/Summary_of_blocks_Overlap_2metrics_introgression_Jan2025.txt", quote = F, row.names = F,
            col.names = T, sep="\t")
