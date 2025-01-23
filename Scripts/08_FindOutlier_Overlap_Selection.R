library(data.table)
library(ggplot2)
library(dplyr)
library(windowscanr)

setwd("~/Desktop/Sarah/Salmon/WGS_Aquagenome/")

#Note that Tajima D and Sweed was run for Canada and Finnmark datasets - although we only looked at selection in Finnmark 
#Can remove Canada analysis here... 

#Read in tajima D for each group 
#Canada is all canadian sites
#Finnmark here is Barents/White sea sites
tajima_canada=read.table("TajimaD/Canada_Results_TajimaD.txt", header=T)
tajima_finmmark=read.table("TajimaD/Finnmark_Results_TajimaD.txt", header=T)

#For Canada dataset - remove NAs
tajima_canada=tajima_canada[!is.nan(tajima_canada$TajimaD),]
order_can_tajima=tajima_canada[order(tajima_canada$TajimaD, decreasing = T),]

#Subset 'outliers' - top 2.5% and lower 2.5% 
outlier_tajima_canada=rbind(order_can_tajima[which(order_can_tajima$TajimaD >= quantile(order_can_tajima$TajimaD, 0.975)),], 
                            order_can_tajima[which(order_can_tajima$TajimaD <= quantile(order_can_tajima$TajimaD, 0.025)),] )

#set values that are not outliers as NA
order_can_tajima_edit=order_can_tajima
order_can_tajima_edit$TajimaD[which(order_can_tajima$TajimaD < quantile(order_can_tajima$TajimaD, 0.975) &
                                      order_can_tajima$TajimaD > quantile(order_can_tajima$TajimaD, 0.025))] <- NA
View(order_can_tajima_edit)

#Repeat same thing for Finnmark dataset
#remove NAs, add SNP name
tajima_finmmark=tajima_finmmark[!is.nan(tajima_finmmark$TajimaD),]
tajima_finmmark$SNP=interaction(tajima_finmmark$CHROM,tajima_finmmark$BIN_START)

#Order data by Tajima D and subset outliers
order_finnmark_tajima=tajima_finmmark[order(tajima_finmmark$TajimaD, decreasing = T),]
outlier_tajima_finnmark=rbind(order_finnmark_tajima[which(order_finnmark_tajima$TajimaD >= quantile(order_finnmark_tajima$TajimaD, 0.975)),], 
                              order_finnmark_tajima[which(order_finnmark_tajima$TajimaD <= quantile(order_finnmark_tajima$TajimaD, 0.025)),] )
#check data
View(outlier_tajima_finnmark[which(outlier_tajima_finnmark$CHROM==5 & outlier_tajima_finnmark$BIN_START > 22000001 ),])
head(outlier_tajima_finnmark)
head(outlier_tajima_canada)

#set non-outliers to NA 
order_finnmark_tajima_edit=order_finnmark_tajima
order_finnmark_tajima_edit$TajimaD[which(order_finnmark_tajima_edit$TajimaD < quantile(order_finnmark_tajima_edit$TajimaD, 0.975) &
                                      order_finnmark_tajima_edit$TajimaD > quantile(order_finnmark_tajima_edit$TajimaD, 0.025))] <- NA
#check data
order_finnmark_tajima_edit

#### Plotting results
library(ggman)
min(outlier_tajima_finnmark$TajimaD)
max(outlier_tajima_finnmark$TajimaD)



plot_tajima_d_finnmark<-ggman(gwas = order_finnmark_tajima,snp="SNP", bp = "BIN_START", title="A",
                             chrom = "Chr2", pvalue = "TajimaD", logTransform = F, pointSize = 1)+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=15), 
                   axis.text = element_text(size=12))+
  geom_hline(yintercept = quantile(order_finnmark_tajima$TajimaD, 0.025), lty=2)+
  geom_hline(yintercept = quantile(order_finnmark_tajima$TajimaD, 0.975), lty=2)+
  ylab("Tajima's D")+xlab("Chromosome") +ylim(min(outlier_tajima_finnmark$TajimaD)-0.5,
                                              max(outlier_tajima_finnmark$TajimaD)+0.5 )

#
####### Next read in SweeD results

#Read in data
sweed_canada=read.table("Sweed/Results/canada/Canada_Results_SweeD.txt", header=T)
sweed_finmmark=read.table("Sweed/Results/finnmark/Finnmark_Results_SweeD.txt", header=T)

#order sweep data for Canada - order by likelihood
order_can_sweed=sweed_canada[order(sweed_canada$Likelihood, decreasing = T),]
nrow(order_can_sweed)*0.05
#Subset outliers for Canada (top 5%)
outlier_sweed_canada=order_can_sweed[which(order_can_sweed$Likelihood >= quantile(order_can_sweed$Likelihood, 0.95)),]

#set non outliers to NAs for Canada - this will be used for sliding window analysis
order_sweed_canada_edit=order_can_sweed
dim(order_sweed_canada_edit)
order_sweed_canada_edit$Likelihood[which(order_sweed_canada_edit$Likelihood < quantile(order_sweed_canada_edit$Likelihood, 0.95))] <-NA

### Repeat for Finnmark dataset 
#Order sweep data for Finnmark  - order by Likelihood
order_finnmark_sweed=sweed_finmmark[order(sweed_finmmark$Likelihood, decreasing = T),]
nrow(order_finnmark_sweed)*0.05
#Subset outliers (top 5%)
outlier_sweed_finnmark=order_finnmark_sweed[which(order_finnmark_sweed$Likelihood >= quantile(order_finnmark_sweed$Likelihood, 0.95)),]

#check data
min(outlier_sweed_finnmark$Likelihood)
#set non outliers to NAs for Finmark
order_finnmark_sweed_edit <- order_finnmark_sweed
order_finnmark_sweed_edit$Likelihood[which(order_finnmark_sweed_edit$Likelihood < quantile(order_finnmark_sweed_edit$Likelihood, 0.95))] <-NA

##Add SNP identifier for finnmark data
head(order_finnmark_sweed)
order_finnmark_sweed$SNP=interaction(order_finnmark_sweed$Chr2, order_finnmark_sweed$Position)
#make sure likelihood is numeric
order_finnmark_sweed$Likelihood=as.numeric(as.character(order_finnmark_sweed$Likelihood))


finnmark_sweed_plot=ggman(gwas = order_finnmark_sweed,snp="SNP", bp = "Position", 
                          chrom = "Chr2",sigLine = NA, title="B",
                          pvalue = "Likelihood", logTransform = F, pointSize = 1)+
  scale_y_continuous(limits = c(0, max(order_finnmark_sweed$Likelihood)+10))+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=15), 
                   axis.text = element_text(size=12))+
  geom_hline(yintercept = quantile(order_finnmark_sweed$Likelihood, 0.95), lty=2)+
  ylab("CLR")+xlab("Chromosome")


#Save as
library(patchwork)
plot_TajD_sweed<- plot_tajima_d_finnmark / finnmark_sweed_plot


ggsave(plot=plot_TajD_sweed,
       filename = "SupplementFig_TajD_CLR.png",
       width = 14.57,
       height = 9.33, units="in",
       device='png', dpi=320,
       bg = NULL)


##Window scanr


###Add "blank" data at ends of chromosomes to ensure all data is evaluated (data at ends of chr)
blank_data = list()

#add blank position at 1MBP for the max. --- for loop gets max value for each chromosome
for(i in 1:29){
  dat_EU=as.data.frame((max(order_finnmark_sweed_edit$Position[which(order_finnmark_sweed_edit$Chr2==i)])))
  dat_EU$i <- i  # maybe you want to keep track of which iteration produced it?
  blank_data[[i]] <- dat_EU # add it to your list  }
}

blank_data_chr = do.call(rbind, blank_data)
colnames(blank_data_chr)=c("Max_Pos", "Chr")
blank_data_chr$blank_Pos=blank_data_chr$Max_Pos+1.1e6

#Create blank data frame to add to include extra windows at the end of chromosome ... blank position have extra 1MBP added to the chromosome
#Make dataframe the same as the one you want to add to
head(order_finnmark_sweed_edit)

blank_data_to_add_FST=as.data.frame(cbind( 
  "Position"= blank_data_chr$blank_Pos,
  "Likelihood"=rep(NA),
  "Alpha"=rep(NA),
  "Chromosome"=blank_data_chr$Chr, 
  "Chr2"=blank_data_chr$Chr
))

#rbind both datasets (blank and real)
head(order_finnmark_sweed_edit)
head(blank_data_to_add_FST)

order_finnmark_sweed_edit_EXTRA <- rbind(as.data.frame(order_finnmark_sweed_edit[,c(1:3,6:7)]), blank_data_to_add_FST)

View(order_finnmark_sweed_edit)
finnmark_Sweed_sliding=winScan(x = order_finnmark_sweed_edit_EXTRA, groups = "Chr2", position = "Position",
                               win_size = 1000000, win_step = 500000,  values = "Likelihood", funs = "mean")


sum(finnmark_Sweed_sliding$Likelihood_n)
nrow(order_finnmark_sweed_edit_EXTRA[!is.na(order_finnmark_sweed_edit_EXTRA$Likelihood),])

#Inlcude blank data in case it gets removed due to missing positions at end of chromosome windows
blank_data = list()
#add blank position at 1MBP for the max. --- for loop gets max value for each chromosome
for(i in 1:29){
  dat_EU=as.data.frame((max(order_finnmark_tajima_edit$BIN_START[which(order_finnmark_tajima_edit$Chr2==i)])))
  dat_EU$i <- i  # maybe you want to keep track of which iteration produced it?
  blank_data[[i]] <- dat_EU # add it to your list  }
}
blank_data_chr = do.call(rbind, blank_data)
colnames(blank_data_chr)=c("Max_Pos", "Chr")
blank_data_chr$blank_Pos=blank_data_chr$Max_Pos+1.1e6

#Create blank data frame to add to include extra windows at the end of chromosome ... blank position have extra 1MBP added to the chromosome
#Make dataframe the same as the one you want to add to
head(order_finnmark_tajima_edit)
blank_data_to_add_FST=as.data.frame(cbind("CHROM"=blank_data_chr$Chr, 
                                         "BIN_START"= blank_data_chr$blank_Pos,
                                         "N_SNPS"=rep(NA),
                                         "TajimaD"=rep(NA),
                                         "Chromosome"=blank_data_chr$Chr, 
                                         "Chr2"=blank_data_chr$Chr,
                                         "SNP"=rep("blank")
                                         
))
#rbind both datasets (blank and real)
order_finnmark_tajima_edit_EXTRA=rbind(as.data.frame(order_finnmark_tajima_edit), blank_data_to_add_FST)
order_finnmark_tajima_edit_EXTRA$BIN_START=as.numeric(as.character(order_finnmark_tajima_edit_EXTRA$BIN_START))

finnmark_Tajima_sliding=winScan(x = order_finnmark_tajima_edit_EXTRA, groups = "Chr2", position = "BIN_START",
                                win_size = 1000000, win_step = 500000,  values = "TajimaD", funs = "mean")


outlier_fin_taj_windows=finnmark_Tajima_sliding[which(finnmark_Tajima_sliding$TajimaD_n>0),]
outlier_fin_sweed_windows=finnmark_Sweed_sliding[which(finnmark_Sweed_sliding$Likelihood_n>0),]

outlier_fin_taj_windows$SNP=interaction(outlier_fin_taj_windows$Chr2, outlier_fin_taj_windows$win_start, sep="_")
outlier_fin_sweed_windows$SNP=interaction(outlier_fin_sweed_windows$Chr2, outlier_fin_sweed_windows$win_start, sep="_")

overlap_selection_finnmark=merge(outlier_fin_taj_windows,outlier_fin_sweed_windows, by="SNP")
nrow(overlap_selection_finnmark)

#############
outliers_introgression=read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Compare_Introgression_Outliers/Overlap_2metrics_introgression_jan2025.txt", header=T)
merge(overlap_selection_finnmark,outliers_introgression, by="SNP" )

##Get regions for Taj D that overlap with introgression tests
finnmark_tajima_outlierIntrogression_Adaptive=merge(outlier_fin_taj_windows,outliers_introgression, by="SNP")

##Get regions for sweed that overlap with introgression tests
finnmark_sweed_outlierIntrogression_Adaptive=merge(outlier_fin_sweed_windows,outliers_introgression, by="SNP")
View(finnmark_tajima_outlierIntrogression_Adaptive)

nrow(finnmark_tajima_outlierIntrogression_Adaptive)
nrow(finnmark_sweed_outlierIntrogression_Adaptive)


###Save results:

write.table(finnmark_sweed_outlierIntrogression_Adaptive, "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Compare_Introgression_Outliers/Finnmark_SelectionSweed_OutliersIntrogression_jan2025.txt", quote = F, row.names = F, col.names = T, sep='\t')
write.table(finnmark_tajima_outlierIntrogression_Adaptive, "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Compare_Introgression_Outliers/Finnmark_SelectionTajima_OutliersIntrogression_jan2025.txt", quote = F, row.names = F, col.names = T, sep='\t')


#
overlap_both_selection_and_intro<- merge(finnmark_sweed_outlierIntrogression_Adaptive, finnmark_tajima_outlierIntrogression_Adaptive, by="SNP")
write.table(overlap_both_selection_and_intro, "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Compare_Introgression_Outliers/Finnmark_SelectionTajima_SweeD_OutliersIntrogression_OVERLAP2_bothMethods_Jan2025.txt", quote = F, row.names = F, col.names = T, sep='\t')
View(overlap_both_selection_and_intro)
