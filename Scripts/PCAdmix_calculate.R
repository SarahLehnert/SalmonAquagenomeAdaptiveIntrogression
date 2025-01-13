library(qqman)
library(raster)
library(ggmap)
library(rgdal)
library(RColorBrewer)
library(broom)
library(maps) # tool for maps
library(mapdata) # all your basemaps are here
library(marmap) # for bathymetry if needed
library(gplots) # for colour range
library(rworldmap)
library(mapplots)
library(rgeos)
library(maptools)
library(lattice)

#Set directory
setwd("~/Desktop/Software/PCAdmix/")

#Read in results from pcadapt (.vit and position of windows .markers)


#For loop here is strickly for plotting results (exploring data) per chromsome - other script below will combine data for all chromosomes and save
#plots were not used in manuscript 

#Do for ssa01-ssa09 first
for(i in 1:9){
  #Read in results (ancestry for windows represented by 1 and 0s for canada and east atlantic)
  data1_conf=read.table(paste0("ssa0",i,"_WGS_FinnmarkIntrogress_jan_2025.vit.txt"), header=F)
  
  #remove first column representing individual names
  data_chr1_toplot=data1_conf[,-1]
  #convert to matrix
  data_chr1_toplot=as.matrix(data_chr1_toplot)
  
  #Read in window positions
  windows=read.table(paste0("ssa0",i,"_WGS_FinnmarkIntrogress_jan_2025.markers.txt"), header=F, fill=T)
  start_poistion=windows$V2
  windows$Start_Position=gsub(start_poistion, pattern = paste0("ssa0",i,":"), replacement = "")
  
  #Get mean value per column (SNP windows) and add position information
  data_chr1=as.data.frame(cbind(colMeans(data_chr1_toplot, na.rm = T), windows$Start_Position))
  
  #ensure values are numeric
  data_chr1$V1=as.numeric(as.character(data_chr1$V1))
  data_chr1$V2=as.numeric(as.character(data_chr1$V2))
  data_chr1[which.max(data_chr1$V1),]
  
  #Plot of values across chromosome - not used in paper - use done for visual representation
  temp_plot=ggplot()+geom_point(data=data_chr1, aes(x=V2, y=V1),size=1)+ylab("Mean Canadian ancestry (PCAdmix)")
  #Save plot per chromosome
  ggsave(temp_plot, file=paste0("Plot_PCAdmix_ssa0", i,".png"), width = 16, height = 10, units = "cm")
  
  #Plot heatmap of results for each chromosme - showing different ancestries in blue and yellow
  mycol<-colorRampPalette(c("navyblue", "yellow"))
  
  png(filename = paste0("Heatmap_ssa0", i,".png"),width = 800, height = 850 )
  heatmap=heatmap.2(data_chr1_toplot,
                    Rowv=FALSE, #rows should be reordered as required
                    Colv = "Rowv", #columns should be treated as rows
                    dendrogram="none", #no trees
                    scale="none",
                    breaks=100, #number of break points used to bin into colours
                    #col=gray.colors(n=99),
                    col=mycol,
                    trace="none", #whether lines should be drawn between cols, rows,
                    margins=c(3,3),#margins for column names and row names
                    labRow= "",
                    labCol= "",
                    #cexCol=0.4, #column label size
                    #cexRow=0.4,#row label size
                    #srtCol = 90,#col label angle, degrees from horizontal
                    #srtRow = 90,
                    key=TRUE,
                    keysize = 1,
                    density.info="none"
                    #lmat = rbind(c(3,1),c(4,2)), #order of elements (1 is heatmap, 4 is key) #makes the plot a 2x2 grid, each "cell" has 1 element - row tree, key, col tree, and heatmap
                    # lhei = c(20,5), #row height for plot elements
                    # lwid = c(8,30)  #column width for plot elements
  )
  dev.off()
  
}

#Same as aboave
#Do for ssa10-ssa29 
for(i in 10:29){
  
  #Read in results (ancestry for windows represented by 1 and 0s for canada and east atlantic)
  data1_conf=read.table(paste0("ssa",i,"_WGS_FinnmarkIntrogress_jan_2025.vit.txt"), header=F)
  
  #remove first column representing individual names
  data_chr1_toplot=data1_conf[,-1]
  #convert to matrix
  data_chr1_toplot=as.matrix(data_chr1_toplot)
  
  #Read in window positions
  windows=read.table(paste0("ssa",i,"_WGS_FinnmarkIntrogress_jan_2025.markers.txt"), header=F, fill=T)
  start_poistion=windows$V2
  windows$Start_Position=gsub(start_poistion, pattern = paste0("ssa",i,":"), replacement = "")
  
  #Get mean value per column (SNP windows) and add position information
  data_chr1=as.data.frame(cbind(colMeans(data_chr1_toplot, na.rm = T), windows$Start_Position))
  
  #ensure values are numeric
  data_chr1$V1=as.numeric(as.character(data_chr1$V1))
  data_chr1$V2=as.numeric(as.character(data_chr1$V2))
  data_chr1[which.max(data_chr1$V1),]
  
  #Plot of values across chromosome - not used in paper - use done for visual representation
  temp_plot=ggplot()+geom_point(data=data_chr1, aes(x=V2, y=V1),size=1)+ylab("Mean Canadian ancestry (PCAdmix)")
  #Save plot per chromosome
  ggsave(temp_plot, file=paste0("Plot_PCAdmix_ssa", i,".png"), width = 16, height = 10, units = "cm")
  
  #Plot heatmap of results for each chromosme - showing different ancestries in blue and yellow
  mycol<-colorRampPalette(c("navyblue", "yellow"))
  
  png(filename = paste0("Heatmap_ssa", i,".png"),width = 800, height = 850 )
  heatmap=heatmap.2(data_chr1_toplot,
                    Rowv=FALSE, #rows should be reordered as required
                    Colv = "Rowv", #columns should be treated as rows
                    dendrogram="none", #no trees
                    scale="none",
                    breaks=100, #number of break points used to bin into colours
                    #col=gray.colors(n=99),
                    col=mycol,
                    trace="none", #whether lines should be drawn between cols, rows,
                    margins=c(3,3),#margins for column names and row names
                    labRow= "",
                    labCol= "",
                    #cexCol=0.4, #column label size
                    #cexRow=0.4,#row label size
                    #srtCol = 90,#col label angle, degrees from horizontal
                    #srtRow = 90,
                    key=TRUE,
                    keysize = 1,
                    density.info="none"
                    #lmat = rbind(c(3,1),c(4,2)), #order of elements (1 is heatmap, 4 is key) #makes the plot a 2x2 grid, each "cell" has 1 element - row tree, key, col tree, and heatmap
                    # lhei = c(20,5), #row height for plot elements
                    # lwid = c(8,30)  #column width for plot elements
  )
                    dev.off()
  rm(data_chr1_toplot) 

}

#######################################################################################
#This script is similar to above - but will save all the results for each chromosome together to use for analyses

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
#write.table(all_PCAdmix_results, "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/All_PCAdmix_results_ProportionCanadian.txt", quote = F, row.names = F, col.names = T, sep="\t")



#Plot manhattan plot for comparisons
all_PCAdmix_results$SNP=interaction(all_PCAdmix_results$Chr, all_PCAdmix_results$Start_Position, sep = "_")

##Selecting outlier regions

#Top 5 % (95%) were chosen as outlier regions
quantile(all_PCAdmix_results$CanadaAncestry, 0.95)
#Save outliers in seperate file
write.table(all_PCAdmix_results[which(all_PCAdmix_results$CanadaAncestry >= quantile(all_PCAdmix_results$CanadaAncestry, 0.95) ),], "~/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/Outliers_PCAdmix_results_ProportionCanadian.txt", quote = F, row.names = F, col.names = T, sep="\t")

#Max position
all_PCAdmix_results[which.max(all_PCAdmix_results$CanadaAncestry),]


manhattan(x=all_PCAdmix_results, 
          snp = "SNP", chr = "Chr2",bp = "Start_Position", 
          p="CanadaAncestry", ylim=c(0, 0.8), logp = F, 
          ylab="Proportion Canada Ancestry" ,
          highlight =all_PCAdmix_results$SNP[which(all_PCAdmix_results$CanadaAncestry >= quantile(all_PCAdmix_results$CanadaAncestry, 0.95) )] )

########################################################################

#Map of site locations and pie chart of proportion ancestry assigned to canada vs south norway
#read coordinates for all sites
coords=read.csv("~/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Site_codes.csv", header=T)

data_table_to_save_ancest=read.table(paste0("ssa0",1,"_WGS_FinnmarkIntrogress_jan_2025.ia.txt"), header=T)[,1]

#read in proportion ancestry for each chromosome
for(i in 1:9){
  
  data1_conf=read.table(paste0("ssa0",i,"_WGS_FinnmarkIntrogress_jan_2025.ia.txt"), header=T)
  
  ances0=as.data.frame(as.numeric(as.character(data1_conf[,2])))
  colnames(ances0)=paste0("ssa0", i)
  data_table_to_save_ancest=cbind(data_table_to_save_ancest,ances0)
  
}

for(i in 10:29){
  
  data1_conf=read.table(paste0("ssa",i,"_WGS_FinnmarkIntrogress_jan_2025.ia.txt"), header=T)
  
  ances0=as.data.frame(as.numeric(as.character(data1_conf[,2])))
  colnames(ances0)=paste0("ssa", i)
  data_table_to_save_ancest=cbind(data_table_to_save_ancest,ances0)
  
}
#Combine results

all_ancestry_results=data_table_to_save_ancest

#Get mean ancestry across all chromosomes per individual
all_ancestry_results$MEAN=rowMeans(all_ancestry_results[,2:30]) #for columns with chromosomes

#Edit IDs to get pop names
all_ancestry_results$pop=gsub("[^_]*_(.*)", "\\1", all_ancestry_results$data_table_to_save_ancest)
all_ancestry_results$pop=gsub("\\_.*", "", all_ancestry_results$pop)

#Get mean ancestry by pop
mean_hap_ancestry=as.data.frame(aggregate(all_ancestry_results$MEAN~ all_ancestry_results$pop, FUN=mean))

colnames(mean_hap_ancestry)=c("Pop", "MeanHapAnces")

#Coordinates for plotting
coords=read.csv("~/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Site_codes.csv", header=T)

sites_ances=merge(x=coords, by.x=2, y=mean_hap_ancestry, by.y=1)

#Plot map of proportion ancestry
Sample.Lat.lim=c(68, 72)
Sample.Long.lim=c(10,35)

par(mar=c(2,2,2,2)) 

sites_ances$MeanHapAnces_1=sites_ances$MeanHapAnces*100
sites_ances$MeanHapAnces_2=100-(sites_ances$MeanHapAnces*100)

maps::map("worldHires", xlim=Sample.Long.lim, ylim=Sample.Lat.lim, mar=rep(1,4),
    col="lightgray", fill=TRUE, lwd=0.001, border=F, resolution=0, lwd=0.2);map.axes()
points(x=sites_ances$Long,  
       y=sites_ances$Lat, pch=25, col="black",bg="yellow", cex=0.75)
text(x=sites_ances$Long,  
     y=sites_ances$Lat, label= sites_ances$Site)
for(i in 1:nrow(sites_ances))
{
  add.pie(as.integer(sites_ances[i,c("MeanHapAnces_1", "MeanHapAnces_2")]),   #Data for making pie proportions
          x=sites_ances$Long[i],y=sites_ances$Lat[i],labels="",radius = 0.1,   #Coordinates for pie and radius is size of pie
          col=c("firebrick2","dodgerblue"))   #Colours for pie
}

#run map script again if looks weird
