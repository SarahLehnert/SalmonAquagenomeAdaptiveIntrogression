
## R script to combine results for introgression into multipanel figure (treemix, dxy, pcadmix)

setwd("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/")
library(ggman)
library(stringr)
library(ggplot2)

#Get dxy outliers/windows
salmondxy=data.table::fread("AllDyx_windows100kbp_Chr1_29_StandardizeDxy_SNor_Finnmark.txt", header=T)
outliers_dyx=read.table("DxyOutliers_windows_Chr1_29_StandardizeDxy_SNor_Finnmark.txt", header=T)
max(outliers_dyx$diff_dxy_mean)

#check data
head(salmondxy)
salmondxy[is.na(salmondxy$diff_dxy_mean),]
max(salmondxy$diff_dxy_mean, na.rm=T)

#Plot dxy results
dxy_plot <- ggman(gwas = salmondxy[!is.na(salmondxy$diff_dxy_mean),],
               snp="SNP",
               bp = "win_start",
               chrom = "Chr2",
               pvalue = "diff_dxy_mean", logTransform = F, pointSize = 1)+
  theme_bw()+theme(panel.grid = element_blank(), title=element_blank(),
        axis.title = element_text(size=15),
                   axis.text = element_text(size=12))+
  ylim(c(min(salmondxy$diff_dxy_mean, na.rm=T)-0.1, max(salmondxy$diff_dxy_mean, na.rm=T)+0.1))+
  ylab("Difference dxy")+xlab("Chromosome")+geom_hline(yintercept = max(outliers_dyx$diff_dxy_mean), lty=2)

#Read in PCAdmix results
salmon_PCAdxmi=data.table::fread("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/Results/All_PCAdmix_results_ProportionCanadian_jan2025.txt", header=T)
head(salmon_PCAdxmi)
salmon_PCAdxmi$SNP=interaction(salmon_PCAdxmi$Chr2, salmon_PCAdxmi$Start_Position, sep = "_")
#Outliers only
outliers_pca=read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/Results/Outliers_PCAdmix_results_ProportionCanadian_jan2025.txt", header=T)
min(outliers_pca$CanadaAncestry)

salmon_PCAdxmi$Mid_position=(salmon_PCAdxmi$Start_Position+salmon_PCAdxmi$End_Position)/2
max(salmon_PCAdxmi$CanadaAncestry)

pcadmix_plot <- ggman(gwas = salmon_PCAdxmi,
                   snp="SNP", bp = "Mid_position", 
                   chrom = "Chr2",sigLine = NA,
                   pvalue = "CanadaAncestry", logTransform = F, pointSize = 1)+
  theme_bw()+theme(axis.title = element_text(size=15), title=element_blank(),
                    axis.text = element_text(size=12),
                   panel.grid = element_blank())+scale_y_continuous(limits=c(0, max(salmon_PCAdxmi$CanadaAncestry)+0.2))+
  ylab("Proportion Canadian ancestry")+xlab("Chromosome")+geom_hline(yintercept = min(outliers_pca$CanadaAncestry), lty=2 )

###Add treemix results to make a single plot:


#read treemix data
treemix<-read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/Results/WGS_salmonpops_results_threepop_final_by_Region_500snp_window_jan2025", header=F)
head(treemix)

#Quick plot of treemix f3 stats
treemix2<-treemix[order(treemix$V2, decreasing = F),]
plot(treemix2[which(treemix2$V2<0.015),]$V2)

#Sepearte columns into groups
column1_2<-as.data.frame(str_split_fixed(treemix2$V1, ";", 2))
extra_columns<-as.data.frame(str_split_fixed(column1_2$V2, ",", 2))

#data frame for plotting
treemix_edited <- as.data.frame(cbind(column1_2$V1, extra_columns, treemix2$V2, treemix2$V3, treemix2$V4))
colnames(treemix_edited)<- c("Admixed", "Source1", "Source2", "F3stat", "SE", "Zscore")

treemix_edited$Popcomparison<- c(1:nrow(treemix_edited))
treemix_edited$significance<- c("nonsign")
treemix_edited$significance[which(treemix_edited$Zscore < -2.5 )]<-"sign"

#Plot 
treemix<- ggplot()+geom_hline(yintercept = 0,  linetype = "longdash")+
  geom_point(data=treemix_edited, aes(y=F3stat, 
  x=Popcomparison, fill=Admixed), size=4, shape=21)+ylim(-0.01, 0.22)+
  scale_fill_manual(values=c("#E21F26","#F69898", "#37A047", "#B5D88A", "#6C3E98", "#FDBD6C"))+
  theme_bw() + 
  xlab("Population comparison") + ylab("F3 statistic")+
  theme(axis.title = element_text(size=15), 
        axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15),
        legend.position = c(0.2, 0.80), 
        legend.box.background = element_rect(colour = "black"))



#Save plot for manuscript Figure. 
figure3 <- treemix + (dxy_plot / pcadmix_plot) + plot_layout(widths = c(1.3,2))


ggsave(plot=figure3,
  filename = "Figure3_treemix_dxy_pcadmix_results.pdf",
  width = 16.70,
  height = 7.32,
  units = "in",
  dpi = 300,
  bg = NULL)

