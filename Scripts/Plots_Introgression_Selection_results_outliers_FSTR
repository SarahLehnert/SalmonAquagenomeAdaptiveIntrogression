library(ggplot2)
setwd("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Figure_introgress_selection/")


# Scripts wiht plots for introgression and selection results with plots showing differences in FST
#One figure from main text, and other figures were in supplement

#Read in FST results for both comparisons
fst <- data.table::fread("~/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_jan2025_updated_ids_for_FST_Finnmark_Canada.fst", header=T)
fst2 <- data.table::fread("~/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST_SouthNor_Canada.fst", header=T)

#merge results by SNP
compare <- merge(fst, fst2, "SNP")
head(compare)

#Get difference in fst at each locus
compare$Diff_Fst=compare$FST.x-compare$FST.y

#Test plot for Chr 21
ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa21"),],
                         aes(y=Diff_Fst,x=POS.x))+xlim(34400001,35500000)+geom_hline(yintercept = 0, col="blue", lty=2)+ylab("Difference FST")

#read in rsults for Sweed and Tajima D
sweed_finn<-read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/SweeD/Results/finnmark/Finnmark_Results_SweeD.txt", header=T)
tajimad_finn<-read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/TajimaD/Finnmark_Results_TajimaD.txt", header=T)

head(sweed_finn)
head(tajimad_finn)

#Read in pcadmix results
pcadmix_order<-read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/Results/All_PCAdmix_results_ProportionCanadian_jan2025.txt", header=T)
dxy_order<-read.table("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/AllDyx_windows100kbp_Chr1_29_StandardizeDxy_SNor_Finnmark.txt", header=T)

min(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==6)], na.rm=T)
max(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==6)], na.rm=T)
max(pcadmix_order$CanadaAncestry[which(pcadmix_order$Chr2==6)], na.rm=T)


introgress_6<- ggplot()+
  geom_rect(aes(xmin=32500001, xmax=33500000, ymin=-1, ymax=1), alpha=0.2)+
  geom_rect(aes(xmin=49000001, xmax=50000000, ymin=-1, ymax=1), alpha=0.2)+
  geom_rect(aes(xmin=63000001, xmax=64000000, ymin=-1, ymax=1), alpha=0.2)+
   geom_smooth(data=pcadmix_order[which(pcadmix_order$Chr2==6),], stat="identity",
              aes(x=Start_Position, y=CanadaAncestry, col="PCAdmix (Proportion Canadian Ancestry)"), se=FALSE )+
  geom_smooth(data=dxy_order[which(dxy_order$Chr2==6),], stat="identity",
              aes(x=win_start, y=diff_dxy_mean, col="Difference Dxy"), se=FALSE )+
  ylab("Metric for introgression")+
  scale_colour_manual(values=c("brown1", "chartreuse4", "cornflowerblue"))+
  theme_classic()+theme(legend.position="top", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.5,0.5), expand = F)+
  xlab("")
NULL


sweed_6 <- ggplot()+  geom_rect(aes(xmin=32500001, xmax=33500000, ymin=0, ymax=70), alpha=0.2)+
  geom_rect(aes(xmin=49000001, xmax=50000000, ymin=0, ymax=70), alpha=0.2)+
  geom_rect(aes(xmin=63000001, xmax=64000000, ymin=0, ymax=70), alpha=0.2)+
 geom_smooth(data=sweed_finn[which(sweed_finn$Chromosome ==6),],stat="identity",
            aes(y=Likelihood ,x= Position), col="goldenrod")+

  theme_classic()+theme(axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(0,(max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==6)])+10)),  
                  expand = F)+xlab("Ssa06 position")+ylab("CLR")

  

min6<- min(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==6)],na.rm=T)
max6<- max(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==6)],na.rm=T)


taji_6<-ggplot()+  geom_rect(aes(xmin=32500001, xmax=33500000, ymin=min6-0.2, ymax=max6+0.2), alpha=0.2)+
  geom_rect(aes(xmin=49000001, xmax=50000000, ymin=min6-0.2, ymax=max6+0.2), alpha=0.2)+
  geom_rect(aes(xmin=63000001, xmax=64000000, ymin=min6-0.2, ymax=max6+0.2), alpha=0.2)+
 geom_smooth(data=tajimad_finn[which(tajimad_finn$Chromosome ==6),], stat="identity",
             aes(y=TajimaD ,x= BIN_START), col="blue")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(min6-0.2, max6+0.2), expand = F)+xlab("")+ylab("Tajima's D")


##SIZE OF FIGURE
## PDF 12.24 x 11.45

  
library(patchwork)  
plotssa6<- introgress_6 /
taji_6 /
sweed_6 +
plot_annotation(tag_levels = "A") &
theme(plot.tag = element_text(size = 25))




ggsave(plot=plotssa6,
       filename = "Ssa06_introgress_sel_plots.pdf",
       width = 12.24,
       height = 11.45,
       units = "in",
       dpi = 300,
       bg = NULL)

ggsave(plot=plotssa6, 
       filename = "Ssa06_introgress_sel_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)



###################### SSA07


min(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==7)], na.rm=T)
max(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==7)], na.rm=T)
max(pcadmix_order$CanadaAncestry[which(pcadmix_order$Chr2==7)], na.rm=T)


introgress_7<- ggplot()+
  geom_rect(aes(xmin=45000001, xmax=46000000, ymin=-1, ymax=1), alpha=0.2)+
  geom_smooth(data=pcadmix_order[which(pcadmix_order$Chr2==7),], stat="identity",
              aes(x=Start_Position, y=CanadaAncestry, col="PCAdmix (Proportion Canadian Ancestry)"), se=FALSE)+
  
  geom_smooth(data=dxy_order[which(dxy_order$Chr2==7),], stat="identity",
              aes(x=win_start, y=diff_dxy_mean, col="Difference Dxy"), se=FALSE)+
  ylab("Metric for introgression")+
  scale_colour_manual(values=c("brown1", "chartreuse4", "cornflowerblue"))+
  theme_classic()+theme(legend.position="top", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.5,0.5), expand = F)+xlab(" ")
NULL


max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==7)])

sweed_7<- ggplot()+ 
  geom_rect(aes(xmin=45000001, xmax=46000000, ymin=0, ymax=75), alpha=0.2)+
  geom_smooth(data=sweed_finn[which(sweed_finn$Chromosome ==7),],stat="identity",
              aes(y=Likelihood ,x= Position), col="goldenrod")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15), 
                        axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(0, max((sweed_finn$Likelihood[which(sweed_finn$Chromosome ==7)]))+10 ),
                           expand = F)+xlab("Ssa07 position")+ylab("CLR")



min7<- min(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==7)],na.rm=T)
max7<- max(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==7)],na.rm=T)


taji_7<-ggplot()+  
  geom_rect(aes(xmin=45000001, xmax=46000000, ymin=min7-0.2, ymax=max7+0.2), alpha=0.2)+
  geom_smooth(data=tajimad_finn[which(tajimad_finn$Chromosome ==7),], stat="identity",
              aes(y=TajimaD ,x= BIN_START), col="blue")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(min7-0.2, max7+0.2), expand = F)+xlab("")+ylab("Tajima's D")

##SIZE OF FIGURE
## PDF 12.24 x 11.45

plotssa7<- introgress_7 /
  taji_7 /
  sweed_7 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))



ggsave(plot=plotssa7,
       filename = "Ssa07_introgress_sel_plots.pdf",
       width = 12.24,
       height = 11.45,
       units = "in",
       dpi = 300,
       bg = NULL)

ggsave(plot=plotssa7, 
       filename = "Ssa07_introgress_sel_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)



######################################################################

min(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==21)], na.rm=T)
max(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==21)], na.rm=T)
max(pcadmix_order$CanadaAncestry[which(pcadmix_order$Chr2==21)], na.rm=T)


introgress_21<- ggplot()+
  geom_rect(aes(xmin=34000001, xmax=35500000, ymin=-1, ymax=1), alpha=0.2)+
  
  geom_smooth(data=pcadmix_order[which(pcadmix_order$Chr2==21),], stat="identity",
              aes(x=Start_Position, y=CanadaAncestry, col="PCAdmix (Proportion Canadian Ancestry)"), se=FALSE)+
  
  geom_smooth(data=dxy_order[which(dxy_order$Chr2==21),], stat="identity",
              aes(x=win_start, y=diff_dxy_mean, col="Difference Dxy"), se=FALSE)+
  ylab("Metric for introgression")+
  scale_colour_manual(values=c("brown1", "chartreuse4", "cornflowerblue"))+
  theme_classic()+theme(legend.position="top", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.65,0.80), expand = F)+xlab("")
NULL


max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==21)])+10

sweed_21<- ggplot()+ 
  geom_rect(aes(xmin=34000001, xmax=35500000, ymin=0, ymax=135), alpha=0.2)+
  geom_smooth(data=sweed_finn[which(sweed_finn$Chromosome ==21),],stat="identity",
              aes(y=Likelihood ,x= Position), col="goldenrod")+
  theme_classic()+theme(legend.position = "none", 
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(0,max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==21)])+10), 
                  expand = F)+xlab("Ssa21 position")+ylab("CLR")



min21<- min(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==21)],na.rm=T)
max21<- max(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==21)],na.rm=T)


taji_21<-ggplot()+  
  geom_rect(aes(xmin=34000001, xmax=35500000, ymin=min21-0.2, ymax=max21+0.2), alpha=0.2)+
  geom_smooth(data=tajimad_finn[which(tajimad_finn$Chromosome ==21),], stat="identity",
              aes(y=TajimaD ,x= BIN_START), col="blue")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(min21-0.2, max21+0.2), expand = F)+xlab("")+ylab("Tajima's D")

##SIZE OF FIGURE
## PDF 12.24 x 11.45


plotssa21<- introgress_21 /
  taji_21 /
  sweed_21+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))



ggsave(plot=plotssa21,
       filename = "Ssa21_introgress_sel_plots.pdf",
       width = 12.24,
       height = 11.45,
       units = "in",
       dpi = 300,
       bg = NULL)

ggsave(plot=plotssa21, 
       filename = "Ssa21_introgress_sel_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)





######################################################################

min(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==27)], na.rm=T)
max(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==27)], na.rm=T)
max(pcadmix_order$CanadaAncestry[which(pcadmix_order$Chr2==27)], na.rm=T)


introgress_27<- ggplot()+
  geom_rect(aes(xmin=41000001, xmax=42000000, ymin=-1, ymax=1), alpha=0.2)+
  geom_rect(aes(xmin=20500001, xmax=21500000, ymin=-1, ymax=1), alpha=0.2)+
  
  geom_smooth(data=pcadmix_order[which(pcadmix_order$Chr2==27),], stat="identity",
              aes(x=Start_Position, y=CanadaAncestry, col="PCAdmix (Proportion Canadian Ancestry)"), se=F)+
  
  geom_smooth(data=dxy_order[which(dxy_order$Chr2==27),], stat="identity",
              aes(x=win_start, y=diff_dxy_mean, col="Difference Dxy"), se=F)+
  ylab("Metric for introgression")+
  scale_colour_manual(values=c("brown1", "chartreuse4", "cornflowerblue"))+
  theme_classic()+theme(legend.position="top", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.8,0.6), expand = F)+xlab("")
NULL


max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==27)])+10

sweed_27<- ggplot()+ 
  geom_rect(aes(xmin=41000001, xmax=42000000, ymin=0, ymax=70), alpha=0.2)+
  geom_rect(aes(xmin=20500001, xmax=21500000, ymin=0, ymax=70), alpha=0.2)+
    geom_smooth(data=sweed_finn[which(sweed_finn$Chromosome ==27),],stat="identity",
              aes(y=Likelihood ,x= Position), col="goldenrod")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(0,max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==27)])+10),
                  expand = F)+xlab("Ssa27 position")+ylab("CLR")



min27<- min(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==27)],na.rm=T)
max27<- max(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==27)],na.rm=T)


taji_27<-ggplot()+  
  geom_rect(aes(xmin=41000001, xmax=42000000, ymin=min27-0.2, ymax=max27+0.2), alpha=0.2)+
  geom_rect(aes(xmin=20500001, xmax=21500000, min=min27-0.2, ymax=max27+0.2), alpha=0.2)+
  
  geom_smooth(data=tajimad_finn[which(tajimad_finn$Chromosome ==27),], stat="identity",
              aes(y=TajimaD ,x= BIN_START), col="blue")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(min27-0.2, max27+0.2), expand = F)+xlab("")+ylab("Tajima's D")


max(compare$Diff_Fst[which(compare$CHR.x=="ssa27")], na.rm=T)
min(compare$Diff_Fst[which(compare$CHR.x=="ssa27")], na.rm=T)

#checking fst
fst_ssa27<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa27"),],
                    aes(y=Diff_Fst,x=POS.x))+
  geom_point(data=compare[which(compare$CHR.x=="ssa27" & compare$POS.x >= 41000001 & compare$POS.x <= 42000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_hline(yintercept = 0, col="blue", lty=2)+ylab("Difference FST")+xlab("Ssa27 position")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.7,0.6), expand = F)+
  geom_smooth(data = compare[which(compare$CHR.x=="ssa27"),],
              aes(y=Diff_Fst,x=POS.x), #stat="identity")+
              method="loess",  span = 0.001, se=FALSE) + 
  
  NULL  


##SIZE OF FIGURE
## PDF 12.24 x 11.45



plotssa27<- introgress_27 /
  taji_27 /
  sweed_27 +
plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))



ggsave(plot=plotssa27,
       filename = "Ssa27_introgress_sel_plots.pdf",
       width = 12.24,
       height = 11.45,
       units = "in",
       dpi = 300,
       bg = NULL)

ggsave(plot=plotssa27, 
       filename = "Ssa27_introgress_sel_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)


######################

##Ssa18
min(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==18)], na.rm=T)
max(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==18)], na.rm=T)
max(pcadmix_order$CanadaAncestry[which(pcadmix_order$Chr2==18)], na.rm=T)


introgress_18<- ggplot()+
  geom_rect(aes(xmin=24500001, xmax=25500000, ymin=-1, ymax=1), alpha=0.2)+

  geom_smooth(data=pcadmix_order[which(pcadmix_order$Chr2==18),], stat="identity",
              aes(x=Start_Position, y=CanadaAncestry, col="PCAdmix (Proportion Canadian Ancestry)"), se=F)+
  
  geom_smooth(data=dxy_order[which(dxy_order$Chr2==18),], stat="identity",
              aes(x=win_start, y=diff_dxy_mean, col="Difference Dxy"), se=F)+
  ylab("Metric for introgression")+
  scale_colour_manual(values=c("brown1", "chartreuse4", "cornflowerblue"))+
  theme_classic()+theme(legend.position="top", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-1,0.6), expand = F)+xlab("")
NULL


max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==18)])+10

sweed_18<- ggplot()+ 
  geom_rect(aes(xmin=24500001, xmax=25500000, ymin=0, ymax=120), alpha=0.2)+
  geom_smooth(data=sweed_finn[which(sweed_finn$Chromosome ==18),],stat="identity",
              aes(y=Likelihood ,x= Position), col="goldenrod")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(0,max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==18)])+10),
                  expand = F)+xlab("Ssa18 position")+ylab("CLR")



min18<- min(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==18)],na.rm=T)
max18<- max(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==18)],na.rm=T)


taji_18<-ggplot()+  
  geom_rect(aes(xmin=24500001, xmax=25500000, ymin=min18-0.2, ymax=max18+0.2), alpha=0.2)+

  geom_smooth(data=tajimad_finn[which(tajimad_finn$Chromosome ==18),], stat="identity",
              aes(y=TajimaD ,x= BIN_START), col="blue")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(min18-0.2, max18+0.2), expand = F)+xlab("")+ylab("Tajima's D")




##SIZE OF FIGURE
## PDF 12.24 x 11.45



plotssa18<- introgress_18 /
  taji_18 /
  sweed_18 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))



ggsave(plot=plotssa18,
       filename = "Ssa18_introgress_sel_plots.pdf",
       width = 12.24,
       height = 11.45,
       units = "in",
       dpi = 300,
       bg = NULL)

ggsave(plot=plotssa18, 
       filename = "Ssa18_introgress_sel_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)




##############

######################

##Ssa17
min(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==17)], na.rm=T)
max(dxy_order$diff_dxy_mean[which(dxy_order$Chr2==17)], na.rm=T)
max(pcadmix_order$CanadaAncestry[which(pcadmix_order$Chr2==17)], na.rm=T)


introgress_17<- ggplot()+
  geom_rect(aes(xmin=48500001, xmax=49500000, ymin=-1, ymax=1), alpha=0.2)+
  
  geom_smooth(data=pcadmix_order[which(pcadmix_order$Chr2==17),], stat="identity",
              aes(x=Start_Position, y=CanadaAncestry, col="PCAdmix (Proportion Canadian Ancestry)"), se=F)+
  
  geom_smooth(data=dxy_order[which(dxy_order$Chr2==17),], stat="identity",
              aes(x=win_start, y=diff_dxy_mean, col="Difference Dxy"), se=F)+
  ylab("Metric for introgression")+
  scale_colour_manual(values=c("brown1", "chartreuse4", "cornflowerblue"))+
  theme_classic()+theme(legend.position="top", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.5,0.4), expand = F)+xlab("")
NULL


max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==17)])+10

sweed_17<- ggplot()+ 
  geom_rect(aes(xmin=48500001, xmax=49500000, ymin=0, ymax=120), alpha=0.2)+
  geom_smooth(data=sweed_finn[which(sweed_finn$Chromosome ==17),],stat="identity",
              aes(y=Likelihood ,x= Position), col="goldenrod")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(0,max(sweed_finn$Likelihood[which(sweed_finn$Chromosome ==17)])+10),
                  expand = F)+xlab("Ssa17 position")+ylab("CLR")



min17<- min(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==17)],na.rm=T)
max17<- max(tajimad_finn$TajimaD[which(tajimad_finn$Chromosome ==17)],na.rm=T)


taji_17<-ggplot()+  
  geom_rect(aes(xmin=48500001, xmax=49500000, ymin=min17-0.2, ymax=max17+0.2), alpha=0.2)+
  
  geom_smooth(data=tajimad_finn[which(tajimad_finn$Chromosome ==17),], stat="identity",
              aes(y=TajimaD ,x= BIN_START), col="blue")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(min17-0.2, max17+0.2), expand = F)+xlab("")+ylab("Tajima's D")




##SIZE OF FIGURE
## PDF 12.24 x 11.45



plotssa17<- introgress_17 /
  taji_17 /
  sweed_17 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))



ggsave(plot=plotssa17,
       filename = "Ssa17_introgress_sel_plots.pdf",
       width = 12.24,
       height = 11.45,
       units = "in",
       dpi = 300,
       bg = NULL)

ggsave(plot=plotssa17, 
       filename = "Ssa17_introgress_sel_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)




################





##FST plots


#This script plots the gene annotation information (downloaded from NCBI) with the q-values from pcadapt
#SNPs that are highly significant (low q-values) were determined in this script
#These SNPs were highlighted in the figure along with the gene that they are found within
#This script could be improved :) 


#All annotations from genome - downloaded all genes for Salmo salar from NCBI in .csv format
annot <- data.table::fread("/Volumes/Accelsior_4TB/Salmon/TransAtlanticPaper/Recombination/ALL_ANNOTAIONS.csv")


####
min(compare$Diff_Fst, na.rm=T)
max(compare$Diff_Fst, na.rm=T)

min(compare$Diff_Fst[which(compare$CHR.x=="ssa06")], na.rm=T)
max(compare$Diff_Fst[which(compare$CHR.x=="ssa06")], na.rm=T)


fst_ssa06_all <-  ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa06"),],
                    aes(y=Diff_Fst,x=POS.x))+
  geom_point(data=compare[which(compare$CHR.x=="ssa06" & compare$POS.x >= 63000001 & compare$POS.x <= 64000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_point(data=compare[which(compare$CHR.x=="ssa06" & compare$POS.x >= 49000001 & compare$POS.x <= 50000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_point(data=compare[which(compare$CHR.x=="ssa06" & compare$POS.x >= 32500001 & compare$POS.x <= 33500000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
              geom_hline(yintercept = 0, col="blue", lty=2)+ylab("Difference FST")+xlab("Ssa06 position")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
NULL  

#get gene annotations
#identify regions in FST plot first and then extract data - visually assess FST regions and then use those regions here for saving/plotting

ssa06_mbp1 <- annot[which(annot$chromosome=="ssa06" & annot$start_position_on_the_genomic_accession>=32800000 & 
                           annot$end_position_on_the_genomic_accession< 33350000),]
write.table(ssa06_mbp1, "Genes_FST_region_Ssa06_1.txt", quote = F, row.names = F, col.names = T, sep="\t")


###

fst_ssa06_zoom1<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa06"),],
                                     aes(y=Diff_Fst,x=POS.x))+
  theme_classic()+theme(legend.position="none", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.7,0.6), xlim=c(31500001,34500000),
                  expand = F)+
  geom_point(data=compare[which(compare$CHR.x=="ssa06" & compare$POS.x >= 32500001 & compare$POS.x <= 33500000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  
  geom_hline(yintercept = 0, col="red", lty=2)+ylab("Difference FST")+xlab("Ssa06 position")+
  
  geom_vline(xintercept =32800000, lty=2 )+
  geom_vline(xintercept =33350000, lty=2)+
  geom_hline(yintercept = -0.6, col="gray70")+
  geom_rect(data=ssa06_mbp1, aes(xmin=ssa06_mbp1$start_position_on_the_genomic_accession,
                           xmax=ssa06_mbp1$end_position_on_the_genomic_accession, 
                          ymin= rep(-0.55, nrow(ssa06_mbp1)), ymax=rep(-0.65, nrow(ssa06_mbp1))), 
     fill="gray90", col="gray70")+
  NULL  


ssa06_mbp2 <- annot[which(annot$chromosome=="ssa06" & annot$start_position_on_the_genomic_accession>=49600001 & 
                            annot$end_position_on_the_genomic_accession< 50000000),]
write.table(ssa06_mbp2, "Genes_FST_region_Ssa06_2.txt", quote = F, row.names = F, col.names = T, sep="\t")


fst_ssa06_zoom2<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa06"),],
                                     aes(y=Diff_Fst,x=POS.x))+
  theme_classic()+theme(legend.position="none", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.7,0.6), xlim=c(48000001,51000000),
                  expand = F)+
  geom_point(data=compare[which(compare$CHR.x=="ssa06" & compare$POS.x >= 49000001 & compare$POS.x <= 50000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  
  geom_hline(yintercept = 0, col="red", lty=2)+ylab("Difference FST")+xlab("Ssa06 position")+
  
  geom_vline(xintercept =49600001, lty=2 )+
  geom_vline(xintercept =50000000, lty=2)+
  geom_hline(yintercept = -0.6, col="gray70")+
  geom_rect(data=ssa06_mbp2, aes(xmin=ssa06_mbp2$start_position_on_the_genomic_accession,
                       xmax=ssa06_mbp2$end_position_on_the_genomic_accession, 
                    ymin= rep(-0.55, nrow(ssa06_mbp2)), ymax=rep(-0.65, nrow(ssa06_mbp2))), 
           fill="gray90", col="gray70")+
  NULL  


ssa06_mbp3 <- annot[which(annot$chromosome=="ssa06" & annot$start_position_on_the_genomic_accession>=62900001 & 
                            annot$end_position_on_the_genomic_accession< 63600000),]
write.table(ssa06_mbp3, "Genes_FST_region_Ssa06_3.txt", quote = F, row.names = F, col.names = T, sep="\t")


fst_ssa06_zoom3 <- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa06"),],
                                     aes(y=Diff_Fst,x=POS.x))+
  theme_classic()+theme(legend.position="none", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.7,0.6), xlim=c(62000001,65000000),
                  expand = F)+
  geom_point(data=compare[which(compare$CHR.x=="ssa06" & compare$POS.x >= 63000001 & compare$POS.x <= 64000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  
  geom_hline(yintercept = 0, col="red", lty=2)+ylab("Difference FST")+xlab("Ssa06 position")+
  
  geom_vline(xintercept =62900001, lty=2 )+
  geom_vline(xintercept =63600000, lty=2)+
  geom_hline(yintercept = -0.6, col="gray70")+
  geom_rect(data=ssa06_mbp3, aes(xmin=ssa06_mbp3$start_position_on_the_genomic_accession,
                                xmax=ssa06_mbp3$end_position_on_the_genomic_accession, 
                                ymin= rep(-0.55, nrow(ssa06_mbp3)), ymax=rep(-0.65, nrow(ssa06_mbp3))), 
            fill="gray90", col="gray70")+
  NULL  




plot_fst_06<- fst_ssa06_all /
  (fst_ssa06_zoom1 +fst_ssa06_zoom2 +fst_ssa06_zoom3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))


ggsave(plot=plot_fst_06, 
       filename = "Ssa06_fst_plots.png",
       width = 16.81,
       height = 10.94, units="in",
       device='png', dpi=320,
       bg = NULL)

####################################################################

##plots for Ssa07

min(compare$Diff_Fst[which(compare$CHR.x=="ssa07")], na.rm=T)
max(compare$Diff_Fst[which(compare$CHR.x=="ssa07")], na.rm=T)


fst_ssa07_all<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa07"),],
                    aes(y=Diff_Fst,x=POS.x))+
  geom_point(data=compare[which(compare$CHR.x=="ssa07" & compare$POS.x >= 45000001 & compare$POS.x <= 46000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_hline(yintercept = 0, col="blue", lty=2)+ylab("Difference FST")+xlab("")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  NULL  


#Next add  gene annotations (using geom_rect)
#And highlight genes of interest (based on looking at row in dataframe that belongs to the gene of interest) - this could be improved :) 
#In this case the gene is in row 5 and 77 in the dataframe ssa13_mbp
#Added gene label manual (see "label" in aes) - looked up gene name to find abbreviation


#Save gene annotations for table
ssa07_mbp <- annot[which(annot$chromosome=="ssa07" & annot$start_position_on_the_genomic_accession>=45600001 & 
                           annot$end_position_on_the_genomic_accession< 45900000),]
write.table(ssa07_mbp, "Genes_FST_region_Ssa07.txt", quote = F, row.names = F, col.names = T, sep="\t")


fst_ssa07_zoom<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa07"),],
                    aes(y=Diff_Fst,x=POS.x))+
  theme_classic()+theme(legend.position="none", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.7,0.6), xlim=c(44500001,46500000),
                  expand = F)+
  geom_point(data=compare[which(compare$CHR.x=="ssa07" & compare$POS.x >= 45000001 & compare$POS.x <= 46000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  
  geom_hline(yintercept = 0, col="red", lty=2)+ylab("Difference FST")+xlab("Ssa07 position")+
  
  geom_vline(xintercept =45600001, lty=2 )+
  geom_vline(xintercept =45900000, lty=2)+
  geom_hline(yintercept = -0.6, col="gray70")+
    geom_rect(data=ssa07_mbp, aes(xmin=ssa07_mbp$start_position_on_the_genomic_accession,
                                xmax=ssa07_mbp$end_position_on_the_genomic_accession, 
                                ymin= rep(-0.55, nrow(ssa07_mbp)), ymax=rep(-0.65, nrow(ssa07_mbp))), 
              fill="gray90", col="gray70")+
  NULL  




plot_fst_07<- fst_ssa07_all /
  fst_ssa07_zoom +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))


ggsave(plot=plot_fst_07, 
       filename = "Ssa07_fst_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)



####################################################################

##plots for Ssa17

min(compare$Diff_Fst[which(compare$CHR.x=="ssa17")], na.rm=T)
max(compare$Diff_Fst[which(compare$CHR.x=="ssa17")], na.rm=T)


fst_ssa17_all<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa17"),],
                                    aes(y=Diff_Fst,x=POS.x))+
  geom_point(data=compare[which(compare$CHR.x=="ssa17" & compare$POS.x >= 48500001 & compare$POS.x <= 49500000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_hline(yintercept = 0, col="blue", lty=2)+ylab("Difference FST")+xlab("")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  NULL  


#Next add  gene annotations (using geom_rect)
#And highlight genes of interest (based on looking at row in dataframe that belongs to the gene of interest) - this could be improved :) 
#In this case the gene is in row 5 and 77 in the dataframe ssa13_mbp
#Added gene label manual (see "label" in aes) - looked up gene name to find abbreviation


#Save gene annotations for table
ssa17_mbp <- annot[which(annot$chromosome=="ssa17" & annot$start_position_on_the_genomic_accession>=48500001 & 
                           annot$end_position_on_the_genomic_accession< 49000000),]
write.table(ssa17_mbp, "Genes_FST_region_Ssa17.txt", quote = F, row.names = F, col.names = T, sep="\t")


fst_ssa17_zoom<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa17"),],
                                     aes(y=Diff_Fst,x=POS.x))+
  theme_classic()+theme(legend.position="none", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.7,0.6), xlim=c(47500001,50500000),
                  expand = F)+
  geom_point(data=compare[which(compare$CHR.x=="ssa17" & compare$POS.x >= 48500001 & compare$POS.x <= 49500000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  
  geom_hline(yintercept = 0, col="red", lty=2)+ylab("Difference FST")+xlab("Ssa17 position")+
  
  geom_vline(xintercept =48500001, lty=2 )+
  geom_vline(xintercept =49000000, lty=2)+
  geom_hline(yintercept = -0.6, col="gray70")+
 geom_rect(data=ssa17_mbp, aes(xmin=ssa17_mbp$start_position_on_the_genomic_accession,
                          xmax=ssa17_mbp$end_position_on_the_genomic_accession, 
                          ymin= rep(-0.55, nrow(ssa17_mbp)), ymax=rep(-0.65, nrow(ssa17_mbp))), 
       fill="gray90", col="gray70")+
  NULL  




plot_fst_17<- fst_ssa17_all /
  fst_ssa17_zoom +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))


ggsave(plot=plot_fst_17, 
       filename = "Ssa17_fst_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)



####################################################################

##plots for Ssa18

min(compare$Diff_Fst[which(compare$CHR.x=="ssa18")], na.rm=T)
max(compare$Diff_Fst[which(compare$CHR.x=="ssa18")], na.rm=T)


fst_ssa18_all<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa18"),],
                                    aes(y=Diff_Fst,x=POS.x))+
  geom_point(data=compare[which(compare$CHR.x=="ssa18" & compare$POS.x >= 24500001 & compare$POS.x <= 25500000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_hline(yintercept = 0, col="blue", lty=2)+ylab("Difference FST")+xlab("")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  NULL  


#Next add  gene annotations (using geom_rect)
#And highlight genes of interest (based on looking at row in dataframe that belongs to the gene of interest) - this could be improved :) 
#In this case the gene is in row 5 and 77 in the dataframe ssa13_mbp
#Added gene label manual (see "label" in aes) - looked up gene name to find abbreviation


#Save gene annotations for table
ssa18_mbp <- annot[which(annot$chromosome=="ssa18" & annot$start_position_on_the_genomic_accession>=24500001 & 
                           annot$end_position_on_the_genomic_accession< 25050000),]
write.table(ssa18_mbp, "Genes_FST_region_Ssa18.txt", quote = F, row.names = F, col.names = T, sep="\t")


fst_ssa18_zoom<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa18"),],
                                     aes(y=Diff_Fst,x=POS.x))+
  theme_classic()+theme(legend.position="none", 
                        legend.background =  element_rect(color=NULL),
                        legend.text = element_text(size=15),
                        legend.title = element_blank(),
                        axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.8,0.6), xlim=c(23500001,26500000),
                  expand = F)+
  geom_point(data=compare[which(compare$CHR.x=="ssa18" & compare$POS.x >= 24500001 & compare$POS.x <= 25500000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  
  geom_hline(yintercept = 0, col="red", lty=2)+ylab("Difference FST")+xlab("Ssa18 position")+
  
  geom_vline(xintercept =24500001, lty=2 )+
  geom_vline(xintercept =25050000, lty=2)+
  geom_hline(yintercept = -0.7, col="gray70")+
  geom_rect(data=ssa18_mbp, aes(xmin=ssa18_mbp$start_position_on_the_genomic_accession,
                        xmax=ssa18_mbp$end_position_on_the_genomic_accession, 
                    ymin= rep(-0.65, nrow(ssa18_mbp)), ymax=rep(-0.75, nrow(ssa18_mbp))), 
            fill="gray90", col="gray70")+
  NULL  




plot_fst_18<- fst_ssa18_all /
  fst_ssa18_zoom +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))


ggsave(plot=plot_fst_18, 
       filename = "Ssa18_fst_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)






######### Plot Ssa27
min(compare$Diff_Fst[which(compare$CHR.x=="ssa27")], na.rm=T)
max(compare$Diff_Fst[which(compare$CHR.x=="ssa27")], na.rm=T)


fst_ssa27_all <- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa27"),],
                    aes(y=Diff_Fst,x=POS.x))+
  geom_point(data=compare[which(compare$CHR.x=="ssa27" & compare$POS.x >= 20500001 & compare$POS.x <= 21500000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_point(data=compare[which(compare$CHR.x=="ssa27" & compare$POS.x >= 41000001 & compare$POS.x <= 42000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  
  geom_hline(yintercept = 0, col="blue", lty=2)+ylab("Difference FST")+xlab("Ssa27 position")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  NULL  




#####ssa27 annoations

ssa27_mbp=annot[which(annot$chromosome=="ssa27" & annot$start_position_on_the_genomic_accession>=41600001 & 
                        annot$end_position_on_the_genomic_accession< 42200000),]

write.table(ssa27_mbp, "Genes_FST_region_Ssa27_1.txt", quote = F, row.names = F, col.names = T, sep="\t")




#Next add  gene annotations (using geom_rect)
#And highlight genes of interest (based on looking at row in dataframe that belongs to the gene of interest) - this could be improved :) 
#In this case the gene is in row 5 and 77 in the dataframe ssa13_mbp
#Added gene label manual (see "label" in aes) - looked up gene name to find abbreviation

fst_ssa27_zoom<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa27"),],
                                     aes(y=Diff_Fst,x=POS.x))+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.8,0.6), xlim=c(40500001,42500000),
                  expand = F)+
  geom_point(data=compare[which(compare$CHR.x=="ssa27" & compare$POS.x >= 41000001 & compare$POS.x <= 42000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
    geom_hline(yintercept = 0, col="red", lty=2)+ylab("Difference FST")+xlab("ssa27 position")+
  
  geom_vline(xintercept =41600001, lty=2 )+
  geom_vline(xintercept =42200000, lty=2)+
  geom_hline(yintercept = -0.7, col="gray70")+
  geom_rect(data=ssa27_mbp, aes(xmin=ssa27_mbp$start_position_on_the_genomic_accession,
                                xmax=ssa27_mbp$end_position_on_the_genomic_accession, 
                                ymin= rep(-0.75, nrow(ssa27_mbp)), ymax=rep(-0.65, nrow(ssa27_mbp))), fill="gray90", col="gray70")+
  
  NULL  


#Annotation for Ssa13
ssa27_mbp2=annot[which(annot$chromosome=="ssa27" & annot$start_position_on_the_genomic_accession>=20550001 & 
                        annot$end_position_on_the_genomic_accession< 21500000),]

write.table(ssa27_mbp2, "Genes_FST_region_Ssa27_2.txt", quote = F, row.names = F, col.names = T, sep="\t")




#Next add  gene annotations (using geom_rect)
#And highlight genes of interest (based on looking at row in dataframe that belongs to the gene of interest) - this could be improved :) 
#In this case the gene is in row 5 and 77 in the dataframe ssa13_mbp
#Added gene label manual (see "label" in aes) - looked up gene name to find abbreviation

fst_ssa27_zoom2 <- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa27"),],
                                     aes(y=Diff_Fst,x=POS.x))+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.8,0.6), xlim=c(19500001,22500000),
                  expand = F)+
  geom_point(data=compare[which(compare$CHR.x=="ssa27" & compare$POS.x >= 20500001 & compare$POS.x <= 21500000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_hline(yintercept = 0, col="red", lty=2)+ylab("Difference FST")+xlab("ssa27 position")+
  
  geom_vline(xintercept =20550001, lty=2 )+
  geom_vline(xintercept =21500000, lty=2)+
  geom_hline(yintercept = -0.7, col="gray70")+
  geom_rect(data=ssa27_mbp2, aes(xmin=ssa27_mbp2$start_position_on_the_genomic_accession,
                                xmax=ssa27_mbp2$end_position_on_the_genomic_accession, 
                                ymin= rep(-0.75, nrow(ssa27_mbp2)), ymax=rep(-0.65, nrow(ssa27_mbp2))), fill="gray90", col="gray70")+
  
  NULL  





plot_fst_27 <- fst_ssa27_all /
 (fst_ssa27_zoom2 + fst_ssa27_zoom ) +
plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))


ggsave(plot=plot_fst_27, 
       filename = "Ssa27_fst_plots.png",
       width = 12.24,
       height = 11.45, units="in",
       device='png', dpi=320,
       bg = NULL)





######### Plots for Ssa21 - for main text
min(compare$Diff_Fst[which(compare$CHR.x=="ssa21")], na.rm=T)
max(compare$Diff_Fst[which(compare$CHR.x=="ssa21")], na.rm=T)



##plots for Ssa21

fst_ssa21_all<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa21"),],
                                    aes(y=Diff_Fst,x=POS.x))+
  geom_point(data=compare[which(compare$CHR.x=="ssa21" & compare$POS.x >= 34000001 & compare$POS.x <= 35000000),],
             aes(y=Diff_Fst,x=POS.x), col="red")+
  geom_hline(yintercept = 0, col="blue", lty=2)+ylab("Difference FST")+xlab("ssa21 position")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  NULL  





#####Ssa21 annoations
#Annotation for Ssa13
ssa21_mbp=annot[which(annot$chromosome=="ssa21" & annot$start_position_on_the_genomic_accession>=34750001 & 
                        annot$end_position_on_the_genomic_accession< 35120000),]

write.table(ssa21_mbp, "Genes_FST_region_Ssa21.txt", quote = F, row.names = F, col.names = T, sep="\t")



#Next add  gene annotations (using geom_rect)
#And highlight genes of interest (based on looking at row in dataframe that belongs to the gene of interest) - this could be improved :) 
#In this case the gene is in row 5 and 77 in the dataframe ssa13_mbp
#Added gene label manual (see "label" in aes) - looked up gene name to find abbreviation

fst_ssa21_zoom<- ggplot()+geom_point(data=compare[which(compare$CHR.x=="ssa21"),],
                                     aes(y=Diff_Fst,x=POS.x), pch=21, fill="black", col="black")+
  theme_classic()+theme(legend.position = "none", axis.title = element_text(size=15),  axis.text = element_text(size=12))+
  coord_cartesian(ylim = c(-0.8,0.6), xlim=c(34400001,36000000),
                  expand = F)+
  
  geom_hline(yintercept = 0, col="red", lty=2)+ ylab(expression('Difference  '*italic(F)[ST]*' '))+
    
 xlab("Ssa21 position")+
  
  geom_vline(xintercept =34750001, lty=2 )+
  geom_vline(xintercept =35120000, lty=2)+
  geom_hline(yintercept = -0.7, col="gray70")+
  geom_rect(data=ssa21_mbp, aes(xmin=ssa21_mbp$start_position_on_the_genomic_accession,
                                xmax=ssa21_mbp$end_position_on_the_genomic_accession, 
                                ymin= rep(-0.75, nrow(ssa21_mbp)), ymax=rep(-0.65, nrow(ssa21_mbp))), fill="gray90", col="gray70")+
  
  NULL  






fst_ssa21_genes<- ggplot() +
  geom_hline(yintercept = 0, col="gray70")+
  geom_rect(data=ssa21_mbp, aes(xmin=ssa21_mbp$start_position_on_the_genomic_accession,
                                xmax=ssa21_mbp$end_position_on_the_genomic_accession, 
                                ymin= rep(-0.15, nrow(ssa21_mbp)), ymax=rep(0.15, nrow(ssa21_mbp))), fill="gray90", col="gray70")+
  ylim(c(-0.6,2))+
  annotate(geom = "text", x = ssa21_mbp$start_position_on_the_genomic_accession, y=0.02, hjust = 0,
           label=ssa21_mbp$description , size=4, angle=45)+theme_bw()+theme(panel.grid = element_blank())


plot_fst_21_main<- introgress_21 /
  sweed_21 /
  fst_ssa21_zoom /
  fst_ssa21_genes +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))


ggsave(plot=plot_fst_21_main, 
       filename = "Ssa21_MainText_plots.png",
       width = 11.65,
       height = 15, units="in",
       device='png', dpi=320,
       bg = NULL)


ggsave(plot=plot_fst_21_main, 
       filename = "Ssa21_MainText_plots.pdf",
       width = 11.65,
       height = 15, units="in",
       device='pdf', dpi=320,
       bg = NULL)


  


