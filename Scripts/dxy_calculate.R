#Libraries
library(ggplot2)
library(PopGenome)
library(stringr)

#Set directory
setwd("~/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/")

#Read map positions
map=data.table::fread("ssa_combined_wgs_aquagenome_finnmark_canada_for_dxy.bim", header=F)
map=map[,1:4]
colnames(map)=c("Chr", "SNP", "cm", "Pos")
map$Pos=as.numeric(as.character(map$Pos))

###Chromosome  Dxy
fam=read.table("ssa_combined_wgs_aquagenome_finnmark_canada_for_dxy.fam", header=F)
fam$intersect=interaction(fam$V1, fam$V2, sep = "_")

EU_POP <-as.character(fam$intersect[which(fam$V1=="Finnmark")])
NA_POP <- as.character(fam$intersect[which(fam$V1=="Newfoundland" | fam$V1=="Labrador"| fam$V1=="QC" | fam$V1=="QC_northeast")])

#loop to read in chromosomes and get dxy
for(i in 1:9){
  
  ALL_CHR_ssa <- readVCF("ssa_combined_wgs_aquagenome_finnmark_canada_for_dxy.vcf.gz", 1000,tid = paste0("ssa0",i),
                           min(map$Pos[which(map$Chr==paste0("ssa0",i))]),
                           max(map$Pos[which(map$Chr==paste0("ssa0",i))]),
                           include.unknown = T)
  
  ALL_CHR_ssa@populations
  ALL_CHR_ssa@n.sites
  ALL_CHR_ssa@n.biallelic.sites
  ALL_CHR_ssa@region.data
  get.individuals(ALL_CHR_ssa)
  
  
  ALL_CHR_ssa <- set.populations(ALL_CHR_ssa, list(NA_POP, EU_POP), diploid = T)
  ALL_CHR_ssa@populations
  get.individuals(ALL_CHR_ssa)
  
  
  ALL_CHR_Slide <- sliding.window.transform(ALL_CHR_ssa, width = 10000,
                                              jump = 10000, type = 2)
  
  
  ALL_CHR_Slide_genome.pos <- sapply(ALL_CHR_Slide@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    8
    val   <- mean(as.numeric(split))
    return(val)
  })
  
  ALL_CHR_Slide <- diversity.stats.between(ALL_CHR_Slide)
  ALL_CHR_Slide@nuc.diversity.between
  
  ALL_CHR_dxy <- as.data.frame(cbind("dxy" = ALL_CHR_Slide@nuc.diversity.between, 
                                       "Position" = as.integer(as.character(ALL_CHR_Slide_genome.pos))))
  colnames(ALL_CHR_dxy) <- c("dxy", "Position")
  ALL_CHR_dxy$Chr=rep(paste0("ssa0",i))
  
  write.table(ALL_CHR_dxy, paste0("Ssa0",i, "_dxy_finnmark_canada.txt"), quote = F, row.names = F, col.names = T, sep="\t")
  print(paste0("Done ssa0", i))
  
}



for(i in 10:29){
  
  ALL_CHR_ssa <- readVCF("ssa_combined_wgs_aquagenome_finnmark_canada_for_dxy.vcf.gz", 1000,tid = paste0("ssa",i),
                           min(map$Pos[which(map$Chr==paste0("ssa",i))]),
                           max(map$Pos[which(map$Chr==paste0("ssa",i))]),
                           include.unknown = T)
  
  ALL_CHR_ssa@populations
  ALL_CHR_ssa@n.sites
  ALL_CHR_ssa@n.biallelic.sites
  ALL_CHR_ssa@region.data
  get.individuals(ALL_CHR_ssa)
  
  
  ALL_CHR_ssa <- set.populations(ALL_CHR_ssa, list(NA_POP, EU_POP), diploid = T)
  ALL_CHR_ssa@populations
  get.individuals(ALL_CHR_ssa)
  
  
  ALL_CHR_Slide <- sliding.window.transform(ALL_CHR_ssa, width = 10000,
                                              jump = 10000, type = 2)
  
  
  ALL_CHR_Slide_genome.pos <- sapply(ALL_CHR_Slide@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    8
    val   <- mean(as.numeric(split))
    return(val)
  })
  
  ALL_CHR_Slide <- diversity.stats.between(ALL_CHR_Slide)
  ALL_CHR_Slide@nuc.diversity.between
  
  ALL_CHR_dxy <- as.data.frame(cbind("dxy" = ALL_CHR_Slide@nuc.diversity.between, 
                                       "Position" = as.integer(as.character(ALL_CHR_Slide_genome.pos))))
  colnames(ALL_CHR_dxy) <- c("dxy", "Position")
  ALL_CHR_dxy$Chr=rep(paste0("ssa",i))
  write.table(ALL_CHR_dxy, paste0("Ssa",i, "_dxy_finnmark_canada.txt"), quote = F, row.names = F, col.names = T, sep="\t")
  
  print(paste0("Done ssa", i))
}





##########################################################################################
##################################


## Southern NORWAY

fam2=read.table("ssa_combined_wgs_aquagenome_southernnorway_canada_for_dxy.fam", header=F)
fam2$intersect=interaction(fam2$V1, fam2$V2, sep = "_")

EU_POP_South <-as.character(fam2$intersect[which(fam2$V1=="SouthernNorway")])
NA_POP <- as.character(fam2$intersect[which(fam2$V1=="Newfoundland" | fam2$V1=="QC" | fam2$V1=="QC_northeast"| fam2$V1=="Labrador")])



for(i in 1:9){
  
  ALL_CHR_ssa_southnor <- readVCF("ssa_combined_wgs_aquagenome_southernnorway_canada_for_dxy.vcf.gz", 1000, tid =  paste0("ssa0",i), 
                                    min(map$Pos[which(map$Chr==paste0("ssa0",i))]),
                                    max(map$Pos[which(map$Chr==paste0("ssa0",i))]),
                                    include.unknown = T)
  
  ALL_CHR_ssa_southnor@populations
  ALL_CHR_ssa_southnor@n.sites
  ALL_CHR_ssa_southnor@n.biallelic.sites
  ALL_CHR_ssa_southnor@region.data
  get.individuals(ALL_CHR_ssa_southnor)
  
  ALL_CHR_ssa_snor <- set.populations(ALL_CHR_ssa_southnor, list(NA_POP, EU_POP_South), diploid = T)
  ALL_CHR_ssa_snor@populations
  get.individuals(ALL_CHR_ssa_snor)
  
  
  ALL_CHR_Slide_southNor <- sliding.window.transform(ALL_CHR_ssa_snor, width = 10000,
                                                       jump = 10000, type = 2)
  
  
  ALL_CHR_Slide_genome.pos_southNor <- sapply(ALL_CHR_Slide_southNor@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    8
    val   <- mean(as.numeric(split))
    return(val)
  })
  
  ALL_CHR_Slide_southNor <- diversity.stats.between(ALL_CHR_Slide_southNor)
  
  ALL_CHR_Slide_southNor@nuc.diversity.between
  
  ALL_CHR_dxy_southNor <- as.data.frame(cbind("dxy" = ALL_CHR_Slide_southNor@nuc.diversity.between, 
                                                "Position" = as.integer(as.character(ALL_CHR_Slide_genome.pos_southNor))))
  colnames(ALL_CHR_dxy_southNor) <- c("dxy", "Position")
  
  
  ALL_CHR_dxy_southNor$Chr=rep(paste0("ssa0",i))
  write.table(ALL_CHR_dxy_southNor, paste0("Ssa0",i, "_dxy_southernnorway_canada.txt"), quote = F, row.names = F, col.names = T, sep="\t")
  
  print(paste0("SNor - Done ssa0", i))
  
  
}




for(i in 10:29){
  ALL_CHR_ssa_southnor <- readVCF("ssa_combined_wgs_aquagenome_southernnorway_canada_for_dxy.vcf.gz", 1000, tid =  paste0("ssa",i), 
                                    min(map$Pos[which(map$Chr==paste0("ssa",i))]),
                                    max(map$Pos[which(map$Chr==paste0("ssa",i))]),
                                    include.unknown = T)
    
  ALL_CHR_ssa_southnor@populations
  ALL_CHR_ssa_southnor@n.sites
  ALL_CHR_ssa_southnor@n.biallelic.sites
  ALL_CHR_ssa_southnor@region.data
  get.individuals(ALL_CHR_ssa_southnor)
  
  ALL_CHR_ssa_snor <- set.populations(ALL_CHR_ssa_southnor, list(NA_POP, EU_POP_South), diploid = T)
  ALL_CHR_ssa_snor@populations
  get.individuals(ALL_CHR_ssa_snor)
  
  
  ALL_CHR_Slide_southNor <- sliding.window.transform(ALL_CHR_ssa_snor, width = 10000,
                                                       jump = 10000, type = 2)
  
  
  ALL_CHR_Slide_genome.pos_southNor <- sapply(ALL_CHR_Slide_southNor@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    8
    val   <- mean(as.numeric(split))
    return(val)
  })
  
  ALL_CHR_Slide_southNor <- diversity.stats.between(ALL_CHR_Slide_southNor)
  
  ALL_CHR_Slide_southNor@nuc.diversity.between
  
  ALL_CHR_dxy_southNor <- as.data.frame(cbind("dxy" = ALL_CHR_Slide_southNor@nuc.diversity.between, 
                                                "Position" = as.integer(as.character(ALL_CHR_Slide_genome.pos_southNor))))
  colnames(ALL_CHR_dxy_southNor) <- c("dxy", "Position")
  
  
  ALL_CHR_dxy_southNor$Chr=rep(paste0("ssa",i))
  write.table(ALL_CHR_dxy_southNor, paste0("Ssa",i, "_dxy_southernnorway_canada.txt"), quote = F, row.names = F, col.names = T, sep="\t")
  
  print(paste0("SNor - Done ssa", i))
  
  
}
