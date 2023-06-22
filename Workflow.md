
Notes for running analyses for WGS aquagenome project
Adaptive Introgression analyses - SL

This workflow includes initial filtering and Treemix threepop analysis

Initial filtering step:
```
#First filter by maf 0.01 -- sites have already been filtered for missingness
#Chose 0.01 as did not want to filter out some rare alleles if present in small populations (White Sea, Landlocked)

vcftools --vcf /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/ssa_combined_wgs_aquagenome_biallelic_PASS.recode.vcf --maf 0.01 --recode --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001.recode.vcf
```

Next use Beagle to phase and impute data

```
#With filtered (maf 0.01) VCF (SNPs  = 2237598)
#Phase and impute in Beagle - use all samples together

java -Xmx300G -jar beagle.18May20.d20.jar gt=/Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001.recode.vcf.recode.vcf out=/Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed  nthreads=30
```

File manipulation in plink to prepare data for Treemix

```
##After phasing/imputing 
#Convert to plink to get allele frequencies for treemix

cd  ~/Desktop/Software/plink_mac_20200219

./plink --vcf /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed.vcf --make-bed --allow-extra-chr --double-id --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink

#Then update names of individuals/populations

./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink --make-bed --allow-extra-chr --update-ids /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Update_IDs.txt --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids

#Then update IDs with geographic regions

./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids --make-bed --allow-extra-chr --update-ids /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Update_IDs2_cluster_region.txt --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION

##Copy files to treemix folder
#Get allele frequency file for treemix

./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION --freq --family --allow-extra-chr --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_freq

##Should do maf>0.01 again after imputing/phasing? (just remove 710 loci)
./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION --maf 0.01 --allow-extra-chr --make-bed --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing


##Copy files to treemix folder
#Get allele frequency file for treemix

./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing --freq --family --allow-extra-chr --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing_allelefreq

```
Preparing data to run treemix

```

#Gzip frequency file
gzip /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_freq.frq.strat

############################
#Move to treemix directory 
cd ..
cd treemix-1.13/src/

############################
#run plink2treemix script to convert format

 ./plink2treemix.py /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_freq.frq.strat.gz WGS_Salmo_for_threepop_REGIONS_may_27_2020.gz

############################
#RUN THREEPOP test
#Use 500 SNP windows
#this will take a while -- no way to multithread??? 
#note that .gz file was saved in treemix/src folder
./threepop -i WGS_Salmo_for_threepop_REGIONS_may_27_2020.gz -k 500 > results_threepop_WGS_salmon_by_REGION_may_27_2020

############################
#save output (remove extra parts that are not needed in text)
cat results_threepop_WGS_salmon_by_REGION_may_27_2020|grep -v Estimating |grep -v nsnp|tr '' ' ' > WGS_salmonpops_results_threepop_final_by_Region_500snp_window

############################ Run with 1000 SNP windows
#RUN THREEPOP test
#Use 1000 SNP windows
#this will take a while -- no way to multithread??? 
#note that .gz file was saved in treemix/src folder
./threepop -i WGS_Salmo_for_threepop_REGIONS_may_27_2020.gz -k 1000 > results_threepop_WGS_salmon_by_REGION_may_28_2020_1000snp

############################
#save output (remove extra parts that are not needed in text)
cat results_threepop_WGS_salmon_by_REGION_may_28_2020_1000snp|grep -v Estimating |grep -v nsnp|tr '' ' ' > WGS_salmonpops_results_threepop_final_by_Region_1000snp_window



#################### TREEMIX  with maf 0.01 after phasing/imputing ######################


#Gzip frequency file
gzip /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing_allelefreq.frq.strat

############################
#Move to treemix directory 
cd ..
cd treemix-1.13/src/

############################
#run plink2treemix script to convert format

 ./plink2treemix.py /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing_allelefreq.frq.strat.gz WGS_Salmo_for_threepop_REGIONS_may_28_2020_MAF001_after.gz

############################
#RUN THREEPOP test
#Use 500 SNP windows
#this will take a while -- no way to multithread??? 
#note that .gz file was saved in treemix/src folder
./threepop -i WGS_Salmo_for_threepop_REGIONS_may_28_2020_MAF001_after.gz -k 500 > results_threepop_WGS_salmon_by_REGION_may_28_2020_maf001_after

############################
#save output (remove extra parts that are not needed in text)
cat results_threepop_WGS_salmon_by_REGION_may_28_2020_maf001_after|grep -v Estimating |grep -v nsnp|tr '' ' ' > WGS_salmonpops_results_threepop_final_by_Region_500snp_window_may28_2020_maf001_after


#################### TREEMIX  with maf 0.01 after phasing/imputing ######################
##By pop (not region) -- do all comparisons to see if there is any evidence for certain pops having contact

./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids --maf 0.01 --allow-extra-chr --make-bed --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_POP_maf_after_phasing

./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_POP_maf_after_phasing --freq --family --allow-extra-chr --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_POP_maf_after_phasing_allelefreq

######################  TREEMIX by POP - maf 0.01 after phasing #########################

#Gzip frequency file
gzip /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_POP_maf_after_phasing_allelefreq.frq.strat

############################
#Move to treemix directory 
cd ..
cd treemix-1.13/src/

############################
#run plink2treemix script to convert format

 ./plink2treemix.py /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_POP_maf_after_phasing_allelefreq.frq.strat.gz WGS_Salmo_for_threepop_POP_may_28_2020_maf001_after.gz

############################
#RUN THREEPOP test
#Use 500 SNP windows
#this will take a while -- no way to multithread??? 
#note that .gz file was saved in treemix/src folder
./threepop -i WGS_Salmo_for_threepop_POP_may_28_2020_maf001_after.gz -k 500 > results_threepop_WGS_salmon_by_POP_may_28_2020_maf001_after

############################
#save output (remove extra parts that are not needed in text)
cat results_threepop_WGS_salmon_by_POP_may_28_2020_maf001_after|grep -v Estimating |grep -v nsnp|tr '' ' ' > WGS_salmonpops_results_threepop_final_by_POP_500snp_window_maf001_after
