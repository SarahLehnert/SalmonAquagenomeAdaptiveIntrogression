
Notes for running analyses for WGS aquagenome project
Adaptive Introgression analyses - SL

This workflow includes initial filtering and Treemix threepop analysis.
Scripts could be updated for efficiency, written awhile ago...

Initial filtering step:
```
#First filter by maf 0.01 -- sites have already been filtered for missingness
#Chose 0.01 as did not want to filter out some rare alleles if present in small populations (White Sea, Landlocked)

vcftools --vcf /filepath/WGS_Aquagenome/vcf_files/Final_Combined_VCF/ssa_combined_wgs_aquagenome_biallelic_PASS.recode.vcf --maf 0.01 --recode --out /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001.recode.vcf
```

Next use Beagle to phase and impute data

```
#With filtered (maf 0.01) VCF (SNPs  = 2237598)
#Phase and impute in Beagle - use all samples together

java -Xmx300G -jar beagle.18May20.d20.jar gt=/filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001.recode.vcf.recode.vcf out=/filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed  nthreads=30
```

File manipulation in plink to prepare data for Treemix - this probably could have been done more efficiently

```
##After phasing/imputing 
#Convert to plink to get allele frequencies for treemix - a few steps here to update pop IDs for geographic regions

cd  ~/Desktop/Software/plink_mac_20200219

./plink --vcf /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed.vcf --make-bed --allow-extra-chr --double-id --out /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink

#Then update names of individuals/populations

./plink --bfile /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink --make-bed --allow-extra-chr --update-ids /filepath/Salmon/WGS_Aquagenome/vcf_files/Update_IDs.txt --out /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids

#Then update IDs with geographic regions

./plink --bfile /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids --make-bed --allow-extra-chr --update-ids /filepath/Salmon/WGS_Aquagenome/vcf_files/Update_IDs3_cluster_region.txt --out /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION

##Copy files to treemix folder
#Get allele frequency file for treemix
##Should do maf>0.01 again after imputing/phasing just in case (just remove 710 loci)
./plink --bfile /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION --maf 0.01 --allow-extra-chr --make-bed --out /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing

./plink --bfile /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing --freq --family --allow-extra-chr --out /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing_allelefreq

```
Prepare file for Treemix, including gzip and running plink2treemix script

```
Gzip frequency file
gzip /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing_allelefreq.frq.strat

############################
#Move to treemix directory 
#run plink2treemix script to convert format

./plink2treemix.py /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing_allelefreq.frq.strat.gz WGS_Salmo_for_threepop_REGIONS_jan_2025.gz
```
Run treepop test in treemix 

```
#Use 500 SNP windows
#this will take a while -- no way to multithread??? 
#note that .gz file was saved in treemix/src folder
./threepop -i WGS_Salmo_for_threepop_REGIONS_jan_2025.gz -k 500 > results_threepop_WGS_salmon_by_REGION_Jan_2025

```
Save Treemix results 
```
############################
#save output (remove extra parts that are not needed in text)
cat results_threepop_WGS_salmon_by_REGION_Jan_2025|grep -v Estimating |grep -v nsnp|tr '' ' ' > WGS_salmonpops_results_threepop_final_by_Region_500snp_window_jan2025
```
