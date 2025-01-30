
Workflow for setting up files for PCAdmix and running analyses (again, this probably could have been done more efficiently)

First use plink to get list of SNPs and IDs in final dataset to run analyses in vcftools:

```
./plink --bfile /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing --allow-extra-chr --make-bed --out /filepath/Salmon/WGS_Aquagenome/Pcadmix/all_snps_all_inds
```
Next - use vcf tools to subset first. Note that data were initially filtered for maf 0.01, then imputed/phased in beagle, then filtered one more time for maf 0.01

For PCAdmix require 2 pure groups and 1 admixed group: subset for these groups now
```
#Subset Canada (pure group)
vcftools --vcf /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed.vcf --snps /filepath/Salmon/WGS_Aquagenome/Pcadmix/SNP_list_all.txt --phased --keep /filepath/Salmon/WGS_Aquagenome/Pcadmix/Canada_IDs_original.txt --recode --out /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_phased.vcf 

#Subset East Atlantic (Southern Norway including landlocked) (pure group)
vcftools --vcf /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed.vcf --snps /filepath/Salmon/WGS_Aquagenome/Pcadmix/SNP_list_all.txt --phased --keep /filepath/Salmon/WGS_Aquagenome/Pcadmix/SouthernNorwayLL_IDs_original.txt --recode --out /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_phased.vcf 

#Subset Western Barents and White Sea group (admixed group)
vcftools --vcf /filepath/Salmon/WGS_Aquagenome/vcf_files/Final_Combined_VCF/Phase_Impute/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed.vcf --snps /filepath/Salmon/WGS_Aquagenome/Pcadmix/SNP_list_all.txt --phased --keep /filepath/Salmon/WGS_Aquagenome/Pcadmix/FinnmarkWhite_IDs_original.txt --recode --out /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_phased.vcf 
```

Next need to get map for each chromosome to specify what SNPs to use as chromosomes are run seperately

```

#generate a map file for pcadamix and do for each chromosome

./plink --bfile /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing --recode --allow-extra-chr --out /filepath/Salmon/WGS_Aquagenome/Pcadmix/all_snps_all_inds

for i in {1..9}; do ./plink --bfile /filepath/Salmon/WGS_Aquagenome/Pcadmix/all_snps_all_inds --chr ssa0${i} --allow-extra-chr --recode --out /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa0${i}all_snps_all_inds; done

for i in {10..29}; do ./plink --bfile /filepath/Salmon/WGS_Aquagenome/Pcadmix/all_snps_all_inds --chr ssa${i} --allow-extra-chr --recode --out /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa${i}all_snps_all_inds; done

```

Convert to format for  pcadmix using beagle

```

cd Desktop/Software/beagle

cat /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_phased.vcf.recode.vcf | java -jar vcf2beagle.jar -1 /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_converted_beagle

cat /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_phased.vcf.recode.vcf | java -jar vcf2beagle.jar -1 /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_converted_beagle

cat /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_phased.vcf.recode.vcf | java -jar vcf2beagle.jar -1 /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_converted_beagle

```
Run PCAdmix, here we used 100 SNP windows and ran each chromosome seperately. 

```
#need to run for each chromosome separately - doesn't delineate between different chromosomes 
#start in different terminal windows to use multiple cores (can't multithread?)
#note might be more efficient to subset vcf for each chr first - but would create a lot more files, so went with this...
#Can change window size to bp size instead - but then lots of variation in number of SNPs used
#no pruning was done here.

cd Desktop/Software/PCAdmix
for i in {1..5}; do ./PCAdmix3_macOSX -anc /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_converted_beagle.bgl /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_converted_beagle.bgl -adm /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_converted_beagle.bgl -map /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa0${i}all_snps_all_inds.map  -w 100 -prune 0 -ld 0 -o ssa0${i}_WGS_FinnmarkIntrogress_jan_2025; done

cd Desktop/Software/PCAdmix
for i in {6..9}; do ./PCAdmix3_macOSX -anc /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_converted_beagle.bgl /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_converted_beagle.bgl -adm /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_converted_beagle.bgl -map /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa0${i}all_snps_all_inds.map  -w 100 -prune 0 -ld 0 -o ssa0${i}_WGS_FinnmarkIntrogress_jan_2025; done

cd Desktop/Software/PCAdmix

for i in {10..15}; do ./PCAdmix3_macOSX -anc /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_converted_beagle.bgl /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_converted_beagle.bgl -adm /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_converted_beagle.bgl -map /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa${i}all_snps_all_inds.map  -w 100 -prune 0 -ld 0 -o ssa${i}_WGS_FinnmarkIntrogress_jan_2025; done 


cd Desktop/Software/PCAdmix

for i in {16..20}; do ./PCAdmix3_macOSX -anc /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_converted_beagle.bgl /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_converted_beagle.bgl -adm /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_converted_beagle.bgl -map /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa${i}all_snps_all_inds.map  -w 100 -prune 0 -ld 0 -o ssa${i}_WGS_FinnmarkIntrogress_jan_2025; done

cd Desktop/Software/PCAdmix

for i in {21..25}; do ./PCAdmix3_macOSX -anc /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_converted_beagle.bgl /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_converted_beagle.bgl -adm /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_converted_beagle.bgl -map /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa${i}all_snps_all_inds.map  -w 100 -prune 0 -ld 0 -o ssa${i}_WGS_FinnmarkIntrogress_jan_2025; done

cd Desktop/Software/PCAdmix

for i in {26..29}; do ./PCAdmix3_macOSX -anc /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_SOUTHERNNORWAY_converted_beagle.bgl /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_CANADA_converted_beagle.bgl -adm /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_converted_beagle.bgl -map /filepath/Salmon/WGS_Aquagenome/Pcadmix/ssa${i}all_snps_all_inds.map  -w 100 -prune 0 -ld 0 -o ssa${i}_WGS_FinnmarkIntrogress_jan_2025; done
```

After running PCadmix - results were imported into R for further analyses - see scripts


