
Workflow for setting up files for Dxy (Dxy script was run in R - see scripts)

First subset files for 1) Canada and Barents White Sea (Finnmark + White sea), and 2) Canada and East Altantic (Southern Norway + Landlocked), seperately
This was done on the same file used for Treemix (maf>0.01 was filtered again after phasing/imputation). First use plink to convert to vcf, then vcf.gz

Create files for Canada and Barents/White Sea:

```
./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing --allow-extra-chr --keep-fam FinCan.txt --make-bed --recode vcf --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/ssa_combined_wgs_aquagenome_finnmark_canada_for_dxy

bgzip -c /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/ssa_combined_wgs_aquagenome_finnmark_canada_for_dxy.vcf > /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/ssa_combined_wgs_aquagenome_finnmark_canada_for_dxy.vcf.gz

tabix -p vcf /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/ssa_combined_wgs_aquagenome_finnmark_canada_for_dxy.vcf.gz

```

Repeat same thing for Canada and East Atlantic:
```
./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing --allow-extra-chr --keep-fam SouthNorCan.txt --make-bed --recode vcf --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/ssa_combined_wgs_aquagenome_southernnorway_canada_for_dxy

bgzip -c /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/ssa_combined_wgs_aquagenome_southernnorway_canada_for_dxy.vcf > /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/ssa_combined_wgs_aquagenome_southernnorway_canada_for_dxy.vcf.gz

tabix -p vcf /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Dxy/ssa_combined_wgs_aquagenome_southernnorway_canada_for_dxy.vcf.gz

```

Use these files for running R script for Dxy (PopGenome)

```
