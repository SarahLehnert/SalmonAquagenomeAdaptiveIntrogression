Workflow for running selection analysis in SweeD 

First - filter data used for PCAdmix by chr in vcftools

```
for i in {1..9}; do vcftools --vcf /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_phased.vcf.recode.vcf  --chr ssa0${i} --recode --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/SweeD/ssa0${i}_wg.phased_finnmark.txt.recode_maf001.vcf --phased --maf 0.01; done
for i in {10..29}; do vcftools --vcf /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Pcadmix/ssa_combined_wgs_for_pcadmix_FINNMARK_phased.vcf.recode.vcf  --chr ssa${i} --recode --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/SweeD/ssa${i}_wg.phased_finnmark.txt.recode_maf001.vcf --phased --maf 0.01; done
```

Use output files (by chromosome) to run SweeD on Finnmark (Western Barents) populations

```
./SweeD-P -input <vcf> -grid <number of grids> -name <output> -folded

for i in {1..9}; do ./SweeD -input /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/SweeD/ssa0${i}_wg.phased_finnmark.txt.recode_maf001.vcf.recode.vcf -grid 500 -name ssa0${i}_wg_sweed_finnmark; done

for i in {10..19}; do ./SweeD -input /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/SweeD/ssa${i}_wg.phased_finnmark.txt.recode_maf001.vcf.recode.vcf -grid 500 -name ssa${i}_wg_sweed_finnmark; done

for i in {20..29}; do ./SweeD -input /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/SweeD/ssa${i}_wg.phased_finnmark.txt.recode_maf001.vcf.recode.vcf -grid 500 -name ssa${i}_wg_sweed_finnmark; done
```

Results were imported into R for further analyses - see scripts
