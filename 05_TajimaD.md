Workflow for running vcftools for Tajima's D
VCFtools - 0.1.15

Analyses were run for Finnmark (Western Barents) populations. Use files created for SweeD analyses (by chromosome seperately)

```
for i in {1..9}; do vcftools --vcf  /filepath/Salmon/WGS_Aquagenome/SweeD/ssa0${i}_wg.phased_finnmark.txt.recode_maf001.vcf.recode.vcf --out /filepath/Salmon/WGS_Aquagenome/TajimaD/ssa0${i}_wg.phased_finnmark_TajimaD --phased --TajimaD 100000; done
for i in {10..29}; do vcftools --vcf  /filepath/Salmon/WGS_Aquagenome/SweeD/ssa${i}_wg.phased_finnmark.txt.recode_maf001.vcf.recode.vcf --out /filepath/Salmon/WGS_Aquagenome/TajimaD/ssa${i}_wg.phased_finnmark_TajimaD --phased --TajimaD 100000; done
```

Results were used in R for further analyses - see scripts
