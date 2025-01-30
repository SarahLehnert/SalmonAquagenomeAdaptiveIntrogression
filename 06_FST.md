Workflow for running FST analyses 

This was all done in plink. An R script was used to create the input files for renaming groups. After plink, files were subset and renamed by specific groups for FST. 

```
#rename individuals - group all canadian samples
  ./plink --bfile /filepath/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing --make-bed --out /filepath/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_jan2025_updated_ids_for_FST --update-ids /filepath/Salmon/WGS_Aquagenome/FST/Groups_for_FST_rename.txt --allow-extra-chr


#extract Barents/White Sea and Canada only, run FST between them
./plink --allow-extra-chr --bfile /filepath/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_jan2025_updated_ids_for_FST --keep /filepath/Salmon/WGS_Aquagenome/FST/Updated_famfile_IDs_Can_Finnmark.txt --make-bed --family --fst --out /filepath/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_jan2025_updated_ids_for_FST_Finnmark_Canada

```
Same thing was done for Canada and East Atlantic: 

```
#extract Southern Norway and Canada only, and run FST betwen them
./plink --allow-extra-chr --bfile /filepath/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_jan2025_updated_ids_for_FST --keep /filepath/Salmon/WGS_Aquagenome/FST/Updated_famfile_IDs_Can_SNor.txt --family --fst --make-bed --out /filepath/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST_SouthNor_Canada
```

Results were used in R for further analyses - see scripts
