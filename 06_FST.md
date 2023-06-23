Workflow for running FST analyses 

This was all done in plink. Files were subset and renamed by specific groups for FST.

```
#rename individuals - group all canadian samples
./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/Treemix/ssa_combined_wgs_aquagenome_biallelic_PASS_maf001_phased_imputed_plink_update_ids_by_REGION_maf_after_phasing --allow-extra-chr --update-ids ~/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/Groups_for_FST_rename.txt --make-bed --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST

#extract Finnmark and Canada only
./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST --allow-extra-chr --keep ~/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/Updated_famfile_IDs_Can_Finnmark.txt --make-bed --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST_Finnmark_Canada

#run FST between Finnmark and Canada
./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST_Finnmark_Canada --fst --family --allow-extra-chr --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST_Finnmark_Canada_results

```
Same thing was done for Canada and Southern Norway:

```
#extract Southern Norway and Canada only
./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST --allow-extra-chr --keep ~/Desktop/Sarah/Salmon/WGS_Aquagenome/FST//Updated_famfile_IDs_Can_SNor.txt --make-bed --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST_SouthNor_Canada

#run FST between Southern Norway and Canada
./plink --bfile /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST_SouthNor_Canada --fst --family --allow-extra-chr --out /Users/ianbradbury/Desktop/Sarah/Salmon/WGS_Aquagenome/FST/ssa_combined_wgs_nov2020_updated_ids_for_FST_SouthNor_Canada_results
```

Results were used in R for further analyses - see scripts
