# Rubin et al - Moon moth bat predator richness

Code used to generate SDM for insectivorous bat predators of moon moths globally (i.e., insectivorous bat species > 10g) to estimate relative richness.

# This repository contains:
## R/

Scripts used to read, tidy, and standardize data and run analyses.
| Script                                 | Description |
| -----------                            | ----------- |
| 00_Setup.R                             | Set up workspace, load widely used libraries |
| 01_loadSpatialData.R                   | Short script to generate country-level maps |
| 02_synthOccurrenceRecords.R            | Integrate GBIF and iDigBio occurrence records |
| 03_findSynonyms.R                      | Manually synonym adjustments and notes for all species with records (see script 2)  |
| 04_cleanCoords.R                       | Automatic cleaning and manual checks and corrections for all coorinates. |
| 06_0_prepToRunSDMPipeline.R            | Load scripts. predictors, species list for SDM pipeline. |
| 06_1_runSDMPipeline.R                  | Run SDM pipeline for all species. |
| 06_2_runSDMPipeline_manualFixes.R      | Manual reruns of SDM pipeline for species that errored out or required manual corrections. |
| 07_checkSDMoutputs.R                   | Application of manual checks of SDM outputs, highlighting pass/fails with expert elicitation. |
| 08_combineContinuousSurfaces.R         | Sum continuous SDM into one richness surface. |
| 09_extractBatRichnessValues.R          | Extract bat richness at sites where moon moth tails had been measured. |
| 10_extractBatSpeciesAndRegions.R       | Define species richness in focal moon moth regions, list bat species present in each region, plot richness map. |


## out/
| File                          | Description |
| -----------                   | ----------- |
| batRegions.csv                | Focal bat species found in moon moth regions. |
| combinedContinuousSurface.tif | Insectivorous bat richness (tif file). |
| includedBatSpecies.csv        | Bat species included in this analysis (successful SDM runs). |


## Session info
Coordinate cleaning and SDM checks were conducted with the following session information:
```
R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 13.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tidyr_1.2.1              ENMeval_2.0.3            dismo_1.3-9             
 [4] raster_3.6-3             rgeos_0.5-9              sp_1.5-1                
 [7] luna_0.3-4               ggpubr_0.5.0             RSQLite_2.2.18          
[10] dbplyr_2.2.1             stringr_1.4.1            CoordinateCleaner_2.0-20
[13] taxotools_0.0.129        rnaturalearth_0.1.0      sf_1.0-9                
[16] data.table_1.14.4        ggplot2_3.4.0            terra_1.7-29            
[19] dplyr_1.0.10             magrittr_2.0.3          

loaded via a namespace (and not attached):
 [1] nlme_3.1-160       bold_1.2.0         oai_0.4.0          bit64_4.0.5       
 [5] httr_1.4.4         rgbif_3.7.3        tools_4.1.2        backports_1.4.1   
 [9] utf8_1.2.2         rgdal_1.6-2        R6_2.5.1           KernSmooth_2.23-20
[13] DBI_1.1.3          lazyeval_0.2.2     colorspace_2.0-3   withr_2.5.0       
[17] tidyselect_1.2.0   bit_4.0.4          curl_4.3.3         compiler_4.1.2    
[21] chron_2.3-58       cli_3.4.1          xml2_1.3.3         scales_1.2.1      
[25] classInt_0.4-8     readr_2.1.3        proxy_0.4-27       digest_0.6.30     
[29] rmarkdown_2.18     stringdist_0.9.10  pkgconfig_2.0.3    htmltools_0.5.3   
[33] fastmap_1.1.0      rlang_1.0.6        rstudioapi_0.14    httpcode_0.3.0    
[37] generics_0.1.3     zoo_1.8-11         jsonlite_1.8.3     car_3.1-1         
[41] geosphere_1.5-14   Rcpp_1.0.10        munsell_0.5.0      fansi_1.0.3       
[45] abind_1.4-5        ape_5.6-2          proto_1.0.0        lifecycle_1.0.3   
[49] sqldf_0.4-11       stringi_1.7.8      whisker_0.4        carData_3.0-5     
[53] plyr_1.8.8         grid_4.1.2         blob_1.2.3         parallel_4.1.2    
[57] crayon_1.5.2       lattice_0.20-45    conditionz_0.1.0   hms_1.1.2         
[61] knitr_1.40         pillar_1.8.1       uuid_1.1-0         taxize_0.9.100    
[65] ggsignif_0.6.4     codetools_0.2-18   crul_1.3           glue_1.6.2        
[69] evaluate_0.18      vctrs_0.5.0        tzdb_0.3.0         foreach_1.5.2     
[73] gtable_0.3.1       purrr_0.3.5        reshape_0.8.9      assertthat_0.2.1  
[77] gsubfn_0.7         cachem_1.0.6       xfun_0.34          broom_1.0.1       
[81] e1071_1.7-12       rstatix_0.7.1      class_7.3-20       tibble_3.1.8      
[85] iterators_1.0.14   memoise_2.0.1      units_0.8-0        ellipsis_0.3.2    
[89] wikitaxa_0.4.0 
```

The SDM pipeline was run for most species under the following configuration:
```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] stringr_1.5.0       ENMeval_2.0.1       dismo_1.3-3         raster_3.5-21      
 [5] rgeos_0.5-2         sp_1.6-0            rnaturalearth_0.1.0 sf_1.0-10          
 [9] data.table_1.12.2   ggplot2_3.4.1       terra_1.6-51        dplyr_1.1.0        
[13] magrittr_1.5       

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.9         pillar_1.7.0       compiler_3.6.3     iterators_1.0.13  
 [5] class_7.3-17       tools_3.6.3        lattice_0.20-38    lifecycle_1.0.3   
 [9] tibble_3.1.8       gtable_0.3.0       pkgconfig_2.0.2    rlang_1.0.6       
[13] foreach_1.5.1      cli_3.6.0          DBI_1.1.1          rstudioapi_0.10   
[17] yaml_2.2.0         e1071_1.7-3        withr_2.5.0        generics_0.0.2    
[21] vctrs_0.5.2        classInt_0.4-2     grid_3.6.3         tidyselect_1.2.0  
[25] glue_1.6.2         R6_2.4.0           fansi_0.4.0        scales_1.2.1      
[29] codetools_0.2-16   ellipsis_0.3.2     units_0.8-0        colorspace_1.4-1  
[33] utf8_1.1.4         KernSmooth_2.23-17 stringi_1.7.12     munsell_0.5.0     
[37] crayon_1.3.4      
```
