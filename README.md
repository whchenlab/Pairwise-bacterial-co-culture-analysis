# Pairwise bacterial co-culture analysis
## Instructions
The related R packages were installed using the functions `install.packages()` and `BiocManager::install()`. It will only take a few minutes on Windows, but compiling the packages on Linux may take longer.

`Figure1.R`: The R script utilized to generate Figure 1 and related supplementary files in the article. The datasets utilized in the code are stored within the directory `Figure1/input`. 

`Figure2.R`: The R script utilized to generate Figure 2 and related supplementary files in the article. The datasets utilized in the code are stored within the directory `Figure1/input`.

Download all data and codes to your local machine, and run the code line by line in the R GUI to reproduce all the results presented in the paper.

## System requirements
```
R version 4.2.0 (2022-04-22 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ComplexHeatmap_2.15.1 ggbreak_0.1.2         gg.gap_1.3            patchwork_1.1.2       ggpmisc_0.5.3        
 [6] ggpp_0.5.2            drc_3.0-1             MASS_7.3-58.1         RColorBrewer_1.1-3    ggpubr_0.6.0         
[11] ggsci_3.0.0           writexl_1.4.2         readxl_1.4.2          lubridate_1.9.2       forcats_1.0.0        
[16] stringr_1.5.0         dplyr_1.1.2           purrr_1.0.1           readr_2.1.4           tidyr_1.3.0          
[21] tibble_3.2.1          ggplot2_3.4.2         tidyverse_2.0.0      

loaded via a namespace (and not attached):
 [1] TH.data_1.1-2          colorspace_2.0-3       rjson_0.2.21           ggsignif_0.6.4         circlize_0.4.15       
 [6] estimability_1.4.1     GlobalOptions_0.1.2    parameters_0.21.1      aplot_0.1.10           clue_0.3-64           
[11] rstudioapi_0.14        MatrixModels_0.5-1     fansi_1.0.4            mvtnorm_1.2-2          codetools_0.2-19      
[16] splines_4.2.0          mnormt_2.1.1           doParallel_1.0.17      zeallot_0.1.0          polynom_1.4-1         
[21] broom_1.0.5            cluster_2.1.4          png_0.1-8              pheatmap_1.0.12        compiler_4.2.0        
[26] emmeans_1.8.7          backports_1.4.1        Matrix_1.5-3           cli_3.6.1              quantreg_5.95         
[31] tools_4.2.0            ggstatsplot_0.11.1     coda_0.19-4            gtable_0.3.3           glue_1.6.2            
[36] Rcpp_1.0.10            carData_3.0-5          cellranger_1.1.0       vctrs_0.6.3            nlme_3.1-160          
[41] iterators_1.0.14       psych_2.3.6            insight_0.19.2         timechange_0.2.0       lifecycle_1.0.3       
[46] gtools_3.9.4           rstatix_0.7.2          zoo_1.8-12             scales_1.2.1           hms_1.1.3             
[51] parallel_4.2.0         sandwich_3.0-2         SparseM_1.81           rematch2_2.1.2         ggfun_0.1.1           
[56] yulab.utils_0.0.6      stringi_1.7.12         paletteer_1.5.0        bayestestR_0.13.1      S4Vectors_0.36.1      
[61] foreach_1.5.2          plotrix_3.8-2          permute_0.9-7          BiocGenerics_0.44.0    shape_1.4.6           
[66] rlang_1.1.1            pkgconfig_2.0.3        matrixStats_0.63.0     lattice_0.20-45        cowplot_1.1.1         
[71] tidyselect_1.2.0       plyr_1.8.8             magrittr_2.0.3         R6_2.5.1               IRanges_2.32.0        
[76] generics_0.1.3         multcomp_1.4-25        pillar_1.9.0           withr_2.5.0            mgcv_1.8-40           
[81] survival_3.4-0         datawizard_0.8.0       abind_1.4-5            crayon_1.5.2           car_3.1-2             
[86] utf8_1.2.3             correlation_0.8.4      tzdb_0.3.0             GetoptLong_1.0.5       vegan_2.6-4           
[91] digest_0.6.31          xtable_1.8-4           gridGraphics_0.5-1     statsExpressions_1.5.1 stats4_4.2.0          
[96] munsell_0.5.0          ggplotify_0.1.0       
```
