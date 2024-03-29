Supplementary information / reproducible research files for the manuscript 
Title: "Sensitivity Analysis of Causal Effects in Observational Studies with
Multivalued Treatments"

Authors: Md Abdul Basit, Mahbub A.H.M. Latif, Abdus S. Wahed
Author of the code: Md Abdul Basit (abasit@isrt.ac.bd)

The code was written/evaluated in R with the following software versions:
R version 4.3.2 (2023-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Asia/Dhaka
tzcode source: system (glibc)

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] kableExtra_1.3.4 lubridate_1.9.3  forcats_1.0.0    stringr_1.5.0    dplyr_1.1.3      purrr_1.0.2     
 [7] readr_2.1.4      tidyr_1.3.0      tibble_3.2.1     ggplot2_3.4.4    tidyverse_2.0.0 

loaded via a namespace (and not attached):
 [1] gtable_0.3.4       xfun_0.41          bench_1.1.3        tzdb_0.4.0         vctrs_0.6.3       
 [6] tools_4.3.2        generics_0.1.3     stats4_4.3.2       fansi_1.0.5        highr_0.10        
[11] R.oo_1.25.0        pkgconfig_2.0.3    data.table_1.14.8  RColorBrewer_1.1-3 webshot_0.5.4     
[16] R.cache_0.16.0     lifecycle_1.0.4    compiler_4.3.2     farver_2.1.1       munsell_0.5.0     
[21] tinytex_0.48       codetools_0.2-19   htmltools_0.5.7    yaml_2.3.7         pillar_1.9.0      
[26] crayon_1.5.2       R.utils_2.12.2     styler_1.10.2      tidyselect_1.2.0   rvest_1.0.3       
[31] digest_0.6.33      stringi_1.7.12     labeling_0.4.3     splines_4.3.2      fastmap_1.1.1     
[36] grid_4.3.2         colorspace_2.1-0   cli_3.6.1          magrittr_2.0.3     fastDummies_1.7.3 
[41] utf8_1.2.4         withr_2.5.2        scales_1.2.1       bit64_4.0.5        timechange_0.2.0  
[46] rmarkdown_2.22     httr_1.4.7         bit_4.0.5          nnet_7.3-19        R.methodsS3_1.8.2 
[51] hms_1.1.3          VGAM_1.1-8         evaluate_0.23      knitr_1.45         viridisLite_0.4.2 
[56] rlang_1.1.2        glue_1.6.2         xml2_1.3.5         svglite_2.1.2      rstudioapi_0.15.0 
[61] vroom_1.6.4        R6_2.5.1           systemfonts_1.0.5 



This folder named "supplementary_files" contains the following sub-folders and files that can be used to reproduce all analysis and figures of the manuscript. Please set the working directory of R to this folder before executing the files given in the folder.


./code/:
    sensitivity_analysis_functions.R: an R script that contains the functions written for conducting
		sensitivity analysis under the proposed framework.
		
    simulation_functions.R: an R script that contains the functions developed to conduct the
		simulation study presented in Section 4 of the submitted manuscript.
    
./data/:
    A subfolder containing the rda file nhanes.fish.rda, which has been extracted from the archived version
	  of the R package CrossScreening available at https://cran.r-project.org/src/contrib/Archive/CrossScreening/
    
./intermediate_results/
    A folder containing the results of simulation_study.Rmd or simulation_study.R 
    
    results_sim.R
    An R script that takes the RData files from the intermediate_results/ subfolder and creates the
    results tables and figures presented in the online supplementary information of the manuscript.

./simulation_study.Rmd:
	 An Rmd file that performs the simulations reported in the manuscript (Section 4). Simulation was performed
   in parallel on 20 cores on an Intel(R) Xeon(R) CPU E5-2650 linux server, but should also work under windows
	 by setting parallel = FALSE in the simulation_final_results() function.

./application.Rmd:
	An Rmd file that reproduces the tables and figures to demonstrate the proposed framework in the manuscript
	(Section 5).

./simulation_study.pdf:
	The rendered pdf output of ./simulation_study.Rmd.

./application.Rmd:
	The rendered pdf output of ./application.Rmd.
	
./simulation_study.R:
	An R script containing all the code chunks and comments of the ./simulation_study.Rmd file.

./application.R:
	An R script containing all the code chunks and comments of the ./application.Rmd file.