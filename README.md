# Overview
The following files allow one to reproduce analyses in the manuscript entitled "Can plant-soil feedbacks scale up from individual species to communities: tests of the richness and biomass-ratio hypotheses".


# Repo Contents
## Code: R code
This script contains all the data analysis and visualization code of Fig.1, Fig.3 & Table S2 S3 S4, Fig.S1,  Fig.S2,  Fig.S3.
## Data  
This script contains all the data of Fig.1, Fig.3 & Table S2 S3 S4, Fig.S1,  Fig.S2,  Fig.S3.
## File overview
1. Fig.2 & Table S1-S3: Fig.2 & TableS1-S3.R, Fig2 & S3 & TableS1-6.xlsx
2. Fig.S3 & TableS4-S6: Fig.S3 & TableS4-S6.R, Fig2 & S3 & TableS1-6.xlsx, Plant_tree.treefile
3. Fig.S1: FigS1_code.R, FigS1.xlsx
4. Fig.S2: FigS2.R, FigS2.xlsx

### Data-specific information for: Fig2 & S3 & TableS1-6.xlsx
1. The "Note" sheet:
   * SpeciesID: conditioning species ID
   * Species_abbrev: abbreviation for conditioning species 
2. The "homeaway_raw" sheet:
   * Tube number: sample ID
   * Pot_res: responding community ID (/pot)
   * Pot_con: conditioning community ID (/pot)
   * Richness_con: species richness of conditioning community
   * AG_con: aboveground biomass of conditioning community
   * TG_res: total biomass of responding community
   * PSFs: the natural logarithm ratio of biomass of the responding communities in conditioned soil and that grown in sterilized soil (average of five)
   * c_con: individual species (c, *Paspalum conjugatum*) presence in both the conditioning and responding phases (binary variable)
   * B_con: individual species (B, *Achyranthes aspera*) presence in both the conditioning and responding phases (binary variable)
   * C_con: individual species (C, *Setaria viridis*) presence in both the conditioning and responding phases (binary variable)
   * F_con: individual species (F, *Senna tora*) presence in both the conditioning and responding phases (binary variable)
   * T_con: individual species (T, *Amaranthus hybridus*) presence in both the conditioning and responding phases (binary variable)
   * V_con: individual species (V, *Senna occidentalis*) presence in both the conditioning and responding phases (binary variable) 
3. The "PSF_raw" sheet:
   * Tube number: sample ID
   * A-b: aboveground biomass of conditioning species plants
   * remark: sample classification: CK, conditioning community ID, conditioning species richness 
   * Lv: sample classification: conditioning species richness 
   * TG: total biomass of responding community
### Data-specific information for: data_FigS1.xlsx
1. The "FigS1" sheet: abundance of conditioning species
   * Tube number: sample ID
### Data-specific information for: data_FigS2.xlsx
1. The "Note" sheet:
   * SpeciesID: conditioning species ID
   * Species_abbrev: abbreviation for conditioning species plants
2. The "AG_con" sheet:
   * Pot_con: conditioning community ID (/pot)
   * SpeciesID: conditioning species ID
   * Species_con: abbreviation for conditioning species plants
   * AG_sp_con: aboveground biomass of conditioning species plants
   * AG_ind_con: mean aboveground biomass of individual conditioning species plants
   * Group: monocultures versus. mixtures
### Data-specific information for: data_FigS3.xlsx
1. The "Respond_AG" sheet:
   * Pot_res: responding community ID (/pot)
   * Paspalum_conjugatum: aboveground biomass of the response species *Paspalum conjugatum* 
   * Achyranthes_aspera: aboveground biomass of the response species *Achyranthes aspera*
   * Setaria_viridis: aboveground biomass of the response species *Setaria viridis*
   * Senna_tora: aboveground biomass of the response species *Senna tora*
   * Amaranthus_hybridus: aboveground biomass of the response species *Amaranthus hybridus*
   * Senna_occidentalis: aboveground biomass of the response species *Senna occidentalis*
2. The "Condition_AGraw" sheet: 
   * Pot_con: conditioning community ID (/pot)
3. The "species_specific_PSF" sheet: aboveground biomass of conditioning species: data in wide format
   * Pot_res: responding community ID (/pot)
   * Pot_con: conditioning community ID (/pot)
   * Richness_con: species richness of conditioning community
   * Pc_AG, Aa_AG, Sv_AG, St_AG, Ah_AG, Soc_AG: aboveground biomass of the responding species *Paspalum conjugatum*,*Achyranthes aspera*,*Setaria viridis*, *Senna tora*, *Amaranthus hybridus*, *Senna occidentalis*
   * Pc_PSF, Aa_PSF, Sv_PSF, St_PSF, Ah_PSF, Soc_PSF: species-specific effect of the responding species *Paspalum conjugatum*,*Achyranthes aspera*,*Setaria viridis*, *Senna tora*, *Amaranthus hybridus*, *Senna occidentalis*
   * Pc_Phylo_Dist, Aa_Phylo_Dist, Sv_Phylo_Dist, St_Phylo_Dist, Ah_Phylo_Dist, Soc_Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for *Paspalum conjugatum*,*Achyranthes aspera*,*Setaria viridis*, *Senna tora*, *Amaranthus hybridus*, *Senna occidentalis*
4. The "all_species" sheet:
   * Tube number: sample ID
   * Species_res: abbreviation for responding species plants 
   * Species_PSF: species-specific effect of each responding species (all responding species were considered together)
   * Species_Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for each response species (all responding species were considered together)
### Data-specific information for: Plant_tree.treefile
Description: Phylogenetic relationship of all experimental species in the plant-soil feedback experiment.

 
# System Requirements
## Hardware Requirements
R package requires only a standard computer with enough RAM to support the in-memory operations.
## Software Requirements
### OS Requirements
The package development version is tested on Windows and Mac systems. 

R version 4.2.2 or higher. 


# Installation Guide
## Intallation of R and Rstudio
R: https://www.r-project.org/

Rstudio: https://www.rstudio.com/

## Note: 

Install from Windows: follow the default path and click "Next".

Rstudio and R associated Settings (Rstudio->Tools->Global options)

The package should take approximately 40-60 seconds to install with vignettes on a recommended computer.

## Install software supporting R: JAVA, Rtools (windows system), Xcode(Mac system)

a) If you do not have Java installed on your computer, 
please download: https://www.java.com/zh-CN/download/ Download and install.

b) If you do not have Rtools installed on your computer, 
please download RTools from https://cran.r-project.org/bin/windows/Rtools/ to download for installation (please choose according to their own version of a R corresponding Rtools version), and check whether the installation is successful.

c) If your computer is a Mac system, 
from https://developer.apple.com/cn/xcode/, please download and install.

The package should take approximately 40-60 seconds to install with vignettes on a recommended computer.

## Package Installation
Users should install the following packages prior to library.
~~~
install.packages(c('readxl', 'openxlsx', 'ggplot2', 'ggpubr', 'ggridges', 'ggepi', 'ggrepel', 'ggspatial', 'ggfortify', 'ggforce', 'ggbeeswarm', 'ggmisc', 'ggsci', 'glmmTMB', 'dplyr', 'tidry', 'tidyverse', 'RColorBrewer', 'randomForest', 'party', 'picante', 'caret', 'export', 'car', 'ape', 'forcats', 'vegan', 'sf', patchwork', 'pheatmap', 'reshape2','cowplot', 'multcomp', 'DHARMa',  'emmeans', 'broom.mixed', 'MuMIn'))
~~~

Note: ggepi & patchwork can be installed via devtools 
~~~
 devtools::install_github("lwjohnst86/ggepi")
 devtools::install_github("thomasp85/patchwork")
~~~

The package should take approximately 50-60 seconds to install with vignettes on a recommended computer.


# Demo
## Instructions to run on data: example-“Fig.3 & Table S2 S3 S4” code
~~~
##################################################################################
#####                                                                        ##### 
#####           Part1---effect of home vs away:PSF                           #####
#####                                                                        #####  
##################################################################################

setwd('D:/2025.10.4-NPH/1-0330/Evan-0624/lu-0705/新建文件夹/Code-0715/Fig.2 & TableS1-S3') 

#-------------------------------------------------------------------------------------------- 
### Table S1 ---   responding community aboveground biomass & responding community PSFs
#-------------------------------------------------------------------------------------------- 
data = read.xlsx("Fig2 & S3 & TableS1-6.xlsx",sheet="pot_data"); data[1:6,1:12]
data <- data %>% slice(-(1:5)); data[1:6,1:12]
data$Richness_con = as.factor(data$Richness_con)

# Fit lm model with TG_res
mod_full1 <- lm(TG_res ~ Richness_con * (Aa_con + Sv_con + St_con + Ah_con + Soc_con), data = data)
anova(mod_full1)-> mod_full1_result; mod_full1_result
p1 <- mod_full1_result$Pr;p1 
p.adjust(p1, "BH")
#p.adjust(p1, "bonferroni")

# Fit lm model with PSF
mod_full2 <- lm(ActualPSF_Pot ~ Richness_con * (Aa_con + Sv_con + St_con + Ah_con + Soc_con), data = data)
anova(mod_full2)-> mod_full2_result; mod_full2_result
p2 <- mod_full2_result$Pr;p2 
p.adjust(p2, "BH")
#p.adjust(p2, "bonferroni")
~~~
## Expected output
~~~
##################################################################################
#####                                                                        ##### 
#####           Part1---effect of home vs away:PSF                           #####
#####                                                                        #####  
##################################################################################
> setwd('D:/2025.10.4-NPH/1-0330/Evan-0624/lu-0705/新建文件夹/Code-0715/Fig.2 & TableS1-S3') 
> #-------------------------------------------------------------------------------------------- 
> ### Table S1 ---   responding community aboveground biomass & responding community PSFs
> #-------------------------------------------------------------------------------------------- 
> data = read.xlsx("Fig2 & S3 & TableS1-6.xlsx",sheet="pot_data"); data[1:6,1:12]
    Pot Alive_all.conditioning.plants Composition_con Richness_con AG_con AG_res BG_res TG_res
1  <NA>                          <NA>            <NA>           NA     NA  74.83  12.71  87.54
2  <NA>                          <NA>            <NA>           NA     NA  77.92  14.22  92.14
3  <NA>                          <NA>            <NA>           NA     NA  85.38  14.94 100.32
4  <NA>                          <NA>            <NA>           NA     NA  69.58  10.60  80.18
5  <NA>                          <NA>            <NA>           NA     NA  63.83  10.69  74.52
6 1_1_1                           YES             Sor            1   8.63  58.04  10.28  68.32
  Pc_average_res Aa_average_res Sv_average_res St_average_res
1          5.815          0.255         17.050          1.150
2          2.265          0.080         24.320          1.175
3          4.950          0.930         26.530          1.410
4          3.230          0.175         14.715          1.130
5          3.915          0.315         10.245          2.070
6          0.410          0.000         26.175          0.495
> data <- data %>% slice(-(1:5)); data[1:6,1:12]
    Pot Alive_all.conditioning.plants Composition_con Richness_con AG_con AG_res BG_res TG_res
1 1_1_1                           YES             Sor            1  8.630  58.04  10.28  68.32
2 1_1_2                           YES             Sor            1 12.215  55.15  10.06  65.21
3 1_1_3                           YES             Sor            1  7.024  63.96  12.37  76.33
4 1_1_4                           YES             Sor            1 11.021  49.75  11.49  61.24
5 1_1_5                           YES             Sor            1 14.738  53.41   6.10  59.51
6 1_2_1                           YES              Aa            1 11.013  57.97  10.18  68.15
  Pc_average_res Aa_average_res Sv_average_res St_average_res
1          0.410          0.000         26.175          0.495
2          1.145          0.055         18.260          4.435
3          1.330          0.050         29.000          0.680
4          2.170          0.025         15.625          1.955
5          0.280          0.075         18.600          1.620
6          0.265          0.125         18.020          7.585
> data$Richness_con = as.factor(data$Richness_con)
> # Fit lm model with TG_res
> mod_full1 <- lm(TG_res ~ Richness_con * (Aa_con + Sv_con + St_con + Ah_con + Soc_con), data = data)
> anova(mod_full1)-> mod_full1_result; mod_full1_result
Analysis of Variance Table

Response: TG_res
                      Df  Sum Sq Mean Sq F value  Pr(>F)  
Richness_con           4   664.9  166.22  1.3891 0.23813  
Aa_con                 1     8.2    8.19  0.0685 0.79378  
Sv_con                 1   523.0  523.04  4.3711 0.03755 *
St_con                 1   197.7  197.73  1.6524 0.19980  
Ah_con                 1   125.2  125.24  1.0466 0.30726  
Soc_con                1     4.9    4.94  0.0413 0.83910  
Richness_con:Aa_con    4   236.2   59.05  0.4935 0.74055  
Richness_con:Sv_con    4   491.6  122.91  1.0272 0.39371  
Richness_con:St_con    4   247.9   61.98  0.5179 0.72261  
Richness_con:Ah_con    4   454.5  113.63  0.9496 0.43586  
Richness_con:Soc_con   4   234.3   58.59  0.4896 0.74337  
Residuals            254 30393.0  119.66                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> p1 <- mod_full1_result$Pr;p1 
 [1] 0.23813017 0.79377615 0.03754687 0.19979947 0.30725677 0.83909649 0.74054825 0.39371391 0.72261087
[10] 0.43586449 0.74337200         NA
> p.adjust(p1, "BH")
 [1] 0.7990849 0.8390965 0.4130156 0.7990849 0.7990849 0.8390965 0.8390965 0.7990849 0.8390965 0.7990849
[11] 0.8390965        NA
> # Fit lm model with PSF
> mod_full2 <- lm(ActualPSF_Pot ~ Richness_con * (Aa_con + Sv_con + St_con + Ah_con + Soc_con), data = data)
> anova(mod_full2)-> mod_full2_result; mod_full2_result
Analysis of Variance Table

Response: ActualPSF_Pot
                      Df Sum Sq  Mean Sq F value  Pr(>F)  
Richness_con           4 0.1320 0.033010  1.2119 0.30620  
Aa_con                 1 0.0008 0.000756  0.0278 0.86779  
Sv_con                 1 0.1275 0.127488  4.6804 0.03144 *
St_con                 1 0.0441 0.044118  1.6197 0.20430  
Ah_con                 1 0.0289 0.028923  1.0618 0.30377  
Soc_con                1 0.0009 0.000859  0.0315 0.85922  
Richness_con:Aa_con    4 0.0471 0.011770  0.4321 0.78537  
Richness_con:Sv_con    4 0.1605 0.040133  1.4734 0.21071  
Richness_con:St_con    4 0.0695 0.017382  0.6381 0.63574  
Richness_con:Ah_con    4 0.1364 0.034097  1.2518 0.28958  
Richness_con:Soc_con   4 0.0439 0.010963  0.4025 0.80678  
Residuals            254 6.9186 0.027239                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> p2 <- mod_full2_result$Pr;p2 
 [1] 0.30620227 0.86778938 0.03144185 0.20429785 0.30377487 0.85922423 0.78537251 0.21070849 0.63574279
[10] 0.28958309 0.80678060         NA
> p.adjust(p2, "BH")
 [1] 0.5613708 0.8677894 0.3458604 0.5613708 0.5613708 0.8677894 0.8677894 0.5613708 0.8677894 0.5613708
[11] 0.8677894        NA
~~~
## Note: 
The dataset stored in this repository is same to the dataset in figshare( ).
Running each line of code on a recommended computer takes about 1-5 seconds. A notable exception is Step 4 in the "#### explanations - Modified to only include the three R² values requested" section of Fig. 3’s code, which takes 5-10 minutes.
### Instructions to run on data: example-“ Step 4” code
~~~
#----------------------------------------------
# Step 4: Run bootstrap iterations (using only RandomForest)
#----------------------------------------------
# Identify indices for stratified sampling
id_lv1 <- which(df.rf$Lv == 1)
id_lvh <- which(df.rf$Lv > 1)
n_lv1 <- length(id_lv1)
n_lvh <- length(id_lvh)

# Progress tracker
cat("Starting bootstrap iterations...\n")
prog_interval <- ceiling(n_iter/10)

for (i in 1:n_iter) {
  # Print progress
  if (i %% prog_interval == 0) {
    cat(sprintf("Completed %d%% (%d/%d iterations)\n", 
                round(i/n_iter*100), i, n_iter))
  }
  
  #----------------------------------------------
  # Standard RandomForest with train/test split
  #----------------------------------------------
  tryCatch({
    # Create a stratified bootstrap sample
    bs_indices <- c(
      sample(id_lv1, ceiling(0.5*n_lv1), replace = TRUE),
      sample(id_lvh, ceiling(0.5*n_lvh), replace = TRUE)
    )
    bs_data <- df.rf[bs_indices, ]
    
    # Train/test split
    train_indices <- sample(nrow(bs_data), round(cv_prop * nrow(bs_data)))
    train_data <- bs_data[train_indices, ]
    test_data <- bs_data[-train_indices, ]
    
    # Fit the three models
    model_Lv <- randomForest(fml_Lv, data = train_data, 
                             ntree = n_tree, mtry = min(2, n_pred_Lv))
    
    model_LvID <- randomForest(fml_LvID, data = train_data, 
                               ntree = n_tree, mtry = min(3, n_pred_LvID))
    
    model_Full <- randomForest(fml_Full, data = train_data, 
                               ntree = n_tree, mtry = min(5, n_pred_Full))
    
    # Make predictions
    pred_Lv <- predict(model_Lv, newdata = test_data)
    pred_LvID <- predict(model_LvID, newdata = test_data)
    pred_Full <- predict(model_Full, newdata = test_data)
    
    # Calculate R² (squared correlation)
    r2_Lv <- cor(pred_Lv, test_data[[response_var]])^2
    r2_LvID <- cor(pred_LvID, test_data[[response_var]])^2
    r2_Full <- cor(pred_Full, test_data[[response_var]])^2
    
    # Store R² values
    idx <- (i-1)*3 + 1:3
    results_rf$R2_RF[idx] <- c(r2_Lv, r2_LvID, r2_Full)
    
  }, error = function(e) {
    cat("Error in RF iteration", i, ":", conditionMessage(e), "\n")
  })
}
~~~
### Expected output

~~~
#----------------------------------------------
> # Step 4: Run bootstrap iterations (using only RandomForest)
> #----------------------------------------------
> # Identify indices for stratified sampling
> id_lv1 <- which(df.rf$Lv == 1)
> id_lvh <- which(df.rf$Lv > 1)
> n_lv1 <- length(id_lv1)
> n_lvh <- length(id_lvh)
> # Progress tracker
> cat("Starting bootstrap iterations...\n")
Starting bootstrap iterations...
> prog_interval <- ceiling(n_iter/10)
> for (i in 1:n_iter) {
+   # Print progress
+   if (i %% prog_interval == 0) {
+     cat(sprintf("Completed %d%% (%d/%d iterations)\n", 
+                 round(i/n_iter*100), i, n_iter))
+   }
+   
+   #----------------------------------------------
+   # Standard RandomForest with train/test split
+   #----------------------------------------------
+   tryCatch({
+     # Create a stratified bootstrap sample
+     bs_indices <- c(
+       sample(id_lv1, ceiling(0.5*n_lv1), replace = TRUE),
+       sample(id_lvh, ceiling(0.5*n_lvh), replace = TRUE)
+     )
+     bs_data <- df.rf[bs_indices, ]
+     
+     # Train/test split
+     train_indices <- sample(nrow(bs_data), round(cv_prop * nrow(bs_data)))
+     train_data <- bs_data[train_indices, ]
+     test_data <- bs_data[-train_indices, ]
+     
+     # Fit the three models
+     model_Lv <- randomForest(fml_Lv, data = train_data, 
+                              ntree = n_tree, mtry = min(2, n_pred_Lv))
+     
+     model_LvID <- randomForest(fml_LvID, data = train_data, 
+                                ntree = n_tree, mtry = min(3, n_pred_LvID))
+     
+     model_Full <- randomForest(fml_Full, data = train_data, 
+                                ntree = n_tree, mtry = min(5, n_pred_Full))
+     
+     # Make predictions
+     pred_Lv <- predict(model_Lv, newdata = test_data)
+     pred_LvID <- predict(model_LvID, newdata = test_data)
+     pred_Full <- predict(model_Full, newdata = test_data)
+     
+     # Calculate R² (squared correlation)
+     r2_Lv <- cor(pred_Lv, test_data[[response_var]])^2
+     r2_LvID <- cor(pred_LvID, test_data[[response_var]])^2
+     r2_Full <- cor(pred_Full, test_data[[response_var]])^2
+     
+     # Store R² values
+     idx <- (i-1)*3 + 1:3
+     results_rf$R2_RF[idx] <- c(r2_Lv, r2_LvID, r2_Full)
+     
+   }, error = function(e) {
+     cat("Error in RF iteration", i, ":", conditionMessage(e), "\n")
+   })
+ }
Completed 10% (100/1000 iterations)
Completed 20% (200/1000 iterations)
Completed 30% (300/1000 iterations)
Completed 40% (400/1000 iterations)
Completed 50% (500/1000 iterations)
Completed 60% (600/1000 iterations)
Completed 70% (700/1000 iterations)
Completed 80% (800/1000 iterations)
Completed 90% (900/1000 iterations)
Completed 100% (1000/1000 iterations)
~~~
# License
This project is covered under the GNU GENERAL PUBLIC LICENSE Version 2 (GPL-2).
