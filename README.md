# Overview
The following files allow one to reproduce analyses in the manuscript entitled "Can plant-soil feedbacks scale up from individual species to communities: tests of the richness and biomass-ratio hypotheses".


# Repo Contents
## Code: R code
This script contains all the data analysis and visualization code of Fig.1, Fig.3 & Table S2 S3 S4, Fig.S1,  Fig.S2,  Fig.S3.
## Data  
This script contains all the data of Fig.1, Fig.3 & Table S2 S3 S4, Fig.S1,  Fig.S2,  Fig.S3.
## File overview
1. Fig.1: Fig1_code.R & data_Fig.1.xlsx
2. Fig.3 & Table S2 S3 S4: Fig.3.R & data_Fig3.xlsx
3. Fig.S1: Fig.S1_code.R & data_Fig.S1.xlsx
4. Fig.S2: Fig.S2.R & data_FigS2.xlsx
5. Fig.S3: Fig.S3.R & data_FigS3.xlsx & Plant_tree.treefile 
### Data-specific information for: data_Fig.1.xlsx
1. The "bar" sheet:
   * MaxNo._Cspp.: the maximum number of conditioning plant species per experimental unit across these studies
   * 1993-2025: the number of PSF studies 
2. The "pie" sheet: 
   * MaxNo._Cspp.: the maximum number of conditioning plant species per experimental unit across these studies
   * MaxNo._Cspp.ID: ID of MaxNo._Cspp.
   * No._study_Cspp._peryear: the number of PSF studies per year for different MaxNo._Cspp.
   * Percentage : the proportion of No._study_Cspp._peryear to total studies (466 PSF experiments)
### Data-specific information for: data_Fig3.xlsx
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
Description: Phylogenetic relationship of all experimental species in Plant-soil feedback experiment.

 
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
###################################################################################
setwd('D:/2025.10.4-NPH/0-Code-alive-1216/code-0108-R1/Fig.3 & Table S2 S3 S4') 

#--------------------------------------
### Table S2
#--------------------------------------
data = read.xlsx("data_Fig3.xlsx",sheet="homeaway_raw"); data[1:6,1:12]
data <- data %>% slice(-(1:5)); data[1:6,1:12]
data$Richness_con = as.factor(data$Richness_con)

# Fit lm model with TG_res
mod_full1 <- lm(TG_res ~ Richness_con * (B_con + C_con + F_con + T_con + V_con), data = data)
anova(mod_full1)-> mod_full1_result; mod_full1_result
p1 <- mod_full1_result$Pr;p1 
p.adjust(p1, "BH")
#p.adjust(p1, "bonferroni")
~~~
## Expected output
~~~
> setwd('D:/2025.10.4-NPH/0-Code-alive-1216/code-0108-R1/Fig.3 & Table S2 S3 S4') 
> #--------------------------------------
> ### Table S2
> #--------------------------------------
> data = read.xlsx("data_Fig3.xlsx",sheet="homeaway_raw"); data[1:6,1:12]
  pot Pot_res      Pot_con Richness_con AG_con TG_res       PSFs c_con B_con C_con F_con T_con
1   1    CK_1 na (sterile)           NA     NA  87.54         NA    NO    NO    NO    NO    NO
2   2    CK_2 na (sterile)           NA     NA  92.14         NA    NO    NO    NO    NO    NO
3   3    CK_3 na (sterile)           NA     NA 100.32         NA    NO    NO    NO    NO    NO
4   4    CK_4 na (sterile)           NA     NA  80.18         NA    NO    NO    NO    NO    NO
5   5    CK_5 na (sterile)           NA     NA  74.82         NA    NO    NO    NO    NO    NO
6   6   1_1_1        1_1_1            1   8.63  68.32 -0.2417056    NO    NO    NO    NO    NO
> data <- data %>% slice(-(1:5)); data[1:6,1:12]
  pot Pot_res Pot_con Richness_con AG_con TG_res       PSFs c_con B_con C_con F_con T_con
1   6   1_1_1   1_1_1            1   8.63  68.32 -0.2417056    NO    NO    NO    NO    NO
2   7   1_1_2   1_1_2            1  12.22  65.21 -0.2882953    NO    NO    NO    NO    NO
3   8   1_1_3   1_1_3            1   7.02  76.33 -0.1308421    NO    NO    NO    NO    NO
4   9   1_1_4   1_1_4            1  11.03  61.24 -0.3511075    NO    NO    NO    NO    NO
5  10   1_1_5   1_1_5            1  14.74  59.51 -0.3797638    NO    NO    NO    NO    NO
6  11   1_2_1   1_2_1            1  11.02  68.15 -0.2441970    NO   YES    NO    NO    NO
> data$Richness_con = as.factor(data$Richness_con)
> # Fit lm model with TG_res
> mod_full1 <- lm(TG_res ~ Richness_con * (B_con + C_con + F_con + T_con + V_con), data = data)
> anova(mod_full1)-> mod_full1_result; mod_full1_result
Analysis of Variance Table

Response: TG_res
                    Df  Sum Sq Mean Sq F value  Pr(>F)  
Richness_con         4   607.6  151.89  1.2643 0.28454  
B_con                1    18.2   18.18  0.1514 0.69757  
C_con                1   449.8  449.76  3.7436 0.05412 .
F_con                1   103.9  103.94  0.8652 0.35317  
T_con                1   123.9  123.92  1.0315 0.31078  
V_con                1     9.2    9.19  0.0765 0.78236  
Richness_con:B_con   4   150.3   37.59  0.3129 0.86926  
Richness_con:C_con   4   506.4  126.60  1.0537 0.38003  
Richness_con:F_con   4   158.8   39.71  0.3305 0.85732  
Richness_con:T_con   4   457.6  114.40  0.9522 0.43438  
Richness_con:V_con   4   257.2   64.31  0.5353 0.70993  
Residuals          254 30515.8  120.14                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> p1 <- mod_full1_result$Pr;p1 
 [1] 0.28454201 0.69756845 0.05412107 0.35317423 0.31078445 0.78236394 0.86925594 0.38002853 0.85732088
[10] 0.43438274 0.70992651         NA
> p.adjust(p1, "BH")
 [1] 0.7963684 0.8692559 0.5953317 0.7963684 0.7963684 0.8692559 0.8692559 0.7963684 0.8692559 0.7963684
[11] 0.8692559        NA
#p.adjust(p1, "bonferroni")
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
