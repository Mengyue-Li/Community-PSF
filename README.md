# Overview
The following files allow one to reproduce analyses in the manuscript entitled "Can plant-soil feedbacks scale up from individual species to communities: tests of the richness and biomass-ratio hypotheses".


# Repo Contents
## Code: R code
This script contains all the data analysis and visualization code of Fig2-3 & TableS1-6, Fig.S1 and Fig.S2.
## Data  
This script contains all the data of Fig2-3 & TableS1-6, Fig.S1 and Fig.S2.
## File overview
1. Fig2-3 & TableS1-6.: Fig2-3 & TableS1-6.R, data-Fig2-3 & TableS1-6.xlsx, Plant_tree.treefile
2. Fig.S1: Fig.S1.R, data-Fig.S1.xlsx
3. Fig.S2: Fig.S2.R, data-Fig.S2.xlsx

### Data-specific information for: data-Fig2-3 & TableS1-6.xlsx
1. The "Fig2_data" sheet:
   * Pot: rsponding species ID
   * Sor-Cs: aboveground biomass of conditioning species plants
   * remark: sample classification: CK, conditioning community ID, conditioning species richness 
   * Richness_con: sample classification: conditioning species richness 
   * TG_Pot_res: total biomass of responding community
   * Pc_averageAG_res: aboveground biomass of responding species Pc *Paspalum conjugatum*
   * Aa_averageAG_res: aboveground biomass of responding species Aa *Achyranthes aspera*
   * Sv_averageAG_res: aboveground biomass of responding species Sv *Setaria viridis*
   * St_averageAG_res: aboveground biomass of responding species St *Senna tora*
   * Ah_averageAG_res: aboveground biomass of responding species Ah *Amaranthus hybridus*
   * Soc_averageAG_res: aboveground biomass of responding species Soc *Senna occidentalis*
2. The "data_regress" sheet:
   * Pot: rsponding species ID
   * Alive_all_conditioning_plants: pots with complete survival of all 12 individuals were scored as 'YES'; pots with missing    biomass (i.e., any mortality) were scored as 'NO'.
   * Composition_con: composition of conditioning community
   * Richness_con: species richness of conditioning community
   * Species: the abbreviation of conditioning plant species
   * Species_full: the Latin name for conditioning plant species
   * ActualPSF: Observed community-level PSFs
   * RichnessPSF: Predicted community-level PSF values based on the richness hypothesis
   * BiomassRatioPSF: predicted community-level PSF values based on the biomass-ratio hypotheses
   * Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for each response species (each responding species were considered together)
   * Species_con: individual species (c, *Paspalum conjugatum*, *Achyranthes aspera*, *Setaria viridis*, *Senna tora*, *Senna occidentalis*) presence in both the conditioning and responding phases (binary variable)  
3. The "Species_AG_con" sheet:
   * SpeciesID: conditioning species ID
   * The Latin name for conditioning plant species 
4. The "pot_data" sheet:
   * Pot: rsponding pot ID
   * Alive_all_conditioning_plants: Pots with complete survival of all 12 individuals were scored as 'YES'; pots with missing    biomass (i.e., any mortality) were scored as 'NO'.
   * Composition_con: Composition of conditioning community
   * Richness_con: species richness of conditioning community
   * AG_con: Aboveground biomass of conditioning community
   * Pc_average_res: aboveground biomass of the response species *Paspalum conjugatum* 
   * Aa_average_res: aboveground biomass of the response species *Achyranthes aspera*
   * Sv_average_res: aboveground biomass of the response species *Setaria viridis*
   * St_average_res: aboveground biomass of the response species *Senna tora*
   * Ah_average_res: aboveground biomass of the response species *Amaranthus hybridus*
   * Soc_average_res: aboveground biomass of the response species *Senna occidentalis*
   * Pc_con：*Paspalum conjugatum* presence in both the conditioning and responding phases (binary variable)  
   * Aa_con：*Achyranthes aspera* presence in both the conditioning and responding phases (binary variable)  
   * Sv_con：*Setaria viridis* presence in both the conditioning and responding phases (binary variable)  
   * St_con：*Senna tora* presence in both the conditioning and responding phases (binary variable)  
   * Ah_con：*Amaranthus hybridus* presence in both the conditioning and responding phases (binary variable)  
   * Soc_con：*Senna occidentalis* presence in both the conditioning and responding phases (binary variable)  
   * ActualPSF: Observed community-level PSFs
   * RichnessPSF: Predicted community-level PSF values based on the richness hypothesis
   * BiomassRatioPSF: predicted community-level PSF values based on the biomass-ratio hypothese
   * ActualPSF_Pc：observed species-specific PSFs of *Paspalum conjugatum*
   * RichnessPSF_Pc：richness predicted species-specific PSFs of *Paspalum conjugatum*
   * BiomassRatioPSF_Pc：biomass-ratio predicted species-specific PSFs of *Paspalum conjugatum*
   * ActualPSF_Aa：observed species-specific PSFs of *Achyranthes aspera*
   * RichnessPSF_Aa：richness predicted species-specific PSFs of *Achyranthes aspera*
   * BiomassRatioPSF_Aa：biomass-ratio predicted species-specific PSFs of *Achyranthes aspera*
   * ActualPSF_Sv：observed species-specific PSFs of *Setaria viridis*
   * RichnessPSF_Sv：richness predicted species-specific PSFs of *Setaria viridis*
   * BiomassRatioPSF_Sv：biomass-ratio predicted species-specific PSFs of *Setaria viridis*
   * ActualPSF_St：observed species-specific PSFs of *Senna tora*
   * RichnessPSF_St：richness predicted species-specific PSFs of *Senna tora*
   * BiomassRatioPSF_St：biomass-ratio predicted species-specific PSFs of *Senna tora*
   * ActualPSF_Ah：observed species-specific PSFs of *Amaranthus hybridus*
   * RichnessPSF_Ah：richness predicted species-specific PSFs of *Amaranthus hybridus*
   * BiomassRatioPSF_Ah：biomass-ratio predicted species-specific PSFs of *Amaranthus hybridus*
   * ActualPSF_Soc：observed species-specific PSFs of *Senna occidentalis*
   * RichnessPSF_Soc：richness predicted species-specific PSFs of *Senna occidentalis*
   * BiomassRatioPSF_Soc：biomass-ratio predicted species-specific PSFs of*Senna occidentalis*
   * Pc_Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for *Paspalum conjugatum*
   * Aa_Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for *Achyranthes aspera*
   * Sv_Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for *Setaria viridis*
   * St_Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for *Senna tora*
   * Ah_Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for *Amaranthus hybridus*
   * Soc_Phylo_Dist: the weighted phylogenetic distance to the conditioning communities for*Senna occidentalis*
### Data-specific information for: data-Fig.S1.xlsx
1. The "FigS1" sheet: abundance of conditioning species
   * Pot_ID: pot IDs under different conditioning species richness
   * Sor-Cs: abundance of conditioning species plants
### Data-specific information for: data-Fig.S2.xlsx
1. The "AG_con" sheet:
   * Pot: conditioning community ID (/pot)
   * SpeciesID: conditioning species ID
   * Species: abbreviation for conditioning species plants
   * AG_sp_con: aboveground biomass of conditioning species plants
   * AG_ind_con: mean aboveground biomass of individual conditioning species plants
   * Group: monocultures versus. mixtures
1. The "mono-data" sheet:
   * Pot: conditioning community ID (/pot)
   * SpeciesID: conditioning species ID
   * Species: abbreviation for conditioning species plants
   * Group: monocultures
   * AG_sp_con: total aboveground biomass for plant monocultures
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
install.packages(c('readxl', 'openxlsx', 'ggplot2', 'ggpubr', 'ggridges', 'ggepi', 'ggrepel', 'ggspatial', 'ggfortify', 'ggforce', 'ggbeeswarm', 'ggmisc', 'ggsci', 'glmmTMB', 'dplyr', 'tidry', 'tidyverse', 'RColorBrewer', 'randomForest', 'party', 'picante', 'caret', 'export', 'car', 'ape', 'forcats', 'vegan', 'sf', patchwork', 'pheatmap', 'reshape2','cowplot', 'multcomp', 'DHARMa',  'emmeans', 'broom.mixed', 'MuMIn', 'multcompView', 'multcomp'))
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

setwd('D:/2025.10.4-NPH/1-0330/1-Revised manuscript 2026/Code-0715/Fig2-3 & TableS1-6')  

#-------------------------------------------------------------------------------------------- 
### Table S1 ---   responding community aboveground biomass & responding community PSFs
#-------------------------------------------------------------------------------------------- 
data = read.xlsx("data-Fig2-3 & TableS1-6.xlsx",sheet="pot_data"); data[1:6,1:12]
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
setwd('D:/2025.10.4-NPH/1-0330/1-Revised manuscript 2026/Code-0715/Fig2-3 & TableS1-6') 
> #-------------------------------------------------------------------------------------------- 
> ### Table S1 ---   responding community aboveground biomass & responding community PSFs
> #-------------------------------------------------------------------------------------------- 
>data = read.xlsx("data-Fig2-3 & TableS1-6.xlsx",sheet="pot_data"); data[1:6,1:12]
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
Running each line of code on a recommended computer takes about 1-5 seconds. A notable exception is Step 4 in the "#### explanations - Modified to only include the three R² values requested" section of Fig. 2’s code, which takes 5-10 minutes.
### Instructions to run on data: example-“Part 4---Random forest analysis (using only actual experimental pot data)” code
~~~
##################################################################################################################
#####                                                                                                        ##### 
#####               Part 4---Random forest analysis (using only actual experimental pot data)                #####
#####                                                                                                        #####  
##################################################################################################################
# ------------------------------------------------------------------------------
#  Step1: – Random forest analysis
# ------------------------------------------------------------------------------
set.seed(123)
 
# filter out the duplicates - each pot should only be used once in the analysis
df.rf <- df[df$remark %in% levels_rich, ] 
df.rf <- merge(df.rf, pot_predictions[, c("Pot", "pred_RichnessPSF_Pot", "pred_BiomassRatioPSF_Pot")], 
               by = "Pot", all.x = TRUE) # Merge predicted values
rownames(df.rf) <- df.rf$Row.names
df.rf$Row.names <- NULL

# Define random forest model formulas
fml_richness <- formula(TG_Pot_res ~ Richness_con)   # Only richness
fml_identity <- formula(paste("TG_Pot_res ~ Richness_con +", paste(stressors, collapse = " + ")))  # Richness + species identity
fml_Effect_size <- formula(paste("TG_Pot_res ~ Richness_con +", paste(stressors, collapse = " + "), 
                                 "+ pred_RichnessPSF_Pot + pred_BiomassRatioPSF_Pot"))  # Richness + identity + effect size
n_tree <- 1000  # Number of trees per forest
n_iter_rf <- 1000  # Bootstrap iterations for RF
cv_prop <- 0.7  # Training set proportion

# Number of predictors in each model
n_pred_richness <- length(all.vars(fml_richness)) - 1
n_pred_identity <- length(all.vars(fml_identity)) - 1
n_pred_Effect_size <- length(all.vars(fml_Effect_size)) - 1

# Data frame to store results
results_rf <- data.frame(Iteration = rep(1:n_iter_rf, each = 3),
                         Model = rep(c("Species_richness", "Species_identity", "Effect_size"), n_iter_rf),
                         R2_RF = NA, stringsAsFactors = FALSE)

# Split samples by richness level (for stratified bootstrap)
id_lv1 <- which(df.rf$Richness_con == 1)
id_lvh <- which(df.rf$Richness_con != 1)
n_lv1 <- length(id_lv1); n_lv1     #********************* 140
n_lvh <- length(id_lvh); n_lvh      #********************* 144

prog_interval <- ceiling(n_iter_rf/10)
for (i in 1:n_iter_rf) {
  if (i %% prog_interval == 0) cat(sprintf("Completed %d%% (%d/%d iterations)\n", round(i/n_iter_rf*100), i, n_iter_rf))
  tryCatch({
    # Stratified bootstrap sampling (preserve proportion of richness 1 and >1)
    bs_indices <- c(sample(id_lv1, ceiling(0.5*n_lv1), replace = TRUE),
                    sample(id_lvh, ceiling(0.5*n_lvh), replace = TRUE))
    bs_data <- df.rf[bs_indices, ]
    train_idx <- sample(nrow(bs_data), round(cv_prop * nrow(bs_data)))
    train_data <- bs_data[train_idx, ]
    test_data <- bs_data[-train_idx, ]
    
    # Train random forest models
    model_richness <- randomForest(fml_richness, data = train_data, ntree = n_tree, mtry = min(2, n_pred_richness))
    model_identity <- randomForest(fml_identity, data = train_data, ntree = n_tree, mtry = min(3, n_pred_identity))
    model_Effect_size <- randomForest(fml_Effect_size, data = train_data, ntree = n_tree, mtry = min(5, n_pred_Effect_size))
    
    # Predictions
    pred_richness <- predict(model_richness, newdata = test_data)
    pred_identity <- predict(model_identity, newdata = test_data)
    pred_Effect_size <- predict(model_Effect_size, newdata = test_data)
    
    # Calculate R² (squared correlation between predicted and observed)
    r2_richness <- cor(pred_richness, test_data$TG_Pot_res)^2
    r2_identity <- cor(pred_identity, test_data$TG_Pot_res)^2
    r2_Effect_size <- cor(pred_Effect_size, test_data$TG_Pot_res)^2
    
    idx <- (i-1)*3 + 1:3
    results_rf$R2_RF[idx] <- c(r2_richness, r2_identity, r2_Effect_size)
  }, error = function(e) cat("Error in iteration", i, ":", e$message, "\n"))
}

results_rf$Model <- factor(results_rf$Model, levels = c("Species_richness", "Species_identity", "Effect_size"))
rf_summary <- results_rf %>%
  group_by(Model) %>%
  summarize(CI.low = quantile(R2_RF, 0.025, na.rm = TRUE),
            Mean = mean(R2_RF, na.rm = TRUE),
            Median = median(R2_RF, na.rm = TRUE),
            CI.high = quantile(R2_RF, 0.975, na.rm = TRUE),
            SD = sd(R2_RF, na.rm = TRUE))
print(rf_summary)
~~~
### Expected output

~~~
> ##################################################################################################################
> #####                                                                                                        ##### 
> #####               Part 4---Random forest analysis (using only actual experimental pot data)                #####
> #####                                                                                                        #####  
> ##################################################################################################################
> # ------------------------------------------------------------------------------
> #  Step1: – Random forest analysis
> # ------------------------------------------------------------------------------
> set.seed(123)
> # filter out the duplicates - each pot should only be used once in the analysis
> df.rf <- df[df$remark %in% levels_rich, ] 
> df.rf <- merge(df.rf, pot_predictions[, c("Pot", "pred_RichnessPSF_Pot", "pred_BiomassRatioPSF_Pot")], 
+                by = "Pot", all.x = TRUE) # Merge predicted values
> rownames(df.rf) <- df.rf$Row.names
> df.rf$Row.names <- NULL
> # Define random forest model formulas
> fml_richness <- formula(TG_Pot_res ~ Richness_con)   # Only richness
> fml_identity <- formula(paste("TG_Pot_res ~ Richness_con +", paste(stressors, collapse = " + ")))  # Richness + species identity
> fml_Effect_size <- formula(paste("TG_Pot_res ~ Richness_con +", paste(stressors, collapse = " + "), 
+                                  "+ pred_RichnessPSF_Pot + pred_BiomassRatioPSF_Pot"))  # Richness + identity + effect size
> n_tree <- 1000  # Number of trees per forest
> n_iter_rf <- 1000  # Bootstrap iterations for RF
> cv_prop <- 0.7  # Training set proportion
> # Number of predictors in each model
> n_pred_richness <- length(all.vars(fml_richness)) - 1
> n_pred_identity <- length(all.vars(fml_identity)) - 1
> n_pred_Effect_size <- length(all.vars(fml_Effect_size)) - 1
> # Data frame to store results
> results_rf <- data.frame(Iteration = rep(1:n_iter_rf, each = 3),
+                          Model = rep(c("Species_richness", "Species_identity", "Effect_size"), n_iter_rf),
+                          R2_RF = NA, stringsAsFactors = FALSE)
> # Split samples by richness level (for stratified bootstrap)
> id_lv1 <- which(df.rf$Richness_con == 1)
> id_lvh <- which(df.rf$Richness_con != 1)
> n_lv1 <- length(id_lv1); n_lv1     #********************* 140
[1] 140
> n_lvh <- length(id_lvh); n_lvh      #********************* 144
[1] 144
> prog_interval <- ceiling(n_iter_rf/10)
> for (i in 1:n_iter_rf) {
+   if (i %% prog_interval == 0) cat(sprintf("Completed %d%% (%d/%d iterations)\n", round(i/n_iter_rf*100), i, n_iter_rf))
+   tryCatch({
+     # Stratified bootstrap sampling (preserve proportion of richness 1 and >1)
+     bs_indices <- c(sample(id_lv1, ceiling(0.5*n_lv1), replace = TRUE),
+                     sample(id_lvh, ceiling(0.5*n_lvh), replace = TRUE))
+     bs_data <- df.rf[bs_indices, ]
+     train_idx <- sample(nrow(bs_data), round(cv_prop * nrow(bs_data)))
+     train_data <- bs_data[train_idx, ]
+     test_data <- bs_data[-train_idx, ]
+     
+     # Train random forest models
+     model_richness <- randomForest(fml_richness, data = train_data, ntree = n_tree, mtry = min(2, n_pred_richness))
+     model_identity <- randomForest(fml_identity, data = train_data, ntree = n_tree, mtry = min(3, n_pred_identity))
+     model_Effect_size <- randomForest(fml_Effect_size, data = train_data, ntree = n_tree, mtry = min(5, n_pred_Effect_size))
+     
+     # Predictions
+     pred_richness <- predict(model_richness, newdata = test_data)
+     pred_identity <- predict(model_identity, newdata = test_data)
+     pred_Effect_size <- predict(model_Effect_size, newdata = test_data)
+     
+     # Calculate R² (squared correlation between predicted and observed)
+     r2_richness <- cor(pred_richness, test_data$TG_Pot_res)^2
+     r2_identity <- cor(pred_identity, test_data$TG_Pot_res)^2
+     r2_Effect_size <- cor(pred_Effect_size, test_data$TG_Pot_res)^2
+     
+     idx <- (i-1)*3 + 1:3
+     results_rf$R2_RF[idx] <- c(r2_richness, r2_identity, r2_Effect_size)
+   }, error = function(e) cat("Error in iteration", i, ":", e$message, "\n"))
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
> results_rf$Model <- factor(results_rf$Model, levels = c("Species_richness", "Species_identity", "Effect_size"))
> rf_summary <- results_rf %>%
+   group_by(Model) %>%
+   summarize(CI.low = quantile(R2_RF, 0.025, na.rm = TRUE),
+             Mean = mean(R2_RF, na.rm = TRUE),
+             Median = median(R2_RF, na.rm = TRUE),
+             CI.high = quantile(R2_RF, 0.975, na.rm = TRUE),
+             SD = sd(R2_RF, na.rm = TRUE))
> print(rf_summary)
# A tibble: 3 × 6
  Model               CI.low   Mean Median CI.high     SD
  <fct>                <dbl>  <dbl>  <dbl>   <dbl>  <dbl>
1 Species_richness 0.0000197 0.0326 0.0176   0.156 0.0420
2 Species_identity 0.00200   0.144  0.123    0.416 0.112 
3 Effect_size      0.00452   0.162  0.144    0.426 0.115 
~~~
# License
This project is covered under the GNU GENERAL PUBLIC LICENSE Version 2 (GPL-2).
