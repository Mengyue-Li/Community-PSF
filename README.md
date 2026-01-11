# Overview
The following files allow one to reproduce analyses in the manuscript entitled "Species-level studies do not upscale to community-wide plant-soil feedbacks".
# Repo Contents
## Code: R code
This script contains all the data analysis and visualization code of Fig.1, Fig.3&Table S2-48, Fig.S1&2.
## Data  
This script contains all the data of Fig.1, Fig.3&Table S2-48, Fig.S1&2.

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
install.packages(c('readxl', 'openxlsx', 'ggfortify', 'ggforce', 'dplyr', 'ggrepel', 'ggpubr', 'export', 'ggspatial', 'ggplot2','cowplot', 'sf', 'multcomp', 'tidyverse',
 'car', 'lme4', 'vegan', 'FUNGuildR', 'lmerTest', 'ggbeeswarm', 'reshape2', 'DHARMa', 'glmmTMB', 'emmeans', 'broom.mixed', 'MuMIn', 'ggridges', 'party', 'caret', 'pheatmap'))
~~~

Note: ggepi & patchwork can be installed via devtools 
~~~
 devtools::install_github("lwjohnst86/ggepi")
 devtools::install_github("thomasp85/patchwork")
~~~

The package should take approximately 50-60 seconds to install with vignettes on a recommended computer.

# Demo
## Instructions to run on data: example-“Fig.3&Table S5-S7 S9” code
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
~~~
## Note: 
The dataset stored in this repository is same to the dataset in figshare( ).
Running each line of code on a recommended computer takes about 1-5 seconds, except for "Step 2: Predictability comparison among random forest models" in Fig.2's code, which takes 5-10 minutes.

# License
This project is covered under the GNU GENERAL PUBLIC LICENSE Version 2 (GPL-2).
