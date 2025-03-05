# Overview
The following files allow one to reproduce analyses in the manuscript entitled "Species-level studies do not upscale to community-wide plant-soil feedbacks".
# Repo Contents
## Code: R code
This script contains all the data analysis and visualization code of Fig.1, Fig.2, Fig.3 and Fig.S2.
## Data  
This script contains all the data of Fig.1, Fig.2, Fig.3 and Fig.S2.

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
install.packages(c('openxlsx', 'ggfortify', 'ggforce', 'dplyr', 'ggrepel', 'ggpubr', 'export', 'ggspatial', 'ggplot2','cowplot', 'sf', 'multcomp', 'tidyverse',
 'car', 'lme4', 'vegan', 'FUNGuildR', 'lmerTest', 'ggbeeswarm', 'reshape2', 'DHARMa', 'glmmTMB', 'emmeans', 'broom.mixed', 'MuMIn', 'ggridges', 'party', 'caret', 'pheatmap'))
~~~

Note: ggepi & patchwork can be installed via devtools 
~~~
 devtools::install_github("lwjohnst86/ggepi")
 devtools::install_github("thomasp85/patchwork")
~~~

The package should take approximately 50-60 seconds to install with vignettes on a recommended computer.

# Demo
## Instructions to run on data: example-Fig.3 code
~~~
##################################################################################
#####                                                                        ##### 
#####           Part0---effect of home vs away:PSF                           #####
#####                                                                        #####  
##################################################################################

setwd('C:/Users/MY/Desktop/CODE/Fig.3')

###Table S9
 
dat <- read_excel("data_Fig.3.xlsx",sheet=1);head(dat)
#Fit lm model with totol biomass
mod_full <- lm(TG ~ Richness_con * (B_con + C_con + F_con + T_con + V_con), data = dat)
qqnorm(resid(mod_full));qqline(resid(mod_full));anova(mod_full)-> mod_full_result; mod_full_result
p <- mod_full_result$Pr;p 
p.adjust(p, "BH")
p.adjust(p, "bonferroni")
~~~
## Expected output
~~~
Response: TG
Df  Sum Sq Mean Sq F value   Pr(>F)   
Richness_con         5  2506.1  501.22  4.1830 0.001133 ** 
B_con                1    18.2   18.18  0.1518 0.697183   
C_con                1   449.8  449.76  3.7535 0.053790 .
F_con                1   103.9  103.94  0.8675 0.352527   
T_con                1   123.9  123.92  1.0342 0.310132   
V_con                1     9.2    9.19  0.0767 0.782079   
Richness_con:B_con   4   150.3   37.59  0.3137 0.868707   
Richness_con:C_con   4   506.4  126.60  1.0565 0.378582   
Richness_con:F_con   4   158.8   39.71  0.3314 0.856731   
Richness_con:T_con   4   457.6  114.40  0.9548 0.432939   
Richness_con:V_con   4   257.2   64.31  0.5367 0.708894   
Residuals          258 30914.8  119.82  

p.adjust(p, "bonferroni")
[1] 0.01246636 1.00000000 0.59168644 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
[9] 1.00000000 1.00000000 1.00000000         NA
~~~
## Note: 
The dataset stored in this repository is same to the dataset in figshare(https://doi.org/10.6084/m9.figshare.28395707).
Running each line of code on a recommended computer takes about 1-5 seconds, except for "Step 2: Predictability comparison among random forest models" in Fig.2's code, which takes 20-30 minutes.

# License
This project is covered under the GNU GENERAL PUBLIC LICENSE Version 2 (GPL-2).
