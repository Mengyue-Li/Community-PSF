
setwd('D:/2025.10.4-NPH/1-0330/2-NEWrawdata-code') 

###########################################################################################
#####                                  Library & mytheme                              #####
###########################################################################################
### Load required packages and data
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(ggepi)
library(party)
library(caret)
library(glmmTMB)
library(dbplyr)
library(tidyr)
library(tidyverse)
library(randomForest)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(ape)
library(dplyr)
library(vegan)
library(picante)
library(RColorBrewer)
library(ggsci)
library(ggpmisc)
library(multcompView)
library(multcomp)


mytheme = theme( panel.background = element_rect(fill='white', colour='black'),legend.position = "none",
                 panel.grid=element_blank(), legend.title = element_blank(),legend.text = element_text(size = 11),
                 legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                 axis.ticks = element_line(color='black'),
                 axis.line = element_line(colour = "black"), 
                 axis.title.x=element_text(colour='black', size=13,vjust = 1),
                 axis.title.y=element_text(colour='black', size=13,vjust = 1),
                 axis.text=element_text(colour='black',size=11),
                 plot.tag = element_text(size = 14, face = "bold"))


## ================================================================================================================================== ##
#                                                                                                                                      #
#                                                 Fig.S3： species-level mixed model                                                        #
#                                                                                                                                      #
## ================================================================================================================================== ##

##################################################################################################################
#####                                                                                                        ##### 
#####              Part 1--- Calculate actual PSF and predicted PSF (Pot-level & Species-level)              #####
#####                                                                                                        #####  
##################################################################################################################

#--------------------------------------
### Step 1: data 
#--------------------------------------

df1 <- read_excel("Fig3 & S3 & TableS1-6.xlsx", sheet = "Fig3_data")

df1 <- as.data.frame(df1); df1

# Define species column names (consistent with Fig3_data)
stressors <- c("Sor", "Aa", "Sv", "Ul", "Vo", "St", "Cal", "Sn", "Ai", "Pf", "Pb", "Md", "Cab", "Mc", 
               "Car", "Sa", "Pa", "Ss", "Sc", "Ah", "Ds", "Soc", "Da", "Av", "At", "Cp", "Cc", "Cs")

# Define list of response variables to compute
responses_vec <- c("TG_Pot_res", 
                   "Pc_averageAG_res", "Aa_averageAG_res", "Sv_averageAG_res", 
                   "St_averageAG_res", "Ah_averageAG_res", "Soc_averageAG_res")

# List to store all results (will be converted to wide format later)
all_results_list <- list()


# ------------------------------------------------------------------------------
#  Step 2: Loop over each response variable
# ------------------------------------------------------------------------------
for (resp in responses_vec) {
  cat("\nProcessing response variable:", resp, "\n")
  
  # 1) Extract response values of the 5 CK controls
  CK_vals <- df1[df1$remark == "CK", resp]
  if (length(CK_vals) != 5) warning(paste("CK_vals length not 5 for", resp))
  
  # 2) Define function to compute PSF for each single species (treat 0 as NA)
  compute_PSF <- function(y, CK_vals) {
    if (is.na(y) || y <= 0) return(NA)   # 0 or negative values are treated as missing
    mean(log(y / CK_vals), na.rm = TRUE)
  }
  
  # 3) Calculate mean PSF for each single species (species_psf)
  species_psf <- sapply(stressors, function(sp) {
    tr_vals <- df1[df1$remark == sp, resp]
    if (length(tr_vals) == 0) return(NA)
    # Calculate PSF for each replicate, ignoring NAs from 0
    psf_vals <- sapply(tr_vals, function(t) compute_PSF(t, CK_vals))
    mean(psf_vals, na.rm = TRUE)
  })
  names(species_psf) <- stressors
  
  # 4) Extract all treated pots (remark in richness levels)
  levels <- c("1", "2", "4", "6", "12")
  pots <- df1[df1$remark %in% levels, ]
  
  # 5) Actual PSF  (also treat 0 as NA)
  pots$actual_PSF <- sapply(pots[, resp], function(t) compute_PSF(t, CK_vals))
  
  # 6)  Richness PSF & Biomass-Ratio PSF
  comp_cols <- which(names(pots) %in% stressors)
  pots$pred_RichnessPSF <- NA          # RichnessPSF prediction (equal-weight average)
  pots$pred_BiomassRatioPSF <- NA      # BiomassRatioPSF prediction (biomass-weighted average)
  
  for (i in 1:nrow(pots)) {
    sp_present <- stressors[ pots[i, comp_cols] > 0 ]# Species present in the pot
    if (length(sp_present) == 0) next
    psfj <- species_psf[sp_present]
    # Remove NAs (species for which PSF could not be calculated)
    valid <- !is.na(psfj)
    if (sum(valid) == 0) next
    psfj <- psfj[valid]
    sp_present <- sp_present[valid]
    
    # RichnessPSF: Simple average of species PSFs
    pots$pred_RichnessPSF[i] <- mean(psfj, na.rm = TRUE)
    
    # BiomassRatioPSF (biomass-weighted average)
    biomass <- unlist(pots[i, sp_present])
    if (sum(biomass, na.rm = TRUE) > 0) {
      weights <- biomass / sum(biomass)
      pots$pred_BiomassRatioPSF[i] <- sum(psfj * weights, na.rm = TRUE)
    } else {
      pots$pred_BiomassRatioPSF[i] <- mean(psfj, na.rm = TRUE)
    }
  }
  
  # 7) Build sub-table for this response variable (long format, to be merged later)
  temp <- data.frame(
    Pot = pots$Pot,
    Remark = pots$remark,
    Richness = pots$Richness_con,
    Response = resp,
    Actual = pots$actual_PSF,
    pred_RichnessPSF = pots$pred_RichnessPSF,
    BiomassRatioPSF = pots$pred_BiomassRatioPSF,
    stringsAsFactors = FALSE
  )
  all_results_list[[resp]] <- temp
}


# ------------------------------------------------------------------------------
#  Step 3: Combine results 
# ------------------------------------------------------------------------------
all_results_long <- bind_rows(all_results_list)

# Convert to wide format using pivot_wider
all_results_wide <- all_results_long %>%
  pivot_wider(
    id_cols = c(Pot, Remark, Richness),
    names_from = Response,
    values_from = c(Actual, pred_RichnessPSF, BiomassRatioPSF),
    names_sep = "_"
  )

# Adjust column order: first Pot, Remark, Richness, then Actual, pred_RichnessPSF, BiomassRatioPSF for each response variable
resp_names <- responses_vec
ordered_cols <- c("Pot", "Remark", "Richness")
for (r in resp_names) {
  ordered_cols <- c(ordered_cols, paste0("Actual_", r), 
                    paste0("pred_RichnessPSF_", r), 
                    paste0("BiomassRatioPSF_", r))
}

# Ensure all columns exist (some may be missing, but theoretically all should be present)
all_results_wide <- all_results_wide[, intersect(ordered_cols, names(all_results_wide))]


# Output results
#write.csv(all_results_wide, "PSF-all-20260518_wide.csv", row.names = FALSE)
head(all_results_wide)



#########################################################################################################################
#####                                                                                                               ##### 
#####    Part 2---Perform regression for each of the two hypotheses against the observed values(Species-level)      #####
#####                                                                                                               #####  
#########################################################################################################################
##--------------------------------------
### Step 1: data 
#--------------------------------------  
data_regress <- read_excel("Fig3 & S3 & TableS1-6.xlsx", sheet = "data_regress")   
glimpse(data_regress)

data_regress <- data_regress %>%
  mutate(
    ActualPSF = as.numeric(ActualPSF),
    RichnessPSF = as.numeric(RichnessPSF),
    BiomassRatioPSF = as.numeric(BiomassRatioPSF),
    Phylo_Dist = as.numeric(Phylo_Dist),
    Richness_con = as.factor(Richness_con),
    Pot = as.factor(Pot),
    Species = as.factor(Species),
    Species_full = as.factor(Species_full),
    Species_con  =  as.factor(Species_con)
  ) %>%
  filter(complete.cases(ActualPSF, RichnessPSF, BiomassRatioPSF, Richness_con, Species, Pot))

# Define species level order and labels
sp_levels <- c("Pc", "Aa", "Sv", "St", "Ah", "Soc")
sp_labels <- c(Pc = "Paspalum_conjugatum", Aa = "Achyranthes_aspera",Sv = "Setaria_viridis",
               St = "Senna_tora", Ah = "Amaranthus_hybridus", Soc = "Senna_occidentalis")
data_regress$Species <- factor(data_regress$Species, levels = sp_levels)
data_regress$Species_full <- factor(data_regress$Species_full, levels = sp_labels)


# ------------------------------------------------------------------------------
#  Step 2: Model 1: RichnessPSF hypothesis --- TableS4
# ------------------------------------------------------------------------------
# Full model
mR_full <- lmer(ActualPSF ~ Species * Richness_con * RichnessPSF + (1 | Pot),data = data_regress, REML = FALSE)
anova_result <- anova(mR_full)
p_table <- anova_result; print(p_table)
p_values <- anova_result$`Pr(>F)`;print(p_values) 
p.adjust(p_values, "BH")


# slopes per species for BiomassRatioPSF
slope_R <- as.data.frame(emtrends(mR_full, ~ 1, var = "RichnessPSF", infer = c(TRUE, TRUE))); slope_R 
slope_R <- slope_R %>% rename(slope = RichnessPSF.trend, slope_SE = SE, slope_df = df,
                              slope_lower = lower.CL, slope_upper = upper.CL,slope_t = t.ratio, slope_p = p.value) ; print(slope_R)

# Figures
xR <- seq(min(data_regress$RichnessPSF, na.rm = TRUE),  max(data_regress$RichnessPSF, na.rm = TRUE), length.out = 100)
grid_R <- emmeans(mR_full, ~ RichnessPSF, at = list(RichnessPSF = xR))
plot_R <- as.data.frame(grid_R)   
dat_R <- data_regress %>% filter(!is.na(ActualPSF), !is.na(RichnessPSF))

p_R <- ggplot() +
  #geom_hline(yintercept = 0, linetype = 2, colour = "grey65") +
  geom_point(data = dat_R, aes(x = RichnessPSF, y = ActualPSF),alpha = 0.3, size = 1.5, colour = "#999999") + 
  geom_line(data = plot_R, aes(x = RichnessPSF, y = emmean), colour = "#737373", linewidth = 1.2) +
  theme_classic() + theme_bw() + mytheme +
  labs(x = "Species-specific predicted PSF\n(Richness hypothesis)",
       y = "Observed species-specific PSF")
print(p_R)


# ------------------------------------------------------------------------------
#  Step 3: Model 2: BiomassRatioPSF hypothesis --- TableS4
# ------------------------------------------------------------------------------
mB_full <- lmer(ActualPSF ~ Species * Richness_con * BiomassRatioPSF + (1 | Pot), data = data_regress, REML = FALSE)
anova_result <- anova(mB_full)
p_table <- anova_result; print(p_table)
p_values <- anova_result$`Pr(>F)`;print(p_values) 
p.adjust(p_values, "BH")

# slopes per species for BiomassRatioPSF
#mB_full_emtrends <- emtrends(mB_full, pairwise ~ Richness_con, var = "BiomassRatioPSF")
#test(mB_full_emtrends, adjust = "BH")
#test(mB_full_emtrends)

slope_B <- as.data.frame(emtrends(mB_full, ~ Richness_con, var = "BiomassRatioPSF", infer = c(TRUE, TRUE)));slope_B
slope_B <- slope_B %>% rename(slope = BiomassRatioPSF.trend, slope_SE = SE, slope_df = df,
                              slope_lower = lower.CL, slope_upper = upper.CL,slope_t = t.ratio, slope_p = p.value) %>%
                              mutate(sig = ifelse(slope_p < 0.05, "significant", "not significant")); slope_B 

# Figures
xB <- seq(min(data_regress$BiomassRatioPSF, na.rm = TRUE),max(data_regress$BiomassRatioPSF, na.rm = TRUE), length.out = 100); xB 
grid_B <- emmeans(mB_full, ~ Richness_con * BiomassRatioPSF, at = list(BiomassRatioPSF = xB, Richness_con = levels(data_regress$Richness_con))); grid_B 
plot_B <- as.data.frame(grid_B) %>% left_join(slope_B %>% dplyr::select(Richness_con, sig), by = "Richness_con"); plot_B
dat_B <- data_regress %>% filter(!is.na(ActualPSF), !is.na(BiomassRatioPSF)); dat_B

p_B <- ggplot() +
  geom_point(data = dat_B, aes(x = BiomassRatioPSF, y = ActualPSF,colour = Richness_con), alpha = 0.4, size = 1.5,) +
  geom_line(data = plot_B, aes(x = BiomassRatioPSF, y = emmean,colour = Richness_con,linetype = sig),linewidth = 1.2) +
  scale_linetype_manual(values = c("significant" = "solid", "not significant" = "dashed")) +
  scale_color_manual(values = c("#3f3e93", "#7ca6e5", "#fedb81",  "#ff783e", "#c04037"),name = "Conditioning species richness") +
  theme_classic() + theme_bw() + mytheme  +
  labs(x = "Species-specific predicted PSF\n(biomass-ratio hypothesis)",
       y = "Observed species-specific PSF") + theme(legend.position = "right")
print(p_B)
combined <- p_R | p_B ; combined 

 



## ================================================================================================================================== ##
#                                                                                                                                      #
#                                              3： phylogenetic signal analysis                                                        #
#                                                                                                                                      #
## ================================================================================================================================== ##
##################################################################################
#####                                                                        ##### 
#####                  Part 1---Calculation of wPDi                          #####
#####                                                                        #####  
##################################################################################
#--------------------------------------
### Step 1: data 
#-------------------------------------- 
### conditioning plant community
Condition_com = read.xlsx("Fig3 & S3 & TableS1-6.xlsx", sheet = "Species_AG_con", colNames = T, rowNames = T)
Condition_com[1:6,1:6]
rel_Condition_com <- decostand(Condition_com, method="total", MARGIN=1) ### Relative biomass
rowSums(rel_Condition_com)
rel_Condition_com [1:6,1:6]


#--------------------------------------
### Step 2: Calculation of wPDi  
#-------------------------------------- 
#0 Plant phylogenetic tree
tree = read.tree("Plant_tree.treefile"); plot(tree)
# Remove outergroup
to_drop <- c("Amborella_trichopoda","")
tree <- drop.tip(as.phylo(tree) , to_drop); plot(tree)

#1 Phylogenetic distance matrix                
plant_dis = cophenetic(tree)  
dis = as.matrix(plant_dis)
#View(dis)

#2 relative biomass table of conditioning community
samp = rel_Condition_com
samp[1:6,1:6]

#3 six species of responding community
Respond_sp <- c(  "Paspalum_conjugatum", "Achyranthes_aspera",	"Setaria_viridis",
                  "Senna_tora",	"Amaranthus_hybridus",	"Senna_occidentalis"); Respond_sp

#4 wPDi: weighted phylogenetic distance between responding species and conditioning communities
N <- dim(samp)[1]
Distance <- NULL

rownames(dis)
colnames(dis)
Respond_sp %in% rownames(dis)

for (i in 1:N) {
  sppInSample <- names(samp)[samp[i, ] > 0]
  sample.dis <- as.data.frame(t(dis[Respond_sp, sppInSample]))
  sample.weights <- as.data.frame(t(samp[i, sppInSample, drop = FALSE]))
  for (A in 1:5) {
    sample.weights <- cbind(sample.weights, sample.weights[ ,1])
  }
  colnames(sample.weights) = colnames(sample.dis)
  distance_group = as.data.frame(t(data.frame(colSums(sample.dis*sample.weights)/colSums(sample.weights)  )))
  rownames(distance_group) = NULL
  distance_group$Sample_ID = rownames(samp[i, ])
  Distance = rbind(Distance, distance_group)
}

colnames(Distance) = c("Pc_dis", "Aa_dis", "Sv_dis", "St_dis", "Ah_dis", "So_dis", "Sample_ID")
head(Distance)
# write.csv(Distance,"Phylo-Dist.csv")



##################################################################################
#####                                                                        ##### 
#####                            Part 2--- MOdel                             #####
#####                                                                        #####  
##################################################################################
#--------------------------------------
### Step 1: data 
#--------------------------------------  
data_regress <- read_excel("Fig3 & S3 & TableS1-6.xlsx", sheet = "data_regress")   
glimpse(data_regress)

data_regress <- data_regress %>%
  mutate(
    ActualPSF = as.numeric(ActualPSF),
    RichnessPSF = as.numeric(RichnessPSF),
    BiomassRatioPSF = as.numeric(BiomassRatioPSF),
    Phylo_Dist = as.numeric(Phylo_Dist),
    Richness_con = as.factor(Richness_con),
    Pot = as.factor(Pot),
    Species = as.factor(Species),
    Species_full = as.factor(Species_full),
    Species_con  = as.factor(Species_con)                 
  ) %>%
  filter(complete.cases(ActualPSF, RichnessPSF, BiomassRatioPSF, Richness_con, Species, Pot))


#--------------------------------------
### Step 2: Model 1  --- TableS6
#-------------------------------------- 
data_regress 

mP_full <- lmer(ActualPSF ~ Species * Richness_con * Phylo_Dist + (1 | Pot), data = data_regress, REML = FALSE)
anova_result <- anova(mP_full)
p_table <- anova_result; print(p_table)
p_values <- anova_result$`Pr(>F)`; print(p_values) 
p.adjust(p_values, "BH")  


#--------------------------------------
### Step 3: Model 2 --- TableS5
#-------------------------------------- 
mC_full <- lmer(ActualPSF ~  Species * Richness_con * Species_con + (1 | Pot), data = data_regress, REML = FALSE)
anova_result <- anova(mC_full)
p_table <- anova_result; print(p_table)
p_values <- anova_result$`Pr(>F)`; print(p_values) 
p.adjust(p_values, "BH")


#--------------------------------------
### Step 4: Model 3 --- TableS6
#--------------------------------------
### general phylogenetic effect among heterospecific-only pots
data_regress_heter <- data_regress %>% filter(Species_con == "0")

mP_full <- lmer(ActualPSF ~ Species * Richness_con * Phylo_Dist + (1 | Pot), data = data_regress_heter, REML = FALSE)
anova_result <- anova(mP_full)
p_table <- anova_result; print(p_table)
p_values <- anova_result$`Pr(>F)`; print(p_values) 
p.adjust(p_values, "BH")  


