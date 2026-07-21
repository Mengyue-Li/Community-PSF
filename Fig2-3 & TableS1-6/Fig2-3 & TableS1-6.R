

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
#                                                   1：Fig.2 & Table S1-S3                                                             #
#                                                                                                                                      #
## ================================================================================================================================== ##
##################################################################################
#####                                                                        ##### 
#####           Part1---effect of home vs away:PSF                           #####
#####                                                                        #####  
##################################################################################

setwd('D:/2025.10.4-NPH/1-0330/1-Revised manuscript 2026/Code-0715/Fig2-3 & TableS1-6') 

#-------------------------------------------------------------------------------------------- 
### Table S1 --- responding community aboveground biomass & responding community PSFs
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
mod_full2 <- lm(ObservedPSF_Pot ~ Richness_con * (Aa_con + Sv_con + St_con + Ah_con + Soc_con), data = data)
anova(mod_full2)-> mod_full2_result; mod_full2_result
p2 <- mod_full2_result$Pr;p2 
p.adjust(p2, "BH")
#p.adjust(p2, "bonferroni")

#--------------------------------------
### Relationship between aboveground biomass of conditioning species or communities and PSFs 
#--------------------------------------
data[1:6,1:12]
### Relationship between species richness and above-ground biomass during the conditioning phase
# F4, 279 = 0.32, P = 0.864
data$Richness_con = as.factor(data$Richness_con)
mod_01 <- lm(AG_con ~ Richness_con, data = data)
anova(mod_01)-> mod_01_result;mod_01_result 
p <- mod_01_result$Pr;p 
p.adjust(p, "BH")


##################################################################################################################
#####                                                                                                        ##### 
#####   Part 1--- Biomass density distributions and means for each treatment (single species and richness)   #####
#####                                                                                                        #####  
##################################################################################################################
#--------------------------------------
### Step 1: data 
#-------------------------------------- 

df <- read_excel("data-Fig2-3 & TableS1-6.xlsx", sheet = "Fig2_data")[, 1:32]

df <- as.data.frame(df)
 
df$TG_Pot_res <- as.numeric(df$TG_Pot_res)
df$Richness_con <- as.factor(df$Richness_con)   
 
stressors <- c("Sor", "Aa", "Sv", "Ul", "Vo", "St", "Cal", "Sn", "Ai", "Pf", "Pb", "Md",  "Cab", "Mc",
               "Car", "Sa", "Pa", "Ss", "Sc", "Ah", "Ds", "Soc", "Da", "Av", "At", "Cp", "Cc", "Cs")
levels_rich <- c("1", "2", "4", "6", "12")
target <- c("CK", stressors, levels_rich)
treatment = as.vector(unique(df$remark)); treatment

n_iter <- 1000
responses <- "TG_Pot_res"

set.seed(123)


#--------------------------------------
### Step 2: Function 
#-------------------------------------- 
# ---------------------------------------------------------------------------------------------   
# Function 1: Estimating mean and its 95% confidence interval
# --------------------------------------------------------------------------------------------- 
BootStrap_mean <- function(response, data = df, n_perm = n_iter) {
  summary <- list()
  for (treatment in target) {
    bs <- numeric(0)
    if (treatment == "1") {
      
      population <- data[data$remark %in% stressors, response]
    } else {
      population <- data[data$remark == treatment, response]
    }
    size <- length(population)
    if (size > 0) {
      for (id in 1:n_perm) {
        k <- mean(sample(population, size, replace = TRUE), na.rm = TRUE)
        bs <- append(bs, k)
      }
      summary[[treatment]] <- c(quantile(bs, 0.025, na.rm = TRUE), 
                                mean(bs, na.rm = TRUE), 
                                quantile(bs, 0.975, na.rm = TRUE))
      names(summary[[treatment]]) <- c("2.5%", "mean", "97.5%")
    } else {
      summary[[treatment]] <- c(NA, NA, NA)
      names(summary[[treatment]]) <- c("2.5%", "mean", "97.5%")
    }
  }
  summary <- t(data.frame(summary))
  summary <- data.frame("target" = target, summary)
  row.names(summary) <- NULL
  return(summary)
}

response_mean <- BootStrap_mean(responses)
 
# Set factor order (for plotting)
response_mean$target <- factor(response_mean$target, levels = c("CK", stressors, levels_rich))
df$remark <- factor(df$remark, levels = levels(response_mean$target))

# ---------------------------------------------------------------------------------------------   
# Function 2: Bootstrap effect size (difference between treatment and control) resampling
# ---------------------------------------------------------------------------------------------  
BootStrap_ES_rep <- function(response, data = df, target_ = unique(df$remark), n_perm = n_iter) {
  resampled <- list()
  population_CK <- data[data$remark == "CK", response]
  for (treatment in target) {
    bs <- numeric(0)
    if (treatment == "1") {
      population_TR <- data[data$remark %in% stressors, response]
    } else if (treatment != "1") {
      population_TR <- data[data$remark == treatment, response]
    }
    size_CK <- length(population_CK)
    size_TR <- length(population_TR)
    for (id in 1:n_perm) {
      k_CK <- mean(sample(population_CK, size_CK, replace = TRUE))
      k_TR <- mean(sample(population_TR, size_TR, replace = TRUE))
      bs <- append(bs, k_TR - k_CK)
    }
    resampled[[treatment]] <- bs
  }
  resampled[["CK"]] <- rep(0, n_perm)
  return(resampled)
}

# ---------------------------------------------------------------------------------------------   
#  Function 3: Summarize bootstrap effect sizes, calculate mean, CI, and p-value
# ---------------------------------------------------------------------------------------------  
BootStrap_ES_summary = function(data){
  summary = list()
  p = 0
  summary[["CK"]] = c(0,0,0,1)
  target = names(data)
  for(treatment in target[-1]){
    bs = data[[treatment]]
    p = length(which(bs>0))/length(bs)
    p = min(p, 1-p)
    summary[[treatment]] = c(quantile(bs, .025), mean(bs), quantile(bs, .975), p)
  }
  summary = t(data.frame(summary))
  colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
  summary = data.frame(target, summary); row.names(summary) = c()
  
  return(summary)
}

# ---------------------------------------------- 
# Step 3: Perform effect size estimation
# ---------------------------------------------- 
set.seed(123)
response_mean_all <- list()
response_ES_all   <- list()
for (i_response in responses) {
  response_mean   <- BootStrap_mean(i_response)
  response_ES_bs  <- BootStrap_ES_rep(i_response)
  response_ES     <- BootStrap_ES_summary(response_ES_bs)
  response_mean_all[[i_response]] <- response_mean
  response_ES_all[[i_response]]   <- response_ES
}


# ----------------------------------------------- 
# Step 4: Table S2
# ----------------------------------------------- 
response_mean 
response_ES_all
#write.csv(response_mean, "TableS3-1-20260518.csv", row.names = FALSE)
#write.csv(response_ES_all, "TableS3-2-20260518.csv", row.names = FALSE)


# ------------------------------------------------------------------------------------------
# Step 5: Figure 3(a) – Biomass density distributions and means for each treatment
# ------------------------------------------------------------------------------------------
mytheme <- theme()   
a <- ggplot() +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  #  "#D3D3D3" region for richness levels 1 to 12
  geom_rect(data = NULL, 
            aes(xmin = -Inf, xmax = Inf, 
                ymin = which(levels(df$remark) == "1") - 0.45, 
                ymax = which(levels(df$remark) == "12") + 1), 
            fill = "#D3D3D3", alpha = 0.4, inherit.aes = FALSE) +
  #  Control mean vertical line
  geom_vline(aes(xintercept = subset(response_mean, target == "CK")$mean),
             linetype = "dashed", color = "black") +
  coord_flip() +
  # Density ridge for control
  stat_density_ridges(data = subset(df, remark == "CK"), 
                      aes(x = TG_Pot_res, y = remark, fill = after_stat(x)), 
                      geom = "density_ridges_gradient", bandwidth = 5.82, rel_min_height = 0.01, 
                      jittered_points = TRUE, color = alpha("#4E307D", 0.5), fill = alpha("#4E307D", 0.3), 
                      alpha = 0.4, position = position_points_jitter(height = 0.15, yoffset = 0.1), scale = 0.6) +
  geom_point(data = subset(response_mean, target == "CK"), 
             aes(x = mean, y = target), shape = 1, size = 2, stroke = 1.2, color = "#4E307D") +
  geom_errorbar(data = subset(response_mean, target == "CK"), 
                aes(xmin = X2.5., xmax = X97.5., y = target), width = 0, color = "#4E307D", linewidth = 0.5) +
  # 1) Single species treatments (excluding richness levels)
  stat_density_ridges(data = subset(df, !remark %in% c("CK", levels_rich)), 
                      aes(x = TG_Pot_res, y = remark, fill = after_stat(x)), 
                      geom = "density_ridges_gradient", bandwidth = 5.02, rel_min_height = 0.01, 
                      jittered_points = TRUE, fill = alpha("#4E307D", 0.2), color = alpha("#4E307D", 0.2), 
                      alpha = 0.4, position = position_points_jitter(height = 0.15, yoffset = 0.1), scale = 0.6) +
  geom_point(data = subset(response_mean, !target %in% c("CK", levels_rich)), 
             aes(x = mean, y = target), shape = 16, size = 3, color = "#4E307D") +
  geom_errorbar(data = subset(response_mean, !target %in% c("CK", levels_rich)), 
                aes(xmin = X2.5., xmax = X97.5., y = target), width = 0, color = "#4E307D", linewidth = 0.5) +
  # 2) Richness levels
  stat_density_ridges(data = subset(df, remark %in% levels_rich), 
                      aes(x = TG_Pot_res, y = remark, fill = after_stat(x)), 
                      geom = "density_ridges_gradient", bandwidth = 4.03, rel_min_height = 0.01, 
                      jittered_points = TRUE, color = alpha("#000000", 0.2), fill = alpha("#000000", 0.2), 
                      alpha = 0.1, position = position_points_jitter(height = 0.15, yoffset = 0.1), scale = 0.6) +
  geom_point(data = subset(response_mean, target %in% levels_rich), 
             aes(x = mean, y = target), shape = 16, size = 3, color = "#000000") +
  geom_errorbar(data = subset(response_mean, target %in% levels_rich), 
                aes(xmin = X2.5., xmax = X97.5., y = target), width = 0, color = "#000000", linewidth = 0.5) +
  labs(x = "Total biomass
of responding community (g)", y = "Single plant species Richness", tag = "(a)") +
  mytheme +  theme_bw() +
  scale_y_discrete(labels = c("CK" = "Control", 
                              "Sor" = "Sor", "Aa" = "Aa", "Sv" = "Sv", "Ul" = "Ul",
                              "Vo" = "Vo", "St" = "St", "Cal" = "Cal", "Sn" = "Sn",
                              "Ai" = "Ai", "Pf" = "Pf", "Pb" = "Pb", "Md" = "Md",
                              "Cab" = "Cab", "Mc" = "Mc", "Car" = "Car", "Sa" = "Sa",
                              "Pa" = "Pa", "Ss" = "Ss", "Sc" = "Sc", "Ah" = "Ah",
                              "Ds" = "Ds", "Soc" = "Soc", "Da" = "Da", "Av" = "Av",
                              "At" = "At", "Cp" = "Cp", "Cc" = "Cc", "Cs" = "Cs",
                              "1" = "1", "2" = "2", "4" = "4", "6" = "6", "12" = "12")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_x_continuous(limits = c(25, 123), expand = c(0, 0), breaks = seq(40, 120, 20))
a



##################################################################################################################
#####                                                                                                        ##### 
#####                              Part 2---PSF prediction for each pot                                      #####
#####                                                                                                        #####  
##################################################################################################################

# --------------------------------------------------------------------------------------
# Step 1: Build bootstrap distribution of PSF for each species (for resampling)
# --------------------------------------------------------------------------------------

CK_vals <- df[df$remark == "CK", "TG_Pot_res"]

species_psf_dist <- lapply(stressors, function(sp) {
  tr_vals <- df[df$remark == sp, "TG_Pot_res"]
  if (length(tr_vals) == 0) return(NA)
  replicate(n_iter, {
    tr_sample <- sample(tr_vals, length(tr_vals), replace = TRUE)
    mean(sapply(tr_sample, function(t) mean(log(t / CK_vals), na.rm = TRUE)))
  })
})
names(species_psf_dist) <- stressors

species_psf <- sapply(stressors, function(sp) {
  tr_vals <- df[df$remark == sp, "TG_Pot_res"]
  if (length(tr_vals) == 0) return(NA)
  pot_psf <- sapply(tr_vals, function(t) mean(log(t / CK_vals), na.rm = TRUE))
  mean(pot_psf, na.rm = TRUE)
})
names(species_psf) <- stressors

 
#--------------------------------------
### Step 2: Function 
#-------------------------------------- 
# ------------------------------------------------------------------------------
# Function 1: Generate null distribution (for each pot at each richness level)
# ------------------------------------------------------------------------------
Null_distribution_new <- function(n_perm = n_iter) {
  output <- list()
  for (Richness_con in levels_rich) {
    resampled <- list()
    if (Richness_con == "1") {
      comp_data <- df[df$remark == "1", c(stressors)] #**************** 140 pot
    } else {
      comp_data <- df[df$remark == Richness_con, c(stressors)] #**************** 144 pot
    }
    for (type in c("Richness_hypothesis", "Biomass_ratio_hypothesis")) {
      resampled[[type]] <- list()
      for (j in 1:nrow(comp_data)) {
        sp_present <- stressors[comp_data[j, ] > 0]
        if (length(sp_present) == 0) {
          resampled[[type]][[j]] <- rep(NA, n_perm)
          next
        }
        if (type == "Biomass_ratio_hypothesis") {
          biomass <- unlist(comp_data[j, sp_present])
          if (sum(biomass) > 0) {
            weights <- biomass / sum(biomass)
          } else {
            weights <- rep(1/length(sp_present), length(sp_present))
          }
        }
        bs <- numeric(n_perm)
        for (id in 1:n_perm) {
          sampled_psf <- sapply(sp_present, function(sp) {
            if (length(species_psf_dist[[sp]]) > 0) sample(species_psf_dist[[sp]], 1) else NA
          })
          if (all(is.na(sampled_psf))) {
            bs[id] <- NA
            next
          }
          if (type == "Richness_hypothesis") {
            bs[id] <- mean(sampled_psf, na.rm = TRUE)
          } else {
            bs[id] <- sum(sampled_psf * weights, na.rm = TRUE)
          }
        }
        resampled[[type]][[j]] <- bs
      }
    }
    output[[Richness_con]] <- resampled
  }
  return(output)
}

# ------------------------------------------------------------------------------
# Function 2: Distribution of Observed PSF (richness levels)
# ------------------------------------------------------------------------------
Observed_PSF_new <- function(n_perm = n_iter) {
  resampled <- list()
  for (Richness_con in levels_rich) {
    if (Richness_con == "1") {
      tr_data <- df[df$remark == "1", "TG_Pot_res"] #**************** 140 pot
    } else {
      tr_data <- df[df$remark == Richness_con, "TG_Pot_res"] #****************144 pot
    }
    if (length(tr_data) == 0) {
      resampled[[Richness_con]] <- rep(NA, n_perm)
      next
    }
    bs <- replicate(n_perm, {
      tr_sample <- sample(tr_data, length(tr_data), replace = TRUE)
      mean(sapply(tr_sample, function(t) mean(log(t / CK_vals), na.rm = TRUE)))
    })
    resampled[[Richness_con]] <- bs
  }
  return(resampled)
}

set.seed(123)
null_dist_new <- Null_distribution_new()
Observed_dist_new <- Observed_PSF_new()

# ------------------------------------------------------------------------------
# Function 3: Hypothesis testing and summary function
# ------------------------------------------------------------------------------
NHST_summary_new <- function(null_data, Observed_data) {
  output <- list()
  for (Richness_con in levels_rich) {
    summary <- list()
    Observed_vec <- Observed_data[[Richness_con]]
    if (length(Observed_vec) == 0 || all(is.na(Observed_vec))) {
      summary[["Observed"]] <- c(NA, NA, NA, NA)
    } else {
      summary[["Observed"]] <- c(quantile(Observed_vec, 0.025, na.rm = TRUE),
                               mean(Observed_vec, na.rm = TRUE),
                               quantile(Observed_vec, 0.975, na.rm = TRUE), NA)
    }
    for (type in c("Richness_hypothesis", "Biomass_ratio_hypothesis")) {
      null_vec <- as.numeric(unlist(null_data[[Richness_con]][[type]]))
      if (length(null_vec) == 0 || all(is.na(null_vec)) || length(Observed_vec) == 0 || all(is.na(Observed_vec))) {
        summary[[type]] <- c(NA, NA, NA, NA)
      } else {
        n_samples <- min(length(Observed_vec), length(null_vec))
        Observed_sample <- Observed_vec[1:n_samples]
        null_sample <- null_vec[1:n_samples]
        diff <- Observed_sample - null_sample
        p <- 2 * min(mean(diff > 0, na.rm = TRUE), mean(diff < 0, na.rm = TRUE))
        summary[[type]] <- c(quantile(null_vec, 0.025, na.rm = TRUE),
                             mean(null_vec, na.rm = TRUE),
                             quantile(null_vec, 0.975, na.rm = TRUE), p)
      }
    }
    summary_df <- as.data.frame(do.call(rbind, summary))
    colnames(summary_df) <- c("2.5%", "mean", "97.5%", "p_value")
    summary_df$ES <- c("Observed", "Richness_hypothesis", "Biomass_ratio_hypothesis")
    output[[Richness_con]] <- summary_df
  }
  return(output)
}

nhst_new <- NHST_summary_new(null_dist_new, Observed_dist_new)


# ------------------------------------------------------------------------------
# Step3: Table S3: Summary results
# ------------------------------------------------------------------------------
ES_plot_data_new <- do.call(rbind, lapply(names(nhst_new), function(Richness_con) {
  df_ <- nhst_new[[Richness_con]]
  data.frame(Richness_con = Richness_con, Low = df_[, "2.5%"], Mean = df_[, "mean"], High = df_[, "97.5%"],
             Model = df_$ES, P_value = df_$p_value, stringsAsFactors = FALSE)
}))

ES_plot_data_new$Richness_con <- factor(ES_plot_data_new$Richness_con, levels = levels_rich)
ES_plot_data_new


# ------------------------------------------------------------------------------
# Step4: Figure 3(b) – Comparison of Observed and predicted PSF
# ------------------------------------------------------------------------------
b <- ggplot() +
  theme_bw() +
  coord_flip() + 
  scale_y_discrete(limits = levels_rich) +
  geom_estci(data = ES_plot_data_new[ES_plot_data_new$Model == "Observed", ], 
             aes(x = Mean, y = Richness_con, xmin = Low, xmax = High, xintercept = 0), 
             color = "#000000", size = 0.6, ci.linesize = 0.5, 
             position = position_nudge(y = 0)) +
  geom_estci(data = ES_plot_data_new[ES_plot_data_new$Model == "Richness_hypothesis", ], 
             aes(x = Mean, y = Richness_con, xmin = Low, xmax = High, xintercept = 0), 
             color = "#7B3D12", shape = 15, size = 0.6, ci.linesize = 0.5, 
             position = position_nudge(y = +0.15)) +
  geom_estci(data = ES_plot_data_new[ES_plot_data_new$Model == "Biomass_ratio_hypothesis", ], 
             aes(x = Mean, y = Richness_con, xmin = Low, xmax = High, xintercept = 0), 
             color = "#DE8125", shape = 17, size = 0.6, ci.linesize = 0.5, 
             position = position_nudge(y = +0.3)) +
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01), 
                     limits = c(-0.5, 0.03), expand = c(0, 0), breaks = seq(-0.5, 0.2, 0.1)) +
  mytheme + labs(x = "Community-level PSFs", y = "Species richness in conditioning community", tag = "(b)")
b



##################################################################################################################
#####                                                                                                        ##### 
#####                    Part 3---Calculate Observed PSF and two predicted PSFs for each pot                   #####
#####                                                                                                        #####  
################################################################################################################## 

# ----------------------------------------------------------------------- 
# Step 1: PSF for each single species (species-level mean)
# ----------------------------------------------------------------------- 
# Extract control values (Sterilized soil treatment)
CK_vals <- df[df$remark == "CK", "TG_Pot_res"]

species_psf <- sapply(stressors, function(sp) {
  tr_vals <- df[df$remark == sp, "TG_Pot_res"]
  if (length(tr_vals) == 0) return(NA)
  pot_psf <- sapply(tr_vals, function(t) mean(log(t / CK_vals), na.rm = TRUE))
  mean(pot_psf, na.rm = TRUE)
})
names(species_psf) <- stressors


# ----------------------------------------------------------------------- 
# Step 2:  Observed PSF (pot-level)
# ----------------------------------------------------------------------- 
pots <- df[df$remark %in% levels_rich, ]
pots$Observed_PSF_Pot <- sapply(pots$TG_Pot_res, function(t) mean(log(t / CK_vals), na.rm = TRUE))


# ----------------------------------------------------------------------- 
# Step 3:  Richness PSF & Biomass-Ratio PSF (pot-level)
# -----------------------------------------------------------------------
# Indices of species composition columns
comp_cols <- which(names(pots) %in% stressors)
pots$pred_RichnessPSF_Pot <- NA          # RichnessPSF
pots$pred_BiomassRatioPSF_Pot <- NA      # BiomassRatioPSF

for (i in 1:nrow(pots)) {
  sp_present <- stressors[pots[i, comp_cols] > 0] # Species present in the pot
  if (length(sp_present) == 0) next
  psfj <- species_psf[sp_present]
  if (all(is.na(psfj))) next
  # RichnessPSF: Simple average of species PSFs
  pots$pred_RichnessPSF_Pot[i] <- mean(psfj, na.rm = TRUE)
  # BiomassRatioPSF: Weighted average by initial biomass
  biomass <- unlist(pots[i, sp_present])
  if (sum(biomass) > 0) {
    weights <- biomass / sum(biomass)
    pots$pred_BiomassRatioPSF_Pot[i] <- sum(psfj * weights, na.rm = TRUE)
  } else {
    pots$pred_BiomassRatioPSF_Pot[i] <- mean(psfj, na.rm = TRUE)
  }
}


# ----------------------------------------------------------------------- 
# Step 4:  Save pot-level predictions (284 pots)
# -----------------------------------------------------------------------
pot_predictions <- pots[, c("Pot", "remark", "Richness_con", "TG_Pot_res", "Observed_PSF_Pot", "pred_RichnessPSF_Pot", "pred_BiomassRatioPSF_Pot")]
print(head(pot_predictions))

 

##################################################################################################################
#####                                                                                                        ##### 
#####               Part 4---Random forest analysis (using only Observed experimental pot data)                #####
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


# ------------------------------------------------------------------------------
# Step2: Figure 3(c) – Distribution of random forest R²
# ------------------------------------------------------------------------------
#plot_data_c <- results_rf %>% select(Model, R2_RF) %>% rename(R2 = R2_RF)
plot_data_c <- results_rf %>% dplyr::select(Model, R2_RF) %>% dplyr::rename(R2 = R2_RF)
model_labels <- c("Species_richness" = "Species richness",
                  "Species_identity" = "+ Species identity",
                  "Effect_size" = "+ Effect size")

c <- ggplot(data = plot_data_c, aes(x = Model, y = R2, fill = Model)) +
  geom_violin(color = "#00000000", alpha = 0.5, position = position_dodge(width = 0.3), trim = TRUE) +
  geom_pointrange(data = rf_summary, 
                  aes(y = Mean, ymax = CI.high, ymin = CI.low, color = Model),
                  position = position_dodge(width = 0.2)) +
  scale_fill_manual(values = c("#CAC5C6", "#DCD6E5", "#f0d99d")) +
  scale_color_manual(values = c("#000000", "#4E307D", "#7B3D12")) +
  theme_bw() + mytheme +
  ylim(c(0, 0.65)) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(x = " ", y = "Variability explained (R²)", tag = "(c)") +
  scale_x_discrete(labels = model_labels) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1))
c


# ------------------------------------------------------------------------------
#  Step3: Final models (for variable importance assessment)
# ------------------------------------------------------------------------------
set.seed(123)
final_model_richness <- randomForest(fml_richness, data = df.rf, ntree = n_tree, mtry = min(2, n_pred_richness))
final_model_identity <- randomForest(fml_identity, data = df.rf, ntree = n_tree, mtry = min(3, n_pred_identity))
final_model_Effect_size <- randomForest(fml_Effect_size, data = df.rf, ntree = n_tree, mtry = min(5, n_pred_Effect_size))

# Calculate variable importance for the full model
var_importance <- importance(final_model_Effect_size)
var_importance_df <- data.frame(
  Variable = rownames(var_importance),
  Importance = var_importance[, 1],
  stringsAsFactors = FALSE
) %>% arrange(desc(Importance))
print("Variable importance in the full model:")
print(var_importance_df)


# Store results in a list (for later use)
results_list <- list(
  summary_stats = rf_summary,
  variable_importance = var_importance_df,
  all_iterations = results_rf,
  panel_c_plot = c,
  final_models = list(
    Species_richness = final_model_richness,
    Species_identity = final_model_identity,
    Effect_size = final_model_Effect_size
  )
)
results_list


# ------------------------------------------------------------------------------
#  Step4: Combine the three subplots and save
# ------------------------------------------------------------------------------
P <- a / (b + c + plot_layout(widths = c(2, 1))) + plot_layout(heights = c(1.2, 1.1)); P
ggsave("Fig.3_modified.pdf", plot = P, width = 12, height = 11, units = "cm", dpi = 600)






## ================================================================================================================================== ##
#                                                                                                                                      #
#                                                     2：Fig.3 & Table S4                                                              #
#                                                                                                                                      #
## ================================================================================================================================== ##

##################################################################################################################
#####                                                                                                        ##### 
#####              Part 1--- Calculate Observed PSF and predicted PSF (Pot-level & Species-level)              #####
#####                                                                                                        #####  
##################################################################################################################

#--------------------------------------
### Step 1: data 
#--------------------------------------

df1 <- read_excel("data-Fig2-3 & TableS1-6.xlsx", sheet = "Fig2_data")

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
  
  # 5) Observed PSF  (also treat 0 as NA)
  pots$Observed_PSF <- sapply(pots[, resp], function(t) compute_PSF(t, CK_vals))
  
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
    Observed = pots$Observed_PSF,
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
    values_from = c(Observed, pred_RichnessPSF, BiomassRatioPSF),
    names_sep = "_"
  )

# Adjust column order: first Pot, Remark, Richness, then Observed, pred_RichnessPSF, BiomassRatioPSF for each response variable
resp_names <- responses_vec
ordered_cols <- c("Pot", "Remark", "Richness")
for (r in resp_names) {
  ordered_cols <- c(ordered_cols, paste0("Observed_", r), 
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
data_regress <- read_excel("data-Fig2-3 & TableS1-6.xlsx", sheet = "data_regress")   
glimpse(data_regress)

data_regress <- data_regress %>%
  mutate(
    ObservedPSF = as.numeric(ObservedPSF),
    RichnessPSF = as.numeric(RichnessPSF),
    BiomassRatioPSF = as.numeric(BiomassRatioPSF),
    Phylo_Dist = as.numeric(Phylo_Dist),
    Richness_con = as.factor(Richness_con),
    Pot = as.factor(Pot),
    Species = as.factor(Species),
    Species_full = as.factor(Species_full),
    Species_con  =  as.factor(Species_con)
  ) %>%
  filter(complete.cases(ObservedPSF, RichnessPSF, BiomassRatioPSF, Richness_con, Species, Pot))

# Define species level order and labels
sp_levels <- c("Pc", "Aa", "Sv", "St", "Ah", "Soc")
sp_labels <- c(Pc = "Paspalum_conjugatum", Aa = "Achyranthes_aspera",Sv = "Setaria_viridis",
               St = "Senna_tora", Ah = "Amaranthus_hybridus", Soc = "Senna_occidentalis")
data_regress$Species <- factor(data_regress$Species, levels = sp_levels)
data_regress$Species_full <- factor(data_regress$Species_full, levels = sp_labels)


# ------------------------------------------------------------------------------
#  Step 2: Model 1: RichnessPSF hypothesis --- Table S4
# ------------------------------------------------------------------------------
# Full model
mR_full <- lmer(ObservedPSF ~ Species * Richness_con * RichnessPSF + (1 | Pot),data = data_regress, REML = FALSE)
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
dat_R <- data_regress %>% filter(!is.na(ObservedPSF), !is.na(RichnessPSF))

p_R <- ggplot() +
  #geom_hline(yintercept = 0, linetype = 2, colour = "grey65") +
  geom_point(data = dat_R, aes(x = RichnessPSF, y = ObservedPSF),alpha = 0.3, size = 1.5, colour = "#999999") + 
  geom_line(data = plot_R, aes(x = RichnessPSF, y = emmean), colour = "#737373", linewidth = 1.2) +
  theme_classic() + theme_bw() + mytheme +
  labs(x = "Species-specific predicted PSF\n(Richness hypothesis)",
       y = "Observed species-specific PSF", tag = "(a)")
print(p_R)


# ------------------------------------------------------------------------------
#  Step 3: Model 2: BiomassRatioPSF hypothesis --- Table S4
# ------------------------------------------------------------------------------
mB_full <- lmer(ObservedPSF ~ Species * Richness_con * BiomassRatioPSF + (1 | Pot), data = data_regress, REML = FALSE)
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
dat_B <- data_regress %>% filter(!is.na(ObservedPSF), !is.na(BiomassRatioPSF)); dat_B

p_B <- ggplot() +
  geom_point(data = dat_B, aes(x = BiomassRatioPSF, y = ObservedPSF,colour = Richness_con), alpha = 0.4, size = 1.5,) +
  geom_line(data = plot_B, aes(x = BiomassRatioPSF, y = emmean,colour = Richness_con,linetype = sig),linewidth = 1.2) +
  scale_linetype_manual(values = c("significant" = "solid", "not significant" = "dashed")) +
  scale_color_manual(values = c("#3f3e93", "#7ca6e5", "#fedb81",  "#ff783e", "#c04037"),name = "Conditioning species richness") +
  theme_classic() + theme_bw() + mytheme  +
  labs(x = "Species-specific predicted PSF\n(biomass-ratio hypothesis)",
       y = "Observed species-specific PSF", tag = "(b)") + theme(legend.position = "right")
print(p_B)
combined <- p_R | p_B ; combined 








## ================================================================================================================================== ##
#                                                                                                                                      #
#                                                         3：Table S5-S6                                                               #
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
Condition_com = read.xlsx("data-Fig2-3 & TableS1-6.xlsx", sheet = "Species_AG_con", colNames = T, rowNames = T)
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
data_regress <- read_excel("data-Fig2-3 & TableS1-6.xlsx", sheet = "data_regress")   
glimpse(data_regress)

data_regress <- data_regress %>%
  mutate(
    ObservedPSF = as.numeric(ObservedPSF),
    RichnessPSF = as.numeric(RichnessPSF),
    BiomassRatioPSF = as.numeric(BiomassRatioPSF),
    Phylo_Dist = as.numeric(Phylo_Dist),
    Richness_con = as.factor(Richness_con),
    Pot = as.factor(Pot),
    Species = as.factor(Species),
    Species_full = as.factor(Species_full),
    Species_con  = as.factor(Species_con)                 
  ) %>%
  filter(complete.cases(ObservedPSF, RichnessPSF, BiomassRatioPSF, Richness_con, Species, Pot))


#--------------------------------------
### Step 2: Model 1  --- Table S5
#-------------------------------------- 
data_regress 

mP_full <- lmer(ObservedPSF ~ Species * Richness_con * Phylo_Dist + (1 | Pot), data = data_regress, REML = FALSE)
anova_result <- anova(mP_full)
p_table <- anova_result; print(p_table)
p_values <- anova_result$`Pr(>F)`; print(p_values) 
p.adjust(p_values, "BH")  


#--------------------------------------
### Step 3: Model 2 --- Table S5
#-------------------------------------- 
### general phylogenetic effect among heterospecific-only pots
data_regress_heter <- data_regress %>% filter(Species_con == "0")

mP_full <- lmer(ObservedPSF ~ Species * Richness_con * Phylo_Dist + (1 | Pot), data = data_regress_heter, REML = FALSE)
anova_result <- anova(mP_full)
p_table <- anova_result; print(p_table)
p_values <- anova_result$`Pr(>F)`; print(p_values) 
p.adjust(p_values, "BH") 


#--------------------------------------
### Step 4: Model 3 --- Table S6
#--------------------------------------
mC_full <- lmer(ObservedPSF ~  Species * Richness_con * Species_con + (1 | Pot), data = data_regress, REML = FALSE)
anova_result <- anova(mC_full)
p_table <- anova_result; print(p_table)
p_values <- anova_result$`Pr(>F)`; print(p_values) 
p.adjust(p_values, "BH")




