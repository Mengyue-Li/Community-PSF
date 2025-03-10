
###########################################################################################
#####                                  Library & mytheme                              #####
###########################################################################################
#chooseCRANmirror()
#devtools::install_github('brendanf/FUNGuildR')
library(readxl)
library(openxlsx)
library(multcomp)
library(tidyverse)
library(car)
library(lme4)
library(vegan)
library(FUNGuildR)
library(dplyr)
library(lmerTest)
library(ggbeeswarm)
library(reshape2)
library(DHARMa)
library(glmmTMB)
library(emmeans)
library(broom.mixed)
library(rdacca.hp)
library(glmm.hp)
library(MuMIn)
library(ggplot2)
library(patchwork)

mytheme= theme(
  legend.position = "none",
  panel.grid=element_blank(), 
  legend.title = element_text(colour='black', size=11), 
  legend.text = element_text(size = 10),  
  legend.background = element_rect(fill = NA),  
  axis.ticks = element_line(color='black'),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(colour='black', size=16),
  axis.title.y = element_text(colour='black', size=16),
  axis.text = element_text(colour='black',size=14),
  plot.tag = element_text(size = 18, face = "bold"))

###########################################################################################
#####                   Part1--- Index calculation & Composition(RDA)                 #####
###########################################################################################

#--------------------------------------------------------------------
# Step 1:  loading data and identify taxonomy levels
#--------------------------------------------------------------------

# raw data
setwd("C:/Users/MY/Desktop/CODE/Fig.2")
  
MetaTab <- read_excel("data_Fig.2.xlsx",sheet=1);
MetaTab <- as.data.frame(MetaTab);rownames(MetaTab) <- MetaTab[,1] ;MetaTab <- MetaTab[MetaTab$Monoculture == "NO",]
Biomass_dat <- read_excel("data_Fig.2.xlsx",sheet=2);Biomass_dat <- as.data.frame(Biomass_dat);rownames(Biomass_dat) <- Biomass_dat[,1] 
seqDat_full <- read_excel("data_Fig.2.xlsx",sheet=3); seqDat_full <-seqDat_full[,1:282]; dim(seqDat_full);seqDat_full <- as.data.frame(seqDat_full);rownames(seqDat_full) <- seqDat_full[,1]
tax <- seqDat_full$Taxonomy; tax <- as.data.frame(tax); rownames(tax) <- rownames(seqDat_full)
seqDat <- seqDat_full[,c(2:281)]

# check the read counts
seqDat <- seqDat[,1:280]
checkReadCounts <- colSums(seqDat)
print(checkReadCounts)
sum(checkReadCounts<100)
range(checkReadCounts[checkReadCounts>50])

# Function to split taxonomic strings
f.split.utax.taxo <- function(x) {
  x <- gsub("\\(.*?\\)", "", x)
  out <- unlist(strsplit(x, ";", TRUE))
  
  if (length(out) > 1) {
    out <- t(sapply(out, function(x) {
      # Properly split each level and ensure the length is always 2
      parts <- unlist(strsplit(x, "__", fixed = TRUE))
      if (length(parts) == 2) {
        return(parts)
      } else {
        return(c(parts[1], "unknown"))
      }
    }))
  } else {
    # Return empty matrix if no valid split
    out <- matrix(ncol = 2, nrow = 0)
  }
  
  return(out)
}

# Function to modify the taxonomic table
f.modify.utax.taxon.table <- function(taxDat, onlyConfident = TRUE) {
  # Ensure we're working with the correct column in the data frame
  colToUse <- 1  
  emptyEntry <- rep("ukn", nrow(taxDat))
  out <- data.frame(
    k = emptyEntry,
    p = emptyEntry,
    c = emptyEntry,
    o = emptyEntry,
    f = emptyEntry,
    g = emptyEntry,
    s = emptyEntry,
    row.names = rownames(taxDat),
    stringsAsFactors = FALSE
  )
  
  for (rn in rownames(taxDat)) {
    temp <- f.split.utax.taxo(taxDat[rn, colToUse])
    
    # Ensure that temp is a matrix and has at least one row with 2 columns
    if (is.matrix(temp) && nrow(temp) > 0 && ncol(temp) == 2) {
      tax_levels <- temp[, 1]
      tax_values <- temp[, 2]
      
      # Filter tax_levels that match valid column names in the output table
      valid_indices <- tax_levels %in% colnames(out)
      
      # Assign values only for valid taxonomic levels
      out[rn, tax_levels[valid_indices]] <- tax_values[valid_indices]
    }
  }
  
  colnames(out) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  # Handle cases where genus is known but higher levels are not
  for (i in ncol(out):2) {
    lowerIsKnown <- out[, i] != "ukn"
    upperIsUnknown <- out[, i - 1] == "ukn"
    out[lowerIsKnown & upperIsUnknown, i - 1] <- paste0("undef_", out[lowerIsKnown & upperIsUnknown, i])
  }
  
  return(out)
}

f.geomeanForDESeq2 <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

f.print.message <- function(x) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste0("=== ", x,"\n")) }

# Applying the function
taxMod <- f.modify.utax.taxon.table(tax, onlyConfident = TRUE)
## taxMod_with_guilds <- funguild_assign(seqDat_full, db = get_funguild_db(), tax_col = "Taxonomy");rownames(taxMod_with_guilds) <- rownames(seqDat) #FUNGuild database:2025.01.27
## If FUNGuild database is updated, please use the following line of code
taxMod_with_guilds <-read_excel("data_Fig.2.xlsx",sheet=3);taxMod_with_guilds<- as.data.frame(taxMod_with_guilds); rownames(taxMod_with_guilds) <- taxMod_with_guilds[,1] 
taxMod_all <- merge(taxMod_with_guilds,taxMod,by=0, all=TRUE);rownames(taxMod_all) <- taxMod_all[,1]

# extract pathogen and AMF 
seqDat_plantPathogen <- taxMod_all[taxMod_all$guild == "Plant Pathogen"| taxMod_all$genus == "Fusarium", c(3:282)]; seqDat_plantPathogen <- na.omit(seqDat_plantPathogen)
seqDat_AMF <- taxMod_all[taxMod_all$guild == "Arbuscular Mycorrhizal",c(3:282)]; seqDat_AMF <- na.omit(seqDat_AMF)


#--------------------------------------------------------------------
# Step 2:  Shannon diversity and inverse Simpson diversity 
#--------------------------------------------------------------------
set.seed(123)
seqDat_Flattening_overall <- as.data.frame(t(rrarefy(t(seqDat), min(colSums(seqDat)))));colSums(seqDat_Flattening_overall)
Shannon_overall <- diversity(seqDat_Flattening_overall, index = "shannon", MARGIN = 2, base = exp(1))
index_overall <- as.data.frame(Shannon_overall)
Invsimp_overall <- diversity(seqDat_Flattening_overall, index = "inv", MARGIN = 2)
index_overall$Invsimp_overall <- Invsimp_overall

set.seed(123)
seqDat_Flattening_pathogen <- as.data.frame(t(rrarefy(t(seqDat_plantPathogen), min(colSums(seqDat_plantPathogen)))));colSums(seqDat_Flattening_pathogen)
Shannon_pathogen <- diversity(seqDat_Flattening_pathogen, index = "shannon", MARGIN = 2, base = exp(1))
index_pathogen <- as.data.frame(Shannon_pathogen)
Invsimp_pathogen <- diversity(seqDat_Flattening_pathogen, index = "inv", MARGIN = 2)
index_pathogen$Invsimp_pathogen <- Invsimp_pathogen

set.seed(123)
seqDat_Flattening_AMF <- as.data.frame(t(rrarefy(t(seqDat_AMF), min(colSums(seqDat_AMF)))));colSums(seqDat_Flattening_AMF)
Shannon_AMF <- diversity(seqDat_Flattening_AMF, index = "shannon", MARGIN = 2, base = exp(1))
index_AMF <- as.data.frame(Shannon_AMF)
Invsimp_AMF <- diversity(seqDat_Flattening_AMF, index = "inv", MARGIN = 2)
index_AMF$Invsimp_AMF <- Invsimp_AMF

# merge data
div_dat_overall <- merge(index_overall,MetaTab,by=0) %>% left_join(Biomass_dat)
div_dat_overall$Treat <- paste(div_dat_overall$FS,div_dat_overall$NR,sep = "_");head(div_dat_overall)

div_dat_pathogen <- merge(index_pathogen,MetaTab,by=0) %>% left_join(Biomass_dat)
div_dat_pathogen$Treat <- paste(div_dat_pathogen$FS,div_dat_pathogen$NR,sep = "_");head(div_dat_pathogen)

div_dat_AMF <- merge(index_AMF,MetaTab,by=0) %>% left_join(Biomass_dat)
div_dat_AMF$Treat <- paste(div_dat_AMF$FS,div_dat_AMF$NR,sep = "_");head(div_dat_AMF)

### Shannon & Invsimp
div_dat_overall <- div_dat_overall[div_dat_overall$Monoculture == "NO",]
# Shannon
model <- glmmTMB(Shannon_overall ~ FS*NR*COT, data = div_dat_overall)
qqnorm(resid(model));qqline(resid(model))
Anova(model, type = "III")-> model_result; model_result
p1 <- model_result$Pr
# Invsimp
model <- glmmTMB(Invsimp_overall ~ FS*NR*COT, data = div_dat_overall) 
qqnorm(resid(model));qqline(resid(model))
Anova(model, type = "III")-> model_result; model_result
p2 <- model_result$Pr

div_dat_pathogen <- div_dat_pathogen[div_dat_pathogen$Monoculture == "NO",]
# Shannon
model <- glmmTMB(Shannon_pathogen ~ FS*NR*COT, data = div_dat_pathogen)
qqnorm(resid(model));qqline(resid(model))
Anova(model, type = "III")-> model_result; model_result
p3 <- model_result$Pr
# Invsimp
model <- glmmTMB(Invsimp_pathogen ~ FS*NR*COT, data = div_dat_pathogen) 
qqnorm(resid(model));qqline(resid(model))
Anova(model, type = "III")-> model_result; model_result
p4 <- model_result$Pr

div_dat_AMF <- div_dat_AMF[div_dat_AMF$Monoculture == "NO",]
# Shannon
model <- glmmTMB(Shannon_AMF ~ FS*NR*COT, data = div_dat_AMF)
qqnorm(resid(model));qqline(resid(model))
Anova(model, type = "III")-> model_result; model_result
p5 <- model_result$Pr
# Invsimp
model <- glmmTMB(Invsimp_AMF ~ FS*NR*COT, data = div_dat_AMF) 
qqnorm(resid(model));qqline(resid(model))
Anova(model, type = "III")-> model_result; model_result
p6 <- model_result$Pr


#--------------------------------------
### Table S8
#--------------------------------------
p <- na.omit(c(p1,p2,p3,p4,p5,p6))
p.adjust(p, "BH")
p.adjust(p, "bonferroni")
 

### since the variable co-occuring species type (COT) is not significant, we indluding it into random faCOTor
#--------------------------------------
### Table S3  
#--------------------------------------
#----------------------Table S3-1: Shannon&Invsimp
# Shannon 
model <- glmmTMB(Shannon_overall ~  FS*NR + (1|COT), data = div_dat_overall)
qqnorm(resid(model));qqline(resid(model))
Anova(model, type = "III")-> model_result; model_result
pp1 <- model_result$Pr
# Invsimp
model <- glmmTMB(Invsimp_overall ~ FS*NR + (1|COT), data = div_dat_overall)
qqnorm(resid(model));qqline(resid(model))
Anova(model, type = "III")-> model_result; model_result
pp2 <- model_result$Pr

# Shannon 
model <- glmmTMB(Shannon_pathogen ~ FS*NR + (1|COT), data = div_dat_pathogen) 
qqnorm(resid(model))
Anova(model, type = "III")-> model_result; model_result
pp3 <- model_result$Pr
# Invsimp  
model <- glmmTMB(Invsimp_pathogen ~  FS*NR + (1|COT), data = div_dat_pathogen)
qqnorm(resid(model))
Anova(model, type = "III")-> model_result; model_result
pp4 <- model_result$Pr

# Shannon
model <- glmmTMB(Shannon_AMF ~ FS*NR + (1|COT), data = div_dat_AMF) 
qqnorm(resid(model))
Anova(model, type = "III")-> model_result; model_result
pp5 <- model_result$Pr
# Invsimp 
model <- glmmTMB(Invsimp_AMF ~  FS*NR + (1|COT), data = div_dat_AMF)
qqnorm(resid(model))
Anova(model, type = "III")-> model_result; model_result
pp6 <- model_result$Pr

pp <- na.omit(c(pp1,pp2,pp3,pp4,pp5,pp6))
p.adjust(pp, "BH")
p.adjust(pp, "bonferroni")


#--------------------------------------------------------------------
# Step 3:  RDA result
#--------------------------------------------------------------------
#----------------------Table S3-2: RDA 
### RDA_result_overall
normData <- read.xlsx("data_Fig.2.xlsx", sheet =4, colNames = T, rowNames = T)
normData_overall <- normData[,rownames(MetaTab)] # removed monoculture

rownames(MetaTab) %in% colnames(normData_overall)## check

selectedOtus <- rownames(normData_overall)
selectedSamples <- rownames(MetaTab)   
subSampleTab <- MetaTab[selectedSamples,]
subSampleTab$Treat <- paste(MetaTab$FS,MetaTab$NR,sep = "_");subSampleTab

RDA2 <- rda(t(normData_overall) ~ FS+ NR + FS:NR +Condition(COT), data = subSampleTab)
RDA_result_overall = anova(RDA2, by = "term", permutations = 999); RDA_result_overall


### RDA_result_pathogen 
normData_pathogen <- normData[,rownames(MetaTab)];normData_pathogen # removed monoculture
normData_pathogen = normData_pathogen[rownames(seqDat_plantPathogen), ]
normData_pathogen = na.omit(normData_pathogen);normData_pathogen

rownames(MetaTab) %in% colnames(normData_pathogen)## check

selectedOtus_pathogen <- rownames(normData_pathogen)
selectedSamples_pathogen <- rownames(MetaTab)
subSampleTab_pathogen <- MetaTab[selectedSamples_pathogen,];subSampleTab_pathogen
subSampleTab_pathogen$Treat <- paste(subSampleTab_pathogen$FS,subSampleTab_pathogen$NR,sep = "_");subSampleTab_pathogen

RDA2 <- rda(t(normData_pathogen) ~ FS+ NR + FS:NR +Condition(COT), data = subSampleTab)
RDA_result_pathogen = anova(RDA2, by = "term", permutations = 999); RDA_result_pathogen

### RDA_result_AMF
normData_AMF <- normData[,rownames(MetaTab)];normData_AMF # removed monoculture
normData_AMF = normData_AMF[rownames(seqDat_AMF), ]
normData_AMF = na.omit(normData_AMF);normData_AMF

rownames(MetaTab) %in% colnames(normData_AMF)## check

selectedOtus_AMF <- rownames(normData_AMF)
selectedSamples_AMF <- rownames(MetaTab)
subSampleTab_AMF <- MetaTab[selectedSamples_AMF,];subSampleTab_AMF
subSampleTab_AMF$Treat <- paste(subSampleTab_AMF$FS,subSampleTab_AMF$NR,sep = "_");subSampleTab_AMF

RDA2 <- rda(t(normData_AMF) ~ FS+ NR + FS:NR +Condition(COT), data = subSampleTab)
RDA_result_AMF = anova(RDA2, by = "term", permutations = 999); RDA_result_AMF


#--------------------------------------
### TableS4
#--------------------------------------
div_dat_long_overall <- melt(div_dat_overall, id= c("Row.names","Pot","FS","NR","COT","CS","Treat","TG_pot","F_survival","Monoculture"));head(div_dat_long_overall)
div_dat_long_pathogen <- melt(div_dat_pathogen, id= c("Row.names","Pot","FS","NR","COT","CS","Treat","TG_pot","F_survival","Monoculture"));head(div_dat_long_overall)
div_dat_long_AMF <- melt(div_dat_AMF, id= c("Row.names","Pot","FS","NR","COT","CS","Treat","TG_pot","F_survival","Monoculture"));head(div_dat_long_overall)

# check the model first
par(mfrow = c(1, 2))
model <- glmmTMB(value ~ NR * FS,family = gaussian,data = filter(div_dat_long_overall, variable == "Shannon_overall"));qqnorm(resid(model));qqline(resid(model))
model <- glmmTMB(value ~ NR * FS,family = gaussian,data = filter(div_dat_long_overall, variable == "Invsimp_overall"));qqnorm(resid(model));qqline(resid(model))

model <- glmmTMB(value ~ NR * FS,family = gaussian,data = filter(div_dat_long_pathogen, variable == "Shannon_pathogen"));qqnorm(resid(model));qqline(resid(model))
model <- glmmTMB(value ~ NR * FS,family = gaussian,data = filter(div_dat_long_pathogen, variable == "Invsimp_pathogen"));qqnorm(resid(model));qqline(resid(model))

model <- glmmTMB(value ~ NR * FS,family = gaussian,data = filter(div_dat_long_AMF, variable == "Shannon_AMF"));qqnorm(resid(model));qqline(resid(model))
model <- glmmTMB(value ~ NR * FS,family = gaussian,data = filter(div_dat_long_AMF, variable == "Invsimp_AMF"));qqnorm(resid(model));qqline(resid(model))

# results of overall  pathogen and AMF
results_overall <- list()
# Loop over each variable
for (var in unique(div_dat_long_overall$variable)) {
  
  # Filter data for the current variable
  data_var <- div_dat_long_overall %>% filter(variable == var)
  
  # Fit the appropriate model using glmmTMB
  if (var %in% c("Shannon_overall", "Invsimp_overall")) {
    model <- glmmTMB(value ~ NR * FS + (1 | COT), family = gaussian, data = data_var,
                     control = glmmTMBControl(optimizer = optim,optArgs = list(method = "BFGS"), profile = TRUE))
  }
  
  # Check for singular fit (glmmTMB models do not have isSingular() directly)
  if (any(is.na(model$fit$par))) {
    warning(paste("Model for variable", var, "failed to converge. Results may not be reliable."))
  }
  # Calculate slopes using emmeans
  emm_slopes <- emtrends(model, specs = ~ FS, var = "NR")
  
  # Extract slope data
  slopes_df <- as.data.frame(summary(emm_slopes))
  
  # Extract confidence intervals dynamically
  ci_df <- confint(emm_slopes, level = 0.95) %>% as.data.frame()
  
  # Standardize CI column names to `LcL` and `UcL`
#  ci_df <- confint(emm_slopes) %>%
#    rename(LcL = lower.CL, UcL = upper.CL) %>%
#    select(FS, LcL, UcL)

  # Combine slope data with standardized confidence intervals
  slopes_df <- slopes_df %>%
    mutate(
      t.ratio = NR.trend / SE, 
      p.value = 2 * pt(-abs(t.ratio), df), 
      std_slope = NR.trend / sd(data_var$NR, na.rm = TRUE),
      variable = var
    )
  
  # Store results
  results_overall[[var]] <- slopes_df
}


results_pathogen <- list()
# Loop over each variable
for (var in unique(div_dat_long_pathogen$variable)) {
  
  # Filter data for the current variable
  data_var <- div_dat_long_pathogen %>% filter(variable == var)
  
  # Fit the appropriate model using glmmTMB
  if (var %in% c("Shannon_pathogen", "Invsimp_pathogen")) {
    model <- glmmTMB(value ~ NR * FS + (1 | COT), family = gaussian, data = data_var,
                     control = glmmTMBControl(optimizer = optim,optArgs = list(method = "BFGS"), profile = TRUE))
  }
  # Check for singular fit (glmmTMB models do not have isSingular() directly)
  if (any(is.na(model$fit$par))) {
    warning(paste("Model for variable", var, "failed to converge. Results may not be reliable."))
  }
  # Calculate slopes using emmeans
  emm_slopes <- emtrends(model, specs = ~ FS, var = "NR")
  
  # Extract slope data
  slopes_df <- as.data.frame(summary(emm_slopes))
  
  # Extract confidence intervals dynamically
  ci_df <- confint(emm_slopes, level = 0.95) %>% as.data.frame()
  
  # Standardize CI column names to `LcL` and `UcL`
#  ci_df <- confint(emm_slopes) %>%
#    rename(LcL = lower.CL, UcL = upper.CL) %>%
#    select(FS, LcL, UcL)

  # Combine slope data with standardized confidence intervals
  slopes_df <- slopes_df %>%
    mutate(
      t.ratio = NR.trend / SE, 
      p.value = 2 * pt(-abs(t.ratio), df), 
      std_slope = NR.trend / sd(data_var$NR, na.rm = TRUE),
      variable = var
    )
  
  # Store results
  results_pathogen[[var]] <- slopes_df
}
 

results_AMF <- list()
# Loop over each variable
for (var in unique(div_dat_long_AMF$variable)) {
  
  # Filter data for the current variable
  data_var <- div_dat_long_AMF %>% filter(variable == var)
  
  # Fit the appropriate model using glmmTMB
  if (var %in% c("Shannon_AMF", "Invsimp_AMF")) {
    model <- glmmTMB(value ~ NR * FS + (1 | COT), family = gaussian, data = data_var,
                     control = glmmTMBControl(optimizer = optim,optArgs = list(method = "BFGS"), profile = TRUE))
  }
  # Check for singular fit (glmmTMB models do not have isSingular() directly)
  if (any(is.na(model$fit$par))) {
    warning(paste("Model for variable", var, "failed to converge. Results may not be reliable."))
  }
  
  # Calculate slopes using emmeans
  emm_slopes <- emtrends(model, specs = ~ FS, var = "NR")
  
  # Extract slope data
  slopes_df <- as.data.frame(summary(emm_slopes))
  
  # Extract confidence intervals dynamically
  ci_df <- confint(emm_slopes, level = 0.95) %>% as.data.frame()
  
  # Standardize CI column names to `LcL` and `UcL`
#  ci_df <- confint(emm_slopes) %>%
#    rename(LcL = lower.CL, UcL = upper.CL) %>%
#    select(FS, LcL, UcL)
  
  # Combine slope data with standardized confidence intervals
  slopes_df <- slopes_df %>%
    mutate(
      t.ratio = NR.trend / SE, 
      p.value = 2 * pt(-abs(t.ratio), df), 
      std_slope = NR.trend / sd(data_var$NR, na.rm = TRUE),
      variable = var
    )
  
  # Store results
  results_AMF[[var]] <- slopes_df
}
results_overall <- bind_rows(results_overall);head(results_overall)
results_pathogen <- bind_rows(results_pathogen);head(results_pathogen)
results_AMF <- bind_rows(results_AMF);head(results_AMF)


#-------------------------------------- 
# index_Visualization: Fig_3A,B,E,F
#--------------------------------------
### Fig_3A B
div_dat_long_overall <- div_dat_long_overall %>% left_join(results_overall %>% dplyr::select(FS, variable, NR.trend, SE, p.value),
                                                           by = c("FS", "variable")) %>%mutate(significant = p.value < 0.05) # Mark significant slopes
div_dat_long_pathogen <- div_dat_long_pathogen %>% left_join(results_pathogen %>% dplyr::select(FS, variable, NR.trend, SE, p.value),
                                                             by = c("FS", "variable")) %>% mutate(significant = p.value < 0.05)
div_dat_long_AMF <- div_dat_long_AMF %>% left_join(results_AMF %>% dplyr::select(FS, variable, NR.trend, SE, p.value),
                                                   by = c("FS", "variable")) %>% mutate(significant = p.value < 0.05)
# Define the critical t-value for 95% confidence intervals
t_critical_pathogen <- qt(0.975, df = max(results_pathogen$df)) # Use the maximum degrees of freedom from the results

# Prepare regression lines with confidence intervals for all slopes
regression_lines_pathogen <- div_dat_long_pathogen %>% group_by(FS, variable) %>% summarize(slope = unique(NR.trend),
                                                                                  intercept = mean(value,na.rm=T) - mean(NR) * unique(NR.trend),SE = unique(SE),
                                                                                  significant = unique(significant),x_min = min(NR),x_max = max(NR),.groups = "drop") %>%
  rowwise() %>% # Ensure row-by-row calculations
  mutate(y_min = intercept + slope * x_min,y_max = intercept + slope * x_max,
         y_min_upper = if (significant) intercept + slope * x_min + t_critical_pathogen * SE else NA,
         y_min_lower = if (significant) intercept + slope * x_min - t_critical_pathogen * SE else NA,
         y_max_upper = if (significant) intercept + slope * x_max + t_critical_pathogen * SE else NA,
         y_max_lower = if (significant) intercept + slope * x_max - t_critical_pathogen * SE else NA ) %>% ungroup() # Remove rowwise grouping for further processing

# Prepare data for significant regression lines in long format
significant_lines_pathogen <- regression_lines_pathogen %>%
  filter(significant == TRUE) %>%
  dplyr::select(FS, variable, x_min, x_max, y_min, y_max) %>%
  pivot_longer(cols = c(x_min, x_max), names_to = "x_type", values_to = "x") %>%
  pivot_longer(cols = c(y_min, y_max), names_to = "y_type", values_to = "y") %>%
  dplyr::filter((x_type == "x_min" & y_type == "y_min") | (x_type == "x_max" & y_type == "y_max")) %>%
  dplyr::select(-x_type, -y_type)

# Prepare data for SE shading and regression lines
line_and_ribbon_data_pathogen <- regression_lines_pathogen %>%dplyr::filter(significant == TRUE) %>%
  dplyr::select(FS, variable, x_min, x_max, y_min, y_max, y_min_lower, y_min_upper, y_max_lower, y_max_upper) %>%
  pivot_longer(cols = c(x_min, x_max), names_to = "x_type",values_to = "x") %>%
  mutate(ymin = ifelse(x_type == "x_min", y_min_lower, y_max_lower),
         ymax = ifelse(x_type == "x_min", y_min_upper, y_max_upper),y = ifelse(x_type == "x_min", y_min, y_max))

# Prepare data for significant regression lines with SE shading
regression_lines_long_pathogen <- regression_lines_pathogen %>%filter(significant == TRUE) %>% 
  mutate(x_min_ymin = intercept + slope * x_min,x_max_ymax = intercept + slope * x_max)

# Subset data for the plot (Richness and Invsimp)
data_richness_invsimp_pathogen <- div_dat_long_pathogen %>% filter(variable %in% c("Shannon_pathogen", "Invsimp_pathogen"))

# add the three monoculture points
monoculture_pathogen <- index_pathogen[c(1, 101, 181), ]

# Create the additional points data frame
additional_points_pathogen <- data.frame(
  Row.names = rownames(monoculture_pathogen),          # Row names from index_overall
  Pot = NA,                                   # Placeholder for Pot column
  FS = NA,                                    # Placeholder for FS column
  NR = 0,                                     # NR = 0 for Monoculture
  COT = NA,                                    # Placeholder for COT column
  CS = NA,                                    # Placeholder for CS column
  Treat = NA,                                 # Placeholder for Treat column
  TG_pot = NA,                                # Placeholder for TG_pot column
  F_survival = NA,                            # Placeholder for F_survival column
  Monoculture = "YES",                        # Mark as Monoculture
  variable = rep(c("Shannon_pathogen", "Invsimp_pathogen"), each = 3),
  value = c(monoculture_pathogen$Shannon_pathogen, monoculture_pathogen$Invsimp_pathogen),  # Add values dynamically
  NR.trend = NA,                              # Placeholder for NR.trend
  SE = NA,                                    # Placeholder for SE
  p.value = NA,                               # Placeholder for p-value
  significant = NA                            # Placeholder for significant flag
)

# Combine the additional points with the existing data
data_richness_invsimp_pathogen <- bind_rows(data_richness_invsimp_pathogen, additional_points_pathogen)
# order the variables 
data_richness_invsimp_pathogen$variable <- factor(data_richness_invsimp_pathogen$variable, levels = c("Shannon_pathogen", "Invsimp_pathogen"))
significant_lines_pathogen$variable <- factor(significant_lines_pathogen$variable, levels = c("Shannon_pathogen", "Invsimp_pathogen"))
line_and_ribbon_data_pathogen$variable <- factor(line_and_ribbon_data_pathogen$variable, levels = c("Shannon_pathogen", "Invsimp_pathogen"))
regression_lines_pathogen$variable <- factor(regression_lines_pathogen$variable, levels = c("Shannon_pathogen", "Invsimp_pathogen"))

pd = position_dodge(.2)
ggplot(data_richness_invsimp_pathogen, aes(x = NR, y = value, color = as.factor(FS))) +
  geom_point(size=2,position=pd,alpha=1,aes(shape = as.factor(NR), fill = FS)) + 
  scale_color_manual(values = c("#6C7C86","#83B8B0","#EAD69D")) +
  ggnewscale::new_scale_color() + ggnewscale::new_scale_fill() + 
  geom_line(data = significant_lines_pathogen %>% filter(variable %in% c("Shannon_pathogen", "Invsimp_pathogen")), aes(x = x, y = y, color = as.factor(FS)), linetype = "solid",alpha = 1) +
  geom_ribbon(data = line_and_ribbon_data_pathogen %>% filter(variable %in% c("Shannon_pathogen", "Invsimp_pathogen")), aes(x = x, ymin = ymin, ymax = ymax, fill = as.factor(FS)), alpha = 0.5, inherit.aes = FALSE) +
  geom_segment(data = regression_lines_pathogen %>% filter(!significant & variable %in% c("Shannon_pathogen", "Invsimp_pathogen")), aes(x = x_min, xend = x_max, y = y_min, yend = y_max, color = as.factor(FS)), linetype = "dashed", linewidth = 0.4,alpha = 0.7) +
  scale_x_continuous(breaks = c(0, 1, 2, 4), labels = c("Monoculture", "1", "2", "4")) +  
  scale_color_manual(values = c("#6C7C86","#83B8B0","#EAD69D")) +
  scale_fill_manual(values = c("#6C7C86","#83B8B0","#EAD69D")) + 
  scale_shape_manual(values = c(1,16,16,16)) +
  facet_wrap(vars(variable), scales = "free", nrow = 2, ncol = 1, labeller = labeller(.multi_line = FALSE)) +
  guides(color = guide_legend(override.aes = list(size = 2)), shape = guide_legend(override.aes = list(size = 2)), alpha = guide_legend(override.aes = list(size = 2))) +
  theme_bw() + mytheme + theme(legend.position = "none") + 
  labs(x = "Neighbour plant species richness", y = NULL, tag = "A") -> Fig_3AB; Fig_3AB


### Fig_3E F
# Define the critical t-value for 95% confidence intervals
t_critical_AMF<- qt(0.975, df = max(results_AMF$df)) # Use the maximum degrees of freedom from the results

# Prepare regression lines with confidence intervals for all slopes
regression_lines_AMF<- div_dat_long_AMF%>% group_by(FS, variable) %>% summarize(slope = unique(NR.trend),
                                                                                            intercept = mean(value,na.rm=T) - mean(NR) * unique(NR.trend),SE = unique(SE),
                                                                                            significant = unique(significant),x_min = min(NR),x_max = max(NR),.groups = "drop") %>%
  rowwise() %>% # Ensure row-by-row calculations
  mutate(y_min = intercept + slope * x_min,y_max = intercept + slope * x_max,
         y_min_upper = if (significant) intercept + slope * x_min + t_critical_AMF * SE else NA,
         y_min_lower = if (significant) intercept + slope * x_min - t_critical_AMF * SE else NA,
         y_max_upper = if (significant) intercept + slope * x_max + t_critical_AMF * SE else NA,
         y_max_lower = if (significant) intercept + slope * x_max - t_critical_AMF * SE else NA ) %>% ungroup() # Remove rowwise grouping for further processing

# Prepare data for significant regression lines in long format
significant_lines_AMF<- regression_lines_AMF%>%
  filter(significant == TRUE) %>%
  dplyr::select(FS, variable, x_min, x_max, y_min, y_max) %>%
  pivot_longer(cols = c(x_min, x_max), names_to = "x_type", values_to = "x") %>%
  pivot_longer(cols = c(y_min, y_max), names_to = "y_type", values_to = "y") %>%
  dplyr::filter((x_type == "x_min" & y_type == "y_min") | (x_type == "x_max" & y_type == "y_max")) %>%
  dplyr::select(-x_type, -y_type)

# Prepare data for SE shading and regression lines
line_and_ribbon_data_AMF<- regression_lines_AMF%>%dplyr::filter(significant == TRUE) %>%
  dplyr::select(FS, variable, x_min, x_max, y_min, y_max, y_min_lower, y_min_upper, y_max_lower, y_max_upper) %>%
  pivot_longer(cols = c(x_min, x_max), names_to = "x_type",values_to = "x") %>%
  mutate(ymin = ifelse(x_type == "x_min", y_min_lower, y_max_lower),
         ymax = ifelse(x_type == "x_min", y_min_upper, y_max_upper),y = ifelse(x_type == "x_min", y_min, y_max))

# Prepare data for significant regression lines with SE shading
regression_lines_long_AMF<- regression_lines_AMF%>%filter(significant == TRUE) %>% 
  mutate(x_min_ymin = intercept + slope * x_min,x_max_ymax = intercept + slope * x_max)

# Subset data for the plot (Richness and Invsimp)
data_richness_invsimp_AMF<- div_dat_long_AMF%>% filter(variable %in% c("Shannon_AMF", "Invsimp_AMF"))
  
# add the three monoculture points
monoculture_AMF<- index_AMF[c(1, 101, 181), ]

# Create the additional points data frame
additional_points_AMF<- data.frame(
  Row.names = rownames(monoculture_AMF),          # Row names from index_overall
  Pot = NA,                                   # Placeholder for Pot column
  FS = NA,                                    # Placeholder for FS column
  NR = 0,                                     # NR = 0 for Monoculture
  COT = NA,                                    # Placeholder for COT column
  CS = NA,                                    # Placeholder for CS column
  Treat = NA,                                 # Placeholder for Treat column
  TG_pot = NA,                                # Placeholder for TG_pot column
  F_survival = NA,                            # Placeholder for F_survival column
  Monoculture = "YES",                        # Mark as Monoculture
  variable = rep(c("Shannon_AMF", "Invsimp_AMF"), each = 3),
  value = c(monoculture_AMF$Shannon_AMF, monoculture_AMF$Invsimp_AMF),  # Add values dynamically
  NR.trend = NA,                              # Placeholder for NR.trend
  SE = NA,                                    # Placeholder for SE
  p.value = NA,                               # Placeholder for p-value
  significant = NA                            # Placeholder for significant flag
)

# Combine the additional points with the existing data
data_richness_invsimp_AMF<- bind_rows(data_richness_invsimp_AMF, additional_points_AMF)
# order the variables 
data_richness_invsimp_AMF$variable <- factor(data_richness_invsimp_AMF$variable, levels = c("Shannon_AMF", "Invsimp_AMF"))
significant_lines_AMF$variable <- factor(significant_lines_AMF$variable, levels = c("Shannon_AMF", "Invsimp_AMF"))
line_and_ribbon_data_AMF$variable <- factor(line_and_ribbon_data_AMF$variable, levels = c("Shannon_AMF", "Invsimp_AMF"))
regression_lines_AMF$variable <- factor(regression_lines_AMF$variable, levels = c("Shannon_AMF", "Invsimp_AMF"))

pd = position_dodge(.2)
ggplot(data_richness_invsimp_AMF, aes(x = NR, y = value, color = as.factor(FS))) +
  geom_point(size=2,position=pd,alpha=1,aes(shape = as.factor(NR), fill = FS)) + 
  scale_color_manual(values = c("#6C7C86","#83B8B0","#EAD69D")) +
  ggnewscale::new_scale_color() + ggnewscale::new_scale_fill() + 
  geom_line(data = significant_lines_AMF%>% filter(variable %in% c("Shannon_AMF", "Invsimp_AMF")), aes(x = x, y = y, color = as.factor(FS)), linetype = "solid",alpha = 1) +
  geom_ribbon(data = line_and_ribbon_data_AMF%>% filter(variable %in% c("Shannon_AMF", "Invsimp_AMF")), aes(x = x, ymin = ymin, ymax = ymax, fill = as.factor(FS)), alpha = 0.5, inherit.aes = FALSE) +
  geom_segment(data = regression_lines_AMF%>% filter(!significant & variable %in% c("Shannon_AMF", "Invsimp_AMF")), aes(x = x_min, xend = x_max, y = y_min, yend = y_max, color = as.factor(FS)), linetype = "dashed", linewidth = 0.4,alpha = 0.7) +
  scale_x_continuous(breaks = c(0, 1, 2, 4), labels = c("Monoculture", "1", "2", "4")) +  
  scale_color_manual(values = c("#6C7C86","#83B8B0","#EAD69D")) +
  scale_fill_manual(values = c("#6C7C86","#83B8B0","#EAD69D")) + 
  scale_shape_manual(values = c(1,16,16,16)) +
  facet_wrap(vars(variable), scales = "free", nrow = 2, ncol = 1, labeller = labeller(.multi_line = FALSE)) +
  guides(color = guide_legend(override.aes = list(size = 2)), shape = guide_legend(override.aes = list(size = 2)), alpha = guide_legend(override.aes = list(size = 2))) +
  theme_bw() + mytheme + theme(legend.position = "none") + 
  labs(x = "Neighbour plant species richness", y = NULL, tag = "E") -> Fig_3EF; Fig_3EF
 
Fig_3AB + Fig_3EF


#--------------------------------------
### RDA_Visualization:Fig.3C & G 
#--------------------------------------
### Fig.3C
RDA2 <- rda(t(normData_pathogen) ~ Treat, subSampleTab_pathogen)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont; temp

pathogen_rda = as.data.frame(axes[,1:2])
pathogen_rda$Pot = rownames(pathogen_rda) 
pathogen_rda = pathogen_rda %>% left_join(subSampleTab_pathogen)  
pathogen_rda$NR = as.factor(pathogen_rda$NR)
colnames(pathogen_rda)  

pathogen_rda_mean = pathogen_rda %>% group_by(FS, NR) %>% 
  summarise(RDA1_mean = mean(RDA1), RDA1_se = sd(RDA1)/(sqrt(length(RDA1))),
            RDA2_mean = mean(RDA2), RDA2_se = sd(RDA2)/(sqrt(length(RDA2))))

pathogen_rda_mean2 = merge(pathogen_rda_mean, pathogen_rda, by = c("FS", "NR"))

pathogen_rda_mean$NR = as.factor(pathogen_rda_mean$NR)
pathogen_rda_mean2$NR = as.factor(pathogen_rda_mean2$NR)

ggplot(pathogen_rda_mean, aes(x = RDA1_mean, y = RDA2_mean)) + 
  geom_point(data = pathogen_rda_mean2,mapping = aes(RDA1, RDA2,color =Treat, 
                                                     size = NR),show.legend = T, pch = 16) + 
  scale_color_manual(values = c(alpha("#2D4452",1),alpha("#2D4452",0.7),alpha("#2D4452",0.4),
                                alpha("#4E9A8E",1),alpha("#4E9A8E",0.7),alpha("#4E9A8E",0.4),
                                alpha("#E2C577",1),alpha("#E2C577",0.7),alpha("#E2C577",0.4))) +
  scale_size_manual(values = c(1,2.5,5.5)) +
  ggnewscale::new_scale_color() + ggnewscale::new_scale_fill() + 
  stat_ellipse(data =pathogen_rda_mean2,mapping = aes(x = RDA1, y = RDA2, 
                                                      group=interaction(NR, FS), color = FS,linetype=NR),
               type="norm",geom = "polygon",fill = NA) + 
  labs(x=paste("RDA1 (", format(100 * temp$importance[2,1], digits=3), "%)", sep=""),
       y=paste("RDA2 (", format(100 * temp$importance[2,2], digits=3), "%)", sep=""), tag = "A") +
  theme_bw() + mytheme +  
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))+
  scale_color_manual(values = c(alpha("#2D4452",0.5),alpha("#4E9A8E",0.5),alpha("#E2C577",0.5))) -> Fig_3C; Fig_3C


### Fig.3G 
RDA2 <- rda(t(normData_AMF) ~ Treat, subSampleTab_AMF)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont; temp

AMF_rda = as.data.frame(axes[,1:2])
AMF_rda$Pot = rownames(AMF_rda) 
AMF_rda = AMF_rda %>% left_join(subSampleTab_AMF)  
AMF_rda$NR = as.factor(AMF_rda$NR)
colnames(AMF_rda)  

AMF_rda_mean = AMF_rda %>% group_by(FS, NR) %>% 
  summarise(RDA1_mean = mean(RDA1), RDA1_se = sd(RDA1)/(sqrt(length(RDA1))),
            RDA2_mean = mean(RDA2), RDA2_se = sd(RDA2)/(sqrt(length(RDA2))))

AMF_rda_mean2 = merge(AMF_rda_mean, AMF_rda, by = c("FS", "NR"))

AMF_rda_mean$NR = as.factor(AMF_rda_mean$NR)
AMF_rda_mean2$NR = as.factor(AMF_rda_mean2$NR)

ggplot(AMF_rda_mean, aes(x = RDA1_mean, y = RDA2_mean)) + 
  geom_point(data = AMF_rda_mean2,mapping = aes(RDA1, RDA2,color =Treat, 
                                                size = NR),show.legend = T, pch = 16) + 
  scale_color_manual(values = c(alpha("#2D4452",1),alpha("#2D4452",0.7),alpha("#2D4452",0.4),
                                alpha("#4E9A8E",1),alpha("#4E9A8E",0.7),alpha("#4E9A8E",0.4),
                                alpha("#E2C577",1),alpha("#E2C577",0.7),alpha("#E2C577",0.4))) +
  scale_size_manual(values = c(1,2.5,5.5)) +
  ggnewscale::new_scale_color() + ggnewscale::new_scale_fill() + 
  stat_ellipse(data =AMF_rda_mean2,mapping = aes(x = RDA1, y = RDA2, 
                                                 group=interaction(NR, FS), color = FS,linetype=NR),
               type="norm",geom = "polygon",fill = NA) + 
  labs(x=paste("RDA1 (", format(100 * temp$importance[2,1], digits=3), "%)", sep=""),
       y=paste("RDA2 (", format(100 * temp$importance[2,2], digits=3), "%)", sep=""), tag = "C") +
  theme_bw() + mytheme +  
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))+
  theme_bw() + mytheme + theme(legend.position.inside = c(0.2,0.76)) +
  scale_color_manual(values = c(alpha("#2D4452",0.5),alpha("#4E9A8E",0.5),alpha("#E2C577",0.5))) -> Fig_3G; Fig_3G



 
###########################################################################################
####              Part2---Explained variation (RDA, Shannon, inverse Simpson)         #####
###########################################################################################

#--------------------------------------
### TableS5
#--------------------------------------
#overall-1------------------- ### Explained variation_overall_composition
variance1 <- as.data.frame(RDA_result_overall$Variance)
rownames(variance1) <- c("FS", "NR", "interaction","Residual");head(variance1)
variance1 <- variance1 %>% mutate_if(is.numeric, function(x){x/sum(x)});variance1
Total_percent1 <- 1 -variance1$`RDA_result_overall$Variance`[4];Total_percent1 # this is the model explained, should be on top of the bar

Total_explained_variance1 <- sum(variance1$`RDA_result_overall$Variance`[1:3]);Total_explained_variance1  # Variance explained by FS, NR, and interaction
variance11 <- variance1 %>% mutate( Relative_Percentage = (`RDA_result_overall$Variance`/ Total_explained_variance1) * 100);variance11
variance_overall_com<- variance11[1:3, c("Relative_Percentage")];variance_overall_com # Print the relative contributions

#overall-2--------------------### Explained variation_overall_Shannon 
model_overall2 <- glmmTMB(Shannon_overall ~ FS * NR + (1 | COT), data = div_dat_overall)
model_overall_no_FS2 <- glmmTMB(Shannon_overall ~ NR + (1 | COT), data = div_dat_overall)
model_overall_no_NR2 <- glmmTMB(Shannon_overall ~ FS + (1 | COT), data = div_dat_overall)
model_overall_no_interaction2 <- glmmTMB(Shannon_overall ~ FS + NR + (1 | COT), data = div_dat_overall)

R2_full2 <- r.squaredGLMM(model_overall2);R2_full2 # this is the model explained, should be on top of the bar
R2_no_FS2 <- r.squaredGLMM(model_overall_no_FS2);R2_no_FS2
R2_no_NR2 <- r.squaredGLMM(model_overall_no_NR2);R2_no_NR2
R2_no_interaction2 <- r.squaredGLMM(model_overall_no_interaction2);R2_no_interaction2

#R2 result
Residual_variance2 <- 1-R2_full2[1];Residual_variance2
R2_interaction2 <- R2_full2[1] - R2_no_interaction2[1];R2_interaction2  # variance explained by interaction
R2_FS2 <- (R2_full2[1] - R2_no_FS2[1]) - R2_interaction2;R2_FS2 # variance explained by FS alone
R2_NR2 <- (R2_full2[1] - R2_no_NR2[1]) - R2_interaction2;R2_NR2 # variance explained by NR alone
R2_overall_Shannon <- na.omit(c(R2_FS2,R2_NR2,R2_interaction2));R2_overall_Shannon 

#Percentage contributions
variance2 <- as.data.frame(R2_overall_Shannon);variance2
Total_explained_variance2 <- sum(variance2$`R2_overall_Shannon`[1:3]);Total_explained_variance2  # Variance explained by FS, NR, and interaction
variance22 <- variance2 %>% mutate( Relative_Percentage = (`R2_overall_Shannon`[1:3]/ Total_explained_variance2) * 100);variance22
variance_overall_Shannon<- variance22[1:3, c("Relative_Percentage")];variance_overall_Shannon # Print the percentage contributions

#overall-3--------------------### Explained variation_overall_Invsimp 
model_overall3 <- glmmTMB(Invsimp_overall ~ FS * NR + (1 | COT), data = div_dat_overall)
model_overall_no_FS3 <- glmmTMB(Invsimp_overall ~ NR + (1 | COT), data = div_dat_overall)
model_overall_no_NR3 <- glmmTMB(Invsimp_overall ~ FS + (1 | COT), data = div_dat_overall)
model_overall_no_interaction3 <- glmmTMB(Invsimp_overall ~ FS + NR + (1 | COT), data = div_dat_overall)

R2_full3 <- r.squaredGLMM(model_overall3);R2_full3# this is the model explained, should be on top of the bar
R2_no_FS3 <- r.squaredGLMM(model_overall_no_FS3);R2_no_FS3
R2_no_NR3 <- r.squaredGLMM(model_overall_no_NR3);R2_no_NR3
R2_no_interaction3 <- r.squaredGLMM(model_overall_no_interaction3);R2_no_interaction3

#R2 result
Residual_variance3 <- 1-R2_full3[1];Residual_variance3
R2_interaction3 <- R2_full3[1] - R2_no_interaction3[1];R2_interaction3 # variance explained by interaction
R2_FS3 <- (R2_full3[1] - R2_no_FS3[1]) - R2_interaction3;R2_FS3  # variance explained by FS alone
R2_NR3 <- (R2_full3[1] - R2_no_NR3[1]) - R2_interaction3;R2_NR3 # variance explained by NR alone
R2_overall_Invsimp <- na.omit(c(R2_FS3,R2_NR3,R2_interaction3));R2_overall_Invsimp 

#Percentage contributions
variance3 <- as.data.frame(R2_overall_Invsimp);variance3
Total_explained_variance3 <- sum(variance3$`R2_overall_Invsimp`[1:3]);Total_explained_variance3  # Variance explained by FS, NR, and interaction
variance33 <- variance3 %>% mutate( Relative_Percentage = (`R2_overall_Invsimp`[1:3]/ Total_explained_variance3) * 100);variance33
variance_overall_Invsimp<- variance33[1:3, c("Relative_Percentage")];variance_overall_Invsimp # Print the relative contributions


### Explained variation_pathogen
#pathogen-1--------------------### Explained variation_pathogen_composition 
variance4 <- as.data.frame(RDA_result_pathogen$Variance)
rownames(variance4) <- c("FS", "NR", "interaction","Residual");head(variance4)
variance4 <- variance4 %>% mutate_if(is.numeric, function(x){x/sum(x)});variance4
Total_percent4 <- 1 -variance4$`RDA_result_pathogen$Variance`[4];Total_percent4 # this is the model explained, should be on top of the bar

Total_explained_variance4 <- sum(variance4$`RDA_result_pathogen$Variance`[1:3]);Total_explained_variance4  # Variance explained by FS, NR, and interaction
variance44 <- variance4 %>% mutate( Relative_Percentage = (`RDA_result_pathogen$Variance`/ Total_explained_variance4) * 100);variance44
variance_pathogen_com<- variance44[1:3, c("Relative_Percentage")];variance_pathogen_com # Print the relative contributions

#pathogen-2--------------------### Explained variation_pathogen_Shannon 
model_pathogen5 <- glmmTMB(Shannon_pathogen ~ FS * NR + (1 | COT), data = div_dat_pathogen);Anova(model_pathogen5,type="III")
model_pathogen_no_FS5 <- glmmTMB(Shannon_pathogen ~ NR + (1 | COT), data = div_dat_pathogen)
model_pathogen_no_NR5 <- glmmTMB(Shannon_pathogen ~ FS + (1 | COT), data = div_dat_pathogen)
model_pathogen_no_interaction5 <- glmmTMB(Shannon_pathogen ~ FS + NR + (1 | COT), data = div_dat_pathogen)

R2_full5 <- r.squaredGLMM(model_pathogen5);R2_full5 # this is the model explained, should be on top of the bar
R2_no_FS5 <- r.squaredGLMM(model_pathogen_no_FS5);R2_no_FS5
R2_no_NR5 <- r.squaredGLMM(model_pathogen_no_NR5);R2_no_NR5
R2_no_interaction5 <- r.squaredGLMM(model_pathogen_no_interaction5);R2_no_interaction5

#R2 result
Residual_variance5 <- 1-R2_full5[1];Residual_variance5
R2_interaction5 <- R2_full5[1] - R2_no_interaction5[1];R2_interaction5  # variance explained by interaction
R2_FS5 <- (R2_full5[1] - R2_no_FS5[1]) - R2_interaction5;R2_FS5 # variance explained by FS alone
R2_NR5 <- (R2_full5[1] - R2_no_NR5[1]) - R2_interaction5;R2_NR5 # variance explained by NR alone
R2_pathogen_Shannon <- na.omit(c(R2_FS5,R2_NR5,R2_interaction5));R2_pathogen_Shannon 

#Percentage contributions
variance5 <- as.data.frame(R2_pathogen_Shannon);variance5
Total_explained_variance5 <- sum(variance5$`R2_pathogen_Shannon`[1:3]);Total_explained_variance5  # Variance explained by FS, NR, and interaction
variance55 <- variance5 %>% mutate( Relative_Percentage = (`R2_pathogen_Shannon`[1:3]/ Total_explained_variance5) * 100);variance55
variance_pathogen_Shannon<- variance55[1:3, c("Relative_Percentage")];variance_pathogen_Shannon # Print the relative contributions

#pathogen-3--------------------### Explained variation_pathogen_Invsimp
model_pathogen6 <- glmmTMB(Invsimp_pathogen ~ FS * NR + (1 | COT), data = div_dat_pathogen);Anova(model_pathogen6,type="III")
model_pathogen_no_FS6 <- glmmTMB(Invsimp_pathogen ~ NR + (1 | COT), data = div_dat_pathogen)
model_pathogen_no_NR6 <- glmmTMB(Invsimp_pathogen ~ FS + (1 | COT), data = div_dat_pathogen)
model_pathogen_no_interaction6 <- glmmTMB(Invsimp_pathogen ~ FS + NR + (1 | COT), data = div_dat_pathogen)

R2_full6 <- r.squaredGLMM(model_pathogen6);R2_full6 # this is the model explained, should be on top of the bar
R2_no_FS6 <- r.squaredGLMM(model_pathogen_no_FS6);R2_no_FS6 
R2_no_NR6 <- r.squaredGLMM(model_pathogen_no_NR6);R2_no_NR6 
R2_no_interaction6 <- r.squaredGLMM(model_pathogen_no_interaction6);R2_no_interaction6  

#R2 result
Residual_variance6 <- 1-R2_full6[1];Residual_variance6
R2_interaction6 <- R2_full6[1] - R2_no_interaction6[1];R2_interaction6  # variance explained by interaction
R2_FS6 <- (R2_full6[1] - R2_no_FS6[1]) - R2_interaction6 # variance explained by FS alone
R2_NR6 <- (R2_full6[1] - R2_no_NR6[1]) - R2_interaction6 # variance explained by NR alone
R2_pathogen_Invsimp <- na.omit(c(R2_FS6,R2_NR6,R2_interaction6));R2_pathogen_Invsimp 

#Percentage contributions
variance6 <- as.data.frame(R2_pathogen_Invsimp);variance6
Total_explained_variance6 <- sum(variance6$`R2_pathogen_Invsimp`[1:3]);Total_explained_variance6  # Variance explained by FS, NR, and interaction
variance66 <- variance6 %>% mutate( Relative_Percentage = (`R2_pathogen_Invsimp`[1:3]/ Total_explained_variance6) * 100);variance66
variance_pathogen_Invsimp<- variance66[1:3, c("Relative_Percentage")];variance_pathogen_Invsimp # Print the relative contributions


### Explained variation_AMF
#AMF-1--------------------### Explained variation_AMF_composition 
variance7 <- as.data.frame(RDA_result_AMF$Variance)
rownames(variance7) <- c("FS", "NR", "interaction","Residual");head(variance7)
variance7 <- variance7 %>% mutate_if(is.numeric, function(x){x/sum(x)});variance7
Total_percent7 <- 1 -variance7$`RDA_result_AMF$Variance`[4];Total_percent7 # this is the model explained, should be on top of the bar

Total_explained_variance7 <- sum(variance7$`RDA_result_AMF$Variance`[1:3]);Total_explained_variance7  # Variance explained by FS, NR, and interaction
variance77 <- variance7 %>% mutate( Relative_Percentage = (`RDA_result_AMF$Variance`/ Total_explained_variance7) * 100);variance77
variance_AMF_com<- variance77[1:3, c("Relative_Percentage")];variance_AMF_com # Print the relative contributions

#AMF-2--------------------### Explained variation_AMF_Shannon 
model_AMF8 <- glmmTMB(Shannon_AMF ~ FS * NR + (1 | COT), data = div_dat_AMF);Anova(model_AMF8,type="III")
model_AMF_no_FS8 <- glmmTMB(Shannon_AMF ~ NR + (1 | COT), data = div_dat_AMF)
model_AMF_no_NR8 <- glmmTMB(Shannon_AMF ~ FS + (1 | COT), data = div_dat_AMF)
model_AMF_no_interaction8 <- glmmTMB(Shannon_AMF ~ FS + NR + (1 | COT), data = div_dat_AMF)

R2_full8 <- r.squaredGLMM(model_AMF8);R2_full8 # this is the model explained, should be on top of the bar
R2_no_FS8 <- r.squaredGLMM(model_AMF_no_FS8);R2_no_FS8
R2_no_NR8 <- r.squaredGLMM(model_AMF_no_NR8);R2_no_NR8
R2_no_interaction8 <- r.squaredGLMM(model_AMF_no_interaction8);R2_no_interaction8

#R2 result
Residual_variance8 <- 1-R2_full8[1];Residual_variance8
R2_interaction8 <- R2_full8[1] - R2_no_interaction8[1];R2_interaction8  # variance explained by interaction
R2_FS8 <- (R2_full8[1] - R2_no_FS8[1]) - R2_interaction8 # variance explained by FS alone
R2_NR8 <- (R2_full8[1] - R2_no_NR8[1]) - R2_interaction8 # variance explained by NR alone
R2_AMF_Shannon <- na.omit(c(R2_FS8,R2_NR8,R2_interaction8));R2_AMF_Shannon 

#Percentage contributions
variance8 <- as.data.frame(R2_AMF_Shannon);variance8
Total_explained_variance8 <- sum(variance8$`R2_AMF_Shannon`[1:3]);Total_explained_variance8  # Variance explained by FS, NR, and interaction
variance88 <- variance8 %>% mutate( Relative_Percentage = (`R2_AMF_Shannon`[1:3]/ Total_explained_variance8) * 100);variance88
variance_AMF_Shannon<- variance88[1:3, c("Relative_Percentage")];variance_AMF_Shannon # Print the relative contributions

#AMF-3--------------------### Explained variation_AMF_Invsimp
model_AMF9 <- glmmTMB(Invsimp_AMF ~ FS * NR + (1 | COT), data = div_dat_AMF)
model_AMF_no_FS9 <- glmmTMB(Invsimp_AMF ~ NR + (1 | COT), data = div_dat_AMF)
model_AMF_no_NR9 <- glmmTMB(Invsimp_AMF ~ FS + (1 | COT), data = div_dat_AMF)
model_AMF_no_interaction9 <- glmmTMB(Invsimp_AMF ~ FS + NR + (1 | COT), data = div_dat_AMF)

R2_full9 <- r.squaredGLMM(model_AMF9);R2_full9 # this is the model explained, should be on top of the bar
R2_no_FS9 <- r.squaredGLMM(model_AMF_no_FS9);R2_no_FS9
R2_no_NR9 <- r.squaredGLMM(model_AMF_no_NR9);R2_no_NR9
R2_no_interaction9 <- r.squaredGLMM(model_AMF_no_interaction9);R2_no_interaction9

#R2 result
Residual_variance9 <- 1-R2_full9[1];Residual_variance9
R2_interaction9 <- R2_full9[1] - R2_no_interaction9[1];R2_interaction9    # variance explained by interaction
R2_FS9 <- (R2_full9[1] - R2_no_FS9[1]) - R2_interaction9 # variance explained by FS alone
R2_NR9 <- (R2_full9[1] - R2_no_NR9[1]) - R2_interaction9 # variance explained by NR alone
R2_AMF_Invsimp <- na.omit(c(R2_FS9,R2_NR9,R2_interaction9));R2_AMF_Invsimp 

#Percentage contributions
variance9 <- as.data.frame(R2_AMF_Invsimp);variance9
Total_explained_variance9 <- sum(variance9$`R2_AMF_Invsimp`[1:3]);Total_explained_variance9  # Variance explained by FS, NR, and interaction
variance99 <- variance9 %>% mutate( Relative_Percentage = (`R2_AMF_Invsimp`[1:3]/ Total_explained_variance9) * 100);variance99
variance_AMF_Invsimp<- variance99[1:3, c("Relative_Percentage")];variance_AMF_Invsimp # Print the relative contributions


### TableS5_results
# overall fungi
variance22 
Residual_variance2
variance33
Residual_variance3
variance11
# Pathogens
variance55 
Residual_variance5
variance66
Residual_variance6
variance44
# AMF
variance88 
Residual_variance8
variance99
Residual_variance9
variance77
 

#model explained: on top of the bartop of the bar (Fig.3G H)
Top_bar <- na.omit(c(#R2_full2[1],R2_full3[1],Total_percent1,
                     R2_full5[1], R2_full6[1],Total_percent4,
                     R2_full8[1], R2_full9[1],Total_percent7));Top_bar


#--------------------------------------
### Visualization:  Fig.3D & H
#--------------------------------------
### Fig.3D 
variance_pathogen = data.frame(factor = rep(c("Focal species", "Neighbouring species richness", "Interaction"), 3),
                              group = c(rep(c("Composition"), 3), rep(c("Shannon"), 3), rep(c("Inverse Simpson"), 3)),
                              value =c(variance_pathogen_com, variance_pathogen_Shannon, variance_pathogen_Invsimp));variance_pathogen
variance_pathogen$value = variance_pathogen$value;variance_pathogen
variance_pathogen$factor = factor(variance_pathogen$factor, levels = c("Focal species", "Neighbouring species richness", "Interaction"))
variance_pathogen$group = factor(variance_pathogen$group, levels = c( "Shannon","Inverse Simpson", "Composition"));variance_pathogen

ggplot(variance_pathogen, aes(fill=factor, y=value, x=group)) + 
  geom_bar(stat="identity", color= "black", width = 0.8) + 
  scale_fill_manual(name = NULL,values = c("#70A7C3","#A67C2A","#D2BEA2"))+
  theme_bw() + mytheme + 
  theme(axis.text.x = element_text(colour='black', size=14, angle = 25, hjust = 1, vjust = 1),
        legend.position = c(0.3,0.91)) +
  scale_y_continuous(expand = c(0,1), limits = c(0,130),breaks = seq(0, 130,50)) + 
  geom_vline(aes(xintercept = 2.5), linetype = "dashed") + 
  labs(x = NULL, y = "Percentage contribution (%)", tag = "B")  -> Fig_3D; Fig_3D


### Fig.3H
variance_AMF = data.frame(factor = rep(c("Focal species", "Neighbouring species richness", "Interaction"), 3),
                          group = c(rep(c("Composition"), 3), rep(c("Shannon"), 3), rep(c("Inverse Simpson"), 3)),
                          value =c(variance_AMF_com, variance_AMF_Shannon, variance_AMF_Invsimp));variance_AMF
variance_AMF$value = variance_AMF$value;variance_AMF
variance_AMF$factor = factor(variance_AMF$factor, levels = c("Focal species", "Neighbouring species richness", "Interaction"))
variance_AMF$group = factor(variance_AMF$group, levels = c( "Shannon","Inverse Simpson", "Composition"));variance_AMF

ggplot(variance_AMF, aes(fill=factor, y=value, x=group)) + 
  geom_bar(stat="identity", color= "black", width = 0.8) + 
  scale_fill_manual(name = NULL,values = c("#70A7C3","#A67C2A","#D2BEA2"))+
  theme_bw() + mytheme + 
  theme(axis.text.x = element_text(colour='black', size=14, angle = 25, hjust = 1, vjust = 1),
        legend.position = c(0.3,0.91)) +
  scale_y_continuous(expand = c(0,1), limits = c(0,130),breaks = seq(0, 130,50)) + 
  geom_vline(aes(xintercept = 2.5), linetype = "dashed") + 
  labs(x = NULL, y = "Percentage contribution (%)", tag = "B")  -> Fig_3H; Fig_3H
 

(Fig_3AB|Fig_3EF)/(Fig_3C|Fig_3G)/(Fig_3D|Fig_3H)->Fig.3;Fig.3


 
