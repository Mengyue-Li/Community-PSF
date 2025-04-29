###########################################################################################
#####                                  Library & mytheme                              #####
###########################################################################################
# Load required packages and data
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggepi)
library(patchwork)
library(party)
library(caret)
library(readxl)
library(glmmTMB)
library(dbplyr)
library(tidyr)
library(tidyverse)
library(randomForest)

##################################################################################
#####                                                                        ##### 
#####           Part1---effect of home vs away:PSF                           #####
#####                                                                        #####  
##################################################################################

setwd('C:/Users/MY/Desktop/CODE/2')

#--------------------------------------
### Table S8
#--------------------------------------
dat <- read_excel("data2.xlsx",sheet=1);head(dat)
# Fit lm model with totol biomass
mod_full <- lm(TG ~ Richness_con * (B_con + C_con + F_con + T_con + V_con), data = dat)
qqnorm(resid(mod_full));qqline(resid(mod_full));anova(mod_full)-> mod_full_result; mod_full_result
p <- mod_full_result$Pr;p 
p.adjust(p, "BH")
p.adjust(p, "bonferroni")


##################################################################################
#####                                                                        ##### 
#####   Part2---Soil-mediated effects of single species and richness         #####
#####                                                                        #####  
##################################################################################
df <- read_excel("data2.xlsx",sheet=2);head(df)
df <- as.data.frame(df); rownames(df) <- df[, 1]

# Define constants
stressors <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "a", "b")
levels <- c("1", "2", "4", "6", "12")
responses <- "TG"
n_iter <- 1000
target <- c("CK", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "a", "b","1", "2", "4", "6", "12")
set.seed(123)
mytheme= theme(legend.position = "none",
               panel.grid=element_blank(), 
               legend.title = element_blank(),
               legend.text = element_text(size = 14), 
               legend.background = element_rect(fill = NA),  
               axis.ticks = element_line(color='black'),
               axis.line = element_line(colour = "black"), 
               axis.title.x = element_text(colour='black', size=14),
               axis.title.y = element_text(colour='black', size=14),
               axis.text = element_text(colour='black',size=11),
               plot.tag = element_text(size = 0))
#----------------------------------------------
# Estimating mean and its 95% confidence interval
#----------------------------------------------
BootStrap_mean <- function(response, data = df, n_perm = n_iter) {  # Removed target default
  summary <- list()
  
  for (treatment in target) {
    bs <- numeric(0)
    if (treatment == "1") {
      population <- data[data$remark %in% stressors, response]
    } else {
      population <- data[data$remark == treatment, response]
    }
    size <- length(population)
    
    if (size > 0) {  # Ensure there's data to sample
      for (id in 1:n_perm) {
        k <- mean(sample(population, size, replace = TRUE), na.rm = TRUE)
        bs <- append(bs, k)
      }
      summary[[treatment]] <- c(quantile(bs, 0.025, na.rm = TRUE), 
                                mean(bs, na.rm = TRUE), 
                                quantile(bs, 0.975, na.rm = TRUE))
      names(summary[[treatment]]) <- c("2.5%", "mean", "97.5%")
    } else {
      summary[[treatment]] <- c(NA, NA, NA)  # Handle empty cases
      names(summary[[treatment]]) <- c("2.5%", "mean", "97.5%")
    }
  }
  summary <- t(data.frame(summary))
  summary <- data.frame("target" = target, summary)
  row.names(summary) <- NULL
  return(summary)
}
response_mean <- BootStrap_mean("TG")


response_mean$target <- factor(response_mean$target,levels = c("CK", LETTERS,"a","b", "1", "2", "4", "6", "12"))
df$remark <- factor(df$remark,levels = levels(response_mean$target)  # Inherit the same order
)
Pa <- ggplot() +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  # Highlight region for levels 1 to 12
  geom_rect(data = NULL, 
            aes(xmin = -Inf, xmax = Inf, 
                ymin = which(levels(df$remark) == "1") - 0.5, 
                ymax = which(levels(df$remark) == "12") + 0.5), 
            fill = "#D3D3D3", alpha = 0.4, inherit.aes = FALSE) +
  # Control mean as vertical dashed line
  geom_vline(aes(xintercept = subset(response_mean, target == "CK")$mean),linetype = "dashed", color = "black") +
  coord_flip() +
  
  # Control as open purple circle
  stat_density_ridges(data = subset(df, remark == "CK"), 
                      aes(x = TG, y = remark, fill = ..x..), 
                      geom = "density_ridges_gradient", bandwidth = 5.82, rel_min_height = 0.01, 
                      jittered_points = TRUE, color = alpha("#4E307D", 0.5), fill = alpha("#4E307D", 0.3), 
                      alpha = 0.4, position = position_points_jitter(height = 0.15, yoffset = 0.1), scale = 0.6) +
  geom_point(data = subset(response_mean, target == "CK"), 
             aes(x = mean, y = target), shape = 1, size = 2, stroke = 1.2, color = "#4E307D") +
  geom_errorbar(data = subset(response_mean, target == "CK"), 
                aes(xmin = X2.5., xmax = X97.5., y = target), width = 0, color = "#4E307D", linewidth = 0.5) +
  
  # Single stressors (excluding levels) as purple filled circles
  stat_density_ridges(data = subset(df, !remark %in% c("CK", "1", "2", "4", "6", "12")), 
                      aes(x = TG, y = remark, fill = ..x..), 
                      geom = "density_ridges_gradient", bandwidth = 5.02, rel_min_height = 0.01, 
                      jittered_points = TRUE, fill = alpha("#4E307D", 0.2), color = alpha("#4E307D", 0.2), 
                      alpha = 0.4, position = position_points_jitter(height = 0.15, yoffset = 0.1), scale = 0.6) +
  geom_point(data = subset(response_mean, !target %in% c("CK", "1", "2", "4", "6", "12")), 
             aes(x = mean, y = target), shape = 16, size = 3, color = "#4E307D") +
  geom_errorbar(data = subset(response_mean, !target %in% c("CK", "1", "2", "4", "6", "12")), 
                aes(xmin = X2.5., xmax = X97.5., y = target), width = 0, color = "#4E307D", linewidth = 0.5) +
  
  # Target levels (1, 2, 4, 6, 12) as black filled circles
  stat_density_ridges(data = subset(df, remark %in% c("1", "2", "4", "6", "12")), 
                      aes(x = TG, y = remark, fill = ..x..), 
                      geom = "density_ridges_gradient", bandwidth = 4.03, rel_min_height = 0.01, 
                      jittered_points = TRUE, color = alpha("#000000", 0.2), fill = alpha("#000000", 0.2), 
                      alpha = 0.1, position = position_points_jitter(height = 0.15, yoffset = 0.1), scale = 0.6) +
  geom_point(data = subset(response_mean, target %in% c("1", "2", "4", "6", "12")), 
             aes(x = mean, y = target), shape = 16, size = 3, color = "#000000") +
  geom_errorbar(data = subset(response_mean, target %in% c("1", "2", "4", "6", "12")), 
                aes(xmin = X2.5., xmax = X97.5., y = target), width = 0, color = "#000000", linewidth = 0.5) +
  
  labs(x = "Total biomass of responding community (g)", y = "Single plant species    Species richness ", tag = "a") +
  mytheme +  # Assuming mytheme is defined; remove if not
  scale_x_continuous(limits = c(25, 123), expand = c(0, 0), breaks = seq(40, 120, 20))
Pa


##################################################################################
#####                                                                        ##### 
#####                           Part3---psf prediction                       #####
#####                                                                        #####  
##################################################################################
 
# Null distribution function (without Dominative)
Null_distribution_rep_log <- function(response, data = df, n_perm = n_iter) {
  output <- list()
  data <- data %>% mutate(across(all_of(response), ~ ifelse(. == 0, 1e-6, .)))
  
  for (Lv in levels) {
    resampled <- list()
    if (Lv == "1") {
      combination <- data[data$remark %in% stressors, 2:29]  # Columns A to b
    } else {
      combination <- data[data$Lv == Lv, 2:29]
    }
    control_log_mean <- log(mean(data[data$remark == "CK", response], na.rm = TRUE))
    
    for (type in c("Additive", "Weighted_Additive")) {  # Removed "Dominative"
      resampled[[type]] <- list()
      single_log_effects <- sapply(stressors, function(s) {
        tr_mean <- mean(data[data$remark == s, response], na.rm = TRUE)
        log(tr_mean) - control_log_mean
      })
      
      for (j in 1:nrow(combination)) {
        bs <- numeric(n_perm)
        selected_stressors_idx <- which(combination[j, ] > 0)
        selected_stressors <- stressors[selected_stressors_idx]
        
        species_weights <- if (type == "Weighted_Additive" && length(selected_stressors) > 0) {
          biomass <- unlist(combination[j, selected_stressors_idx])
          biomass / sum(biomass)
        } else NULL
        
        for (id in 1:n_perm) {
          log_effects <- if (length(selected_stressors) > 0) {
            sapply(selected_stressors, function(s) {
              tr_sample <- sample(data[data$remark == s, response], replace = TRUE)
              log(mean(tr_sample, na.rm = TRUE)) - control_log_mean
            })
          } else {
            numeric(0)
          }
          
          joint_effect <- if (length(log_effects) == 0) {
            NA
          } else {
            switch(
              type,
              "Additive" = sum(log_effects),
              "Weighted_Additive" = sum(log_effects * species_weights)
            )
          }
          
          bs[id] <- joint_effect
        }
        resampled[[type]][[j]] <- bs
      }
    }
    output[[Lv]] <- resampled
  }
  return(output)
}

# Bootstrap actual effect sizes
BootStrap_ES_rep <- function(response, data = df, n_perm = n_iter) {
  data <- data %>% mutate(across(all_of(response), ~ ifelse(. == 0, 1e-6, .)))
  resampled <- list()
  control_mean <- mean(data[data$remark == "CK", response], na.rm = TRUE)
  
  for (Lv in levels) {
    if (Lv == "1") {
      tr_data <- data[data$remark %in% stressors, response]
    } else {
      tr_data <- data[!is.na(data$Lv) & data$Lv == Lv, response]
    }
    
    if (length(tr_data) > 0 && !all(is.na(tr_data))) {
      bs <- replicate(n_perm, {
        sample_mean <- mean(sample(tr_data, replace = TRUE), na.rm = TRUE)
        log(sample_mean / control_mean)
      })
    } else {
      bs <- rep(NA, n_perm)
    }
    
    resampled[[Lv]] <- bs
  }
  return(resampled)
}

# Hypothesis testing and summary
NHST_summary <- function(null_data, Actual_data) {
  output <- list()
  for (Lv in levels) {
    summary <- list()
    
    actual_vec <- as.numeric(Actual_data[[Lv]])
    if (length(actual_vec) == 0 || all(is.na(actual_vec))) {
      summary[["Actual"]] <- c(NA, NA, NA, NA)
    } else {
      summary[["Actual"]] <- c(
        quantile(actual_vec, 0.025, na.rm = TRUE),
        mean(actual_vec, na.rm = TRUE),
        quantile(actual_vec, 0.975, na.rm = TRUE),
        NA
      )
    }
    
    assumptions <- c("Additive", "Weighted_Additive")  # Removed "Dominative"
    for (i_assumption in assumptions) {
      null_vec <- as.numeric(unlist(null_data[[Lv]][[i_assumption]]))
      if (length(null_vec) == 0 || all(is.na(null_vec)) || length(actual_vec) == 0 || all(is.na(actual_vec))) {
        summary[[i_assumption]] <- c(NA, NA, NA, NA)
      } else {
        n_samples <- min(length(actual_vec), length(null_vec))
        actual_sample <- actual_vec[1:n_samples]
        null_sample <- null_vec[1:n_samples]
        
        diff <- actual_sample - null_sample
        p <- 2 * min(mean(diff > 0, na.rm = TRUE), mean(diff < 0, na.rm = TRUE))
        
        summary[[i_assumption]] <- c(
          quantile(null_vec, 0.025, na.rm = TRUE),
          mean(null_vec, na.rm = TRUE),
          quantile(null_vec, 0.975, na.rm = TRUE),
          p
        )
      }
    }
    
    summary_df <- as.data.frame(do.call(rbind, summary))
    colnames(summary_df) <- c("2.5%", "mean", "97.5%", "p_value")
    summary_df$ES <- c("Actual", assumptions)
    output[[Lv]] <- summary_df
  }
  return(output)
}

# Transform for plotting
NHST_summary_transform <- function(nhst_summary) {
  plot_data <- do.call(rbind, lapply(names(nhst_summary), function(Lv) {
    df <- nhst_summary[[Lv]]
    data.frame(
      Lv = Lv,
      Low = df[, "2.5%"],
      Mean = df[, "mean"],
      High = df[, "97.5%"],
      Model = df$ES,
      P_value = df$p_value
    )
  }))
  return(plot_data)
}

# Run the analysis
set.seed(123)
null_dist <- Null_distribution_rep_log(response = "TG", n_perm = 1000)
Actual_data <- BootStrap_ES_rep(response = "TG", n_perm = 1000)
nhst_summary <- NHST_summary(null_dist, Actual_data)
ES_plot_data <- NHST_summary_transform(nhst_summary)

# Print results
print(ES_plot_data)
ES_plot_data$Lv <- factor(ES_plot_data$Lv, levels = as.numeric(levels), labels = levels)
# Plot
Pb <- ggplot() +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(color = "black", size = 12), 
    axis.title.y = element_text(color = "black", size = 12)) +
  labs(x = "Plant-soil feedback effects
Ln(biomass in conditioned soil / biomass in sterilized)", y = "Species richness", tag = "b") +
  coord_flip() +
  scale_y_discrete(limits = factor(c("1", "2", "4", "6", "12"), levels=c("1", "2", "4", "6", "12"))) +
  ### Mean & CI for 2 assumptions + Actual (removed Dominative)
  geom_estci(data = ES_plot_data[ES_plot_data$Model == "Additive" & ES_plot_data$Lv %in% c("1", "2", "4", "6", "12"), ], 
             aes(x = Mean, y = Lv, xmin = Low, xmax = High, xintercept = 0), 
             color = "#7B3D12", size = 0.6, ci.linesize = 0.5, position = position_nudge(y = +0.1)) +
  geom_estci(data = ES_plot_data[ES_plot_data$Model == "Weighted_Additive" & ES_plot_data$Lv %in% c("1", "2", "4", "6", "12"), ], 
             aes(x = Mean, y = Lv, xmin = Low, xmax = High, xintercept = 0), 
             color = "#DE8125", size = 0.6, ci.linesize = 0.5, position = position_nudge(y = +0.2)) +
  geom_estci(data = ES_plot_data[ES_plot_data$Model == "Actual" & ES_plot_data$Lv %in% c("1", "2", "4", "6", "12"), ], 
             aes(x = Mean, y = Lv, xmin = Low, xmax = High, xintercept = 0), 
             color = "#000000", size = 0.6, ci.linesize = 0.5, position = position_nudge(y = 0)) +
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01), 
                     limits = c(-3.5, 0),  # Adjusted limits based on your data
                     expand = c(0, 0), 
                     breaks = seq(-3.5, 0, 1)) + mytheme
Pb


##################################################################################
#####                                                                        ##### 
#####                             Part4---explanations                       #####
#####                                                                        #####  
##################################################################################

#### explanations - Modified to only include the three R² values requested
set.seed(123)
#----------------------------------------------
# Step 1: Data preparation
#----------------------------------------------
# Filter dataframe for treatment levels (excluding control)
df.rf <- df[!is.na(df$Lv) & df$Lv %in% levels, ]

# Extract effect sizes from ES_plot_data (only Additive and Weighted_Additive)
es_data <- ES_plot_data %>%
  filter(Model %in% c("Additive", "Weighted_Additive")) %>% 
  dplyr::select(Lv, Model, Mean) %>%
  pivot_wider(names_from = Model, values_from = Mean)

# Merge effect sizes with df.rf based on Lv
df.rf <- merge(df.rf, es_data, by = "Lv", all.x = TRUE)

# Define constants
n_tree <- 1000       # Number of trees in each forest
n_iter <- 1000       # Number of bootstrap iterations
cv_prop <- 0.7       # Proportion for training in cross-validation

# Define the response variable and stressor variables
response_var <- "TG"  # Target variable

#----------------------------------------------
# Step 2: Define model formulas
#----------------------------------------------
# Model 1: Baseline with species richness only
fml_Lv <- formula(paste(response_var, "~ Lv"))

# Model 2: Refined with species richness + species identity
fml_LvID <- formula(paste(response_var, "~ Lv +", paste(stressors, collapse = " + ")))

# Model 3: Full model with richness, identity, and effect sizes
fml_Full <- formula(paste(response_var, "~ Lv +", 
                          paste(stressors, collapse = " + "), "+",
                          paste(c("Additive", "Weighted_Additive"), collapse = " + ")))

# Count predictors for dynamic mtry parameter
n_pred_Lv <- length(all.vars(fml_Lv)) - 1
n_pred_LvID <- length(all.vars(fml_LvID)) - 1
n_pred_Full <- length(all.vars(fml_Full)) - 1

#----------------------------------------------
# Step 3: Prepare results storage
#----------------------------------------------
# Storage for all iterations
results_rf <- data.frame(
  Iteration = rep(1:n_iter, each = 3),
  Model = rep(c("Species_richness", "Species_identity", "Effect_size"), n_iter),
  R2_RF = NA,            # For standard RandomForest
  stringsAsFactors = FALSE
)

# Storage for additional model metrics
model_metrics <- list()

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
  
  #----------------------
  # Standard RandomForest with train/test split
  #----------------------
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

#----------------------------------------------
# Step 5: Calculate summary statistics
#----------------------------------------------
# Convert Model to a factor with a defined order
results_rf$Model <- factor(results_rf$Model, 
                           levels = c("Species_richness", "Species_identity", "Effect_size"))

# Calculate summary statistics
rf_summary <- results_rf %>%
  group_by(Model) %>%
  summarize(
    CI.low = quantile(R2_RF, 0.025, na.rm = TRUE),
    Mean = mean(R2_RF, na.rm = TRUE),
    Median = median(R2_RF, na.rm = TRUE),
    CI.high = quantile(R2_RF, 0.975, na.rm = TRUE),
    SD = sd(R2_RF, na.rm = TRUE)
  )

# Calculate improvement percentages over baseline
baseline_rf <- rf_summary$Mean[rf_summary$Model == "Species_richness"]

rf_summary <- rf_summary %>%
  mutate(
    Improvement = (Mean - baseline_rf) / baseline_rf * 100
  )

# Print summary results
print("Summary of R² values:")
print(rf_summary)

#----------------------------------------------
# Step 6: Create panel c plot (matching original style)
#----------------------------------------------
# Prepare data for plotting Panel C
plot_data_c <- results_rf %>%
  select(Model, R2_RF) %>%
  rename(R2 = R2_RF)

# Map model names to match original plots
model_labels <- c(
  "Species_richness" = "Species richness",
  "Species_identity" = "+ Species identity",
  "Effect_size" = "+ Effect size"
)

# Create Panel C plot (matching the original style)
Pc <- ggplot(data = plot_data_c, aes(x = Model, y = R2, fill = Model)) +
  # Violin plot
  geom_violin(color = "#00000000", alpha = 0.5, position = position_dodge(width = 0.3), trim = TRUE) +
  # Point range showing mean and CI
  geom_pointrange(data = rf_summary, 
                  aes(y = Mean, ymax = CI.high, ymin = CI.low, color = Model),
                  position = position_dodge(width = 0.2)) +
  # Use color scheme from original plot
  scale_fill_manual(values = c("#CAC5C6", "#DCD6E5", "#f0d99d")) +
  scale_color_manual(values = c("#000000", "#4E307D", "#7B3D12")) +
  # Apply formatting
  theme_bw() +
  mytheme +
  ylim(c(0, 1.0)) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(x = " ", y = "Variability explained (R²)", tag = "c") +
  scale_x_discrete(labels = model_labels) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1))
Pc

#----------------------------------------------
# Step 7: Train final models on complete dataset
#----------------------------------------------
# Using RandomForest for final models
final_model_Lv <- randomForest(fml_Lv, data = df.rf, 
                               ntree = n_tree, mtry = min(2, n_pred_Lv))

final_model_LvID <- randomForest(fml_LvID, data = df.rf, 
                                 ntree = n_tree, mtry = min(3, n_pred_LvID))

final_model_Full <- randomForest(fml_Full, data = df.rf, 
                                 ntree = n_tree, mtry = min(5, n_pred_Full))

# Calculate variable importance for the full model
var_importance <- importance(final_model_Full)
var_importance_df <- data.frame(
  Variable = rownames(var_importance),
  Importance = var_importance[,1],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Importance))

# Print variable importance
print("Variable importance in the full model:")
print(var_importance_df)


#----------------------------------------------
# Step 8: Save outputs
#----------------------------------------------
# Save Panel c to a file (uncomment to save)
# ggsave("panel_c.png", Pc, width = 5, height = 4, dpi = 300)


##################################################################################
#####                                                                        ##### 
#####                        Part5---Resules of Table S3 S4                  #####
#####                                                                        #####  
##################################################################################

#--------------------------------------
### Table S3 S4
#--------------------------------------
# Save results to CSV (uncomment to save)
# write.csv(rf_summary, "R2_summary.csv", row.names = FALSE)
# write.csv(var_importance_df, "variable_importance.csv", row.names = FALSE)

# Return results as a list
results_list <- list(
  summary_stats = rf_summary,
  variable_importance = var_importance_df,
  all_iterations = results_rf,
  panel_c_plot = Pc,
  final_models = list(
    Species_richness = final_model_Lv,
    Species_identity = final_model_LvID,
    Effect_size = final_model_Full
  )
)

##################################################################################
#####                                                                        ##### 
#####                     Part6---Visualization                              #####
#####                                                                        #####  
##################################################################################

# Create the combined plot layout
Pa / (Pb + Pc + plot_layout(widths = c(2, 1))) + plot_layout(heights = c(1, 1))

p <- wrap_plots(Pa, Pb, Pc, design = layout);p 
ggsave("~/Downloads/fig2.tiff", plot = p, width=3500, height=3000,units="px",dpi=300, compression = 'lzw',bg = "white")



