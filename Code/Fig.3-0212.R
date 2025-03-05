
###########################################################################################
#####                                  Library &  mytheme                             #####
###########################################################################################
# NOTE) ggepi & patchwork can be installed via devtools 
# devtools::install_github("lwjohnst86/ggepi")
# devtools::install_github("thomasp85/patchwork")
library(ggplot2)
library(ggepi)
library(ggridges)
library(patchwork)
library(party)
library(caret)
library(dplyr)
library(readxl)
library(glmmTMB)

mytheme= theme(legend.position = "none",
  panel.grid=element_blank(), 
  legend.title = element_blank(),
  legend.text = element_text(size = 14), 
  legend.background = element_rect(fill = NA),  
  axis.ticks = element_line(color='black'),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(colour='black', size=16),
  axis.title.y = element_text(colour='black', size=16),
  axis.text = element_text(colour='black',size=14),
  plot.tag = element_text(size = 18))

##################################################################################
#####                                                                        ##### 
#####           Part0---effect of home vs away:PSF                           #####
#####                                                                        #####  
##################################################################################

setwd('C:/Users/MY/Desktop/CODE/Fig.3')

#--------------------------------------
### Table S9
#--------------------------------------
dat <- read_excel("data_Fig.3.xlsx",sheet=1);head(dat)
# Fit lm model with totol biomass
mod_full <- lm(TG ~ Richness_con * (B_con + C_con + F_con + T_con + V_con), data = dat)
qqnorm(resid(mod_full));qqline(resid(mod_full));anova(mod_full)-> mod_full_result; mod_full_result
p <- mod_full_result$Pr;p 
p.adjust(p, "BH")
p.adjust(p, "bonferroni")

##################################################################################
#####                                                                        ##### 
#####                        Part1---Data info:PSF                           #####
#####                                                                        #####  
##################################################################################

df <- read_excel("data_Fig.3.xlsx",sheet=2)
df <- as.data.frame(df);rownames(df) <- df[,1];head(df)

df$remark = factor(df$remark, levels = unique(df$remark))
levels = c("1", "2", "4", "6", "12")
stressors = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T", "U","V","W","X","Y","Z","a","b")
treatment = as.vector(unique(df$remark))

responses = c("TG")
n_iter = 1000   # permutaion size


#####################################################################
##                      Functions-PSF                              ##
#####################################################################

#----------------------------------------------
# Estimating mean and its 95% confidence interval
#----------------------------------------------
BootStrap_mean = function(response, data=df, target = treatment, n_perm = n_iter){
  summary = list()
  
  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population = data[data$remark%in%stressors, response]
    if(treatment!="1") population = data[data$remark==treatment, response]
    size = length(population)
    
    for(id in c(1:n_perm)){
      k = mean(sample(population, size, replace = T))
      bs = append(bs, k)
    }
    summary[[treatment]] = c(quantile(bs, .025), mean(bs), quantile(bs, .975))
    names(summary[[treatment]]) = c("2.5%", "mean", "97.5%")
  }
  summary = t(data.frame(summary))
  summary = data.frame("target" = target, summary); row.names(summary) = c()
  return(summary)
}

#----------------------------------------------
# Estimating unstandardized effect size and its 95% confidence interval
#----------------------------------------------
BootStrap_ES_rep = function(response, data=df, target = treatment, n_perm = n_iter){
  resampled = list()
  
  population_CT = data[data$remark=="Control", response]
  
  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population_TR = data[data$remark%in%stressors, response]
    if(treatment!="1") population_TR = data[data$remark==treatment, response]
    size_CT = length(population_CT)
    size_TR = length(population_TR)
    
    for(id in c(1:n_perm)){
      k_CT = mean(sample(population_CT, size_CT, replace = T))
      k_TR = mean(sample(population_TR, size_TR, replace = T))
      bs = append(bs, k_TR - k_CT)
    }
    resampled[[treatment]] = bs
  }
  resampled[["Control"]] = rep(0, n_perm)
  return(resampled)
}

BootStrap_ES_summary = function(data){
  summary = list()
  p = 0
  summary[["Control"]] = c(0,0,0,1)
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

#-------------------------------------------
#  Estimating unstandardized ES of joint stressors and its 95% confidence interval
#-------------------------------------------
Null_distribution_rep = function(response, data=df, n_perm=n_iter){
  
  output = list()
  for(Lv in levels){
    resampled = list()
    
    # Checking which stressor combinations were jointly tested
    if(Lv=="1") combination = data[data$remark%in%stressors,c(2:29)]
    if(Lv!="1") combination = data[data["remark"]==Lv,c(2:29)]
    Level = sum(combination[1,])
    
    # Null distributions can be taken based on three different assumptions
    for(type in c("Additive", "Multiplicative", "Dominative")){
      
      population_CT = df[df$remark=="Control", response]
      size_CT = length(population_CT)
      
      # For each combination, bootstrap resampling is conducted
      for(j in c(1:nrow(combination))){
        bs = numeric(0)
        selected_stressors = stressors[which(combination[j,]==1)]
        sub_n_perm = ceiling(n_perm/nrow(combination))
        
        # bootstrap resampling
        for(id in c(1:sub_n_perm)){
          each_effect = numeric(0)
          k_CT = mean(sample(population_CT, size_CT, replace = T))
          
          for(treatment in selected_stressors){
            population_TR = df[df$remark==treatment, response]
            size_TR = length(population_TR)
            k_TR = mean(sample(population_TR, size_TR, replace = T))
            
            # ES estimate depending on the type of null hypotheses
            if(type=="Additive")       each_effect = append(each_effect, (k_TR - k_CT))
            if(type=="Multiplicative") each_effect = append(each_effect, (k_TR - k_CT)/k_CT)
            if(type=="Dominative")     each_effect = append(each_effect, (k_TR - k_CT))
          }
          
          # Calculating an expected ES after collecting the ESs of all relevant single stressors
          if(type=="Additive")       joint_effect = sum(each_effect)
          if(type=="Multiplicative"){
            z = 1
            for(m in c(1:Level)) z = z * (1 + each_effect[m])
            joint_effect = (z - 1)*k_CT
          }
          if(type=="Dominative")      joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
          
          bs = append(bs, joint_effect)
        }
        resampled[[type]][[j]] = bs
      }
      
    }
    output[[Lv]] = resampled
  }  
  return(output)
} 


Null_distribution_rep_transform = function(data){
  output = list()
  for(Lv in levels){
    for(type in c("Additive", "Multiplicative", "Dominative")){
      output[[Lv]][[type]] = sample(unlist(data[[Lv]][[type]]), n_iter, replace=F)
    }
  }
  return(output)
}


NHST_summary = function(null_data, Actual_data){
  output = list()
  for(Lv in levels){
    summary = list()
    summary[["Actual"]] = c(quantile(Actual_data[[Lv]], .025), mean(Actual_data[[Lv]]), quantile(Actual_data[[Lv]], .975), 1)
    p = 0
    assumptions = c("Additive", "Multiplicative", "Dominative")
    
    for(i_assumption in assumptions){
      bs   = (Actual_data[[Lv]] - null_data[[Lv]][[i_assumption]])
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      summary[[i_assumption]] = c(quantile(null_data[[Lv]][[i_assumption]], .025), mean(null_data[[Lv]][[i_assumption]]), quantile(null_data[[Lv]][[i_assumption]], .975), p)
    }
    summary = t(data.frame(summary))
    colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
    summary = data.frame(ES = c("Actual", "Additive","Multiplicative","Dominative"), summary); row.names(summary) = c()
    
    output[[Lv]] = summary
  }
  
  return(output)
}


NHST_summary_transform = function(data){
  output = list()
  for(i in 1:4){
    summary = rbind(data[["1"]][i, 2:4], data[["2"]][i, 2:4], data[["4"]][i, 2:4],
                    data[["6"]][i, 2:4], data[["12"]][i, 2:4])
    summary = cbind(levels, summary)
    colnames(summary) = c("Lv", "Low", "Mean", "High")
    output[[c("Actual", "Additive", "Multiplicative", "Dominative")[i]]] = summary
  }
  return(output)
}


Expected_ES_for_each = function(data){
  output = numeric(0)
  for(type in c("Additive", "Multiplicative", "Dominative")){
    tmp = numeric(0)
    for(Lv in levels){
      n_len = length(data[[Lv]][[type]])
      for(i in 1:n_len){
        tmp = append(tmp, mean(data[[Lv]][[type]][[i]]))
      }
    }
    output = cbind(output,tmp)
  }
  colnames(output)= c("E1", "E2", "E3")
  return(output)
}


#----------------------------------------------
# Step 1: Effect size estimate and null hypothesis significance testing
#----------------------------------------------
set.seed(123)
response_mean_all = list()
response_ES_all   = list()
joint_ES_null_all = list()
ES_for_each_all   = list()
g_rawdata_all     = list()
g_ms_all          = list()

for(i_response in responses){
  # Bootstrap estimate: single stressors
  response_mean   = BootStrap_mean(i_response)
  response_ES_bs  = BootStrap_ES_rep(i_response)
  response_ES     = BootStrap_ES_summary(response_ES_bs)
  
  # Bootstrap estimate: joint stressors
  Null_ES_bs0     = Null_distribution_rep(i_response)
  Null_ES_bs      = Null_distribution_rep_transform(Null_ES_bs0) 
  joint_ES_null   = NHST_summary(Null_ES_bs, response_ES_bs)
  ES_plot         = NHST_summary_transform(joint_ES_null)
  ES_for_each     = Expected_ES_for_each(Null_ES_bs0)
  
  
  #################################################################
  #------ Store the information for each response ----------------
  # Summary table combining the information from all response variables
  response_mean_all[[i_response]] = response_mean
  response_ES_all[[i_response]]   = response_ES
  joint_ES_null_all[[i_response]] = bind_rows(joint_ES_null, .id = "column_label")
  ES_for_each_all[[i_response]]   = ES_for_each
  #################################################################
  
  #### Fig2 A: Ploting the raw data with the mean & its 95% confidence intervals
  g_rawdata_all[[i_response]] =  local({
    i_response = i_response
    response_mean = response_mean
    
    ggplot() +
      theme_bw() + xlab(i_response) + 
      theme(legend.position = 'none', axis.title.x = element_blank()) +
      geom_rect(data = NULL,aes(xmin = -Inf, xmax = Inf, # Full horizontal range
                                ymin = which(levels(df$remark) == "1") - 0.5, # Start just before factor "1"
                                ymax = which(levels(df$remark) == "12") + 0.5), # End just after factor "12"
                fill = "#D3D3D3", alpha = 0.4, inherit.aes = FALSE) + 
      geom_vline(aes(xintercept = subset(response_mean, target == "Control")$mean), linetype = "dashed", color = "black") + 
      coord_flip() + 
      
      # CT as an open purple circle
      stat_density_ridges(data = subset(df, remark =="Control"), 
                          aes(x = .data[[i_response]], y = remark, fill = ..x..), 
                          geom = "density_ridges_gradient", bandwidth = 5.82,rel_min_height = 0.01, 
                          jittered_points = TRUE, color=alpha("#4E307D",0.5),fill=alpha("#4E307D",0.3),alpha = 0.4,  
                          position = position_points_jitter(height = .15, yoffset = .1 ),scale = 0.6) +
      geom_point(data = subset(response_mean, target == "Control"),aes(x = mean, y = target), 
                 shape = 1, size = 2 , stroke = 1.2,color = "#4E307D") +
      geom_errorbar(data = subset(response_mean, target =="Control"),
                    aes(xmin = X2.5., xmax = X97.5., y = target),width = 0, color = "#4E307D", linewidth = 0.5) +
      
      # Rest as purple filled circles
      stat_density_ridges(data = subset(df, !remark %in% c("Control", "1", "2", "4", "6", "12")),
                          aes(x = .data[[i_response]], y = remark, fill = ..x..), 
                          geom = "density_ridges_gradient", bandwidth = 5.02,rel_min_height = 0.01, 
                          jittered_points = TRUE, fill=alpha("#4E307D",0.2),color=alpha("#4E307D",0.2),alpha = 0.4,  
                          position = position_points_jitter(height = .15, yoffset = .1),scale = 0.6) +  
      geom_point(data = subset(response_mean, !target %in% c("Control", "1", "2", "4", "6", "12")),aes(x = mean, y = target), 
                 shape = 16,  size = 3,color = "#4E307D") +
      geom_errorbar(data = subset(response_mean, !target %in% c("Control", "1", "2", "4", "6", "12")),
                    aes(xmin = X2.5., xmax = X97.5., y = target),width = 0 , color = "#4E307D", linewidth = 0.5) + 
      
      # Target 1, 2, 4, 6, 12 as black filled circles
      stat_density_ridges(data = subset(df, remark %in% c("1", "2", "4", "6", "12")),
                          aes(x = .data[[i_response]], y = remark, fill = ..x..), 
                          geom = "density_ridges_gradient", bandwidth = 4.03,rel_min_height = 0.01, 
                          jittered_points = TRUE, color=alpha("#000000",0.2),fill=alpha("#000000",0.2),alpha = 0.1,  
                          position = position_points_jitter(height = .15, yoffset = .1),scale = 0.6) +
      geom_point(data = subset(response_mean, target %in% c("1", "2", "4", "6", "12")),aes(x = mean, y = target), 
                 shape = 16, size = 3,color = "#000000") +
      geom_errorbar(data = subset(response_mean, target %in% c("1", "2", "4", "6", "12")),
                    aes(xmin = X2.5., xmax = X97.5., y = target), width = 0 , color = "#000000", linewidth = 0.5)+
      labs(x ="Total biomass of responding community (g)", y = "Single plant species Richness", tag = "A")+
      mytheme+ scale_x_continuous(limits=c(25, 123), expand = c(0,0),breaks = seq(40, 120,20)) 
  })
  
  #### Fig2 C: Ploting the number of stressors and ES relationship
  g_ms_all[[i_response]] = local({
    i_response = i_response
    ct_value = response_mean[1,"mean"]
    ES_plot = ES_plot
    for(i in 1:4) ES_plot[[i]][,2:4] = ES_plot[[i]][,2:4] + ct_value
    
    ggplot()+
      theme_bw()+
      theme(legend.position = 'none', axis.title.x=element_blank(), axis.title.y=element_blank())+
      theme(legend.position = "none")  +
      coord_flip() +
      scale_y_discrete(limits = factor(c("1", "2", "4", "6", "12"), levels=c("1", "2", "4", "6", "12"))) +
      ### Mean & CI #### 3 assumptions
      geom_estci(data=ES_plot[["Additive"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#7B3D12", size=0.6, ci.linesize = 0.5, position=position_nudge(y = +0.1 )) +
      geom_estci(data=ES_plot[["Multiplicative"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#DE8125", size=0.6, ci.linesize = 0.5, position=position_nudge(y = +0.2)) +
      geom_estci(data=ES_plot[["Dominative"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#F7D8AD", size=0.6, ci.linesize = 0.5, position=position_nudge(y = +0.3)) +
      # Actual ES
      geom_estci(data=ES_plot[["Actual"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#000000", size=0.6, ci.linesize = 0.5, position=position_nudge(y = 0))+
      labs(x =" ", y = "Species richness", tag = "C")+ mytheme 
  })
}

#--------------------------------------
### Table S6 S7
#--------------------------------------
# write.table (response_mean_all, file ="PSF_mean estimate for each condition.csv",sep =",", quote =FALSE)
# write.table (response_ES_all, file ="PSF_Unstandardized effect size for each treatment.csv",sep =",", quote =FALSE)
# write.table (joint_ES_null_all, file ="PSF_null model.csv",sep =",", quote =FALSE)


#----------------------------------------------
# Step 2: Predictability comparison among random forest models
#----------------------------------------------
df.rf = df[df[, "remark"] %in% levels,]

lv_list  = unique(df.rf[,"Lv"])
id_lv1   = which(df.rf[,"Lv"]==1)
id_lvh   = which(df.rf[,"Lv"]>1)

n_data   = nrow(df.rf)
n_lv1    = sum(df.rf[,"Lv"]==10)
n_lvh    = sum(df.rf[,"Lv"]!=1)
n_eachlv = 10
n_tree   = 1000
n_iter2  = 1000

rf.r2 = data.frame(matrix(NA, ncol=3, nrow=3*n_iter2*length(responses)))
rf.r2[,1] = rep(responses,each=3*n_iter2)
rf.r2[,2] = rep(c("Lv", "Lv+ID", "All"), n_iter2*length(responses))
rf.r2[,2] = factor(rf.r2[,2], levels=c("Lv", "Lv+ID", "All"))
colnames(rf.r2) = c("Response", "Model", "R2")

rf.prediction.all = list()

j = 0

for (i_response in responses) {
  
  df.rf.tmp = cbind(df.rf, ES_for_each_all[[i_response]])
  
  # Define formulas
  eval(parse(text = paste("fml = formula(", i_response, " ~ ",paste(c("Lv", stressors, colnames(ES_for_each_all[[i_response]])), collapse = " + "), ")")))
  eval(parse(text = paste("fml.Lv = formula(", i_response, " ~ Lv)")))
  eval(parse(text = paste("fml.LvID = formula(", i_response, " ~ ", 
                          paste(c("Lv", stressors), collapse = " + "), ")")))
  
  # Calculate number of predictors for dynamic mtry
  ninputs_fml = length(all.vars(fml)) - 1
  ninputs_fml_Lv = length(all.vars(fml.Lv)) - 1
  ninputs_fml_LvID = length(all.vars(fml.LvID)) - 1
  
  for (i in 1:n_iter2) {
    j = j + 1
   
    # Bootstrap resampling
    rid = c(sample(id_lv1, n_eachlv, replace = TRUE), sample(id_lvh, n_lvh, replace = TRUE))
    rdf = df.rf.tmp[rid, ]
    
    # Train models with dynamic mtry
    rf_model = tryCatch({
      cforest(fml, data = rdf, control = cforest_control(
        ntree = n_tree, minsplit = 5, minbucket = 2, mtry = min(5, ninputs_fml)
      ))
    }, error = function(e) {
      cforest(fml.Lv, data = rdf, control = cforest_control(
        ntree = n_tree, minsplit = 5, minbucket = 2, mtry = min(5, ninputs_fml_Lv)
      ))
    })
    
    rf_model.Lv = cforest(fml.Lv, data = rdf, control = cforest_control(
      ntree = n_tree, minsplit = 5, minbucket = 2, mtry = min(5, ninputs_fml_Lv)
    ))
    
    rf_model.LvID = cforest(fml.LvID, data = rdf, control = cforest_control(
      ntree = n_tree, minsplit = 5, minbucket = 2, mtry = min(5, ninputs_fml_LvID)
    ))
    
    # Evaluating fitting performance
    rf_prdct.Lv = tryCatch({predict(rf_model.Lv, OOB = TRUE)}, error = function(e) {rep(0, length(rid))})
    rf_prdct.LvID = tryCatch({predict(rf_model.LvID, OOB = TRUE)}, error = function(e) {rep(0, length(rid))})
    rf_prdct.All = tryCatch({predict(rf_model, OOB = TRUE)}, error = function(e) {rep(0, length(rid))})
    
    rf.r2[3 * j - 2, 3] = postResample(rf_prdct.Lv, rdf[, i_response])[2]
    rf.r2[3 * j - 1, 3] = postResample(rf_prdct.LvID, rdf[, i_response])[2]
    rf.r2[3 * j, 3] = postResample(rf_prdct.All, rdf[, i_response])[2]
  }
  
  # Train final models on the full dataset
  rf_model = cforest(fml, data = df.rf.tmp, control = cforest_control(
    ntree = n_tree, minsplit = 5, minbucket = 2, mtry = min(5, ninputs_fml)
  ))
  
  rf_model.Lv = cforest(fml.Lv, data = df.rf.tmp, control = cforest_control(
    ntree = n_tree, minsplit = 5, minbucket = 2, mtry = min(5, ninputs_fml_Lv)
  ))
  
  rf_model.LvID = cforest(fml.LvID, data = df.rf.tmp, control = cforest_control(
    ntree = n_tree, minsplit = 5, minbucket = 2, mtry = min(5, ninputs_fml_LvID)
  ))
  
  # Store predictions
  rf.prediction.all[[i_response]] = data.frame(
    Model = rep(c("Lv", "LvID", "All"), each = nrow(df.rf.tmp)),
    Predicted = c(
      predict(rf_model.Lv, OOB = FALSE), 
      predict(rf_model.LvID, OOB = FALSE), 
      predict(rf_model, OOB = FALSE)
    ),
    Observed = rep(df.rf.tmp[, i_response], 3)
  )
}

rf.r2.summary = data.frame(matrix(NA,ncol=5,nrow=3*length(responses)))
colnames(rf.r2.summary) = c("Response","Model", "CI.low", "Mean", "CI.high")
rf.r2.summary[,1] = rep(responses,each=3)
rf.r2.summary[,2] = rep(c("Lv", "Lv+ID", "All"),length(responses))
rf.r2.summary[,2] = factor(rf.r2.summary[,2],levels=c("Lv", "Lv+ID", "All"))

j  = 0
for(i_response in responses){
  jj = 0
  jjj = 0
  rf.r2.summary[3*j+1,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Lv"),   3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[3*j+2,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Lv+ID"),3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[3*j+3,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="All"),  3],c(.025,.50,.975), na.rm=T)
}

# write.csv(rf.r2.summary, "PSF_randomforest_r2_summary.csv")

df.rf <- df.rf[,-c(37,38)]
g_r2_all = list()

for(i in 1:(length(responses))){
  g_r2_all[[i]] = local({
    i = i
    rf.r2 = rf.r2
    rf.r2.summary = rf.r2.summary 
    
    ggplot(data=rf.r2[rf.r2[,1]==responses[i],],aes(x=Model,y=R2, fill=Model))+
      geom_violin( color="#00000000",alpha=.5,position=position_dodge(width=0.3),trim=T)+ 
      geom_pointrange(data=rf.r2.summary[rf.r2.summary[,1]==responses[i],], 
                      aes(y=Mean, ymax=CI.high, ymin=CI.low,color=Model),
                      position=position_dodge(width=0.2)) +
      scale_fill_manual(values  = c( "#CAC5C6" , "#DCD6E5"  ,  "#F6E8C3" ))+ 
      scale_color_manual(values = c("#000000", "#4E307D", "#7B3D12"))+ 
      theme_bw() +mytheme+ ylim(c(0,1.0))+
      theme(legend.position = 'none', axis.title.x=element_blank())+
      labs(x =" ", y = "Variability explained (R2%)", tag = "E")+
      scale_x_discrete(labels = c("Lv" = "Species richness", "Lv+ID" = "+ Species identity","All" = "+ Effect size" ))+
      theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1))
  })
}


#-------------------------------------------
#  Step 3: Visualization-PSF
#-------------------------------------------
cowplot::plot_grid(g_rawdata_all[[1]])-> Fig_2_A; Fig_2_A
cowplot::plot_grid(g_ms_all[[1]],g_r2_all[[1]], align = "h", ncol = 2,rel_widths= c(2, 1.5))-> Fig_2_BC; Fig_2_BC 
Fig_2_1/ Fig_2_2
