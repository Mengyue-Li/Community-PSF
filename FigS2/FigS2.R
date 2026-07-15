###########################################################################################
#####                                  Library & mytheme                              #####
###########################################################################################
# install.packages("ggpmisc")
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggpmisc)
library(patchwork)
library(FSA)         
library(multcompView) 

mytheme = theme( panel.background = element_rect(fill='white', colour='black'),legend.position = "none",
                 panel.grid=element_blank(), legend.title = element_blank(),legend.text = element_text(size = 11),
                 legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                 axis.ticks = element_line(color='black'),
                 axis.line = element_line(colour = "black"), 
                 axis.title.x=element_text(colour='black', size=13,vjust = 1),
                 axis.title.y=element_text(colour='black', size=13,vjust = 1),
                 axis.text=element_text(colour='black',size=11),
                 plot.tag = element_text(size = 14, face = "bold"))

###########################################################################################
#####                                  main part                                      #####
###########################################################################################

setwd('D:/2025.10.4-NPH/1-0330/2-NEWrawdata-code/FigS2') 

data <- read.xlsx("FigS2.xlsx", sheet="0715"); data[1:5,1:5]

#Extract monoculture / mixture rows and keep only Species_con and AG_ind_con columns
mono_raw <- data %>% filter(Group == "mono") %>% dplyr::select(Species_con, AG_ind_con); mono_raw 
mix_raw <- data %>% filter(Group == "mix")  %>% dplyr::select(Species_con, AG_ind_con); mix_raw 

# Define a function to perform ANOVA + Tukey HSD and return compact letter display
get_letters <- function(df) { aov_res <- aov(AG_ind_con ~ Species_con, data = df)
hsd <- HSD.test(aov_res, "Species_con", group = TRUE) 
letters_df <- data.frame(Species_con = rownames(hsd$groups),letters = as.character(hsd$groups$groups))
              return(letters_df)
}

# Apply the function to monoculture / mixture data
letters_mono <- get_letters(mono_raw); letters_mono
letters_mix  <- get_letters(mix_raw); letters_mix

# Summarise raw data
data_se <- data %>% group_by(Group, Species_con) %>% summarise( mean_AG_ind_con  = mean(AG_ind_con),
                                                                se_AG_ind_con = sd(AG_ind_con) / sqrt(n()),.groups = "drop") ; data_se 
# Join the letters
datamono <- data_se %>% filter(Group == "mono") %>% left_join(letters_mono, by = "Species_con"); datamono
datamix  <- data_se %>% filter(Group == "mix")  %>% left_join(letters_mix, by = "Species_con"); datamix


species_order <- c("Sor","Aa","Sv","Ul","Vo","St","Cal","Sn","Ai","Pf","Pb","Md","Cab","Mc",
                   "Car","Sa","Pa",'Ss',"Sc","Ah","Ds","Soc","Da","Av","At","Cp","Cc","Cs")

#--------------------------------------
### Fig S2a 
#--------------------------------------
pmono <-
  ggplot(datamono , aes(x = Species_con, y = mean_AG_ind_con) ) +
  geom_bar(stat="identity", position="dodge", width = 0.8, color = "#999999", fill = "#999999") +
  geom_errorbar(aes(ymin = mean_AG_ind_con - se_AG_ind_con, ymax = mean_AG_ind_con  + se_AG_ind_con),
                position = position_dodge(0.8), width = 0.2, color = "#333333" ) +
  geom_text(aes(label = letters, y = mean_AG_ind_con + se_AG_ind_con + 0.15),vjust = 0, size = 3.5) +
  scale_y_continuous(limits = c(0, 3.1), expand = c(0, 0), breaks = seq(0, 3, 1)) +
  scale_x_discrete(limits = species_order) +
  labs(y = "Mean aboveground biomass\nof individuals in monocultures (g)", x = "Conditioning plant species", tag = "(a)") +
  mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); pmono

#--------------------------------------
### Fig S2b
#--------------------------------------
pmix <-ggplot(datamix, aes(x = Species_con, y = mean_AG_ind_con)) +
  geom_bar(stat="identity", position="dodge", width = 0.8, color = "#999999", fill = "#999999") +
  geom_errorbar(aes(ymin = mean_AG_ind_con - se_AG_ind_con, ymax = mean_AG_ind_con + se_AG_ind_con),
                position = position_dodge(0.8), width = 0.2, color = "#333333") +
  geom_text(aes(label = letters, y = mean_AG_ind_con + se_AG_ind_con + 0.15), vjust = 0, size = 3.5) +
  scale_y_continuous(limits = c(0, 3.2), expand = c(0, 0), breaks = seq(0, 3, 1)) +
  scale_x_discrete(limits = species_order) +
  labs(y = "Mean aboveground biomass \nof individuals in mixtures (g)", x = "Conditioning plant species", tag = "(b)") +
  mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); pmix

pmono / pmix


#--------------------------------------
### Fig S2c
#--------------------------------------
data_se
df_wide <- data_se %>% pivot_wider(id_cols = Species_con, names_from = Group, values_from = c(mean_AG_ind_con, se_AG_ind_con))
df_wide

mod = lm(mean_AG_ind_con_mono~mean_AG_ind_con_mix,data=df_wide)
anova(mod)-> mod_result; mod_result
p <- mod_result$Pr;p 
p.adjust(p, "BH")

p <-
  ggplot(df_wide, aes(x = mean_AG_ind_con_mono, y = mean_AG_ind_con_mix)) +
  geom_smooth(method = 'lm', se = FALSE, color = "black") +
  geom_errorbarh(aes(xmin = mean_AG_ind_con_mono - se_AG_ind_con_mono,xmax = mean_AG_ind_con_mono + se_AG_ind_con_mono),
                 color = "#999999",height = 0, alpha = 1,size = 0.5 ) +
  geom_errorbar(aes(ymin = mean_AG_ind_con_mix - se_AG_ind_con_mix,ymax = mean_AG_ind_con_mix + se_AG_ind_con_mix ),
                color = "#999999",width = 0, alpha = 1,size = 0.5) +
  geom_point(fill = "#999999", size = 2.2, color = "#333333", pch = 21) +
  stat_cor(method = "pearson",  label.x.npc = 0.6,label.y.npc = 1,size = 4, r.digits = 3, p.digits = 3) +
  labs(x = "Average aboveground biomass of individual 
conditioning species plants grown in monocultures(g)", y = "Average aboveground biomass of individual 
conditioning species plants grown in mixtures (g)", tag = "(c)") +
  #stat_poly_eq(method = 'lm', aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep = "*\", \"*")),size=4 ,alpha = 0.8)+
  mytheme 
p


#--------------------------------------
### Fig S2 Visualization
#--------------------------------------
((pmono / pmix) | p) + plot_layout(widths = c(2, 1.2))

