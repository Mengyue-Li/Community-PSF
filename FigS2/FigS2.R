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

setwd('D:/2025.10.4-NPH/0-Code-alive-1216/FigS2-1222') 

data <- read.xlsx("data_FigS2.xlsx",sheet="AG_con"); data[1:5,1:5]
data_se <- data %>% group_by(Group, Species_con) %>% summarise( mean_AG_ind_con  = mean(AG_ind_con),
                                                                se_AG_ind_con = sd(AG_ind_con) / sqrt(n()),
                                                                .groups = "drop" )  

#--------------------------------------
### Fig S2a 
#--------------------------------------
#datamono <- data [data$Group == "mono", ];  datamono[1:5,1:5]
#datamono_se <- data %>% group_by(Species_con) %>% summarise(mean_AG_ind_con= mean(AG_ind_con),  se_AG_ind_con = sd(AG_ind_con) / sqrt(n()))

datamono_se <- data_se [data_se$Group == "mono", ];  datamono_se[1:4,1:4]

pmono <-
  ggplot(datamono_se, aes(x = Species_con, y = mean_AG_ind_con) ) +
  geom_bar(stat="identity",position="dodge",width = 0.8, color = "#999999",fill= "#999999") +
  # geom_col(position = position_dodge(0.8),, color = "#999999",fill= "#999999") +
  geom_errorbar(aes(ymin = mean_AG_ind_con - se_AG_ind_con, ymax = mean_AG_ind_con  + se_AG_ind_con),
                position = position_dodge(0.8), width = 0.2, color = "#333333" ) +
  scale_y_continuous(limits = c(0, 3.1),expand = c(0, 0), breaks = seq(0, 3, 1))+
  scale_x_discrete(limits = c(c("Sor","Aa","Sv","Ul","Vo","St","Cal","Sn","Ai","Pf","Pb","Md","Cab","Mc",
                                "Car","Sa","Pa",'Ss',"Sc","Ah","Ds","Soc","Da","Av","At","Cp","Cc","Cs")))+ 
  labs(y = "Average aboveground biomass 
of individuals in monocultures (g)", x = " ", tag = "(a)")+
  mytheme+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pmono


#--------------------------------------
### Fig S2b
#--------------------------------------
#datamix <-data [data$Group == "mono", ]; datamix[1:5,1:5]
#datamix_se <- data %>% group_by(Species_con) %>% summarise(mean_AG_ind_con= mean(AG_ind_con),  se_AG_ind_con = sd(AG_ind_con) / sqrt(n()))

datamix_se <- data_se [data_se$Group == "mix", ];  datamono_se[1:4,1:4]
 
pmix <- 
  ggplot(datamix_se, aes(x = Species_con,   y = mean_AG_ind_con)) +
  geom_bar(stat="identity",position="dodge",width = 0.8, color = "#999999",fill= "#999999") +
  #geom_col(position = position_dodge(width = 1 ), color = "#999999",fill= "#999999") +
  geom_errorbar(aes(ymin = mean_AG_ind_con - se_AG_ind_con , ymax = mean_AG_ind_con  + se_AG_ind_con),
                position = position_dodge(0.8), width = 0.2, color = "#333333") +
  scale_y_continuous(limits = c(0, 3.2), expand = c(0, 0), breaks = seq(0, 3, 1))+
  scale_x_discrete(limits = c(c("Sor","Aa","Sv","Ul","Vo","St","Cal","Sn","Ai","Pf","Pb","Md","Cab","Mc",
                                "Car","Sa","Pa",'Ss',"Sc","Ah","Ds","Soc","Da","Av","At","Cp","Cc","Cs")))+ 
 labs(y = "Average aboveground biomass 
of individuals in mixtures (g)", x = "Conditioning plant species", tag = "(b)")+
  mytheme+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pmix
 

#--------------------------------------
### Fig S2c
#--------------------------------------
data_se
df_wide <- data_se %>% pivot_wider(id_cols = Species_con, names_from = Group,
                                 values_from = c(mean_AG_ind_con, se_AG_ind_con))
df_wide
 
p <-
  ggplot(df_wide, aes(x = mean_AG_ind_con_mono, y = mean_AG_ind_con_mix)) +
  geom_smooth(method = 'lm', se = FALSE, color = "black") +
  geom_errorbarh(aes(xmin = mean_AG_ind_con_mono - se_AG_ind_con_mono,
                     xmax = mean_AG_ind_con_mono + se_AG_ind_con_mono),
                    color = "#999999",height = 0, alpha = 1,size = 0.5 ) +
   geom_errorbar(aes(ymin = mean_AG_ind_con_mix - se_AG_ind_con_mix,
                     ymax = mean_AG_ind_con_mix + se_AG_ind_con_mix ),
                    color = "#999999",width = 0, alpha = 1,size = 0.5) +
  geom_point(fill = "#999999", size = 2.2, color = "#333333", pch = 21) +
  stat_cor(method = "pearson",  label.x.npc = 0.6,label.y.npc = 1,size = 4, r.digits = 3, p.digits = 3) +
  labs(x = "Average aboveground biomass of individual 
conditioning species plants grown in monocultures and mixtures (g)", y = "Average aboveground biomass of individual 
conditioning species plants grown in mixtures (g)", tag = "(c)") +
  #stat_poly_eq(method = 'lm', aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep = "*\", \"*")),size=4 ,alpha = 0.8)+
  mytheme 
p


#--------------------------------------
### Fig S2 Visualization
#--------------------------------------
(pmono / pmix) | p



p <-
  ggplot(df_wide, aes(x = mean_AG_ind_con_mono, y = mean_AG_ind_con_mix)) +
  geom_smooth(method = 'lm', se = FALSE, color = "black") +
     geom_point(fill = "#999999", size = 2.2, color = "#333333", pch = 21) +
  stat_cor(method = "pearson",  label.x.npc = 0.0,label.y.npc = 0.1,size = 4) +
  labs(x = "Average aboveground biomass of individual 
conditioning species plants grown in monocultures and mixtures (g)", y = "Average aboveground biomass of individual 
conditioning species plants grown in mixtures (g)", tag = "(c)")+
  stat_poly_eq(method = 'lm', aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep = "*\", \"*")),size=4 ,alpha = 0.8)+
  
  mytheme 

p
