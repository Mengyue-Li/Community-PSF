
###########################################################################################
#####                                  Library & mytheme                              #####
###########################################################################################
# Loading the R packages
library(openxlsx)
library(readxl)
library(ape)
library(dplyr)
library(vegan)
library(ape)
library(picante)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(ggpmisc)
library(picante)


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

setwd('D:/2025.10.4-NPH/0-Code-alive-1216/code-0108-R1/Fig.S3')

#--------------------------------------
### Step 1: data 
#-------------------------------------- 

### conditioning plant community
Condition_com = read.xlsx("data_FigS3.xlsx", sheet = "Condition_AGraw", colNames = T, rowNames = T)
Condition_com[1:6,1:6]
rel_Condition_com <- decostand(Condition_com, method="total", MARGIN=1) ### Relative biomass
rowSums(rel_Condition_com)
rel_Condition_com [1:6,1:6]

### responding plant community
Responde_com = read.xlsx("data_FigS3.xlsx", sheet = "Respond_AG", colNames = T, rowNames = T)
Responde_com[1:6,1:6]


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
Respond_sp = unique(colnames(Responde_com)); Respond_sp

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



#--------------------------------------
### Step 3: lm model   
#--------------------------------------
data0 = read.xlsx("data_FigS3.xlsx", sheet = "species_specific_PSF", colNames = T, rowNames = T); data0[1:6,1:20]
data <- data0 %>% slice(-(1:5));data[1:6,1:20]

### species-specific PSF 
mod_Pc <-lm(Pc_PSF ~ Pc_Phylo_Dist, data = data)
anova(mod_Pc )-> mod_Pc_result;mod_Pc_result 

mod_Sv <- lm(Sv_PSF ~ Sv_Phylo_Dist, data = data)
anova(mod_Sv)-> mod_Sv_result;mod_Sv_result 

mod_Aa <- lm(Aa_PSF ~ Aa_Phylo_Dist, data = data)
anova(mod_Aa)-> mod_Aa_result;mod_Aa_result 

mod_Ah <- lm(Ah_PSF   ~ Ah_Phylo_Dist, data = data)
anova(mod_Ah)-> mod_Ah_result;mod_Ah_result 

mod_St <- lm(St_PSF  ~ St_Phylo_Dist, data = data)
anova(mod_St)-> mod_St_result;mod_St_result 

mod_Soc <-lm(Soc_PSF ~ Soc_Phylo_Dist, data = data)
anova(mod_Soc)-> mod_Soc_result;mod_Soc_result 


### all species were considered together,
data_ALL= read.xlsx("data_FigS3.xlsx", sheet = "all_species" )    
data_ALL[1:6,1:4]

mod <-  lm(Species_PSF ~ Species_Phylo_Dist, data = data_ALL)
anova(mod)-> mod_result;mod_result 


#--------------------------------------
### Step 4: FigS3 Visualization  
#--------------------------------------
P1= ggplot(data,aes(x=(Pc_Phylo_Dist),y=Pc_PSF))+
  #geom_smooth(aes(x = Pc_Phylo_Dist, y = Pc_PSF ), color = "#999999", method = 'lm', level = 0.95, se = FALSE, size = 1.6) +  
  geom_point(color = "#999999",size=1.7,show.legend = T)+ 
  stat_cor(method = "pearson", label.x.npc =0.0,size =4,label.y.npc = 0.1,alpha = 1) +   
  scale_size_manual(values=c(0.5, 0.5 ))+ 
  labs(x = "Weighted phylogenetic distance to the 
conditioning communities for responding species", y = "species-specific 
plant-soil feedback effect", tag = "(a)")+
  scale_x_continuous(limits = c(-0.01, 0.75), expand = c(0, 0), breaks = seq(0, 0.6, 0.2))+
  mytheme+mytheme 
P1

P2=ggplot(data,aes(x=(Sv_Phylo_Dist),y=Sv_PSF))+
  geom_smooth(aes(x = Sv_Phylo_Dist, y = Sv_PSF ), color ="#999999", method = 'lm', level = 0.95, se = FALSE, size = 1.6) +  
  geom_point(color ="#999999",size=1.7 ,show.legend = T)+
  stat_cor(method = "pearson", label.x.npc =0.0,size =4,label.y.npc = 0.1,alpha = 1) +  
  scale_size_manual(values=c(0.5, 0.5 ))+ 
  labs(x = "Weighted phylogenetic distance to the 
conditioning communities for responding species", y = "species-specific 
plant-soil feedback effect", tag = "(b)")+
  mytheme 
P2

P3=ggplot(data,aes(x=(Aa_Phylo_Dist),y=Aa_PSF))+
  #geom_smooth(aes(x = Aa_Phylo_Dist, y = Aa_PSF ), color = "#999999", method = 'lm', level = 0.95, se = FALSE, size = 1.6) +  
  geom_point(color = "#999999",size=1.7 ,show.legend = T)+
  stat_cor(method = "pearson", label.x.npc =0.0,size =4,label.y.npc = 0.1,alpha = 1) + 
  scale_size_manual(values=c(0.5, 0.5 ))+ 
  labs(x = "Weighted phylogenetic distance to the 
conditioning communities for responding species", y = "species-specific 
plant-soil feedback effect",tag = "(c)")+
  mytheme  
P3

P4=ggplot(data,aes(x=(Ah_Phylo_Dist),y=Ah_PSF))+
  #geom_smooth(aes(x = Ah_Phylo_Dist, y = Ah_PSF ), color = "#999999", method = 'lm', level = 0.95, se = FALSE, size = 1.6) +   
  geom_point(color = "#999999",size=1.7 ,show.legend = T)+
  stat_cor(method = "pearson", label.x.npc =0.0,size =4,label.y.npc = 0.1,alpha = 1) + 
  labs(x = "Weighted phylogenetic distance to the 
conditioning communities for responding species", y = "species-specific 
plant-soil feedback effect", tag = "(d)")+
  mytheme 
P4

P5=ggplot(data,aes(x=(St_Phylo_Dist),y=St_PSF))+
  geom_smooth(aes(x = St_Phylo_Dist, y = St_PSF ), color ="#999999", method = 'lm', level = 0.95, se = FALSE, size = 1.6) +   
  geom_point(color ="#999999",size=1.7,show.legend = T)+
  stat_cor(method = "pearson", label.x.npc =0.0,size =4,label.y.npc = ,alpha = 1) + 
  scale_size_manual(values=c(0.5, 0.5 ))+ 
  labs(x = "Weighted phylogenetic distance to the 
conditioning communities for responding species", y = "species-specific 
plant-soil feedback effect", tag = "(e)")+
  mytheme 
P5


P6=ggplot(data,aes(x=(Soc_Phylo_Dist),y=Soc_PSF))+
  #geom_smooth(aes(x = Soc_Phylo_Dist, y = Soc_PSF ), color ="#999999", method = 'lm', level = 0.95, se = FALSE, size = 1.6) +  
  geom_point(color ="#999999",size=1.7,show.legend = T)+
  stat_cor(method = "pearson", label.x.npc =0.0,size =4,label.y.npc = 0.1,alpha = 1) +  
  scale_size_manual(values=c(0.5, 0.5 ))+ 
  labs(x = "Weighted phylogenetic distance to the 
conditioning communities for responding species", y = "species-specific 
plant-soil feedback effect", tag = "(f)")+
  mytheme 
P6

(P1|P2)/(P3|P4)/(P5|P6)

