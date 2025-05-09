###########################################################################################
#####                                  Library &  theme                               #####
###########################################################################################

library(openxlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggforce)
library(RColorBrewer)

mytheme= theme(legend.position = "right",
               panel.grid=element_blank(),
               legend.background = element_rect(fill = NA),
               line = element_line(linewidth = 0.28), 
               axis.line = element_line(colour = "black",linewidth = 0.28),
               axis.ticks = element_line(colour = "black",linewidth = 0.28), 
               legend.title = element_text(colour='black', size=8), 
               legend.text = element_text(size = 7.5, lineheight = 1.33),  
               axis.title = element_text(colour='black', size=8),
               axis.text = element_text(colour='black',size=7),
               plot.tag = element_text(size = 9, face = "bold"))

###########################################################################################
#####                                       Fig.1                                     #####
###########################################################################################
setwd('C:/Users/MY/Desktop/li/Fig.1')                                                                              

data_li = read.xlsx("data_Fig.1.xlsx", sheet="bar" ,colNames = T, rowNames = F);head(data)

long_data <- data_li %>% pivot_longer(cols = -Group, names_to = "Year", values_to = "No._PSF_studies" ) %>% mutate(Year = as.numeric(Year))  
dim(long_data)
long_data$Group = factor(long_data$Group)

long_data$Group = factor(long_data$Group, levels = c("1","2","3","4","5","6","8","9","10","12",">12"))
# bar
ggplot(long_data, aes(x = Year, y = No._PSF_studies, alluvium = Group, fill = Group )) +  
  geom_col(position = "stack") + 
  scale_fill_manual(values = rev(c("#542D20","#b53d01" , "#D75C4D", "#fc8c3b",  "#EFBA55",
                                   "#C5E0D0", "#dbeef0",  "#7FC5D8","#3f95c6", "#2172b4","#1F466F")))+
  theme_bw() + mytheme + 
  theme(legend.position = "right") + 
  theme(legend.text = element_text(size = 14))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0,50)) +
  theme_test()+mytheme+labs(x ="Year", y = "Number of PSF studies")+
  scale_x_continuous( breaks = seq(min(long_data$Year), max(long_data$Year), by = 4))-> P;P



############### pie
data = read.xlsx("data_Fig.1.xlsx", sheet = "pie", rowNames = T, colNames = T);head(data)
total <- sum(data$No._study_Cspp._peryear)
data$start_angle <- NA
data$end_angle <- NA
#start_angle & end_angle
for (i in 1:nrow(data)) {
  if (i == 1) {
    data$start_angle[i] <- 0   
  } else {
    data$start_angle[i] <- data$end_angle[i - 1]  
  }
  data$end_angle[i] <- data$start_angle[i] + (data$No._study_Cspp._peryear[i] / total) * (2 * pi)
};head(data)


df1 = data;head(df1)
df1$MaxNo._Cspp.ID = as.factor(df1$MaxNo._Cspp.ID)
 
ggplot()+
  geom_arc_bar(data=df1,stat ="pie", aes(x0=0,y0=0,r0=0.9,r=2, amount=Percentage,
               fill=MaxNo._Cspp.ID,color=MaxNo._Cspp.ID ), show.legend = F)+
  coord_fixed() + theme_no_axes()+
  scale_fill_manual(values = c(  "#1F466F","#2172b4","#3f95c6","#7FC5D8" ,"#dbeef0" ,"#C5E0D0",
                                 "#EFBA55", "#fc8c3b","#D75C4D","#b53d01","#542D20"))+
  scale_color_manual(values = c( "#1F466F","#2172b4","#3f95c6","#7FC5D8","#dbeef0",  "#C5E0D0",
                                 "#EFBA55", "#fc8c3b","#D75C4D","#b53d01","#542D20"))+
  theme(panel.border = element_blank(),legend.key.width = unit(0.5,"cm"),legend.key.height = unit(0.5,"cm"))->pie;pie

#visualization-Fig.1
P + annotation_custom(grob=ggplotGrob(pie),ymin = 25, ymax =48, xmin=1992, xmax=2008) -> p1;p1
ggsave("Fig.1.pdf", plot = p1,width = 9, height = 7,  units = "cm", dpi = 600)  
 




