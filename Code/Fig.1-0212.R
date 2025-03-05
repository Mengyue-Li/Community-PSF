###########################################################################################
#####                                  Library &  theme                               #####
###########################################################################################

library(openxlsx)
library(ggfortify)
library(ggforce)
library(dplyr)
library(ggrepel)
library(patchwork)
library(ggpubr)
library(export)
library(ggspatial)
library(ggplot2)
library(sf)
library(cowplot)
library(grid)#unit:cmlibrary(raster)

mytheme= theme(legend.position = "none",
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


#####################################################################
#####                        Fig.1B                             #####
#####################################################################

setwd('C:/Users/梦月/Desktop/data-code-li-0113/Fig.1')

data = read.xlsx("data_Fig.1.xlsx", sheet = "pie-460exp", rowNames = T, colNames = T)
head(data)
total <- sum(data$No._literatures_sp)
data$start_angle <- NA
data$end_angle <- NA
#start_angle & end_angle
for (i in 1:nrow(data)) {
  if (i == 1) {
    data$start_angle[i] <- 0   
  } else {
    data$start_angle[i] <- data$end_angle[i - 1]  
  }
 data$end_angle[i] <- data$start_angle[i] + (data$No._literatures_sp[i] / total) * (2 * pi)
}
head(data)

#pie
df1 = data
head(df1)
df1$MaxNo._CspID = as.factor(df1$MaxNo._CspID)

ggplot()+
  geom_arc_bar(data=df1,stat ="pie", aes(x0=0,y0=0,r0=1,r=2, amount=Percentage,
                   fill=MaxNo._CspID,color=MaxNo._CspID,explode=Focus), show.legend = F)+
  geom_arc(data=subset(df1, Group=="1"), size=1,color="#3D3739", 
           aes(x0=0.08,y0=-0.2,r=2,start=start_angle, end=end_angle))+ 
  geom_arc( data= subset(df1, Group =="mid"),size=1,color="#7C777B",
            aes(x0=0,y0=0,r=2.2,start=start_angle, end=end_angle ))+
  geom_arc(data= subset(df1, Group=="high"),size=1,color="#ACACAC",
           aes(x0=0,y0=0,r=2.3,size=index, start=start_angle, end=end_angle))+
  coord_fixed() + theme_no_axes()+
  scale_fill_manual(values = c( "#FCA300","#7E2968","#FDE725","#226312","#5E6D73",
                               "#00204D","#0433FF","#2A83D7","#409FF8","#141414","#636564"))+
  scale_color_manual(values = c("#FCA300","#7E2968","#FDE725","#226312","#5E6D73",
                                "#00204D","#0433FF","#2A83D7","#409FF8","#141414","#636564"))+
  theme(panel.border = element_blank(),legend.key.width = unit(0.5,"cm"),legend.key.height = unit(0.5,"cm"))->pie;pie

#line
df2 = read.xlsx("data_Fig.1.xlsx", sheet = "line-460exp", rowNames = T, colNames = T)
head(df2)
df2$MaxNo._Csp = as.factor(df2$MaxNo._Csp)

ggplot(df2,aes(x=Year,y= No._literatures_sp_peryear,color=MaxNo._CspID,fill=MaxNo._CspID))+
  geom_line(aes(y=No._literatures_sp_peryear), linewidth=1.5)+
  scale_fill_manual(values = c("1"="#FCA300","2"="#7E2968","3"="#FDE725","4"="#226312","5"="#5E6D73",
                               "6"="#00204D","8"= "#0433FF","9"="#2A83D7","10"="#409FF8","12"="#141414",">12"="#636564"))+
  scale_color_manual(values = c("1"="#FCA300","2"="#7E2968","3"="#FDE725","4"="#226312","5"="#5E6D73",
                                "6"="#00204D","8"= "#0433FF","9"="#2A83D7","10"="#409FF8","12"="#141414",">12"="#636564"))+
  scale_x_continuous(limits=c(1993, 2024),breaks = seq(1993, 2024,3) )+
  scale_y_continuous(  limits=c(0, 60), expand = c(0, 0))+
  theme_test()+mytheme+labs(x ="Year", y = "Number of literatures", tag = "B", face = "bold") -> plin;plin

#visualization-Fig.1B
plin + annotation_custom(grob=ggplotGrob(pie),ymin = 25, ymax =58, xmin=1992, xmax=2008) -> p1B;p1B




#####################################################################
#####                        Fig.1C                             #####
#####################################################################
df_map <- read_sf("China_map.json")
site = read.xlsx("data_Fig.1.xlsx", sheet = "map-1977plots", rowNames = F, colNames = T)
site <- site %>%group_by(Longitude, Latitude) %>% mutate(size_value = n()) %>%ungroup()
china_map <- read_sf("China_map.json")
east_asia_map <- read_sf("custom.geo.json")
class(east_asia_map)
east_asia_countries <- c("China", "Japan", "South Korea", "North Korea", "Mongolia", "Taiwan")

ggplot(data = east_asia_map) +    
  geom_sf(data = east_asia_map[east_asia_map$sovereignt %in% east_asia_countries,],color = "black", size = 0.5, fill = NA) + 
  geom_sf(data = china_map, color = "black", size = 0.5, fill = NA) +
  geom_sf(data = china_map[c(3,12,9,10,14,15:20),], color = "black", size = 0.5, fill = "#DCDCDC") + 
  theme_test()+ 
  theme(text = element_blank(), 
        legend.position = "none",  
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),  
        axis.ticks = element_blank(),  
        axis.line = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.3), 
        panel.background = element_rect(fill = NA, color = NA), 
        plot.title = element_blank())+
  annotate("text", x = 85, y = 10, label = "East Asia", size = 6.5, color = "black")-> p1;p1
 
province_labels <- data.frame(name = c("Hebei", "Shandong", "Anhui", "Henan", "Hubei", "Jiangsu", "Jiangxi", "Hunan", "Guangdong", "Guangxi"),  
  Longitude = c(114.5, 117.0, 118, 113, 114, 119, 116, 113, 113, 108),  
  Latitude = c(38.04, 36.6, 32, 34, 31, 33, 29, 28, 23, 23)) 

ggplot() +
  geom_sf(data = china_map[c(3,12,9,10,14,15:20),],size=3.5,fill="#EBEBEB", color = "black")+
  geom_text(data = province_labels, aes(x = Longitude, y = Latitude, label = name), size=4,color = "black") + 
  geom_point(data = site, mapping = aes(x=Longitude,y=Latitude, fill=factor(Color),size= size_value ),pch = 21,alpha=1, color = "black") + 
  scale_color_manual(values=c("#8282B6"))+ scale_fill_manual(values=c("#5A8C20"))+
  theme_test()+ 
  scale_x_continuous(limits = c(105,121),breaks = c(105,110,115,120)) + 
  scale_y_continuous(limits = c(20,40),breaks = c(20,25,30,35,40)) + 
  annotation_scale(location = "bl", style = "ticks",line_width = 0.1,pad_y = unit(0.5, "cm"),text_cex = 1) +############
  theme(legend.position = c(0.9,0.2),plot.margin = margin(0.1,0,0,0, unit = "cm"))+labs(x=NULL, y=NULL, tag = "C")-> p2;p2  
                    
p2 + annotation_custom(grob=ggplotGrob(p1),ymin = 28, ymax =48.561,xmin= 103.95, xmax=113)-> p1C;p1C




#####################################################################
#####                        Fig.1D                             #####
#####################################################################
df4 = read.xlsx("data_Fig.1.xlsx", sheet = "richness-1977plots", rowNames = F, colNames = T)

ggplot(df4, aes(x = Richness, y = Frequency_plot))+
  geom_bar(stat = "identity",position=position_dodge(width = 0.8),fill="#5A8C20",size=0.2,width = 0.8,alpha=1)+
  geom_text(aes(label = Frequency_plot), position = position_dodge(width = 0.8),vjust = -0.45,size = 6.5)+
  scale_color_manual(values=c("black"))+scale_fill_manual(values=c("#8282B6"))+
  scale_x_continuous(limits=c(0, 14), expand = c(0,0),breaks = seq(1, 13, 1)) +
  scale_y_continuous(limits=c(0, 490), expand = c(0,0),breaks = seq(0, 460,100)) +
  theme_test()+ mytheme+ 
  labs(x = expression("Plant richness (# species)"),y = "Frequency (# plot)", tag = "D")-> p1D;p1D
 
