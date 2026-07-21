
###########################################################################################
#####                                    Fig.S1                                       #####
###########################################################################################
library(openxlsx)
library(pheatmap)

setwd('D:/NPH/1-Revised manuscript 2026/Fig.S1')

data1 = read.xlsx("data-Fig.S1.xlsx", sheet = "FigS1", rowNames = T, colNames = T)
head(data1)
dim(data1)
data= data.matrix(data1) 

pheatmap(data1,
         cluster_cols = F, cluster_rows = F, scale = "none", 
         treeheight_col = 0, treeheight_row = 0, 
         display_numbers = F, 
         legend = FALSE,
         #angle_col = 0,
         border_color ="#A2A2A2", 
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 8,  fontsize_row = 8,
         color = colorRampPalette(c( "#FFFFFF" , "grey"))(2),  
         cellwidth =5,  cellheight = 5 ) -> pS1;pS1

 
