
###########################################################################################
#####                                    Fig.S1                                       #####
###########################################################################################

library(pheatmap)
setwd('C:/Users/MY/Desktop/CODE/Fig.S1&2')

data1 = read.xlsx("data_Fig.S1&2.xlsx", sheet = "FigS1", rowNames = T, colNames = T)
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
         fontsize_col = 9,  fontsize_row = 9,
         color = colorRampPalette(c( "#FFFFFF" , "grey"))(2),  
         cellwidth =5,  cellheight = 5) -> pS1;pS1


###########################################################################################
#####                                    Fig.S2                                       #####
###########################################################################################
data2 = read.xlsx("data_Fig.S1&2.xlsx", sheet = "FigS2", rowNames = T, colNames = T)
head(data2)
dim(data2)
data= data.matrix(data2) 

pheatmap(data2,
         cluster_cols = F, cluster_rows = F, scale = "none", 
         treeheight_col = 0, treeheight_row = 0, 
         display_numbers = F, 
         legend = FALSE,
         #angle_col = 0,
         border_color ="#A2A2A2", 
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 12,  fontsize_row = 12,
         color = colorRampPalette(c( "#FFFFFF" , "#d1d1d1"))(2),  
         cellwidth =12,  cellheight = 5) -> pS2;pS2


