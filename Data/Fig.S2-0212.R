
###########################################################################################
#####                                    Fig.S1                                       #####
###########################################################################################

library(pheatmap)

setwd('C:/Users/MY/Desktop/CODE\Fig.S2')

data = read.xlsx("data_Fig.S2.xlsx", sheet = 1, rowNames = T, colNames = T)
head(data)
dim(data)
data= data.matrix(data) 

pheatmap(data,
         cluster_cols = F, cluster_rows = F, scale = "none", 
         treeheight_col = 0, treeheight_row = 0, 
         display_numbers = F, 
         legend = FALSE,
         #angle_col = 0,
         border_color = "#A2A2A2" , 
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 9,  fontsize_row = 9,
         color = colorRampPalette(c("#FFFFFF","grey"))(2),  
         cellwidth =5,  cellheight =5) -> pS1_con;pS1_con
