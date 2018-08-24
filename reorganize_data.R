
library("GenVisR");library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2)
set.seed(426)



x<-read.csv(file = "~/Downloads/7WEX.output.csv", header = T)
table(x$sample)
table(x$variant_class)
samples_to_plot <- c("TCRBOA1_N_WEX", "TCRBOA1_T_WEX", "TCRBOA2_N_WEX", "TCRBOA2_T_WEX", "TCRBOA3_N_WEX", "TCRBOA3_T_WEX", "TCRBOA4_N_WEX", "TCRBOA4_T_WEX", 
                     "TCRBOA5_N_WEX", "TCRBOA5_T_WEX", "TCRBOA6_N_WEX", "TCRBOA6_T_WEX", "TCRBOA7_N_WEX", "TCRBOA7_T_WEX")

pdf(file = "~/Downloads/7WEX.output_reshuffled.pdf",width = 10,height = 15)
waterfall(x = x,fileType = "Custom",variant_class_order = "nonsynonymous SNV",
          plotSamples = samples_to_plot,
          sampOrder = samples_to_plot,
          mainXlabel = T,plotMutBurden = F) 
dev.off()
# x$sample<-as.character(x$sample)
# x$Stage <- ifelse(grepl("_T_", x$sample), "Tumor", "Normal")
# x<-arrange(x,Stage)%>%
#         separate(sample,into="samp",sep = "_")
# 
# #create two new tables; a normal and a tumor table
# normal<-filter(x,Stage=="Normal")
# tumor<-filter(x,Stage=="Tumor")
# 
# #join the two tables based on their gene names
# x<-inner_join(normal,tumor,by="gene")%>%
#         filter(samp.x==samp.y)
# write.csv(x,"7WEX.output_reshuffled.csv",append = F,sep = ",",row.names = F,col.names = T)
# 
# str(x)
# x<-select(x,sample=samp.x,gene,variant_class=variant_class.x)

