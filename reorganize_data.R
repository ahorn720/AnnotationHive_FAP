
library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2)
x<-read.csv(file = "7WEX.output.csv", header = T)
x$sample<-as.character(x$sample)
x$Stage <- ifelse(grepl("_T_", x$sample), "Tumor", "Normal")
arrange(x,Stage)%>%
        separate(sample,into="samp",sep = "_")

#create two new tables; a normal and a tumor table
normal<-filter(x,Stage=="Normal")
tumor<-filter(x,Stage=="Tumor")

#join the two tables based on their gene names
x<-inner_join(normal,tumor,by="gene")%>%
        filter(samp.x==samp.y)
write.csv(x,"7WEX.output_reshuffled.csv",append = F,sep = ",",row.names = F,col.names = T)
