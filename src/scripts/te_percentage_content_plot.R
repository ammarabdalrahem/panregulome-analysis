#Code to install packages if necessary, and read them with library function
required_packages <- c("stringr","ggplot2", "readr","reshape2","ggpubr",
                       "dplyr","tidyverse","dgof","hrbrthemes","ggpol","magrittr")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

#set your data directory
setwd("C:/Users/ammar/Desktop/workR")

# Load dataset
count_te_cds <- read.table(file = "count_te_sm_cds.csv", 
                       sep=",",
                       header = F)      #load cds data

colnames(count_te_cds) <- c("cluster.name","total_length","te_length","TE_cds%")                     #create headers
count_te_cds[, -c(1,4)] <- NULL   #remove  col.

count_te_pr <- read.table(file = "count_te_sm_pr.csv",
                      sep=",",
                      header = F)                                 #load cds data


colnames(count_te_pr) <- c("cluster.name","total_length","te_length","TE_pr%")   
count_te_pr[, -c(1,4)] <- NULL   #remove  col.

count_te_pr$cluster.name<- gsub('BdistachyonBd21v2_283', 'BdistachyonBd21v2_1', count_te_pr$cluster.name) 

te_precentage_all <- merge(count_te_cds, count_te_pr, by="cluster.name", all=TRUE) #merga two tables

#melt cds and promoters in one cone column 
te_precentage_all <- melt(te_precentage_all, measure.vars = c("TE_pr%","TE_cds%"))
colnames(te_precentage_all)<- c("cluster.name","TE_cds_promoter","TEs_precentage_content")
te_precentage_all$TEs_precentage_content <- as.numeric(as.character(te_precentage_all$TEs_precentage_content))
te_precentage_all <- te_precentage_all[!is.na(te_precentage_all$TEs_precentage_content), ]
te_precentage_all <- apply(te_precentage_all, 2, str_remove_all, " ") 
te_precentage_all <-  as.data.frame(te_precentage_all)

#plot
#ggsave(file="ahh.png", width=18, height=7, dpi=300)
ggplot(te_precentage_all, aes(x = cluster.name,  y = TEs_precentage_content,fill= TE_cds_promoter)) +
  geom_histogram(position = "identity", alpha = 0.4)
#dev.off()



ggplot(te_precentage_all, aes(x=cluster.name, y=TEs_precentage_content, fill=TE_cds_promoter)) + 
  geom_bar(stat="identity", position="identity", alpha = 0.4)+
  scale_fill_manual(name="",values = c("#fc1717", "#3422f5"),labels = c("CDS", "Promoters"))+
  scale_colour_manual(name="",values = c("#ff5e5e","#8d85ed"),labels = c("CDS", "Promoters"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  


