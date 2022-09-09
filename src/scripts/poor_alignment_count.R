#Code to install packages if necessary, and read them with library function
required_packages <- c("stringr","ggplot2", "readr","reshape2","ggpubr","dplyr","tidyverse","dgof")
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


occupancy_number <- read.table(file = "occupancy_number.csv",
                               sep=",",
                               header = F)
#define cols names

colnames(occupancy_number)<- c("cluster.name","occupancy")


#remove all spaces from columns 

occupancy_number <- apply(occupancy_number, 2, str_remove_all, " ") 


#convert data to dataframe
occupancy_number <-  as.data.frame(occupancy_number)


# creat levels for x axis
data_levels <- levels(as.factor(as.numeric(as.character(occupancy_number$occupancy))))

#make order to x axis depend on levels orders
occupancy_number$occupancy <- factor(occupancy_number$occupancy, levels=(data_levels)) 

#golbal analysis/promoter
#number of total sequance per cluster
total_seq_pr <- read.table(file = "total_seq_cluster_pr_.csv",
                           sep=",",
                           header = F)

#number of poor alignment for promoter per cluster
poor_number_pr <- read.table(file = "number_of_poor_alignment_pr.csv",
                             sep=",",
                             header = F)

# define coloums headers
colnames(total_seq_pr)<- c("cluster.name","total_number_pr")
colnames(poor_number_pr)<- c("cluster.name","poor_number_pr")


#remove any space
total_seq_pr <- apply(total_seq_pr, 2, str_remove_all, " ") 
total_seq_pr <-  as.data.frame(total_seq_pr)

poor_number_pr <- apply(poor_number_pr, 2, str_remove_all, " ") 
poor_number_pr <-  as.data.frame(poor_number_pr)

#merge side by side by cluster.name
promoter_data <- merge(total_seq_pr,poor_number_pr, by ="cluster.name",all=TRUE)

#remove any space
promoter_data <- apply(promoter_data, 2, str_remove_all, " ") 
promoter_data <-  as.data.frame(promoter_data)



#convert nuber data value as numeric for calculation
promoter_data$total_number_pr <- as.numeric(as.character(promoter_data$total_number_pr))
promoter_data$poor_number_pr <- as.numeric(as.character(promoter_data$poor_number_pr))

#calculate the ratio: (poor number /total nuber )*100
promoter_data$ratio_pr <- (promoter_data$poor_number_pr / promoter_data$total_number_pr)*100

#remove total and poor number keep only ratio col.
promoter_data[2:3] <- NULL   


#golbal analysis/CDS
#number of total sequance per cluster
total_seq_cds <- read.table(file = "total_seq_cluster_cds_.csv",
                            sep=",",
                            header = F)

#number of poor alignment for promoter per cluster
poor_number_cds <- read.table(file = "number_of_poor_alignment_cds.csv",
                              sep=",",
                              header = F)

# define coloums headers
colnames(total_seq_cds)<- c("cluster.name","total_number_cds")
colnames(poor_number_cds)<- c("cluster.name","poor_number_cds")


#remove any space
total_seq_cds <- apply(total_seq_cds, 2, str_remove_all, " ") 
total_seq_cds <-  as.data.frame(total_seq_cds)

poor_number_cds <- apply(poor_number_cds, 2, str_remove_all, " ") 
poor_number_cds <-  as.data.frame(poor_number_cds)

#merge side by side by cluster.name
cds_data <- merge(total_seq_cds,poor_number_cds, by ="cluster.name",all=TRUE)

#remove any space
cds_data <- apply(cds_data, 2, str_remove_all, " ") 
cds_data <-  as.data.frame(cds_data)



#convert nuber data value as numeric for calculation
cds_data$total_number_cds <- as.numeric(as.character(cds_data$total_number_cds))
cds_data$poor_number_cds <- as.numeric(as.character(cds_data$poor_number_cds))

#calculate the ratio: (poor number /total nuber )*100
cds_data$ratio_cds <- (cds_data$poor_number_cds / cds_data$total_number_cds)*100

#remove total and poor number keep only ratio col.
cds_data[2:3] <- NULL  



#merge all data for promoter and CDS
poor_seq_ratio_glob <- merge(promoter_data,cds_data, by ="cluster.name",all=TRUE)

#merge all data with occupancy number

poor_seq_ratio_glob <- merge(poor_seq_ratio_glob,occupancy_number, by ="cluster.name",all=TRUE)
#remove any space
poor_seq_ratio_glob <- apply(poor_seq_ratio_glob, 2, str_remove_all, " ") 
poor_seq_ratio_glob <-  as.data.frame(poor_seq_ratio_glob)

#melt all
poor_seq_ratio_glob <- melt(poor_seq_ratio_glob, measure.vars = c("ratio_pr","ratio_cds"))

# define coloums headers
colnames(poor_seq_ratio_glob)[4] <- "ratio"

#convert ratio value to numeric
poor_seq_ratio_glob$ratio <- as.numeric(as.character(poor_seq_ratio_glob$ratio))

#make order to x axis depend on levels orders
poor_seq_ratio_glob$occupancy<- factor(poor_seq_ratio_glob$occupancy, levels=(data_levels)) 


#plotting data
ggplot(poor_seq_ratio_glob)+(aes(x=occupancy, y=ratio,color=variable,fill=variable)) + 
  geom_boxplot() +
  theme(strip.text.x = element_text(size =15))+
  xlab("Occupancy")+ ylab("Percentage of poor sequance")+
  theme(strip.text.x = element_text(size =15))+
  scale_fill_manual(name="",values = c("#fc1717", "#3422f5"),labels = c("Promoters", "CDS"))+
  scale_colour_manual(name="",values = c("#ff5e5e","#8d85ed"),labels = c("Promoters", "CDS"))+
  stat_compare_means(method="t.test",label ="p.signif")


#clean data
remove(cds_data,promoter_data,poor_number_pr,poor_number_cds,occupancy_number,total_seq_pr,total_seq_cds)
