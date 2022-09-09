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




################################################################################################
#local analysis/promoter

#number of total sequance per local promoter cluster
total_seq_pr_local <- read.table(file = "total_seq_local_pr.csv",
                                 sep=",",
                                 header = F)

# define coloums headers
colnames(total_seq_pr_local)<- c("cluster.name","total_number_local_pr")


#remove any space
total_seq_pr_local <- apply(total_seq_pr_local, 2, str_remove_all, " ") 
total_seq_pr_local <-  as.data.frame(total_seq_pr_local)

#merge side by side by cluster.name
promoter_local_data <- merge(total_seq_pr,total_seq_pr_local, by ="cluster.name",all=TRUE)

#remove any space
promoter_local_data <- apply(promoter_local_data, 2, str_remove_all, " ") 
promoter_local_data <-  as.data.frame(promoter_local_data)



#convert number data value as numeric for calculation
promoter_local_data$total_number_pr <- as.numeric(as.character(promoter_local_data$total_number_pr))
promoter_local_data$total_number_local_pr <- as.numeric(as.character(promoter_local_data$total_number_local_pr))


#calculate the poor local alighment number
promoter_local_data$poor_number_local_pr <- (promoter_local_data$total_number_pr - promoter_local_data$total_number_local_pr )

promoter_local_data[3] <- NULL   





#local analysis/cds

#number of total sequance per local cds cluster
total_seq_cds_local <- read.table(file = "total_seq_local_cds.csv",
                                  sep=",",
                                  header = F)


# define coloums headers
colnames(total_seq_cds_local)<- c("cluster.name","total_number_local_cds")


#remove any space
total_seq_cds_local <- apply(total_seq_cds_local, 2, str_remove_all, " ") 
total_seq_cds_local <-  as.data.frame(total_seq_cds_local)

#merge side by side by cluster.name
cds_local_data <- merge(total_seq_cds,total_seq_cds_local, by ="cluster.name",all=TRUE)

#remove any space
cds_local_data <- apply(cds_local_data, 2, str_remove_all, " ") 
cds_local_data <-  as.data.frame(cds_local_data)



#convert number data value as numeric for calculation
cds_local_data$total_number_cds <- as.numeric(as.character(cds_local_data$total_number_cds))
cds_local_data$total_number_local_cds <- as.numeric(as.character(cds_local_data$total_number_local_cds))
#calculate the poor local alighment number
cds_local_data$poor_number_local_cds <- (cds_local_data$total_number_cds - cds_local_data$total_number_local_cds )

cds_local_data[3] <- NULL  

#summary table
#promoter
#global
promoter_data <- merge(promoter_data,occupancy_number, by ="cluster.name",all=TRUE)
promoter_data$occupancy <- as.numeric(as.character(promoter_data$occupancy ))



table_core_pr <- promoter_data[promoter_data$occupancy >= 51, ]
sum_pr_core <- sum(table_core_pr$poor_number_pr)



table_shell <- promoter_data[promoter_data$occupancy < 51, ]

sum_pr_shell<- sum(table_shell$poor_number_pr)

#local
promoter_local_data <- merge(promoter_local_data,occupancy_number, by ="cluster.name",all=TRUE)
promoter_local_data$occupancy <- as.numeric(as.character(promoter_local_data$occupancy ))



table_core_pr_local <- promoter_local_data[promoter_local_data$occupancy >= 51, ]
sum_pr_core_local <- sum(table_core_pr_local$poor_number_local_pr)



table_pr_shell_local <- promoter_local_data[promoter_local_data$occupancy < 51, ]

sum_pr_shell_local<- sum(table_pr_shell_local$poor_number_local_pr)



#cds
#global
cds_data <- merge(cds_data,occupancy_number, by ="cluster.name",all=TRUE)
cds_data$occupancy <- as.numeric(as.character(cds_data$occupancy ))


table_core_cds <- cds_data[cds_data$occupancy >= 51, ]
sum_cds_core <- sum(table_core_cds$poor_number_cds)



table_shell_cds <- cds_data[cds_data$occupancy < 51, ]

sum_cds_shell<- sum(table_shell_cds$poor_number_cds)


#local
cds_local_data <- merge(cds_local_data,occupancy_number, by ="cluster.name",all=TRUE)
cds_local_data$occupancy <- as.numeric(as.character(cds_local_data$occupancy ))



table_core_cds_local <- cds_local_data[cds_local_data$occupancy >= 51, ]
sum_cds_core_local <- sum(table_core_cds_local$poor_number_local_cds)



table_cds_shell_local <- cds_local_data[cds_local_data$occupancy < 51, ]

sum_pr_shell_local<- sum(table_cds_shell_local$poor_number_local_cds)





sum_pr <- sum(promoter_data$total_number_pr)
sum_pr_poor_global <- sum(promoter_data$poor_number_pr)
sum_pr_poor_local <- sum(promoter_local_data$poor_number_local_pr)


sum_cds <- sum(cds_data$total_number_cds)
sum_cds_poor_global <- sum(cds_data$poor_number_cds)
sum_cds_poor_local <- sum(cds_local_data$poor_number_local_cds)



