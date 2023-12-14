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

# Get the directory path of the currently active script
current_path <- normalizePath(commandArgs(trailingOnly = TRUE)[1])

# Print or use the relative path as needed
print(current_path)

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


te_precentage_all <- merge(count_te_cds, count_te_pr, by="cluster.name", all=TRUE) #merga two tables
#cat pangenome_matrix_t0__shell_list.txt pangenome_matrix_t0__softcore_list.txt > pangenome.csv

pangenome <- read.table(file = "pangenome.csv",
                        sep=",",
                        header = F) #load pangenome data


#define cols names
colnames(pangenome)<- c("cluster.name","pan")#create headers




#remove all spaces from columns 
te_precentage_all <- apply(te_precentage_all, 2, str_remove_all, " ") 
pangenome <- apply(pangenome, 2, str_remove_all, " ") 

#convert data to dataframe
te_precentage_all <-  as.data.frame(te_precentage_all)
pangenome <- as.data.frame(pangenome)

# Merge data according to row names
te_precentage_all <- merge(te_precentage_all,pangenome, by ="cluster.name",all=TRUE)
te_precentage_all <- as.data.frame(te_precentage_all)

#remove all data contain NA
te_precentage_all <- te_precentage_all[!is.na(te_precentage_all$`TE_cds%`), ]
te_precentage_all <- te_precentage_all[!is.na(te_precentage_all$`TE_pr%`), ]


#load occupancy table 
# cat BdistachyonABR2337v1_4taxa_algOMCL_e1_.cluster_list |
# grep "taxa"|cut -d " " -f 4,6|perl -p -e 's/taxa=/,/'|
# perl -lane 'print $F[1],$F[0]'> ../../../tables/occupancy_number.csv

occupancy_number <- read.table(file = "occupancy_number.csv",
                               sep=",",
                               header = F)
#define cols names

colnames(occupancy_number)<- c("cluster.name","occupancy_number")


#remove all spaces from columns 

occupancy_number <- apply(occupancy_number, 2, str_remove_all, " ") 


#convert data to dataframe
occupancy_number <-  as.data.frame(occupancy_number)






# Merge data according to row names
final_te <- merge(te_precentage_all,occupancy_number, by ="cluster.name",all=TRUE)
final_te$`TE_cds%` <- as.numeric(as.character(final_te$`TE_cds%`))
final_te$`TE_pr%` <- as.numeric(as.character(final_te$`TE_pr%`))



#melt cds and promoters in one cone column 
te_data <- melt(final_te, measure.vars = c("TE_cds%","TE_pr%"))
colnames(te_data)<- c("cluster.name","pangenome","occupancy_number","cds_promoter","TEs_precentage")
te_data <- te_data[!is.na(te_data$TEs_precentage), ]

#choose only promoter 
#te_data <- te_data %>% filter(cds_promoter == "TE_pr%")

# creat levels for x axis
data_levels <- levels(as.factor(as.numeric(as.character(te_data$occupancy_number))))

#make order to x axis depend on levels orders
te_data$occupancy_number <- factor(te_data$occupancy_number, levels=(data_levels)) 

# customize the arrangement of pangenome order
te_data$pangenome <- factor(te_data$pangenome, levels=c('shell','core'))

#plot

p <- ggplot(te_data)+(aes(x=occupancy_number, y=TEs_precentage,color=cds_promoter,fill=cds_promoter)) + 
  geom_boxplot() + facet_grid(.~pangenome, scale="free", space="free")+
  theme(strip.text.x = element_text(size =15))+
  scale_fill_manual(name="",values =c("#fc1717", "#3422f5"),labels =  c("CDS", "Promoters"))+
  scale_colour_manual(name="",values = c("#ff5e5e","#8d85ed"),labels = c("CDS", "Promoters"))+
  xlab("Occupancy")+ ylab("TEs_precentage")+
  theme(axis.title=element_text(size=15))+
  theme(legend.text = element_text(size=15)) +
  stat_compare_means(method="t.test",label ="p.signif")

# Save the plot
ggsave(file="figures/ TEs_content.png", plot=p, width=18, height=7, dpi=300)

# Close the graphics device
dev.off()


#Rscript {input.script_mash_plot} /home/ammar/ammar/snakemake_improve
