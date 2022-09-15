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

# Load dataset
mash_cds <- read.table(file = "total_med_mash_dist_cds.csv", 
                     sep=",",
                     header = F)                                 #load cds data
mash_cds[2] <- NULL   #remove second col.


colnames(mash_cds) <- c("cluster.name","mash_cds")                     #create headers


mash_pr <- read.table(file = "total_med_mash_dist_pr.csv",
                    sep=",",
                    header = F)                                 #load cds data
mash_pr[2] <- NULL   #remove second col.


colnames(mash_pr) <- c("cluster.name","mash_pr")                     #create headers



mash_all <- merge(mash_cds, mash_pr, by="cluster.name", all=TRUE) #merga two tables



#cat pangenome_matrix_t0__shell_list.txt pangenome_matrix_t0__softcore_list.txt > pangenome.csv

pangenome <- read.table(file = "pangenome.csv",
                        sep=",",
                        header = F) #load pangenome data


#define cols names
colnames(pangenome)<- c("cluster.name","pan")#create headers




#remove all spaces from columns 
mash_all <- apply(mash_all, 2, str_remove_all, " ") 
pangenome <- apply(pangenome, 2, str_remove_all, " ") 

#convert data to dataframe
mash_all <-  as.data.frame(mash_all)
pangenome <- as.data.frame(pangenome)

# Merge data according to row names
mash_all <- merge(mash_all,pangenome, by ="cluster.name",all=TRUE)
mash_all <- as.data.frame(mash_all)

#remove all data contain NA
mash_all <- mash_all[!is.na(mash_all$mash_cds), ]
mash_all <- mash_all[!is.na(mash_all$mash_pr), ]


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
final_mash <- merge(mash_all,occupancy_number, by ="cluster.name",all=TRUE)


final_mash$mash_cds <- as.numeric(as.character(final_mash$mash_cds))
final_mash$mash_pr <- as.numeric(as.character(final_mash$mash_pr))



#melt cds and promoters in one cone column 
mash_data <- melt(final_mash, measure.vars = c("mash_cds","mash_pr"))
colnames(mash_data)<- c("cluster.name","pangenome","occupancy_number","cds_promoter","mash_distance")
mash_data <- mash_data[!is.na(mash_data$mash_distance), ]

# creat levels for x axis
data_levels <- levels(as.factor(as.numeric(as.character(mash_data$occupancy_number))))

#make order to x axis depend on levels orders
mash_data$occupancy_number <- factor(mash_data$occupancy_number, levels=(data_levels)) 

# customize the arrangement of pangenome order
mash_data$pangenome <- factor(mash_data$pangenome, levels=c('shell','core'))


#plot
#ggsave(file="ahh.png", width=18, height=7, dpi=300)

ggplot(mash_data)+(aes(x=occupancy_number, y=mash_distance,color=cds_promoter,fill=cds_promoter)) + 
  geom_boxplot() + facet_grid(.~pangenome, scale="free", space="free")+geom_boxplot() +
  theme(strip.text.x = element_text(size =15))+
  scale_fill_manual(name="",values = c("#fc1717", "#3422f5"),labels = c("CDS", "Promoters"))+
  scale_colour_manual(name="",values = c("#ff5e5e","#8d85ed"),labels = c("CDS", "Promoters"))+
  xlab("Occupancy")+ ylab("mash distance ")+
  theme(axis.title=element_text(size=15))+
  theme(legend.text = element_text(size=15))+
  stat_compare_means(method="t.test",label ="p.signif")
#dev.off()
