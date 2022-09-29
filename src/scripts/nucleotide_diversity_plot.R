#!/usr/bin/env Rscript

#Code to install packages if necessary, and read them with library function
required_packages <- c("stringr","ggplot2", "readr","reshape2","ggpubr","dplyr","tidyverse")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

#set your data directory
here::set_here()
getwd()

# Load dataset
nd_cds <- read.table(file = "nucleotide_diversity_cds_new.csv",
                     sep=",",
                     header = F)                                 #load cds data

colnames(nd_cds) <- c("cluster.name","piCDS")                    #create headers


nd_pr <- read.table(file = "nucleotide_diversity_pr_new.csv",
                     sep=",",
                     header = F)                                 #load cds data

colnames(nd_pr) <- c("cluster.name","piPromoter")                #create headers



nd_all <- merge(nd_cds, nd_pr, by="cluster.name", all=TRUE) #merga two tables



#cat pangenome_matrix_t0__shell_list.txt pangenome_matrix_t0__softcore_list.txt > pangenome.csv

pangenome <- read.table(file = "pangenome.csv",
                       sep=",",
                       header = F) #load pangenome data


#define cols names
colnames(pangenome)<- c("cluster.name","pan")#create headers

#convert data to dataframe
nd_all <-  as.data.frame(nd_all)
pangenome <- as.data.frame(pangenome)


#remove all spaces from columns
nd_all <- apply(nd_all, 2, str_remove_all, " ")
pangenome <- apply(pangenome, 2, str_remove_all, " ")


# Merge data according to row names
all <- merge(nd_all,pangenome, by ="cluster.name",all=TRUE)
all <- as.data.frame(all)

#remove all data contain NA
all <- all[!is.na(all$piCDS), ]

#remove data

#load occupancy table
# cat BdistachyonABR2337v1_4taxa_algOMCL_e1_.cluster_list |
# grep "taxa"|cut -d " " -f 4,6|perl -p -e 's/taxa=/,/'|
# perl -lane 'print $F[1],$F[0]'> ../../../tables/occupancy_number.csv

occupancy_number <- read.table(file = "occupancy_number.csv",
                               sep=",",
                               header = F)
#define cols names

colnames(occupancy_number)<- c("cluster.name","occupancy_number")

#convert data to dataframe
occupancy_number_new <-  as.data.frame(occupancy_number)


#remove all spaces from columns

occupancy_number_new <- apply(occupancy_number_new, 2, str_remove_all, " ")

# Merge data according to row names
final <- merge(all,occupancy_number_new, by ="cluster.name",all=TRUE)
final_table <- final[!is.na(final$piCDS), ]


#remove data
remove(final,occupancy_number_new,occupancy_number,all,pangenome,nd_all,nd_pr,nd_cds)

final_table$piCDS <- as.numeric(as.character(final_table$piCDS))
final_table$piPromoter <- as.numeric(as.character(final_table$piPromoter))



#melt cds and promoters in one cone column
data <- melt(final_table, measure.vars = c("piCDS","piPromoter"))
colnames(data)<- c("cluster.name","pangenome","Occupancy_number","cds_promoter","Nucleotide_diversity")

# creat levels for x axis
data_levels <- levels(as.factor(as.numeric(as.character(data$Occupancy_number))))

#make order to x axis depend on levels orders
data$Occupancy_number <- factor(data$Occupancy_number, levels=(data_levels))

# customize the arrangement of pangenome order
data$pangenome <- factor(data$pangenome, levels=c('shell','core'))


#plot
#ggsave(file="ahh.png", width=18, height=7, dpi=300)

ggplot(data)+(aes(x=Occupancy_number, y=Nucleotide_diversity,color=cds_promoter,fill=cds_promoter)) +
  geom_boxplot() + facet_grid(.~pangenome, scale="free", space="free")+geom_boxplot() +
  theme(strip.text.x = element_text(size =15))+
  scale_fill_manual(name="",values = c("#fc1717", "#3422f5"),labels = c("CDS", "Promoters"))+
  scale_colour_manual(name="",values = c("#ff5e5e","#8d85ed"),labels = c("CDS", "Promoters"))+
  xlab("Occupancy number")+ ylab("Nucleotide diversity")+
  theme(axis.title=element_text(size=15))+
  theme(legend.text = element_text(size=15))+
  stat_compare_means(method="wilcox.test",label ="p.signif")

#dev.off()

#creat summary table
data_summary <- final_table %>%  #for this table
  group_by(occupancy_number, pan) %>%  #make a group by this variable (tidyverse package)
  summarise(Number = n(), MeanCDS = mean(piCDS), MeanPromoters = mean(piPromoter),medianCDS= median(piCDS),medianPromoters= median(piPromoter)) # create with summaries col you want

#create table per occupancy number
table_10 <- data[data$Occupancy_number ==10,]
table_25 <- data[data$Occupancy_number == 25,]


#histogram to see distribution of this tables
table_10 %>% ggplot( aes(x=Nucleotide_diversity,fill= cds_promoter)) +
      geom_histogram(bins =30, alpha=0.6)


table_25 %>% ggplot( aes(x=Nucleotide_diversity,fill= cds_promoter)) +
  geom_histogram(bins =30, alpha=0.6)
