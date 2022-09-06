#Code to install packages if necessary, and read them with library function
required_packages <- c("car", "tidyverse", "ggplot2", "ggrepel", "readr","here")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {

    library(package, character.only = TRUE)
  } else {

    install.packages(package)

    library(package, character.only = TRUE)
  }
}

#set your in current directory 
here::set_here()
getwd()

#setwd("~/ammar/brachy_snakemake/tables")

# Load dataset os quality table
quality_table <- read.table(file = "tables/quality_table.csv",
                          sep=",", #separator is ,
                          header = T, # first row is a header
                          row.names = 1) # first col. is row names

#ramove number of genes from talbe
df <- subset(quality_table,select = -c(N_genes))

#boxplot for N_count, Gaps, Busco score
tiff("figures/quality_box.eps", units="in", width=5, height=7, res=300) #.eps for illustrator

df_plot <- Boxplot(scale(df[,2:4])) #Boxplot from car package give name outlier

dev.off() #save a plot photo


# define the ecoytypes outlier name without dupilcation 
outlier_ecotypes <- unique(df_plot)

#names of ecotypes to exclude 
outlier_ecotypes


