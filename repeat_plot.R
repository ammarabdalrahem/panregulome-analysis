#!/usr/bin/env Rscript

#Code to install packages if necessary, and read them with library function
required_packages <- c("ggplot2","readr","reshape2","rstudioapi")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}
#set your data directory
#anyway to make the directory what i use automatic?
setwd("~/ammar/work")

# Load dataset
data <- read.table(file = "repeats_analysis.csv",
                   sep=",",
                   header = T)

#lenght of non repeated Seq.
data$Non.repeated <- (data$Genome.size - data$Repeated )

#
data <- subset(data,select = -c(Genome.size))


#data <- melt(data)
data <- melt(data, value.name = "N_bp", variable.name = "Condition")


# Stacked + percent
png("repeats_plot.png")

ggplot(data, aes(fill=Condition,y=N_bp, x=Ecotypes)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  xlab("Ecotypes")+            # for the x axis label
  ylab("length (bp)") 

dev.off()


