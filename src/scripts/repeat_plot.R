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
#set your in current directory 
here::set_here()
getwd()

# Load dataset
data <- read.table(file = "tables/repeats_analysis.csv",
                   sep=",",
                   header = T)

#lenght of non repeated Seq.
data$Non.repeated <- (data$Genome.size - data$Repeated )

#
data <- subset(data,select = -c(Genome.size))


#data <- melt(data)
data <- melt(data, value.name = "N_bp", variable.name = "Condition")


# Stacked + percent

tiff("figures/repeats_plot.eps", units="in", width=7, height=5, res=300) #.eps for illustrator

ggplot(data, aes(fill=Condition,y=N_bp, x=Ecotypes)) + 
  geom_bar(position="stack", stat="identity", width=0.5) +
  theme(axis.text.x=element_text(angle=90,size = rel(0.8), margin = margin(1, unit = "cm"),vjust =1))+ #change size and angle of x axis
  xlab("Ecotypes")+            # for the x axis label
  ylab("length (bp)") + theme(legend.key.height= unit(0.3, 'cm'), #change size of legend.key
        legend.key.width= unit(0.3, 'cm'))+
  theme(legend.title = element_text(size=5))+theme(legend.text = element_text(size=5))+ # change size of tilte and text for legend squre 
theme(axis.title=element_text(size=6,face="bold"))

dev.off() #save a plot photo


