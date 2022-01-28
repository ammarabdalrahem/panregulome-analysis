#Code to install packages if necessary, and read them with library function
required_packages <- c("ggplot2","readr","reshape2","ggforce")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}
#set your data directory
setwd("~/ammar/trial")

# Load dataset
data <- read.table(file = "te_ann.csv",
                      sep=",",
                     header = T,
                     fill = TRUE)


# create TE orders
reorders = factor(data$Order, levels=c("LTR", "DNA", "TIR" ,"LARD", "Helitron" ,    
                                  "Satellite"   ,  "TRIM" , "MITE" ,        
                                  "Unclassified" , "rRNA", "Other"  ,"SINE"  , "DIRS",
                                  "LINE" , "MobileElement","Retroelement",  "LARD|TRIM" ,
                                  "non-LTR(SINE)" ,"SINE|LARD" , "DNAnona" ,    
                                  "nonLTR" , "RC" ,  "DNAauto" ,  "Evirus"))
#creat tool non in
`%nin%` = Negate(`%in%`)

#create data for LTR.DNA.TIR TE order only 
data_LTR.DNA.TIR <- (data[data$Order %in% c("LTR", "DNA","TIR" ),])

#create data for all TE order except (LTR.DNA.TIR) TE order 
data_other <- (data[data$Order %nin% c("LTR", "DNA","TIR" ),])


# plotting specific orders
#ggplot(data , aes(fill=Order,y=Length, x=Ecotype)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  xlab("Ecotypes")+            # for the x axis label
  ylab("length(bp)") +
  facet_zoom(x = Order  %nin% c("LTR", "DNA","TIR" ))

png("TE_orders_plot.png")
#plotting data with zoom the last 0.25% of data 
ggplot(data, aes(fill=reorders,y=Length, x=Ecotype)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  xlab("Ecotypes")+            # for the x axis label
  ylab("length(bp)")+
  facet_zoom(ylim = c(0.0, 0.025))

dev.off()
