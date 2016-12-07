library(shiny)
require(ggplot2)
require(dplyr)
source("/Users/chet/uky/SNP_density_windows/R/HelperScripts.R")
#Load in dataset.
#in the future, rebuild this to allow for switching between multiple datasets

winSize = 1000
Br7raw <- read.table("/Users/chet/uky/SNP_density_windows/inputData/Br7WindowsAvgMin.txt", header=TRUE)
Br7df <- Br7raw %>% group_by(Scaffold, start) %>%
  transform(adjustedSNP= ((winSize-minSNPs)/winSize)*100 )
Br7Lengths<- read.table("/Users/chet/uky/SNP_density_windows/inputData/Br7LengthFile", header=TRUE)

chosenDF <- Br7df
scaffoldList<-as.character(unique(chosenDF$Scaffold))


source("/Users/chet/uky/SNP_density_windows/ShinyViewSNPWindows/Server.R")
source("/Users/chet/uky/SNP_density_windows/ShinyViewSNPWindows/UI.R")

shinyApp(ui = ui, server = server)
