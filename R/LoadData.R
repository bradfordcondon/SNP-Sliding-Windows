##
#Script for loading in and trimming all SNP count data
#Created October 26, 2016
##
#load dependencies etc
library(plyr)
library(ggplot2)
library(knitr)
library(dplyr)
#read in helper functions
source("R/helperFunctions.R")

#Read in SNP reports and scaffold lengths
BR80DF  <- readSNPReports("inputData/br80/")
IA1SDF  <- readSNPReports("inputData/IA1/")
WBKYDF <- readSNPReports("/Users/chet/uky/WBKY_SNP_reports/")


br80AllSnpDF<-  trimFileNamesFromSNPTables(SNPFile = BR80DF, refName = "Br80")
IA1AllSnpDF<-  trimFileNamesFromSNPTables(SNPFile = IA1SDF, refName = "IA1")
WBKYSnpDF<-  trimFileNamesFromSNPTables(SNPFile = WBKYDF, refName = "WBKY")


br80lengths <-read.table("supporting/Br80_length_report.txt")
IA1lengths <- read.table("V1/supporting/IA1_length_report.txt")
colnames(br80lengths)<- c("scaffold", "length")
colnames(IA1lengths)<- c("scaffold", "length")



#read in clade list
clade_list<- read.table("supporting/clade_list_two_class.txt", header = TRUE)
kable(clade_list)
#rename problematic clades in the clade_list
clade_list$ID<- as.character(clade_list$ID)
clade_list[clade_list$ID %in% "G17",2]<- "G17_2M"
clade_list[clade_list$ID %in% "GY11",2]<- "Guy11"
clade_list[clade_list$ID %in% "PY86-1",2]<- "PY86"
clade_list[clade_list$ID %in% "SSFL02-1",2]<- "SSFL02"
clade_list[clade_list$ID %in% "WBKY11",2]<- "WBKY"
#replace lower case oryza with upper case oryza
clade_list[clade_list$WG_clade %in% "oryza",7] <- "Oryza"
#finally, let's merge FLTB with TB
clade_list[clade_list$WG_clade %in% "TB",7] <- "FLTB"

#remove clade list items not found in thedataset
clade_list_br80<- clade_list[clade_list$ID %in% br80AllSnpDF$fileName,]
clade_list_IA1<- clade_list[clade_list$ID %in% IA1AllSnpDF$fileName,]

#Add clade to each SNP 
loop <- c(1:length(clade_list_br80$ID))
br80AllSnpDF$Clade <-""#introduce blank variable

for (i in loop){
  thisTaxon <- as.character(clade_list_br80$ID[i])
  thisClade <- as.character(clade_list_br80[which(clade_list_br80$ID %in% thisTaxon),]$WG_clade)
  br80AllSnpDF[which(br80AllSnpDF$fileName  %in% thisTaxon),]$Clade <-thisClade
}

#drop those missing clades
br80AllSnpDF<- (br80AllSnpDF[which(br80AllSnpDF$Clade != ""),])

#Add clade to each SNP 
loop <- c(1:length(clade_list_IA1$ID))
IA1AllSnpDF$Clade <-""#introduce blank variable

for (i in loop){
  thisTaxon <- as.character(clade_list_IA1$ID[i])
  thisClade <- as.character(clade_list_IA1[which(clade_list_IA1$ID %in% thisTaxon),]$WG_clade)
  IA1AllSnpDF[which(IA1AllSnpDF$fileName  %in% thisTaxon),]$Clade <-thisClade
}

IA1AllSnpDF<- (IA1AllSnpDF[which(IA1AllSnpDF$Clade != ""),])

datasetList<- unique(WBKYAllSnpDF$fileName)
toRemove<- datasetList[!datasetList %in% clade_list$ID]
#remove these from dataset
WBKYAllSnpDF <-filter(WBKYAllSnpDF, fileName!= toRemove)

#remove clade list items not found in thedataset
clade_list_WBKY<- clade_list[clade_list$ID %in% WBKYAllSnpDF$fileName,]


#Add clade to each SNP 
loop <- c(1:length(clade_list_WBKY$ID))
WBKYAllSnpDF$Clade <-""#introduce blank variable

for (i in loop){
  thisTaxon <- as.character(clade_list_WBKY$ID[i])
  thisClade <- as.character(clade_list_WBKY[which(clade_list_WBKY$ID %in% thisTaxon),]$WG_clade)
  WBKYAllSnpDF[which(WBKYAllSnpDF$fileName  %in% thisTaxon),]$Clade <-thisClade
}

#BR7
Br7DF <- readSNPReports("inputData/Br7/")
br7AllSnpDF<-  trimFileNamesFromSNPTables(SNPFile = Br7DF, refName = "Br7")
datasetList<- unique(br7AllSnpDF$fileName)
toRemove<- datasetList[!datasetList %in% clade_list$ID]
#remove these from dataset
br7AllSnpDF <-filter(br7AllSnpDF, fileName!= toRemove)
clade_list_BR7<- clade_list[clade_list$ID %in% br7AllSnpDF$fileName,]
loop <- c(1:length(clade_list_BR7$ID))
br7AllSnpDF$Clade <-""#introduce blank variable
for (i in loop){
  thisTaxon <- as.character(clade_list_BR7$ID[i])
  thisClade <- as.character(clade_list_BR7[which(clade_list_BR7$ID %in% thisTaxon),]$WG_clade)
  br7AllSnpDF[which(br7AllSnpDF$fileName  %in% thisTaxon),]$Clade <-thisClade
}

save.image("LoadDataOut.RDATA")
