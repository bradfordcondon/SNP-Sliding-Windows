####
###Given a SNP file, variable location, and the reference name, trim the BLAST report
####
# trimvariableToTrimsFromSNPTables<- function(SNPFile, variableToTrim= "fileName" , refName) {
#   SNPFile[variabletoTrim]<- gsub(paste(refName, "_merged_masked_v_", sep = ""), "", SNPFile[variableToTrim])
#   SNPFile[variabletoTrim]<- gsub("_masked_out", "", SNPFile[variableToTrim])
#   SNPFile[variabletoTrim]<- gsub("_merged", "", SNPFile[variableToTrim])
#   SNPFile[variabletoTrim]<- gsub("contigs", "", SNPFile[variableToTrim])
#   return(SNPFile)
# }



#for some reason this function kills my r enviornment.  for now will do indiviudally
trimFileNamesFromSNPTables<- function(SNPFile, refName) {
  SNPFile$fileName<- gsub(paste(refName, "_merged_masked_v_", sep = ""), "", SNPFile$fileName)
  SNPFile$fileName<- gsub("_masked_out", "", SNPFile$fileName)
  SNPFile$fileName<- gsub("_merged", "", SNPFile$fileName)
  SNPFile$fileName<- gsub("contigs", "", SNPFile$fileName)
  return(SNPFile)
} 

###
#given a path, read in all files in that path, and create a dataframe assuming all those files are SNP reports
###
readSNPReports <- function(path){
  outDF <- data.frame()
  FileList <- list.files(path)
  for (thisFile in FileList){ 
    thisTable <- read.table(paste(path, thisFile, sep=""), fill = TRUE, col.names=c("scaffoldRef", "ScaffoldQ", "locRef", "locQ", "baseR", "baseQ", "NA", "NA" ))  
    thisTable<- thisTable[,1:6]
    thisTable<- thisTable[complete.cases(thisTable),]
    thisTable$fileName <- thisFile   #add a sixth column with the name of the file
    outDF<- rbind(outDF, thisTable)
  }
  return(outDF)
}


##
#Given: SNP dataframe, (optional: window size and step size), return windows 
#organized by scaffold with the sum of all SNPs for each strain in that window.
##
slidingWindowForSNP <- function(SNPdf, windowSize = 1000, stepSize = 1000){  
  scaffoldLoop<- c(1:length(unique(SNPdf$scaffoldRef)))#Loop for each scaffold
  allScaffoldsTracker <- list()
  for (i in scaffoldLoop) {
    thisScaffoldTracker <- list()
    thisScaffold <- as.character(unique(SNPdf$scaffoldRef)[i])#identify scaffold for this loop
    subsetDF <-  filter(SNPdf, scaffoldRef == thisScaffold)    #subset data to this scaffold
    total <- as.numeric(select(filter(br80lengths, scaffold == thisScaffold), length)) #determine length of this scaffold from scaffold list
    if (total < windowSize){  spots <- seq(from=1, to = 1)  } else {   spots <- seq(from=1, to=(total-windowSize), by=stepSize)} # create sequence of intervals }
    for(j in 1:length(spots)){    #for each interval...
      #filter dataset to this interval
      thisRange <- spots[j]:(spots[j]+windowSize)
      DFforRange<-  filter(subsetDF, locRef %in% thisRange)
      #summarize stats for this interval
      groupedDFrange<- group_by(DFforRange, fileName, Clade)
      thisScaffoldTracker[[j]] <- summarise(groupedDFrange, count = n())
    }
    names(thisScaffoldTracker)<- as.character(spots)  #add names to the sublist
    
    allScaffoldsTracker[[i]]  <-  thisScaffoldTracker  #add this list of summarized windows to the master list
    #names(allScaffoldsTracker[[i]])<- thisScaffold
  }
  names(allScaffoldsTracker)<-  as.character(unique(SNPdf$scaffoldRef))
  return(allScaffoldsTracker)
}

##
#Same as above, except output a dataframe
##

slidingWindowForSNPoutputDF <- function(SNPdf, windowSize = 1000, stepSize = 1000){  
  scaffoldLoop<- c(1:length(unique(SNPdf$scaffoldRef)))#Loop for each scaffold
  allScaffoldsTracker <- data.frame(Scaffold=NA, window=NA, match=NA, clade=NA, value = NA, stringsAsFactors = FALSE)[numeric(0), ]
  for (i in scaffoldLoop) {
    thisScaffold <- as.character(unique(SNPdf$scaffoldRef)[i])#identify scaffold for this loop
    subsetDF <-  filter(SNPdf, scaffoldRef == thisScaffold)    #subset data to this scaffold
    total <- as.numeric(select(filter(br80lengths, scaffold == thisScaffold), length)) #determine length of this scaffold from scaffold list
    if (total < windowSize){  spots <- seq(from=1, to = 1)  } else {   spots <- seq(from=1, to=(total-windowSize), by=stepSize)} # create sequence of intervals }
    for(j in 1:length(spots)){    #for each interval...
      #filter dataset to this interval
      thisRange <- spots[j]:(spots[j]+windowSize)
      DFforRange<-  filter(subsetDF, locRef %in% thisRange)
      #summarize stats for this interval
      groupedDFrange<- summarise(group_by(DFforRange, fileName, Clade), count=n())
      #add back columns for the scaffold and window
      thisWindowDF<- mutate(groupedDFrange, Scaffold = thisScaffold, Window = spots[j])
      #Rbind this to DF
    allScaffoldsTracker<- rbind(allScaffoldsTracker, as.data.frame(thisWindowDF))    
          }
  }
  return(allScaffoldsTracker)
}

###
#Not in
#The inverse in %in%
###
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


###
#Input: Sliding window dataframe
#return: dataframe summary
###

analyzeSlidingWindows<- function(allScaffoldsTracker){
dfTracker <- data.frame(Scaffold=NA, window=NA, match=NA, value = NA, stringsAsFactors = FALSE)[numeric(0), ]
ScaffoldList <- as.character((unique(allScaffoldsTracker$Scaffold)))
loop = c(1:length(ScaffoldList)) #make sure this s valid!
for (i in loop){
  thisScaffold <- ScaffoldList[i]
  thisScaffoldDF <- filter(allScaffoldsTracker, Scaffold == thisScaffold) #filter DF to this scaffold
  windowList <- unique(thisScaffoldDF$Window) #determine how many windows are in this DF
  loopTwo <- c(1:length(windowList))
  for (j in loopTwo){
    thisWindow <-  windowList[j] 
    thisWindowDF<- filter(thisScaffoldDF, Window==thisWindow)
    windowSum<- summarise(group_by(thisWindowDF, Clade), avg = mean(count), min = min(count), max = max(count))
    FLTB <- filter(windowSum, Clade == 'FLTB'  ) #just FLTB
    notFLTB <- arrange(filter(windowSum, Clade != 'FLTB'  ), min)[1,]#lowest minimum, not FLTB
    windowRange <- paste( thisWindow, "-", as.numeric(thisWindow) + windowSize, sep="")
    if(nrow(notFLTB)  ==0) {   #If no FLTB snps, assume FLTB is best match
      dfTracker <- rbind(dfTracker, data.frame(Scaffold=thisScaffoldName, window=windowRange, match="FLTB", value=0 ))
      next } 
    if (is.na(notFLTB[1,1]) | is.na(FLTB[1,1])){next }
    else if (notFLTB[3] < FLTB[2] ){ #if the lowest number of SNPs  not in FLTB is less than the average FLTB snps...
      dfTracker <- rbind(dfTracker, data.frame(Scaffold=thisScaffoldName, window=windowRange, match=notFLTB$Clade[1], value=notFLTB$min[1] ))
    }
    else if (notFLTB[3] > FLTB[2] ){ #
      dfTracker <- rbind(dfTracker, data.frame(Scaffold=thisScaffoldName, window=windowRange, match=FLTB$Clade[1], value=FLTB$avg[1] ))
    }
  }
}
return(dfTracker)
}
##
#Plot SNPs to multiPDF
#Requires: dataframe of SNPs, a name for the PDF.
#Default arguments are for the binsize and the y limit
##

PlotSNPstoMultiPDF <- function( DFtoPrint , pdfName, binSize = 2000, ylim = 50 ){ 
  fullLim = c(1, ylim)
scaffoldList <- row.names(sort(table(DFtoPrint$scaffoldRef), decreasing = TRUE))
pdf(pdfName, 8.5, 11)
for (i in scaffoldList){
  snpsToPlot<- DFtoPrint[DFtoPrint$scaffoldRef == i,]
  title <- i
  
  if (nrow(snpsToPlot) > 0 ) { #dont try and print null
    p <- ggplot(snpsToPlot, aes(x=locRef, fill= Clade)) + geom_histogram(binwidth = binSize)+ ggtitle(title)+xlab("Scaffold position")+
        facet_wrap(~fileName, drop = TRUE, ncol = 2)+coord_cartesian(ylim = fullLim)
    plot(p)
  }
}
dev.off()

}


