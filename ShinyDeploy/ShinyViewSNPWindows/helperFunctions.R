

##
#Converts the raw SNP window input to percent based on a window size
##

convertToPercent <- function(rawSNPinput, windowSize = 1000 ) {
require(dplyr)
  return(rawSNPinput %>% group_by(Scaffold, start) %>%
    transform(adjustedSNP= ((windowSize-minSNPs)/windowSize)*100 ))
}

##
#Simple "not in" function
##
'%!in%' <- function(x,y)!('%in%'(x,y))

###
#Remove outclades from SNP input.
#Then, return the max (IE the most similar) hit for each window
###

outCladeMaxOnly <- function(inputDF, toFilter = c("wheat")){

  #always filter NA's as well
toFilter <- c(toFilter, NA)

percentSNPsfiltered <- inputDF %>% group_by(Scaffold, start, Clade) %>%
  filter(Clade %!in% toFilter)

windowsMaxOnly<-  percentSNPsfiltered %>%group_by(Scaffold, start)%>%
  filter(adjustedSNP == max(adjustedSNP))

return(windowsMaxOnly)
}


#
#Plot
#

plotMaxAncenstorWindow<- function(scaffoldToPlot, title = "add a title!", numberOfBins = 30 ) {
ggplot(scaffoldToPlot, aes(x=start, color = Clade, fill = Clade)) + geom_histogram(bins= numberOfBins)+ 
  ggtitle(title)+
  xlab("Scaffold position")+ylab("number of windows matching clade")
  
}