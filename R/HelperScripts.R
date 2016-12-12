
##
'%!in%' <- function(x,y)!('%in%'(x,y))


##
#Plot selected scaffold
##

plotSelectedScaffold<- function(inputDF, scaffoldLengths, windowSize){  
  #determine scaffold size
  thisScaffoldLength <- scaffoldLengths[which(scaffoldLengths$ID == inputDF$Scaffold[1]),2]
  
binNumber <- thisScaffoldLength/windowSize

  thisTitle <- paste(inputDF$Scaffold[1], "loci identity")
  
  p<- ggplot(inputDF, aes(x=start, color = Clade, fill = Clade)) + 
    geom_histogram(bins=binNumber)+ ggtitle(thisTitle)+xlab("Scaffold position")+
    ylab("number of  windows matching clade")
  return(p)
  }