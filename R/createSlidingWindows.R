#read in helper functions
source("R/helperFunctions.R")

##
#Perform sliding window analyses on genomes and save to an R object.
##

br80WindowsDF<- slidingWindowForSNPoutputDF(SNPdf = br80AllSnpDF, windowSize = 1000, stepSize = 1000)
IA1WindowsDF<- slidingWindowForSNPoutputDF(SNPdf = IA1AllSnpDF, windowSize = 1000, stepSize = 1000)


save(br80WindowsDF, file = "br80Windows_df_1kb.robj")
save(IA1WindowsDF, file = "IA1Windows_df_1kb.robj")
##
#Score sliding windows
##

br80WindowsSummary <- analyzeSlidingWindows(br80WindowsDF)
IA1WindowsSummary <- analyzeSlidingWindows(IA1WindowsDF)

save(br80WindowsSummary, file = "br80Windows_df_1kb_summary.robj")
save(IA1WindowsSummary, file = "IA10Windows_df_1kb_summary.robj")
