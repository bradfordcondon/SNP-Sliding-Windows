##
#Perform sliding window analyses on genomes and save to an R object.
##

br80WindowsDF<- slidingWindowForSNPoutputDF(SNPdf = br80AllSnpDF, windowSize = 1000, stepSize = 1000)

IA1WindowsDF<- slidingWindowForSNPoutputDF(SNPdf = IA1AllSnpDF, windowSize = 1000, stepSize = 1000)


save(br80Windows, file = "br80Windows_df_1kb.robj")
save(IA1WindowsDF, file = "IA1Windows_df_1kb.robj")
