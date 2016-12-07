library(shiny)
require(ggplot2)
require(dplyr)
#setwd("/Users/chet/uky/SNPdensityShinyDraft/")
source("/Users/chet/uky/SNP_density_windows/R/HelperScripts.R")
#Load in dataset.
#in the future, rebuild this to allow for switching between multiple datasets
source("helperFunctions.R")
winSize = 1000
Br7raw <- read.table("inputData//Br7WindowsAvgMin.txt", header=TRUE)
Br7df <- Br7raw %>% group_by(Scaffold, start) %>%
  transform(adjustedSNP= ((winSize-minSNPs)/winSize)*100 )
Br7Lengths<- read.table("inputData/Br7LengthFile", header=TRUE)

chosenDF <- Br7df
scaffoldListDF<-  chosenDF %>%select(Scaffold, start)  %>%distinct(Scaffold, start)%>% group_by(Scaffold) %>% summarise(windowCount= n())


###
#UI
###
library('shiny')
ui <- fluidPage(
  titlePanel("Locus inheritance tracker, v.01"),
  sidebarLayout(
    sidebarPanel(
      uiOutput("scaffoldListOutput"),
    sliderInput(inputId = "minWindows", label= "Minimum number of 1kb windows/scaffold:", min=0, max=350, value=10),
    sliderInput(inputId = "binsToPlot", label= "Number of bins to plot.  (note: more bins = more time)", min=0, max=300, value=30),
        plotOutput("histOfAll")
    ),
    mainPanel(
      verbatimTextOutput("instructionPanel"),
      plotOutput("plotOfScaffold"),
      tableOutput("selectedTable")
    )
  )
)


###
#Server
###


server <- function(input, output) {
  
  ###hist of original
  output$histOfAll<-renderPlot({
    hist(scaffoldListDF$windowCount, main = "number of 1kb windows/scaffold", xlab = "# windows")
  })
  
  #filter list based on slider
  filtered <- reactive({  ###create filtered dataset
  #set the dataset to null if it hasnt loaded yet
    scaffoldListDF %>% filter (windowCount  >input$minWindows) 
      })
  
  ##Create list of scaffolds
  output$scaffoldListOutput <- renderUI({ 
    selectInput("scaffoldInput", "Scaffold to plot", as.character(filtered()$Scaffold), selected = NULL)
  })
  
  #Plot out main window plot
  
  output$plotOfScaffold <- renderPlot({
#    if (is.null(filtered())) {     return()  } #if our dataset is still null, dont try and draw it
   datasetToPlot <-  chosenDF %>% filter (Scaffold ==  input$scaffoldInput ) #filter dataset to selected
   
  toPlot <- outCladeMaxOnly(inputDF = convertToPercent(datasetToPlot, windowSize = 1000), toFilter= "wheat")
  plotMaxAncenstorWindow(toPlot, title = paste(input$scaffoldInput, "tile plot of most similar out taxa clade"),
                         numberOfBins=input$binsToPlot)
   
  })
  ##
  ##end main plot
  
  #make table of selected scaffold
  
  output$selectedTable <- renderTable({
    if (is.null(filtered())) { #if our dataset is still null, dont try and draw it
      return()
    }
    datasetToPlot <-  chosenDF %>% filter (Scaffold ==  input$scaffoldInput ) #filter dataset to selected
    
     toTable<-  outCladeMaxOnly(inputDF = convertToPercent(datasetToPlot, windowSize = 1000), toFilter= "wheat")
    
    toTable %>%  group_by(Clade) %>% summarise( NumberOf1kbWindowsAssigned= n()) #filter dataset to selected
     
  })
  ##
  ##end outputtable
  
 output$instructionPanel<- renderText({paste("Select a scaffold on the dropdown menu to plot it.  
You can filter scaffolds by the number of 1kb windows considered using the slider.
   \n\n Currently plotting:", input$scaffoldInput )})
} 

shinyApp(ui = ui, server = server)
