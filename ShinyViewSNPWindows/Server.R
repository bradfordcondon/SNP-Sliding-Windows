library('shiny')


server <- function(input, output) {
  filtered <- reactive({  ###create filtered dataset
    if (is.null(input$scaffoldInput)) {
      return(NULL)
    }  #set the dataset to null if it hasnt loaded yet
    chosenDF %>%
      filter(
             Scaffold == input$scaffoldInput
      )#
  })
 
  ##
  #Render the UI
  ##
  
  output$countryOutput <- renderUI({  #create country list from selection
    selectInput("scaffoldInput", "Scaffold",
                sort(unique(chosenDF$scaffoldInput)),
                selected = NULL)
  })
  
  ###
  #Render the plot
  ##
  
  output$coolplot <- renderPlot({
    if (is.null(filtered())) { #if our dataset is still null, dont try and draw it
      return()
    }
    plotSelectedScaffold( )# inputDF, scaffoldLengths, windowSize
  })
  output$results <- renderTable({
    filtered()
  })
  
} 



