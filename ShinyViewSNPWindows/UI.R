library('shiny')
ui <- fluidPage(
  titlePanel("Locus inheritance tracker, v.01"),
  sidebarLayout(
    sidebarPanel(
      selectInput("scaffoldInput", "Scaffold to plot", scaffoldList, selected = NULL),
      uiOutput("selectedScaffoldOutput")), #make country a reactive
    
    mainPanel(plotOutput("plotOfScaffold"), tableOutput("results")
              
    )
  )
)




