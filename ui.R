library(shiny)
library(plotly)

final.db <- readRDS("final.db.rds")

fluidPage(
  # App title ----
  titlePanel("Sulk-Seq - Levenshtein distance calculation"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      textInput("pattern", 
                label = "Choose a pattern",
                value = "CAREGFAYW",
                width = NULL,
                placeholder = NULL),
      
      
      textInput("vgene", 
               label = "Choose a V-Gene",
               value = "IGHV1-47",
               width = NULL,
               placeholder = NULL),
      
      
      radioButtons("radio", "Groups",
                   choices = c("All","GF","SPF","MM12"),
                   selected = "All")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotlyOutput(outputId = "distPlot")
      
    )
  )
)