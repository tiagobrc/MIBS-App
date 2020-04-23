library(shiny)
library(plotly)
library(stringdist)
library(ggplot2)

final.db <- readRDS("final.db.rds")
# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  output$distPlot <- renderPlotly({
    
    if(input$radio != "All"){
    db.sub <- final.db %>% filter (STRAIN == input$radio)
    db.sub <- db.sub %>% filter (VGENE == input$vgene)
    }
    else{
      db.sub <- final.db %>% filter (VGENE == input$vgene)
    }
    
    
    
    
    
    dist <- stringdist(c(db.sub$JUNCTION_AA), input$pattern, method = "lv")
    
    
    
    
    
    prop.dist <- data.frame(lv = dist)
    hist <- prop.dist %>% dplyr::group_by(lv) %>% summarise(dplyr::n())
    colnames(hist) <- c("lv", "Counts")
    
    hist$Counts <- hist$Counts/sum(hist$Counts)
    
    x <- reactive(hist$lv)
    y <- reactive(hist$Counts)
    
    
    f <- list(
      family = "Courier New, monospace",
      size = 18,
      color = "#7f7f7f"
    )
    
    
    xlab <- list(
      title = sprintf("CDR3 AA Levenshtein distance From %s", input$pattern),
      titlefont = f
    )
    ylab <- list(
      title = "Ratio (Sum = 1)",
      titlefont = f
    )
    
    
    
    distPlot <- plot_ly(
      x = x(),
      y = y(),
      type = 'scatter',
      mode = 'lines+markers') %>%
      
      
    layout(xaxis = xlab, yaxis = ylab, title = input$pattern)
    
    
    
  })
  
}