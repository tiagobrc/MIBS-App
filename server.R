library(shiny)
library(plotly)
library(stringdist)
library(ggplot2)
library(DT)
library(sulk)
library(dplyr)
library(htmlwidgets)

final.db <- readRDS("final.db.rds")

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {

  db.Filt  <- reactive({
 #Filtering reads
  db.Filt <- final.db[final.db[["READ.COUNTS.SUM"]] >= input$readFilter,]
 #Filtering CloneWells
  cw <- db.Filt %>% count(CLONE_ID) %>% filter(n >= input$cloneWells)
  db.Filt <- db.Filt[db.Filt$CLONE_ID %in% cw[[1]],]
  #Filtering V genes
  #subset options 
  if(input$vgene != "NULL"){
    db.Filt <- db.Filt %>% filter (VGENE == input$vgene)
    }
    
  if(input$norm == TRUE){
dups <- as.data.frame(db.Filt %>% filter(!duplicated(CLONE_ID)))
sums <- db.Filt %>% group_by(CLONE_ID) %>% summarise(sum(READ.COUNTS.SUM))
db.Filt <- full_join(x = sums, y = dups, by = "CLONE_ID")}
  return(db.Filt)
})

  ##Output
hist <- reactive({
hist <- group.dist(x = db.Filt(), summarise.by = "count", split.by = db.Filt()[[input$selec.plotly]], pattern = input$pattern)
 return(hist)
})

  output$distPlot <- renderPlotly({
    #Plotly config
    
    f <- list(
      family = "Arial",
      size = 18,
      color = "#7f7f7f"
    )
    
    xlab <- list(
      title = sprintf("CDR3 AA Levenshtein distance From %s", input$pattern),
      titlefont = f
    )
    
    if(input$radio.fm == "counts"){
      t <- "Number"}
    
    else{t <- "Ratio"}
    
    ylab <- list(
      title = t,
      titlefont = f
    )
    
    distPlot1 <- plot_ly(
      type = 'scatter',
      mode = 'lines+markers')
      
    for(trace in names(hist())){
      distPlot1 <- distPlot1 %>% 
        plotly::add_trace(x = hist()[[trace]][["lv"]], # converted Date as date
                          y = hist()[[trace]][[input$radio.fm]], name = trace)}
    
    distPlot1 <- distPlot1 %>%
    layout(xaxis = xlab, yaxis = ylab, title = input$pattern, showlegend  = TRUE) %>% 
    onRender(js)

})
  
  

  output$text1 <- renderText({ 
    sprintf("You have selected %s clone.wells", nrow(db.Filt()))
  })
  
  output$text2 <- renderText({ 
    if(input$norm == TRUE){
      print("In-mouse clone.wells count as 1")
    }else{
      print("In-mouse clone.wells count as it is")
    }
  })
  
  # download the filtered data
  output$dhist = downloadHandler('histcsv.xlsx', content = function(file) {
      
    openxlsx::write.xlsx(hist(), file, row.names = FALSE)
  })
  
  
  
}