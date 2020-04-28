library(shiny)
library(plotly)
library(stringdist)
library(ggplot2)
library(DT)
library(FSA)
library(sulk)
library(tidyr)
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
    hist <- group.dist(x = db.Filt(), summarise.by = "count", split.by = db.Filt()[[input$selec.plotly]], 
                       pattern = input$pattern)
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
    
    t <- input$radio.fm
    
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
  
 #Table for stats 
  mat <- reactive({
    if(input$table.fmt == "Plain"){
      
      dist <- stringdist(db.Filt()[["JUNCTION_AA"]], input$pattern, method = "lv")
      l <- split(dist, db.Filt()[[input$selec.plotly]])
      l <- lapply(l, sort)
      n.obs <- sapply(l, length)
      seq.max <- seq_len(max(n.obs))
      mat <- sapply(l, "[", i = seq.max)
      
      return(mat)
    }
  })
  
  
  # if(input$table.fmt == "Histogram"){
  #   output$dtable <- renderDataTable({bind_rows(hist(), .id = "Groups")})  
  # }
  
  tab <- reactive({
    if(input$table.fmt == "Plain"){
      return(mat())
    }
    
    if(input$table.fmt == "Histogram"){
      return(
        bind_rows(hist(), .id = "Groups"))
    }
  })
  
  output$dtable <- renderDataTable({
    tab()
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
    
    openxlsx::write.xlsx(tab(), file, row.names = FALSE)
  })
  
  output$boxp <- renderPlot({
    
    if(input$table.fmt == "Plain"){
      dt <- as.data.frame(mat())
      my_data <- gather(dt, "Groups")
      #   
      p1 <- ggpubr::ggboxplot(data = my_data, x = "Groups", y = "value",
                        color = "Groups", palette = scales::hue_pal()(length(unique(my_data$Groups))),
                        order = unique(my_data$Groups),
                        ylab = "Lv Distance", xlab = "Groups")
      return(p1)
    }
    
    # if(input$table.fmt == "Histogram"){
    #   
    #   p2 <- ggpubr::ggboxplot(tab(), x = "Groups", y = "counts",
    #                     color = "Groups", palette = scales::hue_pal()(length(unique(tab()$Groups))),
    #                     order = unique(tab()$Groups),
    #                     ylab = "Lv Distance", xlab = "Groups")
    #   return(p2)
   # }
    
  }) 
  
  # output$stats1 <- renderPrint({
  #   
  #   if(input$table.fmt == "Plain"){
  #     dt <- as.data.frame(mat())
  #     my_data <- gather(dt, "Groups")
  #     
  #     if(unique(my_data$Groups) > 1){
  #       str1 <- kruskal.test(value ~ Groups, data = my_data)
  #       return(str1)
  #       }
      
#       else{str1 <- paste("all observations are in the same group")
#       return(str1)
#       }
#       
#       
#       }
# })
  #Wilcox
  # output$stats2 <- renderDT(rownames = T,{
  #   if(input$table.fmt == "Plain"){
  #   dt <- as.data.frame(mat())
  #   my_data <- gather(dt, "Groups")
  #   
  #   str2 <- pairwise.wilcox.test(my_data$value, my_data$Groups,
  #                                p.adjust.method = "BH")
  #   
  #   return(str2$p.value)
  #   
  #   }
  # })
  
  #Dunn's test
  observeEvent(eventExpr = input$stats,{
    output$stats3 <- renderDT({
      
      if(input$table.fmt == "Plain"){
        dt <- as.data.frame(mat())
        my_data <- gather(dt, "Groups")
        
        str3 <- dunnTest(value ~ Groups, data=my_data, method="bh")
        
        return(as.data.frame(str3$res))
      }
      })
    })
  
output$warn <- renderText({ 
  if(input$table.fmt == "Histogram"){
  return(print("Choose the plain table for testing and make sure to melt data using Plot Group selection"))
}})
  
  # observe({
  #   x <- input$vgene
  # 
  #   # Can use character(0) to remove all choices
  #   if (is.null(x))
  #     x <- character(0)
  # 
  #   # Can also set the label and select items
  #   updateSelectInput(session, "inSelect",
  #                     label = paste("Plot the percentage of: ", "All or ", x),
  #                     choices = c("All", paste(x, "vgene clones")),
  #                     selected = "All"
  #   )
  # })
  # 

}

