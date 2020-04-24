library(shiny)
library(plotly)
library(stringdist)
library(ggplot2)
library(DT)
library(htmlwidgets)

final.db <- readRDS("final.db.rds")
final.db <- final.db %>% dplyr::mutate(TS = paste0(final.db$TISSUE, "_", final.db$STRAIN))
final.db$TS <- gsub("[PMLN]{2,3}_","", final.db$TS)
final.db$TS <- gsub("naive_.*","naive", final.db$TS)

#function to fix xpoints
xdist <- function(x){
  s <- seq(0,20)  
  idx <- s %in% x[["lv"]]
  f <- data.frame(lv = s[!idx], ratio = 0, counts = 0)
  x <- rbind(x,f)
  return(x[order(x[["lv"]]),])
}

js <- c(
  "function(el, x){",
  "  var d3 = Plotly.d3;",
  "  el.on('plotly_restyle', function(evtData) {",
  "    var out = {};",
  "    d3.select('g.legend').selectAll('.traces').each(function(){",
  "      var trace = d3.select(this)[0][0].__data__[0].trace;",
  "      out[trace.name] = trace.visible;",
  "    });",
  "    Shiny.setInputValue('traces', out);",
  "  });",
  "}")



# Define server logic required to draw a histogram ----
server <- function(input, output) {

 
  db.sub  <- reactive({
  #subset options 
  if(input$vgene != "NULL"){
    db.sub <- final.db %>% filter (VGENE == input$vgene)
  }
  else{
    db.sub <- final.db
  }

  #other filters
  
  cw <- db.sub %>% dplyr::count(CLONE_ID) %>% filter(n >= input$cloneWells)
  db.sub <- db.sub[db.sub$CLONE_ID %in% cw[[1]],]
  
  db.sub <- db.sub[db.sub$READ.COUNTS.SUM >= input$readFilter,]

  if(input$norm == TRUE){
  db.sub <- as.data.frame(db.sub %>% dplyr::group_by(CLONE_ID) %>% 
                            summarise(READ.COUNT.SUM = sum(READ.COUNTS.SUM), 
                                      JUNCTION_AA =  JUNCTION_AA[1], 
                                      VGENE = VGENE[1], 
                                      STRAIN = unique(STRAIN),
                                      TS = unique(TS)))
  }
  else{return(db.sub)}

  
})


  ##Output
hist <- reactive({   
  #Calculate dist
  dist <- stringdist(c(db.sub()[["JUNCTION_AA"]]), input$pattern, method = "lv")
  prop.dist <- data.frame(lv = dist)
  hist <- prop.dist %>% dplyr::group_by(lv) %>% summarise(dplyr::n())
  colnames(hist) <- c("lv", "counts")
  hist$ratio <- hist$counts/sum(hist$counts)
  hist <- xdist(hist)
 })  

hists.fn <- reactive({
  xyz <- function(x){
    dist <- stringdist(c(x[["JUNCTION_AA"]]), input$pattern, method = "lv")
    prop.dist <- data.frame(lv = dist)
    hist <- prop.dist %>% dplyr::group_by(lv) %>% summarise(dplyr::n())
    colnames(hist) <- c("lv", "counts")
    hist$ratio <- hist$counts/sum(hist$counts)
    hist <- xdist(hist)
  }
  
  df <- split(db.sub(), db.sub()[["STRAIN"]])
  hists.fn <- lapply(df, FUN = xyz)
  
  return(dplyr::bind_rows(hists.fn, .id = "strain"))
  
})

hists.naive <- reactive({
  xyz <- function(x){
    dist <- stringdist(c(x[["JUNCTION_AA"]]), input$pattern, method = "lv")
    prop.dist <- data.frame(lv = dist)
    hist <- prop.dist %>% dplyr::group_by(lv) %>% summarise(dplyr::n())
    colnames(hist) <- c("lv", "counts")
    hist$ratio <- hist$counts/sum(hist$counts)
    hist <- xdist(hist)
  }
  
  df <- split(db.sub(), db.sub()[["TS"]])
  hists.naive <- lapply(df, FUN = xyz)
  
  return(dplyr::bind_rows(hists.naive, .id = "ts"))
  
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
      t <- "Number"
    }
    else{t <- "Ratio"}
    
    ylab <- list(
      title = t,
      titlefont = f
    )
    
if(input$radio.plotly == "all"){
    distPlot1 <- plot_ly(
      x = hist()[["lv"]],
      y = hist()[[input$radio.fm]],
      name = "all",
      type = 'scatter',
      mode = 'lines+markers') %>%
      layout(xaxis = xlab, yaxis = ylab, title = input$pattern, showlegend  = TRUE)
    distPlot1 %>% onRender(js)
    distPlot1}
    
else if(input$naive1 == FALSE){
  
  distPlot1 <- plot_ly(
  type = 'scatter',
  mode = 'lines+markers') %>%
  layout(xaxis = xlab, yaxis = ylab, title = input$pattern, showlegend  = TRUE)

for(name in c("GF", "MM12", "SPF")){
  df <- dplyr::filter(hists.fn(), strain == name)
  distPlot1 = add_trace(distPlot1,
                        type = "scatter",
                        mode = 'lines+markers',
                        x = as.numeric(df[["lv"]]),
                        y = df[[input$radio.fm]],
                        name = name)
}
distPlot1 %>% onRender(js)
distPlot1}
else if (input$naive1 == TRUE){
  distPlot1 <- plot_ly(
    type = 'scatter',
    mode = 'lines+markers') %>%
    
    layout(xaxis = xlab, yaxis = ylab, title = input$pattern, showlegend  = TRUE)
  for(name in c("GF", "MM12", "SPF", "naive")){
    df <- dplyr::filter(hists.naive(), ts == name)
    distPlot1 = add_trace(distPlot1,
                          type = "scatter",
                          mode = 'lines+markers',
                          x = as.numeric(df[["lv"]]),
                          y = df[[input$radio.fm]],
                          name = name, inherit = T)
  }
  distPlot1 %>% onRender(js)
  
}
    

    
    
    
    
    
    
    })

  

  observeEvent(input$radio.plotly, ignoreNULL = FALSE,{
    if(input$radio.plotly == "dynamic"){
     output$dtable = DT::renderDataTable(hists.fn(), server = FALSE, rownames= FALSE)
    }
    
    if(input$radio.plotly == "all"){
      output$dtable = DT::renderDataTable(hist(), server = FALSE, rownames= FALSE)
    }
    
  })
  
  observeEvent(input$naive1,{
    if(input$naive1 == TRUE){
      output$dtable = DT::renderDataTable(hists.naive(), server = FALSE, rownames= FALSE)  
    }
    else{
      if(input$radio.plotly == "dynamic"){
      output$dtable = DT::renderDataTable(hists.fn(), server = FALSE, rownames= FALSE)}
      if(input$radio.plotly == "all"){
        output$dtable = DT::renderDataTable(hist(), server = FALSE, rownames= FALSE)
      }
      
    }
  })   
  
  
  
  output$text1 <- renderText({ 
    sprintf("You have selected %s clone.wells", nrow(db.sub()))
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