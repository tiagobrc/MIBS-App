library(shiny)
library(plotly)
library(stringdist)
library(ggplot2)
library(DT)
library(FSA)
library(MIBS)
library(tidyr)
library(Biostrings)
library(bios2mds)
library(dplyr)
library(htmlwidgets)
library(msaR)
library(msa)


# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  ##########  Reacts  ##########
 
  #Inputs
  selectedData <- reactiveValues()
  
  #Filter data
  db.Filt  <- reactive({
    #Filtering reads
    db.Filt <- selectedData$db %>% filter(READ.COUNTS.SUM >= input$readFilter)
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
  
  db.seqs <- reactive({
    #Filtering reads
    db.Filt <- selectedData$db %>% filter(READ.COUNTS.SUM >= input$readFilter2)
    #Filtering CloneWells
    cw <- db.Filt %>% count(CLONE_ID) %>% filter(n >= input$cloneWells2)
    db.Filt <- db.Filt[db.Filt$CLONE_ID %in% cw[[1]],]
    
    #Filtering V genes
    #subset options 
    if(input$vgene2 != "NULL"){
      db.Filt <- db.Filt %>% filter(VGENE == input$vgene2)
    }else{
      db.Filt <- db.Filt 
    }
    
    if(input$jgene2 != "NULL"){
      db.Filt <- db.Filt %>% filter(JGENE == input$jgene2)
    }else{
      db.Filt <- db.Filt 
    }
    
    if(input$random != 0){
      db.Filt <- db.Filt %>% sample_n(input$random)
    }else{
      db.Filt <- db.Filt 
    }
    
    selectedData$dbseqs.df <- db.Filt
    seqs <- db.Filt[["SEQUENCE_VDJ"]]
    names(seqs) <- db.Filt[["SEQUENCE_ID"]]
    return(seqs)
    
  })
  
  #Options for groups
  opt1 <- reactive({
    colnames(selectedData$db)[unlist(lapply(selectedData$db, lv.up, min = 2, max = 88))]
  })
  
  #Options for v genes
  opt2 <- reactive({
    unique(selectedData$db[["VGENE"]]) 
  })
  
  opt3 <- reactive({
    unique(selectedData$db[["JGENE"]]) 
  })
 
  
  
  
  hits <- reactive({
    mindex <- vmatchPattern(pattern = input$pattern2, subject = selectedData$translated, min.mismatch = 0)
    hits <- as(mindex, "CompressedIRangesList")
    hits  
  })
  
  
  
  
df.bar  <- eventReactive(input$predict, {
  
  l <- input$lookahead
  
  flanks <- flank(hits(), width = l)  
  
  df <- as.data.frame(extractAt(selectedData$translated, flanks))
  full <- full_join(selectedData$dbseqs.df, df, by = c("SEQUENCE_ID" = "group_name"))
  
  selectedData$missed <- table(is.na(full[["value"]]))[["TRUE"]]
  selectedData$missed.df <- full %>% filter(is.na(full[["value"]])) %>% group_by(CLONE_ID) %>% count() %>% arrange(desc(n))
  
  
  df.bar <- data.frame(Size = c(length(table(full[["JUNCTION_AA"]])), length(table(full[["value"]]))),
                       Groups = c("Truth", "Sequenced"))
  
 df.bar
  
  
 })
  
  
  
  #Histogram
  hist <- reactive({
    hist <- group.dist(x = db.Filt(), summarise.by = "count", split.by = db.Filt()[[input$selec.plotly]], 
                       pattern = input$pattern)
    return(hist)
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
  
  tab <- reactive({
    if(input$table.fmt == "Plain"){
      return(mat())
    }
    
    if(input$table.fmt == "Histogram"){
      return(
        bind_rows(hist(), .id = "Groups"))
    }
  })
  
  
  ##########  Observatory  ##########
  #Update input selections
  observe({
    if (input$radio.db == "sulk1") {
      selectedData$db <- MIBS::sulk1
    }
  })
  
  
  observe({
    
    updateSelectizeInput(session, "selec.plotly",
                         choices = c("All" = "All", opt1()),
                         selected = "All",
                         server = FALSE)
  })
  
  observe({
    
    updateSelectizeInput(session, "vgene",
                         choices = c("All" = 'NULL', opt2()),
                         selected = "NULL",
                         server = FALSE)
  })
  
  observe({
    
    updateSelectizeInput(session, "vgene2",
                         choices = c("All" = 'NULL', opt2()),
                         selected = "NULL",
                         server = FALSE)
  })
  
  observe({
    
    updateSelectizeInput(session, "jgene2",
                         choices = c("All" = 'NULL', opt3()),
                         selected = "NULL",
                         server = FALSE)
  })
  
  proteins <- eventReactive(input$align.seq, {
    set <- Biostrings::DNAStringSet(db.seqs())
    set_subseqs <- lapply(1:3, function(pos) {subseq(c(set), start=pos)})
    
    # for(i in 1:3){
    #   temp <- Biostrings::translate(set_subseqs[[i]], genetic.code = Biostrings::GENETIC_CODE)
    #   st <- grep("\\*", set.jj3)
    #   stops <- c(length(st), stops)
    # }
    suppressWarnings({ 
      translated <- Biostrings::translate(set_subseqs[[2]], genetic.code = Biostrings::GENETIC_CODE)
      })
    
    if(length(grep("\\*", translated)) != 0){
      msa.translated <- msa(inputSeqs = translated[-grep("\\*", translated)], method = "ClustalOmega")
      }
   else{
     msa.translated <- msa(inputSeqs = translated, method = "ClustalOmega")
   }
    
    file1 <- tempfile(fileext = ".fa")
    selectedData$translated <- translated
    export.fasta(msaConvert(msa.translated, "bios2mds::align"), file1)

    proteins <- ape::read.FASTA(file1, type="AA")
    unlink(file1)
    return(proteins)
  })
  
  
  
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
  
  

  ##########  Rendering  ##########
  
  
  output$df.barplot <- renderPlot({
    
    p <- ggplot(data=df.bar(), aes(x=Groups, y=Size, fill = Groups)) +
      geom_bar(stat="identity")+ theme_minimal()
    p
    
  })
  
  output$distPlot <- renderPlotly({
    #Plotly config
    validate(
      need(input$selec.plotly != '', 'Choose one group!'),
      need(input$vgene != '', 'Choose a V gene!')
    )
    f <- list(
      family = "Arial",
      size = 18,
      color = "#7f7f7f"
    )
    
    xlab <- list(
      title = sprintf("CDR3 AA Levenshtein distance from %s", input$pattern),
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
                          y = hist()[[trace]][[input$radio.fm]], name = trace)
    }
    
    distPlot1 <- distPlot1 %>%
      layout(xaxis = xlab, yaxis = ylab, title = input$pattern, showlegend  = TRUE) 
    
  })
  
  output$vcount <- renderTable({
    req(selectedData$translated)
    df <- as.data.frame(table(vcountPattern(pattern = input$pattern2, subject = selectedData$translated)))
    
    colnames(df) <- c("Hits", "Number")
    df
  })
  
  renderText({
    sprintf("The pattern is found n times in the aligned sequences:")
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
    
  })
  
  
  
  
  output$msa <- renderMsaR({
    msaR(msa = proteins(), menu=T, overviewbox = F, seqlogo =T, colorscheme = "clustal")
    })
  
  
  output$dtable <- renderDataTable({
    validate(
      need(input$selec.plotly != '', 'Choose one group!'),
      need(input$vgene != '', 'Choose a V gene!')
    )
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
  
  output$text3 <- renderText({ 
    sprintf("You have selected %s sequences", length(db.seqs()))
  })
  
  output$text4 <- renderText({ 
    x <- sprintf("%s sequences were missed by the pattern priming.", selectedData$missed)
    y <- sprintf("These %s sequences are grouped into %s clones", selectedData$missed, length(table(selectedData$missed.df[["CLONE_ID"]])))
    return(c(x,y))
  })
  
  output$missed.table <- renderDT({
    selectedData$missed.df
  })
  
  output$warn <- renderText({ 
    if(input$table.fmt == "Histogram"){
      return(print("Choose the plain table for testing and make sure to melt data using Plot Group selection"))
    }
    })
  
  
  ##########  Downloads  ##########
  
  output$dhist = downloadHandler('histcsv.xlsx', content = function(file) {
    
    openxlsx::write.xlsx(tab(), file, row.names = FALSE)
  })
  

  
  

  
 

}

