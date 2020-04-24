library(shiny)
library(plotly)
library(DT)

final.db <- readRDS("final.db.rds")

fluidPage(
  # App title ----
  titlePanel("Sulk-Seq analysis"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      sliderInput("readFilter", "Minimum read per sequence:",
                  min = 0, max = 1000, value = 5, step = 5
      ),
      
      sliderInput("cloneWells", "Minimum number of clone.wells",
                  min = 0, max = 32, value = 1
      ),
      
      checkboxInput("norm", "Collapse in-mouse clones?", value = TRUE, width = NULL),
      
      
      h4("Levenshtein distance controls:"),
      textInput("pattern", 
                label = "Choose a pattern",
                value = "CAREGFAYW",
                width = NULL,
                placeholder = NULL),
      
      h6("Choose the V-Gene for subsetting."),
      selectizeInput("vgene", 
                     label = "Choose a V-Gene",
                     choices = c('All V-Genes' = 'NULL', unique(final.db$VGENE)),
                     selected = NULL,
                     width = NULL),
      
      radioButtons("radio.fm", "x-axis format",
                   choices = c("counts","ratio"),
                   selected = "ratio"),
      
      radioButtons("radio.plotly", "Plot type",
                   choices = c("all","dynamic"),
                   selected = "all"),
      
      checkboxInput("naive1", "Detach naive?", value = FALSE, width = NULL),
      
      #Download Plotly graphic
      
      helpText(),
      selectInput("selectplot", "Download plot", choices=list("Plot1.svg")),
      actionButton('download', "Download")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      #n of clone.wells
      textOutput("text1"),
      textOutput("text2"),  
      
      # Output: Histogram ----
      plotlyOutput(outputId = "distPlot"),
      tags$script('
                  document.getElementById("download").onclick = function() {
                  var plot = $("#selectplot").val();
                  var gd = document.getElementById("distPlot");
                  
                  Plotly.Snapshot.toImage(gd, {format: "svg"}).once("success", function(url) {
                    var a = window.document.createElement("a");
                    a.href = url; 
                    a.type = "image/svg";
                    a.download = "plot.svg";
                    document.body.appendChild(a);
                    a.click();
                    document.body.removeChild(a);                      
                  });
                  }
                  '),
      
      DT::dataTableOutput('dtable'),
      
      
      downloadButton('dhist', 'Download hist')
      
      
        

      
      
      
    )
  )
)