library(shiny)
library(plotly)
library(DT)

final.db <- readRDS("final.db.rds")

lv.up <- function(x){length(unique(x)) > 2 && length(unique(x)) < 88}

groupby <- colnames(final.db)[unlist(lapply(final.db, lv.up))]

fluidPage(
  # App title ----
  titlePanel("Sulk-Seq analysis"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
     h5("Filter controls"),
      sliderInput("readFilter", "Minimum read per sequence:",
                  min = 0, max = 1000, value = 5, step = 5
      ),
      
      sliderInput("cloneWells", "Minimum number of clone.wells",
                  min = 0, max = 32, value = 1
      ),
      
      checkboxInput("norm", "Collapse in-mouse clones?", value = TRUE, width = NULL),
      
      h6("V-Gene filter"),
      selectizeInput("vgene", 
                     label = "Choose a V-Gene",
                     choices = c('All V-Genes' = 'NULL', unique(final.db$VGENE)),
                     selected = NULL,
                     width = NULL),
      
      
      h5("Levenshtein distance controls:"),
      textInput("pattern", 
                label = "Choose a pattern",
                value = "CAREGFAYW",
                width = NULL,
                placeholder = NULL),
      
      
      
      radioButtons("radio.fm", "y-axis format",
                   choices = c("counts","ratio","% of all clones in condition"),
                   selected = "ratio"),
      
      selectizeInput("selec.plotly", 
                     label = "Plot group:",
                     choices = c("All" = "All", groupby),
                     selected = "All",
                     width = NULL),
      
    

      
      #Download Plotly graphic
      
      helpText(),
      selectInput("selectplot", "Download plot", choices=list("Plot1.svg")),
      actionButton('download', "Download"),
      
    
      radioButtons("stats", "Statistical tests",
                   choices = c("Kruskal-Wallis"),
                   selected = "Kruskal-Wallis"),
      actionButton('gime.stats', "Give me Stats!"),
      
      
      
      
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