library(shiny)
library(shinythemes)
library(plotly)
library(DT)
library(bslib)
library(MIBS)



fluidPage(
  theme = bs_theme(version = 4, bootswatch = "minty"),
  
  sidebarLayout(
    sidebarPanel(
      style = "overflow-y:scroll; max-height: 800px; position:fixed;width:30%;",
      h2("Sulk App"),
      h5("Filter controls"),
      
      fileInput("dbInput", "Upload Database:", buttonLabel = "Browse...",  placeholder = "No file selected"),
      
      radioButtons("radio.db","Databases:",
                   choices = c("sulk1"),
                   selected = "sulk1"),
      
      sliderInput("readFilter", "Minimum read number per sequence:",
                  min = 0, max = 1000, value = 5, step = 5
      ),
      
      sliderInput("cloneWells", "Spanning n clone.wells:",
                  min = 0, max = 32, value = 1
      ),
      
      checkboxInput("norm", "Collapse in-mouse clones?", value = TRUE, width = NULL),
      
      h6("V-Gene filter"),
      selectizeInput("vgene", 
                     label = "Choose a V-Gene",
                     choices = c('All V-Genes' = 'NULL'),
                     selected = NULL,
                     width = NULL),
      
      
      h5("Levenshtein distance controls:"),
      textInput("pattern", 
                label = "Choose a pattern",
                value = "CAREGFAYW",
                width = NULL,
                placeholder = NULL),

      radioButtons("radio.fm", "y-axis format",
                   choices = c("counts","ratio"),
                   selected = "ratio"),
      
     # checkboxInput("total.perc", "Show percentage of total", value = FALSE),
    
      selectizeInput("selec.plotly", 
                     label = "Plot group:",
                     choices = c("All" = "All"),
                     selected = "All",
                     width = NULL),
      
      
      #Download Plotly graphic
      
      helpText(),
      selectInput("selectplot", "Download plot", choices=list("Plot1.svg")),
      actionButton('download', "Download"),
      
      
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(type = "tabs", id = "tabselected",
                  tabPanel("Levenshtein", value = 1,
                           
                           

                           DTOutput('tbl01'),
                           #n of clone.wells
                           textOutput("text1"),
                           textOutput("text2"),  
      
                           # Output: Histogram ----
                           plotlyOutput(outputId = "distPlot"),
                           tags$script('document.getElementById("download").onclick = function() {
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
                                       }'),
                           
                           DT::dataTableOutput('dtable'),
                           
                           radioButtons("table.fmt", 
                                        label = "Table format",
                                        choices = c("Plain", "Histogram"),
                                        selected = "Histogram",
                                        width = NULL),
                           
                           downloadButton('dhist', 'Download hist'),
                           
                           h5("Dunn (1964) Kruskal-Wallis multiple comparison:"),
                           
                           plotOutput(outputId = "boxp"),
                           radioButtons("stats", "Statistical tests",
                                        choices = c("Kruskal-Wallis"),
                                        selected = "Kruskal-Wallis"),
                           
                           actionButton('stats', "Run Kruskal-Wallis"),
                           DT::dataTableOutput('stats3'),
                           textOutput("warn")
                                       
                           ),
                           tabPanel("Alignment",  value = 2)
                  )
      )
    )
  )
