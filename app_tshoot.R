library(shiny)
library(DT)
library(shinythemes)
library(dplyr)


ui <- navbarPage("Genome quality statics app",
       theme = shinytheme("spacelab"),
       # Include custom CSS in the head tag
       tags$head(
         tags$style(HTML("
           .big-font {
             font-size: 20px;
           }
           .sidebar {
             width: 20%;
           }
           .main {
             width: 80%;
             padding-right: 0;
           }
           .shiny-plot-output, .shiny-data-table-output {
             width: 100% !important;
           }
           .shiny-output-error {
             width: 100%;
           }
         "))
       ),
     tabPanel("Box plots",
                sidebarLayout(
                  sidebarPanel(
                    fileInput("file", "Choose TSV File", accept = c(".tsv")),
                    selectizeInput("box1", "Choose Columns for Boxplot (x-axis):", choices = NULL, multiple = TRUE, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    selectizeInput("box2", "Choose a Column for Boxplot (y-axis):", choices = NULL, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    tags$hr(), # Separator line
                    checkboxInput("unknown", HTML("Filter out genomes with <strong>unknown</strong> metadata annotation")),
                    tags$hr(), # Separator line
                    checkboxInput("filterNA", "Filter NAs", value = FALSE),
                    tags$hr(), # Separator line
                    # Add class to increase font size
                    tags$p(class="big-font", tags$b("To Filter:"), "Select factor, choose level, then hit", tags$b("Filter Button")),
                    tags$p(tags$b("Hitting the filter button is required to filter data")),
                    selectizeInput("factor", "Choose a Factor Column to Filter:", choices = NULL, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    selectizeInput("level", "Choose Level(s) to Filter:", choices = NULL, multiple = TRUE, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    actionButton("filter", strong("Filter")), # Updated line
                    tags$hr(), # Separator line
                    actionButton("reset", "Reset"),
                    width=2
                  ),                          
                  mainPanel(
                    plotOutput("boxPlot", height = "1000px"),
                    downloadButton("downloadBox", "Download Boxplot as PDF"), # add this line
                    dataTableOutput("anovaTable"),
                    downloadButton("downloadAnovaTSV", "Download ANOVA Table as TSV"),
                    dataTableOutput("tukeyTable"),
                    downloadButton("downloadTukeyTSV", "Download Tukey-HSD Table as TSV")  
                  )
                )
                 ),                  
                tabPanel("xy scatter plots",
                 sidebarLayout(
                   sidebarPanel(
                     selectizeInput("scatter_x", "Choose x-axis:", choices = NULL, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                     selectizeInput("scatter_y", "Choose y axis:", choices = NULL, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                     # actionButton("reset2", "Reset")
                     width=2
                   ),
                   mainPanel(
                     # plotOutput("scatterPlot", brush = brushOpts(id = "scatter_brush")), ## causing problems
                     plotOutput("scatterPlot", height = "700px", width = "95%"),
                     downloadButton("downloadScatter", "Download Boxplot as PDF"),
                     dataTableOutput("corrTable"),
                     downloadButton("downloadCorrelationTSV", "Download Correlation Table as TSV")  
                   )
                    )
                 ),
                 tabPanel("Data table",
                          dataTableOutput("dataTable"),
                          downloadButton("downloadDataTableTSV", "Download Data Table")  
                 )
)


server <- function(input, output, session) {
  data <- reactive({
    req(input$file)
    read.delim(input$file$datapath)
  })

  filtered_data <- reactive({
    req(data(), input$factor)
    df <- data()

    # Update the 'level' choices based on the 'factor' chosen
    updateSelectInput(session, "level", "Choose a level", choices = unique(df[[input$factor]]))

    # Filter data
    df %>% filter(.data[[input$factor]] %in% input$level) %>%
      filter(if (input$unknown) (.data[[input$level]] != "unknown") else TRUE)
  })

  output$dataTable <- renderDataTable({
    filtered_data()
  })

  output$downloadDataTableTSV <- downloadHandler(
    filename = function() {
      paste("final_data_table", ".tsv", sep = "")
    },
    content = function(file) {
      write.table(filtered_data(), file, sep = "\t", row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
