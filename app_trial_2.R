library(shiny)
library(ggplot2)
library(broom)
library(dplyr)

ui <- fluidPage(
  fileInput("file", "Choose TSV File", accept = c(".tsv")),  # File input
  selectInput("factor", "Choose a Factor Column to Filter:", choices = NULL),
  selectInput("level", "Choose a Level to Filter:", choices = NULL),
  selectizeInput("box1", "Choose Columns for Boxplot (x-axis):", choices = NULL, multiple = TRUE),
  selectInput("box2", "Choose a Column for Boxplot (y-axis):", choices = NULL),
  plotOutput("boxPlot"),
  tableOutput("anovaTable"),
  tableOutput("tukeyTable")
)

server <- function(input, output, session) {
  
  # Reactive expression to read uploaded file and process the data
  data <- reactive({
    req(input$file)  # ensure the file is available
    
    inFile <- input$file
    input.tmp <- read.delim(inFile$datapath)
    
    # Process the data
    char_cols <- sapply(input.tmp, function(x) if (is.character(x) | is.factor(x)) length(unique(x)) else NA)
    char_cols <- char_cols[!is.na(char_cols)]
    char_to_fac_colNames <- names(char_cols)[char_cols <= 50]
    input.tmp <- input.tmp %>%
      mutate_at(vars(char_to_fac_colNames), as.factor) %>%
      mutate_if(is.integer, as.numeric)
    
    # Update UI based on uploaded file
    updateSelectInput(session, "factor", choices = c("", names(input.tmp[sapply(input.tmp, is.factor)])))
    updateSelectInput(session, "box1", choices = names(input.tmp[sapply(input.tmp, is.factor)]))
    updateSelectInput(session, "box2", choices = names(input.tmp[sapply(input.tmp, is.numeric)]))
    
    return(input.tmp)
  })
  
  observe({
    updateSelectInput(session, "level", choices = c("", levels(data()[[input$factor]])))
  })
  
  output$boxPlot <- renderPlot({
    req(input$box1, input$box2)  
    
    data_filtered <- data()[data()[[input$factor]] == input$level, ]
    
    if(input$factor == "" || input$level == ""){
      data_filtered <- data()
    }
    
    p <- ggplot(data_filtered, aes(x = interaction(data_filtered[, input$box1, drop = FALSE]), 
                                   y = data_filtered[[input$box2]], 
                                   fill = interaction(data_filtered[, input$box1, drop = FALSE]))) +
      geom_boxplot() +
      theme_minimal()
    
    levels_count <- length(unique(interaction(data_filtered[, input$box1, drop = FALSE])))
    
    if (levels_count > 1) {
      p <- p + scale_fill_manual(values = distinctColorPalette(levels_count))
    }
    
    p
  })
  
  output$anovaTable <- renderTable({
    req(input$box1, input$box2)  
    
    data_filtered <- data()[data()[[input$factor]] == input$level, ]
    if(input$factor == "" || input$level == ""){
      data_filtered <- data()
    }
    anova_res <- aov(as.formula(paste(input$box2, "~", paste("interaction(", paste(input$box1, collapse=","), ")", sep=""))), 
                     data = data_filtered)
    tidy(anova_res)  
  })
  
  output$tukeyTable <- renderTable({
    req(input$box1, input$box2)  
    
    data_filtered <- data()[data()[[input$factor]] == input$level, ]
    if(input$factor == "" || input$level == ""){
      data_filtered <- data()
    }
    tukey_res <- TukeyHSD(aov(as.formula(paste(input$box2, "~", paste("interaction(", paste(input$box1, collapse=","), ")", sep=""))), 
                              data = data_filtered))
    tidy(tukey_res)  
  })
}

shinyApp(ui, server)
