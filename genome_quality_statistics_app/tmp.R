## Notes: 
## Currently, turned off brush feature, as not working toegether will the rest of the functions
## Come back to this. 

library(shiny)
library(ggplot2)
library(broom)
library(dplyr)
library(DT)
library(shinythemes)
library(randomcoloR)
library(corrr)
library(tibble) # deframe()

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
       tags$script('
        Shiny.addCustomMessageHandler("resetBrush", function(variable) {
            var brushId = "#" + variable;
            var el = $(brushId);
            var brush = el.data("brush");
            if (brush) {
            brush.clear();
            el.trigger("change");
            }
        });
        '),
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
                     plotOutput("scatterPlot", height = "700px", width = "95%", brush = brushOpts(id = "scatter_brush")), ## causing problems
                    #  plotOutput("scatterPlot", height = "700px", width = "95%"),
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
  
  brushed_data <- reactiveVal()
  filtered_data <- reactiveVal(NULL)  # Initialize filtered_data as NULL
  
  ## READING INPUT FILE AND SETTING UP REACTIVES
  data <- reactive({
     req(input$file)  
     inFile <- input$file
     input.tmp <- read.delim(inFile$datapath)

    ## MAKING ALL COLUMNS W/ GTR 50 UNIQUE VALUES, CHARACTER AND LESS FACTOR
     char_cols <- sapply(input.tmp, function(x) if (is.character(x) | is.factor(x)) length(unique(x)) else NA)
     char_cols <- char_cols[!is.na(char_cols)]
     char_to_fac_colNames <- names(char_cols)[char_cols <= 50]
     input.tmp <- input.tmp %>%
       mutate_at(vars(char_to_fac_colNames), as.factor) %>%
       mutate_if(is.integer, as.numeric)

    ## UPDATE SELECT INPUT SECTION
      updateSelectInput(session, "factor", choices = c("", names(input.tmp[sapply(input.tmp, is.factor)])))
      updateSelectInput(session, "box1", choices = names(input.tmp[sapply(input.tmp, is.factor)]))
      updateSelectInput(session, "box2", choices = c("", names(input.tmp[sapply(input.tmp, is.numeric)])))
      updateSelectInput(session, "scatter_x", choices = c("", names(input.tmp[sapply(input.tmp, is.numeric)]))) 
      updateSelectInput(session, "scatter_y", choices = c("", names(input.tmp[sapply(input.tmp, is.numeric)])))

    return(input.tmp)
  })
  
  # Automatically set scatter_x to the second numeric column if available
  observe({
    req(data())
    num_cols <- names(data()[sapply(data(), is.numeric)])
    if (length(num_cols) > 1) {
      updateSelectInput(session, "scatter_x", selected = num_cols[2])
    }
  })

  # Update level choices when factor changes
  observeEvent(input$factor, {
    if (input$factor != "") {
      updateSelectInput(session, "level", choices = c("", levels(data()[[input$factor]])))
    }
  }, ignoreInit = TRUE)

  # Update scatter_y when box2 changes
  observeEvent(input$box2, {
    if (input$box2 != input$scatter_y) {
      updateSelectInput(session, "scatter_y", selected = input$box2) 
    }
  }, ignoreInit = TRUE) # 'ignoreInit' ensures that this won't trigger at app initialization

  # Update box2 when scatter_y changes
  observeEvent(input$scatter_y, {
    if (input$scatter_y != input$box2) {
      updateSelectInput(session, "box2", selected = input$scatter_y)
    }
  }, ignoreInit = TRUE) # 'ignoreInit' ensures that this won't trigger at app initialization

  reset_triggered <- reactiveVal(FALSE)

  observeEvent(input$reset, {
      # Reset UI components
      updateSelectInput(session, "factor", selected = "")
      updateSelectInput(session, "level", selected = "")
      updateSelectInput(session, "box1", selected = "")
      updateSelectInput(session, "box2", selected = "") 
      updateSelectInput(session, "scatter_y", selected = "") 
      # Reset reactive values
      filtered_data(NULL)  # Reset filtered_data
      brushed_data(NULL)   # Reset brushed_data
      reset_triggered(TRUE)
    # Clear the scatter plot brush
    session$sendCustomMessage(type = 'resetBrush', message = 'scatter_brush')

    })

  observeEvent(input$filter, {
    req(input$factor, input$level)

    # Prioritize brushed_data, then the main data
    data_to_filter <- if (!is.null(brushed_data())) brushed_data() else data()

    filtered <- data_to_filter %>%
      filter(.data[[input$factor]] %in% input$level) %>%
      filter(if(input$unknown) (.data[[input$level]] != "unknown") else TRUE)
       
    # Check if filtered data is not empty
    if (nrow(filtered) > 0) {
      filtered_data(filtered)
    } else {
      filtered_data(NULL)
       brushed_data(NULL) 
    }
  })
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Color palette
fixed_palette <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#DDAA77"))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Boxplot
boxplot_output <- reactive({
  req(input$box1, input$box2)
  # Prioritize brushed_data, then filtered_data, and finally the main data
  data_filtered <- if (!is.null(brushed_data())) brushed_data() else {
    if (!is.null(filtered_data())) filtered_data() else data()
  }

  # Filter out NAs if filterNA is checked
  if (input$filterNA) {
    data_filtered <- data_filtered[!is.na(data_filtered[[input$box2]]), ]
  }
  interaction_values <- interaction(data_filtered[, input$box1, drop = FALSE])
  
  # If input$unknown is true, filter out 'unknown' from interaction_values
  if(input$unknown) {
    data_filtered <- data_filtered[!grepl("unknown", interaction_values),]
    interaction_values <- interaction_values[!grepl("unknown", interaction_values)]
  }

  # Calculate the total count
  total_count <- nrow(data_filtered)
  # Calculate the count for each group
  counts <- data_filtered %>%
    group_by(group = interaction_values) %>%
    summarise(n = n())
  
  p <- ggplot(data_filtered, aes(x = interaction_values, 
                                 y = .data[[input$box2]], 
                                 fill = interaction_values,
                                 color = interaction_values
                                 )) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +
    ylab(input$box2) +
    ggtitle(paste("Filter Column:", input$factor, 
                  "\nFilter Level:", paste(input$level, collapse = ", "))) +
    guides(fill = guide_legend(title = NULL, ncol=1), 
           color = guide_legend(title = NULL, ncol=1)
          )   
  levels_count <- length(unique(interaction_values))

  if (levels_count > 1) {
    p <- p + scale_fill_manual(values = fixed_palette(levels_count)[1:levels_count], 
                              labels = function(x) {
                                sapply(x, function(.x) {
                                  count <- filter(counts, group == .x)$n
                                  paste(.x, "(n =", count, ")")
                                })}) + 
      scale_color_manual(values = fixed_palette(levels_count)[1:levels_count], 
                         labels = function(x) {
                           sapply(x, function(.x) {
                             count <- filter(counts, group == .x)$n
                             paste(.x, "(n =", count, ")")
                           })}) +
      annotate("text", x = Inf, y = Inf,
               label = paste("Total count =", total_count),
               hjust = 1, vjust = 1, size = 6, color = "black")
  }

  # Updating x-axis labels with count
  p <- p + scale_x_discrete(labels = function(x) {
    sapply(x, function(.x) {
      count <- filter(counts, group == .x)$n
      paste(.x, "(n =", count, ")")
    })
  })
  p
})


 output$boxPlot <- renderPlot({boxplot_output()})

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Scatter plot
  # Scatter plot
  scatterplot_output <- reactive({
    req(input$scatter_x, input$scatter_y)

    # Priority changes: if reset is triggered, avoid brushed_data
    data_filtered <- if (!reset_triggered() && !is.null(brushed_data())) {
      brushed_data()
    } else if (!is.null(filtered_data())) {
      filtered_data()
    } else {
      data()
    }
  
# Filter out NAs if filterNA is checked
if (input$filterNA) {
  data_filtered <- data_filtered[!is.na(data_filtered[[input$box2]]), ]
}
  
  if (is.null(input$box1)) {
    total_count <- nrow(data_filtered)
    p <- ggplot(data_filtered, aes_string(x=input$scatter_x, y=input$scatter_y)) +
      geom_point(size = 4, alpha=0.75) +
      ylab(input$scatter_y) +
      ggtitle(paste("Filter Column:", input$factor, 
                "\nFilter Level:", paste(input$level, collapse = ", "))) +
      annotate("text", x = Inf, y = Inf,
                label = paste("Total count =", total_count),
                hjust = 1, vjust = 1, size = 6, color = "black")
  } else {
    interaction_values <- interaction(data_filtered[, input$box1, drop = FALSE])
    
    # If input$unknown is true, filter out 'unknown' from interaction_values
    if(input$unknown) {
      data_filtered <- data_filtered[!grepl("unknown", interaction_values),]
      interaction_values <- interaction_values[!grepl("unknown", interaction_values)]
    }
    
    # Calculate the total count
    total_count <- nrow(data_filtered)
    
    # Calculate the count for each group
    counts <- data_filtered %>%
      group_by(group = interaction_values) %>%
      summarise(n = n()) %>%
      deframe() # Convert to named vector
      
    # Modify counts to include them in legend labels
    new_labels <- paste(names(counts), " (n =", counts, ")")
    
    p <- ggplot(data_filtered, aes_string(x=input$scatter_x, y=input$scatter_y)) +
      geom_point(size = 4, alpha=0.75, aes(color = interaction_values)) +
      ylab(input$scatter_y) +
      ggtitle(paste("Filter Column:", input$factor, 
                    "\nFilter Level:", paste(input$level, collapse = ", "))) +
      guides(color = guide_legend(width = unit(2, "in"), title = paste(input$box1, collapse=", "), ncol=1)) +
      annotate("text", x = Inf, y = Inf,
                label = paste("Total count =", total_count),
                hjust = 1, vjust = 1, size = 6, color = "black")
  }

  if (!is.null(input$box1) && input$box1 %in% names(data_filtered)) {
    interaction_values <- interaction(data_filtered[, input$box1, drop = FALSE])
    levels_count <- length(unique(interaction_values))
    
    if (levels_count > 1) {
      p <- p + scale_color_manual(values = fixed_palette(levels_count)[1:levels_count], labels = new_labels)
    }
  }

  p # using print(p), caused brush reload to fail
})


output$scatterPlot <- renderPlot({scatterplot_output()})

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Added
  # Brush reactive
  observeEvent(input$scatter_brush, {
    reset_triggered(FALSE)  # When brushing happens, reset is not triggered
    brush <- input$scatter_brush
    data_filtered <- if (!is.null(brushed_data())) brushed_data() else {
      if (!is.null(filtered_data())) filtered_data() else data()
    }
    if (!is.null(brush)) {
      brushed_data(brushedPoints(data_filtered, brush, xvar = input$scatter_x, yvar = input$scatter_y))
    } else {
      brushed_data(NULL)  # reset the brushed_data if brush is removed
    }
  })
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
## FILTERED DATA TABLE
output$dataTable <- renderDataTable({
  # Prioritize brushed_data, then filtered_data, and finally the main data
  data_to_display <- if (!is.null(brushed_data())) {
    brushed_data()
  } else if (!is.null(filtered_data())) {
    filtered_data()
  } else {
    data()
  }

  data_to_display
})

output$downloadDataTableTSV <- downloadHandler(
  filename = function() {
    paste("final_data_table", ".tsv", sep = "")
  },
  content = function(file) {
    write.table(if (!is.null(brushed_data())) brushed_data() else if (!is.null(filtered_data())) filtered_data() else data(), 
                file, 
                sep = "\t", 
                row.names = FALSE)
  }
)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
## RESET BUTTON

  # observeEvent(input$reset2, {
  #   updateSelectInput(session, "scatter_x", selected = "")
  #   updateSelectInput(session, "scatter_y", selected = "")
  #   brushed_data(NULL)
  # })

}

shinyApp(ui, server)
