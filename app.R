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
    print(input$factor)
    print(input$level)
    req(input$factor, input$level)

    filtered <- data() %>%
     filter(.data[[input$factor]] %in% input$level) %>%
     filter(if(input$unknown) (.data[[input$level]] != "unknown") else TRUE)
       
   # Check if filtered data is not empty
   if (nrow(filtered) > 0) {
     filtered_data(filtered)
   } else {
     filtered_data(NULL)
   }
  })

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Color palette
fixed_palette <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#DDAA77"))
# fixed_palette <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", 
#   "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", 
#   "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", 
#   "#DDAA77", "#771122", "#AA4455", "#DD7788")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Boxplot
boxplot_output <- reactive({
  req(input$box1, input$box2)
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
          ) +
    theme_classic() + 
    theme(axis.text.x = element_text(size=16, color="black", angle=90),
          axis.text.y = element_text(size=14, color="black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16, color="black"),
          legend.title = element_text(size = 14, color="black"),
          legend.text = element_text(size=12, color="black"))
  
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
 output$downloadBox <- downloadHandler(
   filename = function() {
     paste("boxplot", ".pdf", sep="")
   },
   content = function(file) {
     ggsave(file, plot = boxplot_output(), device = "pdf", width = 10, height = 7)
   }
 )

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
## ANOVA TABLE
generate_anova_table <- reactive({
  req(input$box1, input$box2)  

  # Use filtered_data if available, else use the original data
  data_filtered <- if (!is.null(filtered_data())) filtered_data() else data() 
  
  anova_res <- aov(as.formula(paste(input$box2, "~", 
                                    paste("interaction(", paste(input$box1, collapse=","), ")", sep=""))), 
                   data = data_filtered)
  anova_tidy <- tidy(anova_res)
  anova_tidy <- data.frame(lapply(anova_tidy, function(x) 
    if(is.numeric(x)) format(round(x, 2), nsmall = 2) else x))
  anova_tidy
})

output$anovaTable <- renderDataTable({
  datatable(generate_anova_table(), caption = "ANOVA")
})

output$downloadAnovaTSV <- downloadHandler(
  filename = function() {
    paste("anova_table", ".tsv", sep = "")
  },
  content = function(file) {
    write.table(generate_anova_table(), file, sep = "\t", row.names = FALSE)
  }
)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
## TUKEY HSD TABLE
generate_tukey_table <- reactive({
  req(input$box1, input$box2)  

  # Use filtered_data if available, else use the original data
  data_filtered <- if (!is.null(filtered_data())) filtered_data() else data() 
  tukey_res <- TukeyHSD(aov(as.formula(paste(input$box2, "~", 
                                             paste("interaction(", paste(input$box1, collapse=","), ")", sep=""))), 
                             data = data_filtered))
  tukey_tidy <- tidy(tukey_res)
  tukey_tidy <- data.frame(lapply(tukey_tidy, function(x) 
    if(is.numeric(x)) format(round(x, 2), nsmall = 2) else x))
  tukey_tidy
})

output$tukeyTable <- renderDataTable({
  datatable(generate_tukey_table(), caption = "Tukey-HSD")
})

output$downloadTukeyTSV <- downloadHandler(
  filename = function() {
    paste("tukey_table", ".tsv", sep = "")
  },
  content = function(file) {
    write.table(generate_tukey_table(), file, sep = "\t", row.names = FALSE)
  }
)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
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
      theme_classic() + 
      theme(axis.text.x = element_text(size=16, color="black"),
            axis.text.y = element_text(size=14, color="black"),
            axis.title.x = element_text(size=16, color="black"),
            axis.title.y = element_text(size=16, color="black"),
            legend.title = element_text(size = 14, color="black"),
            legend.text = element_text(size=12, color="black"),
            plot.margin = unit(c(1, 1, 1, 1), "cm")) +
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
      theme_classic() + 
      theme(axis.text.x = element_text(size=16, color="black"),
            axis.text.y = element_text(size=14, color="black"),
            axis.title.x = element_text(size=16, color="black"),
            axis.title.y = element_text(size=16, color="black"),
            legend.title = element_text(size = 14, color="black"),
            legend.text = element_text(size=12, color="black"),
            plot.margin = unit(c(1, 1, 1, 1), "cm")) +
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


  #   guides(color = guide_legend(width = unit(2, "in"), title = paste(input$box1, collapse=", "), ncol=1)) +
  #   theme_classic() + 
  #   theme(axis.text.x = element_text(size=16, color="black"),
  #         axis.text.y = element_text(size=14, color="black"),
  #         axis.title.x = element_text(size=16, color="black"),
  #         axis.title.y = element_text(size=16, color="black"),
  #         legend.title = element_text(size = 14, color="black"),
  #         legend.text = element_text(size=12, color="black"),
  #         plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  #   annotate("text", x = Inf, y = Inf,
  #              label = paste("Total count =", total_count),
  #              hjust = 1, vjust = 1, size = 6, color = "black")



output$scatterPlot <- renderPlot({scatterplot_output()})
output$downloadScatter <- downloadHandler(
  filename = function() {
    paste("scatterplot", ".pdf", sep="")
  },
  content = function(file) {
    ggsave(file, plot = scatterplot_output(), device = "pdf", width = 10, height = 7)
  }
)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CORRELATION TABLE
calculate_pairwise_correlations <- reactive({
  req(input$scatter_x, input$scatter_y)
  
  if (input$scatter_x == input$scatter_y) {
    stop("Correlation Table requires different variables for X and Y axes.")
  }
  
  data_filtered <- if (!is.null(brushed_data())) {
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


  selected_cols <- select(data_filtered, all_of(c(input$scatter_x, input$scatter_y)))
  selected_cols <- as.data.frame(selected_cols)
  
  # Check for NAs in selected columns
  if (any(is.na(selected_cols[[input$scatter_x]])) || any(is.na(selected_cols[[input$scatter_y]]))) {
    return(data.frame(Pair = "Warning", Correlation = "Please select the filterNA filter to remove NAs on the Box plots tab. NAs interfere with the correlation table."))
  }

  # If box1 is not selected, compute a simple correlation
  if (is.null(input$box1) || !(input$box1 %in% names(data_filtered))) {
    corr_value <- round(cor(selected_cols[[input$scatter_x]], selected_cols[[input$scatter_y]], use = "complete.obs"), 2)
    return(data.frame(Pair = "Total", Correlation = corr_value))
  }



  interaction_values <- interaction(data_filtered[, input$box1, drop = FALSE])
  
  # If input$unknown is true, filter out 'unknown' from interaction_values
  if(input$unknown) {
    data_filtered <- data_filtered[!grepl("unknown", interaction_values),]
    interaction_values <- interaction_values[!grepl("unknown", interaction_values)]
  }

  factor_cols <- data_filtered[, input$box1, drop = FALSE]
  combined_factor <- apply(factor_cols, 1, function(row) {
    paste(row, collapse = "_")
  })
  levels <- unique(combined_factor)
  corr_list <- list()
  
  # Calculating correlations for each pair
  for (i in 1:(length(levels) - 1)) {
    for (j in (i + 1):length(levels)) {
      level_i_data <- selected_cols[combined_factor == levels[i], ]
      level_j_data <- selected_cols[combined_factor == levels[j], ]
      combined_data <- rbind(level_i_data, level_j_data)
      corr_value <- round(cor(combined_data[[input$scatter_x]], combined_data[[input$scatter_y]], use = "complete.obs"), 2)
      corr_list[[paste(levels[i], "vs", levels[j])]] <- corr_value
    }
  } 
  
  # Add total correlation (average in this case)
  total_corr <- round(mean(unlist(corr_list)), 2)
  total_df <- data.frame(Pair = "Total", Correlation = total_corr)
  corr_df <- data.frame(Pair = names(corr_list), Correlation = unlist(corr_list))
  
  # Place Total at the top of the data frame
  corr_df <- rbind(total_df, corr_df)
  
  return(corr_df)
})


output$corrTable <- renderDataTable({
  corr_df <- calculate_pairwise_correlations()
  # Use filtered_data if available, else use the original data
  data_filtered <- if (!is.null(filtered_data())) filtered_data() else data() 
  rownames(corr_df) <- NULL  
  datatable(corr_df, options = list(pageLength = 20), 
            caption = "Pairwise Correlations", rownames = FALSE)
})

output$downloadCorrelationTSV <- downloadHandler(
  filename = function() {
    paste("correlation_table", ".tsv", sep = "")
  },
  content = function(file) {
    write.table(calculate_pairwise_correlations(), file, sep = "\t", row.names = FALSE)
  }
)

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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
## FILTERED DATA TABLE
output$dataTable <- renderDataTable({
  # First, check for brushed data. If it's NULL, then revert to filtered_data, 
  # and if that is also NULL, use the main data.
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
    write.table(if (!is.null(brushed_data())) brushed_data() else data(), 
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
