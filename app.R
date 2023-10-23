## Notes: 
## Currently, brush selection tools, flip_xy, filterNA and Reset all (on both tabs) work
## - Only caveat now, Cannot iterively filter after using brush selection, can do this while only on boxplot page, but not after brush selection

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
       id = "main_tabs",  # Setting the ID for the navbarPage
       theme = shinytheme("spacelab"),
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
                    actionButton("sort_box", "Sort Boxplot"),
                    tags$hr(), # Separator line
                    actionButton("flip_coord", "Flip boxplot xy axis"),
                    tags$hr(), # Separator line
                    tags$p(class="big-font", tags$b("To Filter:"), "Select factor, choose level, then hit", tags$b("Filter Button")),
                    tags$p(tags$b("Hitting the filter button is required to filter data")),
                    selectizeInput("filterFactor", "Choose a Factor Column to Filter:", choices = NULL, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    selectizeInput("level", "Choose Level(s) to Filter:", choices = NULL, multiple = TRUE, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    tags$hr(), # Separator line
                    tags$p(class="big-font", tags$b("Secondary Filter:"), "Select factor, choose level"),
                    selectizeInput("filterFactor2", "Choose a Secondary Factor Column to Filter:", choices = NULL, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    selectizeInput("level2", "Choose Level(s) for Secondary Filter:", choices = NULL, multiple = TRUE, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    tags$hr(), # Separator line
                    tags$p(class="big-font", tags$b("Third Filter:"), "Select factor, choose level"),
                    selectizeInput("filterFactor3", "Choose a Third Factor Column to Filter:", choices = NULL, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    selectizeInput("level3", "Choose Level(s) for Third Filter:", choices = NULL, multiple = TRUE, options = list("actions-box" = TRUE, "live-search" = TRUE)),
                    tags$hr(), # Separator line
                    actionButton("filter", "Filter"),
                    tags$hr(), # Separator line
                    actionButton("reset", "Reset"),
                    width=3
                  ),                          
                  mainPanel(
                    plotOutput("boxPlot", height = "1000px"),
                    downloadButton("downloadBox", "Download Boxplot as PDF"),
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
                    tags$hr(), # Separator line
                    actionButton("reset_scatter", "Reset"),

                     width=3
                   ),
                   mainPanel(
                     plotOutput("scatterPlot", height = "700px", width = "95%", brush = brushOpts(id = "scatter_brush")),
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
  
  values <- reactiveValues(brushed_points = NULL)
  brushed_data <- reactiveVal()
  filtered_data <- reactiveVal(NULL)  # Initialize filtered_data as NULL
  sort_boxplot <- reactiveVal(FALSE)
  flip_xy <- reactiveVal(FALSE)

  ## READING INPUT FILE AND SETUP DATA TABLE REACTIVE
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
      updateSelectInput(session, "filterFactor", choices = c("", names(input.tmp[sapply(input.tmp, is.factor)])))
      updateSelectInput(session, "filterFactor2", choices = c("", names(input.tmp[sapply(input.tmp, is.factor)])))
      updateSelectInput(session, "filterFactor3", choices = c("", names(input.tmp[sapply(input.tmp, is.factor)])))
      updateSelectInput(session, "box1", choices = names(input.tmp[sapply(input.tmp, is.factor)]))
      updateSelectInput(session, "box2", choices = c("", names(input.tmp[sapply(input.tmp, is.numeric)])))
      updateSelectInput(session, "scatter_x", choices = c("", names(input.tmp[sapply(input.tmp, is.numeric)]))) 
      updateSelectInput(session, "scatter_y", choices = c("", names(input.tmp[sapply(input.tmp, is.numeric)])))

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#        filtered_data(input.tmp) 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


    return(input.tmp)
  }) 
  # Update the level and level2 dropdowns based on the selected factor
  observe({
      if(!is.null(input$filteredFactor) && input$filteredFactor %in% colnames(data())) {
          updateSelectizeInput(session, "level", choices = unique(data()[[input$filteredFactor]]))
      }
      if(!is.null(input$filteredFactor2) && input$filteredFactor2 %in% colnames(data())) {
          updateSelectizeInput(session, "level2", choices = unique(data()[[input$filteredFactor2]]))
      }
  })
  # Automatically set scatter_x to the second numeric column if available
  observe({
    req(data())
    num_cols <- names(data()[sapply(data(), is.numeric)])
    if (length(num_cols) > 1) {
      updateSelectInput(session, "scatter_x", selected = num_cols[2])
    }
  })
  

  # Update level choices when filterFactor changes
  observeEvent(input$filterFactor, {
    if (input$filterFactor != "") {
      updateSelectInput(session, "level", choices = c("", levels(data()[[input$filterFactor]])))
    }
  }, ignoreInit = TRUE)
  # Update level choices when filterFactor changes
  observeEvent(input$filterFactor2, {
    if (input$filterFactor2 != "") {
      updateSelectInput(session, "level2", choices = c("", levels(data()[[input$filterFactor2]])))
    }
  }, ignoreInit = TRUE)
  # Update level choices when filterFactor changes
  observeEvent(input$filterFactor3, {
    if (input$filterFactor3 != "") {
      updateSelectInput(session, "level3", choices = c("", levels(data()[[input$filterFactor3]])))
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
  # Sorting box plots
  observeEvent(input$sort_box, {
    sort_boxplot(!sort_boxplot())
  })
  # Flip boxplot coords
  observeEvent(input$flip_coord, {
    flip_xy(!flip_xy())
  })

 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  ## COULD POSSIBLE MOVE THIS TO TOP OF SERVER WITH OTHER REACTIVEVALS
  reset_triggered <- reactiveVal(FALSE)
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
 
  observeEvent(input$reset, {
      # Reset UI components
      updateSelectInput(session, "factor", selected = "")
      updateSelectInput(session, "level", selected = "")
      updateSelectInput(session, "level2", selected = "")
      updateSelectInput(session, "box1", selected = "")
      updateSelectInput(session, "box2", selected = "") 
      updateSelectInput(session, "scatter_y", selected = "")
      updateSelectInput(session, "filterFactor", selected = "")
      updateSelectInput(session, "filterFactor2", selected = "")
      updateSelectInput(session, "filterFactor3", selected = "") 
      # Reset reactive values
      filtered_data(NULL)  # Reset filtered_data
      brushed_data(NULL)   # Reset brushed_data
      reset_triggered(TRUE)
      sort_boxplot(FALSE)  # Reset sort_boxplot
      flip_xy(FALSE)       # Reset flip_xy

      # Clear the scatter plot brush
      session$sendCustomMessage(type = 'resetBrush', message = 'scatter_brush')

    })

  observeEvent(input$reset_scatter, {
      # Reset UI components
      updateSelectInput(session, "factor", selected = "")
      updateSelectInput(session, "level", selected = "")
      updateSelectInput(session, "level2", selected = "")
      updateSelectInput(session, "box1", selected = "")
      updateSelectInput(session, "box2", selected = "") 
      updateSelectInput(session, "scatter_y", selected = "")
      updateSelectInput(session, "filterFactor", selected = "")
      updateSelectInput(session, "filterFactor2", selected = "")
      updateSelectInput(session, "filterFactor3", selected = "") 
      # Reset reactive values
      filtered_data(NULL)  # Reset filtered_data
      brushed_data(NULL)   # Reset brushed_data
      reset_triggered(TRUE)
      sort_boxplot(FALSE)  # Reset sort_boxplot
      flip_xy(FALSE)       # Reset flip_xy

      # Clear the scatter plot brush
      session$sendCustomMessage(type = 'resetBrush', message = 'scatter_brush')

    })



  # Clear brush on tab change
observe({
    # Assuming your navbarPage's input ID is 'tabs'
    currentTab <- input$main_tabs  
    if(currentTab != "xy scatter plots") {  # Check the exact name of your scatter plot tab
        session$sendCustomMessage(type = 'resetBrush', message = 'scatter_brush')
        # Also reset any reactive values holding brushed data or state
        values$brushed_points <- NULL
    }
})

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#   observe({
#     brush <- input$scatter_brush
#     if (!is.null(brush)) {
#       # Code to filter data() based on the brush, and assign it to filtered_data
#       brushed_data <- # ... some operation on data() using brush ...
#       filtered_data(brushed_data)
#     }
#   })

#   observe({
#     # If there are additional filters, apply them to filtered_data()
#     # This could depend on various input elements in your UI, not just one
#     if (!is.null(input$some_filter)) {
#       # Code to apply additional filters
#       further_filtered_data <- # ... some operation on filtered_data() using input$some_filter ...
#       filtered_data(further_filtered_data)
#     }
#   })
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 



  observeEvent(input$filter, {
    req(input$filterFactor, input$level)
    filtered <- data()
    if (!is.null(input$filterFactor) && input$filterFactor != "") {
      filtered <- filtered %>%
        filter(.data[[input$filterFactor]] %in% input$level)
    }
    if (!is.null(input$filterFactor2) && input$filterFactor2 != "" && !is.null(input$level2) && length(input$level2) > 0) {
      filtered <- filtered %>%
        filter(.data[[input$filterFactor2]] %in% input$level2)
    }
    if (!is.null(input$filterFactor3) && input$filterFactor3 != "" && !is.null(input$level3) && length(input$level3) > 0) {
      filtered <- filtered %>%
        filter(.data[[input$filterFactor3]] %in% input$level3)
    }
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
  # Initial data filtering
  data_filtered <- if (!is.null(brushed_data())) {
    brushed_data()
  } else if (!is.null(filtered_data())) {
    filtered_data()
  } else {
    data()
  }
  # Handle interaction_values
  if (length(input$box1) > 1) {
    interaction_values <- interaction(data_filtered[, input$box1, drop = FALSE])
  } else {
    interaction_values <- data_filtered[[input$box1]]
  }
  # Filter out NAs if filterNA is checked
  if (input$filterNA) {
    na_indices <- !is.na(data_filtered[[input$box2]])
    data_filtered <- data_filtered[na_indices, ]
    # Make sure interaction_values and data_filtered still align
    # Correcting the subsetting here based on the structure of interaction_values
    if (is.data.frame(interaction_values)) {
      interaction_values <- interaction_values[na_indices, ]
    } else {
      interaction_values <- interaction_values[na_indices]
    }
  }
  # Filter out 'unknown' if input$unknown is true
  if(input$unknown) {
    unknown_indices <- !grepl("unknown", interaction_values)
    data_filtered <- data_filtered[unknown_indices, ]
    interaction_values <- interaction_values[unknown_indices]
  }
  # Sorting boxplot and handling interaction_values
  if (sort_boxplot()) {
    # Handling multiple columns or a single column for interaction_values
    if (length(input$box1) > 1) {
      data_filtered <- data_filtered %>%
        rowwise() %>%
        mutate(interaction_values = paste(c_across(all_of(input$box1)), collapse = "_")) %>%
        ungroup()
    } else {
      data_filtered$interaction_values <- as.character(data_filtered[[input$box1]])
      data_filtered$interaction_values <- factor(data_filtered$interaction_values,
      # original_levels should be defined to represent your initial, unsorted state. 
                                                levels = original_levels)
    }
    # Make interaction_values a factor with proper levels
    data_filtered$interaction_values <- factor(data_filtered$interaction_values, levels = unique(data_filtered$interaction_values))
    # Calculate the mean for input$box2 for each group in 'interaction_values'
    means <- data_filtered %>%
      group_by(interaction_values) %>%
      summarise(mean_val = mean(.data[[input$box2]], na.rm = TRUE)) %>%
      arrange(desc(mean_val))
    # Reorder the levels based on the calculated means
    if (nrow(means) > 0 && !all(is.na(means$mean_val))) {
      data_filtered$interaction_values <- factor(data_filtered$interaction_values, levels = means$interaction_values)
    }
  }

  # Calculate the total count
  total_count <- nrow(data_filtered)
  # Calculate count for each group / interaction values
  counts <- data_filtered %>%
    group_by(group = interaction_values) %>%
    summarise(n = n())

  ## Note: should streamline, throwing all unchanging themes outside ifelse statement
  ## Will reduce lines of code, though had trouble when first attempted
  if(flip_xy()) {
    p <- ggplot(data_filtered, aes(x = forcats::fct_rev(interaction_values), 
                                  y = .data[[input$box2]], 
                                  fill = interaction_values,
                                  color = interaction_values
                                  )) +
       xlab(input$box2) +
       coord_flip() +
       geom_boxplot(alpha = 0.5) +     
       geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +  
       ggtitle(paste("Filter Column:", input$factor, 
                     "\nFilter Level:", paste(input$level, collapse = ", "))) +
       guides(fill = guide_legend(title = NULL, ncol=1), 
             color = guide_legend(title = NULL, ncol=1)
             ) +
       theme_linedraw() + 
       theme(axis.text.x = element_text(size=16, color="black"),
             axis.text.y = element_text(size=14, color="black"),
             axis.title.y = element_blank(),
             axis.title.x = element_text(size=16, color="black", vjust=-1),
             legend.title = element_text(size = 14, color="black"),
             legend.text = element_text(size=12, color="black"))       
  } else {
    p <- ggplot(data_filtered, aes(x = interaction_values, 
                                  y = .data[[input$box2]], 
                                  fill = interaction_values,
                                  color = interaction_values
                                  )) +
        ylab(input$box2) +
        geom_boxplot(alpha = 0.5) +
        geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +  
        ggtitle(paste("Filter Column:", input$factor, 
                      "\nFilter Level:", paste(input$level, collapse = ", "))) +
        guides(fill = guide_legend(title = NULL, ncol=1), 
              color = guide_legend(title = NULL, ncol=1)
              ) +
        theme_linedraw() + 
        theme(axis.text.x = element_text(size=16, color="black", angle=90),
              axis.text.y = element_text(size=14, color="black"),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=16, color="black", vjust=2.5),
              legend.title = element_text(size = 14, color="black"),
              legend.text = element_text(size=12, color="black")) 
  }
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
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
    ## Same as above for boxplot, can probably reduce, but req to maintain colors btwn box and scatter
    # Sort the boxplot based on the mean of each group if sort_boxplot() is TRUE
    if (sort_boxplot()) {
      # If 'input$box1' has more than one column selected, concatenate their values; otherwise, use the selected column directly.
      if (length(input$box1) > 1) {
        data_filtered <- data_filtered %>%
          rowwise() %>%
          mutate(interaction_values = paste(c_across(all_of(input$box1)), collapse = "_")) %>%
          ungroup()
      } else {
        # Directly use the column from 'input$box1' as 'interaction_values'.
        data_filtered$interaction_values <- as.character(data_filtered[[input$box1]])
      }
      # Convert 'interaction_values' to a factor and set its levels.
      data_filtered$interaction_values <- factor(data_filtered$interaction_values, levels = unique(data_filtered$interaction_values))
      # Debugging: Print the levels to check their consistency.
      print(levels(data_filtered$interaction_values))
      # Calculate the mean for input$box2 for each group in 'interaction_values'.
      means <- data_filtered %>%
        group_by(interaction_values) %>%
        summarise(mean_val = mean(.data[[input$box2]], na.rm = TRUE)) %>%
        arrange(desc(mean_val))
      # Reorder the levels based on the calculated means.
      if (nrow(means) > 0 && !all(is.na(means$mean_val))) {
        data_filtered$interaction_values <- factor(data_filtered$interaction_values, levels = means$interaction_values)
      }
    }
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##      
    if (is.null(input$box1)) {
        total_count <- nrow(data_filtered)
        p <- ggplot(data_filtered, aes_string(x=input$scatter_x, y=input$scatter_y)) +
            geom_point(size = 4, alpha=0.75) +
            ylab(input$scatter_y) +
            ggtitle(paste("Filter Column:", input$factor, 
                          "\nFilter Level:", paste(input$level, collapse = ", "))) +
            theme_linedraw() + 
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
        if (input$unknown) {
            data_filtered <- data_filtered[!grepl("unknown", interaction_values),]
            interaction_values <- interaction_values[!grepl("unknown", interaction_values)]
        } 
        total_count <- nrow(data_filtered)
        counts <- data_filtered %>%
            group_by(group = interaction_values) %>%
            summarise(n = n()) %>%
            deframe() # Convert to named vector
        new_labels <- paste(names(counts), " (n =", counts, ")")
        p <- ggplot(data_filtered, aes_string(x=input$scatter_x, y=input$scatter_y)) +
            geom_point(size = 4, alpha=0.75, aes(color = interaction_values)) +
            ylab(input$scatter_y) +
            ggtitle(paste("Filter Column:", input$factor, 
                          "\nFilter Level:", paste(input$level, collapse = ", "))) +
            guides(color = guide_legend(width = unit(2, "in"), title = paste(input$box1, collapse=", "), ncol=1)) +
            theme_linedraw() + 
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
    if (!is.null(input$box1) && any(input$box1 %in% names(data_filtered))) {
        interaction_values <- interaction(data_filtered[, input$box1, drop = FALSE])
        levels_count <- length(unique(interaction_values))
        if (levels_count > 1) {
            p <- p + scale_color_manual(values = fixed_palette(levels_count)[1:levels_count], labels = new_labels)
        }
    }
    p # using print(p), caused brush reload to fail
})

output$scatterPlot <- renderPlot({scatterplot_output()})
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#   output$scatterPlot <- renderPlot({
#     # Check if filtered_data is available, else use data
#     data_to_use <- if (!is.null(filtered_data())) filtered_data() else data()
#     # Code to render the plot using data_to_use
#   })
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 



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
  req(length(input$scatter_x) == 1, length(input$scatter_y) == 1)
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
  if (is.null(input$box1) || length(input$box1) == 0) {
    corr_value <- round(cor(selected_cols[[input$scatter_x]], selected_cols[[input$scatter_y]], use = "complete.obs"), 2)
    return(data.frame(Pair = "Total", Correlation = corr_value))
  }
  # Create interaction values
  interaction_values <- interaction(data_filtered[, input$box1], drop = TRUE, lex.order = TRUE)
  # If input$unknown is true, filter out 'unknown' from interaction_values
  if(input$unknown) {
    data_filtered <- data_filtered[!grepl("unknown", interaction_values),]
    interaction_values <- interaction_values[!grepl("unknown", interaction_values)]
  }
  combined_factor <- as.character(interaction_values)
  levels <- unique(combined_factor)
  # If there's only one unique level in combined_factor, return only the total correlation
  if (length(levels) == 1) {
    corr_value <- round(cor(selected_cols[[input$scatter_x]], selected_cols[[input$scatter_y]], use = "complete.obs"), 2)
    return(data.frame(Pair = "Total", Correlation = corr_value))
  }
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
  # If filterNA is checked, remove rows with NA values in the box2 column
  if (input$filterNA && !is.null(input$box2)) {
    data_to_display <- data_to_display[!is.na(data_to_display[[input$box2]]), ]
  }
  # If unknown is checked, filter out 'unknown' from factor columns
  if (input$unknown) {
    factor_cols <- names(data_to_display[sapply(data_to_display, is.factor)])
    for (col in factor_cols) {
      data_to_display <- data_to_display[!data_to_display[[col]] %in% c("unknown"), ]
    }
  }
  data_to_display
})
output$downloadDataTableTSV <- downloadHandler(
  filename = function() {
    paste("final_data_table", ".tsv", sep = "")
  },
  content = function(file) {
    data_to_save <- if (!is.null(brushed_data())) {
      brushed_data()
    } else if (!is.null(filtered_data())) {
      filtered_data()
    } else {
      data()
    }
    # If filterNA is checked, remove rows with NA values in the box2 column
    if (input$filterNA && !is.null(input$box2)) {
      data_to_save <- data_to_save[!is.na(data_to_save[[input$box2]]), ]
    }
    # If unknown is checked, filter out 'unknown' from factor columns
    if (input$unknown) {
      factor_cols <- names(data_to_save[sapply(data_to_save, is.factor)])
      for (col in factor_cols) {
        data_to_save <- data_to_save[!data_to_save[[col]] %in% c("unknown"), ]
      }
    }
    write.table(data_to_save, file, sep = "\t", row.names = FALSE)
  }
)



}

shinyApp(ui, server)
