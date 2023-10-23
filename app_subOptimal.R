library(shiny)
library(ggplot2)
library(DT)
library(shinyjs)  # Load the shinyjs library
library(broom)

# UI code
ui <- fluidPage(
    useShinyjs(),  # Initialize shinyjs
    
    titlePanel("Genome statistics explorer"),
    
    # Define the sidebar layout
    sidebarLayout(
        sidebarPanel(
            # Add a file input control to upload a TSV file
            fileInput("file", "Upload a TSV file"),
            
            # Add a dropdown menu to select primary grouping variable
            selectInput("group_var_primary", "Select Primary Grouping Variable", choices = character(0)),
            # Add a dropdown menu to select secondary grouping variable
            selectInput("group_var_secondary", "Select Secondary Grouping Variable", choices = character(0)),
            # Add a dropdown menu to select primary numeric variable
            selectInput("numeric_var_primary", "Select Primary Numeric Variable", choices = character(0)),
            # Add a dropdown menu to select secondary numeric variable
            selectInput("numeric_var_secondary", "Select Secondary Numeric Variable", choices = character(0)),
            
            # Add a "Reset All" button
            # actionButton("reset_button", "Reset All")
            
            # Add any additional input controls here if needed
        ),
        
        mainPanel(
            # Create a tabset panel with three tabs
            tabsetPanel(
                tabPanel("Plots",
                    plotOutput("boxplot"),
                    plotOutput("histogram"),
                    plotOutput("scatterplot")  # Add a plotOutput for the scatterplot
                ),
                tabPanel("Statistics",
                    textOutput("statistics"),
                    tableOutput("anova_table"),
                    br(),
                    dataTableOutput("tukey_table")
                ),
                tabPanel("Data Table",
                    DTOutput("datatable")
                )
            )
        )
    )
)

# Server code
server <- function(input, output, session) {
    # Initialize an empty data frame to hold the dataset
    sample_data <- data.frame()
    
    # Function to read the uploaded TSV file
    observe({
        req(input$file)
        
        # Read the uploaded file as a data frame
        df <- read.table(input$file$datapath, header = TRUE, sep = "\t")
        
        # Assign the uploaded data frame to the sample_data variable
        sample_data <<- df
        
        # Get the names of non-integer/numeric variables excluding 'genome_id'
        non_numeric_vars <- names(df)[sapply(names(df), function(x) !is.numeric(df[[x]]) && x != "genome_id")]
        
        # Get the names of numeric variables
        numeric_vars <- names(df)[sapply(df, is.numeric)]
        
        # Populate the dropdown menus
        updateSelectInput(session, "group_var_primary", choices = non_numeric_vars)
        updateSelectInput(session, "group_var_secondary", choices = non_numeric_vars)
        updateSelectInput(session, "numeric_var_primary", choices = numeric_vars)
        updateSelectInput(session, "numeric_var_secondary", choices = numeric_vars)
    })
    
    # Create boxplot
    output$boxplot <- renderPlot({
        # Create a ggplot boxplot based on the selected grouping variable and numeric variable
        group_var_primary <- input$group_var_primary
        numeric_var_primary <- input$numeric_var_primary
        
        # Check if both grouping and numeric variables are selected
        if (!is.null(group_var_primary) && !is.null(numeric_var_primary)) {
            # Load the ggplot2 library if it's not already loaded
            if (!require(ggplot2)) {
                install.packages("ggplot2")
                library(ggplot2)
            }
            
            # Create the boxplot using ggplot
            p_boxplot <- ggplot(sample_data, aes(x = sample_data[[group_var_primary]], y = sample_data[[numeric_var_primary]], fill = sample_data[[group_var_primary]])) +
                geom_boxplot() +
                labs(x = group_var_primary, y = numeric_var_primary, title = "Boxplot by Group") +
                theme_classic() + 
                theme(axis.text.x = element_text(size=16, color="black"),
                        axis.text.y = element_text(size=14, color="black"),
                        axis.title.x = element_text(size=14, color="black"),
                        axis.title.y = element_text(size=16, color="black"),
                        legend.title = element_text(size = 14, color="black"),
                        legend.text = element_text(size=12, color="black"))
            
            # Print the boxplot using print()
            print(p_boxplot)
        } else {
            # If no variables are selected, display a message
            plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), main = "Select grouping and numeric variables.")
        }
    })
    
    # Create histogram
    output$histogram <- renderPlot({
        # Create a ggplot histogram based on the selected grouping variable and numeric variable
        group_var_primary <- input$group_var_primary
        numeric_var_primary <- input$numeric_var_primary
        
        # Check if both grouping and numeric variables are selected
        if (!is.null(group_var_primary) && !is.null(numeric_var_primary)) {
            # Load the ggplot2 library if it's not already loaded
            if (!require(ggplot2)) {
                install.packages("ggplot2")
                library(ggplot2)
            }
            
            # Create the histogram using ggplot
            p_histogram <- ggplot(sample_data, aes(x = sample_data[[numeric_var_primary]], color = sample_data[[group_var_primary]], fill = sample_data[[group_var_primary]])) +
                geom_histogram() +
                labs(x = numeric_var_primary, y = "Frequency", title = "Histogram by Group") +
                theme_classic() + 
                theme(axis.text.x = element_text(size=16, color="black"),
                        axis.text.y = element_text(size=14, color="black"),
                        axis.title.x = element_text(size=14, color="black"),
                        axis.title.y = element_text(size=16, color="black"),
                        legend.title = element_text(size = 14, color="black"),
                        legend.text = element_text(size=12, color="black"))
            
            # Print the histogram using print()
            print(p_histogram)
        } else {
            # If no variables are selected, display a message
            plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), main = "Select grouping and numeric variables.")
        }
    })
    
    output$scatterplot <- renderPlot({
        # Create a ggplot scatterplot based on the selected numeric variables
        numeric_var_primary <- input$numeric_var_primary
        numeric_var_secondary <- input$numeric_var_secondary
        
        # Check if both numeric variables are selected
        if (!is.null(numeric_var_primary) && !is.null(numeric_var_secondary)) {
            # Load the ggplot2 library if it's not already loaded
            if (!require(ggplot2)) {
                install.packages("ggplot2")
                library(ggplot2)
            }
            
            # Create the scatterplot using ggplot
            p_scatterplot <- ggplot(sample_data, aes(x = sample_data[[numeric_var_primary]], y = sample_data[[numeric_var_secondary]])) +
                geom_point(aes(color = sample_data[[input$group_var_primary]])) +  # Add color based on primary grouping variable
                labs(x = numeric_var_primary, y = numeric_var_secondary, title = "Scatterplot") +
                theme_classic() +
                theme(axis.text.x = element_text(size=16, color="black"),
                        axis.text.y = element_text(size=14, color="black"),
                        axis.title.x = element_text(size=14, color="black"),
                        axis.title.y = element_text(size=16, color="black"),
                        legend.title = element_text(size = 14, color="black"),
                        legend.text = element_text(size=12, color="black"))
            
            # Print the scatterplot using print()
            print(p_scatterplot)
        } else {
            # If no numeric variables are selected, display a message
            plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), main = "Select numeric variables for scatterplot.")
        }
    })
    
# Create ANOVA table



output$anova_table <- renderTable({
    # Check if both grouping and numeric variables are selected
    if (!is.null(input$group_var_primary) && !is.null(input$numeric_var_primary)) {
        # Perform ANOVA and store the result
        anova_result <- tidy(aov(sample_data[[input$numeric_var_primary]] ~ sample_data[[input$group_var_primary]], data = sample_data))

        # Modify the column names to include the variable names
        colnames(anova_result) <- c(input$numeric_var_primary, "DF", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
        
        # Replace 'sample_data[[input$group_var_primary]]' in the table with the variable name
        # anova_result$`Source of Variation` <- input$group_var_primary

        anova_result[1,1] <- input$group_var_primary


        # Return the modified ANOVA table
        return(anova_result)
    }
    },
        striped = TRUE, bordered = TRUE, hover = TRUE, spacing = 'm',  width = '100%', 
        align = 'c', rownames = FALSE, caption="ANOVA", 
        caption.placement = getOption("xtable.caption.placement", "top"))

output$tukey_table <- renderDataTable({
    # Check if both grouping and numeric variables are selected
    if (!is.null(input$group_var_primary) && !is.null(input$numeric_var_primary)) {
        # Perform ANOVA and store the result
        tukey_result <- tidy(TukeyHSD(aov(sample_data[[input$numeric_var_primary]] ~ sample_data[[input$group_var_primary]], data = sample_data)))
        # tukey_result[,1] <- input$group_var_primary
        tukey_result <- tukey_result %>%
            select(contrast, adj.p.value)
        # Modify the column names
        colnames(tukey_result) <- c("comparison","adj.p.value")
       # Return the modified Tukey table
        return(tukey_result)
    }
    }, options=list(dom="t"), caption="Tukey-HSD")
        # striped = TRUE, bordered = TRUE, hover = TRUE, spacing = 'm',  width = '100%', 
        # align = 'c', rownames = FALSE, caption="", 
        # caption.placement = getOption("xtable.caption.placement", "top"), options=list(dom="t"))

# Create Tukey pairwise comparison table
# output$tukey_table <- renderDataTable({
#     # Check if both grouping and numeric variables are selected
#     if (!is.null(input$group_var_primary) && !is.null(input$numeric_var_primary)) {
#         # Perform Tukey pairwise comparisons and store the result
#         tukey_result <- TukeyHSD(aov(sample_data[[input$numeric_var_primary]] ~ sample_data[[input$group_var_primary]], data = sample_data))
        
#         # Create a table from Tukey results
#         as.data.frame(tukey_result$`sample_data[[input$group_var_primary]]`)
#     }
# }, 
#     caption="Tukey-HSD", options = list(dom = 't'))


    # Create data table
    output$datatable <- renderDT({
        datatable(sample_data)
    })
}

shinyApp(ui, server)
