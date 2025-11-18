# Outlier Detection Module

# UI function
outlierUI <- function(id) {
  ns <- NS(id)  # namespacing

  tagList(
    fluidRow(
      # Left column: Scrollable explanatory text
      column(
        width = 4,
        div(
          style = "max-height: 500px; overflow-y: auto; padding: 10px; background-color: #f9f9f9;
                   border-radius: 8px; border: 1px solid #ddd;",
          strong("Description"),
          p("Detects outlier genes present in a given tissue of interest, which have a large delta between their RNA expression and protein expression values."),
          p("The function, given a tissue of interest, compares the RNA expression and protein expression of all genes present in the given tissue, and identifies genes which have a large difference in their RNA and protein expression. These genes are identified as outlier genes, and are most likely candidates for post-transcriptional modifications. This function returns the list of outlier genes as well as a plot highlighting the outlier genes."),
          br(),
          strong("Input:"),
          p("The tissue of interest to analyze."),
          strong("Output:"),
          p("Returns a plot highlighting the outlier genes in red and a vector of outlier gene names detected."),
          br(),
          strong("Instructions:"),
          p("1. From the 'Select a tissue' toolbar, choose an input tissue to detect outliers in."),
          p("2. Once selected, click on 'Run Outlier Detection' to view results."),
          p("3. The output plot is under the heading 'Outlier Genes Plot'."),
          p("4. The output vector is under the heading 'Outlier Gene List'."),
          br(),
          strong("Run Examples:"),
          p("To run the example, follow the steps:"),
          p("1. Click on the 'Run Example' button"),
          p("2. Click on the 'Run Outlier Detection' button to view results."),
          br(),
          # Run Example button
          actionButton(ns("run_example"), "Run Example"),
          br(), br(),
          strong("References"),
          p(""),
          p(""),
          p(""),
          p(""),
          p(""),
          br(),
          strong("Explanation of plot"),
          p(""),
          br(),
          strong("Explanation of table"),
          p(""),
          br()

        )
      ),

      # Right column: Input controls and output
      column(
        width = 8,
        selectInput(
          inputId = ns("tissue"),
          label = "Select a tissue:",
          choices = tissue_map$protein_tissue,
          selected = tissue_map$protein_tissue[1]
        ),
        actionButton(ns("run"), "Run Outlier Detection"),
        br(), br(),
        h3("Outlier Genes Plot"),
        plotOutput(ns("outlier_plot")),
        h3("Outlier Gene List"),
        verbatimTextOutput(ns("outlier_vector"))
      )
    )
  )
}

# Server function
outlierServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # EventReactive for normal Run button
    result <- eventReactive(input$run, {
      req(input$tissue)
      detect_outliers(input_tissue = input$tissue)
    })

    # Observe Example button click
    observeEvent(input$run_example, {
      # Update tissue input to "lung"
      updateSelectInput(session, "tissue", selected = "lung")

      # Trigger the same reactive as Run button
      result()
    })

    # Outputs
    output$outlier_plot <- renderPlot({
      req(result())
      result()$outlier_plot
    })

    output$outlier_vector <- renderPrint({
      req(result())
      result()$outlier_vector
    })
  })
}
