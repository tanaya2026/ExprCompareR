# Cancer Protein Module

cancerUI <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId = ns("cancer_type"),
          label = "Select a cancer type:",
          choices = cancer_tissue_map$cancer,
          selected = cancer_tissue_map$cancer[1]
        ),
        actionButton(ns("run"), "Run Cancer Protein Analysis")
      ),
      mainPanel(
        h3("Protein Rank Shift Plot"),
        plotOutput(ns("plot")),

        h3("Protein Expression Table"),
        verbatimTextOutput(ns("table"))
      )
    )
  )
}

cancerServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    result <- eventReactive(input$run, {
      req(input$cancer_type)
      compareCancerProtein(cancer_type = input$cancer_type)
    })

    output$plot <- renderPlot({
      req(result())
      result()$plot
    })

    output$table <- renderPrint({
      req(result())
      result()$table  # Display the full tibble
    })
  })
}
