# Prepare top 1000 genes including example genes
example_genes <- c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
remaining_genes <- setdiff(gene_symbols_list, example_genes)
set.seed(123)
random_genes <- sample(remaining_genes, 995)
top_genes <- c(example_genes, random_genes)


correlationUI <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          inputId = ns("input_type"),
          label = "Select input type:",
          choices = c("Genes", "Tissues"),
          inline = TRUE
        ),
        sliderInput(
          inputId = ns("num_inputs"),
          label = "Number of selections (5-10):",
          min = 5,
          max = 10,
          value = 5,
          step = 1
        ),
        uiOutput(ns("dynamic_selects")),
        actionButton(ns("run"), "Run Correlation")
      ),
      mainPanel(
        h3("Correlation Plot"),
        plotOutput(ns("cor_plot"))
      )
    )
  )
}

correlationServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Dynamically generate dropdowns based on slider input
    output$dynamic_selects <- renderUI({
      ns <- session$ns
      n <- input$num_inputs

      selects <- lapply(1:n, function(i) {
        if (input$input_type == "Genes") {
          # Use limited gene list
          selectizeInput(
            ns(paste0("sel_", i)),
            label = paste("Gene", i),
            choices = top_genes,
            selected = top_genes[i],
            multiple = FALSE,
            options = list(maxOptions = 1000)
          )
        } else {
          # Small list, regular selectInput is fine
          selectInput(
            ns(paste0("sel_", i)),
            label = paste("Tissue", i),
            choices = tissue_map$protein_tissue,
            selected = tissue_map$protein_tissue[i]
          )
        }
      })

      do.call(tagList, selects)
    })

    # Reactive expression to gather selected inputs
    selected_inputs <- reactive({
      n <- input$num_inputs
      vals <- sapply(1:n, function(i) {
        input[[paste0("sel_", i)]]
      })
      vals
    })

    # Run correlation when button is clicked
    result <- eventReactive(input$run, {
      req(selected_inputs())

      if (input$input_type == "Genes") {
        compute_correlation(gene_NAMES = selected_inputs(), tissue_NAMES = NULL)
      } else {
        compute_correlation(gene_NAMES = NULL, tissue_NAMES = selected_inputs())
      }
    })

    output$cor_plot <- renderPlot({
      req(result())
      if (input$input_type == "Genes") {
        result()$per_gene_plot
      } else {
        result()$per_tissue_plot
      }
    })

  })
}
