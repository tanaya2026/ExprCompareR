# Prepare top 1000 genes including example genes
example_genes <- c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
remaining_genes <- setdiff(gene_symbols_list, example_genes)
set.seed(123)
random_genes <- sample(remaining_genes, 995)
top_genes <- c(example_genes, random_genes)


correlationGenesTissuesUI <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarLayout(
      sidebarPanel(
        sliderInput(
          inputId = ns("num_genes"),
          label = "Number of genes (5-10):",
          min = 5,
          max = 10,
          value = 5,
          step = 1
        ),
        uiOutput(ns("gene_selects")),

        sliderInput(
          inputId = ns("num_tissues"),
          label = "Number of tissues (5-10):",
          min = 5,
          max = 10,
          value = 5,
          step = 1
        ),
        uiOutput(ns("tissue_selects")),

        selectInput(
          inputId = ns("plot_choice"),
          label = "Select plot to display:",
          choices = c("per_gene", "per_tissue"),
          selected = "per_gene"
        ),

        actionButton(ns("run"), "Run Correlation")
      ),
      mainPanel(
        h3("Correlation Plot"),
        plotOutput(ns("cor_plot"))
      )
    )
  )
}

correlationGenesTissuesServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # Dynamic gene dropdowns
    output$gene_selects <- renderUI({
      n <- input$num_genes
      selects <- lapply(1:n, function(i) {
        selectizeInput(
          ns(paste0("gene_", i)),
          label = paste("Gene", i),
          choices = top_genes,
          selected = top_genes[i],
          multiple = FALSE,
          options = list(maxOptions = 1000)
        )
      })
      do.call(tagList, selects)
    })

    # Dynamic tissue dropdowns
    output$tissue_selects <- renderUI({
      n <- input$num_tissues
      selects <- lapply(1:n, function(i) {
        selectInput(
          ns(paste0("tissue_", i)),
          label = paste("Tissue", i),
          choices = tissue_map$protein_tissue,
          selected = tissue_map$protein_tissue[i]
        )
      })
      do.call(tagList, selects)
    })

    # Gather selected genes
    selected_genes <- reactive({
      n <- input$num_genes
      sapply(1:n, function(i) input[[paste0("gene_", i)]])
    })

    # Gather selected tissues
    selected_tissues <- reactive({
      n <- input$num_tissues
      sapply(1:n, function(i) input[[paste0("tissue_", i)]])
    })

    # Run correlation when button clicked
    result <- eventReactive(input$run, {
      req(selected_genes(), selected_tissues(), input$plot_choice)
      correlation_genes_tissues(
        gene_NAMES = selected_genes(),
        tissue_NAMES = selected_tissues(),
        plot_choice = input$plot_choice
      )
    })

    # Render selected plot
    output$cor_plot <- renderPlot({
      req(result())
      if (input$plot_choice == "per_gene") {
        result()$per_gene_plot
      } else {
        result()$per_tissue_plot
      }
    })

  })
}
