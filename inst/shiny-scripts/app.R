library(shiny)
library(shinyalert)

# Source modules
source(system.file("shiny-scripts/modules/module-outlier.R", package = "ExprCompareR"), local = TRUE)
source(system.file("shiny-scripts/modules/module-cancer_protein.R", package = "ExprCompareR"), local = TRUE)
source(system.file("shiny-scripts/modules/module-compute_correlation.R", package = "ExprCompareR"), local = TRUE)
source(system.file("shiny-scripts/modules/module-correlation_genes_tissues.R", package = "ExprCompareR"), local = TRUE)

# UI
ui <- navbarPage(
  "ExprCompareR",

  # --- Introduction Tab ---
  tabPanel(
    "Introduction",
    fluidPage(

      # Custom CSS for tables and cards
      tags$head(
        tags$style(HTML("
          body {font-family: Arial, sans-serif;}
          .custom-card {
            background-color: #f9f9f9;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 25px;
            box-shadow: 1px 1px 5px rgba(0,0,0,0.1);
          }
          table.custom-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
          }
          table.custom-table th {
            background-color: #4CAF50;
            color: white;
            padding: 10px;
            text-align: left;
          }
          table.custom-table td {
            border: 1px solid #ddd;
            padding: 8px;
            vertical-align: top;
          }
          table.custom-table tr:nth-child(even){background-color: #f2f2f2;}
          table.custom-table tr:hover {background-color: #e0e0e0;}
          code {background-color: #f1f1f1; padding: 2px 5px; border-radius: 4px;}
        "))
      ),

      # --- Description Card ---
      div(class="custom-card",
          h2("Welcome to ExprCompareR"),
          p("`ExprCompareR` is an R package designed to streamline the exploration and quantification of the relationship between RNA expression and protein expression across human tissues and conditions. Researchers studying post-transcriptional regulation currently face the time-consuming task of manually retrieving RNA-seq and protein-expression data, filtering, normalizing, computing correlations, and creating visualizations for their genes of interest. `ExprCompareR` automates and integrates these steps, providing a user-friendly workflow to quickly retrieve transcriptomic and proteomic data from public repositories such as the Genotype-Tissue Expression Project (GTEx)(Lonsdale et al., 2013) and the Human Protein Atlas (HPA)(Thul et al., 2018), compute correlations, flag outlier genes, and visualize results.")
      ),

      # --- Key Functions Card ---
      div(class="custom-card",
          h3("Summary Table of Key Functions"),
          HTML('
            <table class="custom-table">
              <thead>
                <tr>
                  <th>Function</th>
                  <th>Description</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td><code>compute_correlation()</code></td>
                  <td>Calculates the Spearman correlation between RNA and protein expression for a user-defined list of genes <strong>or</strong> tissues and generates correlation plots. Acts as a wrapper and calls either <code>correlation_genes_only()</code> or <code>correlation_tissues_only()</code> internally.</td>
                </tr>
                <tr>
                  <td><code>correlation_genes_tissues()</code></td>
                  <td>Calculates the Spearman correlation between RNA and protein expression for a user-defined list of genes <strong>and</strong> tissues and generates correlation plots.</td>
                </tr>
                <tr>
                  <td><code>detect_outliers()</code></td>
                  <td>Identifies genes that exhibit significant discrepancies between RNA and protein expression within a selected tissue, and generates outlier plots.</td>
                </tr>
                <tr>
                  <td><code>compareCancerProtein()</code></td>
                  <td>Compares protein expression levels between cancerous and normal tissue samples for a specified cancer type, and plots genes involved.</td>
                </tr>
              </tbody>
            </table>
          ')
      ),

      # --- Navigation Card ---
      div(class="custom-card",
          h3("How to Navigate the Shiny Application"),
          p("Upon launching the ExprCompareR Shiny app, you will see a navigation bar at the top containing four tabs, each corresponding to a core functionality of the app. These tabs allow you to seamlessly switch between different analyses."),
          p("The table below shows which Shiny tab corresponds to each key function in ExprCompareR:"),
          HTML('
            <table class="custom-table">
              <thead>
                <tr>
                  <th>Tab Name</th>
                  <th>Corresponding Function</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>Correlation (Genes OR Tissues)</td>
                  <td><code>compute_correlation()</code></td>
                </tr>
                <tr>
                  <td>Correlation (Genes AND Tissues)</td>
                  <td><code>correlation_genes_tissues()</code></td>
                </tr>
                <tr>
                  <td>Outlier Detection</td>
                  <td><code>detect_outliers()</code></td>
                </tr>
                <tr>
                  <td>Differential Expression in Cancer Protein</td>
                  <td><code>compareCancerProtein()</code></td>
                </tr>
              </tbody>
            </table>
          '),
          p("Simply click on a tab to access its interface, input your data or selections, and view the resulting plots and tables interactively. This layout ensures that all major functions are easily accessible without leaving the app.")
      ),

      # --- Example Datasets Card ---
      div(class="custom-card",
          h3("Example Datasets"),
          p("You can include here a description or links to example datasets for users to explore the functionalities of ExprCompareR.")
      )
    )
  ),

  # --- Functional Tabs ---
  tabPanel("Correlation (Genes OR Tissues)", correlationUI("cor_tab")),
  tabPanel("Correlation (Genes AND Tissues)", correlationGenesTissuesUI("cor_gt_tab")),
  tabPanel("Outlier Detection", outlierUI("outlier_tab")),
  tabPanel("Differential Expression in Cancer Protein", cancerUI("cancer_tab"))
)

# Server
server <- function(input, output, session) {
  correlationServer("cor_tab")
  correlationGenesTissuesServer("cor_gt_tab")
  outlierServer("outlier_tab")
  cancerServer("cancer_tab")
}

# Launch Shiny App
shiny::shinyApp(ui, server)
