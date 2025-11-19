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
            background-color: #00008B;
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
          br(),
          p("ExprCompareR is an R package designed to streamline the exploration and quantification of the relationship between RNA expression and protein expression across human tissues and conditions. Researchers studying post-transcriptional regulation currently face the time-consuming task of manually retrieving RNA-seq and protein-expression data, filtering, normalizing, computing correlations, and creating visualizations for their genes of interest. `ExprCompareR` automates and integrates these steps, providing a user-friendly workflow to quickly retrieve transcriptomic and proteomic data from public repositories such as the Genotype-Tissue Expression Project (GTEx)(Lonsdale et al., 2013) and the Human Protein Atlas (HPA)(Thul et al., 2018), compute correlations, flag outlier genes, and visualize results.")
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
          br(),
          p("On each Shiny tab you will find the following sections:"),
          br(),
          HTML('
            <table class="custom-table">
              <thead>
                <tr>
                  <th>Section</th>
                  <th>Explanation</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>Description</td>
                  <td>Explains the main functionality of the function</td>
                </tr>
                <tr>
                  <td>Input</td>
                  <td>Explains the input for that function</td>
                </tr>
                <tr>
                  <td>Output</td>
                  <td>Explains the output for that function</td>
                </tr>
                <tr>
                  <td>Instructions</td>
                  <td>Explains the UI and how to run the function using shiny</td>
                </tr>
                 <tr>
                  <td>Run Examples </td>
                  <td>Explains how the user can run examples for the function</td>
                </tr>
                 <tr>
                  <td>Interpretation of plot</td>
                  <td>Explains the plot generated and how to interpret the plot</td>
                </tr>
                 <tr>
                  <td>Interpretation of list</td>
                  <td>Explains the list or table generated and how to interpret list or table</td>
                </tr>
                 <tr>
                  <td>References</td>
                  <td>Includes the references for that function</td>
                </tr>
              </tbody>
            </table>
          '),
          br(),
          h4("Simply click on a tab to access its interface, input your selections, and view the resulting plots and tables interactively.")


      ),


      div(class = "custom-card",
          h3("References"),
          p("* Cetinkaya-Rundel M,Cheng J, Grolemund G (2017).Customize your UI with HTML.https://shiny.posit.co/r/articles/build/html-tags/"),
          p("* Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Aden-Buie G, Xie Y,Allen J, McPherson J, Dipert A, Borges B (2025).shiny: Web Application Framework for R. R package version 1.11.1,https://CRAN.R-project.org/package=shiny"),
          p("* OpenAI.(2025). ChatGPT (GPT-5) Large language model.Retrieved November 17, 2025, from https://chatgpt.com/ "),
          p("* Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019).,A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data.BMC Bioinformatics. 2019;20(1):394. https://pubmed.ncbi.nlm.nih.gov/31311497/"),
          p("* Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).HPAanalyze: an R package that facilitates the retrieval and analysis of the Human Protein Atlas data. MC Bioinformatics 20, 463 (2019).https://doi.org/10.1186/s12859-019-3059-z"),
          p("* Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). gtexr: A convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API.Journal of Open Source So ware, 10(109), 8249. ISSN 2475-9066,doi:10.21105/joss.08249}, gigs v0.2.1."),
          p("* Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3 319-24277-4, https://ggplot2.tidyverse.org."),
          p("* Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).dplyr: A Grammar of Data Manipulation. R package version 1.1.4, https://dplyr.tidyverse.org."),
          p("* Wickham H, Henry L (2025). purrr: Functional Programming Tools. R package version 1.1.0, https://purrr.tidyverse.org/."),
          p("* Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data. R package version 1.3.1, https://tidyr.tidyverse.org.")
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
