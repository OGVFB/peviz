#' The application User-Interface
#'
#' @import shiny
#' @importFrom DT DTOutput
#' @noRd
app_ui <- function(request) {

  # Check for local data inside the installed package's inst/extdata directory
  data_path <- system.file("extdata", "peviz_uniprot_data.Rdata", package = "peviz")
  has_swissprot <- file.exists(data_path) && data_path != ""

  # Build the choices vector dynamically
  db_choices <- c('Local UniProt Fasta', 'UniProt API Search')
  if (has_swissprot) {
    db_choices <- c('Swiss-Prot', db_choices)
  }

  fluidPage(
    titlePanel("Protein Evolution Visualization"),

    fluidPage(
      # --- DATA SELECTION BLOCK ---
      fluidRow(
        column(9,
               selectInput('Database', 'Database:',
                           choices = db_choices,
                           selected = 'UniProt API Search'),

               conditionalPanel(
                 condition = "input.Database == 'Local UniProt Fasta'",
                 fileInput('local_fasta', 'Local UniProt Fasta')
               ),

               # --- UniProt API Search Panel ---
               conditionalPanel(
                 condition = "input.Database == 'UniProt API Search'",
                 wellPanel(
                   textInput("api_query", "UniProt Search Query (e.g., gene:CRX AND taxonomy_id:9606):",
                             value = "gene:CRX AND taxonomy_id:9606"),
                   actionButton("search_api", "Search UniProt", style = "background-color: #3269FF; color: white;"),
                   br(), br(),
                   div(style = 'font-size:85%', DT::DTOutput("api_results_table")),
                   br(),
                   actionButton("fetch_api_seqs", "Fetch Selected Sequences", style = "background-color: #28a745; color: white;")
                 )
               ),
               # -------------------------------------

               selectizeInput('uniprot', 'Select proteins:',
                              choices = NULL, multiple = TRUE, width = '100%'),

               selectizeInput('primary_protein', 'Add UniProt domain(s) based on one protein:',
                              choices = NULL, multiple = FALSE, width = '100%',
                              options = list(placeholder = 'Optional')),

               selectizeInput('reference_seq', 'Reference sequence for numbering:',
                              choices = NULL, multiple = FALSE, width = '100%',
                              options = list(placeholder = 'Select protein for axis numbering'))
        )
      ),

      # --- FILTER BLOCK ---
      fluidRow(
        column(3,
               selectizeInput('uniprot_category', 'Filter by Category:',
                              choices = NULL, multiple = TRUE,
                              options = list(placeholder = 'Select categories (e.g., DOMAINS_AND_SITES)'))),
        column(3,
               selectizeInput('uniprot_type', 'Filter by Type:',
                              choices = NULL, multiple = TRUE,
                              options = list(placeholder = 'Select types'))),
        column(3,
               br(),
               checkboxInput('fetch_variants', 'Fetch Variant Data (Slower API Call)', value = FALSE))
      ),

      # --- PLOT SETTINGS & ACTION BUTTON ---
      fluidRow(
        column(3,
               selectizeInput('color', 'Color Scheme',
                              choices = NULL, multiple = FALSE,
                              options = list(placeholder = 'Select Color Scheme'))),
        column(3,
               selectizeInput('method', 'MSA Method',
                              choices = c('ClustalW', 'ClustalOmega', 'Muscle'),
                              multiple = FALSE,
                              options = list(placeholder = 'Select MSA Method')))
      ),

      fluidRow(
        column(9,
               actionButton('Draw', 'Draw Plots',
                            style = 'background-color: #3269FF; color: #ffffff; margin-top: 15px;'))
      ),

      br(), br(),

      # --- VISUALIZATIONS ---
      fluidRow(
        width = 11,
        plotOutput("msa", height = "auto"),
        plotOutput("tree")
      ),

      br(), br(),

      # --- DATA TABLES & DOWNLOADS ---
      fluidRow(
        width = 9,
        div(style = 'font-size:75%', DT::DTOutput('fasta_table'))
      )
    ),

    br(),

    fluidRow(
      column(3, downloadButton('fasta_download', 'Download Fasta',
                               style = 'background-color: #3269FF; color: #ffffff')),
      column(3, downloadButton('fasta_aligned_download', 'Download Aligned Fasta',
                               style = 'background-color: #3269FF; color: #ffffff'))
    ),

    br(), br()
  )
}
