library(shiny)

ui <- fluidPage(
  titlePanel("Protein Evolution Visualization"),

  fluidPage(
    # --- DATA SELECTION BLOCK ---
    fluidRow(
      column(8,
             selectInput('Database', 'Database:',
                         choices = c('Swiss-Prot', 'Local UniProt Fasta'),
                         selected = 'Swiss-Prot'),

             conditionalPanel(
               condition = "input.Database == 'Local UniProt Fasta'",
               fileInput('local_fasta', 'Local UniProt Fasta')
             ),

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
      column(4,
             selectizeInput('uniprot_category', 'Filter by Category:',
                            choices = NULL, multiple = TRUE,
                            options = list(placeholder = 'Select categories (e.g., DOMAINS_AND_SITES)'))),
      column(4,
             selectizeInput('uniprot_type', 'Filter by Type:',
                            choices = NULL, multiple = TRUE,
                            options = list(placeholder = 'Select types'))),
      column(4,
             br(), # Adds spacing so the checkbox aligns vertically with the dropdowns
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
      column(8,
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
      width = 8,
      div(style = 'font-size:75%', DT::DTOutput('fasta_table'))
    )
  ),

  br(),

  fluidRow(
    column(4, downloadButton('fasta_download', 'Download Fasta',
                             style = 'background-color: #3269FF; color: #ffffff')),
    column(4, downloadButton('fasta_aligned_download', 'Download Aligned Fasta',
                             style = 'background-color: #3269FF; color: #ffffff'))
  ),

  br(), br()
)
