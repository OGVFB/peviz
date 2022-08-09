library(shiny)

ui <- fluidPage(

  # App title ----
  titlePanel("peviz"),

  fluidPage(

    fluidRow(width = 8,
      selectInput(inputId = 'Database',
                  label = 'Database:',
                  selected = 'Swiss-Prot',
                  choices = c('Swiss-Prot', 'Local Fasta')),
      conditionalPanel(condition = "input.Database == 'Local Fasta'",
                       fileInput('local_fasta', 'Local Fasta', multiple = FALSE)),
      selectizeInput(inputId = "uniprot",
                     label = "Select proteins:",
                     choices=NULL,
                     multiple=TRUE,
                     width = '100%'
      ),
      selectizeInput(inputId = "primary_protein",
                     label = "Add UniProt domain(s) based on one protein:",
                     choices=NULL,
                     multiple=FALSE,
                     options = list(placeholder = 'Optional'),
                     width = '100%'
      ),
      actionButton('Draw','Draw', alt =  'Draw Plots',
                   style='background-color: #3269FF; color: #ffffff')

    ),
    br(), br(),

    fluidRow(width = 11,
              plotOutput("msa", height = "auto"),
              plotOutput("tree")
    )
  )
)
