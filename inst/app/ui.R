library(shiny)

ui <- fluidPage(

  # App title ----
  titlePanel("peviz"),

  sidebarLayout(

    sidebarPanel(

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

    mainPanel(width = 10,
              plotOutput("msa"),
              plotOutput("tree")
    )
  )
)
