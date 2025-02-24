library(shiny)
library(bslib)
library(DT)
source("pipeline.R")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Shiny GSEA"),
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        "text",
        label = "Enter a list of comma, tab, or newline-separated genes or SNPs:",
        height = '300px',
        width = '200px',
        resize = 'vertical'
      ),
      fileInput("upload", "Upload a gene list (CSV or TSV)"),
      actionButton("submit", "Submit")
    ),
    navset_card_tab(
      title = "Visualizations",
      nav_panel("Gene Table", dataTableOutput("hgnc")),
      #nav_panel("GSEA", dataTableOutput("gsea")),
      nav_panel("OMIM", dataTableOutput("omim"))
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  normalize <- eventReactive(input$submit,
    {
      req(input$text)
      genes <- as.vector(strsplit(input$text, "[,\t\n]+"))
      gene_df <- tibble(normalize_genes(genes))
      return(gene_df)
    }
  )
  
  # gene_set <- eventReactive(input$submit,
  #   {
  #     req(normalize())
  #     type = "CC"
  #     gs <- gsea(normalize()$symbol, type)
  #     set <- parse_gsea(gs)
  #     return(set)
  #   }
  # )
  
  omim_set <- eventReactive(input$submit,
      {
        req(normalize())
        om <- omim_api(normalize()$symbol, api_key)
        dis <- parse_omim(om)
        return(dis)
      }                  
    )
  
  output$hgnc <- renderDataTable({
    normalize()
  })
  
  # output$hierarchy <- renderPlot({
  #   hierarchy(normalize$symbol)  # Assuming hierarchy is a plotting function
  # })
  
  # output$gsea <- renderDataTable({
  #   req(normalize())
  #   gene_set()
  # })
  
  output$omim <- renderDataTable({
    req(normalize())
    omim_set()
  })
  
}

shinyApp(ui = ui, server = server)