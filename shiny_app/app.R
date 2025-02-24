library(shiny)
library(shinyWidgets)
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
    navset_card_pill(
      # title = "Visualizations",
      nav_panel("Gene Table", 
                progressBar(
                  id = "get_hgnc",
                  value = 0,
                  title = "",
                  display_pct = FALSE,
                  striped = TRUE
                ),
                dataTableOutput("hgnc")),
      nav_panel("OMIM", 
                progressBar(
                  id = "get_omim",
                  value = 0,
                  title = "",
                  display_pct = FALSE,
                  striped = TRUE
                ),
                dataTableOutput("omim")),
      nav_panel("GO", 
                progressBar(
                  id = "get_go",
                  value = 0,
                  title = "",
                  display_pct = FALSE,
                  striped = TRUE
                ),
                fluidRow(
                  column(4, actionButton("bp_button", "Biological Process")),
                  column(4, actionButton("mf_button", "Molecular Function")),
                  column(4, actionButton("cc_button", "Cellular Component"))
                ),
                dataTableOutput("go")),
      nav_panel("GSEA", 
                progressBar(
                  id = "get_gsea",
                  value = 0,
                  title = "",
                  display_pct = FALSE,
                  striped = TRUE
                ),
                fluidRow(
                  column(4, actionButton("do_button", "Disease Ontology")),
                  column(4, actionButton("hpo_button", "Human Phenotype Ontology"))
                ),
                dataTableOutput("gsea"))
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  normalize <- eventReactive(input$submit,
    {
      req(input$text)
      genes <- strsplit(input$text, "[,\t\n]+")
      updateProgressBar(session, id = "get_hgnc", value = 25)
      gene_df <- tibble(normalize_genes(genes))
      updateProgressBar(session, id = "get_hgnc", value = 100)
      return(gene_df)
    }
  )
  
  omim_set <- eventReactive(input$submit,
      {
        req(normalize())
        updateProgressBar(session, id = "get_omim", value = 25)
        om <- omim_api(normalize()$symbol, api_key)
        dis <- parse_omim(om)
        updateProgressBar(session, id = "get_omim", value = 100)
        return(dis)
      }                  
    )
  
  selected_annot <- reactiveVal("BP") 
  
  # Update annotation type when buttons are clicked
  observeEvent(input$bp_button, { selected_annot("BP") })
  observeEvent(input$mf_button, { selected_annot("MF") })
  observeEvent(input$cc_button, { selected_annot("CC") })
  
  # Call Panther API based on selected annotation
  go_annot <- reactive({
    req(normalize(), selected_annot())
    genes <- normalize()$symbol
    annot <- selected_annot()  # Use the selected annotation type
    updateProgressBar(session, id = "get_go", value = 25)
    go <- panther_api(genes, annot)
    updateProgressBar(session, id = "get_go", value = 65)
    go_df <- tibble(parse_panther(go))
    updateProgressBar(session, id = "get_go", value = 100)
    return(go_df)
  })
  
  output$hgnc <- renderDataTable({
    normalize()
  })
  
  selected_onto <- reactiveVal("RDO") 
  
  # Update annotation type when buttons are clicked
  observeEvent(input$do_button, { selected_onto("RDO") })
  observeEvent(input$hpo_button, { selected_onto("MP") })

  gene_set <- reactive(
    {
      req(normalize(), selected_onto())
      genes <- normalize()$symbol
      updateProgressBar(session, id = "get_gsea", value = 25)
      onto <- selected_onto()
      gs <- gsea(genes, onto)
      updateProgressBar(session, id = "get_gsea", value = 65)
      set <- parse_gsea(gs)
      updateProgressBar(session, id = "get_gsea", value = 100)
      return(set)
    }
  )
  
  # output$hierarchy <- renderPlot({
  #   hierarchy(normalize$symbol)  # Assuming hierarchy is a plotting function
  # })
  
  output$omim <- renderDataTable({
    req(omim_set())
    omim_set()
  })
  
  output$go <- renderDataTable({
    req(go_annot())
    go_annot()
  })
  
  output$gsea <- renderDataTable({
    req(gene_set())
    gene_set()
  })
  
}

shinyApp(ui = ui, server = server)