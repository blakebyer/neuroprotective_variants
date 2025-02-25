library(shiny)
library(shinyjs)
library(shinyWidgets)
library(bslib)
library(DT)
source("pipeline.R")

options(shiny.maxRequestSize = 10 * 1024^2) # 10 MB file

# Define UI for app that does gene set enrichment
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
      fileInput("upload", "Upload a gene list (CSV, TSV, or TXT)", accept = c(".csv", ".tsv", ".txt")),
      actionButton("submit", "Submit"),
      actionButton("reset", "Reset")
    ),
    navset_card_pill(
      nav_panel("Gene Table",
                progressBar(
                  id = "get_hgnc",
                  value = 0,
                  title = "",
                  display_pct = FALSE,
                  striped = TRUE
                ),
                span(textOutput("hgnc_warnings"), style="color:red"),
                DTOutput("hgnc")
                ),
      nav_panel("OMIM", 
                progressBar(
                  id = "get_omim",
                  value = 0,
                  title = "",
                  display_pct = FALSE,
                  striped = TRUE
                ),
                DTOutput("omim")
                ),
      nav_panel("GO", 
                progressBar(
                  id = "get_go",
                  value = 0,
                  title = "",
                  display_pct = FALSE,
                  striped = TRUE
                ),
                fluidRow(
                  column(4, actionButton("bp_button", "Biological Process", class="btn-primary")),
                  column(4, actionButton("mf_button", "Molecular Function", class="btn-primary")),
                  column(4, actionButton("cc_button", "Cellular Component", class="btn-primary"))
                ),
                h4(textOutput("go_title")),
                navset_card_tab(
                  nav_panel("Gene Set",
                        DTOutput("go")),
                  nav_panel("Hierarchy",
                            fluidRow(
                              column(4, numericInput("go_k", "K Clusters", value = 1, min = 1, max = 100)),
                              column(4, numericInput("go_n_plot", "Plot Cluster K", value = 0, min = 0, max = 100))
                            ),
                            plotOutput("go_hierarchy")),
                  nav_panel("Bar Plot",
                            fluidRow(
                              column(4, sliderInput("go_slice", "Top n:",
                                                    min = 10, max = 100,
                                                    value = 10, step = 10))
                            ),
                            plotOutput("go_plot"))
                            )),
      nav_panel("GSEA",
                progressBar(
                  id = "get_gsea",
                  value = 0,
                  title = "",
                  display_pct = FALSE,
                  striped = TRUE
                ),
                fluidRow(
                  column(3, actionButton("do_button", "Disease Ontology", class="btn-primary")),
                  column(3, actionButton("hpo_button", "Human Phenotype Ontology", class="btn-primary")),
                  column(3, actionButton("pw_button", "Pathway Ontology", class="btn-primary")),
                  column(3, actionButton("chebi_button", "CHEBI", class="btn-primary"))
                ),
                h4(textOutput("gsea_title")),
      navset_card_tab(
        nav_panel("Gene Set",
          DTOutput("gsea")),
        nav_panel("Hierarchy",
                  fluidRow(
                    column(4, numericInput("k", "K Clusters", value = 1, min = 1, max = 100)),
                    column(4, numericInput("n_plot", "Plot Cluster K", value = 0, min = 0, max = 100))
                  ),
                  plotOutput("hierarchy")),
        nav_panel("Bar Plot",
                  fluidRow(
                    column(4, sliderInput("slice", "Top n:",
                                          min = 10, max = 100,
                                          value = 10, step = 10))
                  ),
                  plotOutput("b_plot"))
        )
      )
)
)
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
           txt = vroom::vroom(input$upload$datapath, delim = "\n"),
           validate("Invalid file; Please upload a .csv, .tsv, or .txt file")
    )
  })
  
  normalize <- eventReactive(input$submit,
    {
      req(isTruthy(input$text) || isTruthy(input$upload))
      
      updateProgressBar(session, id = "get_hgnc", value = 25)
      genes <- NULL
      if (!is.null(input$text) && nzchar(input$text)) {
        genes <- strsplit(trimws(input$text), "[,\t\n]+")
      }
      
      # Handle file upload
      if (!is.null(input$upload)) {
        df <- data()
        genes <- unique(unlist(df))  # Flatten and remove duplicates
      }
      
      updateProgressBar(session, id = "get_hgnc", value = 65)
      
      warnings <- NULL
      gene_df <- withCallingHandlers(
        tibble(normalize_genes(genes, session)),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning") # Prevents the warning from printing to the console
        }
      )
      updateProgressBar(session, id = "get_hgnc", value = 100)
      list(data = gene_df, hgnc_warnings = warnings)
      
    }
  )
  
  output$hgnc_warnings <- renderPrint({
    if (!is.null(normalize()$hgnc_warnings)) {
      cat("Warnings:\n", paste(normalize()$hgnc_warnings, collapse = "\n"))
    }
  })
  
  omim_set <- eventReactive(input$submit,
      {
        req(normalize()$data)
        updateProgressBar(session, id = "get_omim", value = 25)
        om <- omim_api(normalize()$data$symbol, api_key)
        
        omim_df <- tibble(parse_omim(om, session))
        updateProgressBar(session, id = "get_omim", value = 65)

        updateProgressBar(session, id = "get_omim", value = 100)
        return(omim_df)
      }                  
  )
  
  selected_annot <- reactiveVal("BP")  # Default annotation is "BP"
  go_title <- reactiveVal("Biological Process")
  
  # Reset annotation to "BP" on submit
  observeEvent(input$submit, {
    selected_annot("BP")
  })
  
  # Update annotation type when buttons are clicked
  observeEvent(input$bp_button, { selected_annot("BP") 
    go_title("Biological Process")
    })
  observeEvent(input$mf_button, { selected_annot("MF") 
    updateProgressBar(session, id = "get_go", value = 0)
    go_title("Molecular Function")
    })
  observeEvent(input$cc_button, 
    { selected_annot("CC")
      updateProgressBar(session, id = "get_go", value = 0)
      go_title("Cellular Component")
    })
  
  # Call Panther API based on selected annotation, triggered by input$submit or selected_annot
  go_annot <- reactive({
    req(isTruthy(normalize()) || isTruthy(selected_annot()))  # Ensure normalized gene data is available
    genes <- normalize()$data$symbol
    annot <- selected_annot()  # Use the selected annotation type
    
    updateProgressBar(session, id = "get_go", value = 25)
    go <- panther_api(genes, annot)
    updateProgressBar(session, id = "get_go", value = 65)
    go_df <- tibble(parse_panther(go, session))
    updateProgressBar(session, id = "get_go", value = 100)
    
    return(go_df)
  })
  
  selected_onto <- reactiveVal("RDO") 
  gsea_title <- reactiveVal("Disease Ontology")
  
  # Reset annotation to "RDO" on submit
  observeEvent(input$submit, {
    selected_onto("RDO")
  })
  
  # Update annotation type when buttons are clicked
  observeEvent(input$do_button, { selected_onto("RDO") 
    gsea_title("Disease Ontology")
    })
  observeEvent(input$hpo_button, { selected_onto("MP") 
    updateProgressBar(session, id = "get_go", value = 0)
    gsea_title("Human Phenotype Ontology")
    })
  
  observeEvent(input$pw_button, { selected_onto("PW") 
    updateProgressBar(session, id = "get_go", value = 0)
    gsea_title("Pathway Ontology")
  })

  observeEvent(input$chebi_button, { selected_onto("CHEBI") 
    updateProgressBar(session, id = "get_go", value = 0)
    gsea_title("Chemical Entities of Biological Interest")
  })
  
  gene_set <- reactive(
    {
      req(isTruthy(normalize()) || isTruthy(selected_onto()))
      genes <- normalize()$data$symbol
      updateProgressBar(session, id = "get_gsea", value = 25)
      onto <- selected_onto()
      gs <- gsea(genes, onto)
      updateProgressBar(session, id = "get_gsea", value = 65)
      set <- parse_gsea(gs, session)
      updateProgressBar(session, id = "get_gsea", value = 100)
      return(set)
    }
  )
  
  observeEvent(input$k, {
    updateNumericInput(session, "n_plot", max = input$k)
  })
  
  observeEvent(input$go_k, {
    updateNumericInput(session, "go_n_plot", max = input$go_k)
  })
  
  observeEvent(gene_set(),
      {
      if (!is.null(gene_set()) && length(gene_set()) > 0) {
        max_clusters <- max(1, length(gene_set()) - 1)
        updateNumericInput(session, "k", max = max_clusters)
      }         
  })
  
  observeEvent(go_annot(),
     {
       if (!is.null(go_annot()) && length(go_annot()) > 0) {
         max_go_clusters <- max(1, length(go_annot()) - 1)
         updateNumericInput(session, "go_k", max = max_go_clusters)
       }         
     })
  
  output$hgnc <- DT::renderDT({
    req(normalize())
    normalize()$data
  })
  
  output$omim <- DT::renderDT({
    req(omim_set())
    omim_set()
  })
  
  output$go_title <- renderText({
    req(go_annot())
    go_title()
  })
  
  output$go <- DT::renderDT({
    req(go_annot())
    go_annot()
  })
  
  output$gsea_title <- renderText({
    req(gene_set())
    gsea_title()
  })
  
  output$gsea <- DT::renderDT({
    req(gene_set())
    gene_set()
  })
  
  output$hierarchy <- renderPlot({
    req(gene_set())
    hierarchy(gene_set(), input$k, input$n_plot) 
  })
  
  output$go_hierarchy <- renderPlot({
    req(go_annot())
    hierarchy(go_annot(), input$go_k, input$go_n_plot) 
  })
  
  output$b_plot <- renderPlot({
    req(gene_set())
    panther_plot(gene_set(), input$slice) 
  })
  
  output$go_plot <- renderPlot({
    req(go_annot())
    panther_plot(go_annot(), input$go_slice) 
  })
  
}

shinyApp(ui = ui, server = server)