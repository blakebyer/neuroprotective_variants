library(tidyverse)
library(tidyr)
library(readxl)
library(httr2)
library(jsonlite)
library(stringr)
library(purrr)

## Read Variant and Gene Data
var_data <- read_xlsx("neuroprotective_variants.xlsx", sheet = "compiled variants")
snps <- var_data %>% select(rsid) %>% distinct() %>% drop_na() %>% as.vector()
genes <- var_data %>% select(gene) %>% distinct() %>% drop_na() %>% as.vector()

## Ensembl

## variant effect prediction
ensembl_api <- function(snps) {
  
  base_url <- "https://rest.ensembl.org"
  ext <- "/vep/human/id"
  snps <- as.list(snps) %>% setNames("ids")
  
  query <- snps
  
  # Make an HTTP request
  post_vep <- request(base_url) |>
    req_url_path(path = ext) |>
    req_headers(accept = "application/json",
                `Content-Type` = "application/json") |>
    req_body_json(query) |>
    req_method("POST")

  resp_vep <- req_perform(post_vep) |>
    resp_body_json()
  
  return(resp_vep)
}

parse_ensembl <- function(ensembl_json) {
  
}

vep <- ensembl_api(snps)

## OMIM

api_key <- "J4eGvqcbQuiOdjflqQQwoQ"

omim_api <- function(genes, api_key) {
  # make request string
  search_query <- paste(unlist(genes), collapse = "+OR+")
  
  # set API key in browser cookies
  base_url <- "https://api.omim.org"
  api_url <-  sprintf("/api/apiKey?apiKey=%s&format=json", api_key)
  query_url <- sprintf(
    "/api/clinicalSynopsis/search?search=%s&format=json&start=0&limit=5000&apiKey=%s",
    search_query, api_key
  )

  # Set API key in browser cookies
  post_api <- request(base_url) |>
    req_url_path(path = api_url) |>
    req_headers(accept = "application/json") |>
    req_method("POST")

  # perform API post
  resp_api <- req_perform(post_api)
  
  # search for clinical terms and gene to phenotype map
  get_ct <- request(base_url) |>
    req_url_path(path = query_url) |>
    req_method("GET")

  resp_ct <- req_perform(get_ct) |>
    resp_body_json()

  return(resp_ct)
}

parse_omim <- function(omim_json) {
  entries <- omim_json$omim$searchResponse$clinicalSynopsisList
  
  omim_df <- map_dfr(entries, ~tibble(
    omim_entry = .x$clinicalSynopsis$mimNumber,
    disease = .x$clinicalSynopsis$preferredTitle,
    matches = .x$clinicalSynopsis$matches
  ))
  
  return(omim_df)
}

ct <- omim_api(genes, api_key)
omim_df <- parse_omim(ct)

# Convert to dataframe in one step
# omim_df <- map_dfr(entries, ~tibble(
#   omim_entry = .x$entry$mimNumber,
#   title = .x$entry$titles$preferredTitle,
#   gene = .x$entry$geneMap$approvedGeneSymbols,
#   gene_name = .x$entry$geneMap$geneName,
#   chromosome = .x$entry$geneMap$chromosomeSymbol,
#   chr_start = .x$entry$geneMap$chromosomeLocationStart,
#   chr_end = .x$entry$geneMap$chromosomeLocationEnd,
#   transcript = .x$entry$geneMap$transcript,
#   matches = .x$entry$matches,
#   phenotype = .x$entry$phenotypeMapList$phenotypeMap
# ))


## Gene Ontology

# "id": "GO%3A0003674",
# "label": "molecular_function"
# 
# "id": "GO%3A0008150",
# "label": "biological_process"

# "id": "GO%3A0005575",
# "label": "cellular_component"

panther_api <- function(genes, annotation) {
  if (str_detect(annotation, ":")) {
    annotation = str_replace(annotation, ":", "%")
  }
  if (annotation == "GO%3A0008150") {
    cat("Returning biological processes...")
  }
  if (annotation == "GO%3A0005575") {
    cat("Returning cellular components...")
  }
  if (annotation == "GO%3A0003674") {
    cat("Returning molecular functions...")
  }
  genes <- paste(unlist(genes), collapse = ",")
  
  # Define base and query url
  base_url <- "https://www.pantherdb.org"
  query_url <- sprintf("/services/oai/pantherdb/enrich/overrep?geneInputList=%s&organism=9606&annotDataSet=%s&enrichmentTestType=FISHER&correction=FDR", genes, annotation)
  
  # Make an HTTP request
  post_panther <- request(base_url) |>
    req_url_path(path = query_url) |>
    req_headers(accept = "application/json") |>
    req_method("POST")

  resp_panther <- req_perform(post_panther) |>
    resp_body_json()

  return(resp_panther)
}

panther_bp <- panther_api(genes, "GO%3A0008150")
panther_cc <- panther_api(genes, "GO%3A0005575")
panther_mf <- panther_api(genes, "GO%3A0003674")

parse_panther <- function(panther_json) {
  
  
}

## KEGG





## 