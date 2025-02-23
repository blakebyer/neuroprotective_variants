library(tidyverse)
library(tidyr)
library(readxl)
library(httr2)
library(jsonlite)
library(stringr)
library(purrr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(DOSE)
library(enrichplot)
library(ggupset)

## Read Variant and Gene Data
var_data <- read_xlsx("neuroprotective_variants.xlsx", sheet = "compiled variants")
snps <- var_data %>% dplyr::select(rsid) %>% distinct() %>% drop_na() %>% as.vector()
genes <- var_data %>% dplyr::select(gene) %>% distinct() %>% drop_na() %>% as.vector()

### Normalize Gene Identifiers
normalize_genes_hgnc <- function(genes) {
  genes_vec <- unlist(genes)
  len <- length(genes_vec)
  genes_string <- paste(genes_vec, collapse = "+OR+")
  
  base_url <- "https://clinicaltables.nlm.nih.gov"
  hgnc_query <- sprintf("/api/genes/v4/search?terms=%s&df=hgnc_id,symbol,name,alias_symbol,location&maxList=%s", genes_string, len*2)
  
  get_hgnc <- request(base_url) |>
      req_url_path(path = hgnc_query) |>
      req_method("GET")

  resp_hgnc <- req_perform(get_hgnc) |>
    resp_body_json()
  
  gene_df <- map_dfr(resp_hgnc[[4]], ~tibble(
          hgnc_id = as.character(.x[1]),
          symbol = as.character(.x[2]),
          name = as.character(.x[3]),
          alias_symbol = as.character(.x[4]),
          location = as.character(.x[5])
        )) %>% mutate(search_term = NA_character_)

  for (i in seq_along(gene_df$symbol)) {
    # iterate over columns
    for (column in c("symbol", "alias_symbol", "name")) {
      for (term in genes_vec) {
        # extract matching term
        match <- str_extract(gene_df[[column]][i], fixed(term, ignore_case = TRUE))
        if (!is.na(match)) {
          gene_df$search_term[i] <- match
          break  # exit loop once match is found
        }
      }
      if (!is.na(gene_df$search_term[i])) break # stop further checks
    }
  }

  gene_df <- gene_df %>%
    filter(!is.na(search_term))

  # mapped genes
  mapped_genes <- unique(gene_df$search_term)
  n_mapped <- length(mapped_genes)
  missing_genes <- genes_vec %>%
    setdiff(mapped_genes)
  n_missing <- length(missing_genes)

  # missing genes
  if (n_missing > 0) {
    missing_genes_str <- paste(missing_genes, collapse = ", ")
    warning(sprintf("%.2f%% of gene symbols are unmapped.\nUnmapped genes are: %s\nCheck spelling; HUGO nomenclature preferred.",
                    (n_missing/len) * 100, missing_genes_str))
  }
  else {
    cat(sprintf("Mapped %.2f%% of gene symbols.", (n_mapped/len) * 100))
  }

  return(gene_df)
}

hgnc <- normalize_genes_hgnc(genes)

## Ensembl


## variant effect prediction
ensembl_api <- function(snps) {
  
  base_url <- "https://rest.ensembl.org"
  ext <- "/vep/human/id"
  snps <- as.list(snps) %>% setNames("ids")
  
  # Make an HTTP request
  post_vep <- request(base_url) |>
    req_url_path(path = ext) |>
    req_headers(accept = "application/json",
                `Content-Type` = "application/json") |>
    req_body_json(snps) |>
    req_method("POST")
  
  # perform request
  resp_vep <- post_vep |>
    req_error(is_error = ~ FALSE) |> # error handling
    req_perform()

  status <- resp_status(resp_vep)
  
  resp <- resp_vep |> resp_body_json()

  if (status >= 400) {
    warning(sprintf("Ensembl server error: HTTP %s", status))
  }
  else {
    return(resp)
  }
}

vep <- ensembl_api(snps)

## OMIM

api_key <- "J4eGvqcbQuiOdjflqQQwoQ"

omim_api <- function(genes, api_key) {
  
  # set API key in browser cookies
  base_url <- "https://api.omim.org"
  api_url <-  sprintf("/api/apiKey?apiKey=%s&format=json", api_key)
  
  # # Set API key in browser cookies
  post_api <- request(base_url) |>
    req_url_path(path = api_url) |>
    req_headers(accept = "application/json") |>
    req_method("POST")
  
  # perform API post
  resp_api <- req_perform(post_api)
  
  fetch_omim <- function(genes_subset) {
    # make request string
    search_query <- paste(genes_subset, collapse = "+OR+")
    
    query_url <- sprintf(
      "/api/clinicalSynopsis/search?search=%s&format=json&start=0&limit=200&apiKey=%s",
      search_query, api_key
    )

    # search for clinical terms
    get_ct <- request(base_url) |>
      req_url_path(path = query_url) |>
      req_method("GET")
    
    resp_ct <- req_perform(get_ct) |>
      resp_body_json()
    
    return(resp_ct)
    
  }
  
  # subset genes into requests of 20 to not exceed search complexity
  subsets <- split(genes, ceiling(seq_along(genes) / 20))
  
  responses <- map(subsets, fetch_omim)

  return(responses)
}

parse_omim <- function(omim_json) {
  # parse json
  omim_df <- map_dfr(names(omim_json), function(i) {
    entries <- omim_json[[i]]$omim$searchResponse$clinicalSynopsisList
    
    map_dfr(entries, function(entry) {
      tibble(
        omim_entry = entry$clinicalSynopsis$mimNumber,
        disease = entry$clinicalSynopsis$preferredTitle,
        matches = paste(entry$clinicalSynopsis$matches, collapse = ", ")
      )
    })
  })
  
  return(omim_df)
}

genes <- hgnc$symbol

ct <- omim_api(genes, api_key)
omim_df <- parse_omim(ct)


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
  else if (annotation == "GO%3A0008150") {
    cat("Returning biological processes...")
  }
  else if (annotation == "GO%3A0005575") {
    cat("Returning cellular components...")
  }
  else if (annotation == "GO%3A0003674") {
    cat("Returning molecular functions...")
  }
  else {
    stop("Improper Gene Ontology annotation type.\nOptions are:\nGO:3A0008150 (biological process)\nGO:3A0005575 (cellular components)\nGO:3A0003674 (molecular functions)")
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
panther_cc <- panther_api(genes, "GO:3A0005575")
panther_mf <- panther_api(genes, "GO%3A0003674")

parse_panther <- function(panther_json) {
  results <- panther_json$results$result
  
  panther_df <- map_dfr(results, ~tibble(
    term_id = .x$term$id,
    term_label = .x$term$label,
    num_genes = .x$number_in_list,
    fold_enrichment = .x$fold_enrichment,
    fdr = .x$fdr,
    pvalue = .x$pValue
  )) %>% filter(num_genes >= 1)
  
  return(panther_df)
}

mf_df <- parse_panther(panther_mf)
bp_df <- parse_panther(panther_bp)
cc_df <- parse_panther(panther_cc)

## KEGG

## step for tomorrow 2/22:
## write a helper function to normalize gene ids to HGNC. First start with rsid matches
## then go to Alias matches to HGNC
## perhaps undo what you did to the xlsx file and start from scratch

## Genes to Entrez IDs for KEGG


## Steps for 2/23:
## Can you write your own code for enrichment analysis of MONDO and HPO?
# Maybe just use this: https://rgd.mcw.edu/wg/new-moet-algorithm/



entrez_ids <-  bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

kk <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

edo <- enrichDO(entrez_ids$ENTREZID)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox,  categorySize="geneNum", layout = "circle")
plot(p1)

edox2 <- pairwise_termsim(edox)
p2 <- treeplot(edox2)
plot(p2)

p3 <- upsetplot(edo)
plot(p3)

## API for HPO and DO enrichment
"https://rest.rgd.mcw.edu/rgdws/enrichment/data"

# MP is phenotype
# MF is molecular function
# RDO is rat disease ontology
# worth a shot to write this code yourself for MONDO and HPO?

# {"species":1,"genes":["BACE1","PTGS2","PACERR","SORL1","TOMM40P3","TOMM40P4","MEF2C-AS2","TOMM40P1","TOMM40P2","MEF2C-AS1","CD33","COX15","EPHA1","SAR1AP1","LRRK2-DT","CELF2-DT","CELF2-AS2","COX10-DT","WASHC5-AS1","SAR1A","BACE1-AS","PIN1-DT","CD36","SOD1-DT","FN1-DT","CELF2-AS1","KANSL1-AS1","DAPK1-IT1","EPHA1-AS1","SORL1-AS1","CIZ1","LRP1-AS","NQO1-DT","DSG2-AS1","ECE1-AS1","IL6-AS1","TRIL","BCCIP","COX10","SAR1AP2","SAR1AP3","SAR1AP4","ADAM10","A2M-AS1","IL10RB-DT","ADAM9","RCOR1","NME8","PANDAR","BNAT1","IL6ST-DT","APP-DT","SOD2-OT1","COX20","ABCA1","HBG2","IL10","IL33","MS4A2","NEDD9","NYAP1","P2RX7","PICALM","PLCG2","PPIAP10","BIN1","TNF","TOMM40","CYP46A1","EGFR","CLU","COMT","ACE","ADAMTS1","ADRA2B","FN1","GAB2","HMGCR","KANSL1","LCORL","LRRK2","MAPT","MTHFR","NQO1","TFCP2","CELF2","ECE1","EXOC3L2","HSPA1A","IL6","LRP11","APOE","ARC","MEF2C","PILRA","PILRB","PON1","PRNP","PTK2B","RBM45","SERPINA3","CCL11","SPI1","TFAM","TREML2","CDK5RAP2","CDKN1A","CHAT","DAPK1","HFE","IDE","IL16","INPP5D","KL","LPL","MS4A6A","NCSTN","NLRP13","NTF3","PLAU","REST","TCN2","TLR4","VAMP1","ZCWPW1","CCR2","CDPF1","CRP","CYP19A1","DSG2","ESR1","ESR2","F5","A2M","ABCA7","GAPDHS","IQCK","LDLR","LRP1","APP","NFIC","NOCT","NRXN3","PIN1","PPP4R3A","PSEN1","RAB10","BCKDK","RELN","RIN3","SLC11A2","SLC24A4","SOD1","CASP7","CASS4","SPPL2A","TMEM106B","TNK1","WASHC5","CETP","CHI3L1","CHRNB2","CYP2C19","KANSL1L-AS1","IL6R-AS1","COX20P2","RCOR2","RCOR3","BACE2","EGFR-AS1","LDLR-AS1","MAPT-IT1","HFE-AS1","MAPT-AS1","ACE2-DT","PIMREG","APTR","IQCF5","SLC35F5","EGILA","TSPOAP1-AS1","ADAMTS19-AS1","LITAF","IGHA2","ADAMTS10","CLUH","TRAP1","ADGRF5","LTA","LTB"],"aspect":"MP","originalSpecies":1}


## Trash pile
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