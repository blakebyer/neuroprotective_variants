library(tidyverse)
library(readxl)
library(httr2)
library(jsonlite)
library(purrr)
library(ggupset)
library(dendextend)
library(shinyhttr)

# ## Read Variant and Gene Data
# var_data <- read_xlsx("neuroprotective_variants.xlsx", sheet = "compiled variants")
# snps <- var_data %>% dplyr::select(rsid) %>% distinct() %>% drop_na() %>% as.vector()
# genes <- var_data %>% dplyr::select(gene) %>% distinct() %>% drop_na() %>% as.vector()

### Normalize Gene Identifiers
normalize_genes <- function(genes) {
  genes_vec <- unlist(genes)
  len <- length(genes_vec)
  genes_string <- paste(genes_vec, collapse = "+OR+")

  base_url <- "https://clinicaltables.nlm.nih.gov"
  hgnc_query <- sprintf("/api/genes/v4/search?terms=%s&df=hgnc_id,hgnc_id_num,symbol,name,alias_symbol,location&maxList=%s&count=%s", genes_string, 500, 500)

  get_hgnc <- request(base_url) |>
      req_url_path(path = hgnc_query) |>
       req_headers(accept = "application/json")|>
      req_method("GET")

  resp_hgnc <- req_perform(get_hgnc) |>
    resp_body_json()

  gene_df <- map_dfr(resp_hgnc[[4]], ~tibble(
          hgnc_id = as.character(.x[1]),
          hgnc_id_num = as.numeric(.x[2]),
          symbol = as.character(.x[3]),
          name = as.character(.x[4]),
          alias_symbol = as.character(.x[5]),
          location = as.character(.x[6])
        )) %>% mutate(search_term = NA_character_)

  if (nrow(gene_df) == 0) {
    stop("No mapped genes found. Try again with valid symbols.")
  }

  # return matches
  for (i in seq_along(gene_df$symbol)) {
    matched_terms <- c()
    for (column in c("symbol", "alias_symbol", "name")) {
      for (term in genes_vec) {
        match <- str_extract(gene_df[[column]][i], fixed(term, ignore_case = TRUE))
        if (!is.na(match)) {
          matched_terms <- c(matched_terms, term)
        }
      }
    }
    if (length(matched_terms) > 0) {
      gene_df$search_term[i] <- paste(unique(matched_terms), collapse = "|")
    }
  }

  hgnc_genes <- paste(gene_df$symbol, collapse = "+OR+")
  ncbi_query <- sprintf("/api/ncbi_genes/v3/search?terms=%s&df=GeneID,HGNC_ID&maxList=%s&count=%s", hgnc_genes, 500, 500)

  # request ncbi gene_ids
  get_ncbi <- request(base_url) |>
    req_url_path(path = ncbi_query) |>
    req_headers(accept = "application/json")|>
    req_method("GET")

  resp_ncbi <- req_perform(get_ncbi) |>
    resp_body_json()

  ncbi_ids <- map_dfr(resp_ncbi[[4]], ~tibble(
     ncbi_id = as.numeric(.x[1]),
     hgnc_id = as.character(.x[2])))

   gene_df <- gene_df %>%
     filter(!is.na(search_term)) %>%
     left_join(ncbi_ids, by = join_by(hgnc_id))

  # # mapped genes
   mapped_genes <- unique(unlist(str_split(gene_df$search_term, "\\|")))
   n_mapped <- length(mapped_genes)
   missing_genes <- setdiff(genes_vec, mapped_genes)
   n_missing <- length(missing_genes)

  # # missing genes
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

#n <- normalize_genes(genes)

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

## Gene Ontology
# "id": "GO%3A0003674",
# "label": "molecular_function"
# 
# "id": "GO%3A0008150",
# "label": "biological_process"

# "id": "GO%3A0005575",
# "label": "cellular_component"

panther_api <- function(genes, annotation) {
  if (annotation == "BP") {
    type <- "GO%3A0008150"
    cat("Returning biological processes...")
  }
  else if (annotation == "CC") {
    type <- "GO%3A0005575"
    cat("Returning cellular components...")
  }
  else if (annotation == "MF") {
    type <- "GO%3A0003674"
    cat("Returning molecular functions...")
  }
  else {
    stop("Improper Gene Ontology annotation type.\nOptions are:\nGO:3A0008150 (biological process)\nGO:3A0005575 (cellular components)\nGO:3A0003674 (molecular functions)")
  }
  genes <- paste(unlist(genes), collapse = ",")
  
  # Define base and query url
  base_url <- "https://www.pantherdb.org"
  query_url <- sprintf("/services/oai/pantherdb/enrich/overrep?geneInputList=%s&organism=9606&annotDataSet=%s&enrichmentTestType=FISHER&correction=FDR", genes, type)
  
  # Make an HTTP request
  post_panther <- request(base_url) |>
    req_url_path(path = query_url) |>
    req_headers(accept = "application/json") |>
    req_method("POST") |>
    req_progress(type = "down")

  resp_panther <- req_perform(post_panther) |>
    resp_body_json()

  return(resp_panther)
}

parse_panther <- function(panther_json) {
  results <- panther_json$results$result
  
  panther_df <- map_dfr(results, ~tibble(
    term_id = .x$term$id,
    term_label = .x$term$label,
    num_genes = .x$number_in_list,
    fold_enrichment = .x$fold_enrichment,
    fdr = .x$fdr,
    pvalue = .x$pValue
  )) %>% filter(num_genes >= 1, pvalue < 0.05)
  
  return(panther_df)
}

## Gene Set Enrichment
gsea <- function(genes, ontology) {
  
  if (!ontology %in% c("MF", "BP", "CC", "MP", "PW", "CHEBI", "RDO")) {
    stop("Improper ontology annotation type.\nOptions are:\nBP (GO biological process)\nCC (GO cellular components)\nMF (GO molecular functions)\nMP (Human Phenotype Ontology)\nRDO (Disease Ontology)\nPW (Pathway Ontology)\nCHEBI (Chemical Entities of Biological Interest)")
  }
  
  query <- list(
    species = 1,
    genes = genes,
    aspect = sprintf("%s", ontology),
    originalSpecies = 1
  )
  
  base_url <- "https://rest.rgd.mcw.edu"
  ext <- "/rgdws/enrichment/data"
  
  post_gsea <- request(base_url) |>
    req_url_path(path = ext) |>
    req_headers(accept = "application/json",
                `Content-Type`="application/json") |>
    req_body_json(query) |>
    req_method("POST")
  
  resp_gsea <- req_perform(post_gsea) |>
    resp_body_json()
  
  return(resp_gsea)
}

parse_gsea <- function(gsea_json) {
  results <- gsea_json$enrichment
  
  gsea_df <- map_dfr(results, ~tibble(
    term_id = as.character(.x$acc),
    term = as.character(.x$term),
    genes = as.numeric(.x$count),
    ref_genes = as.numeric(.x$refCount),
    odds_ratio = as.numeric(.x$oddsratio),
    corrected_pvalue = as.numeric(.x$correctedpvalue),
    pvalue = as.numeric(.x$pvalue)
  )) %>% filter(genes >= 1)
  
  return(gsea_df)
}

hierarchy <- function(gsea_df, k) {
  if ("odds_ratio" %in% colnames(gsea_df)) {
    sets <- gsea_df %>%
      filter(corrected_pvalue < 0.05) %>%
      select(term, odds_ratio) %>%
      mutate(scaled = log(odds_ratio + 1))
  }
  else if ("fold_enrichment" %in% colnames(gsea_df)) {
    sets <- gsea_df %>%
      filter(pvalue < 0.05) %>%
      select(term_label, fold_enrichment) %>%
      mutate(scaled = log(fold_enrichment + 1), 
             term = term_label) %>%
      select(-term_label)
  }
  
  dist_matrix <- dist(sets$scaled, method = "euclidean")
  
  hc <- hclust(dist_matrix, method = "complete")
  
  dend <- as.dendrogram(hc) %>%
    set("labels", sets$term) %>%
    set("branches_k_color", k=k) %>%
    color_labels(k = k) %>% 
    set("labels_cex", 0.5)
  
  # plot entire dendrogram
  dend_plot <- as.ggdend(dend)
  
  j <- ggplot(dend_plot, labels = FALSE) + 
    theme_void()
  
  plot(j)
  
  # assign each path of tree to k clusters
  clusters <- cutree(dend, k=k)
  
  # add clusters to sets
  sets$cluster <- clusters
  
  # height to split graph into k components
  height = as.vector(heights_per_k.dendrogram(dend))[[k]]
  
  # cut trees to expose leaf terms
  sub.trees <- cut(dend, h = height, k = k)
  
  n_cuts <- length(sub.trees$lower)

  # Select lower clusters and plot
  if (length(sub.trees$lower) > 0) {
    for (i in 1:n_cuts) {
      cluster_k <- sub.trees$lower[[i]]
      
      # If the subtree has more than 50 leaves, remove labels
      if (length(labels(cluster_k)) > 50) {
        cluster_k <- cluster_k %>% set("labels", NULL)
      }
      
      # Convert to ggdend object safely
      cluster_k <- tryCatch({
        as.ggdend(cluster_k)
      }, error = function(e) {
        message("Skipping invalid subtree. Try a smaller value of k: ", e$message)
        return(NULL)
      })
      
      # Skip if conversion failed
      if (is.null(cluster_k)) next
      
      p <- ggplot(cluster_k, horiz = TRUE) +
        theme_void()
      
      plot(p)
    }
  }
  return(sets)
  
}