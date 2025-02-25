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
normalize_genes <- function(genes, session) {
  genes_vec <- unlist(genes)
  len <- length(unique(genes_vec))
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
    updateProgressBar(session, id = "get_hgnc", value = 0)
    stop("No mapped diseases found.")
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
    req_headers(accept = "application/json") |>
    req_method("GET")

  resp_ncbi <- req_perform(get_ncbi) |>
    resp_body_json()

  ncbi_ids <- map_dfr(resp_ncbi[[4]], ~tibble(
     ncbi_id = as.numeric(.x[1]),
     hgnc_id = as.character(.x[2])))

   gene_df <- gene_df %>%
     filter(!is.na(search_term)) %>%
     left_join(ncbi_ids)

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
    req_method("POST") |>
    req_throttle(45/60) |>
    req_retry()
  
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
    #make request string
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

parse_omim <- function(omim_json, session) {
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
  
  if (nrow(omim_df) == 0) {
    updateProgressBar(session, id = "get_omim", value = 0)
    stop("No mapped diseases found.")
  }
  
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
    req_throttle(45/60) |>
    req_retry()

  resp_panther <- req_perform(post_panther) |>
    resp_body_json()

  return(resp_panther)
}

parse_panther <- function(panther_json, session) {
  results <- panther_json$results$result
  
  panther_df <- map_dfr(results, ~tibble(
    term_id = .x$term$id,
    term_label = .x$term$label,
    num_genes = .x$number_in_list,
    fold_enrichment = .x$fold_enrichment,
    fdr = .x$fdr,
    pvalue = .x$pValue
  )) %>% filter(num_genes >= 1, pvalue < 0.05, !is.na(fold_enrichment))
  
  if (nrow(panther_df) == 0) {
    updateProgressBar(session, id = "get_go", value = 0)
    stop("No mapped diseases found.")
  }
  
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
    req_method("POST") |>
    req_throttle(45/60) |>
    req_retry()
  
  resp_gsea <- req_perform(post_gsea) |>
    resp_body_json()
  
  return(resp_gsea)
}

parse_gsea <- function(gsea_json, session) {
  results <- gsea_json$enrichment
  
  gsea_df <- map_dfr(results, ~tibble(
    term_id = as.character(.x$acc),
    term = as.character(.x$term),
    genes = as.numeric(.x$count),
    ref_genes = as.numeric(.x$refCount),
    odds_ratio = as.numeric(.x$oddsratio),
    corrected_pvalue = as.numeric(.x$correctedpvalue),
    pvalue = as.numeric(.x$pvalue)
  )) %>% filter(genes >= 1, corrected_pvalue < 0.05, !is.na(odds_ratio))
  
  if (nrow(gsea_df) == 0) {
    updateProgressBar(session, id = "get_gsea", value = 0)
    stop("No mapped annotations found.")
  }
  
  return(gsea_df)
}

hierarchy <- function(gsea_df, k, n_plot) {
  if ("odds_ratio" %in% colnames(gsea_df)) {
    sets <- gsea_df %>%
      select(term, odds_ratio) %>%
      mutate(scaled = log(odds_ratio + 1)) %>%
      select(term, scaled)
  } else if ("fold_enrichment" %in% colnames(gsea_df)) {
    sets <- gsea_df %>%
      select(term_label, fold_enrichment) %>%
      mutate(scaled = log(fold_enrichment + 1), 
             term = term_label) %>%
      select(term, scaled)
  }
  
  # Ensure k is within the valid range
  max_clusters <- max(1, length(gsea_df) - 1)
  if (k > max_clusters) {
    warning(paste("k exceeds the valid range. Resetting to", max_clusters))
    k <- max_clusters
  } else if (k < 1) {
    stop("k must be at least 1")
  }
  
  # Ensure n_plot is between 0 and k
  if (n_plot < 0 || n_plot > k) {
    warning(paste("n_plot out of range [0,", k, "]. Resetting to 0"))
    n_plot <- 0
  }
  
  dist_matrix <- dist(sets$scaled, method = "euclidean")
  hc <- hclust(dist_matrix, method = "complete")
  
  sets$wrapped_term <- str_wrap(sets$term, width = 20)
  
  dend <- as.dendrogram(hc) %>%
    set("labels", sets$wrapped_term) %>%
    set("branches_k_color", k = k) %>%
    color_labels(k = k) %>%
    set("labels_cex", 1)
  
  # Plot entire dendrogram if `n_plot == 0`
  if (n_plot == 0) {
    dend_plot <- as.ggdend(dend)
    j <- ggplot(dend_plot, labels = FALSE) + theme_void()
    return(j)  # Return ggplot object
  }
  
  # Assign clusters
  clusters <- cutree(dend, k = k)
  sets$cluster <- clusters
  
  # Get cutting height
  height <- as.vector(heights_per_k.dendrogram(dend))[[k]]
  sub.trees <- cut(dend, h = height, k = k)
  
  # Check if sub.trees$lower exists and contains trees
  if (!is.null(sub.trees$lower) && length(sub.trees$lower) >= n_plot) {
    cluster_k <- sub.trees$lower[[n_plot]]
    
    if (length(labels(cluster_k)) > 20) {
      cluster_k <- cluster_k %>% set("labels", NULL)
    }
    
    cluster_k <- tryCatch({
      as.ggdend(cluster_k)
    }, error = function(e) {
      warning("Skipping invalid subtree. Try a smaller value of k: ", e$message)
      return(NULL)
    })
    
    if (!is.null(cluster_k)) {
      p <- ggplot(cluster_k, horiz = TRUE) + theme_void() + 
        scale_y_continuous(expand = expansion(mult = c(0,0.2)), transform = "reverse") +
        scale_x_continuous(expand = expansion(mult = 0.1))
      return(p)
    }
  }
  
  warning("No valid subtree found. Try adjusting k or n_plot.")
  
  i <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "Warning: No valid subtree found.\nTry decreasing K.", hjust = 0.5, size = 6) +
    theme_void()
  return(i)
}

panther_plot <- function(panther_df, slice) {
  if ("odds_ratio" %in% colnames(panther_df)) {
    panther_df <- panther_df %>%
      slice_max(n = slice, order_by = odds_ratio)
    i <- ggplot(panther_df) + 
      geom_col(aes(x = odds_ratio, y = reorder(term, -odds_ratio), fill = corrected_pvalue)) + 
      labs(x = "Odds Ratio", y = "Term", fill = "P-value") + 
      scale_fill_viridis_c(option = "plasma") + 
      theme_minimal()
    return(i)
  } else if ("fold_enrichment" %in% colnames(panther_df)) {
    panther_df <- panther_df %>%
      slice_max(n = slice, order_by = fold_enrichment)
    j <- ggplot(panther_df) + 
      geom_col(aes(x = fold_enrichment, y = reorder(term_label, -fold_enrichment), fill = pvalue)) + 
      labs(x = "Fold Enrichment", y = "Term", fill = "P-value") + 
      scale_fill_viridis_c(option = "plasma") + 
      theme_minimal()
    return(j)
  }
  else {return(NULL)}
}
