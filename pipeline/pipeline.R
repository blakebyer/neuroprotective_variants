library(tidyverse)
library(readxl)
library(httr2)
library(jsonlite)
library(purrr)
library(ggupset)
library(dendextend)

## Read Variant and Gene Data
var_data <- read_xlsx("neuroprotective_variants.xlsx", sheet = "compiled variants")
snps <- var_data %>% dplyr::select(rsid) %>% distinct() %>% drop_na() %>% as.vector()
genes <- var_data %>% dplyr::select(gene) %>% distinct() %>% drop_na() %>% as.vector()

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
     HGNC_ID = as.character(.x[2])))

   gene_df <- gene_df %>%
     filter(!is.na(search_term)) %>%
     left_join(ncbi_ids, by = join_by(hgnc_id == HGNC_ID))

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

hgnc <- normalize_genes(genes)
genes <- hgnc$symbol

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

ontology = "RDO"

g <- gsea(genes, ontology)

g_df <- parse_gsea(g)

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

h_plots <- hierarchy(g_df, 4)

h_plot2 <- hierarchy(mf_df, 10)


## Trash pile

## Hierarchical clustering
# sets <- g_df %>%
#   dplyr::filter(corrected_pvalue < 0.05) %>%
#   dplyr::select(term, odds_ratio) %>%
#   mutate(scaled_or = log(odds_ratio + 1))
# 
# rownames(sets) <- sets$term
# 
# dist_matrix <- dist(sets$scaled_or, method = "euclidean")
# 
# hc <- hclust(dist_matrix, method = "complete")
# 
# plot(hc)
# 
# dend <- as.dendrogram(hc) %>%
#   set("labels", sets$term)
# 
# 
# plot(dend, leaflab = "none")
# 
# clusters <- cutree(dend, k=5)
# sets$cluster <- clusters
# 
# tail(sets, 10)
# 
# plot(color_branches(dend, k=5),leaflab="none")
# sub.trees <- cut(dend, h = 1)
# 
# cluster1.tree <- sub.trees$lower[[2]]
# 
# cluster1.tree  %>%
#   set("labels_cex", 0.5) %>%
#   set("labels_col", "tomato") %>%
#   plot(horiz = TRUE)  # plot horizontally


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
# entrez_ids <-  bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# 
# kk <- enrichKEGG(gene         = entrez_ids$ENTREZID,
#                  organism     = 'hsa',
#                  pvalueCutoff = 0.05)
# 
# edo <- enrichDO(entrez_ids$ENTREZID)
# edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
# p1 <- cnetplot(edox,  categorySize="geneNum", layout = "circle")
# plot(p1)
# 
# edox2 <- pairwise_termsim(edox)
# p2 <- treeplot(edox2)
# plot(p2)
# 
# p3 <- upsetplot(edo)
# plot(p3)

# normalize_genes <- function(genes) {
#   genes_vec <- unlist(genes)
#   len <- length(genes_vec)
#   genes_string <- paste(genes_vec, collapse = "+OR+")
# 
#   base_url <- "https://rest.genenames.org"
#   hgnc_query <- sprintf("/search/symbol:%s", genes_string)
# 
#   get_hgnc <- request(base_url) |>
#       req_url_path(path = hgnc_query) |>
#       req_headers(accept = "application/json") |>
#       req_method("GET")
# 
#   resp_hgnc <- req_perform(get_hgnc) |>
#     resp_body_json()
# 
#   gene_df <- map_dfr(resp_hgnc$response$docs, ~tibble(
#           symbol = as.character(.x[1]),
#           hgnc_id = as.character(.x[2]),
#           score = as.numeric(.x[3])))
#   
#   return(gene_df)
# }
# 
# new <- normalize_genes(genes)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(biomaRt)
# library(DOSE)
# library(enrichplot)
# prompt user to pick appropriate k
# if (n_clusters < k) {
#   warning(sprintf("Number of clusters assigned by hclust is less than k\n. Number of clusters: %s\nPlease set k <= %s.", n_clusters, n_clusters))
# }