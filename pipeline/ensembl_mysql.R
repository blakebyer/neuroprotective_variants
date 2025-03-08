# ensembl mysql database
library(DBI)
library(RMySQL)

con <- dbConnect(MySQL(), 
                 user = "anonymous", 
                 password = "", 
                 dbname = "homo_sapiens_variation_110_38", 
                 host = "ensembldb.ensembl.org",
                 port = 3306)

# Check available tables
#dbListTables(con)

snp_string <- paste0("'", unlist(snps), "'", collapse = ",")

query <- paste0("SELECT v.*, vf.*
          FROM variation v
          JOIN variation_feature vf ON v.variation_id = vf.variation_id
          WHERE v.name IN (",snp_string,");")

# Fetch data
snp_data <- dbGetQuery(con, query)


# disconnect when done
dbDisconnect(con)