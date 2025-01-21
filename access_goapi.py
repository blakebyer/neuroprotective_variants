import requests
import json
import pandas as pd

base_url = "https://api.geneontology.org" # API URL
#https://api.geneontology.org/api/bioentity/gene/UniProtKB%3AP34130/function?start=0&rows=100

def search_go(term):
    """Returns a dictionary of all matches to a search term. Use this rarely."""
    url = f"/api/search/entity/autocomplete/{term}"
    response = requests.get(url)
    if response.status_code == 200:
        go_data = response.json()
        return go_data
    else:
        print("Failed to retrieve data.")

def get_subgraph(id):
    """Returns dictionary of descendants and ancestors of a Gene Ontology Term"""
    url = f"{base_url}/api/ontology/term/{id}/subgraph"
    response = requests.get(url) # status code, i.e. 200 'okay' or 404 'Not found'
    
    if response.status_code == 200:
        go_data = response.json()
        return go_data
    else:
        print("Failed to retrieve data.")

def get_gene_functions(id):
    """Returns dictionary of gene functions annotated to the CURIE identifier of the gene. Use UniProtKB:id because it is standard."""
    url = f"{base_url}/api/bioentity/gene/{id}/function"
    response = requests.get(url)

    if response.status_code == 200:
        go_data = response.json()
        df = pd.json_normalize(go_data, 'associations') # flattens json to dataframe
        return df
    else:
        print("Failed to retrieve data.")

## Example usage

# id = "GO:0021675" # neurodevelopment
# go_info = get_subgraph(id)

uniprot_id = "UniProtKB:P34130"

gene_func = get_gene_functions(uniprot_id) # dataframe

human_processes = gene_func[ # filter for human results 
    (gene_func['object.taxon.id'] == "NCBITaxon:9606") & 
    ((gene_func['object.category'].str[0] == "molecular_activity") | # molecular activity or biological process
     (gene_func['object.category'].str[0] == "biological_process"))]

final_df = gene_func[["subject.id", "subject.label", "object.id", "object.label","object.taxon.id", "object.taxon.label", "reference"]]
print(final_df)