import requests
import json
import pandas as pd

base_url = "https://api.geneontology.org" # API URL
#https://api.geneontology.org/api/bioentity/gene/UniProtKB%3AP34130/function?start=0&rows=100

def vep(variants):
    """Returns a dictionary of all variant effect predictions"""
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
        return go_data # returns dictionary
    else:
        print("Failed to retrieve data.")

def get_gene_functions(id):
    """Returns dataframe of gene functions annotated to the CURIE identifier of the gene. Use UniProtKB:id because it is standard."""
    url = f"{base_url}/api/bioentity/gene/{id}/function"
    response = requests.get(url)

    if response.status_code == 200:
        go_data = response.json()
        df = pd.json_normalize(go_data, 'associations') # flattens json to dataframe
        return df
    else:
        print("Failed to retrieve data.")
