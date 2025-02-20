import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

bioP = pd.read_excel("Neuroprotective Genes.xlsx", sheet_name="Biological Process")[6:] # only after header

molF = pd.read_excel("Neuroprotective Genes.xlsx", sheet_name= "Molecular Function")[6:]

celC = pd.read_excel("Neuroprotective Genes.xlsx", sheet_name= "Cellular Component")[6:]

molF = molF.rename(columns = {'Analysis Type:':'Molecular Function', 'PANTHER Overrepresentation Test (Released 20240807)':'Homo Sapiens Reflist', 'Unnamed: 2':'Genes', 'Unnamed: 3':'Expected', 'Unnamed: 4': 'Over/Under', 'Unnamed: 5':'Fold Enrichment', 'Unnamed: 6':'P-value', 
                       'Unnamed: 7': 'FDR'})

bioP = bioP.rename(columns = {'Analysis Type:':'Biological Process', 'PANTHER Overrepresentation Test (Released 20240807)':'Homo Sapiens Reflist', 'Unnamed: 2':'Genes', 'Unnamed: 3':'Expected', 'Unnamed: 4': 'Over/Under', 'Unnamed: 5':'Fold Enrichment', 'Unnamed: 6':'P-value', 
                       'Unnamed: 7': 'FDR'})

celC = celC.rename(columns = {'Analysis Type:':'Cellular Component', 'PANTHER Overrepresentation Test (Released 20240807)':'Homo Sapiens Reflist', 'Unnamed: 2':'Genes', 'Unnamed: 3':'Expected', 'Unnamed: 4': 'Over/Under', 'Unnamed: 5':'Fold Enrichment', 'Unnamed: 6':'P-value', 
                       'Unnamed: 7': 'FDR'})

def clean_fold_enrichment(value):
    if isinstance(value, str) and value.startswith(">"):
        return 100  # Convert ">100" to 100
    return pd.to_numeric(value, errors='coerce')  # Convert other values normally

# Apply function to the column
molF['Fold Enrichment'] = molF['Fold Enrichment'].apply(clean_fold_enrichment).dropna().sort_values(ascending=False)

molF.plot(kind='barh', x='Molecular Function', y='Fold Enrichment', legend=False)
plt.ylabel("Molecular Function")
plt.xlabel("Fold Enrichment")
plt.show()  



