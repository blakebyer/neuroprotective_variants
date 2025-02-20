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
molF['Fold Enrichment'] = molF['Fold Enrichment'].apply(clean_fold_enrichment).sort_values(ascending=False)

molF['Genes'] = pd.to_numeric(molF['Genes'], errors='coerce')

# fig, ax = plt.subplots()
# molF.plot(kind='barh', x='Molecular Function', y='Fold Enrichment', legend=False)
# molF.plot.scatter(x = 'Genes', y = 'Molecular Function')
# plt.subplots_adjust(left=0.6) 
# plt.yticks(fontsize=8)
# plt.ylabel("Molecular Function")
# plt.xlabel("Fold Enrichment")
# plt.show()  

# Create subplots
fig, ax = plt.subplots(figsize=(8, 6))

# Plot horizontal bars
ax.barh(molF["Molecular Function"], molF["Fold Enrichment"], color='skyblue', label="Fold Enrichment")

# Overlay scatter plot (points at end of bars, sized by Gene Count)
ax.scatter(molF["Fold Enrichment"], molF["Molecular Function"], 
           s=molF["Genes"] * 6,  # Scale point size for visibility
           color='red', edgecolors='black', label="Genes")

# Formatting
plt.subplots_adjust(left=0.6)  # Adjust spacing for labels
ax.set_ylabel("Molecular Function")
ax.set_xlabel("Fold Enrichment")
ax.legend()  # Show legend

plt.show()



