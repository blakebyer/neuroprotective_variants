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

def go_enrichment_plot(dataframe, go_type, top_n=None):
    # sort fold enrichment descending 
    dataframe = dataframe.copy()
    dataframe['Fold Enrichment'] = dataframe['Fold Enrichment'].apply(clean_fold_enrichment).sort_values(ascending=True)[:top_n]

    dataframe['Genes'] = pd.to_numeric(dataframe['Genes'], errors='coerce')

    # Create subplots
    fig, ax = plt.subplots(figsize=(10,8))

    # Plot horizontal bars
    ax.barh(dataframe[go_type], dataframe["Fold Enrichment"], color="skyblue", label="Fold Enrichment", height=0.2)

    # Overlay scatter plot (points at end of bars, sized by number of genes
    ax.scatter(dataframe["Fold Enrichment"], dataframe[go_type], 
            s=dataframe["Genes"] * 5,  # Scale point size for visibility
            color='orangered', edgecolors='black', label="Genes")

    # Formatting
    plt.subplots_adjust(left=0.6)  # Adjust spacing for labels
    ax.set_ylabel(go_type)
    ax.set_xlabel("Fold Enrichment")
    ax.legend()  # Show legend
    plt.show()

def go_volcano_plot(dataframe, go_type, color):
    dataframe = dataframe.copy()
    dataframe['Fold Enrichment'] = dataframe['Fold Enrichment'].apply(clean_fold_enrichment).sort_values(ascending=False)
    
    dataframe["FDR"] = pd.to_numeric(dataframe["FDR"], errors='coerce')

    dataframe["FDR"] = np.log10(np.add(dataframe["FDR"], 1))

    dataframe["P-value"] = pd.to_numeric(dataframe["P-value"], errors='coerce')

    dataframe["P-value"] = np.log10(np.add(dataframe["P-value"], 1))

    plt.scatter(dataframe["Fold Enrichment"], dataframe["P-value"], color = color)
    plt.xlabel(f"{go_type} Fold Enrichment")
    plt.ylabel("log(P-value)")
    plt.show()

def volcano_overlay(dataframe, colors):
    dataframe = dataframe.copy()
    dataframe['Fold Enrichment'] = dataframe['Fold Enrichment'].apply(clean_fold_enrichment).sort_values(ascending=False)
    
    dataframe["FDR"] = pd.to_numeric(dataframe["FDR"], errors='coerce')

    dataframe["FDR"] = np.log10(np.add(dataframe["FDR"], 1))

    dataframe["P-value"] = pd.to_numeric(dataframe["P-value"], errors='coerce')

    dataframe["P-value"] = np.log10(np.add(dataframe["P-value"], 1))
    # Create subplots
    fig, ax = plt.subplots(figsize=(5,7))

    # Plot horizontal bars
    ax.barh(dataframe[go_type], dataframe["Fold Enrichment"], color="skyblue", label="Fold Enrichment", height=0.2)

    # Overlay scatter plot (points at end of bars, sized by number of genes
    ax.scatter(dataframe["Fold Enrichment"], dataframe[go_type], 
            s=dataframe["Genes"] * 5,  # Scale point size for visibility
            color='orangered', edgecolors='black', label="Genes")

    # Formatting
    plt.subplots_adjust(left=0.6)  # Adjust spacing for labels
    ax.set_ylabel(go_type)
    ax.set_xlabel("Fold Enrichment")
    ax.legend()  # Show legend
    plt.show()


molF_plot = go_enrichment_plot(molF, "Molecular Function")

bioP_plot = go_enrichment_plot(bioP, "Biological Process", 50) # major overplotting issues, so select bottom 50 enrichment scores (many genes express this process)
bioP_volcano = go_volcano_plot(bioP, "Biological Process", "teal")

celC_plot = go_enrichment_plot(celC, "Cellular Component")

molF_volcano = go_volcano_plot(molF, "Molecular Function", "orangered")

celC_volcano = go_volcano_plot(celC, "Cellular Component", "violet")


