import os
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

def load_SLIDE_res(folder_path):
    """Load SLIDE restulst as a dictionary.
        Args:
            folder_path (str): the path to the folder containing all SLIDE results.
    """
    latent_factors = {}  
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'): 
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, delimiter='\t')
            latent_factors[filename[0:-4]] = df
    return latent_factors

def gene_dict_to_list(gene_dict):
    gene_list = []
    for df in gene_dict.values():
        gene_list.extend(df.iloc[:, -1].tolist())
    gene_list = list(set(gene_list))
    return gene_list

def plot_gene_hist(adata, goi, save_folder):
    sc.get.obs_df(adata, keys=goi, layer="imputed_count").hist()
    plt.show()
    plt.savefig(f"{save_folder}/{goi}_hist.pdf")