import pandas as pd
import os

interaction_table = pd.read_csv("CCC/CellPhoneDB/interaction_table.csv")
protein_table = pd.read_csv("CCC/CellPhoneDB/protein_table.csv")
#%%

# Load in the full SLIDE model for macrophages
all_LFs = [] # a list of all latent factors
folder_path = "ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results"
txt_files = [file for file in os.listdir(folder_path) if file.endswith(".txt")]
for file in txt_files:
    file_path = os.path.join(folder_path, file)
    # read the file using pandas
    df = pd.read_csv(file_path, sep="\t")
    # add a column with the name of the file
    df['LF_name'] = file[:-4]
    df['cell_type'] = 'Macrophage'
    all_LFs.append(df)

LF_df_mac = pd.concat(all_LFs, ignore_index=True)


#%%
# Load in the full SLIDE model for activated fibroblasts
all_LFs = [] # a list of all latent factors
folder_path = "ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final"
txt_files = [file for file in os.listdir(folder_path) if file.endswith(".txt")]
for file in txt_files:
    file_path = os.path.join(folder_path, file)
    # read the file using pandas
    df = pd.read_csv(file_path, sep="\t")
    # add a column with the name of the file
    df['LF_name'] = file[:-4]
    df['cell_type'] = 'AF'
    all_LFs.append(df)

LF_df_AF = pd.concat(all_LFs, ignore_index=True)

#%%
LF_df = pd.concat([LF_df_mac, LF_df_AF], ignore_index=True)
LF_df.to_csv("CCC/All_Mac_AF_CCR2_LFs.txt", sep="\t", index=False)
#%%
LF_gene_names = LF_df['names'].tolist()
# add _HUMAN to the names
LF_gene_names = [name + "_HUMAN" for name in LF_gene_names]
overlap_protein = protein_table[protein_table['protein_name'].isin(LF_gene_names)]


# ccr2_lf = pd.concat([z_37, z_27, z_21, z_17, z_3])

# ccr2_overlap = human_LR[(human_LR['Ligand'].isin(ccr2_lf['names'])) | (human_LR['Receptor'].isin(ccr2_lf['names']))]
