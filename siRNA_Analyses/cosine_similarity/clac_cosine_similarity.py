import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

def split_data_by_TF (data, human_TFs):

    data['Target'] = data['Target'].str[2:].str.upper()
    # find the intersection of the TFs in the data and the TF database
    TF_data = data[data['Target'].isin(human_TFs[0].tolist())]
    nonTF_data = data[~data['Target'].isin(human_TFs[0].tolist())]

    print(f"No. of all targets in data data: {len(data['Target'].unique())}")
    print(f"No. of TFs in the data: {len(TF_data['Target'].unique())}")
    print(f"The TFs in the data are: {TF_data['Target'].unique()}")
    return TF_data, nonTF_data

def cosine_similarity(v1, v2):
    return np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

def calc_pariwsie_cosine_similarity(TF_data, nonTF_data):
    sim_dict = {}
    sim_avg_df = pd.DataFrame(
        index = nonTF_data['Target'].unique(), 
        columns = TF_data['Target'].unique()
    )
    for TF in TF_data['Target'].unique():
        for nonTF in nonTF_data['Target'].unique():
            print(f"TF: {TF}, nonTF: {nonTF}")
            subset_TF = TF_data[TF_data['Target'] == TF]
            subset_nonTF = nonTF_data[nonTF_data['Target'] == nonTF]
            sim_lst = []
            for i in subset_TF.index:
                for j in subset_nonTF.index:
                    # subset the i row and j row from the data
                    sim = cosine_similarity(subset_TF.loc[i][1:].values, subset_nonTF.loc[j][1:].values)
                    sim_lst.append(sim)
            sim_dict[(TF, nonTF)] = sim_lst
            sim_avg_df.loc[nonTF, TF] = np.mean(sim_lst)
    return sim_dict, sim_avg_df


def plot_heatmap(sim_avg_df: pd.DataFrame, output_path: str):
    # Create a custom colormap to match your Plotly colorscale
    colors = [
        'rgb(240, 240, 255)',  # Very light purple
        'rgb(200, 180, 240)',  # Light purple
        'rgb(120, 100, 180)',  # Medium purple
        'rgb(80, 40, 140)',    # Darker purple
        'rgb(40, 0, 80)'       # Very dark purple
    ]

    
    colors_rgb = []
    for color in colors:
        rgb = color.replace('rgb(', '').replace(')', '').split(',')
        rgb_normalized = [float(c)/255 for c in rgb]
        colors_rgb.append(tuple(rgb_normalized))

    
    positions = [0, 0.3, 0.6, 0.8, 1]

    # Create the custom colormap
    custom_cmap = mcolors.LinearSegmentedColormap.from_list(
        'custom_purple', 
        list(zip(positions, colors_rgb))
    )
    cell_size = 0.5
    
    # Calculate figure dimensions that maintain fixed cell size
    width_inches = sim_avg_df.shape[1] * cell_size + 2  # Add space for colorbar
    height_inches = sim_avg_df.shape[0] * cell_size + 1.5  # Add space for title and labels
    
    # Create figure with calculated dimensions
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    # # Use seaborn for better control over heatmap appearance
    ax = sns.heatmap(
        sim_avg_df.astype(float), 
        cmap=custom_cmap,
        vmin=0.7, 
        vmax=1.0,
        annot=False,  # Set to True if you want to see the values
        linewidths=0.3,  # Add grid lines between cells
        linecolor='white',  # Make grid lines white
        cbar_kws={'label': 'Cosine Similarity'},
        ax = ax,
        square = True,
    )

    # Set labels and title
    
    ax.xaxis.set_label_position('top') 
    ax.xaxis.set_ticks_position('top')
    ax.set_aspect('equal')
    plt.title('Cosine Similarity Heatmap between Genes and Transcription Factors', pad=10)
    plt.ylabel('Genes')
    plt.xlabel('Transcription Factors')

    # Save directly to PDF
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

#%% #########################################################################################
'''
qPCR data, TF vs non-TF cosine similarity
'''
# read in the TF database
human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
data = pd.read_excel("siRNA_Analyses/Raw_Data/reformatted_qPCR.xlsx")
data = data.rename(columns = {'Unnamed: 0':'Target'})

# drop all rows with Target = 'NT' or 'Tgfb'
data = data[(data['Target'] != 'NT') & (data['Target'] != 'Tgfb')]
# delete the first 2 elements in data['target'] to match the format of human_TFs
TF_data, nonTF_data = split_data_by_TF(data, human_TFs)

sim_dict, sim_avg_df = calc_pariwsie_cosine_similarity(TF_data, nonTF_data)   
# save the similarity matrix
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/qPCR_TF_nonTF_cosine_similarity.csv')

output_path = 'siRNA_Analyses/cosine_similarity/qPCR_TF_nonTF_cosine_similarity.pdf'
plot_heatmap(sim_avg_df, output_path)

#%% #########################################################################################
'''
qPCR data, TF vs  NT/TGFb cosine similarity
'''
human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
data = pd.read_excel("siRNA_Analyses/Raw_Data/reformatted_qPCR.xlsx")
data = data.rename(columns = {'Unnamed: 0':'Target'})

# drop all rows with Target = 'NT' or 'Tgfb'
control_data = data[data['Target'].isin(['NT', 'Tgfb'])]
# delete the first 2 elements in data['target'] to match the format of human_TFs
TF_data, nonTF_data = split_data_by_TF(data, human_TFs)

sim_dict, sim_avg_df = calc_pariwsie_cosine_similarity(TF_data, control_data)   
# save the similarity matrix
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/qPCR_TF_control_cosine_similarity.csv')

# visualize the similarity with a heatmap
output_path = 'siRNA_Analyses/cosine_similarity/qPCR_TF_control_cosine_similarity.pdf'
plot_heatmap(sim_avg_df, output_path)
# %% ############################################################################################################
'''
flow_freq data, TF vs  nonTF cosine similarity
'''

human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
data = pd.read_excel("siRNA_Analyses/Raw_Data/reformatted_flow_freq.xlsx")
data = data.rename(columns = {'Unnamed: 0':'Target'})

# drop all rows with Target = 'NT' or 'Tgfb'
data = data[(data['Target'] != 'NT') & (data['Target'] != 'Tgfb')]
# delete the first 2 elements in data['target'] to match the format of human_TFs
TF_data, nonTF_data = split_data_by_TF(data, human_TFs)

sim_dict, sim_avg_df = calc_pariwsie_cosine_similarity(TF_data, nonTF_data)   
# save the similarity matrix
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/freq_TF_nonTF_cosine_similarity.csv')

output_path = 'siRNA_Analyses/cosine_similarity/freq_TF_nonTF_cosine_similarity.pdf'
plot_heatmap(sim_avg_df, output_path)


#%% #########################################################################################
'''
freq data, TF vs  NT/TGFb cosine similarity
'''

human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
data = pd.read_excel("siRNA_Analyses/Raw_Data/reformatted_flow_freq.xlsx")
data = data.rename(columns = {'Unnamed: 0':'Target'})

# drop all rows with Target = 'NT' or 'Tgfb'
control_data = data[data['Target'].isin(['NT', 'Tgfb'])]
# delete the first 2 elements in data['target'] to match the format of human_TFs
TF_data, nonTF_data = split_data_by_TF(data, human_TFs)

sim_dict, sim_avg_df = calc_pariwsie_cosine_similarity(TF_data, control_data)   
# save the similarity matrix
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/freq_TF_control_cosine_similarity.csv')

output_path = 'siRNA_Analyses/cosine_similarity/freq_TF_control_cosine_similarity.pdf'
plot_heatmap(sim_avg_df, output_path)


#%% #########################################################################################
'''
MFI data, TF vs non-TF cosine similarity
'''
# read in the TF database
human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
data = pd.read_excel("siRNA_Analyses/Raw_Data/reformatted_flow_MFI.xlsx")
data = data.rename(columns = {'Unnamed: 0':'Target'})

# drop all rows with Target = 'NT' or 'Tgfb'
data = data[(data['Target'] != 'NT') & (data['Target'] != 'Tgfb')]
# delete the first 2 elements in data['target'] to match the format of human_TFs
TF_data, nonTF_data = split_data_by_TF(data, human_TFs)

sim_dict, sim_avg_df = calc_pariwsie_cosine_similarity(TF_data, nonTF_data)   
# save the similarity matrix
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/MFI_TF_nonTF_cosine_similarity.csv')

output_path = 'siRNA_Analyses/cosine_similarity/MFI_TF_nonTF_cosine_similarity.pdf'
plot_heatmap(sim_avg_df, output_path)
#%% ############################################################################################################
'''
MFI data, TF vs  NT/TGFb cosine similarity
'''

human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
data = pd.read_excel("siRNA_Analyses/Raw_Data/reformatted_flow_MFI.xlsx")
data = data.rename(columns = {'Unnamed: 0':'Target'})

# drop all rows with Target = 'NT' or 'Tgfb'
control_data = data[data['Target'].isin(['NT', 'Tgfb'])]
# delete the first 2 elements in data['target'] to match the format of human_TFs
TF_data, nonTF_data = split_data_by_TF(data, human_TFs)

sim_dict, sim_avg_df = calc_pariwsie_cosine_similarity(TF_data, control_data)   
# save the similarity matrix
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/MFI_TF_control_cosine_similarity.csv')


output_path = 'siRNA_Analyses/cosine_similarity/MFI_TF_control_cosine_similarity.pdf'
plot_heatmap(sim_avg_df, output_path)