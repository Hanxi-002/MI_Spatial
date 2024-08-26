import numpy as np
import pandas as pd


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
#%% #########################################################################################
'''
qPCR data, TF vs non-TF cosine similarity
'''
# read in the TF database
human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
data = pd.read_excel("/ix/djishnu/Hanxi/MI_Spatial/siRNA_Analyses/Raw_Data/reformatted_qPCR.xlsx")
data = data.rename(columns = {'Unnamed: 0':'Target'})

# drop all rows with Target = 'NT' or 'Tgfb'
data = data[(data['Target'] != 'NT') & (data['Target'] != 'Tgfb')]
# delete the first 2 elements in data['target'] to match the format of human_TFs
TF_data, nonTF_data = split_data_by_TF(data, human_TFs)

sim_dict, sim_avg_df = calc_pariwsie_cosine_similarity(TF_data, nonTF_data)   
# save the similarity matrix
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/qPCR_TF_nonTF_cosine_similarity.csv')



# visualize the similarity with a heatmap
import seaborn as sns
import matplotlib.pyplot as plt
sim_avg_df = sim_avg_df.astype(float)
sns.heatmap(sim_avg_df, cmap='coolwarm')
# save the plot
plt.savefig('siRNA_Analyses/cosine_similarity/qPCR_TF_nonTF_cosine_similarity.png')
plt.show()

#%% #########################################################################################
'''
qPCR data, TF vs  NT/TGFb cosine similarity
'''

human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
data = pd.read_excel("/ix/djishnu/Hanxi/MI_Spatial/siRNA_Analyses/Raw_Data/reformatted_qPCR.xlsx")
data = data.rename(columns = {'Unnamed: 0':'Target'})

# drop all rows with Target = 'NT' or 'Tgfb'
control_data = data[data['Target'].isin(['NT', 'Tgfb'])]
# delete the first 2 elements in data['target'] to match the format of human_TFs
TF_data, nonTF_data = split_data_by_TF(data, human_TFs)

sim_dict, sim_avg_df = calc_pariwsie_cosine_similarity(TF_data, control_data)   
# save the similarity matrix
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/qPCR_TF_control_cosine_similarity.csv')


# visualize the similarity with a heatmap
import seaborn as sns
import matplotlib.pyplot as plt
sim_avg_df = sim_avg_df.astype(float)
sns.heatmap(sim_avg_df, cmap='coolwarm')
# save the plot
plt.savefig('siRNA_Analyses/cosine_similarity/qPCR_TF_control_cosine_similarity.png')
plt.show()
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


# visualize the similarity with a heatmap
import seaborn as sns
import matplotlib.pyplot as plt
sim_avg_df = sim_avg_df.astype(float)
sns.heatmap(sim_avg_df, cmap='coolwarm')
# save the plot
plt.savefig('siRNA_Analyses/cosine_similarity/freq_TF_nonTF_cosine_similarity.png')
plt.show()

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


# visualize the similarity with a heatmap
import seaborn as sns
import matplotlib.pyplot as plt
sim_avg_df = sim_avg_df.astype(float)
sns.heatmap(sim_avg_df, cmap='coolwarm')
# save the plot
plt.savefig('siRNA_Analyses/cosine_similarity/freq_TF_control_cosine_similarity.png')
plt.show()

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
sim_avg_df.to_csv('siRNA_Analyses/cosine_similarity/MIF_TF_nonTF_cosine_similarity.csv')



# visualize the similarity with a heatmap
import seaborn as sns
import matplotlib.pyplot as plt
sim_avg_df = sim_avg_df.astype(float)
sns.heatmap(sim_avg_df, cmap='coolwarm')
# save the plot
plt.savefig('siRNA_Analyses/cosine_similarity/MFI_TF_nonTF_cosine_similarity.png')
plt.show()

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


# visualize the similarity with a heatmap
import seaborn as sns
import matplotlib.pyplot as plt
sim_avg_df = sim_avg_df.astype(float)
sns.heatmap(sim_avg_df, cmap='coolwarm')
# save the plot
plt.savefig('siRNA_Analyses/cosine_similarity/MFI_TF_control_cosine_similarity.png')
plt.show()