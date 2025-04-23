import sys
sys.path.append("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/COAnalyses")
import pickle as pkl
from adata_oracle import *
from oracle_links import *
from helper_funcs import *
from visualization import *
import plotly.express as px
import pandas as pd

%matplotlib inline

# provide the threshold for magnitude difference
magnitude_thresh = 0.25
all_plot_dfs = []

################################# load the adata object #################################
adata = pkl.load(open("Cell_Oracle/CD68/021424/adata_oracle_CD68_CO_020124.pkl", "rb"))
adata = adata.adata
# plot the umap with leiden cluster colors and status colors
sc.pl.umap(adata, color=['leiden', 'louvain', 'Status'])
status_dict = label_cluster_with_status(adata)

###################################################################################################
# for overlap TFs
################################# load the cell oracle results #################################
folder_path_overlap = "Cell_Oracle/CD68/021424/overlap_figures"
dataframe_dict_overlap = load_magnitude_csvs(folder_path_overlap)

folder_path_linked = "Cell_Oracle/CD68/021424/linked_figures"
dataframe_dict_linked = load_magnitude_csvs(folder_path_linked)

# take the union of the two dictionaries
dataframe_dict = merge_dictionaries(dataframe_dict_overlap, dataframe_dict_linked)


################################# calculation and plotting #################################
all_plot_dfs = []
for TF_name, magnitude_diff in dataframe_dict.items():
    TF_plot_df = calcualte_status_frac_and_diff_per_TF(magnitude_diff, status_dict, magnitude_thresh, TF_name)
    all_plot_dfs.append(TF_plot_df)

df = pd.concat(all_plot_dfs)
# change the datatype of df to float
df['Fraction'] = df['Fraction'].astype(float)
df['Median Differences'] = df['Median Differences'].astype(float)

title = 'CD68_Dist_Bubble'
output_file = 'Cell_Oracle/CD68/021424/' + title + '_V1.pdf'
plot_bubble_plot(df, title, output_file, plot_bubble_size=True)
output_file = 'Cell_Oracle/CD68/021424/size_legend.pdf'
plot_bubble_size_legend([0.05, 0.1, 0.2, 0.3, 0.4, 0.5], 35, "Bubble Size Legend", output_file)

output_file = 'Cell_Oracle/CD68/021424/' + title + '_V2.pdf'
plot_bubble_plot(df, title, output_file, plot_bubble_size=False)


# ###################################################################################################
# # for linked TFs
# ################################# load the cell oracle results #################################
# folder_path = "Cell_Oracle/CD68/021424/linked_figures"
# dataframe_dict= load_magnitude_csvs(folder_path)

# ################################# calculation and plotting #################################
# all_plot_dfs = []
# for TF_name, magnitude_diff in dataframe_dict.items():
#     TF_plot_df = calcualte_status_frac_and_diff_per_TF(magnitude_diff, status_dict, magnitude_thresh, TF_name)
#     all_plot_dfs.append(TF_plot_df)

# df = pd.concat(all_plot_dfs)
# # change the datatype of df to float
# df['Fraction'] = df['Fraction'].astype(float)
# df['Median Differences'] = df['Median Differences'].astype(float)

# title = 'CD68_dist_linked_TFs_Bubble'
# output_file = folder_path + '/' + title + '.pdf'
# plot_bubble_plot(df, title, output_file)