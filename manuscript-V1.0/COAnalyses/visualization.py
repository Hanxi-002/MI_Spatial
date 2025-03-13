import sys
sys.path.append("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/COAnalyses")
import pickle as pkl
from adata_oracle import *
from oracle_links import *
from helper_funcs import *
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

def load_magnitude_csvs(folder_path):
    """load all magnitude difference csvs in a folder into a dictionary of dataframes
    Args:
        folder_path (str): path to the folder containing the csv files
    Returns:
        dict: dictionary with keys as the TF names and values as the corresponding DataFrames.
    """

    dataframes_dict = {}
    # Iterate through all files in the folder
    for file in os.listdir(folder_path):
        if file.endswith('.csv'):  # Process only CSV files
            key = file.split('_')[0]
            df = pd.read_csv(os.path.join(folder_path, file))
            dataframes_dict[key] = df

    
    for key, df in dataframes_dict.items():
        print(f"Key: {key}, DataFrame shape: {df.shape}")
    return dataframes_dict



def label_cluster_with_status(adata, keep_mix = False):
    """
    Label each cluster as HF, Control, or Mix based on the status of the cells in the cluster.

    Args:
        adata (anndata): Anndata object containing the status of each cell.
        keep_mix = False (bool): If True, keep the clusters with a mix of HF and Control cells. If False, exclude the clusters with a mix of HF and Control cells.

    Returns:
        dict: dictionary with keys as the status of the cluster and values as the cluster numbers.
    """
    if keep_mix == True:
        cluster_dict = {'HF': [], 'Control': [], 'Mix': []}
        for cluster in adata.obs['leiden'].unique():
            subset = adata[adata.obs['leiden'] == cluster]
            if (subset.obs['Status'] == 'HF').all():
                cluster_dict['HF'].append(cluster)
            elif (subset.obs['Status']== 'Control').all():
                cluster_dict['Control'].append(cluster)
            else:
                cluster_dict['Mix'].append(cluster)
            cluster_dict['Mix'] = [int(x) for x in cluster_dict['Mix']]

    elif keep_mix == False:
        cluster_dict = {'HF': [], 'Control': []}
        for cluster in adata.obs['leiden'].unique():
            subset = adata[adata.obs['leiden'] == cluster]
            if (subset.obs['Status'] == 'HF').all():
                cluster_dict['HF'].append(cluster)
            elif (subset.obs['Status']== 'Control').all():
                cluster_dict['Control'].append(cluster)

    cluster_dict['HF'] = [int(x) for x in cluster_dict['HF']]
    cluster_dict['Control'] = [int(x) for x in cluster_dict['Control']]
    return cluster_dict

def merge_dictionaries(dict1, dict2):
    """
    Merges two dictionaries into a single dictionary while ensuring that overlapping keys have identical values.

    Args:
        dict1 (dict): The first dictionary to merge.
        dict2 (dict): The second dictionary to merge.

    Returns:
        dict: A new dictionary that is the union of `dict1` and `dict2`.

    Raises:
        AssertionError: If overlapping keys have different values in `dict1` and `dict2`.
    """
    # Ensure that for overlapping keys, the values are the same
    for key in dict1.keys() & dict2.keys():  # Intersection of keys
        print(f"'{key}' is a TF that is both overlapping and linked to SLIDE LFs")
        value1, value2 = dict1[key], dict2[key]
        if isinstance(value1, pd.DataFrame) and isinstance(value2, pd.DataFrame):
            assert value1.equals(value2), f"Value mismatch for key '{key}': DataFrames are not equal"
        else:
            assert value1 == value2, f"Value mismatch for key '{key}': {value1} != {value2}"
    
    # Combine dictionaries
    merged_dict = {**dict1, **dict2}  # dict2 values take precedence if no conflicts
    return merged_dict

def calcualte_cluster_frac_and_diff_per_TF(magnitude_diff, status_dict, magnitude_thresh, TF_name):
# first check if there are any datapoints that has a magnitude difference greater than a certain threshold
    if max(magnitude_diff['magnitude']) > magnitude_thresh:
        # if so, we want to plot this TF. Initialize the plot_dict with required information
        TF_plot_df = pd.DataFrame(columns=['Cluster', 'TF', 'Fraction', 'Median Differences'])
        
        TF_plot_df['Cluster'] = magnitude_diff['cluster_label'].unique()
        TF_plot_df['TF'] = [TF_name] * len(magnitude_diff['cluster_label'].unique())
        
        
        for cluster in magnitude_diff['cluster_label'].unique():
            # extract cluster specific data
            cluster_data = magnitude_diff[magnitude_diff['cluster_label'] == cluster]
            # if this cluster has datapoints with magnitude difference greater than the threshold
            if max(cluster_data['magnitude'] > magnitude_thresh):
                # get the fraction
                cluster_fraction = sum(cluster_data['magnitude'] > magnitude_thresh) / len(cluster_data)
                # get the median differences
                med_diff = np.median(cluster_dat[acluster_data['magnitude'] > magnitude_thresh]['magnitude'])
            else:
                cluster_fraction = 0
                med_diff = 0
            
            TF_plot_df.loc[TF_plot_df['Cluster'] == cluster, 'Fraction'] = cluster_fraction
            TF_plot_df.loc[TF_plot_df['Cluster'] == cluster, 'Median Differences'] = med_diff
        return TF_plot_df      
    else:
        print('no datapoints have a magnitude difference greater than the threshold')
        return None 
    

def calcualte_status_frac_and_diff_per_TF(magnitude_diff, status_dict, magnitude_thresh, TF_name):
# first check if there are any datapoints that has a magnitude difference greater than a certain threshold
    if max(magnitude_diff['magnitude']) > magnitude_thresh:
        magnitude_diff['Status'] = magnitude_diff['cluster_label'].apply(lambda x: 'HF' if x in status_dict['HF'] else 'Control' if x in status_dict['Control'] else 'Mix')
        # if so, we want to plot this TF. Initialize the plot_dict with required information
        TF_plot_df = pd.DataFrame(columns=['Status', 'TF', 'Fraction', 'Median Differences'])
    
        TF_plot_df['Status'] = status_dict.keys()
        TF_plot_df['TF'] = [TF_name] * len(TF_plot_df)
        
        #for cluster in magnitude_diff['cluster_label'].unique():
        for status in status_dict.keys():  
            # extract cluster specific data
            status_data = magnitude_diff[magnitude_diff['cluster_label'].isin(status_dict[status])]
            # if this cluster has datapoints with magnitude difference greater than the threshold
            if max(status_data['magnitude'] > magnitude_thresh):
                # get the fraction
                status_fraction = sum(status_data['magnitude'] > magnitude_thresh) / len(status_data)
                # get the median differences
                med_diff = np.median(status_data[status_data['magnitude'] > magnitude_thresh]['magnitude'])
            else:
                status_fraction = 0
                med_diff = 0
            
            TF_plot_df.loc[TF_plot_df['Status'] == status, 'Fraction'] = status_fraction
            TF_plot_df.loc[TF_plot_df['Status'] == status, 'Median Differences'] = med_diff
        return TF_plot_df      
    else:
        print('no datapoints have a magnitude difference greater than the threshold')
        return None 


def plot_bubble_plot(df, title, output_file):
    # Define fixed bubble size values
    MIN_BUBBLE_SIZE = 0
    MAX_BUBBLE_SIZE = 70  # Fixed bubble size across runs
    
    # Validate that "Fraction" is normalized between [0, 1]
    assert df["Fraction"].max() <= 1 and df["Fraction"].min() >= 0, \
        "Fraction column values must be between 0 and 1."

    # Map "Fraction" directly to bubble sizes
    bubble_sizes = df["Fraction"] * MAX_BUBBLE_SIZE  # Fixed scale mapping

    # Create the scatter plot without automatic size scaling
    fig = px.scatter(
        df,
        x="TF",
        y="Status",
        color="Median Differences",
        color_continuous_scale="Viridis",  # Color scale
        title=title,
        labels={
            "Fraction": "Fraction of Cells",
            "Median Differences": "Median Differences"
        }
    )

    # Update bubble sizes manually
    fig.update_traces(marker=dict(size=bubble_sizes, sizemode='diameter'))

    # Add x-axis and y-axis lines and expand plot dimensions
    fig.update_layout(
        xaxis_title="TF",
        yaxis_title="Status",
        template="plotly_white",
        yaxis=dict(dtick=1),  # Show every integer on y-axis
        shapes=[
            # x-axis line
            dict(
                type="line",
                x0=-0.5,
                y0=-0.5,
                x1=len(df["TF"].unique()) - 0.5,
                y1=-0.5,
                line=dict(color="black", width=2)
            ),
            # y-axis line
            dict(
                type="line",
                x0=-0.5,
                y0=-0.5,
                x1=-0.5,
                y1=len(df["Status"].unique()) - 0.5,
                line=dict(color="black", width=2)
            )
        ],
    )

    # Shrink the colorbar legend
    fig.update_traces(
        marker=dict(
            colorbar=dict(
                thickness=15,  # Adjust thickness (width) of the colorbar
                len=0.5,  # Adjust length of the colorbar (0.0 to 1.0)
                title=dict(
                    text="Median Differences",
                    side="right"  # Position of the colorbar title
                )
            )
        )
    )

    # Save and show the figure
    fig.write_image(output_file, format='pdf')
    fig.show()
# def plot_bubble_plot(df, title, output_file):

#     # Create the bubble plot
#     fig = px.scatter(
#         df,
#         x="TF",
#         y="Status",
#         size="Fraction",
#         color="Median Differences",
#         color_continuous_scale="Viridis",  # Change color scale as needed
#         title=title,
#         labels={
#             "Fraction": "Fraction of Cells",
#             "Median Differences": "Median Differences"
#         },
#         size_max=25,  # Shrinks max bubble size
#     )

#     # Add x-axis and y-axis lines and expand plot dimensions
#     fig.update_layout(
#         xaxis_title="TF",
#         yaxis_title="Status",
#         template="plotly_white",
#         yaxis=dict(dtick=1),  # Show every integer on y-axis
#         shapes=[
#             # x-axis line
#             dict(
#                 type="line",
#                 x0=-0.5,
#                 y0=-0.5,
#                 x1=len(df["TF"].unique()) - 0.5,
#                 y1=-0.5,
#                 line=dict(color="black", width=2)
#             ),
#             # y-axis line
#             dict(
#                 type="line",
#                 x0=-0.5,
#                 y0=-0.5,
#                 x1=-0.5,
#                 y1=len(df["Status"].unique()) - 0.5,
#                 line=dict(color="black", width=2)
#             )
#         ],
#     )
#     # Shrink the colorbar legend
#     fig.update_traces(
#         marker=dict(
#             colorbar=dict(
#                 thickness=15,  # Adjust thickness (width) of the colorbar
#                 len=0.5,  # Adjust length of the colorbar (0.0 to 1.0)
#                 title=dict(
#                     text="Median Differences",
#                     side="right"  # Position of the colorbar title
#                 )
#             )
#         )
#     )

#     fig.write_image(output_file, format='pdf')
#     fig.show()


def plot_bubble_size_legend(size_values, max_size, title, output_file):
    """
    Plot a standalone bubble size legend with percentages for bubble sizes.

    Args:
        size_values (list): List of numeric values representing bubble sizes as decimals.
        max_size (float): Maximum bubble size to scale the bubbles.
        title (str): Title of the bubble size legend.
        output_file (str): Path to save the output image file.
    """
    # Scale the bubble sizes
    bubble_sizes = [v * max_size for v in size_values]
    
    # Create dummy x and y positions to align the bubbles
    x_positions = [1] * len(size_values)  # All bubbles on the same vertical line
    y_positions = list(range(len(size_values), 0, -1))  # Space them vertically

    # Create the figure
    fig = go.Figure()

    # Add bubbles as scatter points
    fig.add_trace(
        go.Scatter(
            x=x_positions,
            y=y_positions,
            mode="markers+text",
            marker=dict(
                size=bubble_sizes,
                color="lightgray",
                line=dict(color="black", width=1)
            ),
            # Convert size values to percentage strings (e.g., 50%)
            text=[f"{v * 100:.0f}%" for v in size_values],
            textposition="middle right",
            showlegend=False
        )
    )

    # Update layout to make it clean and focused
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center"),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        height=300,  # Adjust height based on the number of bubbles
        width=200,   # Narrow width for a compact legend
        margin=dict(l=20, r=20, t=50, b=20),
        plot_bgcolor="white"
    )

    # Save and show the figure
    fig.write_image(output_file, format='pdf')
    fig.show()


#%% Eample code to run the functions
# provide the threshold for magnitude difference
# magnitude_thresh = 0.25
# all_plot_dfs = []

# ################################# load the adata object #################################
# adata = pkl.load(open("Cell_Oracle/Active_Fibro/021324/active_fibroblast_adata_oracle_co_020124.pkl", "rb"))
# adata = adata.adata
# # plot the umap with leiden cluster colors and status colors
# sc.pl.umap(adata, color=['leiden', 'louvain', 'Status'])
# status_dict = label_cluster_with_status(adata)

# ################################# load the adata object #################################
# folder_path = "Cell_Oracle/Active_Fibro/021324/overlap_figures"
# dataframe_dict= load_magnitude_csvs(folder_path)


# all_plot_dfs = []
# for TF_name, magnitude_diff in dataframe_dict.items():
#     TF_plot_df = calcualte_status_frac_and_diff_per_TF(magnitude_diff, status_dict, magnitude_thresh, TF_name)
#     all_plot_dfs.append(TF_plot_df)

# df = pd.concat(all_plot_dfs)
# # change the datatype of df to float
# df['Fraction'] = df['Fraction'].astype(float)
# df['Median Differences'] = df['Median Differences'].astype(float)


# title = 'AF_overlap_TFs'
# output_file = 'Cell_Oracle/Active_Fibro/021324/overlap_figures/AF_overlap_TFs.pdf'
# plot_bubble_plot(df, title, output_file)

