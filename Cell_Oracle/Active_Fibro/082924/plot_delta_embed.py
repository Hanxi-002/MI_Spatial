import sys
sys.path.append("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/COAnalyses")
import dill
import os
from adata_oracle import *
from oracle_links import *
from helper_funcs import *
import plotly.graph_objs as go
import plotly.express as px

# read in the pkl file
with open('/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Active_Fibro/082924/perturbed_oracle/021324_CO/linked_objects/ATF3_adata_oracle.pkl', 'rb') as f:
    data = dill.load(f)

save_folder = "toy"
os.makedirs(save_folder, exist_ok=True)


oracle = data.oracle
adata = data.adata

df_delta = pd.DataFrame(oracle.delta_embedding)
df_delta.index = adata.obs.index
df_delta_rand = pd.DataFrame(oracle.delta_embedding_random)
df_delta_rand.index = adata.obs.index

all_delta = pd.DataFrame()
all_delta_rand = pd.DataFrame()
for j in np.unique(adata.obs[method].astype(int).values):
    cluster_df = adata.obs[adata.obs[method] == str(j)]
    
    # get the rows in df_delta where the index is in cluster_df
    cluster_df_delta = df_delta[df_delta.index.isin(cluster_df.index)]
    cluster_df_delta_rand = df_delta_rand[df_delta_rand.index.isin(cluster_df.index)]

    assert sum(cluster_df_delta.index == cluster_df.index) == cluster_df.shape[0], \
        "index of cluster_df_delta not match with cluster_df"

    # calculate the magnitude of the vectors
    cluster_df_delta['magnitude'] = np.sqrt(cluster_df_delta[0]**2 + cluster_df_delta[1]**2)
    cluster_df_delta_rand['magnitude'] = np.sqrt(cluster_df_delta_rand[0]**2 + cluster_df_delta_rand[1]**2)

    cluster_df_delta['cluster_label'] = j
    cluster_df_delta_rand['cluster_label'] = j

    all_delta = pd.concat([all_delta, cluster_df_delta])
    all_delta_rand = pd.concat([all_delta_rand, cluster_df_delta_rand])

#%% box plot (not differences, just raw magnitues)
all_delta['condition'] = 'Perturbed'
all_delta_rand['condition'] = 'Random'
df_combined = pd.concat([all_delta, all_delta_rand])

fig = go.Figure()

clusters = df_combined['cluster_label'].unique()

# Add pre-treatment box plots
for cluster in clusters:
    fig.add_trace(go.Box(
        y=df_combined[(df_combined['cluster_label'] == cluster) & (df_combined['condition'] == 'Perturbed')]['magnitude'],
        x=[cluster] * len(df_combined[(df_combined['cluster_label'] == cluster) & (df_combined['condition'] == 'Perturbed')]),
        name=f'{cluster} - Perturbed',
        marker_color='blue',
        boxmean='sd',
        legendgroup=f'{cluster}-Pre',
        showlegend=True
    ))

# Add post-treatment box plots
for cluster in clusters:
    fig.add_trace(go.Box(
        y=df_combined[(df_combined['cluster_label'] == cluster) & (df_combined['condition'] == 'Random')]['magnitude'],
        x=[cluster] * len(df_combined[(df_combined['cluster_label'] == cluster) & (df_combined['condition'] == 'Random')]),
        name=f'{cluster} - Random',
        marker_color='red',
        boxmean='sd',
        legendgroup=f'{cluster}-Post',
        showlegend=True
    ))

# Customize the layout
fig.update_layout(
    title='Magnitude by Cluster for Pre and Post Treatment',
    xaxis=dict(title='Cluster'),
    yaxis=dict(title='Magnitude'),
    boxmode='group',  # Group the boxes together for each cluster
    showlegend=True
)

# Show the plot
fig.show()




#%% scatter plot

# Get unique clusters
clusters = all_delta['cluster_label'].unique()
# Generate x positions for each cluster
# Each cluster has a pre- and post- position close to each other
x_positions = np.arange(len(clusters))

# Create figure
fig = go.Figure()

# Adding pre-treatment points
for idx, cluster in enumerate(clusters):
    cluster_data_pre = all_delta[all_delta['cluster_label'] == cluster]
    fig.add_trace(go.Scatter(
        x=[x_positions[idx] - 0.1] * len(cluster_data_pre),  # Slightly left for pre-treatment
        y=cluster_data_pre['magnitude'],
        mode='markers',
        name=f'{cluster} - Perturbed',
        marker=dict(color='blue'),
        showlegend=True
    ))

# Adding post-treatment points
for idx, cluster in enumerate(clusters):
    cluster_data_post = all_delta_rand[all_delta_rand['cluster_label'] == cluster]
    fig.add_trace(go.Scatter(
        x=[x_positions[idx] + 0.1] * len(cluster_data_post),  # Slightly right for post-treatment
        y=cluster_data_post['magnitude'],
        mode='markers',
        name=f'{cluster} - Random',
        marker=dict(color='red'),
        showlegend=True
    ))

# Update x-axis to show cluster names
fig.update_layout(
    title='Magnitude by Cluster for Pre and Post Treatment',
    xaxis=dict(
        tickvals=x_positions,
        ticktext=clusters,
        title='Cluster'
    ),
    yaxis=dict(
        title='Magnitude'
    ),
    showlegend=True
)

# Show the plot
fig.show()