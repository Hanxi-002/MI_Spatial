import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np 

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
    plt.savefig(f"{save_folder}/{goi}_hist.pdf")

def calc_cell_identity_shifts (oracle, goi, save_folder, n_neighbors = 50, scale = 50, perturb_value = 0.0):
    oracle.simulate_shift(perturb_condition={goi: perturb_value},
                        n_propagation=3)

    # Get transition probability
    oracle.estimate_transition_prob(n_neighbors=n_neighbors,
                                    knn_random=True,
                                    sampled_fraction=1)
    # Calculate embedding
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

    scale = scale
    # Show quiver plot
    oracle.plot_quiver(scale=scale, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")
    # Show quiver plot that was calculated with randomized graph.
    oracle.plot_quiver_random(scale=scale, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")
    plt.savefig(f"{save_folder}/{goi}_sim_shift.pdf")
    plt.show()

def get_optimal_min_mass(n_grid, oracle, n_neighbors = 50):
    oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=n_neighbors)
    oracle.suggest_mass_thresholds(n_suggestion=12)

def calc_vector_fields(min_mass, goi, oracle, save_folder, scale_simulation = 10):
    #min_mass = 0.14
    oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
    fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
    #plt.savefig(f"{save_folder}/{goi}_grid_points.pdf")
    # Show quiver plot
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

    # Show quiver plot that was calculated with randomized graph.
    oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")
    plt.savefig(f"{save_folder}/{goi}_vector_field.pdf")
    plt.show()

def plot_vector_filed_on_cluster(oracle, goi, save_folder, scale_simulation=10):
    fig, ax = plt.subplots(figsize=[8, 8])
    oracle.plot_cluster_whole(ax=ax, s=10)
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
    ax.set_title(f"Vector filed on clusters: {goi} KO")
    plt.savefig(f"{save_folder}/{goi}_vector_field_on_cluster.pdf")
    plt.show()

def calc_cluster_vector_diff(oracle, adata, goi, save_folder, method):
    # change a flag for dot plot
    # add the method
    # try out the absolute difference
    df_delta = pd.DataFrame(oracle.delta_embedding)
    df_delta.index = adata.obs.index
    df_delta_rand = pd.DataFrame(oracle.delta_embedding_random)
    df_delta_rand.index = adata.obs.index

    diff = pd.DataFrame()
    method = method
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

        # plot box plot the differences of the magnitudes
        cluster_diff = abs(cluster_df_delta['magnitude'] - cluster_df_delta_rand['magnitude'])
        cluster_diff = pd.DataFrame(cluster_diff)
        cluster_diff['cluster_label'] = j
        #diff.append(cluster_diff)
        #add cluster_diff to diff
        diff = pd.concat([diff, cluster_diff])

    # plot the dot plot with cluster_label as x and magnitude as y
    sns.stripplot(x='cluster_label', y='magnitude', data=diff, jitter=0.2, size=5, linewidth=1, edgecolor='gray')
    sns.boxplot(x='cluster_label', y='magnitude', data=diff, whis=[5, 95], showcaps=True, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':2}, medianprops={'color':'red'})

    #plt.scatter(range(len(diff)), diff)
    plt.xlabel('Clusters')
    plt.ylabel('Differences')
    plt.savefig(f"{save_folder}/{goi}_Vector_Magnitude_Difference.pdf")
    # Display the plot
    plt.show()

    # save diff to csv
    diff.to_csv(f"{save_folder}/{goi}_Vector_Magnitude_Difference.csv")
    
