import os
import dill
import itertools
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def load_all_pickles(directory_path):
    """A function to load all pickle files, in silico perturbation resutls, in a directory.

    Args:
        directory_path (str): the path to the folder

    Returns:
        _type_: a dictionary with the pickle file names as keys and the loaded pickle objects as values.
    """
    # List all files in the directory
    all_files = os.listdir(directory_path)
    
    # Filter out only the pickle files
    pickle_files = [f for f in all_files if f.endswith('.pkl') or f.endswith('.pickle')]
    
    # Dictionary to hold the loaded pickle objects
    loaded_pickles = {}
    
    for pickle_file in pickle_files:
        # Construct the full file path
        file_path = os.path.join(directory_path, pickle_file)
        
        # Load the pickle file
        with open(file_path, 'rb') as f:
            loaded_pickles[pickle_file] = dill.load(f)
    
    return loaded_pickles

def extract_tranition_prob(oracle):
    """Extract the transition probability matrix from the oracle object and cast into a pandas dataframe.

    Args:
        oracle (oracle): an oracle object

    Returns:
        pandas df: a pandas dataframe with the transition probability matrix.
    """
    transition = oracle.transition_prob
    df_transition_prob = pd.DataFrame(transition, index=oracle.adata.obs.index, columns=oracle.adata.obs.index)
    df_transition_prob = df_transition_prob.rename_axis(index="From", columns="To")
    return df_transition_prob


def create_permutation_dict(df):
    """Subroutine for find_max_k.
       Create a dictionary with all possible permutations of two unique values in a column.
       For example, if the column has values "A" and "B", the dictionary will have keys "A_A", "A_B", "B_A", and "B_B".

    Args:
        df (pandas df): a pandas dataframe with only one column. 

    Returns:
        dict: a dictionary with all possible permutations of two unique values in the column as keys and 0 as values.
    """
    # Extract unique values from the column
    unique_values = df.iloc[:, 0].unique()
    
    # Generate all permutations (with repetition)
    permutations = list(itertools.product(unique_values, repeat=2))
    
    # Create a dictionary with permutations as keys
    permutation_dict = {"_".join(permutation): 0 for permutation in permutations}
    
    return permutation_dict



def find_max_k(oracle, k):
    """For the TF of interest perturbation, find the number of times a regions transitions from one status to another status.
       We define the transition of status by using the max k of the transition matrix.

    Args:
        oracle (_type_): an oracle object
        k (int): using top k values to find the most repeated status

    Returns:
        dict: dictionary with all possible permutations of two unique values in the column as keys and the number of times the transition happened as values.
    """
    # transition = oracle.transition_prob
    # df_transition_prob = pd.DataFrame(transition, index=oracle.adata.obs.index, columns=oracle.adata.obs.index)
    # df_transition_prob = df_transition_prob.rename_axis(index="From", columns="To")
    
    df_transition_prob = extract_tranition_prob(oracle)

    status_df = pd.DataFrame(oracle.adata.obs.Status)
    permutation_dict = create_permutation_dict(status_df)
    
    for row in df_transition_prob.iterrows():
        status_from = oracle.adata.obs.Status[row[0]]
        # the the index of max k values in the row
        top_idx =row[1].nlargest(k).index.tolist()
        # get the group of the top_idx
        group = oracle.adata.obs.Status[top_idx]
        # get the sum of all samples in HF_idx
        HF_sum = row[1][group[group == "HF"].index].sum()
        # get the sum of all samples in Conrol_idx
        Control_sum = row[1][group[group == "Control"].index].sum()
        if HF_sum > Control_sum:
            status_to = "HF"
        else:
            status_to = "Control"
        # get the status of the most repeted status in the top k values
        #status_to = oracle.adata.obs.Status[top_idx].mode().iloc[0]
        key = '_'.join([status_from, status_to])

        permutation_dict[key] += 1
    return permutation_dict

def plot_bars_from_dict(data_dict, TF, save_folder):
    """Plot a bar chart from the permutation_dict for each TF.


    Args:
        data_dict (dict): the dictionary with all possible permutations \
                          of two unique values in the column as keys and the \
                          number of times the transition happened as values.

    """
    # Extract keys and values from the dictionary
    labels = list(data_dict.keys())
    values = list(data_dict.values())
    
    # Create a bar chart
    bar_chart = go.Bar(x=labels, y=values)
    
    # Define the layout
    layout = go.Layout(
        title=TF + ' Transition Count',
        xaxis=dict(title="Keys"),
        yaxis=dict(title="Values"),
    )
    
    # Create a figure with the bar chart and layout
    fig = go.Figure(data=[bar_chart], layout=layout)
    output_filename = os.path.join(save_folder + '/' + TF + '_transition_count.pdf')
    fig.write_image(output_filename, format='pdf')
    # Display the figure
    pio.show(fig)


def calc_all_transition_counts(pickles, save_folder, k = 1):
    """Loop through all pickel files, adata_oracle objects, in the directory\
       and find the transition counts from one state to another. .

    Args:
        pickles (dict): a dictionary with the pickle file names as keys and the loaded pickle objects as values.
        save_folder (str): the path to the folder to save the bar charts.
        k (int, optional):using top k values to find the most repeated status. Defaults to 1.
    """
    # initiate transition count dictionary where the TFs are keys and the values are the transition dictionary
    transition_count = {}
    list_adata_oracle = list(pickles.keys())
    for key in list_adata_oracle:
        print(key)
        oracle = pickles[key].oracle
        # get the transition dictionary for the current oracle object
        count = find_max_k(oracle, k)
        # add the transition dictionary to the main dictionary
        TF = key.split('_')[0]
        plot_bars_from_dict(count, TF, save_folder)
        transition_count[TF] = count
    return transition_count


def plot_transition_on_umap(oracle, concat_sum, save_folder, TF):
    """plot the transition probability on UMAP for HF and Control. \
       The gradient of the color represents the transition probability, summed across all HF or Control regions.

    Args:
        oracle (adata oracle):  a perturbed oracle object
        concat_sum (df): a dataframe contains the summation for each region for HF and Control.
        save_folder (str): a path to the folder to save the plots.
        TF (str): the TF name. used for naming the output files.
    """
    for i in range(concat_sum.shape[1]):
        fig = plt.figure(figsize=[12, 10])
        # Show transition probability
        ax = plt.scatter(oracle.embedding[:, 0], # x
                        oracle.embedding[:, 1], # y
                        c=concat_sum.iloc[:, i], # Show transition probability as color
                        norm=colors.LogNorm(),
                        cmap="Blues", lw=0.7, s=38, alpha=0.55,
                        edgecolor="0.7", rasterized=True)

        plt.colorbar(ax)
        plt.title(f"{TF} Transition probability for {concat_sum.columns[i]}")
        plt.show()
        output_filename = os.path.join(save_folder + '/' + TF + '_transition_probability_sum' + concat_sum.columns[i] + '.pdf')
        fig.savefig(output_filename, format='pdf')

def calc_summation_transition_probability(oracle, save_folder, TF, k = 1):
    """For the TF of interest perturbation, for HF and Control, sum up the transition probability to each region.

    Args:
        oracle (oracle): oracle object from the perturbation
        save_folder (str): the path to the folder to save the plots.
        TF (str): the TF of interest. used for naming the output files.
        k (int, optional): the top k values to add. Defaults to 1.

    Returns:
        df: a dataframe contains the summation for each region for HF and Control.
    """
    
    df_transition_prob = extract_tranition_prob(oracle)
    Status =  oracle.adata.obs['Status']
    HF_sum = df_transition_prob[df_transition_prob.index.isin(Status[Status == "HF"].index)]
    # fill HF_sum to all zeros
    HF_sum = pd.DataFrame(0, index=HF_sum.index, columns=HF_sum.columns)
    # fill Control_sum to all zeros
    Control_sum = df_transition_prob[df_transition_prob.index.isin(Status[Status == "Control"].index)]
    Control_sum = pd.DataFrame(0, index=Control_sum.index, columns=Control_sum.columns)

    for row in df_transition_prob.iterrows():
        group = oracle.adata.obs.Status[row[0]]
        if group == "HF":
            top_idx =row[1].nlargest(k).index.tolist()
            HF_sum.loc[row[0], top_idx] = row[1][top_idx]
        elif group == 'Control':
            top_idx =row[1].nlargest(k).index.tolist()
            Control_sum.loc[row[0], top_idx] = row[1][top_idx]
        else:
            print("Error: group not found")
    # concatenate HF_Sum and Control_sum
    concat_sum = pd.concat([HF_sum.sum(axis = 0), Control_sum.sum(axis = 0)], axis = 1)
    concat_sum.columns = ['HF', 'Control']

    plot_transition_on_umap(oracle, concat_sum, save_folder, TF)

    return concat_sum


def all_pickles_transition_summation(pickles, save_folder, k = 1):
    """Loop through all pickel files, adata_oracle objects, in the directory\
       For each perturbation, plot the transition probability on UMAP for HF and Control. 
       The gradient of the color represents the transition probability, summed across all HF or Control regions. 
       For example, consider region1 and extract top k transition probabilities.\
       Then consider region2, and if the top k have overlap with region1, add the transition probability.\

    Args:
        pickles (dict): output for "load_all_pickles". A dictionary contain all pertrubed adata_oracle objects.
        save_folder (str): the path to the folder to save the plots.
        k (int, optional): max k values to consider. Defaults to 1.

    Returns:
        dict: key is the TF and the value is the concatenated sum of transition probability for HF and Control.
    """
    
    transition_sum = dict()
    list_adata_oracle = list(pickles.keys())
    for key in list_adata_oracle:
        print(key)
        oracle = pickles[key].oracle
        TF = key.split('_')[0]
        concat_sum = calc_summation_transition_probability(oracle, save_folder, TF, k)
        transition_sum[key] = concat_sum
    return transition_sum
