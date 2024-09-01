import os
import dill
import itertools
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio

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
    transition = oracle.transition_prob
    df_transition_prob = pd.DataFrame(transition, index=oracle.adata.obs.index, columns=oracle.adata.obs.index)
    df_transition_prob = df_transition_prob.rename_axis(index="From", columns="To")
    
    status_df = pd.DataFrame(oracle.adata.obs.Status)
    permutation_dict = create_permutation_dict(status_df)
    
    for row in df_transition_prob.iterrows():
        # get status of the current row
        status_from = oracle.adata.obs.Status[row[0]]
        # the the index of max k values in the row
        top_idx =row[1].nlargest(k).index.tolist()
        # get the status of the most repeted status in the top k values
        status_to = oracle.adata.obs.Status[top_idx].mode().iloc[0]

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

