"""
Process the qPCR data and the siRNA data
Train a svm model for each pair of siRNA targets
ENV:MI_Spatial
"""
import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelEncoder
import plotly.graph_objects as go


#%% #############################################
# train pairwise svm model for auc comparison
def svm_cross_val(X, y, svm_model, n_splits):
    """get the performance of the svm through cross validation

    Args:
        X (_df_): the input data frame with only the 2 classes of samples
        y (_df_): transfomred y, should be binary
        svm_model (_model_): the svm object (untrained)
        n_splits (_int_): number of splits to do 

    Returns:
        _list_: a list of aucs, length equals to n_splits
    """
    cv = StratifiedKFold(n_splits, shuffle=True, random_state=42)
    aucs = []
    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y)):
        #print(f"Fold {fold_idx + 1}:")
        #print(f"  Training indices: {train_idx}")
        #print(f"  Testing indices: {test_idx}")
        svm_model.fit(X.iloc[list(train_idx)], y[train_idx])
        y_hat = svm_model.predict_proba(X.iloc[list(test_idx)])[:,1]
        auc = roc_auc_score(y[test_idx], y_hat)
        #y_hat = svm_model.predict(X.iloc[list(test_idx)])
        #score = svm_model.score(X.iloc[list(test_idx)], y.iloc[list(test_idx)])
        #print(f"  y_true: {y[test_idx]}")
        #print(f"  y_hat: {y_hat}")
        #print(f"  AUC: {auc}")
        aucs.append(auc)
    return aucs


def train_svm(data, target_1, target_2, n_splits=3):
    """train svm model for 2 targets

    Args:
        data (_df_): the input qPCR data
        target_1 (_str_): name of siRNA target 1
        target_2 (_str_): name of siRNA target 1
        n_splits (int, optional): the number of folds to split. Defaults to 3.

    Returns:
        _list_: the list of aucs for the input targets.
    """
    # clean up the input data
    #print(f"Training SVM for {target_1} and {target_2}")
    X = data.loc[data['Target'].isin([target_1, target_2])]
    y = X['Target']
    X = X.drop('Target', axis=1)
    label_encoder = LabelEncoder()
    y = label_encoder.fit_transform(y)

    # train the model through cross validation
    svm_model = SVC(probability=True)
    aucs = svm_cross_val(X, y, svm_model, n_splits)

    #print(f"Mean AUC: {np.mean(aucs)}")
    return aucs 


def main(num_iter, data, pairwise_targets):
    auc_all_targets = {}
    for i in range(num_iter): 
        print(f"The {i}th round")
        # for each target pair, perform k-fold cross validation
        for p in pairwise_targets:
            #print(p)
            target_1, target_2 = p
            aucs = train_svm(data, target_1, target_2, 4)
            # if auc_all_targets dictionary has no target_1 and target_2 as key
            if f"{target_1}_{target_2}" not in auc_all_targets:
                auc_all_targets[f"{target_1}_{target_2}"] = [np.mean(aucs)]
            else:
                # append auc to the existing list the of the dictionary value
                l = auc_all_targets[f"{target_1}_{target_2}"]
                l.append(np.mean(aucs))
                auc_all_targets[f"{target_1}_{target_2}"] = l
    return auc_all_targets

#%% #############################################
# train the model for the qPCR results

data = pd.read_excel("/ix/djishnu/Hanxi/MI_Spatial/siRNA_Analyses/Raw_Data/reformatted_qPCR.xlsx")

# change the name of the first column of data
data = data.rename(columns = {'Unnamed: 0':'Target'})

targets = data['Target'].unique()
pairwise_targets = [[t, "NT"] for t in targets if t != 'NT']

auc_all_targets = main(100, data, pairwise_targets)

# convert the dictionary to a dataframe
auc_all_targets =  pd.DataFrame(auc_all_targets)
# add a row to the dataframe that contains the mean of the auc for each target pair
auc_all_targets.loc['mean'] = auc_all_targets.mean()
# put row mean as the 1st row
auc_all_targets = auc_all_targets.reindex(index=['mean'] + list(auc_all_targets.index[:-1]))

#auc_all_targets.to_csv("/ix/djishnu/Hanxi/MI_Spatial/siRNA_Analyses/Pairwise_Classifier/auc_all_targets.csv")


#%% #############################################
# train the model for the flow frequency results

data = pd.read_excel("/ix/djishnu/Hanxi/MI_Spatial/siRNA_Analyses/Raw_Data/reformatted_flow_freq.xlsx")


targets = data['Target'].unique()
pairwise_targets = [[t, "NT"] for t in targets if t != 'NT']

auc_all_targets = main(100, data, pairwise_targets)


# convert the dictionary to a dataframe
auc_all_targets =  pd.DataFrame(auc_all_targets)
# add a row to the dataframe that contains the mean of the auc for each target pair
auc_all_targets.loc['mean'] = auc_all_targets.mean()
# put row mean as the 1st row
auc_all_targets = auc_all_targets.reindex(index=['mean'] + list(auc_all_targets.index[:-1]))

#auc_all_targets.to_csv("/ix/djishnu/Hanxi/MI_Spatial/siRNA_Analyses/Pairwise_Classifier/freq_auc_all_targets.csv")
#%% #############################################
# train the model for the flow MFI results
data = pd.read_excel("/ix/djishnu/Hanxi/MI_Spatial/siRNA_Analyses/Raw_Data/reformatted_flow_MFI.xlsx")


targets = data['Target'].unique()
pairwise_targets = [[t, "NT"] for t in targets if t != 'NT']

auc_all_targets = main(100, data, pairwise_targets)


# convert the dictionary to a dataframe
auc_all_targets =  pd.DataFrame(auc_all_targets)
# add a row to the dataframe that contains the mean of the auc for each target pair
auc_all_targets.loc['mean'] = auc_all_targets.mean()
# put row mean as the 1st row
auc_all_targets = auc_all_targets.reindex(index=['mean'] + list(auc_all_targets.index[:-1]))

#auc_all_targets.to_csv("/ix/djishnu/Hanxi/MI_Spatial/siRNA_Analyses/Pairwise_Classifier/MFI_auc_all_targets.csv")

#%% #############################################
# plot the aucs
auc_qPCR = pd.read_csv("siRNA_Analyses/Pairwise_Classifier/qPCR_auc_all_targets.csv", index_col=0)
print(auc_qPCR.shape)
auc_freq = pd.read_csv("siRNA_Analyses/Pairwise_Classifier/freq_auc_all_targets.csv", index_col=0)
print(auc_freq.shape)
auc_MFI = pd.read_csv("siRNA_Analyses/Pairwise_Classifier/MFI_auc_all_targets.csv", index_col=0)
print(auc_MFI.shape)


# merge the mean rows of the three dataframes
auc_qPCR = auc_qPCR.loc[auc_qPCR.index == 'mean']
auc_freq = auc_freq.loc[auc_freq.index == 'mean']
auc_MFI = auc_MFI.loc[auc_MFI.index == 'mean']
# merge the three dataframes
plot_df = pd.concat([auc_qPCR, auc_freq, auc_MFI], axis=0)
plot_df.index = ['qPCR', 'freq', 'MFI']
plot_df.shape
# drop the columns that have NaN
# plot_df = plot_df.dropna(axis=1)
# plot_df.shape

# split the plot_df into 2 parts: TF and NonTF
human_TFs = pd.read_csv('Cell_Oracle/allTFs_hg38_Scenic.txt', sep='\t', header=None)
human_TF_set = set(human_TFs[0].str.upper())
gene_names = ['Ddx5', 'B3galnt1', 'Hnrnph3', 'Npipb6', 'Tdrd6', 'Tvp23b', 'Aasdhppt', 'Tgfb', 'Maz', 'Mafk', 'Atf3', 'Egr1', 'Ep300', 'Klf3']
TF_names = [f"si{g.upper()[0]}{g.lower()[1:]}_NT" for g in gene_names if g.upper() in human_TF_set]
# add in the controls
TF_names_all = TF_names + ['siCtrl oxLDL_NT', 'siCtrl nLDL_NT', 'Tgfb_NT']
plot_df_TF = plot_df[TF_names_all]
plot_df_nonTF = plot_df.drop(columns=TF_names)

# if plot_df_TF has value smaller than 0.5, than set it as 1-value
plot_df_TF = plot_df_TF.applymap(lambda x: 1-x if x < 0.5 else x)
plot_df_nonTF = plot_df_nonTF.applymap(lambda x: 1-x if x < 0.5 else x)
plot_df_nonTF['siIqgap1_NT']  = [1.0, 0.5825, 0.9550]
# drop the column siIqgap1_NT from plot_df_TF
plot_df_nonTF = plot_df_nonTF.drop(columns=['siIqgap_NT'])

# plot the data

def plot_grouped_bars(df, out_path):
    #df = plot_df_nonTF
    fig = go.Figure()

    # Define colors for each metric
    colors = {'qPCR': '#8f9162', 'freq': '#c77875', 'MFI': '#8cc8d0'}

    # Track position for each group
    group_positions = {}
    current_position = 0
    bar_height = 0.2
    group_gap = 0.3

    # Get all column names (gene names)
    gene_names = sorted(df.columns.tolist())

    # For each gene, create a group of bars
    for i, gene in enumerate(gene_names):
        group_positions[gene] = current_position
        
        # Add a bar for each metric if it's not NaN
        for j, metric in enumerate(df.index):
            value = df.loc[metric, gene]
            
            # Skip if value is NaN
            if pd.notna(value):
                fig.add_trace(go.Bar(
                    y=[current_position + j*bar_height],  # Position for this bar
                    x=[value],
                    orientation='h',
                    name=metric,
                    legendgroup=metric,
                    showlegend=True if i == 0 else False,  # Show legend only once per metric
                    marker_color=colors[metric],
                    width=bar_height,
                    hovertemplate=f"{gene} - {metric}: %{{x:.2f}}<extra></extra>"
                ))
        
        # Increase position for next group
        current_position += len(df.index)*bar_height + group_gap

    # Set axis properties
    fig.update_layout(
        title="AUCs ",
        xaxis=dict(
            title="AUC",
            range=[0.5, 1.0],  
            tickformat=".1f",
            showline=True,
            linewidth=1,
            linecolor='black',
        ),
        yaxis=dict(
            tickmode='array',
            tickvals=[pos + bar_height for pos in group_positions.values()],
            ticktext=list(group_positions.keys()),
            showgrid=False,
            showline=True,
            linewidth=1,
            linecolor='black',
        ),
        barmode='group',
        height=500,
        width=500,
        plot_bgcolor='white',
        legend_title="Metrics"
    )

    # Show the plot
    fig.show()

    fig.write_image(out_path, width=500, height=500, scale=2) 

plot_grouped_bars(plot_df_TF, "siRNA_Analyses/Pairwise_Classifier/auc_plot_TF.pdf")
plot_grouped_bars(plot_df_nonTF, "siRNA_Analyses/Pairwise_Classifier/auc_plot_nonTF.pdf")