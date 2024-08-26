"""
Process the qPCR data and the siRNA data
Train a svm model for each pair of siRNA targets
"""
import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelEncoder
from itertools import combinations


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
    for i in range(100): 
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
########## the code below is not what we want ###############

# apply the trained model on each of the siRNA targets

# def predict_target_prob(target, trained_model):
#     """Calculate the probability of the target being a positive sample\
#        based on the trained model.

#     Args:
#         target (_str_): target name
#         trained_model (_model_): a trained model

#     Returns:
#         _ndarray_: a numpy array of shape: # samples by # classes\
#                    each row sums up to be 1. 
#     """
#     target_data = data.loc[data['Target'].str.contains(target)]
#     target_data = target_data.drop('Target', axis=1)
#     prob = trained_model.predict_proba(target_data)
#     return prob


# # train the svm model
# svm_model = SVC(probability=True)
# svm_model.fit(X, y) # train on all nt and tgfb data

# # predict the probability of each target and save in a dictionary
# all_target_prob = {}
# for target in data['Target'].unique():
#     if target == 'NT' or target == 'Tgfb':
#         continue
#     print(f"Target: {target}")
#     target_prob = predict_target_prob(target, svm_model)
#     all_target_prob[target] = target_prob

# all_target_prob["siDdx5"]

 