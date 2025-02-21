import sys
sys.path.append("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/COAnalyses")

from adata_oracle import *
from oracle_links import *
from helper_funcs import *
from transition_matrix import *
#%% ##############################################################
'''
linked_TF (relaxed threshold but includes the original linked_TF),
use k = 1
'''
# loop through all pickel files in the directory and find the transition matrix
directory_path = 'Cell_Oracle/Active_Fibro/082924/perturbed_oracle/082924_CO/linked_objects_k_1' 
pickles = load_all_pickles(directory_path)
len(pickles)

save_folder = directory_path
transition_count_dict = calc_all_transition_counts(pickles, save_folder, k = 1)

#%% ##############################################################
'''
overlap_TF
use k = 1
'''
# loop through all pickel files in the directory and find the transition matrix
directory_path = 'Cell_Oracle/Active_Fibro/082924/perturbed_oracle/021324_CO/overlap_objects_k_1'  
pickles = load_all_pickles(directory_path)
len(pickles)

save_folder = directory_path
transition_count_dict = calc_all_transition_counts(pickles, save_folder, k = 1)


#%%
'''
linked_TF (relaxed threshold but includes the original linked_TF),
use K = 5
'''
# loop through all pickel files in the directory and find the transition matrix
directory_path = 'Cell_Oracle/Active_Fibro/082924/perturbed_oracle/082924_CO/linked_objects_k_1' 
pickles = load_all_pickles(directory_path)
len(pickles)

save_folder = "Cell_Oracle/Active_Fibro/082924/perturbed_oracle/082924_CO/linked_objects_k_5"
transition_count_dict = calc_all_transition_counts(pickles, save_folder, k = 5)

#%%
'''
overlap_TF
use k = 5
'''
# loop through all pickel files in the directory and find the transition matrix
directory_path = 'Cell_Oracle/Active_Fibro/082924/perturbed_oracle/021324_CO/overlap_objects_k_1'  
pickles = load_all_pickles(directory_path)
len(pickles)

save_folder = "Cell_Oracle/Active_Fibro/082924/perturbed_oracle/021324_CO/overlap_objects_k_5"
transition_count_dict = calc_all_transition_counts(pickles, save_folder, k = 5)


#%%
'''
linked_TF
plot the transition matrix in the format of UMAP
'''


directory_path = 'Cell_Oracle/Active_Fibro/082924/perturbed_oracle/082924_CO/linked_objects_k_1'  
pickles = load_all_pickles(directory_path)

save_folder = "Cell_Oracle/Active_Fibro/082924/perturbed_oracle/082924_CO/toy"
transition_count_dict = all_pickles_transition_summation(pickles, save_folder, k = 58)


#%%


