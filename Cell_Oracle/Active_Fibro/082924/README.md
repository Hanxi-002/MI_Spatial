# AF Cell Oracle 082924
## Notes
Rerun Cell Oracle to 
    1. get more relaxed linked TF. 
    2. Save the oracle object with transition matrix. 

## Input
* **active_fibro_norm.csv:** quantil normalized GeoMX measurements for all AF regions.
* **active_fibro.csv:** raw GeoMX measurements for all AF regions.
* **Final_Annotation_Dutta_sample_info.xlsx:** metadata in the excel format. 
* **021324/active_fibroblast_adata_oracle_co_020124.pkl:** adata_oracle object from the previous analyses. Using this object so we can inherit all the pre-process fields such as UMAP coordinates. 
* **021324/active_fibroblast_oracle_links_co_020124.pkl:** oracle_links object from the previous analyses. 

## Scripts
* **AF_082924_CO.py:** code for running in-sillico perturbations and for relaxed linked_TF.
* **calc_transition_mat.py:** code for analyzing the transition matrix and transition counts after the in silico analyes. Using max to identify transitions. 

## Outputs
* **pertrubed_oracle/021324_CO/overlap_objects_k_1:** the folder contains the perturbed adata_oracle object from the original link and adata_oracle objects produced in 021324 analyses, calculating transition counts with k set to 1. 
* **pertrubed_oracle/082924_CO/linked_objects_k_5:** used perturbed oracle from "021324_CO/overlap_objects_k_1". The folder also contains all the bar plots for the transition counts, with k set to 5. 
* **pertrubed_oracle/082924_CO/linked_objects_k_1:** the folder contains perturbed adata_oracle objects from the linked_TFs after relaxing the threshold of source node degrees. The folder also contains all the bar plots for the transition counts, with k set to 1. 
* **pertrubed_oracle/082924_CO/linked_objects_k_1:** used perturbed oracle from "082924_CO/linked_objects_k_1". The folder also contains all the bar plots for the transition counts, with k set to 5. 

* **figres/new_linked_TFs:** contains figures of the linked_TFs after relaxing the source node threshold. This is a super-set of the original linked_TFs in 021324.