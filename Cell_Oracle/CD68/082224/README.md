# CD68 Cell Oracle 082224
## Notes
Rerun Cell Oracle to 
    1. get better UMAP visualization. Results: the visualization looks similar to old version. So continue with older adata-oracle object. 
    2. Get more realxed GRN. 
    3. get more relaxed linked TF. 
    4. Save the oracle object with transition matrix. 

## Input
* **CD68_norm.csv:** quantil normalized GeoMX measurements for all CD68 regions.
* **CD68.csv:** raw GeoMX measurements for all CD68 regions.
* **Final_Annotation_Dutta_sample_info.xlsx:** metadata in the excel format. 
* **020124/adata_oracle_CD68_CO_020124.pkl:** adata_oracle object from the previous analyses. Using this object so we can inherit all the pre-process fields such as UMAP coordinates. 
* **020124/CD68_oracle_links.pkl:** oracle_links object from the previous analyses. 

## Scripts
* **umap_Viz.py:** re-run UMAP visualization. 
* **CD68_082224_CO.py:** code for running in-sillico perturbations.

## Outputs
* **pertrubed_oracle_objects/020124_CO:** the folder contains the perturbed adata_oracle object from the original link and adata_oracle objects produced in 020124 analyses. 
* **pertrubed_oracle_objects/082224_CO/linked_objects:** the folder contains perturbed adata_oracle objects from the linked_TFs after relaxing the threshold of source node degrees. 
* **figres/new_linked_TFs:** contains figures of the linked_TFs after relaxing the source node threshold. This is a super-set of the original linked_TFs in 020124.