# Cell-Cell Communication in MI

## Notes
Running Commot seperately on HF and Control.

## Input Files
    * GeoMxV3: the Geomx V3 cleaned up version. 
    * LR database: the folder that contains the ligand receptor database from cellchat. 
    * ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv: the y that encodes proximity of Macs

## Scripts
    * get_data.R: extract X matrix from GeoMxV3
    * make_adata_obj.py: make adata object and run commot.
    

## Output Files
    * norm_x.txt: the quantiled normalized matrix from GeoMxV3 object. 
    * raw_x.txt: the raaw matrix from GeoMxV3 matrix. 
    * adata_commot_hf.pkl: saved adata object after 
    running spatial_communication
    * adata_commot_control.pkl: saved adata object after running spatial_communication