# Cell-Cell Communication in MI

## Notes
Try to get Spacia to run on the data. 
Couldn't get it run. 

## Input Files
    * GeoMxV3: the Geomx V3 cleaned up version. 
    * LR database : the folder that contains the ligand receptor database from cellchat. 

## Scripts
    * get_data.R: extract X matrix from GeoMxV3
    * get_spacia_input.py: produces the 2 required file to run spacia.


## Output Files
    * norm_x.txt: the quantiled normalized matrix from GeoMxV3 object. 
    * raw_x.txt: the raaw matrix from GeoMxV3 matrix. 
    * spacia_meta.txt: one of the required inputs for spacia
    * spacia_norm_x.txt: the counts required as part of the inputs for spacia. 