# MI_Spatial\spatial network analysis\_hotnet2 \ 051424

## Notes:
* Run hotnet2 network propagation on HF samples. Run each cell type separatly and use the median of the quantile normalized counts as gene score. 


## scripts
* run_hotnet.sh: shell file for sbatch to run the job
* get_gene_score.R: get cell type specific gene score form the GeoMx Data
* analyze_hotnet2_res.py: calls the helper functions to analyze hotnet2 output.
* network_viz.py: combine the SLIDE results with the hotnet2 results (add node attributes for cytoscape visualization)

## result files:
* CD_68.txt: gene score of CD68 regions.
* resting_fibroblast.txt: gene score of resting_fibroblast regions.
* active_fibroblast.txt: gene score of active_fibroblast regions.
* CD68(AF; RF)_hotnet_out: folders containing output of the hotnet2 run. 

## reults:
* CD68 far from AF: no significant modules.