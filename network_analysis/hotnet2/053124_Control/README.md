#MI_Spatial\spatial network analysis\_hotnet2 \ 053124

## Notes:
* Run hotnet2 network propagation on control samples. Run each cell type separatly and use the median of the quantile normalized counts as gene score. We can compare the results of these propagations with HF only hotnet2 analyses (051424).


## scripts
* run_hotnet.sh: shell file for sbatch to run the job
* get_gene_score.R: get cell type specific gene score form the GeoMx Data
* analyze_hotnet2_res.py: calls helper functions to analyze hotnet2 output.

## result files:
* CD_68_control.txt: gene score of CD68 regions.
* resting_fibroblast_control.txt: gene score of resting_fibroblast regions.
* active_fibroblast_control.txt: gene score of active_fibroblast regions.
* CD68(AF; RF)_hotnet_out: folders containing output of the hotnet2 run. 

## results:
* AF control runs do not have any significant modules (p-value below 0.1).
* 