#MI_Spatial\spatial network analysis\_hotnet2 \ 051424

## Notes:
* Run hotnet2 network propagation on HF samples. Run each cell type separatly and use the median of the quantile normalized counts as gene score. 


## scripts
* run_hotnet.sh: shell file for sbatch to run the job
* get_gene_score.R: 

## result files:
* CD_68.txt: gene score of CD68 regions.
* resting_fibroblast.txt: gene score of resting_fibroblast regions.
* active_fibroblast.txt: gene score of active_fibroblast regions.
* CD68(AF; RF)_hotnet_out: folders containing output of the hotnet2 run. 
