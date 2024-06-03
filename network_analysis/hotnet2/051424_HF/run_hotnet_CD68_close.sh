#!/bin/bash

#SBATCH -t 5-00:00 # Runtime in D-HH:MM

#SBATCH --job-name=CD68_close

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --ntasks-per-node=4

#SBATCH --mem=200g

#SBATCH --mail-user=xiaoh@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL


module load python/anaconda2.7-4.2.0


# Set location of HotNet2, number of cores, number of network permutations, and number of heat permutations.
#hotnet2=..
num_cores=-1

python /ix/djishnu/TFH/hotnet2-master/HotNet2.py \
	-nf  /ix/djishnu/Priyamvada/Immport/Network_analysis/Networks/Networks_to_run_hotnet/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5 \
	-pnp /ix/djishnu/Priyamvada/Immport/Network_analysis/Networks/Networks_to_run_hotnet/permuted_networks/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4_##NUM##.h5 \
	-hf  /ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/051424_HF/CD_68_close.txt \
	-np 500 \
	-hp 500 \
	-o  /ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/051424_HF/CD68_close_hotnet_out \
	-c  $num_cores
