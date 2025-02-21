import pandas as pd
import celloracle as co
import os
import glob


human_TF = pd.read_csv('/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/allTFs_hg38_Scenic.txt', header=None)

# Read all .txt files 
txt_files = glob.glob('ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.1_0.5_out_final/corr_networks/*.txt')
dataframes = []
file_names = []

for file in txt_files:
    df = pd.read_csv(file, sep='\t')
    dataframes.append(df)
    file_names.append(os.path.basename(file))  # Save only the file name



human_TF = set(human_TF[0])
# Iterate over each DataFrame and corresponding file name
for df, file_name in zip(dataframes, file_names):
    
    unique_genes =pd.concat([df['target'], df['source']]).drop_duplicates()
    node_attr = pd.DataFrame(index=unique_genes)
    node_attr['TF'] = node_attr.index.map(lambda gene: 1 if gene in human_TF else 0)

    prefix = file_name.split('_')[0]
    output_filename = f"ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.1_0.5_out_final/corr_networks/{prefix}_node_attr.txt"
    
    node_attr.to_csv(output_filename, sep='\t')
