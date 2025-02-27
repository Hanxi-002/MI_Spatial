import pandas as pd

human_LR = pd.read_csv("CCC/LR_database/human_LR.csv", header = None)
human_LR.columns =['Ligand', 'Receptor', 'PathwayNames', 'Method']
z_99 = pd.read_csv("ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/gene_list_Z99.txt", sep = "\t")

z_37 = pd.read_csv("ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/feature_list_Z37.txt", sep = "\t")
z_27 = pd.read_csv("ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/feature_list_Z27.txt", sep = "\t")
z_21 = pd.read_csv("ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/feature_list_Z21.txt", sep = "\t")
z_17 = pd.read_csv("ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/feature_list_Z17.txt", sep = "\t")
z_3 = pd.read_csv("ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/feature_list_Z3.txt", sep = "\t")

z_99_overlap = human_LR[(human_LR['Ligand'].isin(z_99['names'])) | (human_LR['Receptor'].isin(z_99['names']))]
z_37_overlap = human_LR[(human_LR['Ligand'].isin(z_37['names'])) | (human_LR['Receptor'].isin(z_37['names']))]

ccr2_lf = pd.concat([z_37, z_27, z_21, z_17, z_3])

ccr2_overlap = human_LR[(human_LR['Ligand'].isin(ccr2_lf['names'])) | (human_LR['Receptor'].isin(ccr2_lf['names']))]
