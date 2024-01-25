import pandas as pd
import numpy as np


annot = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/raw_data/Final_Annotation_Dutta_sample_info.xlsx")

sum_count = annot.groupby('SegmentLabel').sum()

y = pd.DataFrame(annot["DCCnames"])
y['Labels'] = np.full([len(y)], np.nan)
#%%
for s in annot["SlideName"].unique():
    df_slide = annot[annot["SlideName"] == s]
    for r in df_slide["ROILabel"].unique():
        df_roi = df_slide[df_slide["ROILabel"] == r]
        if len(df_roi) == 3:
            resting_ratio = df_roi[df_roi["SegmentLabel"] == "resting fibroblast"]["AOINucleiCount"].values[0] / sum_count.loc["resting fibroblast"]['AOINucleiCount']
            active_ratio = df_roi[df_roi["SegmentLabel"] == "active fibroblast"]["AOINucleiCount"].values[0] / sum_count.loc["active fibroblast"]['AOINucleiCount']
            dcc = df_roi[df_roi["SegmentLabel"] == "CD 68"]["DCCnames"].values[0]
            print(resting_ratio)
            print(active_ratio)
            
            if resting_ratio > active_ratio:
                y.loc[y["DCCnames"] == dcc, "Labels"] = 0
            elif resting_ratio < active_ratio:    
                y.loc[y["DCCnames"] == dcc, "Labels"] = 1

y = y.dropna()
#y.to_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/052223/Data/y.csv", index = False)
y.to_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/y.csv", index = False)
#%% Get the annotation file including the distance labels
# produce Y for RF and AF concatenated analyses
annot = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/raw_data/Final_Annotation_Dutta_sample_info.xlsx")
sum_count = annot.groupby('SegmentLabel').sum()

# delete the rows where Dccnames is  "DSP-1001660011816-D-B11"
annot = annot[annot["DCCnames"] != "DSP-1001660011816-D-B11"]
# delete the rows where Dccnames is  "DSP-1001660012268-C-H09"
annot = annot[annot["DCCnames"] != "DSP-1001660012268-C-H09"]


annot_hf = annot[annot["Status"] == "HF"]
annot_hf['AFRF_Ratio'] = np.full([len(annot_hf)], np.nan)
# load 
for s in annot_hf["SlideName"].unique():
    df_slide = annot_hf[annot_hf["SlideName"] == s]
    for r in df_slide["ROILabel"].unique():
        df_roi = df_slide[df_slide["ROILabel"] == r]
        if len(df_roi) == 3:
            resting_ratio = df_roi[df_roi["SegmentLabel"] == "resting fibroblast"]["AOINucleiCount"].values[0] / sum_count.loc["resting fibroblast"]['AOINucleiCount']
            active_ratio = df_roi[df_roi["SegmentLabel"] == "active fibroblast"]["AOINucleiCount"].values[0] / sum_count.loc["active fibroblast"]['AOINucleiCount']
            #dcc = df_roi[df_roi["SegmentLabel"] == "CD 68"]["DCCnames"].values[0]
            print(resting_ratio)
            print(active_ratio)
            
            if resting_ratio > active_ratio:
                label = [0] * 3
                #y.loc[y["DCCnames"] == dcc, "Labels"] = 0
            elif resting_ratio < active_ratio:    
                label = [1] * 3
                #y.loc[y["DCCnames"] == dcc, "Labels"] = 1
            annot_hf.loc[(annot_hf["SlideName"] == s) & (annot_hf["ROILabel"] == r), "AFRF_Ratio"]= label

# drop the rows where AFRF_Ratio is nan
annot_hf = annot_hf.dropna(subset = ["AFRF_Ratio"])
annot_hf.to_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/HF_Final_Annotation_RFAF_Ratio.csv", index = False)







