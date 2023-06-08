import pandas as pd
import numpy as np
import rpy2

annot = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/raw_data/Dcc_Initial_Dataset_Dutta02.xlsx")

sum_count = annot.groupby('SegmentLabel').sum()

y = pd.DataFrame(annot["DccNames"])
y['Labels'] = np.full([len(y)], np.nan)
#%%
for s in annot["SlideName"].unique():
    df_slide = annot[annot["SlideName"] == s]
    for r in df_slide["ROILabel"].unique():
        df_roi = df_slide[df_slide["ROILabel"] == r]
        if len(df_roi) == 3:
            resting_ratio = df_roi[df_roi["SegmentLabel"] == "resting fibroblast"]["AOINucleiCount"].values[0] / sum_count.loc["resting fibroblast"]['AOINucleiCount']
            active_ratio = df_roi[df_roi["SegmentLabel"] == "active fibroblast"]["AOINucleiCount"].values[0] / sum_count.loc["active fibroblast"]['AOINucleiCount']
            dcc = df_roi[df_roi["SegmentLabel"] == "CD 68"]["DccNames"].values[0]
            print(resting_ratio)
            print(active_ratio)
            
            if resting_ratio > active_ratio:
                y.loc[y["DccNames"] == dcc, "Labels"] = 0
            elif resting_ratio < active_ratio:    
                y.loc[y["DccNames"] == dcc, "Labels"] = 1

y = y.dropna()
y.to_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/052223/Data/y.csv", index = False)