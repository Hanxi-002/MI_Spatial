#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 13:50:07 2022

@author: xiaoh
"""

#%%
# need to modify DSP file names after the core re-did the DSP files. 
#1001660011816-C goes to 1001660011816-D
#1001660012268-D goes to 1001660012268-C

f = open("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Dutta 01_20220721T1618_LabWorksheet V2.txt", "r")
tmp = f.readlines()
f.close()
new_dsp = list()
for r in tmp:
    l = r.split("\t")
    #print(l)
    if "1001660011816" in l[0]:
       assert l[0][18] == 'C' , "shouldn't be here..."
       tmp_list = list(l[0])
       tmp_list[18] = "D"
       l[0] = "".join(tmp_list)
       
    if "1001660012268" in l[0]:
        assert l[0][18] == 'D' , "shouldn't be here v2..."
        tmp_list = list(l[0])
        tmp_list[18] = "C"
        l[0] = "".join(tmp_list)
    new_dsp.append("\t".join(l))


outfile = open("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Dutta_LabWorksheet_HX.txt", "a")
for r in new_dsp:
    outfile.write(r)
outfile.close()



#%%
#old code for before they fixed the multiplex labeling issue, do not run. 
import pandas as pd
pd.options.mode.chained_assignment = None 
annot = pd.read_excel('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Annotation Dutta sample info V2.xlsx', index_col=0) 
#dcc_name =  pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Dutta 01_20220721T1618_LabWorksheet V2.txt", sep ='\t')
dcc_name =  pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Dutta_LabWorksheet_HX.txt", sep ='\t')



#annot = annot.sort_values("Segment (Name/ Label)", ascending = True)
annot['DCCnames'] = ['na'] * len(annot)

#annot = annot.copy()
# deal with each slide seperatly
group = "collect slide 100 PD78"
subset_annot = annot[annot.index == group]
subset_dcc = dcc_name[dcc_name['scan name'] == group]
assert len(subset_annot) == len(subset_dcc)

# match the roi_id style for the 2 files
roi_id = pd.unique(subset_annot["ROI (label)"])
new_id = [0] * len(roi_id)
for i in range(len(roi_id)):
    if len(str(roi_id[i])) == 1:
        each_id = "00"+ str(roi_id[i])
    else:
        each_id = "0"+ str(roi_id[i])
    new_id[i] = each_id
    
# range through each of the roi_id
for j in range(len(new_id)):
    subset = subset_dcc[subset_dcc['roi'].str.contains(new_id[j])]
    for k in range(len(subset)):
        temp_row = subset.iloc[[k]]`
        seg = temp_row['segment'].values[0]
        annot.loc[(annot.index == group) & (annot["ROI (label)"] == roi_id[j]) & (annot['Segment tags'] == seg), 'DCCnames'] = temp_row['Sample_ID'].values[0]
       

#chec that all dcc files are unique
assert len(pd.unique(annot['DCCnames'])) == len(annot)

#annot.to_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Annotation Dutta sample info V4.xlsx")

#%% 
# check the difference between old and new annotations

old = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Initial Dataset Anno.xlsx", index_col = 0)
new = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Initial Dataset Dutta02.xlsx", index_col = 0)

old_check = old.loc[(old["ScanLabel"] == "collect slide 103 PD81") & (old["ROILabel"] == 1)]
new_check = new.loc[(new["ScanLabel"] == "collect slide 103 PD81") & (new["ROILabel"] == 1)]

#old_dict = old.to_dict("records")
toy = pd.concat([old_check,new_check]).drop_duplicates(keep=False)
toy.to_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/check_difference.xlsx")

#%%
# new code for after they fixed the multiplex issue.
import pandas as pd
pd.options.mode.chained_assignment = None 
#import the new annotation file.
annot = pd.read_excel('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Initial Dataset Dutta02.xlsx', index_col=0) 
dcc_name =  pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Dutta 01_20220721T1618_LabWorksheet V2.txt", sep ='\t')

annot['DCCnames'] = ['na'] * len(annot)


#seperate out the dataframe for each group
group = "colklect slide 98 PD77"
#subset_annot = annot[annot.index == group] # if groups is not slide 98 PD77, run thisline. 
subset_annot = annot[annot["ScanLabel"] == group] # if group is slide 98 PD 77, run this line. There is a typo in the annotation file. 
subset_dcc = dcc_name[dcc_name['scan name'] == group]
assert len(subset_annot) == len(subset_dcc)


roi_id = pd.unique(subset_annot["ROILabel"])
new_id = [0] * len(roi_id)
for i in range(len(roi_id)):
    if len(str(roi_id[i])) == 1:
        each_id = "00"+ str(roi_id[i])
    else:
        each_id = "0"+ str(roi_id[i])
    new_id[i] = each_id
    
for j in range(len(new_id)):
    subset = subset_dcc[subset_dcc['roi'].str.contains(new_id[j])]
    for k in range(len(subset)):
        temp_row = subset.iloc[[k]]
        seg = temp_row['segment'].values[0]
        # if group is not slide 98 PD77, run this line. 
        #annot.loc[(annot.index == group) & (annot["ROILabel"] == roi_id[j]) & (annot['SegmentLabel'] == seg), 'DCCnames'] = temp_row['Sample_ID'].values[0]
        # if group is slide 98 PD 77, run this line. 
        annot.loc[(annot["ScanLabel"] == group) & (annot["ROILabel"] == roi_id[j]) & (annot['SegmentLabel'] == seg), 'DCCnames'] = temp_row['Sample_ID'].values[0]
        
#annot.to_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/Annotation Dutta sample info V3.xlsx")
