#
import pandas as pd

annot = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/raw_data/Final_Annotation_Dutta_sample_info.xlsx", index_col = 0)
x = pd.read_csv("/Users/xiaoh/Desktop/toy/x.csv", index_col = 0)
y = pd.read_csv("/Users/xiaoh/Desktop/toy/y.csv", index_col = 0)
y = y.drop(columns = "Status")


wrong_name = []
for i in range(len(annot["DCCnames"])):
    name = annot["DCCnames"][i]+".dcc"
    # get the row where index is name in x
    x_row = x.loc[name]
    y_row = y.loc[name].values[0]
    y_annot = annot.iloc[i, :]
    assert y_annot["DCCnames"] == annot["DCCnames"][i], print("DCCnames are not the same")
    y_annot = y_annot["Status"]
    if y_annot == "HF":
        if y_row != 1:
            print(name)
            wrong_name.append(name)
    elif y_annot == "Control":
        if y_row != 0:
            print(name)
            wrong_name.append(name)




#%% NEW CELL (Check between 2 excel files)
file1 = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/raw_data/Final_Annotation_Dutta_sample_info.xlsx", index_col = 0)
file2 = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/raw_data/Dcc_Initial_Dataset_Dutta02.xlsx", index_col=0)
file1 = file1.drop(columns = ["QCFlags"])
file2 = file2.drop(columns = ["QCFlags"])

#compare file1 and file2
#iterate through each row in file1
for i in range(len(file1)):
    for j in range(len(file1.iloc[i, :].values)):
        if file1.iloc[i, :].values[j] != file2.iloc[i, :].values[j]:
            print(file1.index[i])
            print(file1.columns[j])
        

