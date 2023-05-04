import numpy as np
from sklearn.mixture import GaussianMixture
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



metadata  = pd.read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/raw_data/Annotation Dutta sample info V3.xlsx")


def get_celltype_data(metadata, celltype):
    celltype_metadata = metadata[metadata['SegmentLabel'] == celltype]
    cord_data = celltype_metadata[['ROICoordinateX', 'ROICoordinateY']]
    return cord_data, celltype_metadata

def perform_gmm(cord_data):
    gm = GaussianMixture(n_components=3, random_state=0).fit(cord_data)
    labels = gm.predict(cord_data)
    return labels

def plot_gmm(cord_data, celltype_metadata, labels):
    cord_data['labels'] = labels
    celltype_metadata['labels'] = labels
    sns.scatterplot(x="ROICoordinateX", y="ROICoordinateY", hue="labels", data=cord_data)
    return cord_data, celltype_metadata

def save_new_metadata(celltype_metadata, path):
    celltype_metadata.to_excel(path)

cord_data, celltype_metadata = get_celltype_data(metadata, "resting fibroblast")
labels = perform_gmm(cord_data)
cord_data, celltype_metadata = plot_gmm(cord_data, celltype_metadata, labels)
save_new_metadata(celltype_metadata, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/RestingFibro/RestingFibroblast_Cluster_Label.xlsx")




