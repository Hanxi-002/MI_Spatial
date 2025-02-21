import pandas as pd
import numpy as np

# spacia needs 2 files: counts.txt and spacia_metadata.txt
# counts
norm_x = pd.read_csv('/ix/djishnu/Hanxi/MI_Spatial/CCC/norm_x.txt', sep = '\t', index_col=0)



# meta data
# read in the excel sheet
meta = pd.read_excel('/ix/djishnu/Hanxi/MI_Spatial/RawData/Final_Annotation_Dutta_sample_info.xlsx')
spacia_meta = meta[['ROICoordinateX', 'ROICoordinateY', 'SegmentLabel', 'DCCnames']]
spacia_meta.set_index('DCCnames', inplace = True)
spacia_meta.columns = ['X', 'Y', 'cell_type']
spacia_meta['cell_type'] = spacia_meta['cell_type'].str.replace(' ', '_')
spacia_meta.to_csv('/ix/djishnu/Hanxi/MI_Spatial/CCC/spacia_meta.txt', sep = '\t')



#%% 
'''
reedit the spacia meta data. Previously, each AOI has the same XY coordinates if they are from the same ROI. 
However, spacia cannot deal with this issue. Hence, adding a small number on each AOI to make them unique.
'''

# spacia_meta = pd.read_csv('CCC/Spacia/spacia_meta.txt', sep = '\t', index_col = 0)

# # for the X column, add a random number drawn from a normal distribution with mean 0 and std 0.5
# spacia_meta['X'] = spacia_meta['X'] + np.random.normal(0, 0.5, spacia_meta.shape[0])
# # do the same for Y
# spacia_meta['Y'] = spacia_meta['Y'] + np.random.normal(0, 0.5, spacia_meta.shape[0])

# spacia_meta.to_csv('/ix/djishnu/Hanxi/MI_Spatial/CCC/new_spacia_meta.txt', sep = '\t')


spacia_meta = pd.read_csv('CCC/Spacia/spacia_meta.txt', sep = '\t', index_col = 0) 

# if it's CD68, add 1 to each X and Y
# if it's active_fibroblast, add 1 to each X and Y
spacia_meta.loc[spacia_meta['cell_type'] == 'CD_68', ['X', 'Y']] += 1
spacia_meta.loc[spacia_meta['cell_type'] == 'active_fibroblast', ['X', 'Y']] -= 1

spacia_meta.to_csv('/ix/djishnu/Hanxi/MI_Spatial/CCC/Spacia/new_spacia_meta.txt', sep = '\t')


count = pd.read_csv('CCC/Spacia/norm_x_log.txt', sep = '\t', index_col = 0) 
count = pd.read_csv('CCC/Spacia_Git/Spacia/test/input/counts.txt', sep = '\t', index_col = 0)
