import pandas as pd
import os
import re
# read in excel 
df = pd.read_excel('Final_Annotation_Dutta_sample_info.xlsx', sheet_name='Sheet1')

# obtain the columns of interest
output_df = df[['SlideName', 'ROILabel', 'SegmentLabel', 'Status', 'DCCnames']]



def map_segment_label(label):
    if label == "CD 68":
        return "CD68"
    elif label == "resting fibroblast":
        return "RF"
    elif label == "active fibroblast":
        return "AF"
    else:
        return label.replace(" ", "_")

# Create a copy of the DataFrame to avoid the SettingWithCopyWarning
output_df = output_df.copy()

# Assuming ROILabel is already a column in your DataFrame
# If it's not, you would need to add it first (e.g., output_df['ROILabel'] = [1, 2, 3, 4, 5])

# Create the Library_Name column using the specified format
output_df.loc[:, 'Library_Name'] = (
    output_df['SlideName'].str.replace(' ', '_') + 
    '_' + 
    output_df['ROILabel'].astype(str) + 
    '_' + 
    output_df['SegmentLabel'].apply(map_segment_label)
)

def transform_dcc_name(name):
    parts = name.split('-')
    if len(parts) >= 3:  # Ensure there are enough parts to work with
        return parts[0] + '-' + parts[1] + '-' + parts[-1] + '.dcc'
    else:
        return name + '.dcc'  # Fallback if format is different
    

output_df.loc[:, 'DCCnames'] = output_df['DCCnames'].apply(transform_dcc_name)


def extract_slide_id(slide_name):
    parts = slide_name.split()
    return parts[-1]  # Get the last part of the slide name

# Create the Titles column by combining the slide ID, Status, ROILabel, and mapped SegmentLabel
output_df.loc[:, 'Titles'] = (
    output_df['SlideName'].apply(extract_slide_id) + 
    '_' + 
    output_df['Status'] + 
    '_' + 
    output_df['ROILabel'].astype(str) + 
    '_' + 
    output_df['SegmentLabel'].apply(map_segment_label)
)



def get_files_from_nested_folders(root_folder):
    """
    Get all filenames from nested folders starting with a root folder.
    
    Args:
        root_folder (str): Path to the root folder (FASTQ in this case)
        
    Returns:
        dict: Dictionary with subfolder names as keys and lists of files as values
    """
    all_files = {}
    
    # Check if the root folder exists
    if not os.path.exists(root_folder):
        print(f"Error: Root folder '{root_folder}' does not exist.")
        return all_files
    
    # Get all subfolders in the root folder
    subfolders = [f for f in os.listdir(root_folder) if os.path.isdir(os.path.join(root_folder, f))]
    
    # Go into each subfolder and collect file names
    for subfolder in subfolders:
        subfolder_path = os.path.join(root_folder, subfolder)
        files_in_subfolder = [f for f in os.listdir(subfolder_path) if os.path.isfile(os.path.join(subfolder_path, f))]
        all_files[subfolder] = files_in_subfolder
        
    return all_files


fastq_folder = "FASTQ"
files_by_subfolder = get_files_from_nested_folders(fastq_folder)
all_files_flat = []
for subfolder, files in files_by_subfolder.items():
    for file in files:
        all_files_flat.append(file)

def modify_file_names(file_list):
    modified_files = []
    
    for filename in file_list:
        parts = filename.split('-')
        if len(parts) >= 4:
            
            modified_parts = [parts[0], parts[1], parts[3]]
            modified_filename = '-'.join(modified_parts)
            modified_files.append(modified_filename)
        else:
            print('should not enter here')
    
    return modified_files

modified_files = modify_file_names(all_files_flat)

file_pairs = {}
for modified_name in modified_files:
    # Extract the base name without R1/R2
    base_name = re.sub(r'_(S\d+_L\d+_)R[12](_\d+\.fastq\.gz)', r'_\1*\2', modified_name)
    
    # Check if it's R1 or R2
    if '_R1_' in modified_name:
        if base_name not in file_pairs:
            file_pairs[base_name] = {'R1': None, 'R2': None}
        file_pairs[base_name]['R1'] = modified_name
    elif '_R2_' in modified_name:
        if base_name not in file_pairs:
            file_pairs[base_name] = {'R1': None, 'R2': None}
        file_pairs[base_name]['R2'] = modified_name

# Function to match DCCnames to Raw1 and Raw2 files
def match_raw_files(dcc_name, file_pairs):
    # Remove the .dcc extension if present
    base_name = dcc_name.replace('.dcc', '')
    
    # Find matching file pairs
    matching_pairs = {}
    for pattern, pair in file_pairs.items():
        pattern_base = pattern.split('_')[0]  # Get the first part of the pattern (e.g., DSP-1001660011816-G03)
        if base_name in pattern_base:
            matching_pairs[pattern] = pair
    
    # Return the matching R1 and R2 files
    if matching_pairs:
        # Just take the first match if there are multiple
        first_match = list(matching_pairs.values())[0]
        return first_match['R1'], first_match['R2']
    else:
        return None, None

# Add Raw1 and Raw2 columns to the DataFrame
output_df['Raw1'] = None
output_df['Raw2'] = None

for idx, row in output_df.iterrows():
    raw1, raw2 = match_raw_files(row['DCCnames'], file_pairs)
    output_df.loc[idx, 'Raw1'] = raw1
    output_df.loc[idx, 'Raw2'] = raw2

output_df.to_excel('editted_library_names.xlsx', index=False, sheet_name='Sheet1', startrow=0, header=True)