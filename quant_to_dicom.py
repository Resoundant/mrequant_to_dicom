import numpy as np
import matplotlib.pyplot as plt
import pydicom
import pydicom.uid
import os
import re
import sys
from copy import deepcopy

def quant_to_dicom(directory_path):
    '''
    quant_to_dicom.py
    Convert MRE quant binary data to DICOM.
    input: directory to process: must contain DICOM and quant files.
    Output: ROIs as dicom

    Theory: Will open one DICOM files to get patient
    and image metadata.  Then will open the quant file
    and stream pixel information into new dicoms.
    '''
    directory_path = os.path.abspath(directory_path)
    print("")
    print(f"Converting MREquant to DICOM in {directory_path}")
    quant_dir = find_quant_directory(directory_path)
    if quant_dir == '':
        print('  ERROR: could not find a valid quant directory (missing dicoms or quant files)')
        return
    print(f"  Found quant directory {quant_dir}")
    # load mmdi-quant dicoms
    all_dicoms = load_dicoms(quant_dir)
    print(f'  Found {len(all_dicoms)} dicoms')
    stiff_dicoms = extract_stiffness_dicoms(all_dicoms)
    if len(stiff_dicoms) == 0:
        print('  ERROR: No stiffness dicoms found.')
        return
    print(f'  Found {len(stiff_dicoms)} stiffness dicoms')
    stiff_dicoms = sort_stiffness_dicoms(stiff_dicoms)
    shape = get_image_shape(stiff_dicoms)

    # find mask files
    mask_filepath, snrmask_filepath = get_mask_filepaths(quant_dir)
    if not check_data_is_related(stiff_dicoms[0], mask_filepath, snrmask_filepath):
        print('  ERROR: Mask files do not appear to be related to DICOM files')
        return
    print(f"  found mask files {os.path.basename(mask_filepath)} and {os.path.basename(snrmask_filepath)}")

    # mask dicom output folder
    # output_folder = os.path.join(quant_dir, 'test_masks')
    output_folder = quant_dir
    # open masks and turn them into DICOMs
    masks = open_mask(mask_filepath, shape)
    mask_dicoms = masks_to_dicom(masks, stiff_dicoms, desc="mask")
    save_dicoms(mask_dicoms, output_folder, desc="mask")
    print(f"  saved {len(mask_dicoms)} mask dicoms to {output_folder}")
    #open snrmasks and save as dicoms
    snrmasks = open_mask(snrmask_filepath, shape)
    snrmask_dicoms = masks_to_dicom(snrmasks, stiff_dicoms, desc="snrmask")
    save_dicoms(snrmask_dicoms, output_folder, desc="snrmask")
    print(f"  saved {len(snrmask_dicoms)} snrmask dicoms to {output_folder}")
    return

def find_quant_directory(directory_path) -> str:
    ''' 
    Return the FIRST directory that contains the mmdi DICOM outputs and MRE quant files, or empty string if not found.
    Specifically, it will look for the existance of one .mask file and one .snrmask (quant files)
    It will also look for MMDI's files by its default file extensions (.sdcopen and .quant)
    This may seem very resctrictive, but we are specifically looking for the "quant" folder
    to exclude possible confusion to other data folders.
    '''
    for root, dirs, files in os.walk(directory_path):
        # skip paths not ending with 'quant'
        if os.path.basename(root) != 'quant':
            continue
        masks_found = 0
        snrmasks_found = 0
        quants_found = 0
        sdcopens_found = 0
        for f in files:
            if f.endswith('.mask'):
                masks_found += 1
            if f.endswith('.snrmask'):
                snrmasks_found += 1
            if f.endswith('.quant'):
                quants_found += 1
            if f.endswith('.sdcopen'):
                sdcopens_found += 1
        if (masks_found == 1) and (snrmasks_found == 1) and (quants_found > 1) and (sdcopens_found > 1):
            return root
    return ''
            

def load_dicoms(filepath) -> list[pydicom.Dataset]:
    ds_list = []
    # assume a flat folder...
    files = [os.path.join(filepath, f)  for f in os.listdir(filepath) if os.path.isfile(os.path.join(filepath, f))]
    for file in files:
        try:
            ds = pydicom.dcmread(file)
            insert_stiffness_based_on_filename(ds, os.path.basename(file))
            ds_list.append(ds)
        except:
            continue
    return ds_list

def insert_stiffness_based_on_filename(ds:pydicom.Dataset, filename:str) -> None:
    if not ds.get('ImageComments', '') and ('s00' in filename):
        #image comments missing, but file is stiffness filename.
        ds.ImageComments = 'stiffness'
    return

def extract_stiffness_dicoms(ds_list:list[pydicom.Dataset]):
    stiff_list = []
    for ds in ds_list:
        if (ds.get('ImageComments', '').strip().lower() == "stiffness") and ds.get('PhotometricInterpretation', '') == "MONOCHROME2":
            stiff_list.append(ds)
    return stiff_list

def sort_stiffness_dicoms(stiff_list:list[pydicom.Dataset]):
    return sorted(stiff_list, key=lambda x: x.SliceLocation)

def get_image_shape(ds_list:list[pydicom.Dataset]) -> tuple[int,int,int]:
    """
    Returns the shape of the 3D image given the list of 2D slices.
    
    Parameters
    ----------
    ds_list : list[pydicom.Dataset]
        A list of 2D slices in the 3D image.
    
    Returns
    -------
    tuple[int,int,int]
        The shape of the 3D image.
    """
    num_slices = len(ds_list)
    return (ds_list[0].Columns, ds_list[0].Rows,num_slices)

def get_mask_filepaths(quant_dir:str) -> tuple[str, str]:
    # return first .mask and .snrmask file 
    mask_filepath = os.path.join(quant_dir, [f for f in os.listdir(quant_dir) if f.endswith('.mask')][0])
    snrmask_filepath = os.path.join(quant_dir, [f for f in os.listdir(quant_dir) if f.endswith('.snrmask')][0])
    return mask_filepath, snrmask_filepath

def check_data_is_related(ds1:pydicom.Dataset, mask_filepath, snrmask_filepath) -> bool:
    mask_study_id = os.path.basename(os.path.splitext(mask_filepath)[0])
    mask_study_id = re.search(r'\d+$', mask_study_id).group()
    snrmask_study_id = os.path.basename(os.path.splitext(snrmask_filepath)[0])
    snrmask_study_id = re.search(r'\d+$', snrmask_study_id).group()
    study_id = ds1.get('StudyID')
    if (mask_study_id == study_id) and (mask_study_id == study_id) and (snrmask_study_id == study_id):
        return True
    return False

def masks_to_dicom(masks:np.ndarray, ds_list:list[pydicom.Dataset], desc='') -> list[pydicom.Dataset]:
    base_series_num = base_series_number(ds_list)
    roi_series_num = int(str(base_series_num) + '97')
    if desc == '':
        desc = 'mrequantfile'
        roi_series_num = int(str(base_series_num) + '97')
    if desc == 'mask':
        roi_series_num = int(str(base_series_num) + '98')
    if desc == 'snrmask':
        roi_series_num = int(str(base_series_num) + '99')
    
    roi_ds_list = []
    roi_series_uid = pydicom.uid.generate_uid()
    for i, ds in enumerate(ds_list):
        roi_ds = deepcopy(ds)
        roi_ds.PixelData = masks[:,:,i].tobytes()
        roi_ds.BitsAllocated = 8
        roi_ds.BitsStored = 8
        roi_ds.HighBit = 7
        roi_ds.PixelRepresentation = 0  # 0 for unsigned integer data, common for 8-bit
        # roi_ds.PixelData.VR = "OB"
        roi_ds.SmallestImagePixelValue = 0
        roi_ds.LargestImagePixelValue = 1
        roi_ds.WindowCenter = 1
        roi_ds.WindowWidth = 2
        roi_ds.SeriesDescription = desc
        roi_ds.ImageComments = "from mrequant"
        roi_ds.SeriesNumber = roi_series_num
        roi_ds.SeriesInstanceUID = roi_series_uid
        roi_ds.SOPInstanceUID = pydicom.uid.generate_uid()

        roi_ds_list.append(roi_ds)
    return roi_ds_list

def base_series_number(ds_list:list[pydicom.Dataset]) -> int:
    series_numbers = [ds.get('SeriesNumber') for ds in ds_list]
    series_numbers = [x for x in series_numbers if x]     # remove all None 
    if len(series_numbers) == 0:
        print('  ERROR: could not find a valid base series number in dicoms')
        return 88888
    series_numbers = list(set(series_numbers))
    if len(series_numbers) > 1:
        print('  ERROR: more than one base series number in dicoms, using the lowest')
    return series_numbers[0]

def save_dicoms(ds_list:list[pydicom.Dataset], mypath:str, desc = ''):
    if desc ==  '':
        desc = 'mrequantfile'
    os.makedirs(mypath, exist_ok=True)
    for i, ds in enumerate(ds_list):
        ds.save_as(os.path.join(mypath, f"{desc}_{i}.dcm"))

def plot_mask(image):
    plt.imshow(image, cmap="gray", vmin=0, vmax=1)
    plt.title("Raw image")
    plt.colorbar()
    plt.show()

def open_mask(mask_filepath, shape:tuple[int,int,int])->np.ndarray:
    ''' shape is (height, width, slices)'''
    with open(mask_filepath, 'rb') as f:
        mask = np.fromfile(f, dtype=np.uint8)
    image = mask.reshape(shape, order='F') * 255 # order Fortran style
    for i in range(image.shape[2]):
        image[:,:,i] = np.transpose(image[:,:,i]) # transpose to match original quant origin
    return image #reshape doesn't create new

def display_help():
    print("Usage: python quant_to_dicom.py <path_to_quant_directory>")

if __name__ == "__main__":
    datapath = sys.argv[1]
    if (len(sys.argv) != 2):
        display_help()
        exit(1)
    if  (not os.path.isdir(datapath)):
        display_help()
        print(f"ERROR: {datapath} is not a directory")
        exit(1)
    quant_to_dicom(datapath)
