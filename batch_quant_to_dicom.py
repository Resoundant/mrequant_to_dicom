import os
import sys
from quant_to_dicom import quant_to_dicom

def find_all_quant_dirs(top_dir):
    top_dir = os.path.abspath(top_dir)
    quant_dirs = []
    for root, _, _ in os.walk(top_dir):
        if os.path.basename(root) == 'quant':
            quant_dirs.append(os.path.join(root))
    return quant_dirs   

def save_quantdir_list(quant_dirs, top_dir):
    filepath = os.path.join(top_dir, 'quantdirs_found.txt')
    print(f"Saving quant directories to {filepath}")
    with open(filepath, 'w') as f:
        for quant_dir in quant_dirs:
            f.write(f"{quant_dir}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python batch_quant_to_dicom.py <path_to_directory>")
        print("Will find all quant directories and run quant_to_dicom on them")
        exit(1)
    top_dir = sys.argv[1]
    quant_dirs = find_all_quant_dirs(top_dir)
    print(f"Found {len(quant_dirs)} quant directories to process in {top_dir} ")
    for dir in quant_dirs:
        quant_to_dicom(dir)