# mrequant_to_dicom
Convert MREquant binary files to DICOM images

Requires python3.   install prerequistes

python3 -m venv venv
venv/Scripts/activate
(venv) python3 -m pip install -r requirements.txt

To run, provide input data folder.

python3 quant_to_dicom.py <data_folder>

Theory of operation:

Will look for 'quant' subfolder of <input_folder>, confirm that mrequant and mmdi files are present, and then create mask and snrmask dicom files from the mrequant data.