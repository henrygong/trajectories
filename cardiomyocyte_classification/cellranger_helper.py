import subprocess
import json
import socket

with open('scRNAseq_5.04.18_180509_K00161_0217_AHVFFLBBXX.json') as f:
    cellranger_command = json.load(f)

host_name = socket.gethostname()
# Working
if host_name == 'tap.soe.ucsc.edu':
    print(cellranger_command['EF10W4D_LV'])
    subprocess.call(cellranger_command['EF10W4D_LV'], shell=True)

# Working
if host_name  == 'citrisdance.soe.ucsc.edu':
    print(cellranger_command['EF10W4D_RA'])
    subprocess.call(cellranger_command['EF10W4D_RA'], shell=True)

# Working
if host_name == 'riverdance.soe.ucsc.edu':
    print(cellranger_command['EF10W4D_RV'])
    subprocess.call(cellranger_command['EF10W4D_RV'], shell=True)

# Working
if host_name == 'waterdance.soe.ucsc.edu':
    print(cellranger_command['EF10W4D_LA'])
    subprocess.call(cellranger_command['EF10W4D_LA'], shell=True)
