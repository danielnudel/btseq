import os
import datetime
import time
import shutil
from sys import stderr

from utils import log
import pandas as pd

cwd = os.getcwd()

# Here I am using logs system with 4 levels,
# from the most informative and verbose (1) to the least informative(4)
ts = time.time()
log(
    datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'),
    4
)
log(f"Working in directory {cwd}", 2)

def exit_and_clean(code=-2):
    log("Exiting and removing temporary files and directories.", 1)
    if os.path.exists(os.path.join(cwd, 'tmp_dir')):
        log("Removing the directory tmp_dir", 1)
        shutil.rmtree(os.path.join(cwd, 'tmp_dir'))
    exit(code)


try:
    ss = pd.read_csv(os.path.join(cwd, 'sample_sheet.txt'), names=["Gene", "Sample", "Barcode", "Seq"], sep='\t')
    log(f"Successfully opened sample_sheet.txt file with {len(ss)} rows.", 1)
except:
    log("Something went wrong when opening sample_sheet.txt file. Exiting with code -2.", 4)
    print("Something wend wrong. Exiting. Check the log.txt file.", file=stderr)
    exit_and_clean()


if not os.path.exists(os.path.join(cwd, 'tmp_dir')):
    log("Creating directory tmp_dir to save working files that will be removed in the end.", 2)
    os.mkdir(os.path.join(cwd, 'tmp_dir'))

sample_counter = 1
sample_ids = []
sample_names_bc_style = []
indices = []
for sample_name, sample_data in ss.groupby('Sample'):
    genes = sample_data['Gene'].values.tolist()
    genes = ['>' + gene for gene in genes]
    sequences = sample_data['Seq'].values.tolist()
    sequences = ['N'*20 + seq + 'N'*20 for seq in sequences]
    if len(genes) != len(sequences):
        log("Number of gene names and number of sequences is not equal. Exiting with code -2", 4)
        print("Something wend wrong. Exiting. Check the log.txt file.", file=stderr)
        exit_and_clean()
    fa = ""
    for gene, seq in zip(genes, sequences):
        fa += gene + '\n' + seq + '\n'
    with open(os.path.join(cwd, 'tmp_dir', f"sample_{sample_counter}.fa"), 'w') as fa_file:
        log(f"Writing fasta file for sample {sample_name} with name sample_{sample_counter}.fa.", 2)
        fa_file.write(fa)
    sample_ids.append(sample_counter)
    if sample_counter > 9:
        sample_names_bc_style.append(f"A0{sample_counter}")
    else:
        sample_names_bc_style.append(f"A00{sample_counter}")
    index = sample_data['Barcode'].unique()
    if len(index) > 1:
        log(f"Sample {sample_name} has more than 1 unique barcode. Exiting.", 4)
        exit_and_clean(-2)
    index_r = ""
    for c in index[0][::-1]:
        if c == 'A':
            index_r += 'T'
        if c == 'C':
            index_r += 'G'
        if c == 'T':
            index_r += 'A'
        if c == 'G':
            index_r += 'C'
    log(f"Created reversed version of index {index[0]} which is {index_r}.", 1)
    indices.append(index_r)
    sample_counter += 1

bc_ss_file_path = os.path.join(cwd, 'tmp_dir', 'bc_ss.csv')
bc_ss = pd.DataFrame([])
bc_ss['Sample_ID'] = sample_ids
bc_ss['Sample_Name'] = sample_names_bc_style
bc_ss['index'] = indices
bc_ss.to_csv(bc_ss_file_path, index=False)

# Read the original content
with open(bc_ss_file_path, "r", newline="") as file:
    content = file.readlines()

# Insert "[DATA]" at the beginning
content.insert(0, "[DATA]\n")

# Write back to the file
with open(bc_ss_file_path, "w", newline="") as file:
    file.writelines(content)

#Create BSTarget_input.txt file
with open(os.path.join(cwd, 'tmp_dir', 'BSTarget_input.txt'), 'w') as f:
    for i in sample_ids:
        path1_name = os.path.join(cwd, 'tmp_dir', 'trim_galore',
                                  f"{sample_names_bc_style[i-1]}_S{i}__R1_001_trimmed.fq.gz")
        path2_name = os.path.join(cwd, 'tmp_dir', f"sample_{i}.fa")
        f.write(f"{i} Sample {path1_name} {path2_name}" + '\n')







