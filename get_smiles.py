from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import random
import pandas as pd
from rdkit.Chem import SaltRemover, MolStandardize
from huggingface_hub import list_repo_files
from datasets import load_dataset
from rdkit.Chem import rdMolDescriptors
import bittensor as bt
import json
from rdkit import rdBase

# Suppress RDKit warnings
rdBase.DisableLog('rdApp.warning')

chunk_size = 300000
dataset_repo = "Metanova/SAVI-2020"

files = list_repo_files(dataset_repo, repo_type='dataset')
files = [file for file in files if file.endswith('.csv')]
random_file = random.choice(files)

similarity_thres = 0.4

file_num = 0

target_index = 5400758
start_index = target_index % chunk_size
end_index = 5400758 + 992
skip_chunks = target_index // chunk_size + 1
bt.logging.info(f"Get chunk at {skip_chunks}, index from {start_index} to {end_index}")

file = files[file_num]
# file_num = 0
# for file in files:
    # bt.logging.info(f"Start file {file_num}")
    # file_num+=1
dataset_dict = load_dataset(
    dataset_repo,
    data_files={'train': file},
    streaming=True,
)
dataset = dataset_dict['train']
batched = dataset.batch(chunk_size)


sample_smiles = None # Example: Q6P6W3
sample_smiles_list = {}

chunk_num = 0
smile_num = 0
total_fetched = 0
for chunk in batched:
    if chunk_num < skip_chunks:
        bt.logging.info(f"Skip {chunk_num} chunk")
        chunk_num += 1
        continue
    df = pd.DataFrame.from_dict(chunk)
    # Clean data
    df['product_name'] = df['product_name'].apply(lambda x: x.replace('"', ''))
    df['product_smiles'] = df['product_smiles'].apply(lambda x: x.replace('"', ''))
    bt.logging.info(f"Fetched: {len(df['product_smiles'].tolist())}")

    # 2. Your database of molecules to search (this is a tiny example list)
    database_smiles_list = df['product_smiles'].tolist()
    
    with open(f"{file_num}_{target_index}_total_chunk.json", 'w') as f:
        json.dump(database_smiles_list, f)
    with open(f"{file_num}_{target_index}_list_results.json", 'w') as f:
        json.dump(database_smiles_list[start_index:end_index], f)