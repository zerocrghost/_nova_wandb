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
skip_chunks = 0

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
# sample_smiles_chunk = {}
results = []


if skip_chunks > 0:
    with open(f"{file_num}_list_results.json", 'r') as f:
        list_results = json.load(f)
        sample_smiles_list = list_results
    # with open(f"{file_num}_chunk_results.json", 'r') as f:
    #     chunk_results = json.load(f)
    #     sample_smiles_chunk = chunk_results

chunk_num = 0
smile_num = 0
total_fetched = 0
for chunk in batched:
    if chunk_num <= skip_chunks:
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
    total_fetched += len(database_smiles_list)
    if sample_smiles == None:
        sample_smiles = database_smiles_list[0]
        sample_smiles_list[sample_smiles] = 0
        # sample_smiles_chunk[sample_smiles] = chunk_num
    # 1. Create a molecule object from SMILES
    mol = Chem.MolFromSmiles(sample_smiles)
    # bt.logging.info(f"Sample: {sample_smiles}")

    # 2. Standardize: Remove salts, neutralize charges, canonicalize tautomers
    remover = SaltRemover.SaltRemover()
    mol = remover.StripMol(mol)  # Remove salts like HCl, Na

    # (Optional but recommended) Use MolStandardize for more robust cleaning
    # normalizer = MolStandardize.normalize.Normalizer()
    # mol = normalizer.normalize(mol)

    # 3. Generate a canonical SMILES string to use as your query
    canonical_smiles_query = Chem.MolToSmiles(mol)
    # bt.logging.info(f"Standardized Query SMILES: {canonical_smiles_query}")

    # 1. Your standardized query molecule
    query_smiles = canonical_smiles_query  # From Step 1
    query_mol = Chem.MolFromSmiles(query_smiles)
    # Generate Morgan Fingerprint for the query (radius=2 is equivalent to ECFP4)
    # query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)
    query_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)

    # 4. Iterate through the database and calculate similarity
    for smi in database_smiles_list:
        db_mol = Chem.MolFromSmiles(smi)
        if db_mol is None: # Skip invalid SMILES
            continue
        db_fp = AllChem.GetMorganFingerprintAsBitVect(db_mol, radius=2, nBits=2048)
        # Calculate Tanimoto similarity
        similarity = DataStructs.TanimotoSimilarity(query_fp, db_fp)
        if similarity >= similarity_thres:
            # Append the result
            results.append({
                'SAMPLE': sample_smiles,
                'SMILES': smi,
                'Similarity': similarity,
                'SMILE_NUM': smile_num
            })
            sample_smiles_list[sample_smiles] += 1
 
        else:
            # bt.logging.info(f"Update smilarity_sample")
            sample_smiles = smi
            sample_smiles_list[sample_smiles] = 0
            # bt.logging.info(f"Sample: {sample_smiles}")
            
            # 1. Create a molecule object from SMILES
            mol = Chem.MolFromSmiles(sample_smiles)

            # 2. Standardize: Remove salts, neutralize charges, canonicalize tautomers
            remover = SaltRemover.SaltRemover()
            mol = remover.StripMol(mol)  # Remove salts like HCl, Na

            # (Optional but recommended) Use MolStandardize for more robust cleaning
            # normalizer = MolStandardize.normalize.Normalizer()
            # mol = normalizer.normalize(mol)

            # 3. Generate a canonical SMILES string to use as your query
            canonical_smiles_query = Chem.MolToSmiles(mol)
            # bt.logging.info(f"Standardized Query SMILES: {canonical_smiles_query}")

            # 1. Your standardized query molecule
            query_smiles = canonical_smiles_query  # From Step 1
            query_mol = Chem.MolFromSmiles(query_smiles)
            # Generate Morgan Fingerprint for the query (radius=2 is equivalent to ECFP4)
            # query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)
            query_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)
            similarity = DataStructs.TanimotoSimilarity(query_fp, db_fp)

            # Append the result
            results.append({
                'SAMPLE': sample_smiles,
                'SMILES': smi,
                'Similarity': similarity,
                'SMILE_NUM': smile_num
            })
            sample_smiles_list[sample_smiles] += 1
        
        smile_num += 1    
            # sample_smiles_chunk[sample_smiles] = chunk_num

    # # 5. Create a DataFrame and sort by similarity (descending order)
    # # results_df = pd.DataFrame(results)
    # # results_df.sort_values(by='Similarity', ascending=False, inplace=True)
    # # print(results_df.tolist())
    # results.sort(key=lambda x: x["Similarity"], reverse=True)
    # temp_top_tens = top_tens + results
    # temp_top_tens.sort(key=lambda x: x["Similarity"], reverse=True)
    # top_tens = temp_top_tens[0:11]
    # bt.logging.info(f"Round {round} top tens list")
    # for top in top_tens:
    #     bt.logging.info(f"{top}")
    # round += 1
    # bt.logging.info(f"sample_smiles_list: {sample_smiles_list}")
    chunk_num += 1
    bt.logging.info(f"Fetching round: {chunk_num}, total: {total_fetched}")
# file_path = os.path.join(os.path.dirname(__file__), file.name)
    with open(f"{file_num}_list_results.json", 'w') as f:
        json.dump(sample_smiles_list, f)
    # with open(f"{file_num}_chunk_results.json", 'w') as f:
    #     json.dump(sample_smiles_chunk, f)