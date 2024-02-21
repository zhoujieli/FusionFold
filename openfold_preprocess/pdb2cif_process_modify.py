import os
import argparse
from Bio.PDB import PDBParser, MMCIFIO
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from tqdm import tqdm
from multiprocessing import Pool

# 忽略 PDB 构建警告
warnings.simplefilter('ignore', PDBConstructionWarning)

def add_virtual_info_to_cif(cif_dict, file_id):
    cif_dict['_entry.id'] = file_id
    cif_dict['_struct.title'] = 'Virtual Structure Title'
    cif_dict['_exptl.method'] = 'X-RAY DIFFRACTION'
    return 

def parse_pdb(pdb_file):
    data = {
        "release_date": "1900-01-01",
        "resolution": None,
        "chain_ids": [],
        "seqs": [],
        "no_chains": 0
    }
    pdb_id = os.path.basename(pdb_file).split('.')[0]

    residues = {}

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("HEADER"):
                pdb_id = line[62:66].strip().lower()  # Extracting PDB ID
            elif line.startswith("REMARK   2 RESOLUTION"):
                try:
                    data["resolution"] = float(line.split()[3])
                except ValueError:
                    pass  # Handle cases where resolution is not a number
            elif line.startswith("ATOM"):
                chain_id = line[21]
                if chain_id not in data["chain_ids"]:
                    data["chain_ids"].append(chain_id)
                    data["no_chains"] += 1
                residue_seq_num = line[22:26].strip()
                residue_name = line[17:20].strip()
                residue_id = (chain_id, residue_seq_num)
                if residue_id not in residues:
                    residues[residue_id] = aa_codes.get(residue_name, 'X')  # 'X' for unknown residues

    # Building sequences for each chain
    for chain_id in data["chain_ids"]:
        chain_sequence = ''
        for res_id in sorted([r for r in residues if r[0] == chain_id], key=lambda x: int(x[1])):
            chain_sequence += residues[res_id]
        data["seqs"].append(chain_sequence)

    return pdb_id, data

def convert_pdb_to_mmcif(pdb_filename, mmcif_filename):
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_filename)
    file_id = pdb_filename.split(".pdb")[0]

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(mmcif_filename)
    
    cif_dict = {} 
    add_virtual_info_to_cif(cif_dict, file_id)
    
    # add release date
    revision_history = """
# 
loop_
_pdbx_audit_revision_history.ordinal
_pdbx_audit_revision_history.data_content_type
_pdbx_audit_revision_history.major_revision
_pdbx_audit_revision_history.minor_revision
_pdbx_audit_revision_history.revision_date 
1 'Structure model' 1 0 1900-01-01
2 'Structure model' 1 1 1900-01-02
3 'Structure model' 1 2 1900-01-03
#
"""

#     # add protein chains
#     protein_chains = """
# # 
# loop_
# _entity_poly_seq.entity_id
# _entity_poly_seq.num
# _entity_poly_seq.mon_id
# _entity_poly_seq.hetero
# 1 1   DC  n
# 1 2   DA  n
# 1 3   DT  n
# 1 4   DT  n
# 1 5   DA  n
# #
# """

    pdb_id, data = parse_pdb(pdb_filename)
    

    with open(mmcif_filename, 'a') as cif_file:
        for key, values in cif_dict.items():
            if isinstance(values, list):
                for value in values:
                    cif_file.write(f'{key}    {value}\n')
            else:
                cif_file.write(f'{key}    {values}\n')
        cif_file.write(revision_history)
    
    
    return 

def process_file(file_info):
    input_folder_path, filename, output_folder_path = file_info
    pdb_file = os.path.join(input_folder_path, filename)
    mmcif_file = os.path.join(output_folder_path, filename.replace('.pdb', '.cif'))
    convert_pdb_to_mmcif(pdb_file, mmcif_file)

def convert_pdb_files_in_folder(input_folder_path, output_folder_path, no_workers, chunksize):
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    pdb_files = [(input_folder_path, f, output_folder_path) for f in os.listdir(input_folder_path) if f.endswith(".pdb")]
    total_files = len(pdb_files)

    with Pool(no_workers) as pool:
        results = list(tqdm(pool.imap(process_file, pdb_files, chunksize=chunksize), total=total_files, desc="Converting files", unit="file"))

def main():
    parser = argparse.ArgumentParser(description="Convert PDB files to mmCIF format.")
    parser.add_argument("input_folder_path", type=str, help="Input folder path containing PDB files.")
    parser.add_argument("output_folder_path", type=str, help="Output folder path for mmCIF files.")
    parser.add_argument("--no_workers", type=int, default=4, help="Number of worker processes to use for conversion.")
    parser.add_argument("--chunksize", type=int, default=1, help="Chunk size for distributing work among processes.")

    args = parser.parse_args()

    convert_pdb_files_in_folder(args.input_folder_path, args.output_folder_path, args.no_workers, args.chunksize)

if __name__ == "__main__":
    main()
