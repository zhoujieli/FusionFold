import os
import json
import argparse
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# Assuming aa_codes and parse_pdb functions are defined here (from previous examples)
aa_codes = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

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

def process_pdb_file(pdb_file_path):
    pdb_id, pdb_data = parse_pdb(pdb_file_path)
    return pdb_id, pdb_data

def main(args):
    pdb_file_path = args.pdb_dir
    output_path = args.output_path
    no_workers = args.no_workers
    chunksize = args.chunksize

    pdb_files = [os.path.join(pdb_file_path, f) for f in os.listdir(pdb_file_path) if f.endswith('.pdb')]
    total_files = len(pdb_files)

    print(f"Processing {total_files} files with {no_workers} workers...")

    with Pool(no_workers) as pool:
        results = list(tqdm(pool.imap(process_pdb_file, pdb_files, chunksize=chunksize), total=total_files))

    final_data = {pdb_id: data for pdb_id, data in results}

    with open(output_path, "w") as outfile:
        json.dump(final_data, outfile, indent=4)

    print(f"Data extracted and saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_dir", type=str, help="Directory containing PDB files")
    parser.add_argument("output_path", type=str, help="Path for JSON output")
    parser.add_argument("--no_workers", type=int, default=cpu_count(), help="Number of workers to use for parsing")
    parser.add_argument("--chunksize", type=int, default=10, help="How many files should be distributed to each worker at a time")

    args = parser.parse_args()
    main(args)
