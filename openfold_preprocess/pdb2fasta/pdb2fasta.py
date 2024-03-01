# Usage : python pdb2fasta.py /data_0/chroma_generation/beta ./input.fasta --num_pdb 100
# Usage2: python pdb2fasta.py /data_0/chroma_generation/beta2 ./input.fasta --all_pdb


import argparse
import os
from Bio import PDB
from tqdm import tqdm
import multiprocessing

def worker(args):
    pdb_path, return_dict, index = args
    fasta_sequences = ""
    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
    pdb_parser = PDB.PDBParser(QUIET=True)
    try:
        structure = pdb_parser.get_structure(pdb_id, pdb_path)
        for model in structure:
            for chain in model:
                sequence = []
                for residue in chain:
                    if PDB.is_aa(residue, standard=True):
                        sequence.append(PDB.Polypeptide.three_to_one(residue.get_resname()))
                fasta_sequence = ''.join(sequence)
                fasta_sequences += f">{pdb_id}_{model.id}_{chain.id}\n{fasta_sequence}\n"
        return_dict[index] = fasta_sequences
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")

def extract_and_save_fasta(pdb_dir, fasta_path, num_pdb=None):
    pdb_files = [os.path.join(pdb_dir, filename) for filename in os.listdir(pdb_dir) if filename.endswith(".pdb")]
    if num_pdb is not None:
        pdb_files = pdb_files[:num_pdb]

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    pool = multiprocessing.Pool(16)  # Specify 16 processes

    tasks = [(pdb_files[i], return_dict, i) for i in range(len(pdb_files))]
    for _ in tqdm(pool.imap_unordered(worker, tasks), total=len(pdb_files), desc="Converting PDB to FASTA"):
        pass

    pool.close()
    pool.join()

    with open(fasta_path, 'w') as file:
        for key in sorted(return_dict.keys()):
            file.write(return_dict[key])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "pdb_dir", type=str,
        help="Path to directory containing PDB files."
    )
    parser.add_argument(
        "fasta_file", type=str,
        help="Output FASTA file."
    )
    parser.add_argument(
        "--num_pdb", type=int, default=None,
        help="Number of PDB files to process."
    )
    parser.add_argument(
        "--all_pdb", action="store_true",
        help="Process all PDB files in the directory."
    )

    args = parser.parse_args()

    if args.all_pdb:
        num_pdb = None
    else:
        num_pdb = args.num_pdb

    extract_and_save_fasta(args.pdb_dir, args.fasta_file, num_pdb)

if __name__ == "__main__":
    main()
