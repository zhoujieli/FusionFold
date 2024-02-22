# Usage : python pdb2fasta.py /data_0/chroma_generation/beta ./input.fasta --num_pdb 100
# Usage2: python pdb2fasta.py /data_0/chroma_generation/beta2 ./input.fasta --all_pdb


import argparse
import os
from Bio import PDB
from tqdm import tqdm

def extract_and_save_fasta(pdb_dir, fasta_path, num_pdb=None):
    with open(fasta_path, 'w') as file:
        pdb_files = [filename for filename in os.listdir(pdb_dir) if filename.endswith(".pdb")]
        if num_pdb is not None:
            pdb_files = pdb_files[:num_pdb]

        for filename in tqdm(pdb_files, desc="Converting PDB to FASTA"):
            pdb_path = os.path.join(pdb_dir, filename)
            pdb_id = os.path.splitext(filename)[0]
            pdb_parser = PDB.PDBParser(QUIET=True)
            structure = pdb_parser.get_structure(pdb_id, pdb_path)

            for model in structure:
                for chain in model:
                    sequence = []
                    for residue in chain:
                        if PDB.is_aa(residue, standard=True):
                            sequence.append(PDB.Polypeptide.three_to_one(residue.get_resname()))
                    fasta_sequence = ''.join(sequence)
                    file.write(f">{pdb_id}_{model.id}_{chain.id}\n{fasta_sequence}\n")

def main():
    # 命令行参数解析
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
        num_pdb = None  # 处理所有PDB文件
    else:
        num_pdb = args.num_pdb

    extract_and_save_fasta(args.pdb_dir, args.fasta_file, num_pdb)

if __name__ == "__main__":
    main()

