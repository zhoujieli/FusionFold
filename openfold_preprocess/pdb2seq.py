from Bio.PDB import PDBParser, PPBuilder
import os
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning


warnings.simplefilter('ignore', PDBConstructionWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="Bio.PDB.DSSP")
def process_pdb_file(pdb_file_path, output_file, checkpoint_file, processed_files):
    if pdb_file_path in processed_files:
        return  # Skip processing if the file has already been processed

    parser = PDBParser()
    structure_id = os.path.splitext(os.path.basename(pdb_file_path))[0]
    structure = parser.get_structure(structure_id, pdb_file_path)
    written_sequences = set()  # Track sequences already written for this file

    with open(output_file, 'a') as file:
        for model in structure:
            chains = list(model.get_chains())
            if len(chains) == 1:  # Check if the structure is a monomer
                sequence = ""
                ppb = PPBuilder()
                print(f"Processing {structure_id}.pdb")
                for pp in ppb.build_peptides(chains[0]):
                    sequence += str(pp.get_sequence())
                if sequence and sequence not in written_sequences:  # Check for non-empty and new sequences
                    file.write(f">{sequence}\n")
                    written_sequences.add(sequence)

    with open(checkpoint_file, 'a') as chk_file:
        chk_file.write(f"{pdb_file_path}\n")

def process_pdb_files_in_directories(directory_paths, output_file, checkpoint_file):
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as chk_file:
            processed_files = {line.strip() for line in chk_file}
    else:
        processed_files = set()

    for directory_path in directory_paths:
        for filename in os.listdir(directory_path):
            if filename.endswith(".pdb"):
                pdb_file_path = os.path.join(directory_path, filename)
                process_pdb_file(pdb_file_path, output_file, checkpoint_file, processed_files)

# Set the list of directory paths and output file path
directory_paths =  ["/media/sampled_esm_pdbs/120000_150000", "/media/sampled_esm_pdbs/150000_180000", "/media/sampled_esm_pdbs/180000_210000","/media/sampled_esm_pdbs/220000_230000","/media/sampled_esm_pdbs/230000_240000","/media/sampled_esm_pdbs/250000_260000","/media/sampled_esm_pdbs/260000_270000"]
output_file = '../esm_pdbs/sequence_3.txt'
checkpoint_file = 'checkpoint_1.txt'

# Process the PDB files in the directory
process_pdb_files_in_directories(directory_paths, output_file, checkpoint_file)
