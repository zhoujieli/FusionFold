import os
import torch
from chroma import Chroma, Protein
from Bio.PDB import PDBParser
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import pickle

# Function to check if a PDB file is a monomer
def is_monomer(pdb_file_path):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file_path)
    chains = list(structure.get_chains())
    return len(chains) == 1

# Function to load the last processed file from the checkpoint
def load_checkpoint(checkpoint_file):
    try:
        with open(checkpoint_file, 'rb') as file:
            return pickle.load(file)
    except FileNotFoundError:
        return set()

# Function to save the current file to the checkpoint
def save_checkpoint(checkpoint_file, processed_files):
    with open(checkpoint_file, 'wb') as file:
        pickle.dump(processed_files, file)

# Initialize Chroma
chroma = Chroma(weights_backbone="/data_0/chroma_weights/chroma_backbone_v1.0.pt",
                weights_design="/data_0/chroma_weights/chroma_design_v1.0.pt").cuda()

# Directories and checkpoint file
input_folders = ['/data_0/input_pdb']  # Add other folders as needed
output_folder = '/data_0/msa_generation/gen2' 
checkpoint_file = 'msa2.pkl'  # Changed to .pkl for pickle format

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Load processed files from checkpoint
processed_files = load_checkpoint(checkpoint_file)

# Process each folder
for input_folder in input_folders:
    for filename in sorted(os.listdir(input_folder)):
        if filename.endswith('.pdb') and filename not in processed_files:
            pdb_filepath = os.path.join(input_folder, filename)
            
            if not is_monomer(pdb_filepath):
                continue  # Skip if not a monomer

            try:
                input_path = pdb_filepath
                output_path = os.path.join(output_folder, f"{os.path.splitext(filename)[0]}_redesign_1.pdb")

                # Load the protein
                protein = Protein(input_path, device='cuda')
                print(f"Processing {filename}...")

                protein = chroma.design(protein, design_selection="resid 1-700 around 5.0")

                protein.to(output_path)

                print(f"Processed {filename} and saved to {output_path}")

                # Update the set of processed files and save to checkpoint
                processed_files.add(filename)
                save_checkpoint(checkpoint_file, processed_files)
                print(f"Total number of files processed: {len(processed_files)}")

            except Exception as e:
                print(f"Error processing {filename}: {e}")

