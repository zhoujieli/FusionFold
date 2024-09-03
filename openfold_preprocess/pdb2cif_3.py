import os
import argparse
from Bio.PDB import PDBParser, PPBuilder, MMCIFIO
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from tqdm import tqdm
from multiprocessing import Pool

warnings.simplefilter('ignore', PDBConstructionWarning)

def extract_sequence_from_structure(structure):
    ppb = PPBuilder()
    sequences = []

    # Iterate through the first model only as we assume homology in structures
    for chain in structure[0]:
        for pp in ppb.build_peptides(chain):
            sequences.append(str(pp.get_sequence()))

    # Combine sequences from all chains
    combined_sequence = "\n".join(sequences)
    return combined_sequence

def extract_sequence(structure):
    sequences_1 = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    sequences_1.append((chain.id, residue.id[1], residue.resname))
    return sequences_1

def convert_pdb_to_mmcif(pdb_filename, mmcif_filename):
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_filename)
    file_id = os.path.basename(pdb_filename).split(".pdb")[0]
    
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(mmcif_filename)
    
    sequence_data = extract_sequence_from_structure(structure)
    sequence_info = extract_sequence(structure)
    # add more info
    add_critical_info =f"""
_entry.id    {file_id}
_struct.title   'Virtual Structure Title'
_exptl.method   'X-RAY DIFFRACTION'
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
"""     
    
    # Example sequence information to be formatted and added
    sequence_section = f"""
# 
_entity_poly.entity_id          1
_entity_poly.type               'polypeptide(L)'
_entity_poly.nstd_linkage       no
_entity_poly.nstd_monomer       yes
_entity_poly.pdbx_seq_one_letter_code      
;{sequence_data}
;
_entity_poly.pdbx_strand_id                 A
_entity_poly.pdbx_target_identifier         ?
#
"""

    with open(mmcif_filename, 'a') as cif_file:
        cif_file.write(add_critical_info)
        cif_file.write(sequence_section)
    
    
    # Preparing the sequence section with correct formatting
    sequence_section_new = """ 
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
"""
    for chain_id, num, mon_id in sequence_info:
        sequence_section_new += f"1 {num} {mon_id} n\n"

    with open(mmcif_filename, 'a') as cif_file:
        cif_file.write(sequence_section_new)
        cif_file.write('#')
    
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
        list(tqdm(pool.imap(process_file, pdb_files, chunksize=chunksize), total=total_files, desc="Converting files", unit="file"))

def main():
    parser = argparse.ArgumentParser(description="Convert PDB files to mmCIF format with sequence information.")
    parser.add_argument("input_folder_path", type=str, help="Input folder path containing PDB files.")
    parser.add_argument("output_folder_path", type=str, help="Output folder path for mmCIF files.")
    parser.add_argument("--no_workers", type=int, default=4, help="Number of worker processes to use for conversion.")
    parser.add_argument("--chunksize", type=int, default=1, help="Chunk size for distributing work among processes.")

    args = parser.parse_args()

    convert_pdb_files_in_folder(args.input_folder_path, args.output_folder_path, args.no_workers, args.chunksize)

if __name__ == "__main__":
    main()

