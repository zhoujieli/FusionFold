import os
from Bio.PDB import PDBParser, MMCIFIO, MMCIF2Dict
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning


warnings.simplefilter('ignore', PDBConstructionWarning)
def convert_pdb_to_mmcif(pdb_filename, mmcif_filename):
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_filename)

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(mmcif_filename)

def add_virtual_info_to_cif(cif_dict, file_id):
    cif_dict['_entry.id'] = file_id
    cif_dict['_struct.title'] = 'Virtual Structure Title'
    cif_dict['_exptl.method'] = 'X-RAY DIFFRACTION'
    return 

def pdb_to_cif_with_virtual_info(pdb_file_path, output_cif_path):
    """
    Convert a PDB file to a CIF file and add virtual information.

    Args:
    pdb_file_path: Path to the input PDB file.
    output_cif_path: Path where the output CIF file will be saved.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_structure", pdb_file_path)
    file_id = pdb_file_path.split(".pdb")[0]

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_cif_path)

    # cif_dict = MMCIF2Dict.MMCIF2Dict(output_cif_path)
    # print(cif_dict)
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
# 
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1   DC  n
1 2   DA  n
1 3   DT  n
1 4   DT  n
1 5   DA  n
#
# """

    with open(output_cif_path, 'a') as cif_file:
        for key, values in cif_dict.items():
            if isinstance(values, list):
                for value in values:
                    cif_file.write(f'{key}    {value}\n')
            else:
                cif_file.write(f'{key}    {values}\n')
        cif_file.write(revision_history)
        # cif_file.write(protein_chains)
        
    return 

# pdb_file_path = 'MGYP002510725049_redesign_1.pdb'  
# output_cif_path = 'test_cif/output.cif'

# The `pdb_to_cif_with_virtual_info` function takes a PDB file as input and converts it to a CIF file.
# It also adds virtual information to the CIF file. The virtual information includes the entry ID,
# structure title, and experimental method. The function uses the `PDBParser` from the `Bio.PDB`
# module to parse the PDB file and create a structure object. Then, it uses the `MMCIFIO` class to
# save the structure object as a CIF file. After saving the CIF file, the function adds the virtual
# information to the CIF file by modifying the CIF dictionary. Finally, the function writes the
# modified CIF dictionary to the CIF file.
# pdb_to_cif_with_virtual_info(pdb_file_path, output_cif_path)

def convert_pdb_files_in_folder(input_folder_path, output_folder_path):
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    for filename in os.listdir(input_folder_path):
        if filename.endswith(".pdb"):
            pdb_file = os.path.join(input_folder_path, filename)
            mmcif_file = os.path.join(output_folder_path, filename.replace('.pdb', '.cif'))
            convert_pdb_to_mmcif(pdb_file, mmcif_file)
            print(f"Converted {pdb_file} to {mmcif_file}")

input_folder_path = 'home/shenyichong/data_process'  
output_folder_path = 'home/shenyichong/data_process' 
convert_pdb_files_in_folder(input_folder_path, output_folder_path)
