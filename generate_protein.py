import os
from chroma import Chroma, Protein

def load_checkpoint(checkpoint_file):
    try:
        with open(checkpoint_file, 'r') as file:
            return file.read().strip()
    except FileNotFoundError:
        return ''

def save_checkpoint(checkpoint_file, filename):
    with open(checkpoint_file, 'w') as file:
        file.write(filename)

chroma = Chroma(weights_backbone="/data_0/chroma_weights/chroma_backbone_v1.0.pt",
                weights_design="/data_0/chroma_weights/chroma_design_v1.0.pt").cuda()

input_folder = '/data_0/chroma_generation/beta'  
output_folder = '/data_0/chroma_generation/beta/round_2' 
checkpoint_file = 'beta_2nd.txt' 

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

last_processed = load_checkpoint(checkpoint_file)
start_processing = False if last_processed else True
processed_files_count = 0

for filename in sorted(os.listdir(input_folder)):
    if filename.endswith('.pdb'):
        if not start_processing:
            if filename == last_processed:
                start_processing = True
            continue

        try:
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, f"{os.path.splitext(filename)[0]}_redesign.pdb")

            protein = Protein(input_path, device='cuda')
            print(f"Processing {filename}...")

            protein = chroma.design(protein, design_selection="resid 350-700 around 5.0")

            protein.to(output_path)

            processed_files_count += 1
            print(f"Processed {filename} and saved to {output_path}")

            save_checkpoint(checkpoint_file, filename)

        except Exception as e:
            print(f"Error processing {filename}: {e}")

        print(f"Total number of files processed: {processed_files_count}")
