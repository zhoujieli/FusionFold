import argparse
import logging
import os
from multiprocessing import Pool

from tqdm import tqdm
from openfold.data import mmcif_parsing
from openfold.np import protein, residue_constants


def process_cif_file(args):
    fname, data_dir, raise_errors = args
    fasta = []
    basename, ext = os.path.splitext(fname)
    # basename = basename.upper()
    fpath = os.path.join(data_dir, fname)
    if ext == ".cif":
        with open(fpath, 'r') as fp:
            mmcif_str = fp.read()
        
        mmcif = mmcif_parsing.parse(
            file_id=basename, mmcif_string=mmcif_str
        )
        if mmcif.mmcif_object is None:
            logging.warning(f'Failed to parse {fname}...')
            if raise_errors:
                raise list(mmcif.errors.values())[0]
            else:
                return None, None  # Return None to indicate failure to parse

        mmcif = mmcif.mmcif_object
        for chain, seq in mmcif.chain_to_seqres.items():
            chain_id = '_'.join([basename, chain])
            fasta.append(f">{chain_id}")
            fasta.append(seq)
        return fasta, os.path.splitext(fname)[0] + '.fasta'
    elif ext == ".core":
        with open(fpath, 'r') as fp:
            core_str = fp.read()

        core_protein = protein.from_proteinnet_string(core_str)
        aatype = core_protein.aatype
        seq = ''.join([
            residue_constants.restypes_with_x[aatype[i]] 
            for i in range(len(aatype))
        ])
        fasta.append(f">{basename}")
        fasta.append(seq)
        return fasta, os.path.splitext(fname)[0] + '.fasta'

    return None, None

def main(args):
    files = [(fname, args.data_dir, args.raise_errors) for fname in os.listdir(args.data_dir)]
    with Pool() as pool:
        results = list(tqdm(pool.imap(process_cif_file, files), total=len(files)))

    for fasta, output_filename in results:
        if fasta and output_filename:
            output_path = os.path.join(args.output_dir, output_filename)
            with open(output_path, "w") as fp:
                fp.write('\n'.join(fasta))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data_dir", type=str,
        help="Path to a directory containing mmCIF or .core files"
    )
    parser.add_argument(
        "output_dir", type=str,
        help="Directory to save output FASTA files"
    )
    parser.add_argument(
        "--raise_errors", type=bool, default=False,
        help="Whether to crash on parsing errors"
    )

    args = parser.parse_args()

    main(args)
