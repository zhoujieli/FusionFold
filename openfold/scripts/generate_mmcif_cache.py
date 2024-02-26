import argparse
from functools import partial
import json
import logging
from multiprocessing import Pool
import os

import sys
sys.path.append(".") # an innocent hack to get this to run from the top level

from tqdm import tqdm

from openfold.data.mmcif_parsing import parse 


def parse_file(f, args):
    try:
        with open(os.path.join(args.mmcif_dir, f), "r") as fp:
            mmcif_string = fp.read()
        file_id = os.path.splitext(f)[0]
        mmcif = parse(file_id=file_id, mmcif_string=mmcif_string)
        # print(str(mmcif))
        if mmcif.mmcif_object is None:
            logging.info(f"Could not parse {f}. Skipping...")
            return {}
        else:
            mmcif = mmcif.mmcif_object

    #     local_data = {}
    #     if "release_date" not in mmcif.header:
    #         local_data["release_date"] = "1990-01-01"
    #     else:
    #         local_data["release_date"] = mmcif.header["release_date"]

    #     chain_ids, seqs = list(zip(*mmcif.chain_to_seqres.items()))
    #     local_data["chain_ids"] = chain_ids
    #     local_data["seqs"] = seqs
    #     local_data["no_chains"] = len(chain_ids)

    #     local_data["resolution"] = mmcif.header["resolution"]

    #     return {file_id: local_data}
    # except KeyError as e:
    #     logging.error(f"KeyError in file {f}: {e}")
    #     return {}
        local_data = {}
        local_data["release_date"] = mmcif.header.get("release_date", "1990-01-01")
        local_data["resolution"] = mmcif.header.get("resolution", 0.0)

        # Handle chain IDs and sequences, assuming a default if missing
        chain_ids = ['A']  # Default chain ID
        seqs = ['Unknown']  # Default sequence
        if mmcif.chain_to_seqres:
            chain_ids, seqs = list(zip(*mmcif.chain_to_seqres.items()))
        local_data["chain_ids"] = chain_ids
        local_data["seqs"] = seqs
        local_data["no_chains"] = len(chain_ids)

        return {file_id: local_data}
    except Exception as e:
        logging.error(f"Error in file, Skipping...")
        return {}

def main(args):
    files = [f for f in os.listdir(args.mmcif_dir) if ".cif" in f]
    fn = partial(parse_file, args=args)
    data = {}
    with Pool(processes=args.no_workers) as p:
        with tqdm(total=len(files)) as pbar:
            for d in p.imap_unordered(fn, files, chunksize=args.chunksize):
                data.update(d)
                pbar.update()

    with open(args.output_path, "w") as fp:
        fp.write(json.dumps(data, indent=4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mmcif_dir", type=str, help="Directory containing mmCIF files"
    )
    parser.add_argument(
        "output_path", type=str, help="Path for .json output"
    )
    parser.add_argument(
        "--no_workers", type=int, default=4,
        help="Number of workers to use for parsing"
    )
    parser.add_argument(
        "--chunksize", type=int, default=10,
        help="How many files should be distributed to each worker at a time"
    )

    args = parser.parse_args()

    main(args)
