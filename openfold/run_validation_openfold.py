import argparse
import logging
import math
import numpy as np
import os
import pandas as pd

from openfold.utils.script_utils import load_models_from_command_line, parse_fasta, run_model, prep_output, \
    update_timings, relax_protein
logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)

import pickle

import random
import time
import torch
torch_versions = torch.__version__.split(".")
torch_major_version = int(torch_versions[0])
torch_minor_version = int(torch_versions[1])
if(
    torch_major_version > 1 or
    (torch_major_version == 1 and torch_minor_version >= 12)
):
    # Gives a large speedup on Ampere-class GPUs
    torch.set_float32_matmul_precision("high")

torch.set_grad_enabled(False)

from openfold.config import model_config
from openfold.data import templates, feature_pipeline, data_pipeline

from openfold.utils.trace_utils import (
    pad_feature_dict_seq,
    trace_model_,
)
from scripts.utils import add_data_args
from openfold.utils.loss import AlphaFoldLoss
from tqdm import tqdm
import torch.multiprocessing as mp
from multiprocessing import Manager


TRACING_INTERVAL = 50

def list_files_with_extensions(dir, extensions):
    return [f for f in os.listdir(dir) if f.endswith(extensions)]

def precompute_alignments(tags, seqs, alignment_dir, args):
    for tag, seq in zip(tags, seqs):
        tmp_fasta_path = os.path.join(args.output_dir, f"tmp_{os.getpid()}.fasta")
        with open(tmp_fasta_path, "w") as fp:
            fp.write(f">{tag}\n{seq}")

        local_alignment_dir = os.path.join(alignment_dir, tag)
        if(args.use_precomputed_alignments is None and not os.path.isdir(local_alignment_dir)):
            logger.info(f"Generating alignments for {tag}...")

            os.makedirs(local_alignment_dir)

            alignment_runner = data_pipeline.AlignmentRunner(
                    jackhmmer_binary_path=args.jackhmmer_binary_path,
                    hhblits_binary_path=args.hhblits_binary_path,
                    hhsearch_binary_path=args.hhsearch_binary_path,
                    uniref90_database_path=args.uniref90_database_path,
                    mgnify_database_path=args.mgnify_database_path,
                    bfd_database_path=args.bfd_database_path,
                    uniclust30_database_path=args.uniclust30_database_path,
                    pdb70_database_path=args.pdb70_database_path,
                    no_cpus=args.cpus,
                )
            alignment_runner.run(
                tmp_fasta_path, local_alignment_dir
            )
        else:
            # logger.info(
            #     f"Using precomputed alignments for {tag} at {alignment_dir}..."
            # )
            pass

        # Remove temporary FASTA file
        os.remove(tmp_fasta_path)

def round_up_seqlen(seqlen):
    return int(math.ceil(seqlen / TRACING_INTERVAL)) * TRACING_INTERVAL

def generate_feature_dict(
    tags,
    seqs,
    alignment_dir,
    data_processor,
    args,
):
    tmp_fasta_path = os.path.join(args.output_dir, f"tmp_{os.getpid()}.fasta")
    if len(seqs) == 1:
        tag = tags[0]
        seq = seqs[0]
        with open(tmp_fasta_path, "w") as fp:
            fp.write(f">{tag}\n{seq}")

        local_alignment_dir = os.path.join(alignment_dir, tag)
        feature_dict = data_processor.process_fasta(
            fasta_path=tmp_fasta_path,
            alignment_dir=local_alignment_dir,
            seqemb_mode=args.use_single_seq_mode,
        )
    else:
        with open(tmp_fasta_path, "w") as fp:
            fp.write(
                '\n'.join([f">{tag}\n{seq}" for tag, seq in zip(tags, seqs)])
            )
        feature_dict = data_processor.process_multiseq_fasta(
            fasta_path=tmp_fasta_path, super_alignment_dir=alignment_dir,
        )

    # Remove temporary FASTA file
    os.remove(tmp_fasta_path)

    return feature_dict


def infer_seqences(infer_model, device, sequences, alignment_dir, data_processor, feature_processor, feature_dicts, args, progress_queue):
    cur_tracing_interval = 0
    infer_result = []
    for (tag, tags), seqs in sequences:
        output_name = f'{tag}_{args.config_preset}'
        if args.output_postfix is not None:
            output_name = f'{output_name}_{args.output_postfix}'

        # Does nothing if the alignments have already been computed
        precompute_alignments(tags, seqs, alignment_dir, args)

        feature_dict = feature_dicts.get(tag, None)
        if(feature_dict is None):
            feature_dict = generate_feature_dict(
                tags,
                seqs,
                alignment_dir,
                data_processor,
                args,
            )

            feature_dicts[tag] = feature_dict

        processed_feature_dict = feature_processor.process_features(
            feature_dict, mode='predict',
        )

        processed_feature_dict = {
            # k:torch.as_tensor(v, device=args.model_device)
            k:torch.as_tensor(v, device=device)
            for k,v in processed_feature_dict.items()
        }

        out = run_model(infer_model, processed_feature_dict, tag, args.output_dir)
        for k,v in out.items():
            if k == "plddt":
                # 计算plddt的均值
                plddt_mean = torch.mean(v)
                # 将tensor值转换为数值
                plddt_mean = plddt_mean.item()
                infer_result.append((tag, tags, seqs, plddt_mean))
            # # 计算loss
                # loss = loss_func(out, processed_feature_dict)
        # update progress bar
        progress_queue.put(1)

    return infer_result

 # 定义worker
def worker(model, device, subset, alignment_dir, data_processor,feature_processor,feature_dicts, args, results_queue, progress_queue):
    # 这里应该是调用您的模型推理函数
    result = infer_seqences(model, device, subset, alignment_dir, data_processor, feature_processor, feature_dicts, args, progress_queue)
    results_queue.put(result)

def main(args):
    # Create the output directory
    os.makedirs(args.output_dir, exist_ok=True)

    config = model_config(args.config_preset, long_sequence_inference=args.long_sequence_inference)
    # loss_function
    loss_func = AlphaFoldLoss(config.loss)

    template_featurizer = templates.TemplateHitFeaturizer(
        mmcif_dir=args.template_mmcif_dir,
        max_template_date=args.max_template_date,
        max_hits=config.data.predict.max_templates,
        kalign_binary_path=args.kalign_binary_path,
        release_dates_path=args.release_dates_path,
        obsolete_pdbs_path=args.obsolete_pdbs_path
    )

    data_processor = data_pipeline.DataPipeline(
        template_featurizer=template_featurizer,
    )

    output_dir_base = args.output_dir
    random_seed = args.data_random_seed
    if random_seed is None:
        random_seed = random.randrange(2**32)

    np.random.seed(random_seed)
    torch.manual_seed(random_seed + 1)

    feature_processor = feature_pipeline.FeaturePipeline(config.data)
    if not os.path.exists(output_dir_base):
        os.makedirs(output_dir_base)
    if args.use_precomputed_alignments is None:
        alignment_dir = os.path.join(output_dir_base, "alignments")
    else:
        alignment_dir = args.use_precomputed_alignments

    validation_num = 16
    tag_list = []
    seq_list = []
    cnt = 0
    for fasta_file in list_files_with_extensions(args.fasta_dir, (".fasta", ".fa")):
        if cnt < validation_num:
            
            # Gather input sequences
            with open(os.path.join(args.fasta_dir, fasta_file), "r") as fp:
                data = fp.read()

            tags, seqs = parse_fasta(data)
            # assert len(tags) == len(set(tags)), "All FASTA tags must be unique"
            tag = '-'.join(tags)

            tag_list.append((tag, tags))
            seq_list.append(seqs)
            cnt += 1
        else:
            break

    seq_sort_fn = lambda target: sum([len(s) for s in target[1]])
    sorted_targets = sorted(zip(tag_list, seq_list), key=seq_sort_fn)
    feature_dicts = {}

    # 解析model_device, model_device = 'cuda:0;cuda:1;cuda:2;cuda:3'
    model_generator = []
    model_devices = args.model_device.split(';')
    num_gpus = len(model_devices)
    for i in range(len(model_devices)):
        model_generator.append(load_models_from_command_line(
            config, model_devices[i], args.openfold_checkpoint_path, args.jax_param_path, args.output_dir
        ))

    # 加入multi-process progress bar
    manager = Manager()
    progress_queue = manager.Queue()

    # 将即将推理的数据按照GPU数量进行分组
    grouped_sequences = np.array_split(sorted_targets, num_gpus)
    
    # 加入multi-process progress bar: set total number of tasks
    total_tasks = sum(len(subset) for subset in grouped_sequences)
    processes = []
    results_queue = mp.Queue()

    # 多进程
    for i, subset in enumerate(grouped_sequences):
        device = f'cuda:{i}' if torch.cuda.is_available() else 'cpu'
        for model, output_directory in model_generator[i]:
            p = mp.Process(target=worker, args=(model, device, subset, alignment_dir, data_processor, feature_processor,feature_dicts, args, results_queue, progress_queue))
            p.start()
            processes.append(p)

    # Progress bar monitoring
    pbar = tqdm(total=total_tasks, desc="Processing Sequences")
    while any(p.is_alive() for p in processes):
        while not progress_queue.empty():
            progress_queue.get()
            pbar.update(1)
    pbar.close()

    # 保存结果
    # 引入infer_result用于保存推理结果，包括tag, tags, seqs, plddt的值
    infer_result = []

    for p in processes:
        p.join()
        infer_result.append(results_queue.get())

    combined_results = [item for sublist in infer_result for item in sublist]
    
    # 将infer_result写入csv文件中，并删选所有plddt小于50的序列写入另外一个csv文件中
    infer_result = pd.DataFrame(combined_results, columns=["tag", "tags", "seqs", "plddt"])
    infer_result.to_csv("infer_result.csv", index=False)
    infer_result = infer_result[infer_result["plddt"] < 50]
    infer_result.to_csv("infer_result_plddt_less_than_50.csv", index=False)


if __name__ == "__main__":
    mp.set_start_method('spawn', force=True) 
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta_dir", type=str,
        help="Path to directory containing FASTA files, one sequence per file"
    )
    parser.add_argument(
        "template_mmcif_dir", type=str,
    )
    parser.add_argument(
        "--use_precomputed_alignments", type=str, default=None,
        help="""Path to alignment directory. If provided, alignment computation 
                is skipped and database path arguments are ignored."""
    )
    parser.add_argument(
        "--use_single_seq_mode", action="store_true", default=False,
        help="""Use single sequence embeddings instead of MSAs."""
    )
    parser.add_argument(
        "--output_dir", type=str, default=os.getcwd(),
        help="""Name of the directory in which to output the prediction""",
    )
    parser.add_argument(
        "--model_device", type=str, default="cpu",
        help="""Name of the device on which to run the model. Any valid torch
             device name is accepted (e.g. "cpu", "cuda:0")"""
    )
    parser.add_argument(
        "--config_preset", type=str, default="model_1",
        help="""Name of a model config preset defined in openfold/config.py"""
    )
    parser.add_argument(
        "--jax_param_path", type=str, default=None,
        help="""Path to JAX model parameters. If None, and openfold_checkpoint_path
             is also None, parameters are selected automatically according to 
             the model name from openfold/resources/params"""
    )
    parser.add_argument(
        "--openfold_checkpoint_path", type=str, default=None,
        help="""Path to OpenFold checkpoint. Can be either a DeepSpeed 
             checkpoint directory or a .pt file"""
    )
    parser.add_argument(
        "--save_outputs", action="store_true", default=False,
        help="Whether to save all model outputs, including embeddings, etc."
    )
    parser.add_argument(
        "--cpus", type=int, default=4,
        help="""Number of CPUs with which to run alignment tools"""
    )
    parser.add_argument(
        "--preset", type=str, default='full_dbs',
        choices=('reduced_dbs', 'full_dbs')
    )
    parser.add_argument(
        "--output_postfix", type=str, default=None,
        help="""Postfix for output prediction filenames"""
    )
    parser.add_argument(
        "--data_random_seed", type=int, default=None
    )
    parser.add_argument(
        "--skip_relaxation", action="store_true", default=False,
    )
    parser.add_argument(
        "--multimer_ri_gap", type=int, default=200,
        help="""Residue index offset between multiple sequences, if provided"""
    )
    parser.add_argument(
        "--subtract_plddt", action="store_true", default=False,
        help=""""Whether to output (100 - pLDDT) in the B-factor column instead
                 of the pLDDT itself"""
    )
    parser.add_argument(
        "--long_sequence_inference", action="store_true", default=False,
        help="""enable options to reduce memory usage at the cost of speed, helps longer sequences fit into GPU memory, see the README for details"""
    )
    parser.add_argument(
        "--cif_output", action="store_true", default=False,
        help="Output predicted models in ModelCIF format instead of PDB format (default)"
    )
    add_data_args(parser)
    args = parser.parse_args()

    if(args.jax_param_path is None and args.openfold_checkpoint_path is None):
        args.jax_param_path = os.path.join(
            "openfold", "resources", "params",
            "params_" + args.config_preset + ".npz"
        )

    if(args.model_device == "cpu" and torch.cuda.is_available()):
        logging.warning(
            """The model is being run on CPU. Consider specifying 
            --model_device for better performance"""
        )

    main(args)