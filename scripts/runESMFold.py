# ================================
#  ESMFold Parallel Research Template
# ================================

import os
import re
import time
import argparse
import numpy as np
import torch
import esm
import torch.distributed as dist
from Bio import SeqIO

# --------------------------------------------------
# Utils
# --------------------------------------------------

def to_numpy(x):
    if torch.is_tensor(x):
        return x.cpu().numpy()
    if isinstance(x, dict):
        return {k: to_numpy(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return type(x)(to_numpy(v) for v in x)
    return x


def sanitize_sequence(seq):
    seq = re.sub("[^A-Z:]", "", seq.replace("/", ":").upper())
    seq = re.sub(":+", ":", seq)
    return seq.strip(":")


# --------------------------------------------------
# Parallel setup
# --------------------------------------------------

def setup_parallel(args):
    """
    Returns:
        parallel_mode
        device (for single / ddp) OR device_ids (for dp)
        rank, world_size (ddp only)
    """
    if args.cpu:
        return "cpu", None, 0, 1

    if args.parallel_mode == "ddp":
        dist.init_process_group(backend="nccl")
        rank = int(os.environ["RANK"])
        world_size = int(os.environ["WORLD_SIZE"])
        local_rank = int(os.environ["LOCAL_RANK"])
        torch.cuda.set_device(local_rank)
        device = torch.device(f"cuda:{local_rank}")
        return "ddp", device, rank, world_size

    # single / dp
    if args.device:
        os.environ["CUDA_VISIBLE_DEVICES"] = args.device

    assert torch.cuda.is_available(), "CUDA not available"
    n = torch.cuda.device_count()
    device_ids = list(range(n))
    return args.parallel_mode, device_ids, 0, 1


# --------------------------------------------------
# Model builder
# --------------------------------------------------

def build_model(parallel_mode, device_or_ids, seq_length):
    base_model = esm.pretrained.esmfold_v1().eval()

    if parallel_mode == "single":
        model = base_model.to("cuda:0")

    elif parallel_mode == "dp":
        assert len(device_or_ids) > 1, "DP requires multiple GPUs"
        base_model = base_model.to("cuda:0")
        model = torch.nn.DataParallel(
            base_model,
            device_ids=device_or_ids,
            output_device=0
        ).module

    elif parallel_mode == "ddp":
        model = base_model.to(device_or_ids)

    else:
        raise ValueError(f"Unknown parallel_mode {parallel_mode}")

    # chunk size strategy (critical for long sequences)
    if seq_length > 1100:
        model.set_chunk_size(64)
        if seq_length > 1250:
            model.set_chunk_size(32)
    else:
        model.set_chunk_size(128)

    return model


# --------------------------------------------------
# ESMFold inference
# --------------------------------------------------

def run_esmfold(jobname, sequence, model, device, out_dir, recycles):
    with torch.no_grad():
        output = model.infer(
            sequence,
            num_recycles=recycles,
            chain_linker="X" * 25,
            residue_index_offset=512
        )

    pdb = model.output_to_pdb(output)[0]
    output = to_numpy(output)

    ptm = output["ptm"][0]
    plddt = output["plddt"][0, ..., 1].mean()

    with open(os.path.join(out_dir, f"{jobname}.pdb"), "w") as f:
        f.write(pdb)

    with open(os.path.join(out_dir, "output.txt"), "a") as f:
        f.write(f"{jobname} ptm:{ptm:.3f} plddt:{plddt:.3f}\n")


# --------------------------------------------------
# Main
# --------------------------------------------------

def main(args):
    parallel_mode, dev, rank, world_size = setup_parallel(args)

    if rank == 0:
        os.makedirs(args.store_folder, exist_ok=True)

    if parallel_mode == "ddp":
        dist.barrier()

    records = [(r.id, sanitize_sequence(str(r.seq)))
               for r in SeqIO.parse(args.fasta_file, "fasta")]
    records.sort(key=lambda x: len(x[1]))

    # DDP sharding
    if parallel_mode == "ddp":
        records = records[rank::world_size]

    recycles = args.recycles if args.recycles is not None else 3

    for jobname, seq in records:
        length = len(seq.replace(":", ""))
        model = build_model(parallel_mode, dev, length)
        run_esmfold(jobname, seq, model,
                    dev if parallel_mode == "ddp" else None,
                    args.store_folder, recycles)

    if parallel_mode == "ddp":
        dist.barrier()
        if rank == 0:
            print("All jobs finished.")


# --------------------------------------------------
# Entry
# --------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser("ESMFold Parallel Template")

    parser.add_argument("-i", dest="fasta_file", required=True)
    parser.add_argument("-o", dest="store_folder", required=True)

    parser.add_argument("-cpu", action="store_true")
    parser.add_argument("-gpu", action="store_true")

    parser.add_argument("--device", type=str,
                        help="Physical GPU ids, e.g. 0,1,2")

    parser.add_argument("--parallel_mode", type=str, default="single",
                        choices=["single", "dp", "ddp"])

    parser.add_argument("-recycles", type=int)

    args = parser.parse_args()
    t0 = time.time()
    main(args)
    print(f"Total time: {time.time() - t0:.1f}s")
