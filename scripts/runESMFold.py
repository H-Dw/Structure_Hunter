# Refer to https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ESMFold.ipynb

from string import ascii_uppercase, ascii_lowercase
import hashlib, re, os
import numpy as np
import esm
from jax.tree_util import tree_map #pip install jax,jaxlib
from scipy.special import softmax
from Bio import SeqIO
import sys
import gc
import zipfile
import subprocess

import torch
import torch.nn as nn
import argparse
import time



def parse_output(output):
  pae = (output["aligned_confidence_probs"][0] * np.arange(64)).mean(-1) * 31
  plddt = output["plddt"][0,:,1]

  bins = np.append(0,np.linspace(2.3125,21.6875,63))
  sm_contacts = softmax(output["distogram_logits"],-1)[0]
  sm_contacts = sm_contacts[...,bins<8].sum(-1)
  xyz = output["positions"][-1,0,:,1]
  mask = output["atom37_atom_exists"][0,:,1] == 1
  o = {"pae":pae[mask,:][:,mask],
       "plddt":plddt[mask],
       "sm_contacts":sm_contacts[mask,:][:,mask],
       "xyz":xyz[mask]}
  return o

def runESMFold(jobname, sequence, deviceC, device_ids, store_folder, recycles):
  jobname = re.sub(r'\W+', '', jobname)[:50]
  sequence = re.sub("[^A-Z:]", "", sequence.replace("/",":").upper())
  sequence = re.sub(":+",":",sequence)
  sequence = re.sub("^[:]+","",sequence)
  sequence = re.sub("[:]+$","",sequence)
  copies = 1 #@param {type:"integer"}
  if copies == "" or copies <= 0: copies = 1
  sequence = ":".join([sequence] * copies)
  num_recycles = recycles
  chain_linker = 25

  ID = jobname
  seqs = sequence.split(":")
  lengths = [len(s) for s in seqs]
  length = sum(lengths)
  print("length",length)

  u_seqs = list(set(seqs))
  if len(seqs) == 1: mode = "mono"
  elif len(u_seqs) == 1: mode = "homo"
  else: mode = "hetero"

  if deviceC == "cpu":
    model = esm.pretrained.esmfold_v1()
    model = model.float()
    model = model.eval().to("cpu")
    actual_model = model
  else: 
    #nvidia-smi
    model = torch.nn.DataParallel(esm.pretrained.esmfold_v1(), device_ids=device_ids)
    model = model.eval().cuda()
    # model = nn.DataParallel(model.eval(), visible_devices).cuda()
    # nn.DataParallel(model.eval().cuda(), visible_devices)
    # nn.DataParallel(model.eval().cuda(), device_ids=[1])

    # model.module.set_chunk_size(128)
    # optimized for different GPU 128, 64, 32
    if length > 1100:
      model.module.set_chunk_size(64) 
      if length > 1250:
        model.module.set_chunk_size(32)
    else:
      model.module.set_chunk_size(128)

    actual_model = model.module

  with torch.no_grad():
    output = actual_model.infer(sequence,
                        num_recycles=num_recycles,
                        chain_linker="X"*chain_linker,
                        residue_index_offset=512)

  pdb_str = actual_model.output_to_pdb(output)[0]
  
  output = tree_map(lambda x: x.cpu().numpy(), output)
  ptm = output["ptm"][0]
  plddt = output["plddt"][0,...,1].mean()
  O = parse_output(output)

  output_text = f'{ID} ptm: {ptm:.3f} plddt: {plddt:.3f}'
  print(output_text)
  outputfile = f"{store_folder}/output.txt"
  
  with open(f"{outputfile}", "a") as file:
    file.write(output_text + "\n")
  # print(f'ptm: {ptm:.3f} plddt: {plddt:.3f}') # before1017


  # os.system(f"mkdir -p esmfoldResult/{ID}")
  # prefix = f"esmfoldResult/{ID}/{ID}"
  prefix = f"{store_folder}/{ID}"
  # np.savetxt(f"{prefix}.pae.txt",O["pae"],"%.3f")
  with open(f"{prefix}.pdb","w") as out:
    out.write(pdb_str)
  # return ID

def main(args):
  if args.cpu:
      device = torch.device('cpu')
      deviceC = "cpu"
      device_ids = []
  elif args.gpu:
      if args.device:
          device_ids = [int(id) for id in args.device.split(',')]
          gpus = torch.cuda.is_available()
          if gpus:
              torch.device('cuda')
              # print("length of input sequence must less than ~1270bp ")
              deviceC = "gpu"
          else:
              print("No GPU devices available.")
      else:
        torch.device('cuda')
        deviceC = "gpu"
        device_ids = [int(1)]
        print("Calculated on GPU 1 (default) when using GPU.")
  else:
    print("Use CPU as default.")
    device = torch.device('cpu')
    deviceC = "cpu"
    device_ids = []

            
  fasta_file_name = args.fasta_file
  print("Fasta file input:\n", fasta_file_name)

  store_folder = args.store_folder
  if store_folder.endswith("/"):
    re.sub(r".$", "",store_folder)

  if os.path.exists(store_folder):
    print("output folder:\n", store_folder)
  else:
    os.makedirs(store_folder)
    print("creat output folder:\n", store_folder)

  sequences = []
  jobnames = []

  for record in SeqIO.parse(fasta_file_name, "fasta"):
    jobnames.append(record.id)
    sequences.append(str(record.seq))

  sorted_sequences = sorted(zip(jobnames, sequences), key=lambda x: len(x[1]))

  if args.recycles:
    recycles = args.recycles
  else:
    print("Use num_recycles = 3 as default.")
    recycles = 3

  for jobname, sequence in sorted_sequences:
    runESMFold(jobname, sequence, deviceC, device_ids, store_folder, recycles)


if __name__ == "__main__":

  time_start = time.time()

  parser = argparse.ArgumentParser(description='ESMFold Script')

  parser.add_argument('-i', dest='fasta_file', type=str, help='Path to FASTA file')
  parser.add_argument('-o', dest='store_folder', type=str, help='Path to store file')  
  parser.add_argument('-cpu', action='store_true', help='Run on CPU only')
  parser.add_argument("-gpu", action="store_true", help="Use GPU for computation")
  parser.add_argument("--device", type=str, help="Comma-separated list of GPU device IDs (command: nvidia-smi)")
  parser.add_argument("-recycles", type=int, help="Set the num_recycles (default: 3)")

  args = parser.parse_args()
  main(args)

  time_end = time.time()
  time_sum = time_end - time_start
  minutes = int(time_sum // 60)
  seconds = int(time_sum % 60)
  milliseconds = int((time_sum - int(time_sum)) * 1000)
  print(f"Runing time：{minutes} min {seconds} s {milliseconds} ms")
