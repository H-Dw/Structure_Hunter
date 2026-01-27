
# Structure Hunter

This pipeline is designed to identify and collect specific-conformation proteins from proteomic data. The workflow integrates signal peptide prediction, structural prediction (ESMFold), and structural similarity search (FoldSeek) to identify specific targets based on secondary structure templates.

NOTE: This pipeline could be applied in other structure-specific protein identification tasks, and follow example shows the Slam-related C-terminal Beta-barrel proteins identification.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Pipeline Overview](#pipeline-overview)
- [Usage](#usage)
  - [1. Signal Peptide Identification](#1-signal-peptide-identification)
  - [2. Structure Prediction (ESMFold)](#2-structure-prediction-esmfold)
  - [3. Structural Search (FoldSeek)](#3-structural-search-foldseek)
  - [4. C-Terminal Filtering](#4-c-terminal-filtering)
- [Scripts Description](#scripts-description)

## Prerequisites

Ensure the following tools and libraries are installed and accessible in your environment:

*   **Linux Environment**
*   **[SignalP 5.0](https://services.healthtech.dtu.dk/service.php?SignalP-5.0)** (Gram-negative organism support)
*   **[SeqKit](https://github.com/shenwei356/seqkit)**
*   **[ESMFold](https://github.com/facebookresearch/esm)** (Requires GPU for efficient processing)
*   **[FoldSeek](https://github.com/steineggerlab/foldseek)**
*   **[DSSP](https://github.com/cmbi/dssp)** (For secondary structure assignment)

---

## Pipeline Overview

1.  **Signal Peptide Screening:** Filter raw proteomes for proteins containing signal peptides suitable for Gram-negative bacteria.
2.  **Sequence Partitioning:** Split sequences by length (current threshold: 1250bp) to optimize computational resources (GPU for short, CPU for long).
3.  **Structure Prediction:** Predict 3D structures using ESMFold.
4.  **Database Creation:** specific structure databases for the target organisms.
5.  **Structural Search:** Query these databases against known Slam-related Beta-barrel templates (TbpB, LbpB, fHbp, HpuA) using FoldSeek.
6.  **Target Filtering:** Filter results to retain only proteins where the match is located at the C-terminus with high confidence.

NOTE: Slam-related proteins refered to publication: [Hooda Y, Lai C C L, Moraes T F. Identification of a Large Family of Slam-Dependent Surface Lipoproteins in Gram-Negative Bacteria [J]. Frontiers in Cellular and Infection Microbiology, 2017, 7.](https://doi.org/10.3389/fcimb.2017.00207).

---

## Installation

Clone the repository and install the required dependencies:
Due to license restrictions, this we cannot provide signalp5 directly. Please download from https://services.healthtech.dtu.dk/service.php?SignalP-5.0.

```bash
git clone https://github.com/H-Dw/Structure_Hunter.git
cd Structure_Hunter
```

```bash
conda env create -f environment.yml
# NOTE: The SignalP need to be installed manually
```

OR you can install following dependencies directly:
```bash
pip install fair-esm
conda install -c conda-forge -c bioconda foldseek seqkit dssp pandas
```

Install SignalP v5.0b after downloading signalp-5.0b.Linux.tar.gz:
```bash
tar -vxf signalp-5.0b.Linux.tar.gz
```
Then, you can write your install path into the ~/.bashrc, so that system can identify signalp directly. 
```bash
vim ~/.bashrc

# Add follow command and save using vim
export PATH="$PATH:/path/to/signalp-5.0b/bin" # /path/to/: Install path of signalp5

# Reactivate system
source ~/.bashrc
```
OR you can modify the running command by adding 'bin' path of signalp, likes:
```bash
/path/to/signalp-5.0b/bin/signalp5 ...
```

## Usage

We used NgSignPDB as an example (processed by ESMFold prediction), the same process can be applied to EcSignPDB, McSignPDB, and NmSignPDB.

### 1. Signal Peptide Identification

Run `SignalP 5.0` to identify proteins with signal peptides, extract their IDs, and separate sequences based on length (1250bp threshold).

```bash
export species="Ng"

# NOTE: 'protein.faa' is the input proteome downloaded from NCBI, it need to be replaced by your own proteome path

# 1. Run SignalP (Gram-negative mode)
signalp -fasta examples/protein.faa -format short -org gram- -prefix ${species}Sign

# 2. Extract IDs with detected signal peptides
awk -F "\t" '$7 ~ /[0-9]+/ { print $1 }' examples/${species}Sign_summary.signalp5 > examples/${species}Sign.txt

# 3. Extract FASTA sequences based on IDs
seqkit grep -f ${species}Sign.txt examples/protein.faa -o examples/${species}Sign.pep.fasta

# 4. Check sequence statistics
cat examples/${species}Sign.pep.fasta | seqkit seq -m 1250 | seqkit stats

# 5. Split sequences by length for resource optimization (optional)
# > 1250bp (Long sequences)
seqkit seq -m 1250 -g examples/${species}Sign.pep.fasta > examples/${species}SignLL.pep.fasta
# < 1250bp (Short sequences)
seqkit seq -M 1250 -g examples/${species}Sign.pep.fasta > examples/${species}Sign1250.pep.fa

```

### 2. Structure Prediction (ESMFold)

Predict structures using ESMFold. Utilize GPU for shorter sequences and CPU for longer ones to avoid memory overflow.

**On GPU Server (Length < 1250bp):**
```bash
# Could add '--device 2' if multiple GPUs are available
nohup python scripts/runESMFold.py -i examples/${species}Sign1250.pep.fasta -o examples/${species}ESMFoldResult -gpu &
```

**On CPU Server (Length > 1250bp):**
```bash
nohup python scripts/runESMFold.py -i examples/${species}LL.pep.fasta -o examples/${species}ESMFoldResult -cpu &
```

**Post-Processing (DSSP):**
Rename PDB files and extract secondary structures.
```bash
python scripts/renamePDB.py -i examples/${species}ESMFoldResult/ -o examples/${species}SignPDB/pdb/
bash scripts/runDSSP.sh examples/${species}SignPDB/pdb/ examples/${species}SignPDB
```

### 3. Structural Search (FoldSeek)

Create FoldSeek databases for the predicted structures and search against the Beta-barrel template library.

**Templates:** *Nm* TbpB, LbpB, fHbp, and HpuA (Derived from DSSP database).

```bash
# 1. Create FoldSeek Databases
mkdir -p examples/foldseekDB
foldseek createdb examples/${species}SignPDB/pdb/ examples/foldseekDB/${species}SignDB

# 2. Run Easy-Search
# Flags: --cov-mode 2 (Query coverage), --format-mode 4 (Custom output)
mkdir -p examples/foldseek_result
foldseek easy-search lipoprotein_template/beta_barrel/ examples/foldseekDB/${species}SignDB examples/foldseek_result/${species}_lipoprotein_beta_barrel.tsv tmp --remove-tmp-files 1 --cov-mode 2 --format-mode 4 --format-output 'query,target,fident,alnlen,nident,mismatch,gapopen,qstart,qend,tstart,tend,evalue,qtmscore,bits,tlen,qaln,taln,tseq'
```

### 4. C-Terminal Filtering

Filter the FoldSeek results to identify "C-terminal" matches.

**Filtering Rules:**
1.  **Position:** The midpoint of the target alignment must be in the last 50% of the target sequence length.
    *   Equation: `(tstart + tend) / 2 > tlen / 2`
2.  **Significance:** `e-value ≤ 0.05`
3.  **Deduplication:** Sort by `bits` score and remove duplicates based on the target column.

```bash
python scripts/extract_cterminal_target.py --input examples/foldseek_result/${species}_lipoprotein_beta_barrel.tsv --output examples/foldseek_result/${species}_lipoprotein_cterminal.tsv
```
