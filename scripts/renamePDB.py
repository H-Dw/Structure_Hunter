import argparse
import os


def pdb_modification(pdb, folder, targetfolder):
    pdb_file = os.path.join(folder, pdb)
    with open(pdb_file, "r") as f:
        pdb_lines = f.readlines()
    new_pdb_lines = []
    skip_header = True

    for line in pdb_lines:
        if skip_header and line.startswith("HEADER"):
            header_line = line
            skip_header = False
            new_pdb_lines.append(header_line)  # Retain HEADER row
        elif "PARENT N/A" in line:
            new_pdb_lines.append("HEADER" + "\n")  # Replace PARENT N/A as HEADER
            new_pdb_lines.append("CRYST1" + "\n")
        else:
            new_pdb_lines.append(line)

    output_file = os.path.join(targetfolder, pdb)
    with open(output_file, "w") as f:
        f.writelines(new_pdb_lines)

    # print("PDB file has been modified with PARENT N/A replaced by HEADER and saved")


parser = argparse.ArgumentParser(description='Rename ESMFold output PDB files')
parser.add_argument('-i', dest='input_fold', type=str, help='Path to input fold')
parser.add_argument('-o', dest='output_fold', type=str, help='Path to output fold')
args = parser.parse_args()

PDB_folder = args.input_fold
new_folder = args.output_fold
os.makedirs(new_folder, exist_ok=True)

for file in os.listdir(PDB_folder):
    if file.endswith(".pdb"):
        print(PDB_folder, file)
        pdb_modification(file,PDB_folder, new_folder)
