import os
import pandas as pd
import subprocess
import shutil

# Load the dataset from the updated TSV file
file_path = "dataDescription.tsv"
data = pd.read_csv(file_path, sep="\t")

# Create the output directory for renamed files
output_dir = "fastq"
os.makedirs(output_dir, exist_ok=True)

# Iterate through the dataset and process files
for _, row in data.iterrows():
    srr_id = row['SRR']
    heat = row['Heat']
    dtag = row['dTAG']
    nelf_line = row['NELF_line']
    replicate = row['Replicate']

    # Construct the base for output filenames
    base_file_name = f"{nelf_line}.{dtag}.{heat}.{replicate}"

    # Prefetch the SRA file
    print(f"Downloading {srr_id}...")
    subprocess.run(["prefetch", srr_id])

    # Path to the subdirectory created by prefetch
    sra_subdir = os.path.join(os.getcwd(), srr_id)
    sra_file_path = os.path.join(sra_subdir, f"{srr_id}.sra")

    if not os.path.exists(sra_file_path):
        print(f"Error: .sra file {sra_file_path} not found for {srr_id}.")
        continue

    # Convert the .sra file to FASTQ in the current working directory
    print(f"Converting {sra_file_path} to FASTQ...")
    subprocess.run(["fastq-dump", "--gzip", "--split-files", sra_file_path])

    # Locate the downloaded FASTQ files in the current working directory
    read1_file = os.path.join(os.getcwd(), f"{srr_id}_1.fastq.gz")
    read2_file = os.path.join(os.getcwd(), f"{srr_id}_2.fastq.gz")

    # Rename and move the first FASTQ file
    if os.path.exists(read1_file):
        read1_new_path = os.path.join(output_dir, f"{base_file_name}_1.fastq.gz")
        print(f"Renaming {read1_file} to {read1_new_path}...")
        shutil.move(read1_file, read1_new_path)
    else:
        print(f"Error: Read 1 file {read1_file} not found for {srr_id}.")

    # Rename and move the second FASTQ file, if it exists
    if os.path.exists(read2_file):
        read2_new_path = os.path.join(output_dir, f"{base_file_name}_2.fastq.gz")
        print(f"Renaming {read2_file} to {read2_new_path}...")
        shutil.move(read2_file, read2_new_path)
    else:
        print(f"Warning: Read 2 file {read2_file} not found for {srr_id} (single-end data?).")

print("All downloads and renaming completed.")

