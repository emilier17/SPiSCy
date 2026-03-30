#==================================================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: differential analysis
# Step: prepare input data for differential analysis
# 
# Author: Émilie Roy
# Date: Feb 2026
# Version: v.1.0
#
# Input : final_samples.csv (marker data)
# Outputs : marker data split into individual csv files according to filename
#==================================================================================================================#


from clustering_utils import (
    snakemake_logs, start_time, end_time, elapsed_time,
    )


snakemake_logs(snakemake)
start = start_time()


# Library
import pandas as pd
import os

# Snakemake variables
# Inputs
samples_csv = snakemake.input[0]
# Outputs
output_dir = snakemake.output[0]
os.makedirs(output_dir, exist_ok=True)

# Adjust depending on memory
CHUNK_SIZE = 1_000_000

# Track which files already have headers written
written_files = set()

# read samples csv chunk by chunk
reader = pd.read_csv(samples_csv, chunksize=CHUNK_SIZE)


for i, chunk in enumerate(reader):

    print(f"Processing chunk #{i}...")

    # Get all filenames present in the chunk
    filenames = chunk["filename"].unique()

    # Iterate over all filenames in the chunk
    for file in filenames:

        # create a mask that indicates if a row is part of the file
        mask = chunk["filename"] == file
        subfile = chunk.loc[mask]

        # Remove the "filename" column
        subfile = subfile.drop(columns=["filename"])

        outpath = os.path.join(output_dir, f"{file}.csv")

        write_header = file not in written_files

        subfile.to_csv(
            outpath,
            mode="w" if write_header else "a",
            header=write_header,
            index=False
        )

        written_files.add(file)


end = end_time()
elapsed_time(start, end)

