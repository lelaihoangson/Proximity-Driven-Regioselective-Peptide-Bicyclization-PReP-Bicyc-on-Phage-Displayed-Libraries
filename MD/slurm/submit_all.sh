#!/bin/bash

# Folder where your SLURM scripts are located
job_folder="selected_2/slurm"

# Submit all .slurm files in the specified folder
for job_file in "$job_folder"/*.slurm; do
    echo "Submitting $job_file ..."
    sbatch "$job_file" &
done

wait
echo "All jobs submitted."