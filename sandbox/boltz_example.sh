#!/bin/sh
  
start_time=$(date +%s)  # Capture the start time

singularity run --nv ./boltz-cuda124-latest.sif  boltz predict --cache /home/ariel/runtime/databases/boltz  --output_format pdb --use_msa_server --diffusion_samples 5 ./examples/ccr1.fasta

end_time=$(date +%s)    # Capture the end time
elapsed=$(( end_time - start_time ))
echo "Elapsed time: ${elapsed} seconds"
