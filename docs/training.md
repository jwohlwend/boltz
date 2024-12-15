# Training

## Download the pre-processed data

Instructions on how to download the processed dataset for training are coming soon, we are currently uploading the data to sharable storage and will update this page when ready.

If you wish to re-run the preprocessing pipeline, or processed your own raw data for training / fine-tuning, see the instructions at the bottom of this page.

## Modify the configuration file

The training script requires a configuration file to run. This file specifies the paths to the data, the output directory, and other parameters of the data, model and training process. 

We provide under `scripts/train/configs` a template configuration file analogous to the one we used for training the structure model (`structure.yaml`) and the confidence model (`confidence.yaml`).

The following are the main parameters that you should modify in the configuration file to get the structure model to train:

```yaml
trainer:
  devices: 1

output: SET_PATH_HERE                 # Path to the output directory  
resume: PATH_TO_CHECKPOINT_FILE       # Path to a checkpoint file to resume training from if any null otherwise

data:
  datasets:
    - _target_: boltz.data.module.training.DatasetConfig
      target_dir: PATH_TO_TARGETS_DIR       # Path to the directory containing the processed structure files
      msa_dir: PATH_TO_MSA_DIR              # Path to the directory containing the processed MSA files

  symmetries: PATH_TO_SYMMETRY_FILE      # Path to the file containing molecule the symmetry information
  max_tokens: 512                        # Maximum number of tokens in the input sequence
  max_atoms: 4608                        # Maximum number of atoms in the input structure
```

`max_tokens` and `max_atoms` are the maximum number of tokens and atoms in the crop. Depending on the size of the GPUs you are using (as well as the training speed desired), you may want to adjust these values. Other recommended values are 256 and 2304, or 384 and 3456 respectively.

## Run the training script

Before running the full training, we recommend using the debug flag. This turns off DDP (sets single device) and sets `num_workers` to 0 so everything is in a single process, as well as disabling wandb:

    python scripts/train/train.py scripts/train/configs/structure.yaml debug=1

Once that seems to run okay, you can kill it and launch the training run:

    python scripts/train/train.py scripts/train/configs/structure.yaml

We also provide a different configuration file to train the confidence model:

    python scripts/train/train.py scripts/train/configs/confidence.yaml


## Processing raw data

There are several steps to perform when processing raw data for boltz training.


#### Step 1: Go to the processing folder

```bash
cd scripts/process
```

#### Step 2: Install requirements

Install the few extra requirements required for processing:

```bash
pip install -r requirements.txt
```

You must also install two external libraries: `mmseqs` and `redis`. Instructions for installation are below:

- `mmseqs`: https://github.com/soedinglab/mmseqs2?tab=readme-ov-file#installation
- `redis`: https://redis.io/docs/latest/operate/oss_and_stack/install/install-redis/

#### Step 3: Preprocess the CCD dictionary


We have already done this for you, the relevant file is here. Unless you wish to do it again yourself, you can skip to the next step!

If you do want to do this your self, you can do so with the following commands:

```bash
wget https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
python ccd.py --components components.cif --ouptut ccd.pkl
```

> Note: runs in parallel by default with as many threads as cpu cores on your machine, can be changed with `--num_processed`

#### Step 4: Create sequence clusters

First, you must create a fasta file containing all the polymer sequences present in your data. You can use any header format you want for the sequences, it will not be used.

For the PDB, this can already be downloaded here:
```bash
wget https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gunzip -d pdb_seqres.txt.gz
```

For the OpenFold data, we provide this for you here:
```
Link to come
``` 

When this is done, you can run the clustering script, which assigns proteins to 40% similarity clusters and rna/dna to a cluster for each unique sequence. For ligands, each CCD code is also assigned to its own cluster.

```bash
python cluster.py --ccd ccd.pkl --sequences pdb_seqres.txt --mmseqs PATH_TO_MMSEQS_EXECUTABLE
```

> Note: you must install mmseqs (see: https://github.com/soedinglab/mmseqs2?tab=readme-ov-file#installation)

#### Step 5: Create MSA's

We have already computed MSA's for all sequences in the PDB at the time of training using the ColabFold `colab_search` tool. You can setup your own local colabfold using instructions provided here: https://github.com/YoshitakaMo/localcolabfold

The raw MSA's for the PDB can be found here:
```
Link to come
```

You can also download the raw OpenFold MSA's here:
```
Link to come
```

If you wish to use your own MSA's, just ensure that their file name is the hash of the query sequence, according to the following function:
```python
import hashlib

def hash_sequence(seq: str) -> str:
    """Hash a sequence."""
    return hashlib.sha256(seq.encode()).hexdigest()
```

#### Step 6: Process MSA's

During MSA processing, we annotate sequences using their taxonomy ID, which is important for MSA pairing during training.

This happens only on MSA sequences with headers that start with the following:

```
>UniRef100_UNIREFID
...
```

This format is the way that MSA's are provided by colabfold. If you use a different MSA pipeline, make sure your Uniref MSA's follow the above format.

Next you should download our provided taxonomy database and place it in this folder:

```bash
Link to come
```

You can now process the raw MSAs. First launch a redis server. We use redis to share the large taxonomy dictionary across workers, so MSA processing can happen in parallel without blowing up the RAM usage.

```bash
redis-server --dbfilename taxonomy.rdb --redis-port 7777
```
> Note: You must have redis installed (see: https://redis.io/docs/latest/operate/oss_and_stack/install/install-redis/)

In a separate shell, run the processing script:
```bash
python msa.py --msadir YOUR_MSA_DIR --outdir YOUR_OUTPUT_DIR --redis-port 7777
```

> Important: the script looks for `.a3m` files in the directory, make sure to match this extension and file format.

#### Step 7: Process structures

Finally, we're ready to process structural data. Here we provide two different scripts for the PDB and for the OpenFold data. In general, we recommend using the `rcsb.py` script for your own data, which is expected in `mmcif` format.

You can download the full RCSB using the instructions here:
https://www.rcsb.org/docs/programmatic-access/file-download-services


```bash
wget Link to come
redis-server --dbfilename ccd.rdb --redis-port 7777
```
> Note: You must have redis installed (see: https://redis.io/docs/latest/operate/oss_and_stack/install/install-redis/)

In a separate shell, run the processing script, make sure to use the `clustering.json` file you previously created.
```bash
python rcsb.py --datadir PATH_TO_MMCIF_DIR --cluster clustering.json --outdir YOUR_OUTPUT_DIR --use-assembly --max-file-size 7000000 --redis-port 7777
```

> Important: the script looks for `.cif` or `cif.gz` files in the directory, make sure to match this extension and file format.

> We skip a few of the very large files, you can modify this using the `--max-file-size` flag, or by removing it.

#### Step 8: Ready!

You're ready to start training the model on your data, make sure to modify the config to assign the paths you created in the previous two steps. If you have any questions, don't hesitate to open an issue or reach out on our community slack channel.