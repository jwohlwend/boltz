# Prediction

Once you have installed `boltz`, you can start making predictions by simply running:

`boltz predict <INPUT_PATH>`

where `<INPUT_PATH>` is a path to the input file or a directory. The input file can either be in fasta (enough for most use cases) or YAML  format (for more complex inputs). If you specify a directory, `boltz` will run predictions on each `.yaml` or `.fasta` file in the directory.

Before diving into more details about the input formats, here are the key differences in what they each support:

| Feature  | Fasta              | YAML    |
| -------- |--------------------| ------- |
| Polymers | :white_check_mark: | :white_check_mark:   |
| Smiles   | :white_check_mark: | :white_check_mark:   |
| CCD code | :white_check_mark: | :white_check_mark:   |
| Custom MSA | :white_check_mark: | :white_check_mark:   |
| Modified Residues | :x:                |  :white_check_mark: |
| Covalent bonds | :x:                | :white_check_mark:   |
| Pocket conditioning | :x:                | :white_check_mark:   |



## YAML format

The YAML format is more flexible and allows for more complex inputs, particularly around covalent bonds. The schema of the YAML is the following:

```yaml
sequences:
    - ENTITY_TYPE:
        id: CHAIN_ID 
        sequence: SEQUENCE    # only for protein, dna, rna
        smiles: SMILES        # only for ligand, exclusive with ccd
        ccd: CCD              # only for ligand, exclusive with smiles
        conformer: PATH       # optional ligand conformer file (SDF, MOL2 or PDB)
        msa: MSA_PATH         # only for protein
        modifications:
          - position: RES_IDX   # index of residue, starting from 1
            ccd: CCD            # CCD code of the modified residue
        
    - ENTITY_TYPE:
        id: [CHAIN_ID, CHAIN_ID]    # multiple ids in case of multiple identical entities
        ...
constraints:
    - bond:
        atom1: [CHAIN_ID, RES_IDX, ATOM_NAME]
        atom2: [CHAIN_ID, RES_IDX, ATOM_NAME]
    - pocket:
        binder: CHAIN_ID
        contacts: [[CHAIN_ID, RES_IDX], [CHAIN_ID, RES_IDX]]
sample_ligand_conformation: true   # optional, set to false to keep ligand fixed
```
`sequences` has one entry for every unique chain/molecule in the input. Each polymer entity as a `ENTITY_TYPE`  either `protein`, `dna` or `rna` and have a `sequence` attribute. Non-polymer entities are indicated by `ENTITY_TYPE` equal to `ligand` and have a `smiles`, `ccd` or `conformer` attribute. `CHAIN_ID` is the unique identifier for each chain/molecule, and it should be set as a list in case of multiple identical entities in the structure. For proteins, the `msa` key is required by default but can be ommited by passing the `--use_msa_server` flag which will auto-generate the MSA using the mmseqs2 server. If you wish to use a precomputed MSA, use the `msa` attribute with `MSA_PATH` indicating the path to the `.a3m` file containing the MSA for that protein. If you wish to explicitly run single sequence mode (which is generally advised against as it will hurt model performance), you may do so by using the special keyword `empty` for that protein (ex: `msa: empty`).

The `modifications` field is an optional field that allows you to specify modified residues in the polymer (`protein`, `dna` or`rna`). The `position` field specifies the index (starting from 1) of the residue, and `ccd` is the CCD code of the modified residue. This field is currently only supported for CCD ligands.

`constraints` is an optional field that allows you to specify additional information about the input structure. Currently, we support just `bond`. The `bond` constraint specifies a covalent bonds between two atoms (`atom1` and `atom2`). It is currently only supported for CCD ligands and canonical residues, `CHAIN_ID` refers to the id of the residue set above, `RES_IDX` is the index (starting from 1) of the residue (1 for ligands), and `ATOM_NAME` is the standardized atom name (can be verified in CIF file of that component on the RCSB website).

As an example:

```yaml
version: 1
sequences:
  - protein:
      id: [A, B]
      sequence: MVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGVFAALAEAPADSGELARRLDCDARAMRVLLDALYAYDVIDRIHDTNGFRYLLSAEARECLLPGTLFSLVGKFMHDINVAWPAWRNLAEVVRHGARDTSGAESPNGIAQEDYESLVGGINFWAPPIVTTLSRKLRASGRSGDATASVLDVGCGTGLYSQLLLREFPRWTATGLDVERIATLANAQALRLGVEERFATRAGDFWRGGWGTGYDLVLFANIFHLQTPASAVRLMRHAAACLAPDGLVAVVDQIVDADREPKTPQDRFALLFAASMTNTGGGDAYTFQEYEEWFTAAGLQRIETLDTPMHRILLARRATEPSAVPEGQASENLYFQ
      msa: ./examples/msa/seq1.a3m
  - ligand:
      id: [C, D]
      ccd: SAH
  - ligand:
      id: [E, F]
      conformer: ./ligand.sdf
sample_ligand_conformation: false
```


## Fasta format

The fasta format is a little simpler, and should contain entries as follows:

```
>CHAIN_ID|ENTITY_TYPE|MSA_PATH
SEQUENCE
```

The `CHAIN_ID` is a unique identifier for each input chain. The `ENTITY_TYPE` can be one of `protein`, `dna`, `rna`, `smiles`, `ccd` (note that we support both smiles and CCD code for ligands). The `MSA_PATH` is only applicable to proteins. By default, MSA's are required, but they can be ommited by passing the `--use_msa_server` flag which will auto-generate the MSA using the mmseqs2 server. If you wish to use a custom MSA, use it toset path to the `.a3m` file containing a pre-computed MSA for this protein. If you wish to explicitly run single sequence mode (which is generally advised against as it will hurt model performance), you may do so by using the special keyword `empty` for that protein (ex: `>A|protein|empty`).

For each of these cases, the corresponding `SEQUENCE` will contain an amino acid sequence (e.g. `EFKEAFSLF`), a sequence of nucleotide bases (e.g. `ATCG`), a smiles string (e.g. `CC1=CC=CC=C1`), or a CCD code (e.g. `ATP`), depending on the entity.

As an example:

```yaml
>A|protein|./examples/msa/seq1.a3m
MVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGVFAALAEAPADSGELARRLDCDARAMRVLLDALYAYDVIDRIHDTNGFRYLLSAEARECLLPGTLFSLVGKFMHDINVAWPAWRNLAEVVRHGARDTSGAESPNGIAQEDYESLVGGINFWAPPIVTTLSRKLRASGRSGDATASVLDVGCGTGLYSQLLLREFPRWTATGLDVERIATLANAQALRLGVEERFATRAGDFWRGGWGTGYDLVLFANIFHLQTPASAVRLMRHAAACLAPDGLVAVVDQIVDADREPKTPQDRFALLFAASMTNTGGGDAYTFQEYEEWFTAAGLQRIETLDTPMHRILLARRATEPSAVPEGQASENLYFQ
>B|protein|./examples/msa/seq1.a3m
MVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGVFAALAEAPADSGELARRLDCDARAMRVLLDALYAYDVIDRIHDTNGFRYLLSAEARECLLPGTLFSLVGKFMHDINVAWPAWRNLAEVVRHGARDTSGAESPNGIAQEDYESLVGGINFWAPPIVTTLSRKLRASGRSGDATASVLDVGCGTGLYSQLLLREFPRWTATGLDVERIATLANAQALRLGVEERFATRAGDFWRGGWGTGYDLVLFANIFHLQTPASAVRLMRHAAACLAPDGLVAVVDQIVDADREPKTPQDRFALLFAASMTNTGGGDAYTFQEYEEWFTAAGLQRIETLDTPMHRILLARRATEPSAVPEGQASENLYFQ
>C|ccd
SAH
>D|ccd
SAH
>E|smiles
N[C@@H](Cc1ccc(O)cc1)C(=O)O
>F|smiles
N[C@@H](Cc1ccc(O)cc1)C(=O)O
```


## Options

The following options are available for the `predict` command:

    boltz predict [OPTIONS] input_path

| **Option**                  | **Type**        | **Default**        | **Description**                                                                 |
|-----------------------------|-----------------|--------------------|---------------------------------------------------------------------------------|
| `--out_dir PATH`            | `PATH`          | `./`             | The path where to save the predictions.                                         |
| `--cache PATH`              | `PATH`          | `~/.boltz`         | The directory where to download the data and model.                             |
| `--checkpoint PATH`         | `PATH`          | None      | An optional checkpoint. Uses the provided Boltz-1 model by default.             |
| `--devices INTEGER`         | `INTEGER`       | `1`                | The number of devices to use for prediction.                                    |
| `--accelerator`             | `[gpu,cpu,tpu]` | `gpu`              | The accelerator to use for prediction.                                          |
| `--recycling_steps INTEGER` | `INTEGER`       | `3`                | The number of recycling steps to use for prediction.                            |
| `--sampling_steps INTEGER`  | `INTEGER`       | `200`              | The number of sampling steps to use for prediction.                             |
| `--diffusion_samples INTEGER` | `INTEGER`       | `1`                | The number of diffusion samples to use for prediction.                          |
| `--output_format`           | `[pdb,mmcif]`   | `mmcif`            | The output format to use for the predictions.                                   |
| `--num_workers INTEGER`     | `INTEGER`       | `2`                | The number of dataloader workers to use for prediction.                         |
| `--override`                | `FLAG`          | `False`            | Whether to override existing predictions if found.                              |
| `--use_msa_server`                | `FLAG`          | `False`            | Whether to use the msa server to generate msa's.

## Output

After running the model, the generated outputs are organized into the output directory following the structure below:
```
out_dir/
├── lightning_logs/                                            # Logs generated during training or evaluation
├── predictions/                                               # Contains the model's predictions
    ├── [input_file1]/
        ├── [input_file1]_model_0.cif                          # The predicted structure in CIF format
        ...
        └── [input_file1]_model_[diffusion_samples-1].cif      # The predicted structure in CIF format
    └── [input_file2]/
        ...
└── processed/                                                 # Processed data used during execution 
```
The `predictions` folder contains a unique folder for each input file. The input folders contain diffusion_samples predictions saved in the output_format. The `processed` folder contains the processed input files that are used by the model during inference.
