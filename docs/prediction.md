# Prediction

Once you have installed `boltz`, you can start making predictions by simply running:

`boltz predict <INPUT_PATH> --use_msa_server`

where `<INPUT_PATH>` is a path to the input file or a directory. The input file can either be in fasta (enough for most use cases) or YAML  format (for more complex inputs). If you specify a directory, `boltz` will run predictions on each `.yaml` or `.fasta` file in the directory. Passing the `--use_msa_server` flag will auto-generate the MSA using the mmseqs2 server, otherwise you can provide a precomputed MSA. 

By default the package will run the Boltz-1x model that includes the inference time potentials that significantly improve the physical quality of the poses. If you find any physical issues with the model predictions, please let us know by opening an issue and including the YAML/FASTA file to replicate, the structure output and a description of the problem. If you want to run the original Boltz-1 model without potentials you can add `--no_potentials` flag. 

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
        smiles: 'SMILES'        # only for ligand, exclusive with ccd
        ccd: CCD              # only for ligand, exclusive with smiles
        msa: MSA_PATH         # only for protein
        modifications:
          - position: RES_IDX   # index of residue, starting from 1
            ccd: CCD            # CCD code of the modified residue
        cyclic: false
     
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
```
`sequences` has one entry for every unique chain/molecule in the input. Each polymer entity as a `ENTITY_TYPE`  either `protein`, `dna` or `rna` and have a `sequence` attribute. Non-polymer entities are indicated by `ENTITY_TYPE` equal to `ligand` and have a `smiles` or `ccd` attribute. `CHAIN_ID` is the unique identifier for each chain/molecule, and it should be set as a list in case of multiple identical entities in the structure. For proteins, the `msa` key is required by default but can be omited by passing the `--use_msa_server` flag which will auto-generate the MSA using the mmseqs2 server. If you wish to use a precomputed MSA, use the `msa` attribute with `MSA_PATH` indicating the path to the `.a3m` file containing the MSA for that protein. If you wish to explicitly run single sequence mode (which is generally advised against as it will hurt model performance), you may do so by using the special keyword `empty` for that protein (ex: `msa: empty`). For custom MSA, you may wish to indicate pairing keys to the model. You can do so by using a CSV format instead of a3m with two columns: `sequence` with the protein sequences and `key` which is a unique identifier indicating matching rows across CSV files of each protein chain.

The `modifications` field is an optional field that allows you to specify modified residues in the polymer (`protein`, `dna` or`rna`). The `position` field specifies the index (starting from 1) of the residue, and `ccd` is the CCD code of the modified residue. This field is currently only supported for CCD ligands. The `cyclic` flag should be used to specify polymer chains (not ligands) that are cyclic. 

`constraints` is an optional field that allows you to specify additional information about the input structure. 


* The `bond` constraint specifies covalent bonds between two atoms (`atom1` and `atom2`). It is currently only supported for CCD ligands and canonical residues, `CHAIN_ID` refers to the id of the residue set above, `RES_IDX` is the index (starting from 1) of the residue (1 for ligands), and `ATOM_NAME` is the standardized atom name (can be verified in CIF file of that component on the RCSB website).

* The `pocket` constraint specifies the residues associated with a ligand, where `binder` refers to the chain binding to the pocket (which can be a molecule, protein, DNA or RNA) and `contacts` is the list of chain and residue indices (starting from 1) associated with the pocket. The model currently only supports the specification of a single `binder` chain (and any number of `contacts` residues in other chains).

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
      smiles: 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'
```


## Fasta format

The fasta format is a little simpler, and should contain entries as follows:

```
>CHAIN_ID|ENTITY_TYPE|MSA_PATH
SEQUENCE
```

The `CHAIN_ID` is a unique identifier for each input chain. The `ENTITY_TYPE` can be one of `protein`, `dna`, `rna`, `smiles`, `ccd` (note that we support both smiles and CCD code for ligands). The `MSA_PATH` is only applicable to proteins. By default, MSA's are required, but they can be omited by passing the `--use_msa_server` flag which will auto-generate the MSA using the mmseqs2 server. If you wish to use a custom MSA, use it to set the path to the `.a3m` file containing a pre-computed MSA for this protein. If you wish to explicitly run single sequence mode (which is generally advised against as it will hurt model performance), you may do so by using the special keyword `empty` for that protein (ex: `>A|protein|empty`). For custom MSA, you may wish to indicate pairing keys to the model. You can do so by using a CSV format instead of a3m with two columns: `sequence` with the protein sequences and `key` which is a unique identifier indicating matching rows across CSV files of each protein chain.

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

    boltz predict input_path [OPTIONS]

As an example, to predict a structure using 10 recycling steps and 25 samples (the default parameters for AlphaFold3) use:

    boltz predict input_path --recycling_steps 10 --diffusion_samples 25

(note however that the prediction will take significantly longer)


| **Option**               | **Type**        | **Default**                 | **Description**                                                                                                                                                                     |
|--------------------------|-----------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--out_dir`              | `PATH`          | `./`                        | The path where to save the predictions.                                                                                                                                             |
| `--cache`                | `PATH`          | `~/.boltz`                  | The directory where to download the data and model. Will use environmnet variable `BOLTZ_CACHE` as an absolute path if set                                                          |
| `--checkpoint`           | `PATH`          | None                        | An optional checkpoint. Uses the provided Boltz-1 model by default.                                                                                                                 |
| `--devices`              | `INTEGER`       | `1`                         | The number of devices to use for prediction.                                                                                                                                        |
| `--accelerator`          | `[gpu,cpu,tpu]` | `gpu`                       | The accelerator to use for prediction.                                                                                                                                              |
| `--recycling_steps`      | `INTEGER`       | `3`                         | The number of recycling steps to use for prediction.                                                                                                                                |
| `--sampling_steps`       | `INTEGER`       | `200`                       | The number of sampling steps to use for prediction.                                                                                                                                 |
| `--diffusion_samples`    | `INTEGER`       | `1`                         | The number of diffusion samples to use for prediction.                                                                                                                              |
| `--step_scale`           | `FLOAT`         | `1.638`                     | The step size is related to the temperature at which the diffusion process samples the distribution. The lower the higher the diversity among samples (recommended between 1 and 2). |
| `--output_format`        | `[pdb,mmcif]`   | `mmcif`                     | The output format to use for the predictions.                                                                                                                                       |
| `--num_workers`          | `INTEGER`       | `2`                         | The number of dataloader workers to use for prediction.                                                                                                                             |
| `--override`             | `FLAG`          | `False`                     | Whether to override existing predictions if found.                                                                                                                                  |
| `--use_msa_server`       | `FLAG`          | `False`                     | Whether to use the msa server to generate msa's.                                                                                                                                    |
| `--msa_server_url`       | str             | `https://api.colabfold.com` | MSA server url. Used only if --use_msa_server is set.                                                                                                                               |
| `--msa_pairing_strategy` | str             | `greedy`                    | Pairing strategy to use. Used only if --use_msa_server is set. Options are 'greedy' and 'complete'                                                                                  |
| `--no_potentials`        | `FLAG`          | `False`                     | Whether to run the original Boltz-1 model without inference time potentials.                                                                                                        |
| `--write_full_pae`       | `FLAG`          | `False`                     | Whether to save the full PAE matrix as a file.                                                                                                                                      |
| `--write_full_pde`       | `FLAG`          | `False`                     | Whether to save the full PDE matrix as a file.                                                                                                                                      |

## Output

After running the model, the generated outputs are organized into the output directory following the structure below:
```
out_dir/
├── lightning_logs/                                            # Logs generated during training or evaluation
├── predictions/                                               # Contains the model's predictions
    ├── [input_file1]/
        ├── [input_file1]_model_0.cif                          # The predicted structure in CIF format, with the inclusion of per token pLDDT scores
        ├── confidence_[input_file1]_model_0.json              # The confidence scores (confidence_score, ptm, iptm, ligand_iptm, protein_iptm, complex_plddt, complex_iplddt, chains_ptm, pair_chains_iptm)
        ├── pae_[input_file1]_model_0.npz                      # The predicted PAE score for every pair of tokens
        ├── pde_[input_file1]_model_0.npz                      # The predicted PDE score for every pair of tokens
        ├── plddt_[input_file1]_model_0.npz                    # The predicted pLDDT score for every token
        ...
        └── [input_file1]_model_[diffusion_samples-1].cif      # The predicted structure in CIF format
        ...
    └── [input_file2]/
        ...
└── processed/                                                 # Processed data used during execution 
```
The `predictions` folder contains a unique folder for each input file. The input folders contain `diffusion_samples` predictions saved in the output_format ordered by confidence score as well as additional files containing the predictions of the confidence model. The `processed` folder contains the processed input files that are used by the model during inference.

The output `.json` file contains various aggregated confidence scores for specific sample. The structure of the file is as follows:
```yaml
{
    "confidence_score": 0.8367,       # Aggregated score used to sort the predictions, corresponds to 0.8 * complex_plddt + 0.2 * iptm (ptm for single chains)
    "ptm": 0.8425,                    # Predicted TM score for the complex
    "iptm": 0.8225,                   # Predicted TM score when aggregating at the interfaces
    "ligand_iptm": 0.0,               # ipTM but only aggregating at protein-ligand interfaces
    "protein_iptm": 0.8225,           # ipTM but only aggregating at protein-protein interfaces
    "complex_plddt": 0.8402,          # Average pLDDT score for the complex
    "complex_iplddt": 0.8241,         # Average pLDDT score when upweighting interface tokens
    "complex_pde": 0.8912,            # Average PDE score for the complex
    "complex_ipde": 5.1650,           # Average PDE score when aggregating at interfaces  
    "chains_ptm": {                   # Predicted TM score within each chain
        "0": 0.8533,
        "1": 0.8330
    },
    "pair_chains_iptm": {             # Predicted (interface) TM score between each pair of chains
        "0": {
            "0": 0.8533,
            "1": 0.8090
        },
        "1": {
            "0": 0.8225,
            "1": 0.8330
        }
    }
}
```
`confidence_score`, `ptm` and `plddt` scores (and their interface and individual chain analogues) have a range of [0, 1], where higher values indicate higher confidence. `pde` scores have a unit of angstroms, where lower values indicate higher confidence.
