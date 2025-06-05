<div align="center">
  <div>&nbsp;</div>
  <img src="docs/boltz_title.png" width="400"/>

[Paper](https://doi.org/10.1101/2024.11.19.624167) |
[Slack](https://join.slack.com/t/boltz-community/shared_invite/zt-34qg8uink-V1LGdRRUf3avAUVaRvv93w) <br> <br>
</div>


![](docs/boltz1_pred_figure.png)


## Introduction

Boltz-1 and Boltz-1x are state-of-the-art open-source models to predict biomolecular structures containing combinations of proteins, RNA, DNA, and other molecules. They also support modified residues, covalent ligands and glycans, as well as conditioning the prediction on specified interaction pockets or contacts. 

All the code and weights are provided under MIT license, making them freely available for both academic and commercial uses. For more information about the model, see our [technical report](https://doi.org/10.1101/2024.11.19.624167). To discuss updates, tools and applications join our [Slack channel](https://join.slack.com/t/boltz-community/shared_invite/zt-34qg8uink-V1LGdRRUf3avAUVaRvv93w).

## Installation
Install boltz with PyPI (recommended):

```
pip install boltz -U
```

or directly from GitHub for daily updates:

```
git clone https://github.com/jwohlwend/boltz.git
cd boltz; pip install -e .
```
> Note: we recommend installing boltz in a fresh python environment

## Inference

You can run inference using Boltz-1x with:

```
boltz predict input_path --use_msa_server
```

Boltz currently accepts three input formats:

1. Fasta file, for most use cases

2. A comprehensive YAML schema, for more complex use cases

3. A directory containing files of the above formats, for batched processing

To see all available options: `boltz predict --help` and for more information on these input formats, see our [prediction instructions](docs/prediction.md).

## Evaluation

To encourage reproducibility and facilitate comparison with other models, we provide the evaluation scripts and predictions for Boltz-1, Boltz-1x, Chai-1 and AlphaFold3 on our test benchmark dataset as well as CASP15. These datasets are created to contain biomolecules different from the training data and to benchmark the performance of these models we run them with the same input MSAs and same number  of recycling and diffusion steps. More details on these evaluations can be found in our [evaluation instructions](docs/evaluation.md).

![Test set evaluations](docs/plot_test.png)
![CASP15 set evaluations](docs/plot_casp.png)


## Training

If you're interested in retraining the model, see our [training instructions](docs/training.md).

## Contributing

We welcome external contributions and are eager to engage with the community. Connect with us on our [Slack channel](https://join.slack.com/t/boltz-community/shared_invite/zt-34qg8uink-V1LGdRRUf3avAUVaRvv93w) to discuss advancements, share insights, and foster collaboration around Boltz-1.

Boltz also runs on Tenstorrent hardware thanks to a [fork](https://github.com/moritztng/tt-boltz) by Moritz Thüning.

## License

Our model and code are released under MIT License, and can be freely used for both academic and commercial purposes.


## Cite

If you use this code or the models in your research, please cite the following paper:

```bibtex
@article{wohlwend2024boltz1,
  author = {Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and Getz, Noah and Reveiz, Mateo and Leidal, Ken and Swiderski, Wojtek and Atkinson, Liam and Portnoi, Tally and Chinn, Itamar and Silterra, Jacob and Jaakkola, Tommi and Barzilay, Regina},
  title = {Boltz-1: Democratizing Biomolecular Interaction Modeling},
  year = {2024},
  doi = {10.1101/2024.11.19.624167},
  journal = {bioRxiv}
}
```

In addition if you use the automatic MSA generation, please cite:

```bibtex
@article{mirdita2022colabfold,
  title={ColabFold: making protein folding accessible to all},
  author={Mirdita, Milot and Sch{\"u}tze, Konstantin and Moriwaki, Yoshitaka and Heo, Lim and Ovchinnikov, Sergey and Steinegger, Martin},
  journal={Nature methods},
  year={2022},
}
```
