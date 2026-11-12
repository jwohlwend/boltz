"""Contracts between clustering, mmCIF preprocessing and training input."""

import importlib
import io
import json
import pickle
from argparse import Namespace
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import gemmi
import numpy as np
import pytest
from rdkit import Chem

from boltz.data import const
from boltz.data.module.training import load_input
from boltz.data.parse.a3m import _parse_a3m
from boltz.data.types import (
    Atom,
    Bond,
    Chain,
    Connection,
    Interface,
    Record,
    Residue,
    Structure,
    StructureInfo,
)


@pytest.fixture
def scripts(monkeypatch):
    pytest.importorskip("redis", reason="Install scripts/process/requirements.txt")
    pytest.importorskip("p_tqdm", reason="Install scripts/process/requirements.txt")
    monkeypatch.syspath_prepend(
        str(Path(__file__).resolve().parents[2] / "scripts/process")
    )
    return importlib.import_module("rcsb"), importlib.import_module("cluster")


@pytest.fixture
def parsed():
    chains = np.zeros(6, dtype=Chain)
    residues = np.zeros(3, dtype=Residue)
    residues["name"] = ["ATP", "NAG", "MAN"]
    sequences = {0: "ACDEFGHIKLMN", 1: "ACDEFGHIKLMN", 2: "AGCT", 3: "ACGU"}
    for i, (name, kind, entity) in enumerate(
        [
            ("A", "PROTEIN", 7),
            ("B", "PROTEIN", 7),
            ("C", "DNA", 8),
            ("D", "RNA", 9),
            ("L", "NONPOLYMER", 10),
            ("M", "NONPOLYMER", 11),
        ]
    ):
        chains[i]["name"] = name
        chains[i]["mol_type"] = const.chain_type_ids[kind]
        chains[i]["entity_id"] = entity
        chains[i]["asym_id"] = i
        chains[i]["res_idx"] = 0 if i < 5 else 1
        chains[i]["res_num"] = (
            len(sequences[i]) if i in sequences else (1 if i == 4 else 2)
        )
    structure = Structure(
        atoms=np.zeros(0, dtype=Atom),
        bonds=np.zeros(0, dtype=Bond),
        residues=residues,
        chains=chains,
        connections=np.zeros(0, dtype=Connection),
        interfaces=np.array([(0, 4), (2, 3)], dtype=Interface),
        mask=np.ones(6, dtype=bool),
    )
    return SimpleNamespace(
        data=structure, info=StructureInfo(num_chains=6), sequences=sequences
    )


def test_cluster_script_output_is_consumed_by_rcsb(
    scripts, parsed, tmp_path, monkeypatch
):
    rcsb, cluster = scripts
    fasta = tmp_path / "seqres.fasta"
    fasta.write_text(
        "".join(f">{name}\n{seq}\n" for name, seq in parsed.sequences.items())
    )
    ccd = tmp_path / "ccd.pkl"
    ccd.write_bytes(pickle.dumps({"ATP": None, "NAG": None, "MAN": None}))
    outdir = tmp_path / "clustering"

    def mmseqs(command, **kwargs):
        proteins = (outdir / "proteins.fasta").read_text().splitlines()
        ids = [line[1:] for line in proteins if line.startswith(">")]
        (outdir / "clust_prot_cluster.tsv").write_text(
            "".join(f"representative\t{seq_id}\n" for seq_id in ids)
        )

    monkeypatch.setattr(cluster.subprocess, "run", mmseqs)
    cluster.main(
        Namespace(
            outdir=str(outdir), sequences=str(fasta), ccd=str(ccd), mmseqs="mmseqs"
        )
    )
    clusters = json.loads((outdir / "clustering.json").read_text())
    monkeypatch.setattr(rcsb, "parse_mmcif", lambda *args: parsed)
    target = rcsb.parse(rcsb.PDB("EXAMPLE", "unused.cif"), {}, clusters)
    assert [c.cluster_id for c in target.record.chains] == [
        "representative",
        "representative",
        cluster.hash_sequence("AGCT"),
        cluster.hash_sequence("ACGU"),
        "ATP",
        -1,
    ]
    assert [c.entity_id for c in target.record.chains] == [7, 7, 8, 9, 10, 11]
    assert [(i.chain_1, i.chain_2) for i in target.record.interfaces] == [
        (0, 4),
        (2, 3),
    ]


@pytest.mark.parametrize("available", [False, True])
def test_only_existing_protein_msas_are_linked_and_loadable(
    scripts, parsed, tmp_path, monkeypatch, available
):
    rcsb, cluster = scripts
    monkeypatch.setattr(rcsb, "parse_mmcif", lambda *args: parsed)
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    seq_id = cluster.hash_sequence(parsed.sequences[0])
    if available:
        msa = _parse_a3m(io.StringIO(">query\nACDEFGHIKLMN\n"), taxonomy=None)
        np.savez_compressed(msa_dir / f"{seq_id}.npz", **asdict(msa))
    # Even similarly named nucleotide files must not become protein MSAs.
    for seq in ["AGCT", "ACGU"]:
        (msa_dir / f"{cluster.hash_sequence(seq)}.npz").write_bytes(
            b"not a protein MSA"
        )
    target = rcsb.parse(rcsb.PDB("example", "unused.cif"), {}, {}, msa_dir=msa_dir)
    assert [c.msa_id for c in target.record.chains] == [
        seq_id if available else ""
    ] * 2 + [""] * 4
    if available:
        structures = tmp_path / "structures"
        structures.mkdir()
        np.savez_compressed(structures / "example.npz", **asdict(target.structure))
        loaded = load_input(target.record, tmp_path, msa_dir)
        assert set(loaded.msa) == {0, 1}
        assert np.array_equal(loaded.msa[0].residues, msa.residues)


def test_legacy_mapping_and_missing_clusters_remain_supported(
    scripts, parsed, monkeypatch
):
    rcsb, cluster = scripts
    monkeypatch.setattr(rcsb, "parse_mmcif", lambda *args: parsed)
    clusters = {
        "example_7": "old",
        cluster.hash_sequence(parsed.sequences[2]): "dna",
        "atp": "ligand",
    }
    target = rcsb.parse(rcsb.PDB("EXAMPLE", "unused.cif"), {}, clusters)
    assert [c.cluster_id for c in target.record.chains] == [
        "old",
        "old",
        "dna",
        -1,
        "ligand",
        -1,
    ]
    assert all(c.msa_id == "" for c in target.record.chains)


def test_sequence_clusters_take_precedence_over_legacy_keys(
    scripts, parsed, monkeypatch
):
    rcsb, cluster = scripts
    monkeypatch.setattr(rcsb, "parse_mmcif", lambda *args: parsed)
    key = cluster.hash_sequence(parsed.sequences[0])
    target = rcsb.parse(
        rcsb.PDB("example", "unused.cif"), {}, {key: "new", "example_7": "old"}
    )
    assert [c.cluster_id for c in target.record.chains[:2]] == ["new", "new"]


def test_legacy_only_cluster_mapping_control(scripts, parsed, monkeypatch):
    rcsb, _ = scripts
    monkeypatch.setattr(rcsb, "parse_mmcif", lambda *args: parsed)
    parsed.sequences = {}
    target = rcsb.parse(rcsb.PDB("EXAMPLE", "unused"), {}, {"example_7": "legacy"})
    assert [c.cluster_id for c in target.record.chains] == [
        "legacy",
        "legacy",
        -1,
        -1,
        -1,
        -1,
    ]


def test_sequence_lookup_uses_chain_index_not_a_truncated_name(
    scripts, parsed, monkeypatch
):
    rcsb, cluster = scripts
    monkeypatch.setattr(rcsb, "parse_mmcif", lambda *args: parsed)
    parsed.data.chains[1]["name"] = parsed.data.chains[0]["name"]
    parsed.sequences[1] = "GGGGGGGGGGGG"
    clusters = {
        cluster.hash_sequence(seq): str(index)
        for index, seq in parsed.sequences.items()
    }
    target = rcsb.parse(rcsb.PDB("example", "unused"), {}, clusters)
    assert [c.cluster_id for c in target.record.chains[:2]] == ["0", "1"]


def test_a_directory_is_not_linked_as_an_msa(scripts, parsed, tmp_path, monkeypatch):
    rcsb, cluster = scripts
    monkeypatch.setattr(rcsb, "parse_mmcif", lambda *args: parsed)
    (tmp_path / f"{cluster.hash_sequence(parsed.sequences[0])}.npz").mkdir()
    target = rcsb.parse(rcsb.PDB("example", "unused"), {}, {}, msa_dir=tmp_path)
    assert all(c.msa_id == "" for c in target.record.chains)


@pytest.mark.parametrize("workers", [1, 2])
def test_msa_directory_is_forwarded_to_serial_and_parallel_workers(
    scripts, tmp_path, monkeypatch, workers
):
    rcsb, _ = scripts
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    clusters = tmp_path / "clusters.json"
    clusters.write_text("{}")
    args = Namespace(
        msa_dir=msa_dir,
        outdir=tmp_path / "output",
        clusters=clusters,
        datadir=tmp_path,
        redis_host="unused",
        redis_port=1,
        num_processes=workers,
    )
    received = []
    monkeypatch.setattr(rcsb, "Resource", lambda **kwargs: {})
    monkeypatch.setattr(
        rcsb, "fetch", lambda *args: [rcsb.PDB("a", "a"), rcsb.PDB("b", "b")]
    )
    monkeypatch.setattr(
        rcsb,
        "process_structure",
        lambda *args, **kwargs: received.append(kwargs["msa_dir"]),
    )
    monkeypatch.setattr(
        rcsb, "p_umap", lambda fn, items, **kwargs: [fn(item) for item in items]
    )
    monkeypatch.setattr(rcsb, "finalize", lambda *args: None)
    monkeypatch.setattr(
        rcsb.rdkit.Chem, "SetDefaultPickleProperties", lambda *args: None
    )
    rcsb.process(args)
    assert received == [msa_dir, msa_dir]


def test_missing_msa_directory_is_refused_before_creating_outputs(scripts, tmp_path):
    rcsb, _ = scripts
    args = Namespace(msa_dir=tmp_path / "missing", outdir=tmp_path / "output")
    with pytest.raises(ValueError, match="Processed MSA directory does not exist"):
        rcsb.process(args)
    assert not args.outdir.exists()


def test_process_structure_writes_linked_record(scripts, parsed, tmp_path, monkeypatch):
    rcsb, cluster = scripts
    monkeypatch.setattr(rcsb, "parse_mmcif", lambda *args: parsed)
    for name in ["structures", "records", "msa"]:
        (tmp_path / name).mkdir()
    seq_id = cluster.hash_sequence(parsed.sequences[0])
    (tmp_path / "msa" / f"{seq_id}.npz").touch()
    rcsb.process_structure(
        rcsb.PDB("example", "unused"),
        {},
        tmp_path,
        [],
        {seq_id: "protein"},
        msa_dir=tmp_path / "msa",
    )
    record = Record.load(tmp_path / "records/example.json")
    assert record.chains[0].msa_id == seq_id
    assert record.chains[0].cluster_id == "protein"
    assert record.chains[0].entity_id == 7


def test_real_mmcif_retains_full_polymer_sequence_and_current_chain_dtype(
    scripts, tmp_path
):
    rcsb, cluster = scripts
    components = {}
    for name in ["ALA", "MET", "GLY"]:
        mol = Chem.RWMol()
        for atom_name in const.ref_atoms[name]:
            atom = Chem.Atom(atom_name[0])
            atom.SetProp("name", atom_name)
            mol.AddAtom(atom)
        mol = mol.GetMol()
        conformer = Chem.Conformer(mol.GetNumAtoms())
        conformer.SetProp("name", "Ideal")
        mol.AddConformer(conformer)
        components[name] = mol
    structure = gemmi.Structure()
    structure.name = "toy"
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for number, name in [(1, "ALA"), (3, "GLY")]:
        residue = gemmi.Residue()
        residue.name, residue.seqid = name, gemmi.SeqId(number, " ")
        residue.subchain, residue.entity_id, residue.label_seq, residue.het_flag = (
            "Axp",
            "1",
            number,
            "A",
        )
        for i, atom_name in enumerate(const.ref_atoms[name]):
            atom = gemmi.Atom()
            atom.name, atom.element = atom_name, gemmi.Element(atom_name[0])
            atom.pos = gemmi.Position(number * 3 + i * 0.1, 0, 0)
            residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    entity = gemmi.Entity("1")
    entity.entity_type, entity.polymer_type = (
        gemmi.EntityType.Polymer,
        gemmi.PolymerType.PeptideL,
    )
    entity.subchains, entity.full_sequence = ["Axp"], ["ALA", "MSE", "GLY"]
    structure.entities.append(entity)
    path = tmp_path / "toy.cif"
    structure.make_mmcif_document().write_file(str(path))
    parsed = rcsb.parse_mmcif(str(path), components)
    assert parsed.sequences == {0: "AMG"}
    assert parsed.data.chains["cyclic_period"].tolist() == [0]
    assert parsed.data.residues["is_present"].tolist() == [True, False, True]
    seq_id = cluster.hash_sequence("AMG")
    target = rcsb.parse(
        rcsb.PDB("toy", str(path)), components, {seq_id: "full-sequence-cluster"}
    )
    assert target.record.chains[0].cluster_id == "full-sequence-cluster"
    assert target.record.chains[0].entity_id == 0
