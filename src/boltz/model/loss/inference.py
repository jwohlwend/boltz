# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# fmt: off
import torch

from boltz.data import const


def compute_chain_clashes(pred_atom_coords, feats, clash_buffer=0.4):
    chain_id = feats["asym_id"]
    with torch.autocast("cuda", enabled=False):
        atom_chain_id = (
            torch.bmm(feats["atom_to_token"].float(), chain_id.unsqueeze(-1).float())
            .squeeze(-1)
            .long()
        )

    vdw_radii = torch.zeros(
        const.num_elements, dtype=torch.float32, device=pred_atom_coords.device
    )
    vdw_radii[1:119] = torch.tensor(
        const.vdw_radii, dtype=torch.float32, device=pred_atom_coords.device
    )
    atom_vdw_radii = (feats["ref_element"].float() @ vdw_radii.unsqueeze(-1)).squeeze(
        -1
    )

    dists = torch.cdist(pred_atom_coords, pred_atom_coords)
    clashes = (
        dists
        < (atom_vdw_radii.unsqueeze(-1) + atom_vdw_radii.unsqueeze(-2)) - clash_buffer
    )

    multiplicity = pred_atom_coords.shape[0]
    num_clashes, num_pairs = {}, {}
    for key in const.out_single_types:
        num_clashes["sym_" + key] = torch.zeros(
            multiplicity, dtype=torch.float32, device=pred_atom_coords.device
        )
        num_pairs["sym_" + key] = torch.zeros(
            multiplicity, dtype=torch.float32, device=pred_atom_coords.device
        )
    for key in const.clash_types:
        num_clashes["asym_" + key] = torch.zeros(
            multiplicity, dtype=torch.float32, device=pred_atom_coords.device
        )
        num_pairs["asym_" + key] = torch.zeros(
            multiplicity, dtype=torch.float32, device=pred_atom_coords.device
        )

    for batch_idx in range(feats["atom_pad_mask"].shape[0]):  # TODO: Batch size > 1
        pair_pad_mask = (
            feats["atom_pad_mask"][batch_idx, :, None]
            * feats["atom_pad_mask"][batch_idx, None, :]
        ).bool()
        if feats["connections_edge_index"][batch_idx].shape[1] > 0:
            pair_pad_mask[
                feats["connections_edge_index"][batch_idx][0],
                feats["connections_edge_index"][batch_idx][1],
            ] = False
            pair_pad_mask[
                feats["connections_edge_index"][batch_idx][1],
                feats["connections_edge_index"][batch_idx][0],
            ] = False
        chain_symmetries = feats["chain_symmetries"][batch_idx]
        chain_id_to_symmetry = {}
        chain_id_to_type = {}
        for idx, symmetry in enumerate(chain_symmetries):
            for chain in symmetry:
                chain_id_to_symmetry[chain[0]] = idx
                chain_id_to_type[chain[0]] = chain[4]
        for i in chain_id_to_symmetry:
            for j in chain_id_to_symmetry:
                type1, type2 = (
                    const.chain_types[chain_id_to_type[i]],
                    const.chain_types[chain_id_to_type[j]],
                )
                if i >= j:
                    continue
                chain_pair_mask = (
                    pair_pad_mask
                    * (atom_chain_id[batch_idx] == i).unsqueeze(-1)
                    * (atom_chain_id[batch_idx] == j).unsqueeze(-2)
                )
                chain_pair_clashes = clashes[:, chain_pair_mask].any(dim=-1)
                if chain_id_to_symmetry[i] == chain_id_to_symmetry[j]:
                    num_clashes[
                        "sym_" + const.chain_type_to_out_single_type[type1]
                    ] += chain_pair_clashes.float()
                    num_pairs["sym_" + const.chain_type_to_out_single_type[type1]] += 1
                else:
                    num_clashes[
                        "asym_"
                        + const.chain_types_to_clash_type[frozenset((type1, type2))]
                    ] += chain_pair_clashes.float()
                    num_pairs[
                        "asym_"
                        + const.chain_types_to_clash_type[frozenset((type1, type2))]
                    ] += 1

    for key in num_clashes:
        if num_pairs[key].sum() > 0:
            num_clashes[key] /= num_pairs[key]

    return num_clashes, num_pairs


def compute_pb_geometry_metrics(
    pred_atom_coords, feats, bond_buffer=0.25, angle_buffer=0.25, clash_buffer=0.3
):
    with torch.autocast("cuda", enabled=False):
        chain_id = feats["asym_id"]
        atom_chain_id = (
            torch.bmm(feats["atom_to_token"].float(), chain_id.unsqueeze(-1).float())
            .squeeze(-1)
            .long()
        )
        _is_ligand_mask = (  # noqa: F841
            torch.bmm(
                feats["atom_to_token"].float(), feats["mol_type"].unsqueeze(-1).float()
            )
            .squeeze(-1)
            .long()
            == const.chain_type_ids["NONPOLYMER"]
        ).float()

    multiplicity = pred_atom_coords.shape[0]
    num_bond_length_failures = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_bond_angle_failures = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_internal_clash_failures = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_ligands = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )

    for index_batch in range(len(feats.get("ligand_edge_index", []))):
        if feats["ligand_edge_index"][index_batch].shape[1] == 0:
            continue
        dists = torch.linalg.norm(
            pred_atom_coords[:, feats["ligand_edge_index"][index_batch][0]]
            - pred_atom_coords[:, feats["ligand_edge_index"][index_batch][1]],
            dim=-1,
        )

        bond_length_violations = (
            (
                dists
                < feats["ligand_edge_lower_bounds"][index_batch] * (1.0 - bond_buffer)
            )
            + (
                dists
                > feats["ligand_edge_upper_bounds"][index_batch] * (1.0 + bond_buffer)
            )
        )[:, feats["ligand_edge_bond_mask"][index_batch]].float()

        bond_angle_violations = (
            (
                dists
                < feats["ligand_edge_lower_bounds"][index_batch] * (1.0 - angle_buffer)
            )
            + (
                dists
                > feats["ligand_edge_upper_bounds"][index_batch] * (1.0 + angle_buffer)
            )
        )[:, feats["ligand_edge_angle_mask"][index_batch]].float()

        internal_clash_violations = (
            dists
            < feats["ligand_edge_lower_bounds"][index_batch] * (1.0 - clash_buffer)
        )[
            :,
            ~(
                feats["ligand_edge_bond_mask"][index_batch]
                + feats["ligand_edge_angle_mask"][index_batch]
            ),
        ].float()

        edge_chain_ids = atom_chain_id[index_batch][
            feats["ligand_edge_index"][index_batch][0]
        ]
        bond_chain_ids = edge_chain_ids[feats["ligand_edge_bond_mask"][index_batch]]
        angle_chain_ids = edge_chain_ids[feats["ligand_edge_angle_mask"][index_batch]]
        internal_clash_chain_ids = edge_chain_ids[
            ~(
                feats["ligand_edge_bond_mask"][index_batch]
                + feats["ligand_edge_angle_mask"][index_batch]
            )
        ]

        num_bond_length_failures += (
            torch.zeros(
                (multiplicity, chain_id.max().item() + 1),
                dtype=torch.float32,
                device=dists.device,
            )
            .scatter_reduce(
                1,
                bond_chain_ids.expand((multiplicity, -1)),
                bond_length_violations,
                reduce="amax",
            )
            .sum(dim=-1)
        )

        num_bond_angle_failures += (
            torch.zeros(
                (multiplicity, chain_id.max().item() + 1),
                dtype=torch.float32,
                device=dists.device,
            )
            .scatter_reduce(
                1,
                angle_chain_ids.expand((multiplicity, -1)),
                bond_angle_violations,
                reduce="amax",
            )
            .sum(dim=-1)
        )

        num_internal_clash_failures += (
            torch.zeros(
                (multiplicity, chain_id.max().item() + 1),
                dtype=torch.float32,
                device=dists.device,
            )
            .scatter_reduce(
                1,
                internal_clash_chain_ids.expand((multiplicity, -1)),
                internal_clash_violations,
                reduce="amax",
            )
            .sum(dim=-1)
        )

        num_ligands += (
            torch.zeros(
                (chain_id.max().item() + 1,), dtype=torch.float32, device=dists.device
            )
            .scatter_reduce(
                0,
                edge_chain_ids,
                torch.ones(
                    edge_chain_ids.shape, dtype=torch.float32, device=dists.device
                ),
                reduce="amax",
            )
            .sum()
        )

    num_bond_length_failures[num_ligands > 0] /= num_ligands[num_ligands > 0]
    num_bond_angle_failures[num_ligands > 0] /= num_ligands[num_ligands > 0]
    num_internal_clash_failures[num_ligands > 0] /= num_ligands[num_ligands > 0]

    return (
        num_bond_length_failures,
        num_bond_angle_failures,
        num_internal_clash_failures,
        num_ligands,
    )


def compute_torsion_angles(coords, torsion_index):
    r_ij = coords.index_select(-2, torsion_index[0]) - coords.index_select(
        -2, torsion_index[1]
    )
    r_kj = coords.index_select(-2, torsion_index[2]) - coords.index_select(
        -2, torsion_index[1]
    )
    r_kl = coords.index_select(-2, torsion_index[2]) - coords.index_select(
        -2, torsion_index[3]
    )

    n_ijk = torch.cross(r_ij, r_kj, dim=-1)
    n_jkl = torch.cross(r_kj, r_kl, dim=-1)

    _r_kj_norm = torch.linalg.norm(r_kj, dim=-1)  # noqa: F841
    n_ijk_norm = torch.linalg.norm(n_ijk, dim=-1)
    n_jkl_norm = torch.linalg.norm(n_jkl, dim=-1)

    sign_phi = torch.sign(
        r_kj.unsqueeze(-2) @ torch.cross(n_ijk, n_jkl, dim=-1).unsqueeze(-1)
    ).squeeze(-1, -2)
    phi = sign_phi * torch.arccos(
        torch.clamp(
            (n_ijk.unsqueeze(-2) @ n_jkl.unsqueeze(-1)).squeeze(-1, -2)
            / (n_ijk_norm * n_jkl_norm),
            -1 + 1e-8,
            1 - 1e-8,
        )
    )
    return phi


def compute_stereo_metrics(pred_atom_coords, feats):
    multiplicity = pred_atom_coords.shape[0]
    num_chiral_atom_violations = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_chiral_atoms = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_stereo_bond_violations = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_stereo_bonds = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )

    for index_batch in range(len(feats.get("ligand_edge_index", []))):
        if feats["ligand_chiral_atom_index"][index_batch].shape[1] > 0:
            pred_chiral_torsion_angles = compute_torsion_angles(
                pred_atom_coords,
                feats["ligand_chiral_atom_index"][index_batch][
                    :, feats["ligand_chiral_check_mask"][index_batch].bool()
                ],
            )
            pred_chiral_atom_orientations = pred_chiral_torsion_angles > 0
            true_chiral_atom_orientations = feats["ligand_chiral_atom_orientations"][
                index_batch
            ][feats["ligand_chiral_check_mask"][index_batch].bool()]
            num_chiral_atom_violations += (
                pred_chiral_atom_orientations != true_chiral_atom_orientations
            ).sum(dim=-1)
            num_chiral_atoms += true_chiral_atom_orientations.shape[0]

        if feats["ligand_stereo_bond_index"][index_batch].shape[1] > 0:
            pred_stereo_torsion_angles = compute_torsion_angles(
                pred_atom_coords,
                feats["ligand_stereo_bond_index"][index_batch][
                    :, feats["ligand_stereo_check_mask"][index_batch].bool()
                ],
            )
            pred_stereo_bond_orientations = (
                torch.abs(pred_stereo_torsion_angles) > torch.pi / 2
            )
            true_stereo_bond_orientations = feats["ligand_stereo_bond_orientations"][
                index_batch
            ][feats["ligand_stereo_check_mask"][index_batch].bool()]
            num_stereo_bond_violations += (
                pred_stereo_bond_orientations != true_stereo_bond_orientations
            ).sum(dim=-1)
            num_stereo_bonds += true_stereo_bond_orientations.shape[0]

    num_chiral_atom_violations[num_chiral_atoms > 0] /= num_chiral_atoms[
        num_chiral_atoms > 0
    ]
    num_stereo_bond_violations[num_stereo_bonds > 0] /= num_stereo_bonds[
        num_stereo_bonds > 0
    ]
    return (
        num_chiral_atom_violations,
        num_chiral_atoms,
        num_stereo_bond_violations,
        num_stereo_bonds,
    )


def compute_pb_flatness_metrics(pred_atom_coords, feats, buffer=0.25):
    multiplicity = pred_atom_coords.shape[0]
    num_aromatic_5_violations = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_aromatic_5_rings = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_aromatic_6_violations = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_aromatic_6_rings = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_double_bond_violations = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )
    num_double_bonds = torch.zeros(
        multiplicity, dtype=torch.float32, device=pred_atom_coords.device
    )

    for index_batch in range(len(feats.get("ligand_aromatic_5_ring_index", []))):
        ring_5_index = feats["ligand_aromatic_5_ring_index"][index_batch].T
        ring_6_index = feats["ligand_aromatic_6_ring_index"][index_batch].T
        double_bond_index = feats["ligand_planar_double_bond_index"][index_batch].T

        ring_5_coords = pred_atom_coords[..., ring_5_index, :]
        ring_6_coords = pred_atom_coords[..., ring_6_index, :]
        double_bond_coords = pred_atom_coords[..., double_bond_index, :]

        centered_ring_5_coords = ring_5_coords - ring_5_coords.mean(
            dim=-2, keepdims=True
        )
        ring_5_vecs = torch.linalg.svd(centered_ring_5_coords)[2][..., -1, :, None]
        ring_5_dists = torch.abs(
            (centered_ring_5_coords @ ring_5_vecs).squeeze(dim=(-1, -2))
        )
        num_aromatic_5_violations += torch.any(ring_5_dists > buffer, dim=-1).sum(
            dim=-1
        )
        num_aromatic_5_rings += ring_5_index.shape[0]

        centered_ring_6_coords = ring_6_coords - ring_6_coords.mean(
            dim=-2, keepdims=True
        )
        ring_6_vecs = torch.linalg.svd(centered_ring_6_coords)[2][..., -1, :, None]
        ring_6_dists = torch.abs(
            (centered_ring_6_coords @ ring_6_vecs).squeeze(dim=(-1, -2))
        )
        num_aromatic_6_violations += torch.any(ring_6_dists > buffer, dim=-1).sum(
            dim=-1
        )
        num_aromatic_6_rings += ring_6_index.shape[0]

        centered_double_bond_coords = double_bond_coords - double_bond_coords.mean(
            dim=-2, keepdims=True
        )
        double_bond_vecs = torch.linalg.svd(centered_double_bond_coords)[2][
            ..., -1, :, None
        ]
        double_bond_dists = torch.abs(
            (centered_double_bond_coords @ double_bond_vecs).squeeze(dim=(-1, -2))
        )
        num_double_bond_violations += torch.any(double_bond_dists > buffer, dim=-1).sum(
            dim=-1
        )
        num_double_bonds += double_bond_index.shape[0]

    num_aromatic_5_violations[num_aromatic_5_rings > 0] /= num_aromatic_5_rings[
        num_aromatic_5_rings > 0
    ]
    num_aromatic_6_violations[num_aromatic_6_rings > 0] /= num_aromatic_6_rings[
        num_aromatic_6_rings > 0
    ]
    num_double_bond_violations[num_double_bonds > 0] /= num_double_bonds[
        num_double_bonds > 0
    ]

    return (
        num_aromatic_5_violations,
        num_aromatic_5_rings,
        num_aromatic_6_violations,
        num_aromatic_6_rings,
        num_double_bond_violations,
        num_double_bonds,
    )
