import itertools
import pickle
import random
from pathlib import Path

import numpy as np
import torch

from boltz.data import const
from boltz.data.pad import pad_dim
from boltz.model.loss.confidence import lddt_dist
from boltz.model.loss.validation import weighted_minimum_rmsd_single


def convert_atom_name(name: str) -> tuple[int, int, int, int]:
    """Convert an atom name to a standard format.

    Parameters
    ----------
    name : str
        The atom name.

    Returns
    -------
    Tuple[int, int, int, int]
        The converted atom name.

    """
    name = name.strip()
    name = [ord(c) - 32 for c in name]
    name = name + [0] * (4 - len(name))
    return tuple(name)


def get_symmetries(path: str) -> dict:
    """Create a dictionary for the ligand symmetries.

    Parameters
    ----------
    path : str
        The path to the ligand symmetries.

    Returns
    -------
    dict
        The ligand symmetries.

    """
    with Path(path).open("rb") as f:
        data: dict = pickle.load(f)  # noqa: S301

    symmetries = {}
    for key, mol in data.items():
        try:
            serialized_sym = bytes.fromhex(mol.GetProp("symmetries"))
            sym = pickle.loads(serialized_sym)  # noqa: S301
            atom_names = []
            for atom in mol.GetAtoms():
                # Get atom name
                atom_name = convert_atom_name(atom.GetProp("name"))
                atom_names.append(atom_name)

            symmetries[key] = (sym, atom_names)
        except Exception:  # noqa: BLE001, PERF203, S110
            pass

    return symmetries


def compute_symmetry_idx_dictionary(data):
    # Compute the symmetry index dictionary
    total_count = 0
    all_coords = []
    for i, chain in enumerate(data.chains):
        chain.start_idx = total_count
        for j, token in enumerate(chain.tokens):
            token.start_idx = total_count - chain.start_idx
            all_coords.extend(
                [[atom.coords.x, atom.coords.y, atom.coords.z] for atom in token.atoms]
            )
            total_count += len(token.atoms)
    return all_coords


def get_current_idx_list(data):
    idx = []
    for chain in data.chains:
        if chain.in_crop:
            for token in chain.tokens:
                if token.in_crop:
                    idx.extend(
                        [
                            chain.start_idx + token.start_idx + i
                            for i in range(len(token.atoms))
                        ]
                    )
    return idx


def all_different_after_swap(l):
    final = [s[-1] for s in l]
    return len(final) == len(set(final))


def minimum_symmetry_coords(
    coords: torch.Tensor,
    feats: dict,
    index_batch: int,
    **args_rmsd,
):
    all_coords = feats["all_coords"][index_batch].unsqueeze(0).to(coords)
    all_resolved_mask = (
        feats["all_resolved_mask"][index_batch].to(coords).to(torch.bool)
    )
    crop_to_all_atom_map = (
        feats["crop_to_all_atom_map"][index_batch].to(coords).to(torch.long)
    )
    chain_symmetries = feats["chain_symmetries"][index_batch]
    amino_acids_symmetries = feats["amino_acids_symmetries"][index_batch]
    ligand_symmetries = feats["ligand_symmetries"][index_batch]

    # Check best symmetry on chain swap
    best_true_coords = None
    best_rmsd = float("inf")
    best_align_weights = None
    for c in chain_symmetries:
        true_all_coords = all_coords.clone()
        true_all_resolved_mask = all_resolved_mask.clone()
        for start1, end1, start2, end2, chainidx1, chainidx2 in c:
            true_all_coords[:, start1:end1] = all_coords[:, start2:end2]
            true_all_resolved_mask[start1:end1] = all_resolved_mask[start2:end2]
        true_coords = true_all_coords[:, crop_to_all_atom_map]
        true_resolved_mask = true_all_resolved_mask[crop_to_all_atom_map]
        true_coords = pad_dim(true_coords, 1, coords.shape[1] - true_coords.shape[1])
        true_resolved_mask = pad_dim(
            true_resolved_mask,
            0,
            coords.shape[1] - true_resolved_mask.shape[0],
        )
        try:
            rmsd, aligned_coords, align_weights = weighted_minimum_rmsd_single(
                coords,
                true_coords,
                atom_mask=true_resolved_mask,
                atom_to_token=feats["atom_to_token"][index_batch : index_batch + 1],
                mol_type=feats["mol_type"][index_batch : index_batch + 1],
                **args_rmsd,
            )
        except:
            print("Warning: error in rmsd computation inside symmetry code")
            continue
        rmsd = rmsd.item()

        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_true_coords = aligned_coords
            best_align_weights = align_weights
            best_true_resolved_mask = true_resolved_mask

    # atom symmetries (nucleic acid and protein residues), resolved greedily without recomputing alignment
    true_coords = best_true_coords.clone()
    true_resolved_mask = best_true_resolved_mask.clone()
    for symmetric_amino in amino_acids_symmetries:
        for c in symmetric_amino:
            # starting from greedy best, try to swap the atoms
            new_true_coords = true_coords.clone()
            new_true_resolved_mask = true_resolved_mask.clone()
            for i, j in c:
                new_true_coords[:, i] = true_coords[:, j]
                new_true_resolved_mask[i] = true_resolved_mask[j]

            # compute squared distance, for efficiency we do not recompute the alignment
            best_mse_loss = torch.sum(
                ((coords - best_true_coords) ** 2).sum(dim=-1)
                * best_align_weights
                * best_true_resolved_mask,
                dim=-1,
            ) / torch.sum(best_align_weights * best_true_resolved_mask, dim=-1)
            new_mse_loss = torch.sum(
                ((coords - new_true_coords) ** 2).sum(dim=-1)
                * best_align_weights
                * new_true_resolved_mask,
                dim=-1,
            ) / torch.sum(best_align_weights * new_true_resolved_mask, dim=-1)

            if best_mse_loss > new_mse_loss:
                best_true_coords = new_true_coords
                best_true_resolved_mask = new_true_resolved_mask

        # greedily update best coordinates after each amino acid
        true_coords = best_true_coords.clone()
        true_resolved_mask = best_true_resolved_mask.clone()

    # Recomputing alignment
    rmsd, true_coords, best_align_weights = weighted_minimum_rmsd_single(
        coords,
        true_coords,
        atom_mask=true_resolved_mask,
        atom_to_token=feats["atom_to_token"][index_batch : index_batch + 1],
        mol_type=feats["mol_type"][index_batch : index_batch + 1],
        **args_rmsd,
    )
    best_rmsd = rmsd.item()

    # atom symmetries (ligand and non-standard), resolved greedily recomputing alignment
    for symmetric_ligand in ligand_symmetries:
        for c in symmetric_ligand:
            new_true_coords = true_coords.clone()
            new_true_resolved_mask = true_resolved_mask.clone()
            for i, j in c:
                new_true_coords[:, j] = true_coords[:, i]
                new_true_resolved_mask[j] = true_resolved_mask[i]
            try:
                # TODO if this is too slow maybe we can get away with not recomputing alignment
                rmsd, aligned_coords, align_weights = weighted_minimum_rmsd_single(
                    coords,
                    new_true_coords,
                    atom_mask=new_true_resolved_mask,
                    atom_to_token=feats["atom_to_token"][index_batch : index_batch + 1],
                    mol_type=feats["mol_type"][index_batch : index_batch + 1],
                    **args_rmsd,
                )
            except Exception as e:
                raise e
                print(e)
                continue
            rmsd = rmsd.item()
            if rmsd < best_rmsd:
                best_true_coords = aligned_coords
                best_rmsd = rmsd
                best_true_resolved_mask = new_true_resolved_mask

        true_coords = best_true_coords.clone()
        true_resolved_mask = best_true_resolved_mask.clone()

    return best_true_coords, best_rmsd, best_true_resolved_mask.unsqueeze(0)


def minimum_lddt_symmetry_coords(
    coords: torch.Tensor,
    feats: dict,
    index_batch: int,
    **args_rmsd,
):
    all_coords = feats["all_coords"][index_batch].unsqueeze(0).to(coords)
    all_resolved_mask = (
        feats["all_resolved_mask"][index_batch].to(coords).to(torch.bool)
    )
    crop_to_all_atom_map = (
        feats["crop_to_all_atom_map"][index_batch].to(coords).to(torch.long)
    )
    chain_symmetries = feats["chain_symmetries"][index_batch]
    amino_acids_symmetries = feats["amino_acids_symmetries"][index_batch]
    ligand_symmetries = feats["ligand_symmetries"][index_batch]

    dmat_predicted = torch.cdist(
        coords[:, : len(crop_to_all_atom_map)], coords[:, : len(crop_to_all_atom_map)]
    )

    # Check best symmetry on chain swap
    best_true_coords = None
    best_lddt = 0
    for c in chain_symmetries:
        true_all_coords = all_coords.clone()
        true_all_resolved_mask = all_resolved_mask.clone()
        for start1, end1, start2, end2, chainidx1, chainidx2 in c:
            true_all_coords[:, start1:end1] = all_coords[:, start2:end2]
            true_all_resolved_mask[start1:end1] = all_resolved_mask[start2:end2]
        true_coords = true_all_coords[:, crop_to_all_atom_map]
        true_resolved_mask = true_all_resolved_mask[crop_to_all_atom_map]
        dmat_true = torch.cdist(true_coords, true_coords)
        pair_mask = (
            true_resolved_mask[:, None]
            * true_resolved_mask[None, :]
            * (1 - torch.eye(len(true_resolved_mask))).to(true_resolved_mask)
        )

        lddt = lddt_dist(
            dmat_predicted, dmat_true, pair_mask, cutoff=15.0, per_atom=False
        )[0]
        lddt = lddt.item()

        if lddt > best_lddt:
            best_lddt = lddt
            best_true_coords = true_coords
            best_true_resolved_mask = true_resolved_mask

    # atom symmetries (nucleic acid and protein residues), resolved greedily without recomputing alignment
    true_coords = best_true_coords.clone()
    true_resolved_mask = best_true_resolved_mask.clone()
    for symmetric_amino_or_lig in amino_acids_symmetries + ligand_symmetries:
        for c in symmetric_amino_or_lig:
            # starting from greedy best, try to swap the atoms
            new_true_coords = true_coords.clone()
            new_true_resolved_mask = true_resolved_mask.clone()
            indices = []
            for i, j in c:
                new_true_coords[:, i] = true_coords[:, j]
                new_true_resolved_mask[i] = true_resolved_mask[j]
                indices.append(i)

            indices = (
                torch.from_numpy(np.asarray(indices)).to(new_true_coords.device).long()
            )

            pred_coords_subset = coords[:, : len(crop_to_all_atom_map)][:, indices]
            true_coords_subset = true_coords[:, indices]
            new_true_coords_subset = new_true_coords[:, indices]

            sub_dmat_pred = torch.cdist(
                coords[:, : len(crop_to_all_atom_map)], pred_coords_subset
            )
            sub_dmat_true = torch.cdist(true_coords, true_coords_subset)
            sub_dmat_new_true = torch.cdist(new_true_coords, new_true_coords_subset)

            sub_true_pair_lddt = (
                true_resolved_mask[:, None] * true_resolved_mask[None, indices]
            )
            sub_true_pair_lddt[indices] = (
                sub_true_pair_lddt[indices]
                * (1 - torch.eye(len(indices))).to(sub_true_pair_lddt).bool()
            )

            sub_new_true_pair_lddt = (
                new_true_resolved_mask[:, None] * new_true_resolved_mask[None, indices]
            )
            sub_new_true_pair_lddt[indices] = (
                sub_new_true_pair_lddt[indices]
                * (1 - torch.eye(len(indices))).to(sub_true_pair_lddt).bool()
            )

            lddt = lddt_dist(
                sub_dmat_pred,
                sub_dmat_true,
                sub_true_pair_lddt,
                cutoff=15.0,
                per_atom=False,
            )[0]
            new_lddt = lddt_dist(
                sub_dmat_pred,
                sub_dmat_new_true,
                sub_new_true_pair_lddt,
                cutoff=15.0,
                per_atom=False,
            )[0]

            if new_lddt > lddt:
                best_true_coords = new_true_coords
                best_true_resolved_mask = new_true_resolved_mask

        # greedily update best coordinates after each amino acid
        true_coords = best_true_coords.clone()
        true_resolved_mask = best_true_resolved_mask.clone()

    # Recomputing alignment
    true_coords = pad_dim(true_coords, 1, coords.shape[1] - true_coords.shape[1])
    true_resolved_mask = pad_dim(
        true_resolved_mask,
        0,
        coords.shape[1] - true_resolved_mask.shape[0],
    )

    try:
        rmsd, true_coords, _ = weighted_minimum_rmsd_single(
            coords,
            true_coords,
            atom_mask=true_resolved_mask,
            atom_to_token=feats["atom_to_token"][index_batch : index_batch + 1],
            mol_type=feats["mol_type"][index_batch : index_batch + 1],
            **args_rmsd,
        )
        best_rmsd = rmsd.item()
    except Exception as e:
        print("Failed proper RMSD computation, returning inf. Error: ", e)
        best_rmsd = 1000

    return true_coords, best_rmsd, true_resolved_mask.unsqueeze(0)


def compute_all_coords_mask(structure):
    # Compute all coords, crop mask and add start_idx to structure
    total_count = 0
    all_coords = []
    all_coords_crop_mask = []
    all_resolved_mask = []
    for i, chain in enumerate(structure.chains):
        chain.start_idx = total_count
        for j, token in enumerate(chain.tokens):
            token.start_idx = total_count - chain.start_idx
            all_coords.extend(
                [[atom.coords.x, atom.coords.y, atom.coords.z] for atom in token.atoms]
            )
            all_coords_crop_mask.extend(
                [token.in_crop for _ in range(len(token.atoms))]
            )
            all_resolved_mask.extend(
                [token.is_present for _ in range(len(token.atoms))]
            )
            total_count += len(token.atoms)
    if len(all_coords_crop_mask) != len(all_resolved_mask):
        pass
    return all_coords, all_coords_crop_mask, all_resolved_mask


def get_chain_symmetries(cropped, max_n_symmetries=100):
    # get all coordinates and resolved mask
    structure = cropped.structure
    all_coords = []
    all_resolved_mask = []
    original_atom_idx = []
    chain_atom_idx = []
    chain_atom_num = []
    chain_in_crop = []
    chain_asym_id = []
    new_atom_idx = 0

    for chain in structure.chains:
        atom_idx, atom_num = (
            chain["atom_idx"],
            chain["atom_num"],
        )

        # compute coordinates and resolved mask
        resolved_mask = structure.atoms["is_present"][atom_idx : atom_idx + atom_num]

        # ensemble_atom_starts = [structure.ensemble[idx]["atom_coord_idx"] for idx in cropped.ensemble_ref_idxs]
        # coords = np.array(
        #    [structure.coords[ensemble_atom_start + atom_idx: ensemble_atom_start + atom_idx + atom_num]["coords"] for
        #     ensemble_atom_start in ensemble_atom_starts])

        coords = structure.atoms["coords"][atom_idx : atom_idx + atom_num]

        in_crop = False
        for token in cropped.tokens:
            if token["asym_id"] == chain["asym_id"]:
                in_crop = True
                break

        all_coords.append(coords)
        all_resolved_mask.append(resolved_mask)
        original_atom_idx.append(atom_idx)
        chain_atom_idx.append(new_atom_idx)
        chain_atom_num.append(atom_num)
        chain_in_crop.append(in_crop)
        chain_asym_id.append(chain["asym_id"])

        new_atom_idx += atom_num

    # Compute backmapping from token to all coords
    crop_to_all_atom_map = []
    for token in cropped.tokens:
        chain_idx = chain_asym_id.index(token["asym_id"])
        start = (
            chain_atom_idx[chain_idx] - original_atom_idx[chain_idx] + token["atom_idx"]
        )
        crop_to_all_atom_map.append(np.arange(start, start + token["atom_num"]))

    # Compute the symmetries between chains
    swaps = []
    for i, chain in enumerate(structure.chains):
        start = chain_atom_idx[i]
        end = start + chain_atom_num[i]
        if chain_in_crop[i]:
            possible_swaps = []
            for j, chain2 in enumerate(structure.chains):
                start2 = chain_atom_idx[j]
                end2 = start2 + chain_atom_num[j]
                if (
                    chain["entity_id"] == chain2["entity_id"]
                    and end - start == end2 - start2
                ):
                    possible_swaps.append((start, end, start2, end2, i, j))
            swaps.append(possible_swaps)
    combinations = itertools.product(*swaps)
    # to avoid combinatorial explosion, bound the number of combinations even considered
    combinations = list(itertools.islice(combinations, max_n_symmetries * 10))
    # filter for all chains getting a different assignment
    combinations = [c for c in combinations if all_different_after_swap(c)]

    if len(combinations) > max_n_symmetries:
        combinations = random.sample(combinations, max_n_symmetries)

    if len(combinations) == 0:
        combinations.append([])

    features = {}
    features["all_coords"] = torch.Tensor(
        np.concatenate(all_coords, axis=0)
    )  # axis=1 with ensemble

    features["all_resolved_mask"] = torch.Tensor(
        np.concatenate(all_resolved_mask, axis=0)
    )
    features["crop_to_all_atom_map"] = torch.Tensor(
        np.concatenate(crop_to_all_atom_map, axis=0)
    )
    features["chain_symmetries"] = combinations

    return features


def get_amino_acids_symmetries(cropped):
    # Compute standard amino-acids symmetries
    swaps = []
    start_index_crop = 0
    for token in cropped.tokens:
        symmetries = const.ref_symmetries.get(const.tokens[token["res_type"]], [])
        if len(symmetries) > 0:
            residue_swaps = []
            for sym in symmetries:
                sym_new_idx = [
                    (i + start_index_crop, j + start_index_crop) for i, j in sym
                ]
                residue_swaps.append(sym_new_idx)
            swaps.append(residue_swaps)
        start_index_crop += token["atom_num"]

    features = {"amino_acids_symmetries": swaps}
    return features


def get_ligand_symmetries(cropped, symmetries):
    # Compute ligand and non-standard amino-acids symmetries
    structure = cropped.structure

    added_molecules = {}
    index_mols = []
    atom_count = 0
    for token in cropped.tokens:
        # check if molecule is already added by identifying it through asym_id and res_idx
        atom_count += token["atom_num"]
        mol_id = (token["asym_id"], token["res_idx"])
        if mol_id in added_molecules.keys():
            added_molecules[mol_id] += token["atom_num"]
            continue
        added_molecules[mol_id] = token["atom_num"]

        # get the molecule type and indices
        residue_idx = token["res_idx"] + structure.chains[token["asym_id"]]["res_idx"]
        mol_name = structure.residues[residue_idx]["name"]
        atom_idx = structure.residues[residue_idx]["atom_idx"]
        mol_atom_names = structure.atoms[
            atom_idx : atom_idx + structure.residues[residue_idx]["atom_num"]
        ]["name"]
        mol_atom_names = [tuple(m) for m in mol_atom_names]
        if mol_name not in const.ref_symmetries.keys():
            index_mols.append(
                (mol_name, atom_count - token["atom_num"], mol_id, mol_atom_names)
            )

    # for each molecule, get the symmetries
    molecule_symmetries = []
    for mol_name, start_mol, mol_id, mol_atom_names in index_mols:
        if not mol_name in symmetries:
            continue
        else:
            swaps = []
            syms_ccd, mol_atom_names_ccd = symmetries[mol_name]
            # Get indices of mol_atom_names_ccd that are in mol_atom_names
            ccd_to_valid_ids = {
                mol_atom_names_ccd.index(name): i
                for i, name in enumerate(mol_atom_names)
            }
            ccd_valid_ids = set(ccd_to_valid_ids.keys())

            syms = []
            # Get syms
            for sym_ccd in syms_ccd:
                sym_dict = {}
                bool_add = True
                for i, j in enumerate(sym_ccd):
                    if i in ccd_valid_ids:
                        if j in ccd_valid_ids:
                            i_true = ccd_to_valid_ids[i]
                            j_true = ccd_to_valid_ids[j]
                            sym_dict[i_true] = j_true
                        else:
                            bool_add = False
                            break
                if bool_add:
                    syms.append([sym_dict[i] for i in range(len(ccd_valid_ids))])

            for sym in syms:
                if len(sym) != added_molecules[mol_id]:
                    raise Exception(
                        f"Symmetry length mismatch {len(sym)} {added_molecules[mol_id]}"
                    )
                # assert (
                #     len(sym) == added_molecules[mol_id]
                # ), f"Symmetry length mismatch {len(sym)} {added_molecules[mol_id]}"
                sym_new_idx = []
                for i, j in enumerate(sym):
                    if i != int(j):
                        sym_new_idx.append((i + start_mol, int(j) + start_mol))
                if len(sym_new_idx) > 0:
                    swaps.append(sym_new_idx)
            if len(swaps) > 0:
                molecule_symmetries.append(swaps)

    features = {"ligand_symmetries": molecule_symmetries}

    return features
