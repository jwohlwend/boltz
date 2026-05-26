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


"""DTensor-based v2 confidence module and heads.

This module provides the distributed implementation of the Boltz-2 ConfidenceModule
and ConfidenceHeads, including support for separate intra/inter-chain heads,
token-level pLDDT confidence, and updated iPLDDT weighting.

Only ``token_level_confidence=True`` is supported.  The atom-level path
(``token_level_confidence=False``) raises ``NotImplementedError``.

Placement conventions (3-D device mesh: [dp, cp_axis_0, cp_axis_1]):
  s:     (Shard(0), Shard(1), Replicate())   — single representation
  z:     (Shard(0), Shard(1), Shard(2))       — pair representation
  d:     (Shard(0), Shard(1), Shard(2))       — distance matrix
  x_pred:(Shard(0), Shard(1), Replicate())   — predicted coords
  scalar metrics: (Shard(0), Replicate(), Replicate())
"""

import warnings
from copy import deepcopy

import torch
from torch import nn
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard

from boltz.data import const
from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.atom_to_token import single_repr_rep_atom_to_token
from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, elementwise_op
from boltz.distributed.model.layers.embedding import EmbeddingParamsReplicated
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.outer_op import OuterOp, replicate_to_shard_outer_op
from boltz.distributed.model.layers.pairformer import PairformerModule
from boltz.distributed.model.layers.redistribute_transpose import redistribute_transpose
from boltz.distributed.model.layers.repeat_interleave import shardwise_repeat_interleave
from boltz.distributed.model.layers.shardwise_op import shardwise_distogram
from boltz.distributed.model.modules.confidence_utils import (
    compute_aggregated_metric,
    compute_ptms,
)
from boltz.distributed.model.modules.encoders import RelativePositionEncoder
from boltz.distributed.model.modules.trunkv2 import ContactConditioning
from boltz.distributed.utils import update_exhaustive_strides
from boltz.model.modules.confidencev2 import (
    IPLDDT_INTERFACE_WEIGHT,
    IPLDDT_LIGAND_WEIGHT,
    IPLDDT_NON_INTERFACE_WEIGHT,
)
from boltz.model.modules.confidencev2 import ConfidenceHeads as SerialConfidenceHeadsV2
from boltz.model.modules.confidencev2 import ConfidenceModule as SerialConfidenceModuleV2


class _ShardwiseWhere(torch.autograd.Function):
    """Select between two 4-D pair DTensors using a 3-D boolean condition, shardwise.

    Computes ``torch.where(cond[..., None], a, b)`` on local shards.
    Gradients flow to *a* where ``cond`` is True and to *b* where False.

    Communication budget: 0 collectives — purely shardwise.

    Parameters (forward)
    ----------
    cond_local : Tensor
        Plain (non-DTensor) bool tensor of shape ``(B_local*mult, N_row, N_col)``.
    a, b : DTensor
        Shape ``(B*mult, N, N, D)`` with identical placements (typically
        ``(Shard(0), Shard(1), Shard(2))``).
    """

    @staticmethod
    def forward(ctx: FunctionCtx, cond_local: torch.Tensor, a: DTensor, b: DTensor) -> DTensor:
        if not isinstance(a, DTensor) or not isinstance(b, DTensor):
            raise TypeError(f"Expected DTensors for a and b, got {type(a)} and {type(b)}")
        if a.device_mesh != b.device_mesh or a.placements != b.placements:
            raise ValueError("a and b must share the same device_mesh and placements")
        for p in a.placements:
            if isinstance(p, Partial):
                raise ValueError("Partial placements are not supported")
        expected_cond_shape = a.to_local().shape[:3]
        if cond_local.shape != expected_cond_shape:
            raise ValueError(
                f"cond_local shape {tuple(cond_local.shape)} does not match "
                f"local a shape prefix {tuple(expected_cond_shape)}"
            )

        a_local = a.to_local()
        b_local = b.to_local()
        cond_expanded = cond_local.unsqueeze(-1)
        result_local = torch.where(cond_expanded, a_local, b_local)

        ctx.save_for_backward(cond_local)
        ctx._a_requires_grad = a.requires_grad
        ctx._b_requires_grad = b.requires_grad
        ctx._device_mesh = a.device_mesh
        ctx._placements = a.placements
        ctx._shape = a.shape
        ctx._stride = a.stride()

        return DTensor.from_local(
            result_local,
            device_mesh=a.device_mesh,
            placements=a.placements,
            shape=a.shape,
            stride=a.stride(),
        )

    @staticmethod
    def backward(ctx: FunctionCtx, grad_output: DTensor):
        (cond_local,) = ctx.saved_tensors
        go_local = grad_output.to_local()
        cond_expanded = cond_local.unsqueeze(-1)
        zero = torch.zeros_like(go_local)

        d_a = (
            DTensor.from_local(
                torch.where(cond_expanded, go_local, zero),
                device_mesh=ctx._device_mesh,
                placements=ctx._placements,
                shape=ctx._shape,
                stride=ctx._stride,
            )
            if ctx._a_requires_grad
            else None
        )
        d_b = (
            DTensor.from_local(
                torch.where(cond_expanded, zero, go_local),
                device_mesh=ctx._device_mesh,
                placements=ctx._placements,
                shape=ctx._shape,
                stride=ctx._stride,
            )
            if ctx._b_requires_grad
            else None
        )
        return None, d_a, d_b


class ConfidenceHeads(nn.Module):
    """DTensor-based v2 confidence heads.

    Wraps the serial ``ConfidenceHeadsV2`` layer, distributing parameters with
    ``LinearParamsReplicated`` and adding sharded metric computation.

    Compared to the v1 distributed ``ConfidenceHeads``:
    * PAE is always computed (no ``compute_pae`` flag).
    * Optional ``use_separate_heads`` splits PAE/PDE into intra/inter-chain projections.
    * iPLDDT weights updated to ``ligand=20, interface=10, non_interface=1``.
    * PTM/iPTM always computed (with try/except fallback).

    Only ``token_level_confidence=True`` is supported.  Constructing with
    ``token_level_confidence=False`` raises ``NotImplementedError``.
    """

    def __init__(
        self,
        layer: SerialConfidenceHeadsV2,
        device_mesh: torch.distributed.device_mesh.DeviceMesh,
        transpose_comm: TransposeComm,
    ):
        super().__init__()

        # token_level_confidence = True is the default setting in the public checkpoint
        if not layer.token_level_confidence:
            raise NotImplementedError(
                "ConfidenceHeads distributed v2 only supports token_level_confidence=True. "
                "The atom-level confidence path is not implemented for DTensor."
            )

        self.device_mesh = device_mesh
        self.transpose_comm = transpose_comm
        self.token_level_confidence = layer.token_level_confidence
        self.use_separate_heads = layer.use_separate_heads

        # --- PAE / PDE heads ---
        if self.use_separate_heads:
            self.to_pae_intra_logits = LinearParamsReplicated(layer.to_pae_intra_logits, device_mesh)
            self.to_pae_inter_logits = LinearParamsReplicated(layer.to_pae_inter_logits, device_mesh)

            self.to_pde_intra_logits = LinearParamsReplicated(layer.to_pde_intra_logits, device_mesh)
            self.to_pde_inter_logits = LinearParamsReplicated(layer.to_pde_inter_logits, device_mesh)
        else:
            self.to_pae_logits = LinearParamsReplicated(layer.to_pae_logits, device_mesh)
            self.to_pde_logits = LinearParamsReplicated(layer.to_pde_logits, device_mesh)

        # --- pLDDT / resolved heads ---
        self.to_plddt_logits = LinearParamsReplicated(layer.to_plddt_logits, device_mesh)
        self.to_resolved_logits = LinearParamsReplicated(layer.to_resolved_logits, device_mesh)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    def forward(
        self,
        s: DTensor,
        z: DTensor,
        x_pred: DTensor,
        d: DTensor,
        feats: dict,
        pred_distogram_logits: DTensor,
        multiplicity: int = 1,
    ) -> dict[str, DTensor]:
        """Compute confidence logits and aggregated metrics.

        Parameters
        ----------
        s : DTensor
            Single representation.  Shape ``(B*mult, N, D_s)``,
            placements ``(Shard(0), Shard(1), Replicate())``.
        z : DTensor
            Pair representation.  Shape ``(B*mult, N, N, D_z)``,
            placements ``(Shard(0), Shard(1), Shard(2))``.
        x_pred : DTensor
            Predicted atom coordinates.  Shape ``(B*mult, N_atoms, 3)``,
            placements ``(Shard(0), Shard(1), Replicate())``.
        d : DTensor
            Token-level distance matrix.  Shape ``(B*mult, N, N)``,
            placements ``(Shard(0), Shard(1), Shard(2))``.
        feats : dict[str, DTensor]
            Feature dictionary.  Required keys: ``token_pad_mask``, ``asym_id``,
            ``mol_type``.
        pred_distogram_logits : DTensor
            Predicted distogram logits.  Shape ``(B, N, N, 64)``,
            placements ``(Shard(0), Shard(1), Shard(2))``.
        multiplicity : int
            Number of diffusion samples per input.

        Returns
        -------
        dict[str, DTensor]
        """
        self._validate_inputs(s, z, x_pred, d, feats, pred_distogram_logits)

        plddt_logits = self.to_plddt_logits(s)
        resolved_logits = self.to_resolved_logits(s)

        # Build same_chain mask once; reused by _ShardwiseWhere (when use_separate_heads)
        # and the no-grad iPLDDT metrics block (always).
        with torch.no_grad():
            same_chain_base = replicate_to_shard_outer_op(
                feats["asym_id"], OuterOp.EQUAL, axis=1, transpose_comm=self.transpose_comm
            ).to_local()  # (B_local, N_row, N_col)

        if self.use_separate_heads:
            # M = same_chain mask, A = intra, B = inter → torch.where(M, A, B)
            same_chain = same_chain_base.repeat_interleave(multiplicity, dim=0) if multiplicity > 1 else same_chain_base

            pae_logits = _ShardwiseWhere.apply(same_chain, self.to_pae_intra_logits(z), self.to_pae_inter_logits(z))

            # proj(z + z^T) = proj(z) + proj(z^T)
            z_pde_intra = self.to_pde_intra_logits(z)
            z_pde_intra_T = redistribute_transpose(
                z_pde_intra, self.transpose_comm, (Shard(0), Shard(1), Shard(2)), 1, 2
            )
            pde_intra = elementwise_op(z_pde_intra, z_pde_intra_T, ElementwiseOp.SUM)

            z_pde_inter = self.to_pde_inter_logits(z)
            z_pde_inter_T = redistribute_transpose(
                z_pde_inter, self.transpose_comm, (Shard(0), Shard(1), Shard(2)), 1, 2
            )
            pde_inter = elementwise_op(z_pde_inter, z_pde_inter_T, ElementwiseOp.SUM)

            pde_logits = _ShardwiseWhere.apply(same_chain, pde_intra, pde_inter)
        else:
            # Original path from boltz1
            pae_logits = self.to_pae_logits(z)

            z_proj = self.to_pde_logits(z)
            z_proj_T = redistribute_transpose(z_proj, self.transpose_comm, (Shard(0), Shard(1), Shard(2)), 1, 2)
            pde_logits = elementwise_op(z_proj, z_proj_T, ElementwiseOp.SUM)

        out_dict: dict[str, DTensor] = {
            "plddt_logits": plddt_logits,
            "pde_logits": pde_logits,
            "resolved_logits": resolved_logits,
            "pae_logits": pae_logits,
        }

        # ==================================================================
        # No-grad aggregated metrics (inference / logging only)
        # ==================================================================
        with torch.no_grad():
            token_pad_mask = feats["token_pad_mask"]
            mask_local = token_pad_mask.to_local()  # (B_local, N_local)
            B_local = mask_local.shape[0]
            N_local = mask_local.shape[1]

            # ---- pLDDT ----
            plddt = compute_aggregated_metric(plddt_logits)
            plddt_local = plddt.to_local()  # (B_local*mult, N_local)
            plddt_reshaped = plddt_local.reshape(B_local, multiplicity, N_local)

            masked_plddt = plddt_reshaped * mask_local.unsqueeze(1)
            num_local = masked_plddt.sum(dim=-1)  # (B_local, mult)
            den_local = mask_local.sum(dim=-1, keepdim=True)  # (B_local, 1)

            group_cp0 = self.device_mesh.get_group("cp_axis_0")
            torch.distributed.all_reduce(num_local, op=torch.distributed.ReduceOp.SUM, group=group_cp0)
            torch.distributed.all_reduce(den_local, op=torch.distributed.ReduceOp.SUM, group=group_cp0)

            complex_plddt_local = (num_local / den_local).reshape(B_local * multiplicity)
            complex_plddt = DTensor.from_local(
                complex_plddt_local,
                device_mesh=self.device_mesh,
                placements=(Shard(0), Replicate(), Replicate()),
                shape=(plddt.shape[0],),
                stride=(1,),
            )

            # ---- iPLDDT  (v2 weights: ligand=20, interface=10, non_interface=1) ----
            ligand_weight = IPLDDT_LIGAND_WEIGHT
            interface_weight = IPLDDT_INTERFACE_WEIGHT
            non_interface_weight = IPLDDT_NON_INTERFACE_WEIGHT

            mol_type_local = feats["mol_type"].to_local()  # (B_local, N_local)
            is_ligand_local = (mol_type_local == const.chain_type_ids["NONPOLYMER"]).float()

            d_local = d.to_local()  # (B_local*mult, N_row, N_col)
            is_contact_local = (d_local < 8).float()

            is_diff_chain_local = (~same_chain_base).float()

            # NOTE: because we use a square grid for now, N_row == N_col
            N_row = d_local.shape[1]
            N_col = d_local.shape[2]
            is_contact_4d = is_contact_local.reshape(B_local, multiplicity, N_row, N_col)
            is_diff_chain_4d = is_diff_chain_local.unsqueeze(1)
            non_ligand_4d = (1 - is_ligand_local).unsqueeze(1).unsqueeze(-1)

            interface_product = is_contact_4d * is_diff_chain_4d * non_ligand_4d
            token_interface_mask_local = interface_product.max(dim=-1).values  # (B_local, mult, N_row)

            group_cp1 = self.device_mesh.get_group("cp_axis_1")
            torch.distributed.all_reduce(
                token_interface_mask_local,
                op=torch.distributed.ReduceOp.MAX,
                group=group_cp1,
            )

            is_ligand_3d = is_ligand_local.unsqueeze(1)  # (B_local, 1, N_local)
            token_non_interface_mask = (1 - token_interface_mask_local) * (1 - is_ligand_3d)
            iplddt_weight_local = (
                is_ligand_3d * ligand_weight
                + token_interface_mask_local * interface_weight
                + token_non_interface_mask * non_interface_weight
            )  # (B_local, mult, N_local)

            masked_iplddt_w = mask_local.unsqueeze(1) * iplddt_weight_local
            num_iplddt = (plddt_reshaped * masked_iplddt_w).sum(dim=-1)
            den_iplddt = masked_iplddt_w.sum(dim=-1)

            torch.distributed.all_reduce(num_iplddt, op=torch.distributed.ReduceOp.SUM, group=group_cp0)
            torch.distributed.all_reduce(den_iplddt, op=torch.distributed.ReduceOp.SUM, group=group_cp0)

            complex_iplddt_local = (num_iplddt / den_iplddt).reshape(B_local * multiplicity)
            complex_iplddt = DTensor.from_local(
                complex_iplddt_local,
                device_mesh=self.device_mesh,
                placements=(Shard(0), Replicate(), Replicate()),
                shape=(plddt.shape[0],),
                stride=(1,),
            )

            # ---- PDE / iPDE ----
            pde = compute_aggregated_metric(pde_logits, end=32)

            pred_disto_local = pred_distogram_logits.to_local()
            if pred_disto_local.ndim == 5:
                if pred_disto_local.shape[-2] != 1:
                    raise ValueError(
                        f"ConfidenceHeads expects num_distograms=1, " f"got shape {pred_disto_local.shape}"
                    )
                pred_disto_local = pred_disto_local.squeeze(-2)
            pred_disto_prob = torch.softmax(pred_disto_local, dim=-1)
            contacts_mask = torch.zeros((1, 1, 1, 64), dtype=pred_disto_prob.dtype, device=pred_disto_prob.device)
            contacts_mask[:, :, :, :20] = 1.0
            prob_contact_local = (pred_disto_prob.unsqueeze(1) * contacts_mask).sum(-1)  # (B_local, 1, N_row, N_col)

            pde_local = pde.to_local().reshape(B_local, multiplicity, N_local, N_local)

            row_mask = mask_local  # (B_local, N_local)
            col_mask = redistribute_transpose(
                token_pad_mask,
                self.transpose_comm,
                (Shard(0), Replicate(), Shard(1)),
                dim0=None,
                dim1=None,
            ).to_local()

            prob_contact_local = prob_contact_local * row_mask[:, None, :, None] * col_mask[:, None, None, :]

            mesh_coord = self.device_mesh.get_coordinate()
            if mesh_coord[1] == mesh_coord[2]:
                diag_idx = torch.arange(0, N_local, device=prob_contact_local.device)
                prob_contact_local[:, :, diag_idx, diag_idx] = 0

            num_pde = (pde_local * prob_contact_local).sum(dim=(2, 3))
            den_pde = prob_contact_local.sum(dim=(2, 3))

            torch.distributed.all_reduce(num_pde, op=torch.distributed.ReduceOp.SUM, group=self.transpose_comm.group)
            torch.distributed.all_reduce(den_pde, op=torch.distributed.ReduceOp.SUM, group=self.transpose_comm.group)

            complex_pde_local = (num_pde / den_pde).reshape(B_local * multiplicity)
            complex_pde = DTensor.from_local(
                complex_pde_local,
                device_mesh=self.device_mesh,
                placements=(Shard(0), Replicate(), Replicate()),
                shape=(pde.shape[0],),
                stride=(1,),
            )

            # iPDE
            token_intf_pair = prob_contact_local * is_diff_chain_local.unsqueeze(1)
            num_ipde = (pde_local * token_intf_pair).sum(dim=(2, 3))
            den_ipde = token_intf_pair.sum(dim=(2, 3))

            torch.distributed.all_reduce(num_ipde, op=torch.distributed.ReduceOp.SUM, group=self.transpose_comm.group)
            torch.distributed.all_reduce(den_ipde, op=torch.distributed.ReduceOp.SUM, group=self.transpose_comm.group)

            complex_ipde_local = (num_ipde / (den_ipde + 1e-5)).reshape(B_local * multiplicity)
            complex_ipde = DTensor.from_local(
                complex_ipde_local,
                device_mesh=self.device_mesh,
                placements=(Shard(0), Replicate(), Replicate()),
                shape=(pde.shape[0],),
                stride=(1,),
            )

            # ---- PAE ----
            pae = compute_aggregated_metric(pae_logits, end=32)

            out_dict["plddt"] = plddt
            out_dict["pde"] = pde
            out_dict["pae"] = pae
            out_dict["complex_plddt"] = complex_plddt
            out_dict["complex_iplddt"] = complex_iplddt
            out_dict["complex_pde"] = complex_pde
            out_dict["complex_ipde"] = complex_ipde

            # --- PTM / iPTM ---
            # No try-except here: the serial v2 code has a broad `except Exception`
            # fallback that silently replaces PTM scores with zeros. We intentionally
            # omit it in the distributed path so that any error surfaces immediately.
            # If the serial path's fallback ever triggers, an equivalence test will
            # catch the mismatch (serial zeros vs distributed crash).
            ptm, iptm, ligand_iptm, protein_iptm, pair_chains_iptm = compute_ptms(
                pae_logits,
                x_pred,
                feats,
                multiplicity,
                self.transpose_comm,
            )
            out_dict["ptm"] = ptm
            out_dict["iptm"] = iptm
            out_dict["ligand_iptm"] = ligand_iptm
            out_dict["protein_iptm"] = protein_iptm
            out_dict["pair_chains_iptm"] = pair_chains_iptm

        return out_dict

    # ------------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------------

    def _validate_inputs(
        self,
        s: DTensor,
        z: DTensor,
        x_pred: DTensor,
        d: DTensor,
        feats: dict,
        pred_distogram_logits: DTensor,
    ) -> None:
        for name, tensor in [("s", s), ("z", z), ("x_pred", x_pred), ("d", d)]:
            if not isinstance(tensor, DTensor):
                raise TypeError(f"Expected DTensor for {name}, got {type(tensor)}")

        expected = {
            "s": (Shard(0), Shard(1), Replicate()),
            "z": (Shard(0), Shard(1), Shard(2)),
            "x_pred": (Shard(0), Shard(1), Replicate()),
            "d": (Shard(0), Shard(1), Shard(2)),
        }
        for name, tensor in [("s", s), ("z", z), ("x_pred", x_pred), ("d", d)]:
            if tensor.placements != expected[name]:
                raise ValueError(f"Expected {name} placements {expected[name]}, got {tensor.placements}")

        if pred_distogram_logits.placements != (Shard(0), Shard(1), Shard(2)):
            raise ValueError(
                f"Expected pred_distogram_logits placements (Shard(0), Shard(1), Shard(2)), "
                f"got {pred_distogram_logits.placements}"
            )

        for key in ("token_pad_mask", "asym_id", "mol_type"):
            feat = feats[key]
            if not isinstance(feat, DTensor):
                raise TypeError(f"Expected DTensor for feats['{key}'], got {type(feat)}")
            if feat.placements != (Shard(0), Shard(1), Replicate()):
                raise ValueError(
                    f"Expected feats['{key}'] placements (Shard(0), Shard(1), Replicate()), " f"got {feat.placements}"
                )

        N_global = feats["token_pad_mask"].shape[1]
        if s.shape[0] != z.shape[0]:
            raise ValueError(f"Batch dim mismatch: s.shape[0]={s.shape[0]} vs z.shape[0]={z.shape[0]}")
        if s.shape[1] != N_global:
            raise ValueError(f"Token dim mismatch: s.shape[1]={s.shape[1]} vs N_global={N_global}")
        if z.shape[1] != N_global or z.shape[2] != N_global:
            raise ValueError(
                f"Pair dims must equal N_global={N_global}, got z.shape[1]={z.shape[1]}, z.shape[2]={z.shape[2]}"
            )


class ConfidenceModule(nn.Module):
    """Distributed ConfidenceModule v2 (Algorithm 31).

    Wraps the serial :class:`~boltz.model.modules.confidencev2.ConfidenceModule`,
    distributing submodule parameters with DTensor-compatible layers and using
    sharded operations for pair computations.

    The forward pass:

    1. Normalize ``s_inputs``, ``s``, ``z``
    2. Optional ``s`` / ``z`` conditioning (``add_s_input_to_s``, ``add_z_input_to_z``)
    3. ``repeat_interleave`` s for multiplicity
    4. Outer-sum ``s_to_z`` pair update
    5. Optional: ``add_s_to_z_prod``
    6. Distogram chain: representative-atom projection → pairwise cdist →
       binning → embedding
    7. ``repeat_interleave`` z for multiplicity
    8. Pairformer stack
    9. :class:`ConfidenceHeads` for logit projections and aggregated metrics

    Only ``token_level_confidence=True`` is supported.

    Communication budget (forward only):

    - Norms / linears: no collectives (params are Replicate)
    - ``replicate_to_shard_outer_op``: 1 all-to-all per call
    - ``single_repr_rep_atom_to_token``: shardwise (0 collectives)
    - ``replicate_to_shard_outer_op(CDIST)``: 1 all-to-all
    - ``shardwise_distogram``: 0 collectives
    - ``PairformerModule``: O(depth) collectives (ring attention + triangle)
    - ``ConfidenceHeads``: O(1) all-reduces for aggregated metrics

    Parameters
    ----------
    module : SerialConfidenceModuleV2
        Initialised serial module whose weights are wrapped / transferred.
    dist_manager : DistributedManager
        Distributed manager with the 3-D device mesh (dp, cp_axis_0, cp_axis_1).
    transpose_comm : TransposeComm
        Base transpose-communication handle.  Deep copies are created
        internally for submodules that store their own handle.
    """

    def __init__(
        self,
        module: SerialConfidenceModuleV2,
        dist_manager: DistributedManager,
        transpose_comm: TransposeComm,
    ) -> None:
        super().__init__()

        if not module.token_level_confidence:
            raise NotImplementedError(
                "ConfidenceModule distributed v2 only supports token_level_confidence=True. "
                "The atom-level confidence path is not implemented for DTensor."
            )

        self.device_mesh = dist_manager.device_mesh_subgroups
        self.transpose_comm = transpose_comm

        self.no_update_s = module.no_update_s
        self.add_s_to_z_prod = module.add_s_to_z_prod
        self.add_s_input_to_s = module.add_s_input_to_s
        self.add_z_input_to_z = module.add_z_input_to_z
        self.return_latent_feats = module.return_latent_feats

        # ---- Buffer (plain tensor, not DTensor) ----
        self.register_buffer("boundaries", module.boundaries)

        # ---- LayerNorms ----
        self.s_inputs_norm = LayerNormParamsReplicated(module.s_inputs_norm, self.device_mesh)
        if not self.no_update_s:
            self.s_norm = LayerNormParamsReplicated(module.s_norm, self.device_mesh)
        self.z_norm = LayerNormParamsReplicated(module.z_norm, self.device_mesh)

        # ---- s → z projections ----
        self.s_to_z = LinearParamsReplicated(module.s_to_z, self.device_mesh)
        self.s_to_z_transpose = LinearParamsReplicated(module.s_to_z_transpose, self.device_mesh)

        if self.add_s_to_z_prod:
            self.s_to_z_prod_in1 = LinearParamsReplicated(module.s_to_z_prod_in1, self.device_mesh)
            self.s_to_z_prod_in2 = LinearParamsReplicated(module.s_to_z_prod_in2, self.device_mesh)
            self.s_to_z_prod_out = LinearParamsReplicated(module.s_to_z_prod_out, self.device_mesh)

        # ---- Optional s_input → s ----
        if self.add_s_input_to_s:
            self.s_input_to_s = LinearParamsReplicated(module.s_input_to_s, self.device_mesh)

        # ---- Optional z-input conditioning (rel_pos, bonds, contacts) ----
        if self.add_z_input_to_z:
            self.rel_pos = RelativePositionEncoder(
                module.rel_pos,
                device_mesh=self.device_mesh,
                transpose_comm=deepcopy(transpose_comm),
            )
            self.token_bonds = LinearParamsReplicated(module.token_bonds, self.device_mesh)
            self.bond_type_feature = getattr(module, "bond_type_feature", False)
            if self.bond_type_feature:
                self.token_bonds_type = EmbeddingParamsReplicated(module.token_bonds_type, self.device_mesh)
            self.contact_conditioning = ContactConditioning(module.contact_conditioning, device_mesh=self.device_mesh)

        # ---- Distogram embedding ----
        self.dist_bin_pairwise_embed = EmbeddingParamsReplicated(module.dist_bin_pairwise_embed, self.device_mesh)

        # ---- Pairformer ----
        self.pairformer_stack = PairformerModule(module.pairformer_stack, dist_manager)

        # ---- Confidence heads ----
        self.confidence_heads = ConfidenceHeads(
            module.confidence_heads,
            self.device_mesh,
            deepcopy(transpose_comm),
        )

    def forward(
        self,
        s_inputs: DTensor,
        s: DTensor,
        z: DTensor,
        x_pred: DTensor,
        feats: dict,
        pred_distogram_logits: DTensor,
        multiplicity: int = 1,
        run_sequentially: bool = False,
    ) -> dict[str, DTensor]:
        """Forward pass through the distributed confidence module.

        Parameters
        ----------
        s_inputs : DTensor
            Input single representation, shape ``(B, N, D_s)``,
            placements ``(Shard(0), Shard(1), Replicate())``.
        s : DTensor
            Trunk single representation (detached), same shape/placements.
        z : DTensor
            Trunk pair representation (detached), shape ``(B, N, N, D_z)``,
            placements ``(Shard(0), Shard(1), Shard(2))``.
        x_pred : DTensor
            Predicted atom coordinates, shape ``(B*mult, N_atoms, 3)``,
            placements ``(Shard(0), Shard(1), Replicate())``.
        feats : dict[str, DTensor]
            Feature dictionary.
        pred_distogram_logits : DTensor
            Predicted distogram logits, shape ``(B, N, N, K, bins)`` or ``(B, N, N, bins)``.
        multiplicity : int
            Number of diffusion samples.
        run_sequentially : bool
            If True and multiplicity > 1, run each multiplicity sample through
            the confidence module one at a time to reduce peak memory usage.

        Returns
        -------
        dict[str, DTensor]
            Confidence outputs including logits and aggregated metrics.
        """
        if run_sequentially and multiplicity > 1:
            return self._forward_sequentially(s_inputs, s, z, x_pred, feats, pred_distogram_logits, multiplicity)

        # ---- 1. Normalize inputs ----
        s_inputs = self.s_inputs_norm(s_inputs)
        if not self.no_update_s:
            s = self.s_norm(s)

        # ---- 2. Optional s_input addition to s ----
        if self.add_s_input_to_s:
            s = elementwise_op(s, self.s_input_to_s(s_inputs), ElementwiseOp.SUM)

        # ---- 3. Normalize z ----
        z = self.z_norm(z)

        # ---- 4. Optional z-input conditioning ----
        if self.add_z_input_to_z:
            z = elementwise_op(z, self.rel_pos(feats), ElementwiseOp.SUM)
            safe_dtype = z.dtype if z.dtype.is_floating_point else torch.float32
            z = elementwise_op(
                z,
                self.token_bonds(feats["token_bonds"].to(dtype=safe_dtype)),
                ElementwiseOp.SUM,
            )
            if self.bond_type_feature:
                z = elementwise_op(
                    z,
                    self.token_bonds_type(feats["type_bonds"].long()),
                    ElementwiseOp.SUM,
                )
            z = elementwise_op(z, self.contact_conditioning(feats), ElementwiseOp.SUM)

        # ---- 5. Repeat s for multiplicity ----
        s = shardwise_repeat_interleave(s, multiplicity, dim=0)

        # ---- 6. Outer-sum s → z ----
        # Serial: z += s_to_z(s_inputs)[:, :, None, :] + s_to_z_T(s_inputs)[:, None, :, :]
        s_to_z_pair = replicate_to_shard_outer_op(
            self.s_to_z(s_inputs),
            OuterOp.SUM,
            axis=1,
            transpose_comm=self.transpose_comm,
            input_t=self.s_to_z_transpose(s_inputs),
        )
        z = elementwise_op(z, s_to_z_pair, ElementwiseOp.SUM)

        # ---- 7. Optional outer-product s → z ----
        if self.add_s_to_z_prod:
            z_prod = replicate_to_shard_outer_op(
                self.s_to_z_prod_in1(s_inputs),
                OuterOp.PROD,
                axis=1,
                transpose_comm=self.transpose_comm,
                input_t=self.s_to_z_prod_in2(s_inputs),
            )
            z = elementwise_op(z, self.s_to_z_prod_out(z_prod), ElementwiseOp.SUM)

        # ---- 8. Distogram: x_pred → representative-atom token repr → cdist → bin → embed ----
        token_to_rep_atom = feats["token_to_rep_atom"]
        x_pred_repr = single_repr_rep_atom_to_token(x_pred, token_to_rep_atom)

        d = replicate_to_shard_outer_op(x_pred_repr, OuterOp.CDIST, axis=1, transpose_comm=self.transpose_comm)
        distogram = shardwise_distogram(d, self.boundaries)
        distogram = self.dist_bin_pairwise_embed(distogram)

        # ---- 9. Repeat z for multiplicity and add distogram ----
        z = shardwise_repeat_interleave(z, multiplicity, dim=0)
        z = elementwise_op(z, distogram, ElementwiseOp.SUM)

        # ---- 10. Masks for pairformer ----
        mask = shardwise_repeat_interleave(feats["token_pad_mask"], multiplicity, dim=0)
        pair_mask = shardwise_repeat_interleave(feats["token_pair_pad_mask"], multiplicity, dim=0)
        mask = mask.to(dtype=s.dtype)
        pair_mask = pair_mask.to(dtype=z.dtype)

        # ---- 11. Pairformer ----
        s, z = self.pairformer_stack(s, z, mask=mask, pair_mask=pair_mask)

        # ---- 12. Output dict ----
        out_dict: dict[str, DTensor] = {}
        if self.return_latent_feats:
            out_dict["s_conf"] = s
            out_dict["z_conf"] = z

        # ---- 13. Confidence heads ----
        out_dict.update(
            self.confidence_heads(
                s=s,
                z=z,
                x_pred=x_pred,
                d=d,
                feats=feats,
                pred_distogram_logits=pred_distogram_logits,
                multiplicity=multiplicity,
            )
        )
        return out_dict

    def _forward_sequentially(
        self,
        s_inputs: DTensor,
        s: DTensor,
        z: DTensor,
        x_pred: DTensor,
        feats: dict,
        pred_distogram_logits: DTensor,
        multiplicity: int,
    ) -> dict[str, DTensor]:
        """Run the confidence module one multiplicity sample at a time.

        This trades throughput for peak memory: instead of processing all
        ``multiplicity`` samples in a single forward pass, each sample is
        processed independently with ``multiplicity=1`` and the results
        are re-assembled at the end.
        """
        x_pred_local = x_pred.to_local()
        assert (
            x_pred_local.shape[0] % multiplicity == 0
        ), f"x_pred.shape[0] must be divisible by multiplicity, got {x_pred.shape[0]} and multiplicity {multiplicity}"
        B_local = x_pred_local.shape[0] // multiplicity
        B_global = x_pred.shape[0] // multiplicity

        if B_local > 1:
            warnings.warn(
                "B_local > 1 could cause deadlocking issues with pair_chains_iptm "
                "when chain counts are different on different dp groups"
            )

        x_pred_single_shape = torch.Size([B_global, *x_pred.shape[1:]])
        x_pred_unflat = x_pred_local.unflatten(0, (B_local, multiplicity))
        x_pred_single_stride = update_exhaustive_strides(x_pred.shape, x_pred.stride(), x_pred_single_shape)

        out_dicts: list[dict] = []
        for mult_idx in range(multiplicity):
            x_pred_sample = DTensor.from_local(
                x_pred_unflat[:, mult_idx : mult_idx + 1].flatten(0, 1),
                device_mesh=x_pred.device_mesh,
                placements=x_pred.placements,
                shape=x_pred_single_shape,
                stride=x_pred_single_stride,
            )
            out_dicts.append(
                self.forward(
                    s_inputs,
                    s,
                    z,
                    x_pred_sample,
                    feats,
                    pred_distogram_logits,
                    multiplicity=1,
                    run_sequentially=False,
                )
            )

        out_dict: dict[str, DTensor] = {}
        B_global_mult = x_pred.shape[0]
        for key in out_dicts[0]:
            if key != "pair_chains_iptm":
                ref = out_dicts[0][key]
                stacked = torch.stack([o[key].to_local() for o in out_dicts], dim=1)
                stacked_flattened = stacked.flatten(0, 1)
                out_shape = torch.Size([B_global_mult, *ref.shape[1:]])
                out_dict[key] = DTensor.from_local(
                    stacked_flattened,
                    device_mesh=ref.device_mesh,
                    placements=ref.placements,
                    shape=out_shape,
                    stride=update_exhaustive_strides(ref.shape, ref.stride(), out_shape),
                )
            else:
                pair_chains_iptm: dict = {}
                for idx1 in out_dicts[0][key]:
                    chain_iptm: dict = {}
                    for idx2 in out_dicts[0][key][idx1]:
                        ref = out_dicts[0][key][idx1][idx2]
                        stacked = torch.stack([o[key][idx1][idx2].to_local() for o in out_dicts], dim=1)
                        stacked_flattened = stacked.flatten(0, 1)
                        ref_shape = torch.Size([B_global_mult, *ref.shape[1:]])
                        chain_iptm[idx2] = DTensor.from_local(
                            stacked_flattened,
                            device_mesh=ref.device_mesh,
                            placements=ref.placements,
                            shape=ref_shape,
                            stride=update_exhaustive_strides(ref.shape, ref.stride(), ref_shape),
                        )
                    pair_chains_iptm[idx1] = chain_iptm
                out_dict[key] = pair_chains_iptm

        return out_dict
