# started from code from https://github.com/lucidrains/alphafold3-pytorch, MIT License, Copyright (c) 2024 Phil Wang
from functools import partial
from math import pi

import torch
from einops import rearrange
from torch import nn
from torch.nn import Module, ModuleList
from torch.nn.functional import one_hot

import boltz.model.layers.initialize as init
from boltz.data import const
from boltz.model.layers.transition import Transition
from boltz.model.modules.transformers import AtomTransformer
from boltz.model.modules.utils import LinearNoBias


class FourierEmbedding(Module):
    """Fourier embedding layer."""

    def __init__(self, dim):
        """Initialize the Fourier Embeddings.

        Parameters
        ----------
        dim : int
            The dimension of the embeddings.

        """
        super().__init__()
        self.proj = nn.Linear(1, dim)
        torch.nn.init.normal_(self.proj.weight, mean=0, std=1)
        torch.nn.init.normal_(self.proj.bias, mean=0, std=1)
        self.proj.requires_grad_(False)

    def forward(
        self,
        times,
    ):
        times = rearrange(times, "b -> b 1")
        rand_proj = self.proj(times)
        return torch.cos(2 * pi * rand_proj)


class RelativePositionEncoder(Module):
    """Relative position encoder."""

    def __init__(self, token_z, r_max=32, s_max=2):
        """Initialize the relative position encoder.

        Parameters
        ----------
        token_z : int
            The pair representation dimension.
        r_max : int, optional
            The maximum index distance, by default 32.
        s_max : int, optional
            The maximum chain distance, by default 2.

        """
        super().__init__()
        self.r_max = r_max
        self.s_max = s_max
        self.linear_layer = LinearNoBias(4 * (r_max + 1) + 2 * (s_max + 1) + 1, token_z)

    def forward(self, feats):
        b_same_chain = torch.eq(
            feats["asym_id"][:, :, None], feats["asym_id"][:, None, :]
        )
        b_same_residue = torch.eq(
            feats["residue_index"][:, :, None], feats["residue_index"][:, None, :]
        )
        b_same_entity = torch.eq(
            feats["entity_id"][:, :, None], feats["entity_id"][:, None, :]
        )
        rel_pos = (
            feats["residue_index"][:, :, None] - feats["residue_index"][:, None, :]
        )
        if torch.any(feats["cyclic_period"] != 0):
            period = torch.where(
                feats["cyclic_period"] > 0,
                feats["cyclic_period"],
                torch.zeros_like(feats["cyclic_period"]) + 10000,
            ).unsqueeze(1)
            rel_pos = (rel_pos - period * torch.round(rel_pos / period)).long()

        d_residue = torch.clip(
            rel_pos + self.r_max,
            0,
            2 * self.r_max,
        )

        d_residue = torch.where(
            b_same_chain, d_residue, torch.zeros_like(d_residue) + 2 * self.r_max + 1
        )
        a_rel_pos = one_hot(d_residue, 2 * self.r_max + 2)

        d_token = torch.clip(
            feats["token_index"][:, :, None]
            - feats["token_index"][:, None, :]
            + self.r_max,
            0,
            2 * self.r_max,
        )
        d_token = torch.where(
            b_same_chain & b_same_residue,
            d_token,
            torch.zeros_like(d_token) + 2 * self.r_max + 1,
        )
        a_rel_token = one_hot(d_token, 2 * self.r_max + 2)

        d_chain = torch.clip(
            feats["sym_id"][:, :, None] - feats["sym_id"][:, None, :] + self.s_max,
            0,
            2 * self.s_max,
        )
        d_chain = torch.where(
            b_same_chain, torch.zeros_like(d_chain) + 2 * self.s_max + 1, d_chain
        )
        a_rel_chain = one_hot(d_chain, 2 * self.s_max + 2)

        p = self.linear_layer(
            torch.cat(
                [
                    a_rel_pos.float(),
                    a_rel_token.float(),
                    b_same_entity.unsqueeze(-1).float(),
                    a_rel_chain.float(),
                ],
                dim=-1,
            )
        )
        return p


class SingleConditioning(Module):
    """Single conditioning layer."""

    def __init__(
        self,
        sigma_data: float,
        token_s=384,
        dim_fourier=256,
        num_transitions=2,
        transition_expansion_factor=2,
        eps=1e-20,
    ):
        """Initialize the single conditioning layer.

        Parameters
        ----------
        sigma_data : float
            The data sigma.
        token_s : int, optional
            The single representation dimension, by default 384.
        dim_fourier : int, optional
            The fourier embeddings dimension, by default 256.
        num_transitions : int, optional
            The number of transitions layers, by default 2.
        transition_expansion_factor : int, optional
            The transition expansion factor, by default 2.
        eps : float, optional
            The epsilon value, by default 1e-20.

        """
        super().__init__()
        self.eps = eps
        self.sigma_data = sigma_data

        input_dim = (
            2 * token_s + 2 * const.num_tokens + 1 + len(const.pocket_contact_info)
        )
        self.norm_single = nn.LayerNorm(input_dim)
        self.single_embed = nn.Linear(input_dim, 2 * token_s)
        self.fourier_embed = FourierEmbedding(dim_fourier)
        self.norm_fourier = nn.LayerNorm(dim_fourier)
        self.fourier_to_single = LinearNoBias(dim_fourier, 2 * token_s)

        transitions = ModuleList([])
        for _ in range(num_transitions):
            transition = Transition(
                dim=2 * token_s, hidden=transition_expansion_factor * 2 * token_s
            )
            transitions.append(transition)

        self.transitions = transitions

    def forward(
        self,
        *,
        times,
        s_trunk,
        s_inputs,
    ):
        s = torch.cat((s_trunk, s_inputs), dim=-1)
        s = self.single_embed(self.norm_single(s))
        fourier_embed = self.fourier_embed(times)
        normed_fourier = self.norm_fourier(fourier_embed)
        fourier_to_single = self.fourier_to_single(normed_fourier)

        s = rearrange(fourier_to_single, "b d -> b 1 d") + s

        for transition in self.transitions:
            s = transition(s) + s

        return s, normed_fourier


class PairwiseConditioning(Module):
    """Pairwise conditioning layer."""

    def __init__(
        self,
        token_z,
        dim_token_rel_pos_feats,
        num_transitions=2,
        transition_expansion_factor=2,
    ):
        """Initialize the pairwise conditioning layer.

        Parameters
        ----------
        token_z : int
            The pair representation dimension.
        dim_token_rel_pos_feats : int
            The token relative position features dimension.
        num_transitions : int, optional
            The number of transitions layers, by default 2.
        transition_expansion_factor : int, optional
            The transition expansion factor, by default 2.

        """
        super().__init__()

        self.dim_pairwise_init_proj = nn.Sequential(
            nn.LayerNorm(token_z + dim_token_rel_pos_feats),
            LinearNoBias(token_z + dim_token_rel_pos_feats, token_z),
        )

        transitions = ModuleList([])
        for _ in range(num_transitions):
            transition = Transition(
                dim=token_z, hidden=transition_expansion_factor * token_z
            )
            transitions.append(transition)

        self.transitions = transitions

    def forward(
        self,
        z_trunk,
        token_rel_pos_feats,
    ):
        z = torch.cat((z_trunk, token_rel_pos_feats), dim=-1)
        z = self.dim_pairwise_init_proj(z)

        for transition in self.transitions:
            z = transition(z) + z

        return z


def get_indexing_matrix(K, W, H, device):
    assert W % 2 == 0
    assert H % (W // 2) == 0

    h = H // (W // 2)
    assert h % 2 == 0

    arange = torch.arange(2 * K, device=device)
    index = ((arange.unsqueeze(0) - arange.unsqueeze(1)) + h // 2).clamp(
        min=0, max=h + 1
    )
    index = index.view(K, 2, 2 * K)[:, 0, :]
    onehot = one_hot(index, num_classes=h + 2)[..., 1:-1].transpose(1, 0)
    return onehot.reshape(2 * K, h * K).float()


def single_to_keys(single, indexing_matrix, W, H):
    B, N, D = single.shape
    K = N // W
    single = single.view(B, 2 * K, W // 2, D)
    return torch.einsum("b j i d, j k -> b k i d", single, indexing_matrix).reshape(
        B, K, H, D
    )


class AtomAttentionEncoder(Module):
    """Atom attention encoder."""

    def __init__(
        self,
        atom_s,
        atom_z,
        token_s,
        token_z,
        atoms_per_window_queries,
        atoms_per_window_keys,
        atom_feature_dim,
        atom_encoder_depth=3,
        atom_encoder_heads=4,
        structure_prediction=True,
        activation_checkpointing=False,
    ):
        """Initialize the atom attention encoder.

        Parameters
        ----------
        atom_s : int
            The atom single representation dimension.
        atom_z : int
            The atom pair representation dimension.
        token_s : int
            The single representation dimension.
        token_z : int
            The pair representation dimension.
        atoms_per_window_queries : int
            The number of atoms per window for queries.
        atoms_per_window_keys : int
            The number of atoms per window for keys.
        atom_feature_dim : int
            The atom feature dimension.
        atom_encoder_depth : int, optional
            The number of transformer layers, by default 3.
        atom_encoder_heads : int, optional
            The number of transformer heads, by default 4.
        structure_prediction : bool, optional
            Whether it is used in the diffusion module, by default True.
        activation_checkpointing : bool, optional
            Whether to use activation checkpointing, by default False.

        """
        super().__init__()

        self.embed_atom_features = LinearNoBias(atom_feature_dim, atom_s)
        self.embed_atompair_ref_pos = LinearNoBias(3, atom_z)
        self.embed_atompair_ref_dist = LinearNoBias(1, atom_z)
        self.embed_atompair_mask = LinearNoBias(1, atom_z)
        self.atoms_per_window_queries = atoms_per_window_queries
        self.atoms_per_window_keys = atoms_per_window_keys

        self.structure_prediction = structure_prediction
        if structure_prediction:
            self.s_to_c_trans = nn.Sequential(
                nn.LayerNorm(token_s), LinearNoBias(token_s, atom_s)
            )
            init.final_init_(self.s_to_c_trans[1].weight)

            self.z_to_p_trans = nn.Sequential(
                nn.LayerNorm(token_z), LinearNoBias(token_z, atom_z)
            )
            init.final_init_(self.z_to_p_trans[1].weight)

            self.r_to_q_trans = LinearNoBias(10, atom_s)
            init.final_init_(self.r_to_q_trans.weight)

        self.c_to_p_trans_k = nn.Sequential(
            nn.ReLU(),
            LinearNoBias(atom_s, atom_z),
        )
        init.final_init_(self.c_to_p_trans_k[1].weight)

        self.c_to_p_trans_q = nn.Sequential(
            nn.ReLU(),
            LinearNoBias(atom_s, atom_z),
        )
        init.final_init_(self.c_to_p_trans_q[1].weight)

        self.p_mlp = nn.Sequential(
            nn.ReLU(),
            LinearNoBias(atom_z, atom_z),
            nn.ReLU(),
            LinearNoBias(atom_z, atom_z),
            nn.ReLU(),
            LinearNoBias(atom_z, atom_z),
        )
        init.final_init_(self.p_mlp[5].weight)

        self.atom_encoder = AtomTransformer(
            dim=atom_s,
            dim_single_cond=atom_s,
            dim_pairwise=atom_z,
            attn_window_queries=atoms_per_window_queries,
            attn_window_keys=atoms_per_window_keys,
            depth=atom_encoder_depth,
            heads=atom_encoder_heads,
            activation_checkpointing=activation_checkpointing,
        )

        self.atom_to_token_trans = nn.Sequential(
            LinearNoBias(atom_s, 2 * token_s if structure_prediction else token_s),
            nn.ReLU(),
        )

    def forward(
        self,
        feats,
        s_trunk=None,
        z=None,
        r=None,
        multiplicity=1,
        model_cache=None,
    ):
        B, N, _ = feats["ref_pos"].shape
        atom_mask = feats["atom_pad_mask"].bool()

        layer_cache = None
        if model_cache is not None:
            cache_prefix = "atomencoder"
            if cache_prefix not in model_cache:
                model_cache[cache_prefix] = {}
            layer_cache = model_cache[cache_prefix]

        if model_cache is None or len(layer_cache) == 0:
            # either model is not using the cache or it is the first time running it

            atom_ref_pos = feats["ref_pos"]
            atom_uid = feats["ref_space_uid"]
            atom_feats = torch.cat(
                [
                    atom_ref_pos,
                    feats["ref_charge"].unsqueeze(-1),
                    feats["atom_pad_mask"].unsqueeze(-1),
                    feats["ref_element"],
                    feats["ref_atom_name_chars"].reshape(B, N, 4 * 64),
                ],
                dim=-1,
            )

            c = self.embed_atom_features(atom_feats)

            # NOTE: we are already creating the windows to make it more efficient
            W, H = self.atoms_per_window_queries, self.atoms_per_window_keys
            B, N = c.shape[:2]
            K = N // W
            keys_indexing_matrix = get_indexing_matrix(K, W, H, c.device)
            to_keys = partial(
                single_to_keys, indexing_matrix=keys_indexing_matrix, W=W, H=H
            )

            atom_ref_pos_queries = atom_ref_pos.view(B, K, W, 1, 3)
            atom_ref_pos_keys = to_keys(atom_ref_pos).view(B, K, 1, H, 3)

            d = atom_ref_pos_keys - atom_ref_pos_queries
            d_norm = torch.sum(d * d, dim=-1, keepdim=True)
            d_norm = 1 / (1 + d_norm)

            atom_mask_queries = atom_mask.view(B, K, W, 1)
            atom_mask_keys = (
                to_keys(atom_mask.unsqueeze(-1).float()).view(B, K, 1, H).bool()
            )
            atom_uid_queries = atom_uid.view(B, K, W, 1)
            atom_uid_keys = (
                to_keys(atom_uid.unsqueeze(-1).float()).view(B, K, 1, H).long()
            )
            v = (
                (
                    atom_mask_queries
                    & atom_mask_keys
                    & (atom_uid_queries == atom_uid_keys)
                )
                .float()
                .unsqueeze(-1)
            )

            p = self.embed_atompair_ref_pos(d) * v
            p = p + self.embed_atompair_ref_dist(d_norm) * v
            p = p + self.embed_atompair_mask(v) * v

            q = c

            if self.structure_prediction:
                # run only in structure model not in initial encoding
                atom_to_token = feats["atom_to_token"].float()

                s_to_c = self.s_to_c_trans(s_trunk)
                s_to_c = torch.bmm(atom_to_token, s_to_c)
                c = c + s_to_c

                atom_to_token_queries = atom_to_token.view(
                    B, K, W, atom_to_token.shape[-1]
                )
                atom_to_token_keys = to_keys(atom_to_token)
                z_to_p = self.z_to_p_trans(z)
                z_to_p = torch.einsum(
                    "bijd,bwki,bwlj->bwkld",
                    z_to_p,
                    atom_to_token_queries,
                    atom_to_token_keys,
                )
                p = p + z_to_p

            p = p + self.c_to_p_trans_q(c.view(B, K, W, 1, c.shape[-1]))
            p = p + self.c_to_p_trans_k(to_keys(c).view(B, K, 1, H, c.shape[-1]))
            p = p + self.p_mlp(p)

            if model_cache is not None:
                layer_cache["q"] = q
                layer_cache["c"] = c
                layer_cache["p"] = p
                layer_cache["to_keys"] = to_keys

        else:
            q = layer_cache["q"]
            c = layer_cache["c"]
            p = layer_cache["p"]
            to_keys = layer_cache["to_keys"]

        if self.structure_prediction:
            # only here the multiplicity kicks in because we use the different positions r
            q = q.repeat_interleave(multiplicity, 0)
            r_input = torch.cat(
                [r, torch.zeros((B * multiplicity, N, 7)).to(r)],
                dim=-1,
            )
            r_to_q = self.r_to_q_trans(r_input)
            q = q + r_to_q

        c = c.repeat_interleave(multiplicity, 0)
        atom_mask = atom_mask.repeat_interleave(multiplicity, 0)

        q = self.atom_encoder(
            q=q,
            mask=atom_mask,
            c=c,
            p=p,
            multiplicity=multiplicity,
            to_keys=to_keys,
            model_cache=layer_cache,
        )

        q_to_a = self.atom_to_token_trans(q)
        atom_to_token = feats["atom_to_token"].float()
        atom_to_token = atom_to_token.repeat_interleave(multiplicity, 0)
        atom_to_token_mean = atom_to_token / (
            atom_to_token.sum(dim=1, keepdim=True) + 1e-6
        )
        a = torch.bmm(atom_to_token_mean.transpose(1, 2), q_to_a)

        return a, q, c, p, to_keys


class AtomAttentionDecoder(Module):
    """Atom attention decoder."""

    def __init__(
        self,
        atom_s,
        atom_z,
        token_s,
        attn_window_queries,
        attn_window_keys,
        atom_decoder_depth=3,
        atom_decoder_heads=4,
        activation_checkpointing=False,
    ):
        """Initialize the atom attention decoder.

        Parameters
        ----------
        atom_s : int
            The atom single representation dimension.
        atom_z : int
            The atom pair representation dimension.
        token_s : int
            The single representation dimension.
        attn_window_queries : int
            The number of atoms per window for queries.
        attn_window_keys : int
            The number of atoms per window for keys.
        atom_decoder_depth : int, optional
            The number of transformer layers, by default 3.
        atom_decoder_heads : int, optional
            The number of transformer heads, by default 4.
        activation_checkpointing : bool, optional
            Whether to use activation checkpointing, by default False.

        """
        super().__init__()

        self.a_to_q_trans = LinearNoBias(2 * token_s, atom_s)
        init.final_init_(self.a_to_q_trans.weight)

        self.atom_decoder = AtomTransformer(
            dim=atom_s,
            dim_single_cond=atom_s,
            dim_pairwise=atom_z,
            attn_window_queries=attn_window_queries,
            attn_window_keys=attn_window_keys,
            depth=atom_decoder_depth,
            heads=atom_decoder_heads,
            activation_checkpointing=activation_checkpointing,
        )

        self.atom_feat_to_atom_pos_update = nn.Sequential(
            nn.LayerNorm(atom_s), LinearNoBias(atom_s, 3)
        )
        init.final_init_(self.atom_feat_to_atom_pos_update[1].weight)

    def forward(
        self,
        a,
        q,
        c,
        p,
        feats,
        to_keys,
        multiplicity=1,
        model_cache=None,
    ):
        atom_mask = feats["atom_pad_mask"]
        atom_mask = atom_mask.repeat_interleave(multiplicity, 0)

        atom_to_token = feats["atom_to_token"].float()
        atom_to_token = atom_to_token.repeat_interleave(multiplicity, 0)

        a_to_q = self.a_to_q_trans(a)
        a_to_q = torch.bmm(atom_to_token, a_to_q)
        q = q + a_to_q

        layer_cache = None
        if model_cache is not None:
            cache_prefix = "atomdecoder"
            if cache_prefix not in model_cache:
                model_cache[cache_prefix] = {}
            layer_cache = model_cache[cache_prefix]

        q = self.atom_decoder(
            q=q,
            mask=atom_mask,
            c=c,
            p=p,
            multiplicity=multiplicity,
            to_keys=to_keys,
            model_cache=layer_cache,
        )

        r_update = self.atom_feat_to_atom_pos_update(q)
        return r_update
