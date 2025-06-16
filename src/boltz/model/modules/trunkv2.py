import torch
from torch import Tensor, nn
from torch.nn.functional import one_hot

from boltz.data import const
from boltz.model.layers.outer_product_mean import OuterProductMean
from boltz.model.layers.pair_averaging import PairWeightedAveraging
from boltz.model.layers.pairformer import (
    PairformerNoSeqLayer,
    PairformerNoSeqModule,
    get_dropout_mask,
)
from boltz.model.layers.transition import Transition
from boltz.model.modules.encodersv2 import (
    AtomAttentionEncoder,
    AtomEncoder,
    FourierEmbedding,
)


class ContactConditioning(nn.Module):
    def __init__(self, token_z: int, cutoff_min: float, cutoff_max: float):
        super().__init__()

        self.fourier_embedding = FourierEmbedding(token_z)
        self.encoder = nn.Linear(
            token_z + len(const.contact_conditioning_info) - 1, token_z
        )
        self.encoding_unspecified = nn.Parameter(torch.zeros(token_z))
        self.encoding_unselected = nn.Parameter(torch.zeros(token_z))
        self.cutoff_min = cutoff_min
        self.cutoff_max = cutoff_max

    def forward(self, feats):
        assert const.contact_conditioning_info["UNSPECIFIED"] == 0
        assert const.contact_conditioning_info["UNSELECTED"] == 1
        contact_conditioning = feats["contact_conditioning"][:, :, :, 2:]
        contact_threshold = feats["contact_threshold"]
        contact_threshold_normalized = (contact_threshold - self.cutoff_min) / (
            self.cutoff_max - self.cutoff_min
        )
        contact_threshold_fourier = self.fourier_embedding(
            contact_threshold_normalized.flatten()
        ).reshape(contact_threshold_normalized.shape + (-1,))

        contact_conditioning = torch.cat(
            [
                contact_conditioning,
                contact_threshold_normalized.unsqueeze(-1),
                contact_threshold_fourier,
            ],
            dim=-1,
        )
        contact_conditioning = self.encoder(contact_conditioning)

        contact_conditioning = (
            contact_conditioning
            * (
                1
                - feats["contact_conditioning"][:, :, :, 0:2].sum(dim=-1, keepdim=True)
            )
            + self.encoding_unspecified * feats["contact_conditioning"][:, :, :, 0:1]
            + self.encoding_unselected * feats["contact_conditioning"][:, :, :, 1:2]
        )
        return contact_conditioning


class InputEmbedder(nn.Module):
    def __init__(
        self,
        atom_s: int,
        atom_z: int,
        token_s: int,
        token_z: int,
        atoms_per_window_queries: int,
        atoms_per_window_keys: int,
        atom_feature_dim: int,
        atom_encoder_depth: int,
        atom_encoder_heads: int,
        activation_checkpointing: bool = False,
        add_method_conditioning: bool = False,
        add_modified_flag: bool = False,
        add_cyclic_flag: bool = False,
        add_mol_type_feat: bool = False,
        use_no_atom_char: bool = False,
        use_atom_backbone_feat: bool = False,
        use_residue_feats_atoms: bool = False,
    ) -> None:
        """Initialize the input embedder.

        Parameters
        ----------
        atom_s : int
            The atom embedding size.
        atom_z : int
            The atom pairwise embedding size.
        token_s : int
            The token embedding size.

        """
        super().__init__()
        self.token_s = token_s
        self.add_method_conditioning = add_method_conditioning
        self.add_modified_flag = add_modified_flag
        self.add_cyclic_flag = add_cyclic_flag
        self.add_mol_type_feat = add_mol_type_feat

        self.atom_encoder = AtomEncoder(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            token_z=token_z,
            atoms_per_window_queries=atoms_per_window_queries,
            atoms_per_window_keys=atoms_per_window_keys,
            atom_feature_dim=atom_feature_dim,
            structure_prediction=False,
            use_no_atom_char=use_no_atom_char,
            use_atom_backbone_feat=use_atom_backbone_feat,
            use_residue_feats_atoms=use_residue_feats_atoms,
        )

        self.atom_enc_proj_z = nn.Sequential(
            nn.LayerNorm(atom_z),
            nn.Linear(atom_z, atom_encoder_depth * atom_encoder_heads, bias=False),
        )

        self.atom_attention_encoder = AtomAttentionEncoder(
            atom_s=atom_s,
            token_s=token_s,
            atoms_per_window_queries=atoms_per_window_queries,
            atoms_per_window_keys=atoms_per_window_keys,
            atom_encoder_depth=atom_encoder_depth,
            atom_encoder_heads=atom_encoder_heads,
            structure_prediction=False,
            activation_checkpointing=activation_checkpointing,
        )

        self.res_type_encoding = nn.Linear(const.num_tokens, token_s, bias=False)
        self.msa_profile_encoding = nn.Linear(const.num_tokens + 1, token_s, bias=False)

        if add_method_conditioning:
            self.method_conditioning_init = nn.Embedding(
                const.num_method_types, token_s
            )
            self.method_conditioning_init.weight.data.fill_(0)
        if add_modified_flag:
            self.modified_conditioning_init = nn.Embedding(2, token_s)
            self.modified_conditioning_init.weight.data.fill_(0)
        if add_cyclic_flag:
            self.cyclic_conditioning_init = nn.Linear(1, token_s, bias=False)
            self.cyclic_conditioning_init.weight.data.fill_(0)
        if add_mol_type_feat:
            self.mol_type_conditioning_init = nn.Embedding(
                len(const.chain_type_ids), token_s
            )
            self.mol_type_conditioning_init.weight.data.fill_(0)

    def forward(self, feats: dict[str, Tensor], affinity: bool = False) -> Tensor:
        """Perform the forward pass.

        Parameters
        ----------
        feats : dict[str, Tensor]
            Input features

        Returns
        -------
        Tensor
            The embedded tokens.

        """
        # Load relevant features
        res_type = feats["res_type"].float()
        if affinity:
            profile = feats["profile_affinity"]
            deletion_mean = feats["deletion_mean_affinity"].unsqueeze(-1)
        else:
            profile = feats["profile"]
            deletion_mean = feats["deletion_mean"].unsqueeze(-1)

        # Compute input embedding
        q, c, p, to_keys = self.atom_encoder(feats)
        atom_enc_bias = self.atom_enc_proj_z(p)
        a, _, _, _ = self.atom_attention_encoder(
            feats=feats,
            q=q,
            c=c,
            atom_enc_bias=atom_enc_bias,
            to_keys=to_keys,
        )

        s = (
            a
            + self.res_type_encoding(res_type)
            + self.msa_profile_encoding(torch.cat([profile, deletion_mean], dim=-1))
        )

        if self.add_method_conditioning:
            s = s + self.method_conditioning_init(feats["method_feature"])
        if self.add_modified_flag:
            s = s + self.modified_conditioning_init(feats["modified"])
        if self.add_cyclic_flag:
            cyclic = feats["cyclic_period"].clamp(max=1.0).unsqueeze(-1)
            s = s + self.cyclic_conditioning_init(cyclic)
        if self.add_mol_type_feat:
            s = s + self.mol_type_conditioning_init(feats["mol_type"])

        return s


class TemplateModule(nn.Module):
    """Template module."""

    def __init__(
        self,
        token_z: int,
        template_dim: int,
        template_blocks: int,
        dropout: float = 0.25,
        pairwise_head_width: int = 32,
        pairwise_num_heads: int = 4,
        post_layer_norm: bool = False,
        activation_checkpointing: bool = False,
        min_dist: float = 3.25,
        max_dist: float = 50.75,
        num_bins: int = 38,
        **kwargs,
    ) -> None:
        """Initialize the template module.

        Parameters
        ----------
        token_z : int
            The token pairwise embedding size.

        """
        super().__init__()
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.num_bins = num_bins
        self.relu = nn.ReLU()
        self.z_norm = nn.LayerNorm(token_z)
        self.v_norm = nn.LayerNorm(template_dim)
        self.z_proj = nn.Linear(token_z, template_dim, bias=False)
        self.a_proj = nn.Linear(
            const.num_tokens * 2 + num_bins + 5,
            template_dim,
            bias=False,
        )
        self.u_proj = nn.Linear(template_dim, token_z, bias=False)
        self.pairformer = PairformerNoSeqModule(
            template_dim,
            num_blocks=template_blocks,
            dropout=dropout,
            pairwise_head_width=pairwise_head_width,
            pairwise_num_heads=pairwise_num_heads,
            post_layer_norm=post_layer_norm,
            activation_checkpointing=activation_checkpointing,
        )

    def forward(
        self,
        z: Tensor,
        feats: dict[str, Tensor],
        pair_mask: Tensor,
        use_kernels: bool = False,
    ) -> Tensor:
        """Perform the forward pass.

        Parameters
        ----------
        z : Tensor
            The pairwise embeddings
        feats : dict[str, Tensor]
            Input features
        pair_mask : Tensor
            The pair mask

        Returns
        -------
        Tensor
            The updated pairwise embeddings.

        """
        # Load relevant features
        asym_id = feats["asym_id"]
        res_type = feats["template_restype"]
        frame_rot = feats["template_frame_rot"]
        frame_t = feats["template_frame_t"]
        frame_mask = feats["template_mask_frame"]
        cb_coords = feats["template_cb"]
        ca_coords = feats["template_ca"]
        cb_mask = feats["template_mask_cb"]
        template_mask = feats["template_mask"].any(dim=2).float()
        num_templates = template_mask.sum(dim=1)
        num_templates = num_templates.clamp(min=1)

        # Compute pairwise masks
        b_cb_mask = cb_mask[:, :, :, None] * cb_mask[:, :, None, :]
        b_frame_mask = frame_mask[:, :, :, None] * frame_mask[:, :, None, :]

        b_cb_mask = b_cb_mask[..., None]
        b_frame_mask = b_frame_mask[..., None]

        # Compute asym mask, template features only attend within the same chain
        B, T = res_type.shape[:2]  # noqa: N806
        asym_mask = (asym_id[:, :, None] == asym_id[:, None, :]).float()
        asym_mask = asym_mask[:, None].expand(-1, T, -1, -1)

        # Compute template features
        with torch.autocast(device_type="cuda", enabled=False):
            # Compute distogram
            cb_dists = torch.cdist(cb_coords, cb_coords)
            boundaries = torch.linspace(self.min_dist, self.max_dist, self.num_bins - 1)
            boundaries = boundaries.to(cb_dists.device)
            distogram = (cb_dists[..., None] > boundaries).sum(dim=-1).long()
            distogram = one_hot(distogram, num_classes=self.num_bins)

            # Compute unit vector in each frame
            frame_rot = frame_rot.unsqueeze(2).transpose(-1, -2)
            frame_t = frame_t.unsqueeze(2).unsqueeze(-1)
            ca_coords = ca_coords.unsqueeze(3).unsqueeze(-1)
            vector = torch.matmul(frame_rot, (ca_coords - frame_t))
            norm = torch.norm(vector, dim=-1, keepdim=True)
            unit_vector = torch.where(norm > 0, vector / norm, torch.zeros_like(vector))
            unit_vector = unit_vector.squeeze(-1)

            # Concatenate input features
            a_tij = [distogram, b_cb_mask, unit_vector, b_frame_mask]
            a_tij = torch.cat(a_tij, dim=-1)
            a_tij = a_tij * asym_mask.unsqueeze(-1)

            res_type_i = res_type[:, :, :, None]
            res_type_j = res_type[:, :, None, :]
            res_type_i = res_type_i.expand(-1, -1, -1, res_type.size(2), -1)
            res_type_j = res_type_j.expand(-1, -1, res_type.size(2), -1, -1)
            a_tij = torch.cat([a_tij, res_type_i, res_type_j], dim=-1)
            a_tij = self.a_proj(a_tij)

        # Expand mask
        pair_mask = pair_mask[:, None].expand(-1, T, -1, -1)
        pair_mask = pair_mask.reshape(B * T, *pair_mask.shape[2:])

        # Compute input projections
        v = self.z_proj(self.z_norm(z[:, None])) + a_tij
        v = v.view(B * T, *v.shape[2:])
        v = v + self.pairformer(v, pair_mask, use_kernels=use_kernels)
        v = self.v_norm(v)
        v = v.view(B, T, *v.shape[1:])

        # Aggregate templates
        template_mask = template_mask[:, :, None, None, None]
        num_templates = num_templates[:, None, None, None]
        u = (v * template_mask).sum(dim=1) / num_templates.to(v)

        # Compute output projection
        u = self.u_proj(self.relu(u))
        return u


class TemplateV2Module(nn.Module):
    """Template module."""

    def __init__(
        self,
        token_z: int,
        template_dim: int,
        template_blocks: int,
        dropout: float = 0.25,
        pairwise_head_width: int = 32,
        pairwise_num_heads: int = 4,
        post_layer_norm: bool = False,
        activation_checkpointing: bool = False,
        min_dist: float = 3.25,
        max_dist: float = 50.75,
        num_bins: int = 38,
        **kwargs,
    ) -> None:
        """Initialize the template module.

        Parameters
        ----------
        token_z : int
            The token pairwise embedding size.

        """
        super().__init__()
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.num_bins = num_bins
        self.relu = nn.ReLU()
        self.z_norm = nn.LayerNorm(token_z)
        self.v_norm = nn.LayerNorm(template_dim)
        self.z_proj = nn.Linear(token_z, template_dim, bias=False)
        self.a_proj = nn.Linear(
            const.num_tokens * 2 + num_bins + 5,
            template_dim,
            bias=False,
        )
        self.u_proj = nn.Linear(template_dim, token_z, bias=False)
        self.pairformer = PairformerNoSeqModule(
            template_dim,
            num_blocks=template_blocks,
            dropout=dropout,
            pairwise_head_width=pairwise_head_width,
            pairwise_num_heads=pairwise_num_heads,
            post_layer_norm=post_layer_norm,
            activation_checkpointing=activation_checkpointing,
        )

    def forward(
        self,
        z: Tensor,
        feats: dict[str, Tensor],
        pair_mask: Tensor,
        use_kernels: bool = False,
    ) -> Tensor:
        """Perform the forward pass.

        Parameters
        ----------
        z : Tensor
            The pairwise embeddings
        feats : dict[str, Tensor]
            Input features
        pair_mask : Tensor
            The pair mask

        Returns
        -------
        Tensor
            The updated pairwise embeddings.

        """
        # Load relevant features
        res_type = feats["template_restype"]
        frame_rot = feats["template_frame_rot"]
        frame_t = feats["template_frame_t"]
        frame_mask = feats["template_mask_frame"]
        cb_coords = feats["template_cb"]
        ca_coords = feats["template_ca"]
        cb_mask = feats["template_mask_cb"]
        visibility_ids = feats["visibility_ids"]
        template_mask = feats["template_mask"].any(dim=2).float()
        num_templates = template_mask.sum(dim=1)
        num_templates = num_templates.clamp(min=1)

        # Compute pairwise masks
        b_cb_mask = cb_mask[:, :, :, None] * cb_mask[:, :, None, :]
        b_frame_mask = frame_mask[:, :, :, None] * frame_mask[:, :, None, :]

        b_cb_mask = b_cb_mask[..., None]
        b_frame_mask = b_frame_mask[..., None]

        # Compute asym mask, template features only attend within the same chain
        B, T = res_type.shape[:2]  # noqa: N806
        tmlp_pair_mask = (
            visibility_ids[:, :, :, None] == visibility_ids[:, :, None, :]
        ).float()

        # Compute template features
        with torch.autocast(device_type="cuda", enabled=False):
            # Compute distogram
            cb_dists = torch.cdist(cb_coords, cb_coords)
            boundaries = torch.linspace(self.min_dist, self.max_dist, self.num_bins - 1)
            boundaries = boundaries.to(cb_dists.device)
            distogram = (cb_dists[..., None] > boundaries).sum(dim=-1).long()
            distogram = one_hot(distogram, num_classes=self.num_bins)

            # Compute unit vector in each frame
            frame_rot = frame_rot.unsqueeze(2).transpose(-1, -2)
            frame_t = frame_t.unsqueeze(2).unsqueeze(-1)
            ca_coords = ca_coords.unsqueeze(3).unsqueeze(-1)
            vector = torch.matmul(frame_rot, (ca_coords - frame_t))
            norm = torch.norm(vector, dim=-1, keepdim=True)
            unit_vector = torch.where(norm > 0, vector / norm, torch.zeros_like(vector))
            unit_vector = unit_vector.squeeze(-1)

            # Concatenate input features
            a_tij = [distogram, b_cb_mask, unit_vector, b_frame_mask]
            a_tij = torch.cat(a_tij, dim=-1)
            a_tij = a_tij * tmlp_pair_mask.unsqueeze(-1)

            res_type_i = res_type[:, :, :, None]
            res_type_j = res_type[:, :, None, :]
            res_type_i = res_type_i.expand(-1, -1, -1, res_type.size(2), -1)
            res_type_j = res_type_j.expand(-1, -1, res_type.size(2), -1, -1)
            a_tij = torch.cat([a_tij, res_type_i, res_type_j], dim=-1)
            a_tij = self.a_proj(a_tij)

        # Expand mask
        pair_mask = pair_mask[:, None].expand(-1, T, -1, -1)
        pair_mask = pair_mask.reshape(B * T, *pair_mask.shape[2:])

        # Compute input projections
        v = self.z_proj(self.z_norm(z[:, None])) + a_tij
        v = v.view(B * T, *v.shape[2:])
        v = v + self.pairformer(v, pair_mask, use_kernels=use_kernels)
        v = self.v_norm(v)
        v = v.view(B, T, *v.shape[1:])

        # Aggregate templates
        template_mask = template_mask[:, :, None, None, None]
        num_templates = num_templates[:, None, None, None]
        u = (v * template_mask).sum(dim=1) / num_templates.to(v)

        # Compute output projection
        u = self.u_proj(self.relu(u))
        return u


class MSAModule(nn.Module):
    """MSA module."""

    def __init__(
        self,
        msa_s: int,
        token_z: int,
        token_s: int,
        msa_blocks: int,
        msa_dropout: float,
        z_dropout: float,
        pairwise_head_width: int = 32,
        pairwise_num_heads: int = 4,
        activation_checkpointing: bool = False,
        use_paired_feature: bool = True,
        subsample_msa: bool = False,
        num_subsampled_msa: int = 1024,
        **kwargs,
    ) -> None:
        """Initialize the MSA module.

        Parameters
        ----------
        token_z : int
            The token pairwise embedding size.

        """
        super().__init__()
        self.msa_blocks = msa_blocks
        self.msa_dropout = msa_dropout
        self.z_dropout = z_dropout
        self.use_paired_feature = use_paired_feature
        self.activation_checkpointing = activation_checkpointing
        self.subsample_msa = subsample_msa
        self.num_subsampled_msa = num_subsampled_msa

        self.s_proj = nn.Linear(token_s, msa_s, bias=False)
        self.msa_proj = nn.Linear(
            const.num_tokens + 2 + int(use_paired_feature),
            msa_s,
            bias=False,
        )
        self.layers = nn.ModuleList()
        for i in range(msa_blocks):
            self.layers.append(
                MSALayer(
                    msa_s,
                    token_z,
                    msa_dropout,
                    z_dropout,
                    pairwise_head_width,
                    pairwise_num_heads,
                )
            )

    def forward(
        self,
        z: Tensor,
        emb: Tensor,
        feats: dict[str, Tensor],
        use_kernels: bool = False,
    ) -> Tensor:
        """Perform the forward pass.

        Parameters
        ----------
        z : Tensor
            The pairwise embeddings
        emb : Tensor
            The input embeddings
        feats : dict[str, Tensor]
            Input features
        use_kernels: bool
            Whether to use kernels for triangular updates

        Returns
        -------
        Tensor
            The output pairwise embeddings.

        """
        # Set chunk sizes
        if not self.training:
            if z.shape[1] > const.chunk_size_threshold:
                chunk_heads_pwa = True
                chunk_size_transition_z = 64
                chunk_size_transition_msa = 32
                chunk_size_outer_product = 4
                chunk_size_tri_attn = 128
            else:
                chunk_heads_pwa = False
                chunk_size_transition_z = None
                chunk_size_transition_msa = None
                chunk_size_outer_product = None
                chunk_size_tri_attn = 512
        else:
            chunk_heads_pwa = False
            chunk_size_transition_z = None
            chunk_size_transition_msa = None
            chunk_size_outer_product = None
            chunk_size_tri_attn = None

        # Load relevant features
        msa = feats["msa"]
        msa = torch.nn.functional.one_hot(msa, num_classes=const.num_tokens)
        has_deletion = feats["has_deletion"].unsqueeze(-1)
        deletion_value = feats["deletion_value"].unsqueeze(-1)
        is_paired = feats["msa_paired"].unsqueeze(-1)
        msa_mask = feats["msa_mask"]
        token_mask = feats["token_pad_mask"].float()
        token_mask = token_mask[:, :, None] * token_mask[:, None, :]

        # Compute MSA embeddings
        if self.use_paired_feature:
            m = torch.cat([msa, has_deletion, deletion_value, is_paired], dim=-1)
        else:
            m = torch.cat([msa, has_deletion, deletion_value], dim=-1)

        # Subsample the MSA
        if self.subsample_msa:
            msa_indices = torch.randperm(msa.shape[1])[: self.num_subsampled_msa]
            m = m[:, msa_indices]
            msa_mask = msa_mask[:, msa_indices]

        # Compute input projections
        m = self.msa_proj(m)
        m = m + self.s_proj(emb).unsqueeze(1)

        # Perform MSA blocks
        for i in range(self.msa_blocks):
            if self.activation_checkpointing and self.training:
                z, m = torch.utils.checkpoint.checkpoint(
                    self.layers[i],
                    z,
                    m,
                    token_mask,
                    msa_mask,
                    chunk_heads_pwa,
                    chunk_size_transition_z,
                    chunk_size_transition_msa,
                    chunk_size_outer_product,
                    chunk_size_tri_attn,
                    use_kernels,
                )
            else:
                z, m = self.layers[i](
                    z,
                    m,
                    token_mask,
                    msa_mask,
                    chunk_heads_pwa,
                    chunk_size_transition_z,
                    chunk_size_transition_msa,
                    chunk_size_outer_product,
                    chunk_size_tri_attn,
                    use_kernels,
                )
        return z


class MSALayer(nn.Module):
    """MSA module."""

    def __init__(
        self,
        msa_s: int,
        token_z: int,
        msa_dropout: float,
        z_dropout: float,
        pairwise_head_width: int = 32,
        pairwise_num_heads: int = 4,
    ) -> None:
        """Initialize the MSA module.

        Parameters
        ----------
        token_z : int
            The token pairwise embedding size.

        """
        super().__init__()
        self.msa_dropout = msa_dropout
        self.msa_transition = Transition(dim=msa_s, hidden=msa_s * 4)
        self.pair_weighted_averaging = PairWeightedAveraging(
            c_m=msa_s,
            c_z=token_z,
            c_h=32,
            num_heads=8,
        )

        self.pairformer_layer = PairformerNoSeqLayer(
            token_z=token_z,
            dropout=z_dropout,
            pairwise_head_width=pairwise_head_width,
            pairwise_num_heads=pairwise_num_heads,
        )
        self.outer_product_mean = OuterProductMean(
            c_in=msa_s,
            c_hidden=32,
            c_out=token_z,
        )

    def forward(
        self,
        z: Tensor,
        m: Tensor,
        token_mask: Tensor,
        msa_mask: Tensor,
        chunk_heads_pwa: bool = False,
        chunk_size_transition_z: int = None,
        chunk_size_transition_msa: int = None,
        chunk_size_outer_product: int = None,
        chunk_size_tri_attn: int = None,
        use_kernels: bool = False,
    ) -> tuple[Tensor, Tensor]:
        """Perform the forward pass.

        Parameters
        ----------
        z : Tensor
            The pairwise embeddings
        emb : Tensor
            The input embeddings
        feats : dict[str, Tensor]
            Input features

        Returns
        -------
        Tensor
            The output pairwise embeddings.

        """
        # Communication to MSA stack
        msa_dropout = get_dropout_mask(self.msa_dropout, m, self.training)
        m = m + msa_dropout * self.pair_weighted_averaging(
            m, z, token_mask, chunk_heads_pwa
        )
        m = m + self.msa_transition(m, chunk_size_transition_msa)

        z = z + self.outer_product_mean(m, msa_mask, chunk_size_outer_product)

        # Compute pairwise stack
        z = self.pairformer_layer(
            z, token_mask, chunk_size_tri_attn, use_kernels=use_kernels
        )

        return z, m


class BFactorModule(nn.Module):
    """BFactor Module."""

    def __init__(self, token_s: int, num_bins: int) -> None:
        """Initialize the bfactor module.

        Parameters
        ----------
        token_s : int
            The token embedding size.

        """
        super().__init__()
        self.bfactor = nn.Linear(token_s, num_bins)
        self.num_bins = num_bins

    def forward(self, s: Tensor) -> Tensor:
        """Perform the forward pass.

        Parameters
        ----------
        s : Tensor
            The sequence embeddings

        Returns
        -------
        Tensor
            The predicted bfactor histogram.

        """
        return self.bfactor(s)


class DistogramModule(nn.Module):
    """Distogram Module."""

    def __init__(self, token_z: int, num_bins: int, num_distograms: int = 1) -> None:
        """Initialize the distogram module.

        Parameters
        ----------
        token_z : int
            The token pairwise embedding size.

        """
        super().__init__()
        self.distogram = nn.Linear(token_z, num_distograms * num_bins)
        self.num_distograms = num_distograms
        self.num_bins = num_bins

    def forward(self, z: Tensor) -> Tensor:
        """Perform the forward pass.

        Parameters
        ----------
        z : Tensor
            The pairwise embeddings

        Returns
        -------
        Tensor
            The predicted distogram.

        """
        z = z + z.transpose(1, 2)
        return self.distogram(z).reshape(
            z.shape[0], z.shape[1], z.shape[2], self.num_distograms, self.num_bins
        )
