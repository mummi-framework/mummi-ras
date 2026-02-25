#!/usr/bin/env python3
#
# Copyright (c) 2025, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# @authors
#           Tim Hsu <hsu16@llnl.gov>
# ------------------------------------------------------------------------------

import torch
import numpy as np
from torch import nn, Tensor
from .basis import GaussianRandomFourierFeatures
from .diffusion import VariancePreservingDiffuser
from .transformer import TransformerEncoderContinuousRoPE, TransformerEncoderLayerContinuousRoPE
from torch.nn.utils.rnn import pad_sequence
from scipy.spatial.transform import Rotation

# TODO: switch to using only 'lightning'
try:
    import lightning as L
    from lightning.pytorch.utilities import grad_norm
except ImportError:
    import pytorch_lightning as L
    from pytorch_lightning.utilities import grad_norm


class TransformerBackbone(nn.Module):
    def __init__(self, pad_dim: int, dim: int, ff_dim: int, num_heads: int, num_layers: int, ucg_idx: list[list[int]], dropout: float):
        super().__init__()
        self.pad_dim    = pad_dim
        self.dim        = dim
        self.ff_dim     = ff_dim
        self.num_heads  = num_heads
        self.num_layers = num_layers
        self.ucg_idx    = ucg_idx
        self.dropout    = dropout
        
        self.embed_cg_patch = nn.Linear(pad_dim*3, dim//2)
        self.embed_ucg_token = nn.Embedding(250, dim//2)
        self.embed_time = GaussianRandomFourierFeatures(dim, input_dim=1)
        self.transformer_encoder = TransformerEncoderContinuousRoPE(
            encoder_layer = TransformerEncoderLayerContinuousRoPE(dim, num_heads, ff_dim, dropout, activation=nn.functional.silu, rope_dim=3),
            num_layers = num_layers,
            norm = nn.LayerNorm(dim),
        )
        self.lin_out = nn.Linear(dim, pad_dim*3)

    def forward(self, ucg_pos: Tensor, cg_disp: Tensor, t: Tensor) -> Tensor:
        B = ucg_pos.size(0)
        device = ucg_pos.device

        # Split cg_disp into num_ucg (250) tokens as nested tensor: [B, num_cg, 3] -> [[B, ?, 3], ...] -> [num_ucg, B, jagged, 3]
        # nested_cg_disp = torch.nested.nested_tensor(
        #     [cg_disp[:, indices, :] for indices in self.ucg_idx],
        #     layout = torch.jagged,
        # )

        # Pad the patches: [num_ucg, B, jagged, 3] -> [num_ucg, B, pad_dim, 3] -> [B, num_ucg, pad_dim, 3] -> [B, num_ucg, pad_dim*3]
        # patch_tokens = nested_cg_disp.to_padded_tensor(padding=0).transpose(0, 1).flatten(2)

        # Split cg_disp into num_ucg (250) patches: [B, num_cg, 3] -> [[B, ?, 3], ...] -> [[?, B, 3], ...]
        nested_cg_disp = [cg_disp[:, indices, :].transpose(0, 1) for indices in self.ucg_idx]

        # Pad the patch: [[?, B, 3], ...] -> [pad_dim, num_ucg, B, 3] -> [B, num_ucg, pad_dim, 3] -> [B, num_ucg, pad_dim*3]
        patch_tokens = pad_sequence(nested_cg_disp).transpose(0, 2).flatten(2)

        # Embed patch tokens: [B, num_ucg, pad_dim*3] -> [B, num_ucg, D/2]
        patch_tokens = self.embed_cg_patch(patch_tokens)

        # Concat UCG bead type information to patches -> [B, num_ucg, D] 
        patch_tokens = torch.cat([patch_tokens, self.embed_ucg_token(torch.arange(250, device=device).expand(B, -1))], dim=-1)

        # Append time token to patch tokens -> [B, num_ucg+1, D]
        all_tokens = torch.cat([patch_tokens, self.embed_time(t)], dim=1)

        # Append all-zero, dummy position for t token to UCG positions -> [B, num_ucg+1, 3]
        pos = torch.cat([ucg_pos, torch.zeros(B, 1, 3, device=device)], dim=1)

        # Feed all tokens and positions into transformer encoder -> [B, num_ucg+1, D]
        output_patches = self.transformer_encoder(all_tokens, pos)

        # Decode transformer outputs: [B, num_ucg+1, D] -> [B, num_ucg, D] -> [B, num_ucg, pad_dim*3] -> [B, num_ucg, pad_dim, 3]
        output_reshaped = self.lin_out(output_patches[:, :-1, :]).view(B, -1, self.pad_dim, 3)

        # Convert to CG-shaped output: [B, num_ucg, pad_dim, 3] -> [B, num_cg, 3]
        cg_disp_out = torch.zeros_like(cg_disp)
        for n, indices in enumerate(self.ucg_idx):
            ucg_size = len(indices)
            cg_disp_out[:, indices, :] = output_reshaped[:, n, :ucg_size, :]
        return cg_disp_out


class UCG2CG_Generator(L.LightningModule):
    def __init__(self, pad_dim: int, dim: int, ff_dim: int, num_heads: int, num_layers: int, ucg_index_file: int, dropout: float, learn_rate: float, T_max: int):
        super().__init__()
        self.save_hyperparameters()

        # Read UCG index file, e.g., ucg_idx: [[2,4,0], [1,3], ...]
        self.ucg_idx = []
        with open(ucg_index_file, 'r') as f:
            for line in f:
                self.ucg_idx.append([int(s) for s in line.rstrip().split(',')])

        # Get scatter indices from ucg_idx, e.g., scatter_idx = [0, 1, 0, 1, 0, ...]
        num_cgs = sum(len(indices) for indices in self.ucg_idx)
        self.scatter_idx = torch.zeros(num_cgs, dtype=torch.long)
        for i, indices in enumerate(self.ucg_idx):
            self.scatter_idx[indices] = i

        # Core model
        self.model = TransformerBackbone(
            pad_dim    = pad_dim,
            dim        = dim,
            ff_dim     = ff_dim,
            num_heads  = num_heads,
            num_layers = num_layers,
            ucg_idx    = self.ucg_idx,
            dropout    = dropout,
        )

        # Data diffuser
        self.diffuser = VariancePreservingDiffuser()

        # Keep track of these parameters
        self.learn_rate = learn_rate
        self.T_max = T_max

    def forward(self, ucg_pos: Tensor, cg_disp: Tensor, t: Tensor) -> Tensor:
        return self.model(ucg_pos, cg_disp, t)

    def _random_rotate(self, batch):
        B = batch['ucg_pos'].size(0)
        theta = np.random.uniform(low=0, high=2*np.pi, size=(B, 1))
        R = Rotation.from_rotvec(np.array([[0., 0., 1.]]) * theta).as_matrix()  # random rotation around z-axis
        R = torch.tensor(R, dtype=torch.float, device=batch['ucg_pos'].device)
        batch['ucg_pos'] @= R  # [B, num_ucg, 3] @ [B, 3, 3] -> [B, num_ucg, 3] 
        batch['cg_disp'] @= R  # [B, num_cg,  3] @ [B, 3, 3] -> [B, nu_ucg,  3] 
        return batch

    def _get_loss(self, batch):
        ucg_pos, cg_disp = batch['ucg_pos'], batch['cg_disp']
        t = torch.empty(ucg_pos.size(0), 1, 1, device=ucg_pos.device).uniform_(self.diffuser.t_min, self.diffuser.t_max)
        noisy_cg_disp, eps = self.diffuser.forward_noise(cg_disp, t)
        pred_eps = self.model(ucg_pos, noisy_cg_disp, t)
        return nn.functional.mse_loss(pred_eps, eps)

    def training_step(self, batch, batch_idx):
        batch = self._random_rotate(batch)
        loss = self._get_loss(batch)
        self.log('train_loss', loss, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        loss = self._get_loss(batch)
        self.log('valid_loss', loss, prog_bar=True, sync_dist=True)
        self.log('hp_metric',  loss, sync_dist=True)

    def configure_optimizers(self):
        optimizer = torch.optim.AdamW(self.parameters(), lr=self.learn_rate)
        scheduler1 = torch.optim.lr_scheduler.LinearLR(optimizer, start_factor=0.1, total_iters=5)
        scheduler2 = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=self.T_max)
        scheduler = torch.optim.lr_scheduler.SequentialLR(
            optimizer, [scheduler1, scheduler2], milestones=[5]
        )
        return {'optimizer': optimizer, 'lr_scheduler': scheduler}
    
    def on_before_optimizer_step(self, optimizer):
        # Compute the 2-norm for each layer
        # If using mixed precision, the gradients are already unscaled here
        norms = grad_norm(self.model, norm_type=2)
        self.log_dict(norms)

    def denoise_by_sde(self, cg_disp: Tensor, ucg_pos: Tensor, t_span: tuple[float, float] = (0.999, 0.001), steps: int = 128, return_last: bool = False) -> Tensor:
        B = ucg_pos.size(0)
        device = ucg_pos.device

        t_start, t_end = t_span
        ts = torch.linspace(t_start, t_end, steps+1, device=device)
        f, g, g2 = self.diffuser.f, self.diffuser.g, self.diffuser.g2
        sigma = self.diffuser.sigma
        
        x = cg_disp
        x_traj = [x.clone()]  # [[B, num_cg, 3]]
        for i, t in enumerate(ts[1:]):
            dt = ts[i+1] - ts[i]
            eps = dt.abs().sqrt() * torch.randn_like(x)
            pred_eps = self.model(ucg_pos, x, t.expand(B, 1, 1))
            dx = (f(t)*x + g2(t)*pred_eps/sigma(t))*dt + g(t)*eps
            x += dx
            x_traj.append(x.clone())

        x_traj = torch.stack(x_traj, dim=1)  # [B, steps, num_cg, 3]
        cg_denoise_traj = ucg_pos[:, None, self.scatter_idx, :] + x_traj
        if return_last:
            return cg_denoise_traj[:, -1, :, :]
        else:
            return cg_denoise_traj

    def generate(self, ucg_pos: Tensor, steps: int = 128, num_cg: int = 4716) -> Tensor:
        B = ucg_pos.size(0)
        cg_disp = torch.randn(B, num_cg, 3, device=ucg_pos.device)
        return self.denoise_by_sde(cg_disp, ucg_pos, steps=steps, return_last=True)
