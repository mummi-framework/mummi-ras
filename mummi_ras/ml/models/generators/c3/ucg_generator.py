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
from .modules.basis import GaussianRandomFourierFeatures
from .modules.diffusion import VariancePreservingDiffuser
from .modules.transformer import TransformerEncoderContinuousRoPE, TransformerEncoderLayerContinuousRoPE
from scipy.spatial.transform import Rotation

try:
    import lightning as L
    from lightning.pytorch.utilities import grad_norm
except ImportError:
    import pytorch_lightning as L
    from pytorch_lightning.utilities import grad_norm


class TransformerBackbone(nn.Module):
    def __init__(self, v_dim: int, dim: int, ff_dim: int, num_heads: int, num_layers: int, dropout: float):
        super().__init__()
        self.v_dim      = v_dim
        self.dim        = dim
        self.ff_dim     = ff_dim
        self.num_heads  = num_heads
        self.num_layers = num_layers
        self.dropout    = dropout

        self.embed_v = nn.Linear(v_dim, dim//2)
        self.embed_time = GaussianRandomFourierFeatures(dim//2, input_dim=1)
        self.embed_token = nn.Embedding(250, dim)
        self.transformer_encoder = TransformerEncoderContinuousRoPE(
            encoder_layer = TransformerEncoderLayerContinuousRoPE(dim, num_heads, ff_dim, dropout, activation=nn.functional.silu, rope_dim=3),
            num_layers = num_layers,
            norm = nn.LayerNorm(dim),
        )
        self.lin_out = nn.Linear(dim, 3)
    
    def forward(self, v: Tensor, ucg: Tensor, t: Tensor) -> Tensor:
        B = v.size(0)
        device = v.device

        v = v.log()

        # Embed UCG bead tokens: [B, num_ucg] -> [B, num_ucg, D]
        ucg_tokens = self.embed_token(torch.arange(250, device=device).expand(B, -1))

        # Embed conditional token containing v and t
        v_emb = self.embed_v(v)[:, None, :]  # [B, 4] -> [B, D/2] -> [B, 1, D/2]
        t_emb = self.embed_time(t)           # [B, 1, 1] -> [B, 1, D/2]
        cond_token = torch.cat([v_emb, t_emb], dim=-1)  # [B, 1, D/2] : [B, 1, D/2] -> [B, 1, D]
        
        # Append conditional v-t token to UCG tokens -> [B, num_ucg+1, D]
        all_tokens = torch.cat([ucg_tokens, cond_token], dim=1)
        
        # Append all-zero, dummy positions for conditional token to UCG positions -> [B, num_ucg+1, 3]
        pos = torch.cat([ucg, torch.zeros(B, 1, 3, device=device)], dim=1)
        
        # Feed all tokens and positions into transformer encoder -> [B, num_ucg+1, D]
        out = self.transformer_encoder(all_tokens, pos)

        # Decode transformer outputs: [B, num_ucg+1, D] -> [B, num_ucg, D] -> [B, num_ucg, 3]
        return self.lin_out(out[:, :-1, :])


class V2UCG_Generator(L.LightningModule):
    def __init__(self, v_dim: int, dim: int, ff_dim: int, num_heads: int, num_layers: int, dropout: float, learn_rate: float, T_max: int):
        super().__init__()
        self.save_hyperparameters()

        # Core model
        self.model = TransformerBackbone(
            v_dim      = v_dim,
            dim        = dim,
            ff_dim     = ff_dim,
            num_heads  = num_heads,
            num_layers = num_layers,
            dropout    = dropout,
        )

        # Data diffuser
        self.diffuser = VariancePreservingDiffuser()

        # Keep track of these parameters
        self.learn_rate = learn_rate
        self.T_max = T_max

    def forward(self, v: Tensor, ucg: Tensor, t: Tensor) -> Tensor:
        return self.model(v, ucg, t)

    def _random_rotate(self, batch):
        B = batch['ucg'].size(0)
        theta = np.random.uniform(low=0, high=2*np.pi, size=(B, 1))
        R = Rotation.from_rotvec(np.array([[0., 0., 1.]]) * theta).as_matrix()  # random rotation around z-axis
        R = torch.tensor(R, dtype=torch.float, device=batch['ucg'].device)
        batch['ucg'] @= R  # [B, num_ucg, 3] @ [B, 3, 3] -> [B, num_ucg, 3] 
        return batch

    def _get_loss(self, batch):
        v, ucg = batch['v'], batch['ucg']
        t = torch.empty(v.size(0), 1, 1, device=v.device).uniform_(self.diffuser.t_min, self.diffuser.t_max)
        eps = torch.randn_like(ucg)
        eps -= eps.mean(dim=1, keepdim=True)
        noisy_ucg = self.diffuser.alpha(t) * ucg + self.diffuser.sigma(t) * eps
        pred_eps = self.model(v, noisy_ucg, t)
        return nn.functional.mse_loss(pred_eps, eps)

    def training_step(self, batch, batch_idx):
        batch = self._random_rotate(batch)
        loss = self._get_loss(batch)
        self.log('train_loss', loss, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        loss = self._get_loss(batch)
        self.log('valid_loss', loss, prog_bar=True, sync_dist=True)
        self.log('hp_metric', loss, sync_dist=True)

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

    def denoise_by_sde(self, ucg: Tensor, v: Tensor, t_span: tuple[float, float] = (0.999, 0.001), steps: int = 128, return_last: bool = False) -> Tensor:
        B = v.size(0)
        device = v.device
        
        t_start, t_end = t_span
        ts = torch.linspace(t_start, t_end, steps+1, device=device)
        f, g, g2 = self.diffuser.f, self.diffuser.g, self.diffuser.g2
        sigma = self.diffuser.sigma

        x = ucg        
        x -= x.mean(dim=1, keepdim=True)
        x_traj = [x.clone()]
        for i, t in enumerate(ts[1:]):
            dt = ts[i+1] - ts[i]
            eps = dt.abs().sqrt() * torch.randn_like(x)
            pred_eps = self.model(v, x, t.expand(B, 1, 1))
            dx = (f(t)*x + g2(t)*pred_eps/sigma(t))*dt + g(t)*eps
            x += dx
            x -= x.mean(dim=1, keepdim=True)
            x_traj.append(x.clone())

        ucg_denoise_traj = torch.stack(x_traj, dim=1)  # [B, steps, num_ucg, 3]
        if return_last:
            return ucg_denoise_traj[:, -1, :, :]
        else:
            return ucg_denoise_traj
    
    def generate(self, v: Tensor, steps: int = 128, num_ucg: int = 250) -> Tensor:
        B = v.size(0)
        noisy_ucg = torch.randn(B, num_ucg, 3, device=v.device)
        return self.denoise_by_sde(noisy_ucg, v, steps=steps, return_last=True)
