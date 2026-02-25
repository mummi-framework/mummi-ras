"""
Adopted from https://github.com/mikaylagawarecki/transformer_tutorial_accompaniment
Explicit differences from nn.TransformerEncoderLayer:
- No need for src_key_padding_mask with nested tensors :)
- Only supports batch_first=True: nested tensors do not support seq_len as the first dimension
- Unnecessary fast path logic is removed
"""

import torch
from torch import nn, Tensor
import copy
from .attention import MultiHeadAttention, MultiHeadAttentionRoPE, MultiHeadAttentionContinuousRoPE


class TransformerEncoderLayerRoPE(nn.Module):
    def __init__(self,
        dim: int, num_heads: int, ff_dim: int, dropout: float, activation: nn.Module,
        layer_norm_eps: float = 1e-5, norm_first: bool = True, bias: bool = True, device=None, dtype=None,
    ):
        factory_kwargs = {"device": device, "dtype": dtype}
        super().__init__()
        self.self_attn = MultiHeadAttentionRoPE(dim, dim, dim, dim, num_heads, dropout, bias=bias, **factory_kwargs)
        
        self.linear1 = nn.Linear(dim, ff_dim, bias=bias, **factory_kwargs)
        self.dropout = nn.Dropout(dropout)
        self.linear2 = nn.Linear(ff_dim, dim, bias=bias, **factory_kwargs)

        self.norm_first = norm_first
        self.norm1 = nn.LayerNorm(dim, eps=layer_norm_eps, bias=bias, **factory_kwargs)
        self.norm2 = nn.LayerNorm(dim, eps=layer_norm_eps, bias=bias, **factory_kwargs)

        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)
        self.activation = activation

    def _sa_block(self, x, attn_mask, is_causal):
        x = self.self_attn(x, x, x, attn_mask=attn_mask, is_causal=is_causal)
        return self.dropout1(x)

    def _ff_block(self, x):
        x = self.linear2(self.dropout(self.activation(self.linear1(x))))
        return self.dropout2(x)

    def forward(self, src: Tensor, src_mask: Tensor | None = None, is_causal: bool = False) -> Tensor:
        """
        Args:
            src (Tensor): input token embeddings (batch_size, seq_len, d_model)
            src_mask (Tensor, optional): attention mask (batch_size, seq_len, seq_len). Default is None
            is_causal (bool, optional): whether to apply causal mask. Default is False 
        """
        x = src
        if self.norm_first:
            x = x + self._sa_block(self.norm1(x), src_mask, is_causal)
            x = x + self._ff_block(self.norm2(x))
        else:
            x = self.norm1(x + self._sa_block(x, src_mask, is_causal))
            x = self.norm2(x + self._ff_block(x))
        return x


class TransformerEncoderLayerContinuousRoPE(nn.Module):
    def __init__(self,
        dim: int, num_heads: int, ff_dim: int, dropout: float, activation: nn.Module, rope_dim : int,
        layer_norm_eps: float = 1e-5, norm_first: bool = True, bias: bool = True, device=None, dtype=None,
    ):
        factory_kwargs = {"device": device, "dtype": dtype}
        super().__init__()
        self.self_attn = MultiHeadAttentionContinuousRoPE(dim, dim, dim, dim, num_heads, dropout, rope_dim=rope_dim, bias=bias, **factory_kwargs)
        
        self.linear1 = nn.Linear(dim, ff_dim, bias=bias, **factory_kwargs)
        self.dropout = nn.Dropout(dropout)
        self.linear2 = nn.Linear(ff_dim, dim, bias=bias, **factory_kwargs)

        self.norm_first = norm_first
        self.norm1 = nn.LayerNorm(dim, eps=layer_norm_eps, bias=bias, **factory_kwargs)
        self.norm2 = nn.LayerNorm(dim, eps=layer_norm_eps, bias=bias, **factory_kwargs)

        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)
        self.activation = activation

    def _sa_block(self, x, pos, attn_mask, is_causal):
        x = self.self_attn(x, x, x, pos, attn_mask=attn_mask, is_causal=is_causal)
        return self.dropout1(x)

    def _ff_block(self, x):
        x = self.linear2(self.dropout(self.activation(self.linear1(x))))
        return self.dropout2(x)

    def forward(self, src: Tensor, pos: Tensor, src_mask: Tensor | None = None, is_causal: bool = False) -> Tensor:
        """
        Args:
            src (Tensor): input token embeddings (batch_size, seq_len, d_model)
            pos (Tensor): continuous token positions (batch_size, seq_len, 3)
            src_mask (Tensor, optional): attention mask (batch_size, seq_len, seq_len). Default is None
            is_causal (bool, optional): whether to apply causal mask. Default is False 
        """
        x = src
        if self.norm_first:
            x = x + self._sa_block(self.norm1(x), pos, src_mask, is_causal)
            x = x + self._ff_block(self.norm2(x))
        else:
            x = self.norm1(x + self._sa_block(x, pos, src_mask, is_causal))
            x = self.norm2(x + self._ff_block(x))
        return x


class TransformerEncoderRoPE(nn.Module):
    def __init__(self, encoder_layer: TransformerEncoderLayerRoPE, num_layers: int, norm: nn.Module | None = None, device=None, dtype=None):
        factory_kwargs = {"device": device, "dtype": dtype}
        super().__init__()
        self.layers = nn.ModuleList([copy.deepcopy(encoder_layer) for _ in range(num_layers)])
        self.num_layers = num_layers
        self.norm = norm

    def forward(self, src: Tensor, mask: Tensor | None = None, is_causal: bool = False) -> Tensor:
        output = src
        for layer in self.layers:
            output = layer(output, mask, is_causal)
        if self.norm is not None:
            output = self.norm(output)
        return output


class TransformerEncoderContinuousRoPE(nn.Module):
    def __init__(self, encoder_layer: TransformerEncoderLayerContinuousRoPE, num_layers: int, norm: nn.Module | None = None, device=None, dtype=None):
        factory_kwargs = {"device": device, "dtype": dtype}
        super().__init__()
        self.layers = nn.ModuleList([copy.deepcopy(encoder_layer) for _ in range(num_layers)])
        self.num_layers = num_layers
        self.norm = norm

    def forward(self, src: Tensor, pos: Tensor, mask: Tensor | None = None, is_causal: bool = False) -> Tensor:
        output = src
        for layer in self.layers:
            output = layer(output, pos, mask, is_causal)
        if self.norm is not None:
            output = self.norm(output)
        return output