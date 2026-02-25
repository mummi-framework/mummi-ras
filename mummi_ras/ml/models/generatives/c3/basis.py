import torch
from torch import nn, Tensor


class BesselEmbeddings(nn.Module):
    def __init__(self, start: float = 0.0, end: float = 1.0, num_basis: int = 8, eps: float = 1e-5) -> None:
        super().__init__()
        self.start     = start
        self.end       = end
        self.num_basis = num_basis
        self.eps       = eps
        self.register_buffer('n', torch.arange(1, num_basis+1, dtype=torch.float), persistent=False)

    def forward(self, x: Tensor) -> Tensor:
        x = x[..., None] - self.start + self.eps
        c = self.end - self.start
        return ((2/c)**0.5) * torch.sin(self.n*torch.pi*x / c) / x

    def extra_repr(self) -> str:
        return f'start={self.start}, end={self.end}, num_basis={self.num_basis}, eps={self.eps}'


class GaussianEmbeddings(nn.Module):
    def __init__(self, start: float = 0.0, end: float = 1.0, num_basis: int = 8) -> None:
        super().__init__()
        self.start     = start
        self.end       = end
        self.num_basis = num_basis
        self.register_buffer('mu', torch.linspace(start, end, num_basis), persistent=False)
    
    def forward(self, x: Tensor) -> Tensor:
        step = self.mu[1] - self.mu[0]
        diff = (x[..., None] - self.mu) / step
        return diff.pow(2).neg().exp().div(1.12) # division by 1.12 so that sum of square is roughly 1

    def extra_repr(self) -> str:
        return f'start={self.start}, end={self.end}, num_basis={self.num_basis}'


class GaussianRandomFourierFeatures(nn.Module):
    """Gaussian random Fourier features.

    Reference: https://arxiv.org/abs/2006.10739
    """
    def __init__(self, embed_dim: int, input_dim: int = 1, sigma: float = 1.0) -> None:
        super().__init__()

        self.embed_dim = embed_dim
        self.input_dim = input_dim
        self.sigma     = sigma

        # Randomly sample weights during initialization. These weights are fixed
        # during optimization and are not trainable.
        self.register_buffer('B', torch.randn(input_dim, embed_dim//2) * sigma)

    def forward(self, v: Tensor) -> Tensor:
        v_proj =  2 * torch.pi * v @ self.B
        return torch.cat([torch.cos(v_proj), torch.sin(v_proj)], dim=-1)

    def extra_repr(self) -> str:
            return f'embed_dim={self.embed_dim}, input_dim={self.input_dim}, sigma={self.sigma}'
