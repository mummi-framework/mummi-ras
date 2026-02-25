import torch
from torch import nn, Tensor


# Rotary PE implementation from https://pytorch.org/torchtune/stable/_modules/torchtune/modules/position_embeddings.html#RotaryPositionalEmbeddings
class RotaryPositionalEmbeddings(nn.Module):
    """
    This class implements Rotary Positional Embeddings (RoPE)
    proposed in https://arxiv.org/abs/2104.09864.

    Reference implementation (used for correctness verfication)
    can be found here:
    https://github.com/meta-llama/llama/blob/main/llama/model.py#L80

    In this implementation we cache the embeddings for each position upto
    ``max_seq_len`` by computing this during init.

    Args:
        dim (int): Embedding dimension. This is usually set to the dim of each
            head in the attention module computed as ``embed_dim // num_heads``
        max_seq_len (int): Maximum expected sequence length for the
            model, if exceeded the cached freqs will be recomputed
        base (int): The base for the geometric progression used to compute
            the rotation angles
    """

    def __init__(self, dim: int, max_seq_len: int = 4096, base: int = 10_000) -> None:
        super().__init__()
        self.dim = dim
        self.base = base
        self.max_seq_len = max_seq_len
        self.rope_init()

    def rope_init(self):
        theta = 1.0 / (
            self.base
            ** (torch.arange(0, self.dim, 2)[: (self.dim // 2)].float() / self.dim)
        )
        self.register_buffer("theta", theta, persistent=False)
        self.build_rope_cache(self.max_seq_len)

    def build_rope_cache(self, max_seq_len: int = 4096) -> None:
        # Create position indexes `[0, 1, ..., max_seq_len - 1]`
        seq_idx = torch.arange(
            max_seq_len, dtype=self.theta.dtype, device=self.theta.device
        )

        # Outer product of theta and position index; output tensor has
        # a shape of [max_seq_len, dim // 2]
        idx_theta = torch.einsum("i, j -> ij", seq_idx, self.theta).float()

        # cache includes both the cos and sin components and so the output shape is
        # [max_seq_len, dim // 2, 2]
        cache = torch.stack([torch.cos(idx_theta), torch.sin(idx_theta)], dim=-1)
        self.register_buffer("cache", cache, persistent=False)

    def forward(self, x: Tensor, *, input_pos: Tensor | None = None) -> Tensor:
        """
        Args:
            x (torch.Tensor): input tensor with shape
                ``[b, s, n_h, h_d]``
            input_pos (Optional[torch.Tensor]): Optional tensor which contains the position ids
                of each token. During training, this is used to indicate the positions
                of each token relative to its sample when packed, shape [b, s].
                During inference, this indicates the position of the current token.
                If none, assume the index of the token is its position id. Default is None.

        Returns:
            torch.Tensor: output tensor with shape ``[b, s, n_h, h_d]``

        Notation used for tensor shapes:
            - b: batch size
            - s: sequence length
            - n_h: num heads
            - h_d: head dim
        """
        # input tensor has shape [b, s, n_h, h_d]
        seq_len = x.size(1)

        # extract the values based on whether input_pos is set or not
        rope_cache = (
            self.cache[:seq_len] if input_pos is None else self.cache[input_pos]
        )

        # reshape input; the last dimension is used for computing the output.
        # Cast to float to match the reference implementation
        # tensor has shape [b, s, n_h, h_d // 2, 2]
        xshaped = x.float().reshape(*x.shape[:-1], -1, 2)

        # reshape the cache for broadcasting
        # tensor has shape [b, s, 1, h_d // 2, 2] if packed samples,
        # otherwise has shape [1, s, 1, h_d // 2, 2]
        rope_cache = rope_cache.view(-1, xshaped.size(1), 1, xshaped.size(3), 2)

        # tensor has shape [b, s, n_h, h_d // 2, 2]
        x_out = torch.stack(
            [
                xshaped[..., 0] * rope_cache[..., 0] - xshaped[..., 1] * rope_cache[..., 1],
                xshaped[..., 1] * rope_cache[..., 0] + xshaped[..., 0] * rope_cache[..., 1],
            ],
            dim = -1,
        )

        # tensor has shape [b, s, n_h, h_d]
        x_out = x_out.flatten(3)
        return x_out.type_as(x)


class ContinuousRotaryPositionalEmbeddings(nn.Module):
    def __init__(self, dim: int, base: int = 10_000):
        """
        Rotary Positional Embeddings for continuous positions
        
        Reference implementation:
        https://pytorch.org/torchtune/stable/_modules/torchtune/modules/position_embeddings.html
        
        Args:
            dim (int): Embedding dimension. This is usually set to the dim of each
                head in the attention module computed as ``embed_dim // num_heads``
            base (int): The base for the geometric progression used to compute
                the rotation angles
        """
        super().__init__()
        self.dim = dim
        self.base = base

        i = torch.arange(1, dim//2 + 1)
        theta = 1.0 / (base ** (2*(i-1) / dim)) # shape [h_d/2]
        self.register_buffer('theta', theta, persistent=False)

    def forward(self, x: Tensor, input_pos: Tensor) -> Tensor:
        """
        Args:
            x (Tensor): input tensor with shape ``[b, s, n_h, h_d]``
            input_pos (Tensor): continuous positions of shape ``[b, s]``

        Returns:
            Tensor: output tensor with shape ``[b, s, n_h, h_d]``
        """
        b, s, n_h, h_d = x.shape

        # Multiply input_pos (aka 'm') with theta: [b, s, 1] * [1, 1, h_d/2] -> [b, s, h_d/2]
        m_theta = input_pos[..., None] * self.theta[None, None, ...]

        # Compute cos(input_pos*theta) and sin(input_pos*theta) -> [b, s, h_d/2]
        cos_m_theta = torch.cos(m_theta)
        sin_m_theta = torch.sin(m_theta)
        
        # Reshape input tensor: [b, s, n_h, h_d] -> [b, s, n_h, h_d/2, 2]
        xshaped = x.reshape(*x.shape[:-1], -1, 2)

        # Reshape cos and sin terms for broadcasting: [b, s, h_d/2] -> [b, s, 1, h_d/2]
        cos_m_theta = cos_m_theta.view(b, s, 1, -1)
        sin_m_theta = sin_m_theta.view(b, s, 1, -1)

        # [b, s, n_h, h_d/2], [b, s, n_h, h_d/2] -> [b, s, n_h, h_d/2, 2]
        x_out = torch.stack(
            [
                xshaped[..., 0] * cos_m_theta - xshaped[..., 1] * sin_m_theta,
                xshaped[..., 1] * cos_m_theta + xshaped[..., 0] * sin_m_theta,
            ],
            dim = -1,
        )

        # [b, s, n_h, h_d/2, 2] -> [b, s, n_h, h_d]
        x_out = x_out.flatten(3)
        return x_out
    
    def extra_repr(self) -> str:
        return f'dim={self.dim}, base={self.base}'


class ContinuousRotaryPositionalEmbeddings3D(nn.Module):
    """
    Rotary Positional Embeddings for continuous 3D positions

    For convenience, the token embedding dim is split evenly into 4 chunks, where:
    - chunk 0 is rotated based on x-axis positions
    - chunk 1 is rotated based on y-axis positions
    - chunk 2 is rotated based on z-axis positions
    - chunk 3 is NOT rotated
    Otherwise, the embedding dim (per-head) has to be divisible by 3, which is annoying
    """
    def __init__(self, dim: int, base: int = 10_000):
        super().__init__()
        assert dim % 4 == 0, 'dim must be divisible by 4'
        self.dim = dim
        self.base = base
        self.rope_1d = ContinuousRotaryPositionalEmbeddings(dim//4, base)
    
    def forward(self, x: Tensor, pos: Tensor) -> Tensor:
        """
        Args:
            x (Tensor): input tensor with shape ``[b, s, n_h, h_d]``
            input_pos (Tensor): continuous positions of shape ``[b, s, 3]``

        Returns:
            Tensor: rotated `x` with shape ``[b, s, n_h, h_d]``
        """
        x0, x1, x2, x3 = x.chunk(4, dim=-1)
        out = torch.cat(
            [
                self.rope_1d(x0, pos[..., 0]),
                self.rope_1d(x1, pos[..., 1]),
                self.rope_1d(x2, pos[..., 2]),
                x3,
            ],
            dim = -1,
        )
        return out


class ContinuousRotaryPositionalEmbeddings4D(nn.Module):
    """
    Rotary Positional Embeddings for continuous 4D positions

    The token embedding dim is split evenly into 4 chunks, where the i-th chunk
    is rotated based on the i-th dimension of the positions. 
    """
    def __init__(self, dim: int, base: int = 10_000):
        super().__init__()
        assert dim % 4 == 0, 'dim must be divisible by 4'
        self.dim = dim
        self.base = base
        self.rope_1d = ContinuousRotaryPositionalEmbeddings(dim//4, base)
    
    def forward(self, x: Tensor, pos: Tensor) -> Tensor:
        """
        Args:
            x (Tensor): input tensor with shape ``[b, s, n_h, h_d]``
            input_pos (Tensor): continuous positions of shape ``[b, s, 4]``

        Returns:
            Tensor: rotated `x` with shape ``[b, s, n_h, h_d]``
        """
        x0, x1, x2, x3 = x.chunk(4, dim=-1)
        out = torch.cat(
            [
                self.rope_1d(x0, pos[..., 0]),
                self.rope_1d(x1, pos[..., 1]),
                self.rope_1d(x2, pos[..., 2]),
                self.rope_1d(x3, pos[..., 3]),
            ],
            dim = -1,
        )
        return out