import torch
from torch import nn, Tensor


class MLP(nn.Module):
    """Multi-layer perceptron.
    """
    def __init__(self, dims: list[int], act: nn.Module | None) -> None:
        """
        Args:
            dims (list of int): Input, hidden, and output dimensions.
            act (activation function, or None): Activation function that
                applies to all but the output layer. For example, 'nn.ReLU()'.
                If None, no activation function is applied.
        """
        super().__init__()
        self.dims = dims
        self.act = act
        
        num_layers = len(dims)

        layers = []
        for i in range(num_layers-1):
            layers += [nn.Linear(dims[i], dims[i+1])]
            if (act is not None) and (i < num_layers-2):
                layers += [act]

        self.mlp = nn.Sequential(*layers)

    def forward(self, x: Tensor) -> Tensor:
        return self.mlp(x)
    
    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(dims={self.dims}, act={self.act})'
