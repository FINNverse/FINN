import torch
import numpy as np
from torch.distributions import Exponential

def sample_poisson_relaxed(lmbd, num_samples=100, temperature = 1e-2):
    z = Exponential(lmbd).rsample([num_samples])
    t = torch.cumsum(z,0)
    relaxed_indicator = torch.sigmoid((1.0 - t) / temperature)
    N = relaxed_indicator.sum(0)
    return N

@torch.jit.script
def pad_tensors_speed_up(value: torch.Tensor, indices: torch.Tensor, org_dim: torch.Tensor) -> torch.Tensor:
    """Pad tensors with speed up.
    
    Args:
        value (torch.Tensor): The input tensor to be padded.
        indices (torch.Tensor): The indices tensor used for padding.
        org_dim (torch.Tensor): The original dimensions of the tensor.
    
    Returns:
        torch.Tensor: The padded tensor.

    dbh, nTree : [sites, patches, cohorts, 1]
    species: [sites, patches, cohorts]
    parGlobal[Species]
    """
    
    KK = torch.tensor_split(value.flatten(0, 1), org_dim[0]*org_dim[1])
    KK_new = []
    for i in range(indices.shape[0]):
        KK_new.append(KK[i][0,indices[i,]])
    return torch.nn.utils.rnn.pad_sequence(KK_new, batch_first=True)