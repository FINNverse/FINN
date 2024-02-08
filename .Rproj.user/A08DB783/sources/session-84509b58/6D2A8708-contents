import torch
import numpy as np
from torch.distributions import Exponential
import matplotlib.pyplot as plt
Bernoulli = torch.distributions.relaxed_bernoulli.RelaxedBernoulli
import pandas as pd
from typing import Union, Tuple, List, Optional, Callable

def sample_poisson_relaxed(lmbd, num_samples=100, temperature = 1e-2):
    z = Exponential(lmbd).rsample([num_samples])
    t = torch.cumsum(z,0)
    relaxed_indicator = torch.sigmoid((1.0 - t) / temperature)
    N = relaxed_indicator.sum(0)
    return N

# Dimensions: sites, populations, cohorts, ...

# Env: sites, time, variables

class CohortMat:
    def __init__(self, dbh=None, nTree=None, Species=None):
        self.nTree = np.random.poisson(10, [13, 33, 10, 1]) if nTree is None else nTree
        self.Species = np.random.randint(5, size=[13, 33, 10]) if Species is None else Species
        self.dbh = np.random.uniform(10, 90, size=[13, 33, 10,1]) if dbh is None else dbh
        
        if not torch.is_tensor(self.nTree):
            self.nTree = torch.tensor(self.nTree, dtype=torch.float32)
        if not torch.is_tensor(self.Species):
            self.Species = torch.tensor(self.Species, dtype=torch.int32)
        if not torch.is_tensor(self.dbh):
            self.dbh = torch.tensor(self.dbh, dtype=torch.float32)
            
@torch.jit.script
def BA_T_P(dbh: torch.Tensor, nTree: torch.Tensor) -> torch.Tensor:
    return torch.pi*(dbh*nTree/100./2.).pow(2)


@torch.jit.script
def BA_P(dbh: torch.Tensor) -> torch.Tensor:
    return torch.pi*(dbh/100./2.).pow(2.0)


@torch.jit.script
def height_P(dbh: torch.Tensor, par: torch.Tensor) -> torch.Tensor:
    height = (dbh*par*0.03).exp()
    return height

def compF_P(dbh: torch.Tensor, Species: torch.Tensor, parGlobal: torch.Tensor, h: Optional[torch.Tensor]=None, minLight: float=50.):
    ba = BA_P(dbh)/0.1
    cohortHeights = height_P(dbh, parGlobal[Species][...,None])
    if h is None:
        h = cohortHeights
        BA_height = (ba*torch.sigmoid((cohortHeights - h.permute(0,1, 3, 2) - 0.1)/1e-3)).sum(-1)
    else:
        BA_height = (ba*torch.sigmoid((cohortHeights - 0.1)/1e-3)).sum(-1) # A[1:3] = 0
    AL = 1.-BA_height/minLight
    AL = torch.clamp(AL, min = 0)
    return AL


def growthFP(dbh: torch.Tensor, Species: torch.Tensor, parGrowth: torch.Tensor, parGlobal: torch.Tensor, env: torch.Tensor, nnGrowthEnv: torch.nn.modules.container.Sequential) -> torch.Tensor:
    AL = compF_P(dbh, Species, parGlobal)
    shade = ((AL**2)*parGlobal[Species]).sigmoid()
    environment = nnGrowthEnv(env)
    environment = torch.gather(environment[:,None,...].repeat(1, Species.shape[1], 1,), 2, Species.type(torch.int64))
    pred = (shade+environment)
    growth = (1.- torch.pow(1.- pred,4.0)) * parShade[Species]
    return torch.nn.functional.softplus(growth)


def mortFP(dbh: torch.Tensor, Species: torch.Tensor, nTree: torch.Tensor, parMort: torch.Tensor,parGlobal: torch.Tensor, envM: torch.Tensor, nnMortEnv: torch.nn.modules.container.Sequential) -> torch.Tensor:
    AL = compF_P(dbh, Species, parGlobal=parGlobal)
    shade = 1-((AL**2)*parGlobal[Species]).sigmoid()
    environment = 1 - nnMortEnv(envM)
    gPSize = 0.1*(torch.clamp(dbh.squeeze(3)/parMort[Species], min = 0.00001) ).pow(2.3) #.reshape([-1,1])
    predM = torch.clamp((shade+environment[:,None,...].repeat(1, Species.shape[1], 1,)+gPSize), 0, 1)
    mort = torch.distributions.Beta(predM*nTree.squeeze(3)+0.00001, nTree.squeeze(3) - predM*nTree.squeeze(3)+0.00001).rsample()*nTree.squeeze(3)
    return mort + mort.round().detach() - mort.detach() 

def regFP(dbh: torch.Tensor, Species: torch.Tensor, parReg: torch.Tensor ,parGlobal: torch.Tensor, envM: torch.Tensor, nnRegEnv: torch.nn.modules.container.Sequential) -> torch.Tensor:
    # TODO: Species specific? Independent of current cohorts?!
    AL = compF_P(dbh, Species, parGlobal=parGlobal, h = torch.zeros([1, 1]))
    regP = torch.sigmoid((parReg - AL)/1e-2)
    environment = nnRegEnv(envM).squeeze()   
    regeneration = sample_poisson_relaxed(regP+environment[:,None,...].repeat(1, Species.shape[1], 1,)*0.9, 20)
    regeneration = regeneration+ regeneration.round().detach() - regeneration.detach() 
    return regeneration


class Parameters:
    def __init__(self, sp, env_n):
      
        self.nnRegEnv = torch.nn.Sequential(
            torch.nn.Linear(env_n, sp),
            torch.nn.Sigmoid()
        )

        self.nnGrowthEnv = torch.nn.Sequential(
            torch.nn.Linear(env_n, sp),
            torch.nn.Sigmoid()
        )

        self.nnMortEnv = torch.nn.Sequential(
            torch.nn.Linear(env_n, sp),
            torch.nn.Sigmoid()
        )

        self.parGlobal = torch.rand([sp], requires_grad=True, dtype=torch.float32)
        self.parGrowth = torch.rand([sp], requires_grad=True, dtype=torch.float32)
        self.parShade = torch.rand([sp], requires_grad=True, dtype=torch.float32)
        self.parMort = torch.rand([sp], requires_grad=True, dtype=torch.float32)
        self.parReg = torch.rand([sp], requires_grad=True, dtype=torch.float32)
        
        
