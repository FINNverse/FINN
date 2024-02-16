import sys
import os
from typing import Union, Tuple, List, Optional, Callable
import itertools
from tqdm import tqdm
import torch
import numpy as np
from torch.distributions import Exponential
import torch_optimizer as optim

os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"

Bernoulli = torch.distributions.relaxed_bernoulli.RelaxedBernoulli

def sample_poisson_relaxed(lmbd, num_samples=100, temperature = 1e-2):
    z = Exponential(lmbd).rsample([num_samples])
    t = torch.cumsum(z,0)
    relaxed_indicator = torch.sigmoid((1.0 - t) / temperature)
    N = relaxed_indicator.sum(0)
    return N


class CohortMat:
    def __init__(self, dbh=None, nTree=None, Species=None, dims=[50, 30, 10,1], sp = 10, device: str='cpu'):
        """Initializes an instance of the class cohortMat.
        
        Args:
            dbh (ndarray, optional): Array of tree diameters at breast height. If not provided, it will be randomly generated using a uniform distribution between 10 and 90. Default is None.
            nTree (ndarray, optional): Array of tree counts. If not provided, it will be randomly generated using a Poisson distribution with a mean of 10. Default is None.
            Species (ndarray, optional): Array of tree species. If not provided, it will be randomly generated using integers between 0 and sp, where sp is the number of species. Default is None.
            dims (list, optional): List of dimensions for the arrays. Default is [50, 30, 10, 1].
            sp (int, optional): Number of species. Default is 10.
            device (str, optional): Device to be used for computation. Default is 'cpu'.
        
        Returns:
            None
        """
        
        device = torch.device(device)
        # TODO initial kohorten
        self.nTree = np.random.poisson(10, size = dims) if nTree is None else nTree 
        self.Species = np.random.randint(0, sp, size=dims[0:3]) if Species is None else Species
        self.dbh = np.random.uniform(10, 90, size=dims) if dbh is None else dbh
        
        if not torch.is_tensor(self.nTree):
            self.nTree = torch.tensor(self.nTree, dtype=torch.float32, device=device)
        if not torch.is_tensor(self.Species):
            self.Species = torch.tensor(self.Species, dtype=torch.int64, device=device)
        if not torch.is_tensor(self.dbh):
            self.dbh = torch.tensor(self.dbh, dtype=torch.float32, device=device)

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


@torch.jit.script
def BA_T_P(dbh: torch.Tensor, nTree: torch.Tensor) -> torch.Tensor:
    """Calculate the basal area of trees.
    
    Args:
        dbh (torch.Tensor): Diameter at breast height of the trees.
        nTree (torch.Tensor): Number of trees.
    
    Returns:
        torch.Tensor: Basal area of the trees.
    
    Example:
        >>> dbh = torch.tensor([10, 15, 20])
        >>> nTree = torch.tensor([5, 10, 15])
        >>> BA_T_P(dbh, nTree)
        tensor([ 0.7854,  1.1781,  1.5708])
    """
    
    return torch.pi*(dbh/100./2.).pow(2)*nTree    
    
@torch.jit.script
def BA_P(dbh: torch.Tensor) -> torch.Tensor:
    """Calculate the basal area of a tree given the diameter at breast height (dbh).
    
    Args:
        dbh (torch.Tensor): The diameter at breast height of the tree.
    
    Returns:
        torch.Tensor: The basal area of the tree.
    
    Example:
        dbh = torch.tensor(50)
        basal_area = BA_P(dbh)
        print(basal_area)  # Output: 1963.4954084936207
    """
    
    return torch.pi*(dbh/100./2.).pow(2.0)

@torch.jit.script
def height_P(dbh: torch.Tensor, parGlobal: torch.Tensor) -> torch.Tensor:
    """Calculate the height of a tree based on its diameter at breast height (dbh) and a parameter (par).
    
    Args:
        dbh (torch.Tensor): The diameter at breast height of the tree.
        par (torch.Tensor): The parameter used to calculate the height.
    
    Returns:
        torch.Tensor: The calculated height of the tree.
    """
    
    height = (dbh*parGlobal*0.03).exp()
    return height

class FINN:
    def __init__(self, 
                 device: str='cpu',  # 'mps'
                 sp: int=5, 
                 env: int=2, 
                 parGlobal: Optional[np.ndarray]=None, # must be dim [species]
                 parGrowth: Optional[np.ndarray]=None, # must be dim [species, 2], first for shade tolerance
                 parMort: Optional[np.ndarray]=None, # must be dim [species, 2], first for shade tolerance,
                 parReg: Optional[np.ndarray]=None, # must be dim [species]
                 hidden_growth: List[int] = [],
                 hidden_mort: List[int]  = [],
                 hidden_reg: List[int]  = []
                 ):
        """Initialize the model.
        
        Args:
            device (str, optional): The device to use for computation. Supported options are 'cpu', 'cuda', and 'mps'. Defaults to 'cpu'.
            sp (int, optional): The number of species. Defaults to 5.
            env (int, optional): The number of environmental covariates. Defaults to 2.
            parGlobal (Optional[np.ndarray], optional): The global parameters. Must be of dimension [species]. Defaults to None.
            parGrowth (Optional[np.ndarray], optional): The growth parameters. Must be of dimension [species, 2], with the first column representing shade tolerance. Defaults to None.
            parMort (Optional[np.ndarray], optional): The mortality parameters. Must be of dimension [species, 2], with the first column representing shade tolerance. Defaults to None.
            parReg (Optional[np.ndarray], optional): The regeneration parameters. Must be of dimension [species]. Defaults to None.
            hidden_growth (List[int], optional): The hidden layers for the growth neural network. Defaults to [].
            hidden_mort (List[int], optional): The hidden layers for the mortality neural network. Defaults to [].
            hidden_reg (List[int], optional): The hidden layers for the regeneration neural network. Defaults to [].
        
        Returns:
            None
        """
        self.device = torch.device(device)
        self.sp = sp
        self.env = env
        self.optimizer = None
        self.dtype = torch.float32
        self.nnRegEnv = self._build_NN(input_shape=env, output_shape=sp, hidden=hidden_reg, activation="selu", bias=[False], dropout=-99)
        self.nnRegEnv.to(self.device)
        self.nnGrowthEnv = self._build_NN(input_shape=env, output_shape=sp, hidden=hidden_growth, activation="selu", bias=[False], dropout=-99)
        self.nnGrowthEnv.to(self.device)
        self.nnMortEnv = self._build_NN(input_shape=env, output_shape=sp, hidden=hidden_mort, activation="selu", bias=[False], dropout=-99)
        self.nnMortEnv.to(self.device)
        
        # TODO: Sind diese Init Parametriesierungen sinnvoll?

        if parGlobal is None:
            self._parGlobal = torch.tensor(np.random.uniform(0.90, 1.0, size = [self.sp]), requires_grad=True, dtype=torch.float32, device=self.device)
        else:
            self._parGlobal = torch.tensor(parGlobal, requires_grad=True, dtype=torch.float32, device=self.device)
            
        if parGrowth is None:
            self._parGrowth = torch.tensor(np.random.uniform(-0.01, 0.01, size = [self.sp,2]), requires_grad=True, dtype=torch.float32, device=self.device)
        else: 
            self._parGrowth = torch.tensor(parGrowth, requires_grad=True, dtype=torch.float32, device=self.device)
        
        if parMort is None:
            self._parMort = torch.tensor(np.random.uniform(0.5, 0.6, size = [self.sp,2]), requires_grad=True, dtype=torch.float32, device=self.device)
        else:
            self._parMort = torch.tensor(parMort, requires_grad=True, dtype=torch.float32, device=self.device)
            
        if parReg is None:
            self._parReg = torch.tensor(np.random.uniform(0.4, 0.5, size = [self.sp]), requires_grad=True, dtype=torch.float32, device=self.device)
        else:
            self._parReg = torch.tensor(parReg, requires_grad=True, dtype=torch.float32, device=self.device)
            
        self.parameters = [[self._parGlobal], [self._parGrowth], [self._parMort], [self._parReg], self.nnRegEnv.parameters(), self.nnGrowthEnv.parameters(), self.nnMortEnv.parameters()]
        
    def _build_NN(self, 
                  input_shape: int, 
                  output_shape: int, 
                  hidden: List, 
                  bias: List[bool], 
                  activation: List[str], 
                  dropout: float) -> torch.nn.modules.container.Sequential:
        """Build neural network

        Args:
            input_shape (int): Number of predictors
            output_shape (int): Number of species
            hidden (List): List of hidden layers
            bias (List[bool]): Biases in hidden layers
            activation (List[str]): List of activation functions
            dropout (float): Dropout rate

        Returns:
            torch.nn.modules.container.Sequential: Sequential neural network object
        """                  
        model_list = torch.nn.ModuleList()
        if len(hidden) != len(activation):
            activation = [activation[0] for _ in range(len(hidden))]

        if len(bias) == 1:
            bias = [bias[0] for _ in range(len(hidden))]
            
        bias.insert(0, False)

        if len(hidden) > 0:
            for i in range(len(hidden)):
                if i == 0:
                    model_list.append(torch.nn.Linear(input_shape, hidden[i], bias=bias[i]).type(self.dtype))
                else:
                    model_list.append(torch.nn.Linear(hidden[i-1], hidden[i], bias=bias[i]).type(self.dtype))

                if activation[i] == "relu":
                     model_list.append(torch.nn.ReLU())
                if activation[i] == "selu":
                     model_list.append(torch.nn.SELU())
                if activation[i] == "leakyrelu":
                     model_list.append(torch.nn.LeakyReLU())
                if activation[i] == "tanh": 
                    model_list.append(torch.nn.Tanh())
                if activation[i] == "sigmoid":
                    model_list.append(torch.nn.Sigmoid())
                if dropout > 0.0:
                    model_list.append(torch.nn.Dropout(p=dropout))

        if len(hidden) > 0:
            model_list.append(torch.nn.Linear(hidden[-1], output_shape, bias=bias[-1]).type(self.dtype))
        else:
            model_list.append(torch.nn.Linear(input_shape, output_shape, bias=False).type(self.dtype))
        model_list.append(torch.nn.Sigmoid())    
        return torch.nn.Sequential(*model_list)       

    # TODO: lighht -> global Parameter, fitbar
    def compF_P(self, dbh: torch.Tensor, Species: torch.Tensor, h: Optional[torch.Tensor]=None, minLight: float=50.) -> torch.Tensor:
        """Competition function
        
        Args:
            dbh (torch.Tensor): Diameter at breast height
            Species (torch.Tensor): Species tensor
            h (Optional[torch.Tensor], optional): Height. Defaults to None.
            minLight (float, optional): Minimal light. Defaults to 50.
        
        Returns:
            torch.Tensor: Allocation of light
        """
        ba = BA_P(dbh)/0.1
        cohortHeights = height_P(dbh, self._parGlobal[Species][...,None])
        if h is None:
            h = cohortHeights
            BA_height = (ba*torch.sigmoid((cohortHeights - h.permute(0,1, 3, 2) - 0.1)/1e-3)).sum(-2)
        else:
            BA_height = (ba*torch.sigmoid((cohortHeights - 0.1)/1e-3)).sum(-2)
        AL = 1.-BA_height/minLight
        AL = torch.clamp(AL, min = 0)
        return AL
    
    def growthFP(self, dbh: torch.Tensor, Species: torch.Tensor, pred: torch.Tensor) -> torch.Tensor:
        """Calculate the growth of a tree based on its diameter at breast height (dbh), species, and environmental predictors.
        
        Args:
            dbh (torch.Tensor): The diameter at breast height of the tree.
            Species (torch.Tensor): The species of the tree.
            pred (torch.Tensor): The environmental predictions nnGrowthEnv(env)[Species].
        
        Returns:
            torch.Tensor: The growth of the tree.
        
        Example:
            growth = growthFP(dbh, Species, pred)
        """
        
        AL = self.compF_P(dbh, Species, self._parGlobal)
        shade = ((AL**2)*self._parGrowth[Species,0]).sigmoid()
        environment = pred[..., Species[3]]
        pred = (shade+environment)
        growth = (1.- torch.pow(1.- pred,4.0)) * self._parGrowth[Species,1]
        return torch.nn.functional.softplus(growth)
    
    def mortFP(self, dbh: torch.Tensor, Species: torch.Tensor, nTree: torch.Tensor, pred: torch.Tensor) -> torch.Tensor:
        """Calculate mortality rate for a given set of parameters.
        
        Args:
            dbh (torch.Tensor): Diameter at breast height.
            Species (torch.Tensor): Species identifier.
            nTree (torch.Tensor): Number of trees.
            pred (torch.Tensor): Prediction.
        
        Returns:
            torch.Tensor: Mortality rate.
        """
        
        AL = self.compF_P(dbh, Species)
        shade = 1-((AL**2)*self._parMort[Species,0]).sigmoid()
        environment = 1 - pred
        gPSize = 0.1*(torch.clamp(dbh.squeeze(3)/(self._parMort[Species,1]*10), min = 0.00001) ).pow(2.3) #.reshape([-1,1])
        predM = torch.sigmoid(shade+environment[..., Species[3]]+gPSize)
        mort = torch.distributions.Beta(predM*nTree.squeeze(3)+0.00001, nTree.squeeze(3) - predM*nTree.squeeze(3)+0.00001).rsample()*nTree.squeeze(3)
        return mort + mort.round().detach() - mort.detach() 
    
    # TODO: wenn alles tot, AL = 1

    def regFP(self, dbh: torch.Tensor, Species: torch.Tensor, pred: torch.Tensor) -> torch.Tensor:
        """Calculate Regeneration rate
        
        Args:
            dbh (torch.Tensor): Diameter at breast height.
            Species (torch.Tensor): Species of the tree.
            pred (torch.Tensor): Prediction.
        
        Returns:
            torch.Tensor: Regeneration of the forest plot.
        """
        
        AL = self.compF_P(dbh, Species, h = torch.zeros([1, 1])) # Was ist wenn alles tot ist?
        regP = torch.sigmoid((self._parReg - AL)/1e-2)
        environment = pred
        regeneration = sample_poisson_relaxed(regP+environment[:,None,...].repeat(1, Species.shape[1], 1,)*0.9, 20)
        regeneration = regeneration + regeneration.round().detach() - regeneration.detach() 
        return regeneration
        
    
    def __aggregate(self, labels, samples, samples_T, Result, Result_T):
        """Aggregate results.
        
        Args:
            labels (ndarray): The labels for each sample.
            samples (ndarray): The samples to be aggregated.
            samples_T (ndarray): The samples for nTrees
            Result (ndarray): The result array to store the aggregated values.
            Result_T (ndarray): The result array to store the aggregated values for nTrees.
        
        Returns:
            list: A list containing the aggregated result arrays [Result, Result_T].
        """
        
        for k in range(labels.shape[0]):
            values, positions = self.__groupby_mean(samples[k,:,:].flatten().view([-1,1]), samples_T[k,:,:].flatten().view([-1,1]), labels[k,:,:].flatten())
            Result[k,positions] = Result[k,positions] + values[0].squeeze()
            Result_T[k,positions] = Result_T[k,positions] + values[1].squeeze()       
        return [Result, Result_T]
    
    def __groupby_mean(self, value:torch.Tensor, value_T: torch.Tensor, labels:torch.LongTensor) -> (torch.Tensor, torch.LongTensor):
        """Groups the given values by their corresponding labels and calculates the mean for each group.
        
        Args:
            value (torch.Tensor): The values to be grouped.
            value_T (torch.Tensor): The transpose of the values to be grouped.
            labels (torch.LongTensor): The labels for grouping.
        
        Returns:
            tuple: A tuple containing two tensors:
                - result (torch.Tensor): The grouped values with mean calculated.
                - result_T (torch.Tensor): The transpose of the grouped values with mean calculated.
            new_labels (torch.LongTensor): The new labels corresponding to the grouped values.
        """
        
        uniques = labels.unique().tolist()
        labels = labels.tolist()
        key_val = {key: val for key, val in zip(uniques, range(len(uniques)))}
        val_key = {val: key for key, val in zip(uniques, range(len(uniques)))}
        labels = torch.tensor(list(map(key_val.get, labels)), device = value.device, dtype=torch.int64)
        labels = labels.view(labels.size(0), 1).expand(-1, value.size(1))
        unique_labels = labels.unique(dim=0, return_counts=False)
        result = torch.zeros_like(unique_labels, dtype = value.dtype).scatter_add_(0, labels, value)
        result_T = torch.zeros_like(unique_labels, dtype = value.dtype).scatter_add_(0, labels, value_T)
        new_labels = torch.LongTensor(list(map(val_key.get, unique_labels[:, 0].tolist())))
        return [result, result_T], new_labels
        
    def __pad_tensors(self, value, indices, org_dim):
        KK = torch.tensor_split(value.flatten(0, 1), org_dim[0]*org_dim[1])
        KK = [KK[i][0,indices[i,]] for i in range(indices.shape[0]) ]
        return torch.nn.utils.rnn.pad_sequence(KK, batch_first=True).unflatten(0, org_dim) 
        

    def predict(self, 
                dbh: Optional[torch.Tensor]=None, 
                nTree: Optional[torch.Tensor]=None, 
                Species: Optional[torch.Tensor]=None, 
                env: Optional[torch.Tensor]=None, 
                record_time: int=0, 
                response: str="dbh",
                pred_growth: Optional[torch.Tensor]=None, 
                pred_morth: Optional[torch.Tensor]=None, 
                pred_reg: Optional[torch.Tensor]=None,
                patches: Optional[float]=50) -> [torch.Tensor, torch.Tensor]:
        """Predicts the growth and mortality of trees based on the given inputs.

        

        
        Args:
            dbh (Optional[torch.Tensor]): The diameter at breast height of the trees. If None, it will be initialized using CohortMat.
            nTree (Optional[torch.Tensor]): The number of trees. If None, it will be initialized using CohortMat.
            Species (Optional[torch.Tensor]): The species of the trees. If None, it will be initialized using CohortMat.
            env (Optional[torch.Tensor]): The environmental data.
            record_time (int): The time at which to start recording the results.
            response (str): The response variable to use for aggregating results. Can be "dbh", "ba", or "ba_p".
            pred_growth (Optional[torch.Tensor]): The predicted growth values.
            pred_morth (Optional[torch.Tensor]): The predicted mortality values.
            pred_reg (Optional[torch.Tensor]): The predicted regeneration values.
            patches (Optional[float]): The number of patches.
        
        Returns:
            Tuple[torch.Tensor, torch.Tensor]: The predicted dbh and number of trees for the recorded times points.



        """
        
        # TODO use start
        if dbh is None:
            cohorts = CohortMat(dims = [env.shape[0], patches, self.sp, 1], sp = self.sp, device=self.device)
            nTree = cohorts.nTree
            Species = cohorts.Species
            dbh = cohorts.dbh

        if type(env) is np.ndarray:
            env = torch.tensor(env, dtype=self.dtype, device=self.device)

        # Predict env niches for all sites and timesteps
        if pred_growth is None:
            pred_growth = self.nnGrowthEnv(env)
        if pred_morth is None:    
            pred_morth = self.nnMortEnv(env)
        if pred_reg is None:    
            pred_reg = self.nnRegEnv(env)
        
        # Result arrays
        ## dbh / ba / ba*nTRee
        Result = torch.zeros([env.shape[0],env.shape[1],  dbh.shape[2]], device=self.device) # sites, time, species
        ## number of trees
        Result_T = torch.zeros([env.shape[0],env.shape[1],  dbh.shape[2]], device=self.device)

        time = env.shape[1]
        patches = dbh.shape[1]
        sp = self.sp

        # Run over timesteps
        for i in range(time):
            
            # aggregate results
            sites = env.shape[0]
            if i > record_time:
                if dbh.shape[2] != 0:
                    if response == "dbh":
                        ba = dbh
                    elif response == "ba":
                        BA_T(dbh, nTree)
                    else:    
                        # ba * nTree
                        ba = BA_T_P(dbh, nTree)
                    labels = Species
                    samples = (ba * torch.sigmoid((nTree - 0.5)/1e-3)).squeeze(3)  # torch.sigmoid((nTree - 0.5)/1e-3) Baeume da oder nicht
                    samples_T = (nTree * torch.sigmoid((nTree - 0.5)/1e-3)).squeeze(3)
                    tmp_res1, tmp_res2 = self.__aggregate(labels, samples, samples_T, Result[:,i,:], Result_T[:,i,:])
                    Result[:,i,:] = Result[:,i,:] + tmp_res1/patches
                    Result_T[:,i,:] = Result_T[:,i,:] + tmp_res2/patches

            # Model
            # envM = env[:,i,:]
            g = self.growthFP(dbh, Species, pred_growth[:,i,:])
            dbh = dbh+g.unsqueeze(3)
            m = self.mortFP(dbh, Species, nTree, pred_morth[:,i,:]).unsqueeze(3)
            nTree = torch.clamp(nTree - m, min = 0)
            r = self.regFP(dbh, Species, pred_reg[:,i,:])

            # New recruits
            new_dbh = ((r-1+0.1)/1e-3).sigmoid() # TODO: check!!! --> when r 0 dann dbh = 0, ansonsten dbh = 1 dbh[r==0] = 0
            new_nTree = r
            new_Species = torch.arange(0, sp, dtype=torch.int64, device = self.device).unsqueeze(0).repeat(r.shape[0], r.shape[1], 1)

            # Combine
            dbh = torch.concat([dbh, new_dbh.unsqueeze(3)], 2)
            nTree = torch.concat([nTree, new_nTree.unsqueeze(3)], 2)
            Species = torch.concat([Species, new_Species], 2)
            
            # Pad tensors, expensive
            if i % 10 == 0:
                indices = (nTree > 0.5).squeeze().flatten(0, 1)
                org_dim = Species.shape[0:2]
                org_dim_t = torch.tensor(org_dim, dtype = torch.long, device = "cpu")
                dbh = pad_tensors_speed_up(dbh.squeeze(3), indices, org_dim_t).unflatten(0, org_dim).unsqueeze(3)
                nTree = pad_tensors_speed_up(nTree.squeeze(3), indices, org_dim_t).unflatten(0, org_dim).unsqueeze(3)
                Species = pad_tensors_speed_up(Species, indices, org_dim_t).unflatten(0, org_dim)
    
            # Quick and dirty padding [sites, patches, cohorts] [,,4] == 0 trees und wenn jam raus damit
            nTree_ind = nTree.squeeze(3) 
            valid_cols = []
            for col_idx in range(nTree_ind.size(2)):
                if not torch.all(nTree_ind[:,:, col_idx] == 0):
                    valid_cols.append(col_idx)
            dbh = dbh[:,:,valid_cols,:]
            nTree = nTree[:,:,valid_cols,:]
            Species = Species[:,:,valid_cols]
        return (Result, Result_T)


    def fit(self, 
            X: Optional[torch.Tensor]=None, 
            Y: Optional[torch.Tensor]=None, 
            initCohort: CohortMat = None,
            epochs: int=2, 
            batch_size: int=20, 
            learning_rate: float=0.1, 
            start_time: float=0.5, 
            patches: int=50, 
            response: str="dbh"):
        """Fits the model to the given data.
        
        Args:
            X (Optional[torch.Tensor]): The input data. Default is None.
            Y (Optional[torch.Tensor]): The target data. Default is None.
            epochs (int): The number of epochs to train the model. Default is 2.
            batch_size (int): The batch size for training. Default is 20.
            learning_rate (float): The learning rate for the optimizer. Default is 0.1.
            start_time (float): The start time for prediction. Default is 0.5.
            patches (int): The number of patches. Default is 50.
            response (str): The response variable. Default is "dbh", other options are 'ba' or 'ba*nT'.
        
        Returns:
            None
        """
        
        if self.optimizer is None:
            self.optimizer = optim.DiffGrad(params = itertools.chain(*self.parameters), lr = learning_rate) # AdaBound was also good
            # TODO scheduler implementieren
            self.scheduler = torch.optim.lr_scheduler.ExponentialLR(self.optimizer, gamma=0.9)
        time = X.shape[1]
        start_time = round(start_time*time)

        desc='loss: Inf'
        stepSize = np.floor(X.shape[0] / batch_size).astype(int) # type: ignore
        
        if self.device.type == 'cuda':
            torch.cuda.set_device(self.device)
            pin_memory = False
        else:
            pin_memory = True
        Counts = np.round(Y[:,:,:,1]).astype(int)
        indices = np.arange(0, X.shape[0]).astype(int)
        data = torch.utils.data.TensorDataset(torch.tensor(X, dtype=self.dtype, device=torch.device('cpu')),
                                              torch.tensor(Y[:,:,:,:], dtype=self.dtype, device=torch.device('cpu')),
                                              torch.tensor(Counts, dtype=torch.int64, device=torch.device('cpu')),
                                              torch.tensor(indices, dtype=torch.int64, device=torch.device('cpu'))
                                              )
        DataLoader = torch.utils.data.DataLoader(data, batch_size=int(batch_size), shuffle=True, num_workers=0, pin_memory=pin_memory, drop_last=True)
        
        batch_loss = np.zeros(stepSize)
        self.history = np.zeros(epochs)
        
        ep_bar = tqdm(range(epochs),bar_format= "Iter: {n_fmt}/{total_fmt} {l_bar}{bar}| [{elapsed}, {rate_fmt}{postfix}]", file=sys.stdout)
        for epoch in ep_bar:
            for step, (x, y, c, ind) in enumerate(DataLoader):
                self.optimizer.zero_grad()
                
                if initCohort is None:
                    cohorts = CohortMat(dims = [x.shape[0], patches, self.sp, 1], sp = self.sp, device=self.device)
                    nTree = cohorts.nTree
                    Species = cohorts.Species
                    dbh = cohorts.dbh
                else:
                    nTree = (initCohort.nTree[ind,...]).to(self.device, non_blocking=True)
                    Species = (initCohort.Species[ind,...]).to(self.device, non_blocking=True)
                    dbh = (initCohort.dbh[ind,...]).to(self.device, non_blocking=True)
                    
                x = x.to(self.device, non_blocking=True)
                y = y.to(self.device, non_blocking=True)
                c = c.to(self.device, non_blocking=True)
                pred = self.predict(dbh, nTree, Species, x, start_time, response)
                loss1 = torch.nn.functional.mse_loss(y[:, start_time:,:,0], pred[0][:,start_time:,:]).mean() # dbh / ba
                #loss2 = torch.nn.functional.mse_loss(y[:, start_time:,:,1], pred[1][:,start_time:,:]).mean() # nTree
                loss2 = torch.distributions.Poisson(pred[1][:,start_time:,:]+0.001).log_prob(c[:, start_time:,:]).mean().negative() # nTree
                loss = (loss1 + loss2)
                loss.backward()
                self.optimizer.step()
                batch_loss[step] = loss.item()
            #self.scheduler.step()    
            bl = np.mean(batch_loss)
            bl = np.round(bl, 3)
            ep_bar.set_postfix(loss=f'{bl}')
            self.history[epoch] = bl
        torch.cuda.empty_cache()

    def continue_fit(self, 
                     X: Optional[torch.Tensor]=None, 
                     Y: Optional[torch.Tensor]=None, 
                     initCohort: CohortMat = None,
                     epochs: int=2, 
                     batch_size: int=20, 
                     learning_rate: float=0.1, 
                     start_time: float=0.5, 
                     patches: int=50, 
                     response: str="dbh"):
        """Continues the training of the model using the provided data.
        
        Args:
            X (Optional[torch.Tensor]): The input data. Defaults to None.
            Y (Optional[torch.Tensor]): The target data. Defaults to None.
            epochs (int): The number of training epochs. Defaults to 2.
            batch_size (int): The batch size for training. Defaults to 20.
            learning_rate (float): The learning rate for training. Defaults to 0.1.
            start_time (float): The starting time for training. Defaults to 0.5.
            patches (int): The number of patches. Defaults to 50.
            response (str): The response type. Defaults to "dbh".
        
        Returns:
            None
        
        Note:
            This function internally calls the `fit` method to continue the training process.
        """
        
        self.fit(X, 
                 Y, 
                 initCohort,
                 epochs, 
                 batch_size, 
                 learning_rate, 
                 start_time, 
                 patches, 
                 response)

    def gradients(self):
        return [p.grad for p in self.optimizer.param_groups[0]["params"]]

    @property
    def parGrowth(self):
        return self._parGrowth.cpu().data.numpy()

    @property
    def parMort(self):
        return self._parMort.cpu().data.numpy()

    @property
    def parReg(self):
        return self._parReg.cpu().data.numpy()

    @property
    def parGlobal(self):
        return self._parGlobal.cpu().data.numpy()

    @property
    def GrowthEnv(self):
        return [(lambda p: p.data.cpu().numpy())(p) for p in self.nnGrowthEnv.parameters()]    

    @property
    def MortEnv(self):
        return [(lambda p: p.data.cpu().numpy())(p) for p in self.nnMortEnv.parameters()]

    @property
    def RegEnv(self):
        return [(lambda p: p.data.cpu().numpy())(p) for p in self.nnRegEnv.parameters()]    
    
    @property
    def weights(self):
        return [(lambda p: p.data.cpu().numpy())(p) for p in self.parameters] 