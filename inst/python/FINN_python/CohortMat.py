import torch
import numpy as np

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