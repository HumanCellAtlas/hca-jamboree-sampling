import warnings

# data processing
with warnings.catch_warnings():
    warnings.simplefilter('ignore', FutureWarning)
    import scanpy.api as sc
import numpy as np
import scanpy.api as sc
import pandas as pd
import os
import scipy.sparse as sp

from sklearn.metrics import pairwise_distances

def _find_cell_indices(adata, cells):
    """Function to find cell indices given index in the AnnData object
    """
    return np.where(adata.obs_names.isin(cells))[0]
        
class SamplingDistanceEstimator:

    def __init__(self, adata: sc.AnnData, ref_set_sigmas: pd.Series, ):
        """Class the computes the distance for the subsample. The score determines the 
        fraction of reference spheres that are unoccupied
        
        :param sc.AnnData adata: Scanpy AnnData object containing the normalized count matrix
        :param pd.Series ref_set_sigmas: Pandas series representing the radii of the 
        reference set
        """
        
        if any(~ref_set_sigmas.index.isin(adata.obs_names)):
            raise ValueError(
                'Some of the cells in the reference set are not in the AnnData object. '
                'Ensure that all the reference cells are in the AnnData object'
            )
        
        self.adata = adata
        self.ref_set = _find_cell_indices(adata, ref_set_sigmas.index)
        self.sigmas = ref_set_sigmas[adata.obs_names[self.ref_set]].values
        
        
    def determine_ref_occupancy(self, test_set: pd.Index, block_size=7500, n_jobs=1):
        """ Function to determine the number of test cells occupying each
        reference sphere
        
        :param pd.Index test_set: Pandas index of test set observation names
        """
        if any(~test_set.isin(self.adata.obs_names)):
            raise ValueError(
                'Some of the cells in the test set are not in the AnnData object. '
                'Ensure that all the test cells are in the AnnData object'
            )
        
        # Test data
        test_data = self.adata[_find_cell_indices(self.adata, test_set),].X
            
        # Compute counts in blocks
        counts = np.zeros(len(self.ref_set))
        blocks = np.linspace(0, len(self.ref_set), 
                             int(len(self.ref_set)/block_size)+1).astype(np.int)
        for b in range(1, len(blocks)):
            test_range = range(blocks[b-1], blocks[b])
            ref_data = self.adata[self.ref_set[test_range],].X
            dists = pairwise_distances(ref_data, test_data, n_jobs=n_jobs)
            counts[test_range] = (dists < self.sigmas[test_range, np.newaxis]).sum(axis=1)
            
        return counts

    
    def determine_distance(self, test_set: pd.Index):
        """ Function to determine the fraction of unoccoupied reference 
        spheres in the test set
        
        :param pd.Index test_set: Pandas index of test set observation names
        """
        if any(~test_set.isin(self.adata.obs_names)):
            raise ValueError(
                'Some of the cells in the test set are not in the AnnData object. '
                'Ensure that all the test cells are in the AnnData object'
            )
        counts = self.determine_ref_occupancy(test_set)
        return np.sum(counts == 0)/len(counts), counts
        
    def determine_distance_from_occupancy(self, ref_sphere_counts):
        """ Function to determine distance from ref spehere counts
        
        :param pd.Index test_set: Pandas index of test set observation names
        """
        return np.sum(ref_sphere_counts==0)/len(ref_sphere_counts)
        
        
        
from typing import List, Tuple
import numpy as np
import h5py

class ImputedSingleCellData:
    
    def __init__(self, US_imputed: np.ndarray, V: np.ndarray, index: np.ndarray, 
                 columns: np.ndarray):
        """Class that contains MAGIC-imputed SVD of a single-cell experiment. Can be used 
        to recover dense imputed data through querying gene names, e.g.: 
        
        >>> im = ImputedSingleCellData.read_h5im('./ica_bone_marrow.h5im')
        >>> im['CD4', 'FOXP3']
        
        :param np.ndarray US_imputed: U_imputed * S, the imputed left eigenvectors of the SVD
          multiplied by the eigenvalues
        :param np.ndarray V: right eigenvectors of the SVD decomposition, contain gene 
          information
        :param np.ndarray index: index (cell ids) of original data matrix
        :param np.ndarray columns: columns (gene ids) of original data matrix
        """
        
        # type checking
        inputs_to_check = {
            'US_imputed': US_imputed,
            'V': V,
            'index': index, 
            'column': columns
        }
        
        for name, array in inputs_to_check.items():
            if not isinstance(array, np.ndarray):
                raise TypeError(
                    'input parameter {param} must be of type np.ndarray not {type}'.format( 
                        param=name, type=type(array))
                )
        
        # verify size of index matches number of cells
        if not US_imputed.shape[0] == index.shape[0]:
            raise ValueError(
                'the number of index values ({index}) must match the number of observations in '
                'US_imputed {US_imputed}'.format(
                    index=index.shape[0], US_imputed=US_imputed.shape[0])
            )
        
        # verify size of columns matches number of genes
        if not V.shape[0] == columns.shape[0]:
            raise ValueError(
                'the number of column values ({columns}) must match the number of variables '
                '(genes) in V\' {V}'.format(columns=columns.shape[0], V=V.shape[0])
            )
        
        self.US_imputed = US_imputed
        self.V = V
        self.index = np.ravel(index.astype('U'))
        self.columns = np.ravel(columns.astype('U'))
    
    @property
    def _columns_to_indices(self) -> pd.Series:
        """private helper function to map column ids to integer indices"""
        return pd.Series(np.arange(len(self.columns)), index=self.columns)
    
    def __getitem__(self, *key: Tuple[str]) -> pd.DataFrame:
        """
        extract imputed gene information
        
        :param Tuple[str] key: any number of string gene ids
        """
        
        # sort out whether a single gene or multiple genes were passed
        if isinstance(key[0], tuple):
            key = key[0]
        
        # translate genes into integer indices
        indices = self._columns_to_indices[list(key)]
        
        # return the imputed data for the requested genes
        return pd.DataFrame(
            data=np.dot(self.US_imputed, self.V[indices, :].T),
            index=self.index,
            columns=list(key)
        )
    
    @classmethod
    def read_h5im(cls, archive_name: str):
        """
        Read an h5im file and return an instance of ImputedSingleCellData
        
        :param str archive_name: filepath to an h5im archive
        """
        with h5py.File(archive_name, mode='r') as f:
            return cls(
                US_imputed=np.array(f['US_imputed']),
                V=np.array(f['V']),
                index=np.array(f['index']).astype('U'),
                columns=np.array(f['columns']).astype('U')
            )
        
    def write_h5im(self, archive_name) -> None:
        """
        write this object to an h5im format archive
        
        :param str archive_name: complete path and filename for the archive
        """
        with h5py.File(archive_name, mode='w') as f:
            f.create_dataset('US_imputed', data=self.US_imputed, compression='lzf')
            f.create_dataset('V', data=self.V, compression='lzf')
            f.create_dataset('index', data=self.index.astype('S'), compression='lzf')
            f.create_dataset('columns', data=self.columns.astype('S'), compression='lzf')


            
def sample_masked(adata, n_cells, cell_mask=None, with_replacement=True):
    """
    :param sc.AnnData adata: scanpy AnnData object to sample from
    :param np.ndarray[Bool] cell_mask: boolean mask over cells. Where false, cells are not 
      considered for sampling
    :param int n_cells: number of cells to draw
    :param bool with_replacement: if True, sample with replacement
    
    :return sc.AnnData: (n_cells x p) AnnData object
    """
    
    # do some usage checking
    if not isinstance(adata, (sp.csr_matrix, sc.AnnData)):
        raise TypeError(
            'required parameter "adata" must be of type scipy.sparse.csr_matrix or '
            'scanpy.AnnData.')
        
    if not isinstance(n_cells, int):
        raise TypeError(
            'required parameter "n_cells" must be of type int and contain '
            'the number of cells to be sampled.')
        
    if cell_mask is not None and not isinstance(cell_mask, np.ndarray):
        raise TypeError(
            'optional parameter "cell_mask" must be of type numpy.ndarray, have the '
            'dtype np.bool, and contain a mask which, when false, disqualifies a cell from '
            'being sampled.')
    
    if cell_mask is not None and not (cell_mask.shape[0] == adata.shape[0]):
        raise ValueError(
            'optional parameter "cell_mask" must have the same length as the number of cells '
            'in "adata".')
    
    if not isinstance(with_replacement, bool):
        raise TypeError(
            'optional parameter "with_replacement" must be of type bool and indicates if '
            'sampling should be done with replacement.')
    
    inds = np.arange(adata.shape[0])
    if cell_mask is not None:
        inds = inds[np.ravel(cell_mask)]
    selected = np.random.choice(inds, size=n_cells, replace=with_replacement)
    
    return adata[selected, :]