from SetEnv import *
from picturedrocks import Rocks, ST
from given_functions import ImputedSingleCellData
from given_functions import SamplingDistanceEstimator
from given_functions import sample_masked

def preproc(adata):
    '''
    standardized preprocessing of adata object
    '''
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)

def subset(adata, cell_list):
    '''
    (for some reason) this way is faster than subsetting using 
    cell_list directly
    '''
    all_index = adata.obs_names.values
    cell_inds = [np.where(all_index == x)[0] for x in cell_list]
    return(adata[cell_inds])
    
def nn(adata, n=300):
    '''
    returns weighted connectivity graph as 2D numpy array
    '''
    sc.pp.pca(adata_initial)
    sc.pp.neighbors(adata_initial, n_neighbors=n, knn=True)
    return(adata.uns['neighbors']['connectivities'])

def dens(adata):
    '''
    computes kNN density based on a weighted connectivity graph
    '''
    k = adata.uns['neighbors']['params']['n_neighbors']
    n = adata.shape[0]
    cell_density = []
    for idx in range(n):
        d = np.max(adata.uns['neighbors']['distances'][idx, :].toarray()[0])
        cell_density.append((k - 1) / (2 * n * d))
    cds = pd.DataFrame(cell_density, index = adata.obs_names)
    return(cds)

def dens_sat(adata, orig_data, num_rand, nn):
    '''
    compares kNN density on orig_data and the union of
    orig_data with another random subsample from adata
    of size num_rand
    '''
    all_index = adata.obs_names.values
    orig_index = orig_data.obs_names.values
    rand_index = np.random.choice(list(set(all_index) - set(orig_index)), num_rand, replace=False)
    aggr_index = np.concatenate([rand_index, orig_index])
    aggr_data = subset(adata, aggr_index)
    nn(aggr_data)
    return(dens_comp(orig_data, aggr_data))

def dens_comp(adata1, adata2):
    '''
    compares kNN density on the intersection of two 
    datasets (without normalization)
    '''
    in_common = adata1.obs_names.intersection(adata2.obs_names)
    dens1 = dens(adata1).loc[in_common]
    dens2 = dens(adata2).loc[in_common]
    plt.plot(np.log(dens1), np.log(dens2), s = 1)
    return(np.corrcoef(np.ravel(dens1), np.ravel(dens2)))

# def _log(x):
#     if (x>0):
#         return(np.log(x))
#     else:
#         return(0)

# 

def query_one_gate(adata, cells, ngenes=6):
    '''
    define gate that separates the cell set from the rest and query
    new data set from interface
    gets gating rule from nicetree()
    gets new data from experiment
    '''
    
    return()

def objective(adata, ref_set_sigmas, test):
    sde = SamplingDistanceEstimator(adata, ref_set_sigmas)
    return(determine_distance(sde, test))

def gate_binary_thres(imputed, gate_genes, gate_list):
    ''' Performs an experiment based on gates
    :param archive_name path to ImputedSingleCellData h5im
    :param np.ndarray columns: gate_genes (str) gene ids to use in gating
    :param np.ndarray columns: gate_list (list) list of gate parameters containing:
        gate_thres (scalar) gate thresholds corresponding to gate_genes
        gate_leq (boolean) chose cells which are 
            larger equal to threshold if True, otherwise cells which are smaller
    '''
    cells_inmask = np.full(len(imputed_mask[gate_genes[0]]), 1, dtype='bool') # initialise all cells to lie in mask
    for i,gene in enumerate(gate_genes):
        if(gate_list[1][i]):
            cells_inmask[np.array(imputed_mask[gene] < gate_list[0][i])] = False
        else:
            cells_inmask[np.array(imputed_mask[gene] >= gate_list[0][i])] = False
    return(cells_inmask)

def gate_linearclass(imputed, gate_genes, gate_list):
    '''
    performs linear classification of imputed cells
    gate list contains weights of the linear classifier corresponding to gate_genes
    and the offset in gate_list[1]
    '''
    cells_inmask = np.full(len(imputed[gate_genes[0]]), 0, dtype='bool') # initialise all cells to lie in mask
    classifier = np.zeros(len(imputed[gate_genes[0]])) + gate_list[1] # initialise to offset
    for i,gene in enumerate(gate_genes):
        imputed_mask = imputed[gene]
        classifier = classifier+np.array(imputed_mask[gene])*gate_list[0][i]
    cells_inmask[classifier > 0] = True
    return(cells_inmask)

def biased_experiment(archive_name, adata_raw, gate_fun, gate_genes, gate_list, n_cells):
    ''' Performs an experiment based on gates
    Wraps ImputedSingleCellData and sample_masked
    :param archive_name path to ImputedSingleCellData h5im
    :param np.ndarray columns: gate_genes (str) gene ids to use in gating
    :param np.ndarray columns: gate_fun (str) gene ids to use in gating
    :param np.ndarray columns: gate_list (list) list of gate parameters to be handed
        to gate function
    '''
    imputed = ImputedSingleCellData.read_h5im(archive_name)
    if(gate_fun == 'binary_thres'):
        cells_inmask = gate_binary_thres(imputed, gate_genes, gate_list)
    if(gate_fun == 'linearclass'):
        cells_inmask = gate_linearclass(imputed, gate_genes, gate_list)
    else:
        print('gate_fun not recognized')
            
    adata_new = sample_masked(
        adata=adata_raw,
        n_cells=n_cells,
        cell_mask=cells_inmask,
        with_replacement=False
    )
    return(adata_new)

def unbiased_experiment(adata_raw, n_cells):
    ''' Performs an experiment based on gates
    Wraps ImputedSingleCellData and sample_masked
    :param archive_name path to ImputedSingleCellData h5im
    '''
    cells_inmask = np.full(adata_raw.X.shape[0], 1, dtype='bool')       
    adata_new = sample_masked(
        adata=adata_raw,
        n_cells=n_cells,
        cell_mask=cells_inmask,
        with_replacement=False
    )
    return(adata_new)

###
# functions for selecting louvain groups based on average score
# that s David s easy fix to avoid neighbourhood detection
def target_clusters(adata, score, n=5):
    '''
    return a list of lists of cells, each list represents one cluster
    that was selected for further investigation which is then to be 
    used to define gates for an experiment
    '''
    louvain(adata)
    score_by_lg = ave_score_bygroup(adata, scores)
    # order groups by ave score
    # return lists of member cells by group
    return()

def louvain(adata):
    '''
    wrapper to run louvain clustering on data set
    '''
    sc.tl.louvain(adata)
    
def ave_score_bygroup(adata, scores):
    '''
    get average bias score by louvain group
    '''
    lg_assigns = adata.obs[''].values
    lgs = np.unique(lg_assigns).tolist()
    score_by_lg = [scores[lg_assigns==lg] for lg in lgs]
    return(score_by_lg)
   
# Find connected components of "low density" cells obtained by 
# restricting the original knn graph to the "low density" cells.
def connected_components(adata, ld_inds):

    ''' 
    ld_inds represents an iterable containing the indices of the cells
    in low density areas.  It should probably be a numpy array.

    Return a list of lists of cells.  Each list represents one connected
    component of the cells in low density areas.  
    
    You should have already computed the knn graph for this to work.
    '''

    num_comps, marks = sp.csgraph.connected_components( 
            adata.uns['neighbors']['connectivities'][np.ix_(ld_inds, ld_inds)]
            )


    return _array2clusts(num_comps, marks)

def _array2clusts(num_comps, marks):
    '''private function which turns an array of the type output by
    sp.csgraph.connected_components into a list of lists.  Each list in the list
    contains the indices of the cells in one of the connected components.

    '''

    clusts = []

    for i in range(num_comps): clusts.append([])

    for cell in range(len(marks)):
        clusts[ marks[cell] ].append(cell)

    return clusts

def auto_prepr(some_scanpy_data, nn=200):
    sc.pp.pca(some_scanpy_data)
    sc.pp.neighbors(some_scanpy_data, n_neighbors = nn, knn = True)


# Need to have found nearest neighbours here first
def genClustsCC(adata, lowdens_cells):
    ''' Generate a vector of cluster assignments - a length p
    vector (where p is the number of cells), entry i being j
    means that cell i is in cluster j.  We start indexing 
    at 0, and put all "high density" cells into cluster 0
    (so cluster 1 is the first low density cluster).

    :param lowdens_cells: the list of indices of cells in 
    areas of low density

    '''
    
    myClusts = connected_components(adata, lowdens_cells)
    clustVec = np.zeros(adata.X.shape[0])

    for ind in range(len(myClusts)):
        clustVec[ np.array(lowdens_cells)[myClusts[ind]] ] = ind + 1
    
    return clustVec

def getMarks(rocksObj, clust, num_marks=6):
    '''Outputs the top num_marks markers for cluster clust.
    Run genConsts before running this function
    '''
    if not rocksObj.cs_generated: 
        print("Please generate the consts object before running this method.")
        return
    
    return np.argsort(np.abs(rocksObj.cs_OvA[clust]))[-num_marks:]

def genConsts(rocksObj, alpha=0.0, lambFac=0.0):
    '''
    Put the relevant data into the Rocks object, compressed in a reasonable way.
    Alpha and lambFac can be chosen arbitrarily as long as they are quite small
    (around 0.1 or less).  You will probably see improvements in the selected
    features in this regime, but it is also safe to leave them at 0.  
    '''
    for ind in range(rocksObj.K):
        marks=rocksObj.findCSmarkersForClust(2.0, alpha=alpha, lambFac=lambFac, clust=ind+1)
        
    return 


def getWeights(rocksObj, clust, myMarks):
    '''
    Get the nonzero entrys of the normal vector to the sparse separating hyperplane 
    that separates the cluster of low density cells from the rest of the cells.

    :param rocksObj: A Rocks object that we have already passed through getMarks
    :param clust: The actual value of the cluster (in the list from the function above.
    :param myMarks:  The list of markers created by the getMarks method.
    '''
    absconsts = np.sort(np.abs(rocksObj.cs_OvA[clust]))[::-1]
    vec = ST( np.asarray(rocksObj.cs_OvA[clust]).flatten(), absconsts[myMarks.shape[0]] )
    vec = vec/np.linalg.norm(vec)

    return (vec*rocksObj.cs_scales[clust])[myMarks]

# weight and scale depend on the cluster that you are looking at
def getOffset(adata, marks, weights, scale):

    cent = adata.X.mean(axis=0)
    offset = (scale*np.asarray(cent).flatten())[marks].dot(weights)

    return -1*offset

