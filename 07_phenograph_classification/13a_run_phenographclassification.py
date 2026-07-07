# By Britney T Forsyth
# Description: Script to run phenograph classification.

# Import packages
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
import re
import os
import sys
from scipy.sparse import csr_matrix
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import matplotlib
import harmony

# Functions 
'''CONSTRUCT AFFINITY MATRIX'''
def _convert_to_affinity(adj, scaling_factors, device, with_self_loops=False):
    """ Convert adjacency matrix to affinity matrix
    """
    N = adj.shape[0]
    rows, cols, dists = find(adj)
    if device == "gpu":
        import cupy as cp
        from cupyx.scipy.sparse import csr_matrix as csr_matrix_gpu
        dists = cp.array(dists) ** 2/(cp.array(scaling_factors.values[rows]) ** 2)
        rows, cols = cp.array(rows), cp.array(cols)
        # Self loops
        if with_self_loops:
            dists = cp.append(dists, cp.zeros(N))
            rows = cp.append(rows, range(N))
            cols = cp.append(cols, range(N))
        aff = csr_matrix_gpu((cp.exp(-dists), (rows, cols)), shape=(N, N)).get()
    elif device == "cpu":
        dists = dists ** 2/(scaling_factors.values[rows] ** 2)

        # Self loops
        if with_self_loops:
            dists = np.append(dists, np.zeros(N))
            rows = np.append(rows, range(N))
            cols = np.append(cols, range(N))
        aff = csr_matrix((np.exp(-dists), (rows, cols)), shape=[N, N])
    return aff

'''CONSTRUCT MUTUAL NEAREST NEIGHBORS GRAPH'''
def _construct_mnn(t1_cells, t2_cells, data_df, n_neighbors,device,n_jobs=-2):
    # FUnction to construct mutually nearest neighbors bewteen two points
    
    if device == "gpu":
        from cuml import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=n_neighbors,
                                metric='cosine')
    elif device == "cpu":
        from sklearn.neighbors import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=n_neighbors,
                                metric='cosine', n_jobs=n_jobs)
    
    print(f't+1 neighbors of t...')
    nbrs.fit(data_df.loc[t1_cells, :].values)
    t1_nbrs = nbrs.kneighbors_graph(
        data_df.loc[t2_cells, :].values, mode='distance')

    print(f't neighbors of t+1...')
    nbrs.fit(data_df.loc[t2_cells, :].values)
    t2_nbrs = nbrs.kneighbors_graph(
        data_df.loc[t1_cells, :].values, mode='distance')

    # Mututally nearest neighbors
    mnn = t2_nbrs.multiply(t1_nbrs.T)
    mnn = mnn.sqrt()
    return mnn

'''COMPUTE SCALING FACTORS'''
def _mnn_scaling_factors(mnn_ka_dists, scaling_factors,device):
    if device == "gpu":
        from cuml import LinearRegression
    else:
        from sklearn.linear_model import LinearRegression
    cells = mnn_ka_dists.index[~mnn_ka_dists.isnull()]

    # Linear model fit
    x = scaling_factors[cells]
    y = mnn_ka_dists[cells]
    lm = LinearRegression()
    lm.fit(x.values.reshape(-1, 1), y.values.reshape(-1, 1))

    # Predict
    x = scaling_factors[mnn_ka_dists.index]
    vals = np.ravel(lm.predict(x.values.reshape(-1, 1)))
    mnn_scaling_factors = pd.Series(vals, index=mnn_ka_dists.index)

    return mnn_scaling_factors

'''CONSTRUCT AFFINITY MATRIX'''
def _mnn_affinity(mnn, mnn_scaling_factors, offset1, offset2, device):
    # Function to convert mnn matrix to affinicty matrix

    # Construct adjacency matrix
    N = len(mnn_scaling_factors)
    x, y, z = find(mnn)
    x = x + offset1
    y = y + offset2
    adj = csr_matrix((z, (x, y)), shape=[N, N])

    # Affinity matrix
    return _convert_to_affinity(adj, mnn_scaling_factors, device, False)

'''
HARMONY AND PALANTIR
'''
def _mnn_ka_distances(mnn, n_neighbors):
    # Function to find distance ka^th neighbor in the mutual nearest neighbor matrix
    ka = int(n_neighbors / 3)
    ka_dists = np.repeat(None, mnn.shape[0])
    x, y, z = find(mnn)
    rows=pd.Series(x).value_counts()
    for r in rows.index[rows >= ka]:
        ka_dists[r] = np.sort(z[x==r])[ka - 1]
    return ka_dists

from harmony import core
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import find, csr_matrix

def harmony_aug_mat_with_pca(projections, timepoints, timepoint_connections, n_neighbors, n_neighbors2):
    # Time point cells and indices
    tp_cells = pd.Series()
    tp_offset = pd.Series()
    offset = 0
    for i in timepoints.unique():
        tp_offset[i] = offset
        tp_cells[i] = list(timepoints.index[timepoints == i])
        offset += len(tp_cells[i])

    # Nearest neighbor graph construction and affinity matrix
    print('Nearest neighbor computation...')
    nbrs = NearestNeighbors(n_neighbors=n_neighbors,
                            metric='cosine', n_jobs=-2).fit(projections.values)

    adj = nbrs.kneighbors_graph(projections.values, mode='distance')
    dists, _ = nbrs.kneighbors(projections.values)
    
    # Scaling factors for affinity matrix construction
    ka = int(n_neighbors / 3)
    scaling_factors = pd.Series(dists[:, ka], index=projections.index)
    
    # Affinity matrix
    nn_aff = _convert_to_affinity(adj, scaling_factors, 'cpu', True)
    n_jobs = -2
    
    # Mututally nearest neighbor affinity matrix
    # Initilze mnn affinity matrix
    N = projections.shape[0]
    full_mnn_aff = csr_matrix(([0], ([0], [0])), [N, N])
    for i in timepoint_connections.index:
        t1, t2 = timepoint_connections.loc[i, :].values
        print(f'Constucting affinities between {t1} and {t2}...')
        # MNN matrix  and distance to ka the distance
        t1_cells = tp_cells[t1]
        t2_cells = tp_cells[t2]
        mnn = _construct_mnn(t1_cells, t2_cells, projections,
                             n_neighbors2, 'cpu', n_jobs)
        
        # MNN Scaling factors
        # Distance to the adaptive neighbor
        ka_dists = pd.Series(0.0, index=t1_cells + t2_cells)
        ka_dists = ka_dists.astype(float)
        # T1 scaling factors
        ka_dists[t1_cells] = _mnn_ka_distances(mnn, n_neighbors2)
        # T2 scaling factors
        ka_dists[t2_cells] = _mnn_ka_distances(mnn.T, n_neighbors2)
        # Scaling factors
        mnn_scaling_factors = pd.Series(0.0, index=projections.index)
#         mnn_scaling_factors[t1_cells] = core._mnn_scaling_factors(
#             ka_dists[t1_cells], scaling_factors,'cpu')
#         mnn_scaling_factors[t2_cells] = core._mnn_scaling_factors(
#             ka_dists[t2_cells], scaling_factors,'cpu')
        mnn_scaling_factors[t1_cells] = _mnn_scaling_factors(
            ka_dists[t1_cells], scaling_factors,'cpu')
        mnn_scaling_factors[t2_cells] = _mnn_scaling_factors(
            ka_dists[t2_cells], scaling_factors,'cpu')
        
        # MNN affinity matrix
#         full_mnn_aff = full_mnn_aff + \
#             core._mnn_affinity(mnn, mnn_scaling_factors,
#                           tp_offset[t1], tp_offset[t2], 'cpu')
        full_mnn_aff = full_mnn_aff + \
            _mnn_affinity(mnn, mnn_scaling_factors,
                          tp_offset[t1], tp_offset[t2], 'cpu')
    # Symmetrize the affinity matrix and return
    aug_aff2 = nn_aff + nn_aff.T + full_mnn_aff + full_mnn_aff.T
    aff2 = nn_aff + nn_aff.T
    return aug_aff2, aff2

'''COMPUTE RANDOM WALK PROBABILITIES'''
def random_walk_probabilities(A, labels):
    D = np.diag(np.sum(A, axis=1))
    L = D - A  # graph laplacian
    seeds = np.array([e != 0 for e in labels], dtype=bool)
    Lu = L[np.invert(seeds),:][:, np.invert(seeds)]  # unlabeled rows, unlabeled cols
    BT = L[np.invert(seeds),:][:, seeds]  # unlabeled rows, labeled cols
    classes = np.unique(labels[labels > 0])
    M = np.zeros((seeds.sum(), len(classes)))
    for k in classes:
        M[labels[seeds] == k, k-1] = 1
    P = np.linalg.lstsq(Lu, np.dot(-BT, M), rcond = None)[0]
    return P

# Set plotting settings
sns.set_style('white')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['image.cmap'] = 'Spectral_r'
matplotlib.rcParams['savefig.dpi'] = 150
matplotlib.style.use("ggplot")
warnings.filterwarnings(action="ignore", module="matplotlib", message="findfont")

import argparse

# Initialize parser
parser = argparse.ArgumentParser(description='Run phenograph classification with adjustable MNN and KNN parameters.')
parser.add_argument('sample', type=str, help='Sample name.')
parser.add_argument('--knn', type=int, default=15, help='Number of neighbors for KNN.')
parser.add_argument('--mnn', type=int, default=30, help='Number of neighbors for MNN.')

# Parse arguments
args = parser.parse_args()

# Assign to variables that can be used globally in this script
sample = args.sample
knn = args.knn
mnn = args.mnn

# Use the arguments
n_neighbors = args.knn
n_neighbors2 = args.mnn

# Set the directory containing the .h5ad files
input_dir = "/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/stratify/"

# Extract sample name from the command line argument
sample = sys.argv[1]

# Construct the file path for the sample
file_path = os.path.join(input_dir, f"{sample}.h5ad")

# Load data
adata_org2 = sc.read_h5ad(file_path)

# Load additional data
adata_patient = sc.read_h5ad('/data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/KG146_Tumor_Mapping_Reference.h5ad')

# Add 'ISC/TA' as a new category if it doesn't exist
# if 'ISC/TA' not in adata_patient.obs['Cell State'].cat.categories:
#     adata_patient.obs['Cell State'] = adata_patient.obs['Cell State'].cat.add_categories(['ISC/TA'])

# Rename the ISC-like cells to ISC/TA cells
# isc_like_cells = adata_patient.obs['Cell State'] == 'ISC-like'
# adata_patient.obs.loc[isc_like_cells, 'Cell State'] = 'ISC/TA'
# ta_like_cells = adata_patient.obs['Cell State'] == 'TA-like'
# adata_patient.obs.loc[ta_like_cells, 'Cell State'] = 'ISC/TA'

# Remove the irrelevant labels
# adata_patient = adata_patient[~adata_patient.obs.loc[:,'Cell State'].isin(['nan','NA','SCC', 'ISC-like', 'TA-like ']),:].copy()
adata_patient = adata_patient[~adata_patient.obs.loc[:,'cell_state_alt1'].isin(['nan', 'NA']),:].copy()

# Prepare matrices
mat_2 = pd.DataFrame(adata_org2.raw.X.todense(), index=adata_org2.raw.obs_names, columns=adata_org2.raw.var_names)
mat_1 = pd.DataFrame(adata_patient.raw.X.todense(), index=adata_patient.raw.obs_names, columns=adata_patient.raw.var_names)

mat_1.index = 'Patient_' + mat_1.index
mat_2.index = 'Org_' + mat_2.index

mat = pd.concat([mat_1, mat_2], join='outer').fillna(0).astype(int)

# Normalize and visualize
adata = sc.AnnData(mat)
adata.raw = adata

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata.obs.loc[:,'log10GenesPerUMI'] = (np.log10(adata.obs.n_genes_by_counts)).div(np.log10(adata.obs.total_counts))

adata.obs['sample'] = [re.sub('_[0-9]+$','',i) for i in adata.obs_names]

sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata, base=2)

adata.layers['counts'] = adata.raw[:,adata.var_names].X

adata.obs['Modality'] = adata.obs.index.str.replace('_.*','', regex=True)

# Feature selection with HVGs
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=500,
    layer="counts",
    flavor="seurat_v3",  # Change to "seurat"
    batch_key="Modality"
)

# Load cell types and markers
cell_types_df = pd.read_excel("/data/chanjlab/CRC_ZFP36L2.092023/ref/GenesSets_MasterSheet.xlsx", sheet_name="HotSpot Modules_Grouped", header=None)

# Create a gene dictionary
genes_dict = {}
for index, row in cell_types_df.iterrows():
    cell_type = row.iloc[0]  # Extract the cell type
    genes = row.iloc[2:].dropna().tolist()
    genes_dict[cell_type] = genes
    
genes_list = []
for index, row in cell_types_df.iterrows():
    genes = row.iloc[2:].dropna().tolist()
    genes_list.extend(genes)

genes_list = [x for x in genes_list if x in adata.var_names]
NE_markers = ['CHGA','CHGB','SYP','ENO2','NCAM1','INSM1']
genes_to_add = list(set(genes_list + NE_markers))

# Set genes_to_add as highly variable in adata.var
adata.var.loc[genes_to_add, 'highly_variable'] = True

# Perform PCA
sc.tl.pca(adata, use_highly_variable=True)

# Define timepoint connections
timepoint_connections = pd.DataFrame(columns=[0, 1])
index = 0
timepoint_connections.loc[index, :] = ['Org', 'Patient']; index += 1

tp = adata.obs.Modality

# Compute augmented and affinity matrix
n_neighbors = 15
n_neighbors2 = 30
pca_merge = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)
aug_mat, aff_mat = harmony_aug_mat_with_pca(pca_merge, tp, timepoint_connections, n_neighbors, n_neighbors2)

# Add matrices to adata
adata.obsm['aug_mat'] = aug_mat.toarray()
adata.obsm['aff_mat'] = aff_mat.toarray()

aug_mat = pd.DataFrame(adata.obsm['aug_mat'], index=adata.obs.index, columns=adata.obs.index)

# Define cell states
bc_intersect = ('Patient_' + adata_patient.obs.index).intersection(adata.obs.index)
adata.obs['cell_state_alt1'] = 'nan'
adata.obs.loc[bc_intersect, 'cell_state'] = adata_patient.obs.loc[bc_intersect.str.replace('Patient_',''), 'cell_state_alt1'].values

cell_types = ['nan'] + adata.obs['cell_state_alt1'].astype('category').cat.categories[:-1].to_list()
ct_labels = adata.obs['cell_state_alt1'].astype('category').cat.reorder_categories(cell_types)
ct_codes = ct_labels.cat.codes

# Compute random walk probabilities
P = random_walk_probabilities(aug_mat.values, ct_codes.values)

# Create ct_prob DataFrame
ct_prob = pd.DataFrame(P, index=aug_mat.index[ct_codes==0], columns=cell_types[1:])
ct_prob.index = ct_prob.index.str.replace('Org_', '')

# Add cd_pred column
cd_pred = ct_prob.idxmax(axis=1)
ct_prob['cell_state_alt1'] = cd_pred

# Save to CSV
# CHANGE THE DATE
ct_prob.to_csv(f"/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/harmony_phenograph/091624/cell_state/adata.probabilities.{sample}_knn{knn}_mnn{mnn}.csv")

# Remove 'Org_' from the row names
ct_prob.index = ct_prob.index.str.replace('Org_', '')
ct_pred = ct_prob.idxmax(axis=1)

# Check if the index matches adata_org2.obs_names
if set(ct_prob.index) == set(adata_org2.obs_names):
    # Add the entire DataFrame to adata_org2.obs
    adata_org2.obs = pd.concat([adata_org2.obs, ct_prob], axis=1)
    adata_org2.obs['Cell State'] = ct_pred
else:
    print("Error: Index mismatch between ct_prob and adata_org2.obs_names")


# Save the adata with predicted cell states and probabilities
file_path = f"/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/harmony_phenograph/adata.probabilities.{sample}.h5ad"
adata_org2.write_h5ad(file_path)

######################################################################################
# # By Britney T Forsyth
# # Description: Script to run phenograph classification.

# # Import packages
# import numpy as np
# import pandas as pd
# import anndata
# import scanpy as sc
# import matplotlib.pyplot as plt
# import re
# import os
# import sys
# from scipy.sparse import csr_matrix
# import seaborn as sns
# import matplotlib.pyplot as plt
# import warnings
# import matplotlib
# import harmony

# # Functions 
# '''CONSTRUCT AFFINITY MATRIX'''
# def _convert_to_affinity(adj, scaling_factors, device, with_self_loops=False):
#     """ Convert adjacency matrix to affinity matrix
#     """
#     N = adj.shape[0]
#     rows, cols, dists = find(adj)
#     if device == "gpu":
#         import cupy as cp
#         from cupyx.scipy.sparse import csr_matrix as csr_matrix_gpu
#         dists = cp.array(dists) ** 2/(cp.array(scaling_factors.values[rows]) ** 2)
#         rows, cols = cp.array(rows), cp.array(cols)
#         # Self loops
#         if with_self_loops:
#             dists = cp.append(dists, cp.zeros(N))
#             rows = cp.append(rows, range(N))
#             cols = cp.append(cols, range(N))
#         aff = csr_matrix_gpu((cp.exp(-dists), (rows, cols)), shape=(N, N)).get()
#     elif device == "cpu":
#         dists = dists ** 2/(scaling_factors.values[rows] ** 2)

#         # Self loops
#         if with_self_loops:
#             dists = np.append(dists, np.zeros(N))
#             rows = np.append(rows, range(N))
#             cols = np.append(cols, range(N))
#         aff = csr_matrix((np.exp(-dists), (rows, cols)), shape=[N, N])
#     return aff

# '''CONSTRUCT MUTUAL NEAREST NEIGHBORS GRAPH'''
# def _construct_mnn(t1_cells, t2_cells, data_df, n_neighbors,device,n_jobs=-2):
#     # FUnction to construct mutually nearest neighbors bewteen two points
    
#     if device == "gpu":
#         from cuml import NearestNeighbors
#         nbrs = NearestNeighbors(n_neighbors=n_neighbors,
#                                 metric='cosine')
#     elif device == "cpu":
#         from sklearn.neighbors import NearestNeighbors
#         nbrs = NearestNeighbors(n_neighbors=n_neighbors,
#                                 metric='cosine', n_jobs=n_jobs)
    
#     print(f't+1 neighbors of t...')
#     nbrs.fit(data_df.loc[t1_cells, :].values)
#     t1_nbrs = nbrs.kneighbors_graph(
#         data_df.loc[t2_cells, :].values, mode='distance')

#     print(f't neighbors of t+1...')
#     nbrs.fit(data_df.loc[t2_cells, :].values)
#     t2_nbrs = nbrs.kneighbors_graph(
#         data_df.loc[t1_cells, :].values, mode='distance')

#     # Mututally nearest neighbors
#     mnn = t2_nbrs.multiply(t1_nbrs.T)
#     mnn = mnn.sqrt()
#     return mnn

# '''COMPUTE SCALING FACTORS'''
# def _mnn_scaling_factors(mnn_ka_dists, scaling_factors,device):
#     if device == "gpu":
#         from cuml import LinearRegression
#     else:
#         from sklearn.linear_model import LinearRegression
#     cells = mnn_ka_dists.index[~mnn_ka_dists.isnull()]

#     # Linear model fit
#     x = scaling_factors[cells]
#     y = mnn_ka_dists[cells]
#     lm = LinearRegression()
#     lm.fit(x.values.reshape(-1, 1), y.values.reshape(-1, 1))

#     # Predict
#     x = scaling_factors[mnn_ka_dists.index]
#     vals = np.ravel(lm.predict(x.values.reshape(-1, 1)))
#     mnn_scaling_factors = pd.Series(vals, index=mnn_ka_dists.index)

#     return mnn_scaling_factors

# '''CONSTRUCT AFFINITY MATRIX'''
# def _mnn_affinity(mnn, mnn_scaling_factors, offset1, offset2, device):
#     # Function to convert mnn matrix to affinicty matrix

#     # Construct adjacency matrix
#     N = len(mnn_scaling_factors)
#     x, y, z = find(mnn)
#     x = x + offset1
#     y = y + offset2
#     adj = csr_matrix((z, (x, y)), shape=[N, N])

#     # Affinity matrix
#     return _convert_to_affinity(adj, mnn_scaling_factors, device, False)

# '''
# HARMONY AND PALANTIR
# '''
# def _mnn_ka_distances(mnn, n_neighbors):
#     # Function to find distance ka^th neighbor in the mutual nearest neighbor matrix
#     ka = int(n_neighbors / 3)
#     ka_dists = np.repeat(None, mnn.shape[0])
#     x, y, z = find(mnn)
#     rows=pd.Series(x).value_counts()
#     for r in rows.index[rows >= ka]:
#         ka_dists[r] = np.sort(z[x==r])[ka - 1]
#     return ka_dists

# from harmony import core
# from sklearn.neighbors import NearestNeighbors
# from scipy.sparse import find, csr_matrix

# def harmony_aug_mat_with_pca(projections, timepoints, timepoint_connections, n_neighbors, n_neighbors2):
#     # Time point cells and indices
#     tp_cells = pd.Series()
#     tp_offset = pd.Series()
#     offset = 0
#     for i in timepoints.unique():
#         tp_offset[i] = offset
#         tp_cells[i] = list(timepoints.index[timepoints == i])
#         offset += len(tp_cells[i])

#     # Nearest neighbor graph construction and affinity matrix
#     print('Nearest neighbor computation...')
#     nbrs = NearestNeighbors(n_neighbors=n_neighbors,
#                             metric='cosine', n_jobs=-2).fit(projections.values)

#     adj = nbrs.kneighbors_graph(projections.values, mode='distance')
#     dists, _ = nbrs.kneighbors(projections.values)
    
#     # Scaling factors for affinity matrix construction
#     ka = int(n_neighbors / 3)
#     scaling_factors = pd.Series(dists[:, ka], index=projections.index)
    
#     # Affinity matrix
#     nn_aff = _convert_to_affinity(adj, scaling_factors, 'cpu', True)
#     n_jobs = -2
    
#     # Mututally nearest neighbor affinity matrix
#     # Initilze mnn affinity matrix
#     N = projections.shape[0]
#     full_mnn_aff = csr_matrix(([0], ([0], [0])), [N, N])
#     for i in timepoint_connections.index:
#         t1, t2 = timepoint_connections.loc[i, :].values
#         print(f'Constucting affinities between {t1} and {t2}...')
#         # MNN matrix  and distance to ka the distance
#         t1_cells = tp_cells[t1]
#         t2_cells = tp_cells[t2]
#         mnn = _construct_mnn(t1_cells, t2_cells, projections,
#                              n_neighbors2, 'cpu', n_jobs)
        
#         # MNN Scaling factors
#         # Distance to the adaptive neighbor
#         ka_dists = pd.Series(0.0, index=t1_cells + t2_cells)
#         ka_dists = ka_dists.astype(float)
#         # T1 scaling factors
#         ka_dists[t1_cells] = _mnn_ka_distances(mnn, n_neighbors2)
#         # T2 scaling factors
#         ka_dists[t2_cells] = _mnn_ka_distances(mnn.T, n_neighbors2)
#         # Scaling factors
#         mnn_scaling_factors = pd.Series(0.0, index=projections.index)
# #         mnn_scaling_factors[t1_cells] = core._mnn_scaling_factors(
# #             ka_dists[t1_cells], scaling_factors,'cpu')
# #         mnn_scaling_factors[t2_cells] = core._mnn_scaling_factors(
# #             ka_dists[t2_cells], scaling_factors,'cpu')
#         mnn_scaling_factors[t1_cells] = _mnn_scaling_factors(
#             ka_dists[t1_cells], scaling_factors,'cpu')
#         mnn_scaling_factors[t2_cells] = _mnn_scaling_factors(
#             ka_dists[t2_cells], scaling_factors,'cpu')
        
#         # MNN affinity matrix
# #         full_mnn_aff = full_mnn_aff + \
# #             core._mnn_affinity(mnn, mnn_scaling_factors,
# #                           tp_offset[t1], tp_offset[t2], 'cpu')
#         full_mnn_aff = full_mnn_aff + \
#             _mnn_affinity(mnn, mnn_scaling_factors,
#                           tp_offset[t1], tp_offset[t2], 'cpu')
#     # Symmetrize the affinity matrix and return
#     aug_aff2 = nn_aff + nn_aff.T + full_mnn_aff + full_mnn_aff.T
#     aff2 = nn_aff + nn_aff.T
#     return aug_aff2, aff2

# '''COMPUTE RANDOM WALK PROBABILITIES'''
# def random_walk_probabilities(A, labels):
#     D = np.diag(np.sum(A, axis=1))
#     L = D - A  # graph laplacian
#     seeds = np.array([e != 0 for e in labels], dtype=bool)
#     Lu = L[np.invert(seeds),:][:, np.invert(seeds)]  # unlabeled rows, unlabeled cols
#     BT = L[np.invert(seeds),:][:, seeds]  # unlabeled rows, labeled cols
#     classes = np.unique(labels[labels > 0])
#     M = np.zeros((seeds.sum(), len(classes)))
#     for k in classes:
#         M[labels[seeds] == k, k-1] = 1
#     P = np.linalg.lstsq(Lu, np.dot(-BT, M), rcond = None)[0]
#     return P

# # Set plotting settings
# sns.set_style('white')
# matplotlib.rcParams['figure.figsize'] = [4, 4]
# matplotlib.rcParams['figure.dpi'] = 100
# matplotlib.rcParams['image.cmap'] = 'Spectral_r'
# matplotlib.rcParams['savefig.dpi'] = 150
# matplotlib.style.use("ggplot")
# warnings.filterwarnings(action="ignore", module="matplotlib", message="findfont")

# # Set the directory containing the .h5ad files
# input_dir = "/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/stratify/"

# # Extract sample name from the command line argument
# sample = sys.argv[1]

# # Construct the file path for the sample
# file_path = os.path.join(input_dir, f"{sample}.h5ad")

# # Load data
# adata_org2 = sc.read_h5ad(file_path)

# # Load additional data
# adata_patient = sc.read_h5ad('/data/chanjlab/CRC_ZFP36L2.092023/ref/KG146_Tumor_Mapping_Reference_alt.h5ad')

# # Add 'ISC/TA' as a new category if it doesn't exist
# # if 'ISC/TA' not in adata_patient.obs['Cell State'].cat.categories:
# #     adata_patient.obs['Cell State'] = adata_patient.obs['Cell State'].cat.add_categories(['ISC/TA'])

# # Rename the ISC-like cells to ISC/TA cells
# # isc_like_cells = adata_patient.obs['Cell State'] == 'ISC-like'
# # adata_patient.obs.loc[isc_like_cells, 'Cell State'] = 'ISC/TA'
# # ta_like_cells = adata_patient.obs['Cell State'] == 'TA-like'
# # adata_patient.obs.loc[ta_like_cells, 'Cell State'] = 'ISC/TA'

# # Remove the irrelevant labels
# # adata_patient = adata_patient[~adata_patient.obs.loc[:,'Cell State'].isin(['nan','NA','SCC', 'ISC-like', 'TA-like ']),:].copy()
# adata_patient = adata_patient[~adata_patient.obs.loc[:,'cell_state_alt1'].isin(['nan', 'NA']),:].copy()

# # Prepare matrices
# mat_2 = pd.DataFrame(adata_org2.raw.X.todense(), index=adata_org2.raw.obs_names, columns=adata_org2.raw.var_names)
# mat_1 = pd.DataFrame(adata_patient.raw.X.todense(), index=adata_patient.raw.obs_names, columns=adata_patient.raw.var_names)

# mat_1.index = 'Patient_' + mat_1.index
# mat_2.index = 'Org_' + mat_2.index

# mat = pd.concat([mat_1, mat_2], join='outer').fillna(0).astype(int)

# # Normalize and visualize
# adata = sc.AnnData(mat)
# adata.raw = adata

# sc.pp.calculate_qc_metrics(adata, inplace=True)
# adata.obs.loc[:,'log10GenesPerUMI'] = (np.log10(adata.obs.n_genes_by_counts)).div(np.log10(adata.obs.total_counts))

# adata.obs['sample'] = [re.sub('_[0-9]+$','',i) for i in adata.obs_names]

# sc.pp.filter_genes(adata, min_cells=1)
# sc.pp.normalize_total(adata)
# sc.pp.log1p(adata, base=2)

# adata.layers['counts'] = adata.raw[:,adata.var_names].X

# adata.obs['Modality'] = adata.obs.index.str.replace('_.*','', regex=True)

# # Feature selection with HVGs
# sc.pp.highly_variable_genes(
#     adata,
#     n_top_genes=500,
#     layer="counts",
#     flavor="seurat_v3",  # Change to "seurat"
#     batch_key="Modality"
# )

# # Load cell types and markers
# cell_types_df = pd.read_excel("/data/chanjlab/CRC_ZFP36L2.092023/ref/GenesSets_MasterSheet.xlsx", sheet_name="HotSpot Modules_Grouped", header=None)

# # Create a gene dictionary
# genes_dict = {}
# for index, row in cell_types_df.iterrows():
#     cell_type = row.iloc[0]  # Extract the cell type
#     genes = row.iloc[2:].dropna().tolist()
#     genes_dict[cell_type] = genes
    
# genes_list = []
# for index, row in cell_types_df.iterrows():
#     genes = row.iloc[2:].dropna().tolist()
#     genes_list.extend(genes)

# genes_list = [x for x in genes_list if x in adata.var_names]
# NE_markers = ['CHGA','CHGB','SYP','ENO2','NCAM1','INSM1']
# genes_to_add = list(set(genes_list + NE_markers))

# # Set genes_to_add as highly variable in adata.var
# adata.var.loc[genes_to_add, 'highly_variable'] = True

# # Perform PCA
# sc.tl.pca(adata, use_highly_variable=True)

# # Define timepoint connections
# timepoint_connections = pd.DataFrame(columns=[0, 1])
# index = 0
# timepoint_connections.loc[index, :] = ['Org', 'Patient']; index += 1

# tp = adata.obs.Modality

# # Compute augmented and affinity matrix
# n_neighbors = 15
# n_neighbors2 = 30
# pca_merge = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)
# aug_mat, aff_mat = harmony_aug_mat_with_pca(pca_merge, tp, timepoint_connections, n_neighbors, n_neighbors2)

# # Add matrices to adata
# adata.obsm['aug_mat'] = aug_mat.toarray()
# adata.obsm['aff_mat'] = aff_mat.toarray()

# aug_mat = pd.DataFrame(adata.obsm['aug_mat'], index=adata.obs.index, columns=adata.obs.index)

# # Define cell states
# bc_intersect = ('Patient_' + adata_patient.obs.index).intersection(adata.obs.index)
# adata.obs['cell_state_alt1'] = 'nan'
# adata.obs.loc[bc_intersect, 'cell_state_alt1'] = adata_patient.obs.loc[bc_intersect.str.replace('Patient_',''), 'cell_state_alt1'].values

# cell_types = ['nan'] + adata.obs['cell_state_alt1'].astype('category').cat.categories[:-1].to_list()
# ct_labels = adata.obs['cell_state_alt1'].astype('category').cat.reorder_categories(cell_types)
# ct_codes = ct_labels.cat.codes

# # Compute random walk probabilities
# P = random_walk_probabilities(aug_mat.values, ct_codes.values)

# # Create ct_prob DataFrame
# ct_prob = pd.DataFrame(P, index=aug_mat.index[ct_codes==0], columns=cell_types[1:])
# ct_prob.index = ct_prob.index.str.replace('Org_', '')

# # Add cd_pred column
# cd_pred = ct_prob.idxmax(axis=1)
# ct_prob['cell_state_alt1'] = cd_pred

# # Save to CSV
# # CHANGE THE DATE
# ct_prob.to_csv(f"/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/harmony_phenograph/091324/cell_state_alt1/adata.probabilities.{sample}.csv")

# # # Remove 'Org_' from the row names
# # ct_prob.index = ct_prob.index.str.replace('Org_', '')
# # ct_pred = ct_prob.idxmax(axis=1)

# # # Check if the index matches adata_org2.obs_names
# # if set(ct_prob.index) == set(adata_org2.obs_names):
# #     # Add the entire DataFrame to adata_org2.obs
# #     adata_org2.obs = pd.concat([adata_org2.obs, ct_prob], axis=1)
# #     adata_org2.obs['Cell State'] = ct_pred
# # else:
# #     print("Error: Index mismatch between ct_prob and adata_org2.obs_names")


# # # Save the adata with predicted cell states and probabilities
# # file_path = f"/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/harmony_phenograph/adata.probabilities.{sample}.h5ad"
# # adata_org2.write_h5ad(file_path)