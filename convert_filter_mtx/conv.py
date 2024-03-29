import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_array
import h5py
import warnings

#adata = sc.read_h5ad("Chapkin_AHR_Clustered_v1.h5ad")

adata = sc.read_mtx('./matrix.mtx')
adata_bc=pd.read_csv('./barcodes.tsv',header=None)
adata_features=pd.read_csv('./features.tsv',header=None)
print("matrix shape:", adata.X.shape)

adata= adata.T
adata.obs['obs_names']= adata_bc[0].tolist()
adata.obs.index= adata.obs['obs_names']
adata.var['var_names']= adata_features[0].tolist()
adata.var.index= adata.var['var_names']

print("barcodes:",str(len(adata.obs_names)))
print("features:",str(len(adata.var_names)))


# reduce cells to 10000

print("intial adata shape:", adata.X.shape)
#sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, min_genes=500)
print("after filtering cells: ", adata.X.shape)
sc.pp.filter_genes(adata, min_cells=500)
print("after filtering genes:", adata.X.shape)
adata.var_names_make_unique()
# Sub-sampling
sc.pp.subsample(adata,n_obs=9000, random_state=0)
print("after subsampling:", adata.X.shape)

print(adata)

#write_10X_h5(adata,'Chapkin_AHR_filtered.h5')
# Write h5
#adata.write("Chapkin_AHR_filtered.h5")

## write csv
#data1 = adata.X.transpose()
#df = pd.DataFrame(data1)
#df.to_csv('matrix.csv')
#print(df)

# write sparse mtx
data1 = adata.X.transpose()
mat_mtx = csr_array(data1)
mmwrite("matrix2.mtx", mat_mtx)

cells = list(adata.obs_names)
genes = list(adata.var_names)
np.savetxt("features2.tsv", genes, fmt="%s", delimiter="\n")
np.savetxt("barcodes2.tsv", cells, fmt="%s", delimiter="\n")

