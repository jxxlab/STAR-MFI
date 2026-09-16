#!/usr/bin/env python
# coding: utf-8

# In[2]:


# import dependencies
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
#from scanpy.plotting._tools.scatterplots import plot_scatter
import seaborn as sns
import matplotlib.pyplot as plt
#import scvelo as scv
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_versions()
sc.settings.njobs = 32;
sc.set_figure_params(dpi=100, fontsize=16, dpi_save=100,figsize=(4,4), frameon=True )
import os
os.chdir('/sdc/xxjiang/rapl');os.getcwd()


# In[3]:


get_ipython().system(' python   --version')
get_ipython().system(' R        --version')


# In[5]:


adata = sc.read("tr79843m-inte.h5ad")
adata.obsm['X_tsne']=adata.obsm['X_umap']
adata


# In[6]:


sc.pl.scatter ( adata, color=['CX1'],basis='umap',legend_loc='on data',title="")


# In[7]:


# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = "tr79843m-inte.loom"
# path to pyscenic output
f_pyscenic_output = "tr79843m-inte_output.loom"
# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = 'tr79843m-inte_scenic_integrated-output.loom'


# In[8]:


adata.X[0:5,0:10].todense()


# ## SCENIC steps

# ### STEP 1: Gene regulatory network inference, and generation of co-expression modules
# #### Phase Ia: GRN inference using the GRNBoost2 algorithm
# 
# For this step the CLI version of SCENIC is used. This step can be deployed on an High Performance Computing system. We use the counts matrix (without log transformation or further processing) from the loom file we wrote earlier.
# _Output:_ List of adjacencies between a TF and its targets stored in `ADJACENCIES_FNAME`.

# In[9]:


# transcription factors list
f_tfs = "/home/soft/scenicDB/allTFs_hg38.txt" # human


# In[10]:


get_ipython().system('pyscenic grn {f_loom_path_scenic} {f_tfs} -o tr79843m-inte_adj.csv --seed=777 --num_workers 40')


# In[ ]:





# read in the adjacencies matrix:

# In[11]:


adjacencies = pd.read_csv("tr79843m-inte_adj.csv", index_col=False, sep=',');adjacencies.head()


# In[12]:


get_ipython().system('md5sum tr79843m-inte_adj.csv')
#a612803b7e925006e4ded77a329e847b  CS78-EPI5718_adj.csv


# In[ ]:





# ### STEP 2-3: Regulon prediction aka cisTarget from CLI
# 
# For this step the CLI version of SCENIC is used. This step can be deployed on an High Performance Computing system.
# 
# _Output:_ List of adjacencies between a TF and its targets stored in `MOTIFS_FNAME`.

# locations for ranking databases, and motif annotations:

# In[13]:


import glob
# ranking databases
f_db_glob = "/home/soft/scenicDB/hg38/*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )

# motif databases
f_motif_path = "/home/soft/scenicDB/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"


# In[14]:


f_db_names


# Here, we use the `--mask_dropouts` option, which affects how the correlation between TF and target genes is calculated during module creation. It is important to note that prior to pySCENIC v0.9.18, the default behavior was to mask dropouts, while in v0.9.18 and later, the correlation is performed using the entire set of cells (including those with zero expression). When using the `modules_from_adjacencies` function directly in python instead of via the command line, the `rho_mask_dropouts` option can be used to control this.

# In[15]:


get_ipython().system('pyscenic -h')


# In[16]:


get_ipython().system('pyscenic ctx tr79843m-inte_adj.csv     {f_db_names}     --annotations_fname {f_motif_path}     --expression_mtx_fname tr79843m-inte.loom     --output tr79843m-inte_reg.csv     --mask_dropouts     --num_workers 12')


# ### STEP 4: Cellular enrichment (aka AUCell) from CLI
# 
# It is important to check that most cells have a substantial fraction of expressed/detected genes in the calculation of the AUC.
# The following histogram gives an idea of the distribution and allows selection of an appropriate threshold.
# In this plot, a few thresholds are highlighted, with the number of genes selected shown in red text and the corresponding percentile in parentheses).
# See [the relevant section in the R tutorial](https://bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html#build-gene-expression-rankings-for-each-cell) for more information.
# 
# By using the default setting for `--auc_threshold` of `0.05`, we see that **1192** genes are selected for the rankings based on the plot below.

# In[17]:


get_ipython().system('pyscenic aucell    tr79843m-inte.loom     tr79843m-inte_reg.csv     --output {f_pyscenic_output}     --num_workers 20')


# ### Visualization of SCENIC's AUC matrix

# First, load the relevant data from the loom we just created

# In[18]:


import json
import zlib
import base64
#import numpy as np
import pandas as pd
import loompy as lp

f_pyscenic_output = "tr79843m-inte_output.loom"
# collect SCENIC AUCell output
lf = lp.connect( f_pyscenic_output, mode='r', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
auc_mtx.shape 


# In[19]:


auc_mtx.head()


# In[ ]:





# In[20]:


import umap

# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
#dr_umap = adata.obsm['X_umap']
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "tr79843m-inte_umap.txt", sep='\t')
# tSNE
tsne = TSNE( n_jobs=36 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "tr79843m-inte_tsne.txt", sep='\t')


# ## Integrate the output
# 
# Here, we combine the results from SCENIC and the Scanpy analysis into a SCope-compatible loom file

# In[21]:


# scenic output
lf = lp.connect( f_pyscenic_output, mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
#exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
dr_umap = pd.read_csv( 'tr79843m-inte_umap.txt', sep='\t', header=0 , index_col=0 )
dr_tsne = pd.read_csv( 'tr79843m-inte_tsne.txt', sep='\t', header=0, index_col=0 )
###


# Fix regulon objects to display properly in SCope:

# In[22]:


auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
# regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )


# Concatenate embeddings (tSNE, UMAP, etc.)

# In[26]:


tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])

Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
Embeddings_X = pd.concat( [
        pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
        #pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
        dr_tsne['X'] ,
        dr_umap['X']
    ], sort=False, axis=1, join='outer' )
Embeddings_X.columns = ['1','2','3']#,'4']

Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
Embeddings_Y = pd.concat( [
        pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
     #   pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
        dr_tsne['Y'] ,
        dr_umap['Y']
    ], sort=False, axis=1, join='outer' )
Embeddings_Y.columns = ['1','2','3']#,'4']


# Metadata:

# In[29]:


### metadata
metaJson = {}

metaJson['embeddings'] = [
    
    {
        "id": -1,
        "name": "tsne"
    },
    {
        "id": 1,
        "name": "umap"
    },
    #{
    #    "id": 2,
    #    "name": "pca"
    #},
        {
        "id": 3,
        "name": "SCENIC AUC t-SNE"
    },
    {
        "id": 4,
        "name": "SCENIC AUC UMAP"
    },
]

metaJson["clusterings"] = [{
            "id": 0,
            "group": "Scanpy",
            "name": "Scanpy louvain default resolution",
            "clusters": [],
        }]

metaJson["metrics"] = [
        {
            "name": "nCount_RNA"
        }, {
            "name": "nFeature_RNA"
        }#, {
         #   "name": "percent_mt"
        #}
]

metaJson["annotations"] = [
    {
        "name": "CX1",
        "values": list(set( adata.obs['CX1'].values ))
    },
    {
        "name": "Typex1",
        "values": list(set( adata.obs['Typex1'].values ))
            
    },

    {
        "name": "sample",
        "values": list(set( adata.obs['sample'].values ))
    },
    {
        "name": "stage",
        "values": list(set( adata.obs['stage'].values ))
    }
    #{
    #    "name": "Genotype",
    #    "values": list(set(adata.obs['Genotype'].values))
    #},
    #{
    #    "name": "Timepoint",
    #    "values": list(set(adata.obs['Timepoint'].values))
    #},
    #{
    #    "name": "Sample",
    #    "values": list(set(adata.obs['Sample'].values))
    #}
]

# SCENIC regulon thresholds:
metaJson["regulonThresholds"] = rt


# ascdfsdcs  Assemble loom file row and column attributesscfds

# In[30]:


def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


# In[31]:


col_attrs = {
    "CellID": np.array(adata.obs.index),
    "nCount_RNA": np.array(adata.obs['nCount_RNA'].values),
    "nFeature_RNA": np.array(adata.obs['nFeature_RNA'].values),
    "CX1": np.array(adata.obs['CX1'].values),
    "Typex1": np.array(adata.obs['Typex1'].values),
    "sample": np.array(adata.obs['sample'].values),
    "stage": np.array(adata.obs['stage'].values),
    #"Genotype": np.array(adata.obs['Genotype'].values),
    #"Timepoint": np.array(adata.obs['Timepoint'].values),
    #"Sample": np.array(adata.obs['Sample'].values),
    #"percent_mt": np.array(adata.obs['percent_mt'].values),
    "Embedding": dfToNamedMatrix(tsneDF),
    "Embeddings_X": dfToNamedMatrix(Embeddings_X),
    "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
    "RegulonsAUC": dfToNamedMatrix(auc_mtx)
}


row_attrs = {
    "Gene": lf.ra.Gene,
    "Regulons": regulons,
}

attrs = {
    "title": "sampleTitle",
    "MetaData": json.dumps(metaJson),
    "Genome": 'hg38',
    "SCopeTreeL1": "",
    "SCopeTreeL2": "",
    "SCopeTreeL3": ""
}

# compress the metadata field:  
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')


# Create a new loom file, copying the expression matrix from the open loom connection:

# In[32]:


lp.create(
    filename = f_final_loom ,
    layers=lf[:,:],
    row_attrs=row_attrs, 
    col_attrs=col_attrs, 
    file_attrs=attrs
)
lf.close() # close original pyscenic loom file


# In[ ]:





# In[ ]:





# In[ ]:




