#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system(' python -V')
import sys
print(sys.version)

import stereo as st
import warnings
from scipy.io import mmwrite
import pandas as pd

warnings.filterwarnings('ignore')

import os
os.chdir("/sdc/xxjiang/rapl-sp/out") 
os.getcwd() 

def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory '{directory}' created.")
    else:
        print(f"Directory '{directory}' already exists.")


# In[2]:


data = st.io.read_gef('GD9B/outs/feature_expression/Y00979K6.tissue.gef', bin_size=50)
adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[3]:


adata.X.T[215:220,200:220].todense()


# In[4]:


gene = adata.var
gene['ID']=adata.var.index
gene['symbol']=adata.var['real_gene_name']
gene.head()


# In[ ]:





# In[ ]:





# In[ ]:





# In[22]:


#os.mkdir('GD11-matrix')
fold="GD9B-matrix"
create_directory(fold)
mmwrite(fold+'/matrix.mtx',adata.X.T ,   field="integer")
pd.DataFrame(adata.var)[["ID" ]].to_csv(fold+'/features.tsv', sep='\t', index=True, header=False)
pd.DataFrame(adata.obs.index).to_csv(fold+'/barcodes.tsv', sep='\t', index=False, header=False)
(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv(fold+'/barcodes_pos.tsv', sep='\t', index=True, header=False)


# In[23]:


get_ipython().system('pigz -p 20 -f GD9B-matrix/*mtx')
get_ipython().system('pigz -p 20 -f  GD9B-matrix/*tsv')


# In[ ]:





# In[5]:


adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[ ]:





# In[6]:


import matplotlib.pyplot as plt
 # 绘制结果
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs[["y"]], s=1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[7]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[35]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]].apply(lambda x: -x)+22300, adata.obs['y'] , s=0.1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[42]:


umap=pd.concat([adata.obs[["x"]].apply(lambda x: -x)+22300+300, adata.obs['y']-2900   ], axis=1)


# In[43]:


print(max(list(umap['x'])))
print(min(list(umap['x'])))

print(max(list(umap['y'])))
print(min(list(umap['y'])))


# In[44]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[45]:


import matplotlib.pyplot as plt

# 绘制结果
plt.figure(figsize=(7, 6))
plt.scatter(umap[["x"]], -umap[["y"]], s=5, edgecolor='k', c='blue')
#plt.title("UMAP")
plt.gca().set_facecolor('none')  
plt.gcf().set_facecolor('none')  
#plt.axis('off')
plt.show()


# In[46]:


(umap/20.5).to_csv('GD9B-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False)  
#(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv('GD11-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False) 
get_ipython().system('pigz -p 20 -f GD9B-matrix/*tsv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[48]:


data = st.io.read_gef('GD11/outs/feature_expression/Y00985CD.tissue.gef', bin_size=50) 
adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[49]:


adata.X.T[215:220,200:220].todense()


# In[50]:


gene = adata.var
gene['ID']=adata.var.index
gene['symbol']=adata.var['real_gene_name']
gene.head()


# In[51]:


print(os.getcwd())
os.listdir() 


# In[ ]:


#os.mkdir('GD11-matrix')
create_directory("GD11-matrix")
mmwrite('GD11-matrix/matrix.mtx',adata.X.T ,   field="integer")
pd.DataFrame(adata.var)[["ID"]].to_csv('GD11-matrix/features.tsv', sep='\t', index=True, header=False)
pd.DataFrame(adata.obs.index).to_csv('GD11-matrix/barcodes.tsv', sep='\t', index=False, header=False)
(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv('GD11-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False)


# In[ ]:


get_ipython().system('pigz -p 20 -f GD11-matrix/*mtx')
get_ipython().system('pigz -p 20 -f GD11-matrix/*tsv')


# In[52]:


adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[ ]:





# In[53]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[ ]:





# In[54]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs[["y"]], s=0.1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[55]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs['y'].apply(lambda x: -x)+21950, s=0.1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[90]:


umap=pd.concat([adata.obs[["x"]]-1600, adata.obs['y'].apply(lambda x: -x)+21950+600], axis=1)


# In[91]:


print(max(list(umap['x'])))
print(min(list(umap['x'])))

print(max(list(umap['y'])))
print(min(list(umap['y'])))


# In[92]:


umap['x']


# In[ ]:





# In[ ]:





# In[93]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[94]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(umap[["x"]], umap[["y"]], s=1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[95]:


(umap/20).to_csv('GD11-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False)  
#(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv('GD11-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False)  
get_ipython().system('pigz -p 20 -f GD11-matrix/*tsv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[96]:


data = st.io.read_gef('GD13/outs/feature_expression/Y00985G2.tissue.gef', bin_size=50)
adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[97]:


adata.X.T[215:220,200:220].todense()


# In[98]:


gene = adata.var
gene['ID']=adata.var.index
gene['symbol']=adata.var['real_gene_name']
gene.head()


# In[64]:


#os.mkdir('GD11-matrix')
fold="GD13-matrix"
create_directory(fold)
mmwrite(fold+'/matrix.mtx',adata.X.T ,   field="integer")
pd.DataFrame(adata.var)[["ID" ]].to_csv(fold+'/features.tsv', sep='\t', index=True, header=False)
pd.DataFrame(adata.obs.index).to_csv(fold+'/barcodes.tsv', sep='\t', index=False, header=False)
(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv(fold+'/barcodes_pos.tsv', sep='\t', index=True, header=False)


# In[65]:


get_ipython().system('pigz -p 20 -f GD13-matrix/*mtx')
get_ipython().system('pigz -p 20 -f  GD13-matrix/*tsv')


# In[ ]:





# In[99]:


adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[ ]:





# In[100]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[101]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs[["y"]], s=1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[126]:


umap=pd.concat([adata.obs[["x"]]-1600, adata.obs['y']-1200], axis=1)


# In[127]:


print(max(list(umap['x'])))
print(min(list(umap['x'])))

print(max(list(umap['y'])))
print(min(list(umap['y'])))


# In[128]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[129]:


import matplotlib.pyplot as plt

# 绘制结果
plt.figure(figsize=(3, 3))
plt.scatter(umap[["x"]], umap[["y"]], s=1, edgecolor='k', c='blue')
#plt.title("UMAP")
plt.gca().set_facecolor('none') 
plt.gcf().set_facecolor('none') 
plt.axis('off')
plt.show()


# In[130]:


(umap/20).to_csv('GD13-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False) 
#(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv('GD11-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False) 
get_ipython().system('pigz -p 20 -f GD13-matrix/*tsv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[131]:


data = st.io.read_gef('GD17/outs/feature_expression/Y00986MA.tissue.gef', bin_size=50)
adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[105]:


adata.X.T[215:220,200:220].todense()


# In[106]:


gene = adata.var
gene['ID']=adata.var.index
gene['symbol']=adata.var['real_gene_name']
gene.head()


# In[107]:


#os.mkdir('GD11-matrix')
fold="GD17-matrix"
create_directory(fold)
mmwrite(fold+'/matrix.mtx',adata.X.T ,   field="integer")
pd.DataFrame(adata.var)[["ID" ]].to_csv(fold+'/features.tsv', sep='\t', index=True, header=False)
pd.DataFrame(adata.obs.index).to_csv(fold+'/barcodes.tsv', sep='\t', index=False, header=False)
(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv(fold+'/barcodes_pos.tsv', sep='\t', index=True, header=False)


# In[108]:


get_ipython().system('pigz -p 20 -f GD17-matrix/*mtx')
get_ipython().system('pigz -p 20 -f  GD17-matrix/*tsv')


# In[ ]:





# In[132]:


adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[ ]:





# In[133]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[134]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs[["y"]], s=1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[112]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs['y'].apply(lambda x: -x)+21400 , s=1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[135]:


umap=pd.concat([adata.obs[["x"]]-2800, -adata.obs['y']+21400+650], axis=1)


# In[136]:


print(max(list(umap['x'])))
print(min(list(umap['x'])))

print(max(list(umap['y'])))
print(min(list(umap['y'])))


# In[137]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[138]:


import matplotlib.pyplot as plt

plt.figure(figsize=(3, 3))
plt.scatter(umap[["x"]], umap[["y"]], s=1, edgecolor='k', c='blue')
#plt.title("UMAP")
plt.gca().set_facecolor('none')  
plt.gcf().set_facecolor('none') 
plt.axis('off')
plt.show()


# In[139]:


(umap/20).to_csv('GD17-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False)
#(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv('GD11-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False)
get_ipython().system('pigz -p 20 -f GD17-matrix/*tsv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[140]:


data = st.io.read_gef('GD20/outs/feature_expression/Y00981E3.tissue.gef', bin_size=50)
adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[124]:


adata.X.T[215:220,200:220].todense()


# In[125]:


gene = adata.var
gene['ID']=adata.var.index
gene['symbol']=adata.var['real_gene_name']
gene.head()


# In[126]:


#os.mkdir('GD11-matrix')
fold="GD20-matrix"
create_directory(fold)
mmwrite(fold+'/matrix.mtx',adata.X.T ,   field="integer")
pd.DataFrame(adata.var)[["ID" ]].to_csv(fold+'/features.tsv', sep='\t', index=True, header=False)
pd.DataFrame(adata.obs.index).to_csv(fold+'/barcodes.tsv', sep='\t', index=False, header=False)
(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv(fold+'/barcodes_pos.tsv', sep='\t', index=True, header=False)


# In[127]:


get_ipython().system('pigz -p 20 -f GD20-matrix/*mtx')
get_ipython().system('pigz -p 20 -f  GD20-matrix/*tsv')


# In[ ]:





# In[128]:


adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[ ]:





# In[129]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[130]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs[["y"]], s=1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[131]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs['y'].apply(lambda x: -x)+23000 , s=1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[151]:


umap=pd.concat([adata.obs[["x"]]-1050, -adata.obs['y']+23000+20], axis=1)


# In[152]:


print(max(list(umap['x'])))
print(min(list(umap['x'])))

print(max(list(umap['y'])))
print(min(list(umap['y'])))


# In[153]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[154]:


import matplotlib.pyplot as plt


plt.figure(figsize=(3, 3))
plt.scatter(umap[["x"]], umap[["y"]], s=1, edgecolor='k', c='blue')
#plt.title("UMAP")
plt.gca().set_facecolor('none')  
plt.gcf().set_facecolor('none')  
plt.axis('off')
plt.show()


# In[155]:


(umap/20).to_csv('GD20-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False)  
#(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv('GD11-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False) 
get_ipython().system('pigz -p 20 -f GD20-matrix/*tsv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[156]:


data = st.io.read_gef('GD26/outs/feature_expression/Y00985C2.tissue.gef', bin_size=50)
adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[142]:


adata.X.T[215:220,200:220].todense()


# In[143]:


gene = adata.var
gene['ID']=adata.var.index
gene['symbol']=adata.var['real_gene_name']
gene.head()


# In[144]:


#os.mkdir('GD11-matrix')
fold="GD26-matrix"
create_directory(fold)
mmwrite(fold+'/matrix.mtx',adata.X.T ,   field="integer")
pd.DataFrame(adata.var)[["ID" ]].to_csv(fold+'/features.tsv', sep='\t', index=True, header=False)
pd.DataFrame(adata.obs.index).to_csv(fold+'/barcodes.tsv', sep='\t', index=False, header=False)
(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv(fold+'/barcodes_pos.tsv', sep='\t', index=True, header=False)


# In[145]:


get_ipython().system('pigz -p 20 -f GD26-matrix/*mtx')
get_ipython().system('pigz -p 20 -f  GD26-matrix/*tsv')


# In[ ]:





# In[157]:


adata=st.io.stereo_to_anndata(data,flavor='seurat')
adata


# In[ ]:





# In[158]:


import matplotlib.pyplot as plt
plt.figure(figsize=(3, 3))
plt.scatter(adata.obs[["x"]], adata.obs[["y"]], s=1, edgecolor='k', c='blue')
plt.title("UMAP")
plt.show()


# In[148]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[169]:


umap=pd.concat([adata.obs[["x"]]-1630, adata.obs['y']-1720], axis=1)


# In[170]:


print(max(list(umap['x'])))
print(min(list(umap['x'])))

print(max(list(umap['y'])))
print(min(list(umap['y'])))


# In[171]:


print(max(list(adata.obs['x'])))
print(min(list(adata.obs['x'])))

print(max(list(adata.obs['y'])))
print(min(list(adata.obs['y'])))


# In[172]:


import matplotlib.pyplot as plt

plt.figure(figsize=(3, 3))
plt.scatter(umap[["x"]], umap[["y"]], s=1, edgecolor='k', c='blue')
#plt.title("UMAP")
plt.gca().set_facecolor('none')  
plt.gcf().set_facecolor('none')  
plt.axis('off')
plt.show()


# In[173]:


(umap/20).to_csv('GD26-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False) 
#(pd.DataFrame(adata.obs)[["x","y"]]/20).to_csv('GD11-matrix/barcodes_pos.tsv', sep='\t', index=True, header=False) 
get_ipython().system('pigz -p 20 -f GD26-matrix/*tsv')


# In[ ]:




