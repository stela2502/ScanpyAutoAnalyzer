---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.12.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

My functions to overcluster and then merge the closest clusters if mean expression corelates better than the threshold.


```python
%load_ext autotime

import scanpy
import anndata
import glob, os
from collections import Counter
import numpy as np
import pandas as pd 

print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))
```

## The functions I need to merge the over clustered

```python
def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out


def mergeClosest( adata, group = 'louvain11_merged.9', cut = 0.9 ):
    ### create mean expression per group, correlate all mean groups and scan all groups for corr > cut
    ### the first group with others cor > cut will merge with the best correlating groups available
    ### only one group will be merged at a time
    ### run as
    ### while ( mergeClosest( adata ) ):
    ###   print ( "-", end="")
    print (str(len(Counter(adata.obs[group])))+ " -", end=" ")
    df = grouped_obs_mean( adata, group )
    df = df.corr()
    for i in range(0,df.shape[0]):
        df.at[df.columns[i],df.columns[i]] = 0
    #print("processing groups")
    #print (df)
    for i in range( 0, len(df.columns)):
        i = df.columns[i]
        col = df[i]
        #print ( "max value="+ str( max( col )) )
        if max( col ) >= cut:
            gr = np.array(adata.obs[group])
            m = [df.columns[a] for a in range(0, len(col)) if col[a] == max(col) ]
            m = m[0]
            #print ("max value="+str(max[col])+ " merge " + str(m) + " to "+ str(i)+".")
            for x in range(0,len(adata.obs[group])):
                    if str(gr[x]) == i:
                        #print ( "changed" + str(x)+ " change " + gr[x]  + " to "+ str(m))
                        gr[x] = m
            #type(gr)
            #print("finished")
            adata.obs[group] = gr
            return True
    print ( "no cluster did match to any other cluster with cor > "+ str(cut) )
    return False


def reIDgroup( adata, group ):
    gr = np.array(adata.obs[group])
    n = 0
    for i in set(adata.obs[group]):
        for x in range(0,len(adata.obs[group])):
            if adata.obs[group][x] == i:
                gr[x] = n
        n = n +1
    adata.obs[group+'_renamed'] = gr

```


## And finally the code that does the trick:

```python
scanpy.tl.louvain( adata, resolution =15,key_added='louvain_clusters15')
adata.obs['louvain15_merged.95'] = [str(n) for n in adata.obs['louvain_clusters15'] ] # make sure it is no numeric!
while ( mergeClosest( adata, group = 'louvain15_merged.95', cut=0.95 ) ):
    print ( "w", end=" - ")

reIDgroup(adata, group='louvain15_merged.95' )
```