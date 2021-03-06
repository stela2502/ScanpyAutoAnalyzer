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

# An Example Analysis

The corner stones here should be flexibility.
2D or 3D analysis - exclude ^MT- and ^RP[LS] or only one of them?
How to deal with this expression if the genes are not excluded?

Make it scriptable - Input either h5ad, loom or CellRanger out files.

```python
CellRangerIn = [ "CELLRANGERDATA" ]
LoomIn = [ "LoomIN" ]
h5file = "./test/test.h5ad"

```

```python
ofile = "h5ad.h5ad"
dimensions = "2"

```

```python
RPexclude = "RPEXCLUDE"
MTexclude = "MTEXCLUDE"

key = "louvain" ## the column to run stats on - leave it like that
key_added = 'KEY_ADDED' ## stats table name and folder name!

```

```python
GOIS = [ "GenesOfInterest" ]
```

```python
import scvelo as scv
import loompy
import scanpy
import igraph
import phate
import glob, os
import pandas as pd
import os
import re
import subprocess
from collections import Counter
import numpy as np
from shutil import rmtree
```

```python
if not CellRangerIn[0] == "CELLRANGERDATA":
    print("reading CellRanger matrix file(s)")
    adata = scanpy.read_10x_mtx( CellRangerIn[0] )
    if len(CellRangerIn) > 1:
        for i in range(1,len(CellRangerIn)):
            tmp = scanpy.read_10x_mtx( CellRangerIn[i] )
            adata = adata0.concatenate(adata, tmp, batch_key='sample')
    adata.var_names_make_unique()
    adata
```

```python
if not LoomIn[0] == "LoomIN":
    print("reading loom file(s)")
    adata = scvelo.read_loom( LoomIn[0] )
    if len(CellRangerIn) > 1:
        for i in range(1,len(LoomIn)):
            tmp = scanpy.read_10x_mtx( LoomIn[i] )
            adata = adata0.concatenate(adata, tmp, batch_key='sample')
    adata.var_names_make_unique()
    adata
```

```python
if not h5file == "./test/test.h5ad":
    print("reading h5 anndata file")
    adata = scanpy.read( h5file )
    adata.var_names_make_unique()
    adata
```

```python
try: adata
except:
    adata = None
if adata is None:
    raise RuntimeError('No input file defined')
```

```python
sampleID = [ str(cellID[-1]) for cellID in adata.obs.index._values ]
Counter(sampleID)
```

```python
#sampleNameDict = {
#    "1": "CD34",
#    "2": "EF_1",
#    "3": "EF_2",
#    "4": "sg_1",
#    "5": "sg_2"
#}
```

```python
#sampleName = [ sampleNameDict[cellID[-1]] for cellID in adata.obs.index._values ]
#Counter(sampleName)
```

```python
#sampleTypeDict = {
#    "1": "CD34",
#    "2": "EF",
#    "3": "EF",
#    "4": "sg",
#    "5": "sg"
#}
```

```python
#sampleType = [ sampleTypeDict[cellID[-1]] for cellID in adata.obs.index._values ]
#Counter(sampleType)
```

```python
df = pd.DataFrame ( {
    "cellID" :   adata.obs.index._values,
    "sampleID" : sampleID  
})
df.index = adata.obs.index
adata.obs= df
```

```python
df
```

```python
RP = re.compile('^RP[SL]')
def testRP(x):
    r= True
    if RP.match(x):
        r = False
    return (r)
RPgenes = [ x  for x in adata.var.index if not testRP(x)]
#print ( Counter(RPgenes) )
adata.obs['RP[LS]sum'] = adata[:,RPgenes].X.sum(axis=1)

if not RPexclude == "RPEXCLUDE":
    OK = [ testRP(x) for x in adata.var.index ]
    adata._inplace_subset_var( np.array(OK) )
    adata
    
```

```python
MT = re.compile('^MT-')
def testMT(x):
    r= True
    if MT.match(x):
        r = False
    return (r)
MTgenes = [ x  for x in adata.var.index if not testMT(x)]
#print ( Counter(RPgenes) )
adata.obs['MTsum'] = adata[:,MTgenes].X.sum(axis=1)

if not MTexclude == "MTEXCLUDE":
    OK = [ testMT(x) for x in adata.var.index ]
    adata._inplace_subset_var( np.array(OK) )
    adata
```

```python
print(adata.n_vars)

```

```python
def testGene(x):
    r= True
    if MT.match(x) or RP.match(x):
        r = False
    return (r)

genes = [ x  for x in adata.var.index if testGene(x)]
adata.obs['geneSum'] = adata[:,genes].X.sum(axis=1)

adata
```

```python
if len(genes) != adata.n_vars:
    ## get rid of cells that have less than 1000 reads from normal genes
    adata = adata[ adata.obs['geneSum'] > 1000]
adata
```

```python
scanpy.pp.filter_cells(adata, min_genes=100 )
adata
```

```python
scanpy.pp.downsample_counts(adata, counts_per_cell= 1000 )
adata
```

```python
scanpy.pp.filter_genes(adata, min_counts=10 )
adata
```

```python
scv.pp.log1p(adata)
scanpy.pp.highly_variable_genes(adata, n_top_genes=3000)
print ( Counter( adata.var['highly_variable']))
adata
```

```python
#scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scanpy.pp.neighbors(adata, n_neighbors=30)
if dimensions == "2":
    dimensions = 2
scanpy.tl.umap(adata,n_components= dimensions)
scanpy.tl.louvain(adata)
```

```python
scv.pl.scatter(adata, color='sampleID', figsize =(15,12), legend_loc='right margin')
```

```python
scv.pl.scatter(adata, color='MTsum', figsize =(7,5), legend_loc='right margin')
```

```python
scv.pl.scatter(adata, color='RP[LS]sum', figsize =(7,5), legend_loc='right margin')
```

```python
scv.pl.scatter(adata, color='geneSum', figsize =(7,5), legend_loc='right margin')
```

```python
if not GOIS[0] == "GenesOfInterest":
    for gene in GOIS:
        scv.pl.scatter(adata, color=gene, cmap="viridis_r", figsize =(15,12), legend_loc='right margin')
```

```python

scanpy.tl.rank_genes_groups(
    adata, 
    groupby   = 'louvain',
    key_added = key_added,
    method    = 'wilcoxon',
)

adata
```

```python
scanpy.pl.rank_genes_groups(adata, key = key_added )
```

```python
scv.pl.scatter(adata, color=['louvain','sampleID'], figsize =(7,5), legend_loc='on data')

scv.pl.scatter(adata, color=['louvain','sampleID'], figsize =(7,5), legend_loc='right margin')
```

```python
diff_results = adata.uns[key_added]
diff_results.keys()
```

```python
columns = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
try: 
    os.mkdir(f"{key_added}") 
except OSError as error: 
    print(error)  
    
for i in range(adata.obs['louvain'].nunique()-1):
    table = {}
    for j in columns:
        table[j] = pd.DataFrame(diff_results[j])[str(i)]
    table = pd.DataFrame(table).to_csv(f"{key_added}/cluster_{i}.csv")
```

```python
if os.path.exists(ofile):
    os.remove(ofile)
adata.write(ofile)
print(ofile)
```

# It has helped tremendousely

to plot the louvain clusters sample specific.
Makes it way more easy to see where the grooups of interest/problem are.

```python
for ID in Counter(sampleID).keys():
    tmp = adata[adata.obs['sampleID'].isin([ID])]
    scv.pl.scatter(tmp,color=['louvain'], figsize =(14,10), legend_loc='right margin', title=f"sampleID {ID}")
```

```python

```
