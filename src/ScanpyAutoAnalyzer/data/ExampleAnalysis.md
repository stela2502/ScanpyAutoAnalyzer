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

The main collection of functions and example usages. You likely want this one ;-)


```python tags=["parameters"]
CellRangerIn = [ "CELLRANGERDATA" ]
CellRangerH5 = [ "CELLRANGERH5" ]
LoomIn = [ "LoomIN" ]
h5file = "H5FILE"
AlevinIN = [ "ALEVIN"]

ofile = "OUTFILE.h5ad"
dimensions = "DIMENSIONS"


RPexclude = "RPEXCLUDE"
MTexclude = "MTEXCLUDE"

key = "leiden" ## the column to run stats on - leave it like that
key_added = 'KEY_ADDED' ## stats table name and folder name!

GOIS = [ "GenesOfInterest" ]
```

```python
%load_ext autotime

import scvelo as scv
import loompy
import scanpy
import igraph
import phate
import glob, os
import pandas as pd
import re
import subprocess
from collections import Counter
import numpy as np
from shutil import rmtree
import anndata

import h5py
from shutil import copyfile

from ScanpyAutoAnalyzer import *


def copyFiles(files, to):
    for f in files:
        name = os.path.basename( f )
        print( f"copy {f} to {to}" )
        copyfile( f, os.path.join(to, name ) )
    print( "all copied" )
```

```python
print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))
```

```python
help("ScanpyAutoAnalyzer")
```


```python        

if not CellRangerIn[0] == "CELLRANGERDATA":
    def sname( path ):
        path, fname = os.path.split(path)
        file, sname = os.path.split(path)
        return(sname)
    print("reading CellRanger matrix file(s)")
    adata = scanpy.read_10x_mtx( CellRangerIn[0] )
    adata.obs['sname'] = sname(CellRangerIn[0])
    print( f"finished file {CellRangerIn[0]}")
    if len(CellRangerIn) > 1:
        for i in range(1,len(CellRangerIn)):
            tmp = scanpy.read_10x_mtx( CellRangerIn[i] )
            tmp.obs['sname'] = sname(CellRangerIn[i])
            adata = adata.concatenate( tmp, batch_key='sample')
            print( f"finished file {CellRangerIn[i]}")
    adata.var_names_make_unique()
    adata
```

```python
if not CellRangerH5[0] == "CELLRANGERH5":
    print("reading CellRanger H5 file(s)")
    adata = scanpy.read_10x_h5( CellRangerH5[0] )
    if len(CellRangerH5) > 1:
        for i in range(1,len(CellRangerH5)):
            tmp = scanpy.read_10x_h5( CellRangerH5[i] )
            adata = adata.concatenate( tmp, batch_key='sample')
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
            adata = adata.concatenate( tmp, batch_key='sample')
    adata.var_names_make_unique()
    adata
```

```python
if not h5file == "H5FILE":
    print("reading h5 anndata file")
    adata = scanpy.read( h5file )
    adata.var_names_make_unique()
    adata
```

```python    

if not AlevinIN[0] == "ALEVIN":
    def sname( path ):
        path, sname = os.path.split(path)
        return(sname)
    print("reading alevin-fry results")
    adata = pyroe.load_fry( AlevinIN[0] )
    adata.obs['sname'] = sname(AlevinIN[0])
    print( f"finished file {AlevinIN[0]}")
    if len(AlevinIN) > 1:
        for i in range(1,len(AlevinIN)):
            tmp = pyroe.load_fry( AlevinIN[i] )
            tmp.obs['sname'] = sname(AlevinIN[i])
            adata = adata.concatenate( tmp, batch_key='sample')
            print( f"finished file {AlevinIN[i]}")
    adata.var_names_make_unique()
    adata
```

```python
adata.write(ofile)
```


```python
try: adata
except:
    adata = None
if adata is None:
    raise RuntimeError('No input file defined')
```

```python
sampleID = [re.sub("[AGCT]*-", "", x) for x in adata.var.index._values ]
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

adata.obs['sampleID'] = sampleID
scanpy.pp.calculate_qc_metrics( adata, inplace=True)

```

```python
MT = re.compile('^[Mm][Tt]-')
RP = re.compile('^R[pP][SsLl]')

genes = [ x  for x in adata.var.index if testGene(x, MT, RP)]
adata.obs['geneSum'] = adata[:,genes].X.sum(axis=1)

adata
```

```python
dropRP( adata, not RPexclude == "RPEXCLUDE" )
dropMT( adata, not MTexclude == "MTEXCLUDE" )
```

```python
print(adata.n_vars)
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
scanpy.pp.filter_genes(adata, min_counts=1 )
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
scanpy.pp.neighbors(adata)
if dimensions == "DIMENSIONS":
    dimensions = 2
scanpy.tl.umap(adata,n_components= dimensions)
scanpy.tl.leiden(adata)
```

```python
# gene names copied from Seurat
s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"] 
g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"]
mouse_s = [ n.capitalize() for n in s_genes]
mouse_g2m = [ n.capitalize() for n in g2m_genes ]

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
        scv.pl.scatter(adata, color=gene, cmap="viridis_r", figsize =(15,12), size = 20, legend_loc='right margin')
```

```python

scanpy.tl.rank_genes_groups(
    adata, 
    groupby   = 'leiden',
    key_added = key_added,
    method    = 'wilcoxon',
)

adata
```

```python
scanpy.pl.rank_genes_groups(adata, key = key_added )
```

```python
scv.pl.scatter(adata, color=['leiden','sampleID'], figsize =(7,5), legend_loc='on data')

scv.pl.scatter(adata, color=['leiden','sampleID'], figsize =(7,5), legend_loc='right margin')
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
    
for i in range(adata.obs['leiden'].nunique()-1):
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

Or a little easier to use - the write_stats_tables function defined earlier:

```python
write_stats_tables( adata, key_added = "leiden_stats" ,cn ="leiden" )
```


# It has helped tremendousely

to plot the leiden clusters sample specific.
Makes it way more easy to see where the grooups of interest/problem are.

```python
for ID in Counter(sampleID).keys():
    tmp = adata[adata.obs['sampleID'].isin([ID])]
    scv.pl.scatter(tmp,color=['leiden'], figsize =(14,10), legend_loc='right margin', title=f"sampleID {ID}")
```

```python
row = "sname"
col = "leiden"
colS = col.replace(" ", "_")
of = "ThisAnalysis"+"_"+row+"_vs_"+colS+".tsv"
test = pd.DataFrame({name : adata[adata.obs[col] ==name].obs[row].value_counts() for name in adata.obs[col].unique() })
test.to_csv(of, sep="\t")
print ( of )
! head {of}
```

```python
adata.obs['sampleID'] = sampleID
adata.obs.pivot_table(values = "louvian", index = "leiden", columns="sampleID", aggfunc='count')
```

One more thing that was usable in a two column analysis:

a more efficien scaling method:

```python
def scale_values(row):
    total_sum = sum(row)
    scaled_values = {key: value * 100 / total_sum for key, value in row.items()}
    return scaled_values
```
And now create and scale the table -> plot it.

```python
import matplotlib.pyplot as plt

sname_tab = adata.obs.pivot_table(values = "n_genes", index = "leiden", columns="sname", aggfunc='count')
sname_tab_scaled = sname_tab.apply( scale_values, axis=1,  result_type='expand')
ax = sname_tab_scaled.plot.bar(stacked=True, figsize=(16, 8))
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('sname_vs_leiden_stacked_bar_plot.svg', format='svg', bbox_inches='tight')
plt.show()
```

A gene dot plot:

```python
markers = ['Nt5e', 'Cd80', 'Pdcd1lg2', 'Aicda', 'Igha', 'Ighd', 
           'Ighg1', 'Ighg2b', 'Bcl6', 'Nr4a1', 'Bhlhe41', 'Irf7', 
           'Tbx21', 'Cd1d1', 'Prdm1', 'Mki67', 'Cd3e' , 'Cd14']
if not GOIS[0] == "GenesOfInterest":
    scanpy.pl.dotplot(adata, GOIS, groupby='leiden', dendrogram=True)
```

