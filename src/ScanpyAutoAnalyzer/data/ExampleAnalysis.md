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
CellRangerH5 = [ "CELLRANGERH5" ]
LoomIn = [ "LoomIN" ]
h5file = "H5FILE"

```

```python
ofile = "OUTFILE.h5ad"
dimensions = "DIMENSIONS"

```

```python
RPexclude = "RPEXCLUDE"
MTexclude = "MTEXCLUDE"

key = "louvain" ## the column to run stats on - leave it like that
key_added = 'KEY_ADDED' ## stats table name and folder name!

```

```python
GOIS = [ "GenesOfInterest" ]

onNode = "ONNODE"
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
import os
import re
import subprocess
from collections import Counter
import numpy as np
from shutil import rmtree
import os

import h5py
from shutil import copyfile


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
def addSampleData( adata, sname, sampleInfo ):
    """Add the sample information obtained by either running demux10x or quantifyRhapsody to a anndata object"""
    adata.obs['sname'] = sname
    adata.obs_names = [ sname +"."+ n for n in adata.obs_names]
    sampleInfo = pd.read_csv( sampleInfo, sep="\t",index_col= 0)
    sampleInfo['CellIDs'] = [ sname +"."+n+"-1" for n in sampleInfo.index.values]
    sampleInfo.set_index('CellIDs', inplace=True)
    adata.obs['CellIDs'] = adata.obs_names
    obs = adata.obs.merge( sampleInfo, left_index=True, right_index=True)
    print(obs.shape)
    obs.drop_duplicates()
    obs.drop( ['CellIDs'], axis=1, inplace=True)
    print(obs.shape)
    adata = adata[ obs.index.values]
    adata.obs = obs
    ADT = pd.read_csv( 'Cell2Sample.ADT'+sname+'.R1_001.fastq.gz.tsv', sep="\t",index_col= 0)
    ADT['CellIDs'] = [ sname+"."+n+"-1" for n in ADT.index.values]
    ADT.set_index( 'CellIDs' , inplace=True)
    
    obs = adata.obs.merge( ADT, left_index=True, right_index=True )
    #print(obs)
    #print (adata)
    adata = adata[ obs.index.values]
    #print(adata)
    adata.obs = obs
    return (adata)

def testRP(x, RP):
    r= True
    if RP.match(x):
        r = False
    return (r)

def dropRP ( adata, drop = True ):
    """ remove all genes matching to /^R[pP][SsLl]/ """
    RP = re.compile('^R[pP][SsLl]')
    RPgenes = [ x  for x in adata.var.index if not testRP(x, RP)]
    print ( Counter(RPgenes) )
    if  len(RPgenes) > 0:
        adata.obs['RP[LS]sum'] = adata[:,RPgenes].X.sum(axis=1)
    else:
        adata.obs['RP[LS]sum'] = 0
    if drop and len(RPgenes) > 0:
        OK = [ testRP(x) for x in adata.var.index ]
        print(Counter(OK))
        adata._inplace_subset_var( np.array(OK) )
    adata

def testMT(x, MT):
    r= True
    if MT.match(x):
        r = False
    return (r)

def dropMT (adata, drop = True ) :
    """ remove all genes matching to /^[Mm][Tt]-/ """
    MT = re.compile('^[Mm][Tt]-')
    MTgenes = [ x  for x in adata.var.index if not testMT(x, MT)]
    print ( Counter(MTgenes) )
    if len(MTgenes) > 0:
        adata.obs['MTsum'] = adata[:,MTgenes].X.sum(axis=1)
    else:
        adata.obs['MTsum'] = 0
    if drop and len(MTgenes) > 0:
        OK = [ testMT(x) for x in adata.var.index ]
        adata._inplace_subset_var( np.array(OK) )


def testGene(x, MT, RP):
    r= True
    if MT.match(x) or RP.match(x):
        r = False
    return (r)
```


```python
if not CellRangerIn[0] == "CELLRANGERDATA":
        def sname( path ):
    path, fname = os.path.split(path)
    file, sname = os.path.split(path)
    return(sname)

if not CellRangerIn[0] == "CELLRANGERDATA":
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
sampleID = [re.sub("[AGCT]*-", "", x) for x in adata.obs.index._values ]
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
scanpy.pp.neighbors(adata, n_neighbors=30)
if dimensions == "DIMENSIONS":
    dimensions = 2
scanpy.tl.umap(adata,n_components= dimensions)
scanpy.tl.louvain(adata)
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
row = "sname"
col = "louvain"
colS = col.replace(" ", "_")
of = "ThisAnalysis"+"_"+row+"_vs_"+colS+".tsv"
test = pd.DataFrame({name : adata[adata.obs[col] ==name].obs[row].value_counts() for name in adata.obs[col].unique() })
test.to_csv(of, sep="\t")
print ( of )
! head {of}
```
