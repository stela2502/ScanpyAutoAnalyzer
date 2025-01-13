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

```python tags=["parameters"]
ifile
ofile
resolution=15
cut=0.95
```

```python
%load_ext autotime

import scanpy
import anndata
import glob, os
from collections import Counter
import numpy as np
import pandas as pd 

from ScanpyAutoAnalyzer import *

print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))
```

## The functions I need to merge the over clustered

```python
help("ScanpyAutoAnalyzer")
```


## And finally the code that does the trick:

```python
adata = anndata.read_h5ad( ifile )
anndata
```

```python
resolution = 10 # way too high!
key = f"leiden_{resolution}"
target = f"{key}_{cut}"
scanpy.tl.leiden( adata, resolution = resolution ,key_added=key)
adata.obs['louvain15_merged.95'] = [str(n) for n in adata.obs[key] ] # make sure it is no numeric!
while ( mergeClosest( adata, group = target, cut=cut ) ):
    print ( "w", end=" - ")

reIDgroup(adata, group=target )
```

```python
adata.write( ofile )
! ls -lh {ofile}
```
