An example of how to impute expression values in a 10X output. Benefits from GPU compute.

```python tags=["parameters"]
data_path
result_path
ofile
```

```python
%load_ext autotime

import anndata
import scanpy
import scvi
import scipy
import os
import gzip
import subprocess
import pandas as pd

from ScanpyAutoAnalyzer import *

```


```python
print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))
```

```python
help("ScanpyAutoAnalyzer")
```


```python
data1 = scanpy.read_10x_mtx( data_path )
data1
```


```python
impute_expression_data( data1, result_path)
```


```python
data1.write( ofile )
```
