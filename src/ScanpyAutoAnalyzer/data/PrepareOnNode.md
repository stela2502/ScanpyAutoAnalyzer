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

# Prepare your Session to be run on an Aurora node

Phthon has problems if the working directory is a network share. 
At least I suppose it is like that as I get problems with larger datasets.

```python
from shutil import rmtree
import shutil
import sys


def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)

try:
    origWD
except NameError:
    print(os.environ['SNIC_TMP'])
    path = os.path.basename( os.environ['SNIC_TMP'] )
    wd = f"{os.environ['SNIC_TMP']}/ananylsis"
    ensure_dir( wd )
    origWD = os.getcwd()
    os.chdir( wd )

```

# And after you have finished this

Please copy your results from the wd to origWD
