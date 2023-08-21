# Integrate several adata objects using scVI

```{python}
import scanpy as sc
import scvi
import torch
from rich import print
from scib_metrics.benchmark import Benchmarker

print("Last run with scvi-tools version:", scvi.__version__)
print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))
```

```{python}
scvi.model.SCVI.setup_anndata(adata, batch_key="Most likely name")
```

```{python}
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
```


```{python}
model.train() # takes a LOOT of time
```


```{python}
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata)

SCVI_MDE_KEY = "X_scVI_MDE"
adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])
```

```{python}
for gene in ["Most likely name", ]:
    sc.pl.embedding(
        adata,
        basis= "X_scVI_MDE" ,
        color=[gene],
        frameon=False,
        ncols=1, color_map="viridis_r",
        save="_"+gene
    )
```