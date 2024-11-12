An example of how to impute expression values in a 10X output.

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
```

    time: 12.6 s (started: 2024-11-12 08:40:39 +01:00)



```python
import tensorflow as tf
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
```

    2024-11-12 08:40:55.562205: E external/local_xla/xla/stream_executor/cuda/cuda_fft.cc:477] Unable to register cuFFT factory: Attempting to register factory for plugin cuFFT when one has already been registered
    WARNING: All log messages before absl::InitializeLog() is called are written to STDERR
    E0000 00:00:1731397255.731887  554669 cuda_dnn.cc:8310] Unable to register cuDNN factory: Attempting to register factory for plugin cuDNN when one has already been registered
    E0000 00:00:1731397255.831343  554669 cuda_blas.cc:1418] Unable to register cuBLAS factory: Attempting to register factory for plugin cuBLAS when one has already been registered


    Num GPUs Available:  1
    time: 6.25 s (started: 2024-11-12 08:40:54 +01:00)



```python
print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))
```

    anndata==0.10.8
    scanpy==1.10.3
    scvi==1.2.0
    scipy==1.14.1
    pandas==2.2.3
    tensorflow==2.18.0
    time: 703 μs (started: 2024-11-12 08:41:00 +01:00)



```python
def remove_gzipped_files(path):
    # Ensure the directory exists
    if os.path.exists(path):
        for file_name in os.listdir(path):
            # Check if the file is a gzipped file
            if file_name.endswith(".gz"):
                gzipped_file_path = os.path.join(path, file_name)
                os.remove(gzipped_file_path)
                print(f"Removed gzipped file: {gzipped_file_path}")
    else:
        print(f"The directory {path} does not exist.")

def gzip_files_in_directory(path):
    # Ensure the directory exists
    if os.path.exists(path):
        # Run the gzip command on all files in the directory
        subprocess.run(f"gzip {path}/*", shell=True, check=True)
        print(f"All files in {path} have been gzipped.")
    else:
        print(f"The directory {path} does not exist.")

def write_data_to_directory(adata, output_dir):
    # Ensure the directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Remove gzipped files before writing new ones
    remove_gzipped_files(output_dir)

    # Write the matrix data to a .mtx file
    from scipy.io import mmwrite
    mmwrite(os.path.join(output_dir, 'matrix.mtx'), adata.X.T )

    # Write the barcodes (cell identifiers) to a .tsv file
    with open(os.path.join(output_dir, 'barcodes.tsv'), 'w') as f:
        for barcode in adata.obs_names.values:
            f.write(f"{barcode}\n")

    # Write the features (gene identifiers) to a .tsv file
    # Create the dataframe for the features.tsv.gz file
    features_df = pd.DataFrame({
        'gene_id': adata.var.index,  # Gene identifiers (e.g., Ensembl IDs)
        'gene_name': adata.var['gene_name'] if 'gene_name' in adata.var.columns else adata.var.index,  # Gene names (or use gene_id if unavailable)
        'gene_type': adata.var['feature_types'] if 'feature_types' in adata.var.columns else  ['Gene Expression'] * len(adata.var)  # Gene type, usually 'Gene Expression'
    })
    features_df.to_csv(os.path.join(output_dir, 'features.tsv'), sep='\t', index=False, header=False )

    # Optionally, gzip the files in the directory
    gzip_files_in_directory(output_dir)


def impute_expression_data( adata, output_dir ):

    # Step 1: Get the number of CPU cores
    num_cores = os.cpu_count()
    
    # Step 2: Set num_workers to the number of cores + 1
    num_workers = num_cores + 1
    
    # Convert to CSR format
    adata.X = scipy.sparse.csr_matrix(adata.X)
    # Set up the model
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    # Train the model
    model.train()
    # Get the imputed expression values for cells
    adata.X = model.posterior_predictive_sample(adata, n_samples=1).to_scipy_sparse()
    write_data_to_directory( adata, output_dir )


```

    time: 1.46 ms (started: 2024-11-12 08:41:04 +01:00)



```python
data1 = scanpy.read_10x_mtx('Rep1_Ram_014-GEX-ADT-CRISPR/filtered_feature_bc_matrix/')
data1
```




    AnnData object with n_obs × n_vars = 19329 × 32738
        var: 'gene_ids', 'feature_types'



    time: 10.9 s (started: 2024-11-12 08:41:07 +01:00)



```python
sub = data1[data1.obs_names[1:100],:].copy()
sub
```




    AnnData object with n_obs × n_vars = 99 × 32738
        obs: '_scvi_batch', '_scvi_labels'
        var: 'gene_ids', 'feature_types'
        uns: '_scvi_uuid', '_scvi_manager_uuid'



    time: 8.27 ms (started: 2024-11-12 08:43:11 +01:00)



```python
impute_expression_data( sub, 'test_IMPUTED/filtered_feature_bc_matrix/')
```

    GPU available: True (cuda), used: True
    TPU available: False, using: 0 TPU cores
    HPU available: False, using: 0 HPUs
    LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0]
    SLURM auto-requeueing enabled. Setting signal handlers.
    /usr/local/lib/python3.12/dist-packages/lightning/pytorch/trainer/connectors/data_connector.py:424: The 'train_dataloader' does not have many workers which may be a bottleneck. Consider increasing the value of the `num_workers` argument` to `num_workers=47` in the `DataLoader` to improve performance.
    /usr/local/lib/python3.12/dist-packages/lightning/pytorch/loops/fit_loop.py:298: The number of training batches (1) is smaller than the logging interval Trainer(log_every_n_steps=10). Set a lower value for log_every_n_steps if you want to see logs for the training epoch.



    Training:   0%|          | 0/400 [00:00<?, ?it/s]


    `Trainer.fit` stopped: `max_epochs=400` reached.


    All files in test_IMPUTED/filtered_feature_bc_matrix/ have been gzipped.
    time: 6.55 s (started: 2024-11-12 08:43:12 +01:00)



```python
impute_expression_data( data1, 'Rep1_Ram_014-GEX-ADT-CRISPR_IMPUTED/filtered_feature_bc_matrix/')
```

    GPU available: True (cuda), used: True
    TPU available: False, using: 0 TPU cores
    HPU available: False, using: 0 HPUs
    LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0]
    SLURM auto-requeueing enabled. Setting signal handlers.
    /usr/local/lib/python3.12/dist-packages/lightning/pytorch/trainer/connectors/data_connector.py:424: The 'train_dataloader' does not have many workers which may be a bottleneck. Consider increasing the value of the `num_workers` argument` to `num_workers=47` in the `DataLoader` to improve performance.



    Training:   0%|          | 0/400 [00:00<?, ?it/s]


    `Trainer.fit` stopped: `max_epochs=400` reached.


    All files in Rep1_Ram_014-GEX-ADT-CRISPR_IMPUTED/filtered_feature_bc_matrix/ have been gzipped.
    time: 12min 40s (started: 2024-11-12 08:43:23 +01:00)



```python
data1.write( "Rep1_Ram_014-GEX-ADT-CRISPR_IMPUTED.h5ad")
```

    time: 364 ms (started: 2024-11-12 08:56:04 +01:00)



```python
data2 = scanpy.read_10x_mtx('Rep2_Ram_015-GEX-ADT-CRISPR/filtered_feature_bc_matrix/')
data2
```




    AnnData object with n_obs × n_vars = 24681 × 32738
        var: 'gene_ids', 'feature_types'



    time: 14.3 s (started: 2024-11-12 08:56:04 +01:00)



```python
impute_expression_data( data2, 'Rep2_Ram_015-GEX-ADT-CRISPR_IMPUTED/filtered_feature_bc_matrix/')
```

    GPU available: True (cuda), used: True
    TPU available: False, using: 0 TPU cores
    HPU available: False, using: 0 HPUs
    LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0]
    SLURM auto-requeueing enabled. Setting signal handlers.
    /usr/local/lib/python3.12/dist-packages/lightning/pytorch/trainer/connectors/data_connector.py:424: The 'train_dataloader' does not have many workers which may be a bottleneck. Consider increasing the value of the `num_workers` argument` to `num_workers=47` in the `DataLoader` to improve performance.



    Training:   0%|          | 0/324 [00:00<?, ?it/s]


    `Trainer.fit` stopped: `max_epochs=324` reached.


    All files in Rep2_Ram_015-GEX-ADT-CRISPR_IMPUTED/filtered_feature_bc_matrix/ have been gzipped.
    time: 14min 43s (started: 2024-11-12 08:56:19 +01:00)



```python
data2.write( "Rep2_Ram_015-GEX-ADT-CRISPR_IMPUTED.h4ad")
```

    time: 537 ms (started: 2024-11-12 09:11:02 +01:00)



```python

```
