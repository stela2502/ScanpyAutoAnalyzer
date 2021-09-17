# ScanpyAutoAnalyzer

ScanpyAutoAnalyzer is a simple Python package that provides scripts, which use jupyter notebooks and the jupytext package to autonaticly run analysis scripts.

The final idea behind this package is that I have now used the same script to analyze the first steps of a single cell 10X experiment four times.
And this analysis has always given publication ready results. 

In addition the data made most of the time sense if analyzed with multiple pseudo timelines. This is not normally the case if you follow the best practice scripts which in my hands always mask developmental steps with cells cycle changes.

As I expect to run more analyes in the future - way more - I am planning to automize this script to ultimately become part of a nextflow pipeline. I plann to prepare this as a singularity image and hence instead of creating a simple script in my scope I plann to install this script globally.


# Aurora users

This library is installed in the singularity image SingSingCell/1.1 on aurora-ls2.
To use this I at the moment use command lines like this:

```
singularity exec -B/projects:/projects -B/local:/mnt /projects/fs1/common/software/SingSingCell/1.1/SingleCells_v1.1.sif ScanpyAnalysis.py  -i /projects/fs3/stefanl/YangLiu/GSE150774 -r -m -o /projects/fs3/stefanl/YangLiu/GSE150774/GSE150774_AutoAnalysis -n NoRiboNoMito
```

This script at the moment still runs on the frontend, but can easily be moved to the blades by putting this command line into a sbatch script.

## more info

The above command line has a logical split into these two elements:

Load the singularity image:
```
singularity exec -B/projects:/projects -B/local:/mnt /projects/fs1/common/software/SingSingCell/1.1/SingleCells_v1.1.sif
```

And run the command
```
ScanpyAnalysis.py  -i /projects/fs3/stefanl/YangLiu/GSE150774 -r -m -o /projects/fs3/stefanl/YangLiu/GSE150774/GSE150774_AutoAnalysis -n NoRiboNoMito
```


Help for this one script:

```
usage: ScanpyAnalysis.py [-h] [-i INPUT] [-o OUTPATH] [-n NAME] [-m] [-r] [-d DIMENSIONS] [-s STATSNAME] [-g GOI [GOI ...]] [-t]

This script runs scanpy analyses based on the lib/ExampleAnalysis.ipynb python notebook. Each run of this script will produce an anndata outfile, a
directory of stats tables and a filled in python notebook. In short the data is read from raw files, mitochondrial and ribosomal gene contens will be
collected and both mitochondrial and ribosomal genes can be excluded from the analysis. Cell with less than 1000 UMIs from not mitochondrial and ribosomal
genes are excluded. Expression is normalized using the scanpy function scanpy.pp.downsample_counts(adata, counts_per_cell= 1000 ). Louvain clusters using
default parameters will be statisically analzed using the scanpy.tl.rank_genes_groups() function (wilcox).

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        the infile or in path
  -o OUTPATH, --outpath OUTPATH
                        the path where the outfiles should be stored
  -n NAME, --name NAME  the name of this analysis
  -m, --mitoEX          exclude mitochondrial transcripts
  -r, --riboEX          exclude ribosomal transcripts
  -d DIMENSIONS, --dimensions DIMENSIONS
                        the umap dimensions to use [2,3], default 2
  -s STATSNAME, --statsName STATSNAME
                        the name of the stats out folder
  -g GOI [GOI ...], --goi GOI [GOI ...]
                        list of genes of interest (will be plotted on the data)
  -t, --test            do not run the script
```


I hope you use this!