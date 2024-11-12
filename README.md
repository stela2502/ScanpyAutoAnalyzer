# ScanpyAutoAnalyzer

This package was planned as an automation package and has turned into a scratchpad manager for me.
It is a collection of my scanpy enhancements. It is nothing fancy, just making ones live a little easier.

It also contains my example analysis scripts that are the basis for quite many published analysies.

How to use it is described [in  this jupyter notebook](./ShortUsageExample.ipynb)


# Install

I am no Python power user and therfore 'only' a small pip command to install from github:

```
pip install git+https://github.com/stela2502/ScanpyAutoAnalyzer
```

I am not using anaconda or anything like that, but I am using a Singularity image to capsulate this functionalite.
The Singularity image is not documented - my bad - and "only" available for Aurora users at the time of this writing.

# Aurora users

This library is installed in the singularity image SingSingCell/1.1 on aurora-ls2.
To use this I use command lines like this at the moment:

```
singularity exec -B/projects:/projects -B/local:/mnt /projects/fs1/common/software/SingSingCell/1.1/SingleCells_v1.1.sif ScanpyAnalysis.py  -i /projects/fs3/stefanl/YangLiu/GSE150774 -r -m -o /projects/fs3/stefanl/YangLiu/GSE150774/GSE150774_AutoAnalysis -n NoRiboNoMito
```

This script runs on the frontend, but can easily be moved to the blades by putting this command line into a sbatch script (untested).

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

Especially important is to mount the /projects folder in the singularity image and use absolute paths for in- and out-put info.


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


The script uses <a href="https://github.com/stela2502/ScanpyAutoAnalyzer/blob/main/src/ScanpyAutoAnalyzer/data/ExampleAnalysis.md">this Markdown document</a> to run the analysis.


I hope you use this!