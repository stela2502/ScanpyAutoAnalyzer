#!python

import argparse, textwrap
import ScanpyAutoAnalyzer.IO as AA

from pathlib import Path
import sys
import re
from os import listdir, makedirs
from os.path import isfile, join, isdir, exists
import pkg_resources
import numpy as np
import os

libFile = pkg_resources.resource_filename('ScanpyAutoAnalyzer',join( 'data', 'VelocytoAnalalysis.md' ))

# Initiate the parser
parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        This script runs scanpy analyses based on the lib/ExampleAnalysis.ipynb python notebook.
        Each run of this script will produce an anndata outfile, a directory of stats tables and a filled in python notebook.

        In short the data is read from raw files, mitochondrial and ribosomal gene contens will be collected and
        both mitochondrial and ribosomal genes can be excluded from the analysis.

        Cell with less than 1000 UMIs from not mitochondrial and ribosomal genes are excluded.
        Expression is normalized using the scanpy function scanpy.pp.downsample_counts(adata, counts_per_cell= 1000 ).

        Louvain clusters using default parameters will be statisically analzed using the scanpy.tl.rank_genes_groups() function (wilcox).

        ''')
        )

# the initial input/output + names
parser.add_argument("-i", "--input", help="the infile or in path")
parser.add_argument("-o", "--outpath", help="the path where the outfiles should be stored")
parser.add_argument("-n", "--name", help="the name of this analysis")

# the analysis parameters

parser.add_argument("-m", "--mitoEX", help="exclude mitochondrial transcripts", action='store_true') # -> MTEXCLUDE
parser.add_argument("-r", "--riboEX", help="exclude ribosomal transcripts", action='store_true')# -> "RPEXCLUDE"
parser.add_argument("-d", "--dimensions", help="the umap dimensions to use [2,3], default 2", action='store')

parser.add_argument("-s", "--statsName", help="the name of the stats out folder", action='store')# -> KEY_ADDED

parser.add_argument("-g", "--goi", help="list of genes of interest (will be plotted on the data)", nargs='+')# -> "GenesOfInterest"

# Run the created script?
parser.add_argument("-t", "--test", help="do not run the script", action='store_true')# -> "GenesOfInterest"

args = parser.parse_args()

if ( args.test ):
    print( vars(args) )

# read the example script
txt = Path( libFile ).read_text()

problems = False


if args.input is None:
    print("No input files detected", file=sys.stderr)
    problems=True

if args.outpath is None:
    print("No output path given", file=sys.stderr)
    problems=True

if args.name is None:
    print("Missing analysis name", file=sys.stderr)
    problems=True

if args.statsName is None:
    args.statsName = f"deg_{args.name}"
    print(f"statsName set to {args.statsName}", file=sys.stderr)
#    problems=True
if args.dimensions is None:
   args.dimensions = "2"
   print(f"dimensions set to {args.dimensions}", file=sys.stderr)
   #    problems=True

if problems:
    print("\n\n", file=sys.stderr)
    parser.print_help(sys.stderr)
    sys.exit()


# input has three checks:
# one h5ad file -> "H5FILE"
# a list of loom files -> "loomIN"
# a path with filtered_feature_bc_matrix subfolders

if not exists(args.input):
    print(f"\ninput is not a path or file: {args.input}\n", file=sys.stderr)
    parser.print_help(sys.stderr)
    sys.exit()
elif isfile(args.input) and  re.search( '.h5ad$', args.input ) :
    txt = txt.replace( "H5FILE", os.path.abspath(args.input), 1 )
elif isfile(args.input) and re.search( '.loom$', args.input) :
    txt = txt.replace( "LoomIN", os.path.abspath( args.input ), 1 )
elif isdir(args.input):
    ## read the whole folder and check
    ## go first for the h5 files

    #H5 = AA.find_files( args.input, '.h5$' )
    #print("len H5: "+ str(len(H5)) )
    #if len(H5) > 0:
    #    H5 = "\", \"".join(str(v) for v in H5)
    #    txt = txt.replace( "CELLRANGERH5", H5, 1 )
    LM = AA.find_files( args.input, '.loom$' )
    print("len LM: "+ str(len(LM)) )
    CR = AA.find_path ( args.input )
    if len(LM) > 0:
        LM = "\", \"".join(str(v) for v in LM)
        txt = txt.replace( "LoomIN", LM, 1 )
    elif len(CR) > 0:
        ## EPIC CRAP - get rid of outs/filtered...
        print(f"\nCellRanger oiutoput is not supported - use kallisto/bustools: {args.input}\n", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit()
        def rem(p):
            p = os.path.split(p)[0]
            p = os.path.split(p)[0]
            return ( p )
        CR = [ rem(p) for p in CR ]
        CR = "\", \"".join(str(v) for v in CR)
        txt = txt.replace( "CELLRANGERDATA", CR, 1 )
    else:
        print(f"\nNo infile could be detected there: {args.input}\n", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit()


txt = txt.replace( "OUTFILE", args.name , 1)
txt = txt.replace( "KEY_ADDED", args.statsName , 1)
txt = txt.replace( "\"DIMENSIONS\"", args.dimensions, 1 )

txt = txt.replace( "\"MTEXCLUDE\"", str(args.mitoEX), 1 )
txt = txt.replace( "\"RPEXCLUDE\"", str(args.riboEX), 1 )

if not exists( args.outpath ):
    makedirs(  args.outpath )

ofile = f"{args.outpath}/{args.name}.md"
with open( ofile, 'w') as f:
    f.write(txt)

jupyter = f"{args.outpath}/{args.name}.ipynb"
## now run these lines:
#f"jupytext --to notebook {ofile}"

if not args.test:
    os.system(f"jupytext --to notebook {ofile}")
    os.system(f"jupytext --set-kernel python3 {ofile}")
    os.system(f"jupytext --to notebook --execute {ofile}")
# name -> "OUTFILE.h5ad"

## infiles - a check would be kind
## One h5ad file
