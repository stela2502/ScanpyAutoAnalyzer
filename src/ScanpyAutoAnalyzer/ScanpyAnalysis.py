#!/usr/bin/env python3.8

import argparse, textwrap
from pathlib import Path
import sys
import re
from os import listdir, makedirs
from os.path import isfile, join, isdir, exists

import os
libFile = join( os.path.dirname(os.path.realpath(__file__)), 'lib', 'ExampleAnalysis.md' )

#>>> parser = argparse.ArgumentParser(prog='PROG')
#>>> parser.add_argument('--foo', nargs='?', help='foo help')
#>>> parser.add_argument('bar', nargs='+', help='bar help')
#>>> parser.print_help()
#usage: PROG [-h] [--foo [FOO]] bar [bar ...]

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

args = parser.parse_args()

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
    txt = txt.replace( "H5FILE", args.input )
elif isfile(args.input) and re.search( '.loom$', args.input ) :
    txt = txt.replace( "LoomIN", args.input )
else:
    ## read the whole folder and check
    mypath = args.input
    CR = [f for f in listdir(mypath) if isdir(join(mypath, f)) and f == "filtered_feature_bc_matrix" ]
    if len(CR) > 0:
        CR = "\", \"".join(CR)
        txt = txt.replace( "CELLRANGERDATA", CR )
    else:
        MT = re.compile('.loom$')
        LM = [f for f in listdir(mypath) if isfile(join(mypath, f)) and MT.match(f) ]
        if len(CR) > 0:
            CR = "\", \"".join(CR)
            txt = txt.replace( "LoomIN", CR )
        else:
            print(f"\nNo infile could be detected there:{args.input}\n", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit()


txt = txt.replace( "OUTFILE", args.name )
txt = txt.replace( "Key_ADDED", args.statsName )
txt = txt.replace( "DIMENSIONS", args.dimensions )

txt = txt.replace( "\\\"MTEXCLUDE\\\"", str(args.mitoEX) )
txt = txt.replace( "\\\"RPEXCLUDE\\\"", str(args.riboEX) )

if not exists( args.outpath ):
    makedirs(  args.outpath )

ofile = f"{args.outpath}/{args.name}.md"
with open( ofile, 'w') as f:
    f.write(txt)

jupyter = f"{args.outpath}/{args.name}.ipynb"
## now run these lines:
#f"jupytext --to notebook {ofile}"
f"jupytext --set-kernel python3 {ofile}"
f"jupytext --to notebook --execute {ofile}"
# name -> "OUTFILE.h5ad"

## infiles - a check would be kind
## One h5ad file
