# ScanpyAutoAnalyzer

ScanpyAutoAnalyzer is a simple Python package that provides scripts, which use jupyter notebooks and the jupytext package to autor run analysis scripts.

The final idea behind this package is that I have now used the same script to analyze the first steps of a single cell 10X experiment in exactly the same way.
The analysis has always given publication ready results. 

Normally the data in addition made most of the time sense if analyzed with multiple pseudo timelines. This is not normally the case if you follow the best practice scripts.

As I expect to run more analyes in the future - way more - I am planning to automize this script to ultimately become part of a nextflow pipeline. I plann to prepare this as a singularity image and hence instead of creating a simple script in my scope I plann to install this script globally.

Fingers crossed. But it already worked as simple script. So this should not be black magic. 
