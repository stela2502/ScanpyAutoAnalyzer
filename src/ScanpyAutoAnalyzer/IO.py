import sys
import re
from os import listdir, makedirs
from os.path import isfile, join, isdir, exists
import pkg_resources
import numpy as np
import os

def checkPath(path):
    OK = None
    try:
        OK = not isfile( path )
    except err:
        print( "got an error from isdir!")
        pass
    return ( OK )


def find_path( path, pattern="filtered_feature_bc_matrix" ):
    paths = listdir(path)
    ret = [ ]
    tmp = None
    for f in paths:
        #print(f + " ispath? "+ str( checkPath(join(path, f)) ) )
        if f == pattern and checkPath(join(path, f)):
            ret.append( os.path.abspath( join( path,f) ) )
        elif checkPath(join(path, f)):
            tmp = find_path( join(path,f) )
            if len(tmp) > 0:
                for i in range(len(tmp)):
                   ret.append( tmp[i] )
    return (np.array(ret).flatten())


def find_files( path, pattern ):
    paths = listdir(path)
    ret = []
    for f in paths:
        if re.search(pattern, f ):
            ret.append( os.path.abspath( join( path,f) ) )
        elif checkPath(join(path, f)):
            tmp = find_files( join(path,f), pattern )
            if (  len(tmp) >0 ):
                for i in range(len(tmp)):
                   ret.append( tmp[i] )    
    return (np.array(ret).flatten())