import pkg_resources
from os.path import isfile, join, isdir, exists

def getFile ( name ):
	libFile = pkg_resources.resource_filename('ScanpyAutoAnalyzer',join( 'data', 'ExampleAnalysis.md' ))
	return ( libFile )
