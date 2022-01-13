import os
import re

# This script should handle all slurm interactions.
# the default options
# the scripts
# the starting to the scripts


storage = "~/surmOptions.txt"


def loadDefault(file=None):
	"""
	loads the options from the default script head
	"""

	if file is None:
		file = storage 

	if not os.file.exists( file ):
		print("please copy a woring SLURM script to {storage} - thank you!")
		sys.exit()

	opts = {}
	with os.open(storage, 'r' ) as f:
		lines = f.read()
		for tup in re.findall( r'^#SBATCH -?-(\w+)\s+(.+)$' , lines):
			opts[tup[0]] = tup[1] 

	opts['J'] = 'NA'
	opts['o'] = None
	opts['e'] = None

	return(opts)

def toScriptHead( opts ):
	"""
	convert an opts array of tuples (option, value) into a slurm script head
	"""
	head = "#!/bin/bash\n"
	for key in opts.keys():
		if len(key) ==1:
			head = head + "#SBATCH -{key} {opts[key]}\n"
		else:
			head = head + "#SBATCH --{key} {opts[key]}\n"

	return (head)


def addOptions( parser ):
	"""
	add the command line options to a script's parser object
	"""
	parser.add_argument( "--slurm_n", help="the amount of cores to use", action='store')
	parser.add_argument( "--slurm_N", help="the amount of nodes to use", action='store')
	parser.add_argument( "--slurm_t", help="the max time to run hh::mm::ss", action='store')
	parser.add_argument( "--slurm_A", help="the allocation to use (lsens2018-3-3)", action='store')
	parser.add_argument( "--slurm_p", help="the partitition to use (dell)", action='store')
	parser.add_argument( "--slurm_J", help="the name of this process - keep it short", action='store')
	#parser.add_argument( "--slurm_n", help="the amount of cores to use", action='store')
	#parser.add_argument( "--slurm_n", help="the amount of cores to use", action='store')
	#parser.add_argument( "--slurm_n", help="the amount of cores to use", action='store')


def collect( args ):
	"""
   	collect the slurm options from the command line args
   	"""
   	opts = loadDefault(None)

   	if args.slurm_n not is None:
   		opts['n'] = args.slurm_n
   	else:
   		print( f"missing slurm option slurm_n - set to {opts['n']}")

   	if args.slurm_N not is None:
   		opts['N'] = args.slurm_N
   	else:
   		print( f"missing slurm option slurm_N - set to {opts['N']}")

   	if args.slurm_t not is None:
   		opts['t'] = args.slurm_t
   	else:
   		print( f"missing slurm option slurm_t - set to {opts['t']}")

	if args.slurm_A not is None:
   		opts['A'] = args.slurm_A
   	else:
   		print( f"missing slurm option slurm_A - set to {opts['A']}")

   	if args.slurm_p not is None:
   		opts['p'] = args.slurm_p
   	else:
   		print( f"missing slurm option slurm_p - set to {opts['p']}")

   	if args.slurm_J not is None:
   		opts['J'] = args.slurm_J
   	else:
   		print( f"missing slurm option slurm_J - set to {opts['J']}")

   	return (opts)


def run( cmd, args, outfile ):
	"""
	run the script - give me the command cmd the args and the outfile!
	"""

	script = toScriptHead ( collect(args) )+ f"\n{cmd}\nexit 0\n"

	if os.file.exists( outfile ):
		os.remove(outfile)
	with open( outfile, "w" ) as f:
		f.write( script )

	os.system ( f"sbatch {outfile}")



