#!/usr/bin/env python
# A script for adding the virtual charged sites for a tip4p/05 model

import MDAnalysis as mda
from add_dummy import Dummys

mda_ver = mda.__version__.split('-')[0] # to remove -dev0 if a development version is being used
print('MDAnalysis version: ' + mda_ver)

# if mda_ver < '2.0.0':
#     Warning('version too low for mpi')
#     print('Warning: version too low for mpi')

# for handling command line arguemnts 

import argparse

parser = argparse.ArgumentParser(description='Generates the positions of dummy atoms of water molecules with a tip4p geometry. This can be done for just a single frame, or an entire trajectory')
parser.add_argument("topfile",help='The topology file name',type=str)
# optional trajectory argument -t
parser.add_argument("-f","--trajfile", default = None,type=str,metavar='Trajectory')
parser.add_argument("-v","--verbose", default = False,action="store_true",help='Increase verbosity')
parser.add_argument("-n","--ncores",default=1,type=int)
parser.add_argument("-q","--charges",default=(-1.1128,0.5564),type=float,
                    help = "Specify the oxygen and then hydrogen charges. Defaults to TIP4P/2005 values",
                    nargs = 2
                    )
parser.add_argument("-d","--Mdist", default= 0.1546, type=float,
                    help = " The distance in box units from the Oxygen atom to the massless virtual site. Defaults to TIP4P/2005 value of 0.1546 Angstrom")
parser.add_argument("--Mpickle", default='m_vects.pkl')
parser.add_argument("-O","--type_O",default=1)
parser.add_argument("-H","--type_H",default=2)
parser.add_argument("-M","--type_M",default='M')
parser.add_argument("--overwrite",default=False, action='store_true')
parser.add_argument("--fname",default=None,help='output file prefix defaults to the top file prefix')
parser.add_argument("-m","--inmem",default=None,help='Trajectory')


# collect the arguments for use in runnign the program
args = parser.parse_args()
print(args.topfile)
if args.trajfile: 
    u_water = mda.Universe(args.topfile,args.trajfile)
else:
    u_water = mda.Universe(args.topfile)

if args.verbose:
    print(f'Loaded a universe: {u_water}')

# defaults to tip4p/2005 parameters. can be overidden with command line argument -q
qM,qH = args.charges

if mda_ver < '2.0.0' and args.ncores > 1:
    print('Warning: version too low for mpi')
    n_cores = 1
else:
    n_cores = args.ncores
    
    
print(qM,qH)

if args.fname:
    fname = args.fname
else:
    fname = 'dummys_added'



# calculate the vectors and store them as a .pkl just in case
w_dummy = Dummys(u_water,type_O=args.type_O,type_H=args.type_H,m_type=args.type_M)
if __name__ == "__main__":
    if not args.overwrite:
        try:
            # try to reload the vectors and initialise the dummy universe
            w_dummy.create_dummy_universe(reload=args.Mpickle)

        except:
            # if cannot reload calculate again
                m_vectors = w_dummy.find_dummy_positions(positioner = 'TIP4P', params = {'m_dist':args.Mdist}, 
                    re_write = args.overwrite, pickle = args.Mpickle, Ncores=n_cores,verbose=args.verbose)
                w_dummy.create_dummy_universe()
    else:
        m_vectors = w_dummy.find_dummy_positions(positioner = 'TIP4P', params = {'m_dist':args.Mdist}, 
                    re_write = args.overwrite, pickle = args.Mpickle, Ncores=n_cores,verbose=args.verbose)
        w_dummy.create_dummy_universe()

    print(w_dummy._M_universe)

    w_dummy.write_Ms(fname, top_format='.pqr', traj_format='.dcd')

    print(w_dummy.merged.dimensions)