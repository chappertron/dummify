#!/usr/bin/env python
# A script for adding the virtual charged sites for a tip4p/05 model



# if mda_ver < '2.0.0':
#     Warning('version too low for mpi')
#     print('Warning: version too low for mpi')

# for handling command line arguemnts 

import argparse
import time

def parser_setup()-> argparse.ArgumentParser: 
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
    parser.add_argument("--fname",default='addedMs',help='output file prefix defaults to the top file prefix')
    parser.add_argument("-m","--inmem",default=False,action="store_true",help='Load Trajectory Into memory')
    parser.add_argument("-r","--reloadM",default=None,help='Reload a previously computed vector of positions')
    parser.add_argument("-o","--trajectory_out",default='.dcd',help='Give a default file name for the trajectory dump, if not specified, tries to do same extension',)
    parser.add_argument("-p","--topfile_out",default='.pdbqt',help='Give a name for outputted top file')

    return parser



def main():
    # collect the arguments for use in runnign the program
    parser = parser_setup()
    args = parser.parse_args()
    print(args.topfile)
    
    import MDAnalysis as mda
    from dummify import Dummys

    mda_ver = mda.__version__.split('-')[0] # to remove -dev0 if a development version is being used
    print('MDAnalysis version: ' + mda_ver)
    
    if args.trajfile: 
        u_water = mda.Universe(args.topfile,args.trajfile,in_memory = args.inmem)
    else:
        u_water = mda.Universe(args.topfile,in_memory = args.inmem)

    if args.verbose:
        print(f'Start time: {time.strftime("%a, %d %b %Y %H:%M:%S")}')

    if args.verbose:
        print(f'Loaded a universe: {u_water} from {args.topfile}')
        if args.trajfile:
            print(f'with a Trajectory of length {u_water.trajectory.n_frames} from {args.trajfile}')

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

    w_dummy = Dummys(u_water,type_O=args.type_O,type_H=args.type_H,m_type=args.type_M)
    
    if False:#args.reloadM: # check if none or notes
        try:
            # try to reload the vectors and initialise the dummy universe
            w_dummy.create_dummy_universe(reload=args.reloadM)

        except:

            # couldn't reload, recalculating
            # if cannot reload, calculate again
                m_vectors = w_dummy.find_dummy_positions(positioner = 'TIP4P', params = {'m_dist':args.Mdist}, 
                    re_write = args.overwrite, pickle = args.Mpickle, Ncores=n_cores,verbose=args.verbose)
                w_dummy.create_dummy_universe()
    else:
        #trying to see if calculating dummys on the fly is tractable
        #m_vectors = w_dummy.find_dummy_positions(positioner = 'TIP4P', params = {'m_dist':args.Mdist}, 
        #             re_write = args.overwrite, pickle = args.Mpickle, Ncores=n_cores,verbose=args.verbose)
        w_dummy.create_dummy_universe()
   
    if args.verbose:
        print('Dummy universe created',w_dummy._M_universe)
    # writing topology
    
    w_dummy.write_Ms(fname, top_format=args.topfile_out, traj_format=args.trajectory_out) # writeing

    print(w_dummy.merged.dimensions)

    if args.verbose:
        print(f'Finish time: {time.strftime("%a, %d %b %Y %H:%M:%S")}')


    # calculate the vectors and store them as a .pkl just in case
if __name__ == "__main__":
    main()

    