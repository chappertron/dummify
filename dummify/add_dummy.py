### TODO add name assignment to the other atoms in the system.
### TODO make this script work with LAMMPS data files too!

import MDAnalysis as mda
import pickle as pkl
import numpy as np
import numpy.linalg as la

import traceback


import multiprocessing
from multiprocessing import Pool
from functools import partial


class Dummys:
    
    def __init__(self,universe, type_O = 1, type_H = 2, m_type='M',load_into_memory = False,verbose=False):
        '''Takes in an mda.Universe object, not the raw files'''
        #self.ftraj = water_traj
        self.type_O = type_O 
        self.type_H = type_H
        
        self._verbose = verbose

        if type(self.type_H) == list or type(self.type_H) == tuple :
            self._select_hs = f"type {' '.join([str(i) for i in self.type_H])}"
        else:
            self._select_hs = f'type {self.type_H}'
        self._select_Os = f'type {self.type_O}'

        # load trajectory with MD analysis Universe
        self.u_water = universe

        # define dummy atom properties
        self.dummy_mass = 0
        self.dummy_charge = self.u_water.select_atoms(self._select_Os)[0].charge ##get the charge of the very first oxygen atom
        self.dummy_type = m_type

    def unit_dip(self, res):
        '''
            Takes in a residue corresponding to a water molecule and returns the unit vector corresponding to the dipole moment direction, relative to the oxygen atom

            res : a MDAnalysis residue object

            returns: a 3 component nd.array of the unit_dipole vector 

            TODO: make this a static method or an independant function??
        '''

        # select hydrogens and oxygen atoms within the molecule(residue)
        O = res.atoms.select_atoms(self._select_Os) # a better approach is to select all Os and then all Hs/every other
        Hs = res.atoms.select_atoms(self._select_hs)


        # check that there are only 2 hydrogens present
        if len(Hs) != 2:
            raise ValueError('There must only be two hydrogens in a water molecule!!!')
        # get the H positions relative to oxygen atom
        # O is indexed because to should be a list of just one atom!
        OH_poses = [h.position - O[0].position for h in Hs]
        
        #can simply take the sum because the OH vectors should lie in the same plane and their x components are opposite in the plane frame, so cancel, leaving only the y compoenent
        
        vect = np.sum(OH_poses,axis=0)

        return vect/la.norm(vect)

    def unit_dips(self):
        '''
            Takes in a residue corresponding to a water molecule and returns the unit vector corresponding to the dipole moment direction, relative to the oxygen atom

            res : a MDAnalysis residue object

            returns: a 3 component nd.array of the unit_dipole vector 

            TODO: Change implementation so this is called directly whilst writing files

        '''

        # select hydrogens and oxygen atoms within the molecule(residue)
        # a better approach is to select all Os and then all Hs/every other
        
        ag = self.u_water.atoms

        Os = ag.select_atoms(self._select_Os)# select all hydrogens!!!
        Hs = ag.select_atoms(self._select_hs)

        # check that there are only 2 hydrogens present
        if len(Hs) != 2*len(Os):
            raise ValueError(
                'There must only be two hydrogens in a water molecule!!!')
        # get the H positions relative to oxygen atom
        # O is indexed because to should be a list of just one atom!
        posO = Os.positions
        posH1 = Hs.positions[::2]
        posH2 = Hs.positions[1::2]
        #OH_poses = [h.position - O[0].position for h in Hs]
        HO1 = posH1-posO
        HO2 = posH2-posO

        #can simply take the sum because the OH vectors should lie in the same plane and their x components are opposite in the plane frame, so cancel, leaving only the y compoenent

        vect = HO1 + HO2
        # was an issue here, because vector of many atoms, nrom is calculated not per atom but of whole vect array 
        return vect/la.norm(vect,axis=1)[:,None] #weird indexing adds another axis

    def tip4p_M_poses(self, M_dist):
        '''
            Sets the virtual site position to be a residue a specified distance away from the oxygen atom, along the direction of the dipole moment vector

            res : MDAnalysis residue - the water molecule/residue to compute the position on
            M_dist : float - the distance along the dipole moment vector from the oxygen atom

            returns np.ndarray - the absolute position of the virtual atom.

        
        '''
        rO = self.u_water.atoms.select_atoms(self._select_Os).positions # slight inefficiency here, having to reselect oxygen positions

        unit_vectors = self.unit_dips()

        return rO + unit_vectors * M_dist

    def _find_dummy_positions_timeframe(self,frame_index, positioner = 'TIP4P', params = {'m_dist': 0.1546}, N_dummy = 1,):
        '''find dummy atom positions for a single timeframe'''
        
        t = frame_index
        # set the trajectory frame
        self.u_water.trajectory[t]
        w_residues = self.u_water.residues
        if positioner == 'TIP4P':
            # use default tip4p_M_pos
            m_dist = params['m_dist']
            return self.tip4p_M_poses(m_dist)
        elif type(positioner) == callable:
            # use this, needs to take in a residue as an input and some parameters
            return [positioner(res, **params) for res in w_residues]
        elif type(positioner) == np.ndarray:
            # just set positions to this!!!! kinda pointless tbh
            return positioner[t]
        else:
            raise TypeError(
                "positioner must either be 'TIP4P', a callable that accepts a residue or a np array ")

    def find_dummy_positions(self, positioner = 'TIP4P', params = {'m_dist':0.1546}, N_dummy = 1, re_write = False, pickle = 'm_vectors.pkl',Ncores=1,verbose=True ):
        """
        docstring

        N_dummy : int - the number of dummy atoms to add to each molecule. the positioner must give the right number of vectors for this
        """
        
        if not re_write:
            try:
                self.m_vectors
                raise ValueError('Dummy atoms already calculated, please set re_write = True to overwrite')
            except AttributeError:
                print('Calculating!')
        # set up the list of residues
        w_residues = self.u_water.residues
        self.m_vectors = np.empty((len(self.u_water.trajectory), len(w_residues), 3))

        frame_indexes  = np.arange(self.u_water.trajectory.n_frames)

        # calculates this for all residues Using partial so nothing needs to be inputted except t)
        run_per_frame = partial(self._find_dummy_positions_timeframe,
                                positioner=positioner, params=params, N_dummy=1)
        if Ncores >1: # Only works for v 2.0.0-dev0 of MDAnalysis or higher. In random tests yields good results, though seemingly parallel mode is single prescision 
            with Pool(Ncores) as worker_pool:
                result = worker_pool.map(run_per_frame,frame_indexes)
            self.m_vectors = np.array(result)
        else:
            for t,ts in enumerate(self.u_water.trajectory):
                if t % 100 == 0 and verbose: print(ts.time)
                self.m_vectors[t] = run_per_frame(t)

        # back up to a file
        if pickle:
            with open(pickle, 'wb') as stream:
                pkl.dump(self.m_vectors, stream)

        return self.m_vectors


    def create_dummy_universe(self, reload = None, M_dist = 0.1546):
        '''Creates empty dummy universe, doesn't need whole array of m_vectors
        
            TODO: create a dt attribute of this class so that the timestep doesn't get set to the weird dcd default.
            have a default of 1 fs instead
            TODO: add option for not computing on the fly for sadists
        
        '''
        # if reload:
        #     #reload specified pickle file of vectors
        #     with open(reload,'rb') as stream:
        #         self.m_vectors = pkl.load(stream)
        
        # 1 m per moleculethe number of residues # may need to account for different numbers for general usage

        # TODO will cause issues if: there are oxygens of the type in other molecules
        # TODO Can't add the non-water molecules back

        n_residues = len(self.u_water.select_atoms(self._select_Os)) ###   rather than u.residues.n_residues, do this for the number of oxygens
        
        print(n_residues)
        resindices = np.arange(0, n_residues)
        segindices = [0]*n_residues

        atom_resindex = resindices
        n_segments = 1
        # create a mda.Universe to merge with the input one and then the trajectory can be written
        self._M_universe = mda.Universe.empty(n_residues, n_residues=n_residues,
                                                atom_resindex=resindices,
                                                residue_segindex=segindices,
                                                trajectory=True)
        # adding various required attributes to the system so that merge occurs properly
        # need dummy atom properties defined in init
        self._M_universe.add_TopologyAttr('resid', list(range(1, n_residues+1)))
        self._M_universe.add_TopologyAttr('type',[self.dummy_type]*n_residues)
        self._M_universe.add_TopologyAttr('name',[self.dummy_type]*n_residues)
        self._M_universe.add_TopologyAttr('segid', ['SYSTEM'])
        self._M_universe.add_TopologyAttr('charge', [self.dummy_charge]*n_residues)
        self._M_universe.add_TopologyAttr('masses', [self.dummy_mass]*n_residues)
        
        # needs to be implemented to have a .data type, for some reason!!! 
        # however need to set bond types
        # self._M_universe.add_TopologyAttr('bonds',[[i,i] for i in range(n_residues)])
        # self._M_universe.add_TopologyAttr('angles',[[i,i,i] for i in range(n_residues)])
        # self._M_universe.add_TopologyAttr('dihedrals',[[i,i,i,i] for i in range(n_residues)])
        # self._M_universe.add_TopologyAttr('impropers',[[i,i,i,i] for i in range(n_residues)])
        self._M_universe.dimensions=self.u_water.dimensions    

        # self._M_universe.bonds.types()
        # self._M_universe.angles.types()
        # self._M_universe.dihedrals.types()
        # self._M_universe.impropers.types()

        #extra depending on the type

        # set positions of the dummy atoms to the first frame

        print(len(self.tip4p_M_poses(M_dist=M_dist)))
        print(self._M_universe.atoms.positions.shape)
        self._M_universe.atoms.positions = self.tip4p_M_poses(M_dist=M_dist) # calculate M_dist positions

        # merged system
        u_water_to_merge = self.u_water.atoms

    

        #M_universe.atoms
        self.merged = mda.Merge(u_water_to_merge,self._M_universe.atoms)

        self.merged.add_TopologyAttr('name')
        #take dimensions from original box
        self.merged.dimensions = self.u_water.dimensions

        # set charges of oxygen in merged system to 0
        # and set names to 'O
        self.merged.select_atoms(self._select_Os).charges=[0]*n_residues
        self.merged.select_atoms(self._select_Os).names=['O']*n_residues
        self.merged.select_atoms(self._select_hs).names=['H']*n_residues*2
        self.merged.select_atoms(f'type {self.dummy_type}').names=['M']*n_residues


    # to do find a decent topology format!!!!
    def write_Ms(self, prefix,top_format = '.pdbqt', traj_format = '.dcd',M_dist = 0.1546,verbose = True): 
        '''Write a trajectory to file that includes the dummy atoms.
            Currently need a topology format that specifies the masses and charges 
         '''

        n_steps = len(self.u_water.trajectory)
        dt = self.u_water.trajectory.dt # NB this is dt between frames not the simulation timestep. Might also be guessed to be 1 fs, need to check other trajectories


        self.u_water.trajectory[0]

        n_atoms = self.merged.atoms.n_atoms
        merged_dummys = self.merged.select_atoms('type M').atoms
        merged_other = self.merged.select_atoms('not (type M)').atoms
        
    
        with mda.Writer(prefix+traj_format,n_atoms,dt=dt) as W:
            ### get atoms from the merged universe
            for t in range(0,n_steps):
                # set orginal water trajectory step
                self.u_water.trajectory[t]
                # set M coordinates to be corresponding from vector
                merged_dummys.positions = self.tip4p_M_poses(M_dist=M_dist) #calc on the fly and dump... probs faster than dumping a pickle
                #set other positons
                merged_other.positions = self.u_water.atoms.positions
                
                if t%100==0 and verbose: print(dt*t,f'step {t}')
                W.write(self.merged)
        

        ## for topology case, wrap the coordinates around in the box, to prevent issues with coordinates being outside the box.
        #write topology first???
        box_vect = self.u_water.dimensions[:3] ### vector of the box coordinates
        try:
            with mda.Writer(prefix+top_format, n_atoms) as W:
                # calc on the fly and dump... probs faster than dumping a pickle
                # set positions of initial dummy atoms
                merged_dummys.positions = self.tip4p_M_poses(M_dist=M_dist) % box_vect
                #set other positons
                merged_other.positions = self.u_water.atoms.positions % box_vect ## wrap coordinates. Ok to do, because topology only definies connectivity
                W.write(self.merged)
        except Exception as e:
            print('Could not write dummed topology file!!!')
            traceback.print_exc(e)  # print the exception
            pass
    
if __name__ == "__main__":
    
    #water_Dummy = Dummys('example_trj/tip4p05.data','example_trj/nemd.dcd')


    # water_Dummy.find_dummy_positions(positioner = 'TIP4P', params = {'m_dist': 0.1546}, N_dummy = 1, re_write = False, pickle = 'm_vectors_para.pkl', Ncores=4)


    #loading the 3 pickles to check parallelisation
    # parrallel mode is 32 bit

    with open('m_vectors.pkl','rb') as stream:
        m_vect_1 = pkl.load(stream)
    
    with open('m_vectors_2.pkl','rb') as stream:
        m_vect_2 = pkl.load(stream)

    with open('m_vectors_para.pkl','rb') as stream:
        m_vect_p = pkl.load(stream)

    m_vects = [m_vect_1,m_vect_2,m_vect_p]

    print([type(i)for i in m_vects])
    print([i.shape for i in m_vects])
    rand_i = np.random.randint(low=0,high = len(m_vect_p))
    rand_j = np.random.randint(low=0,high = len(m_vect_p[0]))
    print(f'random time frame and residues {rand_i} {rand_j}')
    print([i[rand_i][rand_j] for i in m_vects])



    ### testing new implementation,calculating directions of dip vectors

    uni = mda.Universe('example_trj/tip4p05.data', 'example_trj/nemd.dcd')
    water_Dummy = Dummys(uni)
    # this interface doesn't really make sense, should be a static method or just a method with no positional args
    dip_units = water_Dummy.unit_dips()


    water_Dummy.find_dummy_positions(pickle='m_vectors_fast.pkl')

    print(dip_units)
    print(dip_units.shape)
