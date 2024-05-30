# Basics
import os
import numpy as np
from glob import glob

# MDAnalysis
import MDAnalysis

# ASE
from ase import Atoms
from ase import io

# Trajectory source
filename_dcd = 'dyna*.dcd'
split_dcd = ['.', 1] # split a each '.' and take element [1] in filename_dcd
file_psf = 'dyna.psf'

# Data directory ('.' means current directory where .py file is executed)
datadir = '.'

# Detect all dcd files
file_list_dcd = glob(os.path.join(datadir, filename_dcd))

# Get the dcd index numbers
# e.g. if 'glob' returns a list of 
# '["datadir/dyna.0.dcd", "datadir/dyna.1.dcd", "datadir/dyna.1.dcd"]'
# the 'idcds' will become [0, 1, 2]
idcds = np.array([
    int(file_dcd.split('/')[-1].split(split_dcd[0])[split_dcd[1]])
    for file_dcd in file_list_dcd])

# Sort dcd files
# In case the dcd file list is not in correct order 
# (e.g. idcds = [1, 10, 11, ..., 19, 2, 21, ...])
sort_dcd = np.argsort(idcds)
file_list_dcd = np.array(file_list_dcd)[sort_dcd]
idcds = idcds[sort_dcd]

# Prepare PSF file path
file_psf = os.path.join(datadir, file_psf)

# Initialize ASE trajectory
#atoms_traj = io.Trajectory('dyna.traj', 'w')

# Loop over dcd files
for ii, idcd in enumerate(idcds):
    
    # Open dcd file
    file_dcd = file_list_dcd[ii]
    dcd = MDAnalysis.Universe(file_psf, file_dcd)
    
    # Get trajectory parameter
    # Nframes: Number of frames in dcd trajectory file
    # Nskip: equivalent to 'nsavc' in the CHARMM input file
    # dt: equivalent to 'timestep' in the CHARMM input file
    Nframes = len(dcd.trajectory)
    Nskip = int(dcd.trajectory.skip_timestep)
    dt = np.round(
        float(dcd.trajectory._ts_kwargs['dt']), decimals=8)/Nskip
    
    # Get atom types (When passing to ASE, ASE cannot always recognize 
    # the atom element correctly, e.g. atom type 'CA' of alpha carbon
    # will be read as Calcium atom, for a fix see below)
    atoms = np.array([ai for ai in dcd._topology.names.values])
    Natoms = len(atoms)
    
    # Get atom masses
    masses = dcd._topology.masses.values
    
    # Get atom charges
    charges = dcd._topology.charges.values
    
    # Initialize an ASE atoms object with zero positions
    # (Just in preparation that, later on, the positions can be changed and
    # a copy can be added to 'atoms_list'. This is faster than initializing an
    # ASE Atoms object new for every frame)
    #ase_atoms = Atoms(atoms, positions=np.zeros([Natoms, 3], dtype=float))
    
    # If atoms are wrongly recognized, you can correct them manually by:
    # Here, atoms recognized as Calcium (Z_i=20) are changed to carbon (Z_i=6)
    # Z = ase_atoms.get_atomic_numbers()
    # Z[Z==20] = 6
    # ase_atoms = Atoms(Z, positions=np.zeros([Natoms, 3], dtype=float))
    
    # Iterate over frames
    # Eventually, if you only want to add every nth step to the ASE trajectory
    # file, you can change from 'enumerate(dcd.trajectory)' to
    # 'enumerate(dcd.trajectory[::nth])'
    for jj, frame in enumerate(dcd.trajectory):
        
        # Get atom positions in frame
        positions = frame._pos
        
        # (Eventually) Get cell information
        #cell = frame._unitcell
        
        # Set positions to ASE Atoms object (and eventually cell information)
        ase_atoms.set_positions(positions)
        #ase_atoms.set_cell(cell)
        
        # Add a copy to the ASE trajectory file
        atoms_traj.write(ase_atoms.copy())


    
    
