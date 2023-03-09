This package adds virutual massless atoms. The main use for this is to add them in for calculations of the charge density from LAMMPS trajectories, which do not include the position of the water virtual site.

Currently aimed at tip4p water models, really only with TIP4P/2005 parameters as well, flexibility for other 4 body water models is getting there.


The main script for adding the dummy atoms is `tip4p_dummy.py`. Once installing with `pip install .`, this can be accessed from your path. Options for the file can be viewed with the `-h` flag.
