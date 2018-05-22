#setup the gpaw calculation
from ase.build import molecule
struc = molecule('C6H6')
from ase.io import read, write
from gpaw import GPAW
from ase.units import Bohr
from ase.units import Hartree
#from ase.optimize import FIRE #Quasi Newton + friction
from ase.optimize.bfgslinesearch import BFGSLineSearch #Quasi Newton

import os.path
traj_file = 'benzene.traj'

if not os.path.isfile(traj_file):
    struc.set_cell([15,15,15])
    struc.set_pbc([0,0,0])
    struc.center()
    calc  = GPAW(xc='PBE', h=0.2, charge=0,
             spinpol=True, convergence={'energy': 0.001})

    struc.set_calculator(calc)
    #opt   = FIRE(struc, trajectory='benzene.traj', logfile='fire.log')
    dyn = BFGSLineSearch(struc, trajectory=traj_file,
                     restart='bfgs_ls.pckl', logfile='BFGSLinSearch.log')
    dyn.run(fmax=0.05)

# Did not find any documentation on how to directly continue with results
struc = read(traj_file)
calc  = GPAW(xc='PBE', h=0.2, charge=0,
             spinpol=True, convergence={'energy': 0.001})
struc.set_calculator(calc) # is it necessary?
struc.get_potential_energy()

vHtg            = struc.calc.get_electrostatic_potential()
rho4            = struc.calc.get_all_electron_density(gridrefinement=4)
rho_pseudo      = struc.calc.get_pseudo_density()
rho             = struc.calc.get_all_electron_density()

vHtg_hartree                    = vHtg / Hartree
rho4_per_bohr_cube              = rho4 * Bohr**3
rho_pseudo_per_bohr_cube        = rho_pseudo * Bohr**3
rho_per_bohr_cube               = rho * Bohr**3

write('benzene_vHtg.cube', struc, data=vHtg_hartree) 
write('benzene_rho_gridrefinement4.cube', struc, data=rho4 * Bohr**3)
write('benzene_rho.cube', struc, data=rho_per_bohr_cube)
write('benzene_rho_pseudo.cube', struc, data=rho_pseudo_per_bohr_cube)



