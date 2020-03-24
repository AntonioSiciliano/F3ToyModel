from __future__ import print_function
import fforces as ff
import fforces.Calculator

try:
    import cellconstructor as CC
    import cellconstructor.Phonons
except:
    raise ImportError("Error, we need CellConstructor installed!")


import ase
import ase.units
import numpy as np
import matplotlib.pyplot as plt

# Number of random structures
N_RANDOM = 10

def test_SnSe():
    
    # We load the dynamical matrices
    ff_dyn = CC.Phonons.Phonons("ffield_dynq", 3)

    dyn_harm = ff_dyn.Copy()
    dyn_harm.ForcePositiveDefinite()

    # Lets setup the calculator
    calc = ff.Calculator.ToyModelCalculator(ff_dyn, type_cal = "harmx")

    # Now we extract some randomly distributed structures,
    # and we compare the harmonic force with the one from the calculator
    supercell_dyn = ff_dyn.GenerateSupercellDyn(ff_dyn.GetSupercell())
    supercell_dyn_fp = dyn_harm.GenerateSupercellDyn(ff_dyn.GetSupercell())
    structures = supercell_dyn_fp.ExtractRandomStructures(N_RANDOM)

    for i in range(N_RANDOM):
        # Get the harmonic energy forces
        en_harm, f_harm = supercell_dyn.get_energy_forces(structures[i])
        en_harm *= ase.units.Ry
        f_harm *= ase.units.Ry

        # Get the force field energy and forces
        ase_atoms = structures[i].get_ase_atoms()
        ase_atoms.set_calculator(calc)

        energy = ase_atoms.get_total_energy()
        forces = ase_atoms.get_forces()

        # print()
        # print("Structure {}".format(i))
        # print("{:16.6f} {:16.6} eV".format(en_harm, energy))
        # print("Harmonic forces:")
        # print(f_harm)
        # print("ToyModel forces:")
        # print(forces)

        assert np.max (np.abs(en_harm / energy - 1)) < 1e-5
        assert np.max (np.abs((np.real(f_harm) / forces) - 1)) < 1e-5


    s0 = structures[0]
    energies = []
    forces = []
    x_pos = []
    step = 0.05

    for i in range(100):
        s0.coords[2,0] += step

        atoms = s0.get_ase_atoms()
        atoms.set_calculator(calc)
        energies.append(atoms.get_total_energy())
        forces.append(atoms.get_forces()[2,0])
        x_pos.append( s0.coords[2,0]  )
        

    plt.figure()
    plt.plot(x_pos, -np.gradient(energies) / (x_pos[1] - x_pos[0]), label="num  diff")
    plt.plot(x_pos, forces, label="anal diff")
    plt.legend()
    plt.tight_layout()

    
    plt.figure()
    plt.title("Energy")
    plt.plot(x_pos, energies)
    plt.legend()
    plt.tight_layout()
    plt.show()
        
if __name__ == "__main__":
    
    test_SnSe()
