import numpy as np

import ase
import ase.units
import ase.calculators.calculator as calc

import cellconstructor as CC 
import cellconstructor.Structure


# Import the fortran libraries
import fforces.libs



class ToyModelCalculator(calc.Calculator):

    def __init__(self, harmonic_dyn, type_cal = "pdhxx", *args, **kwargs):
        super().__init__(self, *args, **kwargs)

        # Type cal defines the kind of the force field to be used
        self.type_cal = type_cal

        # Setup what properties the calculator can load
        self.implemented_properties = ["energy", "forces"]

        # Prepare the variable specific for this calculator
        # This must be a CellConstructor Phonons
        self.harmonic_dyn = harmonic_dyn

        # The parameters
        # The second order force constant enhancement factor
        self.p2 = 1

        # The third and higher order coefficients
        self.p3 = 0
        self.p4 = 0
        self.p5 = 0
        self.p6 = 0

        self.p3x = 0
        self.p4x = 0
        self.p4f = 0
        self.p4g = 0

        # No clue of what they are
        self.b = 0
        self.c = 0


        # We prepare the data for the faster submission inside the fortran arrays
        superdyn = self.harmonic_dyn.GenerateSupercellDyn(self.harmonic_dyn.GetSupercell())
        phi_sc = superdyn.dynmats[0]
        
        super_structure = self.harmonic_dyn.structure.generate_supercell(self.harmonic_dyn.GetSupercell())
        self.nat_sc = super_structure.N_atoms
        self.phi_sc_harmonic = np.zeros((3,3, self.nat_sc, self.nat_sc), order = "F", dtype = np.double)

        for i in range(self.nat_sc):
            for j in range(self.nat_sc):
                self.phi_sc_harmonic[:, :, i, j] = phi_sc[3*i:3*i+3, 3*j:3*j+3]
        
        self.at_sc = np.zeros((3,3), dtype = np.double, order = "F") 
        self.at_sc[:,:] = super_structure.unit_cell.T

        self.tau_sc = np.zeros((3, self.nat_sc), dtype = np.double, order = "F")
        self.tau_sc[:,:] = super_structure.coords.T


        
    def calculate(self, atoms=None, *args, **kwargs):
        super().calculate(self, atoms, *args, **kwargs)


        # Here we implement the calculation on the atoms object
        structure = CC.Structure.Structure()
        structure.generate_from_ase_atoms(self.atoms)

        assert structure.N_atoms == self.nat_sc, "Error, the structure do not match the harmonic dyn given."

        # Get the vector of the displacements
        u_disp = np.zeros((1, self.nat_sc, 3), order = "F", dtype = np.double)
        u_disp[0, :, :] = structure.coords - self.tau_sc.T 
        # Convert the displacements to bohr
        u_disp *= CC.Units.A_TO_BOHR

        # Get the ityp of the structure
        structure.set_masses(self.harmonic_dyn.masses)
        ityp = structure.get_ityp() + 1

        # Perform the calculation
        forces, v = f3.libs.get_forces_energies(self.phi_sc_harmonic, u_disp, self.type_cal, ityp, \
            self.at_sc, self.tau_sc, self.b, self.c, self.p2, self.p3, self.p4, self.p5, self.p6, \
                self.p4x, self.p3x, self.p4f, self.p4g)

        # The energy is in [Ry] and the force is in Ry/bohr
        # Convert them in [eV] and [eV/A]
        forces /= CC.Units.A_TO_BOHR 
        v *= ase.units.Ry
        forces *= ase.units.Ry

        self.results = {"energy": v[0], "forces": forces[0,:,:]}