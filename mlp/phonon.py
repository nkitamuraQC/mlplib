from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import read_crystal_structure
import numpy as np
from ase import Atoms
from mlplib.mlp.make_calc import get_mlp_calculator, EquiformerWithStress
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

class PhononCalculator:
    def __init__(self, structure_file, force_constants=None):
        """
        Initialize the PhononCalculator with a crystal structure and supercell matrix.

        :param structure_file: Path to the crystal structure file (e.g., POSCAR).
        :param supercell_matrix: Supercell matrix for phonon calculations.
        :param force_constants: Optional force constants for the phonon calculation.
        """
        self.structure = read_crystal_structure(structure_file)
        self.supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        self.force_constants = force_constants
        self.distance = 0.03
        
        self.symbols = None
        self.positions = None
        self.cell_param = None
        self.ase_atoms = None
        self.forces_save = None
        self.phonon = None
        self.ds = None
        self.phonopy_save = "phonopy_params.yaml"

    def arrange_forces(self, forces):
        cumsum = np.cumsum(self.len_supercell)
        for i in range(len(self.len_supercell)):
            if i == 0:
                self.forces_save = forces[:cumsum[0],:]
            else:
                self.forces_save = np.concatenate((self.forces_save, forces[cumsum[i-1]:cumsum[i],:]), axis=0)
        return
    
    def load(self, phonopy_save):
        """
        Load phonon data from a saved file.

        :param phonopy_save: Path to the saved phonon data file.
        """
        self.phonon = Phonopy.load(phonopy_save)
        return

    def calculate_phonons(self):
        """
        Calculate phonons using the provided crystal structure and supercell matrix.
        """
        self.phonon = Phonopy(self.structure, self.supercell_matrix)
        if self.force_constants is not None:
            self.phonon.set_force_constants(self.force_constants)
        self.forces_save = []
        self.phonon.generate_displacements(distance=self.distance)
        supercells = self.phonon.supercells_with_displacements
        self.cell_param = supercells[0].get_cell() * np.array(self.supercell_matrix)
        for supercell in supercells:
            self.symbols = supercell.symbols
            self.positions = supercell.get_positions()
            self.ase_atoms = Atoms(symbols=self.symbols, positions=self.positions, cell=self.cell_param, pbc=True)
            ocp = get_mlp_calculator()
            self.ase_atoms.calc = EquiformerWithStress(ocp)
            e = self.ase_atoms.get_potential_energy()
            forces = self.ase_atoms.get_forces()
            stress = self.ase_atoms.get_stress()
            print("Energy:", e)
            print("Forces:", forces)
            print("Stress:", stress)
            self.forces_save.append(forces)

        # self.arrange_forces(forces)
        self.phonon.forces = self.forces_save
        self.phonon.dataset = self.phonon.displacement_dataset
        self.phonon.produce_force_constants()
        self.phonon.save(filename=self.phonopy_save)
        return

    def get_phonon_band_structure(self, path:list[list[float]], labels:list[str], mesh:list[int]):
        """
        Get the phonon band structure along a specified path in reciprocal space.

        :param qpoints: List of q-points along the path.
        :param path_connections: List of connections between q-points.
        :return: Phonon band structure data.
        """
        qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
        self.phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
        self.phonon.plot_band_structure().show()

        # To plot DOS next to band structure
        self.phonon.run_mesh(mesh)
        self.phonon.run_total_dos()
        self.phonon.plot_band_structure_and_dos().show()

        # To plot PDOS next to band structure
        self.phonon.run_mesh(mesh, with_eigenvectors=True, is_mesh_symmetry=False)
        self.phonon.run_projected_dos()
        self.phonon.plot_band_structure_and_dos(pdos_indices=[[0], [1]]).show()
        return 
