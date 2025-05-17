from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import read_crystal_structure
import numpy as np
from ase import Atoms
from mlplib.mlp.make_calc import get_mlp_calculator, EquiformerWithStress
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from ase.io import read, write
from phono3py import Phono3py

class PhononCalculator:
    def __init__(self, structure_file, force_constants=None):
        """
        Initialize the PhononCalculator with a crystal structure and supercell matrix.

        :param structure_file: Path to the crystal structure file (e.g., POSCAR).
        :param supercell_matrix: Supercell matrix for phonon calculations.
        :param force_constants: Optional force constants for the phonon calculation.
        """
        atoms = read(structure_file)
        write("POSCAR", atoms, format="vasp")
        self.structure = read_crystal_structure("POSCAR")
        print("Structure:", self.structure)
        self.supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        self.force_constants = force_constants
        self.distance = 0.03
        self.use_ph3 = False
        
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

    def load_ph3(self, phonopy_save):
        from phono3py.interface.phono3py_yaml import Phono3pyYaml
        ph3yml = Phono3pyYaml()
        ph3yml.read("phono3py_disp.yaml")
        disp_dataset = ph3yml.dataset
        ph3 = Phono3py(unitcell, supercell_matrix=ph3yml.supercell_matrix, primitive_matrix=ph3yml.primitive_matrix)
        forces = np.loadtxt("FORCES_FC3").reshape(-1, len(ph3.supercell), 3)
        ph3.dataset = disp_dataset
        ph3.forces = forces
        self.phonon = ph3
        return

    def calculate_phonons(self):
        """
        Calculate phonons using the provided crystal structure and supercell matrix.
        """
        if not use_ph3:
            self.phonon = Phonopy(self.structure[0], self.supercell_matrix)
        else:
            self.phonon = Phono3py(
                self.structure[0], 
                supercell_matrix=self.supercell_matrix, 
                primitive_matrix='auto'
            )
        if self.force_constants is not None:
            self.phonon.set_force_constants(self.force_constants)
        self.forces_save = []
        self.phonon.generate_displacements(distance=self.distance)
        supercells = self.phonon.supercells_with_displacements
        self.cell_param = supercells[0].get_cell() * np.array(self.supercell_matrix)
        ocp = get_mlp_calculator()
        eqcalc = EquiformerWithStress(ocp)
        for supercell in supercells:
            self.symbols = supercell.symbols
            self.positions = supercell.get_positions()
            self.ase_atoms = Atoms(symbols=self.symbols, positions=self.positions, cell=self.cell_param, pbc=True)

            self.ase_atoms.calc = eqcalc
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
        if not use_ph3:
            self.phonon.produce_force_constants()
        else:
            self.phonon.produce_fc3()
        self.phonon.save(filename=self.phonopy_save, settings={'force_constants': True})
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
