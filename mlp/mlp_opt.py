from becmlp.mlp.make_calc import get_mlp_calculator, EquiformerWithStress, compute_virial_stress_from_ase
from ase.io import read, write
import numpy as np
from scipy.optimize import minimize
from ase.filters import UnitCellFilter
from ase.optimize import FIRE
from ase import Atoms
from ase.io import Trajectory

class MLPOpt:
    def __init__(
            self, 
            cifname: str, 
            nstep=10, 
            press=0, 
            conv_thr_force=1e-2, 
            conv_thr_stress=1e-5,
            conv_thr=1e-5,
            temp=100
        ):
        self.cifname = cifname
        self.nstep = nstep
        self.press = press / 160.21766208 # GPa->eV/ang^3
        self.press_axis = [[0, 0]]
        self.conv_thr_force = conv_thr_force
        self.conv_thr_stress = conv_thr_stress
        self.conv_thr = conv_thr
        self.atoms = None
        self.temp = temp
        self.bulk_modulus = 1.0 
        self.cell_dyn = True
        self.velocities = None
        self.amu_to_me = 1822.888486209
        self.kB_Hartree = 3.1668114e-6  
        self.bohr2ang = 0.529177210903
        self.use_my_stress = True
        self.output_cif = cifname.replace(".cif", "_opt.cif")
        self.traj = Trajectory('output.traj', 'w')

    def read_cif(self):
        self.atoms = read(self.cifname)
        self.natoms = len(self.atoms)
        return 
    
    def assin_calc(self):
        self.atoms.calc = get_mlp_calculator()
        return
    
    def write_cif(self, atoms):
        atoms.write(self.output_cif, format='cif')
        return

    def get_props(self):
        self.read_cif()
        self.assin_calc()
        energy = atoms.get_potential_energy()
        print(f'Potential energy: {energy:.6f} eV')

        # 原子ごとのフォース（N原子 x 3成分）
        forces = atoms.get_forces()
        print(f'Forces:\n{forces}')

        # 応力テンソル（3x3 行列）単位：eV/Å^3
        stress = atoms.get_stress(voigt=False)
        print(f'Stress tensor (3x3):\n{stress}')
        return energy, forces, stress
    
    def initialize_velocities(self):
        masses_me = self.atoms.get_masses() * self.amu_to_me
        velocities = []
        for m in masses_me:
            sigma = np.sqrt(self.kB_Hartree * self.temp / m)
            v = np.random.normal(0.0, sigma, size=3) * self.bohr2ang
            velocities.append(v)
        velocities = np.array(velocities)
        v_cm = np.sum(velocities.T * masses_me, axis=1) / np.sum(masses_me)
        velocities -= v_cm
        self.velocities = velocities
        return
    
    def update_md(self, timestep_atoms=0.1, timestep_lattice=0.1):
        if self.velocities is None:
            self.initialize_velocities()

        forces = self.atoms.get_forces()
        masses_me = self.atoms.get_masses() * self.amu_to_me
        acc = forces / masses_me[:, None]

        # 位置更新 (Velocity Verlet)
        positions = self.atoms.positions
        positions += self.velocities * timestep_atoms + 0.5 * acc * timestep_atoms **2
        self.atoms.set_positions(positions)

        # 新しい力と加速度（再計算）
        forces_new = self.atoms.get_forces()
        acc_new = forces_new / masses_me[:, None]

        # 速度更新
        self.velocities += 0.5 * (acc + acc_new) * timestep_atoms

        # Stress取得（+加圧補正）
        stress = compute_virial_stress_from_ase(self.atoms, forces_new)
        for paxis in self.press_axis:
            stress[paxis[0], paxis[1]] -= self.press

        # ラティス更新
        if self.cell_dyn:
            strain = timestep_lattice * stress / self.bulk_modulus
            deformation = np.identity(3) + strain
            new_cell = self.atoms.cell @ deformation
            self.atoms.set_cell(new_cell, scale_atoms=True)

        # モニタリング
        self.traj.write(self.atoms)
        f_norm = np.linalg.norm(forces_new)
        s_norm = np.linalg.norm(stress)
        v_norm = np.linalg.norm(self.velocities)
        print(f"Force norm: {f_norm:.3e}, Stress norm: {s_norm:.3e}, Velocity norm: {v_norm:.3e}")
        return f_norm, s_norm
    
    
    def run_md(self, timestep_atoms=0.1, timestep_lattice=0.1):
        self.read_cif()
        self.assin_calc()
        for i in range(self.nstep):
            self.update_md(timestep_atoms=timestep_atoms, timestep_lattice=timestep_lattice)
            if i % 10 == 0:
                print(f"Step {i}: Energy = {self.atoms.get_potential_energy()}")
        self.write_cif(self.atoms)
        return
    
    
    def run_opt(self, fmax=0.05, dt=0.1, maxstep=0.05):
        self.read_cif()
        ocp = get_mlp_calculator()
        self.atoms.calc = EquiformerWithStress(ocp)
        ucf = UnitCellFilter(self.atoms, hydrostatic_strain=True, scalar_pressure=self.press)

        # FIREオプティマイザのセットアップ
        opt = FIRE(ucf, dt=dt, maxstep=maxstep)  # maxstepは1ステップの最大変位（Å）
        
        # 最適化実行（収束条件なども指定可能）
        opt.run(fmax=fmax, steps=self.nstep)

        self.write_cif(self.atoms)
        return
        
