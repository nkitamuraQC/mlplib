from fairchem.core import OCPCalculator
from ase.calculators.calculator import Calculator, all_changes
import numpy as np
from ase import Atoms

def get_mlp_calculator():
    calc = OCPCalculator(
        model_name="EquiformerV2-31M-S2EF-OC20-All+MD",
        local_cache="pretrained_models",
        seed=42,
    )
    return calc

def compute_virial_stress_from_ase(atoms: Atoms, forces: np.ndarray) -> np.ndarray:
    """
    ASEのAtomsオブジェクトとforcesからVirial stress tensorを計算する

    Parameters:
        atoms: ASE Atoms オブジェクト
        forces: (N, 3) numpy array of atomic forces [eV/Å]

    Returns:
        stress_tensor: (3, 3) numpy array of stress tensor [eV/Å^3]
    """
    positions = atoms.get_positions()  # shape: (N, 3)
    volume = atoms.get_volume()

    assert positions.shape == forces.shape, "positions and forces must have the same shape"

    stress_tensor = np.einsum('ni,nj->ij', atoms.positions, forces)

    stress_tensor = -stress_tensor / volume  # minus sign per virial definition
    return stress_tensor

class EquiformerWithStress(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, ocp_calculator, stress_func=None, **kwargs):
        super().__init__(**kwargs)
        self.ocp_calc = ocp_calculator
        self.stress_func = compute_virial_stress_from_ase  # 自前の stress 関数（例: virial stress）

    def calculate(self, atoms=None, properties=['energy', 'forces'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)

        # OCPCalculatorでエネルギーと力を取得
        self.ocp_calc.calculate(atoms, properties, system_changes)
        energy = self.ocp_calc.results['energy']
        forces = self.ocp_calc.results['forces']

        # Stressを自前の関数で計算（例: virial stress）
        if self.stress_func is not None:
            stress = self.stress_func(atoms, forces)
        else:
            stress = np.zeros((3, 3))  # ダミー（要注意）

        # 保存
        self.results['energy'] = energy
        self.results['forces'] = forces
        self.results['stress'] = stress