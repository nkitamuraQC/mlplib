import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

class ConvexHullCalculator:
    def __init__(self):
        self.data = []  # [(fraction_B, energy_per_atom), ...]

    def add_structure(self, a_count, b_count, total_energy):
        """
        構造を追加（原子数と全エネルギー）

        Parameters:
        a_count (int): 元素Aの数
        b_count (int): 元素Bの数
        total_energy (float): 構造の全エネルギー（eV単位）
        """
        total_atoms = a_count + b_count
        if total_atoms == 0:
            raise ValueError("原子数がゼロです")
        fraction_b = b_count / total_atoms
        energy_per_atom = total_energy / total_atoms
        self.data.append((fraction_b, energy_per_atom))

    def compute_hull(self):
        """
        凸包を計算し、プロットする
        """
        if len(self.data) < 3:
            raise ValueError("凸包を描くには最低3点必要です")
        points = np.array(self.data)
        self.hull = ConvexHull(points)
        return self.hull

    def plot_hull(self, show_all_points=True):
        """
        凸包とデータをプロット
        """
        points = np.array(self.data)
        plt.figure(figsize=(6, 5))
        if show_all_points:
            plt.plot(points[:, 0], points[:, 1], 'o', label="Structures")

        for simplex in self.hull.simplices:
            x = points[simplex, 0]
            y = points[simplex, 1]
            plt.plot(x, y, 'k-', linewidth=2)

        plt.xlabel("Fraction of B atoms")
        plt.ylabel("Energy per atom (eV)")
        plt.title("Convex Hull")
        plt.legend()
        plt.tight_layout()
        plt.show()


hull_calc = ConvexHullCalculator()

# サンプル構造を追加（構造ごとの A原子数, B原子数, エネルギー）
hull_calc.add_structure(a_count=1, b_count=1, total_energy=-5.2)   # AB
hull_calc.add_structure(a_count=2, b_count=1, total_energy=-7.7)   # A2B
hull_calc.add_structure(a_count=1, b_count=2, total_energy=-7.3)   # AB2
hull_calc.add_structure(a_count=3, b_count=1, total_energy=-9.8)   # A3B
hull_calc.add_structure(a_count=1, b_count=3, total_energy=-9.0)   # AB3

# 凸包の計算とプロット
hull_calc.compute_hull()
hull_calc.plot_hull()
