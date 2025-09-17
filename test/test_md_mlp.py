from mlp.mlp_opt import MLPOpt
from mlp.phonon import PhononCalculator
from mlp import *


def test_run_mlp_opt():
    cifname = "./input/FeSe_mp-20120_primitive.cif"
    mlp = MLPOpt(cifname)
    mlp.run_opt()

    run_mlp_opt(cifname)
    return

def test_run_mlp_md():
    cifname = "./input/FeSe_mp-20120_primitive.cif"
    mlp = MLPOpt(cifname)
    mlp.run_md()

    run_mlp_md(cifname)
    return

def test_run_mlp_phonon():
    cifname = "./input/FeSe_mp-20120_primitive.cif"
    ph = PhononCalculator(cifname)
    ph.calculate_phonons()
    return

if __name__ == "__main__":
    test_run_mlp_opt()
    test_run_mlp_md()
    test_run_mlp_phonon()
    print("All Passed")
    
