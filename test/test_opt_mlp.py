from becmlp.mlp import run_mlp_opt, run_mlp_opt_simple
import pytest

def test_run_mlp_simple():
    cifname = "../input/FeSe_mp-20120_primitive.cif"
    paxis = [[0, 0], [1, 1], [2, 2]]
    run_mlp_opt_simple(cifname, nstep=100, press=10, paxis=paxis, conv_thr_force=1e-3, conv_thr_stress=1e-4)
    return

def _test_run_mlp_opt():
    cifname = "../input/FeSe_mp-20120_primitive.cif"
    paxis = [[0, 0], [1, 1], [2, 2]]
    run_mlp_opt(cifname, nstep=100, press=10, paxis=paxis)
    return