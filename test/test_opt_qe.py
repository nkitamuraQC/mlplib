from pylattice.pylat.qe_controller import QEController
import pytest
import os

def test_run_qe_opt():
    cifname = "../input/FeSe_mp-20120_primitive.cif"
    pseudo_dict = {
        "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF",
    }
    qe = QEController(cifname, pseudo_dict)
    qe.calculation = "vc-relax"
    qe.pseudo_dir = "../input"
    qe.press = 0.0

    # qe.tefield = True
    # qe.dipfield = True
    # qe.lda_plus_u = True
    # qe.nspin = 2

    inp = qe.make_input()
    qe.write_input(inp)
    qe.exec()

    return
    