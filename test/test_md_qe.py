from pylattice.pylat.qe_md import QEMD
from pylattice.pylat.qe_controller import QEController
import numpy as np
import pytest
import os

def test_qe():
    cifname = "../input/FeSe_mp-20120_primitive.cif"
    pseudo_dict = {
        "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF",
    }
    qe = QEController(cifname, pseudo_dict, supercell=None)
    qe.conv_thr = "1d-5"
    qe.pseudo_dir = "../input"
    qe.ecutwfc = 60.0
    qe.ecutrho = 240.0
    dL = np.zeros((3, 3))
    qemd = QEMD(qe, dL=dL, macro_step=1, temp=100, nkpoints=[8, 8, 8], micro_step=30, press=10)
    qemd.run()
    return