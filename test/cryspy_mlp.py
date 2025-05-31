# from pylattice.pylat.qe_md import QEMD
# from pylattice.pylat.qe_controller import QEController
import numpy as np
import pytest
import os
from mlp.do_cryspy import DoCryspy
from ase import io
from Workflows.workflow_class.gen_psudo_dict import create_upf_dict

def test_cryspy():
    ### https://qiita.com/TomokiYamashit3/items/a07a03a6e4a65537869e

    dc = DoCryspy()
    dc.natot = 3
    dc.atype = "Bi Ti F"
    dc.nat = "1 1 1"
    dc.tot_struc = 10
    dc.njob = 1
    dc.exec_mlp()
    return

# def test_cryspy_qe():
#     ### https://qiita.com/TomokiYamashit3/items/a07a03a6e4a65537869e
#     cif = "./input/FeSe_mp-20120_primitive.cif"
#     pseudo_dict_all = create_upf_dict(qe.pseudo_dir)
#     pseudo_dict = {
#         "Bi": pseudo_dict_all["Bi"],
#         "Ti": pseudo_dict_all["Ti"],
#         "F": pseudo_dict_all["F"],
#     }
#     qe = QEController(cif, pseudo_dict)
#     qe.pseudo_dir = "../../ONCV/paw"
#     # 

#     dc = DoCryspy(qe_ctrl=qe)
#     dc.exec_qe()
#     return

if __name__ == "__main__":
    test_cryspy()
    # test_cryspy_qe()
