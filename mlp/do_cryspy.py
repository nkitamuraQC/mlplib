#### 組成,化学式から3D構造決定

import os
from pymatgen.core import Structure
import periodictable
import numpy as np
import pathlib
from becmlp.mlp.job_cryspy import script_qe, script_ase


def get_section(myclass):
    section_base = {
        "[basic]": [
            "algo",
            "tot_struc",
            "nstage",
            "njob",
            "jobcmd",
            "jobfile",
        ],
        "[structure]": [
            "struc_mode",
            "natot",
            "atype",
            "nat",
        ],
    }
    if myclass.qe_ctrl is not None:
        section_add = {
            "[QE]": [
                "kppvol",
                "qe_infile",
                "qe_outfile",
                "pv_term",
            ],
        }
        section_base.update(section_add)
    elif myclass.calc_code == "ASE":
        section_add = {
            "[ASE]": [
                "ase_python",
            ],
        }
        section_base.update(section_add)
    if myclass.algo == "BO":
        section_add = {
            "[BO]": [
                "nselect_bo",
                "score",
                "dscrpt",
            ]
        }
    elif myclass.algo == "EA":
        section_add = {
            "[EA]": [
                "n_pop",
                "n_crsov",
                "n_perm",
                "n_strain",
                "n_rand",
                "n_elite",
                "slct_func",
            ],
        }
    elif myclass.algo == "RS":
        section_add = {}
    elif myclass.algo == "LAQA":
        section_add = {
            "[LAQA]": ["nselect_laqa"],
        }
    section_base.update(section_add)
    return section_base


class DoCryspy:
    def __init__(self, qe_ctrl=None):
        self.qe_ctrl = qe_ctrl
        self.prefix = None
        if qe_ctrl is not None:
            self.prefix = self.qe_ctrl.prefix
        self.algo = "RS"
        self.calc_code = "ASE"
        self.tot_struc = 100
        self.nstage = 20
        self.njob = 1
        self.jobcmd = "bash"
        self.jobfile = "job_cryspy"

        self.struc_mode = "crystal"
        self.natot = 9
        self.atype = "Bi Ti F"
        self.nat = "1 1 7"
        self.ase_python = "ase_in.py"

        self.kppvol = "40"
        self.qe_infile = f"{self.prefix}.in"
        self.qe_outfile = f"{self.prefix}.out"
        self.pv_term = False

        self.nselect_bo = 5
        self.score = "EI"
        self.dscrpt = "FP"

        self.n_pop = 10
        self.n_crsov = 10
        self.n_perm = 10
        self.n_strain = 10
        self.n_rand = 10
        self.n_elite = 2
        self.slct_func = "TNM"

        self.nselect_laqa = 5
        self.mpi = 1
        
        pathlib.Path("calc_in").mkdir(exist_ok=True)
    

    def make_cryspy(self, txt=""):
        self.section = get_section(self)
        for key in self.section:
            txt += key + "\n"
            for val in self.section[key]:
                txt += f"{val} = {self[val]}\n"
            txt += "\n"
        txt += "[option]\n"
        return txt

    def write_mlp(self):
        from becmlp.mlp.script import script
        fname = f"./calc_in/{self.ase_python}_1"
        wf = open(fname, "w")
        wf.write(script)
        wf.close()

        wf = open(f"./calc_in/{self.jobfile}", "w")
        wf.write(script_ase)
        wf.close()
        return

    def write_cryspy(self, inp):
        wf = open("cryspy.in", "w")
        wf.write(inp)
        wf.close()
        return
    
    def write_job_cryspy(self):
        wf = open(f"./calc_in/{self.jobfile}", "w")
        wf.write(script_qe.format(mpi=self.mpi,prefix=self.prefix))
        wf.close()
        return

    
    def write_qe(self):
        os.chdir("./calc_in")
        self.qe_ctrl.calculation = "relax"
        inp = self.qe_ctrl.make_input_for_cryspy()
        self.qe_ctrl.write_input4cryspy(inp, 1)

        self.qe_ctrl.calculation = "vc-relax"
        inp = self.qe_ctrl.make_input_for_cryspy()
        self.qe_ctrl.write_input4cryspy(inp, 2)

        os.chdir("..")
        return

    def exec_qe(self):
        inp = self.make_cryspy()
        self.write_cryspy(inp)
        self.write_job_cryspy()
        self.write_qe()
        os.system("cryspy")
        return

    def exec_mlp(self):
        inp = self.make_cryspy()
        self.write_cryspy(inp)
        # self.write_job_cryspy()
        self.write_mlp()
        os.system("cryspy")
        return

    def __len__(self):
        return len(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

    def __str__(self):
        return str(self.__dict__)

    def __iter__(self):
        return self.__dict__.iteritems()

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value
