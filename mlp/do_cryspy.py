#### 組成,化学式から3D構造決定

import os
from pymatgen.core import Structure
import periodictable
import numpy as np
import pathlib
from mlp.job_cryspy import script_qe, script_mlp


def get_section(myclass):
    """
    DoCryspyクラスのインスタンスから入力ファイルのセクション情報を生成。
    Args:
        myclass (DoCryspy): 設定情報を持つインスタンス
    Returns:
        dict: セクションごとのパラメータリスト
    """
    section_base = {
        "[basic]": [
            "algo",
            "tot_struc",
            "calc_code",
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
    """
    CrySPY用の入力ファイル生成・ジョブ管理を行うクラス。
    QE/ASE/MLPの各種設定・ファイル出力・実行制御を担当。
    """
    def __init__(self, qe_ctrl=None):
        """
        DoCryspyインスタンスの初期化。
        Args:
            qe_ctrl: Quantum ESPRESSO制御用オブジェクト（任意）
        """
        self.qe_ctrl = qe_ctrl
        self.prefix = None
        if qe_ctrl is not None:
            self.prefix = self.qe_ctrl.prefix
        self.algo = "RS"
        self.calc_code = "ASE"
        self.tot_struc = 5
        self.nstage = 1
        self.njob = 5
        self.jobcmd = "zsh"
        self.jobfile = "job_cryspy"

        self.struc_mode = "crystal"
        self.natot = 2
        self.atype = "Fe Se"
        self.nat = "1 1"
        self.ase_python = "ase_in.py"
        self.qe_infile = "calculation.in"
        self.qe_outfile = "calculation.out"
        self.pv_term = False

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
        self.max_job = 10
    

    def make_cryspy(self, txt=""):
        """
        CrySPY用の入力ファイル内容を生成。
        Args:
            txt (str): 既存テキスト（省略可）
        Returns:
            str: 入力ファイル内容
        """
        self.section = get_section(self)
        for key in self.section:
            txt += key + "\n"
            for val in self.section[key]:
                txt += f"{val} = {self[val]}\n"
            txt += "\n"
        txt += "[option]\n"
        return txt

    def write_mlp(self):
        """
        MLP用のPythonスクリプト・ジョブファイルを出力。
        Returns:
            None
        """
        from mlplib.mlp.script import script
        for i in range(1, self.max_job+1):
            if not os.path.exists(f"./calc_in/{self.ase_python}_{i}"):
                break
        fname = f"./calc_in/{self.ase_python}_{i}"
        wf = open(fname, "w")
        wf.write(script)
        wf.close()

        wf = open(f"./calc_in/{self.jobfile}", "w")
        wf.write(script_mlp)
        wf.close()
        return

    def write_cryspy(self, inp):
        """
        CrySPY用の入力ファイル（cryspy.in）を出力。
        Args:
            inp (str): 入力内容
        Returns:
            None
        """
        wf = open("cryspy.in", "w")
        wf.write(inp)
        wf.close()
        return
    
    def write_job_cryspy(self):
        """
        Quantum ESPRESSO用ジョブファイルを出力。
        Returns:
            None
        """
        wf = open(f"./calc_in/{self.jobfile}", "w")
        wf.write(script_qe.format(mpi=self.mpi,prefix=self.prefix))
        wf.close()
        return

    
    def write_qe(self):
        """
        Quantum ESPRESSO入力ファイル（relax/vc-relax）を出力。
        Returns:
            None
        """
        self.qe_ctrl.calculation = "relax"
        inp = self.qe_ctrl.make_input_for_cryspy()
        self.qe_ctrl.write_input4cryspy(inp, 1)

        self.qe_ctrl.calculation = "vc-relax"
        inp = self.qe_ctrl.make_input_for_cryspy()
        self.qe_ctrl.write_input4cryspy(inp, 2)
        return

    def exec_qe(self):
        """
        QE計算のCrySPYワークフローを実行。
        Returns:
            None
        """
        pathlib.Path("calc_in").mkdir(exist_ok=True)
        inp = self.make_cryspy()
        self.write_cryspy(inp)
        self.write_job_cryspy()
        self.write_qe()
        os.system("cryspy")
        return

    def exec_mlp(self):
        """
        MLP計算のCrySPYワークフローを実行。
        Returns:
            None
        """
        pathlib.Path("calc_in").mkdir(exist_ok=True)
        inp = self.make_cryspy()
        self.write_cryspy(inp)
        # self.write_job_cryspy()
        self.write_mlp()
        os.system("cryspy")
        return
    
    def rm_lock(self):
        """
        CrySPY関連のロック・ログ・データを削除。
        Returns:
            None
        """
        if os.path.exists("lock_cryspy"):
            os.remove("lock_cryspy")
        if os.path.exists("log_cryspy"):
            os.remove("log_cryspy")
        if os.path.exists("cryspy.stat"):
            os.remove("cryspy.stat")
        os.system("rm -r calc_in")
        os.system("rm -r data")
        return

    def __len__(self):
        """
        インスタンスの属性数を返す。
        Returns:
            int: 属性数
        """
        return len(self.__dict__)

    def __repr__(self):
        """
        インスタンスの属性情報を文字列で返す。
        Returns:
            str: 属性情報
        """
        return str(self.__dict__)

    def __str__(self):
        """
        インスタンスの属性情報を文字列で返す。
        Returns:
            str: 属性情報
        """
        return str(self.__dict__)

    def __iter__(self):
        """
        属性辞書のイテレータを返す。
        Returns:
            イテレータ
        """
        return self.__dict__.iteritems()

    def __getitem__(self, key):
        """
        属性辞書からkeyで値を取得。
        Args:
            key: キー
        Returns:
            値
        """
        return self.__dict__[key]

    def __setitem__(self, key, value):
        """
        属性辞書にkeyで値を設定。
        Args:
            key: キー
            value: 設定値
        Returns:
            None
        """
        self.__dict__[key] = value
