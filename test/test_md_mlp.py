from becmlp.mlp import run_mlp_md


def test_run_mlp_opt():
    cifname = "../input/FeSe_mp-20120_primitive.cif"
    run_mlp_md(cifname, nstep=150, press=0, temp=1000000000)
    return

if __name__ == "__main__":
    test_run_mlp_opt()