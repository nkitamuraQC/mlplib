from becmlp.mlp.mlp_opt import MLPOpt

def run_mlp_opt_simple(cifname, nstep=10, press=0, conv_thr_force=1e-2, conv_thr_stress=1e-4, paxis=[[0, 0]]):
    """
    Run MLP optimization on a CIF file.

    Parameters:
    cifname (str): The name of the CIF file.
    nstep (int): Number of optimization steps.
    press (float): Pressure in GPa.
    conv_thr (float): Convergence threshold.

    Returns:
    None
    """
    mlp_opt = MLPOpt(cifname, nstep=nstep, press=press, conv_thr_force=conv_thr_force, conv_thr_stress=conv_thr_stress, temp=0)
    mlp_opt.press_axis = paxis
    mlp_opt.run_md()
    mlp_opt.write_cif(mlp_opt.atoms)
    return

def run_mlp_opt(cifname, nstep=10, press=0, conv_thr=1e-3, paxis=[[0, 0]]):
    """
    Run MLP optimization on a CIF file.

    Parameters:
    cifname (str): The name of the CIF file.
    nstep (int): Number of optimization steps.
    press (float): Pressure in GPa.
    conv_thr (float): Convergence threshold.

    Returns:
    None
    """
    mlp_opt = MLPOpt(cifname, nstep=nstep, press=press, conv_thr=conv_thr)
    mlp_opt.press_axis = paxis
    mlp_opt.run_opt()
    mlp_opt.write_cif(mlp_opt.atoms)
    return

def run_mlp_md(cifname, nstep=10, press=0, conv_thr=1e-5, temp=100):
    """
    Run MLP molecular dynamics on a CIF file.

    Parameters:
    cifname (str): The name of the CIF file.
    nstep (int): Number of MD steps.
    press (float): Pressure in GPa.
    conv_thr (float): Convergence threshold.
    temp (float): Temperature in K.

    Returns:
    None
    """
    mlp_opt = MLPOpt(cifname, nstep=nstep, press=press, conv_thr=conv_thr, temp=temp)
    mlp_opt.run_md()
    mlp_opt.write_cif(mlp_opt.atoms)
    return
