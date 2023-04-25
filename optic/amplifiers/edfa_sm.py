import numpy as np
import os

def signal_args(param_signal):
    return True

def pump_args(param_pump):
    # forward
    param_pump.forPump = getattr(
        param_pump,
        "forPump",
        {"pump_signal": np.array([100e-3]), "pump_lambda": np.array([980e-9])},
    )
    # backward
    param_pump.bckPump = getattr(
        param_pump,
        "bckPump",
        {"pump_signal": np.array([100e-3]), "pump_lambda": np.array([980e-9])},
    )  
    return param_pump

def giles_args(param_giles):
    param_giles.file = getattr(param_giles, "file", "")
    if not (os.path.exists(param_giles.file)):
            raise TypeError(f"{param_giles.file} file doesn't exist.")    
    param_giles.type = getattr(param_giles, "type", "coeff")
    param_giles.fileunit = getattr(param_giles, "fileunit", "nm")
    return param_giles

def edf_args(param_edf):
    # init giles parameters
    param_edf.param_giles = giles_args(param_edf.param_giles)  
    # init geometric and material parameters
    param_edf.a = getattr(param_edf, "a", 1.56e-6)
    param_edf.b = getattr(param_edf, "b", 1.56e-6)
    param_edf.length = getattr(param_edf, "lngth", 8)
    param_edf.rho = getattr(param_edf, "rho", 0.955e25)
    param_edf.na = getattr(param_edf, "na", 0.22)
    # init algorithm parameters
    param_edf.geometric  = getattr(param_edf, "geometric ", "LP01")
    param_edf.algorithm  = getattr(param_edf, "algorithm ", "Giles_spectrum")
    if param_edf.algorithm not in (
        "Giles_spatial",
        "Giles_spectrum",
        "Saleh",
        "Jopson",
        "Inhomogeneous",
    ):
        raise TypeError(
            "edfaSM.algo invalid argument - [Giles_spatial, Giles_spectrum, Saleh, Jopson, Inhomogeneous]."
        )
    param_edf.longitudinal_steps = getattr(param_edf, "longitudinal_steps", 100)
     # init physical paramters
    param_edf.tal = getattr(param_edf, "tal", 10e-3)
    # init signal and pump loss paramters
    param_edf.loss_signal  = getattr(param_edf, "loss_signal", 2.08 * 0.0001 * np.log10(10))
    param_edf.loss_pump  = getattr(param_edf, "loss_pump", 2.08 * 0.0001 * np.log10(10))
    return param_edf

def noise_args(param_noise):
    param_noise.noiseBand = getattr(param_noise, "noiseBand", 125e9)
    return param_noise

def solver_args(param_solver):
    param_solver.longSteps = getattr(param_solver, "longSteps", 100)
    param_solver.tol = getattr(param_solver, "tol", 2 / 100)
    param_solver.tolCtrl = getattr(param_solver, "tolCtrl", 0.5)  # dB
    return param_solver

def control_args(param_control):
    # type of control
    param_control.type = getattr(param_control, "type", "AGC")
    # Verify amplification type
    if param_control.type not in ("AGC", "APC", "none"):
        raise TypeError("param_control.type invalid argument - [AGC, APC, none].")
    # setpoint -> dB (for gain) and dBm (for power)
    param_control.value = getattr(
        param_control, "value", 20
    )  
    param_control.kp = getattr(param_control, "kp", 1e-2)
    param_control.ki = getattr(param_control, "ki", 1e-2)
    param_control.kd = getattr(param_control, "kd", 5e-2)

def edfa_args(param_edfa):     
    param_edfa.param_pump = pump_args(param_edfa.param_pump)
    param_edfa.param_edf = edf_args(param_edfa.param_edf)
    param_edfa.param_noise = noise_args(param_edfa.param_noise)
    param_edfa.param_solver = solver_args(param_edfa.param_solver)
    param_edfa.param_control = control_args(param_edfa.param_control)
    return param_edfa

def edfa_sm(param_edfa):
    param_edfa = edfa_args(param_edfa)
    return True

if __name__ == "__main__":
    from core import parameters
    import os.path as path

    param_edfa = parameters()
    param_edfa.param_pump = parameters()
    param_edfa.param_edf = parameters()
    param_edfa.param_edf.param_giles = parameters()
    param_edfa.param_noise = parameters()
    param_edfa.param_solver = parameters()
    param_edfa.param_control = parameters()

    param_edfa.param_edf.param_giles.file = 'giles_MP980.dat'
    param_edfa.param_edf.param_giles.file = path.join(path.abspath(path.join("../")), 'OptiCommPy', 'optic', 'ampParams', param_edfa.param_edf.param_giles.file)

    edfa_sm(param_edfa)