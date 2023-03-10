import os

import numpy as np
import numpy.matlib as npmat
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from scipy import interpolate
from numpy.fft import fft, ifft, fftfreq
from scipy.integrate import solve_ivp

import pandas as pd
from scipy.signal import find_peaks
from scipy.constants import lambda2nu
from scipy import signal
from scipy.constants import c, Planck

from scipy.constants import c, Planck
from scipy.special import jv, kv

from simple_pid import PID
import logging as logg
import copy

from optic.core import parameters

class signal():
    def __init__(self, frequency, power):
        if len(frequency)!=len(power):
            raise TypeError(
            "signal class invalid argument - invalid lengths."
        )
        self.frequency = frequency
        self.power = power

    def get_frequency_hz(self):
        return self.frequency
    
    def get_frequency_m(self):
        return c/self.frequency
    
    def get_power_w(self):
        return self.power
    
    def get_power_dBm(self):
        return 10*np.log10(1e3*self.power)

    def plot(self):
        plt.plot(1e9*self.get_frequency_m(), self.get_power_dBm())
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Power [dBm]')
        plt.grid(True)

class control():
    def __init__(self, parameter):
        # controller parameters
        self.kp = getattr(parameter, "kp", 1e-2)
        self.ki = getattr(parameter, "ki", 1e-2)
        self.kd = getattr(parameter, "kd", 5e-2)
        # setpoint and control type
        self.type = getattr(parameter, "type", "AGC")
        if self.type in {"AGC", "APC", "none"}:
            self.setpoint = getattr(parameter, "setpoint", 20)  # dB (for gain) and dBm (for power)
        else:
            raise TypeError("control.type invalid argument - [AGC, APC, none].")

class solver():
    def __init__(self, parameter):
        # solver error
        self.tol = getattr(parameter, "tol", 2 / 100)
        self.tolCtrl = getattr(parameter, "tolCtrl", 0.5)     
        # noise parameters
        self.noiseBand = getattr(parameter, "noiseBand", 125e9)

class giles():
    def __init__(self, parameters):
        self.file = getattr(parameters, "file", "")
        if not (os.path.exists(self.file)):
            raise TypeError(f"{self.file} file doesn't exist.")
        self.file_unit = getattr(parameters, "fileunit", "nm")

    def set_giles_frequency_from_file(self):
        file_data = np.loadtxt(self.file)
        # Verify file frequency unit
        if self.file_unit == "nm":
            self.frequency = file_data[:,0] * 1e-9
        elif self.file_unit == "m":
            self.frequency = file_data[:,0]
        elif self.file_unit == "Hz'":
            self.frequency = c/file_data[:,0]
        elif self.file_unit == "THz":
            self.frequency = c/file_data[:,0] * 1e-12
        else:
            raise TypeError('gile.fileunit invalid argument - [nm - m - Hz - THz].')
    
    def set_giles_coefficients_from_file(self):
        file_data = np.loadtxt(self.file)
        self.absorption_coefficient  = 0.1 * np.log(10) * file_data[:,1]
        self.gain_coefficient = 0.1 * np.log(10) * file_data[:,2]

    def set_giles_cross_section_from_file(self):
        file_data = np.loadtxt(self.file)
        self.absorption_cross = file_data[:, 1]
        self.emission_cross = file_data[:, 2]

    def plot(self):
        plt.plot(self.frequency, 1e-25*self.absorption_cross, label='Abs. cross-section')
        plt.plot(self.frequency, 1e-25*self.emission_cross, label='Em. cross-section')
        plt.ylabel(r'Cross-section [10$^{-25}$ m$^2$]')
        plt.grid(True)
        plt.legend()        

class edf(giles):
    def __init__(self,parameter):
        # init giles parameters
        super().__init__(parameter.giles)        
        # init geometric and material parameters
        self.a = getattr(parameter, "a", 1.56e-6)
        self.b = getattr(parameter, "b", 1.56e-6)
        self.length = getattr(parameter, "length", 8)
        self.rho = getattr(parameter, "rho", 0.955e25)
        self.na = getattr(parameter, "na", 0.22)
        # init algorithm parameters
        self.geometric = getattr(parameter, "geometric", "LP01")
        self.algorithm = getattr(parameter, "algorithm", "Giles_spectrum")
        self.longitudinal_steps = getattr(parameter, "longitudinal_steps", 100)   
        # init physical paramters
        self.tal = getattr(parameter, "tal", 10e-3)
        # init signal and pump loss paramters
        self.loss_signal = getattr(parameter, "loss_signal", 2.08 * 0.0001 * np.log10(10))
        self.loss_pump = getattr(parameter, "loss_pump", 2.08 * 0.0001 * np.log10(10))

        if self.algo not in (
        "Giles_spatial",
        "Giles_spectrum",
        "Saleh",
        "Jopson",
        "Inhomogeneous",
        ):
            raise TypeError(
                "edfaSM.algo invalid argument - [Giles_spatial, Giles_spectrum, Saleh, Jopson, Inhomogeneous]."
                )

    def get_radius(self):
        # Logitudinal step
        self.dr = self.a / self.longitudinal_steps
        self.r  = np.arange(0,self.a, self.dr)

    def get_modal_paramters(self)
        self.V = (2 * np.pi / param_edf.lbFl) * param_edfa.a * param_edfa.na
        # u and v calculation for LP01 and Bessel profiles
        self.u = ((1 + np.sqrt(2)) * self.V) / (1 + (4 + self.V ** 4) ** 0.25)
        self.v = np.sqrt(self.V ** 2 - self.u ** 2)

    def get_mode_radius(self, radius, V, v, u):
        if self.gmtc == "Bessel":
            w_gauss = radius * V / u * kv(1, v) / kv(0, v) * jv(0, u)
        elif self.gmtc == "Marcuse":
            w_gauss = radius * (0.650 + 1.619 / V ** 1.5 + 2.879 / V ** 6)
        elif self.gmtc == "Whitley":
            w_gauss = radius * (0.616 + 1.660 / V ** 1.5 + 0.987 / V ** 6)
        elif self.gmtc == "Desurvire":
            w_gauss = radius * (0.759 + 1.289 / V ** 1.5 + 1.041 / V ** 6)
        elif self.gmtc == "Myslinski":
            w_gauss = radius * (0.761 + 1.237 / V ** 1.5 + 1.429 / V ** 6)
        else:
            raise TypeError(
            "model invalid argument - [LP01 - Marcuse - Whitley - Desurvire - Myslinski - Bessel]."
            )


def power_meter(x):
    """
    Calculate the total power of x.

    Parameters
    ----------
    x : np.array
        Signal.

    Returns
    -------
    scalar
        Total signal power of x: P = sum(abs(x)**2).

    """
    return np.sum(np.mean(x * np.conj(x), axis=0).real)

def OSA(x, Fs, Fc=193.1e12):
    """
    Plot the optical spectrum of the signal in X and Y polarizations.

    Parameters
    ----------
    x : np.array
        Signal
    Fs : scalar
        Sampling frequency in Hz.
    Fc : scalar, optional
        Central optical frequency. The default is 193.1e12.

    Returns
    -------
    plot

    """
    freqs, ZX = get_spectrum(x[:,0], Fs, Fc)
    yMin = -70
    yMax = ZX.max() + 10
    _,ax = plt.subplots(1)
    ax.plot( 1e9*freqs, ZX, label="X Pol.")
    if (np.shape(x)[1] == 2):
        freqs, ZY = get_spectrum(x[:,1], Fs, Fc)
        ax.plot( 1e9*freqs, ZY, label="Y Pol.", alpha=0.5)
        yMax = np.array([ZX.max(), ZY.max()]).max() + 10
        ax.legend()
    ax.set_ylim([yMin, yMax])   
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Magnitude [dBm]")
    ax.minorticks_on()
    ax.grid(True)
    
    return ax

def get_spectrum(x, Fs, Fc, xunits = 'm', yunits = 'dBm', window=mlab.window_none, sides="twosided"):
    """
    Calculates the optical spectrum of the signal.

    Parameters
    ----------
    x : np.array
        Signal
    Fs : scalar
        Sampling frequency in Hz.
    Fc : scalar, optional
        Central optical frequency. The default is 193.1e12.

    Returns
    -------
    spectrum : np.array
        Signal's FFT
    frequency: np.array
        Frequency array @ Fc.

    """
    spectrum, frequency = mlab.magnitude_spectrum(
        x, Fs=Fs, window=window, sides=sides
    )
    frequency = c/(frequency + Fc) if (xunits=='m') else frequency + Fc
    spectrum = spectrum*np.conj(spectrum)
    if (yunits=='dBm'):
        spectrum = 10*np.log10(1e3*spectrum)

    return frequency, spectrum

def gilesSpectrum(z, P, properties):
    """
    Routine used to solve the EDFA rate and propagation equations, considering the spectral Giles algorithm.

    Parameters
    ----------
    P : np.array
        Signal power (signal + pump + ASE).
    z : scalar (float)
        Position - erbium doped fiber [0 - edf length].
    properties : object with constants and edfa parameters.

    Returns
    -------
    Eo : np.array
        Increment of the amplified optical signal.

    """
    n2_norm = getN2Pop(P, properties)
    xi_k   = n2_norm * properties.const3 - properties.const4
    tauASE = n2_norm * properties.const5
    return properties.uk * (P * xi_k + properties.ASE * tauASE)

def gilesSpatial(z, P, properties, param_edf):
    """
    Routine used to solve the EDFA rate and propagation equations, considering the spatial Giles algorithm.

    Parameters
    ----------
    P : np.array
        Signal power (signal + pump + ASE).
    z : scalar (float)
        Position - erbium doped fiber [0 - edf length].
    properties : object with constants and edfa parameters.
    param_edf  : object with constants and edf parameters.

    Returns
    -------
    Eo : np.array
        Increment of the amplified optical signal.

    """
    n2_norm = getN2Pop(P, properties)
    intOL = getOverlapInt(n2_norm, properties, param_edf)
    xi_k = intOL * (properties.absCoef + properties.gainCoef) / properties.gamma - (properties.absCoef + properties.lossS)
    tauASE = intOL * (properties.gainCoef / properties.gamma) * Planck * properties.freq * properties.noiseBand
    return properties.uk * (P * xi_k + properties.ASE * tauASE)

def getN2Pop(P, properties):
    """
    Determines the number of carriers at the metastable level, considering the spectral and spatial Giles algorithm.

    Parameters
    ----------
    P : np.array
        Signal power (signal + pump + ASE).
    properties : object with constants and edfa parameters.

    Returns
    -------
    norm2 : np.array
        Number of carriers at the metastable level.

    """
    if (properties.algo == "Giles_spectrum"):
        n2_normT1 = np.dot(P, properties.const1)
        n2_normT2 = np.dot(P, properties.const2) + 1
    elif (properties.algo == "Giles_spatial"):
        n2_normT1 = (properties.tal / Planck) * (properties.i_k @ np.transpose(P * properties.absCross / properties.freq))
        n2_normT2 = (properties.tal / Planck) * (properties.i_k @ np.transpose(P * (properties.absCross + properties.emiCross) / properties.freq)) + 1
    return n2_normT1 / n2_normT2

def getOverlapInt(n2_norm, properties, param_edf):
    """
    Determines the overlap integral between the field envelope and the doping profile.

    Parameters
    ----------
    n2_norm : np.array
        Number of carriers at the metastable level.
    properties : object 
        With constants and edfa parameters.
    param_edf  : object 
        With edf parameters.

    Returns
    -------
    overlapIntegral : np.array
        Overlap integral between the field envelope and the doping profile.

    """     
    dopPrf = npmat.repmat(2 * np.pi * param_edf.r * n2_norm, np.shape(properties.i_k)[1], 1) * param_edf.dr
    return np.trapz(np.transpose(properties.i_k) * dopPrf)



def updtCnst(param):
    xi = np.pi * param.b ** 2 * param.rho / param.tal
    param.const1 = (1 / (Planck * xi)) * (param.absCoef / param.freq)
    param.const2 = (1 / (Planck * xi)) * (param.absCoef + param.gainCoef) / param.freq
    param.const3 = param.absCoef + param.gainCoef
    param.const4 = param.absCoef + param.lossS
    param.const5 = param.gainCoef * Planck * param.freq * param.noiseBand
    return param

def edfaArgs(param_edfa):

    # pump paramters
    param_edfa.forPump = getattr(
        param_edfa,
        "forPump",
        {"pump_signal": np.array([100e-3]), "pump_lambda": np.array([980e-9])},
    )
    param_edfa.bckPump = getattr(
        param_edfa,
        "bckPump",
        {"pump_signal": np.array([100e-3]), "pump_lambda": np.array([980e-9])},
    )    


    return param_edfa

def edfaSM(Ei, Fs, Fc, param_edfa):
    ## Verify arguments
    param_edfa = edfaArgs(param_edfa)

    if (param_edfa.type == "AGC"):
        power_in = power_meter(Ei)
    
    ## Get pump signal frequency points
    freqPmpFor = c / param_edfa.forPump["pump_lambda"]
    freqPmpBck = c / param_edfa.bckPump["pump_lambda"]
    pumpPmpFor = param_edfa.forPump["pump_signal"]
    pumpPmpBck = param_edfa.bckPump["pump_signal"]
    if len(freqPmpFor) != len(pumpPmpFor):
        raise TypeError(
            "edfaSM.forPump invalid argument - number of signals (freq and pump) must be igual."
        )
    if len(freqPmpBck) != len(pumpPmpBck):
        raise TypeError(
            "edfaSM.bckPump invalid argument - number of signals (freq and pump) must be igual."
        )
    lenPmpFor = np.size(freqPmpFor)
    lenPmpBck = np.size(freqPmpBck)

    ## Load Giles file
    # Get EDF cross-sections and coeficients parameters from giles parameters file.
    # absCross, emiCross, absCoef, gainCoef = edfParams(param_edfa, lbFl, fileT)
    param_edf = edfParams(param_edfa)

    ## Format input signal
    # Create second pol, if not exists
    lenFqSg, isy = np.shape(Ei)
    if isy == 1:
        Ei = np.concatenate((Ei, np.zeros(np.shape(Ei))), axis=1)
        isy += 1
    # Get signal in frequency domain
    freqSgn = Fs * fftfreq(len(Ei)) + Fc
    lenFqSg = len(freqSgn)

    ## Create ASE signal components
    # Get optical band and specify frequency points for ASE calculation
    opticalBand = freqSgn.max() - freqSgn.min()
    freqASE = np.arange(-opticalBand / 2, opticalBand / 2, param_edfa.noiseBand) + Fc
    lenASE = np.size(freqASE)

    ## Define the frequency vector used in simulation.
    # SIGNALX + SIGNALY + FASEX + FASEY + FORPUMP
    param_edfa.freq = np.concatenate([freqSgn, freqSgn, freqASE, freqASE, freqPmpFor])
    param_edfa.ASE = np.concatenate(
        [np.zeros(isy * lenFqSg), np.ones(isy * lenASE), np.zeros(lenPmpFor)]
    )
    param_edfa.uk = np.ones(np.size(param_edfa.freq))
    param_edfa.absCoef = np.interp(c / param_edfa.freq, param_edf.lbFl, param_edf.absCoef)
    param_edfa.gainCoef = np.interp(c / param_edfa.freq, param_edf.lbFl, param_edf.gainCoef)
    # BCKPUMP + BASEX + BASEY
    freqAdd = np.concatenate([freqPmpBck, freqASE, freqASE])
    ASEAdd = np.concatenate([np.zeros(lenPmpBck), np.ones(isy * lenASE)])
    ukAdd = np.concatenate([np.ones(lenPmpBck), np.ones(isy * lenASE)])
    absCoefAdd = np.interp(c / freqAdd, param_edf.lbFl, param_edf.absCoef)
    gainCoefAdd = np.interp(c / freqAdd, param_edf.lbFl, param_edf.gainCoef)

    # Indexes of each signal class (signal, ase for, pump for, ase back, pump back)
    idxPS = np.arange(0, isy * lenFqSg)
    idxPAF = np.arange(idxPS[-1] + 1, idxPS[-1] + lenASE * isy + 1)
    idxPPF = np.arange(idxPAF[-1] + 1, idxPAF[-1] + lenPmpFor + 1)
    idxPPB = np.arange(idxPPF[-1] + 1, idxPPF[-1] + lenPmpBck + 1)
    idxPAB = np.arange(idxPPB[-1] + 1, idxPPB[-1] + lenASE * isy + 1)

    ## Create variables us
    # ed in edo solution method
    if param_edfa.algo == "Giles_spatial":
        # SIGNALX + SIGNALY + FASEX + FASEY + FORPUMP
        param_edfa.absCross = np.interp(c / param_edfa.freq, param_edf.lbFl, param_edf.absCross)
        param_edfa.emiCross = np.interp(c / param_edfa.freq, param_edf.lbFl, param_edf.emiCross)
        param_edfa.gamma = np.interp(c / param_edfa.freq, param_edf.lbFl, param_edf.gamma)
        param_edfa.i_k = interpolate.interp2d(
            param_edf.lbFl, param_edf.r, param_edf.i_k, kind="cubic"
        )(c / param_edfa.freq, param_edf.r)
        # BCKPUMP + BASEX + BASEY
        absCrossAdd = np.interp(c / freqAdd, param_edf.lbFl, param_edf.absCross)
        emiCrossAdd = np.interp(c / freqAdd, param_edf.lbFl, param_edf.emiCross)
        gammaAdd = np.interp(c / freqAdd, param_edf.lbFl, param_edf.gamma)
        i_kAdd = interpolate.interp2d(param_edf.lbFl, param_edf.r, param_edf.i_k, kind="cubic")(
            c / freqAdd, param_edf.r
        )
        # Eval string
        evalStr = "solve_ivp(gilesSpatial, zSpan, pInit, method='DOP853', rtol = 5e-4, atol = 5e-7, args=(param_edfa,param_edf))"
    elif (param_edfa.algo == 'Giles_spectrum'):        
        # Update some constants used in rate and propagation equations
        param_edfa = updtCnst(param_edfa)
        # Eval string
        evalStr = "solve_ivp(gilesSpectrum, zSpan, pInit, method='DOP853', rtol = 5e-4, atol = 5e-7, args=(param_edfa,))"

    # Signal power vector and EDF length vector
    EiFt = fft(Ei, axis=0)
    Psgl = np.reshape(np.abs(EiFt / lenFqSg) ** 2, (isy * lenFqSg), order="F")
    zSpan = np.array([0, param_edfa.lngth])

    # Init power conditions
    # SIGNALX + SIGNALY + FASEX + FASEY + FORPUMP
    pInit = np.concatenate([Psgl, np.zeros(isy * lenASE), pumpPmpFor])
    # BCKPUMP + BASEX + BASEY
    pInitAdd = np.concatenate([pumpPmpBck, np.zeros(isy * lenASE)])

    # Solution: 0 -> L without BCKPUMP + BASEX + BASEY
    sol = eval(evalStr)
    Pout = sol["y"]
    zSpan = np.flip(zSpan)

    # Update interpolated parameters
    # SIGNALX + SIGNALY + FASEX + FASEY + FORPUMP
    # BCKPUMP + BASEX + BASEY
    pInit = np.concatenate([pInit, pInitAdd])
    param_edfa.freq = np.concatenate([param_edfa.freq, freqAdd])
    param_edfa.ASE = np.concatenate([param_edfa.ASE, ASEAdd])
    param_edfa.uk = np.concatenate([param_edfa.uk, -ukAdd])
    param_edfa.absCoef = np.concatenate([param_edfa.absCoef, absCoefAdd])
    param_edfa.gainCoef = np.concatenate([param_edfa.gainCoef, gainCoefAdd])

    if (param_edfa.algo == 'Giles_spatial'):
        param_edfa.absCross = np.concatenate([param_edfa.absCross, absCrossAdd])
        param_edfa.emiCross = np.concatenate([param_edfa.emiCross, emiCrossAdd])
        param_edfa.gamma    = np.concatenate([param_edfa.gamma,       gammaAdd])
        param_edfa.i_k      = np.concatenate([param_edfa.i_k,           i_kAdd], axis = 1)
    elif (param_edfa.algo == 'Giles_spectrum'): 
        # Update somes constants used in rate and propagation equations
        param_edfa = updtCnst(param_edfa)

    # update pInit signal with Pout SIGNAL + FASE + FORPUMP values
    pInit[idxPS] = Pout[idxPS, -1]
    pInit[idxPAF] = Pout[idxPAF, -1]
    pInit[idxPPF] = Pout[idxPPF, -1]

    # Variables used in loop
    MaxTry = 15
    tryCtrlLoop = 0
    errorAutoCrtl = 1

    while (np.abs(errorAutoCrtl) > param_edfa.tolCtrl) and (tryCtrlLoop < MaxTry):
        tryLoop = 0
        errorCvg = 1
        ## main loop
        while (np.mean(np.abs(errorCvg)) > param_edfa.tol) and (tryLoop < MaxTry):
            # Solution: L -> 0
            sol = eval(evalStr)
            Pin = sol["y"]
            zSpan = np.flip(zSpan)
            pInit = copy.deepcopy(Pin[:, -1])

            # Reset SIGNAL + FASE + FORPUMP values
            pInit[idxPS] = Psgl
            pInit[idxPAF] = np.zeros(lenASE * isy)
            pInit[idxPPF] = pumpPmpFor

            # Solution: 0 -> L
            sol = eval(evalStr)
            Pout = sol["y"]
            zSpan = np.flip(zSpan)
            pInit = copy.deepcopy(Pout[:, -1])

            # Reset BASE + BCKPUMP values
            pInit[idxPAB] = np.zeros(lenASE * isy)
            pInit[idxPPB] = pumpPmpBck

            # convergence criteria - pump signal power
            if pumpPmpFor == 0:
                errorCvg = 1 - Pout[idxPPB, -1] / pumpPmpBck
            elif pumpPmpBck == 0:
                errorCvg = 1 - Pin[idxPPF, -1] / pumpPmpFor
            else:
                errorCvg = 1 - (
                    np.array([Pout[idxPPB, -1], Pin[idxPPF, -1]])
                ) / np.array([pumpPmpBck, pumpPmpFor])
            logg.info("EDFA SM: loop %2d" % (tryLoop + 1))
            logg.info("Convergence: %5.3f%%.\n" % (100 * np.mean(errorCvg)))

            # Update loop control variable
            tryLoop = tryLoop + 1
            if tryLoop == MaxTry:
                logg.info(
                    "Convergence fail: number of loops greater than max (%d)" % (MaxTry)
                )
        # Automatic gain or power control
        if (param_edfa.type == "AGC") or (param_edfa.type == "APC"):
            power_out = np.sum(
                Pout[np.concatenate([idxPS, idxPAF]), -1]
            )  # Output power - Signal + For. ASE
            if param_edfa.type == "AGC":
                errorAutoCrtl = 10 * np.log10(power_out / power_in)
            else:  # APC
                errorAutoCrtl = 10 * np.log10(1e3 * power_out)
            # PID control - only in forward pumping
            # TODO - and backward pumping? both?
            pid = PID(
                param_edfa.kp, param_edfa.ki, param_edfa.kd, setpoint=param_edfa.value, output_limits=(-pumpPmpFor/2, pumpPmpFor/2)
            )
            pumpPmpFor = pumpPmpFor + pid(errorAutoCrtl)
            errorAutoCrtl = errorAutoCrtl - param_edfa.value

            if np.abs(errorAutoCrtl) > param_edfa.tolCtrl:
                logg.info("EDFA SM: control loop %2d" % (tryCtrlLoop + 1))
                logg.info("Convergence: %5.3f dB" % (errorAutoCrtl))
                logg.info("Pump for.: %5.2f mW\n" % (1e3 * pumpPmpFor))
            tryCtrlLoop = tryCtrlLoop + 1
            if tryCtrlLoop == MaxTry:
                logg.info(
                    "Control fail: number of loops greater than max (%d)" % (MaxTry)
                )
        else:
            errorAutoCrtl = 0
    ## Update pump signal
    PpumpB = Pout[idxPPB, [0, -1]]
    PpumpF = Pout[idxPPF, [0, -1]]

    ## Update signal noise
    # Adjust optical noise level
    freqStep = Fs / lenFqSg
    resolutionOffSet = param_edfa.noiseBand / freqStep
    noiseB = Pout[idxPAB, 0] / resolutionOffSet
    noiseF = Pout[idxPAF, -1] / resolutionOffSet
    # Interpolates the optical noise values ​​and adds phase with normal distribution. It is necessary to divide
    # by sqrt(2) the noise terms, because when adding the phase, the amplitude must be unitary.
    f1_noiseb = interpolate.interp1d(
        freqASE, noiseB[0:lenASE], kind="linear", fill_value="extrapolate"
    )
    f1_noisef = interpolate.interp1d(
        freqASE, noiseF[0:lenASE], kind="linear", fill_value="extrapolate"
    )
    f2_noiseb = interpolate.interp1d(
        freqASE, noiseB[lenASE:], kind="linear", fill_value="extrapolate"
    )
    f2_noisef = interpolate.interp1d(
        freqASE, noiseF[lenASE:], kind="linear", fill_value="extrapolate"
    )
    # Create signal noise
    noiseb = np.concatenate(
        [
            np.sqrt(f1_noiseb(freqSgn), dtype=np.complex),
            np.sqrt(f2_noiseb(freqSgn), dtype=np.complex),
        ]
    )
    noisef = np.concatenate(
        [
            np.sqrt(f1_noisef(freqSgn), dtype=np.complex),
            np.sqrt(f2_noisef(freqSgn), dtype=np.complex),
        ]
    )
    noiseF = (
        noisef
        * (np.random.randn(lenFqSg * isy) + 1j * np.random.randn(lenFqSg * isy))
        / np.sqrt(2)
    )

    ## Update amplified optical signal
    # Update optical signal by adding noise
    Eout = np.reshape(
        np.sqrt(Pout[idxPS, -1], dtype=np.complex), (lenFqSg, isy), order="F"
    )
    Eout = Eout * np.exp(1j * np.angle(EiFt)) + np.reshape(
        noiseF, (lenFqSg, isy), order="F"
    )

    Eout = ifft(Eout * lenFqSg, axis=0)

    return Eout, PpumpF, PpumpB, np.reshape(noisef, (lenFqSg, isy), order="F")