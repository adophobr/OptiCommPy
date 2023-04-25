import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, kv
from giles import giles

def silica_index(wavelength: np.ndarray) -> np.ndarray:
    return np.sqrt((1 + 0.6961663 / (1 - (0.0684043 / wavelength) ** 2) + 
            0.4079426 / (1 - (0.1162414 / wavelength) ** 2) + 
            0.8974794 / (1 - (9.896161 / wavelength) ** 2)))

def normalized_frequency(radius: float,
                         na: float,
                         wavelength: np.ndarray) -> np.ndarray:
    return (2 * np.pi / wavelength) * radius * na
  
class fiber():
    def __init__(self, parameter):
        # init geometric and material parameters
        self.na = getattr(parameter, "na", 0.22)
        self.radius = getattr(parameter, "radius", 1.56e-6)
        self.length = getattr(parameter, "length", 8)
        self.core_index = getattr(parameter, "core_index", 1.45)
    
    @property
    def clad_index(self):
        return np.sqrt(self.core_index**2 - self.na**2)

    def v_parameter(self,
                    wavelength: np.ndarray) -> np.ndarray:
        return normalized_frequency(self.radius, 
                                    self.na, 
                                    wavelength)
    
    def set_radius(self):
        # Logitudinal step
        self.dr = self.a / self.longitudinal_steps
        self.r  = np.arange(0, self.a, self.dr)
   
    def set_gamma(self):
        if (self.geometric == 'LP01'):
            self.gamma = (((self.v * self.b) / (self.a * self.V * jv(1, self.u))) ** 2) * (jv(0, self.u * self.b / self.a) ** 2 + jv(1, self.u * self.b / self.a) ** 2)
        else:
            self.gamma = 1 - np.exp(-2 * (self.b / self.w_gauss) ** 2)

    # Call only if self.algorithm == "Giles_spatial"
    def set_ik(self):
        if (self.geometric == 'LP01'):
            i_k = (lambda r: 1/np.pi*(self.v/(self.a*self.V)*jv(0, self.u*self.r/self.a)/jv(1, self.u))**2)
        else:
            i_k = lambda r: 2 / (np.pi * self.w_gauss ** 2) * np.exp(-2 * (self.r / self.w_gauss) ** 2)
        self.i_k = [i_k(x) for x in self.r]
    
