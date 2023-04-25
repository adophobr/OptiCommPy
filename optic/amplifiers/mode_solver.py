from scipy.special import jn_zeros
import numpy as np
import numpy.typing as npt

#def core_index(wavelength: npt.ArrayLike,
               #na: float) -> npt.ArrayLike:
    #return np.sqrt(silica_index(wavelength) ** 2 + na ** 2)
#
#def u_parameter(radius: float,
                #na: float,
                #wavelength: npt.ArrayLike,
                #eff_rindex: npt.ArrayLike,
                #core_index_function = core_index) -> npt.ArrayLike:
    #core_rindex = core_index_function(wavelength, na)
    #return (2 * np.pi / wavelength) * radius * np.sqrt(core_rindex**2 - eff_rindex**2)
#
#def w_parameter(radius: float,
                #wavelength: npt.ArrayLike,
                #eff_rindex: npt.ArrayLike,
                #clad_index_function = silica_index) -> npt.ArrayLike:
    #clad_rindex = clad_index_function(wavelength)
    #return (2 * np.pi / wavelength) * radius * np.sqrt(eff_rindex**2 - clad_rindex**2)

    """
class fiber_mode():
    def __init__(self, 
                 l_max:float=2, 
                 m_max:float=2,
                 radius:float=1.56e-6, 
                 na:float=0.22, 
                 wavelength:npt.ArrayLike) -> npt.ArrayLike:

    @property
    def propagation_constant(self):
        mode_solver = mode_solver()
        return None

    @property
    def norm_prop_constat(self):
        return None
    
    @property
    def effective_index(self):
        return None
    
    @property
    def phase_velocity(self):
        return None
    
    @property
    def group_velocity(self):
        return None
    
    @property
    def dispersion(self):
        return None
    
    @property
    def confinement_factor(self):
        return None

    def field_intensity(self):
        return None
    
    def plot_field(self):
        return None
"""

class mode_solver():
    """
    Considering LP (lineary polarized) modes in weakly guiding cylindrical step-index fiber.
    """
    def __init__(self, 
                lmax:int,
                mmax:int) -> None:
        self.l = lmax
        self.m = mmax

    def table_modes(self) -> None:
        table = np.zeros([(self.l)*self.m,4])        
        lower = np.zeros([self.l,self.m])
        for l in range(self.l):
            table[l*self.m:(l+1)*self.m,0] = l
            if (l!=0):
                table[l*self.m:(l+1)*self.m,-2] = jn_zeros(l-1, self.m)
            else:
                table[1+l*self.m:(l+1)*self.m,-2] = jn_zeros(1, self.m-1)
            table[l*self.m:(l+1)*self.m,-1] = jn_zeros(l, self.m)
        table[:,1] = 1+np.tile(np.arange(self.m),self.l)       
        self.table = table

    def solve_modes(self, 
                    radius:float, 
                    na:float, 
                    wavelength:npt.ArrayLike) -> npt.ArrayLike:
        v = v_parameter(radius,na,wavelength)
        b = np.zeros([self.l,self.m,len(v)])
        for l in range(self.l):
            for m in range(self.m):
                if (l == 0) and (m == 0):
                    b[l,m,:] = self.solve_he11(v)
                else:
                    b[l,m,:] = self.solve_higher_modes(l, m, v)
        return b

    def solve_he11(self, 
                   v:npt.ArrayLike) -> npt.ArrayLike:
        return 1 - ((1 + np.sqrt(2)) / (1 + (4 + v ** 4) ** 0.25)) ** 2
    
    def solve_higher_modes(self, 
                        l:float,
                        m:float,
                        v:npt.ArrayLike) -> npt.ArrayLike:
        A = np.pi * ((m+1) + 0.5 * (l - 1) - 0.25)
        B = 4 * (l - 1) ** 2
        Uc = A - (B - 1) / (8 * A) - (4 * (B - 1) * (7 * B - 31)) / (3 * (8 * A) ** 3)
        S = np.sqrt(Uc ** 2 - l ** 2 - 1)
        return 1 - (Uc / v) ** 2 * np.exp(2 / S * (np.arcsin(S / Uc) - np.arcsin(S / v)))

if __name__ == '__main__':
    l_max = 4
    m_max = 3
    wavelength = np.array([980e-9, 1550e-9])
    mode = mode_solver(l_max,m_max)
    mode.table_modes()
    print(v_parameter(radius=9e-6,na=0.22,wavelength=wavelength))
    b = mode.solve_modes(radius=9e-6,na=0.22,wavelength=wavelength)    
    print(np.hstack((mode.table,b.reshape(l_max*m_max,len(wavelength)))))