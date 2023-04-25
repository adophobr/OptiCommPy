#%%
import os
import numpy as np
import matplotlib.pyplot as plt

#%%
class giles():
    def __init__(self,
                 file:str,
                 type:str='coeff') -> None:
        self.file = file
        if not (os.path.exists(self.file)):
            raise TypeError(f"{self.file} file doesn't exist.")
        self.type = type
        # set frequency
        self.set_wavelength_from_file()
        if (self.type=='coeff'):
            self.set_coefficients_from_file()
        elif (self.type=='cross'):
            self.set_cross_section_from_file()
        else:
            raise TypeError("gile.type invalid argument - ['coeff', 'cross'].")

    def set_wavelength_from_file(self) -> None:
        # Waiting wavelength in nm
        self.wavelength = 1e-9*np.loadtxt(self.file)[:,0]
    
    def set_coefficients_from_file(self) -> None:
        file_data = np.loadtxt(self.file)
        self.absorption_coefficient  = 0.1 * np.log(10) * file_data[:,1]
        self.gain_coefficient = 0.1 * np.log(10) * file_data[:,2]

    def set_coefficients_from_cross_section(self, 
                                            rho:np.ndarray, 
                                            gamma:np.ndarray) -> None:
        self.absorption_coefficient = self.absorption_cross * rho * gamma
        self.gain_coefficient = self.emission_cross * rho * gamma

    def set_cross_section_from_file(self) -> None:
        file_data = np.loadtxt(self.file)
        self.absorption_cross = file_data[:, 1]
        self.emission_cross = file_data[:, 2]

    def set_cross_section_from_coefficients(self, 
                                            rho:np.ndarray, 
                                            gamma:np.ndarray) -> None:
         self.absorption_cross = self.absorption_coefficient / rho / gamma
         self.emission_cross = self.gain_coefficient / rho / gamma

    def plot(self, 
             type:str) -> None:
        if (type=='coeff'):
            plt.plot(1e9*self.wavelength, self.absorption_coefficient, label='Abs. coeff.')
            plt.plot(1e9*self.wavelength, self.gain_coefficient, label='Gain coeff.')
            plt.ylabel('EDF coeff. [dB/m]')
        elif (type=='cross'):
            plt.plot(1e9*self.wavelength, 1e25*self.absorption_cross, label='Abs. cross-section')
            plt.plot(1e9*self.wavelength, 1e25*self.emission_cross, label='Em. cross-section')
            plt.ylabel(r'Cross-section [10$^{-25}$ m$^2$]')
        else:
            raise TypeError("type invalid argument - ['coeff', 'cross'].")

        plt.xlabel('Wavelenght [nm]')
        plt.title('%d points'%(len(self.wavelength)))
        plt.grid(True)
        plt.legend()
        plt.show()
    
    def giles_interpolates(self,
                           wavelength:np.array,
                           type:str='coeff') -> np.ndarray:
        if (type=='coeff'):
            return self.giles_coeff_interpolates(wavelength)
        elif (type=='cross'):
            return self.giles_cross_interpolates(wavelength)
        else:
            raise TypeError("type invalid argument - ['coeff', 'cross'].")

    def giles_coeff_interpolates(self,
                                 wavelength:np.ndarray) -> np.ndarray:
        intAbsrCoeff = np.interp(wavelength, self.wavelength, self.absorption_coefficient) 
        intGainCoeff = np.interp(wavelength, self.wavelength, self.gain_coefficient) 
        return intAbsrCoeff,intGainCoeff
    
    def giles_cross_interpolates(self,
                                 wavelength:np.ndarray) -> np.ndarray:
        intAbsrCross = np.interp(wavelength, self.wavelength, self.absorption_cross) 
        intEmssCriss = np.interp(wavelength, self.wavelength, self.emission_cross)
        return intAbsrCross,intEmssCriss

#%%
if __name__ == '__main__':
    import os.path as path
    # Set correct filename
    gilesFile = 'giles_MP980.dat'
    gilesFile = path.join(path.abspath(path.join("../")),
                          'OptiCommPy',
                          'optic',
                          'amplifiers', 
                          'ampParams', 
                          gilesFile)
    edf_giles = giles(gilesFile)
    # Plot giles parameters using data from file
    edf_giles.plot('coeff')
    # Interpolates data
    plt.figure()
    wavelength = np.linspace(1500e-9,1600e-9,2000)
    intAbsrCoeff,intGainCoeff = edf_giles.giles_interpolates(wavelength=wavelength)
    plt.plot(1e9*wavelength, intAbsrCoeff, label='Abs. Coeff.')
    plt.plot(1e9*wavelength, intGainCoeff, label='Gain Coeff.')
    plt.xlim([1500,1600])
    plt.xlabel('Wavelenght [nm]')
    plt.ylabel('EDF coeff. [dB/m]')
    plt.title('%d points'%(len(wavelength)))
    plt.grid(True)
    plt.legend()
# %%
