class fiber_doping_profile():
    def __init__(self, parameter):
        self.rho = getattr(parameter, "rho", 0.5e25)
        self.outer_radius = getattr(parameter, "outer_radius", 1.56e-6)
        self.number_sections_radial = getattr(parameter, "number_sections_radial", 100)
        self.number_sections_azimut = getattr(parameter, "number_sections_azimut", 100)
        self.radius_sections = parameter.radius_sections
        self.azimut_sections = parameter.azimut_sections       

        self.nor, self.noa = self.doping_profile()
        self.rho *= np.matmul(self.noa, np.transpose(self.nor))

    def doping_profile(self):
        # Radial profile.
        nor = doping_profile_radial(self.outer_radius,
                                    self.radius_sections, 
                                    self.number_sections_radial)
        # Azimuthal profile.
        noa = doping_profile_azimuthal(self.azimut_sections, 
                                       self.number_sections_azimut)
        # Doping profile
        return nor, noa
    
    def plot(self):
        from matplotlib import cm

        rd = np.linspace(0, self.outer_radius, self.number_sections_radial)
        th = np.linspace(0, 2*np.pi, self.number_sections_azimut)

        fig = plt.figure()
        ax = fig.add_subplot(2,2,1)
        ax.plot(rd, self.nor, 'ro-.')
        ax.set_xlabel(r'Radius [$\mu$m]')
        ax.set_ylabel(r'Er$^{3+}$ density [m$^{-3}$]')
        ax.grid(True)

        ax = fig.add_subplot(2,2,3)
        ax.plot(th, self.noa, 'ro-.')
        ax.set_xlabel(r'Angle [rad]')
        ax.set_ylabel(r'Er$^{3+}$ density [m$^{-3}$]')
        ax.grid(True)

        RD, TH = np.meshgrid(rd, th)
        X, Y = pol2cart(RD, TH)

        ax = fig.add_subplot(2,2,2,projection='3d')
        surf = ax.plot_surface(1e6*X, 1e6*Y, self.rho,  rstride=1, cstride=1, cmap=cm.plasma,
                       linewidth=0, antialiased=False)
        fig.colorbar(surf, shrink=0.5, aspect=10)

        ax = fig.add_subplot(2,2,4,projection='polar')
        pp = ax.contourf(TH, 1e6*RD, self.rho)
        cbar = plt.colorbar(pp, orientation='vertical')
        cbar.ax.set_ylabel('scale label')
        
        plt.tight_layout()
        plt.show()

class doped_fiber(fiber):
    def __init__(self, type, parameter):
        super().__init__(parameter)
        self.tal   = getattr(parameter, "tal", 10e-3)
        self.rho   = getattr(parameter, "rho", 0.95e25)
        self.giles = giles(parameter.gilesFile, 'coeff')
        self.doping= fiber_doping_profile(parameter.doping_profile)

if __name__ == '__main__':
    import os.path as path
    import matplotlib.pyplot as plt
    from optic.core import parameters

    # Set correct filename
    fiber_parameter = parameters()
    fiber_parameter.gilesFile = 'giles_MP980.dat'
    fiber_parameter.gilesFile = path.join(path.abspath(path.join("../")),
                          'OptiCommPy',
                          'optic',
                          'amplifiers', 
                          'ampParams', 
                          fiber_parameter.gilesFile)
    
    # Create fiber object
    passive_fiber = fiber(fiber_parameter)
    # Get normalized frequency
    fig, ax = plt.subplots(1,2)
    wavelength = np.linspace(1500e-9,1600e-9,100)
    ax[0].plot(1e9*wavelength, passive_fiber.v_parameter(wavelength))
    ax[0].set_xlim([1500,1600])
    ax[0].set_xlabel('Wavelenght [nm]')
    ax[0].set_ylabel('Norm. frequency')
    ax[0].grid(True)
    
    # Create doped fiber object
    fiber_parameter.doping_profile = parameters()
    fiber_parameter.doping_profile.radius_sections = np.array([0.85E-06, 1.56e-6])
    fiber_parameter.doping_profile.azimut_sections = np.array([0, 2*np.pi])
    edf_fiber = doped_fiber(fiber_parameter)
   
    # Interpolates data
    wavelength = np.linspace(1500e-9,1600e-9,2000)
    intAbsrCoeff,intGainCoeff = edf_fiber.giles.giles_interpolates(wavelength=wavelength)
    ax[1].plot(1e9*wavelength, intAbsrCoeff, label='Abs. Coeff.')
    ax[1].plot(1e9*wavelength, intGainCoeff, label='Gain Coeff.')
    ax[1].set_xlim([1500,1600])
    ax[1].set_xlabel('Wavelenght [nm]')
    ax[1].set_ylabel('EDF coeff. [dB/m]')
    ax[1].set_title('%d points'%(len(wavelength)))
    ax[1].grid(True)
    ax[1].legend()
    plt.tight_layout()
    plt.show()

    # Plot giles parameters
    edf_fiber.giles.plot('coeff')

    # Plot fiber doping profile
    edf_fiber.doping.plot()