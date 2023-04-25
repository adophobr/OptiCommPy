import numpy as np

def u_parameter(radius: float,
                na: float,
                wavelength: npt.ArrayLike,
                eff_rindex: npt.ArrayLike,
                core_index_function = core_index) -> npt.ArrayLike:
    core_rindex = core_index_function(wavelength, na)
    return (2 * np.pi / wavelength) * radius * np.sqrt(core_rindex**2 - eff_rindex**2)

def w_parameter(radius: float,
                wavelength: npt.ArrayLike,
                eff_rindex: npt.ArrayLike,
                clad_index_function = silica_index) -> npt.ArrayLike:
    clad_rindex = clad_index_function(wavelength)
    return (2 * np.pi / wavelength) * radius * np.sqrt(eff_rindex**2 - clad_rindex**2)

class mode():
    def set_u_modal_parameter(self, wavelength):
        v_parameter = self.set_v_parameter(self, wavelength)
        return ((1 + np.sqrt(2)) * v_parameter) / (1 + (4 + v_parameter ** 4) ** 0.25)

    @staticmethod
    def set_w_modal_parameter(self):
        return np.sqrt(v_parameter ** 2 - u_parameter ** 2)

    @staticmethod
    def mode_radius_profile_bessel(v_parameter, u_parameter, radius):
        w_parameter = fiber.fiber_w_modal_parameter(v_parameter, u_parameter)
        return radius * v_parameter / u_parameter * kv(1, w_parameter) / kv(0, w_parameter) * jv(0, u_parameter)
    
    @staticmethod
    def mode_radius_profile_marcuse(v_parameter, radius):
        return radius * (0.650 + 1.619 / v_parameter ** 1.5 + 2.879 / v_parameter ** 6)
    
    @staticmethod
    def mode_radius_profile_whitley(v_parameter, radius):
        return radius * (0.616 + 1.660 / v_parameter ** 1.5 + 0.987 / v_parameter ** 6)
    
    @staticmethod
    def mode_radius_profile_desurvire(v_parameter, radius):
        return radius * (0.759 + 1.289 / v_parameter ** 1.5 + 1.041 / v_parameter ** 6)
    
    @staticmethod
    def mode_radius_profile_myslinski(v_parameter, radius):
        return radius * (0.761 + 1.237 / v_parameter ** 1.5 + 1.429 / v_parameter ** 6)