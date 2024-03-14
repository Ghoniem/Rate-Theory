import numpy as np


class MaterialProperties:
    def __init__(self, T, nu_v, nu_g, nu_i, a0, b, G, G_He, f, Em_g,
                 Em_v, Em_i, Eb_v_g, Eb_v_2g, Eb_2g, Ef_v, Z_v, Z_i, rho, floor, epsilon):
        self.k_B = 8.617333262145e-5  # Boltzmann constant in eV/K
        self.T = T  # Temperature in Kelvin
        self.nu_v = nu_v  # Vacancy attempt frequency in Hz
        self.nu_g = nu_g  # Helium attempt frequency in Hz
        self.nu_i = nu_i  # SIA attempt frequency in Hz
        self.a0 = a0  # Lattice constant in meters
        self.b = b  # Re-solution parameter (dimensionless)
        self.G = G  # Displacement damage rate (example value)
        self.G_He = G_He
        self.f = f  # Fraction of surviving defects
        self.Em_g = Em_g  # Helium Migration Energy in eV
        self.Em_v = Em_v  # Vacancy Migration Energy in eV
        self.Em_i = Em_i  # SIA Migration Energy in eV
        self.Eb_v_g = Eb_v_g  # Helium Substitutional Binding Energy in eV
        self.Eb_v_2g = Eb_v_2g  # Di-Helium Substitutional Binding Energy in eV
        self.Eb_2g = Eb_2g  # Di-Helium Interstitial Binding Energy in eV
        self.Ef_v = Ef_v  # Vacancy Formation Energy in eV
        self.Z_v = Z_v  # Vacancy Bias Factor
        self.Z_i = Z_i  # SIA Bias Factor
        self.rho = rho  # Dislocation Density
        self.floor = floor  # minimum fractional concentration
        self.epsilon = epsilon  # combinatorial number for bubbles
        self.calculate_properties()

    def calculate_properties(self):
        # Reaction frequencies
        self.alpha = 48 * self.nu_i * np.exp(-self.Em_i / (self.k_B * self.T))
        self.beta = 48 * self.nu_g * np.exp(-self.Em_g / (self.k_B * self.T))
        self.gamma = 48 * self.nu_v * np.exp(-self.Em_v / (self.k_B * self.T))

        # Thermal emission frequencies
        self.e1 = np.exp(-self.Eb_v_g / (self.k_B * self.T))
        self.e2 = np.exp(-self.Eb_v_2g / (self.k_B * self.T))
        self.e3 = np.exp(-self.Eb_2g / (self.k_B * self.T))

        # Additional parameters
        self.delta = self.b * self.G
        self.Omega = self.a0 ** 3 / 4  # Atomic volume

        # Diffusion coefficients and factor
        self.D_i = (self.a0 ** 2 / 48) * self.alpha
        self.D_v = (self.a0 ** 2 / 48) * self.gamma
        self.D_g = (self.a0 ** 2 / 48) * self.beta

        # Sink concentrations
        self.C_s_v = (self.a0 ** 2 / 48) * self.Z_v * self.rho
        self.C_s_i = (self.a0 ** 2 / 48) * self.Z_i * self.rho

        # Equilibrium vacancy concentration
        self.C_v_e = np.exp(-self.Ef_v / (self.k_B * self.T))

