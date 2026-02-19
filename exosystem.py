import numpy as np
from cmf_utility import f_cmf, Abundance_to_CMF # contains SuperEarth code

'''
Utility Class defining exoplanet system to calculate planet mass, radius and CMFs
'''

PREFAC = 28.4329
D2YR = 1/365.15

# All CGS
REARTH = 6.37e8 
MEARTH = 5.97e27
MJUP = 1.898e30
RSUN = 6.957e10

class ExoStar():
    def __init__(self, mass, radius, FeH, MgH, SiH):
        self.mass = mass
        self.radius = radius
        self.FeH = FeH
        self.MgH = MgH
        self.SiH = SiH
        
    # Returns the CMF of the star
    # Methods: 'Fe2Mg' for magnesium or 'Fe2Si' for silicon
    def star_cmf(self, method="Fe2Mg"):
        XH = self.MgH if method=="Fe2Mg" else self.SiH
        return Abundance_to_CMF(Fe_H=self.FeH, X_H=XH, method=method)

class ExoSystem():
    def __init__(self, star:ExoStar, K, period, inclination, rprs, ecc=0):
        self.star = star
        self.K = K
        self.period = period * D2YR # expects period in days
        self.inclination = inclination
        self.rprs = rprs
        self.ecc = ecc
        
    # Returns planet radius in Earth Radii
    def calc_planet_radius(self, rstar=None):
        if rstar == None:
            rstar = self.star.radius
        
        return (self.rprs * rstar) * (RSUN / REARTH)
        
    # Returns planet mass in Earth Masses
    # Assumes mplanet << mstar, and that mstar is in units of solar masses
    def calc_planet_mass(self, mstar=None):
        if mstar == None:
            mstar = self.star.mass
            
        factor = self.K * np.sqrt(1 - self.ecc**2) / PREFAC
        solmass = factor / np.sin(np.deg2rad(self.inclination)) * (mstar**(2/3)) * self.period**(1/3)
        return solmass * (MJUP / MEARTH)
    
    # Returns the CMF of the planet
    def calc_planet_cmf(self, core_silica=0, mantle_iron=0.1):
        pr = self.calc_planet_radius()
        pm = self.calc_planet_mass()
        
        return f_cmf(pm, pr, si=core_silica, fe=mantle_iron)
    
    # Returns the CMF of the star
    # Methods: 'Fe2Mg' for magnesium or 'Fe2Si' for silicon
    def calc_star_cmf(self, method="Fe2Mg"):
        return self.star.star_cmf(method=method)[1][0]