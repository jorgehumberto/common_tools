#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python modules
__author__ = 'jmartins'

import numpy as np

from constants import mJup, mSun, mEarth, rJup, AU


# =======================================================================
#                            CLASSES                                   #
# =======================================================================

# =======================================================================
class clsPlanetParametersOLDOLD:
    """ clsOrbitalParameters:"""

    def __init__(self, params={'Stellar': 1, 'ecc': 0, 'Rp': 0.0}):
        self.period = float(params['p'])  # days
        self.m2sini = float(params['m2sini'])  # Jupiter masses
        self.hostmass = float(params['Stellar'])  # solar masses
        self.ecc = float(params['ecc'])
        self.a = float(params['a'])  # astronomical units
        try:
            self.massratio = float(params['massratio'])
        except:
            self.massratio = self.m2sini / self.hostmass * (mJup / mSun)  # planet/star mass ratio
        self.k1 = float(params['k1']) / 1000  # km/s per second
        self.k2 = float(params['k1']) / self.massratio / 1000  # km/s per second
        self.t0 = float(params['t0'])  # Barycentric Julian Date
        self.w = float(np.radians(float(params['om'])))  # Longitude of periastron
        try:
            self.Sysrv = float(params['Sysrv'])
        except:
            self.Sysrv = float(params['delta3'])
        self.PhaseZero = 0.0
        self.albedo = 0.0
        try:
            self.STD = float(params['STD'])
        except:
            self.STD = 0.0
        try:
            self.Rp = float(params['Rp'])
        except:
            self.Rp = 0.0
        self.I = 90.
        self.FWHMMax = 0.0

        self.PlanetListLenght = 0.0
        self.Significance = 0.
        self.SigErr = 0.

        # Confidence interval
        self.Sig1 = 0.0
        self.Sig2 = 0.0
        self.Sig3 = 0.0

        self.SNFinal = 0.0

    def setFWHMMax(self, FWHMMax):
        self.FWHMMax = float(FWHMMax)

    def GaussFitParams(self, ParamsList, FitErrors):
        """
        stores the fits params as:
        [Amplitude, mean, FWHM, B, TimeStamp]
        """
        self.Amplitude = float(ParamsList[0])
        self.mean = float(ParamsList[1])
        self.FWHM = float(ParamsList[2])
        self.B = float(ParamsList[3])
        self.AmpErr = float(FitErrors[0])
        self.meanErr = float(FitErrors[1])
        self.FWHMErr = float(FitErrors[2])
        self.BErr = float(FitErrors[3])
        self.TimeStamp = float(ParamsList[4])


class clsPlanetParametersOLD:
    """ clsOrbitalParameters:"""

    def __init__(self, params={'massPlanet': 1., 'massStar': 1., 'a': 0.04, 'period': 3., 'ecc': 0., 't0': 0.,\
                               'w': 0, 'I': 90., 'sysRV': 0.0, 'k1': 150.0, 'radiusPlanet':1.0}):
        self.massPlanet = float(params['massPlanet']) * mJup
        self.massStar = float(params['massStar']) * mSun
        self.massRatio = self.massPlanet / self.massStar
        self.a = float(params['a']) * AU
        self.period = float(params['period'])
        self.ecc = float(params['ecc'])
        self.t0 = float(params['t0'])
        # self.phase_offset  = params['phase_offset']
        self.w = float(params['w'])
        self.I = float(params['I'])
        self.k1 = float(params['k1']) / 1000.
        self.sysRV = float(params['sysRV'])
        self.k2 = float(params['k1']) / 1000. * self.massStar / self.massPlanet
        self.radiusPlanet = float(params['radiusPlanet']) * rJup
