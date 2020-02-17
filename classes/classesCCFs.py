#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python modules

import numpy as np
from astropy.io import fits
import os
#from astropy import units as u
from astropy.coordinates import SkyCoord

# =======================================================================
#                            CLASSES                                   #
# =======================================================================


class clsUVESCCF:
    """ clsHARPSCCFLastOrder:"""

    def __init__(self, data, header=None):
        # CCFFile             =     fits.open( FitsFileName )[0].copy()
        self.data = np.array(data).copy()

        if header != None:
            self.CCFStep = header['CDELT1']
            self.CCFWaveIni = header['CRVAL1']
            self.RVC = header['HIERARCH ESO DRS CCF RVC']
            self.RV = header['HIERARCH ESO DRS CCF RV']
            self.RVError = header['HIERARCH ESO DRS CCF NOISE']
            self.FWHM = header['HIERARCH ESO DRS CCF FWHM']
            self.ObsTime = header['MJD-OBS']  # in days!!!
            self.ExpTime = header['EXPTIME']  # in seconds
            self.CCFLines = header['HIERARCH ESO DRS CCF LINES']
            self.MJD = header['MJD-OBS']  # FOR SIMULATED DATA!!!!
            self.contrast = header['HIERARCH ESO DRS CCF CONTRAST']
        else:
            self.CCFStep = 1.0
            self.CCFWaveIni = 0.0
            self.RV = 0.0
            self.RVC = 0.0
            self.RVError = 0.0
            self.FWHM = 0.0
            self.ObsTime = 0.0  # in days!!!
            self.ExpTime = 0.0  # in seconds
            self.CCFLines = 0.0
            self.MJD = 0.0  # FOR SIMULATED DATA!!!!
            self.contrast = (np.nanmax(self.data) - np.nanmin(self.data)) / np.nanmax(self.data)

        self.flux = np.nansum(self.data)
        self.wave = (np.arange(len(self.data), dtype=np.float64) * self.CCFStep + self.CCFWaveIni)
        self.ccfAmplitude = 0.0
        self.ccfMeanPixel = 0.0  # in pixels
        self.ccfFWHMPixels = 0.0  # in pixels
        self.ccfB = 0.0
        self.planetRV = 0.0  # planet RV in km/s
        self.planetRVMeanPixel = 0.0  # predicted position of planet along the array



    def getPhase(self, period, t0):
        self.PlanetPhase = (self.MJD - t0) / period
        self.PlanetPhaseFolded = self.PlanetPhase % 1

    def shiftCCF(self, RV):
        """
        Shift the CCF of a given RV
        """
        self.data = np.interp(self.wave.copy() + RV, self.wave.copy(), self.data.copy())
        # self.wave += RV
        self.RVC += RV

    def normCCF(self, Template):
        """
        Normalize the CCF by a template
        """
        self.data = np.divide(self.data.copy(), Template.copy() * self.flux.copy() / np.nansum(Template.copy()))
        self.flux = np.nansum(self.data.copy())


class clsHARPSCCF:
    """ clsHARPSCCFLastOrder:"""

    def __init__(self, data, header=None):
        # CCFFile             =     fits.open( FitsFileName )[0].copy()
        self.data = np.array(data).copy()
        self.indices = np.indices(np.shape(self.data))[1,:,:]

        if header != None:
            if 'eso' in header['TELESCOP'].lower():
                observatory = 'ESO'
                starCoords = SkyCoord(ra=header['RA'], dec = header['DEC'], unit = 'deg')
            elif 'tng' in header['TELESCOP'].lower():
                observatory = 'TNG'
                starCoords = SkyCoord(ra=header['RA-RAD'], dec = header['DEC-RAD'], unit = 'rad')

            self.CCFStep = header['CDELT1']
            self.CCFWaveIni = header['CRVAL1']
            self.RV = header['HIERARCH {} DRS CCF RV'.format(observatory)]
            self.RVC = header['HIERARCH {} DRS CCF RVC'.format(observatory)]
            #self.RVError = header['HIERARCH {} DRS CCF NOISE'.format(observatory)]
            self.FWHM = header['HIERARCH {} DRS CCF FWHM'.format(observatory)]
            self.ObsDate = header['DATE-OBS'].split('T')[0]  # in days!!!
            self.ExpTime = header['EXPTIME']  # in seconds
            self.CCFLines = header['HIERARCH {} DRS CCF LINES'.format(observatory)]
            self.BJD = header['HIERARCH {} DRS BJD'.format(observatory)]  # FOR SIMULATED DATA!!!!
            self.MJD = header['MJD-OBS'] + 2400000
            self.obsDate = header['DATE-OBS']
            self.contrast = header['HIERARCH {} DRS CCF CONTRAST'.format(observatory)] / 100
            self.BERV = header['HIERARCH {} DRS BERV'.format(observatory)]
            self.SN50 = header['HIERARCH {} DRS SPE EXT SN50'.format(observatory)]
            self.nLines = header['HIERARCH {} DRS CCF LINES'.format(observatory)]
            try:
                self.AirMass = header['HIERARCH {} TEL AIRM START'.format(observatory)]
            except:
                self.AirMass = header['AIRMASS']


            self.ra = starCoords.ra.degree
            self.dec = starCoords.dec.degree

        else:
            self.CCFStep = 1.0
            self.CCFWaveIni = 0.0
            self.RV = 0.0
            self.RVC = 0.0
            self.RVError = 0.0
            self.FWHM = 0.0
            self.ObsDate = 0.0  # in days!!!
            self.ExpTime = 0.0  # in seconds
            self.CCFLines = 0.0
            self.BJD = 0.0  # FOR SIMULATED DATA!!!!
            self.contrast = (np.nanmax(self.data) - np.nanmin(self.data)) / np.nanmax(self.data)
            self.SN50 = 0.0
            self.nLines = 0.0
            self.AirMass = 0.0

        self.flux = np.nansum(self.data)
        self.wave = self.indices.copy()
        self.ccfAmplitude = 0.0
        self.ccfMeanPixel = 0.0  # in pixels
        self.ccfFWHMPixels = 0.0  # in pixels
        self.ccfB = 0.0
        self.planetRV = 0.0  # planet RV in km/s
        self.planetRVMeanPixel = 0.0  # predicted position of planet along the array
        self.planetRVRange = []  # planet RV in km/s
        self.phaseFunction = 0.0

    # def getPhase(self, period, t0, phaseZero = 0):
    #     self.PlanetPhase = (self.BJD - t0) / period + phaseZero
    #     self.PlanetPhaseFolded = self.PlanetPhase % 1

    def shiftCCF(self, RV):
        """
        Shift the CCF of a given RV
        """
        self.data = np.interp(self.wave.copy() + RV, self.wave.copy(), self.data.copy())
        # self.wave += RV
        self.RVC += RV

    def normCCF(self, Template):
        """
        Normalize the CCF by a template
        """
        self.data = np.divide(self.data.copy(), Template.copy() * self.flux.copy() / np.nansum(Template.copy()))
        self.flux = np.nansum(self.data.copy())


    def saveCCF(self,fileName,clobber= False):
        fits.writeto(fileName, self.data, self.header, clobber=clobber)
        return

    def updateGaussParams(self, params, order=None):
        '''
        blah
        :param params:
        :param order:
        :return:
        '''

        if order is None:
            self.ccfAmplitude = params[0]
            self.ccfMeanPixel = params[1]
            self.ccfFWHMPixels = params[2]
            self.ccfB = params[3]
        else:
            self.ccfAmplitude[order] = params[0]
            self.ccfMeanPixel[order] = params[1]
            self.ccfFWHMPixels[order] = params[2]
            self.ccfB[order] = params[3]


class clsFits_CCF:
    """ clsHARPSCCFLastOrder:"""

    def __init__(self, data, header=None, instrument= 'espresso'):
        # CCFFile             =     fits.open( FitsFileName )[0].copy()
        self.data = np.array(data).copy()
        self.indices = np.indices(np.shape(self.data))[1,:,:]

        from ConfigParser import ConfigParser
        config = ConfigParser()
        config.optionxform = str
        config.read(os.path.abspath(os.path.join(os.environ['COMMONPATH'],'keywords/kwList_fits.cfg')))

        kwFits={key:val for key,val in config._sections[instrument].items()}

        # print kwFits

        if header != None:

            # if 'eso' in header['TELESCOP'].lower():
            #     observatory = 'ESO'
            #     starCoords = SkyCoord(ra=header['RA'], dec = header['DEC'], unit = 'deg')
            # elif 'tng' in header['TELESCOP'].lower():
            #     observatory = 'TNG'
            #     starCoords = SkyCoord(ra=header['RA-RAD'], dec = header['DEC-RAD'], unit = 'rad')

            self.CCFStep = header[kwFits['ccf_step']]
            self.CCFWaveIni = header[kwFits['ccf_waveIni']]
            self.RV = header[kwFits['rv']]
            self.FWHM = header[kwFits['ccf_fwhm']]
            self.ObsDate = header[kwFits['obs_date']]  # in days!!!
            self.ExpTime = header[kwFits['exp_time']]  # in seconds
            self.CCFLines = 0#header[kwFits['n_lines']]
            self.MJD = header[kwFits['mjd']] + 2400000
            self.BJD = header[kwFits['bjd']]  # FOR SIMULATED DATA!!!!
            self.contrast = header[kwFits['ccf_contrast']] / 100
            self.BERV = header[kwFits['berv']]
            self.AirMass = header[kwFits['airmass_start']]

        else:
            self.CCFStep = 1.0
            self.CCFWaveIni = 0.0
            self.RV = 0.0
            self.RVC = 0.0
            self.RVError = 0.0
            self.FWHM = 0.0
            self.ObsDate = 0.0  # in days!!!
            self.ExpTime = 0.0  # in seconds
            self.CCFLines = 0.0
            self.BJD = 0.0  # FOR SIMULATED DATA!!!!
            self.contrast = (np.nanmax(self.data) - np.nanmin(self.data)) / np.nanmax(self.data)
            self.SN50 = 0.0
            self.nLines = 0.0
            self.AirMass = 0.0

        self.flux = np.nansum(self.data)
        self.wave = self.indices.copy()
        self.ccfAmplitude = 0.0
        self.ccfMeanPixel = 0.0  # in pixels
        self.ccfFWHMPixels = 0.0  # in pixels
        self.ccfB = 0.0
        self.planetRV = 0.0  # planet RV in km/s
        self.planetRVMeanPixel = 0.0  # predicted position of planet along the array
        self.planetRVRange = []  # planet RV in km/s
        self.phaseFunction = 0.0

    # def getPhase(self, period, t0, phaseZero = 0):
    #     self.PlanetPhase = (self.BJD - t0) / period + phaseZero
    #     self.PlanetPhaseFolded = self.PlanetPhase % 1

    def shiftCCF(self, RV):
        """
        Shift the CCF of a given RV
        """
        self.data = np.interp(self.wave.copy() + RV, self.wave.copy(), self.data.copy())
        # self.wave += RV
        self.RVC += RV

    def normCCF(self, Template):
        """
        Normalize the CCF by a template
        """
        self.data = np.divide(self.data.copy(), Template.copy() * self.flux.copy() / np.nansum(Template.copy()))
        self.flux = np.nansum(self.data.copy())


    def saveCCF(self,fileName,clobber= False):
        fits.writeto(fileName, self.data, self.header, clobber=clobber)
        return

    def updateGaussParams(self, params, order=None):
        '''
        blah
        :param params:
        :param order:
        :return:
        '''

        if order is None:
            self.ccfAmplitude = params[0]
            self.ccfMeanPixel = params[1]
            self.ccfFWHMPixels = params[2]
            self.ccfB = params[3]
        else:
            self.ccfAmplitude[order] = params[0]
            self.ccfMeanPixel[order] = params[1]
            self.ccfFWHMPixels[order] = params[2]
            self.ccfB[order] = params[3]
