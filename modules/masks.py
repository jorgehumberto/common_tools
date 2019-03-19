#!/usr/bin/env python
# -*- coding: utf-8 -*-
########################################################################
#							PYHTON MODULES							   #
########################################################################
import numpy as np


########################################################################
#							INPUT FUNCTIONS							   #
########################################################################

def fnGetMask(filename, waveArray):
    ''' Function that given the  absolute path of a file containing a mask
        will turn it into a numpy array.
        Input parameters:
            filename - absolute path of file containing the mask data
            waveArray - numpy array containing the wavelength ranges data
                        (can have several orders for echelles spectra)

        Input mask data format - 3 columns:
            column 1 - wavelength of beginning of spectral line (lineWavelengthEnd)
            column 2 - wavelength of end of spectral line (lineWavelengthEnd)
            column 3 - depth of spectral line (lineDepth)

        Outpu mask:
            will have the same dimensions of waveArray


    '''
    lineWavelengthIni, lineWavelengthEnd, lineDepth = np.loadtxt(filename, \
                                                                 dtype='float', unpack=True)

    mask = np.ones_like(waveArray)

    for maskElement, waveElement in zip(np.nditer(mask, op_flags=['readwrite']), np.nditer(waveArray)):
        for lineIni, lineEnd, depth in zip(lineWavelengthIni, lineWavelengthEnd, lineDepth):

            if lineIni < waveElement < lineEnd:
                print lineIni, waveElement, lineEnd
                maskElement = 1. - depth
        print waveElement

    return mask
