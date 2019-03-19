# -*- coding: utf-8 -*-
"""
Created on Tue Set 08 2015

@author: Jorge Martins
email: jorge.martins@iastro.pt


"""
###############################################################################
# Python Modules
###############################################################################

import lmfit as lm
import numpy as np
import scipy.optimize as opt


###############################################################################
# Functions
###############################################################################

def fnGauss(data, params):
    '''
    Function to model a Gaussian function with the following parameters:
        params[0] = A               params[1] = mean
        params[2] = FWHM            params[3] = B

        G(data) = B - A * exp(-4*ln(2)*(data - mean)��/ FWHM��)


    It is derived from the common gaussian function:

        G(data) = B + A * exp(-(data - mean)��/ (2*sigma��))

    by defining

        sigma = abs(FWHM)/(2*sqrt(2*ln2))

    '''
    # print 'params',params

    A, mean, FWHM, B = params
    dataAveraged = B - A * np.exp(-4 * np.log(2) * (((data - mean) / FWHM) ** 2))

    return dataAveraged


# -----
def fnGaussFit(data, dataAveraged, GaussParamsInitGuess=[0, 0, 0, 0], fixed=[]):
    # def fnGaussianFit(data, dataAveraged, GaussParamsInitGuess=[0, 0, 0, 0], AmpMax=1.0e99, AmpMin=-1.0e99, \
    #                   MeanMax=1.0e99, MeanMin=-1.0e99, FWHMMax=1.0e99, FWHMMin=-1.0e99, \
    #                   BMax=1.0e99, BMin=-1.0e99, fixed=[]):
    '''
    '''
    resGauss = lambda params, resdata, resdataAveraged: \
        fnGauss(resdata,
                [params['Amplitude'].value, params['mean'].value, params['FWHM'].value, params['B'].value]) - resdataAveraged

    GaussParams = lm.Parameters()
    GaussParams.add('Amplitude', value=GaussParamsInitGuess[0], vary=True)
    GaussParams.add('mean', value=GaussParamsInitGuess[1], vary=True)
    GaussParams.add('FWHM', value=GaussParamsInitGuess[2], vary=True)
    GaussParams.add('B', value=GaussParamsInitGuess[3], vary=True)

    for fix in fixed:
        GaussParams[fix].vary = False

    data = np.array(data)
    GaussParamsMinimize = lm.minimize(resGauss, GaussParams, args=(data, dataAveraged), method='leastsq', \
                                      # maxfev=100000000,
                                      xtol=1.e-13 \
                                      )

    return GaussParamsMinimize.params['Amplitude'].value, GaussParamsMinimize.params['mean'].value, \
           GaussParamsMinimize.params['FWHM'].value, GaussParamsMinimize.params[
               'B'].value  , GaussParams['Amplitude'].stderr, GaussParams['mean'].stderr,GaussParams['FWHM'].stderr, GaussParams['B'].stderr

# ----------------------------------------------------------------------------------------------------------------------
def fnGaussianFitOLD(data, dataAveraged, GaussParamsInitGuess=[0, 0, 1, 0]):
    '''
    '''
    resGauss = lambda params, resdata, resdataAveraged: fnGauss(resdata, params) - resdataAveraged

    data = np.array(data)

    GaussParams, cov, infodict, mesg, ier = opt.leastsq(resGauss, GaussParamsInitGuess, \
                                                        args=(data, dataAveraged), maxfev=1000000000, full_output=True)

    return GaussParams

# ----------------------------------------------------------------------------------------------------------------------
def fnMovingAverage(data, N):
    '''
    :param data: data
    :param N: width of moving average window
    :return: moving average
    adapted from http://mathworld.wolfram.com/MovingAverage.html
    moving average with pixel centered on the moving average window

    '''
    lenData = len(data)
    halfWindow = int(N/2)
    dataAveraged = np.zeros_like(data)

    # initial point (due to divide by zero)
    dataAveraged[0] = data[0]

    # initial adaptive window
    dataAveraged[1:halfWindow] = [np.nansum(data[:2*index -1])/(2*index-1) for index in np.arange(1, halfWindow)]

    # fixed window
    dataAveraged[halfWindow:-halfWindow] = [np.nansum(data[index-halfWindow:index + halfWindow ])/N \
                                   for index in np.arange(halfWindow, lenData - halfWindow )]

    # final adaptive window
    dataAveraged[-halfWindow:] = [np.nansum(data[- 2*index -1 :])/( 2*index +1) \
                         for index in np.arange(halfWindow,0,-1)]

    # import matplotlib.pyplot as mplt
    # mplt.plot(data, 'bo')
    # mplt.plot(dataAveraged, 'r')
    # mplt.show()

    return dataAveraged

# ----------------------------------------------------------------------------------------------------------------------
def fnMovingAverageOld(data, N):
    '''
    :param data: data
    :param N: width of moving average window
    :return: moving average
    adapted from http://mathworld.wolfram.com/MovingAverage.html
    moving average with pixel centered on the moving average window

    '''
    dataAveraged = np.zeros_like(data)
    dataAveraged[int(N/2):-int(N/2)] = [np.nansum(data[index-int(N/2):index + int(N/2) -1])/(N-1) for index in np.arange(int(N/2), len(data) -int(N/2) )]
    return dataAveraged

# region --- Function fnNewGaussianFit
# ----------------------------------------------------------------------------------------------------------------------
def fnGaussianFitWithErr(data, XXXAxis, GaussParamsInitGuess=[0, 0, 1, 0]):
    '''
    '''
    resGauss = lambda params, resdata, resXXXAxis: fnGauss(resdata, params) - resXXXAxis

    data = np.array(data)

    GaussParams, cov, infodict, mesg, ier = opt.leastsq(resGauss, GaussParamsInitGuess, \
                                                        args=(data, XXXAxis), maxfev=1000000000, full_output=True)

    try:
        GaussParamsErrors = np.sqrt(np.diag(cov))
        # print GaussParamsErrors
    except:
        GaussParamsErrors  = (0.0,0,0)

    # print GaussParamsErrors

    return GaussParams, GaussParamsErrors

# endregion