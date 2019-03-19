# -*- coding: utf-8 -*-
"""
Created on Tue Set 08 2015

@author: Jorge Martins
email: jorge.martins@iastro.pt


"""
###############################################################################
# Python Modules
###############################################################################

import numpy as np
from astropy.io import fits
from classes.classesCCFs import clsUVESCCF, clsHARPSCCF
from math_local.mathFunctions import fnGaussianFitOLD
from modules.InOut import fnPrintLine

###############################################################################
# Functions
###############################################################################
def fnOpenUVESFits(FitsFileName):
    ccfFits = fits.open(FitsFileName)[0].copy()

    ccfFile = clsUVESCCF(ccfFits.data[0], ccfFits.header)

    return ccfFile


def fnOpenHARPSFits(FitsFileName, headerOnly=False, allOrders = False):
    if headerOnly == True:
        data = []
        header = fits.getheader(FitsFileName)

    else:
        ccfFits = fits.open(FitsFileName)[0].copy()
        if allOrders:
            data = np.array(ccfFits.data).copy()
        else:
            data = np.array(ccfFits.data)[-1].copy()

        header = ccfFits.header

    ccfFile = clsHARPSCCF(data, header)

    return ccfFile


def fnOpenHARPSFits2(FitsFileName):
    ccfFits = fits.open(FitsFileName)

    ccfClass = clsHARPSCCF(ccfFits[0].data, ccfFits[0].header)

    return ccfClass


def fnCCFunction(lineFraction, lineDepth, fluxPixel):
    '''
    '''
    CCFunction = lineFraction * lineDepth * fluxPixel
    return CCFunction


# ======================================================================================================================


def fnCCFunction(lineFraction, lineDepth, fluxPixel):
    '''
    '''
    CCFunction = lineFraction * lineDepth * fluxPixel
    return CCFunction


# ======================================================================================================================


def fnCCFunction(lineFraction, lineDepth, fluxPixel):
    '''
    '''
    CCFunction = lineFraction * lineDepth * fluxPixel
    return CCFunction


# ======================================================================================================================
def fnRVLambdaShift_OLD(RV, lambdaRest):
    ''' Function to compute a wavelenght shift due to radial velocity
        using RV / c = lambdaShift/lambdaRest
        RV - radial velocity (in m/s)
        lambdaRest - rest wavelenght of the spectral line
        lambdaShift - lambdaFinal - lambdaRest
        '''
    c = 299792.458
    lambdaShift = lambdaRest * (RV / c)

    lambdaFinal = lambdaRest + lambdaShift

    return lambdaFinal


# ======================================================================================================================

def fnRVLambdaShift(lambdaRest, RV):
    ''' Function to compute a wavelenght shift due to radial velocity
        using RV / c = lambdaShift/lambdaRest
        RV - radial velocity (in m/s)
        lambdaRest - rest wavelenght of the spectral line
        lambdaShift - lambdaFinal - lambdaRest
        '''
    c = 299792.458
    lambdaShift = lambdaRest * (RV / c)

    lambdaFinal = lambdaRest + lambdaShift

    return lambdaFinal


# ======================================================================================================================
def fnComputeCCF_OLD(fluxOrder, waveOrder, maskOrder, RV, BERV=0.0):
    '''
    fluxOrder    -   flux of spectrum
    waveOrder    -   wavelength solution
    maskOrder    -   binary mask for the CCF
    '''
    fnPrintLine(None, 'Computing CCF for BERV corrected RV = {:.2f}km/s'.format(RV), sameLine=True)
    waveResolution = waveOrder[1] - waveOrder[0]

    fluxOrderPlusOne = np.roll(fluxOrder, -1)

    waveorderPlusOne = np.roll(waveOrder, -1)

    nLines = 0
    lines = []
    CCF = 0.0
    for lineIni, lineEnd, lineDepth in maskOrder:

        lineIniShifted = fnRVLambdaShift(RV, lineIni)

        lineEndShifted = fnRVLambdaShift(RV, lineEnd)

        if waveOrder[0] < lineIniShifted and lineEndShifted < waveOrder[-1]:
            linePixelIni, linePixelEnd = [np.where(waveOrder > lineIniShifted)[0][0],
                                          np.where(waveOrder < lineEndShifted)[0][-1]]

            lineFractionIni = (waveOrder[linePixelIni] - lineIniShifted) / waveResolution
            lineFractionEnd = (lineEndShifted - waveOrder[linePixelEnd]) / waveResolution

            CCF += lineDepth * (
            np.nansum(fluxOrder[linePixelIni:linePixelEnd]) + (lineFractionIni) * fluxOrder[linePixelIni - 1] + (
            lineFractionEnd) * fluxOrder[linePixelEnd + 1])
            nLines += 1
            lines.append(lineIniShifted)

    return CCF, nLines


# ======================================================================================================================
def fnComputeCCF(fluxOrder, waveOrder, maskOrder, RV, BERV=0.0):
    '''
    fluxOrder    -   flux of spectrum
    waveOrder    -   wavelength solution
    maskOrder    -   binary mask for the CCF
    '''
    fnPrintLine(None, 'Computing CCF for BERV corrected RV = {:.2f}km/s'.format(RV), sameLine=True)
    waveResolution = waveOrder[1] - waveOrder[0]

    fluxOrderPlusOne = np.roll(fluxOrder, -1)

    waveorderPlusOne = np.roll(waveOrder, -1)

    nLines = 0
    lines = []
    CCF = 0.0

    maskOrderShifted = maskOrder.copy()
    maskOrderShifted[:, 0] = fnRVLambdaShift(maskOrder[:, 0], RV)
    maskOrderShifted[:, 1] = fnRVLambdaShift(maskOrder[:, 1], RV)

    maskOrderShiftedOrder = maskOrderShifted[
        (waveOrder[0] < maskOrderShifted[:, 0]) & (maskOrderShifted[:, 1] < waveOrder[-1])]

    for lineIniShifted, lineEndShifted, lineDepth in maskOrderShiftedOrder:
        linePixelIni, linePixelEnd = [np.where(waveOrder > lineIniShifted)[0][0],
                                      np.where(waveOrder < lineEndShifted)[0][-1]]

        lineFractionIni = (waveOrder[linePixelIni] - lineIniShifted) / waveResolution
        lineFractionEnd = (lineEndShifted - waveOrder[linePixelEnd]) / waveResolution

        CCF += lineDepth * (
            np.nansum(fluxOrder[linePixelIni:linePixelEnd]) + (lineFractionIni) * fluxOrder[linePixelIni - 1] + (
                lineFractionEnd) * fluxOrder[linePixelEnd + 1])
        nLines += 1
        lines.append(lineIniShifted)

    return CCF, nLines


# ======================================================================================================================

def fnShiftCCF(CCF, Wave, RV):
    """
    Function to shift a CCF of a given RV
    CCF - CCF to be shifted
    Wave - array with the CCF RV data (l)

    NOTE: the units of RV must match those of Wave

    """

    NewCCF = np.interp(Wave.copy() + RV, Wave.copy(), CCF.copy())

    return NewCCF.copy()


# ======================================================================================================================

def fnBuildStarTemplateOLD(CCFs, Templates=None):
    if Templates == None:
        Templates = CCFs.keys()

    templateMeanPixel = np.nanmean([CCFs[CCFName].ccfMeanPixel for CCFName in Templates])

    templateWave = np.linspace(0, len(CCFs[Templates[0]].data) - 1, num=len(CCFs[Templates[0]].data))

    templateData = np.nansum(
        [fnShiftCCF(CCFs[CCFName].data.copy(), templateWave, CCFs[CCFName].ccfMeanPixel - templateMeanPixel)
         for CCFName in Templates], axis=0)

    templateData = np.divide(templateData.copy(), np.nanmax(templateData.copy()))

    templateCCF = clsUVESCCF(templateData)
    templateCCF.wave = templateWave
    templateCCF.ccfMeanPixel = templateMeanPixel

    templateCCF.ccfAmplitude, templateCCF.ccfMeanPixel, templateCCF.ccfFWHMPixels, templateCCF.ccfB = fnGaussianFitOLD(
        templateCCF.wave, templateCCF.data,
        GaussParamsInitGuess=[max(templateCCF.data) - min(templateCCF.data), len(templateCCF.data) / 2, 10.,
                              max(templateCCF.data)])
    # return templateData, templateWave, templateMeanPixel
    return templateCCF

# ======================================================================================================================

def fnBuildStarTemplate(CCFs, nOrders = 1,Templates=None):
    if Templates == None:
        Templates = sorted(CCFs.keys())

    # defining template CCF class
    # templateCCF = clsHARPSCCF(CCFs[Templates[0]].data)
    templateCCF = clsFits_CCF(CCFs[Templates[0]].data)
    nOrders = len(templateCCF.data)
    templateCCF.wave = np.indices(np.shape(templateCCF.data))[1,:,:]
    templateCCF.ccfAmplitude, templateCCF.ccfMeanPixel, templateCCF.ccfFWHMPixels, templateCCF.ccfB = (np.empty(nOrders),np.empty(nOrders),np.empty(nOrders),np.empty(nOrders))

    # building template orders
    for order, orderWave, orderData in zip(range(nOrders), templateCCF.wave, templateCCF.data):
        meanPixel = np.nanmean([CCFs[CCFName].ccfMeanPixel[order] for CCFName in Templates])
        orderData[:] = np.nansum([fnShiftCCF(CCFs[CCFName].data[order].copy(), orderWave, CCFs[CCFName].ccfMeanPixel[order] - meanPixel)
                                   for CCFName in Templates], axis=0)

        templateCCF.ccfAmplitude[order], templateCCF.ccfMeanPixel[order], templateCCF.ccfFWHMPixels[order], templateCCF.ccfB[order] = fnGaussianFitOLD(
            orderWave, orderData,
        GaussParamsInitGuess=[max(orderData) - min(orderData), len(orderData) / 2, 10.,
                              max(orderData)])

    return templateCCF
#----------------------------------------------------------------------------------------------------------------------

def fnCreateWave(orderHeader, order, numberOrders=1):
    '''
    :param orderHeader:
    :param order:
    :param numberOrders:
    :return:
    '''
    if numberOrders == 1:
        orderWavelengthStart = orderHeader['CRVAL1']
        waveData = np.linspace(orderHeader['CRVAL1'],
                               orderHeader['CRVAL1'] + (orderHeader['CDELT1'] * orderHeader['NAXIS1']),
                               num=orderHeader['NAXIS1'])
    else:
        orderWavelengthLenght = orderHeader['NAXIS1']
        orderWavelengthStart = orderHeader['WSTART{0}'.format(order + 1)]
        orderWavelengthEnd = orderHeader['WEND{0}'.format(order + 1)]
        waveData = np.linspace(orderWavelengthStart, orderWavelengthEnd, num=orderWavelengthLenght)

    return waveData


def fnCreateWave2(orderHeader, order):
    '''
    :param orderHeader:
    :param order:
    :return:
    '''
    orderWavelengthStart = orderHeader['WSTART{0}'.format(order + 1)]
    # orderWavelengthEnd       =  orderHeader['WEND{0}'.format(order+1)]
    orderStep = orderHeader['CDELT1']
    orderLen = orderHeader['NAXIS1']
    orderWavelengthEnd = orderWavelengthStart + orderStep * orderLen
    waveData = np.linspace(orderWavelengthStart, orderWavelengthEnd, num=orderLen)

    return waveData
