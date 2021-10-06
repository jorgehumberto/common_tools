#!/usr/bin/python2.7
'''

main routine
'''

# ======================================================================================================================
#   Python Modules
# ======================================================================================================================
import os
import sys
# from astropy.io import fits
import matplotlib.pyplot as mplt

import scipy.stats as stats

# ======================================================================================================================
#   User Modules
# ======================================================================================================================
# Define Work Paths
WorkPath = os.path.abspath(os.getenv('WorkPath'))
ReductPath = os.path.abspath(os.getenv('ReductPath'))

from defaults.defaultPaths import *

from modules.fileManip import initPaths

from modules.InOut import fnPrintLine, fnGetYOrbit

from classes.classesCCFs import clsUVESCCF
from classes.classesPlanet import clsPlanetParametersOLD
from models.modelsOrbit import fnRVStarOrbitCircular

from math_local.mathFunctions import *

# ======================================================================================================================
#   Settings
# ======================================================================================================================



sys.path.append(WorkPath)
sys.path.append(ReductPath)
settings = __import__(sys.argv[1].split('.')[0])

scienceInputFolder = os.path.abspath('{}/CCFs/{}'.format(WorkPath, settings.ccfFolder))
scienceOuputFolder = os.path.abspath('{}/planetResults/{}/'.format(WorkPath, settings.ccfFolder))
initPaths([scienceOuputFolder])


# ======================================================================================================================
#   Functions
# ======================================================================================================================

def fnGetTemplateList(TemplateType, FitsList, CCFsCopy, PlanetRVs, PlanetCopy):
    if TemplateType == 'All':
        TemplateList = [CCFName for CCFName in FitsList]

    elif TemplateType == 'Opposition':
        TemplateList = [CCFName for CCFName in FitsList if \
                        (abs(PlanetRVs[CCFName] - CCFsCopy[CCFName].RVC) < 2. * PlanetCopy.FWHMMax) and \
                        CCFsCopy[CCFName].PlanetPhaseFolded > 0.5]

    elif TemplateType == 'Transit':
        TemplateList = [CCFName for CCFName in FitsList if \
                        (abs(PlanetRVs[CCFName] - CCFsCopy[CCFName].RVC) < 2. * PlanetCopy.FWHMMax) and \
                        CCFsCopy[CCFName].PlanetPhaseFolded < 0.5]

    elif TemplateType == 'Opp+Trans':
        TemplateList = [CCFName for CCFName in FitsList if \
                        (abs(PlanetRVs[CCFName] - CCFsCopy[CCFName].RVC) < 2. * PlanetCopy.FWHMMax)]
    else:
        print 'NO TEMPLATE HAS BEEN SPECIFIED, ABORTING!'
        sys.exit()
    return TemplateList


def fnShiftCCF(CCF, Wave, RV):
    """
    Function to shift a CCF of a given RV (in m/s)
    Spec - CCF to be shifted\
    Wave - array with the CCF RV data (l)
    """

    NewCCF = np.interp(Wave.copy() + RV, Wave.copy(), CCF.copy())

    return NewCCF.copy()


# =======================================================================
def fnBuildStarTemplate(CCFs, Templates=None, MeanRV=None):
    if Templates == None:
        Templates = CCFs.keys()

    TemplateRVMean = stats.nanmean([CCFs[CCFName].RVC for CCFName in CCFs.keys()])

    Template = np.nansum(
        [fnShiftCCF(CCFs[CCFName].data.copy(), CCFs[CCFName].wave.copy(), CCFs[CCFName].RVC - TemplateRVMean)
         for CCFName in Templates], axis=0)

    TemplateWave = CCFs[Templates[0]].wave - CCFs[Templates[0]].RVC

    Template = np.divide(Template.copy(), np.nanmax(Template.copy()))

    return Template, TemplateWave, TemplateRVMean


# =======================================================================
def fnBuildStarTemplateNew(CCFs, Templates=None):
    if Templates == None:
        Templates = CCFs.keys()

    templateMeanPixel = np.nanmean([CCFs[CCFName].ccfMeanPixel for CCFName in Templates])

    templateWave = np.linspace(0, len(CCFs[Templates[0]].data) - 1, num=len(CCFs[Templates[0]].data))

    templateData = np.nansum(
        [fnShiftCCF(CCFs[CCFName].data.copy(), templateWave, CCFs[CCFName].ccfMeanPixel - templateMeanPixel)
         for CCFName in Templates], axis=0)

    templateData = np.divide(templateData.copy(), np.nanmax(templateData.copy()))

    return templateData, templateWave, templateMeanPixel


# ======================================================================================================================
#   Main
# ======================================================================================================================

def main_function():
    print 'Recover planet CCFs'
    fnPrintLine('CCF', 'Extracting CCFs, please wait')

    fullCCFList = sorted(['{}/{}'.format(scienceInputFolder, fileName) for fileName in os.listdir(scienceInputFolder) if
                          fileName.endswith('REDR_ccf.fits')])

    fullCCFs = {ccfName: clsUVESCCF(ccfName) for ccfName in sorted(fullCCFList)}

    fnPrintLine('Config', 'Extracting planet parameters')
    planetParams = clsPlanetParametersOLD(fnGetYOrbit('{}/{}'.format(scienceInputFolder, settings.orbitParams)))

    print planetParams.__dict__
    # sys.exit()

    phaseZero = 0.5
    # planetRVs = {}

    FIX = True

    # =======================================================================
    fnPrintLine('CCF', 'getting planet RVS')
    rangeRVPixels = np.linspace(0, len(fullCCFs[fullCCFList[0]].data) - 1, len(fullCCFs[fullCCFList[0]].data))

    for ccfName in sorted(fullCCFList):
        if FIX == True:
            fullCCFs[ccfName].getPhase(planetParams.period, planetParams.t0)
            fullCCFs[ccfName].RV = fnRVStarOrbitCircular(planetParams, phaseZero, fullCCFs[ccfName].PlanetPhaseFolded)
            fullCCFs[ccfName].RVC = fullCCFs[ccfName].RV
        else:
            fullCCFs[ccfName].RV -= planetParams.Sysrv + .4
            fullCCFs[ccfName].RVC = fullCCFs[ccfName].RV

        # fullCCFs[ccfName].planetRV = -(fullCCFs[ccfName].RV )/planetParams.massratio
        # fullCCFs[ccfName].planetRVPixels

        fullCCFs[ccfName].ccfAmplitude, fullCCFs[ccfName].ccfMeanPixel, fullCCFs[ccfName].ccfFWHMPixels, fullCCFs[
            ccfName].ccfB = fnGaussianFitOLD(rangeRVPixels, fullCCFs[ccfName].data, GaussParamsInitGuess=[
            max(fullCCFs[ccfName].data) - min(fullCCFs[ccfName].data), len(fullCCFs[ccfName].data) / 2, 10.,
            max(fullCCFs[ccfName].data)])
        fullCCFs[ccfName].planetRV = -(fullCCFs[ccfName].RV) / planetParams.massratio
        fullCCFs[ccfName].planetRVMeanPixel = (fullCCFs[ccfName].planetRV - fullCCFs[ccfName].RVC) / fullCCFs[
            ccfName].CCFStep + fullCCFs[ccfName].ccfMeanPixel

        fullCCFs[ccfName].CCFWaveIni = fullCCFs[ccfName].RV - fullCCFs[ccfName].ccfMeanPixel * fullCCFs[ccfName].CCFStep

        fullCCFs[ccfName].wave = np.arange(len(fullCCFs[ccfName].data), dtype=np.float64) * fullCCFs[ccfName].CCFStep + \
                                 fullCCFs[ccfName].CCFWaveIni
        # print fullCCFs[ccfName].ccfMeanPixel, fullCCFs[ccfName].planetRVMeanPixel/

    # =======================================================================
    fnPrintLine('CCF', 'Building template')
    # starTemplateList = sorted(fullCCFList) #fnGetTemplateList('All', sorted(fullCCFList), fullCCFs, planetRVs.values(),planetParams)

    starTemplate, starTemplateWave, starTemplateMeanPixel = fnBuildStarTemplateNew(fullCCFs)

    templateA, templateMeanPixel, templateFWHM, templateB = fnGaussianFitOLD(rangeRVPixels, starTemplate,
                                                                             GaussParamsInitGuess=[
                                                                                 max(starTemplate) - min(starTemplate),
                                                                                 len(starTemplate) / 2, 10.,
                                                                                 max(starTemplate)])

    # =======================================================================
    planetCCF = []
    for ccfName in sorted(fullCCFList):

        diffPixels = (templateMeanPixel - fullCCFs[ccfName].ccfMeanPixel)  # * fullCCFs[ccfName].CCFStep
        shiftedTemplate = fnShiftCCF(starTemplate, starTemplateWave, diffPixels)

        fullCCFs[ccfName].normCCF(shiftedTemplate)

        # fullCCFs[ccfName].data[int(starMeanPixel-25):int(starMeanPixel+25)] = None
        # fullCCFs[ccfName].data = fnShiftCCF(fullCCFs[ccfName].data, fullCCFs[ccfName].wave, planetRVs[ccfName] )

        planetXXXRange = np.where((fullCCFs[ccfName].wave < fullCCFs[ccfName].planetRV + 40) & (
        fullCCFs[ccfName].wave > fullCCFs[ccfName].planetRV - 40))

        if fullCCFs[ccfName].planetRV - fullCCFs[ccfName].RV > 50:
            planetCCF.append(fullCCFs[ccfName].data[planetXXXRange][:80])
            mplt.plot(fullCCFs[ccfName].wave, fullCCFs[ccfName].data)
            # fullCCFs[ccfName].data[planetXXXRange] = max(fullCCFs[ccfName].data)

    mplt.savefig('{}/{}_normalisedCCFs.png'.format(scienceOuputFolder, settings.ccfFolder))
    mplt.clf()

    mplt.imshow(np.array([fullCCFs[ccfName].data for ccfName in sorted(fullCCFList)]), cmap='Greys')

    mplt.savefig('{}/{}_2D_normCCFs_noOrbit.png'.format(scienceOuputFolder, settings.ccfFolder))
    mplt.clf()
    iii = np.arange(len(fullCCFList))
    for ccfName, ii in zip(sorted(fullCCFList), iii):
        mplt.plot(fullCCFs[ccfName].planetRVMeanPixel, ii, 'ro')

    mplt.ylim([len(fullCCFList) - 1, 0])
    mplt.xlim([0, len(starTemplate) - 1])

    mplt.savefig('{}/{}_2D_normCCFs_withOrbit.png'.format(scienceOuputFolder, settings.ccfFolder))
    mplt.clf()

    finalCCF = np.nansum(planetCCF, axis=0)
    finalCCF = finalCCF / np.nanmedian(finalCCF) - 1
    mplt.plot(np.linspace(-50, 50, num=len(finalCCF)), finalCCF)
    mplt.axhline(1.e-5, color='r')
    mplt.axhline(-1.e-5, color='r')
    axFinalCCF = mplt.gca()
    axFinalCCF.ticklabel_format(useOffset=False)

    mplt.savefig('{}/{}_FinalPlanetCCFs.png'.format(scienceOuputFolder, settings.ccfFolder))
    mplt.clf()


# ======================================================================================================================
#   RunTime
# ======================================================================================================================
if __name__ == '__main__':
    main_function()
