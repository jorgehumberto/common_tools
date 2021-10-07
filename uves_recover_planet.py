#!/usr/bin/env python
'''

main routine
'''

# ======================================================================================================================
#   Initialization
# ======================================================================================================================

from initRecoverPlanet import *

# ======================================================================================================================
#   Settings
# ======================================================================================================================

# fitsFilter = 'REDR_ccf.fits'
fitsFilter = '.fits'
phaseZero = 0.5
FIX = True
planetHalfWidth = 50  # in km/s
distancePlanetStar = 60  # in km/s


# ======================================================================================================================
#   Main
# ======================================================================================================================

def main_function():
    dataFolderFileList = sorted(
        ['{}/{}'.format(scienceInputFolder, fileName) for fileName in os.listdir(scienceInputFolder) if
         fileName.endswith(fitsFilter)])

    try:
        fullCCFList = np.genfromtxt('{}/{}'.format(scienceInputFolder, settings.dataList), dtype=str)
        fullCCFList = ['{}/{}'.format(scienceInputFolder, fileName) for fileName in fullCCFList if
                       fileName.endswith(fitsFilter)]
        fnPrintLine('Config', 'data ccf list file provided: {}'.format(settings.dataList))
    except:
        fullCCFList = dataFolderFileList
        fnPrintLine('Config', 'no data ccf list file provided, using all ccfs.')

    try:
        templateCCFList = np.genfromtxt('{}/{}'.format(scienceInputFolder, settings.templateList), dtype=str)
        templateCCFList = ['{}/{}'.format(scienceInputFolder, fileName) for fileName in templateCCFList if
                           fileName.endswith(fitsFilter)]
        fnPrintLine('Config', 'template ccf list file provided: {}'.format(settings.templateList))
    except:
        templateCCFList = dataFolderFileList
        fnPrintLine('Config', 'no template ccf list file provided, using all ccfs.')

    fnPrintLine('CCF', 'Extracting CCFs, please wait')
    fullCCFs = {ccfName: fnOpenFits(ccfName) for ccfName in sorted(dataFolderFileList)}

    fnPrintLine('Config', 'Extracting planet parameters')
    planetParams = clsPlanetParameters(fnGetYOrbit('{}/{}'.format(scienceInputFolder, settings.orbitParams)))

    # ======================================================================================================================
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


    # ======================================================================================================================
    fnPrintLine('CCF', 'Building template')
    # REDR chip
    starTemplateREDR = fnBuildStarTemplate(fullCCFs,
                                           Templates=[ccfName for ccfName in templateCCFList if 'REDR' in ccfName])

    # REDL side
    starTemplateREDL = fnBuildStarTemplate(fullCCFs,
                                           Templates=[ccfName for ccfName in templateCCFList if 'REDL' in ccfName])

    # ======================================================================================================================
    planetCCFREDR, planetCCFREDL = [], []

    for ccfName in sorted(fullCCFList):
        if 'REDR' in ccfName:
            templateMeanPixel = starTemplateREDR.ccfMeanPixel
            starTemplate = starTemplateREDR.data
            starTemplateWave = starTemplateREDR.wave
        elif 'REDL' in ccfName:
            templateMeanPixel = starTemplateREDL.ccfMeanPixel
            starTemplate = starTemplateREDL.data
            starTemplateWave = starTemplateREDL.wave

        diffPixels = (templateMeanPixel - fullCCFs[ccfName].ccfMeanPixel)
        shiftedTemplate = fnShiftCCF(starTemplate, starTemplateWave, diffPixels)

        fullCCFs[ccfName].normCCF(shiftedTemplate)

        planetXXXRange = np.where((fullCCFs[ccfName].wave < fullCCFs[ccfName].planetRV + planetHalfWidth) & (
        fullCCFs[ccfName].wave > fullCCFs[ccfName].planetRV - planetHalfWidth))

        if abs(fullCCFs[ccfName].planetRV - fullCCFs[ccfName].RV) > distancePlanetStar:
            if 'REDR' in ccfName:
                planetCCFREDR.append(fullCCFs[ccfName].data[planetXXXRange][:int(distancePlanetStar * 2)])
            elif 'REDL' in ccfName:
                planetCCFREDL.append(fullCCFs[ccfName].data[planetXXXRange][:int(distancePlanetStar * 2)])




            # ======================================================================================================================
    fnPrintLine('CCF', 'Generating Plots')
    ccfListREDR = [ccfName for ccfName in fullCCFList if 'REDR' in ccfName]
    ccfListREDL = [ccfName for ccfName in fullCCFList if 'REDL' in ccfName]
    print len(ccfListREDR), len(ccfListREDL)
    # 1D- normalised CCFS
    figNormCCFs1D, (axNormCCFs1DR, axNormCCFs1DL) = mplt.subplots(2, 1, sharex=True)
    axNormCCFs1DR.set_title('REDR')
    axNormCCFs1DL.set_title('REDL')
    for ccfName in ccfListREDR:
        axNormCCFs1DR.plot(fullCCFs[ccfName].wave, fullCCFs[ccfName].data)
    for ccfName in ccfListREDL:
        axNormCCFs1DL.plot(fullCCFs[ccfName].wave, fullCCFs[ccfName].data)

    axNormCCFs1DR.set_xlim([-300, 300])

    figNormCCFs1D.savefig('{}/{}_1DnormalisedCCFs.png'.format(scienceOuputFolder, settings.dataFolder))

    # 2D- normalised CCFS - no orbit
    figNormCCFs2D, (axNormCCFs2DR, axNormCCFs2DL) = mplt.subplots(2, 1, sharex=True)
    axNormCCFs2DR.set_title('REDR')
    axNormCCFs2DL.set_title('REDL')
    axNormCCFs2DR.imshow(np.array([fullCCFs[ccfName].data for ccfName in sorted(ccfListREDR)]), cmap='Greys',
                         extent=(-300, 300, 0, len(ccfListREDR) - 1),
                         aspect='auto')
    axNormCCFs2DL.imshow(np.array([fullCCFs[ccfName].data for ccfName in sorted(ccfListREDL)]), cmap='Greys',
                         extent=(-300, 300, 0, len(ccfListREDL) - 1),
                         aspect='auto')

    axNormCCFs2DR.set_ylim([0, len(ccfListREDR) - 1])
    # axNormCCFs2DR.set_ylim([0, 10])
    axNormCCFs2DR.set_xlim([-300, 300])
    axNormCCFs2DL.set_ylim([0, len(ccfListREDL) - 1])
    # axNormCCFs2DL.set_ylim([0, 200])

    figNormCCFs2D.savefig('{}/{}_2D_normCCFs_noOrbit.png'.format(scienceOuputFolder, settings.dataFolder))

    # 2D- normalised CCFS - with orbit

    for ccfName, ii in zip(ccfListREDR, np.arange(len(ccfListREDR))):
        axNormCCFs2DR.plot(fullCCFs[ccfName].planetRV, ii, 'ro')
    for ccfName, ii in zip(ccfListREDL, np.arange(len(ccfListREDL))):
        axNormCCFs2DL.plot(fullCCFs[ccfName].planetRV, ii, 'ro')

    axNormCCFs2DR.set_ylim([0, len(ccfListREDR) - 1])
    axNormCCFs2DL.set_ylim([0, len(ccfListREDL) - 1])
    # axNormCCFs2DR.set_ylim([0, 10])
    axNormCCFs2DR.set_xlim([-300, 300])
    # axNormCCFs2DL.set_ylim([0, 10])
    # axNormCCFs2DL.set_xlim([0,len(starTemplate)-1])

    mplt.savefig('{}/{}_2D_normCCFs_withOrbit.png'.format(scienceOuputFolder, settings.dataFolder))

    # recovered planet CCF
    figPlanetCCFs1D, (axPlanetCCFs1DR, axPlanetCCFs1DL) = mplt.subplots(2, 1, sharex=True)
    axPlanetCCFs1DR.set_title('REDR')
    axPlanetCCFs1DL.set_title('REDL')

    finalCCFREDR = np.nansum(planetCCFREDR, axis=0)
    finalCCFREDR = finalCCFREDR / np.nanmedian(finalCCFREDR) - 1

    finalCCFREDL = np.nansum(planetCCFREDL, axis=0)
    finalCCFREDL = finalCCFREDL / np.nanmedian(finalCCFREDL) - 1

    axPlanetCCFs1DL.plot(np.linspace(-planetHalfWidth, planetHalfWidth, num=len(finalCCFREDR)), finalCCFREDR)
    axPlanetCCFs1DR.plot(np.linspace(-planetHalfWidth, planetHalfWidth, num=len(finalCCFREDL)), finalCCFREDL)

    axPlanetCCFs1DR.axhline(1.e-5, color='r')
    axPlanetCCFs1DR.axhline(-1.e-5, color='r')

    axPlanetCCFs1DL.axhline(1.e-5, color='r')
    axPlanetCCFs1DL.axhline(-1.e-5, color='r')

    axPlanetCCFs1DR.ticklabel_format(useOffset=False)
    axPlanetCCFs1DR.ticklabel_format(useOffset=False)

    figPlanetCCFs1D.savefig('{}/{}_FinalPlanetCCFs.png'.format(scienceOuputFolder, settings.dataFolder))


# ======================================================================================================================
#   RunTime
# ======================================================================================================================
if __name__ == '__main__':
    main_function()
