#!/usr/bin/env python
'''



'''

__author__ = 'Jorge Martins'
__email__ = 'jorge.martins@iastro.pt'
__version__ = 'alpha 1'
__date__ = '2015/12/02'

# region --- INITIALIZATION
# region --- Python Modules
import ConfigParser
import os
import sys

import matplotlib.pyplot as mplt
import numpy as np
from PyAstronomy.pyTiming import pyPeriod
from PyAstronomy.pyasl.asl.astroTimeLegacy import daycnv
from matplotlib import cm
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter

os.system('cls' if os.name == 'nt' else 'clear')
# endregion

# region --- User Modules
from modules.InOut import fnPrintLine
from classes.classesPlanet import clsPlanetParametersOLD
from recipes.manipCCF import fnOpenHARPSFits, fnBuildStarTemplateOLD, fnShiftCCF
from math_local.mathFunctions import fnGaussianFitOLD, fnGauss
from models.modelsOrbit import fnRVStarOrbitElipse, fnPhaseFunction
from modules.plots import fnPlot2DCCFs
import ephem

# endregion

# region --- Parse Settings
cfgFileName = sys.argv[1]
cfgFile = ConfigParser.RawConfigParser(allow_no_value=True)
cfgFile.optionxform = str
cfgFile.read(cfgFileName)
# endregion

# region --- path, folder and filenames definition
WorkPath = os.path.abspath(os.getenv('WorkPath'))

scienceInputFolder = os.path.join(*[WorkPath, 'CCFs', cfgFile.get('global', 'dataFolder')])

sys.path.append(WorkPath)

resultsFolder = os.path.join(*[WorkPath, 'results', cfgFile.get('global', 'resultsFolder')])
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)




# endregion

# region --- generic settings
starWidth = float(cfgFile.get('detection_settings', 'starWidth'))
# endregion

# endregion


# ======================================================================================================================
#   Main
# ======================================================================================================================

def __main__():
    # region --- Define planet parameters
    planetParams = clsPlanetParametersOLD(params={param: value for param, value in cfgFile.items('orbital_params')})
    # endregion

    # region --- Get files list
    # data files
    try:
        # if a list with fits files from which to recover the planet CCF has been provided
        dataCCFList = sorted([os.path.join(scienceInputFolder, fileName) \
                              for fileName in list(
                np.genfromtxt(os.path.join(scienceInputFolder, cfgFile.get('detection_settings', 'dataFilelist')), \
                              dtype=str))])

        fnPrintLine('Config', 'data ccf list file provided: {}'.format( \
            os.path.join(scienceInputFolder, cfgFile.get('detection_settings', 'dataFilelist'))))
    except:
        dataCCFList = sorted([os.path.join(scienceInputFolder, fileName) for fileName in os.listdir(scienceInputFolder) \
                              if fileName.lower().endswith('fits')])

        fnPrintLine('Config', 'no data ccf list file provided, using all ccfs.')

    # template files
    try:
        # if a list with fits files from which to construct the template has been provided
        templateCCFList = sorted([os.path.join(scienceInputFolder, fileName) \
                                  for fileName in list(
                np.genfromtxt(os.path.join(scienceInputFolder, cfgFile.get('detection_settings', 'templateFilelist')), \
                              dtype=str))])

        fnPrintLine('Config', 'template ccf list file provided: {}'.format( \
            os.path.join(scienceInputFolder, cfgFile.get('detection_settings', 'templateFilelist'))))
    except:
        templateCCFList = sorted(
            [os.path.join(scienceInputFolder, fileName) for fileName in os.listdir(scienceInputFolder) \
             if fileName.lower().endswith('fits')])
        fnPrintLine('Config', 'no template ccf list file provided, using all ccfs.')

    # all files
    fullCCFList = sorted(set(dataCCFList + templateCCFList))

    # endregion

    # region --- load fits files
    fnPrintLine('CCF', 'Extracting CCFs, please wait')
    fullCCFs = {ccfName: fnOpenHARPSFits(ccfName) for ccfName in sorted(fullCCFList)}

    # constrain file lists
    # try:
    #     dataCCFList = [ccfName for ccfName in dataCCFList \
    #                    if fullCCFs[ccfName].SN50 >= float(cfgFile.get('detection_settings', 'minSN50')) \
    #                    ]
    #
    #     templateCCFList = [ccfName for ccfName in templateCCFList \
    #                        if fullCCFs[ccfName].SN50 >= float(cfgFile.get('detection_settings', 'minSN50')) \
    #                        ]
    # except:
    #     pass

    # endregion

    # region --- Finding star CCF parameters
    fnPrintLine('CCF', 'finding star CCF parameters')
    # define xxx axis
    rangeRVPixels = np.linspace(0, len(fullCCFs[fullCCFList[0]].data) - 1, len(fullCCFs[fullCCFList[0]].data))

    ccfPixelTrimMargins = 10.  # in pixels

    for ccfName in sorted(fullCCFList):
        fullCCFs[ccfName].ccfAmplitude, fullCCFs[ccfName].ccfMeanPixel, fullCCFs[ccfName].ccfFWHMPixels, fullCCFs[
            ccfName].ccfB = \
            fnGaussianFitOLD(rangeRVPixels, fullCCFs[ccfName].data, \
                             GaussParamsInitGuess=[max(fullCCFs[ccfName].data) - min(fullCCFs[ccfName].data), \
                                                   len(fullCCFs[ccfName].data) / 2, \
                                                   10., \
                                                   max(fullCCFs[ccfName].data)])

        fullCCFs[ccfName].PlanetPhase = (fullCCFs[ccfName].BJD - planetParams.t0) / planetParams.period + np.radians(
            planetParams.w) / (2 * np.pi)
        fullCCFs[ccfName].PlanetPhaseFolded = (fullCCFs[ccfName].PlanetPhase) % 1
        fullCCFs[ccfName].RV = fullCCFs[ccfName].RVC = fnRVStarOrbitElipse(planetParams, fullCCFs[
            ccfName].PlanetPhaseFolded - np.radians(planetParams.w) / (2 * np.pi))

        fullCCFs[ccfName].planetRV = -fullCCFs[ccfName].RV / planetParams.massRatio

        fullCCFs[ccfName].wave = np.arange(-fullCCFs[ccfName].ccfMeanPixel, \
                                           len(fullCCFs[ccfName].data) - fullCCFs[ccfName].ccfMeanPixel) * \
                                 fullCCFs[ccfName].CCFStep + \
                                 fullCCFs[ccfName].RV

    dateList = set([fullCCFs[ccfName].obsDate for ccfName in sorted(fullCCFList)])

    for date in dateList:
        mplt.plot([fullCCFs[ccfName].PlanetPhaseFolded for ccfName in sorted(fullCCFList)\
                   if date in fullCCFs[ccfName].obsDate], \
                  [fullCCFs[ccfName].contrast for ccfName in sorted(fullCCFList)\
                   if date in fullCCFs[ccfName].obsDate], 'o')
    mplt.show()
    sys.exit()

    # endregion

    # region --- build template

    fnPrintLine('CCF', 'build star template')

    fullCCFs['template'] = fnBuildStarTemplateOLD(fullCCFs, Templates=templateCCFList)
    fullCCFs['template'].CCFStep = fullCCFs[ccfName].CCFStep

    fullCCFs['template'].CCFWaveIni = - fullCCFs['template'].CCFStep * fullCCFs['template'].ccfMeanPixel
    fullCCFs['template'].BJD = 0.0
    fullCCFs['template'].PlanetPhaseFolded = 0.0

    # endregion

    # region --- normalize CCFs by template
    fnPrintLine('CCF', 'normalise by star template')
    dataSelectedCCFList = [ccfName for ccfName in sorted(dataCCFList) if
                           abs(fullCCFs[ccfName].planetRV - fullCCFs[ccfName].RV) > float(
                               cfgFile.get('detection_settings', 'distancePlanetStar'))]

    for ccfName in sorted(fullCCFList):
        templateMeanPixel = fullCCFs['template'].ccfMeanPixel
        starTemplate = fullCCFs['template'].data
        starTemplateWave = fullCCFs['template'].wave

        diffPixels = (templateMeanPixel - fullCCFs[ccfName].ccfMeanPixel)
        shiftedTemplate = fnShiftCCF(starTemplate, starTemplateWave, diffPixels)

        fullCCFs[ccfName].normCCF(shiftedTemplate)

    # endregion

    # region --- searching for planet
    fnPrintLine('CCF', 'searching for planet, please wait')
    rangeRVsPlanet = {}
    planetWindowHalfWidth = float(cfgFile.get('detection_settings', 'planetHalfWidth'))
    for ccfName in sorted(dataCCFList):
        rangeRVsPlanet[ccfName] = np.where( \
            (fullCCFs[ccfName].wave < fullCCFs[ccfName].planetRV + planetWindowHalfWidth) & \
            (fullCCFs[ccfName].wave > fullCCFs[ccfName].planetRV - planetWindowHalfWidth) \
            )

        #fullCCFs[ccfName].planetRVRange = fullCCFs[ccfName].data[rangeRVsPlanet[ccfName]]
        # print fullCCFs[ccfName].planetRVRange[0], len(fullCCFs[ccfName].planetRVRange)

    widthPlanet = np.nanmin([len(rangeRVsPlanet[ccfName][0]) for ccfName in sorted(dataSelectedCCFList)])

    planetCCF = np.nansum(
        [fullCCFs[ccfName].data[rangeRVsPlanet[ccfName]][:widthPlanet] for ccfName in sorted(dataSelectedCCFList)],
        axis=0)

    planetWave = np.linspace(- planetWindowHalfWidth, planetWindowHalfWidth, num=len(planetCCF))

    # linear fit - y = m x + c
    m, c = np.polyfit(planetWave, planetCCF, 1)

    YYY = [m * x + c for x in planetWave]
    planetCCF = planetCCF / YYY
    planetCCF /= np.nanmedian(planetCCF)

    # planet CCF fit
    ccfAmplitude, ccfMean, ccfFWHM, ccfB = fnGaussianFitOLD(planetWave, planetCCF, \
                                                            GaussParamsInitGuess=[max(planetCCF) - min(planetCCF), \
                                                                                  0.0, \
                                                                                  10., \
                                                                                  np.nanmedian(planetCCF)]
                                                            )

    planetFit = fnGauss(planetWave, [ccfAmplitude, ccfMean, ccfFWHM, ccfB])

    # endregion

    # region --- correlating different nights
    fnPrintLine('CCF', 'correlation!!!')

    # delete central region of star CCF for better contrast
    # for ccfName in sorted(fullCCFList):


    nightList = sorted(set([int(fullCCFs[ccfName].BJD) for ccfName in  dataCCFList]))

    ccfListNight = {night:[ccfName for ccfName in sorted(fullCCFList) \
                      if int(fullCCFs[ccfName].BJD) == night] for night in nightList}

    nightCorrelation = {}
    nightTemplates = {night:fnBuildStarTemplateOLD(fullCCFs, Templates=ccfListNight[night]) for night in nightList}

    print len(nightList)

    # remove center of CCF
    for night in nightList:
        nightTemplates[night].data[fullCCFs[ccfName].ccfMeanPixel - starWidth * fullCCFs[ccfName].ccfFWHMPixels:\
            fullCCFs[ccfName].ccfMeanPixel + starWidth *fullCCFs[ccfName].ccfFWHMPixels] = None

        # mplt.plot(nightTemplates[night].data)

    # do correlation
    # for nightFirst, nightSecond in zip(nightList[:-1], nightList[1:] ):
    #     print nightFirst, nightSecond
    #     # # ratio = nightTemplates[nightFirst].data/nightTemplates[nightList[0]].data
    #     # # mplt.plot(ratio)
    #     # # mplt.plot()
    #     # # mplt.show()
    #     # correlation = []#np.correlate(nightTemplates[nightFirst].data, nightTemplates[nightSecond].data,mode='full')
    #     # pixRVmax = 2000
    #     # for pixRV in waveRange(0,pixRVmax):
    #     #     # print pixRV
    #     #     correlation.append(np.corrcoef(nightTemplates[nightFirst].data[pixRV:pixRV-pixRVmax], nightTemplates[nightSecond].data[pixRVmax:])[0,1])
    #     #     print np.corrcoef(nightTemplates[nightFirst].data[pixRV:pixRV-pixRVmax], nightTemplates[nightSecond].data[pixRVmax:])
    #     #     print
    #
    #     print np.corrcoef(nightTemplates[nightFirst].data, nightTemplates[nightSecond].data)
    #
    #     mplt.plot(nightTemplates[nightFirst].data)
    #     mplt.plot(nightTemplates[nightSecond].data)
    #
    #     # mplt.plot(correlation, 'ro')
    #     mplt.show()





    # mplt.show()
    # sys.exit()


    # endregion


    # region --- SN and phase function
    fnPrintLine('CCF', 'SN and phase function')
    expectedSN = np.sqrt(np.nansum(
        [fullCCFs[ccfName].nLines * (fullCCFs[ccfName].SN50 ** 2) for ccfName in sorted(dataSelectedCCFList)]))
    measuredNoise = np.nanstd(planetCCF)

    phaseFunction = {
    ccfName: fnPhaseFunction(planetParams.I, fullCCFs[ccfName].PlanetPhaseFolded - .25, model='Lambert') \
    for ccfName in dataSelectedCCFList}

    medianPhaseFunction = np.nanmedian(phaseFunction.values())
    medianStarCCFContrast = np.nanmedian([fullCCFs[ccfName].contrast for ccfName in dataSelectedCCFList])

    albedoLimitNoise = (planetParams.a / planetParams.radiusPlanet) ** 2 * 1 / (medianPhaseFunction) * (measuredNoise)
    albedoLimitCCF = (planetParams.a / planetParams.radiusPlanet) ** 2 * 1 / (medianPhaseFunction) * (
    ccfAmplitude / medianStarCCFContrast)

    # endregion

    # region --- PERIODOGRAMS
    # fnPrintLine('CCF', 'Periodograms of normalised CCFs')
    # # http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyTimingDoc/pyPeriodDoc/examples.html
    #
    # # periodogram of template
    # fnPrintLine('CCF', 'Periodograms of star template')
    # starMask = np.ones(len(fullCCFs['template'].data), np.bool)
    # starMask[fullCCFs['template'].ccfMeanPixel - 3 * fullCCFs['template'].ccfFWHMPixels: \
    #     fullCCFs['template'].ccfMeanPixel + 3 * fullCCFs['template'].ccfFWHMPixels] = 0
    # fullCCFs['template'].wave = fullCCFs['template'].wave.copy()[starMask]
    # fullCCFs['template'].wave = fullCCFs['template'].wave * fullCCFs['template'].CCFStep + fullCCFs['template'].CCFWaveIni
    #
    #
    # fullCCFs['template'].data = fullCCFs['template'].data.copy()[starMask]
    #
    # periodogramTemplateTimeSeries = pyPeriod.TimeSeries(fullCCFs['template'].wave, fullCCFs['template'].data)
    # periodogramTemplate = pyPeriod.LombScargle(periodogramTemplateTimeSeries, ofac=1, hifac=1)
    #
    #
    #
    # # computes the Lomb Scargle Periodogram of the RV using each frequency as a guess for the planet CCF
    # fnPrintLine('CCF', 'Periodograms of planet CCF')
    #
    # periodogramPlanetCCFTimeSeries = pyPeriod.TimeSeries(planetWave, planetCCF)
    # periodogramPlanetCCF = pyPeriod.LombScargle(periodogramPlanetCCFTimeSeries, ofac=1, hifac=1)
    #
    #
    # powerPeaks = [float(peak) for peak in cfgFile.get('periodograms', 'powerPeaks').strip('[]').split(',')]
    # print powerPeaks
    #
    # for power in powerPeaks:
    #     print power, periodogramPlanetCCF.stats(power)
    #
    # # periodogram of RVs
    # fnPrintLine('CCF', 'Periodograms of RVs')
    #
    # periodogramRVTimeSeries = pyPeriod.TimeSeries(np.array([fullCCFs[ccfName].BJD for ccfName in fullCCFList]), \
    #                                               np.array([fullCCFs[ccfName].RV for ccfName in fullCCFList]))
    # periodogramRV = pyPeriod.LombScargle(periodogramRVTimeSeries, ofac=1, hifac=1)

    # endregion

    # region --- [PLOTS]
    pdfFile = PdfPages(os.path.join(resultsFolder, 'figures_{}.pdf'.format(cfgFile.get('global', 'resultsFolder'))))

    # region [PLOT] orbit points
    fnPrintLine('PLOT', 'orbit points')
    figPhases, axPhases = mplt.subplots(2, 1, sharex=True, figsize=(12, 8))

    phasesFullOrbit = np.arange(0., 1., .01)
    rvFullOrbit = [
        -fnRVStarOrbitElipse(planetParams, phase - np.radians(planetParams.w) / (2 * np.pi)) / planetParams.massRatio
        for phase in phasesFullOrbit]

    for ax in axPhases:
        ax.plot([fullCCFs[ccfName].PlanetPhaseFolded for ccfName in sorted(fullCCFList)], \
                [fullCCFs[ccfName].planetRV for ccfName in sorted(fullCCFList)], 'bd', markersize=8)
        ax.grid(b=True, which='major', color='k', linestyle='--', alpha=.3)
        # ax.set_ylim(-1.1*planetParams.k2,1.1*planetParams.k2)
        # ax.axhline(0.0, label= r'$RV_{CM}$')
        ax.axvline(0.25, label='transit', color='red')
        ax.axvline(0.75, label='opposition', color='green')
        ax.legend(loc='best', title='r$\phi = 0$ - Ascending node')
        ax.set_xlim(0, 1)
        ax.set_xticks(np.arange(0, 1, .05))
        ax.plot(phasesFullOrbit, rvFullOrbit, 'b-', alpha=.3, markersize=8)

    # ax template
    axPhases[0].set_title('template points')
    templatePoints = axPhases[0].plot([fullCCFs[ccfName].PlanetPhaseFolded for ccfName in sorted(templateCCFList)], \
                                      [fullCCFs[ccfName].planetRV for ccfName in sorted(templateCCFList)], \
                                      'r*', markersize=16,
                                      label='{} ({} points)'.format('Template', len(templateCCFList)))
    # ax data
    axPhases[1].set_title('data points')
    dataPoints = axPhases[1].plot([fullCCFs[ccfName].PlanetPhaseFolded for ccfName in sorted(dataCCFList)], \
                                  [fullCCFs[ccfName].planetRV for ccfName in sorted(dataCCFList)], \
                                  'r*', markersize=16, label='{} ({} points)'.format('Data', len(dataCCFList)))
    figPhases.tight_layout()

    figPhases.savefig(pdfFile, format='pdf')

    # endregion

    # region [PLOT] template 2D normalised CCFS
    fnPrintLine('PLOT', '2D normalised CCFs')


    maxFlux = np.nanmax([np.nanmax((fullCCFs[ccfName].data / np.nanmedian(fullCCFs[ccfName].data) - 1) * 1000) \
                         for ccfName in fullCCFList])

    minFlux = np.nanmin([np.nanmin((fullCCFs[ccfName].data / np.nanmedian(fullCCFs[ccfName].data) - 1) * 1000) \
                         for ccfName in fullCCFList])

    # plot template nights 2D norm spectra
    fig2DNormCCFsTemplate = fnPlot2DCCFs(fullCCFs , templateCCFList, stellarCCFWidth = starWidth, title= 'normalised CCFs for template construction')
    fig2DNormCCFsTemplate.savefig(pdfFile, format='pdf')

    # endregion

    # region [PLOT] plot data nights 2D norm spectra

    fig2DNormCCFsPlanet = fnPlot2DCCFs(fullCCFs, dataCCFList, stellarCCFWidth = starWidth, title= 'normalised CCFs for planet recovery',\
                                       RVRanges={ccfName:fullCCFs[ccfName].wave[rangeRVsPlanet[ccfName]] for ccfName in dataCCFList})
    fig2DNormCCFsPlanet.savefig(pdfFile, format='pdf')


    # endregion

    # region [PLOT] phase function
    # fnPrintLine('PLOT', 'Phase function')
    # figPhaseFunction, axPhaseFunction = mplt.subplots()
    # phasefunctionFullOrbit = [fnPhaseFunction(planetParams.I, phase + .75, model='Lambert') for phase in
    #                           phasesFullOrbit]
    #
    # axPhaseFunction.plot(phasesFullOrbit, phasefunctionFullOrbit, 'k', alpha=0.5)
    # axPhaseFunction.plot([fullCCFs[ccfName].PlanetPhaseFolded for ccfName in dataSelectedCCFList], \
    #                      [phaseFunction[ccfName] for ccfName in dataSelectedCCFList], 'ro')
    #
    #
    # axPhaseFunction.set_xlim([0, 1])
    # axPhaseFunction.axvline(.75, color='green', label='Opposition')
    # axPhaseFunction.axvline(.25, color='red', label='Transit')
    # axPhaseFunction.axhline(medianPhaseFunction, color='blue',
    #                         label='average phase function:{:.2f}'.format(medianPhaseFunction))
    # axPhaseFunction.set_title('Phase function')
    #
    # figPhaseFunction.savefig(pdfFile, format='pdf')

    # endregion

    # # region [PLOT] Lomb Scargle Periodogram of template ccf
    fnPrintLine('PLOT', 'Lomb Scargle Periodogram of template ccf')
    figPeriodogramNormCCFs, axPeriodogramNormCCFs = mplt.subplots(2, 1)

    axPeriodogramNormCCFs[0].plot(fullCCFs['template'].wave, fullCCFs['template'].data, 'ko', markersize=2)
    axPeriodogramNormCCFs[0].yaxis.set_major_formatter( \
        FuncFormatter(lambda y, pos: '{:.0f}'.format((y - np.nanmedian(fullCCFs['template'].data)) * 1e3)))
    axPeriodogramNormCCFs[0].set_xlim(min(fullCCFs['template'].wave), max(fullCCFs['template'].wave))

    # axPeriodogramNormCCFs[1].plot([1 / freq for freq in periodogramTemplate.freq], periodogramTemplate.power, 'r')
    # axLombscargleNormCCFs[1].text(0.05,0.95, 'max frequency: {:.2f} km/s'.format(periodogramStarTemplateMaxFrequency*fullCCFs[periodogramFileList[0]].CCFStep) ,\
    #                               fontsize = 12, horizontalalignment='left', verticalalignment='top',transform=axLombscargleNormCCFs[1].transAxes)

    axPeriodogramNormCCFs[0].set_title('Template')
    axPeriodogramNormCCFs[1].set_xlabel(r'RV period [km/s]')
    # axPeriodogramNormCCFs[1].set_xscale('log')

    # try:
    #     axPeriodogramNormCCFs[1].set_xlim(float(cfgFile.get('periodograms', 'frequencyTemplateLow')),
    #                                       float(cfgFile.get('periodograms', 'frequencyTemplateHigh')))
    # except:
    #     pass

    figPeriodogramNormCCFs.savefig(pdfFile, format='pdf')

    # endregion

    # region [PLOT] Lomb Scargle Periodogram of recovered planet CCF
    fnPrintLine('PLOT', 'Lomb Scargle Periodogram of recovered planet CCF')
    figPeriodogramPlanetCCF, axPeriodogramPlanetCCF = mplt.subplots(2, 1)

    axPeriodogramPlanetCCF[0].set_title('Planet CCF')
    axPeriodogramPlanetCCF[0].plot(planetWave, planetCCF, 'ko', markersize=1)
    # axPeriodogramPlanetCCF[0].plot(planetWave, fnGauss(planetWave,[ccfAmplitude, ccfMean, ccfFWHM, ccfB]), 'r')
    axPeriodogramPlanetCCF[0].axhline(1 - 1e-5, color='red')
    axPeriodogramPlanetCCF[0].axhline(1 + 1e-5, color='red')
    axPeriodogramPlanetCCF[0].axhline(1 - 1e-4, color='green')
    axPeriodogramPlanetCCF[0].axhline(1 + 1e-4, color='green')
    axPeriodogramPlanetCCF[0].axvline(-3.5, color='magenta')
    axPeriodogramPlanetCCF[0].axvline(3.5, color='magenta')
    axPeriodogramPlanetCCF[0].yaxis.set_major_formatter( \
        FuncFormatter(lambda y, pos: '{:.0f}'.format((y - np.nanmedian(planetCCF)) * 1e6)))
    axPeriodogramPlanetCCF[0].set_ylabel(r'$1-\frac{F_p}{F/*}$[ppm]')
    axPeriodogramPlanetCCF[0].set_xlabel(r'radial velocity [km/s]')

    # axPeriodogramPlanetCCF[1].plot([1 / freq for freq in periodogramPlanetCCF.freq],
    #                                [power for power in periodogramPlanetCCF.power], 'r', markersize=1)
    axPeriodogramPlanetCCF[1].set_xlabel(r'RV period [km/s]')
    # axPeriodogramPlanetCCF[1].set_xscale('log')
    # try:
    #     axPeriodogramPlanetCCF[1].set_xlim(float(cfgFile.get('periodograms', 'frequencyPlanetLow')),
    #                                        float(cfgFile.get('periodograms', 'frequencyPlanetHigh')))
    # except:
    #     pass

    figPeriodogramPlanetCCF.savefig(pdfFile, format='pdf')

    # endregion

    # region [PLOT] RVs and periodogram
    # figPeriodogramRV, axPeriodogramRV = mplt.subplots(2, 1)
    # axPeriodogramRV[0].set_title('Julian date vs Star RV')
    #
    # axPeriodogramRV[0].plot([fullCCFs[ccfName].BJD - 2.45e6 for ccfName in fullCCFList],
    #                         [fullCCFs[ccfName].RV for ccfName in fullCCFList], 'bo')
    # axPeriodogramRV[0].set_xlabel(r'Julian date - 2450 000 [day]')
    # axPeriodogramRV[1].tick_params(axis='x', which='minor', bottom='on')
    # axPeriodogramRV[1].plot([1 / freq for freq in periodogramRV.freq], [power for power in periodogramRV.power], 'r',
    #                         markersize=1)

    # minor_ticks = np.arange(0, max([1 / freq for freq in periodogramRV.freq]), 10)
    # axPeriodogramRV[1].set_xticks(minor_ticks, minor=True)
    # axPeriodogramRV[1].set_xlabel(r'period [day]')

    # axPeriodogramPlanetPhases[1].set_xscale('log')
    #
    # try:
    #     axPeriodogramRV[1].set_xlim(float(cfgFile.get('periodograms', 'frequencyRVLow')),
    #                                 float(cfgFile.get('periodograms', 'frequencyRVHigh')))
    # except:
    #     pass
    #
    # figPeriodogramRV.savefig(pdfFile, format='pdf')

    # endregion

    # region [PLOT] write config parameters
    fnPrintLine('PLOT', 'write config parameters')
    strFormat = '  {:<20s} = {}\n'
    strSettings = ''
    for section_name in cfgFile.sections():
        strSettings += '[{}]\n'.format(section_name)
        for name, value in cfgFile.items(section_name):
            lineLenght = 60
            value = '\n'.join([value[ini: ini + lineLenght] for ini in np.arange(0, len(value), lineLenght)])
            strSettings += strFormat.format(name, value)
        strSettings += '\n'

    figSettings, axSettings = mplt.subplots(figsize=(8.27, 11.69), tight_layout=True)
    axSettings.text(0.05, 0.95, strSettings, fontsize=12, \
                    horizontalalignment='left', \
                    verticalalignment='top', \
                    )
    axSettings.axis('off')
    figSettings.savefig(pdfFile, format='pdf')

    # endregion

    # region [PLOT] recovered parameters
    fnPrintLine('PLOT', 'recovered parameters')
    strFormat = '  {:<20s} = {}\n'
    strPlanetRecoveredParams = '[Gauss Parameters]\n'
    strPlanetRecoveredParams += strFormat.format('Amplitude ($A_{CCF}$)', '{:.0f} ppm'.format(ccfAmplitude * 1e6))
    strPlanetRecoveredParams += strFormat.format('Mean', '{:.1f} km/s'.format(ccfMean))
    strPlanetRecoveredParams += strFormat.format('FWHM', '{:.1f} km/s'.format(ccfFWHM))

    strPlanetRecoveredParams += '\n[Noise]\n'
    strPlanetRecoveredParams += strFormat.format('expected SN50', '{:.1e}'.format(expectedSN))
    strPlanetRecoveredParams += strFormat.format('expected noise', '{:.1e}'.format(1 / expectedSN))
    strPlanetRecoveredParams += strFormat.format('measured noise', '{:.1e}'.format(measuredNoise))

    strPlanetRecoveredParams += '\n[Albedo]\n'
    strPlanetRecoveredParams += strFormat.format(r'Median Phase function $g(\alpha)$',
                                                 '{:.2f} (maximum at opposition; minimum at transit)'.format(
                                                     medianPhaseFunction))
    strPlanetRecoveredParams += strFormat.format(r'Planet radius', '{:.0f} km'.format(planetParams.radiusPlanet))
    strPlanetRecoveredParams += strFormat.format(r'semi-major axis', '{:.0f} km'.format(planetParams.a))
    strPlanetRecoveredParams += strFormat.format(r'Min albedo [from noise]', '{:.2f}'.format(albedoLimitNoise))
    strPlanetRecoveredParams += '\n\n' + r'$A_g = \left(\frac{a}{R_p}\right)^2 \times\, \frac{1}{avg(g(\alpha))} \times\, noise_{measured}$' + '\n\n'

    strPlanetRecoveredParams += strFormat.format(r'median star contrast', '{:.2f}'.format(medianStarCCFContrast))
    strPlanetRecoveredParams += strFormat.format(r'Min albedo [from detected "ccf"]', '{:.2f}'.format(albedoLimitCCF))

    strPlanetRecoveredParams += '\n\n' + r'$A_g = \left(\frac{a}{R_p}\right)^2 \times\, \frac{1}{avg(g(\alpha))} \times\, \frac{amplitude_{planet}}{median(contrast_{star})}$' + '\n\n'

    strPlanetRecoveredParams += strFormat.format(r'$\left(\frac{R_p}{a} \right)^2$',
                                                 '{:.1e}'.format((planetParams.radiusPlanet / planetParams.a) ** 2))

    figParams, axParams = mplt.subplots(figsize=(8.27, 11.69), tight_layout=True)
    axParams.text(0.05, 0.95, strPlanetRecoveredParams, fontsize=12, \
                  horizontalalignment='left', \
                  verticalalignment='top', \
                  )
    axParams.axis('off')
    figParams.savefig(pdfFile, format='pdf')

    # endregion


    pdfFile.close()
    # endregion

    # region --- [TEXT DATA]
    # region [TEXT DATA] initial configuration parameters
    with open(os.path.join(resultsFolder, 'initConfigParams_{}.txt'.format(cfgFile.get('global', 'resultsFolder'))),
              'w') as textSettings:
        textSettings.write(strSettings)

    # endregion

    # region [TEXT DATA] Planet CCF
    #
    strFormat = '{:.10f}\t{:.10f}\t{:.10f}\n'
    strPlanetData = strFormat.replace('10f', '13s').format('# RV', 'planet_flux', 'ccf_Fit')
    for pixWave, pixPlanet, pixFit in zip(planetWave, planetCCF, planetFit):
        strPlanetData += strFormat.format(pixWave, pixPlanet, pixFit)

    with open(os.path.join(resultsFolder, 'planetCCFData_{}.txt'.format(cfgFile.get('global', 'resultsFolder'))),
              'w') as textPlanetData:
        textPlanetData.write(strPlanetData)

    # endregion

    # region [TEXT DATA] List of CCFs for template
    strFormat = '{:.50s}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\n'
    strData = strFormat.replace('10f', '13s').format('# filename', 'BJD', 'Phase', 'star RV', 'planet RV')
    for ccfName in templateCCFList:
        strData += strFormat.format(ccfName.split('/')[-1], fullCCFs[ccfName].BJD, \
                                    fullCCFs[ccfName].PlanetPhaseFolded, fullCCFs[ccfName].RV, \
                                    fullCCFs[ccfName].planetRV)

    with open(os.path.join(resultsFolder, 'templateListData_{}.txt'.format(cfgFile.get('global', 'resultsFolder'))),
              'w') as textTemplateData:
        textTemplateData.write(strData)

    # endregion

    # region [TEXT DATA] List of CCFs available for planet recovery
    strFormat = '{:.50s}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\n'
    strData = strFormat.replace('.10f', '.13s').format('# filename', 'BJD', 'Phase', 'star RV', 'planet RV')
    for ccfName in dataCCFList:
        strData += strFormat.format(ccfName.split('/')[-1], fullCCFs[ccfName].BJD, \
                                    fullCCFs[ccfName].PlanetPhaseFolded, fullCCFs[ccfName].RV, \
                                    fullCCFs[ccfName].planetRV)

    with open(os.path.join(resultsFolder, 'planetRecoveryData_{}.txt'.format(cfgFile.get('global', 'resultsFolder'))),
              'w') as textTemplateData:
        textTemplateData.write(strData)

    # endregion

    # region [TEXT DATA] List of CCFs selected for planet recovery
    strFormat = '{:.50s}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\n'
    strData = strFormat.replace('10f', '13s').format('# filename', 'BJD', 'Phase', 'star RV', 'planet RV')
    for ccfName in dataSelectedCCFList:
        strData += strFormat.format(ccfName.split('/')[-1], fullCCFs[ccfName].BJD, \
                                    fullCCFs[ccfName].PlanetPhaseFolded, fullCCFs[ccfName].RV, \
                                    fullCCFs[ccfName].planetRV)

    with open(os.path.join(resultsFolder, \
                           'selectedPlanetRecoveryData_{}.txt'.format(cfgFile.get('global', 'resultsFolder'))),
              'w') as textTemplateData:
        textTemplateData.write(strData)

    # endregion

    # region [TEXT DATA] recovered parameters

    with open(
            os.path.join(resultsFolder, 'planetRecoveredParams_{}.txt'.format(cfgFile.get('global', 'resultsFolder'))),
            'w') as textSettings:
        textSettings.write(strPlanetRecoveredParams)

        # endregion





        # endregion


# ======================================================================================================================
#   RunTime
# ======================================================================================================================
if __name__ == '__main__':
    __main__()
