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

import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as mplt
import numpy as np

os.system('cls' if os.name == 'nt' else 'clear')

# endregion


# region --- User Modules
from modules.InOut import fnPrintLine
from classes.classesPlanet import clsPlanetParametersOLD
from recipes.manipCCF import fnOpenHARPSFits
from models.modelsOrbit import fnRVStarOrbitElipse

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


# endregion


# endregion


# ======================================================================================================================
#   Main
# ======================================================================================================================

def __main__():
    # region --- Define planet parameters




    planetParams = clsPlanetParametersOLD(params={param: value for param, value in cfgFile.items('orbital_params')})
    planetParams.sysRV = 0.0
    # endregion




    # region --- Get files list
    # data files
    fnPrintLine('CCF', 'get ccf list')
    fullCCFList = sorted(
        [os.path.join(scienceInputFolder, fileName) for fileName in os.listdir(scienceInputFolder) \
         if fileName.lower().endswith('fits')])

    # region --- load fits files
    fnPrintLine('CCF', 'Extracting CCFs, please wait')
    fullCCFs = {ccfName: fnOpenHARPSFits(ccfName, headerOnly=True) for ccfName in sorted(fullCCFList)}

    # region --- Finding star CCF parameters
    fnPrintLine('CCF', 'Computing star and planet RVs')
    # define xxx axis

    mplt.ion()
    for ccfName in sorted(fullCCFList):
        fullCCFs[ccfName].PlanetPhaseFolded = ((fullCCFs[
                                                    ccfName].BJD - planetParams.t0) / planetParams.period + np.radians(
            planetParams.w) / (2 * np.pi)) % 1
        fullCCFs[ccfName].RV = fullCCFs[ccfName].RVC = fnRVStarOrbitElipse(planetParams, fullCCFs[
            ccfName].PlanetPhaseFolded - np.radians(planetParams.w) / (2 * np.pi))
        # fullCCFs[ccfName].RV = fullCCFs[ccfName].RVC = fnRVStarOrbitElipse(planetParams, fullCCFs[ccfName].PlanetPhaseFolded - np.radians(planetParams.w)/(2*np.pi))

        fullCCFs[ccfName].planetRV = -fullCCFs[ccfName].RV / planetParams.massRatio

    figTemplate, axTemplate = mplt.subplots(figsize=(12, 8))
    maxFWHM = np.nanmax([fullCCFs[ccfName].FWHM for ccfName in fullCCFList])

    allPoints = axTemplate.plot([fullCCFs[ccfName].PlanetPhaseFolded for ccfName in sorted(fullCCFList)], \
                                [fullCCFs[ccfName].planetRV for ccfName in sorted(fullCCFList)], 'bd', markersize=8)

    axTemplate.grid(b=True, which='major', color='k', linestyle='--', alpha=.5)
    axTemplate.set_ylim(-1.1 * planetParams.k2, 1.1 * planetParams.k2)
    axTemplate.set_yticks(
        np.arange(-(np.floor((1.1 * planetParams.k2) / 10) * 10), ((np.ceil(1.1 * planetParams.k2) / 10) * 10), 10))
    axTemplate.axhline(0.0, label=r'$RV_{CM}$')
    axTemplate.axvline(0.25, label='transit', color='red')
    axTemplate.axvline(0.75, label='opposition', color='green')

    axTemplate.axhline(-2 * maxFWHM, color='blue', alpha=.5, label=r'$2 \times FWHM$')
    axTemplate.axhline(2 * maxFWHM, color='blue', alpha=.5, )

    axTemplate.legend(loc='best', title=r'$\phi = 0$ - Ascending node')
    axTemplate.set_xlim(0, 1)
    axTemplate.set_xticks(np.arange(0, 1, .05))
    axTemplate.set_title('template selection - {}'.format(str(cfgFile.get('global', 'dataFolder'))))
    figTemplate.show()
    # endregion

    # region --- get template file list
    # fnPrintLine('template', 'Use all spectra for template?')
    # customTemplate = raw_input('\t\tyes/no: ').lower()
    selectPoints = 'no'
    while selectPoints == 'no':
        templatePhaseRanges = []
        phaseRange = [0, 1]
        print
        fnPrintLine('template', 'Input planet phase ranges for template')
        fnPrintLine(None, '(leave empty when complete)')
        while phaseRange:
            phaseRange = raw_input('\t\tmin_Phase, max_Phase: ').strip('"')
            if phaseRange:
                templatePhaseRanges.append(
                    [float(phaseRange.strip(' ').split(',')[0]), float(phaseRange.strip(' ').split(',')[-1])])
                templateCCFList = [ccfName for ccfName in sorted(fullCCFList) \
                                   for minPhase, maxPhase in templatePhaseRanges \
                                   if minPhase <= float(fullCCFs[ccfName].PlanetPhaseFolded) <= maxPhase \
                                   ]
                selectedPoints, = axTemplate.plot(
                    [fullCCFs[ccfName].PlanetPhaseFolded for ccfName in sorted(templateCCFList)], \
                    [fullCCFs[ccfName].planetRV for ccfName in sorted(templateCCFList)], 'r*', markersize=16)
                figTemplate.show()

                # mplt.plot([fullCCFs[ccfName].PlanetPhaseFolded for ccfName in sorted(fullCCFList)],\
                #   [fullCCFs[ccfName].planetRV for ccfName in sorted(fullCCFList)], 'bd', markersize = 8)

        fnPrintLine(None, 'Is the selection correct?')
        selectPoints = raw_input('\t\tyes/no: ').lower()

        axTemplate.lines.remove(selectedPoints)

    print

    # endregion

    # endregion

    figData, axData = figTemplate, axTemplate
    axData.set_title('data selection - {}'.format(str(cfgFile.get('global', 'dataFolder'))))
    figData.show()

    selectPoints = 'no'
    while selectPoints == 'no':
        dataRVRanges = []
        phaseRange = [0, 1]
        print
        fnPrintLine('data', 'Input planet phase waveRange for planet recovery')
        fnPrintLine(None, '(leave empty when complete)')
        while phaseRange:
            phaseRange = raw_input('\t\tmin_Phase, max_Phase: ').strip('"')
            if phaseRange:
                dataRVRanges.append(
                    [float(phaseRange.strip(' ').split(',')[0]), float(phaseRange.strip(' ').split(',')[-1])])

            dataCCFList = [ccfName for ccfName in sorted(fullCCFList) \
                           for minPhase, maxPhase in dataRVRanges \
                           if minPhase <= float(fullCCFs[ccfName].PlanetPhaseFolded) <= maxPhase \
                           ]
            selectedPoints, = mplt.plot([fullCCFs[ccfName].PlanetPhaseFolded for ccfName in sorted(dataCCFList)], \
                                        [fullCCFs[ccfName].planetRV for ccfName in sorted(dataCCFList)], 'r*',
                                        markersize=16)
            figTemplate.show()

        fnPrintLine(None, 'Is the selection correct?')
        selectPoints = raw_input('\t\tyes/no: ').lower()
        axData.lines.remove(selectedPoints)
        mplt.show()

    # region --- updating file lists and  config
    fnPrintLine('config', 'updating file lists and  config file')
    with open(
            os.path.join(scienceInputFolder, '{}_templateFileList.txt'.format(cfgFile.get('global', 'resultsFolder'))),
            'w') as templateFileList:
        templateFileList.write('\n'.join([fileName.split('/')[-1] for fileName in templateCCFList]))
        cfgFile.set('detection_settings', 'templateFilelist',
                    '{}_templateFileList.txt'.format(cfgFile.get('global', 'resultsFolder')))

    with open(os.path.join(scienceInputFolder, '{}_dataFilelist.txt'.format(cfgFile.get('global', 'resultsFolder'))),
              'w') as dataFileList:
        dataFileList.write('\n'.join([fileName.split('/')[-1] for fileName in dataCCFList]))
        cfgFile.set('detection_settings', 'dataFilelist',
                    '{}_dataFilelist.txt'.format(cfgFile.get('global', 'resultsFolder')))

    with open(os.path.abspath(sys.argv[1]), 'wb') as cfgFileName:
        cfgFile.write(cfgFileName)

    fnPrintLine('config', 'config file {} has been updated with file lists'.format(os.path.abspath(sys.argv[1])))
    # endregion


# ======================================================================================================================
#   RunTime
# ======================================================================================================================
if __name__ == '__main__':
    __main__()
