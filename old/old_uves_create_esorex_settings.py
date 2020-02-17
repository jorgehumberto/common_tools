#!/usr/bin/env python
'''

create settings and scripts to run esorex:

esorex --config=esorex.cfg --recipe-config=51Peg_UVESredchain.cfg recipe sofSettings.sof

'''
__author__ = 'jmartins'

# ======================================================================================================================
#   Python Modules
# ======================================================================================================================
import os
import sys

from astropy.io import fits
# ======================================================================================================================
#   User Modules
# ======================================================================================================================

from bin.classificationRules import *
from bin.manipFiles import fnCreateFolder

# ======================================================================================================================
#   Settings
# ======================================================================================================================
# Define Work Paths
WorkPath = os.path.abspath(os.getenv('WorkPath'))
ReductPath = os.path.abspath(os.getenv('ReductPath'))

sys.path.append(WorkPath)
settings = __import__(sys.argv[1].split('.')[0])

# Defining folders
baseRawDataFolder = os.path.abspath('{}/RawData/{}'.format(WorkPath, settings.rawDataFolder))
fnCreateFolder('{}/reduced/{}/'.format(WorkPath, settings.rawDataFolder))

rawFileList = [fileName for fileName in os.listdir(baseRawDataFolder) if fileName.endswith('.fits')]
stdCalibrationFilelist = [fileName for fileName in os.listdir('{}/RawData/UVESCalibrationFiles/'.format(WorkPath)) if
                          fileName.endswith('.fits')]

# defining filenames
recipeSettings = '{}/scripts/customRecipe_{}.cfg'.format(ReductPath, settings.rawDataFolder)
esorexCommandFilename = '{}/{}.sh'.format(WorkPath, settings.rawDataFolder)


# ======================================================================================================================
#   Main
# ======================================================================================================================

def __main__():
    print 'Create ESOREX scripts'

    # Create settings for esorex
    with open('{}/scripts/defaults/esorex.cfg'.format(ReductPath)) as defaultEsorexFile:
        defaultEsorex = defaultEsorexFile.read()

    # create settings for recipe
    with open('{}/scripts/defaults/{}.cfg'.format(ReductPath, settings.recipeName)) as defaultRecipeFile:
        customRecipe = defaultRecipeFile.read()

    with open(recipeSettings, 'w') as customRecipeFile:
        customRecipeFile.write(customRecipe)

    # Create SOF files for esorex - standard calibration files
    stdCalibsString = ''
    for fileName in stdCalibrationFilelist:
        headerSTDCalib = fits.getheader('{}/RawData/UVESCalibrationFiles/{}'.format(WorkPath, fileName))
        stdCalibsString += '{}/RawData/UVESCalibrationFiles/{}\t{}\n'.format(WorkPath, fileName,
                                                                             headerSTDCalib['HIERARCH ESO PRO CATG'])

    # Create SOF files for esorex - raw night calibration files
    rawCalibrationFileList = [fileName for fileName in rawFileList if
                              'CALIB' in fits.getheader('{}/{}'.format(baseRawDataFolder, fileName))[
                                  'HIERARCH ESO DPR CATG']]
    rawSOFCalibrationString = stdCalibsString

    for fileName in rawCalibrationFileList:
        headerFile = fits.getheader('{}/{}'.format(baseRawDataFolder, fileName))
        header = fnClassificationRules_UVES(headerFile)
        rawSOFCalibrationString += '{}/{}\t{}\n'.format(baseRawDataFolder, fileName, header['D0.CATG'])

    # Create SOF files for esorex -  raw science files
    rawScienceFileList = [fileName for fileName in rawFileList if
                          'SCIENCE' in fits.getheader('{}/{}'.format(baseRawDataFolder, fileName))[
                              'HIERARCH ESO DPR CATG']]

    esorexCommand = ''
    for fileName in rawScienceFileList:
        print fileName

        rawScienceFileKey = 'UVES.{}.{}'.format(header['DATE-OBS'], header['HIERARCH ESO INS PATH '])
        rawScienceString = rawSOFCalibrationString
        headerFile = fits.getheader('{}/{}'.format(baseRawDataFolder, fileName))
        header = fnClassificationRules_UVES(headerFile)
        rawScienceString += '{}/{}\t{}\n'.format(baseRawDataFolder, fileName, header['D0.CATG'])

        esorexSettings = '{}/scripts/customESOREX_{}_{}.cfg'.format(ReductPath, settings.rawDataFolder,
                                                                    rawScienceFileKey)
        customESOREX = defaultEsorex.replace('DEFAULT_OUTPUT_DIR',
                                             '{}/reduced/{}/'.format(WorkPath, settings.rawDataFolder))
        customESOREX = customESOREX.replace('DEFAULT_OUTPUT_PREFIX', rawScienceFileKey)

        with open(esorexSettings, 'w') as customEsorexFile:
            customEsorexFile.write(customESOREX)

        print '\n\ncustom esorex settings file has been created: {}\n\n'.format(esorexSettings)

        sofSettings = '{}/scripts/rawFiles2Reduce_{}_{}.sof'.format(ReductPath, settings.rawDataFolder,
                                                                    rawScienceFileKey)

        with open(sofSettings, 'w') as sofFile:
            sofFile.write(rawScienceString)

        print '\n\nlist of raw data files to be reduced has been created: {}\n\n'.format(rawScienceFileKey)

        # full esorex command
        esorexCommand += 'esorex --config={0} --recipe-config={1} {2} {3} \n'.format(esorexSettings, recipeSettings,
                                                                                     settings.recipeName, sofSettings)

    with open(esorexCommandFilename, 'a') as esorexCommandFile:
        esorexCommandFile.write(esorexCommand)

    print '\n\nto reduce your data please run: {}\n\n{}\n'.format(esorexCommandFilename, esorexCommand)


# ======================================================================================================================
#   RunTime
# ======================================================================================================================
if __name__ == '__main__':
    __main__()
