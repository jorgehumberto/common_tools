#!/usr/bin/env python
'''

create settings and scripts to run esorex:

esorex --config=esorex.cfg --recipe-config=51Peg_UVESredchain.cfg recipe sofSettings.sof

'''

__author__ = 'jmartins'

# region --- Python Modules
import os
import sys

from astropy.io import fits
# endregion


# region --- User Modules

from bin.manipFiles import fnCreateFolder
from modules.InOut import fnPrintLine

# endregion

# region --- Paths, Folders and filenames

# Define Work Paths
WorkPath = os.path.abspath(os.getenv('WorkPath'))
ReductPath = os.path.abspath(os.getenv('ReductPath'))

sys.path.append(WorkPath)
settings = __import__(sys.argv[1].split('.')[0])

# Defining folders
baseRawDataFolder = os.path.abspath('{}/RawData/{}'.format(WorkPath, settings.rawDataFolder))
baseReducedDataFolder = os.path.abspath('{}/reduced/{}/'.format(WorkPath, settings.rawDataFolder))
fnCreateFolder(baseReducedDataFolder)
scriptsFolder = os.path.abspath('{}/scripts/{}/'.format(WorkPath, settings.rawDataFolder))
fnCreateFolder(scriptsFolder)

rawFileList = ['{}/{}'.format(baseRawDataFolder, filename) for filename in os.listdir(baseRawDataFolder) if
               filename.endswith('.fits')]
stdCalibrationFilelist = [filename for filename in os.listdir('{}/RawData/UVESCalibrationFiles/'.format(WorkPath)) if
                          filename.endswith('.fits')]

# defining filenames
# recipeSettings = '{}/scripts/customRecipe_{}.cfg'.format(ReductPath,settings.rawDataFolder)
esorexCommandfilename = '{}/{}.sh'.format(WorkPath, settings.rawDataFolder)

arm = settings.arm

esorexSettings = '{}/scripts/defaults/{}'.format(ReductPath, 'defaults_esorex.cfg')


# endregion

# region --- Main

def __main__():
    '''
        Routine to generate esorex commands to reduce UVES data via ESO's reduction recipes:
        1) create master bias - uves_cal_bias
        2) create guess order and line tables - uves_cal_predict
        3) create order table from order guess table  uves_cal_orderpos
        4) create master flat - uves_cal_mflat
        5) create line calibration table from guess table - uves_cal_wavecal
        6) reduce science frame - uves_obs_scired

    :return:
    '''

    # region - Define File Lists ---------------------------------------------------------------------------------------
    fnPrintLine('FILES', 'Defining file lists')
    # science/data files
    scienceFileList = [filename for filename in rawFileList \
                       if 'OBJECT,POINT'.lower() in fits.getheader(filename)['HIERARCH ESO DPR TYPE'].lower()]

    # bias files
    biasFileList = [filename for filename in rawFileList \
                    if 'BIAS'.lower() in fits.getheader(filename)['HIERARCH ESO DPR TYPE'].lower()]

    # flat files
    flatFileList = [filename for filename in rawFileList \
                    if 'FLAT'.lower in fits.getheader(filename)['HIERARCH ESO DPR TYPE'].lower()]

    esorexBaseString = 'esorex --config={esorexCfg} --recipe-config={recipeSettings} {recipeName} {fileList} \n'

    # endregion --------------------------------------------------------------------------------------------------------

    # Generate reduction scripts


    # region --- 1) create master bias -     uves_cal_mbiasls
    python

    fnCreateFolder('{}/masterBias/'.format(baseReducedDataFolder))
    # writing file list
    masterBiasFileList = '{}/fileList_uves_cal_mbias.sof'.format(scriptsFolder)

    with open(masterBiasFileList, 'w') as writeFile:
        writeFile.write('\n'.join(biasFileList))

    print '\n'.join(biasFileList)

    # writing esorex settings file
    with open(esorexSettings, 'r') as esorexDefaultsFile:
        esorexDefaults = str(esorexDefaultsFile.read())

    esorexSettings_uves_cal_mbias = esorexDefaults.replace('DEFAULT_OUTPUT_DIR',
                                                           '{}/masterBias/'.format(baseReducedDataFolder))
    esorexSettings_uves_cal_mbias = esorexSettings_uves_cal_mbias.replace('DEFAULT_OUTPUT_PREFIX', 'masterBias')
    esorexSettingsFile_uves_cal_mbias = '{}/esorex_uves_cal_mbias.cfg'.format(scriptsFolder)

    with open(esorexSettingsFile_uves_cal_mbias, 'w') as writeFile:
        writeFile.write(esorexSettings_uves_cal_mbias)

    esorexCommands_uves_cal_mbias = esorexBaseString.format( \
        esorexCfg=esorexSettingsFile_uves_cal_mbias, \
        recipeSettings=os.path.abspath(settings.uves_cal_mbias), \
        recipeName='uves_cal_mbias', \
        fileList=masterBiasFileList)

    print esorexCommands_uves_cal_mbias

    # endregion --------------------------------------------------------------------------------------------------------

    # region --- 2) create guess order and line tables - uves_cal_predict




    # esorex_uves_cal_predict = esorexBaseString.format(esorexSettings, recipeSettingsFile, 'uves_cal_predict',
    #                                                  fileListUvesObsScired))




    # endregion --------------------------------------------------------------------------------------------------------

    # region  --- 3) create order table from order guess table      uves_cal_orderpos



    # esorex_uves_cal_orderpos = esorexBaseString.format(esorexSettings, recipeSettingsFile, 'uves_cal_orderpos',
    #                                                  fileListUvesObsScired))




    # endregion --------------------------------------------------------------------------------------------------------

    # region --- 4) create master flat - uves_cal_mflat



    # esorex_uves_cal_mflat = esorexBaseString.format(esorexSettings, recipeSettingsFile, 'uves_cal_mflat',
    #                                                  fileListUvesObsScired))





    # endregion --------------------------------------------------------------------------------------------------------

    # region --- 5) create line calibration table from guess table - uves_cal_wavecal


    # esorex_uves_cal_wavecal = esorexBaseString.format(esorexSettings, recipeSettingsFile, 'uves_cal_wavecal',
    #                                                  fileListUvesObsScired))





    # endregion --------------------------------------------------------------------------------------------------------

    # region --- 6) reduce science frame - uves_obs_scired




    # esorex_uves_obs_scired = esorexBaseString.format(esorexSettings, recipeSettingsFile, 'uves_obs_scired', fileListUvesObsScired))



    # endregion




    # print 'Create ESOREX scripts'
    #
    # # Create settings for esorex
    # with open('{}/scripts/defaults/esorex.cfg'.format(ReductPath)) as defaultEsorexFile:
    #     defaultEsorex = defaultEsorexFile.read()
    #
    # # create settings for recipe
    # with open('{}/scripts/defaults/{}.cfg'.format(ReductPath, settings.recipeName)) as defaultRecipeFile:
    #     customRecipe = defaultRecipeFile.read()
    #
    # with open(recipeSettings,'w') as customRecipeFile:
    #     customRecipeFile.write(customRecipe)
    #
    #
    #
    # # Create SOF files for esorex - standard calibration files
    # stdCalibsString = ''
    # for filename in stdCalibrationFilelist:
    #     headerSTDCalib  =   fits.getheader('{}/RawData/UVESCalibrationFiles/{}'.format(WorkPath,filename))
    #     stdCalibsString   +=  '{}/RawData/UVESCalibrationFiles/{}\t{}\n'.format(WorkPath, filename, headerSTDCalib['HIERARCH ESO PRO CATG'])
    #
    #
    # # Create SOF files for esorex - raw night calibration files
    # rawCalibrationFileList = [filename for filename in rawFileList if 'CALIB' in fits.getheader('{}/{}'.format(baseRawDataFolder, filename))['HIERARCH ESO DPR CATG']]
    # rawSOFCalibrationString = stdCalibsString
    #
    # for filename in rawCalibrationFileList:
    #     headerFile  =   fits.getheader('{}/{}'.format(baseRawDataFolder, filename))
    #     header = fnClassificationRules_UVES(headerFile)
    #     rawSOFCalibrationString   +=  '{}/{}\t{}\n'.format(baseRawDataFolder,filename, header['D0.CATG'])
    #
    #
    # # Create SOF files for esorex -  raw science files
    # rawScienceFileList = [filename for filename in rawFileList if 'SCIENCE' in fits.getheader('{}/{}'.format(baseRawDataFolder, filename))['HIERARCH ESO DPR CATG']]
    #
    # esorexCommand = ''
    # for filename in rawScienceFileList:
    #     print filename
    #
    #     rawScienceFileKey = 'UVES.{}.{}'.format(header['DATE-OBS'], header['HIERARCH ESO INS PATH '])
    #     rawScienceString = rawSOFCalibrationString
    #     headerFile  =   fits.getheader('{}/{}'.format(baseRawDataFolder, filename))
    #     header = fnClassificationRules_UVES(headerFile)
    #     rawScienceString   +=  '{}/{}\t{}\n'.format(baseRawDataFolder,filename, header['D0.CATG'])
    #
    #     esorexSettings = '{}/scripts/customESOREX_{}_{}.cfg'.format( ReductPath,settings.rawDataFolder, rawScienceFileKey)
    #     customESOREX = defaultEsorex.replace('DEFAULT_OUTPUT_DIR', '{}/reduced/{}/'.format(WorkPath, settings.rawDataFolder))
    #     customESOREX = customESOREX.replace('DEFAULT_OUTPUT_PREFIX', rawScienceFileKey)
    #
    #     with open(esorexSettings,'w') as customEsorexFile:
    #         customEsorexFile.write(customESOREX)
    #
    #     print '\n\ncustom esorex settings file has been created: {}\n\n'.format(esorexSettings)
    #
    #
    #     sofSettings = '{}/scripts/rawFiles2Reduce_{}_{}.sof'.format(ReductPath,settings.rawDataFolder,rawScienceFileKey)
    #
    #     with open(sofSettings,'w') as sofFile:
    #         sofFile.write(rawScienceString)
    #
    #     print '\n\nlist of raw data files to be reduced has been created: {}\n\n'.format(rawScienceFileKey)
    #
    #
    #
    #     # full esorex command
    #     esorexCommand += 'esorex --config={0} --recipe-config={1} {2} {3} \n'.format(esorexSettings, recipeSettings,settings.recipeName, sofSettings )
    #
    #
    # with open(esorexCommandfilename,'a') as esorexCommandFile:
    #     esorexCommandFile.write(esorexCommand)
    #
    #
    # print '\n\nto reduce your data please run: {}\n\n{}\n'.format(esorexCommandfilename , esorexCommand)


    return


# endregion

# region --- RunTime

if __name__ == '__main__':
    __main__()


# endregion
