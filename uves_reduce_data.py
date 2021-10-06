#!/usr/bin/env python
'''



'''

__author__ = 'Jorge Martins'
__email__ = 'jorge.martins@iastro.pt'
__version__ = 'alpha 1'
__date__ = '2015/12/02'

# region --- Python Modules
import ConfigParser
import os
import sys

import cpl
from astropy.io import fits
import matplotlib.pyplot as mplt

# endregion


# region --- User Modules

from modules.InOut import fnPrintLine
from modules.recipeTools import fnSetRecipe
from modules.fileManip import fngetRecipeFileList, fnUpdateProductList
from modules.spectraTools import fnInferInstrumentResponse

# endregion

# region --- Init Terminal
os.system("cls")
os.system("clear")
fnPrintLine(None, '')
fnPrintLine(None, '', align='center', flush='=')
fnPrintLine(None, 'reduce UVES data CPL wrapper \t(version: {:<})'.format(__version__), align='center')
# fnPrintLine(None, '', align = 'center', flush = '=')
fnPrintLine(None, '')
fnPrintLine(None, 'Author:{:<} \temail: {:<}'.format(__author__, __email__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, 'Last Update:{:<}'.format(__date__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, '', align='center', flush='=')
fnPrintLine(None, '')

# endregion

# region --- Parse Settings
cfgFileName = sys.argv[1]
cfgFile = ConfigParser.RawConfigParser(allow_no_value=True)
cfgFile.read(cfgFileName)

# endregion


# defining recipes
recipeList = ['uves_cal_mbias', 'uves_cal_predict', 'uves_cal_mflat', 'uves_cal_orderpos', 'uves_cal_wavecal',
              'uves_obs_scired', 'uves_cal_blaze']
calibList = ['uves_cal_mbias', 'uves_cal_predict', 'uves_cal_mflat', 'uves_cal_orderpos', 'uves_cal_wavecal',
             'uves_response']

try:
    recipes2run = cfgFile.get('global', 'recipes2run')
except ConfigParser.NoOptionError:
    recipes2run = recipeList


clpRecipes = {}
clpResults = {}

# region --- Paths, Folders and filenames

# Define Work Paths
WorkPath = os.path.abspath(os.getenv('WorkPath'))
ReductPath = os.path.abspath(os.getenv('ReductPath'))

sys.path.append(WorkPath)
# settings = __import__(sys.argv[1].split('.')[0])

# Defining folders

# data
baseDataFolder_raw = os.path.join(WorkPath, 'RawData', cfgFile.get('global', 'workDir').strip("'"))
baseDataFolder_reduced = os.path.join(WorkPath, 'reduced', cfgFile.get('global', 'workDir').strip("'"))

# standart calibrations
baseCalibrations_raw = os.path.join(WorkPath, 'RawData/UVESCalibrationFiles/')

# recipe folders
outputFolders = {recipe: os.path.join(baseDataFolder_reduced, recipe) for recipe in recipeList}
tmpFolder = os.path.join(baseDataFolder_reduced, 'tmp')

# get filelists
# data
rawData_fileList = [os.path.join(baseDataFolder_raw, filename) for filename in os.listdir(baseDataFolder_raw) if
                    filename.endswith('.fits')]

# Recipe file lists:
recipeFileLists = {}

# endregion

# region --- Recipe Settings

# defining keywords per recipe
# http://www.eso.org/observing/dfo/quality/UVES/pipeline/recipe_calib.html

recipeKeywords = {}
recipeKeywords['uves_cal_mbias'] = {'DPR_CATG': 'CALIB', 'DPR_TYPE': 'BIAS', \
                                    'TAG': 'BIAS_{}'.format(cfgFile.get('global', 'arm')), 'PRO_CATG': ['MASTER_BIAS']}

recipeKeywords['uves_cal_predict'] = {'DPR_CATG': 'CALIB', 'DPR_TYPE': 'FMTCHK', \
                                      'TAG': 'ARC_LAMP_FORM_{}'.format(cfgFile.get('global', 'arm')),
                                      'PRO_CATG': ['LINE_TABLE', 'ORDER_TABLE', 'BACKGROUND_TABLE']}

recipeKeywords['uves_cal_orderpos'] = {'DPR_CATG': 'CALIB', 'DPR_TYPE': 'ORDERDEF', \
                                       'TAG': 'ORDER_FLAT_{}'.format(cfgFile.get('global', 'arm')), 'PRO_CATG': []}

recipeKeywords['uves_cal_mflat'] = {'DPR_CATG': 'CALIB', 'DPR_TYPE': 'FLAT', \
                                    'TAG': 'FLAT_{}'.format(cfgFile.get('global', 'arm')), 'PRO_CATG': []}

recipeKeywords['uves_cal_wavecal'] = {'DPR_CATG': 'CALIB', 'DPR_TYPE': 'WAVE', \
                                      'TAG': 'ARC_LAMP_{}'.format(cfgFile.get('global', 'arm')), 'PRO_CATG': []}

recipeKeywords['uves_obs_scired'] = {'DPR_CATG': 'SCIENCE', 'DPR_TYPE': 'OBJECT', \
                                     'TAG': 'SCIENCE_{}'.format(cfgFile.get('global', 'arm')), 'PRO_CATG': []}

recipeKeywords['uves_cal_blaze'] = {'DPR_CATG': 'CALIB', 'DPR_TYPE': 'FLAT', \
                                    'TAG': 'SCIENCE_{}'.format(cfgFile.get('global', 'arm')), 'PRO_CATG': []}
# standard calibrations
calibrations_fileList = {}
calibrations_fileList['LINE_REFER_TABLE'] = [os.path.join(baseCalibrations_raw, 'thargood_3.fits')]
calibrations_fileList['LINE_INTMON_TABLE'] = [os.path.join(baseCalibrations_raw, 'thar_bright.fits')]


# endregion




# region --- Main



def main_function():
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

    # initialise esorex

    try:
        cpl.esorex.init(source=cfgFile.get('global', 'esorexrc'))
    except:
        cpl.esorex.init()

    # region - Define File Lists ---------------------------------------------------------------------------------------

    # endregion --------------------------------------------------------------------------------------------------------


    # region --- 1) create master bias -     uves_cal_mbias
    fnPrintLine('INIT', 'Recipes to run'.format(','.join([recipe for recipe in recipes2run])))
    recipe = 'uves_cal_mbias'
    if recipe in recipes2run:
        fnPrintLine('CALIB', 'starting {}, please wait...'.format(recipe))
        # get input file list
        recipeFileLists[recipe] = fngetRecipeFileList(recipeName=recipe, \
                                                      dataList=rawData_fileList, \
                                                      keyword=recipeKeywords[recipe]['DPR_TYPE'])

        # setup recipe "uves_cal_mbias"
        clpRecipes[recipe] = fnSetRecipe(recipe, customConfig=cfgFile)

        # run recipe "uves_cal_predict"
        clpResults[recipe] = clpRecipes[recipe](raw={recipeKeywords[recipe]['TAG']: recipeFileLists[recipe]}, \
                                                output_dir=outputFolders[recipe], tmp_dir=tmpFolder)

        clpResults[recipe] = None

        fnPrintLine('CALIB', '{} complete'.format(recipe))

    # endregion --------------------------------------------------------------------------------------------------------

    # region --- 2) create guess order and line tables - uves_cal_predict

    recipe = 'uves_cal_predict'
    if recipe in recipes2run:
        fnPrintLine('CALIB', 'starting {}, please wait...'.format(recipe))
        # get input file list
        recipeFileLists[recipe] = fngetRecipeFileList(recipeName=recipe, \
                                                      dataList=rawData_fileList, \
                                                      keyword=recipeKeywords[recipe]['DPR_TYPE'])

        calibrations_fileList.update(fnUpdateProductList(baseDataFolder_reduced, productList=calibrations_fileList))

        # setup recipe "uves_cal_predict"
        clpRecipes[recipe] = fnSetRecipe(recipe, cfgFile)

        # run recipe "uves_cal_predict"

        clpResults[recipe] = clpRecipes[recipe](raw={recipeKeywords[recipe]['TAG']: recipeFileLists[recipe]}, \
                                                calib=calibrations_fileList, output_dir=outputFolders[recipe],
                                                tmp_dir=tmpFolder)
        clpResults[recipe] = None

        fnPrintLine('CALIB', '{} complete'.format(recipe))

    # endregion --------------------------------------------------------------------------------------------------------

    # region  --- 3) create order table from order guess table      uves_cal_orderpos

    recipe = 'uves_cal_orderpos'
    if recipe in recipes2run:
        fnPrintLine('CALIB', 'starting {}, please wait...'.format(recipe))
        # get input file list
        recipeFileLists[recipe] = fngetRecipeFileList(recipeName=recipe, \
                                                      dataList=rawData_fileList, \
                                                      keyword=recipeKeywords[recipe]['DPR_TYPE'])

        # update calibration lists
        calibrations_fileList.update(fnUpdateProductList(baseDataFolder_reduced, productList=calibrations_fileList))

        # setup recipe "uves_cal_orderpos"
        clpRecipes[recipe] = fnSetRecipe(recipe, cfgFile)

        # run recipe "uves_cal_predict"
        clpResults[recipe] = clpRecipes[recipe](raw={recipeKeywords[recipe]['TAG']: recipeFileLists[recipe]},
                                                calib=calibrations_fileList, \
                                                output_dir=outputFolders[recipe], tmp_dir=tmpFolder)

        clpResults[recipe] = None

        fnPrintLine('CALIB', '{} complete'.format(recipe))

    # endregion --------------------------------------------------------------------------------------------------------

    # region --- 4) create master flat - uves_cal_mflat
    recipe = 'uves_cal_mflat'
    if recipe in recipes2run:
        fnPrintLine('CALIB', 'starting {}, please wait...'.format(recipe))
        # get input file list
        recipeFileLists[recipe] = fngetRecipeFileList(recipeName=recipe, \
                                                      dataList=rawData_fileList, \
                                                      keyword=recipeKeywords[recipe]['DPR_TYPE'])

        # update calibration lists
        calibrations_fileList.update(fnUpdateProductList(baseDataFolder_reduced, productList=calibrations_fileList))

        # setup recipe
        clpRecipes[recipe] = fnSetRecipe(recipe, cfgFile)

        # run recipe "uves_cal_predict"
        clpResults[recipe] = clpRecipes[recipe](raw={recipeKeywords[recipe]['TAG']: recipeFileLists[recipe]},
                                                calib=calibrations_fileList, \
                                                output_dir=outputFolders[recipe], tmp_dir=tmpFolder)

        clpResults[recipe] = None

        fnPrintLine('CALIB', '{} complete'.format(recipe))

    # endregion --------------------------------------------------------------------------------------------------------

    # region --- 5) create line calibration table from guess table - uves_cal_wavecal
    recipe = 'uves_cal_wavecal'
    if recipe in recipes2run:
        fnPrintLine('CALIB', 'startings {}, please wait...'.format(recipe))
        # get input file list
        recipeFileLists[recipe] = fngetRecipeFileList(recipeName=recipe, \
                                                      dataList=rawData_fileList, \
                                                      keyword=recipeKeywords[recipe]['DPR_TYPE'])

        # update calibration lists
        calibrations_fileList.update(fnUpdateProductList(baseDataFolder_reduced, productList=calibrations_fileList))

        # setup recipe "uves_cal_mflat"
        clpRecipes[recipe] = fnSetRecipe(recipe, cfgFile)

        print outputFolders[recipe]

        # run recipe "uves_cal_predict"
        clpResults[recipe] = clpRecipes[recipe](raw={recipeKeywords[recipe]['TAG']: recipeFileLists[recipe]},
                                                calib=calibrations_fileList, \
                                                output_dir=outputFolders[recipe], tmp_dir=tmpFolder)

        clpResults[recipe] = None


        fnPrintLine('CALIB', '{} complete'.format(recipe))

    # endregion --------------------------------------------------------------------------------------------------------

    # region --- 6) reduce science frame - uves_obs_scired

    recipe = 'uves_obs_scired'
    if recipe in recipes2run:
        fnPrintLine('SCIENCE', 'starting {}, please wait...'.format(recipe))

        # get input file list
        recipeFileLists[recipe] = fngetRecipeFileList(recipeName=recipe, \
                                                      dataList=rawData_fileList, \
                                                      keyword=recipeKeywords[recipe]['DPR_TYPE'])

        # update calibration lists
        calibrations_fileList.update(fnUpdateProductList(baseDataFolder_reduced, productList=calibrations_fileList))

        # setup recipe "uves_obs_scired"
        clpRecipes[recipe] = fnSetRecipe(recipe, cfgFile)

        # run recipe "uves_obs_scired"
        for scienceFile in recipeFileLists[recipe]:
            fnPrintLine('SCIENCE', '{}: reducing {}'.format(recipe, scienceFile.split('/')[-1]))
            clpResults[recipe] = clpRecipes[recipe](raw={recipeKeywords[recipe]['TAG']: [scienceFile]},
                                                    calib=calibrations_fileList, \
                                                    output_dir=os.path.join(outputFolders[recipe],
                                                                            fits.getheader(scienceFile)['DATE-OBS']),
                                                    tmp_dir=tmpFolder)

            clpResults[recipe] = None

        fnPrintLine('SCIENCE', '{} complete'.format(recipe))

    # endregion

    # region --- 7) create blaze file
    recipe = 'uves_cal_blaze'
    if recipe in recipes2run:
        fnPrintLine('BLAZE', 'starting {}, please wait...'.format(recipe))

        # get input file list
        reducedFileList = [filenameFull \
                           for root, dirs, files in os.walk(outputFolders['uves_obs_scired'], topdown=False) \
                           for filenameFull in [os.path.join(root, filename) for filename in files] \
                           if filenameFull.endswith('.fits')\
                           and 'WCALIB_FLAT'.lower() in fits.getheader(filenameFull)['HIERARCH ESO PRO CATG'].lower()\
                           ]

        for flatFile in reducedFileList:
            fnPrintLine('BLAZE', '{}: creating blaze from {}'.format(recipe, flatFile.split('/')[-1]))
            flatFits = fits.open(flatFile)

            blazeData = fnInferInstrumentResponse(flatFits[0].data,float(cfgFile.get(recipe, 'window')))
            blazeheader = flatFits[0].header
            # ii=0
            # for orderData, instData in zip(flatFits[0].data,blazeData):
            #     mplt.plot(orderData, 'bo')
            #     mplt.plot(instData, 'r')
            #     mplt.title('order %s' %ii)
            #     mplt.show()
            #     ii+=1

            blazeheader['HIERARCH ESO PRO CATG'] = flatFits[0].header['HIERARCH ESO PRO CATG'].replace('FLAT', 'BLAZE')
            blazeheader['HIERARCH BLAZE WINDOW'] = (float(cfgFile.get(recipe, 'window')), 'blaze moving average window in pixels')


            writeNameBlaze = os.path.join(*[outputFolders['uves_obs_scired'],flatFile.split('/')[-2],\
                                          'blaze_{}.fits'.format(flatFits[0].header['HIERARCH ESO PRO CATG'].split('_')[-1].lower())])


            fits.writeto(writeNameBlaze,blazeData,blazeheader, clobber=True)

        fnPrintLine('BLAZE', '{} complete'.format(recipe))




    # endregion
    return


# endregion

# region --- RunTime

if __name__ == '__main__':
    main_function()


# endregion
