#!/usr/bin/env python
'''

main routine
'''
__author__ = 'Jorge Martins'
__email__ = 'jorge.martins@iastro.pt'
__version__ = 'alpha 1'
__date__ = '2015/09/24'
# ======================================================================================================================
#   Python Modules
# ======================================================================================================================
import sys

import matplotlib.pyplot as mplt
# ======================================================================================================================
#   User Modules
# ======================================================================================================================

from defaults.defaultPaths import *

from recipes.manipCCF import *

from modules.fileManip import initPaths

from modules.InOut import fnPrintLine

from math_local.mathFunctions import *

# Init Terminal
os.system("cls")
os.system("clear")
fnPrintLine(None, '')
fnPrintLine(None, '', align='center', flush='=')
fnPrintLine(None, 'Create star+planet CCFs \t(version: {:<})'.format(__version__), align='center')
# fnPrintLine(None, '', align = 'center', flush = '=')
fnPrintLine(None, '')
fnPrintLine(None, 'Author:{:<} \temail: {:<}'.format(__author__, __email__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, 'Last Update:{:<}'.format(__date__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, '', align='center', flush='=')
fnPrintLine(None, '')

# ======================================================================================================================
#   Settings
# ======================================================================================================================

# Define Work Paths
WorkPath = os.path.abspath(os.getenv('WorkPath'))
ReductPath = os.path.abspath(os.getenv('ReductPath'))

sys.path.append(WorkPath)
sys.path.append(ReductPath)
settings = __import__(sys.argv[1].split('.')[0])

scienceInputFolder = os.path.abspath('{}/reduced/{}/'.format(WorkPath, settings.dataFolder))
CCFFolder = os.path.abspath('{}/CCFs/{}'.format(WorkPath, settings.dataFolder))
initPaths([CCFFolder])
scienceOuputFolder = os.path.abspath('{}/CCFs/{}/'.format(WorkPath, settings.dataFolder))


# ======================================================================================================================
#   Functions
# ======================================================================================================================

def fnCreateWave(orderHeader, order, numberOrders=1):
    '''
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


def fnCreateWave2(orderHeader, order, numberOrders=1):
    '''
    '''
    if numberOrders == 1:
        orderWavelengthStart = orderHeader['CRVAL1']
        waveData = np.linspace(orderWavelengthStart,
                               orderHeader['CRVAL1'] + (orderHeader['CDELT1'] * orderHeader['NAXIS1']),
                               num=orderHeader['NAXIS1'])
    else:
        orderWavelengthStart = orderHeader['WSTART{0}'.format(order + 1)]
        # orderWavelengthEnd       =  orderHeader['WEND{0}'.format(order+1)]
        orderStep = orderHeader['CDELT1']
        orderLen = orderHeader['NAXIS1']
        orderWavelengthEnd = orderWavelengthStart + orderStep * orderLen
        waveData = np.linspace(orderWavelengthStart, orderWavelengthEnd, num=orderLen)

    return waveData


# ======================================================================================================================
#   Main
# ======================================================================================================================

def __main__():
    fnPrintLine('CONFIG', 'Getting file list, please wait')

    if settings.specType.lower() == '1d':
        specType = 'RED_SCI_POINT'
        nOrders = 1
    if settings.specType.lower() == '2d':
        specType = 'WCALIB_SCI_POINT'

    try:
        fnPrintLine('Config', 'list of data files provided: {}'.format(settings.dataList))
        scienceFiles = np.genfromtxt('{}/{}'.format(scienceInputFolder, settings.dataList), dtype=str)

        scienceFiles = [fileName for fileName in scienceFiles \
                        if specType in fits.getheader('{}/{}'.format(scienceInputFolder, fileName))[
                            'HIERARCH ESO PRO CATG']
                        ]

    except:
        fnPrintLine('Config', 'no data file list file provided, creating CCF for all files.')
        scienceFiles = [fileName for fileName in os.listdir(scienceInputFolder) if
                        specType in fits.getheader('{}/{}'.format(scienceInputFolder, fileName))[
                            'HIERARCH ESO PRO CATG']]

    ccfProcessedList = ['{}/{}'.format(scienceOuputFolder, fileName) for fileName in os.listdir(scienceOuputFolder)]


    # getting flats

    if settings.specType.lower() == '1d':
        specTypeFlat = 'RED_FLAT_OBJ'
        nOrders = 1
    if settings.specType.lower() == '2d':
        specTypeFlat = 'WCALIB_FLAT_OBJ'

    fnPrintLine('Config', 'no data file list file provided, creating CCF for all files.')
    masterFlatFiles = [fileName for fileName in os.listdir(scienceInputFolder) if
                        specType in fits.getheader('{}/{}'.format(scienceInputFolder, fileName))[
                            'HIERARCH ESO PRO CATG']]

    ccfProcessedList = ['{}/{}'.format(scienceOuputFolder, fileName) for fileName in os.listdir(scienceOuputFolder)]



    for fitsFile in scienceFiles:
        # Getting science data
        fnPrintLine('CONFIG', 'Processing: {}  '.format(fitsFile))
        scienceFile = fits.open('%s/%s' % (scienceInputFolder, fitsFile))

        rawScienceFileKey = 'UVES.{}.{}{}'.format(scienceFile[0].header['DATE-OBS'],
                                                  scienceFile[0].header['HIERARCH ESO INS PATH'],
                                                  scienceFile[0].header['HIERARCH ESO DET OUT1 NAME'])

        ccfWriteFileName = '{}/{}_ccf.fits'.format(scienceOuputFolder, rawScienceFileKey)

        if ccfWriteFileName in ccfProcessedList:
            fnPrintLine(None, '{}_ccf.fits already been processed, skipping.'.format(rawScienceFileKey))

        else:

            try:
                nOrders = scienceFile[0].header['NAXIS2']
                fnPrintLine('CONFIG', '2D-Spectrum - Number of orders: {:>3} '.format(nOrders))
            except:
                nOrders = 1
                fnPrintLine('CONFIG', '1D-Spectrum')

            # Defining science and wave arrays
            if nOrders == 1:
                scienceData = scienceFile[0].data
                waveData = fnCreateWave2(scienceFile[0].header, 1, numberOrders=1)
                # mplt.plot(waveData,scienceData)

            elif nOrders > 1:
                scienceData = np.empty((nOrders, scienceFile[0].header['NAXIS1']))
                waveData = np.empty((nOrders, scienceFile[0].header['NAXIS1']))
                for order in np.arange(nOrders):
                    scienceData[order] = scienceFile[0].data[order]
                    waveData[order] = fnCreateWave2(scienceFile[0].header, order, numberOrders=nOrders)

            # for order in np.arange(nOrders):
            #     mplt.plot(waveData[order], scienceData[order])
            #
            # mplt.show()
            # sys.exit()

            fnPrintLine('Mask', 'Building mask')
            maskData = np.genfromtxt('{}/masks/{}.mas'.format(ReductPath, settings.maskCCF))

            # Compute CCF
            fnPrintLine('CCF', 'Computing CCF for: {}  '.format(fitsFile))

            # rangeRV = np.arange(settings.guessRV-settings.windowCCF + scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'], settings.guessRV+settings.windowCCF + scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'], settings.stepCCF)

            rangeRV = np.arange(
                settings.guessRV - settings.windowCCF - scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'],
                settings.guessRV + settings.windowCCF - scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'],
                settings.stepCCF)

            if nOrders == 1:
                rangeCCF, nLines = zip(*[fnComputeCCF(scienceData, waveData, maskData, testRV) for testRV in rangeRV])
                fnPrintLine(None, 'Computing CCF for BERV corrected RV = {:.2f}km/s'.format(rangeRV[-1]))
                mplt.plot(rangeRV, rangeCCF / max(rangeCCF))

            elif nOrders > 1:
                rangeCCF = np.empty((nOrders, len(rangeRV)))
                nLinesPerOrder = np.empty(nOrders)
                for order in np.arange(nOrders):
                    fnPrintLine('CCF', 'Computing CCF for order {:>3}/{:<3}         '.format(order + 1, nOrders))

                    rangeCCF[order], nLines = zip(
                        *[fnComputeCCF(scienceData[order], waveData[order], maskData, testRV) for testRV in rangeRV])

                    nLinesPerOrder[order] = np.nanmin(nLines)

                    fnPrintLine(None, 'Number of lines used for order {:>3}/{:<3}: {:>3.0f}'.format(order + 1, nOrders,
                                                                                                    nLinesPerOrder[
                                                                                                        order]))

            # Fit gaussian CCF

            if nOrders == 1:
                A, mean, FWHM, B = fnGaussianFitOLD(rangeRV, rangeCCF,
                                                    GaussParamsInitGuess=[max(rangeCCF) - min(rangeCCF),
                                                                          settings.guessRV - scienceFile[0].header[
                                                                              'HIERARCH ESO QC VRAD BARYCOR'], 10.,
                                                                          max(rangeCCF)])
            elif nOrders > 1:
                A = np.empty(nOrders)
                mean = np.empty(nOrders)
                FWHM = np.empty(nOrders)
                B = np.empty(nOrders)
                AErr = np.empty(nOrders)
                meanErr = np.empty(nOrders)
                FWHMErr = np.empty(nOrders)
                BErr = np.empty(nOrders)

                from math_local.mathFunctions import fnGaussianFit
                for order in np.arange(nOrders):
                    # A[order], mean[order], FWHM[order], B[order] = fnGaussianFitOLD(rangeRV,rangeCCF[order],GaussParamsInitGuess =[max(rangeCCF[order])- min(rangeCCF[order]), settings.guessRV , 10., max(rangeCCF[order])])
                    A[order], mean[order], FWHM[order], B[order], AErr[order], meanErr[order], FWHMErr[order], BErr[
                        order] = fnGaussianFit(rangeRV, rangeCCF[order],
                                               GaussParamsInitGuess=[max(rangeCCF[order]) - min(rangeCCF[order]),
                                                                     settings.guessRV - scienceFile[0].header[
                                                                         'HIERARCH ESO QC VRAD BARYCOR'], 10.,
                                                                     max(rangeCCF[order])])

            # building fits file
            CCFFits = fits.PrimaryHDU()
            CCFFits.header = scienceFile[0].header
            CCFFits.data = rangeCCF

            # NLines per order
            for order in np.arange(nOrders):
                CCFFits.header.set('HIERARCH ESO DRS CCF LINES %2.f' % order, nLinesPerOrder[order],
                                   comment='Number of lines used to compute CCF of order %2.f' % order)

            # CCF Params
            CCFFits.header.set('CDELT1', settings.stepCCF, comment='step of CCF')
            CCFFits.header.set('CRVAL1', rangeRV[0], comment='first RV of CCF / beginning of CCF')

            # Gaussian fit params
            # CCFFits.header.set('HIERARCH ESO DRS CCF RVC', mean+ scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'],comment='Computed Radial Velocity of target / mean of gaussian fit')
            # CCFFits.header.set('HIERARCH ESO DRS CCF RV', mean+ scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'],comment='Computed Radial Velocity of target / mean of gaussian fit')
            # CCFFits.header.set('HIERARCH ESO DRS CCF NOISE', meanErr ,comment='error on radial velocity of target - to be implemented')
            # CCFFits.header.set('HIERARCH ESO DRS CCF FWHM', FWHM ,comment='FWHM of fitted CCF / FWHM of gaussian fit')
            # CCFFits.header.set('HIERARCH ESO DRS CCF CONTRAST', A/B ,comment='Contrast of CCF')

            fits.writeto(ccfWriteFileName, CCFFits.data, CCFFits.header, clobber=True)

            fnPrintLine('PLOT', 'Saving image as pdf in {}'.format('./tmp/{}_ccfs.png'.format(rawScienceFileKey)))
            # mplt.ion()
            fig = mplt.figure(figsize=(100, 80), dpi=50)
            mplt.title(ccfWriteFileName.split('/')[-1])
            for order in np.arange(nOrders):
                axOrder = fig.add_subplot(int(nOrders / 4) + 1, 4, order + 1)
                axOrder.plot(rangeRV, rangeCCF[order] / B[order], 'b')
                gaussFit = fnGauss(rangeRV, [A[order], mean[order], FWHM[order], B[order]])
                axOrder.plot(rangeRV, gaussFit / B[order], 'r', )
                axOrder.annotate('Order: %s \n N_lines: %.0f \n %.2f < lambda < %.2f' % (order, nLinesPerOrder[order],
                                                                                         waveData[order][0],
                                                                                         waveData[order][-1]
                                                                                         ), xy=(0.95, 0.05),
                                 xycoords='axes fraction', va='bottom', ha='right')

            axAll = fig.add_subplot(int(nOrders / 4) + 1, 4, 4 * (int(nOrders / 4) + 1))
            wholeCCF = np.nansum(rangeCCF, axis=0)
            axAll.plot(rangeRV, wholeCCF / max(wholeCCF), 'g')
            axAll.annotate('All orders stacked \n N_lines: %.0f \n %.2f < lambda < %.2f' % (
            np.nansum(nLinesPerOrder), waveData[0][0],
            waveData[-1][-1]), xy=(0.95, 0.05), xycoords='axes fraction', va='bottom', ha='right')

            # mplt.draw()

            fig.savefig('./tmp/{}_ccfs.png'.format(rawScienceFileKey), format='png')
            mplt.close(fig)


# ======================================================================================================================
#   RunTime
# ======================================================================================================================
if __name__ == '__main__':
    __main__()
