#!/usr/bin/env python
'''

main routine
'''
__author__ = 'Jorge Martins'
__email__ = 'jorge.martins@iastro.pt'
__version__ = 'alpha 2'
__date__ = '2016/02/18'

# region --- INITIALIZATION
# region --- Python Modules
import ConfigParser
import os
import sys
import matplotlib.pyplot as mplt
import scipy.optimize as opt

# endregion

# region --- User Modules
from defaults.defaultPaths import *
from recipes.manipCCF import *
from modules.InOut import fnPrintLine
from math_local.mathFunctions import fnGauss, fnGaussFit
from matplotlib.backends.backend_pdf import PdfPages
# endregion

# region --- Parse Settings
cfgFileName = sys.argv[1]
cfgFile = ConfigParser.RawConfigParser(allow_no_value=True)
cfgFile.optionxform = str
cfgFile.read(cfgFileName)
# endregion

# region --- path, folder and filenames definition
WorkPath = os.path.abspath(os.getenv('WorkPath'))
ReductPath = os.path.abspath(os.getenv('ReductPath'))
sys.path.append(WorkPath)
sys.path.append(ReductPath)

scienceInputFolder = os.path.join(*[WorkPath, 'reduced', cfgFile.get('global', 'dataFolder'), 'uves_obs_scired'])
scienceOuputFolder = os.path.join(*[WorkPath, 'CCFs', cfgFile.get('global', 'dataFolder')])
if not os.path.exists(scienceOuputFolder):
    os.makedirs(scienceOuputFolder)


# region --- Function fnGaussian
def fnGaussian(dataXXX, Amp, mean, FWHM, B):
    '''
    Defines an inverted gaussian curve as
            Y = B - Amp * e^((-4 log(2) * ((X - mean))/ FWHM)**2)
    :param dataXXX: X-axis data
    :param Amp: amplitude of gaussian curve
    :param mean: mean of gaussian curve
    :param FWHM: Full-Width Half-Maximum of gaussian curve
    :param B:fnProgress(
    return m * dataXXX + b
    '''
    return B - Amp * np.exp(-4 * np.log(2) * (((dataXXX - mean) / FWHM) ** 2))

# endregion

# endregion


# region --- Init Terminal
os.system('cls' if os.name == 'nt' else 'clear')
os.system("cls")
os.system("clear")
fnPrintLine(None, '')
fnPrintLine(None, '', align='center', flush='=')
fnPrintLine(None, 'Create star+planet CCFs \t(version: {:<})'.format(__version__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, 'Author:{:<} \temail: {:<}'.format(__author__, __email__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, 'Last Update:{:<}'.format(__date__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, '', align='center', flush='=')
fnPrintLine(None, '')
# endregion

# ======================================================================================================================
#   Settings
# ======================================================================================================================

# sys.path.append(WorkPath)
# sys.path.append(ReductPath)
# settings = __import__(sys.argv[1].split('.')[0])




# ======================================================================================================================
#   Functions
# ======================================================================================================================









# ======================================================================================================================
#   Main
# ======================================================================================================================

def __main__():
    # region -- getting file lists
    fnPrintLine('CONFIG', 'Getting file list, please wait')

    specType = 'WCALIB_SCI'

    try:
        fnPrintLine('Config', 'list of data files provided: {}'.format(settings.dataList))
        scienceFiles = np.genfromtxt('{}/{}'.format(scienceInputFolder, settings.dataList), dtype=str)

        scienceFiles = [fileName for fileName in scienceFiles \
                        if specType in fits.getheader('{}/{}'.format(scienceInputFolder, fileName))[
                            'HIERARCH ESO PRO CATG']
                        ]

    except:
        fnPrintLine('Config', 'no data file list file provided, creating CCF for all files.')
        scienceFiles = ['{}/{}'.format(root, fileName) for root, dir, fileNames in
                    os.walk(scienceInputFolder, topdown=False) \
                    for fileName in fileNames \
                    if specType in fits.getheader('{}/{}'.format(root, fileName))['HIERARCH ESO PRO CATG'] \
                    ]

    ccfProcessedList = ['{}/{}'.format(scienceOuputFolder, fileName) for fileName in os.listdir(scienceOuputFolder)]

    # endregion


    # region -- processing files
    for fitsFile in scienceFiles:

        fnPrintLine('CONFIG', 'Processing: {}  '.format(fitsFile))

        # region -- building response curves per order
        blazeFile = fits.open(fitsFile.replace('resampled_science_', 'blaze_'))
        instrumentResponse = blazeFile[0].data

        # for instData in instrumentResponse:
        #     mplt.plot(instData, 'bo')
        #     mplt.show()
        # sys.exit()

        # endregion

        # region -- getting science files
        scienceFile = fits.open(fitsFile)
        # endregion

        # region -- setting write names
        rawScienceFileKey = 'UVES.{}.{}{}'.format(scienceFile[0].header['DATE-OBS'],
                                                  scienceFile[0].header['HIERARCH ESO INS PATH'],
                                                  scienceFile[0].header['HIERARCH ESO DET OUT1 NAME'])
        ccfWriteFileName = '{}/{}_ccf.fits'.format(scienceOuputFolder, rawScienceFileKey)
        # endregion

        # region -- Processing science files
        if ccfWriteFileName in ccfProcessedList:
            fnPrintLine(None, '{}_ccf.fits already been processed, skipping.'.format(rawScienceFileKey))

        else:


            nOrders = scienceFile[0].header['NAXIS2']
            fnPrintLine('CONFIG', '2D-Spectrum - Number of orders: {:>3} '.format(nOrders))

            # region - Defining science and wave arrays
            scienceData = np.empty((nOrders, scienceFile[0].header['NAXIS1']))
            waveData = np.empty((nOrders, scienceFile[0].header['NAXIS1']))
            for order in np.arange(nOrders):
                scienceData[order] = scienceFile[0].data[order]*instrumentResponse[order]
                waveData[order] = fnCreateWave2(scienceFile[0].header, order)

            fnPrintLine('Mask', 'Building mask')
            maskData = np.genfromtxt('{}/masks/{}.mas'.format(ReductPath, cfgFile.get('ccf', 'maskCCF')))
            # endregion

            # region - Compute CCF
            fnPrintLine('CCF', 'Computing CCF for: {}  '.format(fitsFile))
            rangeRV = np.arange(
                float(cfgFile.get('ccf', 'guessRV')) - float(cfgFile.get('ccf', 'windowCCF')) - scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'],
                float(cfgFile.get('ccf', 'guessRV')) + float(cfgFile.get('ccf', 'windowCCF')) - scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'],
                float(cfgFile.get('ccf', 'stepCCF'))
            )

            rangeCCF = np.empty((nOrders, len(rangeRV)))
            nLinesPerOrder = np.empty(nOrders)
            for order in np.arange(nOrders):
                fnPrintLine('CCF', 'Computing CCF for order {:>3}/{:<3}         '.format(order + 1, nOrders))
                nonZeroRange = np.nonzero(scienceData[order])
                rangeCCF[order], nLines = zip(*[fnComputeCCF(scienceData[order][nonZeroRange], waveData[order][nonZeroRange], maskData, testRV) for testRV in rangeRV])
                nLinesPerOrder[order] = np.nanmin(nLines)
                fnPrintLine(None, 'Number of lines used for order {:>3}/{:<3}: {:>3.0f}'.format(\
                    order + 1, nOrders,nLinesPerOrder[order]))

            rangeRV += scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR']

            rangeCCF = np.vstack((rangeCCF, np.nansum(rangeCCF, axis=0)))


            # endregion

            # region - writing CCFs to fits files
            CCFFits = fits.PrimaryHDU()
            CCFFits.header = scienceFile[0].header
            CCFFits.data = rangeCCF

            # NLines per order
            for order in np.arange(nOrders):
                CCFFits.header.set('HIERARCH ESO DRS CCF LINES %2.f' % order, nLinesPerOrder[order],
                                   comment='Num. lines to compute CCF for order %2.f' % order)

            # CCF Params
            CCFFits.header.set('CDELT1', float(cfgFile.get('ccf', 'stepCCF')), comment='step of CCF')
            CCFFits.header.set('CRVAL1', rangeRV[0], comment='first RV of CCF / beginning of CCF / BERV corrected')



            # fit gauss to summed CCF
            fnPrintLine('CONFIG', 'Fitting Gaussian to stacked ordersCCF'.format(fitsFile))
            (A, mean, FWHM, B), pcov = opt.curve_fit(fnGaussian, rangeRV, rangeCCF[-1],
                                                   p0=(np.nanmax(rangeCCF[-1]) - np.nanmin(rangeCCF[-1]),float(cfgFile.get('ccf', 'guessRV')),10.,np.nanmax(rangeCCF[-1])),
                                                #bounds=fitGaussBounds,
                                                   )

            #print cov#
            (A_Err, mean_Err, FWHM_Err, B_Err) = np.sqrt(np.diag(pcov))

            # Gaussian fit params
            # CCFFits.header.set('HIERARCH ESO DRS CCF RVC', mean+ scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'],comment='Computed Radial Velocity of target / mean of gaussian fit')
            CCFFits.header.set('HIERARCH ESO DRS CCF RV', mean + scienceFile[0].header['HIERARCH ESO QC VRAD BARYCOR'],comment='Computed Radial Velocity of target / mean of gaussian fit')
            CCFFits.header.set('HIERARCH ESO DRS CCF RV ERR', mean_Err ,comment='error on radial velocity of target - to be revised')
            CCFFits.header.set('HIERARCH ESO DRS CCF FWHM', FWHM ,comment='FWHM of fitted CCF / FWHM of gaussian fit')
            CCFFits.header.set('HIERARCH ESO DRS CCF FWHM ERR', FWHM_Err ,comment='FWHM of fitted CCF / FWHM of gaussian fit - to be revised')
            CCFFits.header.set('HIERARCH ESO DRS CCF CONTRAST', A/B ,comment='Contrast of CCF')

            fits.writeto(ccfWriteFileName, CCFFits.data, CCFFits.header, clobber=True)
            # endregion

            # region - plotting CCFs as pdf
            fnPrintLine('PLOT', 'Saving image as pdf in {}'.format('./tmp/{}_ccfs.pdf'.format(rawScienceFileKey)))
            # define multipage pdf file
            pdfFile = PdfPages('./tmp/{}_ccfs.pdf'.format(rawScienceFileKey))
            # define figures
            # figCCFs = mplt.figure(figsize=(50, 40), dpi=100)
            # figSpectra = mplt.figure(figsize=(50, 40), dpi=100)
            # # define titles
            # figCCFs.suptitle('CCF name key:{}\n'.format(rawScienceFileKey))
            # figSpectra.suptitle('spectrum name key:{}\n'.format(rawScienceFileKey))

            # setup axes and plot data
            for order in np.arange(nOrders):
                figOrder, [axSpectrum, axCCFs] = mplt.subplots(2,1,figsize=(20, 15), dpi=300)

                # axOrderCCFs = figCCFs.add_subplot(int(nOrders / 4) + 1, 4, order + 1)
                # axOrderSpectra = figSpectra.add_subplot(int(nOrders / 4) + 1, 4, order + 1)
                # CCFs ax
                axCCFs.plot(rangeRV, rangeCCF[order] / np.nanmedian(rangeCCF[order]), 'b')
                axCCFs.set_xlabel('RV [km/s]', fontsize = 20)
                axCCFs.set_xlim(np.nanmin(rangeRV), np.nanmax(rangeRV))

                # spectra ax
                axSpectrum.plot(waveData[order][100:-100],scienceData[order][100:-100]/np.nanmax(scienceData[order][100:-100]), 'b')
                axSpectrum.set_xlabel('wavelength [$\AA$]', fontsize = 20)
                axSpectrum.set_xlim(np.nanmin(waveData[order][100:-100]), np.nanmax(waveData[order][100:-100]))

                # title
                figOrder.suptitle(r'file name: %s' %rawScienceFileKey + '\n' +  \
                             r'Order: %s     N_lines: %.0f     %.2f $\AA < \lambda < %.2f \AA$' \
                                     % (order, nLinesPerOrder[order],waveData[order][0], waveData[order][-1])\
                             , fontsize = 30)
                figOrder.tight_layout()
                figOrder.subplots_adjust(top = .9)
                pdfFile.savefig(figOrder)
                mplt.close(figOrder)

            # axAll = figCCFs.add_subplot(int(nOrders / 4) + 1, 4, 4 * (int(nOrders / 4) + 1))
            wholeCCF = rangeCCF[-1]#np.nansum(rangeCCF, axis=0)
            figAll, axAll = mplt.subplots(1,1,figsize=(20, 15), dpi=300)
            axAll.plot(rangeRV, wholeCCF / max(wholeCCF), 'g')
            axAll.set_xlabel('RV [km/s]', fontsize = 20)
            axAll.set_xlim(np.nanmin(rangeRV), np.nanmax(rangeRV))
            axAll.text(0.05, 0.05,' Contrast: {:.2f}\n RV: {:.2f} +-{:.2f} km/s\n FWHM: {:.2f} +-{:.2f} km/s'.format(A/B, mean, mean_Err, FWHM, FWHM_Err),
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            transform = axAll.transAxes,fontsize = 40)
            figAll.suptitle(r'file name: %s' %rawScienceFileKey + '\n' +  \
                         'All orders stacked: N_lines: %.0f  $%.2f \AA < \lambda < %.2f \AA$' \
                         % (np.nansum(nLinesPerOrder), waveData[0][0], waveData[-1][-1]),\
                         fontsize = 30)
            figAll.tight_layout()
            figAll.subplots_adjust(top = .9)

            pdfFile.savefig(figAll)
            # figCCFs.tight_layout()
            # figSpectra.tight_layout()
            # pdfFile.savefig(figSpectra)
            # pdfFile.savefig(figCCFs)
            pdfFile.close()
            mplt.close("all")


            # endregion

    # endregion


# ======================================================================================================================
#   RunTime
# ======================================================================================================================
if __name__ == '__main__':
    __main__()
