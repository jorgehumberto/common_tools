

import ConfigParser
import os
import sys

import matplotlib.pyplot as mplt
import numpy as np
from PyAstronomy.pyTiming import pyPeriod
from PyAstronomy.pyasl.asl.astroTimeLegacy import daycnv
from matplotlib import cm
from matplotlib import gridspec
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import ephem


# region --- 2D plot of CCFs
def fnPlot2DCCFs(CCFs, plotList, stellarCCFWidth = 2, cols2DCCFsPlot = 32, title = 'None', minFlux = None, maxFlux = None,RVRanges = None):


    nightList = sorted(set([int(CCFs[ccfName].BJD) for ccfName in plotList]))
    # define figure
    verticalSpace = 1
    # axes limits
    ax2DCCFsMin = 5
    axPhasesMin = 16
    axBERVMin = 20
    axAirMassMin = 24
    axSN50Min = 28

    # define figure
    figure2DCCFs = mplt.figure(figsize=(20,int(np.ceil((len(plotList)                   # plot space
                                                        + verticalSpace*(len(nightList)+3)  # vertical spaces between plots
                                                        + 10                            # header
                                                        )/5))), dpi=100)

    # delete central region of star CCF for better contrast
    for ccfName in sorted(plotList):
        CCFs[ccfName].data[CCFs[ccfName].ccfMeanPixel - stellarCCFWidth * CCFs[ccfName].ccfFWHMPixels:\
            CCFs[ccfName].ccfMeanPixel + stellarCCFWidth *CCFs[ccfName].ccfFWHMPixels] = None

    # figure2DCCFs.suptitle(title, fontsize=14)

    # Define gridspec axes
    # rows = int(np.ceil(len(plotList)))
    gsAxes = gridspec.GridSpec(len(plotList)                   # plot space
                                + verticalSpace*(len(nightList)+3)  # vertical spaces between plots
                                + 10, cols2DCCFsPlot)
    ax2DCCFs = range(len(nightList))
    axPhases = range(len(nightList))
    axNotes = range(len(nightList))
    axBERV = range(len(nightList))
    axAirMass = range(len(nightList))
    axSN50 = range(len(nightList))

    # imshow limits
    if maxFlux == None:
        maxFlux = np.nanmax([np.nanmax((CCFs[ccfName].data / np.nanmedian(CCFs[ccfName].data) - 1) * 1000) \
                         for ccfName in plotList])
    if minFlux == None:
        minFlux = np.nanmin([np.nanmin((CCFs[ccfName].data / np.nanmedian(CCFs[ccfName].data) - 1) * 1000) \
                         for ccfName in plotList])

    iniGS = 3


    for n, night in zip(np.arange(len(nightList)), nightList):
        ccfListNight = sorted([ccfName for ccfName in sorted(plotList) \
                      if int(CCFs[ccfName].BJD) == night])


        # Build axes
        axNotes[n]= figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), :ax2DCCFsMin])
        axNotes[n].axis('off')

        if n == 0:
            ax2DCCFs[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), ax2DCCFsMin:axPhasesMin])
            axPhases[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), axPhasesMin:axBERVMin])
            axBERV[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), axBERVMin:axAirMassMin])
            axAirMass[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), axAirMassMin:axSN50Min])
            axSN50[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), axSN50Min:])

        else:
            ax2DCCFs[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), ax2DCCFsMin:axPhasesMin], sharex=ax2DCCFs[0])
            axPhases[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), axPhasesMin:axBERVMin], sharex=axPhases[0])
            axBERV[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), axBERVMin:axAirMassMin], sharex=axBERV[0])
            axAirMass[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), axAirMassMin:axSN50Min], sharex=axAirMass[0])
            axSN50[n] = figure2DCCFs.add_subplot(gsAxes[iniGS :iniGS + len(ccfListNight), axSN50Min:], sharex=axSN50[0])

        iniGS += len(ccfListNight) + verticalSpace

        ax2DCCFs[n].ticklabel_format(useOffset=False, axis='y')

        # configure ticks
        ax2DCCFs[n].yaxis.set_major_locator(MaxNLocator(len(ccfListNight)))
        ax2DCCFs[n].yaxis.grid(b=True, which='major', color='k', linestyle='-')
        ax2DCCFs[n].set_yticklabels(['{:.3f}'.format(CCFs[ccfName].BJD - night) for ccfName in ccfListNight ])

        indexes = np.linspace(0.5,len(ccfListNight)-.5, len(ccfListNight))
        # 2d normalized CCF
        image = ax2DCCFs[n].imshow( \
            np.array([(CCFs[ccfName].data / np.nanmedian(CCFs[ccfName].data) - 1) * 1000 \
                      for ccfName in ccfListNight[::-1]]),\
                     aspect='auto', cmap=cm.Greys_r,
                     extent=(np.nanmin(CCFs[CCFs.keys()[0]].wave), np.nanmax(CCFs[CCFs.keys()[0]].wave), \
                             0,
                             len(ccfListNight)), \
                     vmin=minFlux, vmax=maxFlux \
                     )
        ax2DCCFs[n].set_xlim(np.nanmin(CCFs[CCFs.keys()[0]].wave), np.nanmax(CCFs[CCFs.keys()[0]].wave))
        # RVs

        for ccfName, ii in zip(ccfListNight, indexes):
            ax2DCCFs[n].plot(CCFs[ccfName].planetRV, ii, 'r*')

        if RVRanges != None:

            for ccfName, ii in zip(ccfListNight, indexes):
                ax2DCCFs[n].add_patch(
                    patches.Rectangle(
                            (RVRanges[ccfName][0], ii-.5),  # (x,y)
                            RVRanges[ccfName][-1] - RVRanges[ccfName][0],  # width
                            1.,  # height
                            alpha = 0.3,  # alpha
                            facecolor = 'green',   # color
                            snap = False,        # dunno
                            )
            )


        # phases
        axPhases[n].plot([CCFs[ccfName].PlanetPhaseFolded for ccfName in ccfListNight], indexes, 'r+')
        # BERVs
        axBERV[n].plot([CCFs[ccfName].BERV for ccfName in ccfListNight], indexes, 'r+')
        # AirMAss
        axAirMass[n].plot([CCFs[ccfName].AirMass for ccfName in ccfListNight], indexes, 'r+')
        # SN50
        axSN50[n].plot([CCFs[ccfName].SN50 for ccfName in ccfListNight], indexes, 'r+')
        # Notes

        Moon = ephem.Moon()
        moonDist = {}
        daysToFullMoon = {}
        for ccfName in plotList:
            Moon.compute(CCFs[ccfName].BJD, epoch='2000')
            # compute moon distance
            coordsMoon = (Moon.a_ra, Moon.a_dec)
            coordsStar = (CCFs[ccfName].ra, CCFs[ccfName].dec)

            moonDist[ccfName] = ephem.degrees(ephem.separation(coordsStar , coordsMoon))
            # compute days to next full moon
            daysToFullMoon[ccfName] = ephem.next_full_moon(CCFs[ccfName].BJD) - CCFs[ccfName].BJD


        medianMoonDistance = np.median([np.degrees(moonDist[ccfName]) for ccfName in ccfListNight])
        medianDaysToFullMoon = np.median([daysToFullMoon[ccfName] for ccfName in ccfListNight])

        strNotes = 'Moon Dist:{:.1f} deg \nDays to Full Moon: {:.1f}\n'.format(medianMoonDistance,medianDaysToFullMoon)

        if min([CCFs[ccfName].BJD - night for ccfName in ccfListNight ]) <= .5:
            strNotes += '{} spectra \n date:{} \n JD:{}+'.format(len(ccfListNight),'{}-{}-{}'.format(*daycnv(night)[0:3]), night)
        else:
            strNotes += '{} spectra \n date:{} \n JD:{}+'.format(len(ccfListNight),'{}-{}-{}'.format(*daycnv(night+1)[0:3]), night)



        axNotes[n].text(.75, .5, strNotes,\
                        fontsize=12, horizontalalignment='right',  verticalalignment='center', \
                        transform=axNotes[n].transAxes \
                   )



        # adjust axes
        for ax in [axPhases, axBERV, axAirMass, axSN50]:
            ax[n].set_ylim(0,len(ccfListNight))
            ax[n].yaxis.set_major_locator(MaxNLocator(len(ccfListNight)))
            ax[n].yaxis.grid(b=True, which='major', color='k', linestyle='-')
            ax[n].set_yticklabels([])
            ax[n].yaxis.tick_right()
            ax[n].xaxis.set_major_locator(MaxNLocator(4,prune='both'))


    # color bar
    ax2DCCFsColorbar = figure2DCCFs.add_subplot(gsAxes[2:3,ax2DCCFsMin:axPhasesMin])
    ax2DCCFsColorbar.ticklabel_format(useOffset=False)
    figure2DCCFs.colorbar(image, cax=ax2DCCFsColorbar,orientation="horizontal")

    # Titles
    ax2DCCFsTitle = figure2DCCFs.add_subplot(gsAxes[0:1,4:-4])
    ax2DCCFsTitle.axis('off')
    ax2DCCFsTitle.text(.5, .95, title,fontsize=36, \
                    horizontalalignment='center',  verticalalignment='center'\
                   )

    axPhasesTitle = figure2DCCFs.add_subplot(gsAxes[2:3,axPhasesMin:axBERVMin])
    axPhasesTitle.axis('off')
    axPhasesTitle.text(.5, .5, 'Orbital Phases',fontsize=14, \
                    horizontalalignment='center',  verticalalignment='center'\
                   )

    axBERVTitle = figure2DCCFs.add_subplot(gsAxes[2:3,axBERVMin:axAirMassMin])
    axBERVTitle.axis('off')
    axBERVTitle.text(.5, .5, 'BERV [km/s]',fontsize=14, \
                    horizontalalignment='center',  verticalalignment='center'\
                   )

    axAirMassTitle = figure2DCCFs.add_subplot(gsAxes[2:3,axAirMassMin:axSN50Min])
    axAirMassTitle.axis('off')
    axAirMassTitle.text(.5, .5, 'AirMass',fontsize=14, \
                    horizontalalignment='center',  verticalalignment='center'\
                        )

    axSN50Title = figure2DCCFs.add_subplot(gsAxes[2:3,axSN50Min:])
    axSN50Title.axis('off')
    axSN50Title.text(.5, .5, 'SN50',fontsize=14, \
                    horizontalalignment='center',  verticalalignment='center'\
                        )

    # adjusting figure
    figure2DCCFs.tight_layout()
    figure2DCCFs.subplots_adjust(top=0.85, hspace = 0.5, wspace = 1.)

    return figure2DCCFs



# endregion