#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 13 17:31:16 2014

@author: jorge
"""
import os
import sys
import numpy as np

# path definition
sys.path.append(os.path.abspath("./modules/"))
sys.path.append(os.path.abspath("/home/jmartins/WORK/Python/ReductDRS"))

# get constants
from Constants import *
#from ExoTable import fnFixMolecules, fnGetExoTable
from ExoTable import  fnGetExoTable
from bin.constants import *

# =============================================================================
# FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS
# =============================================================================


def fnComputeFlux2(m1,m2,f1):
    f2 = f1 * 10**((m1-m2)/2.5)
    return f2



# =============================================================================
def main():
    ExoplanetData = fnGetExoTable()

    PlanetList  = [planet for planet in ExoplanetData.keys() \
                    if ~np.isnan(float(ExoplanetData[planet]['mass'])) and \
                    ~np.isnan(float(ExoplanetData[planet]['semi_major_axis'])) and \
                    ~np.isnan(float(ExoplanetData[planet]['orbital_period'])) and\
                    type(ExoplanetData[planet]['mag_v']) != type('')\
                    and (massMin < float(ExoplanetData[planet]['mass']) < massMax) \
                    and (pMin < float(ExoplanetData[planet]['orbital_period']) < pMax) \
                    and (magMin < float(ExoplanetData[planet]['mag_v']) < magMax) \
                    and (decMin < float(ExoplanetData[planet]['dec']) < decMax) \
                    ]

    # Print header
    paramsString = '{0:<20} {1:<8} {2:<8} {3:<8} {4:<7} {5:<7} {6:<7} {7:<15} {8:<12} {9:<12} {10:<12} {11:<12} {12:<12}\n'
    OutputString = paramsString.replace(':<', ':#<').format('','','','','','','','','','','','','')
    OutputString += paramsString.format('name','mass', 'radius','period','* mag_v','* spec','dec', \
                                                'max flux ratio', '3 Sigma SN', '5 Sigma SN','10h SN','3 sig Ratio','5 sig Ratio')
    OutputString += paramsString.replace(':<', ':#<').format('','','','','','','','','','','','','')


    for planet in sorted(PlanetList):
        MaxFluxRatio = ALBEDO * ((float(ExoplanetData[planet]['radius'])*RJup)/(float(ExoplanetData[planet]['semi_major_axis'])*AU))**2 * np.sin(np.radians(float(ExoplanetData[planet]['inclination'])))

        ThreeSigma = 3./MaxFluxRatio
        FiveSigma = 5./MaxFluxRatio
        StarFlux = fnComputeFlux2(magRef,ExoplanetData[planet]['mag_v'],fluxRef)

        StarSN = np.sqrt(StarFlux*2.8)      # 2.8 vem da razao do resultado com o ETC

        StarSN10h = StarSN*np.sqrt(10.*3600/timeRef*nMask)

        # Print data
        OutputString += paramsString.format(planet, '%.3f' %ExoplanetData[planet]['mass'],\
            '%.3f' %ExoplanetData[planet]['radius'],'%.3f' %ExoplanetData[planet]['orbital_period'],\
            '%.2f' %ExoplanetData[planet]['mag_v'],'%s' %ExoplanetData[planet]['star_sp_type'],'%2.0f' %ExoplanetData[planet]['dec'], \
                    '%.2e' %MaxFluxRatio, '%.2e' %ThreeSigma, '%.2e' %FiveSigma, '%.2e' %StarSN10h \
                    ,'%.2f' %(StarSN10h/ThreeSigma),'%.2f' %(StarSN10h/FiveSigma))


    print OutputString
    print paramsString.replace(':<', ':#<').format('','','','','','','','','','','','','','')
    print "sample lenght:", len(PlanetList)

    with open('SelectedTargets_%s.txt' %instrument,'w') as f:
        f.write(OutputString)

# =============================================================================
# MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN
# =============================================================================

if __name__ == '__main__':
    main()



