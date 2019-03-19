#!/usr/bin/env python
# -*- coding: utf-8 -*-
########################################################################
#							PYHTON MODULES							   #
########################################################################
import numpy as np


########################################################################
#							    FUNCTIONS		    				   #
########################################################################

def fnRVStarOrbitCircular(PlanetParams, PhaseZero, Phase):
    RVStar = PlanetParams.k1 * np.sin(2. * np.pi * (Phase + PhaseZero)) * np.sin(np.radians(PlanetParams.I))

    return RVStar


#  -----------------------------------------------------------------------------

def fnRVStarOrbitElipse(PlanetParams, Phase):
    RVPlanet = PlanetParams.k1 * (np.cos(2 * np.pi * (Phase) + np.radians(PlanetParams.w)) \
                                  + PlanetParams.ecc * np.cos(np.radians(PlanetParams.w))) + PlanetParams.sysRV
    return RVPlanet




#  -----------------------------------------------------------------------------
def fnAlpha(Inclination, Phase):
    Alpha = np.arccos(-np.sin(np.radians(Inclination)) * np.cos(2 * np.pi * Phase))
    return Alpha


#  -----------------------------------------------------------------------------
def fnPhaseFunction(Inclination, Phase, model='None'):
    """Inclination in degrees
    phase: 0-1
    """
    Alpha = fnAlpha(Inclination, Phase)

    if model == 'Lambert':
        g = (np.sin(Alpha) + ((np.pi - Alpha) * np.cos(Alpha))) / np.pi
    elif model == 'Venus':
        dm = 0.09 * (np.degrees(Alpha) / 100) + 2.39 * (np.degrees(Alpha) / 100) ** 2 - .65 ** (np.degrees(
            Alpha) / 100) ** 3
        g = 10 ** (-.4 * dm)
    else:
        g = 1.

    return g
