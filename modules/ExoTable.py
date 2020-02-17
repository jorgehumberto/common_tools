#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 13 17:31:16 2014

@author: jorge
"""
import os
import sys
import numpy as np
import gzip
from PyAstronomy import pyasl
from PyAstronomy.pyaC import pyaPermanent
#import matplotlib.pyplot as mplt
import re
# path definition
#sys.path.append(os.path.abspath("./modules/"))

# get constants
from Constants import *

# get planet parameters
#from PlanetDefaults import *
#from PlanetSettings import *


# =============================================================================
# FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS
# =============================================================================
def fnFixMolecules(line):
    '''
        Code from Joao Faria to replace the commas in the molecules
    '''
    r = re.compile(r'''
        \s*                 # Any whitespace.
        (                   # Start capturing here.
          "(?:              # A double-quote followed by a string of characters...
              [^"\\]|\\.    # That are either non-quotes or escaped...
           )*               # ...repeated any number of times.
          "                 # Followed by a closing double-quote.
          |                 # OR
          '(?:[^'\\]|\\.)*' # Same as above, for single quotes.
        )                   # Done capturing.
        \s*                 # Allow arbitrary space before the comma.
        (?:,|$)             # Followed by a comma or the end of a string.
    ''', re.VERBOSE)
    # find all occurences of thingy
    matches = r.findall(line)
    # replace them with sane things
    for match in matches:
        match_ok = match.replace(',', ';')  # replace commas with semi-colons inside each thingy
        line = line.replace(match, match_ok, 1)  # replace in string

    return line


def fnGetExoTable(inclination = 90, **kwargs):
	ExoplanetEU = pyasl.ExoplanetEU()
# 	ExoplanetEU.forceUpdate()
	filename = '%s/%s' %(pyaPermanent.PyAConfig().getDataRoot(),ExoplanetEU.dataFileName)


	with gzip.open(filename,'rb') as File:
	    #Labels = File.readline().replace('\n' , '').replace('\r' , '').replace('#' , '').replace(' ' , '').split(',')
	    Labels = File.readline().replace('\n' , '').replace('\r' , '').replace('#' , '').split(',')
	    Data = {}
	    for line in File.readlines():
	        line = fnFixMolecules(line)
	        #LineData = line.replace('\n' , '').replace('\r' , '').replace('#' , '').replace(' ' , '').split(',')
	        LineData = line.replace('\n' , '').replace('\r' , '').replace('#' , '').split(',')
	        Data[LineData[0]] = {}
	        for label, value in zip(Labels[1:], LineData[1:]):
	            if label == 'radius' and value== '':
	                # if Data[LineData[0]]['mass'] <= 0.035:
	                #     value = .23
	                # elif 0.035 < Data[LineData[0]]['mass'] <= 0.1:
	                #     value = .33
	                # elif 0.1 < Data[LineData[0]]['mass'] <= 0.33:
	                #     value = 1.
	                # else:
	                #     value = RADIUS
					True
	            elif label == 'inclination' and value== '':
	                value = inclination
	            elif value == '':
	                value = np.nan

	            try:
	                Data[LineData[0]][label] = float(value)
	            except:
	                Data[LineData[0]][label] = value

#                if '51' in LineData[0]:
#                    print label, LineData[0],Data[LineData[0]][label], type(Data[LineData[0]][label])
#
#            if '51' in LineData[0]:
#                raw_input()
#

        return Data
