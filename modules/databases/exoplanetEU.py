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

#-------------------------------------------------------------------------
def fnGetExoTable():
	ExoplanetEU = pyasl.ExoplanetEU()
# 	ExoplanetEU.forceUpdate()
	filename = os.path.join(pyaPermanent.PyAConfig().getDataRoot(),ExoplanetEU.dataFileName)


	with gzip.open(filename,'rb') as File:
	    Labels = File.readline().replace('\n' , '').replace('\r' , '').replace('#' , '').split(',')
	    Data = {}
	    for line in File.readlines():
	        line = fnFixMolecules(line)
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
					a = 1
	            elif label == 'inclination' and value== '':
	                value = INCLINATION
	            elif value == '':
	                value = np.nan

	            try:
	                Data[LineData[0]][label] = float(value)
	            except:
	                Data[LineData[0]][label] = value


        return Data
