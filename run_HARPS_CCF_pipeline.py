#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

########################################################################
# 							DESCRIPTION								   #
########################################################################
# python modules
import time
import os
import sys
import ConfigParser

# path definition
#sys.path.append(os.path.abspath("./modules/"))
#settings = __import__(sys.argv[1].split('.')[0])

# user modules
from modules.InOut import fnPrintLine
from modules.manip_data import *
# Initializing terminal
if os.name == 'nt': os.system('cls')
if os.name == 'posix': os.system('clear')
rows, columns = 80,120
print '_'*int(columns)+'\n\n'+("Starting "+os.path.basename(__file__)).center(int(columns)).upper()

start_time = time.time()


# Functions and classes
class clsSettings:
    def __init__(self, plots = False, data = 'data',ccfMask = None, ccfWidth = None, ccfStep = None, theoRV = None ,  **kwargs):
        self.data = data
        self.ccfMask = ccfMask
        self.ccfWidth = ccfWidth
        self.ccfStep = ccfStep
        self.theoRV = theoRV




########################################################################
# 							DEFINITIONS								   #
#########################################################################
strInstruction = 'off_make_ccf_harps.py  {night} {e2ds} {ccfMask} {theoRV} {ccfWidth} {ccfStep} | grep off_make_ccf_harps'

#############
# FUNCTIONS #
#############
def run_bash_script(lista,settingsFile):
    fnPrintLine('CCF','Processing folder {} - {} files '.format(str(settingsFile.data),len(lista) ))
    for fileName, index in zip(lista, range(len(lista))):
        os.system('clear')
        fnPrintLine('{:>5.0f}/{:<5.0f}'.format(index+1,len(lista)),'Processing {}'.format(fileName,len(lista) ))
        instruction = strInstruction.format(night=str(settingsFile.data) ,e2ds=str(fileName) , ccfMask=str(settingsFile.ccfMask) , theoRV=str(settingsFile.theoRV) , ccfWidth=str(settingsFile.ccfWidth) , ccfStep=str(settingsFile.ccfStep) )
        os.system(instruction)
    	
########
# MAIN #
########
def main_function():
	cfgFileName = sys.argv[1]
	cfgFile = ConfigParser.RawConfigParser(allow_no_value=True)
	cfgFile.optionxform = str
	cfgFile.read(cfgFileName)
	print {item[0]: item[1] for section in ['global', 'ccf'] for item in cfgFile.items(section)}
	settings = clsSettings(**{item[0]: item[1] for section in ['global', 'ccf'] for item in cfgFile.items(section)})
	lista = [fileName for fileName in os.listdir('./reduced/{}'.format(settings.data)) if 'e2ds' in fileName]
	create_folder('./CCFs/')
	create_folder('./CCFs/{}'.format(str(settings.data)))
	run_bash_script(lista,settings)
	os.system("mv -vf ./reduced/{0}/*ccf*fits ./CCFs/{0}_CCFs".format(str(settings.data)))
	print 'CCF files have been moved to ./CCFs/{}_CCFs/'.format(str(settings.data))
	print '\n'+'_'*int(columns)+'\n'
	print ("Program ran for "+ str(int (time.time() - start_time)) + " seconds").center(int(columns))
	print_head_text("That's All Folks!!!")



# region --- RunTime

if __name__ == '__main__':
    main_function()
