#!/usr/bin/env python
'''

main routine
'''
__author__ = 'Jorge Martins'
__email__ = 'jorge.martins@iastro.pt'
__version__ = 'alpha 1'
__date__ = '2015/09/24'

# region --- Python Modules
import os
import sys
# endregion

# region --- User Modules
from modules.InOut import fnPrintLine

# endregion

# region --- Init Terminal
os.system("cls")
os.system("clear")
fnPrintLine(None, '')
fnPrintLine(None, '', align='center', flush='=')
fnPrintLine(None, 'Recover planet CCFs \t(version: {:<})'.format(__version__), align='center')
# fnPrintLine(None, '', align = 'center', flush = '=')
fnPrintLine(None, '')
fnPrintLine(None, 'Author:{:<} \temail: {:<}'.format(__author__, __email__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, 'Last Update:{:<}'.format(__date__), align='center')
fnPrintLine(None, '')
fnPrintLine(None, '', align='center', flush='=')
fnPrintLine(None, '')

# endregion

# Paths
WorkPath = os.path.abspath(os.getenv('WorkPath'))
ReductPath = os.path.abspath(os.getenv('ReductPath'))

sys.path.append(WorkPath)
sys.path.append(ReductPath)

# endregion

# settings = __import__(sys.argv[1].split('.')[0])
# scienceInputFolder = os.path.abspath('{}/CCFs/{}'.format(WorkPath,settings.dataFolder))
# scienceOuputFolder = os.path.abspath('{}/planetResults/{}/'.format(WorkPath,settings.dataFolder))
# initPaths([scienceOuputFolder])

# Classes

# Models

# Input/Output

# Math Functions

# CCF Manipulation
