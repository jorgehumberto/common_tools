#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
########################################################################
#							PYHTON MODULES							   #
########################################################################
import re, os, sys, textwrap
from datetime import datetime
from tqdm import tqdm
import subprocess

########################################################################
#							INPUT FUNCTIONS							   #
########################################################################



########################################################################
#							OUTPUT FUNCTIONS						   #
########################################################################

def get_terminal_width():
    '''
    from http://www.brandonrubin.me/2014/03/18/python-snippet-get-terminal-width/
    :return: width
    '''
    try:
        width = int(subprocess.check_output(['tput', 'cols']))
    except OSError as e:
        print("Invalid Command '{0}': exit status ({1})".format(
            command[0], e.errno))
    except subprocess.CalledProcessError as e:
        print("Command '{0}' returned non-zero exit status: ({1})".format(
            command, e.returncode))
    else:
        return width


def fnPrintLine(tag, msg, cols=None, sameLine=False, align='left', flush='', full=False):
    """
    prints a formated line with a tag, message and time to the screen:
    [   TAG    ] This is a message....................................... [ 22:36:39 ]
    """
    if align == 'center':
        halign = '^'
    elif align == 'right':
        halign = '>'
    else:
        halign = '<'

    if cols == None:
        try:
            cols = get_terminal_width()
            if cols < 80:
                raise
        except:
            cols = 100

    if len(msg) > cols - 34:
        msg = textwrap.wrap(msg, width=cols - 34)
        if tag == None:
            string = '{0:^16} {1:{flush}{halign}{w}}'.format('', msg[0], w=cols - 34, halign=halign, flush=flush)
            for line in msg[1:]:
                string += '\n{0:^18} {1:{flush}{halign}{w}}'.format('', line, w=cols - 34, halign=halign, flush=flush)
        else:
            string = '[{0:^16}] {1:{flush}{halign}{w}} [{2:^12}]'.format(tag, msg[0],
                                                                         datetime.now().strftime('%H:%M:%S'),
                                                                         w=cols - 34, halign=halign, flush=flush)
            for line in msg[1:]:
                string += '\n{0:^18} {1:{flush}{halign}{w}} {2:^14}'.format('', line, '', w=cols - 34, halign=halign,
                                                                            flush=flush)

    else:
        if tag == None:
            string = '{0:^18} {1:{flush}{halign}{w}}'.format('', msg, w=cols - 34, halign=halign, flush=flush)
        else:
            string = '[{0:^16}] {1:{flush}{halign}{w}} [{2:^12}]'.format(tag, msg, datetime.now().strftime('%H:%M:%S'),
                                                                         w=cols - 34, halign=halign, flush=flush)

    if sameLine == True:
        sys.stdout.write('{} \r'.format(string))
        sys.stdout.flush()
    elif sameLine == False:
        print string

# ====================================================================================================================
def fnProgress( iterator, description = ' ', cols=None, nested = False):
    if cols == None:
        try:
            cols = get_terminal_width()
            if cols < 80:
                raise
        except:
            cols = 100

    return tqdm(iterator, bar_format= '{n_fmt:>5}/{total_fmt:<5}({percentage:3.0f}%)  {bar}    {remaining:>10}   ', ncols = cols)


# ====================================================================================================================

def fnGetYOrbit(filename):
    """	Function to open a settings file (see template) and store them in
        a dictionary
        FileName - can both be absolute or relative
    """
    with open(filename, "r") as File:
        data = [line.strip("\n") for line in File if not line.startswith("#")
                and not line == "\n"]

    ParamsDic = {'Stellar': line.strip().split(':')[1] for line in data if line.startswith('Stellar')}
    ParamsDic.update(
        {re.sub('\s+', '\t', line).strip().split('\t')[0]: re.sub('\s+', '\t', line).strip().split('\t')[2] for line in
         data if line.startswith(' delta')})

    data = [re.sub('\s+', '\t', line).strip().split('\t')[1:]
            for line in data if line.startswith(" p_1")]

    ParamsDic.update({line[0]: line[2] for line in data})
    ParamsDic.update({line[0]: line[2] for line in data if line[1] == '[mJup]'})

    return ParamsDic


# ====================================================================================================================

def fnInitTerminal(main):
    # os.system('cls' if os.name == 'nt' else 'clear')
    # os.system("cls")
    # os.system("clear")
    fnPrintLine(None, '')
    fnPrintLine(None, '')
    fnPrintLine(None, '')
    fnPrintLine(None, '', align='center', flush='=')
    fnPrintLine(None, '{:<} \t(version: {:<})'.format(main.__title__, main.__version__), align='center')
    fnPrintLine(None, '')
    fnPrintLine(None, 'Author:{:<} \temail: {:<}'.format(main.__author__, main.__email__), align='center')
    fnPrintLine(None, '')
    fnPrintLine(None, 'Last Update:{:<}'.format(main.__date__), align='center')
    fnPrintLine(None, '')
    fnPrintLine(None, '', align='center', flush='=')
    fnPrintLine(None, '')
    # endregion
