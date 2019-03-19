# -*- coding: utf-8 -*-
"""
Created on Tue Set 08 2015

@author: Jorge Martins
email: jorge.martins@iastro.pt


"""
###############################################################################
# Python Modules
###############################################################################
import cpl


###############################################################################
# Functions
###############################################################################

def fnSetRecipe(recipe, customConfig=''):
    '''
    Function that will initialise a cpl recipe instance. If a custom config file has been given, values in this config
    file will be be set on the recipe
    '''
    cplRecipe = cpl.Recipe(recipe)

    # parameters
    for param in cplRecipe.param:
        try:
            param.value = customConfig.get(recipe, param.name)
            # print param.name, param.value, 'new'
        except:

            param.value = param.default
            print param.name, param.value, 'default'
        param.debug = False

    return cplRecipe
