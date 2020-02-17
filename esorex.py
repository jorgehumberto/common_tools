__author__ = 'jmartins'

import cpl


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
        except:
            param.value = param.default

    return cplRecipe
