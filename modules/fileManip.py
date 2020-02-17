# -*- coding: utf-8 -*-

"""
Created on Fri Jan 30 18:06:45 2015

@author: jmartins


"""
import os

from astropy.io import fits


###############################################################################
# Functions
###############################################################################
def initPaths(folders):
    '''
    initPaths() - function that will create the initial file/folder structure
       
    '''
    for folder in folders:
        # Check for folders existance, if not create them
        if not os.path.exists(folder):
            os.makedirs(folder)
    return


# -----------------------------------------------------------------------------------------------------------------------
def fngetRecipeFileList(recipeName='', dataList='', keyword=''):
    return [filename for filename in dataList \
            if keyword.lower() in fits.getheader(filename)['HIERARCH ESO DPR TYPE'].lower()]


# -----------------------------------------------------------------------------------------------------------------------
def fnUpdateProductList(baseFolder, productList={}):
    fileList = [filenameFull \
                for root, dirs, files in os.walk(baseFolder, topdown=False) \
                for filenameFull in [os.path.join(root, filename) for filename in files] \
                if filenameFull.endswith('.fits')]

    fileListHead = [[filenameFull, fits.getheader(filenameFull)['HIERARCH ESO PRO CATG']] \
                    for filenameFull in fileList]

    categoryList = list(set([filename[1] for filename in fileListHead]))

    productList.update({category: [filename[0] for filename in fileListHead if category in filename[1]] \
                        for category in categoryList})

    return productList
