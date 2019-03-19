import numpy as np
from math_local.mathFunctions import fnMovingAverage



#----------------------------------------------------------------------------------------------------------------------
def fnInferInstrumentResponse(data, window = 100):
    '''
    :param data: array containing the reduced flat field
    :return: instrument response defined as a moving average of the flat field array using
    a moving average window of 100 pixels
    '''
    instrumentResponse = np.empty_like(data);
    nOrders = len(data)
    for order in np.arange(nOrders):
        orderData = data[order]
        instrumentResponse[order] = np.zeros_like(orderData)
        dataRange = np.arange(len(orderData))[np.nonzero(orderData)][:-10]
        instrumentResponse[order][dataRange]= fnMovingAverage(orderData[dataRange], window)



    return instrumentResponse

