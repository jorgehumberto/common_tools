__author__ = 'jmartins'
'''

'''


def fnClassificationRules_UVES(fitsHeader):
    '''
    from UVES classification rules (uves.oca) in gasgano
    IN:

    OUT:
        categoryOut =   DO.CATG
        typeOut     =   RAW.TYPE
    '''

    # BIAS FRAMES
    if 'CALIB' in fitsHeader['HIERARCH ESO DPR CATG'] and 'BIAS' in fitsHeader['HIERARCH ESO DPR TYPE']:
        if fitsHeader['HIERARCH ESO DET CHIPS'] == 1:
            fitsHeader.set('D0.CATG', 'BIAS_BLUE')
        elif fitsHeader['HIERARCH ESO DET CHIPS'] == 2:
            fitsHeader.set('D0.CATG', 'BIAS_RED')

    # SCIENCE FRAMES
    if 'SCIENCE' in fitsHeader['HIERARCH ESO DPR CATG'] and 'ECHELLE' in fitsHeader[
        'HIERARCH ESO DPR TECH'] and 'POINT' in fitsHeader[
        'HIERARCH ESO DPR TYPE']:
        if fitsHeader['HIERARCH ESO DET CHIPS'] == 1:
            fitsHeader.set('D0.CATG', 'SCI_POINT_BLUE')
        elif fitsHeader['HIERARCH ESO DET CHIPS'] == 2:
            fitsHeader.set('D0.CATG', 'SCI_POINT_RED')

    # ARC FORM FRAMES
    if 'CALIB' in fitsHeader['HIERARCH ESO DPR CATG'] and 'LAMP' in fitsHeader['HIERARCH ESO DPR TYPE'] and 'FMTCHK' in \
            fitsHeader[
                'HIERARCH ESO DPR TYPE'] and 'ECHELLE' in fitsHeader['HIERARCH ESO DPR TECH']:
        if fitsHeader['HIERARCH ESO DET CHIPS'] == 1:
            fitsHeader.set('D0.CATG', 'ARC_LAMP_FORM_BLUE')
        elif fitsHeader['HIERARCH ESO DET CHIPS'] == 2:
            fitsHeader.set('D0.CATG', 'ARC_LAMP_FORM_RED')

    # ORDER FLAT FRAMES
    if 'CALIB' in fitsHeader['HIERARCH ESO DPR CATG'] and 'LAMP' in fitsHeader[
        'HIERARCH ESO DPR TYPE'] and 'ORDERDEF' in fitsHeader[
        'HIERARCH ESO DPR TYPE'] and 'ECHELLE' in fitsHeader['HIERARCH ESO DPR TECH']:
        if fitsHeader['HIERARCH ESO DET CHIPS'] == 1:
            fitsHeader.set('D0.CATG', 'ORDER_FLAT_BLUE')
        elif fitsHeader['HIERARCH ESO DET CHIPS'] == 2:
            fitsHeader.set('D0.CATG', 'ORDER_FLAT_RED')

    # FLAT FRAMES
    if 'CALIB' in fitsHeader['HIERARCH ESO DPR CATG'] and 'LAMP,FLAT' in fitsHeader[
        'HIERARCH ESO DPR TYPE'] and 'ECHELLE' in fitsHeader[
        'HIERARCH ESO DPR TECH'] and 'FREE' in fitsHeader['HIERARCH ESO INS SLIT1 NAME']:
        if fitsHeader['HIERARCH ESO DET CHIPS'] == 1:
            fitsHeader.set('D0.CATG', 'FLAT_BLUE')
        elif fitsHeader['HIERARCH ESO DET CHIPS'] == 2:
            fitsHeader.set('D0.CATG', 'FLAT_RED')

    # ARC LAMP FRAMES
    if 'CALIB' in fitsHeader['HIERARCH ESO DPR CATG'] and 'LAMP' in fitsHeader['HIERARCH ESO DPR TYPE'] and 'WAVE' in \
            fitsHeader[
                'HIERARCH ESO DPR TYPE'] and 'ECHELLE' in fitsHeader['HIERARCH ESO DPR TECH'] and 'FREE' in fitsHeader[
        'HIERARCH ESO INS SLIT1 NAME']:
        if fitsHeader['HIERARCH ESO DET CHIPS'] == 1:
            fitsHeader.set('D0.CATG', 'ARC_LAMP_BLUE')
        elif fitsHeader['HIERARCH ESO DET CHIPS'] == 2:
            fitsHeader.set('D0.CATG', 'ARC_LAMP_RED')

    return fitsHeader

    # def fnClassificationRules_STDCALIB_UVES(fitsHeader):
    #     # STANDART CALIBRATION FILES
    #     # standardCalibrationKeywords = ['FLUX_STD_TABLE', 'LINE_INTMON_TABLE', 'MASTER_RESPONSE_', 'LINE_REFER_TABLE', 'LINE_INTMON_TABLE', 'CORVEL_MASK', 'EXTCOEFF_TABLE']
    #     # if any(kw in fitsHeader['HIERARCH ESO PRO CATG'] for kw in standardCalibrationKeywords):
    #     fitsHeader.set('D0.CATG',  fitsHeader['HIERARCH ESO PRO CATG'])
    #
    #     return fitsHeader
