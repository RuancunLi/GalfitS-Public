import numpy as np

__all__ = ['get_zp']

def get_zp(filter_name):
    '''
    Get the magnitude zeropoint.

    Parameters
    ----------
    filter_name : string
        Name of the filter.

    Returns
    -------
    zp : float or None
        If found the magnitude zeropoint is returned, otherwise, None is returned.
    '''
    instrument, bandpass = filter_name.split('.')
    if instrument == 'panstarrs':
        bp_dict = filter_panstarrs.get(bandpass, None)
        if bandpass:
            zp = bp_dict['ZP']
        else:
            zp = None
    else:
        zp = None

    return zp

filter_panstarrs = {
    'g_P1' : dict(ZP=24.56),
    'r_P1' : dict(ZP=24.76),
    'i_P1' : dict(ZP=24.74),
    'z_P1' : dict(ZP=24.33),
    'y_P1' : dict(ZP=23.33),
    'w_P1' : dict(ZP=26.04),
    'open' : dict(ZP=26.37),
    'ref' : 'https://ui.adsabs.harvard.edu/abs/2012ApJ...750...99T/abstract'
}
