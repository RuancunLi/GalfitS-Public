import pathlib
import numpy as np
from astropy import units as u

__all__ = ['package_path', 'get_wcs_rotation']

package_path = pathlib.Path(__file__).parent.absolute()

def get_wcs_rotation(wcs):
    '''
    Get east of north rotation of the WCS.

    Parameters
    ----------
    wcs : `astropy.WCS`

    Returns
    -------
    rot : angle (radian)
        Rotation angle, east of north.
    '''
    cdmat = wcs.pixel_scale_matrix
    sgn = np.sign(np.linalg.det(cdmat))
    rot = np.arctan2(sgn * cdmat[1, 0], sgn * cdmat[0, 0]) * u.radian
    return rot
