import numpy as np
from astropy.modeling.functional_models import Sersic1D
from .galaxy import coordinates_transfer,indentify_xy
import cv2

def alpha_tanh(r_3d, r_in, r_out, alpha, theta_out):
    theta_out_rad = theta_out*np.pi/180.
    return theta_out_rad*hypotan(r_3d, r_in, r_out, alpha, theta_out_rad)*np.power(0.5*(r_3d/r_out+1.),alpha)

def hypotan(r_3d, r_in, r_out, alpha, theta_out_rad):
    CDEF = 0.23
    A = 2.*CDEF/(np.abs(theta_out_rad)+CDEF)-1.00001
    B = (2.-np.arctanh(A))*r_out/(r_out-r_in)
    return 0.5*np.tanh(B*(r_3d/r_out-1.)+2.)+0.5

def R_2d(x,y):
    return (np.sqrt(x**2+y**2))

def R_3d(x,y,PA,i):
    phi = PA*np.pi/180.
    dis_proj_square = np.abs(np.cos(phi)*y-np.sin(phi)*x)**2
    return np.sqrt(R_2d(x,y)**2+(1./np.cos(i*np.pi/180.)-1)*dis_proj_square)

def xy_3d(x,y,PA,i):
    phi = PA*np.pi/180.
    dis_proj = np.cos(phi)*y-np.sin(phi)*x
    projvec = (1./np.cos(i*np.pi/180.)-1.)*dis_proj
    return x-projvec*np.sin(phi),y+projvec*np.cos(phi)

def R_fourier(x,y,m,am,theta_m,PA,i):
    lower = y < 0.
    phi_map = np.arccos(x/R_2d(x,y))
    phi_map[lower] = 2.*np.pi-phi_map[lower]
    return R_3d(x,y,PA,i)*(1.+am*np.cos(m*(phi_map+theta_m*np.pi/180.)))

def make_fourier(M, seRe, sen, shape, posi, PA, i, r_in, r_out, alpha, theta_out, m, am, theta_m, i_arm, transpar=None):
    '''
    seRe,sen: control the mass distribution of the arms
    i, PA galaxy projection
    r_in, r_out, alpha, theta_out: rotation curve
    m, am, theta_m, i_arm: arm number, relative strength
    '''
    ny = shape[0]
    nx = shape[1]
    xaxis = np.arange(nx)
    yaxis = np.arange(ny)
    xcen = posi[0]
    ycen = posi[1]
    if transpar is not None:
        xcen,ycen = coordinates_transfer(posi[0],posi[1],transpar)
        seRe = seRe /transpar['pixsc']
        r_in = r_in /transpar['pixsc']
        r_out = r_out /transpar['pixsc']
        testPA = PA+transpar['delta_ang']
        if testPA > 90.:
            testPA -= 180.
        elif testPA < -90.:
            testPA += 180.
        PA = testPA
    else:
        xcen,ycen = indentify_xy(posi[0],posi[1])
    xmesh, ymesh = np.meshgrid(xaxis, yaxis)
    x_p = xmesh + 0.5 - xcen
    y_p = ymesh + 0.5 - ycen
    r_2d = R_2d(x_p,y_p)
    r_3d = R_3d(x_p,y_p,PA,i)
    theta_r = alpha_tanh(r_3d, r_in, r_out, alpha, theta_out)
    x_p,y_p = xy_3d(x_p,y_p,PA,i)
    x_i = x_p*np.cos(theta_r) + y_p*np.sin(theta_r)
    y_i = -x_p*np.sin(theta_r) + y_p*np.cos(theta_r)
    #i_arm = 0.
    r_3d_i = R_fourier(x_i,y_i,m,am,theta_m,PA,i_arm)
    s1 = Sersic1D(amplitude=1, r_eff=seRe, n=sen)
    mod = s1(r_3d_i)
    mass_map = M*mod/np.sum(mod)
    return mass_map

def rotate_CV(image, angel , interpolation):
    '''
        input :
        image           :  image                    : ndarray
        angel           :  rotation angel           : int
        interpolation   :  interpolation mode       : cv2 Interpolation object
        returns :
        rotated image   : ndarray

        '''
    #in OpenCV we need to form the tranformation matrix and apply affine calculations
    #
    h,w = image.shape[:2]
    cX,cY = (w//2,h//2)
    M = cv2.getRotationMatrix2D((cX,cY),angel,1)
    rotated = cv2.warpAffine(image,M , (w,h),flags=interpolation)
    return rotated

def reshape_img(posi, image, shapeF, coordinates_transfer_par0, coordinates_transfer_parF):
    '''
        posi: position of center in [dra, ddec]
    '''
    dtheta = coordinates_transfer_parF['delta_ang']-coordinates_transfer_par0['delta_ang']
    xpix,ypix = coordinates_transfer(posi[0],posi[1],coordinates_transfer_parF)
    xpix0,ypix0 = coordinates_transfer(posi[0],posi[1],coordinates_transfer_par0)
    ny,nx = image.shape
    nyF = shapeF[0]
    nxF = shapeF[1]
    xs = xpix*nx/nxF
    ys = ypix*ny/nyF
    T = np.float32([[1, 0, xs-xpix0], [0, 1, ys-ypix0]])
    img_translation = cv2.warpAffine(image, T, (nx, ny),flags = cv2.INTER_CUBIC)
    M = cv2.getRotationMatrix2D((xs,ys),dtheta,1)
    img_rot = cv2.warpAffine(img_translation, M, (nx, ny),flags = cv2.INTER_CUBIC)
    resized = cv2.resize(img_rot, (nxF, nyF), interpolation = cv2.INTER_CUBIC)
    return resized

def non_para(par, c, plist, image):
    npimg = np.zeros_like(image, dtype=float)
    for loop,p in enumerate(plist):
        npimg[p[1],p[0]] = par['pf_{0}_{1}'.format(c,loop)]
    return npimg
