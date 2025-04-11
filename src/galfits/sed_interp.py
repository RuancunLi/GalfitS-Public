import numpy as np
from scipy.interpolate import RegularGridInterpolator
import extinction
from astropy.convolution import Gaussian1DKernel,convolve_fft
import os
from astropy.table import Table, Column
#from .utils import *

#__all__ = ['get_host_SED']

c=2.9979246e5
C_G=6.67259e-20  #kg km3 s-2
C_year=31556926
C_day=86400.
C_ms=1.9891*1e30   # kg
L_edd=1.25e38
tl2=2*np.sqrt(2*np.log(2))      # fwhm = tl2 * sigma


if "MorphSED_DATA_PATH" not in os.environ:
    raise Exception('You should set environment varialbe `MorphSED_DATA_PATH` in your .bashrc (or rc file for other shells)')
else:
    DATA_PATH=os.environ.get('MorphSED_DATA_PATH')

sed_data = np.load('{0}/templates/thin_disk_lowr.npz'.format(DATA_PATH))
points = (sed_data['spin'], sed_data['logMdot'], sed_data['logM'], sed_data['wave'])
intp_thindisk = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)

sed_data = np.load('{0}/templates/BLRDC.npz'.format(DATA_PATH))
points = (sed_data['spin'], sed_data['logM'], sed_data['logMdot'], sed_data['wave'])
intp_BLRDC = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)

sed_data = np.load('{0}/templates/BLRTOT.npz'.format(DATA_PATH))
points = (sed_data['spin'], sed_data['logM'], sed_data['logMdot'], sed_data['wave'])
intp_BLRTOT = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)

sed_data = np.load('{0}/templates/host_conti.npz'.format(DATA_PATH))
points = (sed_data['Z'], sed_data['age'], sed_data['wave'])
intp_host_cont = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)

sed_data = np.load('{0}/templates/host_inst.npz'.format(DATA_PATH))
points = (sed_data['Z'], sed_data['age'], sed_data['wave'])
intp_host_inst = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)

sed_data = np.load('{0}/templates/dustsb99.npz'.format(DATA_PATH))
points = (sed_data['Z'], sed_data['age'], sed_data['wave'])
intp_nebular99 = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)

sed_data = np.load('{0}/templates/dustsb99_conti.npz'.format(DATA_PATH))
points = (sed_data['Z'], sed_data['age'], sed_data['wave'])
intp_nebular99_conti = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)

sed_data = np.load('{0}/templates/host_conti_hires.npz'.format(DATA_PATH))
points = (sed_data['Z'], sed_data['age'], sed_data['wave'])
intp_host_hres_cont = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)

sed_data = np.load('{0}/templates/host_inst_hires.npz'.format(DATA_PATH))
points = (sed_data['Z'], sed_data['age'], sed_data['wave'])
intp_host_hres_inst = RegularGridInterpolator(points, sed_data['sed'],bounds_error=False,fill_value=0.)
wave_temp = np.loadtxt('{0}/templates/wave_hires.txt'.format(DATA_PATH))

def gaussian(x, amp, cen, wid):
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))

def gaussian3D(x, amp, cen, wid):
    ny,nx = amp.shape
    tlen = len(x)
    IFU = np.zeros((ny,nx,tlen))
    for loop in range(tlen):
        IFU[:,:,loop] = (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x[loop]-cen)**2 / (2*wid**2))
    return IFU


feii_template_op = Table.read('{0}/templates/irontemplate_op_new.ipac'.format(DATA_PATH),format='ascii')
feii_template_uv = Table.read('{0}/templates/irontemplate_uv.ipac'.format(DATA_PATH),format='ascii')
wave_bac = np.genfromtxt('{0}/templates/BLRDC.wavegrid'.format(DATA_PATH), dtype=float, names=None)

def FeII(x,A_uv,A_op,dcen,fwhm):
    wavemodel = np.append(feii_template_uv['Spectral_axis'],feii_template_op['wave'])
    fluxmodel = np.append(A_uv*feii_template_uv['Intensity'],A_op*feii_template_op['flux'])
    if fwhm < 900.:
        fwhm=910.
    fwhm = np.sqrt(fwhm**2-900**2)
    std = fwhm/c/tl2
    dcenlog = dcen/c
    referwave = 2500.
    kersize = np.max([int(std*referwave*5),15])
    kernel = Gaussian1DKernel(referwave*std,x_size=kersize)
    kernel.normalize()
    mask = (wavemodel > (1.-5*std-dcenlog)*x.min())&(wavemodel < (1.+5*std+dcenlog)*x.max())
    if np.sum(mask) == 0:
        return np.zeros_like(x)
    wave = wavemodel[mask]
    flux = fluxmodel[mask]
    return conv_spec(x,wave,flux,kernel.array,dcenlog,referwave)

def BaC(x,cf,logM,logMdot,spin,dcen,fwhm):
    std = fwhm/c/tl2
    dcenlog = dcen/c
    referwave = 2500.
    kersize = np.max([int(std*referwave*5),15])
    kernel = Gaussian1DKernel(referwave*std,x_size=kersize)
    kernel.normalize()
    mask = (wave_bac > (1.-5*std-dcenlog)*x.min())&(wave_bac < (1.+5*std+dcenlog)*x.max())
    if np.sum(mask) == 0:
        return np.zeros_like(x)
    wave = wave_bac[mask]
    flux = 10**intp_BLRDC((spin,logM,logMdot,wave))
    return cf*conv_spec(x,wave,flux,kernel.array,dcenlog,referwave)

def conv_spec(x, wave, flux, kernel, dcenlog, referwave):
    '''
    convolve the give spectra ['wave', 'flux'] with any given 'kernel'
    return the flux respect to x
    '''
    osfactor=2
    logw = np.log(wave)
    wavegrid = np.linspace(np.min(logw),np.max(logw),osfactor*len(wave))
    gridsparse = np.abs(wavegrid[1]-wavegrid[0])
    newkernelsize = int((len(kernel)-1)/(gridsparse*referwave))
    newkernel = np.interp(gridsparse*referwave*np.arange(newkernelsize),np.arange(len(kernel)),kernel)
    #newkernel*=10*np.sqrt(gridsparse*referwave)
    f_new = np.interp(np.exp(wavegrid), wave, flux)
    f_conv = convolve_fft(f_new,kernel=newkernel)
    return np.interp(x, np.exp(wavegrid+dcenlog), f_conv)

def get_host_SED_3D(x, logM, f_cont, age, Z, sigma, Av, C_unit):
    # intp_burst_cont normalized by SFR = 10 Ms/yr
    # intp_burst_inst normalized by M_ast = 10^7 Ms
    # x: wavelength [91, 1.6e6] Angs
    # M: overall stellar mass;
    # f_cont: float, [0, 1.], fraction of the continuum starformation component
    # age: the age of the starburst component [1e-5, 13.69] Gyr
    # Z : metalicity, [0.001,0.04], 0.02 for solar
    # sigma  : velocity dispersion km/s
    # Av : V Attenuation
    # Return: restframe L_lambda in ergs/s
    M=10**logM
    age_yr = age*1e9
    M_cont = M*f_cont
    SFR = M_cont*1e-10
    M_inst = M-M_cont
    sed1 = (10**intp_host_hres_inst((Z, age_yr, wave_temp)))*M_inst*1e-7
    if f_cont == 0. :
        sed2 = np.zeros_like(sed1,dtype=float)
    else:
        sed2 = (10**intp_host_hres_cont((Z, 1e10, wave_temp)))*SFR*1e-1
    wavemodel = wave_temp.copy()
    fluxmodel = C_unit*(sed1+sed2)
    std = sigma/c
    referwave = 5000.
    kersize = np.max([int(std*referwave*5),15])
    kernel = Gaussian1DKernel(referwave*std,x_size=kersize)
    kernel.normalize()
    mask = (wavemodel > (1.-5*std)*x.min())&(wavemodel < (1.+5*std)*x.max())
    wave = wavemodel[mask]
    flux = fluxmodel[mask]
    if mask[0] == 1:
        startx = np.argmax(x > wavemodel[0]*(1.+5*std) )
    else:
        startx = 0
    if mask[-1] == 1:
        stopx = np.argmax(x > wavemodel[-1]*(1.-5*std))-1
    else:
        stopx = -1
    if (stopx == -1)&(startx==0):
        rest_spec = conv_spec(x,wave,flux,kernel.array,0.,referwave)
    else:
        mask = np.zeros_like(x,dtype=bool)
        rest_spec = np.zeros_like(x,dtype=float)
        mask[startx:stopx] = 1
        rest_spec[mask] = conv_spec(x[mask],wave,flux,kernel.array,0.,referwave)
        sed1 = (10**intp_host_inst((Z, age_yr, x[~mask])))*M_inst*1e-7
        sed2 = (10**intp_host_cont((Z, 1e10, x[~mask])))*SFR*1e-1
        rest_spec[~mask]=C_unit*(sed1+sed2)
    cm=extinction.ccm89(x,Av,3.1)/2.5
    return rest_spec/(10**cm)

def get_nebular99(x, logM, f_cont, age, Z, C_unit):
    if f_cont < 0.:
        f_cont = 0.
    elif f_cont > 1. :
        f_cont = 1.
    M=10**logM
    age_yr = age*1e9
    M_cont = M*f_cont
    SFR = M_cont*1e-10
    M_inst = M-M_cont
    sed1 = (10**intp_nebular99((Z, age_yr, x)))*M_inst*1e-7
    if f_cont == 0. :
        sed2 = np.zeros_like(sed1,dtype=float)
    else:
        sed2 = (10**intp_nebular99_conti((Z, 1e10, x)))*SFR*1e-1
    fluxmodel = C_unit*(sed1+sed2)
    return fluxmodel

def get_host_SED(x, logM, f_cont, age, Z, Av, C_unit):
    # intp_burst_cont normalized by SFR = 10 Ms/yr
    # intp_burst_inst normalized by M_ast = 10^7 Ms
    # x: wavelength [91, 1.6e6] Angs
    # M: overall stellar mass;
    # f_cont: float, [0, 1.], fraction of the continuum starformation component
    # age: the age of the starburst component [1e-5, 13.69] Gyr
    # Z : metalicity, [0.001,0.04], 0.02 for solar
    # Av : V Attenuation
    # Return: restframe L_lambda in ergs/s
    if f_cont < 0.:
        f_cont = 0.
    elif f_cont > 1. :
        f_cont = 1.
    M=10**logM
    age_yr = age*1e9
    M_cont = M*f_cont
    SFR = M_cont*1e-10
    M_inst = M-M_cont
    sed1 = (10**intp_host_inst((Z, age_yr, x)))*M_inst*1e-7
    if f_cont == 0. :
        sed2 = np.zeros_like(sed1,dtype=float)
    else:
        sed2 = (10**intp_host_cont((Z, 1e10, x)))*SFR*1e-1
    fluxmodel = C_unit*(sed1+sed2)
    cm=extinction.ccm89(x,Av,3.1)/2.5
    return fluxmodel/(10**cm)

def get_AGN_SED(x, logM,logMdot,spin,C_unit):
    # derive a SED of standard thin disk
    # x:       wavelength, in range[125,24700]
    # logM:    log10(M_BH), in range[5,10]
    # logMdot: log10(dotM), in range[-4,2]
    # spin:    BH spin, in range[0,0.99]
    # Return: restframe L_lambda in ergs/s
    fluxmodel=100.*C_unit*intp_thindisk((spin, logMdot, logM, x))
    return fluxmodel

def sed_to_obse(x,y,z,ebv):
    xnew = x.copy()
    ynew = y.copy()
    av=3.1*ebv
    xnew*=(1+z)
    ynew/=(1+z)**3
    cm=extinction.ccm89(xnew,av,3.1)
    ynew/=np.power(10,cm/(2.5))
    return xnew,ynew

def sed_to_rest(x,y,z,ebv):
    xnew = x.copy()
    ynew = y.copy()
    av=3.1*ebv
    cm=extinction.ccm89(xnew,av,3.1)
    ynew*=np.power(10,cm/(2.5))
    xnew/=(1+z)
    ynew*=(1+z)**3
    return xnew,ynew

def spowerlaw(x,A,C):
    y = A*np.power(x/5100.,-C)/5100
    return y

def brokenpower(x,A,C1,C2):
    y = A*np.power(x/5100.,-C1)
    y[x<5100.] = A*np.power(x[x<5100.]/5100.,-C2)
    return y

def Bbody(x, Teff):
    # Normalized to R_s
    fx = 190095060845894.53/(x**5*(np.exp(143876866./(Teff*x))-1))
    return fx
