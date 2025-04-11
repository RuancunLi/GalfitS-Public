import pyprofit
import lmfit
from  lmfit.models import GaussianModel
from scipy.integrate import trapz
from astropy.convolution import convolve_fft
import numpy as np
import extinction
from astropy.table import Table
from scipy.interpolate import interp1d
import sys
from . import sed_interp as SEDs
from . import emission_lines as EL
from pathlib import Path

'''
"allbands":
'acs_f625w', '4star_m_j2', 'wfc3_f555w', 'wfc3_f139m', 'acs_f475w', 'ukidss_h', 'ndwfs_r', '4star_m_j3', 'acs_f435w',
    'ch1', 'ndwfs_i', '4star_j', 'galex_nuv', 'wfc3_f606w', '4star_m_hlong', 'wfc3_f125w', 'newfirm_j', 'newfirm_ks',
    'sloan_z', 'wfc3_f140w', 'sloan_g', 'sloan_i', 'wfc3_f814w', 'sloan_u', 'wfc3_f775w', 'sloan_r', 'r', 'wise_ch1',
    '4star_m_hshort', 'i', '4star_ks', 'wfpc2_f450w', 'README', 'h', 'wfc3_f275w', '4star_h', 'ch3', 'ukidss_k',
    'wfc3_f218w', 'ch4', 'ukidss_y', '4star_m_j1', 'wfc3_f110w', 'ukidss_j', 'ch2', 'wfc3_f153m', 'acs_f606w', 'ndwfs_bw',
    'wfpc2_f814w', 'galex_fuv', 'wfc3_f225w', 'wfpc2_f675w', 'ks', 'acs_f555w', 'wfc3_f625w', 'wfc3_f127m', 'wfc3_f475w',
    'wfpc2_f555w', 'wfc3_f438w', 'wfc3_f105w', 'newfirm_h', 'wfc3_f160w', 'j', 'v', 'acs_f775w', 'wfpc2_f606w', 'wise_ch2',
    'acs_f814w', 'wfc3_f850lp', 'b', 'wise_ch3', 'wise_ch4', 'acs_f850lp'
'''

filterpath = Path(__file__).parent / 'data/filters'

tl2=2*np.sqrt(2*np.log(2))
c=2.9979246e5

def indentify_xy(x,y):
    '''
    Transform the (x,y) from (\Delat RA, \Delta Dec) space to (pix, pix)
    '''
    return x,y

def coordinates_transfer(x, y, kwargs):
    '''
    Transform the (x,y) from (\Delat RA, \Delta Dec) space to (pix, pix)
    '''
    xp = kwargs['x0'] + kwargs['x0shift'] + x*kwargs['dxra'] + y*kwargs['dxdec']
    yp = kwargs['y0'] + kwargs['y0shift'] + x*kwargs['dyra'] + y*kwargs['dydec']
    return xp,yp

def coordinates_transfer_inverse(xp, yp, kwargs):
    '''
    Transform the (xp,yp) from (pix, pix) space to (\Delat RA, \Delta Dec)
    '''
    x = ((xp-kwargs['x0'])*kwargs['dydec'] - (yp-kwargs['y0'])*kwargs['dxdec'])/(kwargs['dxra']*kwargs['dydec']-kwargs['dyra']*kwargs['dxdec'])
    y = ((xp-kwargs['x0'])*kwargs['dyra'] - (yp-kwargs['y0'])*kwargs['dxra'])/(kwargs['dyra']*kwargs['dxdec']-kwargs['dxra']*kwargs['dydec'])
    return x,y

class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration

    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False

def IFU_to_img(IFU,wave,band,step=0.5):
    '''
    Transform a IFU data cube in to a 2D image

    IFU:   the input 3D array with z as the wavelength dimension
    wave:  1D array shows the sampled wavelength
    band:  choose one from "all bands"
    step:  float, wavelength accuracy to integrate flux
    '''
    resp = Table.read(filterpath / band,format='ascii')
    filter_x=resp['col1']
    filter_y=resp['col2']
    tminx = np.max([np.min(filter_x),np.min(wave)])
    tmaxx = np.min([np.max(filter_x),np.max(wave)])
    interX = np.arange(tminx,tmaxx,step)
    f2=interp1d(filter_x,filter_y,bounds_error=False,fill_value=0.)
    ax=trapz(f2(interX),x=interX)
    nz,ny,nx = IFU.shape
    image = np.zeros((ny,nx))
    for loopy in range(ny):
        for loopx in range(nx):
            f1=interp1d(wave,IFU[:,loopy,loopx],bounds_error=False,fill_value=0.)
            tof=lambda x : f1(x)*f2(x)
            image[loopy][loopx] = trapz(tof(interX),x=interX)
    image /= ax
    return image

def Cal_map(r, type, paradic):
    for case in switch(type):
        if case('linear'):
            return paradic['b'] + paradic['k']*r
            break
        if case('exp'):
            return (paradic['in']-paradic['out'])*np.exp(-r/paradic['k'])+paradic['out']
            break
        if case('const'):
            return np.ones_like(r,dtype=float)*paradic['value']
            break
        if case():
            raise ValueError("Unidentified method for calculate age or Z map")

def Cal_gradient_map(r_grid, r_map, method, type, paradic, flux_band, f2, interX, ax, **kwargs):
    if type =='const':
        return np.ones_like(r_map,dtype=float)
    else:
        fratio_list = []
        targe_para = Cal_map(r_grid,type,paradic)
        for loopr in range(len(r_grid)):
            if method =='age':
                centerSED = SEDs.get_host_SED(interX, 0., kwargs['f_cont_zero'], targe_para[loopr], kwargs['Z_zero'],kwargs['Av_zero'],1.)
            elif method =='f_cont':
                centerSED = SEDs.get_host_SED(interX, 0., targe_para[loopr], kwargs['age_zero'], kwargs['Z_zero'],kwargs['Av_zero'],1.)
            elif method =='Z':
                centerSED = SEDs.get_host_SED(interX, 0., kwargs['f_cont_zero'], kwargs['age_zero'], targe_para[loopr],kwargs['Av_zero'],1.)
            elif method =='Av':
                centerSED = SEDs.get_host_SED(interX, 0., kwargs['f_cont_zero'], kwargs['age_zero'], kwargs['Z_zero'],targe_para[loopr],1.)
            else:
                raise ValueError("Unidentified method for calculate gradient map")
            flux = trapz(centerSED*f2(interX),x=interX)/ax
            fratio_list.append(flux/flux_band)
        fratio_list = np.array(fratio_list)
        targe_map = Cal_map(r_map,type,paradic)
        gradient_map = np.interp(targe_map,np.array(targe_para),fratio_list)
        return gradient_map

class Galaxy(object):
    '''
    the galaxy object
    with physical subcomponents and parameters
    '''
    def __init__(self, mass=1e9, z=0., ebv_G=0.):
        '''
        galaxy object is initialed from a given mass
        '''
        self.mass = mass
        self.Nsub=0
        self.subCs = {}
        self.EL_region = None
        self.ageparams={}
        self.Zparams={}
        self.f_cont={}
        self.Avparams={}
        self.maglist = []
        self.shape = None
        self.mass_map = {}
        self.r_map = None
        self.redshift = z
        self.ebv_G = ebv_G
    def reset_mass(self,mass):
        '''
        reset the mass of a galaxy object
        '''
        self.mass = mass
    def add_subC(self,Pro_names,params,ageparam,Zparam,f_cont,Avparam):
        '''
        To add a subcomponent for a galaxy object

        Pro_names: the name of the profiles
        e.g. = "sersic" "coresersic" "brokenexp" "moffat" "ferrer" "king" "pointsource"

        params: a dictionary of the parameters for this subcomponent
        e.g. for sersic: {'xcen': 50.0, 'ycen': 50.0, 'frac': 0.704, 're': 10.0,
        'nser': 3.0, 'ang': -32.70422048691768, 'axrat': 1.0, 'box': 0.0, 'convolve': False}

        ageparam: a dictionary of the age dsitribution parameters for this subcomponent
        e.g. {'type': 'linear', 'paradic': {'k': -0.05, 'b': 9.0}}
             {'type': 'const', 'paradic': {'value': 5.0}}

        Zparam: a dictionary of the matallicity dsitribution parameters for this subcomponent
        e.g. {'type': 'linear', 'paradic': {'k': 0.0, 'b': 0.02}}

        f_cont: a dictionary of the fraction of starformation dsitribution parameters for this subcomponent
        between [0.,1.], fraction of starformation
        '''
        params['mag']=params['frac']#.pop("frac")
        params['mag'] = 10. - 2.5*np.log10(params['mag'])
        if Pro_names in self.subCs.keys():
            self.subCs[Pro_names].append(params)
            self.ageparams[Pro_names].append(ageparam)
            self.Zparams[Pro_names].append(Zparam)
            self.Avparams[Pro_names].append(Avparam)
            self.f_cont[Pro_names].append(f_cont)
        else:
            self.subCs.update({Pro_names : [params]})
            self.ageparams.update({Pro_names : [ageparam]})
            self.Zparams.update({Pro_names : [Zparam]})
            self.Avparams.update({Pro_names : [Avparam]})
            self.f_cont.update({Pro_names : [f_cont]})
            self.mass_map.update({Pro_names : []})
        #print (self.mass_map)
        self.maglist.append(params['mag'])
        #print (self.maglist)

    def add_arbsubC(self,Pro_names,mass_map,ageparam,Zparam,f_cont,Avparam,sigmaparam={'type': "const", 'paradic':{'value': 100}},arbi_para={'func':'Users'}):
        '''
        To add a subcomponent for a galaxy object

        Pro_names: the name of the profiles
        e.g. = "sersic" "coresersic" "brokenexp" "moffat" "ferrer" "king" "pointsource"

        params: a dictionary of the parameters for this subcomponent
        e.g. for sersic: {'xcen': 50.0, 'ycen': 50.0, 'frac': 0.704, 're': 10.0,
        'nser': 3.0, 'ang': -32.70422048691768, 'axrat': 1.0, 'box': 0.0, 'convolve': False}

        ageparam: a dictionary of the age dsitribution parameters for this subcomponent
        e.g. {'type': 'linear', 'paradic': {'k': -0.05, 'b': 9.0}}
             {'type': 'const', 'paradic': {'value': 5.0}}

        Zparam: a dictionary of the matallicity dsitribution parameters for this subcomponent
        e.g. {'type': 'linear', 'paradic': {'k': 0.0, 'b': 0.02}}

        f_cont: a dictionary of the fraction of starformation dsitribution parameters for this subcomponent
        between [0.,1.], fraction of starformation
        '''
        if Pro_names in self.subCs.keys():
            self.subCs[Pro_names].append(arbi_para)
            self.ageparams[Pro_names].append(ageparam)
            self.Zparams[Pro_names].append(Zparam)
            self.Avparams[Pro_names].append(Avparam)
            self.f_cont[Pro_names].append(f_cont)
        else:
            self.subCs.update({Pro_names : [arbi_para]})
            self.ageparams.update({Pro_names : [ageparam]})
            self.Zparams.update({Pro_names : [Zparam]})
            self.Avparams.update({Pro_names : [Avparam]})
            self.f_cont.update({Pro_names : [f_cont]})
            self.mass_map.update({Pro_names : []})
        self.mass_map[Pro_names].append(mass_map)
        self.shape = mass_map.shape

    def generate_mass_map(self,shape,convolve_func,transpar=None,aperturemask=None):
        '''
        gemerate the mass distribution map for a galaxy object
        shape: return 2D image shape
        convolve_func: a 2D kernel if convolution is needed
        -----
        Caution, in 3D, the r_map calculation can be used as elliptical_func
        pixelscale: in unit pixel/"
        -----
        '''
        mags = np.array(self.maglist)
        magzero = 2.5*np.log10(self.mass/np.sum(np.power(10,mags/(-2.5))))
        ny,nx=shape
        self.shape = shape
        image = np.zeros(shape,dtype=float)
        xaxis = np.arange(nx)
        yaxis = np.arange(ny)
        xmesh, ymesh = np.meshgrid(xaxis, yaxis)
        apertures = []
        self.r_map = None
        for key in self.subCs:
            self.mass_map[key]=[]
            for loop in range(len(self.subCs[key])):
                params = self.subCs[key][loop]
                par_copy = params.copy()
                if transpar is not None:
                    xpix,ypix = coordinates_transfer(params['xcen'],params['ycen'],transpar)
                    par_copy['re'] = par_copy['re'] /transpar['pixsc']
                    testang = par_copy['ang']+transpar['delta_ang']
                    if testang > 90.:
                        testang -= 180.
                    elif testang < -90.:
                        testang += 180.
                    par_copy['ang'] = testang
                else:
                    xpix,ypix = indentify_xy(params['xcen'],params['ycen'])
                par_copy['xcen'] = xpix
                par_copy['ycen'] = ypix
                profit_model = {'width':  nx,
                    'height': ny,
                    'magzero': magzero,
                    'psf': convolve_func,
                    'profiles': {key:[par_copy]}
                   }
                mass_map, _ = pyprofit.make_model(profit_model)
                mass_map = np.array(mass_map)
                mass_map = np.array(mass_map.tolist())
                self.mass_map[key].append(mass_map)
                image += mass_map
                if self.r_map is None:
                    r = np.sqrt( (xmesh+0.5 - xpix)**2. + (ymesh+0.5 - ypix)**2.)
                    if transpar is not None:
                        self.r_map = r*transpar['pixsc']
                    else:
                        self.r_map = r
                if aperturemask is not None:
                    apertures.append(np.sum(mass_map[aperturemask])/ np.sum(mass_map))
        if self.EL_region is not None:
            for loop,eldict in enumerate(self.EL_region):
                params = eldict['params']
                par_copy = params.copy()
                if transpar is not None:
                    xpix,ypix = coordinates_transfer(params['xcen'],params['ycen'],transpar)
                    par_copy['re'] = par_copy['re'] /transpar['pixsc']
                    testang = par_copy['ang']+transpar['delta_ang']
                    if testang > 90.:
                        testang -= 180.
                    elif testang < -90.:
                        testang += 180.
                    par_copy['ang'] = testang
                else:
                    xpix,ypix = indentify_xy(params['xcen'],params['ycen'])
                par_copy['xcen'] = xpix
                par_copy['ycen'] = ypix
                profit_model = {'width':  nx,
                    'height': ny,
                    'magzero': magzero,
                    'psf': convolve_func,
                    'profiles': {eldict['pname']:[par_copy]}
                   }
                mass_map, _ = pyprofit.make_model(profit_model)
                mass_map = np.array(mass_map)
                mass_map = np.array(mass_map.tolist())
                mass_map /= np.sum(mass_map)
                self.EL_region[loop]['normmap'] = mass_map
                if aperturemask is not None:
                    self.EL_region[loop]['speape'] = np.sum(mass_map[aperturemask])/ np.sum(mass_map)
        return image, apertures

    def fiducial_sed(self,wavelength,apertures=None,dustnorm=None,par_tot=None):
        fl = np.zeros_like(wavelength)
        waveintrin = wavelength/(1.+self.redshift)
        fllist = []
        count = 0
        for key in self.subCs:
            for loop in range(len(self.subCs[key])):
                age = Cal_map(0.,self.ageparams[key][loop]['type'],self.ageparams[key][loop]['paradic'])
                Z = Cal_map(0.,self.Zparams[key][loop]['type'],self.Zparams[key][loop]['paradic'])
                f_cont = Cal_map(0.,self.f_cont[key][loop]['type'],self.f_cont[key][loop]['paradic'])
                Av = Cal_map(0.,self.Avparams[key][loop]['type'],self.Avparams[key][loop]['paradic'])
                fc = SEDs.get_host_SED(waveintrin, np.log10(self.mass*self.subCs[key][loop]['frac']/100.), f_cont, age, Z, Av, 1.)
                if dustnorm is not None:
                    fc += SEDs.get_nebular99(waveintrin, np.log10(self.mass*self.subCs[key][loop]['frac']/100.), f_cont, age, Z, dustnorm)
                x,fll =  SEDs.sed_to_obse(waveintrin,fc,self.redshift,self.ebv_G)
                if apertures is not None:
                    fll *= apertures[count]
                fl += fll
                fllist.append(fll)
                count += 1
        if self.EL_region is not None:
            for loop,eldict in enumerate(self.EL_region):
                fc = np.zeros_like(waveintrin)
                prefix = eldict['prefix']
                for line in eldict['Lines']:
                    cwave = par_tot[prefix+line['name']+'_center'].value
                    width = par_tot[prefix+line['name']+'_FWHM'].value*cwave/tl2/c
                    fc += SEDs.gaussian(waveintrin,par_tot[prefix+line['name']+'_Lumin'].value,cwave,width)
                x,fll =  SEDs.sed_to_obse(waveintrin,fc,self.redshift,self.ebv_G)
                if eldict['speape'] is not None:
                    fll *= eldict['speape']
                fl += fll
                fllist.append(fll)
        return fl,fllist

    def generate_image(self,band,convolve_func,inte_step=10,dustnorm=None,par_tot=None):
        resp = Table.read(filterpath / band,format='ascii')
        ny = self.shape[0]
        nx = self.shape[1]
        filter_x=resp['col1']
        filter_y=resp['col2']
        tminx = np.min(filter_x)
        tmaxx = np.max(filter_x)
        interX = np.linspace(tminx,tmaxx,np.max([100,len(filter_x)]))
        f2=interp1d(filter_x,filter_y,bounds_error=False,fill_value=0.)
        ax=trapz(f2(interX),x=interX)
        r_grid = np.logspace(-1,np.log10(np.max(self.r_map)),inte_step)
        totalflux = np.zeros(self.shape,dtype=float)
        interX_intrin = interX/(1.+self.redshift)
        for key in self.subCs:
            for loop in range(len(self.subCs[key])):
                age_zero = Cal_map(0.,self.ageparams[key][loop]['type'],self.ageparams[key][loop]['paradic'])
                Z_zero = Cal_map(0.,self.Zparams[key][loop]['type'],self.Zparams[key][loop]['paradic'])
                f_cont_zero = Cal_map(0.,self.f_cont[key][loop]['type'],self.f_cont[key][loop]['paradic'])
                Av_zero = Cal_map(0.,self.Avparams[key][loop]['type'],self.Avparams[key][loop]['paradic'])
                centerSED = SEDs.get_host_SED(interX_intrin, 0., f_cont_zero, age_zero, Z_zero, Av_zero, 1.)
                if dustnorm is not None:
                    centerSED += SEDs.get_nebular99(interX_intrin, 0., f_cont_zero, age_zero, Z_zero, dustnorm)
                flux_intrin = trapz(centerSED*f2(interX_intrin),x=interX_intrin)/ax
                x, sed_obs = SEDs.sed_to_obse(interX_intrin,centerSED,self.redshift,self.ebv_G)
                flux_band = trapz(sed_obs*f2(interX),x=interX)/ax
                age_gradient = Cal_gradient_map(r_grid, self.r_map, 'age', self.ageparams[key][loop]['type'],self.ageparams[key][loop]['paradic']
                                , flux_intrin, f2, interX_intrin, ax, age_zero=age_zero, Z_zero=Z_zero, f_cont_zero=f_cont_zero, Av_zero=Av_zero)
                Z_gradient = Cal_gradient_map(r_grid, self.r_map, 'Z', self.Zparams[key][loop]['type'],self.Zparams[key][loop]['paradic']
                                , flux_intrin, f2, interX_intrin, ax, age_zero=age_zero, Z_zero=Z_zero, f_cont_zero=f_cont_zero, Av_zero=Av_zero)
                f_cont_gradient = Cal_gradient_map(r_grid, self.r_map, 'f_cont', self.f_cont[key][loop]['type'],self.f_cont[key][loop]['paradic']
                                , flux_intrin, f2, interX_intrin, ax, age_zero=age_zero, Z_zero=Z_zero, f_cont_zero=f_cont_zero, Av_zero=Av_zero)
                Av_gradient = Cal_gradient_map(r_grid, self.r_map, 'Av', self.Avparams[key][loop]['type'],self.Avparams[key][loop]['paradic']
                                , flux_intrin, f2, interX_intrin, ax, age_zero=age_zero, Z_zero=Z_zero, f_cont_zero=f_cont_zero, Av_zero=Av_zero)
                totalflux += flux_band*self.mass_map[key][loop]*age_gradient*Z_gradient*f_cont_gradient*Av_gradient
        if self.EL_region is not None:
            for loop,eldict in enumerate(self.EL_region):
                fc = np.zeros_like(interX)
                prefix = eldict['prefix']
                for line in eldict['Lines']:
                    cwave = par_tot[prefix+line['name']+'_center'].value
                    width = par_tot[prefix+line['name']+'_FWHM'].value*cwave/tl2/c
                    fc += SEDs.gaussian(interX,par_tot[prefix+line['name']+'_Lumin'].value,cwave,width)
                x,fll =  SEDs.sed_to_obse(interX,fc,self.redshift,self.ebv_G)
                flux_band = trapz(fll*f2(interX),x=interX)/ax
                totalflux += flux_band*eldict['normmap']
        return convolve_fft(totalflux,convolve_func)

    def generate_SED_IFU(self,wavelength,resolution=10,sedloopmap = None,dustnorm=None,par_tot=None):
        '''
        gemerate the SED IFU for a galaxy object
        ------
        wavelength: 1D array,
            the wavelength sample
        resolution: int
            number of logr grid to sample SED
        ------
        '''
        ny,nx=self.shape
        tot_IFU = np.zeros((ny,nx,len(wavelength)))
        rmax = np.max(self.r_map)
        r_grid = np.logspace(np.log10(0.5),np.log10(rmax),resolution)
        waveintrin = wavelength/(1.+self.redshift)
        av=3.1*self.ebv_G
        cm=extinction.ccm89(wavelength,av,3.1)
        dimfac = np.power(10,cm/(2.5))*(1+self.redshift)**3
        for key in self.subCs:
            for loop in range(len(self.subCs[key])):
                SED_rgrid = []
                for r in r_grid:
                    age = Cal_map(r,self.ageparams[key][loop]['type'],self.ageparams[key][loop]['paradic'])
                    Z = Cal_map(r,self.Zparams[key][loop]['type'],self.Zparams[key][loop]['paradic'])
                    f_cont = Cal_map(r,self.f_cont[key][loop]['type'],self.f_cont[key][loop]['paradic'])
                    Av = Cal_map(r,self.Avparams[key][loop]['type'],self.Avparams[key][loop]['paradic'])
                    seds_total = SEDs.get_host_SED(waveintrin, 0., f_cont, age, Z, Av, 1.)
                    if dustnorm is not None:
                        seds_total += SEDs.get_nebular99(waveintrin, 0., f_cont, age, Z, dustnorm)
                    SED_rgrid.append(seds_total/dimfac)
                self.subCs[key][loop]['intp']=SED_rgrid
        if sedloopmap is None:
            sedloopmap = np.log10(self.r_map.copy())*resolution/np.log10(rmax)-0.01
            sedloopmap = sedloopmap.astype(int)
        for loopy in range(ny):
            for loopx in range(nx):
                for key in self.subCs:
                    for loop in range(len(self.subCs[key])):
                        tot_IFU[loopy,loopx,:] += self.mass_map[key][loop][loopy][loopx]*self.subCs[key][loop]['intp'][sedloopmap[loopy][loopx]]
        if self.EL_region is not None:
            for loop,eldict in enumerate(self.EL_region):
                fc = np.zeros_like(waveintrin)
                prefix = eldict['prefix']
                for line in eldict['Lines']:
                    cwave = par_tot[prefix+line['name']+'_center'].value
                    width = par_tot[prefix+line['name']+'_FWHM'].value*cwave/tl2/c
                    fc += SEDs.gaussian(waveintrin,par_tot[prefix+line['name']+'_Lumin'].value,cwave,width)
                fll =  fc/dimfac
                prli = np.array(eldict['normmap']).reshape((ny, nx, 1))
                tot_IFU += prli*fll
        return tot_IFU

    def add_extended_EL_region(self, prefix, Pro_names, params, Lines):
        '''
        To add a extended emission line region for a galaxy object

        Pro_names: the name of the profiles
        e.g. = "sersic" "coresersic" "brokenexp" "moffat" "ferrer" "king" "pointsource"

        params: a dictionary of the parameters for this subcomponent
        e.g. for sersic: {'xcen': 50.0, 'ycen': 50.0, 'frac': 0.704, 're': 10.0,
        'nser': 3.0, 'ang': -32.70422048691768, 'axrat': 1.0, 'box': 0.0, 'convolve': False}
        '''
        if self.EL_region is None:
            self.EL_region = []
        eldict = {
            'prefix' : prefix,
            'pname'  : Pro_names,
            'params' : params,
            'normmap': None,
            'Lines'  : Lines,
            'speape' : None,
            }
        self.EL_region.append(eldict)



class AGN(object):
    '''
    the AGN object
    with physical subcomponents and parameters
    '''
    def __init__(self,logM_BH=8.,logLedd=-1.,astar=0., Av=0., z=0., ebv_G=0.,m_BLR=None,m_NLR=None):
        '''
        galaxy object is initialed from a given mass
        '''
        self.logM_BH = logM_BH
        self.logLedd=logLedd
        self.astar = astar
        self.Av = Av
        self.redshift = z
        self.ebv_G = ebv_G
        self.BLR = m_BLR
        self.NLR = m_NLR
        self.sedModel = 'TD'

    def reset_BH(self, logM_BH,logLedd, astar, Av):
        self.logM_BH = logM_BH
        self.logLedd=logLedd
        self.astar = astar
        self.Av = Av
        return

    def change_SED(self, sedmodel):
        self.sedModel = sedmodel

    def set_full_model(self, obsspec, lines_broad, lines_narrow, nbroad=2, nnarrow=1,strict=0.01,broader=12000.,prefix='',**kwargs):
        '''
        set the BLR NLR models for the AGN component
        ------
        obsspec: observed optical spectrum [wave,flux]
            it is used to make an initial guess
        lines_broad: list
            a line list in BLR
        lines_narrow: list
            a line list in NLR
        ------
        '''
        if ('bcentershift' in kwargs.keys()):
            bcentershift = kwargs['bcentershift']
        else:
            bcentershift = 0.
        m_BLR  = lmfit.Model(SEDs.FeII,prefix='FeII{0}'.format(prefix))
        m_BLR += lmfit.Model(SEDs.BaC,prefix='BaC{0}'.format(prefix))
        lab = []
        for line in lines_broad:
            for loop in range(nbroad):
                    m_BLR+=GaussianModel(prefix=line['name']+'b{0}{1}'.format(loop+1,prefix))
            lab.append(np.abs(np.interp(line['wave'],obsspec[0],obsspec[1])*50*np.sqrt(2*np.pi)))
        par_BLR = m_BLR.make_params()
        stdhb = np.median(obsspec[1])
        par_BLR['FeII{0}A_uv'.format(prefix)].set(0.7*stdhb/4000.,min=0.7*stdhb/400000.,max=0.7*stdhb/40.)
        par_BLR['FeII{0}A_op'.format(prefix)].set(stdhb/4000.,min=stdhb/400000.,max=stdhb/40.)
        par_BLR['FeII{0}fwhm'.format(prefix)].set(2000.,min=900.,max=9000.)
        par_BLR['FeII{0}dcen'.format(prefix)].set(0.,min=-3000.,max=3000.)
        par_BLR['BaC{0}cf'.format(prefix)].set(0.2,min=0.,max=1.)
        par_BLR['BaC{0}logM'.format(prefix)].set(expr='1.*agn_logM')
        par_BLR['BaC{0}logMdot'.format(prefix)].set(expr='1.*agn_logLedd')
        par_BLR['BaC{0}spin'.format(prefix)].set(expr='1.*agn_spin')
        par_BLR['BaC{0}dcen'.format(prefix)].set(0.,min=-3000.,max=3000.)
        par_BLR['BaC{0}fwhm'.format(prefix)].set(2000.,min=900.,max=9000.)
        for loopline, line in enumerate(lines_broad):
            la=lab[loopline]
            for loop in range(nbroad):
                if loop ==0:
                    par_BLR['{0}b{1}{2}sigma'.format(line['name'],loop+1,prefix)].set(value=line['wave']*2500./c/tl2,min=line['wave']*1200./c/tl2,max=line['wave']*broader/c/tl2)
                    par_BLR['{0}b{1}{2}amplitude'.format(line['name'],loop+1,prefix)].set(value=la,min=0.,max=100*la)
                    par_BLR['{0}b{1}{2}center'.format(line['name'],loop+1,prefix)].set(value=line['wave'],min=line['wave']-50.,max=line['wave']+50.)
                else:
                    censhi = (bcentershift*loop)/(nbroad-1)
                    par_BLR['{0}b{1}{2}sigma'.format(line['name'],loop+1,prefix)].set(value=line['wave']*5000./c/tl2,min=line['wave']*2000./c/tl2,max=line['wave']*broader/c/tl2)
                    par_BLR['{0}b{1}{2}amplitude'.format(line['name'],loop+1,prefix)].set(value=0.1*la,min=0.,max=10*la)
                    par_BLR['{0}b{1}{2}center'.format(line['name'],loop+1,prefix)].set(value=line['wave']+censhi,min=line['wave']+censhi-200.,max=line['wave']+censhi+200.)
        self.BLR = m_BLR
        nperfix=None
        first=None
        lan=[]
        for line in lines_narrow:
            for loop in range(nnarrow):
                if nperfix is None:
                    nperfix=line
                    first=line
                    m_NLR=GaussianModel(prefix=line['name']+'n{0}{1}'.format(loop+1,prefix))
                else:
                    m_NLR+=GaussianModel(prefix=line['name']+'n{0}{1}'.format(loop+1,prefix))
            lan.append(np.abs(np.interp(line['wave'],obsspec[0],obsspec[1])*5*np.sqrt(2*np.pi)))
        par_NLR = m_NLR.make_params()
        for loopline, line in enumerate(lines_narrow):
            la=lan[loopline]
            if loopline == 0:
                par_NLR['{0}n1{1}center'.format(line['name'],prefix)].set(value=line['wave'],min=(1.-strict)*line['wave'],max=(1.+strict)*line['wave'])
                par_NLR['{0}n1{1}sigma'.format(line['name'],prefix)].set(value=5.,min=line['wave']*100./c/tl2,max=line['wave']*1060./c/tl2)
                par_NLR['{0}n1{1}amplitude'.format(line['name'],prefix)].set(value=la,min=0.,max=100*la)
                for loop in range(nnarrow-1):
                    par_NLR['{0}n{1}{2}center'.format(line['name'],loop+2,prefix)].set(value=line['wave']-5.*(loop+1),min=line['wave']-35.,max=line['wave']+20.)
                    par_NLR['{0}n{1}{2}sigma'.format(line['name'],loop+2,prefix)].set(value=line['wave']*1500./c/tl2,min=line['wave']*200./c/tl2,max=line['wave']*2500./c/tl2)
                    par_NLR['{0}n{1}{2}amplitude'.format(line['name'],loop+2,prefix)].set(value=1.*la/(0.8**(loop+1)),min=0.,max=10*la)
            else:
                for loop in range(nnarrow):
                    par_NLR['{0}n{1}{2}sigma'.format(line['name'],loop+1,prefix)].set(expr='{0}*{1}n{2}{3}sigma'.format(line['wave']/nperfix['wave'],nperfix['name'],loop+1,prefix))
                    par_NLR['{0}n{1}{2}center'.format(line['name'],loop+1,prefix)].set(expr='{0}*{1}n{2}{3}center'.format(line['wave']/nperfix['wave'],nperfix['name'],loop+1,prefix))
                    if loop == 0:
                        if line['name'] == EL.OIII_4959['name']:
                            par_NLR['{0}n{1}{2}amplitude'.format(EL.OIII_4959['name'],loop+1,prefix)].set(expr='0.33557*{0}n{1}{2}amplitude'.format(EL.OIII_5007['name'],loop+1,prefix))
                        elif line['name'] == EL.NII_6549['name']:
                            par_NLR['{0}n{1}{2}amplitude'.format(EL.NII_6549['name'],loop+1,prefix)].set(expr='0.337838*{0}n{1}{2}amplitude'.format(EL.NII_6583['name'],loop+1,prefix))
                        else:
                            par_NLR['{0}n{1}{2}amplitude'.format(line['name'],loop+1,prefix)].set(value=la,min=0.,max=100*la)
                    else:
                        par_NLR['{0}n{1}{2}amplitude'.format(line['name'],loop+1,prefix)].set(expr='1.*{0}n{1}{2}amplitude*{3}n1{2}amplitude/{0}n1{2}amplitude'.format(nperfix['name'],loop+1,prefix,line['name']))
        self.NLR = m_NLR
        return par_BLR,par_NLR

    def generate_image(self, shape, band, convolve_func, psfparams, transpar=None, psftype='psf', par_tot=None):
        '''
        Parameters:
        shape: (y,x) of the output image

        band: band of the output image

        convolve_func: 2D array, the shape of empirical PSF

        {psftype: [psfparams]}: a dict, the point spread function
        eg.  {'psf': [{'xcen':50, 'ycen':50}]}     stands for a point sources which have same shape as the empirical PSF
             {'moffat': [{'xcen':50, 'ycen':50, 'fwhm':3., 'con':'5.'}]}
        '''
        resp = Table.read(filterpath / band,format='ascii')
        ny = shape[0]
        nx = shape[1]
        filter_x=resp['col1']
        filter_y=resp['col2']
        tminx = np.min(filter_x)
        tmaxx = np.max(filter_x)
        interX = np.linspace(tminx,tmaxx,100)
        f2=interp1d(filter_x,filter_y,bounds_error=False,fill_value=0.)
        ax=trapz(f2(interX),x=interX)
        waveintrin = interX/(1.+self.redshift)
        if self.sedModel == 'TD':
            agnsed_rest = SEDs.get_AGN_SED(waveintrin,self.logM_BH,self.logLedd,self.astar,1.)
        elif self.sedModel == 'SP':
            agnsed_rest = SEDs.spowerlaw(waveintrin,10**par_tot['logL5100'].value,par_tot['plC'].value)
        elif self.sedModel == 'BP':
            agnsed_rest = SEDs.brokenpower(waveintrin,10**par_tot['logL5100'].value,par_tot['plC1'].value,par_tot['plC2'].value)
        if self.BLR is None:
            agnsed_rest += 10**intp_BLRTOT((spin,logM,logMdot,waveintrin))
        elif par_tot is not None:
            agnsed_rest += self.BLR.eval(par_tot,x=waveintrin)
        if (self.NLR is not None)&(par_tot is not None):
            agnsed_rest += self.NLR.eval(par_tot,x=waveintrin)
        if self.Av > 0.:
            cm=extinction.ccm89(waveintrin,self.Av,3.1)/2.5
            agnsed_rest /= 10**cm
        x,agnsed = SEDs.sed_to_obse(waveintrin,agnsed_rest,self.redshift,self.ebv_G)
        flux_band = trapz(agnsed*f2(interX),x=interX)/ax
        magzero = 18.
        mag = -2.5*np.log10(flux_band)+magzero
        #print (mag)
        psfparams.update({'mag':mag})
        psfpar_copy=psfparams.copy()
        if transpar is None:
            xpix,ypix = indentify_xy(psfparams['xcen'],psfparams['ycen'])
        else:
            xpix,ypix = coordinates_transfer(psfparams['xcen'],psfparams['ycen'],transpar)
        psfpar_copy['xcen'] = xpix
        psfpar_copy['ycen'] = ypix
        profit_model = {'width':  nx,
            'height': ny,
            'magzero': magzero,
            'psf': convolve_func,
            'profiles': {psftype:[psfpar_copy]}
           }
        agn_map, _ = pyprofit.make_model(profit_model)
        return np.array(agn_map)

    def fiducial_sed(self,wavelength,turnoff=False,par_tot=None):
        waveintrin = wavelength/(1.+self.redshift)
        if not turnoff:
            if self.sedModel == 'TD':
                agnsed_rest = SEDs.get_AGN_SED(waveintrin,self.logM_BH,self.logLedd,self.astar,1.)
            elif self.sedModel == 'SP':
                agnsed_rest = SEDs.spowerlaw(waveintrin,10**par_tot['logL5100'].value,par_tot['plC'].value)
            elif self.sedModel == 'BP':
                agnsed_rest = SEDs.brokenpower(waveintrin,10**par_tot['logL5100'].value,par_tot['plC1'].value,par_tot['plC2'].value)
            if self.BLR is None:
                agnsed_rest += 10**SEDs.intp_BLRTOT((self.astar,self.logM_BH,self.logLedd,waveintrin))
            elif par_tot is not None:
                agnsed_rest += self.BLR.eval(par_tot,x=waveintrin)
            if (self.NLR is not None)&(par_tot is not None):
                agnsed_rest += self.NLR.eval(par_tot,x=waveintrin)
            if self.Av > 0.:
                cm=extinction.ccm89(waveintrin,self.Av,3.1)/2.5
                agnsed_rest /= 10**cm
        else:
            agnsed_rest = np.zeros_like(waveintrin)
            if (self.NLR is not None)&(par_tot is not None):
                agnsed_rest += self.NLR.eval(par_tot,x=waveintrin)
        x,agnsed = SEDs.sed_to_obse(waveintrin,agnsed_rest,self.redshift,self.ebv_G)
        return agnsed

    def generate_SED_IFU(self, wavelength, shape, position, par_BLR=None):
        ny,nx = shape
        tot_IFU = np.zeros((ny,nx,len(wavelength)))
        inty = int(position[1])
        intx = int(position[0])
        waveintrin = wavelength/(1.+self.redshift)
        if self.sedModel == 'TD':
            agnsed_rest = SEDs.get_AGN_SED(waveintrin,self.logM_BH,self.logLedd,self.astar,1.)
        elif self.sedModel == 'SP':
            agnsed_rest = SEDs.spowerlaw(waveintrin,10**par_BLR['logL5100'].value,par_BLR['plC'].value)
        if self.BLR is None:
            agnsed_rest += 10**SEDs.intp_BLRTOT((self.astar,self.logM_BH,self.logLedd,waveintrin))
        elif par_BLR is not None:
            agnsed_rest += self.BLR.eval(par_BLR,x=waveintrin)
        if self.Av > 0.:
            cm=extinction.ccm89(waveintrin,self.Av,3.1)/2.5
            agnsed_rest /= 10**cm
        x,agnsed = SEDs.sed_to_obse(waveintrin,agnsed_rest,self.redshift,self.ebv_G)
        tot_IFU[inty,intx,:]=agnsed
        return tot_IFU

class star(object):
    def __init__(self,xcen,ycen,Teff=5000.,dis=10.):
        self.xcen = xcen
        self.ycen = ycen
        self.Teff = Teff
        self.dis = dis

    def generate_image(self, shape, band, convolve_func, psfparams=None, transpar=None, psftype='psf'):
        '''
        Parameters:
        shape: (y,x) of the output image

        band: band of the output image

        convolve_func: 2D array, the shape of empirical PSF

        {psftype: [psfparams]}: a dict, the point spread function
        eg.  {'psf': [{'xcen':50, 'ycen':50}]}     stands for a point sources which have same shape as the empirical PSF
             {'moffat': [{'xcen':50, 'ycen':50, 'fwhm':3., 'con':'5.'}]}
        '''
        resp = Table.read(filterpath / band,format='ascii')
        ny = shape[0]
        nx = shape[1]
        filter_x=resp['col1']
        filter_y=resp['col2']
        tminx = np.min(filter_x)
        tmaxx = np.max(filter_x)
        interX = np.linspace(tminx,tmaxx,100)
        f2=interp1d(filter_x,filter_y,bounds_error=False,fill_value=0.)
        ax=trapz(f2(interX),x=interX)
        sed_rest = SEDs.Bbody(interX,self.Teff)*(10/self.dis)**2
        flux_band = trapz(sed_rest*f2(interX),x=interX)/ax
        magzero = 18.
        mag = -2.5*np.log10(flux_band)+magzero
        #print (mag)
        psfparams.update({'mag':mag})
        psfpar_copy=psfparams.copy()
        if transpar is None:
            xpix,ypix = indentify_xy(self.xcen,self.ycen)
        else:
            xpix,ypix = coordinates_transfer(self.xcen,self.ycen,transpar)
        psfpar_copy['xcen'] = xpix
        psfpar_copy['ycen'] = ypix
        profit_model = {'width':  nx,
            'height': ny,
            'magzero': magzero,
            'psf': convolve_func,
            'profiles': {psftype:[psfpar_copy]}
           }
        star_map, _ = pyprofit.make_model(profit_model)
        return np.array(star_map)

    def generate_SED_IFU(self, shape, wavelength, convolve_func, psfparams, transpar=None, psftype='psf'):
        ny = shape[0]
        nx = shape[1]
        sed_rest = SEDs.Bbody(wavelength,self.Teff)*(10/self.dis)**2
        magzero = 18.
        mag = magzero
        #print (mag)
        psfparams.update({'mag':mag})
        psfpar_copy=psfparams.copy()
        if transpar is None:
            xpix,ypix = indentify_xy(self.xcen,self.ycen)
        else:
            xpix,ypix = coordinates_transfer(self.xcen,self.ycen,transpar)
        psfpar_copy['xcen'] = xpix
        psfpar_copy['ycen'] = ypix
        profit_model = {'width':  nx,
            'height': ny,
            'magzero': magzero,
            'psf': convolve_func,
            'profiles': {psftype:[psfpar_copy]}
           }
        star_map, _ = pyprofit.make_model(profit_model)
        tot_IFU = np.array(star_map).reshape((ny, nx, 1))
        tot_IFU = sed_rest*tot_IFU
        '''
        tot_IFU = np.zeros((ny,nx,len(wavelength)))
        for loopy in range(ny):
            for loopx in range(nx):
                tot_IFU[loopy,loopx,:] += star_map[loopy][loopx]*sed_rest
        '''
        return tot_IFU

class psfbulge(object):
    def __init__(self,xcen,ycen):
        self.xcen=xcen
        self.ycen=ycen


class PSF(object):
    def __init__(self,xcen,ycen):
        self.xcen=xcen
        self.ycen=ycen

    def generate_image(self,shape,count_rate,convolve_func,psftype='psf',transpar=None):
        '''
        count_rate: I expressed in this way
        '''
        ny,nx=self.shape
        magzero = 18.
        mag = -2.5*np.log10(count_rate)+magzero
        if transpar is None:
            xpix,ypix = indentify_xy(self.xcen,self.ycen)
        else:
            xpix,ypix = coordinates_transfer(self.xcen,self.ycen,transpar)
        psfparams = {'xcen':xpix, 'ycen':ypix}
        psfparams.update({'mag':mag})
        profit_model = {'width':  nx,
            'height': ny,
            'magzero': magzero,
            'psf': convolve_func,
            'profiles': {psftype:[psfparams]}
           }
        psf_map, _ = pyprofit.make_model(profit_model)
        return np.array(psf_map)
