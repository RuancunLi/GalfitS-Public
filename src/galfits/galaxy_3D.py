import pyprofit
from scipy.integrate import trapz
from astropy.convolution import convolve_fft
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
import sys,extinction
from . import sed_interp as SEDs

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
        if case('arctan'):
            R = (r-paradic['r0'])/paradic['rt']
            return paradic['v0'] + 2.*paradic['vc']*np.arctan(R)/np.pi
            break
        if case():
            raise ValueError("Unidentified method for calculate gradient map")


class Galaxy3D(object):
    '''
    the galaxy 3D object
    with physical subcomponents and parameters
    '''
    def __init__(self, mass=1e9, z=0., ebv_G=0.):
        '''
        galaxy object is initialed from a given mass
        '''
        self.mass = mass
        self.Nsub=0
        self.subCs = {}
        self.ageparams={}
        self.Zparams={}
        self.f_cont={}
        self.Avparams={}
        self.sigmaparams={}
        self.maglist = []
        self.imshape = None
        self.mass_map = {}
        self.geometry = None
        self.r_map = None
        self.rotation_map=None
        self.redshift = z
        self.ebv_G = ebv_G
    def reset_mass(self,mass):
        '''
        reset the mass of a galaxy object
        '''
        self.mass = mass
    def add_subC(self,Pro_names,params,ageparam,Zparam,f_cont,Avparam,sigmaparam={'type': "const", 'paradic':{'value': 100}}):
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
        params['mag']=params.pop("frac")
        params['mag'] = 10. - 2.5*np.log10(params['mag'])
        if Pro_names in self.subCs.keys():
            self.subCs[Pro_names].append(params)
            self.ageparams[Pro_names].append(ageparam)
            self.Zparams[Pro_names].append(Zparam)
            self.Avparams[Pro_names].append(Avparam)
            self.f_cont[Pro_names].append(f_cont)
            self.sigmaparams[Pro_names].append(sigmaparam)
        else:
            self.subCs.update({Pro_names : [params]})
            self.ageparams.update({Pro_names : [ageparam]})
            self.Zparams.update({Pro_names : [Zparam]})
            self.Avparams.update({Pro_names : [Avparam]})
            self.f_cont.update({Pro_names : [f_cont]})
            self.sigmaparams.update({Pro_names : [sigmaparam]})
            self.mass_map.update({Pro_names : []})
        self.maglist.append(params['mag'])

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
            self.sigmaparams[Pro_names].append(sigmaparam)
        else:
            self.subCs.update({Pro_names : [arbi_para]})
            self.ageparams.update({Pro_names : [ageparam]})
            self.Zparams.update({Pro_names : [Zparam]})
            self.Avparams.update({Pro_names : [Avparam]})
            self.f_cont.update({Pro_names : [f_cont]})
            self.sigmaparams.update({Pro_names : [sigmaparam]})
            self.mass_map.update({Pro_names : []})
        self.mass_map[Pro_names].append(mass_map)



    def geometry_3D(self,geometry):
        '''
        Construct the 3D geometry of a well-defined galaxy
        ---
        geometry : dict
        {
            'xcen' : center of x, # in pix
            'ycen' : center of x, # in pix
            'i'  : inclination angel,      # in degree
            'PA' : position angle of the major axis (c.c. from x-axis), # in degree
            'rotation_curve': a dict, # {'function_name','params'}
        }
        ---
        '''
        self.geometry = geometry
        ny,nx=self.imshape
        xaxis = np.arange(nx)
        yaxis = np.arange(ny)
        phi = geometry['PA']*np.pi/180.
        xmesh, ymesh = np.meshgrid(xaxis, yaxis)
        r_2d_square = (xmesh+0.5 - geometry['xcen'])**2. + (ymesh+0.5 - geometry['ycen'])**2.
        dis_proj_square = np.abs(np.cos(phi)*(ymesh+0.5 - geometry['ycen'])-np.sin(phi)*(xmesh+0.5 - geometry['xcen']))**2
        r_3d = np.sqrt(r_2d_square+(1./np.cos(geometry['i']*np.pi/180.)-1)*dis_proj_square)
        self.r_map = r_3d
        rotation_v = Cal_map(r_3d,geometry['rotation_curve']['function_name'],geometry['rotation_curve']['params'])
        lower = (ymesh+0.5 - geometry['ycen']) < 0.
        phi_map = np.arccos((xmesh+0.5 - geometry['xcen'])/np.sqrt(r_2d_square))
        phi_map[lower] = 2.*np.pi-phi_map[lower]
        rot_map = rotation_v*np.sin(phi_map-phi)*np.sin(geometry['i']*np.pi/180.)
        self.rotation_map = rot_map
        return rot_map

    def generate_mass_map(self,shape,convolve_func):
        '''
        gemerate the mass distribution map for a galaxy object
        shape: return 2D image shape
        convolve_func: a 2D kernel if convolution is needed
        -----
        Caution, in future, the r_map calculation can be used as elliptical_func
        -----
        '''
        mags = np.array(self.maglist)
        magzero = 2.5*np.log10(self.mass/np.sum(np.power(10,mags/(-2.5))))
        ny,nx=shape
        self.imshape = shape
        image = np.zeros(shape,dtype=float)
        xaxis = np.arange(nx)
        yaxis = np.arange(ny)
        xmesh, ymesh = np.meshgrid(xaxis, yaxis)
        for key in self.subCs:
            self.mass_map[key]=[]
            for loop in range(len(self.subCs[key])):
                params = self.subCs[key][loop]
                profit_model = {'width':  nx,
                    'height': ny,
                    'magzero': magzero,
                    'psf': convolve_func,
                    'profiles': {key:[params]}
                   }
                mass_map, _ = pyprofit.make_model(profit_model)
                mass_map = np.array(mass_map)
                mass_map = np.array(mass_map.tolist())
                self.mass_map[key].append(mass_map)
                image += mass_map
        return image

    def generate_SED_IFU(self,wavelength,resolution=10,highres=True):
        '''
        gemerate the SED IFU for a galaxy object
        ------
        wavelength: 1D array,
            the wavelength sample
        resolution: int
            number of logr grid to sample SED
        ------
        '''
        ny,nx=self.imshape
        tot_IFU = np.zeros((ny,nx,len(wavelength)))
        r_grid = np.logspace(np.log10(0.5),np.log10(np.max(self.r_map)),resolution)
        vmax = np.max(self.rotation_map)
        minwave = np.min(wavelength)
        maxwave = np.max(wavelength)
        meangrid = (maxwave-minwave)/len(wavelength)
        wavelarge = np.append(np.arange(minwave*(1.-1.5*vmax/SEDs.c),minwave,meangrid),wavelength)
        wavelarge = np.append(wavelarge,np.arange(maxwave+meangrid,maxwave*(1.+1.5*vmax/SEDs.c),meangrid))
        waveintrin = wavelarge/(1.+self.redshift)
        av=3.1*self.ebv_G
        cm=extinction.ccm89(wavelarge,av,3.1)
        dimfac = np.power(10,cm/(2.5))*(1+self.redshift)**3
        for key in self.subCs:
            for loop in range(len(self.subCs[key])):
                SED_rgrid = []
                for r in r_grid:
                    age = Cal_map(r,self.ageparams[key][loop]['type'],self.ageparams[key][loop]['paradic'])
                    Z = Cal_map(r,self.Zparams[key][loop]['type'],self.Zparams[key][loop]['paradic'])
                    f_cont = Cal_map(r,self.f_cont[key][loop]['type'],self.f_cont[key][loop]['paradic'])
                    Av = Cal_map(r,self.Avparams[key][loop]['type'],self.Avparams[key][loop]['paradic'])
                    if highres:
                        sigma = Cal_map(r,self.sigmaparams[key][loop]['type'],self.sigmaparams[key][loop]['paradic'])
                        seds_total = SEDs.get_host_SED_3D(waveintrin, 0., f_cont, age, Z, sigma, Av, 1.)
                    else:
                        seds_total = SEDs.get_host_SED(waveintrin, 0., f_cont, age, Z, Av, 1.)
                    SED_rgrid.append(seds_total/dimfac)
                intp_sed = RegularGridInterpolator((r_grid,wavelarge), np.array(SED_rgrid),bounds_error=False,fill_value=0.)
                self.subCs[key][loop]['intp']=intp_sed
        for loopy in range(ny):
            for loopx in range(nx):
                wave_loop = wavelength/(1.+self.rotation_map[loopy][loopx]/SEDs.c)
                for key in self.subCs:
                    for loop in range(len(self.subCs[key])):
                        tot_IFU[loopy,loopx,:] += self.mass_map[key][loop][loopy][loopx]*self.subCs[key][loop]['intp']((self.r_map[loopy][loopx],wave_loop))
        return tot_IFU

    def emission_line(self,wavelength,lines,ampmap,sigmap):
        '''
        generate the emission line IFU
        ------
        wavelength: 1D array,
            the wavelength sample
        lines:  (m,) str array
            line names that included in ALLLINES
        ampmap: (m,) list,
            every element is a (ny,nx) array of the luminosity map
        sigmap: a (ny,nx) array, in km/s
            for each emission line component , same sigma adopted
        ------
        '''
        ny,nx=self.imshape
        tot_IFU = np.zeros((ny,nx,len(wavelength)))
        for loop,line in enumerate(lines):
            sig_pix = sigmap*line['wave']/SEDs.c
            cen_pix = (1.+self.rotation_map/SEDs.c)*line['wave']
            tot_IFU += SEDs.gaussian3D(wavelength,ampmap[loop],cen_pix,sig_pix)
        return tot_IFU
