import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.wcs.utils import proj_plane_pixel_scales
from .plot import plot_image
from .instrument_info import get_zp
from .utils import get_wcs_rotation
from astropy.visualization import simple_norm, make_lupton_rgb
from .math import Maskellipse,polynomialfit,cross_match
from photutils.segmentation import deblend_sources
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_threshold
from photutils import detect_sources
from photutils.segmentation import SourceCatalog
from astropy.table import Table, Column, join, join_skycoord
from astropy.wcs import WCS
from astropy.nddata import NDData
from photutils.psf import extract_stars
import matplotlib.colors as colors
from photutils import EPSFBuilder

__all__ = ['image', 'image_atlas']
tl2=2*np.sqrt(2*np.log(2))

class image(object):
    '''
    A single image object.

    Functions
    ---------
    * Read from fits file use CCDData.
    * get_size : Get the image size.
    * plot : Plot the image.
    * sigma_clipped_stats : Calculate the basic statistics of the image.
    * set_data : Load from numpy array.
    * set_mask : Set image mask.
    * set_pixel_scales : Set the pixel scales along two axes.
    * set_zero_point : Set magnitude zero point.
    '''
    def __init__(self, filename=None, hdu=0, phys_to_image=None, unit=None, zero_point=None, pixel_scales=None, wcs_rotation=None, mask=None, verbose=True,psfFWHM_map=None):
        '''
        Parameters
        ----------
        filename (optional) : string
            FITS file name of the image.
        hdu : int (default: 0)
            The number of extension to load from the FITS file.
        unit (optional) : string
            Unit of the image flux for CCDData.
        zero_point (optional) : float
            Magnitude zero point.
        pixel_scales (optional) : tuple
            Pixel scales along the first and second directions, units: arcsec.
        wcs_rotation (optional) : float
            WCS rotation, east of north, units: radian.
        mask (optional) : 2D bool array
            The image mask.
        verbose : bool (default: True)
            Print out auxiliary data.
        '''
        if filename is None:
            self.data = None
        else:
            self.data = CCDData.read(filename, hdu=hdu, unit=unit, mask=mask)
            if self.data.wcs and (pixel_scales is None):
                pixel_scales = proj_plane_pixel_scales(self.data.wcs) * u.degree.to('arcsec')

        self.zero_point = zero_point
        if pixel_scales is None:
            self.pixel_scales = None
        else:
            self.pixel_scales = (pixel_scales[0]*u.arcsec, pixel_scales[1]*u.arcsec)

        if self.data.wcs and (wcs_rotation is None):
            self.wcs_rotation = get_wcs_rotation(self.data.wcs)
        elif wcs_rotation is not None:
            self.wcs_rotation = wcs_rotation * u.radian
        else:
            self.wcs_rotation = None
        self.sources_catalog = None
        self.sigma_image = None
        self.sources_skycord = None
        self.ss_data = None
        self.PSF = None
        self.psfFWHM_map=psfFWHM_map
        self.psfFWHM_arcsec = None
        self.cut_image = None
        self.cut_sigma_image = None
        self.cut_mask_image = None
        self.coordinates_transfer_para = None
        self.sky_median = 0.
        self.phys_to_counts_rate = 1.
        self.band = None

    def img_cut(self,ra,dec,cutsize,cutposition=None,gain='CELL.GAIN',extime='EXPTIME',sigma_clipped=True):
        #cutsize is arcsec unit
        #delta_ang point to east direction
        srcPstXY = self.data.wcs.all_world2pix([ra], [dec], 1)
        srcXp = srcPstXY[0][0]
        srcYp = srcPstXY[1][0]
        srcPstXY = self.data.wcs.all_world2pix([ra+1./60], [dec], 1)
        srcXpra = srcPstXY[0][0]
        srcYpra = srcPstXY[1][0]
        srcPstXY = self.data.wcs.all_world2pix([ra], [dec+1./60], 1)
        srcXpdec = srcPstXY[0][0]
        srcYpdec = srcPstXY[1][0]
        dxra = (srcXpra-srcXp)/60.
        dyra = (srcYpra-srcYp)/60.
        dxdec = (srcXpdec-srcXp)/60.
        dydec = (srcYpdec-srcYp)/60.
        #print("Calculated pixel difference by increasing 1-arcsec RA: (X, pix)",dxra," (Y, pix)",dyra)
        #print("Calculated pixel difference by increasing 1-arcsec DEC: (X, pix)",dxdec," (Y, pix)",dydec)
        self.sources_skycord = [srcXp,srcYp]
        ny,nx=self.data.data.shape
        if type(cutsize) is float:
            cutsize_int = int(cutsize/self.pixel_scales[0].value)
            minx = np.max([int(srcXp)-cutsize_int,0])
            maxx = np.min([int(srcXp)+cutsize_int,nx])
            miny = np.max([int(srcYp)-cutsize_int,0])
            maxy = np.min([int(srcYp)+cutsize_int,ny])
        else:
            cutsize_intx1 = int(cutsize[0]/self.pixel_scales[0].value)
            cutsize_intx2 = int(cutsize[1]/self.pixel_scales[0].value)
            cutsize_inty1 = int(cutsize[2]/self.pixel_scales[0].value)
            cutsize_inty2 = int(cutsize[3]/self.pixel_scales[0].value)
            minx = np.max([int(srcXp)-cutsize_intx1,0])
            maxx = np.min([int(srcXp)+cutsize_intx2,nx])
            miny = np.max([int(srcYp)-cutsize_inty1,0])
            maxy = np.min([int(srcYp)+cutsize_inty2,ny])
        if cutposition is not None:
            minx = cutposition[0]
            maxx = cutposition[1]
            miny = cutposition[2]
            maxy = cutposition[3]
        else:
            cutposition = [minx,maxx,miny,maxy]
        imcut = self.data.data[miny:maxy,minx:maxx]
        self.cut_image = imcut
        if sigma_clipped:
            sky_mean, sky_median, sky_std = sigma_clipped_stats(imcut, sigma=3.0, maxiters=5)
        else:
            sky_median = np.nanmedian(imcut)
            sky_std = np.nanstd(imcut)
        self.sky_median = sky_median
        if self.data.mask is not None:
            self.cut_mask_image = self.data.mask[miny:maxy,minx:maxx]
        if self.sigma_image is not None:
            self.cut_sigma_image = self.sigma_image[miny:maxy,minx:maxx]
        else:
            if type(gain) is not str:
                GAIN=gain
            else:
                GAIN = self.data.meta[gain]
            if type(extime) is float:
                expt = extime
            else:
                expt = self.data.meta[extime]
            self.cut_sigma_image = np.sqrt(np.abs(imcut)/GAIN+sky_std**2)
        vector = [srcXpdec-srcXp,srcYpdec-srcYp]
        normv = np.sqrt(vector[0]**2+vector[1]**2)
        delta_ang = np.arcsin((srcYpdec-srcYp)/normv)*180./np.pi
        #x/y0shift was added for relative alignment calibration between different instruments
        self.coordinates_transfer_para = {'x0': srcXp-minx,'x0shift':0.0, 'y0': srcYp-miny,'y0shift':0.0, 'dxra':dxra, 'dxdec':dxdec,'dyra':dyra,'dydec':dydec,'pixsc':self.pixel_scales[0].value, 'delta_ang':delta_ang}
        return imcut,cutposition

    def get_size(self, units='pixel'):
        '''
        Get the size of the image.

        Parameters
        ----------
        units : string
            Units of the size (pixel or angular units).

        Returns
        -------
        x, y : float
            Size along X and Y axes.
        '''
        nrow, ncol = self.data.shape
        if units == 'pixel':
            x = ncol
            y = nrow
        else:
            x = ncol * self.pixel_scales[0].to(units).value
            y = nrow * self.pixel_scales[1].to(units).value
        return (x, y)

    def get_size(self, units='pixel'):
        '''
        Get the size of the image.

        Parameters
        ----------
        units : string
            Units of the size (pixel or angular units).

        Returns
        -------
        x, y : float
            Size along X and Y axes.
        '''
        nrow, ncol = self.data.shape
        if units == 'pixel':
            x = ncol
            y = nrow
        else:
            x = ncol * self.pixel_scales[0].to(units).value
            y = nrow * self.pixel_scales[1].to(units).value
        return (x, y)

    def get_data_info(self):
        '''
        Data information to generate model image.

        Returns
        -------
        d : dict
            shape : (ny, nx)
                Image array shape.
            pixel_scale : (pixelscale_x, pixelscale_y), default units: arcsec
                Pixel scales.
            wcs_rotation : angle, default units: radian
                WCS rotation, east of north.
        '''
        d = dict(shape=self.data.shape,
                 pixel_scale=self.pixel_scale,
                 wcs_rotation=self.wcs_rotation)
        return d

    def sigma_clipped_stats(self, **kwargs):
        '''
        Run astropy.stats.sigma_clipped_stats to get the basic statistics of
        the image.

        Parameters
        ----------
        All of the parameters go to astropy.stats.sigma_clipped_stats().

        Returns
        -------
        mean, median, stddev : float
            The mean, median, and standard deviation of the sigma-clipped data.
        '''
        return sigma_clipped_stats(self.data.data, mask=self.data.mask, **kwargs)

    def plot(self, stretch='asinh', units='arcsec', vmin=None, vmax=None, a=None, ax=None, plain=False, **kwargs):
        '''
        Plot an image.

        Parameters
        ----------
        stretch : string (default: 'asinh')
            Choice of stretch: asinh, linear, sqrt, log.
        units : string (default: 'arcsec')
            Units of pixel scale.
        vmin (optional) : float
            Minimal value of imshow.
        vmax (optional) : float
            Maximal value of imshow.
        a (optional) : float
            Scale factor of some stretch function.
        ax (optional) : matplotlib.Axis
            Axis to plot the image.
        plain : bool (default: False)
            If False, tune the image.
        **kwargs : Additional parameters goes into plt.imshow()

        Returns
        -------
        ax : matplotlib.Axis
            Axis to plot the image.
        '''
        assert self.data is not None, 'Set data first!'
        ax = plot_image(self.data, self.pixel_scales, stretch=stretch,
                        units=units, vmin=vmin, vmax=vmax, a=a, ax=ax,
                        plain=plain, **kwargs)
        if plain is False:
            ax.set_xlabel(r'$\Delta X$ ({0})'.format(units), fontsize=24)
            ax.set_ylabel(r'$\Delta Y$ ({0})'.format(units), fontsize=24)
        return ax

    def plot_direction(self, ax, xy=(0, 0), len_E=None, len_N=None, color='k', fontsize=20,linewidth=2, frac_len=0.1, units='arcsec', backextend=0.05):
        '''
        Plot the direction arrow. Only applied to plots using WCS.
        Parameters
        ----------
        ax : Axis
            Axis to plot the direction.
        xy : (x, y)
            Coordinate of the origin of the arrows.
        length : float
            Length of the arrows, units: pixel.
        units: string (default: arcsec)
            Units of xy.
        '''
        xlim = ax.get_xlim()
        len_total = np.abs(xlim[1] - xlim[0])
        pixelscale = self.pixel_scales[0].to('degree').value
        if len_E is None:
            len_E = len_total * frac_len / pixelscale
        if len_N is None:
            len_N = len_total * frac_len / pixelscale

        wcs = self.data.wcs
        header = wcs.to_header()
        d_ra = len_E * pixelscale
        d_dec = len_N * pixelscale
        ra = [header['CRVAL1'], header['CRVAL1']+d_ra, header['CRVAL1']]
        dec = [header['CRVAL2'], header['CRVAL2'], header['CRVAL2']+d_dec]
        ra_pix, dec_pix = wcs.all_world2pix(ra, dec, 1)
        d_arrow1 = [ra_pix[1]-ra_pix[0], dec_pix[1]-dec_pix[0]]
        d_arrow2 = [ra_pix[2]-ra_pix[0], dec_pix[2]-dec_pix[0]]
        l_arrow1 = np.sqrt(d_arrow1[0]**2 + d_arrow1[1]**2)
        l_arrow2 = np.sqrt(d_arrow2[0]**2 + d_arrow2[1]**2)
        d_arrow1 = np.array(d_arrow1) / l_arrow1 * len_E * pixelscale
        d_arrow2 = np.array(d_arrow2) / l_arrow2 * len_N * pixelscale

        def sign_2_align(sign):
            '''
            Determine the alignment of the text.
            '''
            if sign[0] < 0:
                ha = 'right'
            else:
                ha = 'left'
            if sign[1] < 0:
                va = 'top'
            else:
                va = 'bottom'
            return ha, va
        ha1, va1 = sign_2_align(np.sign(d_arrow1))
        ha2, va2 = sign_2_align(np.sign(d_arrow2))

        xy_e = (xy[0] - d_arrow1[0] * backextend, xy[1] - d_arrow1[1] * backextend)
        ax.annotate('E', xy=xy_e, xycoords='data', fontsize=fontsize,
                    xytext=(d_arrow1[0]+xy[0], d_arrow1[1]+xy[1]), color=color,
                    arrowprops=dict(color=color, arrowstyle="<-", lw=linewidth),
                    ha=ha1, va=va1)
        xy_n = (xy[0] - d_arrow2[0] * backextend, xy[1] - d_arrow2[1] * backextend)
        ax.annotate('N', xy=xy_n, xycoords='data', fontsize=fontsize,
                    xytext=(d_arrow2[0]+xy[0], d_arrow2[1]+xy[1]), color=color,
                    arrowprops=dict(color=color, arrowstyle="<-", lw=linewidth),
                    ha=ha2, va=va2)

    def set_data(self, data, unit):
        '''
        Parameters
        ----------
        data : 2D array
            Image data.
        unit : string
            Unit for CCDData.
        '''
        self.data = CCDData(data, unit=unit)

    def source_detection_individual(self, psfFWHM, nsigma=3.0, sc_key=''):
        '''
        Parameters
        ----------
        psfFWHM : float
            FWHM of the imaging point spread function
        nsigma : float
            source detection threshold
        '''
        data = np.array(self.data.copy())
        psfFWHMpix = psfFWHM / self.pixel_scales[0].value
        thresholder = detect_threshold(data, nsigma=nsigma)
        sigma = psfFWHMpix * gaussian_fwhm_to_sigma
        kernel = Gaussian2DKernel(sigma, x_size=5, y_size=5)
        kernel.normalize()
        segm = detect_sources(data, thresholder, npixels=5, kernel=kernel)
        props = SourceCatalog(data, segm)
        tab = Table(props.to_table(['xcentroid','ycentroid','semimajor_sigma','semiminor_sigma','orientation','segment_flux','fwhm','maxval_xindex','maxval_yindex']))
        self.sources_catalog = tab
        srcPstradec = self.data.wcs.all_pix2world(tab['xcentroid'], tab['ycentroid'],1)
        sc = SkyCoord(srcPstradec[0], srcPstradec[1], unit='deg')
        sctab = Table([sc,np.arange(len(sc))],names=['sc','sloop_{0}'.format(sc_key)])
        self.sources_skycord = sctab

    def make_mask(self,sources=None,magnification=3.):
        '''
        make mask for the extension.

        Parameters
        ----------
        sources : a to-be masked source table (can generate from photutils source detection)
                  if None, will use its own source catalog
        magnification : expand factor to generate mask
        '''
        mask=np.zeros_like(self.data, dtype=bool)
        mask[np.isnan(self.data)] = True
        mask[np.isinf(self.data)] = True
        if sources is None:
            sources = self.sources_catalog
        for loop in range(len(sources)):
            position = (sources['xcentroid'][loop],sources['ycentroid'][loop])
            a = sources['semimajor_sigma'][loop]
            b = sources['semiminor_sigma'][loop]
            theta = sources['orientation'][loop]*180./np.pi
            mask=Maskellipse(mask,position,magnification*a,(1-b/a),theta)
        self.data.mask = mask
        if self.ss_data is not None:
            self.ss_data.mask = mask

    def set_mask(self, mask):
        '''
        Set mask for the extension.

        Parameters
        ----------
        mask : 2D array
            The mask.
        '''
        assert self.data.shape == mask.shape, 'Mask shape incorrect!'
        self.data.mask = mask
        if self.ss_data is not Nont:
            self.ss_data.mask = mask

    def set_pixel_scales(self, pixel_scales):
        '''
        Parameters
        ----------
        pixel_scales (optional) : tuple
            Pixel scales along the first and second directions, units: arcsec.
        '''
        self.pixel_scales = (pixel_scales[0]*u.arcsec, pixel_scales[1]*u.arcsec)

    def set_zero_point(self, zp):
        '''
        Set magnitude zero point.
        '''
        self.zero_point = zp

    def sky_subtraction(self, order=3 , filepath = None):
        '''
        Do polynomial-fitting sky subtraction
        Parameters
        ----------
        order (optional) : int
            order of the polynomial
        '''
        data = np.array(self.data.copy())
        maskplus = self.data.mask.copy()
        backR=polynomialfit(data,maskplus.astype(bool),order=order)
        background=backR['bkg']
        self.ss_data = CCDData(data-background, unit=self.data.unit)
        self.ss_data.mask = maskplus
        if filepath is not None:
            hdu_temp = fits.PrimaryHDU(data-background)
            hdu_temp.writeto(filepath, overwrite=True)

    def read_ss_image(self,filepath):
        '''
        read sky subtracted image from "filepath"
        '''
        hdu = fits.open(filepath)
        self.ss_data = CCDData(hdu[0].data, unit=self.data.unit)
        self.ss_data.mask = self.data.mask.copy()

    def cal_sigma_image(self,filepath=None):
        '''
    	Construct sigma map following the same procedure as Galfit (quadruture sum of sigma at each pixel from source and sky background).
    	Note
        ----------
    	'GAIN' keyword must be available in the image header and ADU x GAIN = electron

    	Parameters
        ----------
    		filepath:
    			Whether and where to save sigma map
        '''
        GAIN = self.data.header['CELL.GAIN']
        if self.ss_data is None:
            raise ValueError(" Please do sky subtration first !!!")
        data = np.array(self.ss_data.copy())
        mask = self.ss_data.mask.copy()
        bkgrms = np.nanstd(data[~mask.astype(bool)])
        data[~mask.astype(bool)] = 0.
        sigmap = np.sqrt(data/GAIN+bkgrms**2)
        self.sigma_image = sigmap
        if filepath is not None:
            hdu_temp = fits.PrimaryHDU(sigmap)
            hdu_temp.writeto(filepath, overwrite=True)

    def read_sigmap(self, filepath):
        '''
        read sigma image from "filepath"
        '''
        hdu = fits.open(filepath)
        self.sigma_image = hdu[0].data

    def read_PSF(self, filepath):
        '''
        read PSF image from "filepath"
        '''
        hdu = fits.open(filepath)
        self.PSF = hdu[0].data

    def set_psf_fwhm(self, fwhm=None):
        if fwhm is None:
            self.psfFWHM_arcsec=self.psfFWHM_map[int(self.sources_skycord[1]),int(self.sources_skycord[0])]
        else:
            self.psfFWHM_arcsec = fwhm
        return

    def gaussian_psf(self, size=15):
        '''
        fwhm should in unit of arcsec
        '''
        psfFWHM = self.psfFWHM_arcsec/self.pixel_scales[0].value
        kernel = Gaussian2DKernel(psfFWHM/tl2, x_size=size, y_size=size)
        kernel.normalize()
        self.PSF = np.array(kernel)
        return self.PSF

class image_atlas(object):
    '''
    Many images.
    '''
    def __init__(self, image_list=None, band_list=None):
        '''
        Parameters
        ----------
        image_list (optional) : List
            List of `image`.
        zp_list (optional) : List
            List of magnitude zeropoint.
        band_list (optional) : List
            List of band name.  Check `instrument_info` for band names.
        '''
        if image_list is None:
            self.image_list = []
        else:
            self.image_list = image_list

        if band_list is None:
            self.band_list = []
        else:
            self.band_list = band_list

        '''
        if (zp_list is None) and (band_list is not None):
            zp_list = []
            for b in band_list:
                zp_list.append(get_zp(b))

            for loop, img in enumerate(self.image_list):
                img.set_zero_point(zp_list[loop])
        '''

        self.__length = len(image_list)
        self.common_catalog = None

    def __getitem__(self, key):
        '''
        Get the image data using the filter name or number index.
        '''
        if type(key) is str:
            idx = self.band_list.index(key)
        elif type(key) is int:
            idx = key
        return self.image_list[idx]

    def __len__(self):
        '''
        Get the length of the data list.
        '''
        return self.__length

    def source_detection(self,nsigma=3.0):
        '''
        Do multi-band source detection
        Parameters
        ----------
        nsigma : float, or a array with same size as image_atlas
            source detection threshold
        '''
        if type(nsigma) == float:
            nsigma = nsigma * np.ones(self.__length,dtype=float)
        for loop in range(self.__length):
            fwhm = self.image_list[loop].psfFWHM_arcsec/self.image_list[loop].pixel_scales[0].value
            self.image_list[loop].source_detection_individual(fwhm,nsigma=nsigma[loop],sc_key=loop+1)
        return

    def make_common_catalog(self,CM_separation=2.5,magnification=3.0,applylist=None):
        '''
        Do multi-band source detection
        Parameters
        ----------
        CM_separation : float
        angular separation used to do sky coordinates crossmatching, unit in deg
        magnification : float, or a array with same size as image_atlas
                        magnification for generating mask foe each image
        applylist : [list of index]
        None for all images
        '''
        if type(magnification) == float:
            magnification = magnification * np.ones(self.__length,dtype=float)
        if applylist is None:
            applylist = np.arange(self.__length)
        cats = []
        for loop in applylist:
            cats.append(self.image_list[loop].sources_skycord)
        comc = cross_match(cats,angular_sep = 2.5)
        lencc = len(comc)
        master_a = np.zeros(lencc, dtype = float)
        master_b = np.zeros(lencc, dtype = float)
        for loop in range(len(comc)):
            a = []
            b = []
            for loop2 in applylist:
                a.append(self.image_list[loop2].sources_catalog['semimajor_sigma'][comc['sloop_{0}'.format(loop2+1)][loop]]
                *magnification[loop2]*self.image_list[loop2].pixel_scales[0].value)
                b.append(self.image_list[loop2].sources_catalog['semiminor_sigma'][comc['sloop_{0}'.format(loop2+1)][loop]]
                *magnification[loop2]*self.image_list[loop2].pixel_scales[0].value)
            master_a[loop] = np.max(np.array(a))
            master_b[loop] = np.max(np.array(b))
        comc.add_column(Column(master_a, name = 'master_a'))
        comc.add_column(Column(master_b, name = 'master_b'))
        self.common_catalog = comc

    def sky_subtraction(self,order=3,filepaths=None):
        '''
        Do multi-band sky subtration
        Parameters
        ----------
        order (optional) : int
            order of the polynomial
        filepaths : filepath to store the sky subtracted images
        '''
        if type(order) == int:
            order = order * np.ones(self.__length,dtype=int)
        for loop in range(self.__length):
            if filepaths is None:
                self.image_list[loop].sky_subtraction(order[loop])
            else:
                self.image_list[loop].sky_subtraction(order[loop],filepath=filepaths[loop])
        return

    def master_mask(self, magnification=3.0, applylist=None):
        '''
        Do multi-band source masking
        Parameters
        ----------
        magnification : float, or a array with same size as image_atlas
                        magnification for generating mask foe each image
        applylist : [list of index]
        None for all images
        '''
        if type(magnification) == float:
            magnification = magnification * np.ones(self.__length,dtype=float)
        if applylist is None:
            applylist = np.arange(self.__length)
        comc = self.common_catalog.copy()
        commonsourcelist = []
        for loop2 in applylist:
            newsc = self.image_list[loop2].sources_catalog.copy()
            for loop in range(len(comc)):
                self.image_list[loop2].sources_catalog['semimajor_sigma'][comc['sloop_{0}'.format(loop2+1)][loop]] = comc['master_a'][loop]/(magnification[loop2]*self.image_list[loop2].pixel_scales[0].value)
                self.image_list[loop2].sources_catalog['semiminor_sigma'][comc['sloop_{0}'.format(loop2+1)][loop]] = comc['master_b'][loop]/(magnification[loop2]*self.image_list[loop2].pixel_scales[0].value)
            indexes = np.delete(np.arange(len(self.image_list[loop2].sources_catalog)), comc['sloop_{0}'.format(loop2+1)])
            newsc.remove_rows(indexes)
            commonsourcelist.append(newsc)
        for loop2 in range(self.__length):
            self.image_list[loop2].make_mask(sources=commonsourcelist[loop2],magnification=magnification[loop2])

    def generate_PSFs(self, equivalent_radius=1.5, threshold = 50., axrit=0.75, size = 20.,oversampling=1, plot=None, filepaths=None):
        '''
        Generate effective point spread fuctions (ePSFs) for each image
        Parameters
        ----------
        equivalent_radius : float,
                            radius criteria to indentify star: fwhm of sourve < equivalent_radius*psfFWHM
        size : float, unit pixel
               use what size box to extract stars
        oversampling : int
                       oversample the ePSF
        plot : None for not plot stars & ePSF
               list like [1,2,3] to plot rgb image
        filepaths : filepath to store the ePSFs
        '''
        stars = self.common_catalog.copy()
        remolist = []
        flux_thre = []
        for loop2 in range(self.__length):
            flux_thre.append(np.percentile(self.image_list[loop2].sources_catalog['segment_flux'],threshold))
        for loop in range(len(stars)):
            for loop2 in range(self.__length):
                catloop = stars['sloop_{0}'.format(loop2+1)][loop]
                a = (self.image_list[loop2].sources_catalog['fwhm'][catloop])*self.image_list[loop2].pixel_scales[0].value/self.image_list[loop2].psfFWHM_arcsec
                ar = self.image_list[loop2].sources_catalog['semiminor_sigma'][catloop]/self.image_list[loop2].sources_catalog['semimajor_sigma'][catloop]
                fls = self.image_list[loop2].sources_catalog['segment_flux'][catloop]
                sigvalue = self.image_list[loop2].sigma_image[self.image_list[loop2].sources_catalog['maxval_yindex'][catloop],self.image_list[loop2].sources_catalog['maxval_xindex'][catloop]]
                if (a > equivalent_radius)|(ar<axrit)|(fls<flux_thre[loop2])|np.isinf(sigvalue):
                    remolist.append(loop)
                    break
        stars.remove_rows(remolist)
        star_images = []
        PSFs = []
        for loop2 in range(self.__length):
            newsc = self.image_list[loop2].sources_catalog.copy()
            indexes = np.delete(np.arange(len(self.image_list[loop2].sources_catalog)), stars['sloop_{0}'.format(loop2+1)])
            newsc.remove_rows(indexes)
            stars_tbl = Table()
            stars_tbl['x']=np.array(newsc['maxval_xindex'])
            stars_tbl['y']=np.array(newsc['maxval_yindex'])
            if self.image_list[loop2].ss_data is not None:
                nddata = NDData(data=np.array(self.image_list[loop2].ss_data))
            else:
                nddata = NDData(data=np.array(self.image_list[loop2].data.data-self.image_list[loop2].sky_median))
            Tstar = extract_stars(nddata, stars_tbl, size=size)
            epsf_builder = EPSFBuilder(oversampling=oversampling, maxiters=15,progress_bar=False)
            epsf, fitted_stars = epsf_builder(Tstar)
            self.image_list[loop2].PSF = epsf.data
            if filepaths is not None:
                hdu = fits.PrimaryHDU(epsf.data.astype('float32'))
                After = fits.HDUList([hdu])
                After.writeto(filepaths[loop2],overwrite= True)
            if plot is not None:
                star_images.append(Tstar)
                PSFs.append(epsf.data)
        #print ('len',len(Tstar),len(stars))
        if plot is not None:
            tlens = len(Tstar)
            if (((tlens//5)+1)*5-tlens) < (((tlens//4)+1)*4-tlens):
                ncols = 5
                nrows = (tlens//5)+1
            else:
                ncols = 4
                nrows = (tlens//4)+1
            fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3*ncols, 3*nrows),squeeze=True)
            ax = ax.ravel()
            for i in range(tlens):
                if len(plot) > 2:
                    star_b = star_images[plot[0]][i].data*100./np.sum(star_images[plot[0]][i].data)
                    star_g = star_images[plot[1]][i].data*100./np.sum(star_images[plot[1]][i].data)
                    star_r = star_images[plot[2]][i].data*100./np.sum(star_images[plot[2]][i].data)
                    norm = simple_norm(star_b, 'log', percent=99.)
                    image = make_lupton_rgb(star_r, star_g, star_b, Q=10)
                else:
                    image = star_images[plot[0]][i].data
                    norm = simple_norm(image, 'log', percent=99.)
                ax[i].imshow(image,norm=norm ,origin='lower')
            plt.show()
            fig=plt.figure(figsize=(10,10))
            if len(plot) > 2:
                star_b = PSFs[plot[0]]*100./np.sum(PSFs[plot[0]])
                star_g = PSFs[plot[1]]*100./np.sum(PSFs[plot[1]])
                star_r = PSFs[plot[2]]*100./np.sum(PSFs[plot[2]])
                norm = simple_norm(star_b, 'log', percent=99.)
                image = make_lupton_rgb(star_r, star_g, star_b, Q=10)
            else:
                image = PSFs[plot[0]]
                norm = simple_norm(image, 'log', percent=99.)
            plt.imshow(image,norm=norm ,origin='lower')
            plt.show()
