import numpy as np
import matplotlib.pyplot as plt
#from astropy.io import fits
#from photutils.isophote import (EllipseGeometry, Ellipse, Isophote, IsophoteList)
#from photutils.isophote.sample import EllipseSample, CentralEllipseSample
#from photutils.isophote.fitter import CentralEllipseFitter
#from astropy import units as u
from astropy.visualization import (MinMaxInterval, LinearStretch, SqrtStretch,
                                   LogStretch, AsinhStretch, ImageNormalize)
from astropy.visualization import quantity_support
from matplotlib.widgets import Slider, Button


__all__ = ['plot_image']


def plot_image(image, pixel_scales=None, stretch='asinh', units='arcsec',
               vmin=None, vmax=None, a=None, ax=None, plain=False, **kwargs):
    '''
    Plot an image.

    Parameters
    ----------
    image : 2D array
        The image to be ploted.
    pixel_scales (optional) : tuple
        Pixel scales along the first and second directions, units: arcsec.
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
    if pixel_scales is None:
        extent = None
        units = 'pixel'
    else:
        nrow, ncol = image.shape
        x_len = ncol * pixel_scales[0].to(units).value
        y_len = nrow * pixel_scales[1].to(units).value
        extent = (-x_len/2, x_len/2, -y_len/2, y_len/2)
    stretchDict = {
        'linear': LinearStretch,
        'sqrt': SqrtStretch,
        'log': LogStretch,
        'asinh': AsinhStretch,
    }
    if a is None:
        stretch_use = stretchDict[stretch]()
    else:
        stretch_use = stretchDict[stretch](a=a)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch_use)

    if ax is None:
        plt.figure(figsize=(7, 7))
        ax = plt.gca()

    if 'cmap' not in kwargs:
        kwargs['cmap'] = 'Greys'
    if 'origin' not in kwargs:
        kwargs['origin'] = 'lower'
    if 'norm' not in kwargs:
        kwargs['norm'] = norm
    if 'extent' not in kwargs:
        kwargs['extent'] = extent
    ax.imshow(image, **kwargs)

    if plain is False:
        ax.minorticks_on()
        ax.set_aspect('equal', adjustable='box')
    return ax

def plotIFU(wavelength, IFU, figsize=(10,10)):
    ny,nx,nwave = IFU.shape
    px = int(0.5*nx)
    py = int(0.5*ny)

    quantity_support()

    fig, ax = plt.subplots(figsize=figsize)
    spec1d = IFU[py,px,:]
    ax.set_ylim([0.8*spec1d.min(),1.6*spec1d.max()])
    line, = ax.plot(
      wavelength, spec1d,
      color='black',
      #label="raw spec"
    )


    # adjust the main plot to make room for the sliders
    plt.subplots_adjust(left=0.20, bottom=0.20)

    # slider
    axx = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    x_slider = Slider(
        ax=axx,
        label='x',
        valmin=0,
        valmax=nx,
        valinit=px,
        valstep=1,
    )

    # vertically oriented slider
    axy = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
    y_slider = Slider(
        ax=axy,
        label="y",
        valmin=0,
        valmax=ny,
        valinit=py,
        valstep=1,
        orientation="vertical",
    )

    def update(val):
        inttx = int(x_slider.val)
        intty = int(y_slider.val)
        spec1d = IFU[intty,inttx,:]
        line.set_ydata(spec1d)
        ax.set_ylim([0.8*spec1d.min(),1.6*spec1d.max()])
        fig.canvas.draw_idle()
    slider_updater = update
    x_slider.on_changed(slider_updater)
    y_slider.on_changed(slider_updater)

    # reset button
    resetax  = fig.add_axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')
    def reset(event):
        x_slider.reset()
        y_slider.reset()
        s_slider.reset()
        w_slider.reset()
    button_reseter = reset
    button.on_clicked(button_reseter)


    axw  = fig.add_axes([0.05, 0.25, 0.0225, 0.63])
    w_slider = w_slider = Slider(
        ax=axw,
        label="w",
        valmin= wavelength.min(),
        valmax= wavelength.max(),
        valinit=wavelength[int(0.5*len(wavelength))],
        orientation="vertical"
    )
    axs = fig.add_axes([0.00, 0.25, 0.0225, 0.63])
    s_slider = Slider(
        ax=axs,
        label="s",
        valmin= 1,
        valmax=100,
        valinit=10,
        orientation="vertical",
    )

    aimage = fig.add_axes([0.70,0.70, 0.30,0.30])
    aimage.set_xticks([])
    aimage.set_yticks([])


    image = IFU[:, :, int(0.5*len(wavelength))]
    sky_std = np.nanstd(image)
    norm=ImageNormalize(vmin=0,vmax=10*sky_std,stretch=AsinhStretch())
    oimage = aimage.imshow(image, origin='lower',norm=norm, cmap='Greys',)
    def update_image(val):
        index = np.argmax(wavelength >= w_slider.val)
        print (index)
        image = IFU[:, :, index]
        sky_std = np.nanstd(image)
        norm=ImageNormalize(vmin=0,vmax=s_slider.val*sky_std,stretch=AsinhStretch())
        oimage.set_data(image)
        oimage.set_norm(norm)
    image_updater = update_image
    w_slider.on_changed(image_updater)
    s_slider.on_changed(image_updater)

    marker = aimage.scatter(0, 0, marker='x', s=20, color='red')
    def update_marker(val):
        x = x_slider.val
        y = y_slider.val
        marker.set_offsets([[x, y]])
    update_marker = update_marker
    update_marker(None)
    x_slider.on_changed(update_marker)
    y_slider.on_changed(update_marker)

    wline = ax.axvline(w_slider.val, ls='--', color='red', alpha=0.5)
    def update_wline(val):
        wline.set_xdata(w_slider.val)
        fig.canvas.draw_idle()
    update_wline = update_wline
    w_slider.on_changed(update_wline)

    #ax.legend(loc='lower right')
    return fig
