Examples of Config Files
=========================

Multi-band imaging decomposition
---------------------------------

Pure SED fitting
--------------------


Photometric SED fitting involves modeling the SED of an astronomical object using a set of data points derived from broad-band photometry. In GalfitS, this process begins by transforming photometric data into 'mock' images, followed by the standard fitting routine. The photometric data is typically provided in a table format, including columns for the band, flux, and flux error. Here’s an example:

.. list-table::
   :header-rows: 1

   * - Band
     - Flux
     - Flux_err
   * - NUV
     - xxx
     - xxx
   * - FUV
     - xxx
     - xxx
   * - ...
     - ...
     - ...

**Transforming Photometry Data to Mock Images**

To convert this photometric data into mock images, you can use the ``photometry_to_img`` function from the ``galfits.gsutils`` module. Below is an example Python script to perform this transformation:

.. code-block:: python

    from galfits import gsutils
    object_flux = Table.read('./object.mag', format='ascii')
    z = object['redshift']
    Bands = ['galex_nuv', 'galex_fuv', 'sloan_u', 'sloan_g', 'sloan_r', 'sloan_i', 'sloan_z', 'wise_ch1', 'wise_ch2', 'wise_ch3', 'wise_ch4']
    for loop, band in enumerate(Bands):
        flux_mjy = object_flux[band]
        flux_err = object_flux[band + '_err']
        outputname = './data/' + band + '.fits'
        gsutils.photometry_to_img(band, flux_mjy, flux_err, z, outputname, unit='mJy')

This script generates a set of mock images, one for each photometric band, which are then used in the SED fitting process.

**Configuration File for SED Fitting**

An example configuration file for performing SED fitting from the far-ultraviolet (FUV) to the far-infrared (FIR) is available at `SEDfit.lyric <https://github.com/RuancunLi/GalfitS-Public/tree/main/examples/SEDfit.lyric>`_. In this setup, each photometric data point is treated as an image. Below is an excerpt showing how to configure a photometric data point:

.. code-block:: none

    # Image A
    Ia1) [/home/liruancun/Works/liyang/cutimgs/photometry/2sloan_u.fits,0] 
    Ia2) sloan_u      # band
    Ia3) [/home/liruancun/Works/liyang/cutimgs/photometry/2sloan_u.fits,2] 
    Ia4) [/home/liruancun/Works/liyang/cutimgs/photometry/2sloan_u.fits,3] 
    Ia5) 1                   # PSF fine sampling factor relative to data
    Ia6) [Noimg,0] 
    Ia7) cR                # unit of the image
    Ia8) -1 
    Ia9) 1                # Conversion from image unit to flambda, -1 for default
    Ia10) 28.3291              # Magnitude photometric zeropoint
    Ia11) uniform             # sky model
    Ia12) [[0,-0.5,0.5,0.1,0]]                   # sky parameter, (value, min, max, step)
    Ia13) 0 			# allow relative shifting
    Ia14) [[0,-5,5,0.1,0],[0,-5,5,0.1,0]] # [shiftx, shifty]
    Ia15) 1 # Use SED information

**Important Notes:**

- The ``Ia8)`` parameter is set to -1, indicating that the image size is not used in the conventional sense, as this is a mock image for SED fitting.
- The ``Ia9)`` parameter is set to 1, which is essential for converting image units to flux density units (flambda). Using -1 would apply the default conversion, but specifying 1 ensures the correct transformation for SED fitting.

**Fitting Results**

Once the fitting process is complete, the results are saved in a summary file, accompanied by a figure illustrating the SED fitting outcome. The figure below displays the contributions from different components: stars (blue), nebular emission (green), and dust (red).

.. figure:: ./fig/sedfit.jpg
   :align: center

   Example of SED fitting results, showing contributions from stars (blue), nebular emission (green), and dust (red).

**Performance**

The time required to run this SED fitting example varies depending on the machine. Below is a table summarizing approximate fitting times:

.. list-table::
   :header-rows: 1

   * - Machine
     - Time
   * - RTX 4090
     - 1.5 mins


Pure mophology fitting
-----------------------

Spectrum fitting
---------------------

Imaging + spectrum fitting
----------------------------


Grism imaging fitting
-----------------------------