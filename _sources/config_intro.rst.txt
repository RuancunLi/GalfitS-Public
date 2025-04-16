Introduction to GalfitS Config Files
=====================================

Config the data
---------------------

The following code block shows an example of how to configure the data in the input configuration file for GalfitS:

.. code-block:: none

    # Region information
    R1) PG0050+124                                  # name of the target
    R2) [13.3955833333,12.6933888889]               # sky coordinate of the target [RA, Dec]
    R3) 0.061                                       # redshift of the target

    # Image A
    Ia1) [hst_0_F438W_cut.fits,0]                   # input data image (FITS file)
    Ia2) wfc3_f438w                                 # band of the image
    Ia3) [hst_0_F438W_cut.fits,2]                   # sigma image (automatic calculate from data if blank or "none")
    Ia4) [psf_f438w_cv.fits,0]                      # PSF image
    Ia5) 1                                          # PSF fine sampling factor relative to data
    Ia6) [hst_0_F438W_cut.fits,1]                   # bad pixel mask image (use empty mask if blank or "none")
    Ia7) cR                                         # unit of the image
    Ia8) 22.9602936792869                           # size to make cutout image region for fitting, unit arcsec
    Ia9) 1.5061352645923459e+18                     # conversion factor from erg/s/cm^2/A to image unit, -1 for default
    Ia10) 27.0                                      # magnitude photometric zeropoint
    Ia11) uniform                                   # sky model
    Ia12) [[0,-0.5,0.5,0.1,1]]                      # sky parameter, [[value, min, max, step, vary]]
    Ia13) 0                                         # allow relative shifting
    Ia14) [[0,-0.5,0.5,0.1,0],[0,-0.5,0.5,0.1,0]]   # [shiftx, shifty]
    Ia15) 1  

    # Spectrum 
    Sa1) PG0050+124_spec.txt                        # input spectrum file
    Sa2) 1                                          # conversion from spectrum unit to 1e-17 erg/s/cm^2/A
    Sa3) [4000,8500]                                # fitting wavelength range
    Sa4) 0                                          # use high resolution stellar template 

    # Image atlas
    Aa1) hst                                        # name of the image atlas
    Aa2) ['a', 'b']                                 # images in this atlas
    Aa3) 0                                          # whether the images have same pixel size
    Aa4) 0                                          # link relative shiftings
    Aa5) ['a']                                      # spectra
    Aa6) [[2]]                                      # aperture size, fiber-[radius], long slit-[PA, a, b], PA respect to positive RA 
    Aa7) [0]                                        # references images

**Region Configuration**

In GalfitS, configuration files typically begin by defining a region containing the target of interest. The usual format for defining a region involves directly providing information about the target. For instance, to define a region, the user must specify three key parameters, each starting with the letter ``R``, as shown in the code block above.

- ``R1`` names the target, such as PG0050+124, and this name is then used as a prefix for output file names.
- ``R2`` specifies the sky coordinates of the target in Right Ascension and Declination. This coordinate information serves as a central point for future image cutouts in imaging fitting. Thus, when there are multiple objects of interest, a common center can be used for all.
- ``R3`` indicates the target's redshift, which is essential for transforming the image unit to luminosity based on the luminosity distance at this redshift.

**Image Configuration**

For fittings involving images, the image input section of the config files is indicated by the letter ``I``. Since GalfitS typically handles multiple images, each image should be assigned a distinct character, e.g., ``Ia``, ``Ib``, and ``Ic`` for three input images. For each image, a comprehensive set of fifteen parameters is required:

- ``Ia1``: The input data image file, which is a two-element list: the first element is the path of the FITS file, and the second specifies the FITS extension storing the image array.
- ``Ia2``: The band of the image.
- ``Ia3``: The sigma image, following the same format as the data image. If this list's first element is blank or set to "none", GalfitS automatically calculates a sigma image from the data image using the same method as in GALFIT. Accurate sigma image estimation is crucial for the likelihood value and for balancing weights across different images.
- ``Ia4``: The PSF image, in the same format as the data image.
- ``Ia5``: The PSF fine sampling factor relative to the data.
- ``Ia6``: The bad pixel mask image, which uses '1' for masking and '0' for unmasked pixels. If the first element is blank or "none", an empty mask will be used.
- ``Ia7``: The unit of the image; although it does not affect fitting, it serves as a reminder of the image unit and can be ignored if unclear.
- ``Ia8``: The size for the cutout image region intended for fitting. When a single float, e.g., :math:`s`, is provided, the image is cut into a :math:`2s \times 2s` square box, centered on ``R2``, during fitting. If a list is provided, e.g., :math:`[s_1, s_2, s_3, s_4]`, a rectangular cutout spanning :math:`\Delta \text{RA}` from :math:`-s_1` to :math:`s_2` and :math:`\Delta \text{Dec}` from :math:`-s_3` to :math:`s_4`, centered on ``R2``, is used.
- ``Ia9``: The conversion factor from :math:`\text{erg s}^{-1} \text{cm}^{-2} \mathring{\text{A}}^{-1}` to the image unit. If set to -1, GalfitS will use a default value derived from the magnitude zero point of the standard image filter. However, for images like those from 2MASS where the magnitude zero point is not constant, caution is needed with the value of ``Ia9``.
- ``Ia10``: The magnitude photometric zeropoint, used mainly in single-band imaging fitting without SED information (``Ia15``=0). In multi-band image fittings, this parameter can be ignored.
- ``Ia11``: The sky model, where both "uniform" and "polynomial" models are acceptable. These represent the background emission during fitting using either a constant value or a polynomial function. To avoid unnecessary parameter degeneracy, users are advised to first perform a sky subtraction and then fix the sky model to uniformly zero.
- ``Ia12``: The sky model parameters, which is a list equal in length to one plus the polynomial function's order (a uniform sky model is a zeroth-order polynomial). Each list element is a five-value array: [initial value, minimum value, maximum value, typical variation step, fixed or not]. The typical variation step is useful when applying MCMC optimization, and the last value indicates whether the parameter is free (1) or fixed (0). This five-value array format is standard for initializing parameters in GalfitS.
- ``Ia13``: Settings for allowing relative image shifting, typically used when WCS consistency between images, especially from different instruments, is not exact.
- ``Ia14``: The shift parameters in the x and y directions, in pixel units, defined using the standard five-value array.
- ``Ia15``: When set to 0, shifts GalfitS to a pure photometric tool, akin to GALFIT. By default, it is set to 1, which employs SED information for fitting multi-band images.

This meticulous process ensures that GalfitS receives all the necessary information for accurate modeling and analysis of astronomical data.

**Spectrum Configuration**

GalfitS is capable of performing fittings solely with spectra or combined fittings with both images and spectra. To input a spectrum, the configuration process starts with the letter ``S``. For instance, ``Sa`` and ``Sb`` represent the first and second spectra, respectively, when multiple spectra are used. Each spectrum requires four specific parameters:

- ``Sa1``: The input spectrum file, which is a standard text file containing three columns: wavelength in angstroms (:math:`\mathring{\text{A}}`), flux, and flux error.
- ``Sa2``: The conversion factor from the spectrum unit to :math:`1 \times 10^{-17} \text{erg s}^{-1} \text{cm}^{-2} \mathring{\text{A}}^{-1}`. This conversion is crucial to ensure that the spectrum is in the correct unit for analysis.
- ``Sa3``: The fitting wavelength range, in units of :math:`\mathring{\text{A}}`, specifies the spectrum segment to be used in the fitting process.
- ``Sa4``: Decides whether to use a high-resolution stellar template. If set to 1, the high-resolution template from Starburst99 with a resolution of approximately :math:`30 \, \text{km s}^{-1}` is employed, which is essential for deriving stellar velocity dispersion. However, it’s important to note that the Starburst99 high-resolution template covers only the 3000-7000 :math:`\mathring{\text{A}}` range. Therefore, model spectra outside this range will remain at low resolution. For studies focusing on stellar velocity dispersion, it is recommended to set ``Sa3`` to [3000, 7000].

**Image Atlas Configuration**

In GalfitS, the final step in defining the input data is integrating the images and spectra into an image atlas object. Configuring an image atlas begins with the letter ``A`` and involves seven parameters. Like images and spectra, multiple image atlases are distinguished by identifiers such as ``Aa``, ``Ab``, etc.:

- ``Aa1``: Names the image atlas; for example, "hst" could represent the instrument used for the image set.
- ``Aa2``: Lists the images included in this atlas; for instance, if images ``Ia`` and ``Ib`` are included, ``Aa2`` should be ['a', 'b'].
- ``Aa3``: Indicates whether the images in the atlas have the same pixel size, with a value of 0 signifying inconsistencies.
- ``Aa4``: Links relative shifts among the images, enhancing alignment accuracy within the atlas and avoiding unnecessary parameter degeneracy.
- ``Aa5``: Lists the spectra associated with the atlas, similar to the inclusion of images. Notably, both ``Aa2`` and ``Aa5`` can be empty lists, signifying either spectra-only or images-only fitting.
- ``Aa6``: Specifies the aperture/slit formation of the input spectrum. The format varies: a single float indicates a fiber ([radius]), and three floats represent a long slit ([PA, a, b]), with PA aligned with positive RA, and a and b as the length and width of the slit in arcseconds.
- ``Aa7``: Lists reference images. For example, [0] implies that the first image in this image atlas is used for model spectra calculation. It is advisable to use the image with the highest spatial resolution to ensure the most accurate spectrum integration in the forward modeling process.

This structured method of configuring an image atlas in GalfitS ensures a thorough and systematic integration of both imaging and spectral data for astrophysical analysis.


Config the model
----------------

After configuring the data settings in the config file, it is crucial to define the model within the same text file. GalfitS allows for the inclusion of three types of models: galaxies, AGN, and foreground stars. Users can set up these models following specific guidelines, as detailed below.

**Profile and Galaxy Configuration**

The following code block shows an example of configuring a galaxy model in GalfitS:

.. code-block:: none

    # Profile A - sersic profile
    Pa1) bulge                                      # name of the component
    Pa2) sersic                                     # profile type
    Pa3) [0,-5,5,0.1,1]                             # x-center [arcsec]
    Pa4) [0,-5,5,0.1,1]                             # y-center [arcsec]
    Pa5) [2.296,0.023,4.592,0.01,1]                 # Re [arcsec]
    Pa6) [4,1,6,0.1,1]                              # Sersic index
    Pa7) [0,-90,90,1,1]                             # position angle (PA) [degree]
    Pa8) [0.8,0.6,1,0.01,1]                         # axis ratio (q = b/a)
    Pa9) [[-4,-8,-1,0.1,1]]                         # specific star formation rate of the star forming component
    Pa10) [[5,0.01,11,0.1,1]]                       # burst stellar age of the starburst component [Gyr]
    Pa11) [0.02,0.001,0.04,0.001,1]                 # metallicity [Z=0.02=Solar]
    Pa12) [0.7,0.3,5.1,0.1,1]                       # Av dust extinction [mag]
    Pa13) [100,40,200,1,0]                          # stellar velocity dispersion [km/s]
    Pa14) [10.14,8.5,12,0.1,1]                      # log stellar mass [solar mass]
    Pa15) burst                                     # star formation history type: burst/conti/bins
    Pa16) [-2,-4,-2,0.1,0]                          # logU nebular ionization parameter
    Pa26) [3,0,5,0.1,0]                             # amplitude of the 2175A bump on extinction curve
    Pa27) 0                                         # SED model, 0: full; 1: stellar only; 2: nebular only; 3: dust only
    Pa28) [8.14,4.5,10,0.1,1]                       # log cold dust mass
    Pa29) [1.0, 0.1, 50, 0.1, 0]                    # Umin, minimum radiation field
    Pa30) [1.0, 0.47, 7.32, 0.1, 0]                 # qPAH, mass fraction of PAH
    Pa31) [1.0, 1.0, 3.0, 0.1, 0]                   # alpha, powerlaw slope of U
    Pa32) [0.1, 0, 1.0, 0.1, 1]                     # gamma, fraction illuminated by star forming region

    # Profile B - Fourier sersic profile
    Pb1) arm                                        # name of the component
    Pb2) sersic_f                                   # profile type
    Pb3-16) ......                                  # same as the above
    Pb17) [1,0.3,3,0.1,1]                           # r_in, inner radius of rotation [arcsec]
    Pb18) [2,1,10,0.1,1]                            # r_out, outer radius of rotation [arcsec]
    Pb19) [1,0.5,3,0.1,1]                           # alpha, shape of arctan rotation curve
    Pb20) [70,30,130,1,1]                           # theta_out, the maximum rotation angle [degree]
    Pb21) 1                                         # m, number of Fourier mode
    Pb22) [0.5,0.,1,0.01,1]                         # am, Fourier amplitude
    Pb23) [30,0,180,1,1]                            # theta_m, Fourier position angle
    Pb24) [30,0,90,1,1]                             # i_m, Fourier projection angle
    Pb27-32) ......                                 # same as the above

    # Galaxy A
    Ga1) host                                       # name of the galaxy
    Ga2) ['a','b']                                  # profile component included
    Ga3) [0.061,0.011,0.111,0.01,0]                 # galaxy redshift
    Ga4) 0.055                                      # the EB-V of Galactic dust reddening
    Ga5) [1.,0.8,1.2,0.05,1]                        # normalization of spectrum when images+spec fitting
    Ga6) []                                         # narrow lines in nebular
    Ga7) 1                                          # number of components for narrow lines

For inputting a galaxy model in GalfitS, the process starts with defining individual component profiles and subsequently integrating these into a comprehensive galaxy model. Each component profile begins with the letter ``P`` and can require up to 32 parameters, depending on its type.

- ``Pa1``: Names the component (e.g., bulge, disk, bar).
- ``Pa2``: Specifies the profile type, such as ``sersic``, ``ferrer``, ``edgeondisk``, ``GauRing``, or ``const``. An added ``_f`` indicates the use of Fourier modes for non-axisymmetric profiles.
- ``Pa3`` and ``Pa4``: Define the spatial coordinates for the profile's center in arcseconds, relative to the region's center.
- ``Pa5`` to ``Pa8``: Structural parameters that vary based on the profile type. For a Sersic profile, these include effective radius (Re), Sersic index, position angle (PA) in degrees, and axis ratio (q = b/a).
- ``Pa9`` to ``Pa16``: SED-related parameters, including specific star formation rate (sSFR) of the star-forming component, burst stellar age of the starburst component in Gyr, metallicity (Z, where 0.02 is solar), dust extinction (Av) in magnitudes, stellar velocity dispersion in km/s, log stellar mass in solar masses, star formation history type (e.g., burst, continuous, bins), and the nebular ionization parameter (logU).
- ``Pa26``: Amplitude of the 2175A bump on the extinction curve.
- ``Pa27``: SED model type (0: stellar + nebular, 1: stellar only, 2: nebular only, 3: dust only, 4: stellar + nebular + dust).
- ``Pa28`` to ``Pa32``: Parameters for the DL2014 dust model, including log cold dust mass, minimum radiation field (Umin), mass fraction of PAH (qPAH), power-law slope of the radiation field (alpha), and fraction illuminated by the star-forming region (gamma).

For profiles using Fourier modes (e.g., ``sersic_f``), additional parameters like **Pb17** to **Pb24** define the rotation curve and Fourier mode characteristics, such as inner and outer radii of rotation (r_in, r_out), shape of the arctan rotation curve (alpha), maximum rotation angle (theta_out), number of Fourier modes (m), Fourier amplitude (am), Fourier position angle (theta_m), and Fourier projection angle (i_m).

The galaxy model itself is defined starting with the letter ``G`` and requires seven parameters:

- ``Ga1``: Names the galaxy (e.g., 'host').
- ``Ga2``: Lists the profile components included in the galaxy model (e.g., ['a', 'b']).
- ``Ga3``: Defines the galaxy's redshift.
- ``Ga4``: Specifies the EB-V for Galactic dust reddening.
- ``Ga5``: Sets the normalization of the spectrum for combined image and spectrum fitting.
- ``Ga6`` and ``Ga7``: Configure narrow lines in the nebular component and the number of components for narrow lines, respectively.

This structured approach ensures a detailed and accurate representation of the galaxy for astrophysical analysis.

**AGN Configuration**

The following code block shows an example of configuring an AGN model in GalfitS:

.. code-block:: none

    # Nuclei A
    Na1) AGN                                        # name of the nuclei component
    Na2) [0.061,0.011,0.111,0.01,0]                 # redshift of the nuclei
    Na3) 0.055                                      # the EB-V of Galactic dust reddening
    Na4) [0,-5,5,0.1,1]                             # x-center [arcsec]
    Na5) [0,-5,5,0.1,1]                             # y-center [arcsec]
    Na6) [7,5,10,0.1,0]                             # log black hole mass [solar mass]
    Na7) [-1,-4,2,0.1,0]                            # log L/LEdd
    Na8) [0,0,0.99,0.01,0]                          # black hole spin, a star
    Na9) [0,0,3.1,0.1,0]                            # intrinsic reddening Av [mag]
    Na10) [43,41,47,0.1,1]                          # log L5100 [erg/s]
    Na11) [[1,0,4,0.1,1], [0.6, 0, 5, 0.1,1]]       # power law indexes
    Na12) ['Hg','Hb','Ha']                          # broad emission lines
    Na13) ['Hg','Hb','Ha','OIII_4959','OIII_5007','NII_6549','NII_6583','SII_6716','SII_6731'] # narrow emission lines
    Na14) 2                                         # number of Gaussian components for broad lines
    Na15) 2                                         # number of Gaussian components for narrow lines
    Na16) 0                                         # add Balmer continuum
    Na17) 1                                         # add FeII
    Na18) 0                                         # continuum model, power law/broken power law/thin disk/type 2/arbitrary
    Na19) [1.,0.5,2,0.05,0]                         # normalization of spectrum when images+spec fitting
    Na20) 0                                         # add torus, default emitting at [1, 1000] micro
    Na21) [41,39,44,0.1,0]                          # log torus luminosity [erg/s]
    Na22) [-0.5,-2.5,-0.25,0.05,0]                  # torus a: power-law index of the radial dust cloud distribution
    Na23) [0.5,0.25,1.5,0.05,0]                     # torus h: dimensionless scale height
    Na24) [7,5,10,0.5,0]                            # torus N0: number of clouds along an equatorial line-of-sight
    Na25) [15,0,90,5,0]                             # torus i: inclination angle
    Na26) [1.,0.2,5,0.1,1]                          # normalization of images_atlas
    Na27) [40,38,42,0.1,0]                          # log luminosity of the 1200 K black body for Nuclei hot dust

Defining an AGN model in GalfitS starts with the letter ``N`` and requires 27 parameters, each crucial for capturing the AGN's characteristics:

- ``Na1``: Names the nucleus component (e.g., 'AGN').
- ``Na2``: Specifies the redshift of the nucleus.
- ``Na3``: EB-V for Galactic dust reddening.
- ``Na4`` and ``Na5``: Spatial positioning (x and y centers) in arcseconds.
- ``Na6``: Logarithm of the black hole mass in solar masses.
- ``Na7``: Logarithm of the Eddington ratio (L/LEdd).
- ``Na8``: Black hole spin (a star).
- ``Na9``: Intrinsic reddening (Av) in magnitudes.
- ``Na10``: Logarithm of luminosity at 5100 angstroms in erg/s.
- ``Na11``: Power-law indexes for the continuum.
- ``Na12`` and ``Na13``: Lists of broad and narrow emission lines, respectively.
- ``Na14`` and ``Na15``: Number of Gaussian components for broad and narrow lines.
- ``Na16`` to ``Na18``: Options for adding Balmer continuum, FeII emission, and specifying the continuum model (e.g., power law, broken power law).
- ``Na19``: Normalization of the spectrum for combined image and spectrum fitting.
- ``Na20`` to ``Na25``: Parameters for including and configuring a torus model, including its luminosity, power-law index (a), scale height (h), number of clouds (N0), and inclination angle (i).
- ``Na26``: Normalization of the image atlas.
- ``Na27``: Logarithm of the luminosity of the 1200 K black body for Nuclei hot dust in erg/s.


This detailed configuration ensures an accurate portrayal of the AGN for astrophysical research.

**Foreground Star Configuration**

The following code block shows an example of configuring a foreground star model in GalfitS:

.. code-block:: none

    # Foreground star A
    Fa1) star                                       # name of the foreground star
    Fa2) [0,-5,5,0.1,1]                             # x-center [arcsec]
    Fa3) [0,-5,5,0.1,1]                             # y-center [arcsec]
    Fa4) [5500., 300, 50000, 1., 1]                 # effective temperature of the stellar SED
    Fa5) [0, -1, 3, 0.1, 1]                         # logL [Lsun], assuming distance = 1 kpc
    Fa6) [3, 0, 5, 0.1, 1]                          # logg, [cm s^-2], surface gravity
    Fa7) [0.02, 1e-5, 0.04, 0.01, 1]                # Z, metallicity of the atmosphere
    Fa8) 1                                          # use SED information

To configure a foreground star model in GalfitS, the process starts with the letter ``F`` and requires eight parameters:

- ``Fa1``: Names the foreground star (e.g., 'star').
- ``Fa2`` and ``Fa3``: Define the spatial coordinates (x and y centers) in arcseconds.
- ``Fa4``: Effective temperature of the stellar SED in Kelvin.
- ``Fa5``: Logarithm of luminosity (logL) in solar units (Lsun), assuming a distance of 1 kpc.
- ``Fa6``: Surface gravity (logg) in cm s^-2.
- ``Fa7``: Metallicity (Z) of the atmosphere, where 0.02 is solar.
- ``Fa8``: Flag to use SED information (1 for yes, 0 for no).

This configuration allows for the inclusion of foreground stars in the model, accounting for their spectral energy distribution (SED) and spatial properties in the fitting process.


Constrain parameters
---------------------


After executing a GalfitS run, two files are automatically generated in the specified ``savepath`` directory just prior to the fitting process. These files are named ``targetname.params`` and ``targetname.constrain``, where ``targetname`` corresponds to the designation defined in ``R1`` of the config file.

**Parameter File**

The file ``targetname.params`` serves as the parameter file, summarizing all the parameters utilized in the fitting. It is formatted as a machine-readable table with seven columns, accessible, for instance, via the ``Table.read`` function in the Python package ``astropy``. The columns are structured as follows:

- **Column 1**: Lists the names of the parameters, which users are not permitted to edit.
- **Columns 2-6**: Contain the initial value, minimum value, maximum value, typical variation step, and a flag indicating whether the parameter is variable (1) or fixed (0). Users are free to edit these columns, ensuring that the maximum value is always greater than the minimum value.
- **Column 7**: Details the parameter expressions, with the default value for all parameters set to ``None``. If users wish to link parameters, such as aligning the location of the AGN (e.g., ``agn_x``) with the center of the host galaxy (e.g., ``host_xcen``), they can modify the corresponding row in the parameter file. For example:

  .. code-block:: none

      agn_x 0 -1 1 0.1 False 1*host_xcen

  Here, ``1*host_xcen`` in the last column indicates that the value of parameter ``agn_x`` will be fixed to one times the value of parameter ``host_xcen``. This expression effectively links the AGN's center to that of the host galaxy during the fitting process.

**Constraint File**

While editing the expression column in the parameter file facilitates parameter linking, this linking can also be achieved by directly modifying the constraint file. The file named ``targetname.constrain`` is a text file containing the function ``Update_Constraints``, as illustrated below:

.. code-block:: python

    def Update_Constraints(pardictlc):
        pardictlc['HaAGNb1peak'] = 10**(1.157 * pardictlc['AGNlogL5100'] - 46.19) / (2.5 * pardictlc['HaAGNb1wid'])
        pardictlc['HbAGNb1peak'] = 10**(1.133 * pardictlc['AGNlogL5100'] - 45.70) / (2.5 * pardictlc['HbAGNb1wid'])
        pardictlc['HbAGNn1wid']  = 0.970936 * pardictlc['OIII_5007AGNn1wid']
        pardictlc['HbAGNn1cen']  = 0.970936 * pardictlc['OIII_5007AGNn1cen']

This function, ``Update_Constraints``, is designed to update the local parameter dictionary (``pardictlc``) based on user-defined rules, offering a more flexible approach than editing expressions in the parameter file. This flexibility is particularly advantageous when dealing with complex parameter correlations. For example:

- The initial lines establish a link between the luminosities of the broad Hα and Hβ lines and the AGN continuum monochromatic luminosity at 5100 Å, based on the statistical correlation observed in SDSS AGNs (Greene & Ho, 2005, ApJ).
- Subsequent lines create a correlation between the shape of the narrow component of the AGN Hβ profile and the narrow [OIII] λ5007 emission line, which is crucial for reducing degeneracy in optical spectrum fitting.

**Applying Constraints**

By utilizing these linkages in the parameter file and constraint file, GalfitS allows for a more nuanced and interconnected modeling of astrophysical phenomena, ensuring that the fitting process accounts for the inherent complexities and correlations present in astronomical observations. To apply these constraints, users can pass the file names to the GalfitS run with the following command:

.. code-block:: bash

    PYTHON galfitS.py --config filename --readpar pfile --parconstrain cfile

Here, ``--readpar`` is followed by the full path of the parameter file (``pfile``), and ``--parconstrain`` is followed by the full path of the constraint file (``cfile``).


Applying astrophysical priors
------------------------------

Defining parameter linkages in GalfitS is highly flexible, achieved through direct edits to the parameter or constraint files. This flexibility renders GalfitS a more potent and physically relevant tool compared to other SED fitting codes. The decision to link parameters is informed by the scientific questions at hand. Often, adopting statistical correlations to connect parameters can effectively reduce degeneracies, leading to more robust fitting outcomes. However, considering that many astrophysical correlations exhibit scatter, a rigid fixation on these correlations may not always be desirable. To address this, GalfitS offers the capability to transform these simple linkages into astrophysical priors. The configuration of prior files in GalfitS facilitates this transformation, an example of which is depicted below. By adopting this approach, users can seamlessly integrate empirical relationships into their modeling framework. This integration not only ensures that the fitting process is mathematically sound but also firmly anchored in established astrophysical theories and observations.

**Example Prior File**

.. code-block:: none

    # MSR represent mass size relation
    # we follow the Table 1 in Van Der Wel et tal. (2014), where log(Re/kpc) = logA + \alpha*log(M/Ms) - 10.7 + \sigma
    # we apply the contrain by profile name [logA, alpha, sigma], if sigma = 0, the relation is fixed
    #                                                             if sigma < 0, the constrain is not included
    MSRa1) total 
    MSRa2) [0.6, 0.75, 0.1]

    # MMR represent mass metalicity relation
    # adopt the MMR from Kewley & Ellison with KD02 method
    # [a0, a1, a2, a3, sigma], log(O/H) + 12 = a0 + a1*logM + a2*logM^2 + a3*logM^3 /pm sigma
    MMRa1) total 
    MMRa2) [28.0974, -7.23631, 0.850344, -0.0318315]

    # for exponetial increase, [logSFR0, tau, t0]
    # definition: SFR(t) = SFR0 * exp((t-t0)/tau)
    SFHa1) total
    SFHa2) exponetial
    SFHa3) [[1.,-0.5,2,0.1,1],[0.15,-0.2,0.6,0.1,1],[0.65,0.,1.3,0.1,1]]
    SFHa4) [[1.04,0.16],[0.15,0.04],[0.65,0.18]]

    # contrain on AGN
    # in a2) pro [a, b, sigma], which apply Mbh-Mast relation,  we apply logMbh = a + b*logM(pro) + sigma
    # we can also use Mbh-Mbulge relation, e.g. in a2) bulge [-4.18, 1.17, 0.28], from Kormendy & Ho 2013 Equ. 10
    # a3) and a4) comes from linking broad Halpha/beta flux to L5100, from Greene & Ho 2005
    # a7) - a9) are calculatef from Stern & Laor 2012/2013
    AGNa1) myagn # the nuclei component that want to apply constrain
    AGNa2) total # the component which want to apply M-M correlation
    AGNa3) [-7.0, 1.39, 0.79] # the constrain of black hole mass - stellar mass relation, from Table 6 in Greene et al. (2020)
    AGNa4) [-46.19, 1.157, 0.2] # Halpha flux to L5100
    AGNa5) [-45.70, 1.133, 0.2] # Hbeta flux to L5100
    AGNa6) ['Hg','Pag','Pad','Pae','OVI','Lya','NV','SiIV','CIV','HeII','CIII','MgII','Pab','Brg'] # broad emission line included
    AGNa7) [0.35,0.4,0.225,0.15,1.48148148, 7.40740741, 1.48148148, 1.18518519, 3.7037037 , 1.07407407, 1.66666667, 1.66666667,1,0.8] # broad line ratio to Hbeta
    AGNa8) [-34.05, 0.833, 0.36] # OIII flux to 
    AGNa9) ['Hb','OIII_4959','OI_6302','OI_6365','Ha','NII_6583','NII_6549','SII_6716','SII_6731']
    AGNa10) [0.214,0.336,0.03,0.03,0.685,0.302,0.102,0.13,0.13]
    AGNa11) total # constrain the center of AGN, 'none' for not constrain
    AGNa12) 0.6 # The value of R_FeII, defined in Boroson & Green 1992 

    GP) []     # name of parameter which considering a Gaussian prior rather than uniform

**Mass-Size Relation (MSR)**

In GalfitS, one of the first correlations that can be defined is the Mass-Size Relation (MSR), which is applied to specific profiles as indicated by parameters like ``MSRa1``. For example, ``total`` can be one such profile. The array in ``MSRa2``, consisting of three elements [logA, alpha, sigma], defines the slope, intercept, and intrinsic scatter of the MSR, respectively. This implementation is in line with the findings of van der Wel et al. (2014), where the MSR is expressed as:

.. math::
    \log(r_e) = \log A + \alpha (\log M_\ast - 10.7)

with the effective radius :math:`r_e` in kpc and the stellar mass :math:`M_\ast` in solar mass units. For instance, as illustrated in the example prior file, we might adopt the MSR for early-type galaxies at :math:`z \sim 0.25`, as described in van der Wel et al. (2014), to constrain the ``total`` component. In this scenario, ``MSRa2`` would be set to [0.6, 0.75, 0.1]. The intrinsic scatter, here 0.1, is greater than zero, allowing for adjustments during the fitting process and assuming that the scatter follows a Gaussian distribution with a width of 0.1. Setting sigma to zero would enforce an "exact" prior, strictly adhering to the MSR. Conversely, if sigma < 0, the MSR would not be applied, and one could even remove the relation by deleting ``MSRa1`` and ``MSRa2`` from the file. GalfitS offers the flexibility to apply multiple MSRs to different profiles. This is achieved by appropriately setting parameters like ``MSRa1``, ``MSRb1``, ``MSRc1``, and so on for various profiles, allowing for a diverse and comprehensive application of the Mass-Size Relation in the model fitting.

**Mass-Metallicity Relation (MMR)**

The Mass-Metallicity Relation (MMR) represents the second significant correlation that can be integrated into GalfitS. The parameter ``MMRa1`` is designated to specify the profile to which the MMR is applied, such as ``total``, indicating its application to the total component of a galaxy model. The array ``MMRa2``, comprising five elements [a0, a1, a2, a3, sigma], is utilized to define the coefficients of the MMR as well as its intrinsic scatter. This approach aligns with the methodology outlined in Kewley & Ellison (2008), where the metallicity relationship is expressed as:

.. math::
    \log(\text{O/H}) + 12 = a_0 + a_1 \log M_\ast + a_2 (\log M_\ast)^2 + a_3 (\log M_\ast)^3

For example, adopting metallicity estimates based on the method from Kewley & Dopita (2002), the parameter ``MMRa2`` could be set to [28.0974, -7.23631, 0.850344, -0.0318315] (as shown in the example prior file). Similar to the MSR setting, the sigma value in MMR dictates the handling of scatter within the relation. Allowing for a degree of variability in the relation can be achieved by setting sigma > 0, thereby incorporating a realistic range of metallicity values consistent with observational data. Conversely, setting sigma < 0 effectively removes the correlation.

**Star Formation History (SFH) Constraints**

To constrain the star formation history (SFH) of profiles, GalfitS uses ``SFHx`` priors. For example, ``SFHa1`` specifies the profile (e.g., ``total``), and ``SFHa2`` defines the SFH type, such as ``exponential``. The parameters in ``SFHa3`` and ``SFHa4`` provide the initial values, ranges, and priors for the SFH parameters. In the exponential case, ``SFHa3`` might include [logSFR0, tau, t0] The details in shown in  :doc:`Appendix <appendix>` . 

**AGN Constraints**

In GalfitS, parameter linking for AGN models is informed by statistical studies of nearby AGNs, encompassing a total of 12 parameters. As illustrated in the example prior file, ``AGNa1`` identifies the AGN component, such as 'myagn', for the application of constraints. GalfitS facilitates the linking of black hole mass (:math:`M_\text{BH}`) to the mass of a specific profile, which is particularly useful in reducing parameter degeneracy when employing a thin disk model as the AGN continuum model. For example, to adopt a :math:`M_\text{BH} - M_\text{bulge}` correlation, one can set ``AGNa2`` to 'bulge' and ``AGNa3`` to [-4.18, 1.17, 0.28], following the relation from Kormendy & Ho (2013). The array in ``AGNa3`` is similar to ``MSRa2``, with [logA, alpha, sigma], defining the relationship :math:`\log M_\text{BH} = \log A + \alpha \log M_\ast`. The sigma value, which can be set to 0 or a negative number, dictates the exactness or exclusion of this correlation. In cases where there is a single profile for the galaxy, such as 'total', a :math:`M_\text{BH} - M_\text{total}` correlation can be implemented based on Greene et al. (2020), by setting ``AGNa2`` to 'total' and ``AGNa3`` to [-7.0, 1.39, 0.79]. Parameters ``AGNa4`` and ``AGNa5`` link the broad Hα and Hβ fluxes to the monochromatic luminosity at 5100Å (:math:`L_{5100}`). This linkage is beneficial for imaging-only fittings in GalfitS, where the fluxes of broad Hα and Hβ are unknown but significantly contribute to the broadband flux in certain images. Similarly, the luminosities of other broad emission lines can be linked to the luminosity of Hβ, based on statistics from nearby AGNs, as suggested in Greene & Ho (2005). In ``AGNa6``, other broad emission lines can be specified, and ``AGNa7`` sets their luminosity ratios to broad Hβ. The values in the example prior file for these ratios are derived from SDSS type 1 AGNs (Stern & Laor 2012/2013). Since [OIII] λ5007 is a strong emission line in AGN, ``AGNa8`` allows its luminosity to be linked to the continuum luminosity, akin to ``AGNa4`` and ``AGNa5``. The values for this relation in the example prior file are adopted from Stern & Laor (2012). ``AGNa9`` specifies additional narrow emission lines included in the model, and ``AGNa10`` sets their luminosity ratios to [OIII] λ5007, based on Stern & Laor (2013). The parameter ``AGNa11`` is used for center linking, allowing the AGN's location to be aligned with a specific profile, like ``total``. Setting ``AGNa11`` to 'none' implies no center linking. The final parameter, ``AGNa12``, concerns the line ratio between FeII pseudo-continuum and Hβ, as defined in Boroson & Green (1992). This ratio is crucial, particularly for high accretion rate AGNs where FeII emission can be prominent.

**Gaussian Priors**

To apply Gaussian priors to specific parameters, users can list the parameter names in ``GP``. For example:

.. code-block:: none

    GP) ['param1', 'param2']

This will automatically apply Gaussian priors to the listed parameters, rather than uniform priors, allowing for more constrained fitting based on prior knowledge.

**Energy balance**

To apply Energy balance to specific component, users can list the parameter names in ``EB``. For example:

.. code-block:: none

    EB) ['profilename1', 'profilename2']

Here energy balance means the stellar continuum luminosity absorbed in the UV/optical band is equal to the dust luminosity. This is a common assumption in SED fitting e.g. CIGALE, as it allows for a more accurate representation of the energy balance within the galaxy model. By applying this constraint, users can ensure that the energy output from the stellar component is appropriately accounted for in the dust component, leading to a more physically realistic model.

**Applying Priors**

Once the above astrophysical priors are defined, they can be stored in a text file as the prior file, like the example shown above. To apply these priors, users can pass the file name to the GalfitS run with the following command:

.. code-block:: bash

    PYTHON galfitS.py --config filename --priorpath priorfile

Here, ``--priorpath`` is followed by the full path of the prior file (``priorfile``).