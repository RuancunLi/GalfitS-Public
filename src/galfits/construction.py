from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils import Background2D, MedianBackground
from photutils.segmentation import make_source_mask
from reproject import reproject_exact
from reproject.mosaicking import reproject_and_coadd
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata.utils import Cutout2D
from astropy.visualization import SqrtStretch, LogStretch, AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
norm = ImageNormalize(stretch=LogStretch())



def construct_psf(image='0_i.fits', index=0, star_table='0_stars_table.ipac', output_dir='output/',
				  save_star_cutout=True, save_bkg=True, star_size=33, psf_size=31):
	'''
	construct psf from a list of stars using python package reproject
	:param image: .fits file of an image
	:param index: index of image data in hdu file
	:param star_table: .ipac file containing 'x' and 'y' columns as peaks of stars
	:param output_dir: directory to save results
	:param save_star_cutout: whether to save star cutouts
	:param star_size: size of star cutout to construct psf (should be odd number)
	:param psf_size: size of psf image (should be odd number and <= size of star cutout)
	:return: 2D array of psf
	'''
	print('Start stack stars for: ', image)
	# ---------------------------------------------------------------
	# load the image
	hdu = fits.open(image)
	w = WCS(hdu[index].header)
	data = hdu[index].data
	hdu.close()
	# ---------------------------------------------------------------
	# background properties
	# ---------------------------------------------------------------
	# mask sources
	mask = make_source_mask(data, 5, 5, dilate_size=11)

	# subtract 2D background
	sigma_clip = SigmaClip(sigma=3.)
	bkg_estimator = MedianBackground()
	bkg = Background2D(data, (60, 60), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=mask)
	data = data - bkg.background

	if save_bkg:
	# plot source mask and bkg
		plt.imshow(bkg.background, origin='lower', cmap='Greys_r',
			   interpolation='nearest')
		plt.colorbar()
		plt.savefig(output_dir+image+'_bkg.pdf', bbox_inches='tight')
		plt.close()

	# ---------------------------------------------------------------
	# load psf star table
	stars_tbl = Table.read(star_table, format='ascii.ipac')

	hdu_List = []
	for each in range(len(stars_tbl)):
		new_wcs = w
		new_wcs.wcs.crval = (0.0000, 0.0000)
		new_wcs.wcs.crpix = [stars_tbl[each]['x']+1, stars_tbl[each]['y']+1]

		cut = Cutout2D(data, (stars_tbl[each]['x'], stars_tbl[each]['y']), star_size, wcs=new_wcs)
		hdu_temp = fits.PrimaryHDU(cut.data, header=cut.wcs.to_header())
		hdu_temp.header['skystd'] = bkg.background_rms_median
		hdu_temp.header['x'] = stars_tbl[each]['x']
		hdu_temp.header['y'] = stars_tbl[each]['y']
		# save star images to fits files
		if save_star_cutout:
			hdu_temp.writeto(output_dir+'{0}_stars_{1}.fits'.format(image.split('.fits')[0], each), overwrite=True)
		hdu_List.append(hdu_temp)

	AGN_wcs = WCS(hdu_List[0].header)
	AGN_wcs.wcs.crval = (0.000000, 0.000000)
	AGN_wcs.wcs.crpix = ((psf_size+1)/2, (psf_size+1)/2)
	#
	# Use reproject package to generate psf
	array, footprint = reproject_and_coadd(hdu_List, AGN_wcs, shape_out=(psf_size, psf_size), reproject_function=reproject_exact,
										   combine_function='sum')
	# save psf to fits file
	hdu_temp = fits.PrimaryHDU(array)
	hdu_temp.writeto(output_dir+'{}_psf.fits'.format(image.split('.fits')[0]), overwrite=True)
	return array

def construct_sigma_map(image='0_i.fits', index=0, output_dir='output/', dilate_size=11, box_size=(60, 60), plot_sigma_map=True, save_sigma_map=True):
	'''
		Construct sigma map following the same procedure as Galfit (quadruture sum of sigma at each pixel from source and sky background).

	Note
    ----------
	'GAIN' keyword must be available in the image header and ADU x GAIN = electron

	Parameters
    ----------
		image: string
			.fits file of an image
		index: int
			index of image data in hdu file
		output_dir: string
			directory to save results
		dilate_size: int
			The size of the square array used to dilate the segmentation image for the detected source.
		box_size: int or array_like (int)
			The box size along each axis.  If ``box_size`` is a scalar then
			a square box of size ``box_size`` will be used.  If ``box_size``
			has two elements, they should be in ``(ny, nx)`` order.  For
			best results, the box shape should be chosen such that the
			``data`` are covered by an integer number of boxes in both
			dimensions.  When this is not the case, see the ``edge_method``
			keyword for more options.
		plot_sigma_map: bool
			Whether to plot sigma map
		save_sigma_map: bool
			Whether to save sigma map
	Return
	------------
	np.ndarray
		sigma image with the same size of input image
	'''

	# Load image
	hdu = fits.open(image)
	GAIN = hdu[index].header['GAIN']
	data = hdu[0].data
	hdu.close()

	# Generate source mask
	mask = make_source_mask(data, 5, 5, dilate_size=dilate_size)
	sigma_clip = SigmaClip(sigma=3.)
	bkg_estimator = MedianBackground()
	# Generate 2D background
	bkg = Background2D(data, box_size, filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
					   mask=mask)
	# Subtract background from image
	data = data - bkg.background

	# mask background
	data[~mask] = 0
	# background rms
	bkgrms = bkg.background_rms

	# sigma map
	sigmap = np.sqrt(data/GAIN+bkgrms**2)

	# plot sigma map
	if plot_sigma_map:
		plt.imshow(sigmap, origin='lower', cmap='Greys_r', norm=norm,
				   interpolation='nearest')
		plt.colorbar()
		plt.savefig(output_dir + image + '_sigma_map.pdf', bbox_inches='tight')
		plt.close()

	if save_sigma_map:
		hdu_temp = fits.PrimaryHDU(sigmap)
		hdu_temp.writeto(output_dir + '{}_sigma_map.fits'.format(image.split('.fits')[0]), overwrite=True)

	return sigmap
