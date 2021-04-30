
import os
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy import stats
from astropy.nddata import CCDData, Cutout2D
from astropy.convolution import Gaussian2DKernel
from astropy.modeling import models, fitting

import photutils
import ccdproc
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astroscrappy import detect_cosmics

def ARTNreduce(filename='', to_fits=True):
    """
    Take a raw FITS image from an ARTN imager (currently only mont4k is supported), perform basic overscan subtraction
    and trimming, and then stitch the images into a single image with the correct sky orientation (E left, N up).
    """
    reduced = []
    hdus = (1, 2)
    fullim = CCDData.read(filename)
    xbin = fullim.header['CCDBIN1']
    ybin = fullim.header['CCDBIN2']
    airm = fullim.header['AIRMASS']

    for h in hdus:
        im = CCDData.read(filename, hdu=h)
        oscansec = im.header['BIASSEC']
        trimsec = im.header['TRIMSEC']
        poly_model = models.Polynomial1D(1)
        im = ccdproc.subtract_overscan(im, fits_section=oscansec, overscan_axis=None,model=poly_model)
        im = ccdproc.trim_image(im, fits_section=trimsec)
        reduced.append(im)

    if xbin != ybin:
        raise Exception("ARTNreduce requires equal binning in both axes.")

    # hacky, hard-coded way to make combined image for mont4k. result is E to left, N up.
    w = reduced[1].wcs.copy()
    ysize, xsize = reduced[1].shape
    w.array_shape = (ysize, xsize*2)
    blank = np.zeros((ysize, xsize*2))
    blank[:, :xsize] = reduced[1].data
    blank[:, xsize:] = np.fliplr(reduced[0].data)
    # line added by PND to do the correct flip!
    blank[:, :] = np.fliplr(blank)
    cr_mask, clean_data = detect_cosmics(blank, sigclip=5., niter=5, cleantype='medmask', psffwhm=30./xbin)
    stitched = CCDData(clean_data, wcs=w, unit=reduced[0].unit)
    #stitched.header['BINNING'] = xbin
    #stitched.header['AIRMASS'] = airm
    stitched.header = fullim.header


    # added by PND
    #if to_fits:
    #    _base = os.path.basename(filename)
    #    _fits = f"{_base.split('.')[0]}_sitched.fits"
    #i    stitched.write(_fits)

    return stitched


def sub_background(image: CCDData, filter_size: int = 27, box_size: int = 150):
    """
    Perform background subtraction using photutils' median background estimator over a 2D mesh.
    """
    binning = image.header['BINNING']
    filter_size = int(filter_size/binning)
    box_size = int(box_size/binning)
    bkg_estimator = photutils.MedianBackground()
    bkg = photutils.Background2D(
        image,
        (box_size, box_size),
        filter_size=(filter_size, filter_size),
        bkg_estimator=bkg_estimator
    )
    sub = image.data - bkg.background
    return sub


def find_sources(
    image: CCDData,
    snr: float = 3.,    # Threshold SNR for segmentation
    fwhm: float = 5.,  # Kernel FWHM for segmentation
    ksize: int = 5,    # Kernel size
    npixels: int = 10   # Number of connected pixels required to be considered a source
):
    """
    Find extended sources in image with default parameters tuned for expected donut size.
    """
    binning = image.header['BINNING']
    fwhm = int(fwhm/binning)
    ksize = int(ksize/binning)
    npixels = int(npixels/binning)
    threshold = photutils.detect_threshold(image, nsigma=snr)
    sigma = fwhm * stats.gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=ksize, y_size=ksize)
    kernel.normalize()
    segm = photutils.detect_sources(image.data, threshold, npixels=npixels, filter_kernel=kernel)
    cat = photutils.source_properties(image.data, segm, wcs=image.wcs)
    return segm, cat

def measure_fwhm(array,displ):
    """Fit a Gaussian2D model to a PSF and return the FWHM

    Parameters
    ----------
    array : numpy.ndarray
        Array containing PSF

    Returns
    -------
    x_fwhm : float
        FWHM in x direction in units of pixels

    y_fwhm : float
        FWHM in y direction in units of pixels
    """
    yp, xp = array.shape
    y, x, = np.mgrid[:yp, :xp]
    p_init = models.Moffat2D(amplitude=np.max(array),x_0=xp/2,y_0=yp/2)
    fit_p = fitting.LevMarLSQFitter()
#    print(fit_p.fit_info)
    fitted_psf = fit_p(p_init, x, y, array,maxiter=50)
#    print(fit_p.fit_info)
    if displ==True:
        plt.figure(figsize=(8, 2.5))
        plt.subplot(1, 3, 1)
        plt.imshow(array, origin='lower', interpolation='nearest')
        plt.title("Data")
        plt.subplot(1, 3, 2)
        plt.imshow(fitted_psf(x, y), origin='lower', interpolation='nearest')
        plt.title("Model")
        plt.subplot(1, 3, 3)
        plt.imshow(array - fitted_psf(x, y), origin='lower', interpolation='nearest')
        plt.title("Residual")
        plt.show()
    return fitted_psf.fwhm

def cutout_sources(
    image: CCDData,
    cat: photutils.segmentation.properties.SourceCatalog,
    size: int = 150,                # cutout size
    buffer: int = 10,                # edge buffer
    saturation: int = 60000,        # when saturation is reached
    max_ellipticity: float = 0.25,  # maximum ellipticity
    nsigma: float = 10.0            # noise threshold for max value
):
    """
    Create cleaned list of donuts and return list of cutouts for the donuts and their mean width
    """
    cutouts = []
    good = []
    fwhms = []
    binning = image.header['BINNING']
    size = int(size/binning)
    buffer = int(buffer/binning)
    minpos = u.pix * (int(size / 2.) + buffer)  # give a bit of a buffer at the edge
    maxpos = image.shape[0] * u.pix - minpos
    mean, median, stddev = stats.sigma_clipped_stats(image, sigma=2.0, maxiters=None)
    count=0
    displ=False
    for s in cat:
        valid_pos = s.xcentroid > minpos and s.xcentroid < maxpos and s.ycentroid > minpos and s.ycentroid < maxpos
        unsaturated = s.max_value < saturation
        is_round = s.ellipticity < max_ellipticity
        if valid_pos and unsaturated  and is_round:
            cutout = Cutout2D(image, (s.xcentroid.value, s.ycentroid.value), (size, size), wcs=image.wcs)
            y_sum = np.sum(cutout.data, axis=0)
            x_sum = np.sum(cutout.data, axis=1)
            y_fwhm = len(np.where(y_sum > np.max(y_sum)/2.)[0])
            x_fwhm = len(np.where(x_sum > np.max(x_sum)/2.)[0])
            fwhm = 0.5 * (y_fwhm + x_fwhm)
            if count > 3.0:displ=False 
            moff_fwhm=measure_fwhm(cutout.data,displ)
#            print(fwhm,moff_fwhm)
            cutouts.append(cutout)
            good.append(s)
            if np.isnan(moff_fwhm)==False:fwhms.append(moff_fwhm)
            count=count+1
     
    mean_fwhm,med_fwhm,std_fwhm = sigma_clipped_stats(fwhms)
    clean_cat = photutils.SourceCatalog(good)
    return (clean_cat, cutouts, med_fwhm,std_fwhm)
