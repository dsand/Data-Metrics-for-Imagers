from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import CCDData, Cutout2D
from imutils import ARTNreduce,sub_background,find_sources,cutout_sources

ast = AstrometryNet()
ast.api_key = 'guisvvvtcwovgney'



def solve(img):

    for i in range(10):print('')
    print('Subtracting overscan and stitching images together.\n\n')
    im = ARTNreduce(img)
    print('Subtracting remaining background with a 2D background mesh\n\n')
    im.data = sub_background(im)
    print('Finding sources in image using photutils segmentation\n\n')
    segm, cat = find_sources(im)
    print('Select individual sources with some filters\n\n')
    print('Fit 2d Moffat to determine FWHM')
    print('-------------------------------')
    clean_cat, cutouts, fwhm, std = cutout_sources(im, cat)
    x,y,ssum=[],[],[]
    for s in clean_cat:
        x.append(s.xcentroid.value)
        y.append(s.ycentroid.value)
        ssum.append(s.source_sum)
    zp = zip(x,y,ssum)
    zp_sort=sorted(zp,key=lambda t: t[2],reverse=True)
    x,y,ssum=zip(*zp_sort)
    print(ssum)
    
#    hdus = (1, 2)
    fullim = CCDData.read(img)
    zra = fullim.header['RA']
    zdec = fullim.header['DEC']
    binning = fullim.header['CCDBIN1']
    c = SkyCoord(zra,zdec,unit=(u.hourangle, u.deg))
    ra,dec=float(c.ra.deg),float(c.dec.deg)

#    for h in hdus:
#        im = CCDData.read(filename, hdu=h)
#        oscansec = im.header['BIASSEC']
#        trimsec = im.header['TRIMSEC']
#        im = ccdproc.subtract_overscan(im, fits_section=oscansec, overscan_axis=None)
#        im = ccdproc.trim_image(im, fits_section=trimsec)
#        reduced.append(im)



    image_width = len(im.data[:,0])
    image_height = len(im.data[0,:])
    print(image_width,image_height,ra,dec)
    print((float(x[0])))
    wcs_header = ast.solve_from_source_list(x, y, image_width, image_height,solve_timeout=120,scale_units='arcsecperpix',scale_type='ev',scale_est=0.14*binning,scale_err=20,parity=2,center_ra=ra,center_dec=dec,radius=1,publicly_visible='n')
    print(wcs_header)




solve('/home/mlundquist/Projects/QA/data/20200529071618-485-RA.fits')
