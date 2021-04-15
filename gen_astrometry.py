from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import CCDData, Cutout2D
from astro_imutils import imreduce,sub_background,find_all_sources,clean_sources

ast = AstrometryNet()
ast.api_key = 'guisvvvtcwovgney'



def solve(img,clean_cat,ra,dec,binning):
    x,y,ssum=[],[],[]
    for s in clean_cat:
        x.append(s.xcentroid.value)
        y.append(s.ycentroid.value)
        ssum.append(s.source_sum)
    zp = zip(x,y,ssum)
    zp_sort=sorted(zp,key=lambda t: t[2],reverse=True)
    x,y,ssum=zip(*zp_sort)
    
    image_width = len(im.data[:,0])
    image_height = len(im.data[0,:])
    wcs_header = ast.solve_from_source_list(x, y, image_width, image_height,solve_timeout=180,scale_units='arcsecperpix',scale_type='ev',scale_est=0.14*binning,scale_err=20,parity=2,center_ra=ra,center_dec=dec,radius=1,publicly_visible='n')
    print(wcs_header)
    return wcs_header



if True:

    infile='/home/mlundquist/Projects/QA/data/20200529071618-485-RA.fits'

    fullim = CCDData.read(infile)
    xbin = fullim.header['CCDBIN1']
    ybin = fullim.header['CCDBIN2']
    airm = fullim.header['AIRMASS']
    n_ext = fullim.header['NEXTEND']
    c = SkyCoord(fullim.header['RA'],fullim.header['DEC'],unit=(u.hourangle, u.deg))
    ra,dec=float(c.ra.deg),float(c.dec.deg)


    for h in range(1,1+n_ext):
        print('HDU',h)
        im = CCDData.read(infile, hdu=h)
        im.header['BINNING'] = xbin
        red = imreduce(im)
        segm,cat = find_all_sources(red,xbin)
        clean_cat = clean_sources(im, cat)
        wcs_header=solve(red,clean_cat,ra,dec,xbin)

