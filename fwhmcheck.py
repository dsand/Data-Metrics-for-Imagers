import argparse
import multiprocessing
import warnings
from multiprocessing import Pool, Manager
from multiprocessing import Queue, Process, cpu_count
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time,os
from imutils import ARTNreduce,sub_background,find_sources,cutout_sources
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch,PercentileInterval
from astropy.visualization.mpl_normalize import ImageNormalize





#####   TODO needs to check that file as finished downloading before loading.






def get_source_file(file):
    """
    Perform data reduction and calculate focus offset from pupil size of detected stars
    """
    for i in range(10):print('')
    print('Subtracting overscan and stitching images together.\n\n')
    im = ARTNreduce(file)
    print('Subtracting remaining background with a 2D background mesh\n\n')
    im.data = sub_background(im)
    print('Finding sources in image using photutils segmentation\n\n')
    segm, cat = find_sources(im)
    print('Select individual sources with some filters\n\n')
    print('Fit 2d Moffat to determine FWHM')
    print('-------------------------------')
    clean_cat, cutouts, fwhm, std = cutout_sources(im, cat)
    print(clean_cat)
    print('FWHM pix:',fwhm,'err:',std)
    bin = im.header['BINNING']
    print('FWHM arcsec:',fwhm*float(bin)*.14,'err:',std*float(bin)*.14)
    print('Zenith Corrected FWHM',fwhm*float(bin)*.14*(float(im.header['AIRMASS'])**(-0.6)),'err:',std*float(bin)*.14*(float(im.header['AIRMASS'])**(-0.6)))
#    print(im.header)
    view_cutouts=False
    if view_cutouts==True:
        for i in range(len(cutouts)):
            fix,ax=plt.subplots(1,1,figsize=(10,10))
            ax.imshow(cutouts[i].data,origin='lower',cmap='gray',vmin=-50,vmax=100)
            plt.show()

    view_image=False
    if view_image==True:
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(16,8))
 
        ax1.imshow(im,origin='lower',cmap='gray',vmin=-50,vmax=100)
        for s in clean_cat:
            ax1.scatter(s.xcentroid,s.ycentroid,marker=0,color='cyan',s=100)
        ax2.imshow(segm,origin='lower')
        plt.tight_layout()
        plt.show()

def get_source_path(event):
    """
    Perform data reduction and calculate focus offset from pupil size of detected stars
    """
    try:
        try:
            file = str(event.dest_path)
        except AttributeError:
            file = str(event.src_path) #get name of new file
    except AttributeError: #if event is a file
            file = event
    print(file)
    if file.split('.')[-1] == 'fits':
        for i in range(10):print('')
        print('Subtracting overscan and stitching images together.\n\n')
        im = ARTNreduce(file)
        print('Subtracting remaining background with a 2D background mesh\n\n')
        im.data = sub_background(im)
        print('Finding sources in image using photutils segmentation\n\n')
        segm, cat = find_sources(im)
        print('Select individual sources with some filters\n\n')
        print('Fit 2d Moffat to determine FWHM')
        print('-------------------------------')
        clean_cat, cutouts, fwhm, std = cutout_sources(im, cat)
        print(clean_cat)
        print('FWHM pix:',fwhm,'err:',std)
        bin = im.header['BINNING']
        print('FWHM arcsec:',fwhm*float(bin)*.14,'err:',std*float(bin)*.14)
        print('Zenith Corrected FWHM',fwhm*float(bin)*.14*(float(im.header['AIRMASS'])**(-0.6)),'err:',std*float(bin)*.14*(float(im.header['AIRMASS'])**(-0.6)))
    #    print(im.header)
        view_cutouts=False
        if view_cutouts==True:
            for i in range(len(cutouts)):
                fix,ax=plt.subplots(1,1,figsize=(10,10))
                ax.imshow(cutouts[i].data,origin='lower',cmap='gray',vmin=-50,vmax=100)
                plt.show()
 
        view_image=False
        if view_image==True:
            fig,(ax1,ax2)=plt.subplots(1,2,figsize=(16,8))
    
            ax1.imshow(im,origin='lower',cmap='gray',vmin=-50,vmax=100)
            for s in clean_cat:
                ax1.scatter(s.xcentroid,s.ycentroid,marker=0,color='cyan',s=100)
            ax2.imshow(segm,origin='lower')
            plt.tight_layout()
            plt.show()




class FileWatcher(FileSystemEventHandler,object):

    def __init__(self, queue): #parameters needed for action
        self._queue = queue

    def on_created(self, event):
        '''Action to take for new files.

        :param event: new event found
        :type event: event
        '''
        file_size = -1
        while file_size != os.path.getsize(event.src_path):
            file_size = os.path.getsize(event.src_path)
            time.sleep(1)

        self._queue.apply_async(get_source_path,[event])


warnings.filterwarnings('ignore', category=UserWarning, append=True)
parser = argparse.ArgumentParser()
parser.add_argument("-f", dest='file', help="Run on a file")
parser.add_argument("-w", dest='watchpath', help="Run on new files in directory")
args = parser.parse_args()

if args.file:
    get_source_file(args.file)
elif args.watchpath:
    pool = Pool(20) #create pool with given CPUs and queue feeding into action function 
    observer = Observer() #create observer
    observer.schedule(FileWatcher(pool), args.watchpath, recursive=True) #setup observer
    observer.start() #start observe
    while True:
        time.sleep(1)

else:
    print('Not Valid')
    print('e.g. python fwhmcheck.py -f /home/mlundquist/Projects/QA/data/20200529071618-485-RA.fits')
    print('e.g. python fwhmcheck.py -w /home/mlundquist/Projects/QA/data/')


