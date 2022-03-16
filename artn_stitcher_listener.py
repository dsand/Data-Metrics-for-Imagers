import sys
import os
import datetime

sys.path.insert(0, '/home/mtnops/git-clones/Data-Metrics-for-Imagers')
from imutils import ARTNreduce

def main(isot='20210227', runtry = True):

    '''
    Loop over the /rts2data/Kuiper/Mont4k/ directories
    Get all the images in each respective directory
    Get all the images in the stitched directory
    Check if the raw image is in the stitched directory
        if it isn't 
            run the imutils.ARTNreduce method
            and save it in the stitched directory
        if it is
            pass
    '''
    
    data_dir = '/rts2data/Kuiper/Mont4k/{}'.format(isot)
    image_dirs = os.listdir(data_dir)
    for i_dir in image_dirs:

        image_dir = data_dir + '/{}'.format(i_dir)
        stitched_dir = image_dir + '/stitched'

        raw_images = [x for x in os.listdir(image_dir) if '.fits' in x]
        stitched_images = [x for x in os.listdir(stitched_dir) if '.fits' in x]

        for r_image in raw_images:
            if r_image not in stitched_images:
                print('Found unstitched image: {}'.format(image_dir + '/{}'.format(r_image)))
                '''
                run the ARTNreduce()
                '''
                if runtry:
                    try:
                        stitched_image = ARTNreduce(image_dir + '/{}'.format(r_image))
                        stitched_image.write(stitched_dir +'/{}'.format(r_image))
                    except:
                        print('Error stitching {}'.format(stitched_dir +'/{}'.format(r_image)))
                else:
                    stitched_image = ARTNreduce(image_dir + '/{}'.format(r_image))
                    stitched_image.write(stitched_dir +'/{}'.format(r_image))


def history_stitch():
    d = datetime.datetime(2021, 1, 1)
    today = datetime.datetime.now()
    while d < today:
        y = d.year
        m = d.month
        if m < 10:
            m = '0{}'.format(m)
        day = d.day
        if day < 10:
            day = '0{}'.format(day)
        isot = '{}{}{}'.format(y,m,day)
        print('Working through {}'.format(isot))
        main(isot=isot)
        d = d + datetime.timedelta(days=1)

if __name__ == '__main__':
    d = datetime.datetime.now()
    if d.hour < 12:
        d = d - datetime.timedelta(days=1)

    year = d.year
    month = str(d.month) if d.month >= 10 else '0{}'.format(d.month)
    day = str(d.day) if d.day >= 10 else '0{}'.format(d.day)
    isot = '{}{}{}'.format(year, month, day)

    print('Stitching images for {}\n'.format(isot))

    main(isot=isot, runtry=False)
