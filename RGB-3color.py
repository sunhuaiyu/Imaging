from pylab import *
from glob import glob

#RGB: splitting in the order of R-G-B
masks = [(1, 0, 0), (0, 1, 0), (0, 0, 1)] 
#merge = {'RG':(1, 1, 0), 'RB':(1, 0, 1), 'GB':(0, 1, 1)} #in case

image_types = ['.tif', '.jpg', '.png', '.bmp', '.gif']

for t in image_types:
    image_files = glob('/Users/sun_huaiyu/Documents/Projects/RGBsplitter/*' + t) #change path
    for f in image_files:
        im = imread(f)
        if im.ndim == 3 and im.shape[2] == 3:
            for i in range(3):
                imsave(f[:-len(t)] + '_' + 'RGB'[i] + t, im * masks[i]) 
        else:
            print(f, im.shape, '......not processed')
