# modified 2014-10-31

from skimage.filter import threshold_otsu as otsu
from scipy.ndimage import label
from pylab import *
from glob import glob

RGB = {'r': 0, 'g': 1, 'b': 2}

def nuclear_dots(file, t=0, c='g'):
    '''
    t: threshold
    c: color, 'r', 'g', 'b'
    '''
    im = imread(file)[:, :, RGB[c]]
    if t == 0:
        t = otsu(im)
    return label(im >= t)

def nb(file, t=240, c='g'):
    im = imread(file)[:, :, RGB[c]]  
    label_dots, n_dots = label(im >= t)
    return n_dots

nb_labels, n = nuclear_dots('20141023U2OS_EGFP-ArkN865_8a.tif')
nb_sizes = array([sum(nb_labels==i) for i in arange(1, n+1)])


# 2015-04-14 try to compare the size of HLB foci formed by FLASHwt and FLASHsim
wt = glob('
sim = 
wt_array = array([]).astype(int)
sim_array = array([]).astype(int)

for i in wt:
    nb_labels, n = nuclear_dots(i)
    nb_sizes = array([sum(nb_labels==i) for i in arange(1, n+1)])
    wt_array = append(wt_array, nb_sizes)

hist([wt_array, sim_array],
      bins=30, normed=0, range=[5, 80], color=['b', 'y'], linewidth=0)