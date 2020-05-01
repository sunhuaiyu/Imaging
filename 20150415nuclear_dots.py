# 2015-04-14 try to compare the size of HLB foci formed by FLASHwt and FLASHsim
from scipy.ndimage import label
from pylab import *
from glob import glob

# 2015-04-14 try to compare the size of HLB foci formed by FLASHwt and FLASHsim
def HLB_foci(files, t):
    '''
    files: file name array from glob
    t: threshold
    '''
    sizes = array([]).astype(int)
    intensity = array([]).astype(float)
    for i in files:
        image = imread(i)
        nb_labels, n = label(image >= t)
        nb_sizes = array([sum(nb_labels==i) for i in arange(1, n+1)])
        sizes = append(sizes, nb_sizes)
        intensity = append(intensity, image[image > t])
    return sizes, intensity.mean()

# read in file names/data sets
wt = glob('*FLASHwt*ax.tif')
sim = glob('*FLASHsim*ax.tif')

t = 30000; min_size = 5                #threshold
wt_sizes, wt_mean = HLB_foci(wt, t)
sim_sizes, sim_mean = HLB_foci(sim, t)

print wt_mean, sim_mean   
hist([wt_sizes, sim_sizes],
      bins=30, normed=0, range=[5, 60], color=['b', 'y'], linewidth=0)
savefig('20150416nb_sizes.pdf')

# maximum intensity projection and noise reduction
wt_stacks = [glob('z*'.join(i.split('Max'))) for i in wt]
sim_stacks = [glob('z*'.join(i.split('Max'))) for i in sim]

wt1a = array([imread(j) for j in wt_stacks[0]])

wt1a[wt1a <=100] = 0  #initial thresholding
wt1a[:, ((wt1a>0).sum(0) <= 2)] = 0   # reduce to zero if only 2 pixels above threshold 

# save as maximum intensity projection
imsave('wt1a.tif', wt1a.max(0), cmap='Greens_r')

