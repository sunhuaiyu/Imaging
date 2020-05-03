# Python scripts to normalize and plot FRAP result 
#

from glob import glob
import pandas as pd; from pandas import Series, DataFrame
import matplotlib.pylab as plt

# intensity normalized to percentile of recovery relative to pre-bleaching
# correct time point 0 to the lowest (time zero post-bleaching)
files = [ i for i in glob('*.txt') ]
for f in files:
    m = pd.read_csv(f, sep='\t').dropna(0)
    intensity_max = m.ix[0, 1]
    intensity_min = m.ix[:50, 1].min()
    m['region1_%'] = m.ix[:, 1]/intensity_max
    zero_time = m.ix[:50, 0][(m.ix[:50, 1] - intensity_min) <= 0.00001]
    m['time_corrected'] = m.ix[:, 0].apply(lambda x: x - zero_time)
    m.to_csv(f, sep='\t', index=False) 
 
# plot time series 
files = [ i for i in glob('*.txt') ]
for f in files:
    m = pd.read_csv(f, sep='\t')
    plt.plot(m.ix[:800, 'time_corrected'], m.ix[:800, 'region1_%'], '.')
    plt.ylim(0, 1)
    plt.savefig(f[:-3] + 'jpg', dpi=600)
    close()

# 20150324: summarize multiple experiments
wt_files = [ i for i in glob('*FLASHwt*.txt') ]
sim_files = [ i for i in glob('*FLASHsim*.txt') ]

wt_all = DataFrame({
     f[:-4]: pd.read_csv(f, sep='\t')['region1_%']
     for f in wt_files}) 

sim_all = DataFrame({
     f[:-4]: pd.read_csv(f, sep='\t')['region1_%']
     for f in sim_files}) 

fig, ax = plt.subplots()
plt.errorbar(arange(750), wt_all.ix[:750, :].mean(1)[:750], 
             xerr=0, yerr=wt_all.ix[:750, :].std(1)[:750], 
             alpha=0.8, linewidth=0.5, capsize=0, fmt='b')
plt.errorbar(arange(750), sim_all.ix[:750, :].mean(1)[:750], 
             xerr=0, yerr=sim_all.ix[:750, :].std(1)[:750], 
             alpha=0.5, linewidth=0.5, capsize=0, fmt='y')             
#plt.plot(wt_all.ix[:750, :], 'b.')
#plt.plot(sim_all.ix[:750, :], 'y.')
plt.xlim(-10, 750); plt.ylim(0, 1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
savefig('20150730_FLASH_summary.pdf')


# curve fitting 
from scipy.optimize import curve_fit

constructs = [['WT', 750], ['sim', 750]]
colors = {'WT': 'b', 'sim': 'y'}   

processed_data = {'WT': wt_all, 'sim': sim_all}
fig, ax = plt.subplots()

print 'construct' + '\t' + 'y0' + '\t' + 't1/2' + '\t' + 'mobile%' +'\tn'

frap = lambda x, plateau, t_half: \
       plateau - (plateau - Y[0]) * 2.0 ** (- x / t_half) 
    
# regression to calculate the half recovery time and mobile fraction       
for i in constructs:
    Y = np.array(processed_data[i[0]].mean(axis=1, skipna=True)[2:i[1]])
    x = np.arange(len(Y))
    p, q = curve_fit(frap, x, Y, p0=[0.2, 40])
    plt.plot(x, frap(x, *p), label=i[0], color=colors[i[0]], linewidth=2)
    print i[0], '\t', np.round(Y[0], 3), '\t', np.round(p[1], 3), '\t', \
         np.round((p[0] - Y[0])/(1 - Y[0]), 3), '\t', processed_data[i[0]].shape[1]

plt.xlim(-10, 400); plt.ylim(0, 1.0)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig('20150730curve_fit.pdf')

'''
construct	y0	t1/2	mobile%	n
WT 		0.094 	80.108 	0.279 	26
sim 	0.097 	60.492 	0.3 	23
'''
