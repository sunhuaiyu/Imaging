from glob import glob
import pandas as pd; from pandas import Series, DataFrame
import matplotlib.pylab as plt
from pylab import *

# intensity normalized to percentile of recovery relative to pre-bleaching
# correct time point 0 to the lowest (time zero post-bleaching)
for f in glob('201505*FRAP*.txt'):
    m = pd.read_csv(f, sep='\t').dropna(0)
    intensity_max = m.ix[:20, 1].max()
    m['region1_%'] = m.ix[:, 1]/intensity_max
    zero_timepoint = m.ix[m.ix[:20, 1].argmin(), 0]
    m['time_corrected'] = m.ix[:, 0].apply(lambda x: x - zero_timepoint)
    m.to_csv(f, sep='\t', index=False) 

# plot the time series
for f in glob('201505*FRAP*.txt'):
    m = pd.read_csv(f, sep='\t')
    plt.plot(m.ix[:300, 'time_corrected'], m.ix[:300, 'region1_%'], '.')
    plt.ylim(0, 1)
    plt.savefig(f[:-3] + 'jpg')
    plt.close()

# summarize all ArkN FRAP data got so far (2015-05-19)
# read normalized numbers from multiple experiments
base_path = '~/Documents/Projects/RING_SUMO/Raw_Data/'

raw_data = { 
'N725': glob(base_path + '20150518FRAP_ArkNsim_data/*N725_A4*.txt')
      + glob(base_path + '20150407U2OS_EGFP-ArkN_FRAP_summary/*N725*.txt'), 
'N725sim': glob(base_path + '20150518FRAP_ArkNsim_data/*N725sim_A5*.txt'),
'N865': glob(base_path + '20150518FRAP_ArkNsim_data/*N865_A1*.txt') 
      + glob(base_path + '20150407U2OS_EGFP-ArkN_FRAP_summary/*N865*.txt'),
'N865sim': glob(base_path + '20150518FRAP_ArkNsim_data/*N865sim_A2*.txt'),
'N865simhis': glob(base_path + '20150518FRAP_ArkNsim_data/*N865simhis_A3*.txt') }

# make sure the zero time point is at the minimum value
processed_data = dict()
for construct, files in raw_data.items():
    container = dict()
    for f in files:
        a = np.array(pd.read_csv(f, sep='\t')['region1_%'])
        container[f] = Series(a[a.argmin() - 2:])
    processed_data[construct] = DataFrame(container)

# color scheme for plotting
colors = {'N865': 'y', 'N865sim': 'c', 'N865simhis': 'g', 'N725': 'b', 'N725sim': 'r'}   

# plot N865, N865sim, N865simhis data with error bars
fig, ax = plt.subplots()
for i in ['N865', 'N865sim', 'N865simhis']:
    plt.errorbar(arange(-2, 398), processed_data[i].mean(1)[:400], 
             xerr=0, yerr=processed_data[i].std(1)[:400], 
             alpha=0.8, linewidth=0.5, capsize=0, fmt=colors[i])
plt.xlim(-10, 398)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig('20150526errorbars1.pdf')
plt.close()

# plot N865, N725, and N725sim (with N725sim only up to 200 sec)
fig, ax = plt.subplots()
i = 'N725'
plt.errorbar(arange(-2, 398), processed_data[i].mean(axis=1, skipna=True)[:400],
             xerr=0, yerr=processed_data[i].std(1)[:400], 
             alpha=0.8, linewidth=0.5, capsize=0, fmt=colors[i])
i = 'N865'
plt.errorbar(arange(-2, 398), processed_data[i].mean(axis=1, skipna=True)[:400],
             xerr=0, yerr=processed_data[i].std(1)[:400], 
             alpha=0.8, linewidth=0.5, capsize=0, fmt=colors[i])
i = 'N725sim'
plt.errorbar(arange(-2, 198), processed_data[i].mean(axis=1, skipna=True)[:200],
             xerr=0, yerr=processed_data[i].std(1)[:200], 
             alpha=0.8, linewidth=0.5, capsize=0, fmt=colors[i])
plt.xlim(-10, 398)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig('20150526errorbars2.pdf')
plt.close()

# curve fitting 
from scipy.optimize import curve_fit

constructs = [['N865', 750], ['N865sim', 400], ['N865simhis', 400],
              ['N725', 400], ['N725sim', 200]]
fig, ax = plt.subplots()

print 'construct' + '\t' + 'y0' + '\t' + 't1/2' + '\t' + 'mobile%' +'\tn'
for i in constructs:
    Y = np.array(processed_data[i[0]].mean(axis=1, skipna=True)[2:i[1]])
    x = np.arange(len(Y))
    
    frap = lambda x, plateau, t_half, decay: \
           plateau - (plateau - Y[0]) * 2.0 ** (- x / t_half) + decay * x
    
    # regression to calculate the half recovery time and mobile fraction       
    p, q = curve_fit(frap, x, Y, p0=[0.6, 30, 0])
    
    plt.plot(x, frap(x, *p), label=i[0], color=colors[i[0]], linewidth=2)
    print i[0], '\t', np.round(Y[0], 3), '\t', np.round(p[1], 3), '\t', \
          np.round((p[0] - Y[0])/(1 - Y[0]), 3), '\t', len(raw_data[i[0]])

plt.xlim(-10, 400); plt.ylim(0, 1.0)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig('20150526curve_fit.pdf')

'''
construct	y0		t1/2	mobile%	n
N865 		0.082 	76.847 	0.315 	18
N865sim 	0.093 	40.635 	0.426 	23
N865simhis 	0.113 	16.912 	0.673 	30
N725 		0.093 	17.555 	0.601 	22
N725sim 	0.163 	3.121 	0.69 	28

'''

# curve fitting for each sample 20150911
base_path = '/Users/sun_huaiyu/Documents/Projects/RING_SUMO/Raw_Data/'

raw_data = { 
'N725': glob(base_path + '20150518FRAP_ArkNsim_data/*N725_A4*.txt')
      + glob(base_path + '20150407U2OS_EGFP-ArkN_FRAP_summary/*N725*.txt'), 
'N725sim': glob(base_path + '20150518FRAP_ArkNsim_data/*N725sim_A5*.txt'),
'N865': glob(base_path + '20150518FRAP_ArkNsim_data/*N865_A1*.txt') 
      + glob(base_path + '20150407U2OS_EGFP-ArkN_FRAP_summary/*N865*.txt'),
'N865sim': glob(base_path + '20150518FRAP_ArkNsim_data/*N865sim_A2*.txt'),
'N865simhis': glob(base_path + '20150518FRAP_ArkNsim_data/*N865simhis_A3*.txt') }

c = ['N725', 'N725sim', 'N865', 'N865sim', 'N865simhis']

all_data = DataFrame([(i, j) for i in c for j in raw_data[i]],
                     columns=['construct', 'file_name'])

from scipy.optimize import curve_fit

def fit(f):
    df = pd.read_csv(f, sep='\t').dropna(0) 
    start = df[df['time_corrected'] == 0.0].index[0]
    end = len(df) - 2
    x = np.array(df['time_corrected'][start:end])
    y = np.array(df['region1_%'][start:end])
    Y0 = y[0]
    
    frap_model = lambda x, plateau, t_half, decay: \
               plateau - (plateau - Y0) * (2.0 ** (- x / t_half)) + decay * x
    
    p, q = curve_fit(frap_model, x, y, p0=[0.6, 30, 0])
    
    return (p[1], (p[0]- Y0) / (1 - Y0), p[2])

result = pd.merge(all_data, DataFrame([fit(i) for i in all_data.file_name], 
                columns=['t_half', 'mobile_frac', 'decay']), 
                left_index=True, right_index=True)

result.boxplot(by='construct', column='t_half', grid=False)
result.boxplot(by='construct', column='mobile_frac', grid=0)

from seaborn import stripplot
stripplot(x='construct', y='t_half', data=result, jitter=1)
