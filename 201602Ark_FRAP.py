# Grant Figures 2016-02-26

from glob import glob
import pandas as pd; from pandas import Series, DataFrame
import matplotlib.pylab as plt
from pylab import *

base_path = '/Users/sun_huaiyu/Documents/Projects/RING_SUMO/Raw_Data/'

raw_data = { 
  'N725': glob(base_path + '20150518FRAP_ArkNsim_data/*N725_A4*.txt') + \
        glob(base_path + '20150407U2OS_EGFP-ArkN_FRAP_summary/*N725*.txt'),
  'N725sim': glob(base_path + '20150518FRAP_ArkNsim_data/*N725sim_A5*.txt'),
  'N725dHis': glob(base_path + '20160218FRAP_EGFP-ArkN_delHis/*N725dHis*.txt'),
  'N865': glob(base_path + '20150518FRAP_ArkNsim_data/*N865_A1*.txt') + \
          glob(base_path + '20150407U2OS_EGFP-ArkN_FRAP_summary/*N865*.txt'),
  'N865sim': glob(base_path + '20150518FRAP_ArkNsim_data/*N865sim_A2*.txt'),
  'N865dHis': glob(base_path + '20160218FRAP_EGFP-ArkN_delHis/*N865dHis*.txt'),
  'N865simdHis': glob(base_path + '20150518FRAP_ArkNsim_data/*N865simhis_A3*.txt')
}


# make sure the zero time point is at the minimum value
processed_data = dict()
for construct, files in raw_data.items():
    container = dict()
    for f in files:
        a = np.array(pd.read_csv(f, sep='\t')['region1_%'])
        container[f] = Series(a[a.argmin() - 2:])
    processed_data[construct] = DataFrame(container)

# color scheme for plotting
colors = {'N725': 'b', 'N725dHis': 'c', 'N725sim': 'r',
          'N865': 'y', 'N865dHis': 'c', 'N865sim': 'r', 'N865simdHis': 'g'}   

# plot all N865 data with error bars
fig, ax = plt.subplots()
for i in ['N865', 'N865dHis', 'N865sim', 'N865simdHis']:
    plt.errorbar(arange(-2, 198), processed_data[i].mean(1)[:200], 
             xerr=0, yerr=processed_data[i].std(1)[:200], 
             alpha=0.8, linewidth=0.5, capsize=0, fmt=colors[i])
plt.xlim(-10, 198)
plt.ylim(0, 1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig('/Users/sun_huaiyu/Documents/Projects/Grant/PolyHis/PolyHis_Figures/'
             + '20160226errorbars1.pdf')
plt.close()



# plot N865, N725, N725dHis and N725sim (up to 200 sec)
fig, ax = plt.subplots()
for i in ['N725', 'N725sim', 'N725dHis']:
    plt.errorbar(arange(-2, 198), processed_data[i].mean(axis=1, skipna=True)[:200],
             xerr=0, yerr=processed_data[i].std(1)[:200], 
             alpha=0.8, linewidth=0.5, capsize=0, fmt=colors[i])
plt.xlim(-10, 198)
plt.ylim(0, 1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig('/Users/sun_huaiyu/Documents/Projects/Grant/PolyHis/PolyHis_Figures/' +
            '20160226errorbars2.pdf')
plt.close()




