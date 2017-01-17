#!/home/birac/anaconda2/bin/python
#
# python script to plot any Gromacsdata file
# 

### Version control

import scipy as sc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit
from scipy import optimize
from scipy import signal
import os
import re
import sys
import csv
import shutil
import pandas as pd
import subprocess
import os.path
from scipy.interpolate import spline

def loaddata ( str ):
        print "%s" % y
	return

filename = sys.argv[1]
jin = sys.argv[2]
j = int(jin)
#legend = sys.argv[3]
xlabel = "time ($\mu$s)"
ylabel = "residue to CL distance (nm)"
xdata, ydata = np.loadtxt(fname='%s' % filename, delimiter=' ', usecols=(0,1), unpack=True)	

params = {'legend.fontsize': 'large',
	'axes.labelsize': 'x-large',
	'xtick.labelsize': 'x-large',
	'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['font.sans-serif'] = "cmss10"
plt.rcParams['axes.unicode_minus']=False

plt.plot(xdata, ydata, color='gray', marker='.', markersize=0)
smooth = sc.signal.savgol_filter(ydata, 11, 1, deriv=0, delta=1, axis=-1, mode='interp', cval=0.0)
colors = cm.rainbow(np.linspace(0, 1, 16))
plt.plot(xdata, smooth, color=colors[j], linewidth=0.5)

plt.title('%s' % filename )
plt.ylim([0,5])
plt.xlabel('%s' % xlabel, fontname="cmss10", fontsize=25 )
plt.ylabel('%s' % ylabel, fontname="cmss10", fontsize=25 )
#plt.legend()
plt.savefig('%s.png' % filename, bbox_inches='tight')
