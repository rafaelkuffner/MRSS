#!/usr/bin/env python3

import sys, os
from numpy import *
from matplotlib.pyplot import *

datapath = sys.argv[1]
print('loading', datapath, '...')
data = loadtxt(datapath, delimiter=',')
print('loaded', data.shape[0], 'vertices')

area = data[:,0]
meancurv = data[:,1]
doncurv = data[:,2]

def plot_curv(curv, bins, curvlbl):
	figure()
	hist(curv, weights=area, bins=bins)
	xlabel(curvlbl)
	ylabel('Area')
	title('Vertex Area by {}'.format(curvlbl))
	figure()
	plot(sqrt(area), curv, '.')
	xlabel('√ vertex area')
	ylabel(curvlbl)
	title('{} vs Vertex Area'.format(curvlbl))
# }

plot_curv(abs(meancurv), 256, 'Absolute mean curvature')
plot_curv(doncurv, linspace(0, 1, 257), 'DoN curvature')

show(block=True)
