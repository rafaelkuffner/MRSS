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

figure()
plot(sqrt(area), abs(meancurv), '.')
title('Mean Curvature vs Vertex Area')
xlabel('√ vertex area')
ylabel('Absolute mean curvature')

figure()
plot(sqrt(area), doncurv, '.')
title('DoN curvature vs Vertex Area')
xlabel('√ vertex area')
ylabel('DoN curvature')

show(block=True)
