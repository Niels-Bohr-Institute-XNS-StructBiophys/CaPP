##       TESTING                ##

## Libraries
import re
import os
import sys
import math
import linecache
import re
import platform
import ntpath

try:
    import wx
except:
    print("")
    print("******************************************************************")
    print(" CaPP failed to import wxPython - is it correctly installed...?   ")
    print("******************************************************************")
    print("")
    sys.exit(1)

import wx.lib.scrolledpanel

programpath = os.getcwd()

## Import plotting libraries
import matplotlib
matplotlib.interactive(True)
matplotlib.use('WXAgg')
import pylab
import numpy as np

## import fitting libraries
from scipy.optimize import curve_fit

# testing
print("TESTING:")

a = np.array([1., 2., 3., 4., 5.])

b = np.array([4., 2., 3., 4., 5.])

c = np.sin(a)/b

print(c)


#a = np.array([1., 2., 3., 4., 5.])
#print a
#da = 0.1
# print da
#sigma = np.array([-3,-2,-1,0,1,2,3])
# print sigma
#a_MAT = np.matrix([a,a,a,a,a,a,a])
#for i in range(0,len(sigma)):
#     a_MAT[i] = sigma[i] * da + a
#print a_MAT
#print a_MAT[0]
#print a_MAT.item(0,1)
#print a_MAT.item((0,1))
#print a_MAT.item(0)
#print np.array(a_MAT)[0]
#print np.array(a_MAT)[0][1]

# print len(sigma)


#pcov = [[  4.41083708e-06  -7.19939497e-09  -1.55289185e-05][ -7.19939497e-09   4.94555459e-09  -1.15417531e-07][ -1.55289185e-05  -1.15417531e-07   1.00461684e-04]]


# q = np.genfromtxt("/Users/anlarsen/Seafile/GitHub/SASDAG2.dat", skiprows=2,usecols=[0], unpack=True) #import q vector
# if isinstance(list, np.genfromtxt("/Users/anlarsen/Seafile/GitHub/SASDAG2.dat", skiprows=2,usecols=[3], unpack=True)):
#     print "yes"
# else:
#     print "no"







# dq = q*0.1
# sigma = np.array([-3,-2,-1,0,1,2,3])
# q_long = np.array([])
# for i in range(0,3):
#     x = sigma * dq[i] + q[i]
#     if i == 0:
#         q_long = x
#     else:
        #print "i:"
        #print i
        #print "q[i]:"
        #print q[i]
        #print "q_long:"
        #print q_long
        #print "x[i]:"
        #print x
        #q_long = np.concatenate((q_long,x),axis=1)
        #print "q_long:"
        #print q_long


# PqRES   = [0]*2*2
# PqMAT   = np.matrix([PqRES , PqRES])
# print "PqRES:"
# print PqRES
# print "PqMAT:"
# print PqMAT
# print "PqMAT[1]:"
# print PqMAT[1]
# print "np.array(PqMAT)[1]:"
# print np.array(PqMAT)[1]

# for i in range(-3,3):
#     print "i = %d" % i 