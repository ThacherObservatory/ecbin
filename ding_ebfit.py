#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:59:12 2016

@author: Ding
"""

from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import emcee
import constants as c
#import occultquad as oq
#import scipy as sp
#import matplotlib.pyplot as plt
#import robust as rb
#import sys
#import math
#from scipy.ndimage.filters import uniform_filter
#from mpl_toolkits.mplot3d import Axes3D
#import time
#import constants as c
#import os
import pdb
from PyAstronomy import pyasl
from scipy import stats

def run():
    
    ndim= 7
    nwalkers = 250
    #check dimensions
    m1 = np.random.normal(0.6, 0.01, nwalkers)
    m2 = np.random.normal(0.4, 0.01, nwalkers)
    t0 = np.random.normal(0, 100, nwalkers)
    #inc = np.random.normal(90, ,2 nwalkers)
    ecosw = np.random.normal(0, 0.1, nwalkers)
    esinw = np.random.normal(0, 0.1, nwalkers)
    offset = np.random.normal(0, 1, nwalkers)
    period = np.random.normal(4.2, .1, nwalkers)
    
    
    
    p0 = np.array([m1, m2, t0,  ecosw, esinw, offset, period]).T

    #p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim)) #random locations for 250 walkers in 8 dimension
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads = 4) #instantiate emcee 
    #pdb.set_trace()
    
    pos, prob, state = sampler.run_mcmc(p0, 100) #burn-in
    sampler.reset()


    sampler.run_mcmc(pos, 1000)
    #samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    
    for i in range(ndim):
        plt.figure()
        #y, x, _ = plt.hist(sampler.flatchain[:, i])
        y, x, _ = plt.hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
        #max = x.max()  # Find the maximum y value
        #max_x = x[sampler.flatchain[:, i].index(max)]  # Find the x value corresponding to the maximum y value
        #max_x = stats.mode(sampler.flatchain[:, i])
        print(x)
        
        print(y)
        #print(max)
        plt.title("Dimension {0:d}".format(i))
        plt.show()
        
    return sampler
    
def read_data():
    #read files
    time1, rv1, rv1error = np.loadtxt("Data_large.txt", unpack=True)
    time2, rv2, rv2error = np.loadtxt("Data_small.txt", unpack=True)

    data_dict = {"TimeOne" : time1, "rv1" : rv1, "rv1error": rv1error, "TimeTwo": time2, "rv2": rv2, "rv2error": rv2error}
    return data_dict

def lnlike(params, data_dictionary):

   
    rvModel = get_model_rv(data_dictionary["TimeOne"], params=params) 
    
    returnval =  -np.sum(((rvModel - data_dictionary["rv1"]) ** 2) / (2 * (data_dictionary["rv1error"] ** 2)))
 
    return returnval
    
def lnprior(params):
    p = params
    period = p[6] * 24.0 * 3600.0
    per    = period * 24.0 * 3600.0
    M1     = p[0] * c.Msun
    M2     = p[1] * c.Msun
    t0     = p[2]
   # inc    = p[3]
    ecosw  = p[3]
    esinw  = p[4]
    offset = p[5]
    ecc = np.sqrt(ecosw**2+esinw**2)
    if ecc > 1:
        return -np.inf
        
    if per < 0:
        return -np.inf
    return 1
    
def lnprob(params):
    data_dictionary = read_data()
    lp = lnprior(params)
    if not np.isfinite(lp):
        return -np.inf

    hold = lp + lnlike(params, data_dictionary)
    
    return hold 
    
    
    

def get_model_rv(t,params=False):

    
#    if len(params) == 9:
#    p = np.double(params)
#    else:
#    p = p0
    #pdb.set_trace()
    p = params
    period = p[6] * 24.0 * 3600.0
    per    = period * 24.0 * 3600.0
    M1     = p[0] * c.Msun
    M2     = p[1] * c.Msun
    t0     = p[2]
    #inc    = p[3]
    ecosw  = p[3]
    esinw  = p[4]
    offset = p[5]
#        drift  = p[8]
    drift  = 0.0
    inc = 90
    ecc  = np.sqrt(ecosw**2+esinw**2)
 
    omega    = ((np.degrees(np.arctan2(esinw,ecosw))) + 360.0) % 360.0

#    modelt = np.linspace(0,period/(24.0*3600.0),1000,endpoint=True)

    num = float(1)/3
    holder = (2.0 * np.pi * c.G/per)
    KpartOne = np.power(holder, (num))
    KpartTwo = (M2 * np.sin(np.radians(inc)) / ((1.0 - ecc**2) ** 0.5 * (M1+M2) ** (2.0/3.0)))
    K1 = KpartOne * KpartTwo
    #K1  = ((2.0 * np.pi * c.G * M2**3.0 * np.sin(np.radians(inc))**3.0 / \
                #(per * (1.0 - ecc**2)**(1.5) * (M1 + M2)**2))**(1.0/3.0))/ 1e5

    n = 2.0 * np.pi / period # in days

    # Ephemeris given w.r.t. eclipse, not periastron passage!!!
    M = n * (t - t0) % (2.0 * np.pi)

    if np.isnan(K1):
        pdb.set_trace()
    
# Instantiate the solver
    ks = pyasl.MarkleyKESolver()
    if len(M) > 0:
        E = np.array([ks.getE(m, ecc) for m in M])
    else:
        E = ks.getE(M, ecc)

    f = 2.0 * np.arctan(np.sqrt( (1.0 + ecc) / (1.0 - ecc) ) * np.tan(E/2.0))


    rv1_model = K1*(np.cos(np.radians(omega) + f) + ecc*np.cos(np.radians(omega))) + \
        offset + drift*(t - t0)

    rv2_model = -1.0*K1*(np.cos(np.radians(omega) + f) + \
                             ecc*np.cos(np.radians(omega)))*M1/M2 + \
                             offset + drift*(t - t0)

    return rv1_model,rv2_model