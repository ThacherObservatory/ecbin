#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:59:12 2016

@author: Ding
"""

from __future__ import print_function
import numpy as np
import emcee
import occultquad as oq
import scipy as sp
import matplotlib.pyplot as plt
import robust as rb
import sys
import math
from scipy.ndimage.filters import uniform_filter
from mpl_toolkits.mplot3d import Axes3D
import time
import constants as c
import os
import pdb
from PyAstronomy import pyasl

def run():
    read_data()
    ndim, nwalkers = 3, 100
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    import corner
    fig = corner.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],truths=[m_true, b_true, np.log(f_true)])
    fig.savefig("triangle.png")
    
def read_data():
    #read files
    time1, rv1, rv1error = np.loadtxt("Data_large.txt", unpack=True)
    time2, rv2, rv2error = np.loadtxt("Data_small.txt", unpack=True)

    data_dict = {"TimeOne" : time1, "rv1" : rv1, "rv1error": rv1error, "TimeTwo": time2, "rv2": rv2, "rv2error": rv2error}
    return data_dict

def lnlike(params, data_dictionary):
    rvModel = get_model_rv(t, params)
    return -np.sum(((rvModel - rvActual) ** 2) / (2 * (error ** 2)))
    

def get_model_rv(t,params=False):
#    if len(params) == 9:
#    p = np.double(params)
#    else:
#    p = p0
    p = params
#        period = p[0] * 24.0 * 3600.0
    per    = period * 24.0 * 3600.0
    M1     = p[0] * c.Msun
    M2     = p[1] * c.Msun
    t0     = p[2]
    inc    = p[3]
    ecosw  = p[4]
    esinw  = p[5]
    offset = p[6]
#        drift  = p[8]
    drift  = 0.0

    ecc  = np.sqrt(ecosw**2+esinw**2)

    omega    = ((np.degrees(np.arctan2(esinw,ecosw))) + 360.0) % 360.0

#    modelt = np.linspace(0,period/(24.0*3600.0),1000,endpoint=True)


    K1  = ((2.0 * np.pi * c.G * M2**3.0 * np.sin(np.radians(inc))**3.0 / \
                (per * (1.0 - ecc**2)**(1.5) * (M1 + M2)**2))**(1.0/3.0))/ 1e5

    n = 2.0 * np.pi / period # in days

    # Ephemeris given w.r.t. eclipse, not periastron passage!!!
    M = n * (t - t0) % (2.0 * np.pi)

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