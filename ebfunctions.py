# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 22:24:50 2016

@author: Yousef L, Katie O'Neill, Jeffrey Ding, Joe Hardewicke, Jon Swift
"""

import math
import numpy as np



def markovChain(data = np.random.normal(5,s,100), plotDataError = False, stepSize = .01,nSteps = 10000, C0 = 5):
    if plotDataError:
        plt.errorbar(np.arrange(len(y)),data,yerr= np.ones(len(x)), font ='o', color = 'black')
    
    cvec = []
    lnlike = []
    
    q = 1/(Math.sqrt(2*Math.pi)*s)*math.e**(-())
    
    for i in range(nSteps):
        L0 = lnLike(C0,data)
        cnew = C0 + random
        Lnew = lnLike(cnew,data)
        if Lnew < L0:
            #acceptance probability: metropolis hasting
            
    
def lnLike(x,d,sigma):
    #x is a constant, d is an array, sigma is an array; returns p(D|m)
    return np.sum(-(x-d)**2/(2*sigma**2))


#r1 is radius of primary (bigger) star
#r2 is radius of secondary (smaller) star
#b is impact parameter (projected distance b/w their centers)
#area of chord tangent to bottom of smaller star
#A = (r1^2/2)(theta-sintheta)
# theta = angle of circular segment


def getImpactParamTra(a, i, r1, e, argofperi):
    i = math.radians(i)
    return (a*math.cos(i)/r1)*((1-e**2)/(1+e*math.sin(argofperi)))
    

def getImpactParamOcc(a, i, r1, e, argofperi):
    i = math.radians(i)
    return (a*math.cos(i)/r1)*((1-e**2)/(1-e*math.sin(argofperi))) 


def primary_total_area(r1,r2,b):
    #x = chord length/2
    adjacentSide = (b * r1) - r2
    theta = (2 * math.acos(adjacentSide / r1))
    total_area = ((r1 ** 2) / 2) * (theta - math.sin(theta))
    return total_area
    
def fractional_area(r1, r2, b, verbose=False):
    adjacentSide = (b * r1) - r2
    theta = (2 * math.acos(adjacentSide / r1))
    total_area = ((r1 ** 2) / 2) * (theta - math.sin(theta))
    frac_area = 100. * (total_area / (r1 ** 2 * math.pi))
    if verbose:
        print("Fractional area covered = %.2f percent" % (frac_area))
    return frac_area
