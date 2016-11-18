# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 22:24:50 2016

@author: Yousef L, Katie O'Neill, Jeffrey Ding, Joe Hardewicke, Jon Swift
"""

import math,pdb
import numpy as np
import matplotlib.pyplot as plt



def markovChain(c0 = 10, mean = 5, sigma = 1, nPts = 100,burnMax = 10000, 
                plotDataError = False, stepSize = .1,nSteps = 20000):
   # if plotDataError:
      #  plt.errorbar(np.arrange(len(y)),data,yerr= np.ones(len(x)), font ='o', color = 'black')
    
    data = np.random.normal(mean,sigma,nPts)
    l0 = lnLike(c0, data,sigma)
    cvec = np.array([c0])
    lnlike = np.array([l0])
    burnCounter = 0
    
    cfinal = []
    for i in range(nSteps):
        #if i % 2000 == 0:        
            #pdb.set_trace() 
        cprev = cvec[-1]
        lprev = lnLike(cprev,data,sigma)
       
        random = np.random.normal(0,stepSize)
        ccur = cprev + random
        lcur = lnLike(ccur,data,sigma)

        if lcur < lprev:
            lnratio = lcur - lprev
            ratio = np.exp(lnratio)
            prop = np.random.uniform(0,1)
            if ratio > prop:
                cvec = np.append(cvec,ccur)
                lnlike = np.append(lnlike,lcur)
                if burnCounter > burnMax:
                    cfinal = np.append(cfinal,ccur)
            else: 
                cvec = np.append(cvec,cprev)
                lnlike = np.append(lnlike,lprev)
                if burnCounter > burnMax:
                    cfinal = np.append(cfinal,cprev)
        else:
            cvec = np.append(cvec,ccur)
            lnlike = np.append(lnlike,lcur)            
            if burnCounter > burnMax:
                cfinal = np.append(cfinal,ccur)
            
        burnCounter = burnCounter + 1
            
    plt.ion()
    plt.figure(3)
    plt.clf()    
    plt.errorbar(np.arange(len(data)),data,yerr=np.ones(len(data)),fmt='o')
    
    plt.figure(4)
    plt.clf()
    plt.hist(cfinal,bins=100)
    
    return cfinal
        
            
    
def lnLike(x,data,sigma):
    #x is a constant, d is an array, sigma is an array; returns p(D|m)
    return np.sum(-1*(x-data)**2/(2*sigma**2))

#def like()


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
