# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 22:24:50 2016

@author: Yousef L, Katie O'Neill, Jeffrey Ding, Joe Hardewicke, Jon Swift
"""

import math

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


def primary_frac_area(r1,r2,b):
    #x = chord length/2
    x = math.sqrt(abs(r1**2 - (b-r2)**2))
    coveredArea = ((r1**2)/2)*(2*math.asin(math.radians(x/r1))-math.sin(math.radians(2*math.asin(math.radians(x/r1)))))
    fracAreaOfCoveredPrimary = coveredArea / (math.pi*r1**2)
    return fracAreaOfCoveredPrimary
   
   

   
   
    
