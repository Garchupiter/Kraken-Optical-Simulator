#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 18:59:22 2021

@author: joelherreravazquez
"""
import numpy as np
import scipy

def R_RMS_delta(delta_Z,L,M,N,X,Y):

    X = ((L/N)*(delta_Z)) + X
    Y = ((M/N)*(delta_Z)) + Y
    cenX=np.mean(X)
    cenY=np.mean(Y)
    x1 = X - cenX
    y1 = Y - cenY
    R2 = (( x1 * x1 ) + ( y1 * y1 ))
    R_RMS = np.sqrt(np.mean( R2 ))
    return R_RMS

def RMS(X,Y,Z,L,M,N):
    cenX=np.mean(X)
    cenY=np.mean(Y)
    x1 = X - cenX
    y1 = Y - cenY
    R2 = (( x1 * x1 ) + ( y1 * y1 ))
    R_RMS = np.sqrt(np.mean( R2 ))
    return R_RMS

def BestFocus(X,Y,Z,L,M,N):
    delta_Z=0
    ZZ=L,M,N,X,Y
    v=scipy.optimize.fsolve(R_RMS_delta, delta_Z,  args=(ZZ))
    return v[0]
