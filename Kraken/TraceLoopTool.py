#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 14:05:41 2021

@author: joelherreravazquez
"""


def TraceLoop(x,y,z,L,M,N,W,Container):
    """TraceLoop.

    :param x:
    :param y:
    :param z:
    :param L:
    :param M:
    :param N:
    :param W:
    :param Container:
    """
    System=Container.SYSTEM
    for i in range(0,len(x)):
        pSource_0 = [x[i], y[i], z[i]]
        dCos=[L[i], M[i], N[i]]
        System.Trace(pSource_0, dCos, W)
        Container.push()
    return 0

def NsTraceLoop(x,y,z,L,M,N,W,Container):
    """NsTraceLoop.

    :param x:
    :param y:
    :param z:
    :param L:
    :param M:
    :param N:
    :param W:
    :param Container:
    """
    System=Container.SYSTEM
    for i in range(0,len(x)):
        pSource_0 = [x[i], y[i], z[i]]
        dCos=[L[i], M[i], N[i]]
        System.Trace(pSource_0, dCos, W)
        Container.push()
    return 0
