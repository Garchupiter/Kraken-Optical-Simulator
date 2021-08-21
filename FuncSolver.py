#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 22:12:12 2021

@author: joelherreravazquez
"""

from scipy.optimize import fsolve
import numpy as np

def func(x):
    y=x[0]
    z=x[1]

    v1 = y * np.cos(z) - 4
    v2 = z * y - z - 5
    return [v1, v2]

root = fsolve(func, [1, 1])
print(root)

aa=np.isclose(func(root), [0.0, 0.0])  # func(root) should be almost 0.0.
print(aa)
