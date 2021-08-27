# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 16:14:04 2021

@author: JOELHERRERAVAZQUEZ
"""


# from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import axes3d
# import matplotlib.pyplot as plt

from scipy.stats import uniform
import random
import numpy as np


class SourceRnd:
    def __init__(self):
        
        self.fun=0
        self.field=89.9
        self.type=0
        self.dim=10
        self.num=100
        

    def rays(self):
        def ff(x):
            return x/x
        
        f=self.fun
        if f==0:
            f=ff
         
        r_source=self.dim
        lim = np.deg2rad(self.field)
    
        A = np.linspace(lim/self.num, lim-(lim/self.num), self.num)
        prob = f(A)*(2*np.pi*np.sin(A))
        prob=prob-np.min(prob)
        prob=prob/np.max(prob)
        
        Rand_num = random.choices(A, prob, k=self.num)
        Delta = (lim/ self.num)*(uniform.rvs(size=self.num, loc = 0, scale = 2)-1)    
        
        Theta=Rand_num+Delta
        Phi=uniform.rvs(size=self.num, loc = 0, scale = 2.0*np.pi)
        
        L = np.sin(Theta)*np.cos(Phi)
        M = np.sin(Theta)*np.sin(Phi)
        N = np.cos(Theta)
        
        
        if self.type==0: # Circle
            AR = np.linspace(r_source/self.num, r_source-(r_source/self.num), self.num)
            prob_R = 2*np.pi*AR
            lt0=np.argwhere(prob_R<0)
            
            prob_R = prob_R-np.min(prob_R)
            
            
            
            prob_R = prob_R/np.max(prob_R)
            
            
            Rand_num_R = random.choices(AR, prob_R, k=self.num)
            Delta_R = (lim/ self.num)*(uniform.rvs(size=self.num, loc = 0, scale = 2)-1)    
            R = Rand_num_R+Delta_R
                
            Phi_R=uniform.rvs(size=self.num, loc = 0, scale = 2.0*np.pi)
            X = R*np.cos(Phi_R)
            Y = R*np.sin(Phi_R)
            
        
        if self.type==1: #Square
            X = r_source*(uniform.rvs(size=self.num, loc = 0, scale = 2)-1)
            Y = r_source*(uniform.rvs(size=self.num, loc = 0, scale = 2)-1)
        
        Z=np.zeros_like(X)
        
        return L,M,N,X,Y,Z




# phi = np.linspace(0, np.pi, 20)
# theta = np.linspace(0, 2 * np.pi, 40)
# x = np.outer(np.sin(theta), np.cos(phi))
# y = np.outer(np.sin(theta), np.sin(phi))
# z = np.outer(np.cos(theta), np.ones_like(phi))

# Sol=SourceRand()

# def f(x):
#     res=np.cos(x*2*10)
#     return res

# Sol.fun=f
# Sol.dim=1
# Sol.num=10000
# L, M, N, X, Y, Z = Sol.rays()

# fig, ax = plt.subplots(1, 1, subplot_kw={'projection':'3d', 'aspect':'auto'})

# ax.plot_wireframe(x, y, z, color='k', rstride=1, cstride=1)

# ax.scatter(L, M, N, s=0.1, c='r', zorder=10)
# plt.plot(X, Y, ',', color='black');
        