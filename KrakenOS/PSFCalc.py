#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 07:54:16 2021

@author: joelherreravazquez & Antonio A.
"""

import numpy as np
import matplotlib.pyplot as plt
from .MathShapesClass import *



def psf(COEF, Focal, Diameter, Wave, pixels=600, PupilSample=4,  plot=0, sqr=0):

    ''' Generando los polinomios de Zernike (Nomenclarura Noll)con:
                   zernike_expand de Jherrera
    '''
    L=38 # Numero de terminos para la expanción polinomial


    ########################################################################

    '''Utilizando un solo valor de la pupila (x,y) con
                  Wavefront_Zernike_Phase
    '''
    Zern_pol, z_pow = zernike_expand(len(COEF))

    # ####################################################################

    """
            Generando un mapa de pupila aberrado
    """

    #Número de Elementos en la Matriz
    N=pixels
    #Variable de Muestreo
    Q=PupilSample
    #Diámetro de la Apertura [m]
    D=Diameter/1000.0
    #Diámetro de la Obstrucción[m]
    Do = 0.5
    #Número F
    FocalD=Focal/1000.0
    #Longitud de Onda [m]
    wvl=Wave*1e-6 # Cambiando wave a metros

    TamImag = int(N)
    r = (TamImag / (Q*2.0))

    center = int((TamImag / 2.0))
    xy=np.arange(0,TamImag)
    [X,Y]=np.meshgrid(xy, xy)
    x = ((X - center) / r)
    y = ((Y - center) / r)
    R=np.sqrt((x**2)+(y**2))

    W = Wavefront_Zernike_Phase(x, y, COEF)

    f=R>1
    W[f]=0.0

    T=np.copy(W)
    f=R<=1
    T[f]=1.0
    #Area del Círculo
    a=np.sum(np.sum(T))

    #Función de Pupila
    ##############################################################
    #Mapa del Frente de Onda con Máscara [wvl rms]
    # Wt=W*T
    #Campo complejo
    U=T*np.exp(-1j*2*np.pi*W)

    #PSF-Fraunhofer
    ##############################################################
    #TF de la distribución de Amplitud en Pupila
    F0=np.fft.fftshift(np.fft.fft2(np.fft.fftshift(U)))/a
    #Irradiancia en el Plano de Observación
    I=np.abs(F0)**2

    #Vector de Muestras con Cero en el Centro
    ##############################################################
    #Vector de Muestras
    v=np.arange(0,N)
    #Elemento central [bin]
    c=np.floor(N/2.)
    c=int(c)
    #Corriendo el cero al centro
    vx=v-c
    vy=c-v

    #Coordenadas del Plano de la Apertura
    ##############################################################
    #Largo del Plano
    L=Q*D
    ##############################################################
    #Coordenadas del Plano de Imagen
    #Tamaño del paso en frecuencia [m^-1]
    f=1/L
    # Vector de Coordenadas fx y fy
    fx=f*vx
    fy=f*vy

    # Vector de Coordenadas x y y
    u=fx*wvl*FocalD
    v=fy*wvl*FocalD

    umx=u[N-1]*1e6
    umn=u[0]*1e6
    vmx=v[N-1]*1e6
    vmn=v[0]*1e6

    I = np.rot90(I)
    # Saturando la Imagen

    if plot !=0:
        plt.figure(1)
        if sqr == 0:
            plt.imshow(I,extent=[umn,umx,vmx,vmn], cmap= plt.cm.bone)
            plt.colorbar()
            plt.ylabel('V[μm]')
            plt.xlabel('U[μm]')
            plt.title('Fraunhofer Prop - PSF')

        if sqr == 1:
            II=np.sqrt(I/np.max(I))*255
            plt.imshow(II,extent=[umn,umx,vmx,vmn], cmap= plt.cm.bone)
            plt.colorbar()
            plt.ylabel('V[μm]')
            plt.xlabel('U[μm]')
            plt.title('Fraunhofer Prop - PSF ( note: sqrt(I) )')

        if sqr == 2:
            III = np.log(I)
            III = III-np.min(III)
            III = 255*III/np.max(III)

            plt.imshow(III,extent=[umn,umx,vmx,vmn], cmap= plt.cm.bone)
            plt.colorbar()
            plt.ylabel('V[μm]')
            plt.xlabel('U[μm]')
            plt.title('Fraunhofer Prop - PSF ( note: sqrt(I) )')

    return I




