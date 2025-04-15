#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 07:54:16 2021

@author: joelherreravazquez & Antonio A.
"""

import numpy as np
import matplotlib.pyplot as plt
from .MathShapesClass import *

from scipy.ndimage import rotate


import numpy as np
import matplotlib.pyplot as plt

###############################################################################


import numpy as np
import matplotlib.pyplot as plt

def psf4mtf(COEF, Focal, Diameter, Wave, pixels=600, PupilSample=4):
    """Función que calcula la PSF a partir de los coeficientes de Zernike."""
    N = pixels  # Número de elementos en la matriz
    Q = PupilSample  # Muestreo de la pupila
    # D = Diameter / 1000.0  # Diámetro de la apertura en metros
    # FocalD = Focal / 1000.0  # Longitud focal en metros
    # wvl = Wave * 1e-6  # Longitud de onda en metros

    # Crear la cuadrícula de coordenadas
    TamImag = int(N)
    r = (TamImag / (Q * 2.0))
    center = int((TamImag / 2.0))
    xy = np.arange(0, TamImag)
    X, Y = np.meshgrid(xy, xy)
    x = (X - center) / r
    y = (Y - center) / r
    R = np.sqrt(x**2 + y**2)

    # Calcular el frente de onda usando los polinomios de Zernike
    W = Wavefront_Zernike_Phase(x, y, COEF)
    W[R > 1] = 0.0  # Máscara para limitar a la apertura

    # Crear la función pupila
    T = np.ones_like(W)
    T[R > 1] = 0.0  # Máscara circular

    # Campo complejo de la pupila
    U = T * np.exp(-1j * 2 * np.pi * W)

    # Cálculo de la PSF mediante la transformada de Fourier
    F0 = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(U))) / np.sum(T)
    I = np.abs(F0)**2  # PSF: Intensidad

    return I

def calculate_mtf(COEF, Focal, Diameter, w, pixels=600, PupilSample=4):
    """Calcula la MTF a partir de la PSF."""
    psf = psf4mtf(COEF, Focal, Diameter, w, pixels=600, PupilSample=4)
    psf = psf / np.sum(psf)  # Normalizar la PSF por su área total
    F0 = np.fft.fft2(psf)  # Transformada de Fourier de la PSF
    F0 = np.fft.fftshift(F0)  # Centrar la transformada de Fourier
    mtf = np.abs(F0)  # Módulo de la OTF
    mtf = mtf / np.max(mtf)  # Normalización para que MTF(0) = 1
    return mtf

def plot_mtf(mtf, Diameter, w, freq_limit = 1100):
    """Función para graficar las MTF tangencial y sagital."""
    N = mtf.shape[0]

    # Calcular la frecuencia máxima en cycles per mm
    freq_max = (Diameter / (w * 1e-3))  # Frecuencia máxima
    freq = np.linspace(0, freq_max, N // 2)  # Frecuencia espacial positiva

    # Curvas tangencial y sagital (direcciones X y Y de la MTF)
    mtf_tangential = mtf[N // 2, N // 2:]  # Corte horizontal (tangencial)
    mtf_sagittal = mtf[N // 2:, N // 2]  # Corte vertical (sagital)

    # Graficar ambas curvas en la misma figura
    plt.figure(figsize=(10, 6))
    plt.plot(freq/10, mtf_tangential[:len(freq)], label='Tangential MTF', color='blue')
    plt.plot(freq/10, mtf_sagittal[:len(freq)], label='Sagittal MTF', color='red')
    plt.xlabel('Spatial Frequency in cycles per mm')
    plt.ylabel('Modulus of the OTF')
    plt.title('Polychromatic Diffraction MTF')
    plt.grid(True)
    plt.legend()

    plt.xlim(0, freq_limit)
    plt.ylim(0, 1)
    plt.show()



#########################################################################333

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

    # Número de Elementos en el arreglo
    N = pixels

    # Variable de Muestreo
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

    # f=R<0.5
    # W[f]=1.0


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

    # I = W
    # Saturando la Imagen
    if plot !=0:
        plt.figure("PSF")
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




def PsfPlus(COEF, Focal, Diameter, Wave, pixels=600, PupilSample=4,  plot=0, sqr=0):

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

    # Número de Elementos en el arreglo
    N = pixels

    # Variable de Muestreo
    Q = PupilSample

    #Diámetro de la Apertura [m]
    D = Diameter / 1000.0

    #Diámetro de la Obstrucción[m]
    Do = 0.5
    #Número F
    FocalD = Focal / 1000.0
    #Longitud de Onda [m]
    wvl = Wave*1e-6 # Cambiando wave a metros

    TamImag = int(N)
    r = (TamImag / (Q * 2.0))

    center = int((TamImag / 2.0))
    xy = np.arange(0,TamImag)
    [X,Y] = np.meshgrid(xy, xy)
    x = ((X - center) / r)
    y = ((Y - center) / r)
    R = np.sqrt((x**2)+(y**2))

    W = Wavefront_Zernike_Phase(x, y, COEF)

    f = R>1
    W[f] = 0.0

    # # Telescope secundary mirror and spider
    # # f=R<0.3
    # # T[f]=0.0

    # # rx = np.abs(x)
    # # f=rx<0.001
    # # T[f]=0.0

    # # ry = np.abs(y)
    # # f=ry<0.001
    # # T[f]=0.0
    # # T = rotate(T, 45, reshape=False)


    T = np.copy(W)
    f = R <= 1
    T[f] = 1.0
    #Area del Círculo
    a = np.sum(np.sum(T))

    #Función de Pupila
    ##############################################################
    #Mapa del Frente de Onda con Máscara [wvl rms]
    # Wt=W*T

    #Campo complejo
    U = T * np.exp( -1j * 2 * np.pi * W )

    #PSF-Fraunhofer
    ##############################################################
    #TF de la distribución de Amplitud en Pupila
    F0 = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(U))) / a

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

    # I = W

    # plt.imshow(W, cmap= plt.cm.bone)
    # plt.colorbar()

    # Saturando la Imagen
    if plot !=0:
        plt.figure("PSF")
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

    # return I, I

    return I, vmn*2




