#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulation of a Fresnel lens using KrakenOS

This script loads the profile of a Fresnel lens from a text file,
builds an optical system in KrakenOS, and traces rays through the system to
visualize the lens behavior.
"""

import pkg_resources
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# -----------------------------------------------------------------------------
# Check that KrakenOS is available
# -----------------------------------------------------------------------------

required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("KrakenOS is not installed. Adding relative path for local environment.")
    sys.path.append("../..")

import KrakenOS as Kos

# -----------------------------------------------------------------------------
# Object surface (input of the system)
# -----------------------------------------------------------------------------

P_Obj = Kos.surf()
P_Obj.Rc = 0.0  # Radius of curvature (0 = flat)
P_Obj.Thickness = 150  # Thickness to the next surface
P_Obj.Glass = "AIR"
P_Obj.Diameter = 50.0  # Effective diameter

# -----------------------------------------------------------------------------
# Class to process the Fresnel lens profile
# -----------------------------------------------------------------------------

class FresnelPrepare:
    def __init__(self, file):
        """Initialize profile from a CSV file with x, y columns."""
        data = np.loadtxt(file, delimiter=',')
        x = data[:, 0]
        y = data[:, 1]

        # Discrete derivatives (slopes between points)
        ny = np.roll(y, shift=1)
        nny = ny - y
        nx = np.roll(x, shift=1)
        nnx = nx - x
        m = nny / nnx

        # Thresholding to detect sharp slope transitions (zone edges)
        threshold = 0
        m[m < threshold] = 0
        m[m > threshold] = 1

        # Detect decreasing transitions (zone edges)
        vi = 0
        ARG = []
        for i in range(1, len(m)):
            v = m[i]
            if v < vi:
                vi = v
                ARG.append(i)
            vi = v
        ARG = np.asarray(ARG)

        # Adjust indices to center segments
        ARG1 = np.roll(ARG, shift=1)
        ARG1[0] = 0
        ARG = ARG - ((ARG - ARG1) / 2.0)
        ARG = ARG.astype(int)

        # Use points 10 positions before and after to define facets
        self.X0 = x[ARG - 10]
        self.Y0 = y[ARG - 10]
        self.X1 = x[ARG + 10]
        self.Y1 = y[ARG + 10]

        # Compute slope and intercept for each facet
        self.M = (self.Y1 - self.Y0) / (self.X1 - self.X0)
        self.b = self.Y0 - (self.M * self.X0)

    def find_closest_indices(self, x, r):
        """Find the index in r closest to each value in x."""
        closest_indices = np.empty(len(x), dtype=int)
        for i, value_x in enumerate(x):
            differences = np.abs(r - value_x)
            closest_indices[i] = np.argmin(differences)
        return closest_indices

    def calculate(self, x, y, E):
        """Compute surface height z at coordinates x, y."""
        r = np.sqrt(x*x + y*y)
        ARG_C = self.find_closest_indices(r, self.X1)
        Mp = self.M[ARG_C]
        bp = self.b[ARG_C]
        z = Mp * r + bp
        return z

# -----------------------------------------------------------------------------
# Define the Fresnel lens from the profile file
# -----------------------------------------------------------------------------

file = "R1064_F1800.txt"  # Lens profile file
fresnel = FresnelPrepare(file)
E = []  # Placeholder for additional data if needed

# Structured surface (zones)
L1a = Kos.surf()
L1a.Rc = 0.0
L1a.Thickness = 5.0
L1a.Glass = "PMMA"
L1a.Diameter = 2128.0
L1a.ExtraData = [fresnel.calculate, E]  # Defined with external function
L1a.Res = 1
L1a.Name = 'zone side'
L1a.Nm_Pos = (-500, 200)
# DerPres controls the numerical precision of derivative calculations.
# It should only be modified in cases where the surface includes extremely fine structures
# such as those present in Fresnel optical elements.
L1a.DerPres = 0.00001

# Flat back surface
L1b = Kos.surf()
L1b.Thickness = 1795.0
L1b.Glass = "AIR"
L1b.Diameter = 2128.0
L1b.Name = 'flat side'
L1b.Nm_Pos = (500, 200)

# Image surface
P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 8000.0
P_Ima.Name = "Image plane"

# -----------------------------------------------------------------------------
# Build the complete optical system
# -----------------------------------------------------------------------------

A = [P_Obj, L1a, L1b, P_Ima]
Config_1 = Kos.Setup()
Lens = Kos.system(A, Config_1)
Rays = Kos.raykeeper(Lens)

# -----------------------------------------------------------------------------
# Generate and trace rays
# -----------------------------------------------------------------------------

SemiDiameter = 2128.0 / 2.0
num = 300
I = np.linspace(-SemiDiameter, SemiDiameter, num=num)  # Y-axis positions
Wav = 0.55  # Wavelength (microns)

for i in I:
    pSource = [0.0, i, 0.0]  # Ray origin
    dCos = [0.0, 0.0, 1.0]   # Direction (Z axis)
    Lens.Trace(pSource, dCos, Wav)
    Rays.push()

# -----------------------------------------------------------------------------
# Visualization of system and rays
# -----------------------------------------------------------------------------

Kos.display3d(Lens, Rays, 0)  # 3D view
Kos.display2d(Lens, Rays, 0)  # 2D view
