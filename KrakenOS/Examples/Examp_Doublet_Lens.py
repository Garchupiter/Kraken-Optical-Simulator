#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: Doublet Lens with Focus Optimization
This script demonstrates how to:
  1. Build a doublet lens system using KrakenOS.
  2. Trace rays for multiple wavelengths.
  3. Optimize the image plane focus using a Newton–Raphson method.
  4. Adjust the image plane position, re-trace the rays, and display an updated 2D plot.

Author: Joel Herrera V.
Date: 10/03/2025
"""

import numpy as np
import pkg_resources

# =============================================================================
# Check if KrakenOS is installed. If not, assume that the code is run from a 
# downloaded GitHub folder and add the relative path.
# =============================================================================
required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("KrakenOS is not installed. Using local GitHub folder.")
    import sys
    sys.path.append("../..")  # Adjust this path if necessary

import KrakenOS as Kos  # Using KrakenOS for optical simulation

# =============================================================================
# Helper functions for calculating the RMS spot size and its derivative.
# =============================================================================
def R_RMS(L, M, N, X, Y, delta_Z):
    """
    Calculate the RMS radius (spot size) of rays on the image plane after a 
    shift in the Z direction (delta_Z).
    """
    cenX = np.mean(X)
    cenY = np.mean(Y)
    nX = X - cenX
    nY = Y - cenY
    # Propagate the deviations in X and Y using the ray direction cosines and delta_Z.
    x1 = ((L / N) * delta_Z) + nX
    y1 = ((M / N) * delta_Z) + nY
    R2 = (x1 * x1) + (y1 * y1)
    return np.sqrt(np.mean(R2))

def DER_R_RMS(L, M, N, X, Y, delta_Z):
    """
    Compute the numerical derivative of R_RMS with respect to delta_Z using 
    a central difference method.
    """
    h = 0.001  # Small step for numerical differentiation
    f1 = R_RMS(L, M, N, X, Y, delta_Z + h)
    f2 = R_RMS(L, M, N, X, Y, delta_Z - h)
    return (f1 - f2) / (2.0 * h)

# =============================================================================
# Define the optical surfaces for the doublet lens system.
# The system is composed of:
#   - Object plane (flat, air)
#   - First lens surface (convex, BK7)
#   - Second lens surface (concave, F2)
#   - Air gap before the image plane
#   - Image plane (detector)
# =============================================================================
# Object Plane
P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 10
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# First Lens Surface (convex, BK7)
L1a = Kos.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

# Second Lens Surface (concave, F2)
L1b = Kos.surf()
L1b.Rc = -3.071608670000159E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# Air gap before image plane
L1c = Kos.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001
L1c.Glass = "AIR"
L1c.Diameter = 30

# Image Plane
P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 3.0
P_Ima.Name = "Image Plane"

# Build the system
# Note: The order of surfaces in list A determines the optical sequence.
A = [P_Obj, L1a, L1b, L1c, P_Ima]
config_1 = Kos.Setup()
Doublet = Kos.system(A, config_1)  # 'Doublet' is the system (formerly 'Doblete')

# =============================================================================
# Create ray containers for storing the traced rays.
# We'll use a global container for all wavelengths.
# =============================================================================
raysWavelength1 = Kos.raykeeper(Doublet)  # For wavelength 0.4
raysWavelength2 = Kos.raykeeper(Doublet)  # For wavelength 0.5
raysWavelength3 = Kos.raykeeper(Doublet)  # For wavelength 0.6
raysTotal = Kos.raykeeper(Doublet)        # Global container

# =============================================================================
# Generate rays across a circular aperture and trace them through the system.
# =============================================================================
grid_resolution = 10   # Grid resolution parameter (number of steps)
aperture_radius = 10.0  # Maximum radius of the circular aperture (mm)

for j in range(-grid_resolution, grid_resolution + 1):
    for i in range(-grid_resolution, grid_resolution + 1):
        # Compute the ray origin (x, y) in the object plane.
        x0 = (i / grid_resolution) * aperture_radius
        y0 = (j / grid_resolution) * aperture_radius
        if np.sqrt(x0**2 + y0**2) < aperture_radius:
            # For simplicity, rays are launched parallel to the optical axis (0° tilt).
            angle_deg = 0.0
            pSource = [x0, y0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(angle_deg)), np.cos(np.deg2rad(angle_deg))]
            
            # Trace ray at wavelength 0.4
            wavelength = 0.4
            Doublet.Trace(pSource, dCos, wavelength)
            raysWavelength1.push()
            raysTotal.push()
            
            # Trace ray at wavelength 0.5
            wavelength = 0.5
            Doublet.Trace(pSource, dCos, wavelength)
            raysWavelength2.push()
            raysTotal.push()
            
            # Trace ray at wavelength 0.6
            wavelength = 0.6
            Doublet.Trace(pSource, dCos, wavelength)
            raysWavelength3.push()
            raysTotal.push()

# =============================================================================
# Display the initial 2D plot of the ray paths.
# =============================================================================
Kos.display2d(Doublet, raysTotal, 0)

# =============================================================================
# Focus optimization using the Newton–Raphson method.
#
# Extract the ray data from the last surface (image plane) to compute the RMS
# spot size and its derivative.
# =============================================================================
X, Y, Z, L, M, N = raysTotal.pick(-1)

damping = 0.5       # Damping factor to smooth the update
tolerance = 1e-6    # Tolerance to consider convergence
max_iterations = 100
dz = 0.001          # Initial guess for delta_Z
iteration = 0

while iteration < max_iterations:
    rms_value = R_RMS(L, M, N, X, Y, dz)
    rms_derivative = DER_R_RMS(L, M, N, X, Y, dz)
    new_dz = dz - damping * (rms_value / rms_derivative)
    print("Iteration:", iteration,
          "delta_Z:", dz,
          "R_RMS:", rms_value,
          "dR_RMS/dZ:", rms_derivative,
          "new delta_Z:", new_dz)
    
    # If the change is smaller than the tolerance, we consider the algorithm converged.
    if abs(new_dz - dz) < tolerance:
        dz = new_dz
        break
    
    dz = new_dz
    iteration += 1

print("Optimized delta_Z:", dz)

# =============================================================================
# Adjust the image plane position using the computed dz.
#
# Here we assume that the image plane is located at index 3 in the system's SDT.
# (Ensure that this index correctly corresponds to your image plane.)
# =============================================================================
print("Original image plane thickness:", Doublet.SDT[3].Thickness)
Doublet.SDT[3].Thickness = Doublet.SDT[3].Thickness + dz
Doublet.SetData()  # Update the system data after modification
Doublet.SetSolid()  # Update the 3D solid representation if necessary
print("Adjusted image plane thickness:", Doublet.SDT[3].Thickness)

# =============================================================================
# Re-trace the rays after adjusting the image plane.
#
# Create a new ray container and repeat the ray tracing procedure.
# =============================================================================
newRaysTotal = Kos.raykeeper(Doublet)

for j in range(-grid_resolution, grid_resolution + 1):
    for i in range(-grid_resolution, grid_resolution + 1):
        x0 = (i / grid_resolution) * aperture_radius
        y0 = (j / grid_resolution) * aperture_radius
        if np.sqrt(x0**2 + y0**2) < aperture_radius:
            angle_deg = 0.0
            pSource = [x0, y0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(angle_deg)), np.cos(np.deg2rad(angle_deg))]
            
            # Re-trace for wavelength 0.4
            wavelength = 0.4
            Doublet.Trace(pSource, dCos, wavelength)
            newRaysTotal.push()
            
            # Re-trace for wavelength 0.5
            wavelength = 0.5
            Doublet.Trace(pSource, dCos, wavelength)
            newRaysTotal.push()
            
            # Re-trace for wavelength 0.6
            wavelength = 0.6
            Doublet.Trace(pSource, dCos, wavelength)
            newRaysTotal.push()

# =============================================================================
# Display the updated 2D plot after the focus adjustment.
# =============================================================================
Kos.display2d(Doublet, newRaysTotal, 0)
