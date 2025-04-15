#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens"""


import numpy as np
import pkg_resources
import os
""" Looking for if KrakenOS is installed, if not, it assumes that
an folder downloaded from github is run"""

required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("Not installed")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos

LENSCAT_PATH = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + "/KrakenOS/LensCat"


ED_CAT_PATH = os.path.join(LENSCAT_PATH, 'Edmund Optics 2019.ZMF')

THOR_CAT_PATH = os.path.join(LENSCAT_PATH, 'THORLABS.ZMF')

# load lens catalog
PATH1 = os.path.join(LENSCAT_PATH, ED_CAT_PATH)

PATH2 = os.path.join(LENSCAT_PATH, THOR_CAT_PATH)

cat = Kos.zmf2dict([PATH1, PATH2])

for key in cat:
    print(key)
    print(key, '->', cat[key])

print(8)

'''
convert catalog to surface lis7
cat2surf convert the dictionary form of the surfaces to list of surf() class
surflist2dict re-convert the list of surf() class to the catalog dictionary form
'''
surf_list = Kos.cat2surf(cat['84383'], inverse=True, Thickness=0, DespY=5, TiltX=30, AxisMove=0)
cat_dictionary = Kos.surflist2dict(surf_list)








Zf = os.path.join(LENSCAT_PATH, 'zmax_84383.zmx')

Z = Kos.zmx_read(Zf)


config = Kos.Setup()



def rays_2D(system, angle: float=0, src: list=[-12.5, 12.5], num_rays: int=20, wavelength: float=0.55):
    # convert angle degree to radian
    angle = np.deg2rad(angle)

    rays = Kos.raykeeper(system)
    for src_position in src:
        for ang in np.linspace(-angle, angle, num_rays): # in radian
            dCos = [0.0, np.sin(ang), -np.cos(ang)]
            source = [0, src_position, 0]
            system.Trace(source, dCos, wavelength)
            rays.push()
    return rays



angle = np.rad2deg(np.arcsin(0.6))
lens_list = [Kos.surf(Diameter=7, Thickness=0),
              Kos.SurfBlock(Z, inverse=True, name='asphere1', **{'DespY' : 3, 'AxisMove' : 0}),
            Kos.SurfBlock(Z, inverse=False, name='lens0', Thickness=0),
            Kos.SurfBlock(cat['67252'], inverse=True, name='lens1', Thickness=0),
            Kos.SurfBlock(cat['49104'], inverse=False, name='lens2', Thickness=0),
            Kos.surf(Diameter=17.6, Thickness=1.05*np.sqrt(2), Rc=0, Glass='F_SILICA', TiltX=45, AxisMove=0),
            Kos.surf(Diameter=17.6, Thickness=10, TiltX=45, AxisMove=0, ShiftY=0, Name='dichroic'),
            Kos.surf(Diameter=25.4, Thickness=0)]


magnification = 37.5/15
sensor = 1.4
print(f'source position {sensor/magnification/2:.3f} | source size {sensor/magnification:.3f} | sensor size {sensor} | magnification {magnification:.2f} |(mm)')
wavelength = 0.55

'''
Alig
'''
Distances = {0 : 33.72, 'asphere1' : 10.14282592, 'lens0' : 41.4911102, 'lens1' : 10.59216003, 'lens2' : 16.82, 'dichroic' : 18.13360847}
align = Kos.system(Kos.alignment(lens_list, Distances), config)

# ray tracing
src = [-sensor/magnification/2, 0, sensor/magnification/2]
src = [0]
ray = rays_2D(align, angle, src, num_rays=200, wavelength=wavelength)


fig = Kos.display2d(align, ray, figsize=(20, 8))
xlim = [0, 180]
ylim = [-15, 15]

# plt.show()


