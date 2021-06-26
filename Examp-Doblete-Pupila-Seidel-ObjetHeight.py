#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""
import Kraken as kn

##############################################################    

P_Obj = kn.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 50
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0
P_Obj.Name = "P Obj"

L1a = kn.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

L1b = kn.surf()
L1b.Rc = -3.071608670000159E+001

L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

L1c = kn.surf()
L1c.Rc = -78.184019
L1c.Thickness = 97.375174 - 40
L1c.Glass = "AIR"
L1c.Diameter = 30

pupila = kn.surf()
pupila.Rc = 0
pupila.Thickness = 40.
pupila.Glass = "AIR"
pupila.Diameter = 10
pupila.Name = "Ap Stop"
pupila.DespY = 0.

P_Ima = kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 20.0
P_Ima.Name = "Plano imagen"

A = [P_Obj, L1a, L1b, L1c, P_Ima]

config_1 = kn.Kraken_setup()

Doblete = kn.system(A, config_1)

sup = 1
W = .6
AperType = "EPD"
AperVal = 20  # L1a.Diameter
field = 4.0
fieldType = "height"

AB = kn.Seidel(Doblete, sup, W, AperType, AperVal, field, fieldType)
print("Aberraciones: ", AB)
