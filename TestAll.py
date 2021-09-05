#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import glob

myFiles = glob.glob('*.py')

for i in myFiles:
    if i == "TestAll.py" or i == "Examp_MultiCore.py" or i == "Examp_RonchiTest.py":
        print("---------",i, "No se ejecuta")
    else:
        print("-------------------")
        print(i)
        exec(open(i).read())

