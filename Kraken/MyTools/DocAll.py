#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import glob
import os

# https://github.com/heavenshell/py-doq

myFiles = glob.glob('*.py')

for i in range(0, len(myFiles)):
    if myFiles[i] != "DocAll.py" or myFiles[i] !="__init__.py":
        filePy = myFiles[i]
        ToDoc = "(cat "+filePy+" | doq)>New/"+filePy
        os.system(ToDoc)
        print(ToDoc)


