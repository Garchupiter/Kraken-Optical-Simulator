#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import glob
import os



myFiles = glob.glob('*.py')

for i in range(0, len(myFiles)):
    if myFiles[i] != "DocAll.py" or myFiles[i] !="__init__.py":
        filePy = myFiles[i]
        ToDoc = "(cat "+filePy+" | doq --formatter=google)>New_"+filePy
        os.system(ToDoc)
        print(ToDoc)


