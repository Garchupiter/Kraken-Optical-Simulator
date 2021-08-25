#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import glob

myFiles = glob.glob('*.py')

for i in range(0, len(myFiles)):
    if myFiles[i] != "TestAll.py" or myFiles[i] != "Examp-MultiCore.py":
        filePy = myFiles[i]
        print(filePy)
        exec(open(filePy).read())

