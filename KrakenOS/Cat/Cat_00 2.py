# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 22:23:29 2020

@author: JOELHERRERAVAZQUEZ
"""

import glob

txtfiles = []
for file in glob.glob("*.AGF"):
    txtfiles.append(file)

print(txtfiles)
