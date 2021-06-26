#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 08:46:43 2021

@author: joelherreravazquez
"""

from multiprocessing import Process
 
 
def square(x):
 
    for x in range(0,10):
        print('%s squared  is  %s' % (x, x**2))
 
 
def is_even(x):
    
    for x in range(0,10):
        if x % 2 == 0:
            print('%s is an even number ' % x)
 
 
if __name__ == '__main__':
    
    
    p1 = Process(target=square, args=('x',))
    p2 = Process(target=is_even, args=('x',))
 
    p1.start()
    p2.start()
 
    p1.join()
    p2.join()
 
    print ("Done")
     