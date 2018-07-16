# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 00:31:23 2018

@author: Admaren02
"""

from itertools import groupby
import csv

def duppoints(filepath, xloc = 0):
    with open(filepath, 'r') as _file:
        x = filepath
        linelist = _file.readlines()
        out_file = open('%smodified'%x, 'w')
        

        keyfunc = lambda p: p[:2]
        mypoints = [(float(x.strip()), float(y.strip()),float(z.strip()) )\
           for(x,y,z) in [ v.strip().split(",")  for v  in linelist ]]
        for k, g in groupby(sorted(zip(x, y, z), key=keyfunc), keyfunc):
            out_file.append(list(g)[0])
        out_file.close()
    _file.close()
    

if __name__== "__main__":
    
    duppoints(r"C:\Users\Admaren02\Desktop\ship sections\aft\xlocation900.00.csv")