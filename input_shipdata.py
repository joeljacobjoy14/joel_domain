#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 11:13:39 2018

@author: admaren01
"""

def readinputs( inputfilename = "holtropinput.txt"):
    """Read input data from the file
    
    inputfilename :
        
        output : [dictionay] {lwl
                              lbp 
                              b 
                              tf
                              ta
                              vol
                              lcb as percentage
                              abt tranverse bulb area
                              center of bulb area above keel
                              cm
                              cwp
                              at
                              Sapp
                              stern shape parameter
                              D
                              Z
                              clearance from keel line
                              V
                              }
    
    """
    output_dic = {}
    _file = open(inputfilename, 'r')
    lines = _file.readlines()
    output_dic ['lwl'] =   float(lines[0].strip()) 
    output_dic['lpp'] = float(lines[1].strip())
    output_dic['b'] = float(lines[2].strip())
    output_dic['tf'] = float(lines[3].strip())
    output_dic['ta'] = float(lines[4].strip())
    output_dic['vol'] = float(lines[5].strip())
    output_dic['lcb'] = float(lines[6].strip())
    output_dic['abt'] = float(lines[7].strip())
    output_dic['hb'] = float(lines[8].strip())
    output_dic['cm'] = float(lines[9].strip())
    output_dic['cwp'] = float(lines[10].strip())
    output_dic['at'] = float(lines[11].strip())
    output_dic['sapp'] = float(lines[12].strip())
    output_dic['cstern'] = float(lines[13].strip())
    output_dic['d'] = float(lines[14].strip())
    output_dic['z'] = float(lines[15].strip())
    output_dic['prpclr'] = float(lines[16].strip())
    output_dic['v'] = float(lines[17].strip())
    
    
    return output_dic

#print("Current file -->:" + __name__)
    

if __name__ == "__main__":
    holtropinput = readinputs()
    print(holtropinput)
    # {"l": 205.0, "b" : 30 . ....}
    #lpp = float(holtropinput["lpp"])