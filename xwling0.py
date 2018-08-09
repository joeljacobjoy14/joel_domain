# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 14:43:56 2018

@author: Admaren02
"""

import xlwings as xw
import os
import re

def parselocation(locstr, xref = 0):
    patt = re.compile(r"(\d*\.?\d*)\s?([fF|aA|sS|pP|uU])")

    match = patt.match(locstr)
    if match:
        val = match.group(1)
        flg =  match.group(2)
        if flg == "A" or flg =="a":
            return xref - float(val)
        elif flg =="f" or flg =="F":

            return xref + float(val)
        elif flg =="s"or flg =="S":
            return  float(val)
        elif flg == "p"or flg =="P":
            return float(val)*-1.0
        elif flg == "u"or flg == "U":
            return float(val)
    else:
        try:
            return float(locstr)
        except:
            print( "Error in parsing value : %s"%locstr)
            return locstr

def conv(arg):
    cmd = r'C:\Users\Admaren02\Desktop\tank_soundings\RTF2Text.exe %s %s'%(arg, arg[:-3]+"txt")
    print(cmd)
    os.system(cmd)
    return arg[:-3]+"txt"




def processrtfout(rtfpath):    
    """Convert rtf text into dataframe
    """
    hdr = ['SOUNDING(m)', 'VOLUME','WEIGHT','LCG','TCG','VCG', 'FSM',\
                     'ULLAGE(m)']
    txtpath = conv('"'+rtfpath+'"')
    print( txtpath)
    _f = open(txtpath.replace('"',''), 'r')
    tab = []
    _flg = 0
    for line in _f:
        if "FULL" in line:
            
            break
        if _flg:
            l = line.strip().replace("%","").split("$")

            lvals = [parselocation(v.replace(" ","").replace(",","").strip()) \
                     for v in l if v.replace(" ","").replace(",","").strip()]
            if len(lvals) == 8:
                tab.append(lvals)

        if line.startswith("Sounding"): # start point
            _f.readline()
            _flg  = 1
#    print tab
   
    return tab
    
def excelsheet(rtfpath):
    
    tab = processrtfout(rtfpath)

#    print(tab)
    
    xw.Book()
    for i,l in enumerate(tab):
        for j, v in enumerate(l):
            xw.Range((i+1, j+1)).value = v



excelsheet(r"C:\Users\Admaren02\Desktop\tank_soundings\BTK6S.RTF")



