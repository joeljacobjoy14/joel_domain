#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 11:39:16 2018

@author: admaren01
"""

from input_shipdata import readinputs
import math
from matplotlib import pyplot as plt
inputdict = readinputs()

print(inputdict)

def CB(vol, l, b, t):
    '''to find block coefficient'''
    
    return vol/(l*b*t)  

def reynolds_no(v):
    """ calculates kinematic viscosity
    
    v : velocity in m/s
    l : characteristic length in m
    mu : kinematic viscosity in si 
    """
    l = inputdict['lwl']
    nu = 1.158* (10**-6)
    return (v * l) / nu

def fric_coeff(rn):
    '''
    calculated the coefficient of frictional resistance
    rn : reynolds no. in SI
    '''
    return 0.075/((math.log10(rn) - 2)**2)

def wetarea(vol, l, t, b, cm, cwp, abt ):
    
    cb = CB(vol, l, b, t)
    bt=  b/t# B/T
    ac = abt/ cb#Abt/Cb
    
    S =  l*((2*t )+ b)*math.sqrt(cm)*\
    ((0.453+(0.4425*cb))-(0.2862*cm)-(0.003467*bt)\
            +(0.3696*cwp))+(2.38*ac)
    return S
    
        

def RF(vol, l, t, b, cm, cwp, abt, v):
    """ Calculates frictional resistance
    
    """
    rho = 1.025
    S = wetarea(vol, l, t, b, cm, cwp, abt)
    #print(round(S, 2))
    rn = reynolds_no(v)
    cf = fric_coeff(rn)
   # print('the coeff of friction is :',cf)
    return ((0.5 * cf * rho * S * (v**2)))

def AM(b,t,cm):
    '''
    calculates area of midship at lwl
    b: breadth
    t: draft
    cm: coeff of midship
    '''
    return (b*t*cm)

def CP(vol,l,b,t,cm):
    ''' 
    calculates prismatic coeff
    vol : volume of displacement in m3
    l : lwl m
     b: breadth
    t: draft
    cm: coeff of midship
    '''
    am= AM(b,t,cm)
    return(vol/(l*am))
    
def LR(vol,l,b,t,cm,lcb):
    ''' 
    calculated length of run
    vol : volume of displacement in m3
    l : lwl m
     b: breadth
    t: draft
    cm: coeff of midship
    lcb : longitudanal; position of center of boyancy frwd of 0.5l as 
          percentage of l
    '''      
    
    cp= CP(vol,l,b,t,cm)
    return ((1-cp+0.06*cp*lcb/(4*cp-1))*l)

 
def C12(t,l):
    '''
    calculates coeff c12 
    t : mean draft (m)
    l : lwl
    '''
    
    if((t/l)>0.05):
        return (t/l)**0.2228446
    if ((t/l)>0.02 and (t/l)<0.05):
        return 48.20*(((t/l)-0.02)**2.078 )+0.479948
    if ((t/l)<0.02):
        return 0.479948   

    
def FF(vol,l,b,t,cm,lcb):
    '''
    calculates the form factor
    vol : volume of displacement in m3
    l : lwl m
     b: breadth
    t: draft
    cm: coeff of midship
    lcb : longitudanal; position of center of boyancy frwd of 0.5l as 
          percentage of 
          
    '''
 
    Cstern = float (input('enter value for stern coefficient(10, 0, -10)'))
    c13 = (1 + (0.003 * Cstern))
    c12 = C12(t,l)
    lr = LR(vol,l,b,t,cm,lcb)
    cp = CP(vol,l,b,t,cm)
    
    return c13*(0.93 + (c12*((b/lr)**0.92497))*((0.95-cp)**(-0.521448))*\
                (1-cp+0.0225*lcb)**0.6906)
    
#def K2(sapp) :
   # r1 = float (input('k2 value for rudder behind skeg  1.5...2.0'))
    #r2 = float (input('k2 value for rudder behind stern 1.3....1.5'))
   #t1 = float (input('k2 value for twin screw balance rudders 2.8'))
    #s1 = float (input('k2 value for shaft brackets 3.0'))
    #s2 = float (input('k2 value for skeg  1.5...2.0'))
def RAPP(v, sapp, rho= 1.025):
    
    '''
    calculates the appendage resistance
    v : speed in m/s
    sapp : wetted surface area of appendages
    rho : density of sea water in t/m3
    '''

    k2 = float (input('enter value for 1+k2 (1.5)'))
    rn = reynolds_no(v)
    cf = fric_coeff(rn)
    
    return 0.5*rho*(v**2)*sapp*k2*cf

def IE(vol,l,b,t,cm,lcb,cwp):
    
    '''
    calculates the half angle of entrance
    
    vol : volume of displacement in m3
    l : lwl m
     b: breadth
    t: draft
    cm: coeff of midship
    lcb : longitudanal; position of center of boyancy frwd of 0.5l as 
          percentage of 
    cwp : coeff of water plane area

    '''  
    
    
    cp = CP(vol,l,b,t,cm)
    lr = LR(vol,l,b,t,cm,lcb)
    
    a1 = ((-(l/b)**0.80856 ) * ((1-cwp)**0.30484)*((1-cp-(0.0225*lcb))**0.6367)\
          * ((lr/b)**0.34574)*(((100 * vol)/(l**3))**0.16302))
    
    
    return 1 + (89 * math.exp(a1))


def C7(l,b):
    
    if (b/l < 0.11):
        return (0.229577*((b/l)**0.33333))
    if ((b/l >0.11) and (b/l< 0.25)):
        return b/l
    if (b/l > 0.25):
        return (0.5 - (0.0625*(l/b)))
    
        

def C1(vol,l,b,t,cm,lcb,cwp):
    
    ie = IE(vol,l,b,t,cm,lcb,cwp)
    c7 = C7(l,b)
    
    return 2223105 * (c7**3.78613) * ((t/b)**1.07961)*((90-ie)**(-1.37565))

def C3(b,t,tf,abt,hb):
    
    return (0.56*(abt**1.5))/(b*t*(0.31*math.sqrt(abt)+tf-hb))

def C2(b,t,tf,abt,hb):
    c3 = C3(b,t,tf,abt,hb)
    
    return math.exp(-1.89*math.sqrt(c3))


def C16(vol, l, b, t, cm):
    cp = CP(vol, l, b, t, cm)
    
    if cp < 0.8:
        return (8.07981*cp )- (13.8673 * (cp**2))+ (6.984388*(cp**3))
    else:
        return (1.73014-(0.7067*cp))
    
def M1(vol, l, b, t, cm):
    c16 = C16(vol, l, b, t, cm)
    
    return (0.0140407*(l/t))-(1.75254*(vol**0.333333)/l)-(4.79323*b/l)-c16 

def C15(vol,l):
    
    if (l**3)/vol < 512:
        return -1.69385
    if (l**3)/vol > 1727:
        return 0.0
    if ((l**3)/vol> 512) and ((l**3)/vol<1727):
        return -1.69385 + (l/(vol**0.33333) - 8.0)/2.36
    
def FN(v,l):
    '''
    v : velocity of ship in m/s
    l : lwl m
    '''
    return v/math.sqrt(9.81*l)

def M2(vol,l,b,t,cm,v):
    c15 = C15(vol,l)
    cp = CP(vol,l,b,t,cm)
    fn = FN(v,l)
    
    return (c15*(cp**2)*(math.exp(-0.1*(fn**-2))))

def LAMBDA(vol, l, b, t, cm):
    cp = CP(vol, l, b, t, cm)
    if (l/b < 12):
        return (1.446*cp - 0.03*l/b)
    else:
        return 1.446*cp - 0.36
    


    
    
def RW(vol,l,b,t,cm,lcb,cwp,tf,abt,hb,at ,v= 12.86, rho = 1.025):
    
    c1 = C1(vol,l,b,t,cm,lcb,cwp)
    c2 = C2(b,t,tf,abt,hb)
    c5 = 1-(0.8*at/(b*t*cm))
    m1 = M1(vol, l, b, t, cm)
    m2 = M2(vol,l,b,t,cm,v)
    fn = FN(v,l)
    lam = LAMBDA(vol, l, b, t, cm)
    return c1*c2*c5*vol*rho*9.81*math.exp(m1*(fn**-0.9)+ \
                                          m2*math.cos(lam*(fn ** -2)))
    

def PB(abt, tf, hb):
    
    return 0.56*math.sqrt(abt)/(tf - 1.5*hb)

def FNI(v, tf, hb, abt):
    return v/math.sqrt(9.81*(tf-hb-(0.25*math.sqrt(abt)))+0.15*(v**2))

def RB(v, tf, hb, abt,rho = 1.025):
    pb = PB(abt, tf, hb)
    fni = FNI(v, tf, hb, abt)
    
    return 0.11*math.exp(-3*(pb**-2))*(fni**3)*\
                                         (abt**1.5)*rho*9.81/(1+ (fni**2))
                                         
def FNT(v,at,b,cwp):
    
    return v/math.sqrt(2*9.81*at/(b+b*cwp))

def C6(v,at,b,cwp):
    
    fnt = FNT(v,at,b,cwp)
    if fnt >= 5:
        return 0
    else:
        return 0.2*(1-0.2*fnt)

def RTR(v,b,cwp,at,rho=1.025):
    
    c6 = C6(v,at,b,cwp)
    return 0.5*rho*(v**2)*at*c6

def C4(l, tf):
    
    if tf/l <=0.04:
        return tf/l
    else:
        return 0.04
    
def CA(vol, l, b, t, tf, abt, hb):
    
    cb = CB(vol, l, b, t)
    c2 = C2(b, t, tf, abt, hb)
    c4 = C4(l, tf)
    
    return 0.006*((l+100)**-0.16)- 0.00205 + 0.003*math.sqrt(l/7.5)*(cb**4)*c2\
                                                *(0.04-c4)
 
def RA(vol, l, b, t, tf, cm, cwp, abt, hb, v, rho=1.025):
    
    ca = CA(vol, l, b, t, tf, abt, hb)
    s = wetarea(vol, l, t, b, cm, cwp, abt )
    
    return 0.5*rho*(v**2)*s*ca

def RTOT(vol, l, b, t, cm, lcb, cwp, tf, abt, hb, at, v, rho = 1.025):
    rf = RF(vol, l, t, b, cm, cwp, abt, v)
    k1 = FF(vol,l,b,t,cm,lcb)
    rapp = RAPP(v, sapp, rho)
    rw = RW(vol,l,b,t,cm,lcb,cwp,tf,abt,hb,at ,v, rho)
    rb = RB(v, tf, hb, abt,rho)
    rtr = RTR(v,b,cwp,at,rho)
    ra = RA(vol, l, b, t, tf, cm, cwp, abt, hb, v, rho)
    
    return rf*k1 + rapp + rw + rb + rtr + ra

    
    
                                              
        
    
    



    
    



        
    


    

                
    
    
    

    

if __name__ == "__main__":                           
    #S= 7381.45
    
    pass
   # vknots = [v for v in range(10,35)]
    #vms = [0.5144 * v  for v in vknots]
   # rfvals = [RF(v, l, S) for v in vms]
   # plt.plot(vknots, rfvals, "+-")
  # Rf= RF(12.86)
  # print(Rf)
   
   

    
    
    