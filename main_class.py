#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 11:39:16 2018

@author: admaren01
"""

from input_shipdata import readinputs
import math
from matplotlib import pyplot as plt

class resistance_calc(object):
    
    def __init__(self):
        print("init")
        self.vol = None
        self.l = None
        self.b = None
        self.tf = None
        self.ta = None
        self.t = None
        self.lcb = None
        self.abt = None
        self.hb = None
        self.cm = None
        self.cwp = None
        self.at = None
        self.sapp = None
        self.cstren = None
        self.v = None
    def CB(self, vol, l, b, t):
        
        
        """to find block coefficient"""
        return self.vol/(self.l*self.b*self.t) 
 
    def reynolds_no(self):
        
        """ calculates kinematic viscosity
        v : velocity in m/s
        l : characteristic length in m
        mu : kinematic viscosity in si """
        nu = 1.188* (10**-6)
        return (self.v * self.l) / nu

    
    def fric_coeff(self, rn):
            
           """
           calculated the coefficient of frictional resistance
           rn : reynolds no. in SI """
           return 0.075/((math.log10(self.rn) - 2)**2)
 
    def wetarea(self, vol, l, t, b, cm, cwp, abt):
           """calculates the total wetted surface area"""  
           self.cb = CB(vol, l, b, t)
           self.bt = b/t# B/T
           self.ac = abt/ cb#Abt/Cb
           self.S = l*((2*t)+ b)*math.sqrt(cm)*((0.453+(0.4425*cb))-(0.2862*cm)-\
     (0.003467*bt)+(0.3696*cwp))+(2.38*ac)
           return S    
             
    def RF(self, vol, l, t, b, cm, cwp, abt, v):
         """ Calculates frictional resistance"""    
         self.rho = 1.025
         self.S = wetarea(vol, l, t, b, cm, cwp, abt)
         #print(round(S, 2))
         self.rn = reynolds_no(l, v)
         self.cf = fric_coeff(rn)
        # print('the coeff of friction is :',cf)
         print (0.5 * cf * rho * S * (v**2))
         return 0.5 * cf * rho * S * (v**2)
     
    def AM(self, b, t, cm):
         """
         calculates area of midship at lwl
         b: breadth
         t: draft
         cm: coeff of midship
         """
         return b*t*cm
     
    def CP(self, vol, l, b, t, cm):
         """ 
         calculates prismatic coeff
         vol : volume of displacement in m3
         l : lwl m
          b: breadth
         t: draft
         cm: coeff of midship
         """
         self.am = AM(b, t, cm)
         return vol/(l*am)
         
    def LR(self, vol, l, b, t, cm, lcb):
         """ 
         calculated length of run
         vol : volume of displacement in m3
         l : lwl m
          b: breadth
         t: draft
         cm: coeff of midship
         lcb : longitudanal; position of center of boyancy frwd of 0.5l as 
               percentage of l
         """     
         self.cp = CP(vol, l, b, t, cm)
         return (1-cp+0.06*cp*lcb/(4*cp-1))*l
     
    def C12(self, t, l):
         """
         calculates coeff c12 
         t : mean draft (m)
         l : lwl
         """   
         if(t/l) > 0.05:
             return (t/l)**0.2228446
         if(t/l) > 0.02 and(t/l) < 0.05:
             return 48.20*(((t/l)-0.02)**2.078)+0.479948
         if(t/l) < 0.02:
             return 0.479948   
         
    def FF(self, vol, l, b, t, cm, lcb, Cstern=10):
         """
         calculates the form factor
         vol : volume of displacement in m3
         l : lwl m
          b: breadth
         t: draft
         cm: coeff of midship
         lcb : longitudanal; position of center of boyancy frwd of 0.5l as 
               percentage of l
         Cstern : stern shape parameter      
         """
         self.c13 = (1 + (0.003 * Cstern))
         self.c12 = C12(t, l)
         self.lr = LR(vol, l, b, t, cm, lcb)
         self.cp = CP(vol, l, b, t, cm)
         
         return c13*(0.93 + (c12*((b/lr)**0.92497))*((0.95-cp)**(-0.521448))*\
                     (1-cp+0.0225*lcb)**0.6906)
         
     #def K2(sapp) :
        # r1 = float (input('k2 value for rudder behind skeg  1.5...2.0'))
         #r2 = float (input('k2 value for rudder behind stern 1.3....1.5'))
        #t1 = float (input('k2 value for twin screw balance rudders 2.8'))
         #s1 = float (input('k2 value for shaft brackets 3.0'))
         #s2 = float (input('k2 value for skeg  1.5...2.0'))
    def RAPP(self, l, v, sapp, rho=1.025):
         """
         calculates the appendage resistance
         v : speed in m/s
         sapp : wetted surface area of appendages
         rho : density of sea water in t/m3
         """
         self.k2 = float(input('enter value for 1+k2 (1.5)'))
         self.rn = reynolds_no(l, v)
         self.cf = fric_coeff(rn)
         return 0.5*rho*(v**2)*sapp*k2*cf
     
    def IE(self, vol, l, b, t, cm, lcb, cwp):
         
         """
         calculates the half angle of entrance
         vol : volume of displacement in m3
         l : lwl m
          b: breadth
         t: draft
         cm: coeff of midship
         lcb : longitudanal; position of center of boyancy frwd of 0.5l as 
               percentage of l
         cwp : coeff of water plane area
     
         """   
         self.cp = CP(vol, l, b, t, cm)
         self.lr = LR(vol, l, b, t, cm, lcb)
         self.a1 = -(l/b)**0.80856 *((1-cwp)**0.30484)*((1-cp-(0.0225*lcb))**0.6367)\
               * ((lr/b)**0.34574)*(((100 * vol)/(l**3))**0.16302)   
         return 1 +(89 * math.exp(a1))
     
    def C7(self, l, b):
         if b/l < 0.11:
             return 0.229577*((b/l)**0.33333)
         if(b/l > 0.11) and(b/l < 0.25):
             return b/l
         if b/l > 0.25:
             return 0.5 -(0.0625*(l/b))
         
    def C1(self, vol, l, b, t, cm, lcb, cwp):
         
         self.ie = IE(vol, l, b, t, cm, lcb, cwp)
         self.c7 = C7(l, b)
         return 2223105 *(c7**3.78613) *((t/b)**1.07961)*((90-ie)**(-1.37565))
     
    def C3(self, b, t, tf, abt, hb):
         return (0.56*(abt**1.5))/(b*t*(0.31*math.sqrt(abt)+tf-hb))
     
    def C2(self, b, t, tf, abt, hb):
         self.c3 = C3(b, t, tf, abt, hb)
         return math.exp(-1.89*math.sqrt(c3))
     
    def C16(self, vol, l, b, t, cm):
         self.cp = CP(vol, l, b, t, cm) 
         if cp < 0.8:
             return(8.07981*cp)-(13.8673 *(cp**2))+(6.984388*(cp**3))
         else:
             return 1.73014-(0.7067*cp)
         
    def M1(self, vol, l, b, t, cm):
         self.c16 = C16(vol, l, b, t, cm)    
         return (0.0140407*(l/t))-(1.75254*(vol**0.333333)/l)-(4.79323*b/l)-c16 
     
    def C15(self, vol, l):    
         if (l**3)/vol < 512:
             return -1.69385
         if (l**3)/vol > 1727:
             return 0.0
         if ((l**3)/vol > 512) and ((l**3)/vol < 1727):
             return -1.69385 + (l/(vol**0.33333) - 8.0)/2.36
         
    def FN(self, v, l):
         """
         calculates the froude number 
         v : velocity of ship in m/s
         l : lwl m
         """
         return v/math.sqrt(9.81*l)
     
    def M2(self, vol, l, b, t, cm, v):
         self.c15 = C15(vol, l)
         self.cp = CP(vol, l, b, t, cm)
         self.fn = FN(v, l) 
         return c15*(cp**2)*(math.exp(-0.1*(fn**-2)))
     
    def LAMBDA(self, vol, l, b, t, cm):
         self.cp = CP(vol, l, b, t, cm)
         if l/b < 12:
             return 1.446*cp - 0.03*l/b
         else:
             return 1.446*cp - 0.36
             
    def RW(self, vol, l, b, t, cm, lcb, cwp, tf, abt, hb, at, v, rho=1.025): 
         """
         claculates the wave making resistance 
         calls in the required coeff functions
         vol : volume of displacement in m3
         l : lwl m
         b: breadth m
         t: draft m
         cm: coeff of midship
         lcb : longitudanal; position of center of boyancy frwd of 0.5l as 
               percentage of l
         cwp : coeff of water plane area
         tf : the fprward draft
         abt : transverse bulb area
         hb : center of bulb area above keel line
         at : transon area
         v : velocity of ship in m/s
         rho is taken as 1.025 t/m3
         """
         self.c1 = C1(vol, l, b, t, cm, lcb, cwp)
         self.c2 = C2(b, t, tf, abt, hb)
         self.c5 = 1-(0.8*at/(b*t*cm))
         self.m1 = M1(vol, l, b, t, cm)
         self.m2 = M2(vol, l, b, t, cm, v)
         self.fn = FN(v, l)
         self.lam = LAMBDA(vol, l, b, t, cm)
         return c1*c2*c5*vol*rho*9.81*math.exp(m1*(fn**-0.9)+ \
                                               m2*math.cos(lam*(fn ** -2)))
     
    def PB(self, abt, tf, hb):
         return 0.56*math.sqrt(abt)/(tf - 1.5*hb)
     
    def FNI(self, v, tf, hb, abt):
         return v/math.sqrt(9.81*(tf-hb-(0.25*math.sqrt(abt)))+0.15*(v**2))
     
    def RB(self, v, tf, hb, abt, rho=1.025):
         """ 
         calculates additonal resistance of bulbous bow at the surface
         v : velocity of ship in m/s
         tf : the fprward draft
         hb : center of bulb area above keel line
         abt : transverse bulb area
         rho is taken as 1.025 t/m3
         """
         self.pb = self.PB(abt, tf, hb)
         self.fni = self.FNI(v, tf, hb, abt)    
         return 0.11*math.exp(-3*(pb**-2))*(fni**3)*\
                                              (abt**1.5)*rho*9.81/(1+ (fni**2))
                                              
    def FNT(self):
         return self.v/math.sqrt(2*9.81*self.at/(self.b+self.b*self.cwp))
     
    def C6(self):
         fnt = self.FNT()
         if fnt >= 5:
             return 0
         else:
             return 0.2*(1-0.2*fnt)
     
    def RTR(self, v, b, cwp, at, rho=1.025):
         self.c6 = C6(v, at, b, cwp)
         return 0.5*rho*(v**2)*at*c6
     
    def C4(self, l, tf):
         if tf/l <= 0.04:
             return tf/l
         else:
             return 0.04
         
    def CA(self, vol, l, b, t, tf, abt, hb):    
         self.cb = CB(vol, l, b, t)
         self.c2 = C2(b, t, tf, abt, hb)
         self.c4 = C4(l, tf)
         return 0.006*((l+100)**-0.16)- 0.00205 + 0.003*\
                math.sqrt(l/7.5)*(cb**4)*c2*(0.04-c4)
      
    def RA(self, vol, l, b, t, tf, cm, cwp, abt, hb, v, rho=1.025):
         self.ca = CA(vol, l, b, t, tf, abt, hb)
         self.s = wetarea(vol, l, t, b, cm, cwp, abt)
         return 0.5*rho*(v**2)*s*ca
     
    def RTOT(self, vol, l, b, t, cm, lcb, cwp, tf, abt, hb, sapp, at, v, Cstern, \
                                                                    rho=1.025):
         self.rf = RF(vol, l, t, b, cm, cwp, abt, v)
         self.k1 = FF(vol, l, b, t, cm, lcb, Cstern)
         self.rapp = RAPP(l, v, sapp, rho)
         self.rw = RW(vol, l, b, t, cm, lcb, cwp, tf, abt, hb, at, v, rho)
         self.rb = RB(v, tf, hb, abt, rho)
         self.rtr = RTR(v, b, cwp, at, rho)
         self.ra = RA(vol, l, b, t, tf, cm, cwp, abt, hb, v, rho) 
# =============================================================================
# =============================================================================
         return rf*k1 + rapp + rw + rb + rtr + ra
   
    def Calculate(self):
           self.Rt = self.RTOT(vol, l, b, t, cm, lcb, cwp, tf, abt, \
                               hb, sapp, at, v, Cstern, rho=1.025)
           print("the total resistance is : ",Rt)
   
if __name__ == "__main__":
    
    
     
                   
     res1 = resistance_calc()
     #inputdict = readinputs()
     res1.vol = 37500
     res1.l = 205
     res1.b = 32
     res1.tf = 10
     res1.ta = 10
     res1.t = (res1.ta + res1.tf)/2
     res1.lcb = -0.75
     res1.abt = 20
     res1.hb = 4
     res1.cm = 0.98
     res1.cwp = 0.75
     res1.at = 16
     res1.sapp = 50
     res1.cstren = 10
     res1.v = 12.86

     #res1.Calculate()
     print("finish")


#     
#    vknots = [v for v in range(5,35)]
#    vms = [0.5144 * v  for v in vknots]
#    #rfvals = [RF(37500, 205, 10, 32, 0.98, 0.75, 20, v)for v in vms]
#    #plt.plot(vknots, rfvals, "+-")
#    rwvals = [RW(37500, 205, 32, 10, 0.98, -0.75, 0.75, 10, 20, 4, 16, v)for\
#                                                                       v in vms]
#    plt.plot(vknots, rwvals, "+-")
#   # Rf= RF(12.86)
#   # print(Rf)
# =============================================================================
   
   

    
    
    