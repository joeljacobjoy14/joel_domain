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
    def CB(self):
        
        
        """to find block coefficient"""
        return self.vol/(self.l*self.b*self.t) 
 
    def reynolds_no(self):
        
        """ calculates kinematic viscosity
        v : velocity in m/s
        l : characteristic length in m
        mu : kinematic viscosity in si """
        nu = 1.188* (10**-6)
        return (self.v * self.l) / nu

    
    def fric_coeff(self):
            
           """
           calculated the coefficient of frictional resistance
           rn : reynolds no. in SI """
           rn = self.reynolds_no()
           return 0.075/((math.log10(rn) - 2)**2)
 
    def wetarea(self):
           """calculates the total wetted surface area"""  
           cb = self.CB()
           bt = self.b/self.t# B/T
           ac = self.abt/ cb#Abt/Cb
           S = self.l*((2*self.t)+ self.b)*math.sqrt(self.cm)*((0.453+(0.4425*\
cb))-(0.2862*self.cm)-(0.003467*bt)+(0.3696*self.cwp))+(2.38*ac)
           return S    
             
    def RF(self):
         """ Calculates frictional resistance"""    
         rho = 1.025
         S = self.wetarea()
         #print(round(S, 2))
         cf = self.fric_coeff()
        # print('the coeff of friction is :',cf)
        # print (0.5 * cf * rho * S * (self.v**2))
         return 0.5 * cf * rho * S * (self.v**2)
     
    def AM(self):
         """
         calculates area of midship at lwl
         b: breadth
         t: draft
         cm: coeff of midship
         """
         return self.b*self.t*self.cm
     
    def CP(self):
         """ 
         calculates prismatic coeff
         vol : volume of displacement in m3
         l : lwl m
          b: breadth
         t: draft
         cm: coeff of midship
         """
         am = self.AM()
         return self.vol/(self.l*am)
         
    def LR(self):
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
         cp = self.CP()
         return (1-cp+0.06*cp*self.lcb/(4*cp-1))*self.l
     
    def C12(self):
         """
         calculates coeff c12 
         t : mean draft (m)
         l : lwl
         """   
         if(self.t/self.l) > 0.05:
             return (self.t/self.l)**0.2228446
         if(self.t/self.l) > 0.02 and(self.t/self.l) < 0.05:
             return 48.20*(((self.t/self.l)-0.02)**2.078)+0.479948
         if(self.t/self.l) < 0.02:
             return 0.479948   
         
    def FF(self):
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
         c13 = (1 + (0.003 * self.cstern))
         c12 = self.C12()
         lr = self.LR()
         cp = self.CP()
         
         return c13*(0.93 + (c12*((self.b/lr)**0.92497))*((0.95-cp)**\
                             (-0.521448))*(1-cp+0.0225*self.lcb)**0.6906)
         
     #def K2(sapp) :
        # r1 = float (input('k2 value for rudder behind skeg  1.5...2.0'))
         #r2 = float (input('k2 value for rudder behind stern 1.3....1.5'))
        #t1 = float (input('k2 value for twin screw balance rudders 2.8'))
         #s1 = float (input('k2 value for shaft brackets 3.0'))
         #s2 = float (input('k2 value for skeg  1.5...2.0'))
    def RAPP(self):
         """
         calculates the appendage resistance
         v : speed in m/s
         sapp : wetted surface area of appendages
         rho : density of sea water in t/m3
         """
         rho = 1.025
         k2 = float(input('enter value for 1+k2 (1.5)'))
         rn = self.reynolds_no()
         cf = self.fric_coeff()
         return 0.5*rho*(self.v**2)*self.sapp*k2*cf
     
    def IE(self):
         
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
         cp = self.CP()
         lr = self.LR()
         a1 = -(self.l/self.b)**0.80856 *((1-self.cwp)**0.30484)*((1-\
               cp-(0.0225*self.lcb))**0.6367)* ((lr/self.b)**0.34574)*\
               (((100 * self.vol)/(self.l**3))**0.16302)   
         return 1 +(89 * math.exp(a1))
     
    def C7(self):
         if self.b/self.l < 0.11:
             return 0.229577*((self.b/self.l)**0.33333)
         if(self.b/self.l > 0.11) and(self.b/self.l < 0.25):
             return self.b/self.l
         if self.b/self.l > 0.25:
             return 0.5 -(0.0625*(self.l/self.b))
         
    def C1(self):
         
         ie = self.IE()
         c7 = self.C7()
         return 2223105 *(c7**3.78613) *((self.t/self.b)**1.07961)*((90-ie)\
                                                                **(-1.37565))
     
    def C3(self):
         return (0.56*(self.abt**1.5))/(self.b*self.t*(0.31*\
                                       math.sqrt(self.abt)+self.tf-self.hb))
     
    def C2(self):
         c3 = self.C3()
         return math.exp(-1.89*math.sqrt(c3))
     
    def C16(self):
         cp = self.CP() 
         if cp < 0.8:
             return(8.07981*cp)-(13.8673 *(cp**2))+(6.984388*(cp**3))
         else:
             return 1.73014-(0.7067*cp)
         
    def M1(self):
         c16 = self.C16()    
         return (0.0140407*(self.l/self.t))-(1.75254\
                *(self.vol**0.333333)/self.l)-(4.79323*self.b/self.l)-c16 
     
    def C15(self):    
         if (self.l**3)/self.vol < 512:
             return -1.69385
         if (self.l**3)/self.vol > 1727:
             return 0.0
         if ((self.l**3)/self.vol > 512) and ((self.l**3)/self.vol < 1727):
             return -1.69385 + (self.l/(self.vol**0.33333) - 8.0)/2.36
         
    def FN(self):
         """
         calculates the froude number 
         v : velocity of ship in m/s
         l : lwl m
         """
         return self.v/math.sqrt(9.81*self.l)
     
    def M2(self):
         c15 = self.C15()
         cp = self.CP()
         fn = self.FN() 
         return c15*(cp**2)*(math.exp(-0.1*(fn**-2)))
     
    def LAMBDA(self):
         cp = self.CP()
         if self.l/self.b < 12:
             return 1.446*cp - 0.03*self.l/self.b
         else:
             return 1.446*cp - 0.36
             
    def RW(self): 
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
         rho = 1.025
         c1 = self.C1()
         c2 = self.C2()
         c5 = 1-(0.8*self.at/(self.b*self.t*self.cm))
         m1 = self.M1()
         m2 = self.M2()
         fn = self.FN()
         lam = self.LAMBDA()
         return c1*c2*c5*self.vol*rho*9.81*math.exp(m1*(fn**-0.9)+ \
                                               m2*math.cos(lam*(fn ** -2)))
     
    def PB(self):
         return 0.56*math.sqrt(self.abt)/(self.tf - 1.5*self.hb)
     
    def FNI(self):
         return self.v/math.sqrt(9.81*(self.tf-self.hb-(0.25*\
                                        math.sqrt(self.abt)))+0.15*(self.v**2))
     
    def RB(self):
         """ 
         calculates additonal resistance of bulbous bow at the surface
         v : velocity of ship in m/s
         tf : the fprward draft
         hb : center of bulb area above keel line
         abt : transverse bulb area
         rho is taken as 1.025 t/m3
         """
         rho = 1.025
         pb = self.PB()
         fni = self.FNI()    
         return 0.11*math.exp(-3*(pb**-2))*(fni**3)*\
                        (self.abt**1.5)*rho*9.81/(1+ (fni**2))
                                              
    def FNT(self):
         return self.v/math.sqrt(2*9.81*self.at/(self.b+self.b*self.cwp))
     
    def C6(self):
         fnt = self.FNT()
         if fnt >= 5:
             return 0
         else:
             return 0.2*(1-0.2*fnt)
     
    def RTR(self):
        rho = 1.025
        c6 = self.C6()
        return 0.5*rho*(self.v**2)*self.at*c6
     
    def C4(self):
         if self.tf/self.l <= 0.04:
             return self.tf/self.l
         else:
             return 0.04
         
    def CA(self):    
         cb = self.CB()
         c2 = self.C2()
         c4 = self.C4()
         return 0.006*((self.l+100)**-0.16)- 0.00205 + 0.003*\
                math.sqrt(self.l/7.5)*(cb**4)*c2*(0.04-c4)
      
    def RA(self):
        rho = 1.025
        ca = self.CA()
        s = self.wetarea()
        return 0.5*rho*(self.v**2)*s*ca
     
    def RTOT(self):
        rho = 1.025
        rf = self.RF()
        k1 = self.FF()
        rapp = self.RAPP()
        rw = self.RW()
        rb = self.RB()
        rtr = self.RTR()
        ra = self.RA() 
        return rf*k1 + rapp + rw + rb + rtr + ra
   
    def Calculate(self):
           Rt = self.RTOT()
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
     res1.cstern = 10
     res1.v = 12.86

     res1.Calculate()
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
   
   

    
    
    