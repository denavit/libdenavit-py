import dataclasses
from math import sqrt,atan,radians,tan,pi
from ..design import available_strength
from . import database


class WideFlangeDB:

    def __init__(self,name):
        self.name = name
        self.data = database.wide_flange_database[name.upper()]

    @property
    def d(self):
        return self.data['d']

    @property
    def bf(self):
        return self.data['bf']

    @property
    def tf(self):
        return self.data['tf']

    @property
    def tw(self):
        return self.data['tw']

    @property
    def A(self):
        return self.data['A']

    @property
    def Ix(self):
        return self.data['Ix']

    @property
    def Zx(self):
        return self.data['Zx']

    @property
    def Sx(self):
        return self.data['Sx']

    @property
    def rx(self):
        return self.data['rx']

    @property
    def Iy(self):
        return self.data['Iy']

    @property
    def Zy(self):
        return self.data['Zy']
        
    @property
    def Sy(self):
        return self.data['Sy']

    @property
    def ry(self):
        return self.data['ry']

    @property
    def J(self):
        return self.data['J']
        
    @property
    def Cw(self):
        return self.data['Cw']
    
    @property
    def rts(self):
        return self.data['rts']
    
    @property
    def ho(self):
        return self.data['ho']

    @property
    def h_over_tw(self):
        return self.data['h/tw']
        
    @property
    def bf_over_2tf(self):
        return self.data['bf/2tf']
        
 
class WideFlangeMember_AISC2016:
    def __init__(self,section,Fy,E,G,strength_type):
        self.section = section
        self.Fy = Fy
        self.E = E
        self.G = G
        self.strength_type = strength_type

    def Ae(self,Fcr):
        Ae = self.section.A
    
        # Effective flange width
        lam = self.section.bf_over_2tf
        lam_r = 0.56*sqrt(self.E/self.Fy)
        if lam > lam_r*sqrt(self.Fy/Fcr):
            c1 = 0.22
            c2 = 1.49
            Fel = (c2*lam_r/lam)**2*self.Fy
            be = self.section.bf*(1-c1*sqrt(Fel/Fcr))*sqrt(Fel/Fcr)
            Ae = Ae - 2*(self.section.bf-be)*self.section.tf
        
        # Effective web height
        lam = self.section.h_over_tw
        lam_r = 1.49*sqrt(self.E/self.Fy)
        if lam > lam_r*sqrt(self.Fy/Fcr):
            c1 = 0.18
            c2 = 1.31
            Fel = (c2*lam_r/lam)**2*self.Fy
            h = lam*self.section.tw
            he = h*(1-c1*sqrt(Fel/Fcr))*sqrt(Fel/Fcr)
            Ae = Ae - (h-he)*self.section.tw
        
        return Ae
           
    def Pnt(self):
        Pn = self.Fy*self.section.A
        return  available_strength(Pn,self.strength_type,0.9,1.67)

    def Pnc(self,Lcx,Lcy):
        # Does not compute torsional buckling
        
        # Major-axis flexural buckling
        if Lcx == 0.0:
            Fcrx = self.Fy
        else:
            Fe = pi**2*self.E/(Lcx/self.section.rx)**2
            if self.Fy/Fe <= 2.25:
                Fcrx = 0.658**(self.Fy/Fe)*self.Fy
            else:
                Fcrx = 0.877*Fe            
        
        # Minor-axis flexural buckling
        if Lcy == 0.0:
            Fcry = self.Fy
        else:
            Fe = pi**2*self.E/(Lcy/self.section.ry)**2
            if self.Fy/Fe <= 2.25:
                Fcry = 0.658**(self.Fy/Fe)*self.Fy
            else:
                Fcry = 0.877*Fe  
                
        # Compute strength
        Fcr = min(Fcrx,Fcry)
        Pn = Fcr*self.Ae(Fcr)
        
        return available_strength(Pn,self.strength_type,0.9,1.67)
      
    def Mn(self,Lb,Cb):
        
        # Check width-to-thickness ratios
        if self.section.bf_over_2tf > 0.38*sqrt(self.E/self.Fy):
            raise Exception('Mn not yet implemented for noncompact or slender flanges') 
        if self.section.h_over_tw > 3.76*sqrt(self.E/self.Fy) :
            raise Exception('Mn not yet implemented for noncompact or slender webs')           
        
        # Yielding 
        Mn_Y = self.Fy*self.section.Zx
        
        # Lateral-Torsional Buckling
        Lp = 1.76*self.section.ry*sqrt(self.E/self.Fy)
        Lr = 1.95*self.section.rts*(self.E/(0.7*self.Fy))*sqrt((self.section.J/(self.section.Sx*self.section.ho))+sqrt(((self.section.J/(self.section.Sx*self.section.ho))**2)+(6.76*((0.7*self.Fy)/self.E)**2)))
        
        if Lb <= Lp:
            Mn = Mn_Y
        if Lb <= Lr:
            Mn_LTB = Cb*((self.Fy*self.section.Zx)-(((self.Fy*self.section.Zx)-(0.7*self.Fy*self.section.Sx))*((Lb-Lp)/(Lr-Lp))))
            Mn = min(Mn_LTB,Mn_Y)
        else:
            Fcr= Cb*(pi**2)*((self.E)/((Lb/self.section.rts)**2))*sqrt(1+(0.078*((self.section.J/(self.section.Sx*self.section.ho))*((Lb/self.section.rts)**2))))
            Mn_LTB = Fcr*self.section.Sx
            Mn = min(Mn_LTB,Mn_Y)
                
        return available_strength(Mn,self.strength_type,0.9,1.67)
    
