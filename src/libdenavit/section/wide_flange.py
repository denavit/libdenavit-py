import dataclasses
import warnings
from math import sqrt,pi,ceil
from libdenavit import opensees as ops
from libdenavit.design import available_strength
from libdenavit.section import GeometricShape,database

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
        



class I_shape(GeometricShape):
        
    def __init__(self, d, tw, bf, tf, Fy, E,
                 A=None, Ix=None, Zx=None, Sx=None, rx=None,
                 Iy=None, Zy=None, Sy=None, ry=None,
                 J=None, Cw=None, rts=None, ho=None, h_over_tw=None):
        self.d = d
        self.tw = tw
        self.bf = bf
        self.tf = tf

        self.Fy = Fy
        self.E = E

        self._A = A
        self._Ix = Ix
        self._Zx = Zx
        self._Sx = Sx
        self._rx = rx
        self._Iy = Iy
        self._Zy = Zy
        self._Sy = Sy
        self._ry = ry
        self._J = J
        self._Cw = Cw
        self._rts = rts
        self._ho = ho
        self._h_over_tw = h_over_tw

        self.has_steel=True
        self.has_concrete=False

    @classmethod
    def from_database(cls, section_name,Fy,E):

        db = WideFlangeDB(section_name)
        return cls(
            d=db.d,
            tw=db.tw,
            bf=db.bf,
            tf=db.tf,
            Fy=Fy,
            E=E,
            A=db.A,
            Ix=db.Ix,
            Zx=db.Zx,
            Sx=db.Sx,
            rx=db.rx,
            Iy=db.Iy,
            Zy=db.Zy,
            Sy=db.Sy,
            ry=db.ry,
            J=db.J,
            Cw=db.Cw,
            rts=db.rts,
            ho=db.ho,
            h_over_tw=db.h_over_tw,
        )


    @property
    def A(self):
        if self._A is not None:
            return self._A
        else:
            return 2 * self.bf * self.tf + (self.d - 2 * self.tf) * self.tw

    @A.setter
    def A(self, x):
        self._A = x

    @property
    def Ix(self):
        if self._Ix is not None:
            return self._Ix
        else:
            flange = (self.bf * self.tf**3) / 6
            web = (self.tw * (self.d - 2 * self.tf)**3) / 12
            return 2 * flange + web

    @Ix.setter
    def Ix(self, x):
        self._Ix = x

    @property
    def Zx(self):
        if self._Zx is not None:
            return self._Zx
        else:
            raise ValueError("Zx formula not set")

    @Zx.setter
    def Zx(self, x):
        self._Zx = x

    @property
    def Sx(self):
        if self._Sx is not None:
            return self._Sx
        else:
            raise ValueError("Sx formula not set")

    @Sx.setter
    def Sx(self, x):
        self._Sx = x

    @property
    def rx(self):
        if self._rx is not None:
            return self._rx
        else:
            raise ValueError("rx formula not set")


    @rx.setter
    def rx(self, x):
        self._rx = x
        
        
    @property
    def Iy(self):
        if self._Iy is not None:
            return self._Iy
        else:
            flange = (self.tf * self.bf**3) / 12
            web = ((self.d - 2 * self.tf) * self.tw**3) / 12
            return 2 * flange + web

    @Iy.setter
    def Iy(self, x):
        self._Iy = x

    @property
    def Zy(self):
        if self._Zy is not None:
            return self._Zy
        else:
            raise ValueError("Zy formula not set")


    @Zy.setter
    def Zy(self, x):
        self._Zy = x

    @property
    def Sy(self):
        if self._Sy is not None:
            return self._Sy
        else:
            raise ValueError("Sy formula not set")


    @Sy.setter
    def Sy(self, x):
        self._Sy = x
        

    @property
    def ry(self):
        if self._ry is not None:
            return self._ry
        else:
            raise ValueError("ry formula not set")

    @ry.setter
    def ry(self, x):
        self._ry = x

    @property
    def J(self):
        if self._J is not None:
            return self._J
        else:
            raise ValueError("J formula not set")

    @J.setter
    def J(self, x):
        self._J = x

    @property
    def Cw(self):
        if self._Cw is not None:
            return self._Cw
        else:
            raise ValueError("Cw formula not set")


    @Cw.setter
    def Cw(self, x):
        self._Cw = x

    @property
    def rts(self):
        if self._rts is not None:
            return self._rts
        else:
            raise ValueError("rts formula not set")


    @rts.setter
    def rts(self, x):
        self._rts = x

    @property
    def ho(self):
        if self._ho is not None:
            return self._ho
        else:
            raise ValueError("ho formula not set")


    @ho.setter
    def ho(self, x):
        self._ho = x

    @property
    def bf_over_2tf(self):
        return self.bf / (2 * self.tf)
    
    @property
    def h_over_tw(self):
        if self._h_over_tw is not None:
            return self._h_over_tw
        else:
            return self.dw /self.tw

    @h_over_tw.setter
    def h_over_tw(self, x):
        self._h_over_tw = x

    @property
    def dw(self):
        return self.d - 2 * self.tf


    @property
    def p0(self):
        # Nominal axial yield strength (short steel section)
        return self.Fy * self.A
    
    def depth(self, axis):
        if axis=='x' or axis==None:
            return self.d
        else:
            return self.bf
        
    def width(self, axis):
        if axis=='x' or axis==None:
            return self.bf
        else:
            return self.d
    

    def build_ops_fiber_section(self, section_id, start_material_id=1, mat_type=None, 
        nfy=20, nfx=20, frc=0, num_regions=10, 
        stiffness_reduction=1.0,strength_reduction=1.0, 
        GJ=1.0e6, axis=None, 
        hardening_ratio=0.001, MultiLinearPts=None):

        E_reduced=self.E*stiffness_reduction
        Fy_reduced=self.Fy*strength_reduction
        
        ## Define Base Material
        base_material_id = start_material_id
        if mat_type == 'Elastic':
            ops.uniaxialMaterial('Elastic', base_material_id, E_reduced)
        elif mat_type == 'ElasticPP':
            ops.uniaxialMaterial('ElasticPP', base_material_id, E_reduced, Fy_reduced/E_reduced)
        elif mat_type == 'Steel01':
            b = hardening_ratio / (1 + hardening_ratio)
            ops.uniaxialMaterial('Steel01', base_material_id, Fy_reduced, E_reduced, b)
        elif mat_type == 'Hardening':
            ops.uniaxialMaterial('Hardening', base_material_id, E_reduced, Fy_reduced, 0.0, hardening_ratio*E_reduced)
        elif mat_type == 'MultiLinear':
            if MultiLinearPts == None:
                Fy_reduced = Fy_reduced
                ey = Fy_reduced/(E_reduced)
                MultiLinearPts = [
                    0.5000*ey,0.5*Fy_reduced,
                    0.7222*ey,0.7*Fy_reduced,
                    0.8556*ey,0.8*Fy_reduced,
                    1.0556*ey,0.9*Fy_reduced,
                    1.7222*ey,Fy_reduced,                            
                    101.7222*ey,1.1*Fy_reduced,
                ]
            ops.uniaxialMaterial('MultiLinear', base_material_id, *MultiLinearPts)
        else:
            raise Exception(
                'Input Error - unknown material type (%s)' % mat_type)

        ops.section('Fiber', section_id, '-GJ', GJ)

        if axis=='x':
            Nfw = ceil(self.dw * (nfy / self.d))
            Nff = ceil(self.tf * (nfy / self.d))

            if frc == 0 or mat_type == 'Elastic':
                # Web patch
                ops.patch('rect', base_material_id, Nfw, 1, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)
                # Flange patches
                ops.patch('rect', base_material_id, Nff, 1, self.dw / 2, -self.bf / 2, self.d / 2, self.bf / 2)
                ops.patch('rect', base_material_id, Nff, 1, -self.d / 2, -self.bf / 2, -self.dw / 2, self.bf / 2)

            else:   
                frt = -frc * (self.bf * self.tf) / (self.bf * self.tf + self.tw * self.dw)
                material_idi = start_material_id+1
                ops.uniaxialMaterial('InitStressMaterial',material_idi, base_material_id, frt)
                # Web patch
                ops.patch('rect', material_idi, Nfw, 1, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)

                region_width = self.bf / num_regions
                for i in range(num_regions):
                    fri = frc + ((i + 0.5) / num_regions) * (frt - frc)
                    material_idi = start_material_id + 2 + i 
                    ops.uniaxialMaterial('InitStressMaterial', material_idi, base_material_id, fri)
                    # Top flange
                    ops.patch('rect', material_idi, Nff, 1, self.dw / 2, -region_width / 2, self.d / 2, region_width / 2)
                    # Bottom flange
                    ops.patch('rect', material_idi, Nff, 1, -self.d / 2, -region_width / 2, -self.dw / 2, region_width / 2)

        elif axis=='y':
            Nfw = ceil(self.tw * (nfx / self.bf))
            Nff = ceil(self.bf * (nfx / self.bf))

            if frc == 0 or mat_type == 'Elastic':
                # Web patch
                ops.patch('rect', base_material_id, Nfw, 1, -self.tw/2, -self.dw/2, self.tw/2, self.dw/2)
                # Flange patches
                ops.patch('rect', base_material_id, Nff, 1, -self.bf/2, -self.d/2, self.bf/2, -self.dw/2)
                ops.patch('rect', base_material_id, Nff, 1, -self.bf/2, self.dw/2, self.bf/2, self.d/2)

            else:
                frt = -frc * (self.bf * self.tf) / (self.bf * self.tf + self.tw * self.dw)
                material_idi = start_material_id+1
                ops.uniaxialMaterial('InitStressMaterial',material_idi, base_material_id, frt)
                # Web patch
                ops.patch('rect', material_idi, Nfw, 1, -self.tw/2, -self.dw/2, self.tw/2, self.dw/2)

                half_flange_width = self.bf / 2.0
                region_width = half_flange_width / num_regions
                Nff_region = ceil(region_width * (nfx / self.bf))
                for i in range(num_regions):
                    fri = frc + ((i + 0.5) / num_regions) * (frt - frc)
                    material_idi = start_material_id + 2 + i 
                    ops.uniaxialMaterial('InitStressMaterial', material_idi, base_material_id, fri)

                    y_start_right = half_flange_width - (i + 1) * region_width
                    y_end_right   = half_flange_width - i * region_width
                    y_start_left = -y_end_right
                    y_end_left   = -y_start_right
                    # Flange patches
                    ops.patch('rect', material_idi, Nff_region, 1, y_start_right, self.dw/2, y_end_right, self.d/2)
                    ops.patch('rect', material_idi, Nff_region, 1, y_start_left,  self.dw/2, y_end_left,  self.d/2)
                    ops.patch('rect', material_idi, Nff_region, 1, y_start_right, -self.d/2, y_end_right, -self.dw/2)
                    ops.patch('rect', material_idi, Nff_region, 1, y_start_left,  -self.d/2, y_end_left,  -self.dw/2)
                    
        elif axis is None:
            Nfw_x = ceil(self.tw * (nfx / self.bf))
            Nff_x = ceil(self.bf * (nfx / self.bf))
            Nfw_y = ceil(self.dw * (nfy / self.d))
            Nff_y = ceil(self.tf * (nfy / self.d))

            if frc == 0 or mat_type == 'Elastic':
                # Web patch
                ops.patch('rect', base_material_id, Nfw_y, Nfw_x, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)
                # Flange patches
                ops.patch('rect', base_material_id, Nff_y, Nff_x, self.dw / 2, -self.bf / 2, self.d / 2, self.bf / 2)
                ops.patch('rect', base_material_id, Nff_y, Nff_x, -self.d / 2, -self.bf / 2, -self.dw / 2, self.bf / 2)

            else:
                frt = -frc * (self.bf * self.tf) / (self.bf * self.tf + self.tw * self.dw)
                material_idi = start_material_id+1
                ops.uniaxialMaterial('InitStressMaterial',material_idi, base_material_id, frt)
                # Web patch
                ops.patch('rect', material_idi, Nfw_y, Nfw_x, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)

                region_width = self.bf / num_regions
                Nfw_x = ceil(self.tw * (nfx / region_width))
                Nff_x = ceil(region_width * (nfx / region_width))
                for i in range(num_regions):

                    fri = frc + ((i + 0.5) / num_regions) * (frt - frc)
                    material_idi = start_material_id + 2 + i 
                    ops.uniaxialMaterial('InitStressMaterial', material_idi, base_material_id, fri)
                    # Top flange
                    ops.patch('rect', material_idi, Nff_y, Nff_x, self.dw / 2, -region_width / 2, self.d / 2, region_width / 2)
                    # Bottom flange
                    ops.patch('rect', material_idi, Nff_y, Nff_x, -self.d / 2, -region_width / 2, -self.dw / 2, region_width / 2)

        else:
            raise ValueError("Please give valid axis, 'x' or 'y' or set it to None to use a 3d fiber section")



    def maximum_compression_strain(self, axial_strain, curvatureX=0, curvatureY=0):
        extreme_strain = axial_strain - self.d/2 * abs(curvatureX) \
                            - self.bf/2 * abs(curvatureY)
        return extreme_strain



    def maximum_tensile_strain(self, axial_strain, curvatureX=0.0, curvatureY=0.0):
        extreme_strain = axial_strain + self.d/2 * abs(curvatureX) \
                                    + self.bf/2 * abs(curvatureY)
        return extreme_strain

class WideFlangeMember_AISC2022:
    def __init__(self,section,Fy,E,L,**kwargs):
        self.section = section
        self.Fy = Fy
        self.E = E
        self.L=L

        defaults={'strength_type':'lrfd'
                  }
        
        for key,value in defaults.items():
            setattr(self,key,kwargs.get(key,value))

    def check_compactness_of_flange(self):
        if self.section.bf_over_2tf> 0.38*sqrt(self.E/self.Fy) :
            raise Exception('Mn not yet implemented for noncompact or slender flanges')
    def check_compactness_of_web(self):
        if self.section.h_over_tw > 3.76*sqrt(self.E/self.Fy) :
            raise Exception('Mn not yet implemented for noncompact or slender webs')
    def check_slenderness_of_flange(self):
        limit = 0.56 * sqrt(self.E / self.Fy)
        return self.section.bf_over_2tf > limit
    def check_slenderness_of_web(self):
        limit = 1.49 *sqrt(self.E / self.Fy)
        return self.section.h_over_tw > limit  
          
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

    def Mnx(self,Lb,Cb):
        """Moment strength of member for major-axis bending.

        Parameters
        ----------
        Lb : float
            Unbraced length of the member.
        Cb : float
            Lateral-torsional buckling modification factor.

        Notes
        -----
        - Not yet implemented for noncompact or slender webs.

        Reference: AISC Specification Chapter F; Sections F1 -- F3
        """

        # Check width-to-thickness ratios
        self.check_compactness_of_web()           
               
        # Yielding 
        Mp = self.Fy*self.section.Zx
        Mn = Mp

        # Compression Flange Local Buckling
        λ  = self.section.bf_over_2tf
        λp = 0.38*sqrt(self.E/self.Fy)
        λr = 1.0*sqrt(self.E/self.Fy)
        if λ <= λp:
            pass
        elif λ <= λr:
            Mn_CFLB = Mp - (Mp-0.7*self.Fy*self.section.Sx)*(λ-λp)/(λr-λp)
            Mn = min(Mn,Mn_CFLB)
        else:
            kc = 4/sqrt(self.section.h_over_tw)
            if kc > 0.76:
                kc = 0.76
            if kc < 0.35:
                kc = 0.35
            Mn_CFLB = (0.9*self.E*kc*self.section.Sx)/(λ**2)
            Mn = min(Mn,Mn_CFLB)

        # Lateral-Torsional Buckling
        Lp = 1.76*self.section.ry*sqrt(self.E/self.Fy)
        Lr = 1.95*self.section.rts*(self.E/(0.7*self.Fy))*sqrt((self.section.J/(self.section.Sx*self.section.ho))+sqrt(((self.section.J/(self.section.Sx*self.section.ho))**2)+(6.76*((0.7*self.Fy)/self.E)**2)))
        
        if Lb <= Lp:
            pass
        if Lb <= Lr:
            Mn_LTB = Cb*((self.Fy*self.section.Zx)-(((self.Fy*self.section.Zx)-(0.7*self.Fy*self.section.Sx))*((Lb-Lp)/(Lr-Lp))))
            Mn = min(Mn,Mn_LTB)
        else:
            Fcr= Cb*(pi**2)*((self.E)/((Lb/self.section.rts)**2))*sqrt(1+(0.078*((self.section.J/(self.section.Sx*self.section.ho))*((Lb/self.section.rts)**2))))
            Mn_LTB = Fcr*self.section.Sx
            Mn = min(Mn,Mn_LTB)
                
        return available_strength(Mn,self.strength_type,0.9,1.67)

    def Mny(self):
        """Moment strength of member for minor-axis bending.

        Reference: AISC Specification Sections F1, F6
        """
        # Yielding
        Mp = min(self.Fy*self.section.Zy, 1.6*self.Fy*self.section.Sy)
        Mn = Mp

        # Flange local buckling
        λ  = self.section.bf_over_2tf
        λp = 0.38*sqrt(self.E/self.Fy)
        λr = 1.0*sqrt(self.E/self.Fy)

        if λ <= λp:
            pass
        elif λ <= λr:
            Mn_CFLB = Mp - (Mp-0.7*self.Fy*self.section.Sy)*(λ-λp)/(λr-λp)
            Mn = min(Mn,Mn_CFLB)
        else:
            Fcr = 0.7*self.E/λ**2
            Mn_CFLB = Fcr*self.section.Sy
            Mn = min(Mn,Mn_CFLB)
        
        return available_strength(Mn,self.strength_type,0.9,1.67)

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
        slender=self.check_slenderness_of_flange() or self.check_slenderness_of_web()
        if not slender:
            A=self.section.A
        else:
            A=self.Ae(Fcr)
        Pn = Fcr*A
        
        return available_strength(Pn,self.strength_type,0.9,1.67)



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
        return available_strength(Pn,self.strength_type,0.9,1.67)

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
        warnings.warn("Mn() is deprecated; use Mnx() instead.", stacklevel=2)
        return self.Mnx(Lb,Cb)

    def Mnx(self,Lb,Cb):
        """Moment strength of member for major-axis bending.

        Parameters
        ----------
        Lb : float
            Unbraced length of the member.
        Cb : float
            Lateral-torsional buckling modification factor.

        Notes
        -----
        - Not yet implemented for noncompact or slender webs.

        Reference: AISC Specification Chapter F; Sections F1 -- F3
        """

        # Check width-to-thickness ratios
        if self.section.h_over_tw > 3.76*sqrt(self.E/self.Fy) :
            raise Exception('Mn not yet implemented for noncompact or slender webs')           
               
        # Yielding 
        Mp = self.Fy*self.section.Zx
        Mn = Mp

        # Compression Flange Local Buckling
        λ  = self.section.bf_over_2tf
        λp = 0.38*sqrt(self.E/self.Fy)
        λr = 1.0*sqrt(self.E/self.Fy)
        if λ <= λp:
            pass
        elif λ <= λr:
            Mn_CFLB = Mp - (Mp-0.7*self.Fy*self.section.Sx)*(λ-λp)/(λr-λp)
            Mn = min(Mn,Mn_CFLB)
        else:
            kc = 4/sqrt(self.section.h_over_tw)
            if kc > 0.76:
                kc = 0.76
            if kc < 0.35:
                kc = 0.35
            Mn_CFLB = (0.9*self.E*kc*self.section.Sx)/(λ**2)
            Mn = min(Mn,Mn_CFLB)

        # Lateral-Torsional Buckling
        Lp = 1.76*self.section.ry*sqrt(self.E/self.Fy)
        Lr = 1.95*self.section.rts*(self.E/(0.7*self.Fy))*sqrt((self.section.J/(self.section.Sx*self.section.ho))+sqrt(((self.section.J/(self.section.Sx*self.section.ho))**2)+(6.76*((0.7*self.Fy)/self.E)**2)))
        
        if Lb <= Lp:
            pass
        if Lb <= Lr:
            Mn_LTB = Cb*((self.Fy*self.section.Zx)-(((self.Fy*self.section.Zx)-(0.7*self.Fy*self.section.Sx))*((Lb-Lp)/(Lr-Lp))))
            Mn = min(Mn,Mn_LTB)
        else:
            Fcr= Cb*(pi**2)*((self.E)/((Lb/self.section.rts)**2))*sqrt(1+(0.078*((self.section.J/(self.section.Sx*self.section.ho))*((Lb/self.section.rts)**2))))
            Mn_LTB = Fcr*self.section.Sx
            Mn = min(Mn,Mn_LTB)
                
        return available_strength(Mn,self.strength_type,0.9,1.67)
    
    def Mny(self):
        """Moment strength of member for minor-axis bending.

        Reference: AISC Specification Sections F1, F6
        """
        # Yielding
        Mp = min(self.Fy*self.section.Zy, 1.6*self.Fy*self.section.Sy)
        Mn = Mp

        # Flange local buckling
        λ  = self.section.bf_over_2tf
        λp = 0.38*sqrt(self.E/self.Fy)
        λr = 1.0*sqrt(self.E/self.Fy)

        if λ <= λp:
            pass
        elif λ <= λr:
            Mn_CFLB = Mp - (Mp-0.7*self.Fy*self.section.Sy)*(λ-λp)/(λr-λp)
            Mn = min(Mn,Mn_CFLB)
        else:
            Fcr = 0.69*self.E/λ**2 ## Has this changed to 0.7 in AISC 2022
            Mn_CFLB = Fcr*self.section.Sy
            Mn = min(Mn,Mn_CFLB)
        
        return available_strength(Mn,self.strength_type,0.9,1.67)

    def Vn(self):
        Aw = self.section.d*self.section.tw
        kv = 5.34 # For webs without transverse stiffeners
        if self.section.h_over_tw <= 2.24*sqrt(self.E/self.Fy):
            # Cv1 = 1.0
            Vn = 0.6*self.Fy*Aw
            return available_strength(Vn,self.strength_type,1.00,1.50)
        elif self.section.h_over_tw <= 1.10*sqrt(kv*self.E/self.Fy):
            # Cv1 = 1.0
            Vn = 0.6*self.Fy*Aw
            return available_strength(Vn,self.strength_type,0.90,1.67)
        else:
            raise Exception('Vn for shear buckling not yet implemented')
