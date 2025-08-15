import libdenavit.section.database.aisc as section
import libdenavit.section.wide_flange as database
from libdenavit.section.geometric_shape import *
from libdenavit.OpenSees.get_fiber_data import *
from math import pi, ceil
import openseespy.opensees as ops
from Units import *
import matplotlib.pyplot as plt
import opsvis as opsv

class WF_Database:
    def __init__(self,Section_name,unit=inch):

        self.section=Section_name
        self.d=section.wide_flange_database[self.section]['d']*unit
        self.tw=section.wide_flange_database[self.section]['tw']*unit
        self.bf=section.wide_flange_database[self.section]['bf']*unit
        self.tf=section.wide_flange_database[self.section]['tf']*unit
        self.A=section.wide_flange_database[self.section]['A']*(unit**2)
        self.Ix=section.wide_flange_database[self.section]['Ix']*(unit**4)
        self.Iy=section.wide_flange_database[self.section]['Iy']*(unit**4)

class wf_Database:
        def __init__(self,Section_name,unit=inch):
        
            db = database.WideFlangeDB(Section_name)
            self.d   = db.d   * unit
            self.tw  = db.tw  * unit
            self.bf  = db.bf  * unit
            self.tf  = db.tf  * unit
            self.A   = db.A   * unit**2
            self.Ix  = db.Ix  * unit**4
            self.Zx  = db.Zx  * unit**3
            self.Sx  = db.Sx  * unit**3
            self.rx  = db.rx  * unit
            self.Iy  = db.Iy  * unit**4
            self.Zy  = db.Zy  * unit**3
            self.Sy  = db.Sy  * unit**3
            self.ry  = db.ry  * unit
            self.J   = db.J   * unit**4
            self.Cw  = db.Cw  * unit**6
            self.rts = db.rts * unit
            self.ho  = db.ho  * unit


class I_shape(GeometricShape):
        
    def __init__(self, d, tw, bf, tf, fy, E, Hk,
                 A=None, Ix=None, Zx=None, Sx=None, rx=None,
                 Iy=None, Zy=None, Sy=None, ry=None,
                 J=None, Cw=None, rts=None, ho=None):
        self.d = d
        self.tw = tw
        self.bf = bf
        self.tf = tf

        self.fy = fy
        self.E = E
        self.Hk = Hk
        self.b=Hk / (E + Hk)

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


        self.num_regions = 10   # Number of regions for the discretization of residual stress
        self.num_elements = 20
        # self.num_fiber = 20
        self.num_steps = 100


    @classmethod
    def from_database(cls, section_name,fy,E,Hk):
        db = database.WideFlangeDB(section_name)
        return cls(
            d=db.d,
            tw=db.tw,
            bf=db.bf,
            tf=db.tf,
            fy=fy,
            E=E,
            Hk=Hk,
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
            ho=db.ho
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

    @property
    def Sx(self):
        if self._Sx is not None:
            return self._Sx
        else:
            raise ValueError("Sx formula not set")

    @property
    def rx(self):
        if self._rx is not None:
            return self._rx
        else:
            raise ValueError("rx formula not set")

    
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

    @property
    def Sy(self):
        if self._Sy is not None:
            return self._Sy
        else:
            raise ValueError("Sy formula not set")

    @property
    def ry(self):
        if self._ry is not None:
            return self._ry
        else:
            raise ValueError("ry formula not set")

    @property
    def J(self):
        if self._J is not None:
            return self._J
        else:
            raise ValueError("J formula not set")

    @property
    def Cw(self):
        if self._Cw is not None:
            return self._Cw
        else:
            raise ValueError("Cw formula not set")

    @property
    def rts(self):
        if self._rts is not None:
            return self._rts
        else:
            raise ValueError("rts formula not set")

    @property
    def ho(self):
        if self._ho is not None:
            return self._ho
        else:
            raise ValueError("ho formula not set")


    @property
    def bf_over_2tf(self):
        if hasattr(self, "bf") and hasattr(self, "tf") and self.tf != 0:
            return self.bf / (2 * self.tf)
        else:
            raise ValueError("bf_over_2tf formula not set or invalid (tf = 0)")

    @property
    def dw(self):
        if hasattr(self, "d") and hasattr(self, "tf"):
            return self.d - 2 * self.tf
        else:
            raise ValueError("dw formula not set")

    @property
    def p0(self):
        # Nominal axial yield strength (short steel section)
        return self.fy * self.A
    

    def build_ops_fiber_section(self, section_id, start_material_id, mat_type, nfy, nfx, frc, GJ=1.0e6,axis=None):
                                
        self.num_fiber=max(nfx,nfy)    ### used only for 2d fiber sections

        self.frc = frc



        if axis=='x':
            Nfw = ceil(self.dw * (self.num_fiber / self.d))
            Nff = ceil(self.tf * (self.num_fiber / self.d))

            # if not Residual_Stress:

            if self.frc == 0 or mat_type == 'Elastic':
                if mat_type == 'Elastic':
                    ops.uniaxialMaterial('Elastic', start_material_id, self.E)
                elif mat_type == 'ElasticPP':
                    ops.uniaxialMaterial('ElasticPP', start_material_id,
                                        self.E, self.fy/self.E)
                elif mat_type == 'Steel01':
                    ops.uniaxialMaterial('Steel01', start_material_id, self.fy, self.E, self.b)
                elif mat_type == 'Hardening':
                    ops.uniaxialMaterial('Hardening', start_material_id,
                                        self.E, self.fy, 0.0, self.Hk)
                else:
                    raise Exception(
                        'Input Error - unknown material type (%s)' % mat_type)
                
                ops.section('Fiber', section_id, '-GJ', GJ)

                # Web patch
                ops.patch('rect', start_material_id, Nfw, 1, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)
                # Flange patches
                ops.patch('rect', start_material_id, Nff, 1, self.dw / 2, -self.bf / 2, self.d / 2, self.bf / 2)
                ops.patch('rect', start_material_id, Nff, 1, -self.d / 2, -self.bf / 2, -self.dw / 2, self.bf / 2)

                # fib_sec = [

                #         ['section', 'Fiber', start_material_id, '-GJ', 1.0e4],

                #         ["patch","rect", start_material_id, Nfw,1,-self.dw/2, -self.tw/2, self.dw/2, self.tw/2],
                #         ["patch","rect",start_material_id, Nff,1,self.dw/2, -self.bf/2,self.d/2, self.bf/2],
                #         ["patch","rect", start_material_id, Nff,1,-self.d/2, -self.bf/2, -self.dw/2, self.bf/2]

                #         ]
                
                # matcolor = [ 'gold']*10000
                # opsv.plot_fiber_section(fib_sec, matcolor=matcolor)
                # plt.axis('equal')
                # plt.title('section_name')
                # plt.show()
            
            
            else:
                ops.section('Fiber', section_id, '-GJ', GJ)
                

                frt = -self.frc * (self.bf * self.tf) / (self.bf * self.tf + self.tw * self.dw)

                ## web patch
                if mat_type == 'ElasticPP':
                    ops.uniaxialMaterial(
                        'ElasticPP', start_material_id, self.E, self.fy/self.E, -self.fy/self.E, frt/self.E)
                elif mat_type == 'Steel01':
                    ops.uniaxialMaterial('Steel01', start_material_id+1, self.fy, self.E, self.b)
                    ops.uniaxialMaterial('InitStressMaterial',
                                        start_material_id, start_material_id+1, frt)
                elif mat_type == 'Hardening':
                    ops.uniaxialMaterial('Hardening', start_material_id+1,
                                        self.E, self.fy, 0.0, self.Hk)
                    ops.uniaxialMaterial('InitStressMaterial',
                                        start_material_id, start_material_id+1, frt)
                else:
                    raise Exception(
                        'Input Error - unknown material type (%s)' % mat_type)


                ops.patch('rect', start_material_id, Nfw, 1, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)



                region_width = self.bf / self.num_regions

                for i in range(self.num_regions):
                    fri = self.frc + ((i + 0.5) / self.num_regions) * (frt - self.frc)
                    start_material_idi = start_material_id + 2 * (i + 1)

                    if mat_type == 'ElasticPP':
                        ops.uniaxialMaterial(
                            'ElasticPP', start_material_idi, self.E, self.fy/self.E, -self.fy/self.E, fri/self.E)
                    elif mat_type == 'Steel01':
                        ops.uniaxialMaterial(
                            'Steel01', start_material_idi+1, self.fy, self.E, self.b)
                        ops.uniaxialMaterial(
                            'InitStressMaterial', start_material_idi, start_material_idi+1, fri)
                    elif mat_type == 'Hardening':
                        ops.uniaxialMaterial(
                            'Hardening', start_material_idi+1, self.E, self.fy, 0.0, self.Hk)
                        ops.uniaxialMaterial(
                            'InitStressMaterial', start_material_idi, start_material_idi+1, fri)
                    else:
                        raise Exception(
                            'Input Error - unknown material type (%s)' % mat_type)

                    # Top flange
                    ops.patch('rect', start_material_idi, Nff, 1, self.dw / 2, -region_width / 2, self.d / 2, region_width / 2)
                    # Bottom flange
                    ops.patch('rect', start_material_idi, Nff, 1, -self.d / 2, -region_width / 2, -self.dw / 2, region_width / 2)

        



        elif axis=='y':
            Nfw = ceil(self.tw * (self.num_fiber / self.bf))
            Nff = ceil(self.bf * (self.num_fiber / self.bf))

            # if not Residual_Stress:

            if self.frc == 0 or mat_type == 'Elastic':
                if mat_type == 'Elastic':
                    ops.uniaxialMaterial('Elastic', start_material_id, self.E)
                elif mat_type == 'ElasticPP':
                    ops.uniaxialMaterial('ElasticPP', start_material_id,
                                        self.E, self.fy/self.E)
                elif mat_type == 'Steel01':
                    ops.uniaxialMaterial('Steel01', start_material_id, self.fy, self.E, self.b)
                elif mat_type == 'Hardening':
                    ops.uniaxialMaterial('Hardening', start_material_id,
                                        self.E, self.fy, 0.0, self.Hk)
                else:
                    raise Exception(
                        'Input Error - unknown material type (%s)' % mat_type)
                
                ops.section('Fiber', section_id, '-GJ', GJ)


                # Web patch
                ops.patch('rect', start_material_id, Nfw, 1, -self.tw/2, -self.dw/2, self.tw/2, self.dw/2)

                # Flange patches
                ops.patch('rect', start_material_id, Nff, 1, -self.bf/2, -self.d/2, self.bf/2, -self.dw/2)
                ops.patch('rect', start_material_id, Nff, 1, -self.bf/2, self.dw/2, self.bf/2, self.d/2)


            else:


                frt = -self.frc * (self.bf * self.tf) / (self.bf * self.tf + self.tw * self.dw)

                ## web patch
                if mat_type == 'ElasticPP':
                    ops.uniaxialMaterial(
                        'ElasticPP', start_material_id, self.E, self.fy/self.E, -self.fy/self.E, frt/self.E)
                elif mat_type == 'Steel01':
                    ops.uniaxialMaterial('Steel01', start_material_id+1, self.fy, self.E, self.b)
                    ops.uniaxialMaterial('InitStressMaterial',
                                        start_material_id, start_material_id+1, frt)
                elif mat_type == 'Hardening':
                    ops.uniaxialMaterial('Hardening', start_material_id+1,
                                        self.E, self.fy, 0.0, self.Hk)
                    ops.uniaxialMaterial('InitStressMaterial',
                                        start_material_id, start_material_id+1, frt)
                else:
                    raise Exception(
                        'Input Error - unknown material type (%s)' % mat_type)
                
                ops.section('Fiber', section_id, '-GJ', GJ)                
                # Web patch
                ops.patch('rect', start_material_id, Nfw, 1, -self.tw/2, -self.dw/2, self.tw/2, self.dw/2)

                region_width = self.bf / self.num_regions
                Nff = ceil(region_width * (self.num_fiber / self.bf))

                for i in range(self.num_regions):
                    fri = self.frc + ((i + 0.5) / self.num_regions) * (frt - self.frc)
                    start_material_idi = start_material_id + 2 * (i + 1)

                    if mat_type == 'ElasticPP':
                        ops.uniaxialMaterial(
                            'ElasticPP', start_material_idi, self.E, self.fy/self.E, -self.fy/self.E, fri/self.E)
                    elif mat_type == 'Steel01':
                        ops.uniaxialMaterial(
                            'Steel01', start_material_idi+1, self.fy, self.E, self.b)
                        ops.uniaxialMaterial(
                            'InitStressMaterial', start_material_idi, start_material_idi+1, fri)
                    elif mat_type == 'Hardening':
                        ops.uniaxialMaterial(
                            'Hardening', start_material_idi+1, self.E, self.fy, 0.0, self.Hk)
                        ops.uniaxialMaterial(
                            'InitStressMaterial', start_material_idi, start_material_idi+1, fri)
                    else:
                        raise Exception(
                            'Input Error - unknown material type (%s)' % mat_type)

                    # Top flange segment
                    ops.patch('rect', start_material_idi, Nff, 1, -region_width / 2, self.dw/2, region_width / 2, self.d/2)
                    # Bottom flange segment
                    ops.patch('rect', start_material_idi, Nff, 1, -region_width / 2, -self.d/2, region_width / 2, -self.dw/2)

        elif axis is None:

            # raise ValueError(' The code has not been tested yet for the case when axis is not set to "x" or "y"')

            Nfw_x = ceil(self.tw * (nfx / self.bf))
            Nff_x = ceil(self.bf * (nfx / self.bf))
            Nfw_y = ceil(self.dw * (nfy / self.d))
            Nff_y = ceil(self.tf * (nfy / self.d))

            # if not Residual_Stress:

            if self.frc == 0 or mat_type == 'Elastic':
                if mat_type == 'Elastic':
                    ops.uniaxialMaterial('Elastic', start_material_id, self.E)
                elif mat_type == 'ElasticPP':
                    ops.uniaxialMaterial('ElasticPP', start_material_id,
                                        self.E, self.fy/self.E)
                elif mat_type == 'Steel01':
                    ops.uniaxialMaterial('Steel01', start_material_id, self.fy, self.E, self.b)
                elif mat_type == 'Hardening':
                    ops.uniaxialMaterial('Hardening', start_material_id,
                                        self.E, self.fy, 0.0, self.Hk)
                else:
                    raise Exception(
                        'Input Error - unknown material type (%s)' % mat_type)
                
                ops.section('Fiber', section_id, '-GJ', GJ)

                # Web patch
                ops.patch('rect', start_material_id, Nfw_y, Nfw_x, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)
                # Flange patches
                ops.patch('rect', start_material_id, Nff_y, Nff_x, self.dw / 2, -self.bf / 2, self.d / 2, self.bf / 2)
                ops.patch('rect', start_material_id, Nff_y, Nff_x, -self.d / 2, -self.bf / 2, -self.dw / 2, self.bf / 2)
                # fib_sec = [

                #         ['section', 'Fiber', section_id, '-GJ', 1.0e4],

                #         ["patch","rect",start_material_id, Nfw_y, Nfw_x, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2],
                #         ["patch","rect",start_material_id, Nff_y, Nff_x, self.dw / 2, -self.bf / 2, self.d / 2, self.bf / 2],
                #         ["patch","rect", start_material_id, Nff_y, Nff_x, -self.d / 2, -self.bf / 2, -self.dw / 2, self.bf / 2]

                #         ]
                
                # matcolor = [ 'gold']*1000000
                # opsv.plot_fiber_section(fib_sec,matcolor=matcolor)
                # plt.axis('equal')
                # plt.title('section_name')
                # plt.show()

            else:

                ops.section('Fiber', section_id, '-GJ', GJ)
                

                frt = -self.frc * (self.bf * self.tf) / (self.bf * self.tf + self.tw * self.dw)

                ## web patch
                if mat_type == 'ElasticPP':
                    ops.uniaxialMaterial(
                        'ElasticPP', start_material_id, self.E, self.fy/self.E, -self.fy/self.E, frt/self.E)
                elif mat_type == 'Steel01':
                    ops.uniaxialMaterial('Steel01', start_material_id+1, self.fy, self.E, self.b)
                    ops.uniaxialMaterial('InitStressMaterial',
                                        start_material_id, start_material_id+1, frt)
                elif mat_type == 'Hardening':
                    ops.uniaxialMaterial('Hardening', start_material_id+1,
                                        self.E, self.fy, 0.0, self.Hk)
                    ops.uniaxialMaterial('InitStressMaterial',
                                        start_material_id, start_material_id+1, frt)
                else:
                    raise Exception(
                        'Input Error - unknown material type (%s)' % mat_type)


                ops.patch('rect', start_material_id, Nfw_y, Nfw_x, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)



                region_width = self.bf / self.num_regions
                Nfw_x = ceil(self.tw * (nfx / region_width))
                Nff_x = ceil(region_width * (nfx / region_width))
                for i in range(self.num_regions):
                    fri = self.frc + ((i + 0.5) / self.num_regions) * (frt - self.frc)
                    start_material_idi = start_material_id + 2 * (i + 1)

                    if mat_type == 'ElasticPP':
                        ops.uniaxialMaterial(
                            'ElasticPP', start_material_idi, self.E, self.fy/self.E, -self.fy/self.E, fri/self.E)
                    elif mat_type == 'Steel01':
                        ops.uniaxialMaterial(
                            'Steel01', start_material_idi+1, self.fy, self.E, self.b)
                        ops.uniaxialMaterial(
                            'InitStressMaterial', start_material_idi, start_material_idi+1, fri)
                    elif mat_type == 'Hardening':
                        ops.uniaxialMaterial(
                            'Hardening', start_material_idi+1, self.E, self.fy, 0.0, self.Hk)
                        ops.uniaxialMaterial(
                            'InitStressMaterial', start_material_idi, start_material_idi+1, fri)
                    else:
                        raise Exception(
                            'Input Error - unknown material type (%s)' % mat_type)

                    # Top flange
                    ops.patch('rect', start_material_idi, Nff_y, Nff_x, self.dw / 2, -region_width / 2, self.d / 2, region_width / 2)
                    # Bottom flange
                    ops.patch('rect', start_material_idi, Nff_y, Nff_x, -self.d / 2, -region_width / 2, -self.dw / 2, region_width / 2)

        
        else:
            raise ValueError("Please give valid axis, 'x' or 'y' or set it to None to use a 3d fiber section")


    def maximum_compression_strain(self, axial_strain, curvatureX=0, curvatureY=0):
       
        extreme_strain = axial_strain - self.d/2 * abs(curvatureX) \
                            - self.bf/2 * abs(curvatureY)
        return extreme_strain


    def maximum_tensile_strain(self, axial_strain, curvatureX=0.0, curvatureY=0.0):
        bf = self.bf
        d  = self.d
        xh = bf / 2.0
        yh = d  / 2.0

        # strains at four extreme flange corners
        s1 = axial_strain - (+yh)*curvatureX - (+xh)*curvatureY  # top-right
        s2 = axial_strain - (+yh)*curvatureX - (-xh)*curvatureY  # top-left
        s3 = axial_strain - (-yh)*curvatureX - (+xh)*curvatureY  # bottom-right
        s4 = axial_strain - (-yh)*curvatureX - (-xh)*curvatureY  # bottom-left

        return max(s1, s2, s3, s4)

    def plot_fiber_section(section_id):
        get_fiber_data(section_tag=section_id,plot_fibers=True,keep_json=True)
     