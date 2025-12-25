from math import pi, sqrt, ceil
from libdenavit import opensees as ops
from libdenavit.OpenSees import circ_patch_2d
from libdenavit.section import AciStrainCompatibility, FiberSingle, FiberCirclePatch, FiberSection
from libdenavit.section.database import reinforcing_bar_database
import numpy as np


class CCFT:
    """
     Class representing a Concrete-Filled Circular Tube (CFCT) cross-section.

     Parameters
     ----------
     D : float
         Outer diameter of the steel tube.
     t : float
         Wall thickness of the steel tube.
     Fy : float
         Yield stress of the steel tube.
     fc : float
         Compressive strength of the concrete.
     units : str
         Units system, either 'US' or 'SI'.
     num_bars : int, optional
         Number of longitudinal reinforcing bars. Default is 0 (no bars).
     bar_size : str, optional
         Standard size designation of the reinforcing bars (e.g., '#8').
         If provided, `Ab` and `db` will be taken from the database.
     Ab : float, optional
         Cross-sectional area of one reinforcing bar (if bar_size not provided).
     db : float, optional
         Diameter of one reinforcing bar (if bar_size not provided).
     Fylr : float, optional
         Yield stress of the longitudinal reinforcing bars.
     Dp : float, optional
         Radial distance from the outside surface of the steel tube to the center of the longitudinal reinforcement.
     neglect_local_buckling : bool, optional
         Whether to neglect local buckling effects. Default is False.
     """
    _Es = None
    _Ec = None
    
    def __init__(self, D, t, Fy, fc, units, num_bars=0, bar_size=None, Ab=None, db=None, Fylr=None, Dp=None,
                 reinforcement=None, neglect_local_buckling=False):
        # Main Parameters
        self.D = D
        self.t = t
        self.Fy = Fy
        self.fc = fc
        self.units = units
        self._eps_c = None
        
        # Reinforcment
        self.reinforcement = reinforcement
        self.num_bars = num_bars
        if bar_size is None:
            self.Ab = Ab
            self.db = db
        else:
            bar_data = reinforcing_bar_database[bar_size]
            self.Ab = bar_data['area']
            self.db = bar_data['diameter']
        self.Fylr = Fylr
        self.Dp = Dp
        
        # Options
        self.neglect_local_buckling = neglect_local_buckling


    @property
    def Es(self):
        if self._Es is not None:
            return self._Es

        if self.units.lower() == 'us':
            return 29000

        if self.units.lower() == 'si':
            return 200000

        raise ValueError(f'Es is not set and default is not implemented for {self.units = }')

    @Es.setter
    def Es(self, x):
        self._Es = x
         
    @property
    def Ec(self):
        if self._Ec is not None:
            return self._Ec

        if self.units.lower() == 'us':
            return 57 * sqrt(self.fc * 1000)

        if self.units.lower() == 'si':
            return 4700 * sqrt(self.fc)

        raise ValueError(f'Ec is not set and default is not implemented for {self.units = }')

    @Ec.setter
    def Ec(self, x):
        self._Ec = x
         
    @property
    def As(self):
        return 0.25*pi*(self.D**2 - (self.D-2*self.t)**2)

    @property
    def reinforcement(self):
        return self._reinforcement

    @reinforcement.setter
    def reinforcement(self, x):
        if x is None:
            self._reinforcement = None
        elif type(x) == list:
            self._reinforcement = x
        else:
            self._reinforcement = [x]

    @property
    def Asr(self):
        if self.num_bars == 0:
            return 0.0
        else:
            return self.num_bars*self.Ab
    
    @property
    def Ac(self):
        return 0.25*pi*(self.D-2*self.t)**2 - self.Asr
    
    def Is(self,axis='x'):
        return pi/64*(self.D**4 - (self.D-2*self.t)**4)

    def Isr(self, axis='x'):
        I_bars = 0.0
        if self.reinforcement is None:
            if self.num_bars > 0:
                x, y = self.reinforcing_coordinates()
                if axis.lower() == 'x':
                    # Iy = sum(Ab * y^2)
                    for i in range(self.num_bars):
                        I_bars += self.Ab * (y[i] ** 2)
                elif axis.lower() == 'y':
                    # Ix = sum(Ab * x^2)
                    for i in range(self.num_bars):
                        I_bars += self.Ab * (x[i] ** 2)
                else:
                    raise ValueError("axis must be 'x' or 'y'")
        else:
            i = 0
            for j in self.reinforcement:
                i += j.I(axis)
            return i

        return I_bars
            
    def Ic(self,axis='x'):
        return pi/64*(self.D-2*self.t)**4 - self.Isr(axis)

    @property
    def eps_c(self):
        if self._eps_c is not None:
            return self._eps_c
        if self.units.lower() == 'us':
            return (self.fc * 1000) ** (1 / 4) / 4000
        if self.units.lower() == 'si':
            return (self.fc * 145.038) ** (1 / 4) / 4000

        raise ValueError(f'eps_c is not set and default is not impleted for {self.units = }')

    @property
    def Ag(self):
        return 0.25*pi*self.D**2

    def Ig(self, axis='x'):
        return pi/64*(self.D**4)

    @eps_c.setter
    def eps_c(self, x):
        self._eps_c = x

    def reinforcing_coordinates(self):
        angles = np.linspace(0,2*pi,self.num_bars,endpoint=False)
        r = 0.5*self.D - self.Dp
        x = r*np.sin(angles)
        y = r*np.cos(angles)
        return x,y

    def maximum_concrete_compression_strain(self, axial_strain, curvatureX=0, curvatureY=0):
        extreme_strain = axial_strain - (self.D / 2 - self.t) * sqrt(curvatureX ** 2 + curvatureY ** 2)
        return extreme_strain

    def maximum_tensile_steel_strain(self, axial_strain, curvatureX=0, curvatureY=0):
        max_strain = float('-inf')
        for i in self.reinforcement[0].coordinates:
            strain = axial_strain - (i[1] - self.reinforcement[0].yc) * curvatureX \
                     - (i[0] - self.reinforcement[0].xc) * curvatureY
            if strain > max_strain:
                max_strain = strain
        return max_strain


    def aci_strain_compatibility_object(self):
        id_steel = 1
        id_conc = 2
        id_reinf = 3       
        fs = self.fiber_section_object(id_steel, id_conc, id_reinf)
        scACI = AciStrainCompatibility(fs)
        scACI.add_concrete_boundary(0.0, 0.0, 0.5*self.D-self.t)
        scACI.add_steel_boundary(0.0, 0.0, 0.5*self.D)
        scACI.add_material(id_steel, 'steel', self.Fy, self.Es)
        scACI.add_material(id_conc, 'concrete', self.fc, self.units)
        scACI.add_material(id_reinf, 'steel', self.Fylr, self.Es)
        return scACI    
    
    def fiber_section_object(self, id_steel, id_conc, id_reinf, nfx=200, nfy=200):
        fs = FiberSection(nfx, nfy)
        # Steel
        fs.add_fibers(FiberCirclePatch(0,0,0.5*self.D-self.t,0.5*self.D,id_steel))
        # Concrete 
        fs.add_fibers(FiberCirclePatch(0,0,0,0.5*self.D-self.t,id_conc))
        # Reinforcement
        if self.num_bars > 0:
            x,y = self.reinforcing_coordinates()
            for i in range(self.num_bars):
                fs.add_fibers(FiberSingle(self.Ab,x[i],y[i],id_reinf,id_conc))
        return fs

    def build_ops_fiber_section(self, section_id, start_material_id, steel_mat_type, conc_mat_type, tube_mat_type,
                                nfx, nfy, GJ=1.0e6, axis=None):

        # region Define Material IDs
        steel_material_id = start_material_id
        concrete_material_id = start_material_id + 1
        tube_material_id = start_material_id + 2
        # endregion

        # region Define ops Section
        ops.section('Fiber', section_id, '-GJ', GJ)
        # endregion

        # region Define Steel Material
        if self.num_bars != 0:
            if steel_mat_type == "ElasticPP":
                ops.uniaxialMaterial("ElasticPP", steel_material_id, self.Es, self.Fylr / self.Es)

            elif steel_mat_type == "Hardening":
                ops.uniaxialMaterial("Hardening", steel_material_id, self.Es, self.Fylr, 0.001*self.Es, 0)

            elif steel_mat_type == "ReinforcingSteel":
                ops.uniaxialMaterial("ReinforcingSteel", steel_material_id, self.Fylr, self.Fylr * 1.5, self.Es, self.Es / 2,
                                     0.002, 0.008)

            elif steel_mat_type == "Elastic":
                ops.uniaxialMaterial("Elastic", steel_material_id, self.Es)

            else:
                raise ValueError(f"Steel material {steel_mat_type} not supported")
        # endregion

        # region Define Concrete Material
        if conc_mat_type == "Concrete04_no_confinement":
            ops.uniaxialMaterial("Concrete04", concrete_material_id, -self.fc, -self.eps_c, -1.0, self.Ec)

        elif conc_mat_type == "Concrete01_no_confinement":
            ops.uniaxialMaterial("Concrete01", concrete_material_id, -self.fc, -2 * self.fc / self.Ec, 0, 0.006)

        elif conc_mat_type == "ENT":
            ops.uniaxialMaterial('ENT', concrete_material_id, self.Ec)

        elif conc_mat_type == "Elastic":
            ops.uniaxialMaterial('Elastic', concrete_material_id, self.Ec)

        else:
            raise ValueError(f"Concrete material {conc_mat_type} not supported")
        # endregion

        # region Define Tube Material
        if tube_mat_type == "ElasticPP":
            ops.uniaxialMaterial("ElasticPP", tube_material_id, self.Es, self.Fy / self.Es)
        # endregion

        # region Define Fibers and Patches
        if (axis is None) or (axis == '3d'):
            raise NotImplementedError('3D fiber section not implemented')

        elif axis in ['x', 'y']:
            if self.num_bars != 0:
                for i in self.reinforcement:
                    if axis == 'x':
                        rebar_coords = i.coordinates[1]
                    else:
                        rebar_coords = i.coordinates[0]

                    for index, value in enumerate(i.coordinates[1]):
                        ops.fiber(value, 0, i.Ab, steel_material_id)
                        negative_area_material_id = concrete_material_id
                        ops.fiber(value, 0, -i.Ab, negative_area_material_id)

            d = self.D - 2 * self.t
            if axis == 'x':
                max_fiber_size = d / nfy
            else:
                max_fiber_size = d / nfx

            nf = ceil(0.5*d/max_fiber_size)
            circ_patch_2d(concrete_material_id, nf, d)
            circ_patch_2d(tube_material_id, ceil(nf * (self.D-d)/d), self.D, Di=d)
        else:
            raise NotImplementedError(f'Axis {axis} not implemented')

        # endregion


def run_example():
    #section = CCFT(24,0.25,50,6,'US')
    section = CCFT(24,0.25,50,6,'US',num_bars=8, bar_size='#8', Fylr=60, Dp=3.0)
    
    fs = section.fiber_section_object(1,2,3)
    fs.print_section_properties()
    
    print(f'As  = {section.As}')
    print(f'Ac  = {section.Ac}')
    print(f'Asr = {section.Asr}')
    
    
if __name__ == '__main__':
    run_example()   
    
    