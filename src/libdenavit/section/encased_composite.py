from math import ceil, pi, sqrt
import matplotlib.pyplot as plt
import numpy as np
from libdenavit import opensees as ops
from libdenavit.section import AciStrainCompatibility, FiberSingle, FiberCirclePatch, FiberSection, Rectangle, FiberQuadPatch, RectangularTube
from libdenavit.section.database import reinforcing_bar_database


class SRC:

    def __init__(self, H, B, d, bf, tf, tw, Fy, fc, units, nbH, nbB, Fylr, Dp, bar_size=None, Ab=None, db=None):
        # Main Parameters
        self.H = H
        self.B = B
        self.d = d
        self.bf = bf
        self.tf = tf
        self.tw = tw
        self.Fy = Fy
        self.fc = fc
        self.units = units

        # Reinforcing
        self.nbH = nbH
        self.nbB = nbB
        if bar_size is None:
            self.Ab = Ab
            self._db = db
        else:
            bar_data = reinforcing_bar_database[bar_size]
            self.Ab = bar_data['area']
            self._db = bar_data['diameter']
        self.Fylr = Fylr  # Yield strength of longitudinal reinforcing
        self.Dp = Dp  # Distance from outside surface of concrete to center of longitudinal reinforcing

        # Options
        self.C2 = 0.85

        # Internal variables
        self._Es = None
        self._Ec = None


    @property
    def num_bars(self):
        if self.nbH != 0 and self.nbB != 0:
            return 2 * (self.nbH + self.nbB) - 4
        else:
            return 0

    @property
    def db(self):
        if self._db is not None:
            return self._db
        else:
            return sqrt(4 * self.Ab / pi)

    @db.setter
    def db(self, x):
        self._db = x

    @property
    def Asr(self):
        if self.num_bars == 0:
            return 0.0
        else:
            return self.num_bars * self.Ab

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
    def h(self):
        return self.d - 2*self.tf

    @property
    def As(self):
        return self.h*self.tw + 2*self.bf*self.tf

    @property
    def Ac(self):
        return self.H * self.B - self.As - self.Asr

    @property
    def Ag(self):
        return self.B * self.H

    @property
    def p0(self):
        # Nominal axial strength of section
        p0 = self.Fy * self.As + self.Fylr * self.Asr + 0.85 * self.fc * self.Ac
        return p0

    def depth(self, axis):
        if axis.lower() == 'x':
            return self.H 
        elif axis.lower() == 'y':
            return self.B
        else:
            raise ValueError(f'Invalid axis: {axis}')

    def Is(self, axis):
        if axis.lower() == 'x':
            return (1/12)*self.bf*self.d**3 - (1/12)*(self.bf-self.tw)*self.h**3 
        elif axis.lower() == 'y':
            return 2*(1/12)*self.tf*self.bf**3 - (1/12)*self.h*self.tw**3
        else:
            raise ValueError(f'Invalid axis: {axis}')

    def Ic(self, axis):
        return self.Ig(axis) - self.Is(axis) - self.Isr(axis)

    def Isr(self, axis):
        I = 0
        (x,y) = self.reinforcing_coordinates()
        if axis == 'x':
            for yi in y:
                I += yi*yi
        elif axis == 'y':
            for xi in x:
                I += xi*xi
        else:
            raise ValueError(f'Invalid axis: {axis}') 
        I = I*self.Ab
        return I

    def Ig(self, axis):
        shp = Rectangle(self.H, self.B)
        return shp.I(axis)

    def reinforcing_coordinates(self):
        if self.nbH != 0 and self.nbB != 0:
            x = []
            y = []
            for i in range(self.nbH):
                for j in range(self.nbB):
                    if i in [0, self.nbH - 1] or j in [0, self.nbB - 1]:
                        x.append(j * (self.B - 2 * self.Dp) / (self.nbB - 1) - self.B / 2 + self.Dp)
                        y.append(i * (self.H - 2 * self.Dp) / (self.nbH - 1) - self.H / 2 + self.Dp)
        return x, y

    def aci_strain_compatibility_object(self):
        id_steel = 1
        id_conc = 2
        id_reinf = 3
        fs = self.fiber_section_object(id_steel, id_conc, id_reinf)
        scACI = AciStrainCompatibility(fs, "RCFT")
        scACI.add_concrete_boundary(-self.B/2, -self.H/2, 0)
        scACI.add_concrete_boundary(-self.B/2, +self.H/2, 0)
        scACI.add_concrete_boundary(+self.B/2, -self.H/2, 0)
        scACI.add_concrete_boundary(+self.B/2, +self.H/2, 0)

        assert (self.B - 2*self.Dp) > self.bf, 'reinforcing bars must be placed wider than the steel shape'
        assert (self.H - 2*self.Dp) > self.h,  'reinforcing bars must be placed deeper than the steel shape'
        scACI.add_steel_boundary(-self.B/2 + self.Dp, -self.H/2 + self.Dp, 0)
        scACI.add_steel_boundary(-self.B/2 + self.Dp, +self.H/2 - self.Dp, 0)
        scACI.add_steel_boundary(+self.B/2 - self.Dp, -self.H/2 + self.Dp, 0)
        scACI.add_steel_boundary(+self.B/2 - self.Dp, +self.H/2 - self.Dp, 0)

        scACI.add_material(id_steel, 'steel', self.Fy, self.Es)
        scACI.add_material(id_conc, 'concrete', self.fc, self.units)
        scACI.add_material(id_reinf, 'steel', self.Fylr, self.Es)

        return scACI

    def fiber_section_object(self, id_steel=1, id_conc=2, id_reinf=3, nfx=200, nfy=200):
        fs = FiberSection(nfx, nfy)

        tw2 = self.tw/2
        bf2 = self.bf/2
        B2 = self.B/2
        h2 = self.h/2
        d2 = self.d/2
        H2 = self.H/2
        
        # Steel 
        fs.add_fibers(FiberQuadPatch(-bf2,  h2, -bf2,  d2, bf2,  d2, bf2,  h2, id_steel))
        fs.add_fibers(FiberQuadPatch(-tw2, -h2, -tw2,  h2, tw2,  h2, tw2, -h2, id_steel))      
        fs.add_fibers(FiberQuadPatch(-bf2, -d2, -bf2, -h2, bf2, -h2, bf2, -d2, id_steel))

        # Concrete
        fs.add_fibers(FiberQuadPatch( -B2,  d2,  -B2,  H2,  B2,  H2,  B2,  d2, id_conc))
        fs.add_fibers(FiberQuadPatch( -B2,  h2,  -B2,  d2,-bf2,  d2,-bf2,  h2, id_conc))
        fs.add_fibers(FiberQuadPatch( bf2,  h2,  bf2,  d2,  B2,  d2,  B2,  h2, id_conc))
        fs.add_fibers(FiberQuadPatch( -B2, -h2,  -B2,  h2,-tw2,  h2,-tw2, -h2, id_conc))
        fs.add_fibers(FiberQuadPatch( tw2, -h2,  tw2,  h2,  B2,  h2,  B2, -h2, id_conc))
        fs.add_fibers(FiberQuadPatch( -B2, -d2,  -B2, -h2,-bf2, -h2,-bf2, -d2, id_conc))
        fs.add_fibers(FiberQuadPatch( bf2, -d2,  bf2, -h2,  B2, -h2,  B2, -d2, id_conc))
        fs.add_fibers(FiberQuadPatch( -B2, -H2,  -B2, -d2,  B2, -d2,  B2, -H2, id_conc))
        
        # Reinforcement
        if self.nbH > 0:
            x, y = self.reinforcing_coordinates()
            for i in range(self.num_bars):
                fs.add_fibers(FiberSingle(self.Ab, x[i], y[i], id_reinf, id_conc))
        return fs

    def build_ops_fiber_section(self, section_id, steel_mat_type, reinf_mat_type, conc_mat_type, nfy, nfx, 
                                GJ=1.0e6, axis=None, start_material_id=1):
        """ Builds the fiber section object

        Parameters
        ----------
        section_id : int
            The id of the section
        steel_mat_type : str
            The type of the structural steel material
        reinf_mat_type : str
            The type of the reinforcing steel material
        conc_mat_type : str
            The type of the concrete material
        nfy : int
            The minimum number of fibers in the y direction
        nfx : int
            The minimum number of fibers in the x direction
        GJ : float (default 1.0e6)
            The torsional rigidity of the cross section
        axis : str (default None)
            Optional argument for defining a fiber section for 2-dimensional analysis
              - If "None", then a 3-dimensional fiber section will be defined
              - If "x", then a 2-dimensional fiber section will be defined based on 
                bending about the x-axis (note that the value nfx will be ignored) 
              - If "y", then a 2-dimensional fiber section will be defined based on 
                bending about the y-axis (note that the value nfy will be ignored)
        start_material_id : int (default 1)
            The id of the first uniaxial material to be defined (others will be defined sequentially)
        """

        # Three or four uniaxial materials are defined in this function
        steel_material_id = start_material_id
        reinforcing_material_id = start_material_id+1
        concrete_material_id = start_material_id+2        # Used if no confinement
        cover_concrete_material_id = start_material_id+2  # Used if confinement
        core_concrete_material_id  = start_material_id+3  # Used if confinement
        
        # Define section      
        ops.section('Fiber', section_id, '-GJ', GJ)

        # region Define Structural Steel Material
        if steel_mat_type == "Elastic":
            ops.uniaxialMaterial("Elastic", steel_material_id, self.Es)
        
        elif steel_mat_type == "ElasticPP":
            ops.uniaxialMaterial("ElasticPP", steel_material_id, self.Es, self.Fy / self.Es)

        elif steel_mat_type == "Hardening":
            ops.uniaxialMaterial("Hardening", steel_material_id, self.Es, self.Fy, 0.001*self.Es, 0)

        else:
            raise ValueError(f"Steel material {steel_mat_type} not supported")
        # endregion


        # region Define Reinforcing Steel Material
        if reinf_mat_type == "Elastic":
            ops.uniaxialMaterial("Elastic", reinforcing_material_id, self.Es)
        
        elif reinf_mat_type == "ElasticPP":
            ops.uniaxialMaterial("ElasticPP", reinforcing_material_id, self.Es, self.Fylr / self.Es)

        elif reinf_mat_type == "Hardening":
            ops.uniaxialMaterial("Hardening", reinforcing_material_id, self.Es, self.Fylr, 0.001*self.Es, 0)

        else:
            raise ValueError(f"Steel material {reinf_mat_type} not supported")
        # endregion

            
        # region Define Concrete Material
        if conc_mat_type == "Elastic":
            ops.uniaxialMaterial('Elastic', concrete_material_id, self.Ec)
            confinement = False
        
        elif conc_mat_type == "ENT":
            ops.uniaxialMaterial('ENT', concrete_material_id, self.Ec)
            confinement = False

        elif conc_mat_type == "Concrete01_no_confinement":
            ops.uniaxialMaterial("Concrete01", concrete_material_id, -self.fc, -2*self.fc/self.Ec, 0, -0.006)
            confinement = False       

        elif conc_mat_type == "Concrete04_no_confinement":
            ops.uniaxialMaterial("Concrete04", concrete_material_id, -self.fc, -self.eps_c, -1.0, self.Ec)
            confinement = False

        else:
            raise ValueError(f"Concrete material {conc_mat_type} not supported")
        # endregion


        # region Define fibers and patches
        tw2 = self.tw/2
        bf2 = self.bf/2
        B2 = self.B/2
        h2 = self.h/2
        d2 = self.d/2
        H2 = self.H/2 
        
        if (axis is None) or (axis == '3d') or (axis == 'x'): 

            # Structural Steel 
            nfIJ = ceil(self.tf*nfy/self.H)
            nfJK = ceil(self.bf*nfx/self.B)
            ops.patch('quad', steel_material_id, nfIJ, nfJK,  h2, -bf2,  d2, -bf2,  d2, bf2,  h2, bf2)
            ops.patch('quad', steel_material_id, nfIJ, nfJK, -d2, -bf2, -h2, -bf2, -h2, bf2, -d2, bf2)
            
            nfIJ = ceil(self.h*nfy/self.H)
            nfJK = ceil(self.tw*nfx/self.B)
            ops.patch('quad', steel_material_id, nfIJ, nfJK, -h2, -tw2,  h2, -tw2,  h2, tw2, -h2, tw2)      
        
            # Reinforcing Steel
            (x_reinf,y_reinf) = self.reinforcing_coordinates()
            for x, y in zip(x_reinf,y_reinf):
                ops.fiber(y, x, self.Ab, reinforcing_material_id)
                if confinement:
                    negative_area_material_id = core_concrete_material_id
                else:
                    negative_area_material_id = concrete_material_id
                ops.fiber(y, x, -self.Ab, negative_area_material_id)
                    
            # Concrete
            if confinement:
                raise ValueError('Not yet implemented for concrete confinement')
            else:
                nfIJ = ceil((H2-d2)*nfy/self.H)
                nfJK = ceil(self.B*nfx/self.B)
                ops.patch('quad', concrete_material_id, nfIJ, nfJK,  d2, -B2,  H2,  -B2,  H2,  B2,  d2,  B2)
                ops.patch('quad', concrete_material_id, nfIJ, nfJK, -H2, -B2, -d2,  -B2, -d2,  B2, -H2,  B2)
                                                                                                      
                nfIJ = ceil(self.tf*nfy/self.H)
                nfJK = ceil((B2-bf2)*nfx/self.B)
                ops.patch('quad', concrete_material_id, nfIJ, nfJK,  h2, -B2,  d2,  -B2,  d2,-bf2,  h2,-bf2)
                ops.patch('quad', concrete_material_id, nfIJ, nfJK,  h2, bf2,  d2,  bf2,  d2,  B2,  h2,  B2)
                ops.patch('quad', concrete_material_id, nfIJ, nfJK, -d2, -B2, -h2,  -B2, -h2,-bf2, -d2,-bf2)
                ops.patch('quad', concrete_material_id, nfIJ, nfJK, -d2, bf2, -h2,  bf2, -h2,  B2, -d2,  B2)
                                                                                                      
                nfIJ = ceil(self.h*nfy/self.H)
                nfJK = ceil((B2-tw2)*nfx/self.B)
                ops.patch('quad', concrete_material_id, nfIJ, nfJK, -h2, -B2,  h2,  -B2,  h2,-tw2, -h2,-tw2)
                ops.patch('quad', concrete_material_id, nfIJ, nfJK, -h2, tw2,  h2,  tw2,  h2,  B2, -h2,  B2)

        elif axis in ['x', 'y']:
            raise ValueError('build_ops_fiber_section not yet implemented for two-dimensional analyses')

        else:
            raise ValueError(f'Invalid axis: {axis}')
        # endregion

    def maximum_concrete_compression_strain(self, axial_strain, curvatureX=0, curvatureY=0):
        extreme_strain = axial_strain - self.H/2 * abs(curvatureX) - self.B/2 * abs(curvatureY)
        return extreme_strain

    def maximum_tensile_steel_strain(self, axial_strain, curvatureX=0, curvatureY=0):
        max_strain = float('-inf')
        (x_reinf,y_reinf) = self.reinforcing_coordinates()
        for x, y in zip(x_reinf,y_reinf):
            strain = axial_strain - y * curvatureX - x * curvatureY
            if strain > max_strain:
                max_strain = strain
        return max_strain
        
    def plot_section(self, show=True, **kwargs):
        plt.figure()

        # Concrete
        x = [self.B / 2, - self.B / 2, - self.B / 2, self.B / 2, self.B / 2]
        y = [self.H / 2, self.H / 2, - self.H / 2, - self.H / 2, self.H / 2]
        plt.fill(x, y, edgecolor='k', facecolor=[0.9,0.9,0.9], **kwargs)

        # Steel Shape
        x = [self.bf / 2, -self.bf / 2, - self.bf / 2, -self.tw / 2, -self.tw / 2, -self.bf / 2, - self.bf / 2, 
             self.bf / 2,  self.bf / 2,   self.tw / 2,  self.tw / 2,  self.bf / 2,  self.bf / 2]
        y = [self.d / 2, self.d / 2, self.h / 2, self.h / 2, -self.h / 2, -self.h / 2, -self.d / 2, 
             -self.d / 2, -self.h / 2, -self.h / 2, self.h / 2, self.h / 2, self.d / 2]
        plt.fill(x, y, edgecolor='k', facecolor=[0.5,0.5,0.5], **kwargs)
        
        # Reinforcement
        (x_reinf,y_reinf) = self.reinforcing_coordinates()
        angles = np.linspace(0,2*pi,16)
        x_bar = (self.db/2)*np.cos(angles)
        y_bar = (self.db/2)*np.sin(angles)
        for x, y in zip(x_reinf,y_reinf):
            plt.fill(x_bar+x, y_bar+y, edgecolor='k', facecolor=[0.5,0.5,0.5], **kwargs)
        
        plt.box(False)
        plt.axis('scaled')
        if show:
            plt.show()
        
def run_example():
    H = 28
    B = 24
    d = 14
    bf = 14
    tf = 1.5
    tw = 0.75
    Fy = 50
    fc = 4
    units = 'US'   
    nbH = 2
    nbB = 2
    Fylr = 60
    Dp = 3
    bar_size = '#8'

    section = SRC(H, B, d, bf, tf, tw, Fy, fc, units, nbH, nbB, Fylr, Dp, bar_size)

    fs = section.fiber_section_object()
    fs.plot_fibers()
    fs.print_section_properties()

    print(f'As  = {section.As}')
    print(f'Ac  = {section.Ac}')
    print(f'Asr = {section.Asr}')
    
    axis = 'x'
    Is  = section.Is(axis)
    Ic  = section.Ic(axis)
    Isr = section.Isr(axis)
    print(f'Is  = {Is}')
    print(f'Ic  = {Ic}')
    print(f'Isr = {Isr}')

if __name__ == '__main__':
    run_example()

