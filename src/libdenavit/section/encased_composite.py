from math import pi, sqrt
from libdenavit.section import AciStrainCompatibility, FiberSingle, FiberCirclePatch, FiberSection, Rectangle, FiberQuadPatch, RectangularTube
from libdenavit.section.database import reinforcing_bar_database


class SRC:
    _Es = None
    _Ec = None

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
            self.db = db
        else:
            bar_data = reinforcing_bar_database[bar_size]
            self.Ab = bar_data['area']
            self.db = bar_data['diameter']
        self.Fylr = Fylr  # Yield strength of longitudinal reinforcing
        self.Dp = Dp  # Distance from outside surface of concrete to center of longitudinal reinforcing

        # Options
        self.C2 = 0.85

    @property
    def num_bars(self):
        if self.nbH != 0 and self.nbB != 0:
            return 2 * (self.nbH + self.nbB) - 4
        else:
            return 0

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

    def Is(self, axis):
        if axis.lower() == 'x':
            return (1/12)*self.bf*self.d**3 - (1/12)*(self.bf-self.tw)*self.h**3 
        elif axis.lower() == 'y':
            return 2*(1/12)*self.tf*self.bf**3 - (1/12)*self.h*self.tw**3
        else:
            raise Exception(f'Unknown axis: {axis}')

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
            raise ValueError(f'Unknown option for axis: {axis}') 
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

