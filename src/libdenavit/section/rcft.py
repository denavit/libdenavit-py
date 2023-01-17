from math import pi, sqrt
from libdenavit.section import AciStrainCompatibility, FiberSingle, FiberCirclePatch, FiberSection, Rectangle, FiberQuadPatch, RectangularTube
from libdenavit.section.database import reinforcing_bar_database


class RCFT:
    _Es = None
    _Ec = None

    def __init__(self, H, B, t, Fy, fc, units, ri=None, nbH=0, nbB=0, bar_size=None, Ab=None, db=None, Fylr=None, Dp=None,
                 neglect_local_buckling=False):
        # Main Parameters
        self.H = H
        self.B = B
        self.t = t
        self.Fy = Fy
        self.fc = fc
        self.units = units

        # Internal radius of tube
        if ri is None:
            self.ri = 0.
        else:
            self.ri = ri

        # Reinforcing
        if type(nbH) is not int or type(nbB) is not int:
            raise ValueError('nbH and nbB must be integers')
        if 1 in (nbH, nbB) or (nbH < 0 or nbB < 0):
            raise ValueError('nbH or nbB value is not valid')
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
        self.Dp = Dp  # Distance from outside surface of steel tube to center of longitudinal reinforcing

        # Options
        self.neglect_local_buckling = neglect_local_buckling
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
    def ro(self):
        if self.ri == 0:
            return 0
        else:
            return self.ri + self.t

    @property
    def Hc(self):
        return self.H - 2 * self.t

    @property
    def Bc(self):
        return self.B - 2 * self.t

    @property
    def As(self):
        return self.Ag - self.Ac - self.Asr

    @property
    def Ac(self):
        return self.Bc * self.Hc - (4 - pi) * self.ri ** 2 - self.Asr

    @property
    def Ag(self):
        return self.B * self.H - (4 - pi) * self.ri ** 2

    def Is(self, axis):
        shp = RectangularTube(self.H, self.B, self.t, self.ro)
        return shp.I(axis)

    def Ic(self, axis):
        shp = Rectangle(self.Hc, self.Bc, self.ri)
        return shp.I(axis) - self.Isr(axis)

    def Isr(self, axis):
        if self.nbH == 0 or self.nbB == 0:
            return 0

        else:
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
        shp = Rectangle(self.H, self.B, self.ri)
        return shp.I(axis)

    def Js(self, axis):
        shp = RectangularTube(self.H, self.B, self.t, self.ro)
        return shp.J(axis)

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
        else:
            return 0, 0

    def aci_strain_compatibility_object(self):
        id_steel = 1
        id_conc = 2
        id_reinf = 3
        fs = self.fiber_section_object(id_steel, id_conc, id_reinf)
        scACI = AciStrainCompatibility(fs, "RCFT")
        scACI.add_concrete_boundary(-self.Bc/2 + self.ri, -self.Hc/2 + self.ri, self.ri)
        scACI.add_concrete_boundary(-self.Bc/2 + self.ri, +self.Hc/2 - self.ri, self.ri)
        scACI.add_concrete_boundary(+self.Bc/2 - self.ri, -self.Hc/2 + self.ri, self.ri)
        scACI.add_concrete_boundary(+self.Bc/2 - self.ri, +self.Hc/2 - self.ri, self.ri)

        scACI.add_steel_boundary(-self.B/2 + self.ro, -self.H/2 + self.ro, self.ro)
        scACI.add_steel_boundary(-self.B/2 + self.ro, +self.H/2 - self.ro, self.ro)
        scACI.add_steel_boundary(+self.B/2 - self.ro, -self.H/2 + self.ro, self.ro)
        scACI.add_steel_boundary(+self.B/2 - self.ro, +self.H/2 - self.ro, self.ro)

        scACI.add_material(id_steel, 'steel', self.Fy, self.Es)
        scACI.add_material(id_conc, 'concrete', self.fc, self.units)
        scACI.add_material(id_reinf, 'steel', self.Fylr, self.Es)

        return scACI

    def fiber_section_object(self, id_steel, id_conc, id_reinf, nfx=200, nfy=200):
        fs = FiberSection(nfx, nfy)
        # Steel and Concrete
        if self.ri == 0:
            b1 = self.B/2 - self.t
            b2 = self.B/2
            d1 = self.H/2 - self.t
            d2 = self.H/2

            # Steel Section
            fs.add_fibers(FiberQuadPatch(-b1, d1,-b1, d2, b1, d2, b1, d1, id_steel))
            fs.add_fibers(FiberQuadPatch(-b1,-d2,-b1,-d1, b1,-d1, b1,-d2, id_steel))
            fs.add_fibers(FiberQuadPatch(-b2,-d2,-b2, d2,-b1, d2,-b1,-d2, id_steel))
            fs.add_fibers(FiberQuadPatch( b1,-d2, b1, d2, b2, d2, b2,-d2, id_steel))

            # Concrete Section
            fs.add_fibers(FiberQuadPatch(-b1,-d1,-b1, d1, b1, d1, b1,-d1, id_conc))

        else:
            b1 = self.B / 2 - self.ro
            b2 = self.B / 2 - (self.ro - self.ri)
            b3 = self.B / 2
            d1 = self.H / 2 - self.ro
            d2 = self.H / 2 - (self.ro - self.ri)
            d3 = self.H / 2

            # Steel Section
            fs.add_fibers(FiberQuadPatch(-b1, d2, -b1, d3, b1, d3, b1, d2, id_steel))
            fs.add_fibers(FiberQuadPatch(-b1, -d3, -b1, -d2, b1, -d2, b1, -d3, id_steel))
            fs.add_fibers(FiberQuadPatch(b2, -d1, b2, d1, b3, d1, b3, -d1, id_steel))
            fs.add_fibers(FiberQuadPatch(-b3, -d1, -b3, d1, -b2, d1, -b2, -d1, id_steel))
            fs.add_fibers(FiberCirclePatch(b1, d1, self.ri, self.ro, id_steel, a1=0, a2=pi / 2))
            fs.add_fibers(FiberCirclePatch(-b1, d1, self.ri, self.ro, id_steel, a1=pi / 2, a2=pi))
            fs.add_fibers(FiberCirclePatch(-b1, -d1, self.ri, self.ro, id_steel, a1=pi, a2=3*pi / 2))
            fs.add_fibers(FiberCirclePatch(b1, -d1, self.ri, self.ro, id_steel, a1=3 * pi / 2, a2=2 * pi))

            # Concrete Section
            fs.add_fibers(FiberQuadPatch(-b2, -d1, -b2, d1, b2, d1, b2, -d1, id_conc))
            fs.add_fibers(FiberQuadPatch(-b1, d1, -b1, d2, b1, d2, b1, d1, id_conc))
            fs.add_fibers(FiberQuadPatch(-b1, -d2, -b1, -d1, b1, -d1, b1, -d2, id_conc))
            fs.add_fibers(FiberCirclePatch(b1, d1, 0.0, self.ri, id_conc, a1=0, a2=pi/2))
            fs.add_fibers(FiberCirclePatch(-b1, d1, 0.0, self.ri, id_conc, a1=pi/2, a2=pi))
            fs.add_fibers(FiberCirclePatch(-b1, -d1, 0.0, self.ri, id_conc, a1=pi, a2=3*pi / 2))
            fs.add_fibers(FiberCirclePatch(b1, -d1, 0.0, self.ri, id_conc, a1=3*pi / 2, a2=2*pi))

        # Reinforcement
        if self.nbH > 0:
            x, y = self.reinforcing_coordinates()
            for i in range(self.num_bars):
                fs.add_fibers(FiberSingle(self.Ab, x[i], y[i], id_reinf, id_conc))
        return fs


def run_example():
    H = 40
    B = 20
    t = 1
    fy = 60
    fc = 4
    ri = 0

    section = RCFT(H, B, t, fy, fc, 'US', ri, 1, 4, Ab=0.2, Fylr=60, Dp=3)

    fs = section.fiber_section_object(1, 2, 3)
    fs.plot_fibers()
    fs.print_section_properties()

    print(f'As  = {section.As}')
    print(f'Ac  = {section.Ac}')
    print(f'Asr = {section.Asr}')
    
    axis = 'y'
    Is  = section.Is(axis)
    Ic  = section.Ic(axis)
    Isr = section.Isr(axis)
    print(f'Is  = {Is}')
    print(f'Ic  = {Ic}')
    print(f'Isr = {Isr}')

    scACI = section.aci_strain_compatibility_object()


if __name__ == '__main__':
    run_example()

