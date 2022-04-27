from math import sqrt, e, pi, ceil
from libdenavit.section import AciStrainCompatibility, FiberSingle, FiberSection, ACI_phi
import matplotlib.pyplot as plt
import numpy as np
import openseespy.opensees as ops
import openseespy.postprocessing.ops_vis as opsv
import openseespy.postprocessing.Get_Rendering as opsplt


class RC:
    confinement = True
    transverse_reinf_type = ''
    treat_reinforcement_as_point = True
    _reinforcement = None
    _Ec = None
    _Es = None
    _eps_c = None

    def __init__(self, conc_cross_section, reinforcement, fc, fy, units, dbt=None, s=None, fyt=None):
        self.conc_cross_section = conc_cross_section
        self.reinforcement = reinforcement
        self.fc = fc
        self.fy = fy
        self.units = units
        self.dbt = dbt
        self.s = s
        self.fyt = fyt

    @property
    def Ec(self):
        if self._Ec is not None:
            return self._Ec

        if self.units.lower() == 'us':
            return 57 * sqrt(self.fc * 1000)

        if self.units.lower() == 'si':
            return 4700 * sqrt(fc)

        raise ValueError(f'Ec is not set and default is not impleted for {units = }')

    @Ec.setter
    def Ec(self, x):
        self._Ec = x

    @property
    def Es(self):
        if self._Es is not None:
            return self._Es

        if self.units.lower() == 'us':
            return 29000

        if self.units.lower() == 'si':
            return 2e5

        raise ValueError(f'Es is not set and default is not impleted for {units = }')

    @Es.setter
    def Es(self, x):
        self._Es = x

    @property
    def reinforcement(self):
        return self._reinforcement

    @reinforcement.setter
    def reinforcement(self, x):
        if type(x) == list:
            self._reinforcement = x
        else:
            self._reinforcement = [x]

    def depth(self, axis):
        d = self.conc_cross_section.depth(axis)
        return d

    @property
    def Ag(self):
        return self.conc_cross_section.A

    @property
    def Ac(self):
        return self.Ag - self.Asr

    @property
    def Asr(self):
        a = 0
        for i in self.reinforcement:
            a += i.num_bars * i.Ab
        return a

    def Ig(self, axis):
        return self.conc_cross_section.I(axis)

    def Ic(self, axis):
        return self.Ig(axis) - self.Isr(axis)

    def Isr(self, axis):
        i = 0
        for j in self.reinforcement:
            i += j.I(axis)
        return i

    @property
    def p0(self):
        p0 = 0.85 * self.fc * (self.Ag - self.Asr) + self.fy * self.Asr
        return p0

    @property
    def pnco(self):
        # See Section 22.4.2 of ACI318-19
        if self.transverse_reinf_type.lower() == 'ties':
            pnco = 0.80 * self.p0
        elif self.transverse_reinf_type.lower() == 'spiral' or 'spirals':
            pnco = 0.85 * self.p0
        else:
            raise ValueError("Unknown transverse_reinf_type")
        return pnco

    def phi(self, et):
        f = ACI_phi(self.transverse_reinf_type, et, self.fy / self.Es)
        return f

    def plot_section(self, line_width=2, save_name=None):
        plt.figure()
        self.conc_cross_section.plot_section(line_width)
        for i in range(len(self.reinforcement)):
            self.reinforcement[i].plot_section()
        plt.box(False)
        plt.axis('scaled')

        if save_name:
            plt.savefig(save_name)
        else:
            plt.show()

        plt.clf()
        plt.close()

    def aci_strain_compatibility_object(self):
        id_conc = 1
        id_reinf = 2
        fs = self.fiber_section_object(id_conc, id_reinf)
        fs.print_section_properties()
        scACI = AciStrainCompatibility(fs)
        # add concrete boundaries
        x, y, r = self.conc_cross_section.boundary_points
        for i in range(len(x)):
            scACI.add_concrete_boundary(x[i], y[i], r[i])
        # add steel boundary
        for i in self.reinforcement:
            x, y = i.coordinates

            for j in range(len(x)):
                if self.treat_reinforcement_as_point:
                    scACI.add_steel_boundary(x[j], y[j], 0)
                else:
                    raise ValueError("Not implemented yet")

        scACI.maxCompressiveStrength = -self.pnco
        scACI.add_material(id_conc, 'concrete', self.fc, self.units)
        scACI.add_material(id_reinf, 'steel', self.fy, self.Es)
        return scACI

    def fiber_section_object(self, id_conc, id_reinf, nfx=200, nfy=200):
        fs = FiberSection(nfx, nfy)
        self.conc_cross_section.add_to_fiber_section(fs, id_conc)
        for i in self.reinforcement:
            i.add_to_fiber_section(fs, id_reinf, id_conc)
        return fs

    def section_interaction_2d(self):
        scACI = self.aci_strain_compatibility_object()

        plot_P = []
        plot_Mx = []
        plot_et = []
        scACI.build_data()

        upper = self.conc_cross_section.H * 1.51
        lower = -self.conc_cross_section.H * 2
        r1 = np.array([-1e99])
        r2 = np.arange(lower, lower / 2, 0.2)
        r3 = np.arange(lower / 2, upper, 0.2)
        r4 = np.array([1e99])
        d = np.concatenate((r1, r2, r3, r4))

        for i in d:
            P, Mx, My, et = scACI.compute_point(0, i, 0)
            plot_P.append(-P)
            plot_Mx.append(Mx)
            plot_et.append(et)

        return plot_P, plot_Mx, plot_et

    def plot_interaction_diagram(self, save_name=None):
        plot_P, plot_Mx, plot_et = self.section_interaction_2d()
        phi = []
        for i in plot_et:
            phi.append(self.phi(i))
        plt.plot(plot_Mx, plot_P)
        phi = np.array(phi)
        plot_Mx = np.array(plot_Mx)
        plot_P = np.array(plot_P)
        plt.plot(phi * plot_Mx, phi * plot_P)
        plt.xlabel('M')
        plt.ylabel('P')
        plt.title("Interaction Diagram")

        if save_name:
            plt.savefig(f'{save_name}')
        else:
            plt.show()

        plt.clf()
        plt.close()


def run_example():
    from libdenavit.section import Rectangle, ReinfRect

    # Define concrete cross section
    H = 40
    B = 20
    conc_cross_section = Rectangle(H, B)

    # Define reinforcement
    rhosr = 0.06
    nbB = 2
    nbH = 2
    cover = 0.15 * H
    Ab = H * B * rhosr / (2 * nbB + 2 * nbH - 4)
    reinforcement = ReinfRect(B - 2 * cover, H - 2 * cover, nbB, nbH, Ab)

    # Define materials
    fc = 4
    fy = 60
    units = 'US'

    # Define RC object
    section = RC(conc_cross_section, reinforcement, fc, fy, units)
    section.transverse_reinf_type = 'ties'
    section.plot_section()
    section.plot_interaction_diagram()

if __name__ == '__main__':
    run_example()
