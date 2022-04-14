import matplotlib.pyplot as plt
from math import inf
import numpy as np
import pandas as pd


class FiberSection:
    def __init__(self, nfx=1, nfy=1):
        self.nfx = nfx
        self.nfy = nfy
        self.default_axes_origin = "gross_section_centroid"
        self._fibers_object_list = []

    def add_fibers(self, *new_fibers):
        for i in new_fibers:
            self._fibers_object_list.append(i)

    def get_bounds(self):
        xmin = inf
        xmax = -inf
        ymin = inf
        ymax = -inf

        for fibers in self._fibers_object_list:
            ixmin, ixmax, iymin, iymax = fibers.get_bounds()
            if ixmin < xmin:
                xmin = ixmin

            if ixmax > xmax:
                xmax = ixmax

            if iymin < ymin:
                ymin = iymin

            if iymax > ymax:
                ymax = iymax

        return xmin, xmax, ymin, ymax

    def get_fiber_data(self):
        # Initialize arrays
        A = np.array([])
        x = np.array([])
        y = np.array([])
        m = np.array([], dtype=int)

        # Compute maximum fiber size
        xmin, xmax, ymin, ymax = self.get_bounds()
        sfx = (xmax - xmin) / self.nfx
        sfy = (ymax - ymin) / self.nfy

        # Loop through fibers dictionary and add to total list.
        for fibers in self._fibers_object_list:
            iA, ix, iy, im = fibers.get_fiber_data(sfx, sfy)
            A = np.append(A, iA)
            x = np.append(x, ix)
            y = np.append(y, iy)
            m = np.append(m, im)

        return A, x, y, m

    def plot_fibers(self, scale=1):
        # Does not plot fibers with negative area
        A, x, y, m = self.get_fiber_data()
        color_map = {1: 'bo', 2: 'ro', 3: 'ko', 4: 'yo', 5: 'co', 6: 'go'}
        for i, j in enumerate(m):
            if A[i] > 0:
                plt.plot(x[i], y[i], color_map.get(j, 'ko'), markersize=A[i] * 1 * scale)

        plt.axis('equal')
        plt.show()

    def unique_mat_ids(self):
        _, _, _, m = self.get_fiber_data()
        mat = np.unique(m)
        return mat

    def print_section_properties(self):
        '''
        This function spits out a table like this.
        
        Mat   num_fibers  A      Ix      Iy
         1     ###       ###    ####    ####
         2     ###       ###    ####    ####
        Total  ###       ###    ####    ####
        '''

        A, x, y, m = self.get_fiber_data()
        m_unique = np.unique(m)

        Ix = A * y ** 2
        Iy = A * x ** 2
        num_fibers = {i: 0 for i in m_unique}
        A_print = {i: 0 for i in m_unique}
        Ix_print = {i: 0 for i in m_unique}
        Iy_print = {i: 0 for i in m_unique}

        for i, j in enumerate(m):
            num_fibers[j] += 1
            A_print[j] += A[i]
            Ix_print[j] += Ix[i]
            Iy_print[j] += Iy[i]

        m = np.append(m, "Total")
        num_fibers[-1] = sum(num_fibers.values())
        A_print[-1] = sum(A_print.values())
        Ix_print[-1] = sum(Ix_print.values())
        Iy_print[-1] = sum(Iy_print.values())

        print_list = pd.DataFrame({"Material": np.unique(m).tolist(), "num_fibers": num_fibers.values(),\
                                   "Area": A_print.values(), "Ix": Ix_print.values(), "Iy": Iy_print.values()})

        print(print_list.to_string(index=False))


def run_example():
    from libdenavit.section import FiberSingle, FiberQuadPatch

    a1 = FiberSingle(2.0, 1, 2, 2, 1)  # A x y m m_neg
    a2 = FiberSingle(2.0, 1, 8, 2, 1)
    a3 = FiberSingle(2.0, 4, 8, 2, 1)
    a4 = FiberSingle(2.0, 4, 2, 2, 1)
    b1 = FiberQuadPatch(0, 0, 0, 2, 5, 2, 5, 0, 1)  # xI yI ... m
    b2 = FiberQuadPatch(0, 2, 0, 8, 1, 8, 1, 2, 1)
    b3 = FiberQuadPatch(0, 8, 0, 10, 5, 10, 5, 8, 1)
    b4 = FiberQuadPatch(4, 2, 4, 8, 5, 8, 5, 2, 1)
    b5 = FiberQuadPatch(1, 2, 1, 8, 4, 8, 4, 2, 1)
    c = FiberSection(10, 10)
    c.add_fibers(a1, a2, a3, a4, b1, b2, b3, b4, b5)
    c.plot_fibers(scale=5)
    c.print_section_properties()


if __name__ == '__main__':
    run_example()
