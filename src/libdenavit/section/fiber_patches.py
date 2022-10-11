from dataclasses import dataclass
import numpy as np
import math

@dataclass
class FiberSingle:
    """
    FiberSingle is a class that represents a single fiber.

    Parameters
    ----------

    A : float
     the area of the fiber
    x : float
     the x-coordinate of the fiber
    y : float
     the y-coordinate of the fiber
    m : float
     the material of the fiber
    m_neg : float
     the material of the negative fiber (optional)

    """
    A: float
    x: float
    y: float
    m: int
    m_neg: int = None  # default to none


    def get_bounds(self, ):
        return self.x, self.x, self.y, self.y

    def get_fiber_data(self, sfx, sfy):

        A = [self.A]
        x = [self.x]
        y = [self.y]
        m = [self.m]

        if self.m_neg is not None:

            A.append(-self.A)
            x.append(self.x)
            y.append(self.y)
            m.append(self.m_neg)

        return np.array(A), np.array(x), np.array(y), np.array(m)

@dataclass
class FiberQuadPatch:

    """
    # Creates a set of fibers to describe a quadrilateral patch
    #
    # x and y parameters are the coordiantes of the four corners of the patch as labeled below:
    #
    #   y |
    #     |    J  o-----o K
    #     |      /     /
    #     |     /     /
    #     |  I o-----o L
    #     |______________
    #                   x
    """


    # region input arguments
    xI: float
    yI: float
    xJ: float
    yJ: float
    xK: float
    yK: float
    xL: float
    yL: float
    mat_id: int

    # endregion

    def get_bounds(self, ):
        xmin = min(self.xI, self.xJ, self.xK, self.xL)
        xmax = max(self.xI, self.xJ, self.xK, self.xL)
        ymin = min(self.yI, self.yJ, self.yK, self.yL)
        ymax = max(self.yI, self.yJ, self.yK, self.yL)
        return xmin, xmax, ymin, ymax

    def get_fiber_data(self, sfx, sfy):

        A = []  # Make these numpy arrays
        x = []
        y = []
        m = []

        nfx = math.ceil(max(1, abs(self.xJ - self.xK) / sfx))  # number of fibers in the x
        nfy = math.ceil(max(1, abs(self.yI - self.yJ) / sfy))  # number of fibers in the y

        for i in range(nfx):
            for j in range(nfy):
                B = abs(self.xJ - self.xK)
                H = abs(self.yK - self.yL)

                A.append(B / nfx * H / nfy)  # Area
                x.append(self.xI + (i * 2 + 1) * B / (2 * nfx))  # Centroid z
                y.append(self.yI + (j * 2 + 1) * H / (2 * nfy))  # Centroid y
                m.append(self.mat_id)  # Material

        return np.array(A), np.array(x), np.array(y), np.array(m)


@dataclass
class FiberCirclePatch:
    """
    FiberCirclePatch is a class that represents a Circle Patch.

    Parameters
    ----------
    xc : float
     the x-coordinate of the center of the circle
    yc : float
     the y-coordinate of the center of the circle
    ri : float
     the radius of the inner circle
    ro : float
     the radius of the outer circle
    mat_id : int
     the material of the fiber
    is_neg : float
     the material of the negative fiber (optional)
    a1 : float
     the starting angle of the circle
    a2 : float
     the ending angle of the circle
    """
    xc: float
    yc: float
    ri: float
    ro: float
    mat_id: int
    is_neg: bool = None  # default to none
    a1: float = 0
    a2: float = 2 * np.pi
    
    def get_bounds(self):
        xmin = self.xc - self.ro
        xmax = self.xc + self.ro
        ymin = self.yc - self.ro
        ymax = self.yc + self.ro
        return xmin, xmax, ymin, ymax

    def get_fiber_data(self, sfx, sfy):
        A = []
        x = []
        y = []

        # Compute fiber data
        sf = min(sfx,sfy)
        nf_rad = math.ceil((self.ro-self.ri) / sf)
        x_rad = np.linspace(self.ri, self.ro, nf_rad + 1)

        for i in range(nf_rad):
            ri_ith = x_rad[i]
            ro_ith = x_rad[i + 1]
            nf_arc = max(1, math.ceil(ro_ith * (self.a2 - self.a1) / sf))
            d_ang = (self.a2 - self.a1) / nf_arc
            x_ang = np.linspace(self.a1, self.a2, nf_arc + 1)
            Ai = []
            xi = []
            yi = []
            for j in range(nf_arc):
                Ai_1 = (d_ang / 2) * ro_ith ** 2
                Ai_2 = -(d_ang / 2) * ri_ith ** 2
                Ai_th = Ai_1 + Ai_2
                yi_1 = 4 * ro_ith * np.sin(d_ang / 2) / 3 / d_ang
                yi_2 = 4 * ri_ith * np.sin(d_ang / 2) / 3 / d_ang
                yi_th = ((Ai_1 * yi_1 + Ai_2 * yi_2) / Ai_th)
                ang_i = (x_ang[j + 1] + x_ang[j]) / 2

                # Assign fiber data
                Ai.append(Ai_th)
                xi.append(self.xc + yi_th * np.cos(ang_i))
                yi.append(self.yc + yi_th * np.sin(ang_i))

            A.extend(Ai)
            x.extend(xi)
            y.extend(yi)

        A = np.array(A)
        x = np.array(x)
        y = np.array(y)
        m = np.array([self.mat_id] * len(A))
        if self.is_neg:
            A = -A
        return A, x, y, m
