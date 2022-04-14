from dataclasses import dataclass
import numpy as np
import math

@dataclass
class FiberSingle:
    # region input arguments
    A: float
    x: float
    y: float
    m: int
    m_neg: int = None  # default to none

    # endregion

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
