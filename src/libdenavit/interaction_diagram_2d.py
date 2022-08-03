import numpy as np
import math
from libdenavit import find_limit_point_in_list, interpolate_list
from find_intersection_between_two_lines import find_intersection_between_two_lines
import operator

def cart2pol(x, y):
    assert type(x) == type(y), "x and y must be of the same type"
    if type(x) in [int, float]:
        x = [x]
        y = [y]

    assert len(x) == len(y), "x and y must be of the same length"

    rho = []
    phi = []
    for i, _ in enumerate(x):
        rho.append(np.sqrt(x[i]**2 + y[i]**2))
        phi.append(np.arctan2(y[i], x[i]))

    return rho, phi


class InteractionDiagram2D():

    def __init__(self, idx: list, idy: list):
        _, q = cart2pol(idx, idy)
        q = [i%(2*np.pi) for i in q]
        ind = np.argsort(q)
        q = np.sort(q)
        idx = [idx[i] for i in ind]
        idy = [idy[i] for i in ind]

        self.idx = idx
        self.idy = idy
        self.q = q


    def radial_distance(self, angles, degrees=False):
        if degrees:
            angles = np.deg2rad(angles)

        angles = [i%(2*np.pi) for i in angles]
        d = [None] * len(angles)

        for i in range(len(angles)):
            ind = None

            for j, k in enumerate(self.q):
                if k >= angles[i]:
                    ind = j
                    break

            if ind is None:
                d[i] = None

            elif ind == 0:
                if self.q[0] == angles[i]:
                    d[i] = math.hypot(self.idx[0], self.idy[0])
                else:
                    d[i] = 0
            else:
                Ax = self.idx[ind]
                Ay = self.idy[ind]
                Bx = self.idx[ind-1]
                By = self.idy[ind-1]
                Ix, Iy = find_intersection_between_two_lines(Ax, Ay, Bx, By, 0, 0, np.cos(angles[i]), np.sin(angles[i]))
                try:
                    if len(Ix) == 0:
                        raise ValueError('Bad intersection diagram')
                except:
                    d[i] = math.hypot(Ix, Iy)

        return d


    def compare_two(self, test_id, angles):
        d_base = self.radial_distance(angles)
        d_test = test_id.radial_distance(angles)

        errors = [(i-j) / i for i, j in zip(d_base, d_test)]
        return errors


    def check_points(self, pointsX, pointsY):
        d_pts, q_pts = cart2pol(pointsX, pointsY)
        d = self.radial_distance(q_pts)
        max_ratio = np.max(np.divide(d_pts, d))
        return max_ratio


    def findIntersection(self, pathX, pathY):
        d_path, q_path = cart2pol(pathX, pathY)
        d = self.radial_distance(q_path)
        ind, x = find_limit_point_in_list(list(map(operator.sub, d_path, d)), 0)
        if ind == None:
            X = None
            Y = None
        else:
            X = interpolate_list(pathX, ind, x)
            Y = interpolate_list(pathY, ind, x)

        return X, Y, ind, x


    def findXgivenY(self, Y, signX):
        if signX.lower() in ['+', 'positive', 'pos']:
            peakX = 1.1 * np.max(self.idx)
        elif signY.lower() in ['-', 'negative', 'neg']:
            peakX = 1.1 * np.min(self.idx)
        else:
            raise ValueError('signX must be positive or negative')

        npts = 1000
        pathX = np.linspace(0, peakX, npts)
        pathY = Y * np.ones(npts)
        X, _, _, _ = self.findIntersection(pathX, pathY)
        return X


    def findYgivenX(self, X, signY):
        if signY.lower() in ['+', 'positive', 'pos']:
            peakY = 1.1 * max(self.idy)
        elif signY.lower() in ['-', 'negative', 'neg']:
            peakY = 1.1 * min(self.idy)
        else:
            raise ValueError('signY must be positive or negative')

        npts = 1000
        pathX = X * np.ones(npts)
        pathY = np.linspace(0, peakY, npts)
        _, Y, _, _ = self.findIntersection(pathX, pathY)
        return Y


    def plot(self, *args):
        plt.plot(self.idx, self.idy, *args)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    a1 = [0.9, 1.0, 0.8, 0.0]
    b1 = [0.0, 0.4, 0.9, 1.0]
    c1 = InteractionDiagram2D(a1, b1)

    a2 = [0.8, 1.0, 0.8, 0.0]
    b2 = [0.0, 0.4, 1.0, 1.1]
    c2 = InteractionDiagram2D(a2, b2)

    angles = np.linspace(0, np.pi/2, 50)
    distance = c1.radial_distance(angles)

    print("Compare errors= \n",
          c1.compare_two(c2, angles))
    print("Check point (0.9, 0.0): \n",
          c1.check_points(0.9, 0.0))
    print("Find Intersection: \n",
          c1.findIntersection(a2, b2))
    print("Find X given Y: \n",
          c1.findXgivenY(0.9, '+'))
    print("Find Y given X: \n",
          c1.findYgivenX(0.9, '+'))

    plt.plot(distance*np.cos(angles), distance*np.sin(angles), 'bo', label='Angle Points')
    plt.plot(a1, b1, '-r', label='Interaction Diagram 1')
    plt.plot(a2, b2, '-g', label='Interaction Diagram 2')
    plt.legend()

    plt.show()

